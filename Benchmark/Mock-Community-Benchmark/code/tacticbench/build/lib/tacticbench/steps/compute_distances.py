from __future__ import annotations

from tempfile import TemporaryDirectory
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

from ..log import get_logger
from ..utils.paths import distances_dir, expected_dir, method_dir, qc_dir
from ..utils.phylo import UniFracTree
from ..utils.blast import BlastTools, read_hits, run_blast
from ..utils.io import canonical_species_name


CLUSTERING = {"TIC", "TAC", "OTU"}
METHOD_KEY_ALIASES = {
    "ASV": ["asv", "dada2"],
    "zOTU": ["zotu", "unoise3"],
    "OTU": ["otu", "uparse"],
    "TIC": ["tic"],
    "TAC": ["tac"],
}


def _resolve_target_level(cfg, method: str) -> str:
    default_level = "species" if method in CLUSTERING else "tip"
    raw = getattr(cfg.mapping, "target_level_by_method", {})
    if not isinstance(raw, dict) or not raw:
        return default_level
    for k in METHOD_KEY_ALIASES.get(method, []):
        v = str(raw.get(k, "")).strip().lower()
        if v in {"tip", "species"}:
            return v
    vdef = str(raw.get("default", "")).strip().lower()
    if vdef in {"tip", "species"}:
        return vdef
    return default_level


def _jaccard(a: np.ndarray, b: np.ndarray, extra_inferred_only: int = 0) -> float:
    a1 = a > 0
    b1 = b > 0
    inter = np.logical_and(a1, b1).sum()
    union = np.logical_or(a1, b1).sum()
    if extra_inferred_only > 0:
        union += int(extra_inferred_only)
    return 0.0 if union == 0 else float(1.0 - (inter / union))


def _weighted_jaccard(a: np.ndarray, b: np.ndarray) -> float:
    den = np.maximum(a, b).sum()
    return 0.0 if den == 0 else float(1.0 - (np.minimum(a, b).sum() / den))


def _braycurtis(a: np.ndarray, b: np.ndarray) -> float:
    denom = a.sum() + b.sum()
    return 0.0 if denom == 0 else float(np.abs(a - b).sum() / denom)


def _aitchison(a: np.ndarray, b: np.ndarray, pseudocount: float = 1e-6) -> float:
    """
    Aitchison distance between two compositions.

    Steps:
      1) add pseudocount to handle zeros
      2) closure to unit sum
      3) CLR transform
      4) Euclidean distance in CLR space
    """
    pc = float(pseudocount)
    if not np.isfinite(pc) or pc <= 0.0:
        pc = 1e-6

    xa = np.asarray(a, dtype=float)
    xb = np.asarray(b, dtype=float)
    if xa.size == 0 or xb.size == 0:
        return float("nan")

    xa = np.clip(xa, 0.0, None) + pc
    xb = np.clip(xb, 0.0, None) + pc

    sa = float(xa.sum())
    sb = float(xb.sum())
    if sa <= 0.0 or sb <= 0.0:
        return float("nan")

    pa = xa / sa
    pb = xb / sb

    la = np.log(pa)
    lb = np.log(pb)
    clra = la - float(la.mean())
    clrb = lb - float(lb.mean())
    return float(np.linalg.norm(clra - clrb))


def _safe_unifrac(tree: Optional[UniFracTree], vec_a: np.ndarray, vec_b: np.ndarray) -> float:
    if tree is None:
        return float("nan")
    return tree.unifrac_unweighted(vec_a, vec_b)


def _safe_gunifrac(tree: Optional[UniFracTree], vec_a: np.ndarray, vec_b: np.ndarray, alpha: float) -> float:
    if tree is None:
        return float("nan")
    return tree.gunifrac(vec_a, vec_b, alpha=alpha)


def _load_tree(path: Path) -> Optional[UniFracTree]:
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return UniFracTree(path)
    except Exception:
        return None


def _load_unmapped_fraction_map(qc_path: Path) -> Dict[str, float]:
    """
    Load sample -> unmapped fraction from mapping_qc.tsv.
    Prefers mapping_rate when available; otherwise uses mapped/unmapped mass.
    """
    if not qc_path.exists() or qc_path.stat().st_size == 0:
        return {}
    try:
        qcdf = pd.read_csv(qc_path, sep="\t")
    except Exception:
        return {}

    out: Dict[str, float] = {}
    for _, row in qcdf.iterrows():
        sid = str(row.get("sample_id", ""))
        if not sid:
            continue

        val = np.nan
        mr = row.get("mapping_rate", np.nan)
        if pd.notna(mr):
            val = 1.0 - float(mr)
        else:
            mm = row.get("mapped_mass", np.nan)
            um = row.get("unmapped_mass", np.nan)
            if pd.notna(mm) and pd.notna(um) and (float(mm) + float(um)) > 0:
                val = float(um) / (float(mm) + float(um))

        if pd.isna(val):
            continue
        out[sid] = float(np.clip(val, 0.0, 1.0))
    return out


def _load_unmapped_mass_map(qc_path: Path) -> Dict[str, float]:
    """
    Load sample -> unmapped mass from mapping_qc.tsv.
    """
    if not qc_path.exists() or qc_path.stat().st_size == 0:
        return {}
    try:
        qcdf = pd.read_csv(qc_path, sep="\t")
    except Exception:
        return {}

    out: Dict[str, float] = {}
    for _, row in qcdf.iterrows():
        sid = str(row.get("sample_id", ""))
        if not sid:
            continue
        um = row.get("unmapped_mass", np.nan)
        if pd.isna(um):
            continue
        out[sid] = float(max(0.0, float(um)))
    return out


def _load_unmapped_feature_count_map(unmapped_path: Path, table_std_path: Path) -> Dict[str, int]:
    """
    Load sample -> count of unmapped inferred features with positive abundance.

    This is used by Jaccard so that each unmapped inferred feature contributes
    one presence/absence element instead of all unmapped mass being collapsed
    into a single pseudo-feature.
    """
    if (not unmapped_path.exists()) or (not table_std_path.exists()):
        return {}
    try:
        unmapped_df = pd.read_csv(unmapped_path, sep="\t")
        tab = pd.read_csv(table_std_path, sep="\t", index_col=0)
    except Exception:
        return {}

    if "inferred_feature" not in unmapped_df.columns or tab.empty:
        return {}

    unmapped = [str(x) for x in unmapped_df["inferred_feature"].dropna().astype(str).tolist()]
    if not unmapped:
        return {}

    tab = tab.copy()
    tab.index = tab.index.astype(str)
    common = tab.index.intersection(pd.Index(unmapped))
    if len(common) == 0:
        return {}

    sub = tab.loc[common, :].astype(float)
    out: Dict[str, int] = {}
    for s in sub.columns.astype(str):
        abund = sub[s].to_numpy(dtype=float)
        out[s] = int(np.count_nonzero(abund > 0.0))
    return out


def _load_mapping_mismatch_map(mapping_path: Path, table_std_path: Path) -> Dict[str, float]:
    """
    Load sample -> abundance-weighted mapping mismatch from mapping.tsv and table_std_tss1000.tsv.
    mismatch(feature) = 1 - pident/100 for mapped inferred features.
    """
    if not mapping_path.exists() or not table_std_path.exists():
        return {}
    try:
        mp = pd.read_csv(mapping_path, sep="\t")
        tab = pd.read_csv(table_std_path, sep="\t", index_col=0)
    except Exception:
        return {}

    if "inferred_feature" not in mp.columns or "pident" not in mp.columns:
        return {}
    if tab.empty:
        return {}

    feat_to_pid = (
        mp[["inferred_feature", "pident"]]
        .dropna(subset=["inferred_feature", "pident"])
        .copy()
    )
    # soft mapping mode may contain multiple rows per inferred_feature with weights.
    if "weight" in mp.columns:
        feat_to_pid["weight"] = pd.to_numeric(mp["weight"], errors="coerce").fillna(0.0).clip(lower=0.0)
        # compute weighted mean pident per inferred feature; fallback to simple mean if total weight <= 0
        rows = []
        for feat, g in feat_to_pid.groupby("inferred_feature", sort=False):
            p = pd.to_numeric(g["pident"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
            w = g["weight"].to_numpy(dtype=float)
            ws = float(w.sum())
            if ws > 0:
                pid = float((p * w).sum() / ws)
            else:
                pid = float(p.mean()) if len(p) > 0 else float("nan")
            rows.append({"inferred_feature": str(feat), "pident": pid})
        feat_to_pid = pd.DataFrame(rows)
    else:
        feat_to_pid = feat_to_pid.drop_duplicates(subset=["inferred_feature"], keep="first")
    if feat_to_pid.empty:
        return {}

    tab = tab.copy()
    tab.index = tab.index.astype(str)
    pid = pd.Series(feat_to_pid["pident"].astype(float).values, index=feat_to_pid["inferred_feature"].astype(str))

    common = tab.index.intersection(pid.index)
    if len(common) == 0:
        return {}

    mismatch = (1.0 - (pid.loc[common] / 100.0)).clip(lower=0.0, upper=1.0)
    sub = tab.loc[common, :].astype(float)

    out: Dict[str, float] = {}
    for s in sub.columns.astype(str):
        abund = sub[s].to_numpy(dtype=float)
        total = float(abund.sum())
        if total <= 0:
            out[s] = 0.0
            continue
        m = float((abund * mismatch.to_numpy(dtype=float)).sum() / total)
        out[s] = float(np.clip(m, 0.0, 1.0))
    return out


def _apply_penalties(
    distance: float,
    unmapped_fraction: float,
    mapping_mismatch: float,
    coverage_deficit: float,
    penalize_unmapped: bool,
    unmapped_penalty_lambda: float,
    penalize_mapping_mismatch: bool,
    mapping_mismatch_lambda: float,
    penalize_coverage_deficit: bool,
    coverage_deficit_lambda: float,
    penalty_transform: str,
    penalty_exp_k: float,
    bounded_0_1: bool = True,
) -> float:
    def _transform(x: float) -> float:
        xv = float(np.clip(x, 0.0, 1.0))
        if penalty_transform == "exp":
            k = max(0.0, float(penalty_exp_k))
            return float(1.0 - np.exp(-k * xv))
        return xv

    if np.isnan(distance):
        return float(distance)
    penalty = 0.0
    if penalize_unmapped:
        penalty += float(unmapped_penalty_lambda) * _transform(unmapped_fraction)
    if penalize_mapping_mismatch:
        penalty += float(mapping_mismatch_lambda) * _transform(mapping_mismatch)
    if penalize_coverage_deficit:
        penalty += float(coverage_deficit_lambda) * _transform(coverage_deficit)
    penalized = float(distance) + penalty
    if bounded_0_1:
        return float(np.clip(penalized, 0.0, 1.0))
    return penalized


def _coverage_deficit(a: np.ndarray, b: np.ndarray) -> float:
    """
    Proxy for over-merging: fraction of expected-present targets not covered by inferred.
    """
    exp_present = a > 0
    den = int(exp_present.sum())
    if den == 0:
        return 0.0
    covered = int(np.logical_and(exp_present, b > 0).sum())
    return float(np.clip(1.0 - (covered / den), 0.0, 1.0))


def _vectorize_on_tree(
    table: pd.DataFrame,
    sample: str,
    tree: UniFracTree,
) -> np.ndarray:
    """
    Build a vector aligned to tree.tips from a feature x sample table.
    Missing features get 0.
    """
    idx = {name: i for i, name in enumerate(tree.tips)}
    v = np.zeros(len(tree.tips), dtype=float)
    if sample not in table.columns:
        return v
    s = table[sample]
    # only iterate over non-zero to be safe/fast
    nz = s[s != 0]
    for feat, val in nz.items():
        i = idx.get(str(feat))
        if i is not None:
            v[i] = float(val)
    return v


def _write_fasta_subset(src: Path, wanted_ids: set[str], dest: Path) -> int:
    """
    Write a FASTA containing only records whose record.id is in wanted_ids.
    Returns the number of records written.
    """
    if not src.exists() or not wanted_ids:
        return 0

    remaining = set(str(x) for x in wanted_ids)
    written = 0
    with dest.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(src), "fasta"):
            rid = str(rec.id)
            if rid not in remaining:
                continue
            SeqIO.write([rec], out, "fasta")
            written += 1
            remaining.discard(rid)
            if not remaining:
                break
    return written


def _build_augmented_phylo(
    base_tree: Optional[UniFracTree],
    mapping_path: Path,
    blast_hits_path: Path,
    std_table_path: Path,
    query_fasta_path: Optional[Path],
    blast_db_prefix: Optional[Path],
    blast_tools: Optional[BlastTools],
    tip_to_species_path: Optional[Path],
    target_level: str,
    beta: float,
    eps: float,
    include_unmapped: bool,
) -> Tuple[Optional[UniFracTree], Dict[str, int]]:
    """
    Build augmented tree and inferred_feature -> synthetic tip index map for phylo distances.
    """
    if base_tree is None:
        return None, {}
    if (not mapping_path.exists()) or (not std_table_path.exists()):
        return base_tree, {}
    try:
        mp = pd.read_csv(mapping_path, sep="\t")
        tab = pd.read_csv(std_table_path, sep="\t", index_col=0)
    except Exception:
        return base_tree, {}

    if "inferred_feature" not in mp.columns or "pident" not in mp.columns:
        return base_tree, {}

    species_target = str(target_level).strip().lower() == "species"
    if species_target:
        target_col = "target_species"
        if target_col not in mp.columns:
            return base_tree, {}
    else:
        # Backward compatibility for older mapping.tsv files that stored only best_tip.
        if "target_tip" in mp.columns:
            target_col = "target_tip"
        elif "best_tip" in mp.columns:
            target_col = "best_tip"
        else:
            return base_tree, {}

    tip_to_species: Dict[str, str] = {}
    if species_target and tip_to_species_path is not None and tip_to_species_path.exists():
        try:
            t2s = pd.read_csv(tip_to_species_path, sep="\t")
            if {"tip", "species"}.issubset(set(t2s.columns)):
                tip_to_species = dict(zip(t2s["tip"].astype(str), t2s["species"].astype(str)))
        except Exception:
            tip_to_species = {}

    def _species_fallback(tip: str) -> str:
        return canonical_species_name(str(tip))

    def _target_from_sseqid(sseqid: str) -> str:
        if species_target:
            s = str(sseqid)
            return canonical_species_name(str(tip_to_species.get(s, _species_fallback(s))))
        return str(sseqid)

    attached: List[Tuple[str, str, float]] = []
    feat_to_syn_name: Dict[str, str] = {}
    seen_syn = set()

    for _, row in mp.iterrows():
        feat = str(row.get("inferred_feature", "")).strip()
        target = str(row.get(target_col, "")).strip()
        if not feat or not target:
            continue
        pid = row.get("pident", np.nan)
        if pd.isna(pid):
            continue
        mismatch = float(np.clip(1.0 - (float(pid) / 100.0), 0.0, 1.0))
        bl = float(beta * mismatch + eps)
        syn = f"__INF__{feat}"
        if syn in seen_syn:
            continue
        seen_syn.add(syn)
        feat_to_syn_name[feat] = syn
        attached.append((target, syn, bl))

    unattached: List[Tuple[str, float]] = []
    if include_unmapped and not tab.empty:
        mapped = set(feat_to_syn_name.keys())
        unmapped_feats = [f for f in tab.index.astype(str) if f not in mapped]
        placed = set()

        if blast_hits_path.exists() and unmapped_feats:
            try:
                hits = read_hits(blast_hits_path)
            except Exception:
                hits = pd.DataFrame()

            if not hits.empty:
                hs = hits[hits["qseqid"].astype(str).isin(unmapped_feats)].copy()
                if not hs.empty:
                    hs = hs.sort_values(
                        ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                        ascending=[True, False, False, False, True, False],
                    )
                    best = hs.groupby("qseqid", sort=False).first().reset_index()
                    for _, row in best.iterrows():
                        feat = str(row["qseqid"]).strip()
                        if not feat:
                            continue
                        target = _target_from_sseqid(str(row["sseqid"]).strip())
                        if not target:
                            continue
                        pid = row.get("pident", np.nan)
                        if pd.isna(pid):
                            continue
                        mismatch = float(np.clip(1.0 - (float(pid) / 100.0), 0.0, 1.0))
                        bl = float(beta * mismatch + eps)
                        syn = f"__UNMAPPED__{feat}"
                        if syn in seen_syn:
                            continue
                        seen_syn.add(syn)
                        feat_to_syn_name[feat] = syn
                        attached.append((target, syn, bl))
                        placed.add(feat)

        fallback_feats = [feat for feat in unmapped_feats if feat not in placed]
        if (
            fallback_feats
            and query_fasta_path is not None
            and query_fasta_path.exists()
            and blast_db_prefix is not None
            and blast_tools is not None
        ):
            try:
                with TemporaryDirectory(prefix="tacticbench_augmented_mismatch_") as tmpdir:
                    tmpdir_p = Path(tmpdir)
                    subset_fasta = tmpdir_p / "unmapped_subset.fasta"
                    n_written = _write_fasta_subset(query_fasta_path, set(fallback_feats), subset_fasta)
                    if n_written > 0:
                        fallback_hits_path = tmpdir_p / "fallback_hits.tsv"
                        fallback_tools = BlastTools(
                            makeblastdb=blast_tools.makeblastdb,
                            blastn=blast_tools.blastn,
                        )
                        ok = run_blast(
                            tools=fallback_tools,
                            query_fasta=subset_fasta,
                            db_prefix=blast_db_prefix,
                            out_tsv=fallback_hits_path,
                            min_ident=0.0,
                        )
                        if ok and fallback_hits_path.exists() and fallback_hits_path.stat().st_size > 0:
                            try:
                                fallback_hits = read_hits(fallback_hits_path)
                            except Exception:
                                fallback_hits = pd.DataFrame()
                            if not fallback_hits.empty:
                                fh = fallback_hits[fallback_hits["qseqid"].astype(str).isin(fallback_feats)].copy()
                                if not fh.empty:
                                    fh = fh.sort_values(
                                        ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                                        ascending=[True, False, False, False, True, False],
                                    )
                                    best_fb = fh.groupby("qseqid", sort=False).first().reset_index()
                                    for _, row in best_fb.iterrows():
                                        feat = str(row["qseqid"]).strip()
                                        if not feat:
                                            continue
                                        target = _target_from_sseqid(str(row["sseqid"]).strip())
                                        if not target:
                                            continue
                                        pid = row.get("pident", np.nan)
                                        if pd.isna(pid):
                                            continue
                                        mismatch = float(np.clip(1.0 - (float(pid) / 100.0), 0.0, 1.0))
                                        bl = float(beta * mismatch + eps)
                                        syn = f"__UNMAPPED__{feat}"
                                        if syn in seen_syn:
                                            continue
                                        seen_syn.add(syn)
                                        feat_to_syn_name[feat] = syn
                                        attached.append((target, syn, bl))
                                        placed.add(feat)
            except Exception:
                # Best-effort fallback only.
                pass

        # fallback for unmapped features with no usable BLAST hit:
        for feat in unmapped_feats:
            if feat in placed:
                continue
            syn = f"__UNMAPPED__{feat}"
            if syn in seen_syn:
                continue
            seen_syn.add(syn)
            feat_to_syn_name[feat] = syn
            unattached.append((syn, float(beta + eps)))

    aug = base_tree.with_synthetic_tips(attached=attached, unattached=unattached)
    feat_to_syn_idx = {f: aug.tip_index[syn] for f, syn in feat_to_syn_name.items() if syn in aug.tip_index}
    return aug, feat_to_syn_idx


def _vectorize_augmented_inferred(
    table_std: pd.DataFrame,
    sample: str,
    n_tips: int,
    feat_to_syn_idx: Dict[str, int],
) -> np.ndarray:
    v = np.zeros(n_tips, dtype=float)
    if sample not in table_std.columns:
        return v
    s = table_std[sample]
    nz = s[s != 0]
    for feat, val in nz.items():
        i = feat_to_syn_idx.get(str(feat))
        if i is not None:
            v[i] = float(val)
    return v


def run(cfg) -> int:
    """
    Compute distances between expected composition and each method per sample.

    Inputs per community (produced by previous steps):
      Expected:
        - expected/expected_tip_abundance_tss1000.tsv
        - expected/expected_species_abundance_tss1000.tsv
        - expected/expected_tree_full.nwk (for ASV/zOTU UniFrac/gUniFrac)
        - expected/centroid_species_tree.nwk (for TIC/TAC/OTU UniFrac/gUniFrac)

      Methods (resolved by mapping.target_level_by_method):
        - tip target: methods/<m>/mapped_tip_table_tss1000.tsv
        - species target: methods/<m>/mapped_species_table_tss1000.tsv

    Outputs:
      Per-community:
        - communities/<cid>/distances/distances_long.tsv
        - communities/<cid>/distances/distances_wide.tsv

      Global:
        - outdir/distances_all_long.tsv
        - outdir/distances_all_wide.tsv
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    manifest = pd.read_csv(outdir / "manifest.tsv", sep="\t")
    blast_tools = BlastTools(makeblastdb=cfg.tree.makeblastdb_exe, blastn=cfg.tree.blastn_exe)

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    # metrics registry
    metric_fns: Dict[str, Callable[[np.ndarray, np.ndarray], float]] = {
        "jaccard": _jaccard,
        "weighted_jaccard": _weighted_jaccard,
        "braycurtis": _braycurtis,
    }
    # phylo metrics handled separately
    requested_methods: List[str] = list(getattr(cfg.distances, "methods", ["zOTU", "ASV", "TIC", "TAC", "OTU"]))
    requested_metrics: List[str] = list(
        getattr(
            cfg.distances,
            "metrics",
            [
                "jaccard",
                "weighted_jaccard",
                "braycurtis",
                "unifrac",
                "gunifrac_a0.0",
                "gunifrac_a0.5",
                "gunifrac_a1.0",
            ],
        )
    )
    aitchison_pseudocount = float(getattr(cfg.distances, "aitchison_pseudocount", 1e-6))
    penalize_unmapped = bool(getattr(cfg.distances, "penalize_unmapped", False))
    unmapped_penalty_lambda = float(getattr(cfg.distances, "unmapped_penalty_lambda", 1.0))
    penalize_mapping_mismatch = bool(getattr(cfg.distances, "penalize_mapping_mismatch", False))
    mapping_mismatch_lambda = float(getattr(cfg.distances, "mapping_mismatch_lambda", 1.0))
    penalize_coverage_deficit = bool(getattr(cfg.distances, "penalize_coverage_deficit", False))
    coverage_deficit_lambda = float(getattr(cfg.distances, "coverage_deficit_lambda", 1.0))
    penalty_transform = str(getattr(cfg.distances, "penalty_transform", "linear")).strip().lower()
    penalty_transform = "exp" if penalty_transform == "exp" else "linear"
    penalty_exp_k = float(getattr(cfg.distances, "penalty_exp_k", 3.0))
    phylo_mapping_mode = str(getattr(cfg.distances, "phylo_mapping_mode", "mapped_projection")).strip().lower()
    if phylo_mapping_mode == "classic":
        phylo_mapping_mode = "mapped_projection"  # legacy alias
    phylo_mapping_mode = (
        "augmented_mismatch" if phylo_mapping_mode == "augmented_mismatch" else "mapped_projection"
    )
    augmented_beta = float(getattr(cfg.distances, "augmented_mismatch_beta", 1.0))
    augmented_eps = float(getattr(cfg.distances, "augmented_mismatch_eps", 1e-6))
    augmented_include_unmapped = bool(getattr(cfg.distances, "augmented_mismatch_include_unmapped", False))

    # Mode-specific semantics for unmapped handling:
    # - mapped_projection: penalize_unmapped => additive post-hoc penalty by unmapped fraction
    # - augmented_mismatch: augmented_mismatch_include_unmapped => include unmapped synthetic tips
    use_additive_unmapped_penalty = penalize_unmapped and (phylo_mapping_mode == "mapped_projection")
    use_additive_mapping_mismatch_penalty = (
        penalize_mapping_mismatch and (phylo_mapping_mode == "mapped_projection")
    )
    use_additive_coverage_deficit_penalty = (
        penalize_coverage_deficit and (phylo_mapping_mode == "mapped_projection")
    )
    use_augmented_unmapped_tips = augmented_include_unmapped and (phylo_mapping_mode == "augmented_mismatch")

    all_rows: List[dict] = []

    for _, r in manifest.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue

        ed = expected_dir(outdir, cid)
        exp_tip_p = ed / "expected_tip_abundance_tss1000.tsv"
        exp_sp_p = ed / "expected_species_abundance_tss1000.tsv"
        if not exp_tip_p.exists() or not exp_sp_p.exists():
            continue

        exp_tip = pd.read_csv(exp_tip_p, sep="\t", index_col=0)
        exp_sp = pd.read_csv(exp_sp_p, sep="\t", index_col=0)

        tree_full = _load_tree(ed / "expected_tree_full.nwk")
        tree_cent = _load_tree(ed / "centroid_species_tree.nwk")

        comm_rows: List[dict] = []

        for method in requested_methods:
            md = method_dir(outdir, cid, method)
            target_level = _resolve_target_level(cfg, method)
            species_target = target_level == "species"
            unmapped_fraction_by_sample = _load_unmapped_fraction_map(md / "mapping_qc.tsv")
            unmapped_mass_by_sample = _load_unmapped_mass_map(md / "mapping_qc.tsv")
            unmapped_feature_count_by_sample = _load_unmapped_feature_count_map(
                md / "unmapped_features.tsv",
                md / "table_std_tss1000.tsv",
            )
            mapping_mismatch_by_sample = _load_mapping_mismatch_map(md / "mapping.tsv", md / "table_std_tss1000.tsv")
            std_table = None
            phylo_tree = None
            feat_to_syn_idx: Dict[str, int] = {}

            # choose expected space + inferred table + tree
            if species_target:
                expected = exp_sp
                inferred_p = md / "mapped_species_table_tss1000.tsv"
                if not inferred_p.exists():
                    continue
                inferred = pd.read_csv(inferred_p, sep="\t", index_col=0)

                # align inferred to expected index (keep all expected features!)
                inferred = inferred.reindex(expected.index).fillna(0.0)

                # for phylo distances, use centroid tree (species tips)
                tree = tree_cent
            else:
                expected = exp_tip
                inferred_p = md / "mapped_tip_table_tss1000.tsv"
                if not inferred_p.exists():
                    continue
                inferred = pd.read_csv(inferred_p, sep="\t", index_col=0)

                inferred = inferred.reindex(expected.index).fillna(0.0)

                # for phylo distances, use full tree (tip tips)
                tree = tree_full

            if phylo_mapping_mode == "augmented_mismatch":
                std_p = md / "table_std_tss1000.tsv"
                if std_p.exists():
                    try:
                        std_table = pd.read_csv(std_p, sep="\t", index_col=0)
                    except Exception:
                        std_table = None
                phylo_tree, feat_to_syn_idx = _build_augmented_phylo(
                    base_tree=tree,
                    mapping_path=md / "mapping.tsv",
                    blast_hits_path=md / "blast_hits.tsv",
                    std_table_path=std_p,
                    query_fasta_path=md / "seqs.fasta",
                    blast_db_prefix=ed / "blastdb_expected",
                    blast_tools=blast_tools,
                    tip_to_species_path=(ed / "tip_to_species.tsv") if species_target else None,
                    target_level=target_level,
                    beta=augmented_beta,
                    eps=augmented_eps,
                    include_unmapped=use_augmented_unmapped_tips,
                )
            else:
                phylo_tree = tree

            # samples to compute (should match after standardization)
            samples = [s for s in expected.columns if s in inferred.columns]
            if not samples:
                continue

            for s in samples:
                a = expected[s].to_numpy(dtype=float)
                b = inferred[s].to_numpy(dtype=float)
                unmapped_fraction = float(unmapped_fraction_by_sample.get(s, 0.0))
                unmapped_mass = float(unmapped_mass_by_sample.get(s, 0.0))
                mapping_mismatch = float(mapping_mismatch_by_sample.get(s, 0.0))
                coverage_deficit = _coverage_deficit(a, b)

                # non-phylo metrics
                for met in ("jaccard", "weighted_jaccard", "braycurtis", "aitchison"):
                    if met not in requested_metrics:
                        continue
                    is_bounded_metric = met != "aitchison"
                    if met == "jaccard":
                        # Treat every unmapped inferred feature as its own
                        # present/absent element so Jaccard does not collapse
                        # all unmapped inferred mass into a single pseudo-feature.
                        unmapped_feature_count = int(unmapped_feature_count_by_sample.get(s, 0))
                        if unmapped_feature_count <= 0 and unmapped_mass > 0.0 and not unmapped_feature_count_by_sample:
                            unmapped_feature_count = 1
                        d = float(_jaccard(a, b, extra_inferred_only=unmapped_feature_count))
                        d = _apply_penalties(
                            d,
                            0.0,
                            mapping_mismatch,
                            coverage_deficit,
                            use_additive_unmapped_penalty,
                            unmapped_penalty_lambda,
                            use_additive_mapping_mismatch_penalty,
                            mapping_mismatch_lambda,
                            use_additive_coverage_deficit_penalty,
                            coverage_deficit_lambda,
                            penalty_transform,
                            penalty_exp_k,
                            bounded_0_1=is_bounded_metric,
                        )
                    else:
                        # For the abundance-sensitive metrics, keep the current
                        # behavior: collapse all unmapped mass into one
                        # inferred-only pseudo-feature.
                        if unmapped_mass > 0.0:
                            a_np = np.append(a, 0.0)
                            b_np = np.append(b, unmapped_mass)
                            if met == "aitchison":
                                d = float(_aitchison(a_np, b_np, pseudocount=aitchison_pseudocount))
                            else:
                                d = float(metric_fns[met](a_np, b_np))
                        else:
                            if met == "aitchison":
                                d = float(_aitchison(a, b, pseudocount=aitchison_pseudocount))
                            else:
                                d = float(metric_fns[met](a, b))
                        d = _apply_penalties(
                            d,
                            unmapped_fraction if unmapped_mass <= 0.0 else 0.0,
                            mapping_mismatch,
                            coverage_deficit,
                            use_additive_unmapped_penalty,
                            unmapped_penalty_lambda,
                            use_additive_mapping_mismatch_penalty,
                            mapping_mismatch_lambda,
                            use_additive_coverage_deficit_penalty,
                            coverage_deficit_lambda,
                            penalty_transform,
                            penalty_exp_k,
                            bounded_0_1=is_bounded_metric,
                        )
                    comm_rows.append(
                        {
                            "community_id": cid,
                            "sample_id": s,
                            "method": method,
                            "metric": met,
                            "distance": d,
                        }
                    )

                # phylo metrics (vectorize on tree tips; restrict to what tree contains)
                if any(m in requested_metrics for m in ("unifrac", "gunifrac_a0.0", "gunifrac_a0.5", "gunifrac_a1.0")):
                    if phylo_tree is None:
                        if "unifrac" in requested_metrics:
                            d = _apply_penalties(
                                float("nan"),
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "unifrac",
                                    "distance": d,
                                }
                            )
                        if "gunifrac_a0.0" in requested_metrics:
                            d = _apply_penalties(
                                float("nan"),
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "gunifrac_a0.0",
                                    "distance": d,
                                }
                            )
                        if "gunifrac_a0.5" in requested_metrics:
                            d = _apply_penalties(
                                float("nan"),
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "gunifrac_a0.5",
                                    "distance": d,
                                }
                            )
                        if "gunifrac_a1.0" in requested_metrics:
                            d = _apply_penalties(
                                float("nan"),
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "gunifrac_a1.0",
                                    "distance": d,
                                }
                            )
                    else:
                        if phylo_mapping_mode == "augmented_mismatch" and std_table is not None:
                            va = _vectorize_on_tree(expected, s, phylo_tree)
                            vb = _vectorize_augmented_inferred(
                                table_std=std_table,
                                sample=s,
                                n_tips=len(phylo_tree.tips),
                                feat_to_syn_idx=feat_to_syn_idx,
                            )
                        else:
                            va = _vectorize_on_tree(expected, s, phylo_tree)
                            vb = _vectorize_on_tree(inferred, s, phylo_tree)

                        if "unifrac" in requested_metrics:
                            d = _safe_unifrac(phylo_tree, va, vb)
                            d = _apply_penalties(
                                d,
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "unifrac",
                                    "distance": d,
                                }
                            )
                        if "gunifrac_a0.0" in requested_metrics:
                            d = _safe_gunifrac(phylo_tree, va, vb, alpha=0.0)
                            d = _apply_penalties(
                                d,
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "gunifrac_a0.0",
                                    "distance": d,
                                }
                            )
                        if "gunifrac_a0.5" in requested_metrics:
                            d = _safe_gunifrac(phylo_tree, va, vb, alpha=0.5)
                            d = _apply_penalties(
                                d,
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "gunifrac_a0.5",
                                    "distance": d,
                                }
                            )
                        if "gunifrac_a1.0" in requested_metrics:
                            d = _safe_gunifrac(phylo_tree, va, vb, alpha=1.0)
                            d = _apply_penalties(
                                d,
                                unmapped_fraction,
                                mapping_mismatch,
                                coverage_deficit,
                                use_additive_unmapped_penalty,
                                unmapped_penalty_lambda,
                                use_additive_mapping_mismatch_penalty,
                                mapping_mismatch_lambda,
                                use_additive_coverage_deficit_penalty,
                                coverage_deficit_lambda,
                                penalty_transform,
                                penalty_exp_k,
                                bounded_0_1=True,
                            )
                            comm_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "method": method,
                                    "metric": "gunifrac_a1.0",
                                    "distance": d,
                                }
                            )

        # write per-community
        if comm_rows:
            dd = distances_dir(outdir, cid)
            dd.mkdir(parents=True, exist_ok=True)

            long_df = pd.DataFrame(comm_rows)
            long_df.to_csv(dd / "distances_long.tsv", sep="\t", index=False)

            wide = long_df.pivot_table(
                index=["community_id", "sample_id"],
                columns=["method", "metric"],
                values="distance",
                aggfunc="first",
            )
            wide.columns = [f"{m}|{k}" for (m, k) in wide.columns.to_list()]
            wide.reset_index().to_csv(dd / "distances_wide.tsv", sep="\t", index=False)

            log.info(f"[{cid}] wrote distances to {dd}")

        all_rows.extend(comm_rows)

    # write global
    if all_rows:
        all_long = pd.DataFrame(all_rows)
        all_long.to_csv(outdir / "distances_all_long.tsv", sep="\t", index=False)

        all_wide = all_long.pivot_table(
            index=["community_id", "sample_id"],
            columns=["method", "metric"],
            values="distance",
            aggfunc="first",
        )
        all_wide.columns = [f"{m}|{k}" for (m, k) in all_wide.columns.to_list()]
        all_wide.reset_index().to_csv(outdir / "distances_all_wide.tsv", sep="\t", index=False)

        log.info("Wrote global distances TSVs")

    return 0
