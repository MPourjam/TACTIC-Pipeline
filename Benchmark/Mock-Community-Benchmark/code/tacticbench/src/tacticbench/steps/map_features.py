from __future__ import annotations

import json
from pathlib import Path
import random

import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.blast import BlastTools, filter_hits, make_db, read_hits, run_blast, stable_one_to_one
from ..utils.io import fasta_ids, genus_species, read_expected_taxonomy, tip_to_species_map
from ..utils.paths import expected_dir, method_dir, qc_dir


METHODS = ["ASV", "zOTU", "TIC", "TAC", "OTU"]
CLUSTERING = {"TIC", "TAC", "OTU"}  # map to species (genus_species) not tips
METHOD_KEY_ALIASES = {
    "ASV": ["asv", "dada2"],
    "zOTU": ["zotu", "unoise3"],
    "OTU": ["otu", "uparse"],
    "TIC": ["tic"],
    "TAC": ["tac"],
}


def _normalize_soft_score_field(raw: str) -> str:
    v = str(raw or "").strip().lower()
    if v in {"bitscore", "pident", "length", "uniform"}:
        return v
    return "bitscore"


def _as_percent_identity(raw: float) -> float:
    v = float(raw)
    if 0.0 <= v <= 1.0:
        v *= 100.0
    return float(v)


def _resolve_target_level(method: str, target_cfg: object) -> str:
    default_level = "species" if method in CLUSTERING else "tip"
    if not isinstance(target_cfg, dict) or not target_cfg:
        return default_level
    for k in METHOD_KEY_ALIASES.get(method, []):
        v = str(target_cfg.get(k, "")).strip().lower()
        if v in {"tip", "species"}:
            return v
    vdef = str(target_cfg.get("default", "")).strip().lower()
    if vdef in {"tip", "species"}:
        return vdef
    return default_level


def _resolve_min_effective(cfg, method: str) -> float:
    by_method = getattr(cfg.mapping, "min_effective_pident_by_method", {})
    if isinstance(by_method, dict) and by_method:
        for k in METHOD_KEY_ALIASES.get(method, []):
            if k in by_method:
                return _as_percent_identity(float(by_method[k]))
        if "default" in by_method:
            return _as_percent_identity(float(by_method["default"]))

    min_effective = _as_percent_identity(float(cfg.mapping.min_effective_pident))
    if method in {"TIC", "TAC"} and getattr(cfg.mapping, "min_effective_pident_tic_tac", None) is not None:
        min_effective = _as_percent_identity(float(cfg.mapping.min_effective_pident_tic_tac))
    elif method in CLUSTERING and getattr(cfg.mapping, "min_effective_pident_clustering", None) is not None:
        min_effective = _as_percent_identity(float(cfg.mapping.min_effective_pident_clustering))
    return min_effective


def _mapping_signature(
    *,
    method: str,
    target_level: str,
    min_effective: float,
    blast_min_ident: float,
    min_qcov: float,
    assignment_mode: str,
    soft_score_field: str,
    soft_min_rel: float,
    seed: int,
) -> dict:
    tie_policy = "shuffle_equal_score_candidates"
    if str(assignment_mode).strip().lower() == "tie_aware_split":
        tie_policy = "explicit_top_ties_equal_weight"
    return {
        "version": 4,
        "method": str(method),
        "target_level": str(target_level),
        "min_effective_pident": round(float(min_effective), 12),
        "blast_record_min_raw_pident": round(float(blast_min_ident), 12),
        "min_query_coverage": round(float(min_qcov), 12),
        "assignment_mode": str(assignment_mode),
        "soft_species_score_field": str(soft_score_field),
        "soft_species_min_relative_score": round(float(soft_min_rel), 12),
        "exact_tie_policy": tie_policy,
        "seed": int(seed),
    }


def _score_key(rec: dict) -> tuple[float, float, float, float, float]:
    """Ranking key used for one-to-one candidate comparisons."""
    return (
        float(rec.get("effective_identity", rec.get("pident", 0.0))),
        float(rec.get("pident", 0.0)),
        float(rec.get("bitscore", 0.0)),
        -float(rec.get("evalue", 0.0)),
        float(rec.get("length", 0.0)),
    )


def _shuffle_equal_score_blocks(opts: list[dict], rng: random.Random) -> list[dict]:
    """
    Remove first-row bias by shuffling candidates that are exactly tied by score.

    We keep block order by score (best->worst), but randomize within each
    equal-score block using the seeded RNG for reproducibility.
    """
    if len(opts) <= 1:
        return opts

    out: list[dict] = []
    i = 0
    n = len(opts)
    while i < n:
        j = i + 1
        key_i = _score_key(opts[i])
        while j < n and _score_key(opts[j]) == key_i:
            j += 1
        block = opts[i:j]
        if len(block) > 1:
            rng.shuffle(block)
        out.extend(block)
        i = j
    return out


def _top_tied_options(opts: list[dict]) -> list[dict]:
    """Return the best-score tie block from a best->worst candidate list."""
    if not opts:
        return []
    best = _score_key(opts[0])
    out: list[dict] = []
    for rec in opts:
        if _score_key(rec) != best:
            break
        out.append(rec)
    return out


def _weights_from_group(g: pd.DataFrame, score_field: str, min_relative_score: float) -> pd.DataFrame:
    out = g.copy()
    if score_field == "uniform":
        out["__score"] = 1.0
    else:
        out["__score"] = pd.to_numeric(out[score_field], errors="coerce").fillna(0.0).clip(lower=0.0)

    rel = float(max(0.0, min(1.0, float(min_relative_score))))
    max_s = float(out["__score"].max()) if not out.empty else 0.0
    if rel > 0.0 and max_s > 0.0:
        out = out[out["__score"] >= (max_s * rel)].copy()
        if out.empty:
            return out

    total = float(out["__score"].sum())
    if total <= 0:
        out["weight"] = 1.0 / float(len(out))
    else:
        out["weight"] = out["__score"] / total
    return out


def run(cfg) -> int:
    """
    Create mapping files for inferred features -> expected references using BLAST.

    Reference:
      - expected sequences: outdir/communities/<cid>/expected/expected_sequences.fasta
      - BLAST DB:           outdir/communities/<cid>/expected/blastdb_expected.*

    For each method:
      - query sequences: outdir/communities/<cid>/methods/<method>/seqs.fasta
      - run BLAST and store hits: blast_hits.tsv
        (BLAST recording threshold configurable by mapping.blast_record_min_raw_pident)

    Mapping logic:
      - Accept hits with effective_identity >= cfg.mapping.min_effective_pident
        where effective_identity = pident * qcov_exact
        and qcov_exact = (abs(qend - qstart) + 1) / qlen
      - Optionally require minimum query coverage.
      - ASV/zOTU map to expected TIP IDs (feature-level)
      - TIC/TAC/OTU map to SPECIES IDs (genus_species):
          - convert each hit's sseqid (tip id) -> species via expected tip->species map
          - for each (query, species), keep the best hit
            (highest effective_identity, then pident, then bitscore, then lower evalue, then longer)
      - assignment_mode controls final assignment:
          - one_to_one:
              enforce unique targets with stable_one_to_one()
              and seeded shuffling inside exact score-tie blocks
          - tie_aware_split:
              keep only the top score-tie block per inferred feature and
              split weight equally across tied targets (many-to-many allowed)
          - soft_species:
              ASV/zOTU + species target only, weighted soft assignment across
              species candidates

    Outputs per method:
      - mapping.tsv
      - unmapped_features.tsv

    mapping.tsv columns:
      - for ASV/zOTU:
          one_to_one mode:
            inferred_feature, target_tip, effective_identity, pident, bitscore, evalue, length
          tie_aware_split mode:
            tip target:
              inferred_feature, target_tip, weight, tie_size,
              effective_identity, pident, bitscore, evalue, length
            species target:
              inferred_feature, target_species, best_tip, target_tip, weight, tie_size,
              effective_identity, pident, bitscore, evalue, length
          soft_species mode:
            inferred_feature, target_species, best_tip, target_tip, weight, effective_identity, pident, bitscore, evalue, length
      - for TIC/TAC/OTU:
          one_to_one mode:
            inferred_feature, target_species, best_tip, target_tip, effective_identity, pident, bitscore, evalue, length
          tie_aware_split mode:
            inferred_feature, target_species, best_tip, target_tip, weight, tie_size,
            effective_identity, pident, bitscore, evalue, length
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    manifest = pd.read_csv(outdir / "manifest.tsv", sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    tools = BlastTools(makeblastdb=cfg.tree.makeblastdb_exe, blastn=cfg.tree.blastn_exe)
    rng = random.Random(int(cfg.mapping.seed))
    assignment_mode = str(getattr(cfg.mapping, "assignment_mode", "one_to_one")).strip().lower()
    if assignment_mode not in {"one_to_one", "soft_species", "tie_aware_split"}:
        assignment_mode = "one_to_one"
    soft_score_field = _normalize_soft_score_field(getattr(cfg.mapping, "soft_species_score_field", "bitscore"))
    soft_min_rel = float(getattr(cfg.mapping, "soft_species_min_relative_score", 0.0))
    mapping_min_qcov = float(getattr(cfg.mapping, "min_query_coverage", 0.0))

    for _, r in manifest.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue

        ed = expected_dir(outdir, cid)
        ref_fasta = ed / "expected_sequences.fasta"
        if not ref_fasta.exists():
            continue

        db_prefix = ed / "blastdb_expected"

        # build DB if missing
        if not (ed / "blastdb_expected.nsq").exists() and not (ed / "blastdb_expected.00.nsq").exists():
            ok = make_db(tools, ref_fasta, db_prefix)
            if not ok:
                log.info(f"[{cid}] BLAST DB build failed (BLAST missing?); skipping mapping.")
                continue

        # tip->species mapping
        exp_tax_path = Path(str(r.get("expected_taxonomy", "")))
        if exp_tax_path.exists():
            raw = read_expected_taxonomy(exp_tax_path)
            tip_to_sp = tip_to_species_map(ref_fasta, raw)
        else:
            tip_to_sp = {}

        for method in METHODS:
            md = method_dir(outdir, cid, method)
            q_fasta = md / "seqs.fasta"
            if not q_fasta.exists():
                continue

            try:
                all_query_features = [str(x) for x in fasta_ids(q_fasta)]
            except Exception:
                all_query_features = []

            target_level = _resolve_target_level(method, getattr(cfg.mapping, "target_level_by_method", {}))
            species_target = target_level == "species"
            min_effective = _resolve_min_effective(cfg, method)

            # Decouple BLAST output threshold from mapping threshold:
            # blast_hits.tsv can retain lower-identity hits for diagnostics/exploration,
            # while final mapping still applies method-specific min_ident below.
            blast_record_min_raw = getattr(cfg.mapping, "blast_record_min_raw_pident", None)
            if blast_record_min_raw is None:
                blast_min_ident = float(min_effective)
            else:
                blast_min_ident = float(_as_percent_identity(float(blast_record_min_raw)))
                if blast_min_ident > float(min_effective):
                    log.info(
                        f"[{cid}] {method}: blast_record_min_raw_pident ({blast_min_ident}) > "
                        f"mapping threshold ({min_effective}); this may reduce retained hits."
                    )

            hits_path = md / "blast_hits.tsv"
            map_path = md / "mapping.tsv"
            params_path = md / "mapping_params.json"
            desired_sig = _mapping_signature(
                method=method,
                target_level=target_level,
                min_effective=min_effective,
                blast_min_ident=blast_min_ident,
                min_qcov=mapping_min_qcov,
                assignment_mode=assignment_mode,
                soft_score_field=soft_score_field,
                soft_min_rel=soft_min_rel,
                seed=int(cfg.mapping.seed),
            )
            if map_path.exists() and map_path.stat().st_size > 0:
                regen_reason = "missing required columns"
                mapping_ok = False
                try:
                    existing_map = pd.read_csv(map_path, sep="\t")
                    existing_cols = set(existing_map.columns.astype(str).tolist())
                except Exception:
                    existing_map = pd.DataFrame()
                    existing_cols = set()
                if {"target_tip", "effective_identity"}.issubset(existing_cols):
                    # Guard against legacy/corrupted effective_identity values
                    # where effective_identity may exceed pident.
                    if {"effective_identity", "pident"}.issubset(existing_cols):
                        eff = pd.to_numeric(existing_map["effective_identity"], errors="coerce")
                        pid = pd.to_numeric(existing_map["pident"], errors="coerce")
                        if np.isfinite(eff).any() and np.isfinite(pid).any():
                            if ((eff - pid) > 1e-9).any():
                                regen_reason = "legacy effective_identity (> pident)"
                            else:
                                mapping_ok = True
                        else:
                            mapping_ok = True
                    else:
                        mapping_ok = True

                if mapping_ok:
                    hits_ok = False
                    if hits_path.exists() and hits_path.stat().st_size > 0:
                        try:
                            n_cols = int(pd.read_csv(hits_path, sep="\t", header=None, nrows=1).shape[1])
                        except Exception:
                            n_cols = 0
                        if n_cols >= 10:
                            hits_ok = True
                        else:
                            regen_reason = "legacy blast_hits.tsv missing qstart/qend columns"
                    else:
                        regen_reason = "missing blast_hits.tsv required for strict effective_identity"
                    if hits_ok:
                        try:
                            existing_sig = {}
                            if params_path.exists() and params_path.stat().st_size > 0:
                                existing_sig = json.loads(params_path.read_text(encoding="utf-8"))
                            if existing_sig == desired_sig:
                                continue
                            regen_reason = "mapping parameters changed"
                        except Exception:
                            regen_reason = "mapping params missing/invalid"
                log.info(f"[{cid}] {method}: existing mapping.tsv {regen_reason}; regenerating.")

            ok = run_blast(
                tools=tools,
                query_fasta=q_fasta,
                db_prefix=db_prefix,
                out_tsv=hits_path,
                min_ident=blast_min_ident,
            )
            if not ok or (not hits_path.exists()) or hits_path.stat().st_size == 0:
                log.info(f"[{cid}] {method}: no BLAST hits / blast failed.")
                pd.DataFrame(columns=["inferred_feature"]).to_csv(md / "unmapped_features.tsv", sep="\t", index=False)
                try:
                    params_path.write_text(json.dumps(desired_sig, sort_keys=True, indent=2) + "\n", encoding="utf-8")
                except Exception:
                    pass
                continue

            hits = read_hits(hits_path)
            hits = filter_hits(
                hits,
                min_effective_pident=min_effective,
                min_qcovhsp=mapping_min_qcov,
                strict_pident_gt=False,
            )
            if hits.empty:
                log.info(
                    f"[{cid}] {method}: no BLAST hits passed filters "
                    f"(effective_pident>={min_effective}, qcov>={mapping_min_qcov})."
                )
                if (assignment_mode == "soft_species") and (method not in CLUSTERING):
                    empty_cols = [
                        "inferred_feature",
                        "target_species",
                        "best_tip",
                        "target_tip",
                        "weight",
                        "effective_identity",
                        "pident",
                        "bitscore",
                        "evalue",
                        "length",
                    ]
                elif assignment_mode == "tie_aware_split":
                    if species_target:
                        empty_cols = [
                            "inferred_feature",
                            "target_species",
                            "best_tip",
                            "target_tip",
                            "weight",
                            "tie_size",
                            "effective_identity",
                            "pident",
                            "bitscore",
                            "evalue",
                            "length",
                        ]
                    else:
                        empty_cols = [
                            "inferred_feature",
                            "target_tip",
                            "weight",
                            "tie_size",
                            "effective_identity",
                            "pident",
                            "bitscore",
                            "evalue",
                            "length",
                        ]
                elif method in CLUSTERING:
                    empty_cols = [
                        "inferred_feature",
                        "target_species",
                        "best_tip",
                        "target_tip",
                        "effective_identity",
                        "pident",
                        "bitscore",
                        "evalue",
                        "length",
                    ]
                else:
                    empty_cols = [
                        "inferred_feature",
                        "target_tip",
                        "effective_identity",
                        "pident",
                        "bitscore",
                        "evalue",
                        "length",
                    ]
                pd.DataFrame(columns=empty_cols).to_csv(map_path, sep="\t", index=False)
                pd.DataFrame({"inferred_feature": all_query_features}).to_csv(
                    md / "unmapped_features.tsv", sep="\t", index=False
                )
                try:
                    params_path.write_text(json.dumps(desired_sig, sort_keys=True, indent=2) + "\n", encoding="utf-8")
                except Exception:
                    pass
                continue

            # sort best-first for each query
            hits = hits.sort_values(
                ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                ascending=[True, False, False, False, True, False],
            )

            rows = []
            unmapped: list[str] = []

            use_soft_species = (assignment_mode == "soft_species") and (method not in CLUSTERING) and species_target
            use_tie_aware_split = assignment_mode == "tie_aware_split"
            if use_soft_species:
                # Soft species-level assignment for ASV/zOTU:
                # keep candidate species hits and distribute each inferred feature by normalized weights.
                hs = hits.copy()
                hs["species"] = hs["sseqid"].map(lambda x: tip_to_sp.get(str(x), genus_species(str(x))))
                best_q_sp = hs.groupby(["qseqid", "species"], sort=False).first().reset_index()
                best_q_sp = best_q_sp.sort_values(
                    ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                    ascending=[True, False, False, False, True, False],
                )

                mapped_features: set[str] = set()
                for q, g in best_q_sp.groupby("qseqid", sort=False):
                    gw = _weights_from_group(g, score_field=soft_score_field, min_relative_score=soft_min_rel)
                    if gw.empty:
                        continue
                    qid = str(q)
                    mapped_features.add(qid)
                    for _, rec in gw.iterrows():
                        rows.append(
                            {
                                "inferred_feature": qid,
                                "target_species": str(rec["species"]),
                                "best_tip": str(rec["sseqid"]),
                                "target_tip": str(rec["sseqid"]),
                                "weight": float(rec["weight"]),
                                "effective_identity": float(rec.get("effective_identity", rec["pident"])),
                                "pident": float(rec["pident"]),
                                "bitscore": float(rec["bitscore"]),
                                "evalue": float(rec["evalue"]),
                                "length": float(rec["length"]),
                            }
                        )
                all_q = set(all_query_features) if all_query_features else set(hs["qseqid"].astype(str).tolist())
                unmapped = sorted([q for q in all_q if q not in mapped_features])
            elif use_tie_aware_split:
                mapped_features: set[str] = set()
                if species_target:
                    hs = hits.copy()
                    hs["species"] = hs["sseqid"].map(
                        lambda x: tip_to_sp.get(str(x), genus_species(str(x)))
                    )
                    best_q_sp = hs.groupby(["qseqid", "species"], sort=False).first().reset_index()
                    best_q_sp = best_q_sp.sort_values(
                        ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                        ascending=[True, False, False, False, True, False],
                    )
                    for q, g in best_q_sp.groupby("qseqid", sort=False):
                        opts = []
                        for _, rec in g.iterrows():
                            opts.append(
                                {
                                    "qseqid": str(rec["qseqid"]),
                                    "target": str(rec["species"]),
                                    "best_tip": str(rec["sseqid"]),
                                    "effective_identity": float(rec.get("effective_identity", rec["pident"])),
                                    "pident": float(rec["pident"]),
                                    "length": float(rec["length"]),
                                    "bitscore": float(rec["bitscore"]),
                                    "evalue": float(rec["evalue"]),
                                }
                            )
                        ties = _top_tied_options(opts)
                        if not ties:
                            continue
                        mapped_features.add(str(q))
                        tie_size = int(len(ties))
                        w = 1.0 / float(tie_size)
                        for rec in ties:
                            rows.append(
                                {
                                    "inferred_feature": str(q),
                                    "target_species": rec["target"],
                                    "best_tip": rec.get("best_tip", ""),
                                    "target_tip": rec.get("best_tip", ""),
                                    "weight": w,
                                    "tie_size": tie_size,
                                    "effective_identity": rec.get("effective_identity", rec["pident"]),
                                    "pident": rec["pident"],
                                    "bitscore": rec["bitscore"],
                                    "evalue": rec["evalue"],
                                    "length": rec["length"],
                                }
                            )
                    all_q = set(all_query_features) if all_query_features else set(best_q_sp["qseqid"].astype(str).tolist())
                else:
                    for q, g in hits.groupby("qseqid", sort=False):
                        opts = []
                        for _, rec in g.iterrows():
                            opts.append(
                                {
                                    "qseqid": str(rec["qseqid"]),
                                    "target": str(rec["sseqid"]),
                                    "effective_identity": float(rec.get("effective_identity", rec["pident"])),
                                    "pident": float(rec["pident"]),
                                    "length": float(rec["length"]),
                                    "bitscore": float(rec["bitscore"]),
                                    "evalue": float(rec["evalue"]),
                                }
                            )
                        ties = _top_tied_options(opts)
                        if not ties:
                            continue
                        mapped_features.add(str(q))
                        tie_size = int(len(ties))
                        w = 1.0 / float(tie_size)
                        for rec in ties:
                            rows.append(
                                {
                                    "inferred_feature": str(q),
                                    "target_tip": rec["target"],
                                    "weight": w,
                                    "tie_size": tie_size,
                                    "effective_identity": rec.get("effective_identity", rec["pident"]),
                                    "pident": rec["pident"],
                                    "bitscore": rec["bitscore"],
                                    "evalue": rec["evalue"],
                                    "length": rec["length"],
                                }
                            )
                    all_q = set(all_query_features) if all_query_features else set(hits["qseqid"].astype(str).tolist())
                unmapped = sorted([q for q in all_q if q not in mapped_features])
            else:
                candidates: dict[str, list[dict]] = {}
                if species_target:
                    # map hits to species targets (and keep best tip per species)
                    hits = hits.copy()
                    hits["species"] = hits["sseqid"].map(
                        lambda x: tip_to_sp.get(str(x), genus_species(str(x)))
                    )

                    # pick best hit per (q, species) by taking first row after sorting
                    best_q_sp = hits.groupby(["qseqid", "species"], sort=False).first().reset_index()
                    best_q_sp = best_q_sp.sort_values(
                        ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                        ascending=[True, False, False, False, True, False],
                    )

                    for q, g in best_q_sp.groupby("qseqid", sort=False):
                        opts = []
                        for _, rec in g.iterrows():
                            opts.append(
                                {
                                    "qseqid": str(rec["qseqid"]),
                                    "target": str(rec["species"]),  # target species
                                    "best_tip": str(rec["sseqid"]),
                                    "effective_identity": float(rec.get("effective_identity", rec["pident"])),
                                    "pident": float(rec["pident"]),
                                    "length": float(rec["length"]),
                                    "bitscore": float(rec["bitscore"]),
                                    "evalue": float(rec["evalue"]),
                                }
                            )
                        candidates[str(q)] = _shuffle_equal_score_blocks(opts, rng)
                else:
                    # map to tip targets directly
                    for q, g in hits.groupby("qseqid", sort=False):
                        opts = []
                        for _, rec in g.iterrows():
                            opts.append(
                                {
                                    "qseqid": str(rec["qseqid"]),
                                    "target": str(rec["sseqid"]),  # target tip
                                    "effective_identity": float(rec.get("effective_identity", rec["pident"])),
                                    "pident": float(rec["pident"]),
                                    "length": float(rec["length"]),
                                    "bitscore": float(rec["bitscore"]),
                                    "evalue": float(rec["evalue"]),
                                }
                            )
                        candidates[str(q)] = _shuffle_equal_score_blocks(opts, rng)

                # 1-to-1 assignment on targets
                assigned = stable_one_to_one(candidates, rng=rng)

                # write mapping
                for q, rec in assigned.items():
                    if species_target:
                        rows.append(
                            {
                                "inferred_feature": q,
                                "target_species": rec["target"],
                                "best_tip": rec.get("best_tip", ""),
                                "target_tip": rec.get("best_tip", ""),
                                "effective_identity": rec.get("effective_identity", rec["pident"]),
                                "pident": rec["pident"],
                                "bitscore": rec["bitscore"],
                                "evalue": rec["evalue"],
                                "length": rec["length"],
                            }
                        )
                    else:
                        rows.append(
                            {
                                "inferred_feature": q,
                                "target_tip": rec["target"],
                                "effective_identity": rec.get("effective_identity", rec["pident"]),
                                "pident": rec["pident"],
                                "bitscore": rec["bitscore"],
                                "evalue": rec["evalue"],
                                "length": rec["length"],
                            }
                        )

                all_q = all_query_features if all_query_features else list(candidates.keys())
                unmapped = [q for q in all_q if q not in assigned]

            pd.DataFrame(rows).to_csv(map_path, sep="\t", index=False)
            pd.DataFrame({"inferred_feature": unmapped}).to_csv(
                md / "unmapped_features.tsv", sep="\t", index=False
            )
            try:
                params_path.write_text(json.dumps(desired_sig, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            except Exception:
                pass

            if use_soft_species:
                n_features = hits["qseqid"].astype(str).nunique()
                log.info(f"[{cid}] {method}: soft-species mapped {n_features - len(unmapped)}/{n_features} inferred features")
            else:
                log.info(f"[{cid}] {method}: mapped {len(rows)}/{len(rows) + len(unmapped)}")

    return 0
