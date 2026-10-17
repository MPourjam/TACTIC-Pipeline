from __future__ import annotations

import copy
from fnmatch import fnmatch
from itertools import combinations
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.blast import filter_hits, read_hits
from ..utils.labels import display_method
from ..utils.normalize import tss
from ..utils.paths import figures_root
from ..utils.io import canonical_species_name

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}
METRIC_LABEL = {
    "jaccard": "Jaccard",
    "weighted_jaccard": "Weighted Jaccard",
    "braycurtis": "Bray-Curtis",
    "unifrac": "UniFrac",
    "gunifrac_a0.0": "gUniFrac (alpha=0.0)",
    "gunifrac_a0.5": "gUniFrac (alpha=0.5)",
    "gunifrac_a1.0": "gUniFrac (alpha=1.0)",
}
AUTO_BETA_METRICS = {"weighted_jaccard", "braycurtis", "gunifrac_a0.0", "gunifrac_a0.5", "gunifrac_a1.0"}
CLUSTERING = {"TIC", "TAC", "OTU"}
METHOD_KEY_ALIASES = {
    "ASV": ["asv", "dada2"],
    "zOTU": ["zotu", "unoise3"],
    "OTU": ["otu", "uparse"],
    "TIC": ["tic"],
    "TAC": ["tac"],
}


def _normalize_plot_type(raw: str, default: str = "violin") -> str:
    v = str(raw or "").strip().lower()
    if v in {"box", "boxplot"}:
        return "boxplot"
    if v in {"violin", "violinplot"}:
        return "violin"
    return default


def _normalize_pairwise_test(raw: str, default: str = "sign_test_exact") -> str:
    v = str(raw or "").strip().lower()
    if v in {"signed_rank_permutation", "sign_test_exact"}:
        return v
    return default


def _normalize_stats_model(raw: str, default: str = "mixed_effect_glmmtmb_auto") -> str:
    v = str(raw or "").strip().lower()
    if v in {"mixed_effect_glmmtmb_auto", "mixed_effect_glmmtmb_beta", "paired_permutation"}:
        return v
    return default


def _normalize_projection_mapping_mode(raw: str, default: str = "standard") -> str:
    v = str(raw or "").strip().lower()
    if v in {"best_hit_aggregate", "best-hit-aggregate", "best_hit", "best-hit"}:
        return "best_hit_aggregate"
    return default


def _normalize_emmeans_scale(raw: str, default: str = "link") -> str:
    v = str(raw or "").strip().lower()
    if v in {"link", "response"}:
        return v
    return default


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _as_percent_identity(raw: float) -> float:
    v = float(raw)
    if 0.0 <= v <= 1.0:
        v *= 100.0
    return float(v)


def _resolve_target_level_override(cfg, method: str) -> str:
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


def _resolve_min_effective_override(cfg, method: str) -> float:
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


def _best_hit_aggregate_projection(cfg, outdir: Path, log, figure_label: str = "Figure 11") -> None:
    """
    Figure 11-only projection mode:
      - tip target: choose best hit(s) per inferred feature above threshold;
        split abundance equally across ties; unmapped inferred features are
        excluded from projected mapped tables.
      - species target: choose best hit(s) per inferred feature above threshold,
        collapse tied hits to unique species targets, then split abundance
        equally across tied species; unmapped inferred features are excluded.
    """
    manifest_p = outdir / "manifest.tsv"
    if not manifest_p.exists():
        return
    try:
        manifest = pd.read_csv(manifest_p, sep="\t")
    except Exception:
        return

    requested_methods = list(getattr(cfg.distances, "methods", METHOD_ORDER))
    qcov = float(getattr(cfg.mapping, "min_query_coverage", 0.0))

    for cid in manifest.get("community_id", pd.Series(dtype=str)).astype(str).tolist():
        ed = outdir / "communities" / cid / "expected"
        exp_tip_p = ed / "expected_tip_abundance_tss1000.tsv"
        exp_sp_p = ed / "expected_species_abundance_tss1000.tsv"
        tip_to_species_p = ed / "tip_to_species.tsv"
        if not exp_tip_p.exists() or not exp_sp_p.exists():
            continue

        try:
            exp_tip = pd.read_csv(exp_tip_p, sep="\t", index_col=0)
            exp_sp = pd.read_csv(exp_sp_p, sep="\t", index_col=0)
        except Exception:
            continue
        exp_tip.index = exp_tip.index.astype(str)
        exp_sp.index = exp_sp.index.astype(str)

        tip_to_species: dict[str, str] = {}
        if tip_to_species_p.exists():
            try:
                t2s = pd.read_csv(tip_to_species_p, sep="\t")
                if {"tip", "species"}.issubset(set(t2s.columns)):
                    tip_to_species = dict(zip(t2s["tip"].astype(str), t2s["species"].astype(str)))
            except Exception:
                tip_to_species = {}

        for method in requested_methods:
            md = outdir / "communities" / cid / "methods" / method
            table_std_p = md / "table_std_tss1000.tsv"
            blast_hits_p = md / "blast_hits.tsv"
            if not table_std_p.exists():
                continue
            try:
                tab = pd.read_csv(table_std_p, sep="\t", index_col=0)
            except Exception:
                continue
            tab.index = tab.index.astype(str)
            tab = tss(tab, total=cfg.normalization.tss_total, axis=0)

            try:
                hits = read_hits(blast_hits_p) if blast_hits_p.exists() else pd.DataFrame()
            except Exception:
                hits = pd.DataFrame()
            min_effective = _resolve_min_effective_override(cfg, method)
            if not hits.empty:
                hits = filter_hits(
                    hits,
                    min_effective_pident=min_effective,
                    min_qcovhsp=qcov,
                    strict_pident_gt=False,
                )
                hits = hits.sort_values(
                    ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                    ascending=[True, False, False, False, True, False],
                )

            target_level = _resolve_target_level_override(cfg, method)
            species_target = target_level == "species"

            mapped_feats: set[str] = set()
            if species_target:
                out = pd.DataFrame(0.0, index=exp_sp.index.astype(str), columns=tab.columns)
                hs = hits.copy() if not hits.empty else pd.DataFrame()
                if not hs.empty:
                    hs["target_species"] = hs["sseqid"].map(
                        lambda x: canonical_species_name(str(tip_to_species.get(str(x), str(x))))
                    )
                for feat in tab.index.astype(str):
                    if hs.empty:
                        continue
                    g = hs[hs["qseqid"].astype(str) == feat]
                    if g.empty:
                        continue
                    best_eff = float(pd.to_numeric(g["effective_identity"], errors="coerce").max())
                    best = g[np.isclose(pd.to_numeric(g["effective_identity"], errors="coerce"), best_eff, atol=1e-12)]
                    species = sorted(set(best["target_species"].astype(str).tolist()))
                    species = [s for s in species if s in out.index]
                    if not species:
                        continue
                    mapped_feats.add(feat)
                    vals = tab.loc[feat, :].to_numpy(dtype=float) / float(len(species))
                    for spx in species:
                        out.loc[spx, :] += vals
                if cfg.normalization.renormalize_after_mapping:
                    out = tss(out, total=cfg.normalization.tss_total, axis=0)
                out.to_csv(md / "mapped_species_table_tss1000.tsv", sep="\t")
            else:
                out = pd.DataFrame(0.0, index=exp_tip.index.astype(str), columns=tab.columns)
                hs = hits.copy() if not hits.empty else pd.DataFrame()
                for feat in tab.index.astype(str):
                    g = hs[hs["qseqid"].astype(str) == feat] if not hs.empty else pd.DataFrame()
                    if g.empty:
                        continue
                    best_eff = float(pd.to_numeric(g["effective_identity"], errors="coerce").max())
                    best = g[np.isclose(pd.to_numeric(g["effective_identity"], errors="coerce"), best_eff, atol=1e-12)]
                    tips = sorted(set(best["sseqid"].astype(str).tolist()))
                    tips = [t for t in tips if t in out.index]
                    if not tips:
                        continue
                    mapped_feats.add(feat)
                    share = tab.loc[feat, :].to_numpy(dtype=float) / float(len(tips))
                    for t in tips:
                        out.loc[t, :] += share
                if cfg.normalization.renormalize_after_mapping:
                    out = tss(out, total=cfg.normalization.tss_total, axis=0)
                out.to_csv(md / "mapped_tip_table_tss1000.tsv", sep="\t")

            unmapped_feats = [f for f in tab.index.astype(str).tolist() if f not in mapped_feats]
            unmapped_mass = tab.loc[unmapped_feats].sum(axis=0) if unmapped_feats else pd.Series(
                0.0, index=tab.columns
            )
            mapped_mass = tab.drop(index=unmapped_feats, errors="ignore").sum(axis=0)
            qcdf = pd.DataFrame(
                {
                    "sample_id": tab.columns.astype(str),
                    "mapped_mass": mapped_mass.values,
                    "unmapped_mass": unmapped_mass.values,
                    "mapping_rate": np.where(
                        (mapped_mass + unmapped_mass).values > 0,
                        mapped_mass.values / (mapped_mass + unmapped_mass).values,
                        np.nan,
                    ),
                }
            )
            qcdf.to_csv(md / "mapping_qc.tsv", sep="\t", index=False)

        log.info(f"[{cid}] {figure_label} best-hit aggregate projection built")


def _prepare_fig11_override_workspace(src_outdir: Path, tmp_outdir: Path) -> bool:
    manifest_src = src_outdir / "manifest.tsv"
    if not manifest_src.exists():
        return False

    tmp_outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(manifest_src, tmp_outdir / "manifest.tsv")

    cmap_src = src_outdir / "community_name_mapping.tsv"
    if cmap_src.exists():
        shutil.copy2(cmap_src, tmp_outdir / "community_name_mapping.tsv")

    qc_src = src_outdir / "qc"
    if qc_src.exists():
        shutil.copytree(qc_src, tmp_outdir / "qc", dirs_exist_ok=True)

    man = pd.read_csv(manifest_src, sep="\t")
    for cid in man["community_id"].astype(str).tolist():
        src_comm = src_outdir / "communities" / cid
        if not src_comm.exists():
            continue
        dst_comm = tmp_outdir / "communities" / cid
        exp_src = src_comm / "expected"
        meth_src = src_comm / "methods"
        if exp_src.exists():
            shutil.copytree(exp_src, dst_comm / "expected", dirs_exist_ok=True)
        if meth_src.exists():
            shutil.copytree(meth_src, dst_comm / "methods", dirs_exist_ok=True)
    return True


def _load_fig11_distances(cfg, log) -> pd.DataFrame:
    outdir = Path(cfg.paths.outdir)
    long_p = outdir / "distances_all_long.tsv"

    target_overrides = getattr(cfg.figures, "fig11_target_level_by_method", {})
    min_ident_overrides = getattr(cfg.figures, "fig11_min_effective_pident_by_method", {})
    projection_mode = _normalize_projection_mapping_mode(
        getattr(cfg.figures, "fig11_projection_mapping_mode", "standard")
    )
    use_overrides = (
        isinstance(target_overrides, dict)
        and bool(target_overrides)
    ) or (
        isinstance(min_ident_overrides, dict)
        and bool(min_ident_overrides)
    ) or (projection_mode != "standard")

    if not use_overrides:
        if not long_p.exists():
            log.info("Figure 11: distances_all_long.tsv not found; run compute-distances first.")
            return pd.DataFrame()
        return pd.read_csv(long_p, sep="\t")

    log.info("Figure 11: using figure-specific overrides; recomputing temporary distances.")
    with tempfile.TemporaryDirectory(prefix="tacticbench_fig11_override_") as td:
        tmp_outdir = Path(td) / "out"
        ok = _prepare_fig11_override_workspace(src_outdir=outdir, tmp_outdir=tmp_outdir)
        if not ok:
            log.info("Figure 11: failed to prepare temporary workspace for figure-specific mapping overrides.")
            return pd.DataFrame()

        tmp_cfg = copy.deepcopy(cfg)
        tmp_cfg.paths.outdir = tmp_outdir
        if isinstance(target_overrides, dict) and target_overrides:
            tmp_cfg.mapping.target_level_by_method = dict(target_overrides)
        if isinstance(min_ident_overrides, dict) and min_ident_overrides:
            tmp_cfg.mapping.min_effective_pident_by_method = dict(min_ident_overrides)

        from . import compute_distances, map_features, project_abundances

        map_features.run(tmp_cfg)
        if projection_mode == "best_hit_aggregate":
            _best_hit_aggregate_projection(tmp_cfg, tmp_outdir, log)
        else:
            project_abundances.run(tmp_cfg)
        compute_distances.run(tmp_cfg)

        tmp_long_p = tmp_outdir / "distances_all_long.tsv"
        if not tmp_long_p.exists():
            log.info("Figure 11: temporary distances were not generated for figure-specific mapping overrides.")
            return pd.DataFrame()
        return pd.read_csv(tmp_long_p, sep="\t")


def _bh_adjust(pvals: list[float]) -> list[float]:
    out = [float("nan")] * len(pvals)
    finite_idx = [i for i, p in enumerate(pvals) if np.isfinite(p)]
    if not finite_idx:
        return out

    order = sorted(finite_idx, key=lambda i: float(pvals[i]))
    m = len(order)
    raw = [0.0] * m
    for rank1, idx in enumerate(order, start=1):
        raw[rank1 - 1] = float(pvals[idx]) * m / rank1

    mono = [0.0] * m
    running_min = 1.0
    for j in range(m - 1, -1, -1):
        running_min = min(running_min, raw[j])
        mono[j] = min(1.0, running_min)

    for j, idx in enumerate(order):
        out[idx] = mono[j]
    return out


def _pvalue_stars(p: float, alpha: float = 0.05) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < alpha:
        return "*"
    return "ns"


def _friedman_statistic(values: np.ndarray) -> float:
    if values.ndim != 2:
        return float("nan")
    n, k = values.shape
    if n < 2 or k < 2:
        return float("nan")

    ranks = pd.DataFrame(values).rank(axis=1, method="average").to_numpy(dtype=float)
    rank_sums = ranks.sum(axis=0)
    q = (12.0 / (n * k * (k + 1.0))) * float(np.sum(rank_sums**2)) - 3.0 * n * (k + 1.0)

    tie_num = 0.0
    for i in range(n):
        _, counts = np.unique(values[i, :], return_counts=True)
        tie_num += float(np.sum(counts**3 - counts))
    tie_den = float(n * (k**3 - k))
    if tie_den > 0:
        c = 1.0 - (tie_num / tie_den)
        if c > 0:
            q = q / c
    return float(q)


def _friedman_permutation_pvalue(values: np.ndarray, n_perm: int = 4000, seed: int = 1) -> tuple[float, float]:
    obs = _friedman_statistic(values)
    if not np.isfinite(obs):
        return float("nan"), float("nan")

    rng = np.random.default_rng(seed)
    ge = 0
    n, k = values.shape
    n_perm = max(1, int(n_perm))

    for _ in range(n_perm):
        permuted = np.empty_like(values)
        for i in range(n):
            permuted[i, :] = values[i, rng.permutation(k)]
        h = _friedman_statistic(permuted)
        if np.isfinite(h) and h >= obs:
            ge += 1
    p = (ge + 1.0) / (n_perm + 1.0)
    return float(obs), float(p)


def _paired_signed_rank_permutation(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = 4000,
    seed: int = 1,
) -> tuple[float, float, int, float]:
    d = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    d = d[np.isfinite(d)]
    d = d[d != 0]
    n = int(d.size)
    if n == 0:
        return float("nan"), 1.0, 0, 0.0

    ranks = pd.Series(np.abs(d)).rank(method="average").to_numpy(dtype=float)
    w_plus = float(ranks[d > 0].sum())
    total = float(ranks.sum())
    center = 0.5 * total
    obs_dev = abs(w_plus - center)

    rng = np.random.default_rng(seed)
    ge = 0
    n_perm = max(1, int(n_perm))
    for _ in range(n_perm):
        signs = rng.choice(np.array([-1.0, 1.0]), size=n, replace=True)
        w_perm = float(ranks[signs > 0].sum())
        if abs(w_perm - center) >= (obs_dev - 1e-12):
            ge += 1
    p = (ge + 1.0) / (n_perm + 1.0)
    return w_plus, float(p), n, float(np.median(d))


def _paired_sign_test_exact(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[float, float, int, int, int, int, float]:
    d_all = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    d_all = d_all[np.isfinite(d_all)]
    n_zero = int(np.sum(d_all == 0))
    d = d_all[d_all != 0]

    n = int(d.size)
    if n == 0:
        return float("nan"), 1.0, 0, 0, 0, n_zero, 0.0

    n_pos = int(np.sum(d > 0))
    n_neg = int(np.sum(d < 0))
    k = min(n_pos, n_neg)

    tail_count = 0
    for i in range(0, k + 1):
        tail_count += math.comb(n, i)
    p = float((2.0 * tail_count) / (2.0**n))
    p = float(min(1.0, max(0.0, p)))
    return float(n_pos), p, n, n_pos, n_neg, n_zero, float(np.median(d))


def _apply_pairwise_adjustments(
    pairs: list[dict], alpha: float, multiple_test_correction: bool
) -> tuple[list[dict], list[dict]]:
    p_raw_vals = [float(r.get("p_raw", np.nan)) for r in pairs]
    p_adj_bh_vals = _bh_adjust(p_raw_vals)
    p_used_vals = p_adj_bh_vals if multiple_test_correction else p_raw_vals
    for i, p_used in enumerate(p_used_vals):
        p_adj = p_adj_bh_vals[i]
        pairs[i]["p_adj_bh"] = p_adj
        pairs[i]["p_value_used"] = p_used
        pairs[i]["multiple_test_correction"] = bool(multiple_test_correction)
        pairs[i]["p_adjustment"] = "bh" if multiple_test_correction else "none"
        pairs[i]["significant_bh"] = bool(np.isfinite(p_adj) and p_adj < alpha)
        pairs[i]["significant_alpha"] = bool(np.isfinite(p_used) and p_used < alpha)
    sig_pairs = [r for r in pairs if bool(r["significant_alpha"])]
    return pairs, sig_pairs


def _compute_panel_stats_permutation(
    df: pd.DataFrame,
    method_order: list[str],
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
) -> dict:
    wide = df.pivot_table(index="sample_key", columns="method", values="distance", aggfunc="first")
    methods = [m for m in method_order if m in wide.columns]
    if not methods:
        return {
            "n_blocks": 0,
            "methods": [],
            "omnibus_stat": float("nan"),
            "omnibus_p": float("nan"),
            "multiple_test_correction": bool(multiple_test_correction),
            "stats_model": "paired_permutation",
            "pairwise": [],
            "pairwise_sig": [],
        }

    wide = wide[methods].dropna(axis=0, how="any")
    n_blocks = int(wide.shape[0])
    if n_blocks < 2 or len(methods) < 2:
        return {
            "n_blocks": n_blocks,
            "methods": methods,
            "omnibus_stat": float("nan"),
            "omnibus_p": float("nan"),
            "multiple_test_correction": bool(multiple_test_correction),
            "stats_model": "paired_permutation",
            "pairwise": [],
            "pairwise_sig": [],
        }

    values = wide.to_numpy(dtype=float)
    omnibus_stat, omnibus_p = _friedman_permutation_pvalue(values, n_perm=n_perm, seed=seed)

    pairs: list[dict] = []
    for pair_idx, (a, b) in enumerate(combinations(methods, 2), start=1):
        xa = wide[a].to_numpy(dtype=float)
        xb = wide[b].to_numpy(dtype=float)
        if pairwise_test == "sign_test_exact":
            stat, p_raw, n_nonzero, n_pos, n_neg, n_zero, median_diff = _paired_sign_test_exact(xa, xb)
            pairs.append(
                {
                    "method_a": a,
                    "method_b": b,
                    "pairwise_test": "sign_test_exact",
                    "statistic": stat,
                    "w_plus": float("nan"),
                    "p_raw": p_raw,
                    "n_nonzero_pairs": n_nonzero,
                    "n_pos_pairs": n_pos,
                    "n_neg_pairs": n_neg,
                    "n_zero_pairs": n_zero,
                    "median_diff_a_minus_b": median_diff,
                    "estimate": median_diff,
                    "std_error": float("nan"),
                    "model_family": "permutation",
                }
            )
        else:
            w_plus, p_raw, n_nonzero, median_diff = _paired_signed_rank_permutation(
                xa,
                xb,
                n_perm=n_perm,
                seed=seed + 1000 + pair_idx,
            )
            pairs.append(
                {
                    "method_a": a,
                    "method_b": b,
                    "pairwise_test": "signed_rank_permutation",
                    "statistic": w_plus,
                    "w_plus": w_plus,
                    "p_raw": p_raw,
                    "n_nonzero_pairs": n_nonzero,
                    "n_pos_pairs": float("nan"),
                    "n_neg_pairs": float("nan"),
                    "n_zero_pairs": float("nan"),
                    "median_diff_a_minus_b": median_diff,
                    "estimate": median_diff,
                    "std_error": float("nan"),
                    "model_family": "permutation",
                }
            )

    pairs, sig_pairs = _apply_pairwise_adjustments(
        pairs=pairs, alpha=alpha, multiple_test_correction=multiple_test_correction
    )
    return {
        "n_blocks": n_blocks,
        "methods": methods,
        "omnibus_stat": omnibus_stat,
        "omnibus_p": omnibus_p,
        "multiple_test_correction": bool(multiple_test_correction),
        "stats_model": "paired_permutation",
        "pairwise": pairs,
        "pairwise_sig": sig_pairs,
    }


def _run_beta_mixed_r(
    model_input: pd.DataFrame,
    method_order: list[str],
    emmeans_scale: str,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    r_code = textwrap.dedent(
        """
        suppressWarnings(suppressMessages({
          has_glmmTMB <- requireNamespace("glmmTMB", quietly = TRUE)
          has_emmeans <- requireNamespace("emmeans", quietly = TRUE)
        }))

        args <- commandArgs(trailingOnly = TRUE)
        in_tsv <- args[[1]]
        out_omni <- args[[2]]
        out_pair <- args[[3]]
        methods <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
        emmeans_scale <- args[[5]]

        d <- read.delim(in_tsv, sep = "\\t", stringsAsFactors = FALSE, check.names = FALSE)
        d$community_id <- factor(d$community_id)
        d$sample_key <- factor(d$sample_key)
        d$method <- factor(d$method, levels = methods)
        d <- d[!is.na(d$method), , drop = FALSE]
        d$distance <- as.numeric(d$distance)
        d <- d[is.finite(d$distance), , drop = FALSE]
        n_blocks <- length(unique(as.character(d$sample_key)))
        n_methods <- length(unique(as.character(d$method)))

        write_empty <- function(note, model_fit = FALSE) {
          write.table(data.frame(n_blocks=n_blocks,
                                 n_methods=n_methods,
                                 omnibus_stat=NA_real_,
                                 omnibus_df=NA_real_,
                                 omnibus_p=NA_real_,
                                 model_family="beta",
                                 model_fit=model_fit,
                                 note=note),
                      file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
          write.table(data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                                 se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character()),
                      file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
        }

        if (!has_glmmTMB || !has_emmeans) {
          missing <- c()
          if (!has_glmmTMB) missing <- c(missing, "glmmTMB")
          if (!has_emmeans) missing <- c(missing, "emmeans")
          write_empty(paste0("missing_R_packages:", paste(missing, collapse=",")), FALSE)
          quit(save="no", status=0)
        }

        if (nrow(d) < 2 || n_methods < 2) {
          write_empty("insufficient_data", FALSE)
          quit(save="no", status=0)
        }

        n_obs <- nrow(d)
        d$y_beta <- (d$distance * (n_obs - 1.0) + 0.5) / n_obs
        d$y_beta <- pmax(pmin(d$y_beta, 1.0 - 1e-6), 1e-6)

        fit <- tryCatch(
          glmmTMB::glmmTMB(y_beta ~ method + (1|community_id) + (1|sample_key),
                           family = glmmTMB::beta_family(link = "logit"),
                           data = d),
          error = function(e) e
        )
        fit0 <- tryCatch(
          glmmTMB::glmmTMB(y_beta ~ 1 + (1|community_id) + (1|sample_key),
                           family = glmmTMB::beta_family(link = "logit"),
                           data = d),
          error = function(e) e
        )

        if (inherits(fit, "error") || inherits(fit0, "error")) {
          msg <- if (inherits(fit, "error")) conditionMessage(fit) else conditionMessage(fit0)
          write_empty(msg, FALSE)
          quit(save="no", status=0)
        }

        fit_ok <- TRUE
        fit_note <- "ok"
        omnibus_stat <- NA_real_
        omnibus_df <- NA_real_
        omnibus_p <- NA_real_
        lrt <- tryCatch(stats::anova(fit0, fit, test = "Chisq"), error = function(e) e)
        if (inherits(lrt, "error")) {
          fit_ok <- FALSE
          fit_note <- conditionMessage(lrt)
        } else if (nrow(lrt) < 2) {
          fit_ok <- FALSE
          fit_note <- "anova_output_unexpected"
        } else {
          lrt_row <- lrt[2, , drop = FALSE]
          chisq_col <- intersect(c("Chisq", "LRT", "Deviance"), colnames(lrt_row))
          df_col <- intersect(c("Df", "df"), colnames(lrt_row))
          p_col <- grep("^Pr\\\\(", colnames(lrt_row), value = TRUE)
          if (length(chisq_col) > 0) omnibus_stat <- as.numeric(lrt_row[[chisq_col[1]]])
          if (length(df_col) > 0) omnibus_df <- as.numeric(lrt_row[[df_col[1]]])
          if (length(p_col) > 0) omnibus_p <- as.numeric(lrt_row[[p_col[1]]])
          if (!is.finite(omnibus_p)) {
            fit_note <- "anova_pvalue_missing"
          }
        }

        method_levels <- methods[methods %in% unique(as.character(d$method))]
        pair_df <- data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                              se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character(),
                              stringsAsFactors=FALSE)
        if (fit_ok && length(method_levels) >= 2) {
          em <- tryCatch(emmeans::emmeans(fit, specs = ~ method, type = emmeans_scale), error = function(e) e)
          if (inherits(em, "error")) {
            fit_ok <- FALSE
            fit_note <- conditionMessage(em)
          } else {
            pw <- tryCatch(
              as.data.frame(emmeans::contrast(em, method = "pairwise", adjust = "none")),
              error = function(e) e
            )
            if (inherits(pw, "error")) {
              fit_ok <- FALSE
              fit_note <- conditionMessage(pw)
            } else if (nrow(pw) > 0) {
              if (!("contrast" %in% colnames(pw))) {
                fit_ok <- FALSE
                fit_note <- "emmeans_pairs_missing_contrast"
              } else {
                pw$contrast <- as.character(pw$contrast)
                pw$method_a <- trimws(sub(" - .*", "", pw$contrast))
                pw$method_b <- trimws(sub(".* - ", "", pw$contrast))
                stat_col <- if ("z.ratio" %in% colnames(pw)) "z.ratio" else if ("t.ratio" %in% colnames(pw)) "t.ratio" else NA_character_
                se_col <- if ("SE" %in% colnames(pw)) "SE" else if ("se" %in% colnames(pw)) "se" else NA_character_
                est_col <- if ("estimate" %in% colnames(pw)) "estimate" else NA_character_
                p_col <- if ("p.value" %in% colnames(pw)) "p.value" else NA_character_
                pair_df <- data.frame(
                  method_a = pw$method_a,
                  method_b = pw$method_b,
                  estimate = if (!is.na(est_col)) as.numeric(pw[[est_col]]) else NA_real_,
                  se = if (!is.na(se_col)) as.numeric(pw[[se_col]]) else NA_real_,
                  statistic = if (!is.na(stat_col)) as.numeric(pw[[stat_col]]) else NA_real_,
                  p_raw = if (!is.na(p_col)) as.numeric(pw[[p_col]]) else NA_real_,
                  pairwise_test = "mixed_beta_glmmtmb_emmeans",
                  stringsAsFactors = FALSE
                )
                pair_df <- pair_df[
                  pair_df$method_a %in% method_levels & pair_df$method_b %in% method_levels,
                  ,
                  drop = FALSE
                ]
              }
            }
          }
        }

        write.table(data.frame(n_blocks=n_blocks,
                               n_methods=n_methods,
                               omnibus_stat=omnibus_stat,
                               omnibus_df=omnibus_df,
                               omnibus_p=omnibus_p,
                               model_family="beta",
                               model_fit=fit_ok,
                               note=fit_note),
                    file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
        if (!fit_ok) {
          pair_df <- data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                                se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character(),
                                stringsAsFactors=FALSE)
        }
        write.table(pair_df, file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
        """
    )

    with tempfile.TemporaryDirectory(prefix="tacticbench_fig11_beta_") as td:
        td_path = Path(td)
        in_tsv = td_path / "fig11_input.tsv"
        out_omni = td_path / "fig11_omnibus.tsv"
        out_pair = td_path / "fig11_pairwise.tsv"
        r_script = td_path / "fit_fig11_beta.R"
        model_input.to_csv(in_tsv, sep="\t", index=False)
        r_script.write_text(r_code, encoding="utf-8")
        methods_csv = ",".join([str(m) for m in method_order])
        cp = subprocess.run(
            [
                "Rscript",
                str(r_script),
                str(in_tsv),
                str(out_omni),
                str(out_pair),
                methods_csv,
                str(emmeans_scale),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if cp.returncode != 0:
            msg = cp.stderr.strip() or cp.stdout.strip() or f"Rscript exited with code {cp.returncode}"
            return pd.DataFrame(), pd.DataFrame(), msg
        if (not out_omni.exists()) or (not out_pair.exists()):
            return pd.DataFrame(), pd.DataFrame(), "Rscript did not write expected outputs."
        return pd.read_csv(out_omni, sep="\t"), pd.read_csv(out_pair, sep="\t"), ""


def _compute_panel_stats_beta(
    df: pd.DataFrame,
    method_order: list[str],
    alpha: float,
    multiple_test_correction: bool,
    emmeans_scale: str,
    log,
) -> dict | None:
    model_input = df[["community_id", "sample_key", "method", "distance"]].copy()
    model_input["community_id"] = model_input["community_id"].astype(str)
    model_input["sample_key"] = model_input["sample_key"].astype(str)
    model_input["method"] = model_input["method"].astype(str)
    model_input["distance"] = pd.to_numeric(model_input["distance"], errors="coerce")
    model_input = model_input[np.isfinite(model_input["distance"])].copy()

    if model_input.empty:
        return None

    methods = [m for m in method_order if m in set(model_input["method"].astype(str).tolist())]
    n_blocks = int(model_input["sample_key"].astype(str).nunique())
    if len(methods) < 2 or n_blocks < 2:
        return {
            "n_blocks": n_blocks,
            "methods": methods,
            "omnibus_stat": float("nan"),
            "omnibus_p": float("nan"),
            "multiple_test_correction": bool(multiple_test_correction),
            "stats_model": "mixed_effect_glmmtmb_beta",
            "pairwise": [],
            "pairwise_sig": [],
        }

    omni_df, pair_df, err = _run_beta_mixed_r(model_input=model_input, method_order=methods, emmeans_scale=emmeans_scale)
    if err:
        if log is not None:
            log.info(f"Figure 11 mixed beta stats failed: {err}")
        return None
    if omni_df.empty:
        return None

    omni = omni_df.iloc[0]
    model_fit = str(omni.get("model_fit", "False")).strip().lower()
    fit_ok = model_fit in {"true", "t", "1", "yes", "y"}
    if not fit_ok:
        note = str(omni.get("note", "model_fit_false"))
        if log is not None:
            log.info(f"Figure 11 mixed beta model not fitted: {note}")
        return None

    pairs: list[dict] = []
    if not pair_df.empty:
        for _, row in pair_df.iterrows():
            pairs.append(
                {
                    "method_a": str(row.get("method_a", "")),
                    "method_b": str(row.get("method_b", "")),
                    "pairwise_test": str(row.get("pairwise_test", "mixed_beta_glmmtmb_emmeans")),
                    "statistic": float(row.get("statistic", np.nan)),
                    "w_plus": float("nan"),
                    "p_raw": float(row.get("p_raw", np.nan)),
                    "n_nonzero_pairs": float("nan"),
                    "n_pos_pairs": float("nan"),
                    "n_neg_pairs": float("nan"),
                    "n_zero_pairs": float("nan"),
                    "median_diff_a_minus_b": float(row.get("estimate", np.nan)),
                    "estimate": float(row.get("estimate", np.nan)),
                    "std_error": float(row.get("se", np.nan)),
                    "model_family": str(omni.get("model_family", "beta")),
                }
            )

    pairs, sig_pairs = _apply_pairwise_adjustments(
        pairs=pairs, alpha=alpha, multiple_test_correction=multiple_test_correction
    )
    return {
        "n_blocks": int(omni.get("n_blocks", n_blocks)),
        "methods": methods,
        "omnibus_stat": float(omni.get("omnibus_stat", np.nan)),
        "omnibus_p": float(omni.get("omnibus_p", np.nan)),
        "multiple_test_correction": bool(multiple_test_correction),
        "stats_model": "mixed_effect_glmmtmb_beta",
        "pairwise": pairs,
        "pairwise_sig": sig_pairs,
    }


def _compute_metric_stats(
    df: pd.DataFrame,
    metric: str,
    method_order: list[str],
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
    stats_model: str,
    emmeans_scale: str,
    log,
) -> dict:
    if stats_model == "mixed_effect_glmmtmb_auto":
        if metric in AUTO_BETA_METRICS:
            mixed = _compute_panel_stats_beta(
                df=df,
                method_order=method_order,
                alpha=alpha,
                multiple_test_correction=multiple_test_correction,
                emmeans_scale=emmeans_scale,
                log=log,
            )
            if mixed is not None:
                return mixed
            if log is not None:
                log.info(f"Figure 11: beta mixed model failed for '{metric}', falling back to paired_permutation.")
        return _compute_panel_stats_permutation(
            df=df,
            method_order=method_order,
            n_perm=n_perm,
            seed=seed,
            alpha=alpha,
            pairwise_test=pairwise_test,
            multiple_test_correction=multiple_test_correction,
        )

    if stats_model == "mixed_effect_glmmtmb_beta":
        mixed = _compute_panel_stats_beta(
            df=df,
            method_order=method_order,
            alpha=alpha,
            multiple_test_correction=multiple_test_correction,
            emmeans_scale=emmeans_scale,
            log=log,
        )
        if mixed is not None:
            return mixed
        if log is not None:
            log.info(f"Figure 11: beta mixed model failed for '{metric}', falling back to paired_permutation.")

    return _compute_panel_stats_permutation(
        df=df,
        method_order=method_order,
        n_perm=n_perm,
        seed=seed,
        alpha=alpha,
        pairwise_test=pairwise_test,
        multiple_test_correction=multiple_test_correction,
    )


def _draw_pairwise_brackets(
    ax,
    pairwise_rows: list[dict],
    method_order: list[str],
    alpha: float = 0.05,
    label_fontsize: float = 6.5,
    star_fontsize: float = 7.2,
) -> None:
    if not pairwise_rows:
        return

    method_to_x = {m: i for i, m in enumerate(method_order, start=1)}
    intervals: list[tuple[int, int, float]] = []
    for row in pairwise_rows:
        a = str(row.get("method_a", ""))
        b = str(row.get("method_b", ""))
        p = float(row.get("p_value_used", row.get("p_adj_bh", row.get("p_raw", np.nan))))
        if (a not in method_to_x) or (b not in method_to_x) or (not np.isfinite(p)):
            continue
        x1 = method_to_x[a]
        x2 = method_to_x[b]
        if x1 == x2:
            continue
        if x1 > x2:
            x1, x2 = x2, x1
        intervals.append((x1, x2, p))
    if not intervals:
        return

    intervals.sort(key=lambda t: ((t[1] - t[0]), t[2]))
    levels: list[list[tuple[int, int]]] = []
    placed: list[tuple[int, int, int, float]] = []
    for x1, x2, p in intervals:
        placed_level = None
        for lv, spans in enumerate(levels):
            if all((x2 < a) or (x1 > b) for a, b in spans):
                spans.append((x1, x2))
                placed_level = lv
                break
        if placed_level is None:
            levels.append([(x1, x2)])
            placed_level = len(levels) - 1
        placed.append((x1, x2, placed_level, p))

    y_min, y_max = ax.get_ylim()
    y_range = float(y_max - y_min)
    if y_range <= 0:
        y_range = 1.0

    base = y_min + (0.58 * y_range)
    step = 0.08 * y_range
    cap_h = 0.014 * y_range
    text_off = 0.007 * y_range
    star_off = 0.015 * y_range

    for x1, x2, level, p in placed:
        y = min(base + (level * step), y_max - 0.12 * y_range)
        y_cap = min(y + cap_h, y_max - 0.11 * y_range)
        ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], color="#111111", linewidth=0.9, zorder=5)
        label_y = min(y_cap + text_off, y_max - 0.09 * y_range)
        star_y = min(label_y + star_off, y_max - 0.06 * y_range)
        ax.text(
            (x1 + x2) / 2.0,
            label_y,
            f"p={p:.3g}",
            ha="center",
            va="bottom",
            fontsize=label_fontsize,
            color="#111111",
            zorder=6,
        )
        ax.text(
            (x1 + x2) / 2.0,
            star_y,
            _pvalue_stars(p, alpha=alpha),
            ha="center",
            va="bottom",
            fontsize=star_fontsize,
            fontweight="bold",
            color="#111111",
            zorder=6,
        )


def run(cfg) -> int:
    """
    Figure 11 (SVG):
      One violin/box plot per distance metric.
      - x-axis: analysis methods
      - y-axis: distance values across all included samples
      - optional repeated-measures statistics per metric
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig11", []) or [])

    long = _load_fig11_distances(cfg, log)
    if long.empty:
        log.info("Figure 11: no distance rows to plot.")
        return 0

    if exclude_patterns:
        keep = ~long["community_id"].astype(str).map(lambda cid: _is_excluded(cid, exclude_patterns))
        long = long[keep].copy()
        if long.empty:
            log.info("Figure 11: no rows left after figure-level community exclusions.")
            return 0

    configured_methods = list(getattr(cfg.distances, "methods", METHOD_ORDER) or METHOD_ORDER)
    method_order = [m for m in METHOD_ORDER if m in configured_methods] + [m for m in configured_methods if m not in METHOD_ORDER]
    long = long[long["method"].isin(method_order)].copy()
    if long.empty:
        log.info("Figure 11: no rows for configured methods.")
        return 0

    fig11_metrics_override = list(getattr(cfg.figures, "fig11_metrics", []) or [])
    configured_metrics = (
        fig11_metrics_override
        if fig11_metrics_override
        else list(getattr(cfg.distances, "metrics", list(METRIC_LABEL.keys())) or list(METRIC_LABEL.keys()))
    )
    metrics = [
        m
        for m in configured_metrics
        if m in long["metric"].astype(str).unique().tolist() and m in METRIC_LABEL
    ]
    if not metrics:
        log.info("Figure 11: no metrics found in distances table.")
        return 0

    long["sample_key"] = long["community_id"].astype(str) + "|" + long["sample_id"].astype(str)

    out_dir = figures_root(outdir) / "fig11_distance_distributions"
    out_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(out_dir / "distance_distribution_points.tsv", sep="\t", index=False)

    title = getattr(cfg.figures, "fig11_title", None) or "Distance Distributions by Method"
    panel_title_by_metric = getattr(cfg.figures, "fig11_panel_title_by_metric", {}) or {}
    combined_2x2_title = getattr(cfg.figures, "fig11_combined_2x2_title", None)
    combined_3x2_title = getattr(cfg.figures, "fig11_combined_3x2_title", None)
    plot_type = _normalize_plot_type(getattr(cfg.figures, "fig11_plot_type", "violin"), default="violin")
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))

    stats_enabled = bool(getattr(cfg.figures, "fig11_stats_enabled", True))
    stats_permutations = max(100, int(getattr(cfg.figures, "fig11_stats_permutations", 4000)))
    stats_seed = int(getattr(cfg.figures, "fig11_stats_seed", 1))
    stats_alpha = float(getattr(cfg.figures, "fig11_stats_alpha", 0.05))
    stats_show_pairwise = bool(getattr(cfg.figures, "fig11_stats_show_pairwise", True))
    stats_model = _normalize_stats_model(getattr(cfg.figures, "fig11_stats_model", "mixed_effect_glmmtmb_auto"))
    stats_pairwise_test = _normalize_pairwise_test(getattr(cfg.figures, "fig11_stats_pairwise_test", "sign_test_exact"))
    stats_multiple_test_correction = bool(getattr(cfg.figures, "fig11_stats_multiple_test_correction", True))
    stats_glmm_emmeans_scale = _normalize_emmeans_scale(getattr(cfg.figures, "fig11_stats_glmm_emmeans_scale", "link"))
    fig11_font_size = max(6.0, float(getattr(cfg.figures, "fig11_font_size", 10.0)))
    fig11_tick_size = max(6.0, fig11_font_size - 2.0)
    fig11_suptitle_size = fig11_font_size + 2.0
    fig11_stats_text_size = max(6.0, fig11_font_size - 3.0)
    fig11_combined_stats_text_size = max(6.0, fig11_stats_text_size - 0.6)
    fig11_pairwise_label_size = max(6.0, fig11_font_size - 3.5)
    fig11_pairwise_star_size = max(6.0, fig11_font_size - 2.8)
    fig11_combined_note_size = max(6.0, fig11_font_size - 2.0)

    def _metric_panel_title(metric: str) -> str:
        default = METRIC_LABEL.get(metric, metric)
        if isinstance(panel_title_by_metric, dict):
            for key in (metric, str(metric).lower(), str(metric).upper()):
                v = panel_title_by_metric.get(key)
                if v is not None:
                    txt = str(v).strip()
                    if txt:
                        return txt
        return default

    metric_stats: dict[str, dict] = {}
    omnibus_rows: list[dict] = []
    pairwise_rows: list[dict] = []

    if stats_enabled:
        for i, metric in enumerate(metrics, start=1):
            sub = long[long["metric"] == metric].copy()
            st = _compute_metric_stats(
                df=sub,
                metric=metric,
                method_order=method_order,
                n_perm=stats_permutations,
                seed=stats_seed + (i * 10000),
                alpha=stats_alpha,
                pairwise_test=stats_pairwise_test,
                multiple_test_correction=stats_multiple_test_correction,
                stats_model=stats_model,
                emmeans_scale=stats_glmm_emmeans_scale,
                log=log,
            )
            metric_stats[metric] = st
            omnibus_rows.append(
                {
                    "metric": metric,
                    "label": METRIC_LABEL.get(metric, metric),
                    "stats_model": st.get("stats_model", stats_model),
                    "pairwise_test": stats_pairwise_test,
                    "multiple_test_correction": bool(st.get("multiple_test_correction", stats_multiple_test_correction)),
                    "n_blocks": st.get("n_blocks", 0),
                    "omnibus_stat": st.get("omnibus_stat", np.nan),
                    "omnibus_p": st.get("omnibus_p", np.nan),
                    "omnibus_significant_alpha": bool(
                        np.isfinite(float(st.get("omnibus_p", np.nan))) and float(st.get("omnibus_p", np.nan)) < stats_alpha
                    ),
                    "alpha": stats_alpha,
                    "n_permutations": stats_permutations,
                }
            )
            for row in st.get("pairwise", []):
                pairwise_rows.append(
                    {
                        "metric": metric,
                        "label": METRIC_LABEL.get(metric, metric),
                        "method_a": row.get("method_a", ""),
                        "method_b": row.get("method_b", ""),
                        "stats_model": st.get("stats_model", stats_model),
                        "pairwise_test": row.get("pairwise_test", stats_pairwise_test),
                        "statistic": row.get("statistic", row.get("w_plus", np.nan)),
                        "w_plus": row.get("w_plus", np.nan),
                        "n_nonzero_pairs": row.get("n_nonzero_pairs", np.nan),
                        "n_pos_pairs": row.get("n_pos_pairs", np.nan),
                        "n_neg_pairs": row.get("n_neg_pairs", np.nan),
                        "n_zero_pairs": row.get("n_zero_pairs", np.nan),
                        "median_diff_a_minus_b": row.get("median_diff_a_minus_b", np.nan),
                        "estimate": row.get("estimate", row.get("median_diff_a_minus_b", np.nan)),
                        "std_error": row.get("std_error", np.nan),
                        "model_family": row.get("model_family", ""),
                        "p_raw": row.get("p_raw", np.nan),
                        "p_adj_bh": row.get("p_adj_bh", np.nan),
                        "p_value_used": row.get("p_value_used", row.get("p_adj_bh", np.nan)),
                        "p_adjustment": row.get("p_adjustment", "bh"),
                        "multiple_test_correction": bool(
                            row.get("multiple_test_correction", stats_multiple_test_correction)
                        ),
                        "significant_alpha": row.get("significant_alpha", False),
                        "significant_bh_alpha": row.get("significant_bh", False),
                        "alpha": stats_alpha,
                        "n_permutations": stats_permutations,
                    }
                )

        pd.DataFrame(omnibus_rows).to_csv(out_dir / "distance_distribution_stats_omnibus.tsv", sep="\t", index=False)
        pd.DataFrame(pairwise_rows).to_csv(out_dir / "distance_distribution_stats_pairwise.tsv", sep="\t", index=False)

    def _draw_metric_panel(
        ax,
        metric: str,
        panel_title: str,
        title_fontsize: float,
        tick_fontsize: float,
        text_fontsize: float,
        show_threshold_note: bool,
    ) -> bool:
        sub = long[long["metric"] == metric].copy()
        if sub.empty:
            return False
        positions = np.arange(1, len(method_order) + 1, dtype=float)
        method_data: list[tuple[int, str, np.ndarray]] = []

        for i, method in enumerate(method_order):
            vals = sub[sub["method"] == method]["distance"].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size > 0:
                method_data.append((i, method, vals))

        if plot_type == "boxplot":
            for i, method, vals in method_data:
                ax.boxplot(
                    [vals],
                    positions=[positions[i]],
                    widths=0.55,
                    patch_artist=True,
                    boxprops=dict(facecolor=METHOD_COLORS.get(method, "#999999"), alpha=0.40, edgecolor="#333333"),
                    medianprops=dict(color="#111111", linewidth=1.3),
                    whiskerprops=dict(color="#333333", linewidth=0.9),
                    capprops=dict(color="#333333", linewidth=0.9),
                )
        else:
            violin_vals: list[np.ndarray] = []
            violin_pos: list[float] = []
            violin_methods: list[str] = []
            for i, method, vals in method_data:
                if vals.size >= 2:
                    violin_vals.append(vals)
                    violin_pos.append(float(positions[i]))
                    violin_methods.append(method)
            if violin_vals:
                vp = ax.violinplot(
                    violin_vals,
                    positions=violin_pos,
                    widths=0.78,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                )
                for body, method in zip(vp["bodies"], violin_methods):
                    body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                    body.set_edgecolor("#333333")
                    body.set_alpha(0.5)
                if "cmedians" in vp:
                    vp["cmedians"].set_color("#111111")
                    vp["cmedians"].set_linewidth(1.2)

        rng = np.random.default_rng(1)
        for i, method, vals in method_data:
            xj = positions[i] + rng.normal(0.0, 0.04, size=vals.size)
            ax.scatter(
                xj,
                vals,
                s=point_size,
                alpha=0.40,
                color=METHOD_COLORS.get(method, "#666666"),
                edgecolors="none",
                zorder=3,
            )

        ax.set_xticks(positions)
        ax.set_xticklabels([display_method(m) for m in method_order], rotation=0, fontsize=tick_fontsize)
        ax.set_xlabel("Method", fontsize=fig11_font_size)
        ax.set_ylabel("Distance", fontsize=fig11_font_size)
        ax.tick_params(axis="y", labelsize=tick_fontsize)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        ax.set_title(panel_title, fontsize=title_fontsize)

        if stats_enabled:
            st = metric_stats.get(metric, None)
            if isinstance(st, dict):
                if stats_show_pairwise:
                    _draw_pairwise_brackets(
                        ax=ax,
                        pairwise_rows=list(st.get("pairwise_sig", [])),
                        method_order=method_order,
                        alpha=stats_alpha,
                        label_fontsize=fig11_pairwise_label_size,
                        star_fontsize=fig11_pairwise_star_size,
                    )
                    thresh_text = (
                        f"Pairwise BH q-threshold: alpha={stats_alpha:.3g}\n*** q<0.001, ** q<0.01, * q<{stats_alpha:.3g}"
                        if stats_multiple_test_correction
                        else f"Pairwise p-threshold: alpha={stats_alpha:.3g}\n*** p<0.001, ** p<0.01, * p<{stats_alpha:.3g}"
                    )
                    if show_threshold_note:
                        ax.text(
                            0.02,
                            0.03,
                            thresh_text,
                            transform=ax.transAxes,
                            ha="left",
                            va="bottom",
                            fontsize=text_fontsize,
                            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
                        )

                op = float(st.get("omnibus_p", np.nan))
                os = float(st.get("omnibus_stat", np.nan))
                model_used = str(st.get("stats_model", stats_model))
                ax.text(
                    0.98,
                    0.03,
                    f"{model_used}\nomnibus={os:.3g}, p={op:.3g}" if np.isfinite(op) else f"{model_used}\nomnibus=NA, p=NA",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=text_fontsize,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
                )
        return True

    for metric in metrics:
        fig = plt.figure(figsize=(8.8, 5.2))
        ax = fig.add_subplot(1, 1, 1)
        ok = _draw_metric_panel(
            ax=ax,
            metric=metric,
            panel_title=f"{title}: {_metric_panel_title(metric)}",
            title_fontsize=fig11_font_size,
            tick_fontsize=fig11_tick_size,
            text_fontsize=fig11_stats_text_size,
            show_threshold_note=True,
        )
        if not ok:
            plt.close(fig)
            continue

        suffix = "violin" if plot_type == "violin" else "boxplot"
        out_svg = out_dir / f"{metric}_{suffix}.svg"
        fig.tight_layout()
        fig.savefig(out_svg, format="svg")
        plt.close(fig)

    # Additional combined Figure 11 panel set (fixed 2x2 layout):
    # Jaccard, Weighted Jaccard, Bray-Curtis, and one gUniFrac panel.
    gunifrac_metric = None
    if "gunifrac_a1.0" in metrics:
        gunifrac_metric = "gunifrac_a1.0"
    elif "gunifrac_a0.5" in metrics:
        gunifrac_metric = "gunifrac_a0.5"
    elif "gunifrac_a0.0" in metrics:
        gunifrac_metric = "gunifrac_a0.0"

    panel_order = ["jaccard", "weighted_jaccard", "braycurtis"]
    if gunifrac_metric is not None:
        panel_order.append(gunifrac_metric)

    panel_metrics = [m for m in panel_order if m in metrics]
    if len(panel_metrics) < 4:
        log.info(
            "Figure 11: combined 2x2 requested metrics are partially missing; "
            f"rendering {len(panel_metrics)}/4 available panels."
        )

    combined = plt.figure(figsize=(13.2, 10.4))
    for idx in range(4):
        ax = combined.add_subplot(2, 2, idx + 1)
        if idx >= len(panel_metrics):
            ax.axis("off")
            continue
        metric = panel_metrics[idx]
        _draw_metric_panel(
            ax=ax,
            metric=metric,
            panel_title=_metric_panel_title(metric),
            title_fontsize=fig11_font_size,
            tick_fontsize=fig11_tick_size,
            text_fontsize=fig11_combined_stats_text_size,
            show_threshold_note=False,
        )

    if stats_enabled and stats_show_pairwise:
        combined_thresh_text = (
            f"Pairwise BH q-threshold: alpha={stats_alpha:.3g}\n*** q<0.001, ** q<0.01, * q<{stats_alpha:.3g}"
            if stats_multiple_test_correction
            else f"Pairwise p-threshold: alpha={stats_alpha:.3g}\n*** p<0.001, ** p<0.01, * p<{stats_alpha:.3g}"
        )
        combined.text(
            0.01,
            0.01,
            combined_thresh_text,
            ha="left",
            va="bottom",
            fontsize=fig11_combined_note_size,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
        )
    combined.suptitle(
        str(combined_2x2_title).strip() if combined_2x2_title is not None else str(title),
        fontsize=fig11_suptitle_size,
    )
    combined.tight_layout(rect=[0, 0.03, 1, 0.96])
    suffix = "violin" if plot_type == "violin" else "boxplot"
    combined_out_svg = out_dir / f"distance_distributions_combined_2x2_{suffix}.svg"
    combined.savefig(combined_out_svg, format="svg")
    plt.close(combined)

    # Additional combined Figure 11 panel set (fixed 3x2 layout):
    # includes up to first 6 configured metrics.
    panel_metrics_3x2 = metrics[:6]
    if len(metrics) > 6:
        log.info("Figure 11: more than 6 metrics configured; combined 3x2 figure includes first 6 metrics only.")

    combined_3x2 = plt.figure(figsize=(16.8, 15.0))
    for idx in range(6):
        ax = combined_3x2.add_subplot(3, 2, idx + 1)
        if idx >= len(panel_metrics_3x2):
            ax.axis("off")
            continue
        metric = panel_metrics_3x2[idx]
        _draw_metric_panel(
            ax=ax,
            metric=metric,
            panel_title=_metric_panel_title(metric),
            title_fontsize=fig11_font_size,
            tick_fontsize=fig11_tick_size,
            text_fontsize=fig11_combined_stats_text_size,
            show_threshold_note=False,
        )

    if stats_enabled and stats_show_pairwise:
        combined_thresh_text = (
            f"Pairwise BH q-threshold: alpha={stats_alpha:.3g}\n*** q<0.001, ** q<0.01, * q<{stats_alpha:.3g}"
            if stats_multiple_test_correction
            else f"Pairwise p-threshold: alpha={stats_alpha:.3g}\n*** p<0.001, ** p<0.01, * p<{stats_alpha:.3g}"
        )
        combined_3x2.text(
            0.01,
            0.01,
            combined_thresh_text,
            ha="left",
            va="bottom",
            fontsize=fig11_combined_note_size,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
        )
    combined_3x2.suptitle(
        str(combined_3x2_title).strip() if combined_3x2_title is not None else str(title),
        fontsize=fig11_suptitle_size,
    )
    combined_3x2.tight_layout(rect=[0, 0.03, 1, 0.96])
    combined_3x2_out_svg = out_dir / f"distance_distributions_combined_3x2_{suffix}.svg"
    combined_3x2.savefig(combined_3x2_out_svg, format="svg")
    plt.close(combined_3x2)

    log.info(f"Wrote Figure 11 to {out_dir}")
    return 0
