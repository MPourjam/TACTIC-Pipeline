from __future__ import annotations

from fnmatch import fnmatch
from itertools import combinations
from pathlib import Path
import subprocess
import tempfile
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from ..utils.labels import display_methods
from .plot_fig11 import _normalize_pairwise_test, _paired_sign_test_exact

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
CLUSTERING = {"OTU", "TIC", "TAC"}
METHOD_KEY_ALIASES = {
    "ASV": ["asv", "dada2"],
    "zOTU": ["zotu", "unoise3"],
    "OTU": ["otu", "uparse"],
    "TIC": ["tic"],
    "TAC": ["tac"],
}
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}
METRICS = [
    ("observed_species", "Observed Species"),
    ("shannon", "Shannon"),
    ("effective_shannon", "Effective Shannon"),
    ("simpson", "Simpson (1 - sum p^2)"),
    ("effective_simpson", "Effective Simpson (1 / sum p^2)"),
]
METRICS_RICHNESS_EFFECTIVE = [
    ("observed_species", "Observed Species"),
    ("effective_shannon", "Effective Shannon"),
    ("effective_simpson", "Effective Simpson (1 / sum p^2)"),
]


def _normalize_plot_type(raw: str, default: str = "boxplot") -> str:
    v = str(raw or "").strip().lower()
    if v in {"box", "boxplot"}:
        return "boxplot"
    if v in {"violin", "violinplot"}:
        return "violin"
    return default


def _normalize_stats_model(raw: str, default: str = "mixed_effect_glmmtmb") -> str:
    v = str(raw or "").strip().lower()
    if v in {"mixed_effect_glmmtmb", "paired_permutation"}:
        return v
    return default


def _normalize_emmeans_scale(raw: str, default: str = "response") -> str:
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


def _resolve_target_level(cfg, method: str) -> str:
    """
    Resolve expected-composition level used for Figure 07 per method.
      - ASV/zOTU: tip level
      - OTU/TIC/TAC: species level
    """
    default_level = "species" if method in CLUSTERING else "tip"
    raw = getattr(cfg.figures, "fig07_target_level_by_method", {})
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


def _paired_signed_rank_permutation(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = 2000,
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


def _draw_pairwise_brackets(
    ax,
    pairwise_rows: list[dict],
    method_order: list[str],
    alpha: float = 0.05,
    placement: str = "top",
    spacing_scale: float = 1.2,
) -> None:
    if not pairwise_rows:
        return
    place_bottom = str(placement).strip().lower() == "bottom"
    spacing_scale = max(0.8, float(spacing_scale))

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
    n_levels = max(1, len(levels))

    y_min, y_max = ax.get_ylim()
    y_range = float(y_max - y_min)
    if y_range <= 0:
        y_range = max(1.0, abs(y_max))

    if place_bottom:
        base = y_min - (0.06 * y_range)
    else:
        base = y_max + (0.08 * y_range)
    step = 0.12 * y_range * spacing_scale
    cap_h = 0.026 * y_range
    star_off = 0.030 * y_range

    bracket_lw = 1.55
    star_fs = 10.0

    top = y_max
    bottom = y_min
    for x1, x2, level, p in placed:
        y = base + (level * step)
        if place_bottom:
            y_cap = y - cap_h
            ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], color="#111111", linewidth=bracket_lw, zorder=5)
            star_y = y_cap - star_off
            ax.text(
                (x1 + x2) / 2.0,
                star_y,
                _pvalue_stars(p, alpha=alpha),
                ha="center",
                va="top",
                fontsize=star_fs,
                fontweight="bold",
                color="#111111",
                zorder=6,
            )
            bottom = min(bottom, star_y)
        else:
            y_cap = y + cap_h
            ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], color="#111111", linewidth=bracket_lw, zorder=5)
            star_y = y_cap + star_off
            ax.text(
                (x1 + x2) / 2.0,
                star_y,
                _pvalue_stars(p, alpha=alpha),
                ha="center",
                va="bottom",
                fontsize=star_fs,
                fontweight="bold",
                color="#111111",
                zorder=6,
            )
            top = max(top, star_y)

    if place_bottom:
        ax.set_ylim(bottom - (0.05 * y_range), y_max)
    else:
        ax.set_ylim(y_min, top + (0.08 * y_range))


def _run_gaussian_mixed_r(
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
        d$value <- as.numeric(d$value)
        d <- d[is.finite(d$value), , drop = FALSE]

        n_obs <- nrow(d)
        n_blocks <- length(unique(as.character(d$sample_key)))
        n_methods <- length(unique(as.character(d$method)))

        fit_ok <- FALSE
        fit_note <- "ok"
        omnibus_stat <- NA_real_
        omnibus_df <- NA_real_
        omnibus_p <- NA_real_
        pair_df <- data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                              se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character(),
                              stringsAsFactors=FALSE)

        if (!has_glmmTMB || !has_emmeans) {
          missing <- character()
          if (!has_glmmTMB) missing <- c(missing, "glmmTMB")
          if (!has_emmeans) missing <- c(missing, "emmeans")
          fit_note <- paste0("missing_r_packages:", paste(missing, collapse=","))
        } else if (n_obs < 2 || n_methods < 2 || n_blocks < 2) {
          fit_note <- "insufficient_data"
        } else {
          fit <- tryCatch(
            glmmTMB::glmmTMB(value ~ method + (1|community_id) + (1|sample_key),
                             family = gaussian(link = "identity"),
                             data = d),
            error = function(e) e
          )
          fit0 <- tryCatch(
            glmmTMB::glmmTMB(value ~ 1 + (1|community_id) + (1|sample_key),
                             family = gaussian(link = "identity"),
                             data = d),
            error = function(e) e
          )
          if (inherits(fit, "error")) {
            fit_note <- paste0("fit_error:", conditionMessage(fit))
          } else if (inherits(fit0, "error")) {
            fit_note <- paste0("null_fit_error:", conditionMessage(fit0))
          } else {
            fit_ok <- TRUE
            an <- tryCatch(as.data.frame(anova(fit, fit0, test = "Chisq")), error = function(e) e)
            if (!inherits(an, "error") && nrow(an) >= 2) {
              cn <- colnames(an)
              stat_col <- if ("Chisq" %in% cn) "Chisq" else if ("LRT" %in% cn) "LRT" else NA_character_
              p_col <- if ("Pr(>Chisq)" %in% cn) "Pr(>Chisq)" else if ("p.value" %in% cn) "p.value" else NA_character_
              df_col <- if ("Df" %in% cn) "Df" else if ("df" %in% cn) "df" else NA_character_
              if (!is.na(stat_col)) omnibus_stat <- suppressWarnings(as.numeric(an[2, stat_col]))
              if (!is.na(df_col)) omnibus_df <- suppressWarnings(as.numeric(an[2, df_col]))
              if (!is.na(p_col)) omnibus_p <- suppressWarnings(as.numeric(an[2, p_col]))
            } else {
              fit_note <- "anova_failed"
            }

            em <- tryCatch(emmeans::emmeans(fit, specs = ~ method, type = emmeans_scale), error = function(e) e)
            if (!inherits(em, "error")) {
              pw <- tryCatch(
                as.data.frame(emmeans::contrast(em, method = "pairwise", adjust = "none")),
                error = function(e) e
              )
              if (!inherits(pw, "error") && nrow(pw) > 0 && ("contrast" %in% colnames(pw))) {
                parts <- strsplit(as.character(pw$contrast), " - ", fixed = TRUE)
                method_a <- vapply(parts, function(x) if (length(x) >= 1) x[[1]] else NA_character_, "")
                method_b <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, "")
                stat_col <- if ("z.ratio" %in% colnames(pw)) "z.ratio" else if ("t.ratio" %in% colnames(pw)) "t.ratio" else NA_character_
                se_col <- if ("SE" %in% colnames(pw)) "SE" else if ("se" %in% colnames(pw)) "se" else NA_character_
                est_col <- if ("estimate" %in% colnames(pw)) "estimate" else NA_character_
                p_col <- if ("p.value" %in% colnames(pw)) "p.value" else NA_character_
                pair_df <- data.frame(
                  method_a = method_a,
                  method_b = method_b,
                  estimate = if (!is.na(est_col)) as.numeric(pw[[est_col]]) else NA_real_,
                  se = if (!is.na(se_col)) as.numeric(pw[[se_col]]) else NA_real_,
                  statistic = if (!is.na(stat_col)) as.numeric(pw[[stat_col]]) else NA_real_,
                  p_raw = if (!is.na(p_col)) as.numeric(pw[[p_col]]) else NA_real_,
                  pairwise_test = "mixed_glmmtmb_emmeans",
                  stringsAsFactors = FALSE
                )
                pair_df <- pair_df[
                  pair_df$method_a %in% methods & pair_df$method_b %in% methods,
                  ,
                  drop = FALSE
                ]
              }
            } else {
              fit_note <- "emmeans_failed"
            }
          }
        }

        write.table(data.frame(n_blocks=n_blocks,
                               n_methods=n_methods,
                               omnibus_stat=omnibus_stat,
                               omnibus_df=omnibus_df,
                               omnibus_p=omnibus_p,
                               model_family="gaussian",
                               model_fit=fit_ok,
                               note=fit_note),
                    file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
        write.table(pair_df, file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
        """
    )

    with tempfile.TemporaryDirectory(prefix="tacticbench_fig07_glmmtmb_") as td:
        td_path = Path(td)
        in_tsv = td_path / "fig07_input.tsv"
        out_omni = td_path / "fig07_omnibus.tsv"
        out_pair = td_path / "fig07_pairwise.tsv"
        r_script = td_path / "fit_fig07_glmmtmb.R"
        model_input.to_csv(in_tsv, sep="\t", index=False)
        r_script.write_text(r_code, encoding="utf-8")
        methods_csv = ",".join([str(m) for m in method_order])
        cp = subprocess.run(
            ["Rscript", str(r_script), str(in_tsv), str(out_omni), str(out_pair), methods_csv, emmeans_scale],
            capture_output=True,
            text=True,
            check=False,
        )
        if cp.returncode != 0:
            msg = cp.stderr.strip() or cp.stdout.strip() or f"Rscript exited with code {cp.returncode}"
            return pd.DataFrame(), pd.DataFrame(), msg
        if (not out_omni.exists()) or (not out_pair.exists()):
            return pd.DataFrame(), pd.DataFrame(), "Rscript did not write expected outputs."
        omni_df = pd.read_csv(out_omni, sep="\t")
        pair_df = pd.read_csv(out_pair, sep="\t")
        return omni_df, pair_df, ""


def _compute_metric_stats(
    panel_df: pd.DataFrame,
    method_order: list[str],
    n_perm: int,
    seed: int,
    stats_model: str,
    pairwise_test: str,
    emmeans_scale: str,
    alpha: float,
    multiple_test_correction: bool,
    log,
) -> dict:
    d = panel_df.copy()
    d["community_id"] = d["community_id"].astype(str)
    d["sample_key"] = d["sample_key"].astype(str)
    d["method"] = d["method"].astype(str)
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d = d[np.isfinite(d["value"])].copy()

    groups = []
    for method in method_order:
        vals = d.loc[d["method"] == method, "value"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size > 0:
            groups.append(vals)

    pairwise_test = _normalize_pairwise_test(pairwise_test, default="signed_rank_permutation")

    if stats_model == "mixed_effect_glmmtmb":
        methods = [m for m in method_order if m in set(d["method"].tolist())]
        if len(methods) >= 2 and int(d["sample_key"].nunique()) >= 2:
            model_input = d[["community_id", "sample_key", "method", "value"]].copy()
            omni_df, pair_df, err = _run_gaussian_mixed_r(
                model_input=model_input, method_order=methods, emmeans_scale=emmeans_scale
            )
            if err:
                if log is not None:
                    log.info(f"Figure 07 GLMM failed; falling back to permutation stats: {err}")
            elif not omni_df.empty:
                omni = omni_df.iloc[0]
                fit_ok = bool(str(omni.get("model_fit", "False")).lower() in {"true", "1", "yes"})
                if fit_ok:
                    pairs: list[dict] = []
                    if not pair_df.empty:
                        for _, row in pair_df.iterrows():
                            pairs.append(
                                {
                                    "method_a": str(row.get("method_a", "")),
                                    "method_b": str(row.get("method_b", "")),
                                    "estimate": float(row.get("estimate", np.nan)),
                                    "std_error": float(row.get("se", np.nan)),
                                    "statistic": float(row.get("statistic", np.nan)),
                                    "p_raw": float(row.get("p_raw", np.nan)),
                                    "pairwise_test": str(row.get("pairwise_test", "mixed_glmmtmb_emmeans")),
                                }
                            )
                        pvals = [float(x.get("p_raw", np.nan)) for x in pairs]
                        pads = _bh_adjust(pvals)
                        for i in range(len(pairs)):
                            pairs[i]["p_adj_bh"] = float(pads[i])
                            p_used = float(pairs[i]["p_adj_bh"] if multiple_test_correction else pairs[i]["p_raw"])
                            pairs[i]["p_value_used"] = p_used
                            pairs[i]["significant_alpha"] = bool(np.isfinite(p_used) and (p_used < alpha))
                    return {
                        "stats_model": "mixed_effect_glmmtmb",
                        "omnibus_stat": float(omni.get("omnibus_stat", np.nan)),
                        "omnibus_df": float(omni.get("omnibus_df", np.nan)),
                        "omnibus_p": float(omni.get("omnibus_p", np.nan)),
                        "model_family": str(omni.get("model_family", "gaussian")),
                        "note": str(omni.get("note", "")),
                        "n_blocks": int(omni.get("n_blocks", d["sample_key"].nunique())),
                        "pairwise_test": "mixed_glmmtmb_emmeans",
                        "pairwise": pairs,
                    }
                if log is not None:
                    log.info(
                        f"Figure 07 GLMM fit was not usable (note={str(omni.get('note', ''))}); falling back to permutation stats."
                    )

    h, p = _permutation_kw_pvalue(groups, n_perm=n_perm, seed=seed)
    pairs: list[dict] = []
    wide = d.pivot_table(index="sample_key", columns="method", values="value", aggfunc="first")
    methods = [m for m in method_order if m in wide.columns]
    wide = wide[methods].dropna(axis=0, how="any") if methods else wide
    if len(methods) >= 2 and wide.shape[0] >= 2:
        for pair_idx, (a, b) in enumerate(combinations(methods, 2), start=1):
            xa = wide[a].to_numpy(dtype=float)
            xb = wide[b].to_numpy(dtype=float)
            if pairwise_test == "sign_test_exact":
                stat, p_raw, n_nonzero, n_pos, n_neg, n_zero, median_diff = _paired_sign_test_exact(xa, xb)
                pairs.append(
                    {
                        "method_a": a,
                        "method_b": b,
                        "estimate": float(median_diff),
                        "std_error": float("nan"),
                        "statistic": float(stat),
                        "w_plus": float("nan"),
                        "p_raw": float(p_raw),
                        "pairwise_test": "sign_test_exact",
                        "n_nonzero_pairs": int(n_nonzero),
                        "n_pos_pairs": int(n_pos),
                        "n_neg_pairs": int(n_neg),
                        "n_zero_pairs": int(n_zero),
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
                        "estimate": float(median_diff),
                        "std_error": float("nan"),
                        "statistic": float(w_plus),
                        "w_plus": float(w_plus),
                        "p_raw": float(p_raw),
                        "pairwise_test": "signed_rank_permutation",
                        "n_nonzero_pairs": int(n_nonzero),
                        "n_pos_pairs": float("nan"),
                        "n_neg_pairs": float("nan"),
                        "n_zero_pairs": float("nan"),
                    }
                )
        pvals = [float(x.get("p_raw", np.nan)) for x in pairs]
        pads = _bh_adjust(pvals)
        for i in range(len(pairs)):
            pairs[i]["p_adj_bh"] = float(pads[i])
            p_used = float(pairs[i]["p_adj_bh"] if multiple_test_correction else pairs[i]["p_raw"])
            pairs[i]["p_value_used"] = p_used
            pairs[i]["significant_alpha"] = bool(np.isfinite(p_used) and (p_used < alpha))

    return {
        "stats_model": "paired_permutation",
        "omnibus_stat": float(h),
        "omnibus_df": float("nan"),
        "omnibus_p": float(p),
        "model_family": "",
        "note": "",
        "n_blocks": int(d["sample_key"].nunique()),
        "pairwise_test": pairwise_test,
        "pairwise": pairs,
    }



def _kruskal_h_stat(groups: list[np.ndarray]) -> float:
    """
    Kruskal-Wallis H statistic (with tie correction) computed from pooled ranks.
    """
    groups = [g[np.isfinite(g)] for g in groups if g.size > 0]
    if len(groups) < 2:
        return float("nan")
    n_total = int(sum(g.size for g in groups))
    if n_total <= 1:
        return float("nan")

    x = np.concatenate(groups)
    ranks = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    sizes = [g.size for g in groups]

    # Sum of ranks by group
    r_sums = []
    i0 = 0
    for n in sizes:
        r_sums.append(float(ranks[i0 : i0 + n].sum()))
        i0 += n

    h = (12.0 / (n_total * (n_total + 1.0))) * sum((rs * rs) / n for rs, n in zip(r_sums, sizes))
    h -= 3.0 * (n_total + 1.0)

    # Tie correction
    _, counts = np.unique(x, return_counts=True)
    t_num = float(np.sum(counts**3 - counts))
    t_den = float(n_total**3 - n_total)
    if t_den > 0:
        c = 1.0 - (t_num / t_den)
        if c > 0:
            h = h / c
    return float(h)


def _permutation_kw_pvalue(
    groups: list[np.ndarray],
    n_perm: int = 2000,
    seed: int = 1,
) -> tuple[float, float]:
    """
    Distribution-free permutation p-value for group difference using Kruskal H statistic.
    """
    groups = [g[np.isfinite(g)] for g in groups if g.size > 0]
    if len(groups) < 2:
        return float("nan"), float("nan")

    obs = _kruskal_h_stat(groups)
    if not np.isfinite(obs):
        return float("nan"), float("nan")

    sizes = [g.size for g in groups]
    all_vals = np.concatenate(groups).copy()
    rng = np.random.default_rng(seed)

    ge = 0
    for _ in range(max(1, int(n_perm))):
        rng.shuffle(all_vals)
        pos = 0
        perm_groups = []
        for n in sizes:
            perm_groups.append(all_vals[pos : pos + n])
            pos += n
        h = _kruskal_h_stat(perm_groups)
        if np.isfinite(h) and h >= obs:
            ge += 1
    p = (ge + 1.0) / (max(1, int(n_perm)) + 1.0)
    return float(obs), float(p)


def _alpha_metrics(tab: pd.DataFrame, sample: str) -> dict[str, float]:
    if sample not in tab.columns:
        return {k: float("nan") for k, _ in METRICS}
    vals = tab[sample].astype(float).to_numpy()
    vals = vals[vals > 0]
    if vals.size == 0:
        return {
            "observed_species": 0.0,
            "shannon": 0.0,
            "effective_shannon": 1.0,
            "simpson": 0.0,
            "effective_simpson": 0.0,
        }
    p = vals / float(vals.sum())
    observed = float(p.size)
    shannon = float(-(p * np.log(p)).sum())
    eff_shannon = float(np.exp(shannon))
    lambda2 = float((p * p).sum())
    simpson = float(1.0 - lambda2)
    eff_simpson = float(1.0 / lambda2) if lambda2 > 0 else float("nan")
    return {
        "observed_species": observed,
        "shannon": shannon,
        "effective_shannon": eff_shannon,
        "simpson": simpson,
        "effective_simpson": eff_simpson,
    }


def run(cfg) -> int:
    """
    Figure 07 (SVG):
      Alpha diversity error plots, comparing methods against expected:
        - ASV/zOTU against expected tip-level composition
        - OTU/TIC/TAC against expected species-level composition
      Expected tables and sample labels are read from tacticbench-prepared outputs
      under outdir/communities/<cid>/expected to keep sample mapping consistent
      with other figures (e.g., Figure 09).
        - Observed species
        - Shannon and Effective Shannon
        - Simpson and Effective Simpson
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig07", []) or [])

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 07: manifest.tsv not found; run discover first.")
        return 0
    mf = pd.read_csv(mf_p, sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    rows: list[dict] = []
    for _, r in mf.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue
        if _is_excluded(cid, exclude_patterns):
            continue

        ed = expected_dir(outdir, cid)
        exp_sp_p = ed / "expected_species_abundance_tss1000.tsv"
        if not exp_sp_p.exists():
            continue
        exp_sp = pd.read_csv(exp_sp_p, sep="\t", index_col=0)
        exp_sp.columns = exp_sp.columns.astype(str)
        exp_tip = pd.DataFrame()
        exp_tip_p = ed / "expected_tip_abundance_tss1000.tsv"
        if exp_tip_p.exists():
            try:
                exp_tip = pd.read_csv(exp_tip_p, sep="\t", index_col=0)
                exp_tip.columns = exp_tip.columns.astype(str)
            except Exception:
                exp_tip = pd.DataFrame()

        method_tabs: dict[str, pd.DataFrame] = {}
        method_expected: dict[str, pd.DataFrame] = {}
        for method in METHOD_ORDER:
            md = method_dir(outdir, cid, method)
            mp = md / "table_std_tss1000.tsv"
            if mp.exists():
                tab = pd.read_csv(mp, sep="\t", index_col=0)
                tab.columns = tab.columns.astype(str)
                target_level = _resolve_target_level(cfg, method)
                if target_level == "tip":
                    if exp_tip.empty:
                        log.info(f"[{cid}] Figure 07: expected tip table missing; skipping {method}.")
                        continue
                    exp_ref = exp_tip
                else:
                    exp_ref = exp_sp
                method_tabs[method] = tab
                method_expected[method] = exp_ref

        for method in METHOD_ORDER:
            tab = method_tabs.get(method)
            exp_ref = method_expected.get(method)
            if not isinstance(tab, pd.DataFrame) or not isinstance(exp_ref, pd.DataFrame):
                continue
            tab_cols = set(tab.columns.astype(str).tolist())
            shared_samples = [str(s) for s in exp_ref.columns if str(s) in tab_cols]
            if not shared_samples:
                log.info(f"[{cid}] Figure 07: no shared samples for {method}; skipping.")
                continue
            for s in shared_samples:
                exp_m = _alpha_metrics(exp_ref, s)
                inf_m = _alpha_metrics(tab, s)
                rec = {
                    "community_id": cid,
                    "sample_id": s,
                    "sample_key": f"{cid}|{s}",
                    "method": method,
                }
                for k, _ in METRICS:
                    rec[f"expected_{k}"] = float(exp_m[k])
                    rec[f"inferred_{k}"] = float(inf_m[k])
                    rec[f"error_{k}"] = float(inf_m[k] - exp_m[k])
                    rec[f"abs_error_{k}"] = float(abs(inf_m[k] - exp_m[k]))
                rows.append(rec)

    if not rows:
        log.info("Figure 07: no rows to plot.")
        return 0

    df = pd.DataFrame(rows)
    out_dir = figures_root(outdir) / "fig07_alpha_diversity_error"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "alpha_diversity_errors.tsv", sep="\t", index=False)
    write_long_tsv = bool(getattr(cfg.figures, "fig07_write_long_tsv", True))
    if write_long_tsv:
        long_blocks: list[pd.DataFrame] = []
        for metric_key, metric_label in METRICS:
            block = df[
                [
                    "community_id",
                    "sample_id",
                    "sample_key",
                    "method",
                    f"expected_{metric_key}",
                    f"inferred_{metric_key}",
                    f"error_{metric_key}",
                    f"abs_error_{metric_key}",
                ]
            ].copy()
            block = block.rename(
                columns={
                    f"expected_{metric_key}": "expected",
                    f"inferred_{metric_key}": "inferred",
                    f"error_{metric_key}": "error",
                    f"abs_error_{metric_key}": "abs_error",
                }
            )
            block["metric"] = metric_key
            block["metric_label"] = metric_label
            long_blocks.append(block)
        if long_blocks:
            long_df = pd.concat(long_blocks, ignore_index=True)
            long_df = long_df[
                [
                    "community_id",
                    "sample_id",
                    "sample_key",
                    "method",
                    "metric",
                    "metric_label",
                    "expected",
                    "inferred",
                    "error",
                    "abs_error",
                ]
            ]
            long_df.to_csv(out_dir / "alpha_diversity_errors_long.tsv", sep="\t", index=False)

    n_perm = int(getattr(cfg.figures, "fig07_stats_permutations", 2000))
    stats_seed = int(getattr(cfg.figures, "fig07_stats_seed", 1))
    show_pvalue_label = bool(getattr(cfg.figures, "fig07_show_pvalue_label", True))
    stats_model = _normalize_stats_model(getattr(cfg.figures, "fig07_stats_model", "mixed_effect_glmmtmb"))
    stats_pairwise_test = _normalize_pairwise_test(
        getattr(cfg.figures, "fig07_stats_pairwise_test", "signed_rank_permutation"),
        default="signed_rank_permutation",
    )
    stats_glmm_emmeans_scale = _normalize_emmeans_scale(
        getattr(cfg.figures, "fig07_stats_glmm_emmeans_scale", "response")
    )
    stats_alpha = float(getattr(cfg.figures, "fig07_stats_alpha", 0.05))
    stats_multiple_test_correction = bool(getattr(cfg.figures, "fig07_stats_multiple_test_correction", True))
    stats_show_pairwise = bool(getattr(cfg.figures, "fig07_stats_show_pairwise", True))
    plot_type = _normalize_plot_type(getattr(cfg.figures, "fig07_plot_type", "boxplot"), default="boxplot")
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))
    def _render_metric_set(
        metric_defs: list[tuple[str, str]],
        fig_name: str,
        stats_name: str,
        title_text: str,
    ) -> None:
        fig = plt.figure(figsize=(max(10.8, 3.6 * len(metric_defs)), 4.6))
        n = len(metric_defs)
        rng = np.random.default_rng(1)
        stats_rows: list[dict] = []
        pairwise_rows: list[dict] = []

        for i, (metric_key, metric_label) in enumerate(metric_defs, start=1):
            ax = fig.add_subplot(1, n, i)
            y_col = f"abs_error_{metric_key}"
            groups_for_test: list[np.ndarray] = []
            method_data: list[tuple[int, str, np.ndarray]] = []
            for j, method in enumerate(METHOD_ORDER, start=1):
                sub = df[df["method"] == method]
                vals = sub[y_col].to_numpy(dtype=float)
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    continue
                groups_for_test.append(vals)
                method_data.append((j, method, vals))

            if plot_type == "violin":
                violin_vals: list[np.ndarray] = []
                violin_pos: list[float] = []
                violin_methods: list[str] = []
                for j, method, vals in method_data:
                    if vals.size >= 2:
                        violin_vals.append(vals)
                        violin_pos.append(float(j))
                        violin_methods.append(method)
                if violin_vals:
                    vp = ax.violinplot(
                        violin_vals,
                        positions=violin_pos,
                        widths=0.64,
                        showmeans=False,
                        showmedians=True,
                        showextrema=False,
                    )
                    for body, method in zip(vp["bodies"], violin_methods):
                        body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                        body.set_edgecolor("#333333")
                        body.set_alpha(0.32)
                    if "cmedians" in vp:
                        vp["cmedians"].set_color("#111111")
                        vp["cmedians"].set_linewidth(1.2)
            else:
                for j, method, vals in method_data:
                    ax.boxplot(
                        [vals],
                        positions=[j],
                        widths=0.5,
                        patch_artist=True,
                        boxprops=dict(facecolor=METHOD_COLORS.get(method, "#999999"), alpha=0.32, edgecolor="#333333"),
                        medianprops=dict(color="#111111", linewidth=1.2),
                        whiskerprops=dict(color="#333333", linewidth=0.9),
                        capprops=dict(color="#333333", linewidth=0.9),
                    )

            for j, method, vals in method_data:
                xj = j + rng.normal(0, 0.045, size=vals.size)
                ax.scatter(
                    xj,
                    vals,
                    s=point_size,
                    alpha=0.45,
                    color=METHOD_COLORS.get(method, "#666666"),
                    edgecolors="none",
                )

            ax.set_title(metric_label, fontsize=10)
            ax.set_xticks(range(1, len(METHOD_ORDER) + 1))
            ax.set_xticklabels(display_methods(METHOD_ORDER), rotation=45, ha="right", fontsize=8)
            ax.set_ylabel("Absolute Error", fontsize=8)
            ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)

            panel_df = df[["community_id", "sample_key", "method", y_col]].copy()
            panel_df = panel_df.rename(columns={y_col: "value"})
            st = _compute_metric_stats(
                panel_df=panel_df,
                method_order=METHOD_ORDER,
                n_perm=n_perm,
                seed=stats_seed + i,
                stats_model=stats_model,
                pairwise_test=stats_pairwise_test,
                emmeans_scale=stats_glmm_emmeans_scale,
                alpha=stats_alpha,
                multiple_test_correction=stats_multiple_test_correction,
                log=log,
            )
            p = float(st.get("omnibus_p", np.nan))
            stat = float(st.get("omnibus_stat", np.nan))
            p_txt = f"{p:.3g}" if np.isfinite(p) else "NA"
            stat_txt = f"{stat:.3f}" if np.isfinite(stat) else "NA"
            model_used = str(st.get("stats_model", stats_model))
            pairwise_used = str(st.get("pairwise_test", stats_pairwise_test))
            if model_used == "mixed_effect_glmmtmb":
                ann = f"GLMM LRT={stat_txt}\np={p_txt}"
            else:
                ann = f"KW-perm H={stat_txt}\np={p_txt}"
            if show_pvalue_label:
                ax.text(
                    0.02,
                    0.98,
                    ann,
                    transform=ax.transAxes,
                    va="top",
                    ha="left",
                    fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
                )
            stats_rows.append(
                {
                    "metric": metric_key,
                    "label": metric_label,
                    "stats_model": model_used,
                    "pairwise_test": pairwise_used,
                    "omnibus_stat": stat,
                    "omnibus_df": float(st.get("omnibus_df", np.nan)),
                    "omnibus_p": p,
                    "model_family": str(st.get("model_family", "")),
                    "note": str(st.get("note", "")),
                    "n_blocks": int(st.get("n_blocks", 0)),
                    # Backward-compatible fields:
                    "kruskal_h": stat if model_used == "paired_permutation" else float("nan"),
                    "p_value_perm": p if model_used == "paired_permutation" else float("nan"),
                    "n_permutations": n_perm,
                }
            )
            for row in list(st.get("pairwise", [])):
                pairwise_rows.append(
                    {
                        "metric": metric_key,
                        "label": metric_label,
                        "stats_model": model_used,
                        "multiple_test_correction": bool(stats_multiple_test_correction),
                        "method_a": str(row.get("method_a", "")),
                        "method_b": str(row.get("method_b", "")),
                        "estimate": float(row.get("estimate", np.nan)),
                        "std_error": float(row.get("std_error", np.nan)),
                        "statistic": float(row.get("statistic", np.nan)),
                        "p_raw": float(row.get("p_raw", np.nan)),
                        "p_adj_bh": float(row.get("p_adj_bh", np.nan)),
                        "p_value_used": float(row.get("p_value_used", np.nan)),
                        "pairwise_test": str(row.get("pairwise_test", pairwise_used)),
                        "significant_alpha": bool(row.get("significant_alpha", False)),
                    }
                )
            if stats_show_pairwise:
                sig_pairs = [r for r in list(st.get("pairwise", [])) if bool(r.get("significant_alpha", False))]
                _draw_pairwise_brackets(
                    ax=ax,
                    pairwise_rows=sig_pairs,
                    method_order=METHOD_ORDER,
                    alpha=stats_alpha,
                    placement="top",
                    spacing_scale=1.2,
                )

        fig.suptitle(str(title_text), fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(out_dir / fig_name, format="svg")
        plt.close(fig)
        if stats_rows:
            pd.DataFrame(stats_rows).to_csv(out_dir / stats_name, sep="\t", index=False)
        if pairwise_rows:
            pairwise_name = stats_name.replace(".tsv", "_pairwise.tsv")
            pd.DataFrame(pairwise_rows).to_csv(out_dir / pairwise_name, sep="\t", index=False)

    title = getattr(cfg.figures, "fig07_title", None) or "Alpha Diversity Error Plots"
    _render_metric_set(
        metric_defs=METRICS,
        fig_name="alpha_diversity_error_plots.svg",
        stats_name="alpha_diversity_error_stats.tsv",
        title_text=str(title),
    )
    include_richness_effective_variant = bool(getattr(cfg.figures, "fig07_include_richness_effective_variant", True))
    if include_richness_effective_variant:
        _render_metric_set(
            metric_defs=METRICS_RICHNESS_EFFECTIVE,
            fig_name="alpha_diversity_error_plots_richness_effective.svg",
            stats_name="alpha_diversity_error_stats_richness_effective.tsv",
            title_text=f"{title} (Richness + Effective Diversity)",
        )

    log.info(f"Wrote Figure 07 to {out_dir}")
    return 0
