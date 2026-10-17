from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from ..utils.labels import display_methods

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}


def _normalize_plot_type(raw: str, default: str = "boxplot") -> str:
    v = str(raw or "").strip().lower()
    if v in {"box", "boxplot"}:
        return "boxplot"
    if v in {"violin", "violinplot"}:
        return "violin"
    return default


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _tip_to_species_table(tab: pd.DataFrame, tip_to_species: dict[str, str]) -> pd.DataFrame:
    tmp = tab.copy()
    tmp["__species__"] = [str(tip_to_species.get(str(t), str(t))) for t in tmp.index]
    out = tmp.groupby("__species__", sort=False).sum(numeric_only=True)
    return out


def _braycurtis(u: np.ndarray, v: np.ndarray) -> float:
    den = float(np.sum(u + v))
    if den <= 0:
        return 0.0
    return float(np.sum(np.abs(u - v)) / den)


def _pairwise_braycurtis(mat: np.ndarray) -> np.ndarray:
    n = mat.shape[0]
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = _braycurtis(mat[i], mat[j])
            out[i, j] = d
            out[j, i] = d
    return out


def _upper_triangle_values(m: np.ndarray) -> np.ndarray:
    if m.shape[0] < 2:
        return np.array([], dtype=float)
    iu = np.triu_indices(m.shape[0], k=1)
    return m[iu].astype(float)


def _spearman_from_vectors(x: np.ndarray, y: np.ndarray) -> float:
    if x.size != y.size or x.size < 2:
        return float("nan")
    xr = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    yr = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if float(np.std(xr)) == 0.0 or float(np.std(yr)) == 0.0:
        return float("nan")
    return float(np.corrcoef(xr, yr)[0, 1])


def run(cfg) -> int:
    """
    Figure 08 (SVG):
      Beta-diversity concordance at species level, by method.

      For each community with >=3 samples:
        - compute expected pairwise Bray-Curtis distances between samples
        - compute inferred pairwise Bray-Curtis distances between samples
        - compute Spearman correlation between the two upper-triangle vectors
      Then summarize per-method concordance across communities.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig08", []) or [])

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 08: manifest.tsv not found; run discover first.")
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
        exp_p = ed / "expected_species_abundance_tss1000.tsv"
        if not exp_p.exists():
            continue
        exp = pd.read_csv(exp_p, sep="\t", index_col=0)
        samples = [str(s) for s in exp.columns]
        if len(samples) < 3:
            continue

        tip_to_species: dict[str, str] = {}
        tip2sp_p = ed / "tip_to_species.tsv"
        if tip2sp_p.exists():
            try:
                m = pd.read_csv(tip2sp_p, sep="\t")
                if {"tip", "species"}.issubset(set(m.columns)):
                    tip_to_species = dict(zip(m["tip"].astype(str), m["species"].astype(str)))
            except Exception:
                tip_to_species = {}

        # expected matrix: samples x species
        exp_s = exp.T.copy()
        exp_s.index = exp_s.index.astype(str)
        exp_species = [str(c) for c in exp_s.columns]
        exp_mat = exp_s.to_numpy(dtype=float)
        exp_d = _pairwise_braycurtis(exp_mat)
        exp_vec = _upper_triangle_values(exp_d)
        if exp_vec.size < 3:
            continue

        for method in METHOD_ORDER:
            md = method_dir(outdir, cid, method)
            if method in {"OTU", "TIC", "TAC"}:
                mp = md / "mapped_species_table_tss1000.tsv"
                if not mp.exists():
                    continue
                inf = pd.read_csv(mp, sep="\t", index_col=0)
            else:
                mp = md / "mapped_tip_table_tss1000.tsv"
                if not mp.exists():
                    continue
                tip_tab = pd.read_csv(mp, sep="\t", index_col=0)
                inf = _tip_to_species_table(tip_tab, tip_to_species)

            # align inferred to expected sample/species space
            inf = inf.reindex(index=exp_species, columns=samples).fillna(0.0)
            inf_s = inf.T
            inf_s.index = inf_s.index.astype(str)
            inf_mat = inf_s.to_numpy(dtype=float)
            inf_d = _pairwise_braycurtis(inf_mat)
            inf_vec = _upper_triangle_values(inf_d)
            rho = _spearman_from_vectors(exp_vec, inf_vec)
            rows.append(
                {
                    "community_id": cid,
                    "method": method,
                    "n_samples": len(samples),
                    "n_pairs": int(exp_vec.size),
                    "spearman_rho": rho,
                }
            )

    if not rows:
        log.info("Figure 08: no communities with enough samples to compute concordance.")
        return 0

    df = pd.DataFrame(rows)
    out_dir = figures_root(outdir) / "fig08_beta_diversity_concordance"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "beta_diversity_concordance.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(8.8, 5.8))
    ax = fig.add_subplot(1, 1, 1)
    rng = np.random.default_rng(1)
    plot_type = _normalize_plot_type(getattr(cfg.figures, "fig08_plot_type", "boxplot"), default="boxplot")
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))
    method_data: list[tuple[int, str, np.ndarray]] = []

    for i, method in enumerate(METHOD_ORDER, start=1):
        sub = df[df["method"] == method]
        vals = sub["spearman_rho"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        method_data.append((i, method, vals))

    if plot_type == "violin":
        violin_vals: list[np.ndarray] = []
        violin_pos: list[float] = []
        violin_methods: list[str] = []
        for i, method, vals in method_data:
            if vals.size >= 2:
                violin_vals.append(vals)
                violin_pos.append(float(i))
                violin_methods.append(method)
        if violin_vals:
            vp = ax.violinplot(
                violin_vals,
                positions=violin_pos,
                widths=0.72,
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
        for i, method, vals in method_data:
            ax.boxplot(
                [vals],
                positions=[i],
                widths=0.55,
                patch_artist=True,
                boxprops=dict(facecolor=METHOD_COLORS.get(method, "#999999"), alpha=0.32, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.3),
                whiskerprops=dict(color="#333333", linewidth=0.9),
                capprops=dict(color="#333333", linewidth=0.9),
            )

    for i, method, vals in method_data:
        xj = i + rng.normal(0, 0.045, size=vals.size)
        ax.scatter(xj, vals, s=point_size, alpha=0.55, color=METHOD_COLORS.get(method, "#666666"), edgecolors="none")

    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_ylim(-1.02, 1.02)
    ax.set_xticks(range(1, len(METHOD_ORDER) + 1))
    ax.set_xticklabels(display_methods(METHOD_ORDER), rotation=0, fontsize=9)
    ax.set_ylabel("Spearman Concordance (Expected vs Inferred Bray-Curtis)")
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    title = getattr(cfg.figures, "fig08_title", None) or "Beta-Diversity Concordance"
    ax.set_title(str(title), fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "beta_diversity_concordance.svg", format="svg")
    plt.close(fig)

    log.info(f"Wrote Figure 08 to {out_dir}")
    return 0
