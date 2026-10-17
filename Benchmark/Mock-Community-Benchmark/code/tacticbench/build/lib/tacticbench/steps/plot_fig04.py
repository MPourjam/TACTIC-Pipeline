from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from ..utils.labels import display_method

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_COLORS = {
    # Color-blind-friendly palette (Okabe-Ito family)
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _species_richness_from_table(tab: pd.DataFrame, sample: str) -> int:
    if sample not in tab.columns:
        return 0
    return int((tab[sample].astype(float) > 0).sum())


def _tip_to_species_table(tab: pd.DataFrame, tip_to_species: dict[str, str]) -> pd.DataFrame:
    tmp = tab.copy()
    tmp["__species__"] = [str(tip_to_species.get(str(t), str(t))) for t in tmp.index]
    out = tmp.groupby("__species__", sort=False).sum(numeric_only=True)
    return out


def run(cfg) -> int:
    """
    Figure 04 (SVG):
      Species Richness Recovery scatter plot:
        - x-axis: expected species richness per sample
        - y-axis: recovered species richness per sample
        - one point per sample per method
        - color by method
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig04", []) or [])

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 04: manifest.tsv not found; run discover first.")
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

        tip_to_species: dict[str, str] = {}
        tip2sp_p = ed / "tip_to_species.tsv"
        if tip2sp_p.exists():
            try:
                m = pd.read_csv(tip2sp_p, sep="\t")
                if {"tip", "species"}.issubset(set(m.columns)):
                    tip_to_species = dict(zip(m["tip"].astype(str), m["species"].astype(str)))
            except Exception:
                tip_to_species = {}

        for method in METHOD_ORDER:
            md = method_dir(outdir, cid, method)
            if method in {"OTU", "TIC", "TAC"}:
                m_p = md / "mapped_species_table_tss1000.tsv"
                if not m_p.exists():
                    continue
                inf = pd.read_csv(m_p, sep="\t", index_col=0)
            else:
                m_p = md / "mapped_tip_table_tss1000.tsv"
                if not m_p.exists():
                    continue
                inf_tip = pd.read_csv(m_p, sep="\t", index_col=0)
                inf = _tip_to_species_table(inf_tip, tip_to_species)

            for s in samples:
                exp_r = _species_richness_from_table(exp, s)
                inf_r = _species_richness_from_table(inf, s)
                rows.append(
                    {
                        "community_id": cid,
                        "sample_id": s,
                        "method": method,
                        "expected_richness": exp_r,
                        "recovered_richness": inf_r,
                    }
                )

    if not rows:
        log.info("Figure 04: no rows to plot.")
        return 0

    df = pd.DataFrame(rows)
    max_v = int(max(df["expected_richness"].max(), df["recovered_richness"].max()))

    fig = plt.figure(figsize=(8.8, 6.6))
    ax = fig.add_subplot(1, 1, 1)

    for method in METHOD_ORDER:
        sub = df[df["method"] == method]
        if sub.empty:
            continue
        ax.scatter(
            sub["expected_richness"],
            sub["recovered_richness"],
            s=28,
            alpha=0.78,
            color=METHOD_COLORS.get(method, "#444444"),
            label=display_method(method),
            edgecolors="none",
        )

    ax.plot([0, max_v], [0, max_v], linestyle="--", color="#222222", linewidth=1.4, label="1:1")
    ax.set_xlim(0, max_v * 1.02 if max_v > 0 else 1)
    ax.set_ylim(0, max_v * 1.02 if max_v > 0 else 1)
    ax.set_xlabel("Expected Species Richness")
    ax.set_ylabel("Recovered Species Richness")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", frameon=False, fontsize=8)

    title = getattr(cfg.figures, "fig04_title", None) or "Species Richness Recovery"
    ax.set_title(str(title), fontsize=11)
    fig.tight_layout()

    out_dir = figures_root(outdir) / "fig04_species_richness_recovery"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_svg = out_dir / "species_richness_recovery.svg"
    fig.savefig(out_svg, format="svg")
    plt.close(fig)

    # Also export source points for downstream custom plotting/reporting
    df.to_csv(out_dir / "species_richness_recovery_points.tsv", sep="\t", index=False)

    log.info(f"Wrote Figure 04 to {out_dir}")
    return 0
