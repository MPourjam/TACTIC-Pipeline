from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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


def _presence_set(tab: pd.DataFrame, sample: str) -> set[str]:
    if sample not in tab.columns:
        return set()
    s = tab[sample].astype(float)
    return set(s.index[s > 0].astype(str))


def _tip_to_species_table(tab: pd.DataFrame, tip_to_species: dict[str, str]) -> pd.DataFrame:
    tmp = tab.copy()
    tmp["__species__"] = [str(tip_to_species.get(str(t), str(t))) for t in tmp.index]
    out = tmp.groupby("__species__", sort=False).sum(numeric_only=True)
    return out


def run(cfg) -> int:
    """
    Figure 05 (SVG):
      Species-presence Precision-Recall scatter:
        - one point per sample per method
        - x-axis: Recall
        - y-axis: Precision
        - color by method
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig05", []) or [])

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 05: manifest.tsv not found; run discover first.")
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
                exp_set = _presence_set(exp, s)
                inf_set = _presence_set(inf, s)
                tp = len(exp_set & inf_set)
                fp = len(inf_set - exp_set)
                fn = len(exp_set - inf_set)
                precision = float(tp / (tp + fp)) if (tp + fp) > 0 else np.nan
                recall = float(tp / (tp + fn)) if (tp + fn) > 0 else np.nan
                f1 = float((2 * precision * recall) / (precision + recall)) if (precision + recall) > 0 else np.nan
                rows.append(
                    {
                        "community_id": cid,
                        "sample_id": s,
                        "method": method,
                        "tp": tp,
                        "fp": fp,
                        "fn": fn,
                        "precision": precision,
                        "recall": recall,
                        "f1": f1,
                    }
                )

    if not rows:
        log.info("Figure 05: no rows to plot.")
        return 0

    df = pd.DataFrame(rows)
    out_dir = figures_root(outdir) / "fig05_species_presence_precision_recall"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "species_presence_precision_recall_points.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(8.8, 6.6))
    ax = fig.add_subplot(1, 1, 1)

    for method in METHOD_ORDER:
        sub = df[df["method"] == method].dropna(subset=["precision", "recall"])
        if sub.empty:
            continue
        ax.scatter(
            sub["recall"],
            sub["precision"],
            s=30,
            alpha=0.78,
            color=METHOD_COLORS.get(method, "#444444"),
            label=display_method(method),
            edgecolors="none",
        )
        # method centroid marker
        ax.scatter(
            [sub["recall"].mean()],
            [sub["precision"].mean()],
            s=95,
            marker="X",
            color=METHOD_COLORS.get(method, "#444444"),
            edgecolors="#111111",
            linewidths=0.6,
        )

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Recall (Species Presence)")
    ax.set_ylabel("Precision (Species Presence)")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    title = getattr(cfg.figures, "fig05_title", None) or "Precision-Recall (Species Presence)"
    ax.set_title(str(title), fontsize=11)
    fig.tight_layout()

    out_svg = out_dir / "precision_recall_species_presence.svg"
    fig.savefig(out_svg, format="svg")
    plt.close(fig)

    log.info(f"Wrote Figure 05 to {out_dir}")
    return 0
