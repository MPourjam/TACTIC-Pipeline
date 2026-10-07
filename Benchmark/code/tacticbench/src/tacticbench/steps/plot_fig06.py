from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from ..utils.labels import display_method, display_sample_key

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
LINE_COLORS = {
    "EXPECTED": "#111111",
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


def _tip_to_species_table(tab: pd.DataFrame, tip_to_species: dict[str, str]) -> pd.DataFrame:
    tmp = tab.copy()
    tmp["__species__"] = [str(tip_to_species.get(str(t), str(t))) for t in tmp.index]
    out = tmp.groupby("__species__", sort=False).sum(numeric_only=True)
    return out


def _rank_abundance(tab: pd.DataFrame, sample: str) -> np.ndarray:
    if sample not in tab.columns:
        return np.array([], dtype=float)
    vals = tab[sample].astype(float).to_numpy()
    vals = vals[vals > 0]
    if vals.size == 0:
        return np.array([], dtype=float)
    vals = vals / float(vals.sum())
    vals = np.sort(vals)[::-1]
    return vals


def run(cfg) -> int:
    """
    Figure 06 (SVG):
      Rank-Abundance curves (species):
        - one curve for expected species composition
        - one curve per method (species-space mapped composition)
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig06", []) or [])

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 06: manifest.tsv not found; run discover first.")
        return 0
    mf = pd.read_csv(mf_p, sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    sample_curves: dict[str, dict[str, np.ndarray]] = {}
    point_rows: list[dict] = []

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

        method_tabs: dict[str, pd.DataFrame] = {}
        for method in METHOD_ORDER:
            md = method_dir(outdir, cid, method)
            if method in {"OTU", "TIC", "TAC"}:
                mp = md / "mapped_species_table_tss1000.tsv"
                if mp.exists():
                    method_tabs[method] = pd.read_csv(mp, sep="\t", index_col=0)
            else:
                mp = md / "mapped_tip_table_tss1000.tsv"
                if mp.exists():
                    tip_tab = pd.read_csv(mp, sep="\t", index_col=0)
                    method_tabs[method] = _tip_to_species_table(tip_tab, tip_to_species)

        for s in samples:
            sk = f"{cid}|{s}"
            curves: dict[str, np.ndarray] = {"EXPECTED": _rank_abundance(exp, s)}
            for method in METHOD_ORDER:
                tab = method_tabs.get(method)
                curves[method] = _rank_abundance(tab, s) if isinstance(tab, pd.DataFrame) else np.array([], dtype=float)
            sample_curves[sk] = curves

            for label, arr in curves.items():
                for i, v in enumerate(arr, start=1):
                    point_rows.append(
                        {
                            "community_id": cid,
                            "sample_id": s,
                            "sample_key": sk,
                            "source": label,
                            "rank": i,
                            "relative_abundance": float(v),
                        }
                    )

    if not sample_curves:
        log.info("Figure 06: no sample curves to plot.")
        return 0

    out_dir = figures_root(outdir) / "fig06_rank_abundance_species"
    out_dir.mkdir(parents=True, exist_ok=True)
    if point_rows:
        pd.DataFrame(point_rows).to_csv(out_dir / "rank_abundance_species_points.tsv", sep="\t", index=False)

    title = getattr(cfg.figures, "fig06_title", None) or "Rank-Abundance Curves (Species)"

    def _draw_on_ax(ax, curves: dict[str, np.ndarray], sub_title: str) -> None:
        max_rank = max((len(a) for a in curves.values()), default=0)
        for label in ["EXPECTED"] + METHOD_ORDER:
            arr = curves.get(label, np.array([], dtype=float))
            if arr.size == 0:
                continue
            x = np.arange(1, arr.size + 1)
            lw = 2.1 if label == "EXPECTED" else 1.6
            shown_label = "Expected" if label == "EXPECTED" else display_method(label)
            ax.plot(x, arr, color=LINE_COLORS.get(label, "#444444"), linewidth=lw, label=shown_label)
        ax.set_xlim(1, max(1, max_rank))
        ax.set_yscale("log")
        ax.set_xlabel("Rank", fontsize=7)
        ax.set_ylabel("Relative Abundance", fontsize=7)
        ax.grid(True, which="both", alpha=0.22, linewidth=0.5)
        ax.set_title(display_sample_key(sub_title), fontsize=8)

    # One SVG per sample
    for sk in sorted(sample_curves.keys()):
        fig = plt.figure(figsize=(6.8, 4.8))
        ax = fig.add_subplot(1, 1, 1)
        _draw_on_ax(ax, sample_curves[sk], sk)
        ax.legend(loc="upper right", frameon=False, fontsize=8)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        safe = sk.replace("|", "__").replace("/", "_")
        fig.savefig(out_dir / f"{safe}.svg", format="svg")
        plt.close(fig)

    # Paged overview grid
    keys = sorted(sample_curves.keys())
    R = int(getattr(cfg.figures, "fig06_grid_rows", None) or getattr(cfg.figures, "grid_rows", 5))
    C = int(getattr(cfg.figures, "fig06_grid_cols", None) or getattr(cfg.figures, "grid_cols", 6))
    per_page = max(1, R * C)

    for p0 in range(0, len(keys), per_page):
        chunk = keys[p0 : p0 + per_page]
        page = (p0 // per_page) + 1
        fig = plt.figure(figsize=(C * 3.1, R * 2.7))
        for i in range(per_page):
            ax = fig.add_subplot(R, C, i + 1)
            if i >= len(chunk):
                ax.axis("off")
                continue
            sk = chunk[i]
            _draw_on_ax(ax, sample_curves[sk], sk)
            if i == 0:
                handles, labels = ax.get_legend_handles_labels()
                if handles:
                    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=6, frameon=False, fontsize=8)

        fig.suptitle(f"{title} (page {page})", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        suffix = f"_p{page:02d}" if len(keys) > per_page else ""
        fig.savefig(out_dir / f"all_samples{suffix}.svg", format="svg")
        plt.close(fig)

    log.info(f"Wrote Figure 06 to {out_dir}")
    return 0
