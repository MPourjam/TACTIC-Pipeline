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


def _load_tip_to_species(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    try:
        m = pd.read_csv(path, sep="\t")
        if {"tip", "species"}.issubset(set(m.columns)):
            return dict(zip(m["tip"].astype(str), m["species"].astype(str)))
    except Exception:
        return {}
    return {}


def _vec_with_other(series: pd.Series, species_order: list[str]) -> np.ndarray:
    vals = np.array([float(series.get(sp, 0.0)) for sp in species_order], dtype=float)
    return vals


def run(cfg) -> int:
    """
    Figure 10 (SVG):
      Species confusion heatmaps by method.

      Construction (abundance-based, sample-aggregated):
        - diagonal: shared mass min(expected_i, inferred_i)
        - off-diagonal: expected deficits redistributed to inferred excess species
          proportionally to excess abundance.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig10", []) or [])
    top_n = int(getattr(cfg.figures, "fig10_top_n_species", 30))

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 10: manifest.tsv not found; run discover first.")
        return 0
    mf = pd.read_csv(mf_p, sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    # First pass: global expected species mass to select top-N species
    expected_mass: dict[str, float] = {}
    communities: list[tuple[str, Path, pd.DataFrame, list[str], dict[str, str]]] = []
    for _, r in mf.iterrows():
        cid = str(r["community_id"])
        if cid in skipped or _is_excluded(cid, exclude_patterns):
            continue
        ed = expected_dir(outdir, cid)
        exp_p = ed / "expected_species_abundance_tss1000.tsv"
        if not exp_p.exists():
            continue
        exp = pd.read_csv(exp_p, sep="\t", index_col=0)
        samples = [str(s) for s in exp.columns]
        tip_to_species = _load_tip_to_species(ed / "tip_to_species.tsv")
        communities.append((cid, ed, exp, samples, tip_to_species))
        for sp_name, v in exp.sum(axis=1).items():
            k = str(sp_name)
            expected_mass[k] = expected_mass.get(k, 0.0) + float(v)

    if not communities:
        log.info("Figure 10: no community data to plot.")
        return 0

    top_species = [k for k, _ in sorted(expected_mass.items(), key=lambda kv: kv[1], reverse=True)[: max(1, top_n)]]
    species_order = top_species + ["Other"]
    n = len(species_order)

    out_dir = figures_root(outdir) / "fig10_species_confusion_heatmaps"
    out_dir.mkdir(parents=True, exist_ok=True)

    title = getattr(cfg.figures, "fig10_title", None) or "Species Confusion Heatmaps"

    for method in METHOD_ORDER:
        M = np.zeros((n, n), dtype=float)
        for cid, ed, exp, samples, tip_to_species in communities:
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

            for s in samples:
                if s not in exp.columns:
                    continue
                exp_s = exp[s].astype(float)
                if s in inf.columns:
                    inf_s = inf[s].astype(float)
                else:
                    inf_s = pd.Series(dtype=float)

                # top species + other bucket
                e_top = _vec_with_other(exp_s, top_species)
                i_top = _vec_with_other(inf_s, top_species)
                e_other = max(0.0, float(exp_s.sum() - e_top.sum()))
                i_other = max(0.0, float(inf_s.sum() - i_top.sum()))
                e = np.append(e_top, e_other)
                q = np.append(i_top, i_other)

                shared = np.minimum(e, q)
                M += np.diag(shared)
                deficit = e - shared
                excess = q - shared
                ex_sum = float(excess.sum())
                if ex_sum > 0:
                    for i in range(n):
                        if deficit[i] > 0:
                            M[i, :] += deficit[i] * (excess / ex_sum)

        # Normalize by expected row mass for interpretability
        row_sums = M.sum(axis=1, keepdims=True)
        M_norm = np.divide(M, row_sums, out=np.zeros_like(M), where=row_sums > 0)

        raw_df = pd.DataFrame(M, index=species_order, columns=species_order)
        norm_df = pd.DataFrame(M_norm, index=species_order, columns=species_order)
        raw_df.to_csv(out_dir / f"{method}_confusion_raw.tsv", sep="\t")
        norm_df.to_csv(out_dir / f"{method}_confusion_row_normalized.tsv", sep="\t")

        label_fs = 9 if n <= 35 else 8
        fig_w = max(9.0, n * 0.34)
        fig_h = max(7.5, n * 0.30)
        fig = plt.figure(figsize=(fig_w, fig_h))
        ax = fig.add_subplot(1, 1, 1)
        im = ax.imshow(M_norm, aspect="auto", cmap="cividis", vmin=0.0, vmax=1.0)
        ax.set_xticks(range(n))
        ax.set_xticklabels(species_order, rotation=90, fontsize=label_fs)
        ax.set_yticks(range(n))
        ax.set_yticklabels(species_order, fontsize=label_fs)
        ax.set_xlabel("Inferred Species", fontsize=10)
        ax.set_ylabel("Expected Species", fontsize=10)
        ax.set_title(f"{title}: {display_method(method)}", fontsize=12)
        cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
        cbar.ax.tick_params(labelsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / f"{method}_species_confusion_heatmap.svg", format="svg")
        plt.close(fig)

    log.info(f"Wrote Figure 10 to {out_dir}")
    return 0
