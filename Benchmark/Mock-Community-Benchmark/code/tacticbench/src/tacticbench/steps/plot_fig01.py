from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from ..utils.io import load_fasta
from ..utils.labels import display_community_id, display_method, display_sample_key

BAR_ORDER = ["zOTU", "ASV", "EXPECTED", "TIC", "TAC", "OTU"]
BAR_COLORS = {
    # Color-blind-safe palette (Okabe-Ito family)
    "zOTU": "#56B4E9",
    "ASV": "#0072B2",
    "EXPECTED": "#F0E442",
    "EXPECTED_FEATURES": "#999999",
    "TIC": "#009E73",
    "TAC": "#E69F00",
    "OTU": "#D55E00",
}


def _count_nonzero(table: pd.DataFrame, sample: str) -> int:
    if sample not in table.columns:
        return 0
    return int((table[sample] > 0).sum())


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _community_sort_key(cid: str) -> tuple[int, int, str]:
    """
    Order communities for publication-style Figure 01 panels:
      1) Reindexed mocks ("Mock N" display names) in numeric order
      2) Named controls (ZIEL I/II, ZYMO)
      3) Any remaining communities in lexical order
    """
    disp = display_community_id(str(cid))
    if disp.startswith("Mock "):
        try:
            return (0, int(disp.split(" ", 1)[1]), "")
        except Exception:
            return (0, 10**9, str(cid))
    if disp == "ZIEL I":
        return (1, 1, "")
    if disp == "ZIEL II":
        return (1, 2, "")
    if disp == "ZYMO":
        return (1, 3, "")
    return (2, 10**9, str(cid))


def run(cfg) -> int:
    """
    Figure 01 (SVG):
      For each sample (across all communities), draw bars (left->right):
        1) zOTU count (TACTIC de-novo zOTU)
        2) ASV count (DADA2)
        3) Expected species count (dashed)   [= #rows in taxonomy_collapsed_species.tsv]
        4) TIC SOTU count
        5) TAC OTU count
        6) OTU count (de-novo UPARSE OTUs)

      One subplot per sample in a grid of cfg.figures.grid_rows x cfg.figures.grid_cols.

    Source tables for counts:
      - methods/<method>/table_std_tss1000.tsv   (features x expected samples)
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig01", []) or [])

    mf = pd.read_csv(outdir / "manifest.tsv", sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    rows = []

    for _, r in mf.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue
        if _is_excluded(cid, exclude_patterns):
            continue

        ed = expected_dir(outdir, cid)
        exp_sp_p = ed / "taxonomy_collapsed_species.tsv"
        exp_seq_p = ed / "expected_sequences.fasta"
        if not exp_sp_p.exists():
            continue

        exp_sp = pd.read_csv(exp_sp_p, sep="\t", index_col=0)
        expected_n = int(exp_sp.shape[0])
        samples = list(exp_sp.columns.astype(str))
        expected_features_n = 0
        if exp_seq_p.exists():
            try:
                seqs = load_fasta(exp_seq_p)
                expected_features_n = int(len(set(seqs.values())))
            except Exception:
                expected_features_n = 0

        # load method tables once per community (if present)
        tables = {}
        for method in ["zOTU", "ASV", "TIC", "TAC", "OTU"]:
            tp = method_dir(outdir, cid, method) / "table_std_tss1000.tsv"
            if tp.exists():
                tables[method] = pd.read_csv(tp, sep="\t", index_col=0)
            else:
                tables[method] = None

        for s in samples:
            rec = {
                "community_id": cid,
                "sample_id": s,
                "sample_key": f"{cid}|{s}",
                "EXPECTED": expected_n,
                "EXPECTED_FEATURES": expected_features_n,
            }
            for method in ["zOTU", "ASV", "TIC", "TAC", "OTU"]:
                t = tables.get(method)
                rec[method] = _count_nonzero(t, s) if isinstance(t, pd.DataFrame) else 0
            rows.append(rec)

    if not rows:
        log.info("Figure 01: no samples found; nothing to plot.")
        return 0

    df = pd.DataFrame(rows)
    ordered_communities = sorted(df["community_id"].astype(str).unique().tolist(), key=_community_sort_key)
    community_rank = {cid: i for i, cid in enumerate(ordered_communities)}
    community_sample_counts = (
        df.groupby("community_id")["sample_id"].nunique().astype(int).to_dict()
    )
    df["__community_rank"] = df["community_id"].astype(str).map(community_rank).fillna(10**9).astype(int)
    df = df.sort_values(["__community_rank", "sample_id"]).drop(columns=["__community_rank"]).reset_index(drop=True)

    R = int(getattr(cfg.figures, "grid_rows", 5))
    C = int(getattr(cfg.figures, "grid_cols", 6))

    def _plot_grid(order: list[str], labels: list[str], out_svg: Path, expected_idx: int, fig_title: str | None) -> None:
        fig = plt.figure(figsize=(C * 4.5, R * 4))

        for i in range(R * C):
            ax = fig.add_subplot(R, C, i + 1)
            if i >= len(df):
                ax.axis("off")
                continue

            row = df.iloc[i]
            vals = [row[k] for k in order]
            x = list(range(len(vals)))

            ax.bar(x, vals, color=[BAR_COLORS[k] for k in order])
            ax.bar(expected_idx, row["EXPECTED"], fill=False, edgecolor="#111111", linestyle="--", linewidth=2)

            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=18)
            ax.tick_params(axis="y", labelsize=20)
            cid = str(row["community_id"])
            if int(community_sample_counts.get(cid, 0)) <= 1:
                panel_title = display_community_id(cid)
            else:
                panel_title = display_sample_key(str(row["sample_key"]))
            ax.set_title(panel_title, fontsize=20)

        if fig_title:
            fig.suptitle(fig_title, fontsize=22)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
        else:
            fig.tight_layout()
        fig.savefig(out_svg, format="svg")
        plt.close(fig)
        log.info(f"Wrote {out_svg}")

    fr = figures_root(outdir)
    fr.mkdir(parents=True, exist_ok=True)
    out_dir = fr / "fig01_feature_counts"
    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_grid(
        order=["zOTU", "ASV", "EXPECTED", "TIC", "TAC", "OTU"],
        labels=[display_method("zOTU"), display_method("ASV"), "Exp", display_method("TIC"), display_method("TAC"), display_method("OTU")],
        out_svg=out_dir / "fig01_feature_counts_grid.svg",
        expected_idx=2,
        fig_title=getattr(cfg.figures, "fig01_title", None),
    )
    _plot_grid(
        order=["zOTU", "ASV", "EXPECTED_FEATURES", "EXPECTED", "TIC", "TAC", "OTU"],
        labels=[display_method("zOTU"), display_method("ASV"), "ExpFeat", "Exp", display_method("TIC"), display_method("TAC"), display_method("OTU")],
        out_svg=out_dir / "fig01_feature_counts_grid_with_expected_features.svg",
        expected_idx=3,
        fig_title=getattr(cfg.figures, "fig01_title_with_expected_features", None),
    )

    return 0
