from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.io import genus_species
from ..utils.labels import (
    configure_community_display_mapping_from_manifest,
    display_community_id,
    display_method,
)
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
BAR_ORDER = ["Expected"] + METHOD_ORDER
METHOD_ALIASES = {
    "asv": "ASV",
    "dada2": "ASV",
    "dada": "ASV",
    "zotu": "zOTU",
    "unoise3": "zOTU",
    "unoise": "zOTU",
    "otu": "OTU",
    "uparse": "OTU",
    "tic": "TIC",
    "tac": "TAC",
}

BAR_BASE_COLORS = {
    "Expected": "#6E6E6E",
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}


def _normalize_method_token(raw: object) -> str | None:
    key = str(raw).strip()
    if not key:
        return None
    return METHOD_ALIASES.get(key.lower(), key if key in METHOD_ORDER else None)


def _normalize_method_order(raw: object, default: list[str] | None = None) -> list[str]:
    fallback = list(default or METHOD_ORDER)
    if raw is None:
        return fallback
    if isinstance(raw, str):
        parts = [p.strip() for p in re.split(r"[,;]", raw) if p.strip()]
    else:
        try:
            parts = list(raw)  # type: ignore[arg-type]
        except TypeError:
            parts = [raw]
    out: list[str] = []
    for part in parts:
        method = _normalize_method_token(part)
        if method and method not in out:
            out.append(method)
    return out or fallback


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _community_grid_sort_key(community_id: str) -> tuple[int, int | str]:
    s = str(community_id)
    m = re.match(r"^mock-(\d+)$", s, flags=re.IGNORECASE)
    if m:
        return (0, int(m.group(1)))
    return (1, display_community_id(s).lower())


def _community_species_label_style(community_id: str) -> tuple[float, bool]:
    special_multiline = str(community_id).lower() in {"mock-12", "mock-16"}
    return (45.0 if special_multiline else 0.0, special_multiline)


def _read_numeric_table(path: Path) -> pd.DataFrame:
    try:
        t = pd.read_csv(path, sep="	", index_col=0)
    except Exception:
        return pd.DataFrame()
    if t.empty:
        return pd.DataFrame()
    t = t.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    t.index = t.index.astype(str)
    t.columns = t.columns.astype(str)
    return t


def _load_tip_to_species(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    try:
        t = pd.read_csv(path, sep="	")
    except Exception:
        return {}
    if t.empty or "tip" not in t.columns or "species" not in t.columns:
        return {}
    return dict(zip(t["tip"].astype(str), t["species"].astype(str)))


def _species_from_tip(tip: str, tip_to_species: dict[str, str]) -> str:
    raw = str(tip)
    sp = tip_to_species.get(raw)
    if sp:
        return genus_species(str(sp))
    base = re.sub(r"__Copy\d+$", "", raw)
    return genus_species(base)


def _species_display_name(raw: str) -> str:
    gs = genus_species(str(raw))
    parts = [p for p in gs.split("_") if p]
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    if len(parts) == 1:
        return parts[0]
    return gs.replace("_", " ")


def _species_display_name_multiline(raw: str) -> str:
    label = _species_display_name(raw)
    parts = label.split(" ", 1)
    if len(parts) == 2:
        return f"{parts[0]}\n{parts[1]}"
    return label


def _shade_palette(base_hex: str, n: int) -> list[tuple[float, float, float]]:
    if n <= 0:
        return []
    base = np.array(to_rgb(base_hex), dtype=float)
    if n == 1:
        return [tuple(base)]
    vals = np.linspace(0.40, 1.00, n)
    out: list[tuple[float, float, float]] = []
    for v in vals:
        col = 1.0 - (1.0 - base) * float(v)
        col = np.clip(col, 0.0, 1.0)
        out.append((float(col[0]), float(col[1]), float(col[2])))
    return out


def _normalize_sample_aggregation(raw: object) -> str:
    v = str(raw or "mean").strip().lower()
    if v in {"mean", "median"}:
        return v
    return "mean"


def _aggregate_component_relative(
    component_table: pd.DataFrame,
    full_table: pd.DataFrame,
    sample_cols: list[str],
    sample_aggregation: str,
) -> pd.Series:
    if component_table.empty or full_table.empty or not sample_cols:
        return pd.Series(dtype=float)

    comp = component_table.reindex(columns=sample_cols).fillna(0.0)
    full = full_table.reindex(columns=sample_cols).fillna(0.0)
    denom = pd.to_numeric(full.sum(axis=0), errors="coerce")
    denom = denom.where(denom > 0.0, np.nan)
    rel = comp.div(denom, axis=1).replace([np.inf, -np.inf], np.nan).fillna(0.0)

    if sample_aggregation == "median":
        out = rel.median(axis=1)
    else:
        out = rel.mean(axis=1)
    out = pd.to_numeric(out, errors="coerce").fillna(0.0)
    out = out[out > 0.0]
    return out.sort_values(ascending=False)


def _community_legend_handles(bar_order: list[str] | None = None) -> tuple[list[Patch], list[str]]:
    groups = list(bar_order or BAR_ORDER)
    handles: list[Patch] = []
    labels: list[str] = []
    for group in groups:
        base = BAR_BASE_COLORS.get(group, "#777777")
        handles.append(Patch(facecolor=base, edgecolor=base))
        labels.append("Expected" if group == "Expected" else display_method(group))
    return handles, labels


def _bundle_ymax(bundle: dict[str, object], bar_order: list[str] | None = None) -> float:
    groups = list(bar_order or BAR_ORDER)
    species: list[str] = list(bundle.get("multi_species", []))
    comp_map: dict[str, dict[str, pd.Series]] = dict(bundle.get("species_components", {}))
    ymax = 0.0
    for sp in species:
        bars = comp_map.get(sp, {})
        for group in groups:
            ymax = max(ymax, float(bars.get(group, pd.Series(dtype=float)).sum()))
    return float(ymax)


def _draw_community_panel(
    ax,
    community_id: str,
    bundle: dict[str, object],
    title_prefix: str,
    font_size: float,
    *,
    show_ylabel: bool = True,
    show_yticklabels: bool = True,
    title_text: str | None = None,
    y_top_override: float | None = None,
    species_label_rotation: float = 0.0,
    species_label_multiline: bool = False,
    bar_order: list[str] | None = None,
) -> bool:
    groups = list(bar_order or BAR_ORDER)
    species: list[str] = list(bundle.get("multi_species", []))
    comp_map: dict[str, dict[str, pd.Series]] = dict(bundle.get("species_components", {}))

    if not species:
        return False

    x = np.arange(len(species), dtype=float)
    n_groups = len(groups)
    cluster_width = 0.86
    bar_w = cluster_width / max(1, n_groups)
    offsets = (np.arange(n_groups, dtype=float) - (n_groups - 1) / 2.0) * bar_w

    ymax = _bundle_ymax(bundle, bar_order=groups)

    for gi, group in enumerate(groups):
        base = BAR_BASE_COLORS.get(group, "#777777")
        for si, sp in enumerate(species):
            s = comp_map.get(sp, {}).get(group, pd.Series(dtype=float))
            xpos = x[si] + offsets[gi]
            if s.empty:
                ax.bar(
                    xpos,
                    0.0,
                    width=bar_w * 0.92,
                    facecolor="none",
                    edgecolor=base,
                    linewidth=0.9,
                )
                continue

            left = 0.0
            colors = _shade_palette(base, len(s))
            for ci, (_comp, v) in enumerate(s.items()):
                vv = float(v)
                if vv <= 0.0:
                    continue
                ax.bar(
                    xpos,
                    vv,
                    width=bar_w * 0.92,
                    bottom=left,
                    color=colors[ci],
                    edgecolor=base,
                    linewidth=0.45,
                )
                left += vv

    species_labels = [
        _species_display_name_multiline(sp) if species_label_multiline else _species_display_name(sp)
        for sp in species
    ]
    ax.set_xticks(x)
    if abs(float(species_label_rotation)) < 1e-9:
        x_ha, x_va = ("center", "top")
    else:
        x_ha, x_va = ("right" if float(species_label_rotation) > 0 else "left", "top")
    ax.set_xticklabels(
        species_labels,
        rotation=float(species_label_rotation),
        ha=x_ha,
        va=x_va,
        fontsize=max(7.0, font_size - 2.5),
        rotation_mode="anchor",
    )
    ax.set_ylabel(
        "Relative abundance (%)" if show_ylabel else "",
        fontsize=font_size if show_ylabel else max(7.0, font_size - 1.0),
    )
    ax.yaxis.set_major_formatter(FuncFormatter(lambda yy, _pos: f"{yy * 100:.1f}"))
    ax.tick_params(axis="y", labelsize=max(7.5, font_size - 2.0), labelleft=show_yticklabels)
    ax.grid(True, axis="y", alpha=0.28, linewidth=0.65)
    ax.set_axisbelow(True)

    y_top = float(y_top_override) if y_top_override is not None else max(0.01, min(1.0, ymax) * 1.14)
    ax.set_ylim(0.0, y_top)

    community_disp = display_community_id(community_id)
    title = str(title_text) if title_text is not None else f"{title_prefix} - {community_disp}"
    ax.set_title(title, fontsize=max(9.0, font_size + 0.5))
    return True


def _feature_to_species_mapping(
    mapping_tsv: Path,
    tip_to_species: dict[str, str],
) -> dict[str, str]:
    if not mapping_tsv.exists():
        return {}

    try:
        m = pd.read_csv(mapping_tsv, sep="	")
    except Exception:
        return {}
    if m.empty or "inferred_feature" not in m.columns:
        return {}

    feature = m["inferred_feature"].astype(str)

    if "target_species" in m.columns:
        sp = m["target_species"].astype(str).map(genus_species)
    elif "target_tip" in m.columns:
        sp = m["target_tip"].astype(str).map(lambda t: _species_from_tip(t, tip_to_species))
    elif "best_tip" in m.columns:
        sp = m["best_tip"].astype(str).map(lambda t: _species_from_tip(t, tip_to_species))
    else:
        return {}

    out: dict[str, str] = {}
    for f, s in zip(feature, sp):
        sf = str(f).strip()
        ss = str(s).strip()
        if not sf or not ss or ss.lower() == "nan":
            continue
        out[sf] = genus_species(ss)
    return out


def _collect_multistrain_species_data_for_community(
    outdir: Path,
    community_id: str,
    sample_aggregation: str,
    method_order: list[str] | None = None,
) -> dict[str, object]:
    ed = expected_dir(outdir, community_id)
    exp_tip_p = ed / "expected_tip_abundance_tss1000.tsv"
    tip_to_species_p = ed / "tip_to_species.tsv"

    exp_tip = _read_numeric_table(exp_tip_p)
    tip_to_species = _load_tip_to_species(tip_to_species_p)

    if exp_tip.empty or not tip_to_species:
        return {
            "community_id": community_id,
            "sample_cols": [],
            "multi_species": [],
            "species_components": {},
            "status": "missing_expected",
        }

    species_to_tips: dict[str, list[str]] = {}
    for tip, sp in tip_to_species.items():
        tt = str(tip)
        if tt not in exp_tip.index:
            continue
        gs = genus_species(str(sp))
        species_to_tips.setdefault(gs, []).append(tt)

    multi_species = sorted([sp for sp, tips in species_to_tips.items() if len(set(tips)) > 1])
    if not multi_species:
        return {
            "community_id": community_id,
            "sample_cols": list(exp_tip.columns.astype(str)),
            "multi_species": [],
            "species_components": {},
            "status": "no_multistrain_species",
        }

    expected_cols = [str(c) for c in exp_tip.columns.astype(str).tolist()]
    shared_cols = set(expected_cols)
    method_tables: dict[str, pd.DataFrame] = {}
    feature_species: dict[str, dict[str, str]] = {}
    methods = list(method_order or METHOD_ORDER)

    for method in methods:
        md = method_dir(outdir, community_id, method)
        t = _read_numeric_table(md / "table_std_tss1000.tsv")
        method_tables[method] = t
        if not t.empty:
            shared_cols &= set(t.columns.astype(str).tolist())
        feature_species[method] = _feature_to_species_mapping(md / "mapping.tsv", tip_to_species)

    sample_cols = [c for c in expected_cols if c in shared_cols]
    if not sample_cols:
        sample_cols = expected_cols

    species_components: dict[str, dict[str, pd.Series]] = {}
    for species in multi_species:
        tips = [t for t in species_to_tips.get(species, []) if t in exp_tip.index]
        if not tips:
            continue

        exp_comp = _aggregate_component_relative(
            component_table=exp_tip.loc[tips],
            full_table=exp_tip,
            sample_cols=sample_cols,
            sample_aggregation=sample_aggregation,
        )

        bars: dict[str, pd.Series] = {"Expected": exp_comp}

        for method in methods:
            tab = method_tables.get(method, pd.DataFrame())
            fmap = feature_species.get(method, {})
            if tab.empty or not fmap:
                bars[method] = pd.Series(dtype=float)
                continue

            method_cols = [c for c in sample_cols if c in tab.columns]
            if not method_cols:
                bars[method] = pd.Series(dtype=float)
                continue

            features = [f for f, sp in fmap.items() if sp == species and f in tab.index]
            if not features:
                bars[method] = pd.Series(dtype=float)
                continue

            comp = _aggregate_component_relative(
                component_table=tab.loc[features],
                full_table=tab,
                sample_cols=method_cols,
                sample_aggregation=sample_aggregation,
            )
            bars[method] = comp

        species_components[species] = bars

    multi_species_sorted = sorted(
        species_components.keys(),
        key=lambda sp: float(species_components[sp].get("Expected", pd.Series(dtype=float)).sum()),
        reverse=True,
    )

    return {
        "community_id": community_id,
        "sample_cols": sample_cols,
        "multi_species": multi_species_sorted,
        "species_components": species_components,
        "status": "ok",
    }


def _plot_community_figure(
    community_id: str,
    bundle: dict[str, object],
    out_path: Path,
    title_prefix: str,
    font_size: float,
    sample_aggregation: str,
    bar_order: list[str] | None = None,
) -> None:
    if not bundle.get("multi_species", []):
        return

    species: list[str] = list(bundle.get("multi_species", []))
    fig_w = max(12.0, 3.6 + 1.08 * len(species))
    fig_h = 8.2
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    species_label_rotation, species_label_multiline = _community_species_label_style(community_id)
    _draw_community_panel(
        ax=ax,
        community_id=community_id,
        bundle=bundle,
        title_prefix=title_prefix,
        font_size=font_size,
        show_ylabel=True,
        show_yticklabels=True,
        species_label_rotation=species_label_rotation,
        species_label_multiline=species_label_multiline,
        bar_order=bar_order,
    )

    handles, labels = _community_legend_handles(bar_order=bar_order)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=min(len(labels), 6),
        frameon=False,
        fontsize=max(7.0, font_size - 2.0),
    )

    fig.tight_layout(rect=(0.01, 0.01, 0.995, 0.94))
    fig.savefig(out_path, format="svg", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def _plot_community_grid(
    bundles: list[tuple[str, dict[str, object]]],
    out_path: Path,
    title_prefix: str,
    font_size: float,
    sample_aggregation: str,
    shared_ylim: bool = True,
    bar_order: list[str] | None = None,
) -> None:
    if not bundles:
        return

    n = len(bundles)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / float(ncols)))
    fig_w = max(13.5, 6.8 * ncols)
    fig_h = max(4.2, 3.95 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), sharey=shared_ylim)
    axes_arr = np.atleast_1d(axes).reshape(nrows, ncols)
    grid_ymax = max((_bundle_ymax(bundle, bar_order=bar_order) for _community_id, bundle in bundles), default=0.0)
    grid_y_top = max(0.01, min(1.0, grid_ymax) * 1.05)

    for idx, (community_id, bundle) in enumerate(bundles):
        r = idx // ncols
        c = idx % ncols
        ax = axes_arr[r, c]
        local_y_top = max(0.01, min(1.0, _bundle_ymax(bundle, bar_order=bar_order)) * 1.05)
        species_label_rotation, species_label_multiline = _community_species_label_style(community_id)
        ok = _draw_community_panel(
            ax=ax,
            community_id=community_id,
            bundle=bundle,
            title_prefix=title_prefix,
            font_size=max(8.0, font_size - 1.0),
            show_ylabel=(c == 0),
            show_yticklabels=(c == 0) if shared_ylim else True,
            title_text=display_community_id(community_id),
            y_top_override=grid_y_top if shared_ylim else local_y_top,
            species_label_rotation=species_label_rotation,
            species_label_multiline=species_label_multiline,
            bar_order=bar_order,
        )
        if not ok:
            ax.set_visible(False)

    for idx in range(n, nrows * ncols):
        r = idx // ncols
        c = idx % ncols
        axes_arr[r, c].set_visible(False)

    handles, labels = _community_legend_handles(bar_order=bar_order)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.97),
        ncol=min(len(labels), 6),
        frameon=False,
        fontsize=max(7.0, font_size - 2.0),
    )
    fig.suptitle(title_prefix, fontsize=max(10.0, font_size + 1.0), y=0.985)
    fig.tight_layout(rect=(0.01, 0.01, 0.995, 0.90))
    fig.subplots_adjust(hspace=0.32, wspace=0.18)
    fig.savefig(out_path, format="svg", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def run(cfg) -> int:
    """
    Figure 20:
      Per-community multi-strain species abundance bars.

      For each community containing at least one species with >1 expected strain,
      render grouped stacked bars where each species has:
        - one expected bar (stacked by expected strains/tips)
        - one bar per method (stacked by inferred mapped features for that species)

      Relative abundances are computed sample-wise and aggregated by mean/median
      across shared expected/method sample columns.
    """
    log = get_logger()
    if not bool(getattr(cfg.figures, "fig20", True)):
        log.info("Figure 20 disabled in config; skipping.")
        return 0

    outdir = Path(cfg.paths.outdir)
    manifest_p = outdir / "manifest.tsv"
    if not manifest_p.exists():
        log.info("Figure 20: manifest.tsv not found; run discover first.")
        return 0

    try:
        mf = pd.read_csv(manifest_p, sep="	")
    except Exception:
        log.info("Figure 20: failed reading manifest.tsv.")
        return 0
    if mf.empty or "community_id" not in mf.columns:
        log.info("Figure 20: manifest.tsv has no community_id rows.")
        return 0

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="	")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig20", []) or [])
    configure_community_display_mapping_from_manifest(
        outdir,
        figure_exclude_patterns=exclude_patterns,
    )

    sample_aggregation = _normalize_sample_aggregation(getattr(cfg.figures, "fig20_sample_aggregation", "mean"))
    font_size = float(getattr(cfg.figures, "fig20_font_size", 10.0))
    method_order = _normalize_method_order(getattr(cfg.figures, "fig20_methods", None))
    bar_order = ["Expected"] + method_order
    title_prefix = str(
        getattr(cfg.figures, "fig20_title", None)
        or "Figure 20: Multi-strain species abundance recovery"
    )

    out_fig = figures_root(outdir) / "fig20_multistrain_species_abundance"
    out_fig.mkdir(parents=True, exist_ok=True)
    out_plots = out_fig / "community_plots"
    out_plots.mkdir(parents=True, exist_ok=True)

    component_rows: list[dict[str, object]] = []
    total_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    grid_bundles: list[tuple[str, dict[str, object]]] = []

    community_ids = list(dict.fromkeys(mf["community_id"].astype(str).tolist()))

    n_plotted = 0
    for cid in community_ids:
        if cid in skipped or _is_excluded(cid, exclude_patterns):
            continue

        bundle = _collect_multistrain_species_data_for_community(
            outdir=outdir,
            community_id=cid,
            sample_aggregation=sample_aggregation,
            method_order=method_order,
        )

        species = list(bundle.get("multi_species", []))
        comp_map = dict(bundle.get("species_components", {}))
        n_samples = len(list(bundle.get("sample_cols", [])))
        status = str(bundle.get("status", "ok"))

        plot_path = ""
        if species:
            safe_cid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(cid))
            p = out_plots / f"figure20_multistrain_species_abundance__community-{safe_cid}.svg"
            _plot_community_figure(
                community_id=cid,
                bundle=bundle,
                out_path=p,
                title_prefix=title_prefix,
                font_size=font_size,
                sample_aggregation=sample_aggregation,
                bar_order=bar_order,
            )
            plot_path = str(p)
            n_plotted += 1
            grid_bundles.append((str(cid), bundle))

            for spc in species:
                bars = comp_map.get(spc, {})
                for group in bar_order:
                    s = bars.get(group, pd.Series(dtype=float))
                    total = float(s.sum()) if not s.empty else 0.0
                    n_components = int(len(s)) if not s.empty else 0
                    total_rows.append(
                        {
                            "community_id": str(cid),
                            "community_display": str(display_community_id(cid)),
                            "species": str(spc),
                            "species_display": _species_display_name(spc),
                            "bar_group": str(group),
                            "bar_group_display": "Expected" if group == "Expected" else str(display_method(group)),
                            "total_relative_abundance": float(total),
                            "n_components": int(n_components),
                            "n_samples_aggregated": int(n_samples),
                            "sample_aggregation": str(sample_aggregation),
                            "methods_included": ",".join(method_order),
                        }
                    )

                    if s.empty:
                        component_rows.append(
                            {
                                "community_id": str(cid),
                                "community_display": str(display_community_id(cid)),
                                "species": str(spc),
                                "species_display": _species_display_name(spc),
                                "bar_group": str(group),
                                "bar_group_display": "Expected" if group == "Expected" else str(display_method(group)),
                                "component_id": "",
                                "component_type": "missing",
                                "component_relative_abundance": 0.0,
                                "total_relative_abundance": float(total),
                                "n_components": int(n_components),
                                "n_samples_aggregated": int(n_samples),
                                "sample_aggregation": str(sample_aggregation),
                                "methods_included": ",".join(method_order),
                            }
                        )
                        continue

                    ctype = "expected_tip" if group == "Expected" else "inferred_feature"
                    for comp_id, val in s.items():
                        component_rows.append(
                            {
                                "community_id": str(cid),
                                "community_display": str(display_community_id(cid)),
                                "species": str(spc),
                                "species_display": _species_display_name(spc),
                                "bar_group": str(group),
                                "bar_group_display": "Expected" if group == "Expected" else str(display_method(group)),
                                "component_id": str(comp_id),
                                "component_type": str(ctype),
                                "component_relative_abundance": float(val),
                                "total_relative_abundance": float(total),
                                "n_components": int(n_components),
                                "n_samples_aggregated": int(n_samples),
                                "sample_aggregation": str(sample_aggregation),
                                "methods_included": ",".join(method_order),
                            }
                        )

        summary_rows.append(
            {
                "community_id": str(cid),
                "community_display": str(display_community_id(cid)),
                "n_multi_strain_species": int(len(species)),
                "n_samples_aggregated": int(n_samples),
                "sample_aggregation": str(sample_aggregation),
                "methods_included": ",".join(method_order),
                "status": str(status),
                "plot_path": str(plot_path),
            }
        )

    grid_bundles.sort(key=lambda item: _community_grid_sort_key(item[0]))

    summary_df = pd.DataFrame(summary_rows)
    total_df = pd.DataFrame(total_rows)
    comp_df = pd.DataFrame(component_rows)

    summary_df.to_csv(out_fig / "fig20_multistrain_species_plot_summary.tsv", sep="	", index=False)
    total_df.to_csv(out_fig / "fig20_multistrain_species_bar_totals.tsv", sep="	", index=False)
    comp_df.to_csv(out_fig / "fig20_multistrain_species_components_long.tsv", sep="	", index=False)

    grid_path = out_fig / "fig20_multistrain_species_abundance_grid.svg"
    if grid_bundles:
        _plot_community_grid(
            bundles=grid_bundles,
            out_path=grid_path,
            title_prefix=title_prefix,
            font_size=font_size,
            sample_aggregation=sample_aggregation,
            shared_ylim=True,
            bar_order=bar_order,
        )
        _plot_community_grid(
            bundles=grid_bundles,
            out_path=out_fig / "fig20_multistrain_species_abundance_grid_local_ylim.svg",
            title_prefix=title_prefix,
            font_size=font_size,
            sample_aggregation=sample_aggregation,
            shared_ylim=False,
            bar_order=bar_order,
        )

    log.info(
        f"Wrote Figure 20 to {out_fig} "
        f"(communities_considered={len(summary_rows)}, community_plots={n_plotted}, grid={'yes' if grid_bundles else 'no'})"
    )
    return 0
