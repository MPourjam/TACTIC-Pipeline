from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator, MultipleLocator
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.io import fasta_ids, genus_species
from ..utils.labels import display_community_id, display_method
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from . import plot_fig09

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
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

METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}

FIG16_RELATIVE_XLABEL = "Relative Abundance Deviation (Inferred - Expected)"
FIG16_RAW_XLABEL = "Abundance Deviation (Inferred - Expected)"


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


def _normalize_projection_abundance_mode(raw: object, default: str = "relative") -> str:
    v = str(raw).strip().lower()
    if v in {"source_relative", "original_relative"}:
        return "original_relative"
    if v in {"relative", "raw", "original_relative"}:
        return v
    return default


def _is_relative_projection_mode(projection_mode: str) -> bool:
    return str(projection_mode).strip().lower() in {"relative", "original_relative"}


def _resolve_fig16_xlabel(configured: object, projection_mode: str) -> str:
    if projection_mode == "raw":
        fallback = FIG16_RAW_XLABEL
    else:
        fallback = FIG16_RELATIVE_XLABEL
    txt = str(configured).strip() if configured is not None else ""
    if not txt:
        return fallback
    lower = txt.lower()
    if projection_mode == "raw" and ("relative" in lower or "%" in txt):
        return FIG16_RAW_XLABEL
    if _is_relative_projection_mode(projection_mode) and ("raw" in lower and "relative" not in lower):
        return fallback
    return txt


def _parse_xlim_mode(raw: object) -> tuple[str, float | None]:
    """
    Parse Figure 16 x-axis limit setting.
    Returns:
      - ("auto", None): one common limit from max range across communities
      - ("unbound", None): per-community local limit
      - ("fixed", frac): fixed symmetric limit in fraction scale (0..1)
    """
    if raw is None:
        return ("unbound", None)
    if isinstance(raw, (int, float)):
        val = float(raw)
        if 0.0 <= val <= 1.0:
            return ("fixed", max(0.0, min(1.0, val)))
        return ("fixed", max(0.0, min(1.0, val / 100.0)))
    txt = str(raw).strip().lower()
    if txt in {"auto", "unbound"}:
        return (txt, None)
    try:
        val = float(txt)
        if 0.0 <= val <= 1.0:
            return ("fixed", max(0.0, min(1.0, val)))
        return ("fixed", max(0.0, min(1.0, val / 100.0)))
    except Exception:
        return ("unbound", None)


def _x_tick_alignment(rotation: float) -> tuple[str, str]:
    rot = float(rotation)
    if abs(rot) < 1e-9:
        return ("center", "top")
    return ("right" if rot > 0 else "left", "top")


def _y_tick_alignment(rotation: float) -> tuple[str, str]:
    rot = float(rotation)
    if abs(rot) < 1e-9:
        return ("right", "center")
    return ("right" if rot > 0 else "left", "center")


def _driver_deviation_matrix(
    exp_sub: pd.DataFrame,
    proj_sub: pd.DataFrame,
    projection_mode: str,
    projection_totals: pd.Series | None = None,
    expected_totals: pd.Series | None = None,
) -> pd.DataFrame:
    """
    Compute the deviation matrix in either relative or raw abundance space.
    """
    if projection_mode == "raw":
        return proj_sub - exp_sub

    def _to_relative(df: pd.DataFrame) -> pd.DataFrame:
        col_sum = pd.to_numeric(df.sum(axis=0), errors="coerce")
        col_sum = col_sum.where(col_sum > 0.0, np.nan)
        rel = df.div(col_sum, axis=1)
        return rel.fillna(0.0)

    def _to_relative_with_totals(df: pd.DataFrame, totals: pd.Series | None) -> pd.DataFrame:
        if totals is None:
            return _to_relative(df)
        denom = pd.to_numeric(pd.Series(totals), errors="coerce").reindex(df.columns)
        denom = denom.where(denom > 0.0, np.nan)
        rel = df.div(denom, axis=1)
        return rel.fillna(0.0)

    if projection_mode == "original_relative":
        exp_rel = _to_relative_with_totals(exp_sub, expected_totals)
        proj_rel = _to_relative_with_totals(proj_sub, projection_totals)
        return proj_rel - exp_rel

    exp_rel = _to_relative(exp_sub)
    proj_rel = _to_relative(proj_sub)
    return proj_rel - exp_rel


def _top_species_local_xmax(
    matrix: pd.DataFrame,
    top_n: int,
    projection_mode: str,
    method_order: list[str] | None = None,
) -> float:
    """
    Compute local x-extent from plotted top-N species only.
    """
    if matrix.empty:
        return 0.0
    max_abs = matrix.abs().max(axis=1)
    nonzero = max_abs[max_abs > 0.0]
    if nonzero.empty:
        return 0.0
    top_idx = nonzero.sort_values(ascending=False).head(int(max(1, top_n))).index.tolist()
    methods = [m for m in (method_order or METHOD_ORDER) if m in matrix.columns]
    if not methods:
        return 0.0
    top = matrix.loc[top_idx, methods]
    try:
        x_max = float(np.nanmax(np.abs(top.to_numpy(dtype=float))))
    except Exception:
        return 0.0
    if _is_relative_projection_mode(projection_mode):
        return min(1.0, max(0.02, x_max))
    return max(0.02, x_max)


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _species_display_name(raw: str) -> str:
    gs = genus_species(str(raw))
    parts = [p for p in gs.split("_") if p]
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    if len(parts) == 1:
        return parts[0]
    return gs.replace("_", " ")


def _split_tip_copy(raw: str) -> tuple[str, int | None]:
    tip = str(raw)
    m = re.match(r"^(.*)__Copy(\d+)$", tip)
    if m:
        return str(m.group(1)), int(m.group(2))
    return tip, None


def _tip_display_name(raw: str, reset_copy_number: int | None = None) -> str:
    base, parsed_copy = _split_tip_copy(raw)
    sp = _species_display_name(base)
    if reset_copy_number is not None:
        return f"{sp} | Copy{int(reset_copy_number)}"
    if parsed_copy is not None:
        return f"{sp} | Copy{int(parsed_copy)}"
    return sp


def _build_tip_copy_reset_mapping(
    outdir: Path,
    community_id: str,
    tips: list[str],
) -> tuple[dict[str, int], pd.DataFrame]:
    """
    Build a deterministic per-community copy-number remapping table for Figure 16
    V4 drivers, so displayed copy labels start at 1 within each species.
    """
    ed = expected_dir(outdir, community_id)
    expected_fasta = ed / "expected_sequences.fasta"
    fasta_tip_order: list[str] = []
    if expected_fasta.exists():
        try:
            fasta_tip_order = [str(t) for t in fasta_ids(expected_fasta)]
        except Exception:
            fasta_tip_order = []

    dedup_tips: list[str] = []
    seen_tips: set[str] = set()
    for tip in fasta_tip_order + [str(t) for t in tips]:
        t = str(tip)
        if t in seen_tips:
            continue
        seen_tips.add(t)
        dedup_tips.append(t)

    order_lookup: dict[str, int] = {t: i for i, t in enumerate(dedup_tips, start=1)}

    by_species: dict[str, list[dict[str, object]]] = {}
    for tip in dedup_tips:
        base, parsed_copy = _split_tip_copy(tip)
        species = genus_species(base)
        by_species.setdefault(species, []).append(
            {
                "tip": tip,
                "base": base,
                "parsed_copy": parsed_copy,
                "expected_tip_order": int(order_lookup.get(tip, 0)),
            }
        )

    community_disp = display_community_id(community_id)
    tip_to_reset_copy: dict[str, int] = {}
    mapping_rows: list[dict[str, object]] = []
    for species in sorted(by_species):
        recs = by_species[species]
        recs_sorted = sorted(
            recs,
            key=lambda r: (
                int(r["expected_tip_order"]),
                int(r["parsed_copy"]) if r["parsed_copy"] is not None else -1,
                str(r["tip"]),
            ),
        )
        for reset_copy, rec in enumerate(recs_sorted, start=1):
            tip = str(rec["tip"])
            parsed_copy = rec["parsed_copy"]
            actual_copy = int(parsed_copy) if parsed_copy is not None else 0
            tip_to_reset_copy[tip] = int(reset_copy)
            mapping_rows.append(
                {
                    "community_id": str(community_id),
                    "community_display": str(community_disp),
                    "species": str(species),
                    "species_display": str(_species_display_name(species)),
                    "tip": str(tip),
                    "tip_display_original": str(_tip_display_name(tip)),
                    "tip_display_reset": str(_tip_display_name(tip, reset_copy_number=reset_copy)),
                    "expected_tip_order": int(rec["expected_tip_order"]),
                    "actual_copy_number": int(actual_copy),
                    "actual_copy_label": str(f"Copy{actual_copy}"),
                    "reset_copy_number": int(reset_copy),
                    "reset_copy_label": str(f"Copy{int(reset_copy)}"),
                }
            )

    mapping_df = pd.DataFrame(mapping_rows)
    if not mapping_df.empty:
        mapping_df = mapping_df.sort_values(
            ["community_id", "species_display", "reset_copy_number", "tip"],
            kind="mergesort",
        ).reset_index(drop=True)
    return tip_to_reset_copy, mapping_df


def _read_numeric_table(path: Path) -> pd.DataFrame:
    try:
        t = pd.read_csv(path, sep="\t", index_col=0)
    except Exception:
        return pd.DataFrame()
    if t.empty:
        return pd.DataFrame()
    t = t.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    t.index = t.index.astype(str)
    t.columns = t.columns.astype(str)
    return t


def _load_method_source_table_and_mapping(
    outdir: Path,
    community_id: str,
    method: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    md = method_dir(outdir, community_id, method)
    table_p = md / "table_std_tss1000.tsv"
    map_p = md / "mapping.tsv"
    if not table_p.exists() or not map_p.exists():
        return pd.DataFrame(), pd.DataFrame()
    tab = _read_numeric_table(table_p)
    if tab.empty:
        return pd.DataFrame(), pd.DataFrame()
    try:
        mp = pd.read_csv(map_p, sep="\t")
    except Exception:
        return pd.DataFrame(), pd.DataFrame()
    if mp.empty or "inferred_feature" not in mp.columns:
        return pd.DataFrame(), pd.DataFrame()
    mp = mp.copy()
    mp["inferred_feature"] = mp["inferred_feature"].astype(str)
    mp = mp.drop_duplicates(subset=["inferred_feature"], keep="first")
    return tab, mp


def _load_method_species_projection_from_source(
    outdir: Path,
    community_id: str,
    method: str,
    tip_to_species: dict[str, str],
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Build species x samples projection directly from the original inferred table,
    keeping the original inferred sample totals as the denominator source.
    """
    tab, mp = _load_method_source_table_and_mapping(outdir, community_id, method)
    if tab.empty or mp.empty:
        return pd.DataFrame(), pd.Series(dtype=float)

    source_totals = pd.to_numeric(tab.sum(axis=0), errors="coerce").fillna(0.0)
    species_targets: list[str] = []
    feature_targets: dict[str, str] = {}
    has_target_species = "target_species" in mp.columns
    has_target_tip = "target_tip" in mp.columns
    if not has_target_species and not has_target_tip:
        return pd.DataFrame(), source_totals

    for rec in mp.itertuples(index=False):
        feat = str(getattr(rec, "inferred_feature", "") or "")
        if feat not in tab.index:
            continue
        target_species = str(getattr(rec, "target_species", "") or "") if has_target_species else ""
        target_tip = str(getattr(rec, "target_tip", "") or "") if has_target_tip else ""
        if target_species:
            species = genus_species(target_species)
        elif target_tip:
            species = genus_species(tip_to_species.get(target_tip) or re.sub(r"__Copy\d+$", "", target_tip))
        else:
            continue
        if not species:
            continue
        feature_targets[feat] = str(species)
        species_targets.append(str(species))

    if not species_targets:
        return pd.DataFrame(), source_totals

    out = pd.DataFrame(0.0, index=list(dict.fromkeys(species_targets)), columns=tab.columns)
    for feat, species in feature_targets.items():
        if species in out.index:
            out.loc[species, :] += tab.loc[feat, :].values
    out.index = [genus_species(i) for i in out.index.astype(str)]
    out = out.groupby(level=0, sort=False).sum(numeric_only=True)
    return out, source_totals


def _load_method_tip_projection_from_source(
    outdir: Path,
    community_id: str,
    method: str,
    tip_to_species: dict[str, str],
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Build tip(V4-copy) x samples projection from the original inferred table.
    If the method has no tip-level projected table, preserve current Figure 16
    behavior by aggregating to species first and then splitting equally across
    the expected tips of that species.
    """
    md = method_dir(outdir, community_id, method)
    tip_p = md / "mapped_tip_table_tss1000.tsv"

    tab, mp = _load_method_source_table_and_mapping(outdir, community_id, method)
    if tab.empty or mp.empty:
        return pd.DataFrame(), pd.Series(dtype=float)
    source_totals = pd.to_numeric(tab.sum(axis=0), errors="coerce").fillna(0.0)

    if tip_p.exists() and "target_tip" in mp.columns:
        feature_targets: dict[str, str] = {}
        tip_targets: list[str] = []
        for rec in mp.itertuples(index=False):
            feat = str(getattr(rec, "inferred_feature", "") or "")
            target_tip = str(getattr(rec, "target_tip", "") or "")
            if feat not in tab.index or not target_tip:
                continue
            feature_targets[feat] = target_tip
            tip_targets.append(target_tip)
        if not tip_targets:
            return pd.DataFrame(), source_totals
        out = pd.DataFrame(0.0, index=list(dict.fromkeys(tip_targets)), columns=tab.columns)
        for feat, target_tip in feature_targets.items():
            if target_tip in out.index:
                out.loc[target_tip, :] += tab.loc[feat, :].values
        return out, source_totals

    spp, source_totals = _load_method_species_projection_from_source(
        outdir=outdir,
        community_id=community_id,
        method=method,
        tip_to_species=tip_to_species,
    )
    if spp.empty:
        return pd.DataFrame(), source_totals

    sp_to_tips = _species_to_tips_map(tip_to_species)
    expected_tips = list(dict.fromkeys([str(t) for t in tip_to_species.keys()]))
    if not expected_tips:
        return pd.DataFrame(), source_totals
    out = pd.DataFrame(0.0, index=expected_tips, columns=spp.columns)
    for sp, row in spp.iterrows():
        tips = sp_to_tips.get(genus_species(str(sp)), [])
        if not tips:
            continue
        share = row.to_numpy(dtype=float) / float(len(tips))
        for tip in tips:
            if tip in out.index:
                out.loc[tip, :] += share
    return out, source_totals


def _load_method_species_projection(
    outdir: Path,
    community_id: str,
    method: str,
    tip_to_species: dict[str, str],
) -> pd.DataFrame:
    """
    Return species x samples projected abundance table for one method.
    For denoising methods we collapse mapped-tip projection to species.
    """
    md = method_dir(outdir, community_id, method)
    species_p = md / "mapped_species_table_tss1000.tsv"
    tip_p = md / "mapped_tip_table_tss1000.tsv"

    if species_p.exists():
        t = _read_numeric_table(species_p)
        if t.empty:
            return pd.DataFrame()
        t.index = [genus_species(i) for i in t.index.astype(str)]
        t = t.groupby(level=0, sort=False).sum(numeric_only=True)
        return t

    if not tip_p.exists():
        return pd.DataFrame()

    t = _read_numeric_table(tip_p)
    if t.empty:
        return pd.DataFrame()

    species_idx = []
    for tip in t.index.astype(str):
        sp = tip_to_species.get(tip)
        if not sp:
            sp = genus_species(re.sub(r"__Copy\d+$", "", tip))
        species_idx.append(genus_species(str(sp)))
    t = t.copy()
    t.index = species_idx
    t = t.groupby(level=0, sort=False).sum(numeric_only=True)
    return t


def _species_to_tips_map(tip_to_species: dict[str, str]) -> dict[str, list[str]]:
    sp_to_tips: dict[str, list[str]] = {}
    for tip, sp in tip_to_species.items():
        sp_key = genus_species(str(sp))
        if sp_key not in sp_to_tips:
            sp_to_tips[sp_key] = []
        sp_to_tips[sp_key].append(str(tip))
    return sp_to_tips


def _load_method_tip_projection(
    outdir: Path,
    community_id: str,
    method: str,
    tip_to_species: dict[str, str],
) -> pd.DataFrame:
    """
    Return tip(V4-copy) x samples projected abundance table for one method.
    If only species-level projection exists, split species mass equally across
    expected tips within that species.
    """
    md = method_dir(outdir, community_id, method)
    tip_p = md / "mapped_tip_table_tss1000.tsv"
    species_p = md / "mapped_species_table_tss1000.tsv"

    if tip_p.exists():
        t = _read_numeric_table(tip_p)
        if t.empty:
            return pd.DataFrame()
        return t

    if not species_p.exists():
        return pd.DataFrame()

    spp = _read_numeric_table(species_p)
    if spp.empty:
        return pd.DataFrame()
    spp.index = [genus_species(i) for i in spp.index.astype(str)]
    spp = spp.groupby(level=0, sort=False).sum(numeric_only=True)

    sp_to_tips = _species_to_tips_map(tip_to_species)
    expected_tips = list(dict.fromkeys([str(t) for t in tip_to_species.keys()]))
    if not expected_tips:
        return pd.DataFrame()
    out = pd.DataFrame(0.0, index=expected_tips, columns=spp.columns)
    for sp, row in spp.iterrows():
        tips = sp_to_tips.get(genus_species(str(sp)), [])
        if not tips:
            continue
        share = row.to_numpy(dtype=float) / float(len(tips))
        for tip in tips:
            if tip in out.index:
                out.loc[tip, :] += share
    return out


def _build_species_driver_matrix_for_community(
    outdir: Path,
    community_id: str,
    sample_aggregation: str = "mean",
    projection_mode: str = "relative",
    method_order: list[str] | None = None,
) -> tuple[pd.DataFrame, int]:
    """
    Build species x methods signed deviation matrix where:
      deviation = projected_abundance - expected_abundance
    in either relative or raw abundance space.
    Deviations are aggregated across samples by mean/median.
    """
    exp_p = expected_dir(outdir, community_id) / "expected_species_abundance_tss1000.tsv"
    expected = _read_numeric_table(exp_p)
    if expected.empty:
        return pd.DataFrame(), 0
    expected.index = [genus_species(i) for i in expected.index.astype(str)]
    expected = expected.groupby(level=0, sort=False).sum(numeric_only=True)

    tip_to_species = plot_fig09._load_tip_to_species(expected_dir(outdir, community_id) / "tip_to_species.tsv")
    samples_used: set[str] = set()
    per_method: dict[str, pd.Series] = {}
    all_species: set[str] = set(expected.index.astype(str).tolist())
    expected_totals = pd.to_numeric(expected.sum(axis=0), errors="coerce").fillna(0.0)

    agg_mode = str(sample_aggregation or "mean").strip().lower()
    if agg_mode not in {"mean", "median"}:
        agg_mode = "mean"

    methods = list(method_order or METHOD_ORDER)

    for method in methods:
        source_totals = pd.Series(dtype=float)
        if projection_mode == "original_relative":
            proj, source_totals = _load_method_species_projection_from_source(
                outdir, community_id, method, tip_to_species
            )
        else:
            proj = _load_method_species_projection(outdir, community_id, method, tip_to_species)
        if proj.empty:
            per_method[method] = pd.Series(dtype=float)
            continue
        all_species.update(proj.index.astype(str).tolist())
        shared_samples = [s for s in expected.columns.astype(str).tolist() if s in proj.columns.astype(str).tolist()]
        if not shared_samples:
            per_method[method] = pd.Series(dtype=float)
            continue
        samples_used.update(shared_samples)
        sp_all = sorted(set(expected.index.astype(str).tolist()) | set(proj.index.astype(str).tolist()))
        exp_sub = expected.reindex(sp_all).fillna(0.0)[shared_samples]
        proj_sub = proj.reindex(sp_all).fillna(0.0)[shared_samples]
        dev = _driver_deviation_matrix(
            exp_sub=exp_sub,
            proj_sub=proj_sub,
            projection_mode=projection_mode,
            projection_totals=source_totals.reindex(shared_samples) if projection_mode == "original_relative" else None,
            expected_totals=expected_totals.reindex(shared_samples) if projection_mode == "original_relative" else None,
        )
        if agg_mode == "median":
            dev_s = dev.median(axis=1)
        else:
            dev_s = dev.mean(axis=1)
        per_method[method] = dev_s.astype(float)

    all_species_sorted = sorted(all_species)
    mat = pd.DataFrame(index=all_species_sorted, columns=methods, dtype=float)
    for method in methods:
        s = per_method.get(method)
        if s is None or s.empty:
            mat[method] = 0.0
            continue
        mat[method] = pd.to_numeric(s.reindex(all_species_sorted), errors="coerce").fillna(0.0).astype(float)
    return mat, len(samples_used)


def _build_tip_driver_matrix_for_community(
    outdir: Path,
    community_id: str,
    sample_aggregation: str = "mean",
    projection_mode: str = "relative",
    method_order: list[str] | None = None,
) -> tuple[pd.DataFrame, int]:
    """
    Build tip(V4-copy) x methods signed deviation matrix where:
      deviation = projected_abundance - expected_abundance
    in either relative or raw abundance space.
    Deviations are aggregated across samples by mean/median.
    """
    exp_p = expected_dir(outdir, community_id) / "expected_tip_abundance_tss1000.tsv"
    expected = _read_numeric_table(exp_p)
    if expected.empty:
        return pd.DataFrame(), 0

    tip_to_species = plot_fig09._load_tip_to_species(expected_dir(outdir, community_id) / "tip_to_species.tsv")
    samples_used: set[str] = set()
    per_method: dict[str, pd.Series] = {}
    all_tips: set[str] = set(expected.index.astype(str).tolist())
    expected_totals = pd.to_numeric(expected.sum(axis=0), errors="coerce").fillna(0.0)

    agg_mode = str(sample_aggregation or "mean").strip().lower()
    if agg_mode not in {"mean", "median"}:
        agg_mode = "mean"

    methods = list(method_order or METHOD_ORDER)

    for method in methods:
        source_totals = pd.Series(dtype=float)
        if projection_mode == "original_relative":
            proj, source_totals = _load_method_tip_projection_from_source(
                outdir, community_id, method, tip_to_species
            )
        else:
            proj = _load_method_tip_projection(outdir, community_id, method, tip_to_species)
        if proj.empty:
            per_method[method] = pd.Series(dtype=float)
            continue
        all_tips.update(proj.index.astype(str).tolist())
        shared_samples = [s for s in expected.columns.astype(str).tolist() if s in proj.columns.astype(str).tolist()]
        if not shared_samples:
            per_method[method] = pd.Series(dtype=float)
            continue
        samples_used.update(shared_samples)
        tips_all = sorted(set(expected.index.astype(str).tolist()) | set(proj.index.astype(str).tolist()))
        exp_sub = expected.reindex(tips_all).fillna(0.0)[shared_samples]
        proj_sub = proj.reindex(tips_all).fillna(0.0)[shared_samples]
        dev = _driver_deviation_matrix(
            exp_sub=exp_sub,
            proj_sub=proj_sub,
            projection_mode=projection_mode,
            projection_totals=source_totals.reindex(shared_samples) if projection_mode == "original_relative" else None,
            expected_totals=expected_totals.reindex(shared_samples) if projection_mode == "original_relative" else None,
        )
        if agg_mode == "median":
            dev_s = dev.median(axis=1)
        else:
            dev_s = dev.mean(axis=1)
        per_method[method] = dev_s.astype(float)

    all_tips_sorted = sorted(all_tips)
    mat = pd.DataFrame(index=all_tips_sorted, columns=methods, dtype=float)
    for method in methods:
        s = per_method.get(method)
        if s is None or s.empty:
            mat[method] = 0.0
            continue
        mat[method] = pd.to_numeric(s.reindex(all_tips_sorted), errors="coerce").fillna(0.0).astype(float)
    return mat, len(samples_used)


def _plot_species_driver_per_community(
    matrix: pd.DataFrame,
    community_id: str,
    n_samples: int,
    out_path: Path,
    font_size: float,
    projection_mode: str = "relative",
    top_n: int = 25,
    x_limit_override: float | None = None,
    x_tick_rotation: float = 45.0,
    y_tick_rotation: float = 45.0,
    title_prefix: str = "Species Driver Plot",
    xlabel: str = "Relative Abundance Deviation (Inferred - Expected)",
    ylabel_template: str = "Top {n} Species by |Deviation|",
    method_order: list[str] | None = None,
) -> pd.DataFrame:
    methods = [m for m in (method_order or METHOD_ORDER) if m in matrix.columns]
    if matrix.empty:
        return pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "species",
                "species_display",
                "method",
                "method_display",
                "deviation_relative_abundance",
                "abs_deviation_relative_abundance",
                "deviation_abundance",
                "abs_deviation_abundance",
                "projection_abundance_mode",
                "rank_by_max_abs_deviation",
                "n_samples_aggregated",
            ]
        )

    max_abs = matrix.abs().max(axis=1)
    nonzero = max_abs[max_abs > 0.0]
    if nonzero.empty:
        return pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "species",
                "species_display",
                "method",
                "method_display",
                "deviation_relative_abundance",
                "abs_deviation_relative_abundance",
                "deviation_abundance",
                "abs_deviation_abundance",
                "projection_abundance_mode",
                "rank_by_max_abs_deviation",
                "n_samples_aggregated",
            ]
        )

    top_idx = nonzero.sort_values(ascending=False).head(int(max(1, top_n))).index.tolist()
    if not methods:
        return pd.DataFrame()

    top = matrix.loc[top_idx, methods].copy()
    top["__max_abs__"] = top.abs().max(axis=1)
    top = top.sort_values("__max_abs__", ascending=False, kind="mergesort")
    top = top.drop(columns=["__max_abs__"])

    plotted_rows: list[dict[str, object]] = []
    community_disp = display_community_id(community_id)
    for rank, species in enumerate(top.index.astype(str).tolist(), start=1):
        sp_disp = _species_display_name(species)
        for method in methods:
            val = float(top.at[species, method])
            plotted_rows.append(
                {
                    "community_id": str(community_id),
                    "community_display": str(community_disp),
                    "species": str(species),
                    "species_display": str(sp_disp),
                    "method": str(method),
                    "method_display": str(display_method(method)),
                    "deviation_relative_abundance": float(val),
                    "abs_deviation_relative_abundance": float(abs(val)),
                    # Backward-compatible aliases; values are relative-abundance deviations.
                    "deviation_abundance": float(val),
                    "abs_deviation_abundance": float(abs(val)),
                    "projection_abundance_mode": str(projection_mode),
                    "rank_by_max_abs_deviation": int(rank),
                    "n_samples_aggregated": int(n_samples),
                }
            )

    species_labels = [_species_display_name(s) for s in top.index.astype(str).tolist()]
    y = np.arange(len(species_labels), dtype=float)
    n_methods = len(methods)
    group_height = 0.82
    bar_h = group_height / max(1, n_methods)
    offsets = (np.arange(n_methods, dtype=float) - (n_methods - 1) / 2.0) * bar_h

    fig_h = max(8.0, 2.4 + 0.36 * len(species_labels))
    fig, ax = plt.subplots(1, 1, figsize=(12.8, fig_h))

    for i, method in enumerate(methods):
        vals = top[method].to_numpy(dtype=float)
        ax.barh(
            y + offsets[i],
            vals,
            height=bar_h * 0.92,
            color=METHOD_COLORS.get(method, "#777777"),
            edgecolor="#333333",
            linewidth=0.35,
            label=display_method(method),
        )

    if x_limit_override is not None:
        x_max = float(x_limit_override)
    else:
        x_max = float(np.nanmax(np.abs(top.to_numpy(dtype=float))))
        if _is_relative_projection_mode(projection_mode):
            x_max = min(1.0, max(0.02, x_max))
        else:
            x_max = max(0.02, x_max)
    pad = x_max * 0.08
    ax.set_xlim(-(x_max + pad), x_max + pad)
    ax.axvline(0.0, color="#222222", linewidth=1.0)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)

    tick_size = max(7.0, font_size - 2.0)
    y_tick_size = max(6.5, font_size - 3.0)
    ax.set_yticks(y)
    ax.set_yticklabels(species_labels, fontsize=y_tick_size)
    ax.invert_yaxis()
    ax.tick_params(axis="y", pad=6)
    y_ha, y_va = _y_tick_alignment(y_tick_rotation)
    plt.setp(
        ax.get_yticklabels(),
        rotation=float(y_tick_rotation),
        ha=y_ha,
        va=y_va,
        rotation_mode="anchor",
    )
    ax.tick_params(axis="x", labelsize=tick_size)
    ax.tick_params(axis="x", pad=6)
    if _is_relative_projection_mode(projection_mode):
        ax.xaxis.set_major_locator(MultipleLocator(0.10))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x:.2f}"))
    else:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x:.1f}"))
    x_ha, x_va = _x_tick_alignment(x_tick_rotation)
    plt.setp(
        ax.get_xticklabels(),
        rotation=float(x_tick_rotation),
        ha=x_ha,
        va=x_va,
        rotation_mode="anchor",
    )
    ax.set_xlabel(str(xlabel), fontsize=font_size)
    try:
        ylab = str(ylabel_template).format(n=len(species_labels), n_species=len(species_labels))
    except Exception:
        ylab = str(ylabel_template)
    ax.set_ylabel(ylab, fontsize=font_size)
    ax.set_title(f"{str(title_prefix)} - {community_disp}", fontsize=max(font_size + 0.5, 9.5))

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.01),
            ncol=5,
            frameon=False,
            fontsize=tick_size,
        )
    fig.tight_layout(rect=(0.03, 0.08, 0.99, 0.98))
    fig.savefig(out_path, format="svg", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    return pd.DataFrame(plotted_rows)


def _plot_tip_driver_per_community(
    matrix: pd.DataFrame,
    community_id: str,
    n_samples: int,
    out_path: Path,
    font_size: float,
    tip_copy_reset_map: dict[str, int] | None = None,
    projection_mode: str = "relative",
    top_n: int = 25,
    x_limit_override: float | None = None,
    x_tick_rotation: float = 45.0,
    y_tick_rotation: float = 45.0,
    title_prefix: str = "V4-Copy Driver Plot",
    xlabel: str = "Relative Abundance Deviation (Inferred - Expected)",
    ylabel_template: str = "Top {n} V4 Copies by |Deviation|",
    method_order: list[str] | None = None,
) -> pd.DataFrame:
    methods = [m for m in (method_order or METHOD_ORDER) if m in matrix.columns]
    if matrix.empty:
        return pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "tip",
                "tip_display",
                "species",
                "species_display",
                "method",
                "method_display",
                "deviation_relative_abundance",
                "abs_deviation_relative_abundance",
                "deviation_abundance",
                "abs_deviation_abundance",
                "projection_abundance_mode",
                "actual_copy_number",
                "reset_copy_number",
                "rank_by_max_abs_deviation",
                "n_samples_aggregated",
            ]
        )

    max_abs = matrix.abs().max(axis=1)
    nonzero = max_abs[max_abs > 0.0]
    if nonzero.empty:
        return pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "tip",
                "tip_display",
                "species",
                "species_display",
                "method",
                "method_display",
                "deviation_relative_abundance",
                "abs_deviation_relative_abundance",
                "deviation_abundance",
                "abs_deviation_abundance",
                "projection_abundance_mode",
                "actual_copy_number",
                "reset_copy_number",
                "rank_by_max_abs_deviation",
                "n_samples_aggregated",
            ]
        )

    top_idx = nonzero.sort_values(ascending=False).head(int(max(1, top_n))).index.tolist()
    if not methods:
        return pd.DataFrame()

    top = matrix.loc[top_idx, methods].copy()
    top["__max_abs__"] = top.abs().max(axis=1)
    top = top.sort_values("__max_abs__", ascending=False, kind="mergesort")
    rank_by_max = {str(t): i for i, t in enumerate(top.index.astype(str).tolist(), start=1)}
    top["__species__"] = [genus_species(_split_tip_copy(t)[0]) for t in top.index.astype(str).tolist()]
    top["__species_display__"] = [_species_display_name(sp) for sp in top["__species__"].astype(str).tolist()]
    top["__reset_copy__"] = [
        int((tip_copy_reset_map or {}).get(str(t), 0)) for t in top.index.astype(str).tolist()
    ]
    top["__tip_display__"] = [
        _tip_display_name(str(t), reset_copy_number=(tip_copy_reset_map or {}).get(str(t)))
        for t in top.index.astype(str).tolist()
    ]
    # Keep top-N selection by deviation, then alphabetize within selected tips so copies
    # from the same species are adjacent and easy to scan.
    top = top.sort_values(
        ["__species_display__", "__reset_copy__", "__tip_display__", "__max_abs__"],
        ascending=[True, True, True, False],
        kind="mergesort",
    )
    top = top.drop(columns=["__max_abs__"])
    top_numeric = top[methods].copy()

    plotted_rows: list[dict[str, object]] = []
    community_disp = display_community_id(community_id)
    for tip in top.index.astype(str).tolist():
        base, parsed_copy = _split_tip_copy(tip)
        tip_reset_copy = (tip_copy_reset_map or {}).get(str(tip))
        tip_disp = _tip_display_name(tip, reset_copy_number=tip_reset_copy)
        species = genus_species(base)
        sp_disp = _species_display_name(species)
        actual_copy_number = int(parsed_copy) if parsed_copy is not None else 0
        reset_copy_number = int(tip_reset_copy) if tip_reset_copy is not None else int(actual_copy_number + 1)
        rank = int(rank_by_max.get(str(tip), 0))
        for method in methods:
            val = float(top_numeric.at[tip, method])
            plotted_rows.append(
                {
                    "community_id": str(community_id),
                    "community_display": str(community_disp),
                    "tip": str(tip),
                    "tip_display": str(tip_disp),
                    "species": str(species),
                    "species_display": str(sp_disp),
                    "method": str(method),
                    "method_display": str(display_method(method)),
                    "deviation_relative_abundance": float(val),
                    "abs_deviation_relative_abundance": float(abs(val)),
                    "deviation_abundance": float(val),
                    "abs_deviation_abundance": float(abs(val)),
                    "projection_abundance_mode": str(projection_mode),
                    "actual_copy_number": int(actual_copy_number),
                    "reset_copy_number": int(reset_copy_number),
                    "rank_by_max_abs_deviation": int(rank),
                    "n_samples_aggregated": int(n_samples),
                }
            )

    tip_labels = [
        _tip_display_name(str(t), reset_copy_number=(tip_copy_reset_map or {}).get(str(t)))
        for t in top.index.astype(str).tolist()
    ]
    y = np.arange(len(tip_labels), dtype=float)
    n_methods = len(methods)
    group_height = 0.82
    bar_h = group_height / max(1, n_methods)
    offsets = (np.arange(n_methods, dtype=float) - (n_methods - 1) / 2.0) * bar_h

    fig_h = max(8.0, 2.4 + 0.36 * len(tip_labels))
    fig, ax = plt.subplots(1, 1, figsize=(12.8, fig_h))

    for i, method in enumerate(methods):
        vals = top_numeric[method].to_numpy(dtype=float)
        ax.barh(
            y + offsets[i],
            vals,
            height=bar_h * 0.92,
            color=METHOD_COLORS.get(method, "#777777"),
            edgecolor="#333333",
            linewidth=0.35,
            label=display_method(method),
        )

    if x_limit_override is not None:
        x_max = float(x_limit_override)
    else:
        x_max = float(np.nanmax(np.abs(top_numeric.to_numpy(dtype=float))))
        if _is_relative_projection_mode(projection_mode):
            x_max = min(1.0, max(0.02, x_max))
        else:
            x_max = max(0.02, x_max)
    pad = x_max * 0.08
    ax.set_xlim(-(x_max + pad), x_max + pad)
    ax.axvline(0.0, color="#222222", linewidth=1.0)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)

    tick_size = max(7.0, font_size - 2.0)
    y_tick_size = max(6.5, font_size - 3.0)
    ax.set_yticks(y)
    ax.set_yticklabels(tip_labels, fontsize=y_tick_size)
    ax.invert_yaxis()
    ax.tick_params(axis="y", pad=6)
    y_ha, y_va = _y_tick_alignment(y_tick_rotation)
    plt.setp(
        ax.get_yticklabels(),
        rotation=float(y_tick_rotation),
        ha=y_ha,
        va=y_va,
        rotation_mode="anchor",
    )
    ax.tick_params(axis="x", labelsize=tick_size)
    ax.tick_params(axis="x", pad=6)
    if _is_relative_projection_mode(projection_mode):
        ax.xaxis.set_major_locator(MultipleLocator(0.10))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x:.2f}"))
    else:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x:.1f}"))
    x_ha, x_va = _x_tick_alignment(x_tick_rotation)
    plt.setp(
        ax.get_xticklabels(),
        rotation=float(x_tick_rotation),
        ha=x_ha,
        va=x_va,
        rotation_mode="anchor",
    )
    ax.set_xlabel(str(xlabel), fontsize=font_size)
    try:
        ylab = str(ylabel_template).format(n=len(tip_labels), n_v4=len(tip_labels))
    except Exception:
        ylab = str(ylabel_template)
    ax.set_ylabel(ylab, fontsize=font_size)
    ax.set_title(f"{str(title_prefix)} - {community_disp}", fontsize=max(font_size + 0.5, 9.5))

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.01),
            ncol=5,
            frameon=False,
            fontsize=tick_size,
        )
    fig.tight_layout(rect=(0.03, 0.08, 0.99, 0.98))
    fig.savefig(out_path, format="svg", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    return pd.DataFrame(plotted_rows)


def _write_species_driver_outputs(
    cfg,
    outdir: Path,
    out_dir: Path,
    included_communities: list[str],
    font_size: float,
    method_order: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    projection_mode = _normalize_projection_abundance_mode(
        getattr(cfg.figures, "fig16_projection_abundance_mode", "relative")
    )
    top_n = int(max(1, getattr(cfg.figures, "fig16_species_driver_top_n", 25)))
    agg_mode = str(getattr(cfg.figures, "fig16_species_driver_sample_aggregation", "mean") or "mean")
    agg_mode = agg_mode.strip().lower()
    if agg_mode not in {"mean", "median"}:
        agg_mode = "mean"
    title_prefix = str(
        getattr(cfg.figures, "fig16_species_driver_title_prefix", None)
        or getattr(cfg.figures, "fig16_title", None)
        or "Species Driver Plot"
    )
    xlabel = str(
        _resolve_fig16_xlabel(
            getattr(cfg.figures, "fig16_species_driver_xlabel", None),
            projection_mode,
        )
    )
    ylabel_template = str(
        getattr(
            cfg.figures,
            "fig16_species_driver_ylabel_template",
            "Top {n} Species by |Deviation|",
        )
        or "Top {n} Species by |Deviation|"
    )
    x_tick_rotation = float(getattr(cfg.figures, "fig16_xtick_rotation", 45.0))
    y_tick_rotation = float(getattr(cfg.figures, "fig16_ytick_rotation", 45.0))
    xlim_mode, xlim_fixed = _parse_xlim_mode(getattr(cfg.figures, "fig16_species_driver_xlim", "unbound"))

    driver_dir = out_dir / "species_driver_plots_by_community"
    driver_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    matrices_by_community: dict[str, tuple[pd.DataFrame, int]] = {}
    for cid in included_communities:
        mat, n_samples = _build_species_driver_matrix_for_community(
            outdir=outdir,
            community_id=cid,
            sample_aggregation=agg_mode,
            projection_mode=projection_mode,
            method_order=method_order,
        )
        matrices_by_community[str(cid)] = (mat, int(n_samples))

    global_xmax = 0.0
    if xlim_mode == "auto":
        for cid in included_communities:
            mat, n_samples = matrices_by_community.get(str(cid), (pd.DataFrame(), 0))
            if mat.empty or n_samples <= 0:
                continue
            global_xmax = max(
                global_xmax,
                _top_species_local_xmax(
                    mat,
                    top_n=top_n,
                    projection_mode=projection_mode,
                    method_order=method_order,
                ),
            )
        global_xmax = float(max(0.02, float(global_xmax)))

    for cid in included_communities:
        mat, n_samples = matrices_by_community.get(str(cid), (pd.DataFrame(), 0))
        if mat.empty or n_samples <= 0:
            summary_rows.append(
                {
                    "community_id": str(cid),
                    "community_display": display_community_id(cid),
                    "n_samples_aggregated": int(n_samples),
                    "n_species_plotted": 0,
                    "sample_aggregation": agg_mode,
                    "projection_abundance_mode": projection_mode,
                    "methods_included": ",".join(method_order),
                    "plot_path": "",
                }
            )
            continue

        safe_cid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(cid))
        out_path = driver_dir / f"figure16_species_driver__community-{safe_cid}.svg"
        if xlim_mode == "fixed":
            x_limit_override = xlim_fixed
        elif xlim_mode == "auto":
            x_limit_override = global_xmax
        else:
            x_limit_override = None
        rows_df = _plot_species_driver_per_community(
            matrix=mat,
            community_id=cid,
            n_samples=n_samples,
            out_path=out_path,
            font_size=font_size,
            projection_mode=projection_mode,
            top_n=top_n,
            x_limit_override=x_limit_override,
            x_tick_rotation=x_tick_rotation,
            y_tick_rotation=y_tick_rotation,
            title_prefix=title_prefix,
            xlabel=xlabel,
            ylabel_template=ylabel_template,
            method_order=method_order,
        )
        if rows_df.empty:
            summary_rows.append(
                {
                    "community_id": str(cid),
                    "community_display": display_community_id(cid),
                    "n_samples_aggregated": int(n_samples),
                    "n_species_plotted": 0,
                    "sample_aggregation": agg_mode,
                    "projection_abundance_mode": projection_mode,
                    "methods_included": ",".join(method_order),
                    "plot_path": "",
                }
            )
            continue

        all_rows.append(rows_df)
        summary_rows.append(
            {
                "community_id": str(cid),
                "community_display": display_community_id(cid),
                "n_samples_aggregated": int(n_samples),
                "n_species_plotted": int(rows_df["species"].astype(str).nunique()),
                "sample_aggregation": agg_mode,
                "projection_abundance_mode": projection_mode,
                "methods_included": ",".join(method_order),
                "plot_path": str(out_path),
            }
        )

    if all_rows:
        driver_df = pd.concat(all_rows, ignore_index=True)
    else:
        driver_df = pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "species",
                "species_display",
                "method",
                "method_display",
                "deviation_relative_abundance",
                "abs_deviation_relative_abundance",
                "deviation_abundance",
                "abs_deviation_abundance",
                "projection_abundance_mode",
                "rank_by_max_abs_deviation",
                "n_samples_aggregated",
            ]
        )
    summary_df = pd.DataFrame(summary_rows)
    return driver_df, summary_df


def _write_tip_driver_outputs(
    cfg,
    outdir: Path,
    out_dir: Path,
    included_communities: list[str],
    font_size: float,
    method_order: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    projection_mode = _normalize_projection_abundance_mode(
        getattr(cfg.figures, "fig16_projection_abundance_mode", "relative")
    )
    top_n = int(max(1, getattr(cfg.figures, "fig16_species_driver_top_n", 25)))
    agg_mode = str(getattr(cfg.figures, "fig16_species_driver_sample_aggregation", "mean") or "mean")
    agg_mode = agg_mode.strip().lower()
    if agg_mode not in {"mean", "median"}:
        agg_mode = "mean"
    title_prefix = str(getattr(cfg.figures, "fig16_v4_driver_title_prefix", "V4-Copy Driver Plot") or "V4-Copy Driver Plot")
    xlabel = str(
        _resolve_fig16_xlabel(
            getattr(cfg.figures, "fig16_v4_driver_xlabel", None),
            projection_mode,
        )
    )
    ylabel_template = str(
        getattr(
            cfg.figures,
            "fig16_v4_driver_ylabel_template",
            "Top {n} V4 Copies by |Deviation|",
        )
        or "Top {n} V4 Copies by |Deviation|"
    )
    x_tick_rotation = float(getattr(cfg.figures, "fig16_xtick_rotation", 45.0))
    y_tick_rotation = float(getattr(cfg.figures, "fig16_ytick_rotation", 45.0))
    xlim_mode, xlim_fixed = _parse_xlim_mode(getattr(cfg.figures, "fig16_species_driver_xlim", "unbound"))

    driver_dir = out_dir / "v4_copy_driver_plots_by_community"
    driver_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    mapping_rows: list[pd.DataFrame] = []
    matrices_by_community: dict[str, tuple[pd.DataFrame, int]] = {}
    for cid in included_communities:
        mat, n_samples = _build_tip_driver_matrix_for_community(
            outdir=outdir,
            community_id=cid,
            sample_aggregation=agg_mode,
            projection_mode=projection_mode,
            method_order=method_order,
        )
        matrices_by_community[str(cid)] = (mat, int(n_samples))

    global_xmax = 0.0
    if xlim_mode == "auto":
        for cid in included_communities:
            mat, n_samples = matrices_by_community.get(str(cid), (pd.DataFrame(), 0))
            if mat.empty or n_samples <= 0:
                continue
            global_xmax = max(
                global_xmax,
                _top_species_local_xmax(
                    mat,
                    top_n=top_n,
                    projection_mode=projection_mode,
                    method_order=method_order,
                ),
            )
        global_xmax = float(max(0.02, float(global_xmax)))

    for cid in included_communities:
        mat, n_samples = matrices_by_community.get(str(cid), (pd.DataFrame(), 0))
        if mat.empty or n_samples <= 0:
            summary_rows.append(
                {
                    "community_id": str(cid),
                    "community_display": display_community_id(cid),
                    "n_samples_aggregated": int(n_samples),
                    "n_v4_copies_plotted": 0,
                    "sample_aggregation": agg_mode,
                    "projection_abundance_mode": projection_mode,
                    "methods_included": ",".join(method_order),
                    "plot_path": "",
                }
            )
            continue

        safe_cid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(cid))
        out_path = driver_dir / f"figure16_v4_copy_driver__community-{safe_cid}.svg"
        tip_copy_reset_map, tip_copy_map_df = _build_tip_copy_reset_mapping(
            outdir=outdir,
            community_id=str(cid),
            tips=mat.index.astype(str).tolist(),
        )
        if not tip_copy_map_df.empty:
            mapping_rows.append(tip_copy_map_df)
        if xlim_mode == "fixed":
            x_limit_override = xlim_fixed
        elif xlim_mode == "auto":
            x_limit_override = global_xmax
        else:
            x_limit_override = None
        rows_df = _plot_tip_driver_per_community(
            matrix=mat,
            community_id=cid,
            n_samples=n_samples,
            out_path=out_path,
            font_size=font_size,
            tip_copy_reset_map=tip_copy_reset_map,
            projection_mode=projection_mode,
            top_n=top_n,
            x_limit_override=x_limit_override,
            x_tick_rotation=x_tick_rotation,
            y_tick_rotation=y_tick_rotation,
            title_prefix=title_prefix,
            xlabel=xlabel,
            ylabel_template=ylabel_template,
            method_order=method_order,
        )
        if rows_df.empty:
            summary_rows.append(
                {
                    "community_id": str(cid),
                    "community_display": display_community_id(cid),
                    "n_samples_aggregated": int(n_samples),
                    "n_v4_copies_plotted": 0,
                    "sample_aggregation": agg_mode,
                    "projection_abundance_mode": projection_mode,
                    "methods_included": ",".join(method_order),
                    "plot_path": "",
                }
            )
            continue

        all_rows.append(rows_df)
        summary_rows.append(
            {
                "community_id": str(cid),
                "community_display": display_community_id(cid),
                "n_samples_aggregated": int(n_samples),
                "n_v4_copies_plotted": int(rows_df["tip"].astype(str).nunique()),
                "sample_aggregation": agg_mode,
                "projection_abundance_mode": projection_mode,
                "methods_included": ",".join(method_order),
                "plot_path": str(out_path),
            }
        )

    if all_rows:
        driver_df = pd.concat(all_rows, ignore_index=True)
    else:
        driver_df = pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "tip",
                "tip_display",
                "species",
                "species_display",
                "method",
                "method_display",
                "deviation_relative_abundance",
                "abs_deviation_relative_abundance",
                "deviation_abundance",
                "abs_deviation_abundance",
                "projection_abundance_mode",
                "actual_copy_number",
                "reset_copy_number",
                "rank_by_max_abs_deviation",
                "n_samples_aggregated",
            ]
        )
    summary_df = pd.DataFrame(summary_rows)
    if mapping_rows:
        mapping_df = pd.concat(mapping_rows, ignore_index=True)
    else:
        mapping_df = pd.DataFrame(
            columns=[
                "community_id",
                "community_display",
                "species",
                "species_display",
                "tip",
                "tip_display_original",
                "tip_display_reset",
                "expected_tip_order",
                "actual_copy_number",
                "actual_copy_label",
                "reset_copy_number",
                "reset_copy_label",
            ]
        )
    return driver_df, summary_df, mapping_df


def run(cfg) -> int:
    """
    Figure 16:
      Per-community species and V4-copy driver plots of abundance deviation
      (projected - expected), with one horizontal grouped bar plot per community.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 16: manifest.tsv not found; run discover first.")
        return 0
    mf = pd.read_csv(mf_p, sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig16", []) or [])

    included_communities: list[str] = []
    community_ids = list(dict.fromkeys(mf["community_id"].astype(str).tolist()))
    for cid in community_ids:
        if cid in skipped or _is_excluded(cid, exclude_patterns):
            continue
        included_communities.append(str(cid))

    out_dir = figures_root(outdir) / "fig16_species_abundance_drivers"
    out_dir.mkdir(parents=True, exist_ok=True)

    font_size = float(getattr(cfg.figures, "fig16_font_size", 11.0))
    method_order = _normalize_method_order(getattr(cfg.figures, "fig16_methods", None))
    driver_df, driver_summary_df = _write_species_driver_outputs(
        cfg=cfg,
        outdir=outdir,
        out_dir=out_dir,
        included_communities=included_communities,
        font_size=font_size,
        method_order=method_order,
    )
    driver_df.to_csv(out_dir / "species_abundance_driver_by_community.tsv", sep="\t", index=False)
    driver_summary_df.to_csv(out_dir / "species_abundance_driver_plot_summary.tsv", sep="\t", index=False)
    tip_driver_df, tip_driver_summary_df, tip_copy_mapping_df = _write_tip_driver_outputs(
        cfg=cfg,
        outdir=outdir,
        out_dir=out_dir,
        included_communities=included_communities,
        font_size=font_size,
        method_order=method_order,
    )
    tip_driver_df.to_csv(out_dir / "v4_copy_abundance_driver_by_community.tsv", sep="\t", index=False)
    tip_driver_summary_df.to_csv(out_dir / "v4_copy_abundance_driver_plot_summary.tsv", sep="\t", index=False)
    tip_copy_mapping_df.to_csv(out_dir / "v4_copy_number_reset_mapping.tsv", sep="\t", index=False)

    n_plots = 0
    if not driver_summary_df.empty and "n_species_plotted" in driver_summary_df.columns:
        n_plots = int(pd.to_numeric(driver_summary_df["n_species_plotted"], errors="coerce").fillna(0).gt(0).sum())
    n_tip_plots = 0
    if not tip_driver_summary_df.empty and "n_v4_copies_plotted" in tip_driver_summary_df.columns:
        n_tip_plots = int(pd.to_numeric(tip_driver_summary_df["n_v4_copies_plotted"], errors="coerce").fillna(0).gt(0).sum())
    log.info(
        f"Wrote Figure 16 to {out_dir} "
        f"(communities={len(included_communities)}, species_plots={n_plots}, v4_copy_plots={n_tip_plots})"
    )
    return 0
