from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.normalize import tss
from ..utils.paths import expected_dir, method_dir, qc_dir


METHODS = ["ASV", "zOTU", "TIC", "TAC", "OTU"]
CLUSTERING = {"TIC", "TAC", "OTU"}
METHOD_KEY_ALIASES = {
    "ASV": ["asv", "dada2"],
    "zOTU": ["zotu", "unoise3"],
    "OTU": ["otu", "uparse"],
    "TIC": ["tic"],
    "TAC": ["tac"],
}


def _resolve_target_level(cfg, method: str) -> str:
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


def run(cfg) -> int:
    """
    Project each method's standardized table onto expected references using mapping.tsv.

    Inputs per community:
      - expected tip table:    expected/expected_tip_abundance_tss1000.tsv
      - expected species table expected/expected_species_abundance_tss1000.tsv
      - method standardized:   methods/<method>/table_std_tss1000.tsv
      - method mapping:        methods/<method>/mapping.tsv

    Outputs per method:
      - ASV/zOTU: mapped_tip_table_tss1000.tsv (index=expected tips)
      - TIC/TAC/OTU: mapped_species_table_tss1000.tsv (index=species)
      - mapping_qc.tsv (optional): per-sample mapped/unmapped mass and mapping_rate

    Notes:
      - Input method tables are already TSS-normalized to 1000 in prepare_methods.
        We keep that, and if cfg.normalization.renormalize_after_mapping=True,
        re-normalize the mapped table back to 1000 after discarding unmapped features.
      - Unmapped inferred features are not included in mapped tables.
        Their mass is recorded in mapping_qc.tsv if enabled.
      - If mapping.tsv contains weighted rows (target_species/target_tip + weight),
        feature abundance is distributed by those weights.
      - For tip-level outputs with species-weighted rows, abundance is split
        across expected tips within each species.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    manifest = pd.read_csv(outdir / "manifest.tsv", sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    for _, r in manifest.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue

        ed = expected_dir(outdir, cid)
        exp_tip_p = ed / "expected_tip_abundance_tss1000.tsv"
        exp_sp_p = ed / "expected_species_abundance_tss1000.tsv"
        if not exp_tip_p.exists() or not exp_sp_p.exists():
            continue

        exp_tip = pd.read_csv(exp_tip_p, sep="\t", index_col=0)
        exp_sp = pd.read_csv(exp_sp_p, sep="\t", index_col=0)
        exp_tip.index = exp_tip.index.astype(str)
        exp_sp.index = exp_sp.index.astype(str)

        # species -> expected tips and tip -> species
        sp_to_tips: dict[str, list[str]] = {}
        tip_to_species: dict[str, str] = {}
        tipmap_p = ed / "tip_to_species.tsv"
        if tipmap_p.exists():
            try:
                tipmap = pd.read_csv(tipmap_p, sep="\t")
                if {"tip", "species"}.issubset(set(tipmap.columns)):
                    for _, tr in tipmap.iterrows():
                        tip = str(tr["tip"])
                        spx = str(tr["species"])
                        if tip in exp_tip.index:
                            tip_to_species[tip] = spx
                            sp_to_tips.setdefault(spx, []).append(tip)
            except Exception:
                sp_to_tips = {}
                tip_to_species = {}

        for method in METHODS:
            md = method_dir(outdir, cid, method)
            table_p = md / "table_std_tss1000.tsv"
            map_p = md / "mapping.tsv"
            if not table_p.exists() or not map_p.exists():
                continue
            target_level = _resolve_target_level(cfg, method)
            species_target = target_level == "species"

            tab = pd.read_csv(table_p, sep="\t", index_col=0)
            mp = pd.read_csv(map_p, sep="\t")
            tab.index = tab.index.astype(str)

            # Safety: ensure TSS 1000 (in case someone runs this step standalone)
            tab = tss(tab, total=cfg.normalization.tss_total, axis=0)

            all_feats = tab.index.astype(str).tolist()
            mapped_feats: set[str] = set()
            weighted_species_rows: Optional[pd.DataFrame] = None
            weighted_tip_rows: Optional[pd.DataFrame] = None

            has_weight = {"inferred_feature", "weight"}.issubset(set(mp.columns))
            if has_weight:
                if "target_species" in mp.columns:
                    mps = mp[["inferred_feature", "target_species", "weight"]].copy()
                    mps["inferred_feature"] = mps["inferred_feature"].astype(str)
                    mps["target_species"] = mps["target_species"].astype(str)
                    mps["weight"] = pd.to_numeric(mps["weight"], errors="coerce").fillna(0.0).clip(lower=0.0)
                    mps = mps[mps["weight"] > 0].copy()
                    mps = mps[mps["target_species"].isin(set(exp_sp.index.astype(str)))].copy()
                    if not mps.empty:
                        mps["weight"] = mps.groupby("inferred_feature", sort=False)["weight"].transform(
                            lambda s: s / float(s.sum()) if float(s.sum()) > 0 else 0.0
                        )
                        weighted_species_rows = mps[mps["weight"] > 0].copy()
                        mapped_feats |= set(weighted_species_rows["inferred_feature"].astype(str).tolist())

                if "target_tip" in mp.columns:
                    mpt = mp[["inferred_feature", "target_tip", "weight"]].copy()
                    mpt["inferred_feature"] = mpt["inferred_feature"].astype(str)
                    mpt["target_tip"] = mpt["target_tip"].astype(str)
                    mpt["weight"] = pd.to_numeric(mpt["weight"], errors="coerce").fillna(0.0).clip(lower=0.0)
                    mpt = mpt[mpt["weight"] > 0].copy()
                    mpt = mpt[mpt["target_tip"].isin(set(exp_tip.index.astype(str)))].copy()
                    if not mpt.empty:
                        mpt["weight"] = mpt.groupby("inferred_feature", sort=False)["weight"].transform(
                            lambda s: s / float(s.sum()) if float(s.sum()) > 0 else 0.0
                        )
                        weighted_tip_rows = mpt[mpt["weight"] > 0].copy()
                        mapped_feats |= set(weighted_tip_rows["inferred_feature"].astype(str).tolist())

                # Fallback: derive species weights from tip weights if species rows were not provided.
                if (weighted_species_rows is None) and isinstance(weighted_tip_rows, pd.DataFrame):
                    if not weighted_tip_rows.empty:
                        mps = weighted_tip_rows.copy()
                        mps["target_species"] = mps["target_tip"].map(lambda t: str(tip_to_species.get(str(t), "")))
                        mps = mps[mps["target_species"].astype(str).str.strip() != ""].copy()
                        if not mps.empty:
                            mps = (
                                mps.groupby(["inferred_feature", "target_species"], sort=False, as_index=False)["weight"]
                                .sum()
                            )
                            mps["weight"] = mps.groupby("inferred_feature", sort=False)["weight"].transform(
                                lambda s: s / float(s.sum()) if float(s.sum()) > 0 else 0.0
                            )
                            weighted_species_rows = mps[mps["weight"] > 0].copy()
                            mapped_feats |= set(weighted_species_rows["inferred_feature"].astype(str).tolist())

            if not mapped_feats:
                mapped_feats = set(mp["inferred_feature"].astype(str).tolist())

            unmapped_feats = [f for f in all_feats if f not in mapped_feats]

            unmapped_mass = tab.loc[unmapped_feats].sum(axis=0) if unmapped_feats else pd.Series(
                0.0, index=tab.columns
            )
            mapped_mass = tab.drop(index=unmapped_feats, errors="ignore").sum(axis=0)

            if species_target:
                out = pd.DataFrame(0.0, index=exp_sp.index.astype(str), columns=tab.columns)

                if isinstance(weighted_species_rows, pd.DataFrame):
                    for feat, g in weighted_species_rows.groupby("inferred_feature", sort=False):
                        if feat not in tab.index:
                            continue
                        vals = tab.loc[feat, :].to_numpy(dtype=float)
                        for _, rec in g.iterrows():
                            spx = str(rec["target_species"])
                            w = float(rec["weight"])
                            if w <= 0:
                                continue
                            if spx in out.index:
                                out.loc[spx, :] += (vals * w)
                else:
                    fmap = dict(
                        zip(mp["inferred_feature"].astype(str), mp["target_species"].astype(str))
                    )

                    for feat in all_feats:
                        tgt = fmap.get(feat)
                        if tgt is None:
                            continue
                        if tgt in out.index:
                            out.loc[tgt, :] += tab.loc[feat, :].values

                if cfg.normalization.renormalize_after_mapping:
                    out = tss(out, total=cfg.normalization.tss_total, axis=0)

                out.to_csv(md / "mapped_species_table_tss1000.tsv", sep="\t")

            else:
                out = pd.DataFrame(0.0, index=exp_tip.index.astype(str), columns=tab.columns)

                if isinstance(weighted_tip_rows, pd.DataFrame):
                    for feat, g in weighted_tip_rows.groupby("inferred_feature", sort=False):
                        if feat not in tab.index:
                            continue
                        vals = tab.loc[feat, :].to_numpy(dtype=float)
                        for _, rec in g.iterrows():
                            tgt = str(rec["target_tip"])
                            w = float(rec["weight"])
                            if w <= 0:
                                continue
                            if tgt in out.index:
                                out.loc[tgt, :] += (vals * w)
                elif isinstance(weighted_species_rows, pd.DataFrame):
                    for feat, g in weighted_species_rows.groupby("inferred_feature", sort=False):
                        if feat not in tab.index:
                            continue
                        vals = tab.loc[feat, :].to_numpy(dtype=float)
                        for _, rec in g.iterrows():
                            spx = str(rec["target_species"])
                            w = float(rec["weight"])
                            if w <= 0:
                                continue
                            tips = sp_to_tips.get(spx, [])
                            if not tips:
                                continue
                            share = vals * w / float(len(tips))
                            for tip in tips:
                                out.loc[tip, :] += share
                else:
                    fmap = dict(zip(mp["inferred_feature"].astype(str), mp["target_tip"].astype(str)))

                    for feat in all_feats:
                        tgt = fmap.get(feat)
                        if tgt is None:
                            continue
                        if tgt in out.index:
                            out.loc[tgt, :] += tab.loc[feat, :].values

                if cfg.normalization.renormalize_after_mapping:
                    out = tss(out, total=cfg.normalization.tss_total, axis=0)

                out.to_csv(md / "mapped_tip_table_tss1000.tsv", sep="\t")

            if cfg.mapping.store_unmapped_mass:
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

        log.info(f"[{cid}] projected abundances")

    return 0
