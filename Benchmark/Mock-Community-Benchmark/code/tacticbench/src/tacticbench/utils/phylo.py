from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
from Bio import Phylo


@dataclass
class Edge:
    length: float
    descendants: List[int]  # indices of descendant tips in the tip vector


class UniFracTree:
    """
    Minimal tree wrapper for computing:
      - unweighted UniFrac
      - generalized UniFrac (gUniFrac) with alpha

    Notes:
      - Input vectors must be aligned to self.tips order.
      - Branch lengths are taken from the Newick; missing => 0.0.
      - We midpoint-root best-effort for stability (not required for these formulas).
    """

    def __init__(self, newick_path: Path):
        self.tree = Phylo.read(str(newick_path), "newick")
        try:
            self.tree.root_at_midpoint()
        except Exception:
            pass

        self.tips = [t.name for t in self.tree.get_terminals()]
        self.tip_index = {n: i for i, n in enumerate(self.tips)}
        self.edges: List[Edge] = []
        self._build_edges()

    def _build_edges(self) -> None:
        """
        Build a list of edges, each storing branch length and descendant tip indices.
        """
        def postorder(clade):
            if clade.is_terminal():
                d = [self.tip_index[clade.name]] if clade.name in self.tip_index else []
                clade.__dict__["_desc"] = d
                return d
            d: List[int] = []
            for c in clade.clades:
                d.extend(postorder(c))
            clade.__dict__["_desc"] = d
            return d

        postorder(self.tree.root)

        for parent in self.tree.find_clades(order="preorder"):
            for child in parent.clades:
                desc = child.__dict__.get("_desc", [])
                bl = child.branch_length if child.branch_length is not None else 0.0
                self.edges.append(Edge(length=float(bl), descendants=desc))

    def with_synthetic_tips(
        self,
        attached: List[Tuple[str, str, float]],
        unattached: List[Tuple[str, float]] | None = None,
    ) -> "UniFracTree":
        """
        Return a tree copy augmented with synthetic tips.

        attached: (target_tip_name, synthetic_tip_name, branch_length)
          Synthetic tip contributes to all ancestral edges of target_tip and has
          its own pendant edge with branch_length.
        unattached: (synthetic_tip_name, branch_length)
          Synthetic tips added as standalone pendant edges (no ancestral path).
        """
        aug = object.__new__(UniFracTree)
        aug.tree = self.tree
        aug.tips = list(self.tips)
        aug.tip_index = dict(self.tip_index)
        aug.edges = [Edge(length=e.length, descendants=list(e.descendants)) for e in self.edges]

        unattached = unattached or []

        def add_tip(name: str) -> int:
            if name in aug.tip_index:
                return aug.tip_index[name]
            idx = len(aug.tips)
            aug.tips.append(name)
            aug.tip_index[name] = idx
            return idx

        for target_tip, syn_tip, bl in attached:
            tgt_idx = self.tip_index.get(str(target_tip))
            if tgt_idx is None:
                continue
            syn_idx = add_tip(str(syn_tip))
            for e in aug.edges:
                if tgt_idx in e.descendants:
                    e.descendants.append(syn_idx)
            aug.edges.append(Edge(length=float(max(0.0, bl)), descendants=[syn_idx]))

        for syn_tip, bl in unattached:
            syn_idx = add_tip(str(syn_tip))
            aug.edges.append(Edge(length=float(max(0.0, bl)), descendants=[syn_idx]))

        return aug

    @staticmethod
    def _subtree_sum(vec: np.ndarray, desc: List[int]) -> float:
        return float(vec[desc].sum()) if desc else 0.0

    def unifrac_unweighted(self, a: np.ndarray, b: np.ndarray) -> float:
        """
        Unweighted UniFrac distance:
          sum_{edges} L_e * I(present in one but not the other) / sum_{edges} L_e * I(present in either)
        Presence is defined by subtree sum > 0.
        """
        uniq = 0.0
        tot = 0.0
        for e in self.edges:
            pa = self._subtree_sum(a, e.descendants) > 0
            pb = self._subtree_sum(b, e.descendants) > 0
            if pa or pb:
                tot += e.length
                if pa != pb:
                    uniq += e.length
        return 0.0 if tot == 0.0 else float(uniq / tot)

    def gunifrac(self, a: np.ndarray, b: np.ndarray, alpha: float) -> float:
        """
        Generalized UniFrac (gUniFrac) distance.

        Implementation uses the ratio form:
          d = (sum L_e |x_e - y_e|) / (sum L_e (x_e + y_e))
        where x_e and y_e are transformed subtree relative abundances:
          - normalize to relative abundances (sum=1) per sample (if sum>0)
          - x_e = (subtree_sum)^alpha

        alpha=1.0 corresponds to weighted UniFrac-like behavior in this formulation.
        alpha=0.5 corresponds to gUniFrac(0.5).
        """
        ta = float(a.sum())
        tb = float(b.sum())
        if ta <= 0 and tb <= 0:
            return 0.0

        pa = a / ta if ta > 0 else np.zeros_like(a, dtype=float)
        pb = b / tb if tb > 0 else np.zeros_like(b, dtype=float)

        num = 0.0
        den = 0.0
        for e in self.edges:
            xa = self._subtree_sum(pa, e.descendants)
            xb = self._subtree_sum(pb, e.descendants)

            # generalized power transform
            if alpha == 0.0:
                xa = 1.0 if xa > 0 else 0.0
                xb = 1.0 if xb > 0 else 0.0
            else:
                xa = float(xa**alpha)
                xb = float(xb**alpha)

            num += e.length * abs(xa - xb)
            den += e.length * (xa + xb)

        return 0.0 if den == 0.0 else float(num / den)
