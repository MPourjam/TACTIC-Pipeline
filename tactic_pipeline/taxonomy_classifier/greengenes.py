"""
This module holds GreenGenes taxonomy classifier related code.
It contains functions and classes necessary to classify sequences based on
the GreenGenes taxonomy database.
"""
from __future__ import annotations
from pathlib import Path
from tactic_pipeline.taxonomy_classifier.classifier import TaxonomyClassifier
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union



class GreenGenesTaxonomyClassifier(TaxonomyClassifier):
    """
    A concrete implementation of TaxonomyClassifier for the GreenGenes database.
    """
    pass
