from __future__ import annotations

import numpy as np
import pandas as pd


def tss(df: pd.DataFrame, total: float = 1000.0, axis: int = 0) -> pd.DataFrame:
    """
    Total Sum Scaling (TSS) normalization.

    Parameters
    ----------
    df : pd.DataFrame
        Input table.
    total : float
        Target total per sample.
    axis : int
        axis=0 => normalize each column to sum(total)  (features x samples tables)
        axis=1 => normalize each row to sum(total)

    Returns
    -------
    pd.DataFrame
        Normalized table (float).
    """
    if axis == 0:
        s = df.sum(axis=0).replace(0, np.nan)
        return (df / s).fillna(0.0) * float(total)

    s = df.sum(axis=1).replace(0, np.nan)
    return df.div(s, axis=0).fillna(0.0) * float(total)
