"""Shared helpers for loading the KRAS lipidomics workbook.

The workbook has a three-row header block above the data:
    row 1: column names (feature, class, ... , sample IDs)
    row 2: short sample names
    row 3: cohort label for each sample (e.g. "WT-Control", "KG12V-UGCGi")
The first 8 columns are lipid metadata; every column after that is one sample.
"""

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_FILE = REPO_ROOT / "data" / "lipidomics" / "lipidomics.xlsx"
RESULTS_DIR = REPO_ROOT / "results" / "lipidomics"
FIGURES_DIR = REPO_ROOT / "figures" / "lipidomics"

N_META_COLS = 8

# Cohort labels used in the workbook
WT_CONTROL = "WT-Control"
KRAS_CONTROL = "KG12V-Control"
WT_UGCGI = "WT-UGCGi"
KRAS_UGCGI = "KG12V-UGCGi"


def load_lipidomics(path=DATA_FILE):
    """Return (meta, abundance) where abundance columns are cohort labels.

    Missing values are filled with 0, matching how the original analysis
    summed lipid classes. Columns are not unique: each cohort label appears
    once per replicate (n=3).
    """
    raw = pd.read_excel(path, header=None)

    # Locate the header row (the one whose first cell is "feature")
    header_idx = raw.index[raw.iloc[:, 0] == "feature"][0]
    header = raw.iloc[header_idx]
    cohorts = raw.iloc[header_idx + 2, N_META_COLS:].astype(str).values
    body = raw.iloc[header_idx + 3:].reset_index(drop=True)

    meta = body.iloc[:, :N_META_COLS].copy()
    meta.columns = header.iloc[:N_META_COLS].astype(str).values

    abundance = body.iloc[:, N_META_COLS:].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    abundance.columns = cohorts
    return meta, abundance


def significance_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def bh_adjust(pvalues):
    """Benjamini-Hochberg FDR adjustment (no extra dependencies)."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out
