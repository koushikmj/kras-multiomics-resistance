#!/usr/bin/env python3
"""Step 1: summarise the raw lipidomics workbook into the tables behind Figures 1-3.

Outputs (results/lipidomics/):
    figure1_control_comparison.csv   lipid-class totals, WT vs KRAS G12V (control)
    figure2_ugcgi_effects.csv        glycosphingolipid totals, control vs UGCGi, per genotype
    figure3_ps_species.csv           individual phosphatidylserine species, WT vs KRAS G12V

Statistics: class total = sum of all species in that class per sample; mean and
sample SD (ddof=1) across n=3 replicates; two-sided Student's t-test (equal variance).
"""

import numpy as np
import pandas as pd
from scipy import stats

from lipid_io import (KRAS_CONTROL, KRAS_UGCGI, RESULTS_DIR, WT_CONTROL, WT_UGCGI,
                      load_lipidomics, significance_stars)

FIG1_CLASSES = ["HexCer", "Cer", "SM", "PS", "CE", "PC", "PE"]
FIG2_CLASSES = ["HexCer", "SM", "Cer"]


def ttest_p(a, b):
    return float(stats.ttest_ind(a, b).pvalue)


def figure1(class_totals):
    rows = []
    for cls in FIG1_CLASSES:
        wt = class_totals.loc[cls, WT_CONTROL].values
        kr = class_totals.loc[cls, KRAS_CONTROL].values
        p = ttest_p(wt, kr)
        rows.append({
            "Lipid_Class": cls,
            "WT_mean": wt.mean(), "WT_std": wt.std(ddof=1),
            "KRAS_mean": kr.mean(), "KRAS_std": kr.std(ddof=1),
            "Fold_Change": kr.mean() / wt.mean(),
            "p_value": p, "significance": significance_stars(p),
        })
    return pd.DataFrame(rows)


def figure2(class_totals):
    rows = []
    for cls in FIG2_CLASSES:
        for genotype, ctrl, treat in [("WT", WT_CONTROL, WT_UGCGI),
                                      ("KRAS G12V", KRAS_CONTROL, KRAS_UGCGI)]:
            c = class_totals.loc[cls, ctrl].values
            t = class_totals.loc[cls, treat].values
            p = ttest_p(c, t)
            rows.append({
                "Lipid_Class": cls, "Genotype": genotype,
                "Control_mean": c.mean(), "Control_std": c.std(ddof=1),
                "UGCGi_mean": t.mean(), "UGCGi_std": t.std(ddof=1),
                "Percent_Change": (t.mean() - c.mean()) / c.mean() * 100,
                "p_value": p, "significance": significance_stars(p),
            })
    return pd.DataFrame(rows)


def figure3(meta, abundance):
    is_ps = (meta["class"] == "PS").values
    ps = abundance[is_ps]
    names = meta.loc[is_ps, "feature"].values
    rows = []
    for name, (_, row) in zip(names, ps.iterrows()):
        wt = row[WT_CONTROL].values.astype(float)
        kr = row[KRAS_CONTROL].values.astype(float)
        if wt.mean() == 0 and kr.mean() == 0:
            continue  # not detected in either group
        # NOTE: species detected in only one group get Log2_FC = 0 (as in the
        # original analysis) because the ratio is undefined. See README limitations.
        log2fc = np.log2(kr.mean() / wt.mean()) if wt.mean() > 0 and kr.mean() > 0 else 0.0
        p = ttest_p(wt, kr)
        rows.append({"PS_Species": name, "WT_mean": wt.mean(), "KRAS_mean": kr.mean(),
                     "Log2_FC": log2fc, "p_value": p, "significance": significance_stars(p)})
    return pd.DataFrame(rows).sort_values("Log2_FC", ascending=False, kind="stable")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    meta, abundance = load_lipidomics()
    class_totals = abundance.groupby(meta["class"].values).sum()

    outputs = {
        "figure1_control_comparison.csv": figure1(class_totals),
        "figure2_ugcgi_effects.csv": figure2(class_totals),
        "figure3_ps_species.csv": figure3(meta, abundance),
    }
    for fname, df in outputs.items():
        df.to_csv(RESULTS_DIR / fname, index=False)
        print(f"Saved {RESULTS_DIR / fname} ({len(df)} rows)")


if __name__ == "__main__":
    main()
