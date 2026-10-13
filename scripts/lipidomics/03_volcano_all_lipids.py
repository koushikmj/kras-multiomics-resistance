#!/usr/bin/env python3
"""Step 3: species-level volcano plot, KRAS G12V vs WT (control cells).

Cleaned-up version of the original lipidomics2.py. Zeros are treated as
"not detected" (NaN); a t-test is run only when both groups have >= 2
detected replicates, otherwise p is set to 1.

Outputs:
    results/lipidomics/volcano_all_lipids_results.csv
    figures/lipidomics/volcano_all_lipids_kras_g12v_vs_wt.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from lipid_io import (FIGURES_DIR, KRAS_CONTROL, RESULTS_DIR, WT_CONTROL,
                      bh_adjust, load_lipidomics)

FC_THRESHOLD = 1.0      # |log2FC| > 1  (2-fold)
PVAL_THRESHOLD = 0.05   # raw p-value, as in the original analysis
N_LABELS = 8


def main():
    sns.set_style("whitegrid")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    meta, abundance = load_lipidomics()
    abundance = abundance.replace(0, np.nan)
    kras = abundance[KRAS_CONTROL]
    wt = abundance[WT_CONTROL]

    res = pd.DataFrame({"feature": meta["feature"].values, "class": meta["class"].values})
    res["KG12V_mean"] = kras.mean(axis=1).values
    res["WT_mean"] = wt.mean(axis=1).values
    res["log2FC"] = np.log2((res["KG12V_mean"] + 1e-10) / (res["WT_mean"] + 1e-10))

    pvals = []
    for i in range(len(res)):
        k = kras.iloc[i].dropna().values
        w = wt.iloc[i].dropna().values
        pvals.append(stats.ttest_ind(k, w).pvalue if len(k) >= 2 and len(w) >= 2 else 1.0)
    res["pvalue"] = np.nan_to_num(np.asarray(pvals, dtype=float), nan=1.0)
    res["-log10_pvalue"] = -np.log10(res["pvalue"])

    res = res.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC", "-log10_pvalue"])
    res["padj_BH"] = bh_adjust(res["pvalue"].values)  # reported, not used for calls

    res["category"] = "Not Significant"
    sig = res["pvalue"] < PVAL_THRESHOLD
    res.loc[sig & (res["log2FC"] > FC_THRESHOLD), "category"] = "Upregulated"
    res.loc[sig & (res["log2FC"] < -FC_THRESHOLD), "category"] = "Downregulated"

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(12, 8))
    colors = {"Upregulated": "red", "Downregulated": "blue", "Not Significant": "gray"}
    for cat in ["Not Significant", "Downregulated", "Upregulated"]:
        sub = res[res["category"] == cat]
        ax.scatter(sub["log2FC"], sub["-log10_pvalue"], c=colors[cat], alpha=0.6, s=30,
                   label=f"{cat} (n={len(sub)})", edgecolors="none")

    ax.axhline(-np.log10(PVAL_THRESHOLD), color="black", ls="--", lw=1, alpha=0.5)
    for x in (FC_THRESHOLD, -FC_THRESHOLD):
        ax.axvline(x, color="black", ls="--", lw=1, alpha=0.5)

    top = res[res["category"] != "Not Significant"].nlargest(N_LABELS, "-log10_pvalue")
    for _, row in top.iterrows():
        ax.annotate(row["feature"], xy=(row["log2FC"], row["-log10_pvalue"]),
                    xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.8,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.3),
                    arrowprops=dict(arrowstyle="->", lw=0.5))

    ax.set_xlabel("log₂(Fold Change)", fontsize=14, fontweight="bold")
    ax.set_ylabel("-log₁₀(p-value)", fontsize=14, fontweight="bold")
    ax.set_title("Volcano Plot: KRAS G12V vs WT Control (all lipid species)",
                 fontsize=16, fontweight="bold", pad=20)
    ax.legend(loc="upper left", frameon=True, fontsize=10)
    ax.set_xlim(res["log2FC"].min() - 0.5, res["log2FC"].max() + 0.5)
    ax.set_ylim(0, res["-log10_pvalue"].max() * 1.1)
    plt.tight_layout()
    out_png = FIGURES_DIR / "volcano_all_lipids_kras_g12v_vs_wt.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    out_csv = RESULTS_DIR / "volcano_all_lipids_results.csv"
    res.to_csv(out_csv, index=False)

    counts = res["category"].value_counts()
    print(f"Lipids tested: {len(res)}")
    for cat in ["Upregulated", "Downregulated", "Not Significant"]:
        print(f"  {cat}: {counts.get(cat, 0)}")
    print(f"  Significant after BH-FDR < 0.05: {(res['padj_BH'] < 0.05).sum()}")
    print(f"Saved {out_png}\nSaved {out_csv}")


if __name__ == "__main__":
    main()
