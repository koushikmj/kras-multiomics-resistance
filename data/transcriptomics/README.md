# Transcriptomics input (not included)

`scripts/transcriptomics/transcriptomics_analysis.py` expects two files here:

| File | Contents |
|---|---|
| `counts_kras_a146.csv` | Gene × sample count matrix. First column = gene ID, then 5 `CTRL_*` and 5 `A146_*` sample columns. |
| `meta_kras_a146.csv` | Sample metadata table (one row per sample). |

These were not in the project archive. The dataset is **synthetic** (gene IDs are
`Gene_1` … `Gene_9999` plus `KRAS`), so it does not describe real biology.

The outputs it produced are already committed in `results/transcriptomics/`
and `figures/transcriptomics/`, so the figures in the report can be inspected
without rerunning the script.
