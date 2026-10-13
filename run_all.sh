#!/usr/bin/env bash
# Regenerate every result table and figure that can be rebuilt from data in this repo.
set -euo pipefail
cd "$(dirname "$0")"

echo "== Lipidomics =="
python scripts/lipidomics/01_summarize_lipid_classes.py
python scripts/lipidomics/02_plot_lipid_figures.py
python scripts/lipidomics/03_volcano_all_lipids.py

echo "== Transcriptomics =="
if [[ -f data/transcriptomics/counts_kras_a146.csv ]]; then
  python scripts/transcriptomics/transcriptomics_analysis.py
else
  echo "Skipping: data/transcriptomics/counts_kras_a146.csv not present (see data/transcriptomics/README.md)"
fi
