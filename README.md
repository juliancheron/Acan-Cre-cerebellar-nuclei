# Acan-Cre-cerebellar-nuclei
## Data availability

The BARseq3 dataset associated with this repository is publicly available on Zenodo:

DOI: 10.5281/zenodo.20214221

The dataset contains the processed BARseq3 AnnData (.h5ad) object used for projection mapping and cell-type classification analyses.

## Repository contents

This repository contains MATLAB scripts used for BARseq3 cell-type classification and downstream analysis of cerebellar nuclei neurons from wild-type and Acan-P2A-Cre mice.

### Main scripts

- `classify_BARseq3_cells.m`  
  Classification pipeline used to categorize cells into A cells, B cells, inhibitory neurons, non-neuronal cells, and unassigned cells based on marker gene expression thresholds.
