#CytokineLink Pipeline

A command line R pipeline for cytokine network inference and validation from single-cell RNA-seq data, as introduced in Olbei et al., 2026. 

**Given a single-cell dataset, the pipeline:**

- Predicts cytokine-cytokine interactions across all pairwise cell type combinations using NicheNet
- Validates interactions via 5-fold cross-validated random forest classification
- Outputs a final combined network with cell-cell communication context



Example:
```bash
# Process UC inflamed with naive cells
Rscript unified_pipeline.R ../rds/uc_inflamed_naive.rds minor_cluster uc_inflamed_naive
```

The final cytokine networks will be stored as `<prefix>_combined.rds`.

## Requirements

- R packages: `tidyverse`, `Seurat`, `nichenetr`
- NicheNet prior model files (download from [Zenodo](https://zenodo.org/records/7074291)), placed in `rds/`
- A cytokine list (`cytokine_list.tsv`) and OmniPath network snapshot (`OP_nichenet_networks.RDS`) — provided in this repository

See [PIPELINE_USAGE.md](PIPELINE_USAGE.md) for a more detailed documentation including setup, parameter details, expected outputs, and troubleshooting.
