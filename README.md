# CytokineLink Pipeline

A command line R pipeline for cytokine network inference and validation from single-cell RNA-seq data, as introduced in [Olbei et al., 2026.](https://www.science.org/doi/10.1126/scisignal.adt0986) 

**Given a single-cell dataset, the pipeline:**

- Predicts cytokine-cytokine interactions across all pairwise cell type combinations using NicheNet
- Validates interactions via 5-fold cross-validated random forest classification
- Outputs a final combined network with cell-cell communication context

**Basic usage**

```bash
Rscript unified_pipeline.R <seurat_rds_file> <idents_slot> [output_prefix]
```

**Arguments**

1. **seurat_rds_file** (required): Path to your input Seurat object as an `.rds` file
2. **idents_slot** (required): Name of the metadata column containing cell type annotations
3. **output_prefix** (optional): Prefix for all output files (default: `cytokine_network`)

Note: in case of multi-condition Seurat objects (i.e. health and disease) split the conditions into separate objects, and run them seperately.

Example:
```bash
# Process UC inflamed naive dataset as an example
Rscript unified_pipeline.R ../rds/uc_inflamed_naive.rds minor_cluster uc_inflamed_naive
```

The final cytokine networks will be stored as `<prefix>_combined.rds`.

## Requirements

- R packages: `tidyverse`, `Seurat`, `nichenetr`
- NicheNet prior model files (download from [Zenodo](https://zenodo.org/records/7074291)), placed in `rds/`
- A cytokine list (`cytokine_list.tsv`) and OmniPath network snapshot (`OP_nichenet_networks.RDS`), provided in this repository

See [scripts/PIPELINE_USAGE.md](scripts/PIPELINE_USAGE.md) for a more detailed documentation including setup, parameter details, expected outputs, and troubleshooting.
