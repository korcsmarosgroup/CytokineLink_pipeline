# Unified NicheNet Pipeline Usage Guide

## Overview

The `unified_pipeline.R` is a more user friendly version of the cytokinelink pipeline introduced in Olbei et al., 2026. It automates the entire process of:

1. Generating preliminary cytokine networks using NicheNet
2. Validating the networks with random forest classification
3. Combining validated and original networks into a final output

## Requirements

The script expects the following folder / file structure:
```md
project_folder/
├─ scripts/
│  ├─ unified_pipeline.R
│  ├─ iteratecells.R
├─ cytokine_list.tsv
├─ rds/
│  ├─ OP_nichenet_networks.RDS
│  ├─ **weighted_networks_nsga2r_final.rds**
│  ├─ **ligand_tf_matrix_nsga2r_final.rds**
│  ├─ **ligand_target_matrix_nsga2r_final.rds**
```
The NicheNet prior models highlighted in **bold** need to be downloaded from: https://zenodo.org/records/7074291

- R packages: `tidyverse`, `Seurat`, `nichenetr`
- NicheNet resource files in `../rds/`:
  - `ligand_target_matrix_nsga2r_final.rds`
  - `weighted_networks_nsga2r_final.rds`
  - `ligand_tf_matrix_nsga2r_final.rds`
  - `OP_nichenet_networks.RDS`
- Cytokine list: `../cytokine_list.tsv`
- Helper script: `iteratecells.R` (in the same directory)

### Notes on usage and default settings:
- cytokine_list.tsv is a list of cytokines to query - feel free to modify this
- OP_nichenet_networks.RDS is a snapshot of OmniPath interactions from 2024 that was used in the analysis, kept here for reproducibility. Feel free to update it.
- The script uses a few default parameters to be aware of, that were used with the original publication. These will be added as modifiable optional arguments in the future: gene expression percentage is set at 10%, ligand activity cutoff is set to pearson correlation >= 0.1.

## Usage

### Basic Command

```bash
Rscript unified_pipeline.R <seurat_rds_file> <idents_slot> [output_prefix]
```

### Arguments

1. **seurat_rds_file** (required): Path to your input Seurat object as an `.rds` file
2. **idents_slot** (required): Name of the metadata column containing cell type annotations
3. **output_prefix** (optional): Prefix for all output files (default: `cytokine_network`)

### Examples

#### Example 1: Using default output prefix
```bash
cd data_processing
Rscript unified_pipeline.R ../rds/my_seurat_object.rds minor_cluster
```

This will create:
- `cytokine_network_preliminary.tsv`
- `cytokine_network_validated.rds`
- `cytokine_network_combined.rds`

#### Example 2: Using custom output prefix
```bash
cd data_processing
Rscript unified_pipeline.R ../rds/uc_inflamed_naive.rds minor_cluster uc_inflamed_naive
```

This will create:
- `uc_inflamed_naive_preliminary.tsv`
- `uc_inflamed_naive_validated.rds`
- `uc_inflamed_naive_combined.rds`

#### Example 3: Processing multiple datasets
```bash
cd data_processing

# Process UC inflamed with naive cells
Rscript unified_pipeline.R ../rds/uc_inflamed_naive.rds minor_cluster uc_inflamed_naive

# Process UC inflamed without naive cells
Rscript unified_pipeline.R ../rds/uc_inflamed_no_naive.rds minor_cluster uc_inflamed_no_naive

# Process healthy samples
Rscript unified_pipeline.R ../rds/scibd_healthy.rds minor_cluster healthy
```

## Output Files

The pipeline generates three output files:

1. **`<prefix>_preliminary.tsv`**: Initial cytokine network with all predicted interactions
   - Columns: `ligand`, `target`, `weight`, `source_cell`, `target_cell`

2. **`<prefix>_validated.rds`**: Random forest validated interactions
   - RDS list with validated ligand-target pairs per cell-cell interaction

3. **`<prefix>_combined.rds`**: Final combined network
   - Columns: `source`, `target`, `cells`, `n`
   - Contains validated interactions with cell-cell communication context

## Workflow Details

### Step 1: Load NicheNet Resources
- Loads ligand-target matrices, weighted networks, and cytokine gene sets
- Prepares ligand-receptor network and regulatory networks

### Step 2: Generate Preliminary Network
- Normalizes input Seurat object
- Extracts cell types from specified Idents slot
- Runs pairwise cell-cell interaction analysis using `iterateCells()`
- Saves preliminary network as TSV

### Step 3: Validate Network
- Uses random forest classification to validate ligand-target relationships
- Performs 5-fold cross-validation with 5 iterations
- Filters for true cytokine targets
- Saves validated results as RDS

### Step 4: Combine Networks
- Joins validated and original networks
- Collapses by ligand-target pairs across cell-cell contexts
- Generates final combined network with communication metadata

## Troubleshooting

### Common Issues

1. **"Error: object 'ligand_target_matrix' not found"**
   - Ensure you're running the script from the `scripts/` directory
   - Check that NicheNet resource files exist in `../rds/`

2. **"Error: Invalid Idents slot"**
   - Verify the Idents slot name matches a column in your Seurat object metadata
   - Use `colnames(seurat_obj@meta.data)` to check available columns

3. **"No significant cytokine ligands found"**
   - This is expected for some cell type pairs with no meaningful interactions
   - The script will continue processing other pairs

## Notes

- The script maintains minimal changes from the original workflows
- Progress is printed to console for each step
- Processing time depends on the number of cell types (N×N pairwise comparisons)
- Memory usage scales with dataset size and number of cell types
