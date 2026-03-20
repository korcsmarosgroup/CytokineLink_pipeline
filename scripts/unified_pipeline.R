#!/usr/bin/env Rscript
# Unified NicheNet Pipeline
# Combines model generation, validation, and network combination
# Author: Marton Olbei

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(nichenetr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript unified_pipeline.R <seurat_rds_file> <idents_slot> [output_prefix]\n")
  cat("\nArguments:\n")
  cat("  seurat_rds_file : Path to input Seurat object (.rds file)\n")
  cat("  idents_slot     : Name of the Idents slot containing celltypes\n")
  cat("  output_prefix   : Optional output file prefix (default: 'cytokine_network')\n")
  quit(status = 1)
}

seurat_file <- args[1]
idents_slot <- args[2]
output_prefix <- ifelse(length(args) >= 3, args[3], "cytokine_network")

cat("=======================================================\n")
cat("UNIFIED NICHENET PIPELINE\n")
cat("=======================================================\n")
cat("Input Seurat file:", seurat_file, "\n")
cat("Idents slot:", idents_slot, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("=======================================================\n\n")

# Source the iterateCells function
source('iteratecells.R')

# ============================================================
# STEP 1: LOAD NICHENET RESOURCES
# ============================================================
cat("[STEP 1] Loading NicheNet resources...\n")

ligand_target_matrix <- readRDS("../rds/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("../rds/weighted_networks_nsga2r_final.rds")
ligand_tf_matrix <- readRDS('../rds/ligand_tf_matrix_nsga2r_final.rds')
networks <- readRDS("../rds/OP_nichenet_networks.RDS")

lr_network <- networks$lr_network %>%
  dplyr::filter(from %in% colnames(ligand_target_matrix))
lr_network <- dplyr::rename(lr_network, ligand = from, receptor = to)

sig_network <- networks$signaling_network
gr_network <- networks$gr_network

geneset_oi <- read_tsv("../cytokine_list.tsv") %>%
  pull(Name) %>%
  unique() %>%
  .[. %in% rownames(ligand_target_matrix)] %>%
  .[. %in% colnames(ligand_target_matrix)]

cat("NicheNet resources loaded successfully.\n\n")

# ============================================================
# STEP 2: GENERATE PRELIMINARY NETWORK
# ============================================================
cat("[STEP 2] Generating preliminary cytokine network...\n")

# Load Seurat object
seurat_obj <- readRDS(seurat_file)
Idents(seurat_obj) <- idents_slot

# Normalize data
seurat_obj <- NormalizeData(seurat_obj)

# Extract celltypes
celltypes <- seurat_obj %>%
  pull(!!sym(idents_slot)) %>%
  as.character() %>%
  unique()

cat("Found", length(celltypes), "cell types\n")
cat("Cell types:", paste(celltypes, collapse = ", "), "\n\n")

# Run pairwise cell-cell interactions
cat("Running pairwise cell-cell interaction analysis...\n")
results <- list()
total_pairs <- length(celltypes) * length(celltypes)
current_pair <- 0

for (i in celltypes) {
  for (j in celltypes) {
    current_pair <- current_pair + 1
    cat(sprintf("  Processing pair %d/%d: %s -> %s\n", current_pair, total_pairs, i, j))
    hits <- iterateCells(i, j, seurat_obj)
    results[[paste0(i, ":", j)]] <- hits
  }
}

# Combine results and save
preliminary_network <- keep(results, is.data.frame) %>% bind_rows()
preliminary_file <- paste0(output_prefix, "_preliminary.tsv")
write_tsv(preliminary_network, preliminary_file)

cat("Preliminary network saved to:", preliminary_file, "\n")
cat("Total interactions found:", nrow(preliminary_network), "\n\n")

# ============================================================
# STEP 3: VALIDATE NETWORK WITH RANDOM FOREST
# ============================================================
cat("[STEP 3] Validating network with random forest...\n")

# Define validation function
rf_validation <- function(seurat_obj, tsv_file, output_file) {
  Idents(seurat_obj) <- idents_slot

  ligand_df <- read_tsv(tsv_file) %>%
    dplyr::select(target_cell, ligand) %>%
    distinct()

  cat("Validating", nrow(ligand_df), "unique ligand-receiver pairs...\n")

  tryCatch(
    expr = {},
    error = function(e) {
      cat("An error has occurred:\n")
      cat(e$message, "\n")
    },
    warning = function(w) {
      cat("Warning:\n")
      cat(w$message, "\n")
    },
    finally = {
      validationRes <- list()

      for (idx in 1:nrow(ligand_df)) {
        row <- ligand_df[idx, ]
        upstreamCytokine <- row$ligand %>% as.character()
        receiver <- row$target_cell

        cat(sprintf("  Validating %d/%d: %s -> %s\n", idx, nrow(ligand_df), upstreamCytokine, receiver))

        expressed_genes_receiver <- get_expressed_genes(receiver, seurat_obj, pct = 0.10)
        background_expressed_genes <- expressed_genes_receiver %>%
          .[. %in% rownames(ligand_target_matrix)]

        ligands_of_interest <- c(upstreamCytokine, upstreamCytokine)

        gene_predictions_top20_list <- seq(5) %>%
          lapply(assess_rf_class_probabilities,
                 folds = 5,
                 geneset = geneset_oi,
                 background_expressed_genes = background_expressed_genes,
                 ligands_oi = ligands_of_interest,
                 ligand_target_matrix = ligand_target_matrix)

        # Get top predicted genes
        top_predicted_genes <- seq(length(gene_predictions_top20_list)) %>%
          lapply(get_top_predicted_genes, gene_predictions_top20_list) %>%
          reduce(full_join, by = c("gene", "true_target"))

        top_predicted_cytokine_genes <- top_predicted_genes %>%
          filter(true_target) %>%
          dplyr::filter(gene %in% geneset_oi) %>%
          drop_na()

        top_predicted_cytokine_genes$upstrm <- upstreamCytokine
        top_predicted_cytokine_genes$target_cell <- receiver
        top_predicted_cytokine_genes <- top_predicted_cytokine_genes %>%
          dplyr::select('upstrm', 'gene', 'target_cell')

        validationRes[[paste0(upstreamCytokine, ':', receiver)]] <- as_tibble(top_predicted_cytokine_genes)
      }

      saveRDS(validationRes, output_file)
      cat("Validation results saved to:", output_file, "\n")
    }
  )
}

# Run validation
validated_file <- paste0(output_prefix, "_validated.rds")
rf_validation(seurat_obj, preliminary_file, validated_file)

cat("Validation complete.\n\n")

# ============================================================
# STEP 4: COMBINE VALIDATED AND ORIGINAL NETWORKS
# ============================================================
cat("[STEP 4] Combining validated and original networks...\n")

# Read original and validated networks
orig <- read_tsv(preliminary_file) %>%
  rename(source = ligand)

valid <- readRDS(validated_file) %>%
  bind_rows() %>%
  rename(source = upstrm, target = gene)

# Join and process
net <- left_join(valid, orig) %>%
  drop_na()

net_collapse <- net %>%
  tidyr::unite(c('source_cell', 'target_cell'), col = 'CCC', sep = ':') %>%
  group_by(source, target) %>%
  summarise(cells = paste(CCC, collapse = ","), n = n(), .groups = 'drop')

# Save combined network
combined_file <- paste0(output_prefix, "_combined.rds")
saveRDS(net_collapse, combined_file)

cat("Combined network saved to:", combined_file, "\n")
cat("Total validated interactions:", nrow(net_collapse), "\n\n")

# ============================================================
# PIPELINE COMPLETE
# ============================================================
cat("=======================================================\n")
cat("PIPELINE COMPLETE!\n")
cat("=======================================================\n")
cat("Output files:\n")
cat("  1. Preliminary network:", preliminary_file, "\n")
cat("  2. Validated network:", validated_file, "\n")
cat("  3. Combined network:", combined_file, "\n")
cat("=======================================================\n")
