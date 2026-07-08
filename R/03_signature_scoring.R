#!/usr/bin/env Rscript
# ============================================================
# Step 3: Signature scoring analysis using GSVA
# ============================================================

library(GEOquery)
library(limma)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(reshape2)

out_dir <- "D:/MedResearch/METRN_IVDD"
geo_cache <- file.path(out_dir, "geo_cache")

# ---- Define gene sets ----
cat("=== Defining gene sets ===\n")

mechanotransduction_genes <- c("PIEZO1", "PIEZO2", "TRPV4", "YAP1", "WWTR1", 
  "CTGF", "CYR61", "ITGA1", "ITGA2", "ITGA3", "ITGA5", "ITGA7", "ITGA10", "ITGA11",
  "ITGB1", "ITGB3", "ITGB5", "ITGB6", "ITGB8",
  "PTK2", "RHOA", "ROCK1", "ROCK2", "LIMK1", "CFL1", 
  "VCL", "TLN1", "FOSB", "CREB3L2", "TMEM102")

ecm_genes <- c("COL2A1", "ACAN", "SOX9", "COL1A1", "FN1", "COMP", 
  "MMP1", "MMP3", "MMP13", "ADAMTS4", "ADAMTS5", "COL10A1", "COL3A1")

autophagy_genes <- c("ATG5", "ATG7", "ATG12", "ULK1", "BECN1", 
  "MAP1LC3B", "SQSTM1", "ATG3", "ATG4B", "ATG4D", "ATG10", "ATG16L1")

senescence_genes <- c("CDKN1A", "CDKN2A", "TP53", "IL6", "CXCL8", 
  "MMP1", "MMP3", "MMP13", "SERPINE1", "RB1", "CCND1")

er_stress_genes <- c("CALR", "HSPA5", "ATF4", "DDIT3", "XBP1", 
  "ERN1", "CREB3L2", "HSP90B1", "PDIA3", "EIF2AK3")

inflammation_genes <- c("IL1B", "IL6", "TNF", "NFKB1", "CXCL8", 
  "CCL2", "TLR2", "TLR4", "IL10", "TNFAIP3")

ferroptosis_genes <- c("GPX4", "SLC7A11", "ACSL4", "FTH1", "TFRC", 
  "NCOA4", "ALOX15", "POR", "FDFT1", "HSPB1")

# Combine into gene set list
gene_sets <- list(
  Mechanotransduction = mechanotransduction_genes,
  ECM_Homeostasis = ecm_genes,
  Autophagy = autophagy_genes,
  Senescence = senescence_genes,
  ER_Stress = er_stress_genes,
  Inflammation = inflammation_genes,
  Ferroptosis = ferroptosis_genes
)

cat("Gene set sizes:\n")
for (nm in names(gene_sets)) {
  cat(nm, ": ", length(gene_sets[nm]), " genes\n")
}

# ---- Load expression data ----
cat("\n=== Loading expression data ===\n")

# GSE124272 (IVDD peripheral blood) - main analysis dataset
gse1 <- getGEO(filename = file.path(geo_cache, "GSE124272_series_matrix.txt.gz"))
expr1 <- exprs(gse1)
pdata1 <- pData(gse1)
fd1 <- fData(gse1)

# GSE23130 (IVDD disc tissue)
gse3 <- getGEO(filename = file.path(geo_cache, "GSE23130_series_matrix.txt.gz"))
expr3 <- exprs(gse3)
pdata3 <- pData(gse3)
fd3 <- fData(gse3)

cat("GSE124272:", nrow(expr1), "probes,", ncol(expr1), "samples\n")
cat("GSE23130:", nrow(expr3), "probes,", ncol(expr3), "samples\n")

# ---- Map probes to gene symbols ----
cat("\n=== Mapping probes to genes ===\n")

map_to_genes <- function(expr_mat, fdata) {
  # Use GENE_SYMBOL column
  gene_col <- "GENE_SYMBOL"
  if (!gene_col %in% colnames(fdata)) {
    for (col in colnames(fdata)) {
      if (grepl("gene|symbol", col, ignore.case = TRUE)) {
        gene_col <- col
        break
      }
    }
  }
  cat("Using gene column:", gene_col, "\n")
  
  # Create gene-to-probe mapping
  fdata$gene_sym <- fdata[[gene_col]]
  # Remove empty gene symbols
  fdata_filtered <- fdata[fdata$gene_sym != "" & !is.na(fdata$gene_sym), ]
  
  # Collapse multiple probes per gene by median
  expr_filtered <- expr_mat[fdata_filtered$ID, ]
  
  # Group by gene symbol
  gene_expr <- t(sapply(unique(fdata_filtered$gene_sym), function(gene) {
    probes <- fdata_filtered$ID[fdata_filtered$gene_sym == gene]
    if (length(probes) == 1) {
      return(expr_filtered[probes, ])
    } else {
      return(apply(expr_filtered[probes, ], 2, median))
    }
  }))
  
  rownames(gene_expr) <- unique(fdata_filtered$gene_sym)
  return(gene_expr)
}

gene_expr1 <- map_to_genes(expr1, fd1)
cat("GSE124272 gene-level matrix:", nrow(gene_expr1), "genes,", ncol(gene_expr1), "samples\n")

gene_expr3 <- map_to_genes(expr3, fd3)
cat("GSE23130 gene-level matrix:", nrow(gene_expr3), "genes,", ncol(gene_expr3), "samples\n")

# ---- Find METRN in gene-level matrix ----
cat("\n=== Finding METRN ===\n")
metrn1 <- gene_expr1["METRN", ]
cat("GSE124272 METRN found:", !is.null(metrn1), "\n")
if (!is.null(metrn1)) cat("METRN values:", round(metrn1, 3), "\n")

metrn3 <- gene_expr3["METRN", ]
cat("GSE23130 METRN found:", !is.null(metrn3), "\n")

# ---- Filter gene sets for available genes ----
cat("\n=== Checking gene set coverage ===\n")

filter_geneset <- function(gene_set, available_genes) {
  found <- intersect(gene_set, available_genes)
  missing <- setdiff(gene_set, available_genes)
  cat("  Found:", length(found), "/", length(gene_set), "genes\n")
  if (length(missing) > 0) cat("  Missing:", missing, "\n")
  return(found)
}

# GSE124272
filtered_sets1 <- lapply(gene_sets, filter_geneset, rownames(gene_expr1))
# Remove empty sets
filtered_sets1 <- filtered_sets1[sapply(filtered_sets1, length) > 0]
cat("Usable gene sets for GSE124272:", names(filtered_sets1), "\n")

# GSE23130
filtered_sets3 <- lapply(gene_sets, filter_geneset, rownames(gene_expr3))
filtered_sets3 <- filtered_sets3[sapply(filtered_sets3, length) > 0]
cat("Usable gene sets for GSE23130:", names(filtered_sets3), "\n")

# ---- GSVA analysis ----
cat("\n=== Running GSVA ===\n")

# Convert to GeneSetCollection
make_gsc <- function(filtered_sets) {
  gsc <- GeneSetCollection()
  for (nm in names(filtered_sets)) {
    gsc <- c(gsc, GeneSet(filtered_sets[[nm]], geneIdType = SymbolIdentifier(), setName = nm))
  }
  return(gsc)
}

run_gsva <- function(gene_expr, filtered_sets, dataset_name) {
  gsc <- make_gsc(filtered_sets)
  cat("Running GSVA for", dataset_name, "...\n")
  
  scores <- gsva(gene_expr, gsc, method = "gsva", kcdf = "Gaussian",
                 rnaseq = FALSE, verbose = FALSE)
  
  cat("GSVA scores dimensions:", dim(scores), "\n")
  return(scores)
}

scores1 <- run_gsva(gene_expr1, filtered_sets1, "GSE124272")
scores3 <- run_gsva(gene_expr3, filtered_sets3, "GSE23130")

# ---- Correlate METRN with signature scores ----
cat("\n=== Correlating METRN with signature scores ===\n")

# GSE124272
cat("\nGSE124272 METRN-Signature correlations:\n")
cor_results1 <- data.frame()
if (!is.null(metrn1)) {
  for (sig_name in rownames(scores1)) {
    cor_val <- cor(metrn1, scores1[sig_name, ], method = "spearman")
    cor_test <- cor.test(metrn1, scores1[sig_name, ], method = "spearman")
    cor_results1 <- rbind(cor_results1, data.frame(
      Dataset = "GSE124272",
      Signature = sig_name,
      Spearman_r = cor_val,
      P_value = cor_test$p.value,
      n_genes_in_set = length(filtered_sets1[[sig_name]])
    ))
    cat(sig_name, ": r =", round(cor_val, 3), "P =", format(cor_test$p.value, digits = 3), "\n")
  }
  
  # FDR correction
  cor_results1$FDR <- p.adjust(cor_results1$P_value, method = "BH")
  cat("\nWith FDR correction:\n")
  print(cor_results1)
}

# GSE23130
cat("\nGSE23130 METRN-Signature correlations:\n")
cor_results3 <- data.frame()
if (!is.null(metrn3)) {
  for (sig_name in rownames(scores3)) {
    cor_val <- cor(metrn3, scores3[sig_name, ], method = "spearman")
    cor_test <- cor.test(metrn3, scores3[sig_name, ], method = "spearman")
    cor_results3 <- rbind(cor_results3, data.frame(
      Dataset = "GSE23130",
      Signature = sig_name,
      Spearman_r = cor_val,
      P_value = cor_test$p.value,
      n_genes_in_set = length(filtered_sets3[[sig_name]])
    ))
    cat(sig_name, ": r =", round(cor_val, 3), "P =", format(cor_test$p.value, digits = 3), "\n")
  }
  
  # FDR correction
  cor_results3$FDR <- p.adjust(cor_results3$P_value, method = "BH")
  cat("\nWith FDR correction:\n")
  print(cor_results3)
}

# ---- Save results ----
cat("\n=== Saving results ===\n")

# Combine correlation results
all_cor <- rbind(cor_results1, cor_results3)
write.csv(all_cor, file.path(out_dir, "signature_correlation_results.csv"), row.names = FALSE)

# Save GSVA scores
write.csv(scores1, file.path(out_dir, "GSVA_scores_GSE124272.csv"))
write.csv(scores3, file.path(out_dir, "GSVA_scores_GSE23130.csv"))

# ---- Visualization ----
cat("\n=== Generating plots ===\n")

# Heatmap of METRN vs signature scores
if (!is.null(metrn1) && nrow(cor_results1) > 0) {
  # Bar plot of correlations
  p <- ggplot(cor_results1, aes(x = reorder(Signature, -Spearman_r), y = Spearman_r, fill = FDR < 0.05)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    scale_fill_manual(values = c("TRUE" = "#F44336", "FALSE" = "#9E9E9E"),
                      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR >= 0.05")) +
    labs(title = "METRN-Signature Correlations (GSE124272)",
         x = "Signature", y = "Spearman r",
         fill = "Significance") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "METRN_signature_correlation_GSE124272.png"), p, width = 8, height = 5)
  ggsave(file.path(out_dir, "METRN_signature_correlation_GSE124272.pdf"), p, width = 8, height = 5)
  cat("GSE124272 correlation bar plot saved\n")
}

if (!is.null(metrn3) && nrow(cor_results3) > 0) {
  p <- ggplot(cor_results3, aes(x = reorder(Signature, -Spearman_r), y = Spearman_r, fill = FDR < 0.05)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    scale_fill_manual(values = c("TRUE" = "#F44336", "FALSE" = "#9E9E9E"),
                      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR >= 0.05")) +
    labs(title = "METRN-Signature Correlations (GSE23130)",
         x = "Signature", y = "Spearman r",
         fill = "Significance") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "METRN_signature_correlation_GSE23130.png"), p, width = 8, height = 5)
  ggsave(file.path(out_dir, "METRN_signature_correlation_GSE23130.pdf"), p, width = 8, height = 5)
  cat("GSE23130 correlation bar plot saved\n")
}

# ---- METRN-PIEZO1 scatter plot (GSE124272) ----
cat("\n=== Generating METRN-PIEZO1 scatter ===\n")
piezo1_vals <- gene_expr1["PIEZO1", ]
if (!is.null(metrn1) && !is.null(piezo1_vals)) {
  labels1_plot <- rep("Control", ncol(gene_expr1))
  for (i in 1:ncol(gene_expr1)) {
    ti <- pdata1$title[i]
    if (grepl("patient|IDD", ti, ignore.case = TRUE)) labels1_plot[i] <- "IDD"
  }
  
  df_scatter <- data.frame(METRN = metrn1, PIEZO1 = piezo1_vals, Group = labels1_plot)
  
  cor_test <- cor.test(df_scatter$METRN, df_scatter$PIEZO1, method = "spearman")
  
  p <- ggplot(df_scatter, aes(x = METRN, y = PIEZO1, color = Group)) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    scale_color_manual(values = c("Control" = "#4CAF50", "IDD" = "#F44336")) +
    labs(title = paste0("METRN-PIEZO1 Correlation (GSE124272)\nSpearman r = ", 
                        round(cor_test$estimate, 3), ", P = ", format(cor_test$p.value, digits = 3)),
         x = "METRN Expression", y = "PIEZO1 Expression") +
    theme_bw() + theme(plot.title = element_text(size = 13, face = "bold"))
  
  ggsave(file.path(out_dir, "METRN_PIEZO1_scatter_GSE124272.png"), p, width = 7, height = 5)
  ggsave(file.path(out_dir, "METRN_PIEZO1_scatter_GSE124272.pdf"), p, width = 7, height = 5)
  cat("METRN-PIEZO1 scatter saved\n")
}

cat("\n=== All Step 3 analyses completed ===\n")
cat("Output saved to:", out_dir, "\n")
