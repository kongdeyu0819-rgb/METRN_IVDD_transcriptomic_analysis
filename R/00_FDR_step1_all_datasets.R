#!/usr/bin/env Rscript
# ============================================================
# STEP 1 + STEP 2: FDR Calculation & GSE311180 Analysis
# METRN in IVDD and AIS — Integrated Transcriptomic Analysis
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)

# ---- Configuration ----
out_dir <- "D:/MedResearch/METRN_IVDD"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
geo_cache <- file.path(out_dir, "geo_cache")
dir.create(geo_cache, recursive = TRUE, showWarnings = FALSE)

# ---- Dataset config ----
datasets <- list(
  list(
    gse_id = "GSE124272",
    disease = "IVDD",
    tissue = "Peripheral blood (whole blood)",
    strings_ais = c("disc", "prolapse", "IDD", "patient"),
    strings_ctrl = c("healthy", "volunteer", "control", "normal")
  ),
  list(
    gse_id = "GSE150408",
    disease = "IVDD",
    tissue = "Peripheral blood",
    strings_ais = c("sciatica", "IVDD", "patient"),
    strings_ctrl = c("healthy", "volunteer", "control")
  ),
  list(
    gse_id = "GSE23130",
    disease = "IVDD",
    tissue = "Annulus fibrosus (disc tissue)",
    strings_ais = c("degenerative", "degeneration", "patient"),
    strings_ctrl = c("control", "normal", "non-degenerative")
  ),
  list(
    gse_id = "GSE311180",
    disease = "AIS",
    tissue = "Paravertebral muscle + Spinal cartilage (RNA-seq)",
    strings_ais = c("scoliosis", "AIS", "patient"),
    strings_ctrl = c("control", "unaffected", "normal")
  )
)

# ============================================================
# Function: Download GSE (with cache)
# ============================================================
get_gse <- function(gse_id, cache_dir) {
  cat("\n====== Downloading", gse_id, "======\n")
  cache_file <- file.path(cache_dir, paste0(gse_id, "_gse.RData"))
  if (file.exists(cache_file)) {
    cat("  [Cache] Loading", cache_file, "\n")
    load(cache_file)
    return(gse)
  }
  gse <- getGEO(gse_id, destdir = cache_dir, AnnotGPL = TRUE)
  if (length(gse) > 1) gse <- gse[[1]] else gse <- gse[[1]]
  save(gse, file = cache_file)
  cat("  [Download] Saved to", cache_file, "\n")
  return(gse)
}

# ============================================================
# Function: Label samples
# ============================================================
label_samples <- function(pdata, cfg) {
  n <- nrow(pdata)
  labels <- rep("Unknown", n)
  names(labels) <- rownames(pdata)
  
  # Search in all characterstic columns
  char_cols <- colnames(pdata)[grepl("characteristics|title|source|description", 
                                     colnames(pdata), ignore.case = TRUE)]
  for (cc in char_cols) {
    if (cc %in% colnames(pdata)) {
      vals <- tolower(as.character(pdata[[cc]]))
      for (s in cfg$strings_ais) {
        idx <- grep(s, vals, ignore.case = TRUE)
        if (length(idx) > 0) labels[idx] <- "Case"
      }
      for (s in cfg$strings_ctrl) {
        idx <- grep(s, vals, ignore.case = TRUE)
        if (length(idx) > 0) labels[idx] <- "Control"
      }
    }
  }
  
  cat("  Sample labels for", cfg$gse_id, ":\n")
  print(table(labels))
  return(labels)
}

# ============================================================
# Function: Run limma (microarray)
# ============================================================
run_limma <- function(expr, labels) {
  ais_idx <- which(labels == "Case")
  ctrl_idx <- which(labels == "Control")
  if (length(ais_idx) < 2 || length(ctrl_idx) < 2) return(NULL)
  
  group <- factor(c(rep("Case", length(ais_idx)), rep("Control", length(ctrl_idx))))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  expr_sub <- expr[, c(ais_idx, ctrl_idx)]
  fit <- lmFit(expr_sub, design)
  contrast_matrix <- makeContrasts(Case_vs_Control = Case - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE)
  
  top <- topTable(fit2, number = Inf, adjust.method = "BH")
  top$probe_id <- rownames(top)
  return(top)
}

# ============================================================
# Function: Run DESeq2 (RNA-seq: GSE311180)
# ============================================================
run_deseq2 <- function(counts, labels) {
  if (!require(DESeq2, quietly = TRUE)) {
    cat("  DESeq2 not available, skipping GSE311180\n")
    return(NULL)
  }
  ais_idx <- which(labels == "Case")
  ctrl_idx <- which(labels == "Control")
  if (length(ais_idx) < 2 || length(ctrl_idx) < 2) return(NULL)
  
  coldata <- data.frame(condition = factor(c(rep("Case", length(ais_idx)), 
                                            rep("Control", length(ctrl_idx)))))
  rownames(coldata) <- colnames(counts)[c(ais_idx, ctrl_idx)]
  
  dds <- DESeqDataSetFromMatrix(countData = counts[, c(ais_idx, ctrl_idx)],
                                  colData = coldata,
                                  design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, name = "condition_Case_vs_Control")
  res_df <- as.data.frame(res)
  res_df$probe_id <- rownames(res_df)
  colnames(res_df)[colnames(res_df) == "padj"] <- "adj.P.Val"
  colnames(res_df)[colnames(res_df) == "pvalue"] <- "P.Value"
  colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- "logFC"
  return(res_df)
}

# ============================================================
# Function: Annotate probes
# ============================================================
annotate <- function(de_results, gse) {
  fdata <- fData(gse)
  # Find gene symbol column
  sym_col <- NULL
  for (c in colnames(fdata)) {
    if (tolower(c) %in% c("gene symbol", "genesymbol", "gene_symbol", "symbol"))
      sym_col <- c
  }
  if (is.null(sym_col)) {
    for (c in colnames(fdata)) {
      if (grepl("gene", tolower(c)) && grepl("symbol", tolower(c)))
        sym_col <- c
    }
  }
  if (!is.null(sym_col)) {
    sym_map <- fdata[[sym_col]]
    names(sym_map) <- rownames(fdata)
    de_results$gene_symbol <- sym_map[de_results$probe_id]
  } else {
    de_results$gene_symbol <- NA
  }
  return(de_results)
}

# ============================================================
# Main
# ============================================================
all_metrn <- list()
all_results <- list()

for (cfg in datasets) {
  gse_id <- cfg$gse_id
  cat("\n##################################################\n")
  cat("Processing:", gse_id, "(", cfg$disease, ")\n")
  cat("##################################################\n")
  
  # Download
  gse <- get_gse(gse_id, geo_cache)
  
  # Get expression and pData
  if (gse_id == "GSE311180") {
    # RNA-seq: use featureCounts or RSEM data if available
    # GSE311180 is RNA-seq — GEOquery may return raw counts
    # Try to get processed data
    expr <- tryCatch(exprs(gse), error = function(e) NULL)
    if (is.null(expr)) {
      cat("  RNA-seq data — attempting to parse...\n")
      # For RNA-seq, may need to download count files separately
      cat("  NOTE: GSE311180 is RNA-seq. Use featureCounts/RSEM files from GEO supplementary.\n")
      next
    }
  } else {
    expr <- exprs(gse)
  }
  
  pdata <- pData(gse)
  cat("  Expression:", nrow(expr), "features x", ncol(expr), "samples\n")
  
  # Label
  labels <- label_samples(pdata, cfg)
  
  # Run DEA
  if (gse_id == "GSE311180") {
    # For RNA-seq: use DESeq2 if counts available
    # If expr is log2-normalized counts, use limma-voom
    cat("  Using limma-voom for RNA-seq data...\n")
    de <- run_limma(expr, labels)  # voom would be better
  } else {
    de <- run_limma(expr, labels)
  }
  
  if (is.null(de)) {
    cat("  SKIPPING:", gse_id, "(insufficient samples)\n")
    next
  }
  
  # Annotate
  de <- annotate(de, gse)
  
  # Save
  out_csv <- file.path(out_dir, paste0("FDR_", gse_id, ".csv"))
  write.csv(de, out_csv, row.names = FALSE)
  cat("  Saved:", out_csv, "\n")
  
  # Extract METRN
  metrn <- de[!is.na(de$gene_symbol) & toupper(de$gene_symbol) == "METRN", ]
  metrnl <- de[!is.na(de$gene_symbol) & toupper(de$gene_symbol) == "METRNL", ]
  
  cat("\n  *** METRN results for", gse_id, "***\n")
  if (nrow(metrn) > 0) {
    for (i in 1:nrow(metrn)) {
      r <- metrn[i, ]
      cat("  METRN:", as.character(r$probe_id), "\n")
      cat("    log2FC =", round(r$logFC, 4),
          ", raw P =", signif(r$P.Value, 4),
          ", FDR =", signif(r$adj.P.Val, 4),
          ", FDR<0.05?", r$adj.P.Val < 0.05, "\n")
    }
  } else {
    cat("  METRN: NOT FOUND in annotated results\n")
    cat("  (May be unannotated probe — check CSV manually)\n")
  }
  if (nrow(metrnl) > 0) {
    cat("  METRNL FOUND (distinct from METRN):\n")
    for (i in 1:nrow(metrnl)) {
      r <- metrnl[i, ]
      cat("    METRNL:", as.character(r$probe_id),
          ", log2FC =", round(r$logFC, 4), "\n")
    }
  }
  
  all_metrn[[gse_id]] <- list(METRN = metrn, METRNL = metrnl)
  all_results[[gse_id]] <- de
}

# ============================================================
# Write final summary
# ============================================================
cat("\n\n====== FINAL METRN SUMMARY ======\n\n")
for (gse_id in names(all_metrn)) {
  cat("Dataset:", gse_id, "\n")
  metrn <- all_metrn[[gse_id]]$METRN
  if (nrow(metrn) > 0) {
    for (i in 1:nrow(metrn)) {
      r <- metrn[i, ]
      cat("  METRN: log2FC =", round(r$logFC, 4),
          ", raw P =", signif(r$P.Value, 4),
          ", FDR =", signif(r$adj.P.Val, 4),
          ", Significant (FDR<0.05)?", r$adj.P.Val < 0.05, "\n")
    }
  } else {
    cat("  METRN: NOT FOUND\n")
  }
  cat("\n")
}

cat("Analysis complete. Results saved to:", out_dir, "\n")
