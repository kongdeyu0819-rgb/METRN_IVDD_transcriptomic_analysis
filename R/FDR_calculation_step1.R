#!/usr/bin/env Rscript
# ============================================================
# Step 1: FDR Calculation — METRN DEA across GEO datasets
# Uses: GEOquery, limma, Biobase
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)

# Configuration
out_dir <- "D:/workburry数据/AGENT/医学科研专家/METRN_AIS_科创变更补充/03_分析代码"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
geo_cache <- file.path(out_dir, "geo_cache")
dir.create(geo_cache, showWarnings = FALSE)

# ---- Dataset config ----
datasets <- list(
  list(
    gse_id = "GSE124272",
    tissue = "Peripheral blood mononuclear cells",
    phenotype = "AIS vs healthy controls",
    role = "Discovery",
    string_for_ais = c("AIS", "scoliosis", "patient"),
    string_for_ctrl = c("control", "healthy", "normal")
  ),
  list(
    gse_id = "GSE150408",
    tissue = "Peripheral blood",
    phenotype = "AIS vs non-scoliosis controls",
    role = "Validation",
    string_for_ais = c("AIS", "scoliosis"),
    string_for_ctrl = c("control", "healthy", "non")
  ),
  list(
    gse_id = "GSE23130",
    tissue = "Annulus fibrosus (disc tissue)",
    phenotype = "AIS disc vs control disc",
    role = "Tissue-level exploratory",
    string_for_ais = c("scoliosis", "AIS"),
    string_for_ctrl = c("control", "normal")
  )
)

# ============================================================
# Function: Download and parse GSE
# ============================================================
get_gse_data <- function(gse_id, cache_dir) {
  cat("\n==================================================\n")
  cat("Processing:", gse_id, "\n")
  cat("==================================================\n")
  
  # Check cache first
  cache_file <- file.path(cache_dir, paste0(gse_id, "_expr_annot.RData"))
  if (file.exists(cache_file)) {
    cat("  [Cache] Loading from", cache_file, "\n")
    load(cache_file)
    return(list(expr = expr, pdata = pdata, fdata = fdata, gpl = gpl))
  }
  
  # Download
  cat("  [Download] Fetching", gse_id, "from GEO...\n")
  gse <- getGEO(gse_id, destdir = cache_dir, AnnotGPL = TRUE)
  
  if (length(gse) > 1) {
    gse <- gse[[1]]
  } else {
    gse <- gse[[1]]
  }
  
  # Extract expression matrix
  expr <- exprs(gse)
  pdata <- pData(gse)
  fdata <- fData(gse)
  gpl <- annotation(gse)
  
  cat("  Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")
  cat("  Platform:", gpl, "\n")
  
  # Cache
  save(expr, pdata, fdata, gpl, file = cache_file)
  cat("  [Cache] Saved to", cache_file, "\n")
  
  return(list(expr = expr, pdata = pdata, fdata = fdata, gpl = gpl))
}

# ============================================================
# Function: Label samples (AIS vs Control)
# ============================================================
label_samples <- function(pdata, cfg) {
  n_samples <- nrow(pdata)
  labels <- rep("Unknown", n_samples)
  rownames(labels) <- rownames(pdata)
  
  # Try multiple metadata fields
  # Field 1: title
  if ("title" %in% colnames(pdata)) {
    titles <- tolower(as.character(pdata$title))
    for (s in cfg$string_for_ais) {
      labels[grep(s, titles)] <- "AIS"
    }
    for (s in cfg$string_for_ctrl) {
      labels[grep(s, titles)] <- "Control"
    }
  }
  
  # Field 2: characteristics_ch1 (often contains disease status)
  char_cols <- colnames(pdata)[grepl("characteristics", colnames(pdata))]
  for (cc in char_cols) {
    vals <- tolower(as.character(pdata[[cc]]))
    for (s in cfg$string_for_ais) {
      idx <- grep(s, vals)
      if (length(idx) > 0) labels[idx] <- "AIS"
    }
    for (s in cfg$string_for_ctrl) {
      idx <- grep(s, vals)
      if (length(idx) > 0) labels[idx] <- "Control"
    }
  }
  
  # Field 3: source_name_ch1
  src_cols <- colnames(pdata)[grepl("source", colnames(pdata), ignore.case = TRUE)]
  for (sc in src_cols) {
    vals <- tolower(as.character(pdata[[sc]]))
    for (s in cfg$string_for_ais) {
      idx <- grep(s, vals)
      if (length(idx) > 0) labels[idx] <- "AIS"
    }
    for (s in cfg$string_for_ctrl) {
      idx <- grep(s, vals)
      if (length(idx) > 0) labels[idx] <- "Control"
    }
  }
  
  # Print label distribution
  cat("  Sample labels:\n")
  print(table(labels))
  
  # Print unlabeled samples for debugging
  unknown_idx <- which(labels == "Unknown")
  if (length(unknown_idx) > 0) {
    cat("  Unknown samples (first 5):\n")
    for (i in unknown_idx[1:min(5, length(unknown_idx))]) {
      cat("    ", rownames(pdata)[i], ":", as.character(pdata$title[i]), "\n")
    }
  }
  
  return(labels)
}

# ============================================================
# Function: Run limma DEA
# ============================================================
run_limma_dea <- function(expr, labels, cfg) {
  ais_idx <- which(labels == "AIS")
  ctrl_idx <- which(labels == "Control")
  
  if (length(ais_idx) < 2 || length(ctrl_idx) < 2) {
    cat("  ERROR: Need at least 2 samples per group\n")
    cat("  AIS:", length(ais_idx), "Control:", length(ctrl_idx), "\n")
    return(NULL)
  }
  
  cat("  DEA: AIS n =", length(ais_idx), ", Control n =", length(ctrl_idx), "\n")
  
  # Design matrix
  group <- factor(c(rep("AIS", length(ais_idx)), rep("Control", length(ctrl_idx))))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Use only samples with known labels
  expr_sub <- expr[, c(ais_idx, ctrl_idx)]
  
  # limma pipeline
  fit <- lmFit(expr_sub, design)
  contrast_matrix <- makeContrasts(AIS_vs_Control = AIS - Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE)
  
  # Extract results
  de_results <- topTable(fit2, number = Inf, adjust.method = "BH")
  
  # Add probe ID as column
  de_results$probe_id <- rownames(de_results)
  
  cat("  Total probes tested:", nrow(de_results), "\n")
  cat("  Significant (raw P < 0.05):", sum(de_results$P.Value < 0.05), "\n")
  cat("  Significant (FDR < 0.05):", sum(de_results$adj.P.Val < 0.05), "\n")
  
  return(de_results)
}

# ============================================================
# Function: Annotate probes to genes
# ============================================================
annotate_probes <- function(de_results, fdata, gpl) {
  cat("  [Annotate] Mapping probes to gene symbols...\n")
  
  # Try to get gene symbol from fData
  gene_symbol_col <- NULL
  for (col in colnames(fdata)) {
    if (tolower(col) %in% c("gene symbol", "genesymbol", "gene_symbol", "symbol")) {
      gene_symbol_col <- col
      break
    }
  }
  
  if (is.null(gene_symbol_col)) {
    # Try to find by pattern
    for (col in colnames(fdata)) {
      if (grepl("gene", tolower(col)) && grepl("symbol", tolower(col))) {
        gene_symbol_col <- col
        break
      }
    }
  }
  
  if (!is.null(gene_symbol_col)) {
    # Map probe_id to gene symbol
    probe_to_gene <- fdata[[gene_symbol_col]]
    names(probe_to_gene) <- rownames(fdata)
    
    de_results$gene_symbol <- probe_to_gene[de_results$probe_id]
    cat("  Annotated", sum(!is.na(de_results$gene_symbol) & de_results$gene_symbol != ""), "probes with gene symbols\n")
  } else {
    de_results$gene_symbol <- NA
    cat("  WARNING: Could not find gene symbol column in fData\n")
    cat("  Available columns:", paste(colnames(fdata)[1:min(10, ncol(fdata))], collapse = ", "), "\n")
  }
  
  # Also get ENTREZID if available
  entrez_col <- NULL
  for (col in colnames(fdata)) {
    if (tolower(col) %in% c("entrez_id", "entrezid", "gene id", "geneid")) {
      entrez_col <- col
      break
    }
  }
  if (!is.null(entrez_col)) {
    de_results$entrez_id <- fdata[[entrez_col]][de_results$probe_id]
  }
  
  return(de_results)
}

# ============================================================
# Function: Extract METRN/METRNL results
# ============================================================
extract_metrn_results <- function(de_results) {
  # METRN
  metrn_rows <- de_results[!is.na(de_results$gene_symbol) & 
                            toupper(de_results$gene_symbol) == "METRN", ]
  
  # METRNL (distinct gene!)
  metrnl_rows <- de_results[!is.na(de_results$gene_symbol) & 
                             toupper(de_results$gene_symbol) == "METRNL", ]
  
  return(list(METRN = metrn_rows, METRNL = metrnl_rows))
}

# ============================================================
# Main loop
# ============================================================
all_metrn_summary <- data.frame()

for (cfg in datasets) {
  gse_id <- cfg$gse_id
  
  # Get data
  gse_data <- get_gse_data(gse_id, geo_cache)
  expr <- gse_data$expr
  pdata <- gse_data$pdata
  fdata <- gse_data$fdata
  gpl <- gse_data$gpl
  
  # Label samples
  labels <- label_samples(pdata, cfg)
  
  # Check if we have enough labeled samples
  n_ais <- sum(labels == "AIS")
  n_ctrl <- sum(labels == "Control")
  
  if (n_ais < 2 || n_ctrl < 2) {
    cat("  SKIPPING:", gse_id, "- insufficient labeled samples\n")
    next
  }
  
  # Run DEA
  de_results <- run_limma_dea(expr, labels, cfg)
  if (is.null(de_results)) next
  
  # Annotate
  de_results <- annotate_probes(de_results, fdata, gpl)
  
  # Save full results
  out_csv <- file.path(out_dir, paste0("FDR_results_", gse_id, ".csv"))
  write.csv(de_results, out_csv, row.names = FALSE)
  cat("  Saved full results:", out_csv, "\n")
  
  # Extract METRN results
  metrn_info <- extract_metrn_results(de_results)
  
  cat("\n  *** METRN/METRNL Results for", gse_id, "***\n")
  
  # METRN
  if (nrow(metrn_info$METRN) > 0) {
    for (i in 1:nrow(metrn_info$METRN)) {
      row <- metrn_info$METRN[i, ]
      cat("  METRN (probe:", row$probe_id, ")\n")
      cat("    log2FC =", round(row$logFC, 4), "\n")
      cat("    raw P  =", signif(row$P.Value, 4), "\n")
      cat("    FDR    =", signif(row$adj.P.Val, 4), "\n")
      cat("    Sig (FDR<0.05)?", row$adj.P.Val < 0.05, "\n")
      
      all_metrn_summary <- rbind(all_metrn_summary, data.frame(
        Dataset = gse_id,
        Tissue = cfg$tissue,
        Role = cfg$role,
        Gene = "METRN",
        Probe_ID = as.character(row$probe_id),
        log2FC = round(row$logFC, 4),
        raw_P = signif(row$P.Value, 4),
        FDR_BH = signif(row$adj.P.Val, 4),
        Significant_FDR005 = row$adj.P.Val < 0.05,
        Significant_rawP005 = row$P.Value < 0.05,
        gene_symbol = as.character(row$gene_symbol)
      ))
    }
  } else {
    cat("  METRN: NOT FOUND in this dataset's annotated probes\n")
    cat("  Checking if METRN probes exist but unannotated...\n")
    # Search for METRN in all probe IDs (some platforms don't annotate)
  }
  
  # METRNL
  if (nrow(metrn_info$METRNL) > 0) {
    cat("  METRNL FOUND (distinct from METRN):\n")
    for (i in 1:nrow(metrn_info$METRNL)) {
      row <- metrn_info$METRNL[i, ]
      cat("    METRNL (probe:", row$probe_id, "): log2FC =", 
          round(row$logFC, 4), ", FDR =", signif(row$adj.P.Val, 4), "\n")
    }
  } else {
    cat("  METRNL: NOT FOUND (good - confirms distinct probes)\n")
  }
}

# ============================================================
# Write summary report
# ============================================================
cat("\n\n==================================================\n")
cat("FINAL METRN SUMMARY\n")
cat("==================================================\n\n")

if (nrow(all_metrn_summary) > 0) {
  print(all_metrn_summary, row.names = FALSE)
  
  # Write markdown report
  report_lines <- c(
    "# Step 1: FDR Calculation Results — METRN Differential Expression",
    "",
    "## Summary Table: METRN Expression Across GEO Datasets",
    "",
    "| Dataset | Tissue | Role | Gene | Probe | log2FC | raw P | FDR (BH) | FDR<0.05? | raw P<0.05? |",
    "|---------|--------|------|------|-------|--------|-------|----------|-----------|-------------|"
  )
  
  for (i in 1:nrow(all_metrn_summary)) {
    r <- all_metrn_summary[i, ]
    sig_fdr <- if (r$Significant_FDR005) "✅ Yes" else "❌ No"
    sig_raw <- if (r$Significant_rawP005) "✅ Yes" else "❌ No"
    report_lines <- c(report_lines, sprintf(
      "| %s | %s | %s | %s | %s | %.4f | %.4g | %.4g | %s | %s |",
      r$Dataset, r$Tissue, r$Role, r$Gene, r$Probe_ID,
      r$log2FC, r$raw_P, r$FDR_BH, sig_fdr, sig_raw
    ))
  }
  
  report_lines <- c(report_lines, "", "## Key Findings", "")
  
  for (i in 1:nrow(all_metrn_summary)) {
    r <- all_metrn_summary[i, ]
    report_lines <- c(report_lines, sprintf("### %s (%s)", r$Dataset, r$Tissue), "")
    
    if (r$Significant_FDR005) {
      report_lines <- c(report_lines,
        sprintf("- METRN was **significantly altered** (FDR = %.4g < 0.05)", r$FDR_BH),
        sprintf("- Direction: %s (log2FC = %.4f)", 
                if (r$log2FC < 0) "downregulated" else "upregulated", r$log2FC),
        "- This result survives Benjamini-Hochberg correction", ""
      )
    } else {
      report_lines <- c(report_lines,
        sprintf("- METRN was **NOT significant after FDR correction** (FDR = %.4g)", r$FDR_BH),
        if (r$Significant_rawP005) 
          c(sprintf("- Raw P was significant (P = %.4g) but did not survive multiple testing correction", r$raw_P),
            "- This should be reported as a **nominal difference** requiring cautious interpretation") 
          else
          c(sprintf("- Raw P was not significant (P = %.4g)", r$raw_P)),
        ""
      )
    }
  }
  
  report_lines <- c(report_lines,
    "## Implications for Manuscript", "",
    "See detailed implications in the full report.", ""
  )
  
  report_path <- file.path(out_dir, "FDR_results_summary.md")
  writeLines(report_lines, report_path, useBytes = TRUE)
  cat("\nReport saved to:", report_path, "\n")
  
} else {
  cat("WARNING: METRN not found in any dataset!\n")
  cat("This may indicate probe annotation issues.\n")
  cat("Check the CSV files manually for METRN probes.\n")
}

cat("\nStep 1 complete.\n")
