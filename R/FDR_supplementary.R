#!/usr/bin/env Rscript
# ============================================================
# Supplementary: PIEZO1-METRN correlation + GSE56081 fix
# ============================================================

library(GEOquery)
library(limma)

out_dir <- "D:/MedResearch/METRN_IVDD"
geo_cache <- file.path(out_dir, "geo_cache")

# ---- GSE124272 METRN-PIEZO1 correlation ----
cat("=== METRN-PIEZO1 correlation in GSE124272 ===\n")
load(file.path(geo_cache, "GSE124272_gse.RData"))
expr1 <- exprs(gse)
pd1 <- pData(gse)

# Re-read FDR results to get annotation
de1 <- read.csv(file.path(out_dir, "FDR_GSE124272.csv"))

# Find METRN probe(s)
metrn_probes <- de1$probe_id[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "METRN"]
piezo1_probes <- de1$probe_id[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "PIEZO1"]

cat("  METRN probes:", as.character(metrn_probes), "\n")
cat("  PIEZO1 probes:", as.character(piezo1_probes), "\n")

if (length(metrn_probes) > 0 && length(piezo1_probes) > 0) {
  # Use all 16 samples for correlation
  metrn_vals <- as.numeric(expr1[metrn_probes[1], ])
  piezo1_vals <- as.numeric(expr1[piezo1_probes[1], ])
  
  cor_pearson <- cor.test(metrn_vals, piezo1_vals, method = "pearson")
  cor_spearman <- cor.test(metrn_vals, piezo1_vals, method = "spearman")
  
  cat("  Pearson r =", round(cor_pearson$estimate, 4), ", P =", signif(cor_pearson$p.value, 4), "\n")
  cat("  Spearman r =", round(cor_spearman$estimate, 4), ", P =", signif(cor_spearman$p.value, 4), "\n")
  
  # Also in just cases
  case_idx <- 1:8  # First 8 are cases based on earlier labeling
  ctrl_idx <- 9:16
  cat("\n  In cases only:\n")
  cor_s_case <- cor.test(as.numeric(expr1[metrn_probes[1], case_idx]), 
                          as.numeric(expr1[piezo1_probes[1], case_idx]), method = "spearman")
  cat("  Spearman r =", round(cor_s_case$estimate, 4), ", P =", signif(cor_s_case$p.value, 4), "\n")
  
  cat("\n  In controls only:\n")
  cor_s_ctrl <- cor.test(as.numeric(expr1[metrn_probes[1], ctrl_idx]),
                          as.numeric(expr1[piezo1_probes[1], ctrl_idx]), method = "spearman")
  cat("  Spearman r =", round(cor_s_ctrl$estimate, 4), ", P =", signif(cor_s_ctrl$p.value, 4), "\n")
}

# ---- GSE56081 ----
cat("\n\n=== GSE56081 (IVDD nucleus pulposus) ===\n")
load(file.path(geo_cache, "GSE56081_gse.RData"))
expr5 <- exprs(gse)
pd5 <- pData(gse)

# Print all columns
cat("  pData columns:\n")
print(colnames(pd5))

# Print title and source for each sample
cat("\n  Sample details:\n")
for (i in 1:nrow(pd5)) {
  cat("  Sample", i, ": title =", as.character(pd5$title[i]), "\n")
  cat("    source =", as.character(pd5$source_name_ch1[i]), "\n")
  for (cc in colnames(pd5)[grepl("characteristics", colnames(pd5))]) {
    cat("    ", cc, "=", as.character(pd5[[cc]][i]), "\n")
  }
}

# Label based on characteristics
labels5 <- rep("Unknown", nrow(pd5))
for (cc in colnames(pd5)[grepl("characteristics", colnames(pd5))]) {
  vals <- tolower(as.character(pd5[[cc]]))
  idx_degenerat <- grep("degenerat|patient|ivdd|herniated|prolapse", vals)
  idx_normal <- grep("normal|control|healthy|non-degenerat|non-degenerate", vals)
  if (length(idx_degenerat) > 0) labels5[idx_degenerat] <- "Case"
  if (length(idx_normal) > 0) labels5[idx_normal] <- "Control"
}
cat("\n  Labels:\n")
print(table(labels5))

# If still all Unknown, try title
if (all(labels5 == "Unknown")) {
  for (i in 1:nrow(pd5)) {
    title <- tolower(as.character(pd5$title[i]))
    if (grepl("degenerat|herniated|prolapse|patient", title)) labels5[i] <- "Case"
    if (grepl("normal|control|healthy|non-degenerat", title)) labels5[i] <- "Control"
  }
  cat("  Labels (from title):\n")
  print(table(labels5))
}

case_idx5 <- which(labels5 == "Case")
ctrl_idx5 <- which(labels5 == "Control")
cat("  Case:", length(case_idx5), "Control:", length(ctrl_idx5), "\n")

if (length(case_idx5) >= 2 && length(ctrl_idx5) >= 2) {
  group5 <- factor(c(rep("Case", length(case_idx5)), rep("Control", length(ctrl_idx5))))
  design5 <- model.matrix(~0 + group5)
  colnames(design5) <- levels(group5)
  expr_sub5 <- expr5[, c(case_idx5, ctrl_idx5)]
  fit5 <- lmFit(expr_sub5, design5)
  contrast5 <- makeContrasts(Case_vs_Control = Case - Control, levels = design5)
  fit5b <- contrasts.fit(fit5, contrast5)
  fit5b <- eBayes(fit5b, trend = TRUE)
  de5 <- topTable(fit5b, number = Inf, adjust.method = "BH")
  de5$probe_id <- rownames(de5)
  
  fdata5 <- fData(gse)
  sym_col5 <- colnames(fdata5)[grepl("gene.symbol|Gene Symbol", colnames(fdata5), ignore.case = TRUE)]
  if (length(sym_col5) > 0) {
    sym_map5 <- fdata5[[sym_col5[1]]]
    names(sym_map5) <- rownames(fdata5)
    de5$gene_symbol <- sym_map5[de5$probe_id]
  }
  
  write.csv(de5, file.path(out_dir, "FDR_GSE56081.csv"), row.names = FALSE)
  
  metrn5 <- de5[!is.na(de5$gene_symbol) & toupper(de5$gene_symbol) == "METRN", ]
  cat("\n  *** METRN results for GSE56081 ***\n")
  if (nrow(metrn5) > 0) {
    for (i in 1:nrow(metrn5)) {
      cat("  METRN: log2FC =", round(metrn5$logFC[i], 4),
          ", raw P =", signif(metrn5$P.Value[i], 4),
          ", FDR =", signif(metrn5$adj.P.Val[i], 4), "\n")
    }
  } else {
    cat("  METRN: NOT FOUND\n")
  }
  
  key_genes <- c("PIEZO1", "ATG12", "ULK1", "COL2A1", "ACAN", "SOX9", "MMP13", "CALR")
  cat("\n  Key genes:\n")
  for (g in key_genes) {
    rows <- de5[!is.na(de5$gene_symbol) & toupper(de5$gene_symbol) == g, ]
    if (nrow(rows) > 0) {
      cat("  ", g, ": log2FC =", round(rows$logFC[1], 4),
          ", raw P =", signif(rows$P.Value[1], 4),
          ", FDR =", signif(rows$adj.P.Val[1], 4), "\n")
    }
  }
} else {
  cat("  Insufficient labeled samples. Manual labeling needed.\n")
}

cat("\n\nDone.\n")
