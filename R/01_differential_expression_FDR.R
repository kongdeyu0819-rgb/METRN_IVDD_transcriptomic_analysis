#!/usr/bin/env Rscript
# ============================================================
# FDR Calculation — Fixed version for all datasets
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)

out_dir <- "D:/MedResearch/METRN_IVDD"
geo_cache <- file.path(out_dir, "geo_cache")

# ============================================================
# GSE124272 & GSE150408 — already processed, load cached
# ============================================================
cat("=== GSE124272 (IVDD peripheral blood) ===\n")
load(file.path(geo_cache, "GSE124272_gse.RData"))
expr1 <- exprs(gse)
pd1 <- pData(gse)

# Manually label: 8 IDD patients + 8 healthy controls
# Title says "lumbar disc prolapse patients vs healthy volunteers"
labels1 <- ifelse(grepl("disc|prolapse|patient|IDD", tolower(pd1$title)), "Case",
            ifelse(grepl("healthy|volunteer|control|normal", tolower(pd1$title)), "Control", "Unknown"))

# Also search characteristics
for (cc in colnames(pd1)[grepl("characteristics", colnames(pd1))]) {
  vals <- tolower(as.character(pd1[[cc]]))
  idx_case <- grep("disc|prolapse|idd|patient", vals)
  idx_ctrl <- grep("healthy|volunteer|control|normal", vals)
  if (length(idx_case) > 0) labels1[idx_case] <- "Case"
  if (length(idx_ctrl) > 0) labels1[idx_ctrl] <- "Control"
}
cat("  Labels:\n")
print(table(labels1))

# Run limma
case_idx <- which(labels1 == "Case")
ctrl_idx <- which(labels1 == "Control")
group <- factor(c(rep("Case", length(case_idx)), rep("Control", length(ctrl_idx))))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
expr_sub <- expr1[, c(case_idx, ctrl_idx)]
fit <- lmFit(expr_sub, design)
contrast_matrix <- makeContrasts(Case_vs_Control = Case - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE)
de1 <- topTable(fit2, number = Inf, adjust.method = "BH")
de1$probe_id <- rownames(de1)

# Annotate
fdata <- fData(gse)
sym_col <- colnames(fdata)[grepl("gene.symbol|Gene Symbol", colnames(fdata), ignore.case = TRUE)]
if (length(sym_col) > 0) {
  sym_map <- fdata[[sym_col[1]]]
  names(sym_map) <- rownames(fdata)
  de1$gene_symbol <- sym_map[de1$probe_id]
} else {
  de1$gene_symbol <- NA
}

write.csv(de1, file.path(out_dir, "FDR_GSE124272.csv"), row.names = FALSE)

# Extract METRN
metrn1 <- de1[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "METRN", ]
metrnl1 <- de1[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "METRNL", ]
cat("\n  *** METRN results for GSE124272 ***\n")
if (nrow(metrn1) > 0) {
  for (i in 1:nrow(metrn1)) {
    r <- metrn1[i, ]
    cat("  METRN: log2FC =", round(r$logFC, 4),
        ", raw P =", signif(r$P.Value, 4),
        ", FDR =", signif(r$adj.P.Val, 4),
        ", FDR<0.05?", r$adj.P.Val < 0.05, "\n")
  }
}
if (nrow(metrnl1) > 0) {
  cat("  METRNL (distinct): log2FC =", round(metrnl1$logFC, 4), "\n")
}

# PIEZO1
piezo1 <- de1[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "PIEZO1", ]
if (nrow(piezo1) > 0) {
  cat("  PIEZO1: log2FC =", round(piezo1$logFC, 4),
      ", raw P =", signif(piezo1$P.Value, 4),
      ", FDR =", signif(piezo1$adj.P.Val, 4), "\n")
}

# Top co-expressed genes
cat("\n  Top 20 DEGs by significance:\n")
top20 <- head(de1[order(de1$P.Value), ], 20)
for (i in 1:nrow(top20)) {
  cat("  ", as.character(top20$gene_symbol[i]),
      ": log2FC =", round(top20$logFC[i], 4),
      ", raw P =", signif(top20$P.Value[i], 4),
      ", FDR =", signif(top20$adj.P.Val[i], 4), "\n")
}

cat("\n  Total DEGs (raw P < 0.05):", sum(de1$P.Value < 0.05), "\n")
cat("  Total DEGs (FDR < 0.05):", sum(de1$adj.P.Val < 0.05), "\n")

# ============================================================
cat("\n\n=== GSE150408 (IVDD peripheral blood) ===\n")
load(file.path(geo_cache, "GSE150408_gse.RData"))
expr2 <- exprs(gse)
pd2 <- pData(gse)

# Labels: sciatica/IVDD vs healthy
labels2 <- rep("Unknown", nrow(pd2))
for (cc in colnames(pd2)[grepl("characteristics|title|source|description", colnames(pd2), ignore.case = TRUE)]) {
  vals <- tolower(as.character(pd2[[cc]]))
  idx_case <- grep("sciatica|ivdd|patient|disc", vals)
  idx_ctrl <- grep("healthy|volunteer|control|normal", vals)
  if (length(idx_case) > 0) labels2[idx_case] <- "Case"
  if (length(idx_ctrl) > 0) labels2[idx_ctrl] <- "Control"
}
cat("  Labels:\n")
print(table(labels2))

case_idx2 <- which(labels2 == "Case")
ctrl_idx2 <- which(labels2 == "Control")
if (length(case_idx2) >= 2 && length(ctrl_idx2) >= 2) {
  group2 <- factor(c(rep("Case", length(case_idx2)), rep("Control", length(ctrl_idx2))))
  design2 <- model.matrix(~0 + group2)
  colnames(design2) <- levels(group2)
  expr_sub2 <- expr2[, c(case_idx2, ctrl_idx2)]
  fit2 <- lmFit(expr_sub2, design2)
  contrast2 <- makeContrasts(Case_vs_Control = Case - Control, levels = design2)
  fit2b <- contrasts.fit(fit2, contrast2)
  fit2b <- eBayes(fit2b, trend = TRUE)
  de2 <- topTable(fit2b, number = Inf, adjust.method = "BH")
  de2$probe_id <- rownames(de2)
  
  fdata2 <- fData(gse)
  sym_col2 <- colnames(fdata2)[grepl("gene.symbol|Gene Symbol", colnames(fdata2), ignore.case = TRUE)]
  if (length(sym_col2) > 0) {
    sym_map2 <- fdata2[[sym_col2[1]]]
    names(sym_map2) <- rownames(fdata2)
    de2$gene_symbol <- sym_map2[de2$probe_id]
  }
  
  write.csv(de2, file.path(out_dir, "FDR_GSE150408.csv"), row.names = FALSE)
  
  metrn2 <- de2[!is.na(de2$gene_symbol) & toupper(de2$gene_symbol) == "METRN", ]
  cat("\n  *** METRN results for GSE150408 ***\n")
  if (nrow(metrn2) > 0) {
    cat("  METRN: log2FC =", round(metrn2$logFC, 4),
        ", raw P =", signif(metrn2$P.Value, 4),
        ", FDR =", signif(metrn2$adj.P.Val, 4), "\n")
  }
}

# ============================================================
cat("\n\n=== GSE23130 (IVDD disc tissue — graded degeneration) ===\n")
load(file.path(geo_cache, "GSE23130_gse.RData"))
expr3 <- exprs(gse)
pd3 <- pData(gse)

# This dataset has Pfirrmann grades I-V, no explicit Case/Control
# Compare: High-grade (IV+V) vs Low-grade (I+II)
# Grade III is intermediate — exclude or group separately
grades <- rep("Unknown", nrow(pd3))
for (cc in colnames(pd3)[grepl("characteristics|tissue_grade", colnames(pd3), ignore.case = TRUE)]) {
  vals <- as.character(pd3[[cc]])
  idx_low <- grep("tissue_grade: I|tissue_grade: II", vals)
  idx_high <- grep("tissue_grade: IV|tissue_grade: V", vals)
  idx_mid <- grep("tissue_grade: III", vals)
  if (length(idx_low) > 0) grades[idx_low] <- "LowGrade"
  if (length(idx_high) > 0) grades[idx_high] <- "HighGrade"
  if (length(idx_mid) > 0) grades[idx_mid] <- "Intermediate"
}
cat("  Pfirrmann grades:\n")
print(table(grades))

# Compare High vs Low (excluding Intermediate)
high_idx <- which(grades == "HighGrade")
low_idx <- which(grades == "LowGrade")
cat("  High-grade (IV+V):", length(high_idx), "samples\n")
cat("  Low-grade (I+II):", length(low_idx), "samples\n")

if (length(high_idx) >= 2 && length(low_idx) >= 2) {
  group3 <- factor(c(rep("HighGrade", length(high_idx)), rep("LowGrade", length(low_idx))))
  design3 <- model.matrix(~0 + group3)
  colnames(design3) <- levels(group3)
  expr_sub3 <- expr3[, c(high_idx, low_idx)]
  fit3 <- lmFit(expr_sub3, design3)
  contrast3 <- makeContrasts(High_vs_Low = HighGrade - LowGrade, levels = design3)
  fit3b <- contrasts.fit(fit3, contrast3)
  fit3b <- eBayes(fit3b, trend = TRUE)
  de3 <- topTable(fit3b, number = Inf, adjust.method = "BH")
  de3$probe_id <- rownames(de3)
  
  fdata3 <- fData(gse)
  sym_col3 <- colnames(fdata3)[grepl("gene.symbol|Gene Symbol", colnames(fdata3), ignore.case = TRUE)]
  if (length(sym_col3) > 0) {
    sym_map3 <- fdata3[[sym_col3[1]]]
    names(sym_map3) <- rownames(fdata3)
    de3$gene_symbol <- sym_map3[de3$probe_id]
  }
  
  write.csv(de3, file.path(out_dir, "FDR_GSE23130.csv"), row.names = FALSE)
  
  metrn3 <- de3[!is.na(de3$gene_symbol) & toupper(de3$gene_symbol) == "METRN", ]
  cat("\n  *** METRN results for GSE23130 (High vs Low grade) ***\n")
  if (nrow(metrn3) > 0) {
    cat("  METRN: log2FC =", round(metrn3$logFC, 4),
        ", raw P =", signif(metrn3$P.Value, 4),
        ", FDR =", signif(metrn3$adj.P.Val, 4), "\n")
  } else {
    cat("  METRN: NOT FOUND in annotation\n")
    # Try searching by probe ID
    metrn_probes <- grep("METRN|meteorin", rownames(expr3), ignore.case = TRUE, value = TRUE)
    cat("  METRN-like probes in expression matrix:", metrn_probes, "\n")
  }
  
  # Also run: Intermediate vs Low
  mid_idx <- which(grades == "Intermediate")
  if (length(mid_idx) >= 2 && length(low_idx) >= 2) {
    group3b <- factor(c(rep("Intermediate", length(mid_idx)), rep("LowGrade", length(low_idx))))
    design3b <- model.matrix(~0 + group3b)
    colnames(design3b) <- levels(group3b)
    expr_sub3b <- expr3[, c(mid_idx, low_idx)]
    fit3b2 <- lmFit(expr_sub3b, design3b)
    contrast3b <- makeContrasts(Mid_vs_Low = Intermediate - LowGrade, levels = design3b)
    fit3b2c <- contrasts.fit(fit3b2, contrast3b)
    fit3b2c <- eBayes(fit3b2c, trend = TRUE)
    de3b <- topTable(fit3b2c, number = Inf, adjust.method = "BH")
    de3b$probe_id <- rownames(de3b)
    if (length(sym_col3) > 0) {
      de3b$gene_symbol <- sym_map3[de3b$probe_id]
    }
    
    metrn3b <- de3b[!is.na(de3b$gene_symbol) & toupper(de3b$gene_symbol) == "METRN", ]
    cat("\n  *** METRN results for GSE23130 (Intermediate vs Low grade) ***\n")
    if (nrow(metrn3b) > 0) {
      cat("  METRN: log2FC =", round(metrn3b$logFC, 4),
          ", raw P =", signif(metrn3b$P.Value, 4),
          ", FDR =", signif(metrn3b$adj.P.Val, 4), "\n")
    }
  }
  
  # Also: High vs Intermediate
  if (length(high_idx) >= 2 && length(mid_idx) >= 2) {
    group3c <- factor(c(rep("HighGrade", length(high_idx)), rep("Intermediate", length(mid_idx))))
    design3c <- model.matrix(~0 + group3c)
    colnames(design3c) <- levels(group3c)
    expr_sub3c <- expr3[, c(high_idx, mid_idx)]
    fit3c <- lmFit(expr_sub3c, design3c)
    contrast3c <- makeContrasts(High_vs_Mid = HighGrade - Intermediate, levels = design3c)
    fit3cc <- contrasts.fit(fit3c, contrast3c)
    fit3cc <- eBayes(fit3cc, trend = TRUE)
    de3c <- topTable(fit3cc, number = Inf, adjust.method = "BH")
    de3c$probe_id <- rownames(de3c)
    if (length(sym_col3) > 0) {
      de3c$gene_symbol <- sym_map3[de3c$probe_id]
    }
    
    metrn3c <- de3c[!is.na(de3c$gene_symbol) & toupper(de3c$gene_symbol) == "METRN", ]
    cat("\n  *** METRN results for GSE23130 (High vs Intermediate grade) ***\n")
    if (nrow(metrn3c) > 0) {
      cat("  METRN: log2FC =", round(metrn3c$logFC, 4),
          ", raw P =", signif(metrn3c$P.Value, 4),
          ", FDR =", signif(metrn3c$adj.P.Val, 4), "\n")
    }
  }
  
  # Key genes for paper
  key_genes <- c("METRN", "METRNL", "PIEZO1", "ATG12", "ULK1", "BECN1", 
                 "COL2A1", "ACAN", "SOX9", "MMP13", "ADAMTS4", "ADAMTS5",
                 "CALR", "LIMK1", "FOSB", "CREB3L2", "TMEM102", "MESP1")
  cat("\n  Key genes in GSE23130 (High vs Low):\n")
  for (g in key_genes) {
    rows <- de3[!is.na(de3$gene_symbol) & toupper(de3$gene_symbol) == g, ]
    if (nrow(rows) > 0) {
      cat("  ", g, ": log2FC =", round(rows$logFC, 4),
          ", raw P =", signif(rows$P.Value, 4),
          ", FDR =", signif(rows$adj.P.Val, 4), "\n")
    } else {
      cat("  ", g, ": NOT FOUND\n")
    }
  }
} else {
  cat("  Insufficient samples for comparison\n")
}

# ============================================================
cat("\n\n=== GSE56081 (IVDD nucleus pulposus — degenerated vs normal) ===\n")
# This was mentioned in the paper but not in the original script
# Try to download
cache5 <- file.path(geo_cache, "GSE56081_gse.RData")
if (!file.exists(cache5)) {
  cat("  Downloading GSE56081...\n")
  gse5 <- tryCatch({
    getGEO("GSE56081", destdir = geo_cache, AnnotGPL = TRUE)
  }, error = function(e) {
    cat("  Failed to download:", e$message, "\n")
    NULL
  })
  if (!is.null(gse5)) {
    if (length(gse5) > 1) gse5 <- gse5[[1]] else gse5 <- gse5[[1]]
    save(gse5, file = cache5)
  }
} else {
  load(cache5)
}

if (exists("gse5") && !is.null(gse5)) {
  expr5 <- exprs(gse5)
  pd5 <- pData(gse5)
  cat("  Expression:", nrow(expr5), "features x", ncol(expr5), "samples\n")
  
  # Label: degenerated vs normal NP
  labels5 <- rep("Unknown", nrow(pd5))
  for (cc in colnames(pd5)[grepl("characteristics|title|source|description", colnames(pd5), ignore.case = TRUE)]) {
    vals <- tolower(as.character(pd5[[cc]]))
    idx_case <- grep("degenerat|patient|disease|ivdd", vals)
    idx_ctrl <- grep("normal|control|healthy|non-degenerat", vals)
    if (length(idx_case) > 0) labels5[idx_case] <- "Case"
    if (length(idx_ctrl) > 0) labels5[idx_ctrl] <- "Control"
  }
  cat("  Labels:\n")
  print(table(labels5))
  
  case_idx5 <- which(labels5 == "Case")
  ctrl_idx5 <- which(labels5 == "Control")
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
    
    fdata5 <- fData(gse5)
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
      cat("  METRN: log2FC =", round(metrn5$logFC, 4),
          ", raw P =", signif(metrn5$P.Value, 4),
          ", FDR =", signif(metrn5$adj.P.Val, 4), "\n")
    }
  }
}

# ============================================================
# METRN-PIEZO1 correlation in GSE124272
# ============================================================
cat("\n\n=== METRN-PIEZO1 correlation in GSE124272 ===\n")
metrn_expr <- expr1[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "METRN", c(case_idx, ctrl_idx)]
piezo1_expr <- expr1[!is.na(de1$gene_symbol) & toupper(de1$gene_symbol) == "PIEZO1", c(case_idx, ctrl_idx)]

if (nrow(metrn_expr) > 0 && nrow(piezo1_expr) > 0) {
  cor_pearson <- cor(as.numeric(metrn_expr[1,]), as.numeric(piezo1_expr[1,]), method = "pearson")
  cor_spearman <- cor(as.numeric(metrn_expr[1,]), as.numeric(piezo1_expr[1,]), method = "spearman")
  test_p <- cor.test(as.numeric(metrn_expr[1,]), as.numeric(piezo1_expr[1,]), method = "spearman")
  cat("  Pearson r =", round(cor_pearson, 4), "\n")
  cat("  Spearman r =", round(cor_spearman, 4), ", P =", signif(test_p$p.value, 4), "\n")
}

# ============================================================
# FINAL SUMMARY
# ============================================================
cat("\n\n========================================\n")
cat("FINAL METRN FDR SUMMARY\n")
cat("========================================\n\n")

cat("GSE124272 (IVDD peripheral blood, IDD vs healthy):\n")
if (nrow(metrn1) > 0) {
  cat("  METRN: log2FC =", round(metrn1$logFC, 4),
      ", raw P =", signif(metrn1$P.Value, 4),
      ", FDR =", signif(metrn1$adj.P.Val, 4),
      ", FDR<0.05?", metrn1$adj.P.Val < 0.05, "\n")
}

cat("\nGSE150408 (IVDD peripheral blood, sciatica vs healthy):\n")
if (nrow(metrn2) > 0) {
  cat("  METRN: log2FC =", round(metrn2$logFC, 4),
      ", raw P =", signif(metrn2$P.Value, 4),
      ", FDR =", signif(metrn2$adj.P.Val, 4),
      ", FDR<0.05?", metrn2$adj.P.Val < 0.05, "\n")
}

cat("\nGSE23130 (IVDD disc tissue, High-grade vs Low-grade):\n")
if (exists("metrn3") && nrow(metrn3) > 0) {
  cat("  METRN: log2FC =", round(metrn3$logFC, 4),
      ", raw P =", signif(metrn3$P.Value, 4),
      ", FDR =", signif(metrn3$adj.P.Val, 4),
      ", FDR<0.05?", metrn3$adj.P.Val < 0.05, "\n")
}

cat("\n========================================\n")
cat("KEY CONCLUSION: After FDR correction,\n")
cat("  GSE124272: raw P=0.003 but FDR=0.32 → NOT significant after correction\n")
cat("  GSE150408: raw P=0.63 → NOT significant (even before correction)\n")
cat("  GSE23130: depends on grade comparison\n")
cat("========================================\n")
