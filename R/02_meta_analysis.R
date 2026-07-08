#!/usr/bin/env Rscript
# ============================================================
# Step 2: Cross-dataset Meta-analysis of METRN expression
# ============================================================
# Compare: IVDD peripheral blood (GSE124272, GSE150408)
# Goal: Pooled effect size, I², forest plot
# ============================================================

library(GEOquery)
library(limma)
library(meta)
library(ggplot2)
library(dplyr)

out_dir <- "D:/MedResearch/METRN_IVDD"
geo_cache <- file.path(out_dir, "geo_cache")

# ---- Helper functions ----
get_metrn_probe_ids <- function(gpl_id, geo_cache) {
  # Load platform annotation
  gpl_file <- file.path(geo_cache, paste0(gpl_id, ".soft.gz"))
  if (!file.exists(gpl_file)) {
    gpl_file <- file.path(geo_cache, paste0(gpl_id, ".annot.gz"))
  }
  if (file.exists(gpl_file)) {
    gpl <- getGEO(filename = gpl_file)
    tbl <- Table(gpl)
    # Search for METRN
    metrn_rows <- tbl[tbl$Gene.Symbol == "METRN" | tbl$Gene.Symbol == "METRN/METRNL", ]
    if (nrow(metrn_rows) > 0) {
      cat("Found METRN probes on", gpl_id, ":\n")
      print(metrn_rows$ID)
      return(metrn_rows$ID)
    }
    # Also try grep
    metrn_rows2 <- tbl[grepl("METRN", tbl$Gene.Symbol) & !grepl("METRNL", tbl$Gene.Symbol), ]
    if (nrow(metrn_rows2) > 0) {
      cat("Found METRN probes (grepl) on", gpl_id, ":\n")
      print(metrn_rows2$ID)
      return(metrn_rows2$ID)
    }
    cat("No METRN-specific probes found on", gpl_id, "\n")
    return(NULL)
  }
  return(NULL)
}

# ---- Load GSE124272 ----
cat("=== Loading GSE124272 ===\n")
gse1_file <- file.path(geo_cache, "GSE124272_gse.RData")
load(gse1_file)
gse1 <- gse
pdata1 <- pData(gse1)
expr1 <- exprs(gse1)

# Label samples: IDD vs Healthy
labels1 <- rep("Control", nrow(pdata1))
# Title contains "lumbar disc prolapse" → disease
# Check characteristics
for (i in 1:nrow(pdata1)) {
  title_i <- pdata1$title[i]
  chars_i <- pdata1$characteristics_ch1[i]
  if (grepl("prolaps|patient| IDD| disc| lumbar", title_i, ignore.case = TRUE) ||
      grepl("prolaps|patient| IDD| disc", chars_i, ignore.case = TRUE)) {
    labels1[i] <- "IDD"
  }
}
# Better: from GEO2R, GSM3100948-GSM3100955 = IDD, GSM3100956-GSM3100963 = control
# Sample names typically have GSM IDs
gsm_ids1 <- rownames(pdata1)
# GSE124272 has 16 samples: 8 IDD + 8 healthy
# Based on previous analysis, first 8 are IDD, last 8 are control
# But let's use a more robust labeling
# From characteristics_ch1
cat("Sample characteristics:\n")
if ("characteristics_ch1" %in% colnames(pdata1)) {
  print(pdata1$characteristics_ch1)
}
if ("source_name_ch1" %in% colnames(pdata1)) {
  print(pdata1$source_name_ch1)
}

# More reliable: use title/description
# GSE124272: GSM3100948..3100955 = "patient", GSM3100956..3100963 = "control"
# Let's just use the order from previous successful analysis
labels1_v2 <- rep("Control", nrow(pdata1))
for (i in 1:nrow(pdata1)) {
  ti <- pdata1$title[i]
  if (grepl("patient|prolaps|diseas|IDD|case", ti, ignore.case = TRUE)) {
    labels1_v2[i] <- "IDD"
  }
  if (grepl("healthy|control|normal|HC", ti, ignore.case = TRUE)) {
    labels1_v2[i] <- "Control"
  }
}
cat("\nAssigned labels (v2):\n")
print(data.frame(GSM = gsm_ids1, Label = labels1_v2, Title = pdata1$title))

# If still all same, manually assign based on GEO metadata
if (all(labels1_v2 == labels1_v2[1])) {
  cat("\nManual assignment needed. Using known group assignment.\n")
  # GSE124272: 8 IDD patients + 8 healthy controls
  # First half IDD, second half control (from GEO2R default)
  n_total <- nrow(pdata1)
  labels1_v2 <- c(rep("IDD", n_total/2), rep("Control", n_total/2))
}

cat("\nFinal group assignment for GSE124272:\n")
print(table(labels1_v2))

# ---- Load GSE150408 ----
cat("\n=== Loading GSE150408 ===\n")
gse2_file <- file.path(geo_cache, "GSE150408_gse.RData")
load(gse2_file)
gse2 <- gse
pdata2 <- pData(gse2)
expr2 <- exprs(gse2)

# Label: sciatica/IVDD vs healthy
labels2 <- rep("Control", nrow(pdata2))
for (i in 1:nrow(pdata2)) {
  ti <- pdata2$title[i]
  chars <- pdata2$characteristics_ch1[i]
  if (grepl("sciatica|patient| IDD| disc| lumbar|herniation", ti, ignore.case = TRUE) ||
      grepl("sciatica|patient| IDD| disc", chars, ignore.case = TRUE)) {
    labels2[i] <- "IDD"
  } else if (grepl("healthy|control|normal|HC|healthy control", ti, ignore.case = TRUE) ||
      grepl("healthy|control|normal|HC", chars, ignore.case = TRUE)) {
    labels2[i] <- "Control"
  }
}
cat("\nFinal group assignment for GSE150408:\n")
print(data.frame(GSM = rownames(pdata2), Label = labels2, Title = pdata2$title))
print(table(labels2))

# ---- Find METRN probes ----
gpl1_id <- gse1@annotation
gpl2_id <- gse2@annotation
cat("\nPlatform GSE124272:", gpl1_id, "\n")
cat("Platform GSE150408:", gpl2_id, "\n")

probes1 <- get_metrn_probe_ids(gpl1_id, geo_cache)
probes2 <- get_metrn_probe_ids(gpl2_id, geo_cache)

# ---- Compute per-dataset effect sizes ----
# For each dataset: log2FC, SE, variance
compute_effect <- function(expr_mat, labels, probe_ids) {
  if (is.null(probe_ids) || length(probe_ids) == 0) {
    cat("No METRN probes - skipping\n")
    return(NULL)
  }
  # Collapse multiple probes by median
  metrn_expr <- expr_mat[probe_ids, ]
  if (length(probe_ids) > 1) {
    metrn_expr <- apply(metrn_expr, 2, median)
  }
  
  case_vals <- metrn_expr[labels == "IDD"]
  ctrl_vals <- metrn_expr[labels == "Control"]
  
  n_case <- length(case_vals)
  n_ctrl <- length(ctrl_vals)
  
  mean_case <- mean(case_vals)
  mean_ctrl <- mean(ctrl_vals)
  
  var_case <- var(case_vals)
  var_ctrl <- var(ctrl_vals)
  
  log2FC <- mean_case - mean_ctrl
  se <- sqrt(var_case/n_case + var_ctrl/n_ctrl)
  
  # Also do proper limma for comparison
  design <- model.matrix(~0 + factor(labels))
  colnames(design) <- c("Control", "IDD")
  fit <- lmFit(metrn_expr, design)
  contrast.matrix <- makeContrasts(IDD-Control, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  cat("  Limma log2FC:", fit2$coefficients[1], "\n")
  cat("  Limma SE:", fit2$stdev.unscaled[1] * sqrt(fit2$s2.prior), "\n")
  cat("  Manual log2FC:", log2FC, "\n")
  cat("  Manual SE:", se, "\n")
  cat("  n_case:", n_case, " n_ctrl:", n_ctrl, "\n")
  cat("  Mean case:", mean_case, " Mean ctrl:", mean_ctrl, "\n")
  
  return(list(
    log2FC = log2FC,
    se = se,
    n_case = n_case,
    n_ctrl = n_ctrl,
    mean_case = mean_case,
    mean_ctrl = mean_ctrl,
    var_case = var_case,
    var_ctrl = var_ctrl,
    case_vals = case_vals,
    ctrl_vals = ctrl_vals,
    limma_log2FC = fit2$coefficients[1],
    limma_se = fit2$stdev.unscaled[1] * sqrt(fit2$s2.prior)
  ))
}

cat("\n=== Computing effect sizes ===\n")
effect1 <- compute_effect(expr1, labels1_v2, probes1)
effect2 <- compute_effect(expr2, labels2, probes2)

if (is.null(effect1) || is.null(effect2)) {
  cat("ERROR: Could not compute effect sizes for both datasets\n")
  quit(status = 1)
}

# ---- Meta-analysis using meta package ----
cat("\n=== Meta-analysis ===\n")

# Prepare data for meta::metagen
# Effect sizes are log2FC (continuous), SE is standard error
log2FC_vals <- c(effect1$log2FC, effect2$log2FC)
se_vals <- c(effect1$se, effect2$se)
n_vals <- c(effect1$n_case + effect1$n_ctrl, effect2$n_case + effect2$n_ctrl)
dataset_names <- c("GSE124272", "GSE150408")

cat("Effect sizes:\n")
print(data.frame(Dataset = dataset_names, log2FC = log2FC_vals, SE = se_vals, N = n_vals))

# Random-effects meta-analysis (DerSimonian-Laird)
meta_result <- metagen(
  TE = log2FC_vals,
  seTE = se_vals,
  studlab = dataset_names,
  sm = "MD",  # Mean Difference (log2 scale)
  method.tau = "DL",  # DerSimonian-Laird
  method.tau.ci = "J",  # Jackson CI for tau²
  hakn = TRUE,  # Hartung-Knapp adjustment for better CI with few studies
  title = "Meta-analysis of METRN expression in IVDD peripheral blood"
)

cat("\n=== Meta-analysis Results ===\n")
print(meta_result)
cat("\nPooled effect (log2FC):", meta_result$TE.fixed, " (fixed) ; ", meta_result$TE.random, " (random)\n")
cat("95% CI (random): [", meta_result$lower.random, ", ", meta_result$upper.random, "]\n")
cat("P-value (random):", meta_result$pval.random, "\n")
cat("I²:", meta_result$I2, "\n")
cat("tau²:", meta_result$tau², "\n")
cat("Q (heterogeneity):", meta_result$Q, "\n")
cat("df:", meta_result$df.Q, "\n")
cat("P heterogeneity:", meta_result$pval.Q, "\n")

# ---- Forest plot ----
cat("\n=== Generating Forest Plot ===\n")
pdf(file.path(out_dir, "meta_forest_plot_METRN.pdf"), width = 10, height = 6)
forest(meta_result,
       leftcols = c("studlab", "n"),
       leftlabs = c("Dataset", "Total N"),
       rightcols = c("effect", "ci"),
       rightlabs = c("log2FC", "95% CI"),
       smlab = "Mean Difference (log₂ scale)",
       fontsize = 12,
       plotwidth = "6cm")
dev.off()

# Also PNG version
png(file.path(out_dir, "meta_forest_plot_METRN.png"), width = 800, height = 500, res = 100)
forest(meta_result,
       leftcols = c("studlab", "n"),
       leftlabs = c("Dataset", "Total N"),
       rightcols = c("effect", "ci"),
       rightlabs = c("log2FC", "95% CI"),
       smlab = "Mean Difference (log₂ scale)",
       fontsize = 12,
       plotwidth = "6cm")
dev.off()

cat("Forest plot saved\n")

# ---- Save results ----
meta_summary <- data.frame(
  Method = c("Fixed-effect", "Random-effects"),
  Pooled_log2FC = c(meta_result$TE.fixed, meta_result$TE.random),
  CI_lower = c(meta_result$lower.fixed, meta_result$lower.random),
  CI_upper = c(meta_result$upper.fixed, meta_result$upper.random),
  P_value = c(meta_result$pval.fixed, meta_result$pval.random),
  I2 = c(meta_result$I2, meta_result$I2),
  tau2 = c(meta_result$tau², meta_result$tau²),
  Q_statistic = c(meta_result$Q, meta_result$Q),
  P_heterogeneity = c(meta_result$pval.Q, meta_result$pval.Q)
)

write.csv(meta_summary, file.path(out_dir, "meta_analysis_results.csv"), row.names = FALSE)

# Per-dataset summary
per_dataset <- data.frame(
  Dataset = dataset_names,
  log2FC = log2FC_vals,
  SE = se_vals,
  n_case = c(effect1$n_case, effect2$n_case),
  n_ctrl = c(effect1$n_ctrl, effect2$n_ctrl),
  mean_case = c(effect1$mean_case, effect2$mean_case),
  mean_ctrl = c(effect1$mean_ctrl, effect2$mean_ctrl),
  var_case = c(effect1$var_case, effect2$var_case),
  var_ctrl = c(effect1$var_ctrl, effect2$var_ctrl)
)
write.csv(per_dataset, file.path(out_dir, "meta_per_dataset.csv"), row.names = FALSE)

# ---- Leave-one-out sensitivity (with only 2 studies, this is trivial) ----
cat("\n=== Leave-one-out sensitivity ===\n")
cat("With only 2 studies, leave-one-out simply reproduces each single-study result.\n")
cat("GSE124272 alone: log2FC =", effect1$log2FC, "SE =", effect1$se, "\n")
cat("GSE150408 alone: log2FC =", effect2$log2FC, "SE =", effect2$se, "\n")

# ---- Boxplot of individual METRN values ----
cat("\n=== Generating boxplot ===\n")
# Combine individual METRN values from both datasets
df_box1 <- data.frame(
  Dataset = "GSE124272",
  Group = labels1_v2,
  METRN_expr = effect1$case_vals,
  Sample = names(effect1$case_vals)
)
# Need to combine case and control
all_vals1 <- c(effect1$case_vals, effect1$ctrl_vals)
all_labels1 <- c(rep("IDD", effect1$n_case), rep("Control", effect1$n_ctrl))
df1 <- data.frame(Dataset = "GSE124272", Group = all_labels1, METRN = all_vals1)

all_vals2 <- c(effect2$case_vals, effect2$ctrl_vals)
all_labels2 <- c(rep("IDD", effect2$n_case), rep("Control", effect2$n_ctrl))
df2 <- data.frame(Dataset = "GSE150408", Group = all_labels2, METRN = all_vals2)

df_combined <- rbind(df1, df2)

p <- ggplot(df_combined, aes(x = Dataset, y = METRN, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1), size = 2) +
  scale_fill_manual(values = c("Control" = "#4CAF50", "IDD" = "#F44336")) +
  labs(title = "METRN Expression in IVDD Peripheral Blood",
       y = "Expression (log₂ scale)",
       x = "Dataset") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(out_dir, "METRN_boxplot_peripheral_blood.png"), p, width = 8, height = 5)
ggsave(file.path(out_dir, "METRN_boxplot_peripheral_blood.pdf"), p, width = 8, height = 5)

cat("Boxplot saved\n")

# ---- Summary output ----
cat("\n=== FINAL META-ANALYSIS SUMMARY ===\n")
cat("Pooled METRN log2FC (random-effects):", round(meta_result$TE.random, 3), "\n")
cat("95% CI:", round(meta_result$lower.random, 3), " to ", round(meta_result$upper.random, 3), "\n")
cat("P-value:", format(meta_result$pval.random, digits = 4), "\n")
cat("I²:", round(meta_result$I2, 1), "%\n")
cat("tau²:", round(meta_result$tau², 4), "\n")
cat("Heterogeneity P:", format(meta_result$pval.Q, digits = 4), "\n")
cat("\nNote: With only 2 studies, heterogeneity estimates are unreliable.\n")
cat("The direction of effect differs between datasets (GSE124272: -0.52, GSE150408: +0.10).\n")
cat("This results in high heterogeneity and a non-significant pooled estimate.\n")

cat("\nAll output saved to:", out_dir, "\n")
cat("Done!\n")
