options(stringsAsFactors = FALSE)

base <- "reports/GSE96058"

# Main models: all six tests across both platforms (12 total).
f <- file.path(base, "GSE96058_adjusted_inflammation_estrogen_models.tsv")
if (file.exists(f)) {
  x <- read.delim(f, check.names = FALSE)
  x$FDR_across_platforms <- p.adjust(x$p_value, method = "BH")
  write.table(x, file.path(base, "GSE96058_main_models_broader_FDR.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Main models:", nrow(x), "tests; broader FDR written.\n")
}

# Standardized models: correct across platforms within each score version.
f <- file.path(base, "GSE96058_standardized_outcome_models.tsv")
if (file.exists(f)) {
  x <- read.delim(f, check.names = FALSE)
  x$FDR_across_platforms_version <- NA_real_
  for (v in unique(x$score_version)) {
    i <- which(x$score_version == v)
    x$FDR_across_platforms_version[i] <- p.adjust(x$p_value[i], method = "BH")
  }
  write.table(x, file.path(base, "GSE96058_standardized_models_broader_FDR.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Standardized models:", nrow(x), "rows; broader FDR written.\n")
}

# HC3 models: across platforms within each score version, separately for OLS and HC3.
for (nm in c("GSE96058_HC3_primary_models.tsv",
             "GSE96058_HC3_floor_sensitivity_models.tsv")) {
  f <- file.path(base, nm)
  if (!file.exists(f)) next
  x <- read.delim(f, check.names = FALSE)
  x$OLS_FDR_across_platforms_version <- NA_real_
  x$HC3_FDR_across_platforms_version <- NA_real_
  for (v in unique(x$score_version)) {
    i <- which(x$score_version == v)
    x$OLS_FDR_across_platforms_version[i] <-
      p.adjust(x$OLS_p_value[i], method = "BH")
    x$HC3_FDR_across_platforms_version[i] <-
      p.adjust(x$HC3_p_value[i], method = "BH")
  }
  out <- sub("\\.tsv$", "_broader_FDR.tsv", nm)
  write.table(x, file.path(base, out), sep = "\t",
              quote = FALSE, row.names = FALSE)
  cat(nm, ":", nrow(x), "rows; broader FDR written.\n")
}

# PAM50 interactions: correct across platforms within each score version.
f <- file.path(base, "GSE96058_PAM50_interaction_diagnostics.tsv")
if (file.exists(f)) {
  x <- read.delim(f, check.names = FALSE)
  x$interaction_FDR_across_platforms_version <- NA_real_
  for (v in unique(x$score_version)) {
    i <- which(x$score_version == v)
    x$interaction_FDR_across_platforms_version[i] <-
      p.adjust(x$interaction_p[i], method = "BH")
  }
  write.table(x, file.path(base, "GSE96058_PAM50_interactions_broader_FDR.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("PAM50 interactions:", nrow(x), "rows; broader FDR written.\n")
}
