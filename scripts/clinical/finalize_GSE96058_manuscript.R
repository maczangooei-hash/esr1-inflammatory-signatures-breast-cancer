options(stringsAsFactors = FALSE)

base <- "reports/GSE96058"
out_file <- file.path(base, "GSE96058_final_manuscript_audit.txt")

read_result <- function(file) {
  path <- file.path(base, file)
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }
  read.delim(path, check.names = FALSE)
}

lines <- c(
  "GSE96058 MANUSCRIPT RESULTS AUDIT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Interpretation note:",
  paste(
    "GPL11154 and GPL18573 are platform-defined subsets of GSE96058,",
    "not independent external patient cohorts."
  ),
  ""
)

# 1. Main adjusted models
x <- read_result("GSE96058_main_models_broader_FDR.tsv")
if (!is.null(x)) {
  lines <- c(lines, "=== PRIMARY ADJUSTED MODELS ===")
  lines <- c(lines, capture.output(
    print(x[, c(
      "platform", "outcome", "inflammatory_score", "n",
      "beta_per_SD", "CI95_low", "CI95_high",
      "p_value", "FDR_across_platforms"
    )], row.names = FALSE)
  ))
  lines <- c(lines, "")
}

# 2. HC3 primary models
x <- read_result("GSE96058_HC3_primary_models_broader_FDR.tsv")
if (!is.null(x)) {
  lines <- c(lines, "=== HC3 ROBUSTNESS: PRIMARY SCORES ===")
  lines <- c(lines, capture.output(
    print(x[, c(
      "platform", "estrogen_outcome", "inflammatory_predictor",
      "score_version", "n", "beta_standardized",
      "HC3_CI95_low", "HC3_CI95_high",
      "HC3_p_value", "HC3_FDR_across_platforms_version"
    )], row.names = FALSE)
  ))
  lines <- c(lines, "")
}

# 3. Shared-gene-removal sensitivity
x <- read_result("GSE96058_HC3_floor_sensitivity_models_broader_FDR.tsv")
if (!is.null(x)) {
  x <- x[x$score_version == "overlap_removed", ]
  lines <- c(lines, "=== HC3 SENSITIVITY: SHARED GENES REMOVED ===")
  lines <- c(lines, capture.output(
    print(x[, c(
      "platform", "estrogen_outcome", "inflammatory_predictor",
      "n", "beta_standardized", "HC3_CI95_low",
      "HC3_CI95_high", "HC3_p_value",
      "HC3_FDR_across_platforms_version"
    )], row.names = FALSE)
  ))
  lines <- c(lines, "")
}

# 4. PAM50 interaction tests
x <- read_result("GSE96058_PAM50_interactions_broader_FDR.tsv")
if (!is.null(x)) {
  lines <- c(lines, "=== PAM50 INTERACTION TESTS ===")
  keep <- intersect(c(
    "platform", "score_version", "outcome",
    "inflammatory_predictor", "interaction_p",
    "interaction_FDR", "interaction_FDR_across_platforms_version"
  ), names(x))
  lines <- c(lines, capture.output(print(x[, keep], row.names = FALSE)))
  lines <- c(lines, "")
}

# 5. Subtype-specific slope table
x <- read_result("GSE96058_manuscript_PAM50_slopes.tsv")
if (!is.null(x)) {
  lines <- c(lines, "=== PAM50 SUBTYPE-SPECIFIC SLOPES ===")
  lines <- c(lines, capture.output(
    print(utils::head(x, 40), row.names = FALSE)
  ))
  lines <- c(lines, "")
}

lines <- c(
  lines,
  "=== INTERPRETATION GUARDRAILS ===",
  "- These are associations, not causal effects.",
  "- Gene-set scores are transcriptomic proxies, not direct ER protein activity measurements.",
  "- The two platforms are not independent patient cohorts.",
  "- Inflammatory scores are highly correlated; interpret them as related signatures.",
  "- Shared-gene-removal sensitivity supports robustness in GPL11154,",
  "  but replication is weaker in GPL18573.",
  "- PAM50 interaction results do not alone identify the subtype driving heterogeneity.",
  "- Report the multiple-testing family explicitly in Methods.",
  ""
)

writeLines(lines, out_file)
cat("Audit report written to:\n", out_file, "\n")
cat("Lines:", length(lines), "\n")
