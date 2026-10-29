options(stringsAsFactors = FALSE, warn = 1)
expr_file  <- "data/processed/GSE96058/GSE96058_transformed_expression_primary.rds"
meta_file  <- "data/processed/GSE96058/GSE96058_primary_metadata.tsv"
scores_file <- "data/processed/GSE96058/GSE96058_Hallmark_scores.tsv"
gmt_file   <- "data/references/h.all.v2025.1.Hs.symbols.gmt"
out_dir    <- "reports/GSE96058"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(expr_file, meta_file, scores_file, gmt_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

expr <- as.matrix(readRDS(expr_file))
storage.mode(expr) <- "numeric"
meta <- read.delim(meta_file, check.names = FALSE, na.strings = c("", "NA"))
scores <- read.delim(scores_file, check.names = FALSE, na.strings = c("", "NA"))

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have gene-symbol row names and sample-title column names.")
}
meta <- meta[match(colnames(expr), meta$title), , drop = FALSE]
scores <- scores[match(meta$title, scores$title), , drop = FALSE]
if (anyNA(meta$title) || anyNA(scores$title) ||
    !all(meta$title == colnames(expr)) || !all(scores$title == meta$title)) {
  stop("Sample-title alignment failed among expression, metadata, and score files.")
}

# Read MSigDB Hallmark gene sets.
gmt_lines <- readLines(gmt_file, warn = FALSE)
gmt <- list()
for (line in gmt_lines) {
  f <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(f) >= 3) gmt[[f[1]]] <- unique(f[-c(1, 2)])
}

inflammation_sets <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
estrogen_sets <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)
all_sets <- c(inflammation_sets, estrogen_sets)
if (!all(all_sets %in% names(gmt))) {
  stop("Missing required Hallmark sets in GMT: ",
       paste(setdiff(all_sets, names(gmt)), collapse = ", "))
}
if (!"ESR1" %in% rownames(expr)) stop("ESR1 not found in expression matrix.")
if (!all(inflammation_sets %in% names(scores))) {
  stop("One or more inflammatory score columns are missing from the existing score file.")
}

# Identify the known lower measurement bound from the transformed expression matrix.
# The matrix is already log2(FPKM + 0.1); do not log-transform it again.
floor_value <- log2(0.1)
floor_tolerance <- 1e-7
floor_fraction <- list()
score_reconstruction <- list()
score_floor_audit <- list()

# Platform-relative per-gene z-scoring, matching the project's original scoring method.
score_genes <- function(genes, sample_idx) {
  genes <- intersect(unique(genes), rownames(expr))
  if (length(genes) < 10) stop("Fewer than 10 mapped genes remain in a signature.")
  x <- expr[genes, sample_idx, drop = FALSE]
  gene_z <- t(scale(t(x), center = TRUE, scale = TRUE))
  gene_z[!is.finite(gene_z)] <- NA_real_
  n_variable <- sum(apply(x, 1, function(v) {
    v <- v[is.finite(v)]
    length(v) >= 2 && is.finite(sd(v)) && sd(v) > 0
  }))
  sc <- colMeans(gene_z, na.rm = TRUE)
  sc[!is.finite(sc)] <- NA_real_
  list(score = sc, n_genes = length(genes), n_variable = n_variable)
}

platforms <- unique(meta$platform)
if (!all(c("GPL11154", "GPL18573") %in% platforms)) {
  warning("Expected platforms GPL11154 and GPL18573; found: ", paste(platforms, collapse = ", "))
}

# Outcomes are scored six ways; model_versions also include analyses in which
# inflammatory predictors are floor-filtered as well as the estrogen outcome.
outcome_scores <- list()
predictor_scores <- list()
outcome_versions <- c(
  "original", "overlap_removed", "floor50_outcome_only",
  "floor75_outcome_only", "overlap_removed_floor50",
  "overlap_removed_floor75"
)
model_versions <- c(
  outcome_versions,
  "floor50_all_signatures", "floor75_all_signatures",
  "overlap_removed_floor50_all_signatures",
  "overlap_removed_floor75_all_signatures"
)
shared_inflammation_genes <- unique(unlist(gmt[inflammation_sets], use.names = FALSE))

for (plat in platforms) {
  idx <- which(meta$platform == plat)
  floor_hit <- abs(expr[, idx, drop = FALSE] - floor_value) <= floor_tolerance
  floor_prop <- rowMeans(floor_hit, na.rm = TRUE)
  floor_fraction[[plat]] <- floor_prop

  # Store the original inflammatory scores plus floor-filtered variants.
  for (set_name in inflammation_sets) {
    base_inf_genes <- intersect(gmt[[set_name]], rownames(expr))
    predictor_scores[[paste(plat, set_name, "original", sep = "|")]] <-
      scores[[set_name]][idx]
    for (threshold in c(0.50, 0.75)) {
      floor_genes <- names(floor_prop)[floor_prop >= threshold]
      retained <- setdiff(base_inf_genes, floor_genes)
      sc_inf <- score_genes(retained, idx)
      tag <- if (threshold == 0.50) "floor50" else "floor75"
      predictor_scores[[paste(plat, set_name, tag, sep = "|")]] <- sc_inf$score
      score_floor_audit[[length(score_floor_audit) + 1]] <- data.frame(
        platform = plat, signature_type = "floor_filtered_inflammatory_set",
        gene_set = set_name, floor_threshold = threshold,
        mapped_genes = length(base_inf_genes),
        genes_excluded_for_floor = length(base_inf_genes) - length(retained),
        genes_retained = length(retained),
        variable_retained_genes = sc_inf$n_variable,
        n_samples = length(sc_inf$score),
        n_nonfinite_scores = sum(!is.finite(sc_inf$score))
      )
    }
  }

  # Check whether the provided baseline scores reconstruct from the same expression data.
  for (set_name in all_sets) {
    genes <- intersect(gmt[[set_name]], rownames(expr))
    rebuilt <- score_genes(genes, idx)$score
    provided <- scores[[set_name]][idx]
    finite <- is.finite(rebuilt) & is.finite(provided)
    if (sum(finite) < 10) {
      corr <- NA_real_
      max_abs <- NA_real_
    } else {
      corr <- suppressWarnings(cor(rebuilt[finite], provided[finite]))
      max_abs <- max(abs(rebuilt[finite] - provided[finite]))
    }
    score_reconstruction[[length(score_reconstruction) + 1]] <- data.frame(
      platform = plat, gene_set = set_name,
      n_samples_compared = sum(finite), correlation = corr,
      max_absolute_difference = max_abs
    )

    # Composition audit for all five original gene sets.
    for (threshold in c(0.50, 0.75)) {
      is_floor_dominated <- floor_prop[genes] >= threshold
      remain <- genes[!is_floor_dominated]
      stats <- score_genes(remain, idx)
      score_floor_audit[[length(score_floor_audit) + 1]] <- data.frame(
        platform = plat, signature_type = "original_set",
        gene_set = set_name, floor_threshold = threshold,
        mapped_genes = length(genes),
        genes_excluded_for_floor = sum(is_floor_dominated),
        genes_retained = length(remain),
        variable_retained_genes = stats$n_variable,
        n_samples = length(stats$score),
        n_nonfinite_scores = sum(!is.finite(stats$score))
      )
    }
  }

  for (estrogen_set in estrogen_sets) {
    base_genes <- intersect(gmt[[estrogen_set]], rownames(expr))
    no_overlap_genes <- setdiff(base_genes, shared_inflammation_genes)
    outcome_scores[[paste(plat, estrogen_set, "original", sep = "|")]] <-
      scores[[estrogen_set]][idx]

    base_no_overlap <- score_genes(no_overlap_genes, idx)$score
    outcome_scores[[paste(plat, estrogen_set, "overlap_removed", sep = "|")]] <- base_no_overlap

    for (threshold in c(0.50, 0.75)) {
      # Floor filter is platform-specific because the matrix's floor prevalence can differ by platform.
      floor_genes <- names(floor_prop)[floor_prop >= threshold]
      floor_filtered <- setdiff(base_genes, floor_genes)
      floor_filtered_no_overlap <- setdiff(no_overlap_genes, floor_genes)

      label_floor <- if (threshold == 0.50) "floor50_outcome_only" else "floor75_outcome_only"
      label_combo <- if (threshold == 0.50) "overlap_removed_floor50" else "overlap_removed_floor75"

      s_floor <- score_genes(floor_filtered, idx)
      s_combo <- score_genes(floor_filtered_no_overlap, idx)
      outcome_scores[[paste(plat, estrogen_set, label_floor, sep = "|")]] <- s_floor$score
      outcome_scores[[paste(plat, estrogen_set, label_combo, sep = "|")]] <- s_combo$score

      score_floor_audit[[length(score_floor_audit) + 1]] <- data.frame(
        platform = plat, signature_type = "overlap_removed_estrogen_set",
        gene_set = estrogen_set, floor_threshold = threshold,
        mapped_genes = length(no_overlap_genes),
        genes_excluded_for_floor = length(no_overlap_genes) - length(floor_filtered_no_overlap),
        genes_retained = length(floor_filtered_no_overlap),
        variable_retained_genes = s_combo$n_variable,
        n_samples = length(s_combo$score),
        n_nonfinite_scores = sum(!is.finite(s_combo$score))
      )
    }
  }
}

reconstruction_out <- do.call(rbind, score_reconstruction)
floor_audit_out <- unique(do.call(rbind, score_floor_audit))
write.table(reconstruction_out,
  file.path(out_dir, "GSE96058_score_reconstruction_check.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(floor_audit_out,
  file.path(out_dir, "GSE96058_score_floor_filter_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

reconstruction_bad <- !is.finite(reconstruction_out$correlation) |
  reconstruction_out$correlation < 0.999999 |
  !is.finite(reconstruction_out$max_absolute_difference) |
  reconstruction_out$max_absolute_difference > 1e-6
if (any(reconstruction_bad)) {
  warning("At least one precomputed score did not reconstruct closely. Inspect GSE96058_score_reconstruction_check.tsv before interpreting model results.")
}

# HC3 heteroskedasticity-consistent covariance estimator; implemented directly
# so no additional package is required. This does not correct for dependence
# between observations or unmeasured confounding.
hc3_stats <- function(fit, term = "inflammation_z") {
  X <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  df <- df.residual(fit)
  if (any(1 - h <= .Machine$double.eps^0.5)) {
    stop("Leverage too close to 1 for HC3 adjustment.")
  }
  bread <- solve(crossprod(X))
  adjusted_resid <- e / (1 - h)
  Xw <- X * as.numeric(adjusted_resid)
  meat <- crossprod(Xw)
  V <- bread %*% meat %*% bread
  V <- (V + t(V)) / 2
  if (!term %in% names(coef(fit))) stop("Term not found in model: ", term)
  j <- match(term, names(coef(fit)))
  est <- unname(coef(fit)[j])
  se <- sqrt(max(0, V[j, j]))
  tval <- est / se
  p <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
  crit <- qt(0.975, df = df)
  c(beta = est, HC3_SE = se, HC3_CI_low = est - crit * se,
    HC3_CI_high = est + crit * se, HC3_p_value = p)
}

result_rows <- list()
kk <- 0
for (plat in platforms) {
  platform_idx <- which(meta$platform == plat)
  for (estrogen_set in estrogen_sets) {
    for (predictor in inflammation_sets) {
      # Create one matched complete-case set for all score versions for this outcome/predictor.
      pred_vec <- scores[[predictor]][platform_idx]
      er <- suppressWarnings(as.numeric(as.character(meta$er.status[platform_idx])))
      esr1 <- as.numeric(expr["ESR1", platform_idx])
      age <- suppressWarnings(as.numeric(as.character(meta$age.at.diagnosis[platform_idx])))
      tumor_size <- suppressWarnings(as.numeric(as.character(meta$tumor.size[platform_idx])))
      pam50 <- meta$pam50.subtype[platform_idx]
      all_outcomes <- lapply(outcome_versions, function(v) {
        outcome_scores[[paste(plat, estrogen_set, v, sep = "|")]]
      })
      names(all_outcomes) <- outcome_versions
      all_predictors <- list(
        original = predictor_scores[[paste(plat, predictor, "original", sep = "|")]],
        floor50 = predictor_scores[[paste(plat, predictor, "floor50", sep = "|")]],
        floor75 = predictor_scores[[paste(plat, predictor, "floor75", sep = "|")]]
      )
      common_keep <- !is.na(er) & er == 1 &
        is.finite(esr1) & is.finite(age) & is.finite(tumor_size) & !is.na(pam50)
      for (v in outcome_versions) common_keep <- common_keep & is.finite(all_outcomes[[v]])
      for (v in names(all_predictors)) common_keep <- common_keep & is.finite(all_predictors[[v]])
      common_keep[is.na(common_keep)] <- FALSE

      n_common <- sum(common_keep)
      if (n_common < 30) stop("Too few complete cases for ", plat, " / ", estrogen_set, " / ", predictor)

      outcome_map <- c(
        original = "original",
        overlap_removed = "overlap_removed",
        floor50_outcome_only = "floor50_outcome_only",
        floor75_outcome_only = "floor75_outcome_only",
        overlap_removed_floor50 = "overlap_removed_floor50",
        overlap_removed_floor75 = "overlap_removed_floor75",
        floor50_all_signatures = "floor50_outcome_only",
        floor75_all_signatures = "floor75_outcome_only",
        overlap_removed_floor50_all_signatures = "overlap_removed_floor50",
        overlap_removed_floor75_all_signatures = "overlap_removed_floor75"
      )
      predictor_map <- c(
        original = "original",
        overlap_removed = "original",
        floor50_outcome_only = "original",
        floor75_outcome_only = "original",
        overlap_removed_floor50 = "original",
        overlap_removed_floor75 = "original",
        floor50_all_signatures = "floor50",
        floor75_all_signatures = "floor75",
        overlap_removed_floor50_all_signatures = "floor50",
        overlap_removed_floor75_all_signatures = "floor75"
      )

      for (version in model_versions) {
        out_version <- unname(outcome_map[[version]])
        pred_version <- unname(predictor_map[[version]])
        dat <- data.frame(
          y = all_outcomes[[out_version]][common_keep],
          inflammation = all_predictors[[pred_version]][common_keep],
          ESR1 = esr1[common_keep], age = age[common_keep],
          tumor_size = tumor_size[common_keep], PAM50 = factor(pam50[common_keep])
        )
        dat$PAM50 <- droplevels(dat$PAM50)
        if (nlevels(dat$PAM50) < 2) stop("PAM50 has fewer than two levels in the analysis subset.")
        dat$y_z <- as.numeric(scale(dat$y))
        dat$inflammation_z <- as.numeric(scale(dat$inflammation))
        dat$ESR1_z <- as.numeric(scale(dat$ESR1))
        dat$age_z <- as.numeric(scale(dat$age))
        dat$tumor_size_z <- as.numeric(scale(dat$tumor_size))

        fit <- lm(y_z ~ inflammation_z + ESR1_z + age_z + tumor_size_z + PAM50, data = dat)
        co <- summary(fit)$coefficients
        classic_ci <- confint(fit, "inflammation_z", level = 0.95)
        robust <- hc3_stats(fit, "inflammation_z")

        kk <- kk + 1
        result_rows[[kk]] <- data.frame(
          platform = plat, estrogen_outcome = estrogen_set,
          inflammatory_predictor = predictor, score_version = version,
          n = nrow(dat), beta_standardized = unname(co["inflammation_z", "Estimate"]),
          OLS_SE = unname(co["inflammation_z", "Std. Error"]),
          OLS_CI95_low = unname(classic_ci[1]), OLS_CI95_high = unname(classic_ci[2]),
          OLS_p_value = unname(co["inflammation_z", "Pr(>|t|)"]),
          HC3_SE = unname(robust["HC3_SE"]),
          HC3_CI95_low = unname(robust["HC3_CI_low"]),
          HC3_CI95_high = unname(robust["HC3_CI_high"]),
          HC3_p_value = unname(robust["HC3_p_value"]),
          R_squared = summary(fit)$r.squared,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

out <- do.call(rbind, result_rows)
out$OLS_FDR_within_platform_version <- NA_real_
out$HC3_FDR_within_platform_version <- NA_real_
for (plat in unique(out$platform)) {
  for (version in unique(out$score_version)) {
    ix <- which(out$platform == plat & out$score_version == version)
    out$OLS_FDR_within_platform_version[ix] <- p.adjust(out$OLS_p_value[ix], method = "BH")
    out$HC3_FDR_within_platform_version[ix] <- p.adjust(out$HC3_p_value[ix], method = "BH")
  }
}

out <- out[order(out$platform, out$score_version, out$estrogen_outcome,
                 match(out$inflammatory_predictor, inflammation_sets)), ]
write.table(out,
  file.path(out_dir, "GSE96058_HC3_floor_sensitivity_models.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Concise comparison for the six main models per platform, original vs overlap-removed,
# focusing on whether significance and direction persist under HC3 inference.
main <- out[out$score_version %in% c("original", "overlap_removed"), ]
write.table(main,
  file.path(out_dir, "GSE96058_HC3_primary_models.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\n=== Score reconstruction check ===\n")
print(reconstruction_out, row.names = FALSE)
cat("\n=== Floor-filter score audit ===\n")
print(floor_audit_out, row.names = FALSE)
cat("\n=== HC3 and floor-sensitivity models (selected columns) ===\n")
print(out[, c("platform", "estrogen_outcome", "inflammatory_predictor",
              "score_version", "n", "beta_standardized",
              "OLS_FDR_within_platform_version", "HC3_CI95_low",
              "HC3_CI95_high", "HC3_p_value", "HC3_FDR_within_platform_version")],
      row.names = FALSE, digits = 4)
cat("\nSaved files:\n")
cat(file.path(out_dir, "GSE96058_score_reconstruction_check.tsv"), "\n")
cat(file.path(out_dir, "GSE96058_score_floor_filter_audit.tsv"), "\n")
cat(file.path(out_dir, "GSE96058_HC3_floor_sensitivity_models.tsv"), "\n")
cat(file.path(out_dir, "GSE96058_HC3_primary_models.tsv"), "\n")
