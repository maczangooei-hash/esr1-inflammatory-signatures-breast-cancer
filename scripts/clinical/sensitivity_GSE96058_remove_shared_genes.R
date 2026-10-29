options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(limma))

expr <- readRDS(
  "data/processed/GSE96058/GSE96058_transformed_expression_primary.rds"
)

meta <- read.delim(
  "data/processed/GSE96058/GSE96058_primary_metadata.tsv",
  check.names = FALSE
)

scores <- read.delim(
  "data/processed/GSE96058/GSE96058_Hallmark_scores.tsv",
  check.names = FALSE
)

gmt_lines <- readLines(
  "data/references/h.all.v2025.1.Hs.symbols.gmt",
  warn = FALSE
)

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

shared_inflammation_genes <- unique(
  unlist(gmt[inflammation_sets], use.names = FALSE)
)

# Remove all genes belonging to any inflammatory signature from
# each estrogen signature, then calculate platform-specific gene z-scores.
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

if (is.null(rownames(expr)) || is.null(colnames(expr))) {
  stop("Expression matrix must have gene row names and sample column names.")
}

meta <- meta[match(colnames(expr), meta$title), , drop = FALSE]
if (anyNA(meta$title) || any(meta$title != colnames(expr))) {
  stop("Metadata titles do not align with expression columns.")
}

scores <- scores[match(meta$title, scores$title), , drop = FALSE]
if (anyNA(scores$title) || any(scores$title != meta$title)) {
  stop("Hallmark score rows do not align with expression samples.")
}

results <- list()
overlap_audit <- list()
k <- 0
a <- 0

for (estrogen_set in estrogen_sets) {
  original_genes <- intersect(gmt[[estrogen_set]], rownames(expr))
  retained_genes <- setdiff(original_genes, shared_inflammation_genes)

  if (length(retained_genes) < 20) {
    stop("Too few genes remain for ", estrogen_set)
  }

  a <- a + 1
  overlap_audit[[a]] <- data.frame(
    estrogen_set = estrogen_set,
    original_genes_present = length(original_genes),
    inflammatory_genes_removed = length(setdiff(original_genes, retained_genes)),
    retained_genes = length(retained_genes)
  )

  for (platform in unique(meta$platform)) {
    idx <- which(meta$platform == platform)
    x <- expr[retained_genes, idx, drop = FALSE]

    # Match the original scoring strategy: standardize each gene
    # across samples within each platform, then average z-scores.
    gene_z <- t(scale(t(x)))
    gene_z[!is.finite(gene_z)] <- NA_real_
    new_score <- colMeans(gene_z, na.rm = TRUE)

    if (!paste0(estrogen_set, "_noInflammatoryOverlap") %in% names(scores)) {
      scores[[paste0(estrogen_set, "_noInflammatoryOverlap")]] <- NA_real_
    }

    # Assign by sample title, not by assuming metadata and score row order.
    score_rows <- match(colnames(expr)[idx], scores$title)
    if (anyNA(score_rows) || length(score_rows) != length(new_score)) {
      stop("Sample alignment failure while assigning revised estrogen scores.")
    }

    scores[[paste0(estrogen_set, "_noInflammatoryOverlap")]][score_rows] <-
      new_score
  }
}

write.table(
  do.call(rbind, overlap_audit),
  "reports/GSE96058/GSE96058_no_overlap_gene_audit.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Standardize ESR1 within each platform, consistent with the original models.
esr1 <- expr["ESR1", ]
meta$ESR1_z <- NA_real_
for (platform in unique(meta$platform)) {
  idx <- which(meta$platform == platform)
  meta$ESR1_z[idx] <- as.numeric(scale(esr1[idx]))
}

# Model each inflammatory score against the revised estrogen outcomes.
# Keep the original covariates and complete-case selection.
for (platform in unique(meta$platform)) {
  for (outcome_set in estrogen_sets) {
    outcome <- paste0(outcome_set, "_noInflammatoryOverlap")

    for (inflammation_set in inflammation_sets) {
      inflammation_col <- inflammation_set
      if (!inflammation_col %in% names(scores)) {
        stop("Missing inflammatory score column: ", inflammation_col)
      }

      dat <- data.frame(
        outcome = scores[[outcome]][meta$platform == platform],
        inflammation = scores[[inflammation_col]][meta$platform == platform],
        ESR1_z = meta$ESR1_z[meta$platform == platform],
        age = meta$age.at.diagnosis[meta$platform == platform],
        tumor_size = meta$tumor.size[meta$platform == platform],
        PAM50 = meta$pam50.subtype[meta$platform == platform],
        er_status = meta$er.status[meta$platform == platform]
      )

      # Handle numeric fields that may have been imported as text.
      dat$age <- suppressWarnings(as.numeric(as.character(dat$age)))
      dat$tumor_size <- suppressWarnings(as.numeric(as.character(dat$tumor_size)))

      # Reproduce the original regression cohort:
      # clinically ER-positive, with complete ESR1, age, tumor size,
      # PAM50, predictor and outcome data.
      platform_idx <- which(meta$platform == platform)

      keep <- (
        !is.na(meta$er.status[platform_idx]) &
        meta$er.status[platform_idx] == 1 &
        !is.na(meta$age.at.diagnosis[platform_idx]) &
        !is.na(meta$tumor.size[platform_idx]) &
        !is.na(meta$pam50.subtype[platform_idx]) &
        is.finite(scores[[outcome]][platform_idx]) &
        is.finite(scores[[inflammation_col]][platform_idx]) &
        is.finite(meta$ESR1_z[platform_idx])
      )

      # Never allow NA values in a logical row-selection vector.
      keep[is.na(keep)] <- FALSE

      # Match the same platform rows in dat, which was built in platform order.
      dat <- dat[keep, , drop = FALSE]

      # Match the original model: standardize the inflammatory
      # predictor within the platform-specific complete-case cohort.
      dat$inflammation_z <- as.numeric(scale(dat$inflammation))

      dat$PAM50 <- droplevels(factor(dat$PAM50))

      if (nrow(dat) == 0 || nlevels(dat$PAM50) < 2) {
        stop(
          "Insufficient complete cases or PAM50 levels for ",
          platform, " / ", outcome_set, " / ", inflammation_set
        )
      }

      fit <- lm(
        outcome ~ inflammation_z + ESR1_z + age + tumor_size + PAM50,
        data = dat
      )

      sm <- summary(fit)$coefficients
      ci <- confint(fit, "inflammation_z", level = 0.95)

      k <- k + 1
      results[[k]] <- data.frame(
        platform = platform,
        estrogen_outcome = outcome_set,
        inflammatory_predictor = inflammation_set,
        n = nrow(dat),
        beta = unname(sm["inflammation_z", "Estimate"]),
        SE = unname(sm["inflammation_z", "Std. Error"]),
        CI_low = ci[1],
        CI_high = ci[2],
        p_value = unname(sm["inflammation_z", "Pr(>|t|)"])
      )
    }
  }
}

out <- do.call(rbind, results)
out$FDR_within_platform <- NA_real_

for (platform in unique(out$platform)) {
  idx <- which(out$platform == platform)
  out$FDR_within_platform[idx] <- p.adjust(out$p_value[idx], method = "BH")
}

write.table(
  out,
  "reports/GSE96058/GSE96058_adjusted_models_no_inflammatory_overlap.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(do.call(rbind, overlap_audit))
cat("\nSensitivity model results:\n")
print(out, row.names = FALSE)

cat("\nSaved sensitivity results to:\n")
cat("reports/GSE96058/GSE96058_adjusted_models_no_inflammatory_overlap.tsv\n")
