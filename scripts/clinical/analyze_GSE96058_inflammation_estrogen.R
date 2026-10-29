options(stringsAsFactors = FALSE)
input_file <- "data/processed/GSE96058/GSE96058_Hallmark_scores.tsv"
out_dir <- "reports/GSE96058"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  input_file, check.names = FALSE, na.strings = c("", "NA")
)

# Restrict to clinically ER-positive tumors.
dat <- dat[!is.na(dat$er.status) & dat$er.status == 1, , drop = FALSE]

inflammation_scores <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

estrogen_scores <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)

required <- c(
  "platform", "ESR1_expression", "age.at.diagnosis",
  "tumor.size", "pam50.subtype",
  inflammation_scores, estrogen_scores
)
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Convert clinical numeric variables explicitly.
dat$age_num <- suppressWarnings(as.numeric(dat$age.at.diagnosis))
dat$tumor_size_num <- suppressWarnings(as.numeric(dat$tumor.size))

# Fit separate models within each platform.
platforms <- c("GPL11154", "GPL18573")
results <- list()
sample_flow <- list()

for (plat in platforms) {
  d <- dat[dat$platform == plat, , drop = FALSE]

  # Standardize continuous predictors within this platform.
  d$ESR1_z <- as.numeric(scale(d$ESR1_expression))
  d$age_z <- as.numeric(scale(d$age_num))
  d$tumor_size_z <- as.numeric(scale(d$tumor_size_num))

  d$pam50_subtype <- factor(d$pam50.subtype)

  sample_flow[[plat]] <- data.frame(
    platform = plat,
    clinically_ER_positive = nrow(d),
    complete_covariates = sum(
      complete.cases(d[, c(
        "ESR1_z", "age_z", "tumor_size_z", "pam50_subtype"
      )])
    )
  )

  for (outcome in estrogen_scores) {
    for (predictor in inflammation_scores) {
      vars <- c(
        outcome, predictor, "ESR1_z", "age_z",
        "tumor_size_z", "pam50_subtype"
      )
      dd <- d[complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]

      # Scale the inflammatory predictor so the coefficient is per SD.
      dd$inflammation_z <- as.numeric(scale(dd[[predictor]]))

      if (nrow(dd) < 30) {
        warning("Too few complete cases for ", plat, " / ", predictor)
        next
      }

      formula <- as.formula(paste(
        "`", outcome, "` ~ inflammation_z + ESR1_z + age_z + ",
        "tumor_size_z + pam50_subtype", sep = ""
      ))

      fit <- lm(formula, data = dd)
      sm <- summary(fit)$coefficients
      ci <- confint(fit, "inflammation_z", level = 0.95)

      results[[length(results) + 1]] <- data.frame(
        platform = plat,
        outcome = outcome,
        inflammatory_score = predictor,
        n = nrow(dd),
        beta_per_SD = unname(sm["inflammation_z", "Estimate"]),
        SE = unname(sm["inflammation_z", "Std. Error"]),
        CI95_low = unname(ci[1]),
        CI95_high = unname(ci[2]),
        p_value = unname(sm["inflammation_z", "Pr(>|t|)"]),
        ESR1_beta = unname(sm["ESR1_z", "Estimate"]),
        R_squared = summary(fit)$r.squared,
        adjusted_R_squared = summary(fit)$adj.r.squared
      )
    }
  }
}

res <- do.call(rbind, results)

# Adjust for the six inflammation-by-estrogen-score tests within each platform.
res$FDR_within_platform <- NA_real_
for (plat in unique(res$platform)) {
  ix <- which(res$platform == plat)
  res$FDR_within_platform[ix] <- p.adjust(res$p_value[ix], method = "BH")
}

res <- res[order(res$platform, res$outcome, res$p_value), ]

write.table(
  res,
  file.path(out_dir, "GSE96058_adjusted_inflammation_estrogen_models.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

flow <- do.call(rbind, sample_flow)
write.table(
  flow,
  file.path(out_dir, "GSE96058_regression_sample_flow.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nER-positive regression sample flow:\n")
print(flow, row.names = FALSE)

cat("\nAdjusted regression results:\n")
print(res, row.names = FALSE, digits = 4)

cat("\nSaved results to:\n")
cat(file.path(out_dir, "GSE96058_adjusted_inflammation_estrogen_models.tsv"), "\n")
cat(file.path(out_dir, "GSE96058_regression_sample_flow.tsv"), "\n")
