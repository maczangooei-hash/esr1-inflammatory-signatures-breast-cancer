options(stringsAsFactors = FALSE)
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

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"
meta <- meta[match(colnames(expr), meta$title), , drop = FALSE]
scores <- scores[match(meta$title, scores$title), , drop = FALSE]

stopifnot(
  !anyNA(meta$title), !anyNA(scores$title),
  all(meta$title == colnames(expr)),
  all(scores$title == meta$title)
)

inflammation_sets <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
estrogen_sets <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)

shared_genes <- unique(unlist(gmt[inflammation_sets], use.names = FALSE))

# Recalculate overlap-removed estrogen scores.
for (set in estrogen_sets) {
  genes <- setdiff(intersect(gmt[[set]], rownames(expr)), shared_genes)
  new_col <- paste0(set, "_noInflammatoryOverlap")
  scores[[new_col]] <- NA_real_

  for (plat in unique(meta$platform)) {
    idx <- which(meta$platform == plat)
    gz <- t(scale(t(expr[genes, idx, drop = FALSE])))
    gz[!is.finite(gz)] <- NA_real_
    scores[[new_col]][idx] <- colMeans(gz, na.rm = TRUE)
  }
}

score_versions <- c("original", "overlap_removed")
diagnostics <- list()
interactions <- list()
k <- 0
j <- 0

for (plat in unique(meta$platform)) {
  pidx <- which(meta$platform == plat)

  for (estrogen_set in estrogen_sets) {
    revised_set <- paste0(estrogen_set, "_noInflammatoryOverlap")

    for (version in score_versions) {
      outcome_col <- if (version == "original") {
        estrogen_set
      } else {
        revised_set
      }

      for (inflammation_set in inflammation_sets) {
        dat <- data.frame(
          er_status = suppressWarnings(as.numeric(as.character(
            meta$er.status[pidx]
          ))),
          y = scores[[outcome_col]][pidx],
          inflammation = scores[[inflammation_set]][pidx],
          ESR1 = as.numeric(expr["ESR1", pidx]),
          age = suppressWarnings(as.numeric(as.character(
            meta$age.at.diagnosis[pidx]
          ))),
          tumor_size = suppressWarnings(as.numeric(as.character(
            meta$tumor.size[pidx]
          ))),
          PAM50 = meta$pam50.subtype[pidx]
        )

        keep <- with(
          dat,
          !is.na(er_status) & er_status == 1 &
          is.finite(y) & is.finite(inflammation) &
          is.finite(ESR1) & is.finite(age) &
          is.finite(tumor_size) & !is.na(PAM50)
        )
        keep[is.na(keep)] <- FALSE
        dat <- dat[keep, , drop = FALSE]
        dat$PAM50 <- droplevels(factor(dat$PAM50))

        if (nrow(dat) < 30 || nlevels(dat$PAM50) < 2) next

        # Standardize within this exact analysis cohort.
        dat$y_z <- as.numeric(scale(dat$y))
        dat$inflammation_z <- as.numeric(scale(dat$inflammation))
        dat$ESR1_z <- as.numeric(scale(dat$ESR1))
        dat$age_z <- as.numeric(scale(dat$age))
        dat$tumor_size_z <- as.numeric(scale(dat$tumor_size))

        fit <- lm(
          y_z ~ inflammation_z + ESR1_z + age_z +
            tumor_size_z + PAM50,
          data = dat
        )
        b_full <- unname(coef(fit)["inflammation_z"])

        # Remove the most influential 1% by Cook's distance.
        cooks <- cooks.distance(fit)
        n_remove <- max(1, ceiling(0.01 * nrow(dat)))
        remove_idx <- order(cooks, decreasing = TRUE)[seq_len(n_remove)]
        reduced <- dat[-remove_idx, , drop = FALSE]
        reduced$PAM50 <- droplevels(factor(reduced$PAM50))

        if (nlevels(reduced$PAM50) >= 2) {
          fit_reduced <- lm(
            y_z ~ inflammation_z + ESR1_z + age_z +
              tumor_size_z + PAM50,
            data = reduced
          )
          b_reduced <- unname(coef(fit_reduced)["inflammation_z"])
        } else {
          b_reduced <- NA_real_
        }

        k <- k + 1
        diagnostics[[k]] <- data.frame(
          platform = plat,
          estrogen_outcome = estrogen_set,
          score_version = version,
          inflammatory_predictor = inflammation_set,
          n = nrow(dat),
          beta_full = b_full,
          beta_after_removing_top_1pct_Cooks = b_reduced,
          absolute_beta_change = abs(b_reduced - b_full),
          direction_preserved = sign(b_full) == sign(b_reduced),
          samples_removed = n_remove,
          max_Cooks_distance = max(cooks, na.rm = TRUE)
        )

        # Test whether the inflammation slope differs by PAM50 subtype.
        fit_interaction <- lm(
          y_z ~ inflammation_z * PAM50 + ESR1_z +
            age_z + tumor_size_z,
          data = dat
        )
        fit_additive <- lm(
          y_z ~ inflammation_z + PAM50 + ESR1_z +
            age_z + tumor_size_z,
          data = dat
        )
        comparison <- anova(fit_additive, fit_interaction)

        j <- j + 1
        interactions[[j]] <- data.frame(
          platform = plat,
          estrogen_outcome = estrogen_set,
          score_version = version,
          inflammatory_predictor = inflammation_set,
          n = nrow(dat),
          interaction_df = comparison$Df[2],
          interaction_F = comparison$F[2],
          interaction_p = comparison$`Pr(>F)`[2]
        )
      }
    }
  }
}

diag_out <- do.call(rbind, diagnostics)
int_out <- do.call(rbind, interactions)

# Correct the interaction tests within each platform and score version.
int_out$interaction_FDR <- NA_real_
for (plat in unique(int_out$platform)) {
  for (version in unique(int_out$score_version)) {
    idx <- which(
      int_out$platform == plat &
      int_out$score_version == version
    )
    int_out$interaction_FDR[idx] <-
      p.adjust(int_out$interaction_p[idx], method = "BH")
  }
}

write.table(
  diag_out,
  "reports/GSE96058/GSE96058_influence_diagnostics.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  int_out,
  "reports/GSE96058/GSE96058_PAM50_interaction_diagnostics.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

cat("\n=== Influence diagnostics ===\n")
print(diag_out, row.names = FALSE)

cat("\n=== PAM50 interaction diagnostics ===\n")
print(int_out, row.names = FALSE)

cat("\nSaved reports:\n")
cat("reports/GSE96058/GSE96058_influence_diagnostics.tsv\n")
cat("reports/GSE96058/GSE96058_PAM50_interaction_diagnostics.tsv\n")
