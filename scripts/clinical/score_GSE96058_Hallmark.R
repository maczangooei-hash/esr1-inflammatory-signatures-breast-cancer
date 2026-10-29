options(stringsAsFactors = FALSE)

# Make the existing user-installed R packages visible.
expr_file <- "data/processed/GSE96058/GSE96058_transformed_expression_primary.rds"
meta_file <- "data/processed/GSE96058/GSE96058_primary_metadata.tsv"
gmt_file <- "data/references/h.all.v2025.1.Hs.symbols.gmt"

out_dir <- "data/processed/GSE96058"
report_dir <- "reports/GSE96058"

stopifnot(file.exists(expr_file), file.exists(meta_file), file.exists(gmt_file))

expr <- readRDS(expr_file)
meta <- read.delim(meta_file, check.names = FALSE, na.strings = c("", "NA"))

stopifnot(identical(colnames(expr), meta$title))

# Read the Hallmark GMT file.
gmt_lines <- readLines(gmt_file, warn = FALSE)
gmt <- list()

for (line in gmt_lines) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) >= 3) {
    gmt[[fields[1]]] <- unique(fields[-c(1, 2)])
  }
}

set_names <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)

missing_sets <- setdiff(set_names, names(gmt))
if (length(missing_sets)) {
  stop("Missing Hallmark sets: ", paste(missing_sets, collapse = ", "))
}

# Standardize each gene across samples WITHIN each platform, then average
# the standardized values for the genes in each set.
# This produces platform-relative scores, not absolute cross-platform values.
score_matrix <- matrix(
  NA_real_,
  nrow = length(set_names),
  ncol = ncol(expr),
  dimnames = list(set_names, colnames(expr))
)

coverage <- data.frame(
  gene_set = set_names,
  set_size = integer(length(set_names)),
  genes_found = integer(length(set_names)),
  coverage_fraction = numeric(length(set_names)),
  stringsAsFactors = FALSE
)

for (i in seq_along(set_names)) {
  genes <- gmt[[set_names[i]]]
  found <- intersect(genes, rownames(expr))
  coverage$set_size[i] <- length(genes)
  coverage$genes_found[i] <- length(found)
  coverage$coverage_fraction[i] <- length(found) / length(genes)
  if (length(found) < 10) {
    stop("Too few genes found for ", set_names[i], ": ", length(found))
  }
}

# Score each platform separately to avoid mixing platform distributions.
platforms <- unique(meta$platform)

for (plat in platforms) {
  sample_idx <- which(meta$platform == plat)
  x <- expr[, sample_idx, drop = FALSE]

  for (set_name in set_names) {
    genes <- intersect(gmt[[set_name]], rownames(x))
    gene_expr <- x[genes, , drop = FALSE]

    # Z-score each gene across samples in this platform.
    gene_z <- t(scale(t(gene_expr), center = TRUE, scale = TRUE))
    gene_z[!is.finite(gene_z)] <- NA_real_

    score_matrix[set_name, sample_idx] <- colMeans(gene_z, na.rm = TRUE)
  }
}

# Include clinical metadata and ESR1 expression as a separate covariate.
if (!"ESR1" %in% rownames(expr)) stop("ESR1 not found in expression matrix.")

scores <- data.frame(
  meta,
  ESR1_expression = as.numeric(expr["ESR1", meta$title]),
  t(score_matrix[, meta$title, drop = FALSE]),
  check.names = FALSE
)

write.table(
  scores,
  file.path(out_dir, "GSE96058_Hallmark_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

write.table(
  coverage,
  file.path(report_dir, "GSE96058_Hallmark_gene_set_coverage.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Report the distribution of scores separately by platform.
summary_rows <- list()
for (plat in platforms) {
  ix <- which(scores$platform == plat)
  for (set_name in set_names) {
    v <- scores[[set_name]][ix]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      platform = plat,
      gene_set = set_name,
      n = sum(is.finite(v)),
      mean = mean(v, na.rm = TRUE),
      sd = sd(v, na.rm = TRUE),
      min = min(v, na.rm = TRUE),
      median = median(v, na.rm = TRUE),
      max = max(v, na.rm = TRUE)
    )
  }
}
score_summary <- do.call(rbind, summary_rows)

write.table(
  score_summary,
  file.path(report_dir, "GSE96058_Hallmark_score_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nGene-set coverage:\n")
print(coverage, row.names = FALSE)

cat("\nScore summaries by platform:\n")
print(score_summary, row.names = FALSE)

cat("\nOutput files:\n")
cat(file.path(out_dir, "GSE96058_Hallmark_scores.tsv"), "\n")
cat(file.path(report_dir, "GSE96058_Hallmark_gene_set_coverage.tsv"), "\n")
cat(file.path(report_dir, "GSE96058_Hallmark_score_summary.tsv"), "\n")
