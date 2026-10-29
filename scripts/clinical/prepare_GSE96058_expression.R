options(stringsAsFactors = FALSE)

expr_file <- "data/raw/GSE96058/GSE96058_transformed_expression.csv.gz"
meta_file <- "data/metadata/GSE96058/GSE96058_clinical_key_fields.tsv"
out_dir <- "data/processed/GSE96058"
report_dir <- "reports/GSE96058"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_file)) stop("Expression matrix not found: ", expr_file)
if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)

message("Reading clinical metadata...")
meta <- read.delim(
  meta_file, check.names = FALSE, na.strings = c("", "NA")
)

required_cols <- c("GSM", "title", "platform", "er.status",
                   "pam50.subtype", "age.at.diagnosis",
                   "tumor.size")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols)) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(meta$title)) {
  stop("Metadata titles are not unique; refusing ambiguous mapping.")
}

message("Reading transformed expression matrix...")
expr_df <- read.csv(
  gzfile(expr_file),
  header = TRUE,
  check.names = FALSE,
  row.names = NULL,
  na.strings = c("", "NA"),
  quote = "\""
)

if (ncol(expr_df) != 3410L) {
  warning("Expected 3410 columns including gene IDs; found ", ncol(expr_df))
}

gene_ids <- as.character(expr_df[[1]])
if (anyDuplicated(gene_ids)) {
  warning("Duplicate gene identifiers found; they will be audited, not silently collapsed.")
}

expr <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
colnames(expr) <- names(expr_df)[-1]
rownames(expr) <- gene_ids
rm(expr_df)
gc()

message("Checking sample mapping...")
if (anyDuplicated(colnames(expr))) stop("Duplicate expression sample labels.")
if (anyDuplicated(rownames(expr))) {
  message("Duplicate gene IDs: ", sum(duplicated(rownames(expr))))
}

idx <- match(colnames(expr), meta$title)
if (anyNA(idx)) {
  stop("Expression labels without metadata: ",
       paste(colnames(expr)[is.na(idx)][1:min(10, sum(is.na(idx)))],
             collapse = ", "))
}

meta_aligned <- meta[idx, , drop = FALSE]
stopifnot(identical(colnames(expr), meta_aligned$title))

# Replicates have titles ending in "repl"; preserve originals only.
is_replicate <- grepl("repl$", meta_aligned$title, ignore.case = TRUE)
message("Total expression samples: ", ncol(expr))
message("Replicate samples excluded from primary set: ", sum(is_replicate))
message("Original samples retained: ", sum(!is_replicate))

expr_primary <- expr[, !is_replicate, drop = FALSE]
meta_primary <- meta_aligned[!is_replicate, , drop = FALSE]

stopifnot(identical(colnames(expr_primary), meta_primary$title))

# Matrix-level quality checks.
n_nonfinite <- sum(!is.finite(expr_primary))
n_missing <- sum(is.na(expr_primary))
finite_vals <- expr_primary[is.finite(expr_primary)]

qc <- data.frame(
  metric = c(
    "genes_rows", "all_samples_columns", "primary_samples_columns",
    "replicate_samples_excluded", "duplicate_gene_ids",
    "missing_values", "nonfinite_values",
    "minimum_finite_value", "maximum_finite_value",
    "median_finite_value", "fraction_values_at_minus_3.321928"
  ),
  value = c(
    nrow(expr), ncol(expr), ncol(expr_primary),
    sum(is_replicate), sum(duplicated(rownames(expr))),
    n_missing, n_nonfinite,
    if (length(finite_vals)) min(finite_vals) else NA,
    if (length(finite_vals)) max(finite_vals) else NA,
    if (length(finite_vals)) median(finite_vals) else NA,
    if (length(finite_vals)) mean(abs(finite_vals + 3.32192809488736) < 1e-8) else NA
  )
)

write.table(
  qc, file.path(report_dir, "GSE96058_expression_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save aligned metadata and the primary expression matrix.
write.table(
  meta_primary,
  file.path(out_dir, "GSE96058_primary_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

saveRDS(
  expr_primary,
  file.path(out_dir, "GSE96058_transformed_expression_primary.rds"),
  compress = "gzip"
)

# Audit gene identifiers and platform/sample composition.
gene_audit <- data.frame(
  gene_id = rownames(expr),
  duplicated_id = duplicated(rownames(expr)) |
                  duplicated(rownames(expr), fromLast = TRUE)
)
write.table(
  gene_audit, file.path(report_dir, "GSE96058_gene_ID_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

platform_counts <- as.data.frame(table(
  platform = meta_primary$platform,
  useNA = "ifany"
))
names(platform_counts)[2] <- "n_samples"
write.table(
  platform_counts,
  file.path(report_dir, "GSE96058_primary_platform_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

er_counts <- as.data.frame(table(
  ER_status = meta_primary$er.status,
  useNA = "ifany"
))
names(er_counts)[2] <- "n_samples"
write.table(
  er_counts,
  file.path(report_dir, "GSE96058_primary_ER_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("\nQC summary:")
print(qc, row.names = FALSE)
message("\nPlatform counts:")
print(platform_counts, row.names = FALSE)
message("\nER status counts:")
print(er_counts, row.names = FALSE)
message("\nPreparation complete.")
message("Expression RDS: ", file.path(out_dir, "GSE96058_transformed_expression_primary.rds"))
message("Metadata TSV: ", file.path(out_dir, "GSE96058_primary_metadata.tsv"))
message("QC report: ", file.path(report_dir, "GSE96058_expression_QC.tsv"))
