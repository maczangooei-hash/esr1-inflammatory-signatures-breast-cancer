input <- "data/metadata/GSE96058/GSE96058_family.soft.gz"
outdir <- "data/metadata/GSE96058"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

con <- gzfile(input, open = "rt")
lines <- readLines(con, warn = FALSE)
close(con)

starts <- which(grepl("^\\^SAMPLE = ", lines))
if (!length(starts)) stop("No SAMPLE records found in SOFT file.")

ends <- c(starts[-1] - 1L, length(lines))
records <- vector("list", length(starts))

for (i in seq_along(starts)) {
  block <- lines[starts[i]:ends[i]]
  gsm <- sub("^\\^SAMPLE = ", "", block[1])
  title <- sub("^!Sample_title = ", "", 
               block[grepl("^!Sample_title = ", block)][1])
  source <- sub("^!Sample_source_name_ch1 = ", "",
                block[grepl("^!Sample_source_name_ch1 = ", block)][1])
  platform <- sub("^!Sample_platform_id = ", "",
                  block[grepl("^!Sample_platform_id = ", block)][1])

  chars <- block[grepl("^!Sample_characteristics_ch1 = ", block)]
  chars <- sub("^!Sample_characteristics_ch1 = ", "", chars)

  one <- data.frame(
    GSM = gsm,
    title = if (length(title)) title else NA_character_,
    source = if (length(source)) source else NA_character_,
    platform = if (length(platform)) platform else NA_character_,
    characteristic = NA_character_,
    value = NA_character_,
    stringsAsFactors = FALSE
  )

  if (length(chars)) {
    parts <- strsplit(chars, ": ", fixed = TRUE)
    keys <- vapply(parts, function(x) x[1], character(1))
    vals <- vapply(parts, function(x) {
      if (length(x) >= 2) paste(x[-1], collapse = ": ") else NA_character_
    }, character(1))

    records[[i]] <- data.frame(
      GSM = gsm, title = one$title, source = one$source,
      platform = one$platform,
      characteristic = keys, value = vals,
      stringsAsFactors = FALSE
    )
  } else {
    records[[i]] <- one
  }
}

long <- do.call(rbind, records)
write.table(long,
            file = file.path(outdir, "GSE96058_metadata_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Wide table for key variables; preserve raw values and labels.
wanted <- c(
  "scan-b external id", "er status", "pgr status", "her2 status",
  "pam50 subtype", "er prediction mgc", "er prediction sgc",
  "overall survival days", "overall survival event",
  "endocrine treated", "chemo treated", "age at diagnosis",
  "tumor size", "instrument model"
)

base <- unique(long[, c("GSM", "title", "source", "platform")])
for (key in wanted) {
  vals <- long[long$characteristic == key, c("GSM", "value")]
  vals <- vals[!duplicated(vals$GSM), ]
  names(vals)[2] <- make.names(key)
  base <- merge(base, vals, by = "GSM", all.x = TRUE, sort = FALSE)
}

write.table(base,
            file = file.path(outdir, "GSE96058_clinical_key_fields.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat("SOFT sample records:", length(starts), "\n")
cat("Unique GSM IDs:", length(unique(long$GSM)), "\n")
cat("Clinical annotation table:", nrow(base), "rows\n\n")

cat("Source types:\n")
print(sort(table(base$source, useNA = "ifany"), decreasing = TRUE))

cat("\nClinical ER status:\n")
print(table(base$er.status, useNA = "ifany"))

cat("\nPAM50 subtype:\n")
print(sort(table(base$pam50.subtype, useNA = "ifany"), decreasing = TRUE))

cat("\nPlatform:\n")
print(table(base$platform, useNA = "ifany"))

cat("\nSample titles containing 'repl':\n")
print(sum(grepl("repl", base$title, ignore.case = TRUE), na.rm = TRUE))

cat("\nDuplicate external IDs:\n")
external_col <- "scan.b.external.id"
if (external_col %in% names(base)) {
  ids <- base[[external_col]]
  print(sum(duplicated(ids[!is.na(ids) & ids != ""])))
} else {
  cat("External-ID column not found; inspect the output table.\n")
}

cat("\nSaved files:\n")
cat(file.path(outdir, "GSE96058_metadata_long.tsv"), "\n")
cat(file.path(outdir, "GSE96058_clinical_key_fields.tsv"), "\n")
