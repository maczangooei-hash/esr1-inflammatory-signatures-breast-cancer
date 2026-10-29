suppressPackageStartupMessages(library(limma))

base <- "data/processed/GSE126003"
out <- "reports/GSE126003_t47d_siRNA_analysis"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

matrix_file <- file.path(
  base, "GSE126003_t47d_siRNA_RPKM_deduplicated.tsv.gz"
)
mapping_file <- "reports/dataset_audit/ensembl115_gene_mapping.tsv"

stopifnot(file.exists(matrix_file), file.exists(mapping_file))

# Read expression matrix and gene mapping
x <- read.delim(matrix_file, check.names=FALSE,
                stringsAsFactors=FALSE)
ids <- sub("\\..*$", "", x[[1]])
expr <- as.matrix(x[, -1, drop=FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- ids

mapping <- read.delim(mapping_file, stringsAsFactors=FALSE)
mapping <- mapping[!duplicated(mapping$Ensembl_ID_without_version), ]

idx <- match(rownames(expr), mapping$Ensembl_ID_without_version)
symbols <- mapping$Gene_symbol[idx]

# Parse sample names
samples <- colnames(expr)
sirna <- sub(".*_t47d_(si[^_]+)_.*", "\\1", samples)
treatment <- sub(".*_t47d_si[^_]+_(il6|pbs)_.*", "\\1", samples)
treatment <- toupper(treatment)

if (any(!sirna %in% c("sictl", "sistat3", "sier", "sifoxa1")) ||
    any(!treatment %in% c("IL6", "PBS"))) {
  stop("Could not parse all sample names. Inspect colnames(expr).")
}

meta <- data.frame(
  sample=samples,
  siRNA=sirna,
  treatment=treatment,
  stringsAsFactors=FALSE
)
meta$siRNA <- factor(meta$siRNA,
                     levels=c("sictl", "sistat3", "sier", "sifoxa1"))
meta$treatment <- factor(meta$treatment, levels=c("PBS", "IL6"))

write.table(meta, file.path(out, "sample_metadata.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# Verify balanced replicate groups
group_counts <- as.data.frame.matrix(
  table(meta$siRNA, meta$treatment)
)
write.table(group_counts, file.path(out, "group_counts.tsv"),
            sep="\t", quote=FALSE, col.names=NA)
if (any(table(meta$siRNA, meta$treatment) == 0))
  stop("At least one siRNA-by-treatment group is empty.")

# Processed RPKM: use log2(RPKM + 1), not a raw-count model
y <- log2(expr + 1)

# Remove features with no variation across samples
keep <- apply(y, 1, function(z) all(is.finite(z)) && var(z) > 0)
y <- y[keep, , drop=FALSE]
symbols <- symbols[keep]
ids <- rownames(y)

design <- model.matrix(~ siRNA * treatment, data=meta)
colnames(design) <- gsub(":", ".", colnames(design), fixed=TRUE)
fit <- lmFit(y, design)

# IL6 versus PBS within each siRNA group
contrast_matrix <- makeContrasts(
  IL6_vs_PBS_sictl =
    treatmentIL6,
  IL6_vs_PBS_sistat3 =
    treatmentIL6 + siRNAsistat3.treatmentIL6,
  IL6_vs_PBS_sier =
    treatmentIL6 + siRNAsier.treatmentIL6,
  IL6_vs_PBS_sifoxa1 =
    treatmentIL6 + siRNAsifoxa1.treatmentIL6,
  Interaction_sistat3_vs_sictl =
    siRNAsistat3.treatmentIL6,
  Interaction_sier_vs_sictl =
    siRNAsier.treatmentIL6,
  Interaction_sifoxa1_vs_sictl =
    siRNAsifoxa1.treatmentIL6,
  levels=design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust=TRUE)

# Save complete results for each contrast
for (nm in colnames(contrast_matrix)) {
  tab <- topTable(fit2, coef=nm, number=Inf, sort.by="P")
  tab$Ensembl_ID <- rownames(tab)
  tab$Gene_symbol <- symbols[match(rownames(tab), ids)]
  tab <- tab[, c("Ensembl_ID", "Gene_symbol",
                 setdiff(names(tab), c("Ensembl_ID", "Gene_symbol")))]
  write.table(tab, file.path(out, paste0(nm, ".tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
}

# Predefined biological target panel
targets <- c(
  "ESR1", "ESR2", "TFF1", "PGR", "CCND1", "IL6", "STAT3",
  "NFKBIA", "CXCL8", "JUN", "FOS", "IL1R1", "TNFAIP3",
  "PTGS2", "CCL2", "CXCL2", "ICAM1", "BIRC3", "RELA"
)

for (nm in colnames(contrast_matrix)) {
  tab <- topTable(fit2, coef=nm, number=Inf, sort.by="P")
  tab$Ensembl_ID <- rownames(tab)
  tab$Gene_symbol <- symbols[match(rownames(tab), ids)]
  target_tab <- tab[tab$Gene_symbol %in% targets, ]
  target_tab <- target_tab[, c(
    "Ensembl_ID", "Gene_symbol",
    setdiff(names(target_tab), c("Ensembl_ID", "Gene_symbol"))
  )]
  write.table(target_tab,
              file.path(out, paste0(nm, "_target_panel.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
}

# Record model information and feature handling
writeLines(c(
  "Dataset: GSE126003, T47D siRNA experiment",
  paste("Input:", matrix_file),
  paste("Samples:", ncol(y)),
  paste("Features tested:", nrow(y)),
  "Expression transformation: log2(RPKM + 1)",
  "Method: limma linear model and empirical Bayes moderation",
  "Multiple testing: Benjamini-Hochberg FDR in topTable",
  "Main contrasts: IL6 versus PBS within each siRNA background",
  "Interaction contrasts: change in IL6 response versus control siRNA",
  "Duplicate handling: one representative per stable Ensembl ID;",
  "the prior audit preferred rows without _PAR_Y; no summation.",
  "Gene symbols joined using Ensembl IDs without version suffixes.",
  "Interpretation is exploratory; processed RPKM is not raw counts.",
  "Interaction contrasts, not differences in individual p-values,",
  "test whether IL6 responses differ between siRNA backgrounds."
), file.path(out, "analysis_notes.txt"))

cat("Analysis complete.\n")
cat("Samples:", ncol(y), "\n")
cat("Features tested:", nrow(y), "\n")
cat("Results directory:", out, "\n")
