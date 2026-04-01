#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("GMMAT", quietly = TRUE)) {
    stop("Package 'GMMAT' is required. Install it in R before running this script.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript run_gmmat_gwas.R <phenotype_tsv> <vcf_or_vcfgz> <grm_tsv> [output_prefix]"
    )
  )
}

phenotype_path <- normalizePath(args[[1]], mustWork = TRUE)
vcf_path <- normalizePath(args[[2]], mustWork = TRUE)
grm_path <- normalizePath(args[[3]], mustWork = TRUE)
output_prefix <- if (length(args) >= 4) args[[4]] else "simulated_data/gmmat"

parse_gt_to_dosage <- function(gt) {
  gt_main <- sub(":.*$", "", gt)
  if (gt_main %in% c(".", "./.", ".|.")) {
    return(NA_real_)
  }
  alleles <- unlist(strsplit(gt_main, "[/|]"))
  if (length(alleles) != 2 || any(alleles == ".")) {
    return(NA_real_)
  }
  sum(as.numeric(alleles))
}

read_vcf_dosage_matrix <- function(vcf_file) {
  open_fun <- if (grepl("\\.gz$", vcf_file, ignore.case = TRUE)) gzfile else file
  con <- open_fun(vcf_file, open = "rt")
  on.exit(close(con), add = TRUE)

  sample_ids <- NULL
  variant_rows <- list()
  line_index <- 0L

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    if (startsWith(line, "##")) {
      next
    }
    if (startsWith(line, "#CHROM")) {
      header_fields <- strsplit(line, "\t", fixed = FALSE)[[1]]
      sample_ids <- header_fields[-(1:9)]
      next
    }

    fields <- strsplit(line, "\t", fixed = FALSE)[[1]]
    line_index <- line_index + 1L
    dosages <- vapply(fields[-(1:9)], parse_gt_to_dosage, numeric(1))
    variant_rows[[line_index]] <- c(fields[1], fields[3], fields[2], dosages)
  }

  if (is.null(sample_ids)) {
    stop("VCF header not found in ", vcf_file)
  }

  dosage_df <- as.data.frame(do.call(rbind, variant_rows), stringsAsFactors = FALSE)
  colnames(dosage_df) <- c("CHR", "SNP", "POS", sample_ids)
  dosage_df$POS <- as.integer(dosage_df$POS)
  for (sample_id in sample_ids) {
    dosage_df[[sample_id]] <- as.numeric(dosage_df[[sample_id]])
  }

  list(sample_ids = sample_ids, dosage_df = dosage_df)
}

normalize_phenotype <- function(phenotype_df) {
  real_cols <- c("queen_id", "indv", "offspring_index", "phenotype")
  toy_cols <- c("queen_id", "offspring_id", "offspring_index", "recombination_breaks")

  if (all(real_cols %in% colnames(phenotype_df))) {
    out <- phenotype_df[, real_cols, drop = FALSE]
    colnames(out) <- c("queen_id", "offspring_id", "offspring_index", "recombination_breaks")
    return(out)
  }

  if (all(toy_cols %in% colnames(phenotype_df))) {
    return(phenotype_df[, toy_cols, drop = FALSE])
  }

  stop(
    paste(
      "Phenotype file is missing a supported column set.",
      "Expected either:",
      paste(real_cols, collapse = ", "),
      "or",
      paste(toy_cols, collapse = ", ")
    )
  )
}

phenotype_raw <- read.table(
  phenotype_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
phenotype <- normalize_phenotype(phenotype_raw)

vcf_data <- read_vcf_dosage_matrix(vcf_path)
sample_ids <- vcf_data$sample_ids
dosage_df <- vcf_data$dosage_df

grm_df <- read.table(grm_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
grm <- as.matrix(grm_df)
if (!all(sample_ids %in% rownames(grm))) {
  stop("Some VCF sample IDs are missing from the GRM.")
}
grm <- grm[sample_ids, sample_ids, drop = FALSE]

if (!all(unique(phenotype$queen_id) %in% sample_ids)) {
  stop("Phenotype queens do not match VCF sample IDs.")
}
phenotype$queen_id <- as.character(phenotype$queen_id)
phenotype$offspring_id <- as.character(phenotype$offspring_id)
phenotype$offspring_index <- as.integer(phenotype$offspring_index)
phenotype$recombination_breaks <- as.numeric(phenotype$recombination_breaks)
phenotype <- phenotype[order(phenotype$offspring_index, phenotype$queen_id), ]

null_model <- GMMAT::glmmkin(
  fixed = recombination_breaks ~ 1,
  data = phenotype,
  kins = grm,
  id = "queen_id",
  family = gaussian(link = "identity")
)

genotype_text_path <- paste0(output_prefix, ".geno.tsv")
null_model_path <- paste0(output_prefix, ".null_model.rds")
assoc_path <- paste0(output_prefix, ".assoc.tsv")

write.table(
  dosage_df,
  file = genotype_text_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(null_model, file = null_model_path)

select_idx <- match(sample_ids, unique(null_model$id_include))
if (any(is.na(select_idx))) {
  stop("Failed to match genotype sample order to null-model IDs.")
}

assoc_result <- GMMAT::glmm.score(
  obj = null_model,
  infile = genotype_text_path,
  outfile = assoc_path,
  infile.sep = "\t",
  infile.ncol.skip = 3,
  infile.ncol.print = 1:3,
  infile.header.print = c("CHR", "SNP", "POS"),
  select = select_idx,
  is.dosage = FALSE
)

if (file.exists(assoc_path)) {
  assoc_table <- read.table(
    assoc_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  assoc_table <- assoc_table[order(assoc_table$PVAL), ]
  top_table <- utils::head(assoc_table, 10)
  print(top_table)
} else {
  print(assoc_result)
}

cat("Phenotype file:", phenotype_path, "\n")
cat("VCF file:", vcf_path, "\n")
cat("GRM file:", grm_path, "\n")
cat("Genotype matrix:", genotype_text_path, "\n")
cat("Null model RDS:", null_model_path, "\n")
cat("Association results:", assoc_path, "\n")
