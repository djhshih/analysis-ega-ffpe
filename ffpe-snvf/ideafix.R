#!/usr/bin/env Rscript
library(ideafix)
library(io)
library(argparser)

# Parse Command-Line Arguments
p <- arg_parser("Perform FFPE SNVF with IdeaFix")
p <- add_argument(p, "--vcf", help = "path to the VCF file")
p <- add_argument(p, "--ref", help = "path to the reference genome")
p <- add_argument(p, "--outdir", help = "Output directory for the results", default = ".")
p <- add_argument(p, "--sample_name", help = "Name of Sample being filtered", default = NA)
argv <- parse_args(p)

# Assign Variables
vcf_path <- argv$vcf
ref_genome <- argv$ref
sample_name <- ifelse(is.na(argv$sample_name), unlist(strsplit(basename(vcf_path), "\\."))[1], argv$sample_name)

outdir <- file.path(argv$outdir, "ideafix", sample_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (endsWith(basename(vcf_path), ".vcf")){
	vcf_filename <- vcf_path
} else if (endsWith(basename(vcf_path), ".vcf.gz")) {
	message(cat(sprintf("\nVCF is compressed. Decompressing to %s \n", outdir)))
	# run bcftools and write vcf to outdir
	vcf_filename <- file.path(outdir, paste0(sample_name, ".vcf"))
	system2("bcftools", args = c("view", "-Oz", "-o", vcf_filename, vcf_path))
}

message(cat(sprintf("\nProcessing %s \n", sample_name)))
message(cat("\tGetting descriptors from VCF \n"))
descriptors <- get_descriptors(vcf_filename = vcf_filename, fasta_filename = ref_genome)

message(cat("\n\tRunning IdeaFix XGBoost Model \n"))
predictions_xgboost <- classify_variants(variant_descriptors = descriptors, algorithm = "XGBoost")
names(predictions_xgboost) <- tolower(names(predictions_xgboost))

# message(cat("\n\tRunning IdeaFix RF Model \n"))
# predictions_rf <- classify_variants(variant_descriptors = descriptors, algorithm = "RF")
# names(predictions_rf) <- tolower(names(predictions_rf))

message(cat(sprintf("\n\tWriting outputs to: %s \n", outdir)))
qwrite(descriptors, file.path(outdir, sprintf("%s.ideafix.descriptors.tsv", sample_name)))
qwrite(predictions_xgboost, file.path(outdir, sprintf("%s.ideafix.tsv", sample_name)))
# qwrite(predictions_rf, file.path(outdir, sprintf("%s.ideafix-rf.tsv", sample_name)))

message(cat("\n\tComplete"))
