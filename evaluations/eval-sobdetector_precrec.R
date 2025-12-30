#!/usr/bin/env Rscript
library(io)
library(precrec)
source("../common-ffpe-snvf/R/eval.R")

#######################################################

#' Load or Build Truth Set for a Case and Variant Caller
#'
#' This function loads a previously cached truth set for a given case and variant caller,
#' or constructs a new truth set if one does not already exist. The truth set is created
#' by taking the union of variants from matched fresh-frozen (FF) VCF files and cached by
#' assigning unique identifiers to each variant.
#'
#' @param case_id Character. The identifier for the case/tumor being analyzed.
#' @param variant_caller Character. The name of the variant caller used to generate
#'   the VCF files (e.g., "mutect2", "varscan2").
#' @param matched_ff_vcf_paths Character vector. File paths to the matched fresh-frozen
#'   VCF files to be used for creating the truth set.
#' @param outdir Character. The output directory where truth sets will be cached.
#'   Defaults to \code{main.outdir}.
#'
#' @return A data frame or tibble containing the truth set with unique variant identifiers.
#'   If cached, returns the previously saved truth set. Otherwise, returns a newly
#'   constructed truth set created from the union of variants in matched_ff_vcf_paths.
#'
#' @details
#' Truth sets are cached in the \code{truth_sets} subdirectory of \code{outdir}
#' with the naming convention: \code{<case_id>_<variant_caller>.tsv}
#'
#' @examples
#' \dontrun{
#' truth <- get_truth_set("CASE001", "mutect2", c("file1.vcf", "file2.vcf"))
get_truth_set <- function(case_id, matched_ff_vcf_paths, outdir = main.outdir) {
	truth_set_dir <- file.path(outdir, "truth_sets")
	dir.create(truth_set_dir, showWarnings = FALSE, recursive = TRUE)
	truth_set_path <- file.path(truth_set_dir, paste(basename(matched_ff_vcf_paths), collapse = "_"), sprintf("%s.tsv", case_id))

	if (file.exists(truth_set_path)){
		message(cat("\tGround truth set exists, reading from file"))
		truth <- qread(truth_set_path)
	} else {
		message(cat("\tCached ground truth set does not exists, generating ground truth set from:"))
		message(cat("\t", matched_ff_vcf_paths))
		truth <- snv_union(matched_ff_vcf_paths)
		truth <- add_id(truth)
		qwrite(truth, truth_set_path)
	}
	return(truth)
}

#' Evaluate a Set of Tumor Samples
#'
#' This function evaluates the performance of a variant calling model on a set of FFPE (Formalin-Fixed Paraffin-Embedded)
#' tumor samples by comparing model predictions against a ground truth set derived from matched frozen tumor samples.
#' It processes each sample individually, then performs an overall evaluation across all samples.
#'
#' @param ffpe_tumoral A data frame containing FFPE tumor sample metadata. Must include columns:
#'   - case_id: character, unique identifier for the case
#'   - tumor_bam_fid: character, file ID for the tumor BAM file
#'   - workflow_type: character, name of the variant caller used
#'
#' @param frozen_tumoral A data frame containing frozen tumor sample metadata. Must include columns:
#'   - case_id: character, unique identifier matching cases in ffpe_tumoral
#'   - vcf_fid: character, file ID for the VCF file
#'
#' @param model_name A character string specifying the name of the model being evaluated.
#'   Used for file naming and directory organization.
#'
#' @param outdir_root A character string specifying the root output directory where evaluation
#'   results will be written.
#'
#' @param ff_vcf_dir A character string specifying the directory containing frozen tumor VCF files.
#'
#' @param ffpe_snvf_dir A character string specifying the directory containing FFPE model score files.
#'   Expected structure: ffpe_snvf_dir/model_name/case_id/
#'
#' @return Invisibly returns NULL. The function writes evaluation results to disk:
#'   - Per-sample evaluation results
#'   - Overall evaluation across all samples
#'
#' @details
#' For each sample, the function:
#' 1. Retrieves FFPE and matched frozen tumor metadata
#' 2. Reads model predictions and ground truth variants
#' 3. Preprocesses the data using model-specific logic
#' 4. Evaluates model performance if both true and false labels exist
#' 5. Writes results to disk
#'
#' Samples are skipped if no overlap exists between FFPE and frozen variants or if
#' either true or false truth labels are absent.
evaluate_sample_set <- function(
	ffpe_tumoral,
	frozen_tumoral,
	model_name,
	ff_vcf_dir,
	ffpe_snvf_dir,
	case_id_col = "tissue_type",
	sample_id_col = "sample_name"
) {

	for (i in seq_len(nrow(ffpe_tumoral))){

		## Get FFPE sample metadata
		case_id <- ffpe_tumoral[i, case_id_col]
		sample_name <- ffpe_tumoral[i, sample_id_col]
		variant_caller <- "Mutect2"
		dataset <- basename(dirname(ffpe_snvf_dir))
		variant_set <- basename(ffpe_snvf_dir)
		outdir_root <- file.path(dataset, variant_set)

		## Getting matched FF metadata by matching case ID and workflow type i.e variant caller
		matched_ff_annot <- frozen_tumoral[(frozen_tumoral[[case_id_col]] == case_id), ]
		matched_ff_sample_names <- matched_ff_annot[[sample_id_col]]
		matched_ff_vcf_paths <- file.path(ff_vcf_dir, matched_ff_sample_names, sprintf("%s.vcf.gz", matched_ff_sample_names))

		snvf_path <- file.path(ffpe_snvf_dir, model_name, sample_name, sprintf("%s.%s.snv", sample_name, model_name))

		if (!file.exists(snvf_path)){
			message(sprintf("	Warning: %s SNVF does not exist at %s . Skipping", model_name, snvf_path))
			next
		}

		# read model's score for each sample
		d <- read.delim(snvf_path)
		# Get ground truth for that sample
		truth <- get_truth_set(case_id, matched_ff_vcf_paths, outdir = outdir_root)
		# Apply model specific processing
		d <- preprocess_sobdetector(d, truth)
		
		## Check if true labels exist in the variant_score_truth table (d)
		## If not this means there's no overlap between FFPE and FF variants
		## Cases like these are skipped as evaluation is not supported by precrec
		if((nrow(d[d$truth, ]) == 0)){
			message(sprintf("	no true labels exist for %s", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}
		if((nrow(d[!d$truth, ]) == 0)){
			message(sprintf("	no false labels exist for %s", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}

		message(sprintf("	%s", snvf_path))

		# Evaluate the filter's performance
		sobdetector_res <- evaluate_filter(d, model_name, sample_name)

		# write results
		write_sample_eval(d, sobdetector_res, outdir_root, sample_name, model_name)
		
	}

	# Overall Evaluation
	## The scores annotated with ground truth is combined into a single dataframe
	message("	performing Evaluation across all samples")

	sobdetector_all_score_truth <- do.call(
		rbind,
		lapply(ffpe_tumoral[[sample_id_col]], function(sample_name) {
			path <- file.path(outdir_root, "model-scores_truths", sample_name, sprintf("%s_%s-scores_truths.tsv", sample_name, model_name))
			if (!file.exists(path)){
				message(sprintf("	Warning: %s was not was not found. SKIPPING", path))
			} else {
				d <- read.delim(path)
				d$sample_name <- sample_name
				d
			}
		})
	)

	# Evaluate across all samples
	sobdetector_overall_res <- evaluate_filter(sobdetector_all_score_truth, model_name, "all-samples")
	write_overall_eval(sobdetector_all_score_truth, sobdetector_overall_res, outdir_root, "all-samples", model_name)

}


########################################################

# Read Annotation Table linking FF sample to FFPE samples by case
annot_table <- read.delim("../annot/EGAD00001004066/sample_annotations.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
annot_table <- annot_table[annot_table$sample_type == "Tumoral", ]

ffpe_tumoral <- annot_table[(annot_table$preservation == "FFPE"), ]
frozen_tumoral <- annot_table[(annot_table$preservation == "Frozen"), ]


# Evaluate the ffpe filter
message("Evaluating sobdetector:")
model_name <- "sobdetector"

## Evaluate the tumor-only dataset
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/EGAD00001004066/somatic_vcf",
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_vcf"
)

## Evaluate the tumor-only dataset with DP>=10
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/EGAD00001004066/somatic_vcf-dp10",
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_vcf-dp10"
)

## Evaluate the tumor-only dataset with DP>=10 and MICR artifacts removed [MicroSEC Filter 1234]
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ff_vcf_dir = "../vcf/EGAD00001004066/somatic_vcf-dp10",
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_vcf-dp10_msec-1234-excluded"
)