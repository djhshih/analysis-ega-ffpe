#!/usr/bin/env Rscript
library(io)
library(precrec)
source("../common-ffpe-snvf/R/eval.R")

#######################################################

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
	ffpe_snvf_dir,
	ground_truth_dir,
	case_id_col = "tissue_type",
	sample_id_col = "sample_name",
	snvf_ext="tsv"
) {

	message(sprintf("Evaluating %s | Using FFPE VCFs from: %s, Ground Truths from: %s", model_name, ffpe_snvf_dir, ground_truth_dir))

	for (i in seq_len(nrow(ffpe_tumoral))){

		## Get FFPE sample metadata
		case_id <- ffpe_tumoral[i, case_id_col]
		sample_name <- ffpe_tumoral[i, sample_id_col]
		variant_caller <- "Mutect2"
		dataset <- basename(dirname(ffpe_snvf_dir))
		variant_set <- basename(ffpe_snvf_dir)
		outdir_root <- file.path(dataset, variant_set)

		message(sprintf("	%d. Processing Case: %s Sample: %s", i, case_id, sample_name))

		snvf_path <- file.path(ffpe_snvf_dir, model_name, sample_name, sprintf("%s.%s.%s", sample_name, model_name, snvf_ext))

		if (!file.exists(snvf_path)){
			message(sprintf("		Warning: %s SNVF does not exist at %s . Skipping", model_name, snvf_path))
			next
		}

		## Read in the ground truth
		ground_truth <- read.delim(file.path(ground_truth_dir, sample_name, sprintf("%s.ground-truth.tsv", sample_name)))
		ground_truth$truth <- as.logical(ground_truth$truth)

		# read model's score for each sample
		d <- read.delim(snvf_path)
		# Apply model specific processing
		d <- preprocess_filter(d, model_name)
		d <- merge(ground_truth, d, by = c("chrom", "pos", "ref", "alt"))
		
		## Check if true labels exist in the variant_score_truth table (d)
		## If not this means there's no overlap between FFPE and FF variants
		## Cases like these are skipped as evaluation is not supported by precrec
		if((nrow(d[d$truth, ]) == 0)){
			message(sprintf("		No true labels exist for %s | Skipping...", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}
		if((nrow(d[!d$truth, ]) == 0)){
			message(sprintf("		No false labels exist for %s | Skipping...", snvf_path))
			write_sample_eval(d, NULL, outdir_root, sample_name, model_name)
			next
		}

		# Evaluate the filter's performance
		res <- evaluate_filter(d, model_name, sample_name)

		# write results
		write_sample_eval(d, res, outdir_root, sample_name, model_name)
		
	}

	# Overall Evaluation
	## The scores annotated with ground truth is combined into a single dataframe
	message("	Performing Evaluation across all samples")

	all_score_truth <- do.call(
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
	overall_res <- evaluate_filter(all_score_truth, model_name, "all-samples")
	write_overall_eval(all_score_truth, overall_res, outdir_root, "all-samples", model_name)

	# Evaluate across colon samples
	colon_score_truth <- all_score_truth[grepl("Colon", all_score_truth$sample_name), ]
	colon_res <- evaluate_filter(colon_score_truth, model_name, "all-colon-samples")
	write_overall_eval(colon_score_truth, colon_res, outdir_root, "all-colon-samples", model_name)

	# Evaluate across liver samples
	liver_score_truth <- all_score_truth[grepl("Liver", all_score_truth$sample_name), ]
	liver_res <- evaluate_filter(liver_score_truth, model_name, "all-liver-samples")
	write_overall_eval(liver_score_truth, liver_res, outdir_root, "all-liver-samples", model_name)

}


########################################################

# Read Annotation Table linking FF sample to FFPE samples by case
annot_table <- read.delim("../annot/EGAD00001004066/sample_annotations.tsv")

# Stratify annotation table based on FFPE and FF Somatic Variants
annot_table <- annot_table[annot_table$sample_type == "Tumoral", ]

ffpe_tumoral <- annot_table[(annot_table$preservation == "FFPE"), ]
frozen_tumoral <- annot_table[(annot_table$preservation == "Frozen"), ]


# Evaluate the ffpe filter
model_name <- "ideafix"

## Evaluate the tumor-only variants with DP>=10, blacklist and MICR artifacts removed [MicroSEC Filter 1234]
evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered",
	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered"
)

# ## Evaluate the tumor-only variants with DP>=10, blacklist and MICR artifacts removed [MicroSEC Filter 1234]
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered-dp20",
# 	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered-dp20"
# )

# ## Evaluate the tumor-only variants with DP>=10, blacklist and MICR artifacts removed [MicroSEC Filter 1234]
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered-dp20-blacklist",
# 	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered-dp20-blacklist"
# )

# ## Evaluate the tumor-only variants with DP>=10, blacklist and MICR artifacts removed [MicroSEC Filter 1234]
# evaluate_sample_set(
# 	ffpe_tumoral = ffpe_tumoral,
# 	frozen_tumoral = frozen_tumoral,
# 	model_name = model_name,
# 	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered-dp20-blacklist-micr1234",
# 	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered-dp20-blacklist"
# )

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered-exome",
	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered-exome"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered-exome-blacklist",
	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered-exome-blacklist"
)

evaluate_sample_set(
	ffpe_tumoral = ffpe_tumoral,
	frozen_tumoral = frozen_tumoral,
	model_name = model_name,
	ffpe_snvf_dir = "../ffpe-snvf/EGAD00001004066/somatic_filtered-exome-blacklist-micr1234",
	ground_truth_dir = "../ground-truth/EGAD00001004066/somatic_filtered-exome-blacklist"
)

