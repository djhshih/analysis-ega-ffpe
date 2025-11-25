library(io)
library(precrec)
library(argparser)

# TODO move common functions into module
# TODO refactoring in progress

source("../common-ffpe-snvf/R/eval.R")

# Setup Directories

## specific for each dataset

dataset_id <- "EGAD00001004066"

# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/somatic_vcf", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/somatic_vcf", dataset_id)
# Output directory
main.outdir <- sprintf("%s/somatic_vcf", dataset_id)
score_truth_outdir <- sprintf("%s/model-scores_truths/", main.outdir)
eval_outdir <- sprintf("%s/roc-prc-auc/precrec", main.outdir)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample_annotations.tsv", dataset_id))

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE" & lookup_table$sample_type == "Tumoral"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen" & lookup_table$sample_type == "Tumoral"), ]


#######################################

# Evaluate ideafix
## Per Sample
## Each sample is evaluated first due to the necessity of independently annotating the scores with ground truth

## Workflow:
# The SNVs with the score annotations are first read
# Then the ground truth is constructed by reading SNVs from the match FF sample. 
# An SNV union is made from multiple FF samples if they exist
# The truths are annotated via intersection of the FF SNVs and the FFPE SNVs
# Evaluation is then carried out using precrec
# The per sample evaluations result (AUC, ROC coordinates and PRC coordinates) and the SNV score truth annotation table is then saved
message("Evaluating ideafix:")
model_name <- "ideafix"
for (index in seq_len(nrow(ffpe_tumoral))){
	
	metadata <- set_up(ffpe_tumoral, index)
	message(sprintf("	%s", metadata$sample_name))
	
	ideafix_processed <- process_sample(
		read_snv,
		construct_ground_truth,
		preprocess_ideafix,
		evaluate_filter,
		sample_name = metadata$sample_name,
		tissue = metadata$tissue,
		filter_name = model_name,
		snvf_dir = ffpe_snvf.dir,
		gt_vcf_dir = vcf.dir,
		gt_annot_d = frozen_tumoral
	)

	write_sample_eval(ideafix_processed$d, ideafix_processed$res, main.outdir, metadata$sample_name, model_name)

}


# Overall Evaluation
## The scores annotated with ground truth is combined into a single dataframe
message("	performing Evaluation across all samples")
ideafix_all_score_truth <- do.call(
	rbind,
	lapply(seq_len(nrow(ffpe_tumoral)), function(i) {
		meta <- set_up(ffpe_tumoral, i)
		path <- file.path(score_truth_outdir, meta$sample_name, sprintf("%s_%s-scores_truths.tsv", meta$sample_name, model_name))
		d <- read.delim(path)
		d$sample_name <- meta$sample_name
		d
	})
)


# Evaluate across all samples
ideafix_overall_res <- evaluate_filter(ideafix_all_score_truth, model_name, "all-samples")
write_overall_eval(ideafix_all_score_truth, ideafix_overall_res, score_truth_outdir, eval_outdir, "all-samples", model_name)


# Evaluate across colon samples
message("	performing Evaluation across all colon samples")
ideafix_colon_score_truth <- ideafix_all_score_truth[grepl("Colon", ideafix_all_score_truth$sample_name), ]
ideafix_colon_res <- evaluate_filter(ideafix_colon_score_truth, model_name, "all-colon-samples")
write_overall_eval(ideafix_colon_score_truth, ideafix_colon_res, score_truth_outdir, eval_outdir, "all-colon-samples", model_name)


## Evaluate across liver samples
message("	performing Evaluation across all liver samples")
ideafix_liver_score_truth <- ideafix_all_score_truth[grepl("Liver", ideafix_all_score_truth$sample_name), ]
ideafix_liver_res <- evaluate_filter(ideafix_liver_score_truth, model_name, "all-liver-samples")
write_overall_eval(ideafix_liver_score_truth, ideafix_liver_res, score_truth_outdir, eval_outdir, "all-liver-samples", model_name)

message("Done.")
