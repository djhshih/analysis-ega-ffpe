library(io)
library(precrec)
library(argparser)

# TODO split up filteres into modules
# TODO refactoring in progress

source("../R/eval.R")

# Setup Directories

## specific for each dataset

dataset_id <- "EGAD00001004066"

# Directory for FFPE SNVF inputs
ffpe_snvf.dir <- sprintf("../ffpe-snvf/%s/somatic_vcf", dataset_id)
# Directory for the Somatic VCFs
vcf.dir <- sprintf("../vcf/%s/somatic_vcf", dataset_id)
# Output directory
main.outdir <- sprintf("%s/somatic_vcf", dataset_id)


# Read Annotation Table
lookup_table <- read.delim(sprintf("../annot/%s/sample_annotations.tsv", dataset_id))
# create a sample name column which are the file names for the samples
lookup_table$sample_name <- paste0(gsub(" ", "-", lookup_table$title), "_", lookup_table$sample_alias)

# Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE" & lookup_table$sample_type == "Tumoral"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen" & lookup_table$sample_type == "Tumoral"), ]

#################################


read_snv <- function(sample_name, filter_name) {
	path <- file.path(ffpe_snvf.dir, filter_name, sample_name, sprintf("%s.mobsnvf.snv", sample_name))
	read.delim(path)
}

# Construct the Ground Truth SNVs
## To do this we make a union set for all the FF samples matched to this FFPE sample
## This dataset has matched FFPE and FF from the same patient.
## All samples within the same tissue type are matched
## Hence, we select the frozen variants from the same tissue type as the FFPE sample
## @params d		data.frame containing sample annotation for frozen samples
## @params tissue_type		string describing the type of tissue		
construct_ground_truth <- function(d, tissue){

	d <- d[d$tissue_type == tissue, ]

	## Obtain the sample names for the frozen samples within the same tissue
	truth_samples <- d$sample_name
	truth_sample_paths <- file.path(vcf.dir, truth_samples, sprintf("%s.vcf.gz", truth_samples))
	truths <- snv_union(truth_sample_paths)
	truths <- add_id(truths)
	truths
}

# @param d  data.frame of variant annotation by mobsnvf
# @param truths  data.frame of ground-truth variants
preprocess_mobsnvf <- function(d, truths) {
	# mobsnvf sets FOBP to NA for variants that are not C>T
	d <- d[!is.na(d$FOBP), ];
	# lower score signifies real mutation:  
	# hence, we flip scores to make higher score signify a real mutation
	d$score <- -d$FOBP;
	d <- add_id(d);
	d <- annotate_truth(d, truths)
	d
}

# @param d  data.frame of variant annotation by vafsnvf
# @param truths  data.frame of ground-truth variants
preprocess_vafsnvf <- function(d, truths) {
	d <- d[!is.na(d$VAFF), ]
	# vafsnvf sets VAFF to NA for variants that are not C>T
	d$score <- d$VAFF;
	d <- add_id(d);
	d <- annotate_truth(d, truths)
	d
}

# @param d  data.frame of variant annotation by sobdetector
# @param truths  data.frame of ground-truth variants
preprocess_sobdetect <- function(d, truths) {
	# SOBDetector output column explanations:
	# 		artiStatus: Binary classification made by SOBDetector. Values are "snv" or "artifact"
	# 		SOB: This is the strand oreintation bias score column which ranges from 0 and 1. Exception values: "." or NaN. 
	# variants ignored by SOBdetector have score of "."
	d <- d[!(d$SOB == "."), ]
	# now, it's safe to convert to numeric
	d$SOB <- as.numeric(d$SOB)
	# SOBdetector score = 0 indicates that it is not artifact
	# variants classified by SOBdetect as "snv" have score of NaN
	d$SOB <- ifelse(is.nan(d$SOB), 0, d$SOB)
	### Precaution to remove the any NA Scores if present
	d <- d[!is.na(d$SOB), ]
	# in this set the TRUE label indicates real mutation and FALSE indicates artifacts
	# Hence, scores needs to be adjusted so that higher score represents a real mutation.
	d$score <- -d$SOB;
	# Keep only C>T variants
	d <- ct_filter(d)
	d <- add_id(d);
	d <- annotate_truth(d, truths)
	d
}

# @param d  data.frame of variant annotation by GATK Orientation Bias Mixture Model
# @param truths  data.frame of ground-truth variants
preprocess_gatk_obmm <- function(d, truths){
	### GATK Orientation Bias mixture model makes binary classification. This is casted into scores 0 and 1
	d$score <- ifelse(grepl("orientation", d$filter), 0, 1)
	d <- d[!is.na(d$score), ]
	# Keep only C>T variants
	d <- ct_filter(d)
	d <- add_id(d)
	d <- annotate_truth(d, truths)
	d
}


# @param d       data.frame of variant annotation
evaluate_filter <- function(d, name) {;
	eval_obj <- with(d, evalmod(scores = score, labels = truth, modnames = name));

	eval_auc <- auc(eval_obj)
	eval_auc$dsids <- NULL

	# Collect ROC PRC coordinates
	## Columns -> (x, y, modname, type):
	## x and y are coordinates for making ROC and PRC plots,
	## plot type values are ROC or PRC
	roc_prc <- as.data.frame(eval_obj)
	roc_prc$dsid <- NULL
	# Stratify ROC PRC
	roc <- roc_prc[roc_prc$type == "ROC", ]
	prc <- roc_prc[roc_prc$type == "PRC", ]

	list(
		eval = eval_obj,
		auc = eval_auc,
		roc = roc,
		prc = prc
	)
}

# Remark: evalmod will be replaced with another function in the future!

# higher-order function
# custom functions specific for each filter: read_f, preprocess_f, evaluate_f
# @return list consisting of data.frame "d" and evaluation result object "res".
process_sample <- function(sample_name, read_f, preprocess_f, evaluate_f) {
	d <- read_f(sample_name);
	d <- preprocess_f(d);
	res <- evaluate_f(d);
	
	list(
		d = d,
		res = res
	)
}

# @params d		data.frame containing sample annotation
# @params index		the index of the sample to process
set_up <- function(d, index) {
	list(
		sample_name = d[index, "sample_name"],
		tissue = d[index, "tissue_type"]
	)
}


## Create a write res function

# process across all samples
# THIS IS IMPORTANT PART
# process mobsnvf filter
# process vafsnvf filter

# process each sample individually


# process each sample individually

# Evaluate mobsnvf
message("Evaluating Mobsnvf:")
for (index in seq_len(nrow(ffpe_tumoral))){
	
	sample_name <- ffpe_tumoral[index, "sample_name"]
	tissue <- ffpe_tumoral[index, "tissue_type"]
	message(sprintf("	%s", sample_name))
	
	## Obtain ground truth for this sample
	truths <- construct_ground_truth(frozen_tumoral, tissue)
	
	mobsnvf_d <- read.delim(file.path(ffpe_snvf.dir, "mobsnvf", sample_name, sprintf("%s.mobsnvf.snv", sample_name)))
	mobsnvf_d <- preprocess_mobsnvf(mobsnvf_d, truths);
	mobsnvf_res <- evaluate_filter(mobsnvf_d, "mobsnvf");


	# Saving the variant set with scores and ground truth labels for each sample
	out_dir <- file.path(main.outdir, "model-scores_truths", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(mobsnvf_d, file.path(out_dir, sprintf("%s_mobsnvf-scores_truths.tsv", sample_name)))


	# Create an output directory for saving the ground truth for each sample
	out_dir <- file.path(main.outdir, "ground_truth", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(truths, file.path(out_dir, sprintf("%s_ground_truth_variants.tsv", sample_name)))


	# Create output directory for roc prc and auc evaluation
	out_dir <- file.path(main.outdir, "roc-prc-auc", "precrec", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

	qwrite(mobsnvf_res$eval, file.path(out_dir, sprintf("%s_precrec_eval.rds", sample_name)))
	qwrite(mobsnvf_res$auc, file.path(out_dir, sprintf("%s_auc_table.tsv", sample_name)))
	qwrite(mobsnvf_res$roc, file.path(out_dir, sprintf("%s_roc_coordinates.tsv", sample_name)))
	qwrite(mobsnvf_res$prc, file.path(out_dir, sprintf("%s_prc_coordinates.tsv", sample_name)))

}


# for one filter
for (index in seq_len(nrow(ffpe_tumoral))){
	dobj <- set_up(index, ...);
	res <- process_sample(dobj$sample_name);		
	# write res to file
	qwrite(...)
}

# for another filter
for (index in seq_len(nrow(ffpe_tumoral))){
	dobj <- set_up(index, ...);
	res <- process_sample(dobj$sample_name);		
	# write res to file
	qwrite(...)
}


# Perform per sample evaluation
message("Evaluating: ")
## Initialize dataframe to collect AUC across all samples
auc_compilation <- data.frame()

for (index in seq_len(nrow(ffpe_tumoral))){

	# Prepare data for evaluation

	# obtain metadata for current sample
	sample_name <- ffpe_tumoral[index, "sample_name"]
	tissue <- ffpe_tumoral[index, "tissue_type"]
	variant_caller <- "MuTect2"

	message(sprintf("%d. %s", index, sample_name))

	# Construct the Ground Truth SNVs
	## To do this we make a union set for all the FF samples matched to this FFPE sample
	## This dataset has matched FFPE and FF from the same patient.
	## All samples within the same tissue type are matched
	## Hence, we select the frozen variants from the same tissue type as the FFPE sample
	frozen_tumoral_subset <- frozen_tumoral[frozen_tumoral$tissue_type == tissue, ]

	## Obtain the sample names for the frozen samples within the same tissue
	truth_samples <- frozen_tumoral_subset$sample_name
	truth_sample_paths <- file.path(vcf.dir, truth_samples, sprintf("%s.vcf.gz", truth_samples))
	truths <- snv_union(truth_sample_paths)
	truths <- add_id(truths)


	# Read in the FFPE-Filtered output for current sample
	mobsnvf_d <- read.delim(file.path(ffpe_snvf.dir, "mobsnvf", sample_name, sprintf("%s.mobsnvf.snv", sample_name)))
	mobsnvf_d <- preprocess_mobsnvf(mobsnvf_d);
	mobsnvf_res <- evaluate_filter(mobsnvf_d, truths, "mobsnvf");

	vafsnvf_d <- read.delim(file.path(ffpe_snvf.dir, "vafsnvf", sample_name, sprintf("%s.vafsnvf.snv", sample_name)))
	vafsnvf_d <- preprocess_vafsnvf(vafsnvf_d);
	vafsnvf_res <- evaluate_filter(vafsnvf_d, truths, "vafsnvf");

	# SOBDetector output column explanations:
	# 		artiStatus: Binary classification made by SOBDetector. Values are "snv" or "artifact"
	# 		SOB: This is the strand oreintation bias score column which ranges from 0 and 1. Exception values: "." or NaN. 
	sobdetector_d <- read.delim(file.path(ffpe_snvf.dir, "sobdetector", sample_name, sprintf("%s.sobdetector.snv", sample_name)))
	sobdetector_d <- preprocess_sobdetect(sobdetector_d);
	sobdetector_res <- evaluate_filter(sobdetector_d, truths, "sobdetector");


	# FIXME put into function
	## GATK Orientation Bias Mixture Model outputs
	gatk_obmm_output <- read_vcf(file.path(vcf.dir, sample_name, sprintf("%s.vcf.gz", sample_name)), columns = c("chrom", "pos", "ref", "alt", "filter"))
	### GATK Orientation Bias mixture model makes binary classification. This is casted into scores 0 and 1
	gatk_obmm_output$obmm <- ifelse(grepl("orientation", gatk_obmm_output$filter), 0, 1)
	gatk_obmm_output <- gatk_obmm_output[complete.cases(gatk_obmm_output$obmm), ]
	gatk_obmm_eval <- with(gatk_obmm_output, evalmod(scores = obmm, labels = truth, modnames = "gatk-obmm"))


	# Combine all evaluation dataframes into one
	# Get evaluation dataframe. 
	## Columns consists of:
	## x and y coordinates for making ROC and PRC plots,
	## model name, and plot type (ROC or PRC)
	all_models_eval_df <- rbind(
		as.data.frame(mobsnvf_res$eval),
		as.data.frame(vafsnvf_res$eval),
		as.data.frame(sobdetector_res$eval),
		as.data.frame(gatk_obmm_eval)
	);
	all_models_eval_df$dsid <- NULL

	# Stratify ROC PRC coordinates
	all_models_roc <- all_models_eval_df[all_models_eval_df$type == "ROC", ]
	all_models_prc <- all_models_eval_df[all_models_eval_df$type == "PRC", ]

	# Get AUC for each model
	vafsnvf_auc <- auc(vafsnvf_eval)
	sobdetector_auc <- auc(sobdetector_eval)
	gatk_obmm_auc <- auc(gatk_obmm_eval)


	# Compile model's AUC
	all_auc <- rbind(mobsnvf_auc, vafsnvf_auc, sobdetector_auc, gatk_obmm_auc)
	all_auc$dsids <- NULL
	all_auc <- reshape(all_auc, idvar = c("modnames"), timevar = "curvetypes", direction = "wide")
	colnames(all_auc) <- gsub("aucs\\.ROC", "auroc", colnames(all_auc))
	colnames(all_auc) <- gsub("aucs\\.PRC", "auprc", colnames(all_auc))
	colnames(all_auc) <- gsub("modnames", "model", colnames(all_auc))
	all_auc$sample_id <- sample_name
	all_auc <- all_auc[, c("sample_id", "model", "auroc", "auprc")]
	auc_compilation <- rbind(auc_compilation, all_auc)

	# Combine precrec eval objects, for saving as RDS
	all_models_eval <- list(
		mobsnvf = mobsnvf_eval, 
		vafsnf = vafsnvf_eval,
		sobdetector = sobdetector_eval,
		gatk_obmm = gatk_obmm_eval
	)

	# Create output directory for saving the variant set with scores and ground truth labels for each sample
	out_dir <- file.path(main.outdir, "model-scores_truths", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(mobsnvf_output, file.path(out_dir, sprintf("%s_mobsnvf-scores_truths.tsv", sample_name)))
	qwrite(vafsnvf_output, file.path(out_dir, sprintf("%s_vafsnvf-scores_truths.tsv", sample_name)))
	qwrite(sobdetector_output, file.path(out_dir, sprintf("%s_sobdetector-scores_truths.tsv", sample_name)))
	qwrite(gatk_obmm_output, file.path(out_dir, sprintf("%s_gatk-obmm-scores_truths.tsv", sample_name)))

	# Create an output directory for saving the ground truth for each sample
	out_dir <- file.path(main.outdir, "ground_truth", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(truths, file.path(out_dir, sprintf("%s_ground_truth_variants.tsv", sample_name)))

	# Create output directory for roc prc and auc evaluation
	out_dir <- file.path(main.outdir, "roc-prc-auc", "precrec", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

	qwrite(all_models_eval, file.path(out_dir, sprintf("%s_precrec_eval.rds", sample_name)))
	qwrite(all_auc, file.path(out_dir, sprintf("%s_auc_table.tsv", sample_name)))
	qwrite(all_models_roc, file.path(out_dir, sprintf("%s_roc_coordinates.tsv", sample_name)))
	qwrite(all_models_prc, file.path(out_dir, sprintf("%s_prc_coordinates.tsv", sample_name)))

}

# Create output directory for roc prc and auc evaluation
message("Compiling individual AUCs from each sample")
out_dir <- file.path(main.outdir, "roc-prc-auc", "precrec")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
qwrite(auc_compilation, file.path(out_dir, sprintf("compiled_per_sample_auc_table.tsv")))


# Perform overall evaluation

## Read in paths for the evaluation data for each sample samples
scores_truths_datasets_paths <- Sys.glob("EGAD00001004066/somatic_vcf/model-scores_truths/*/*.tsv")

model_score_columns <- c(
	"gatk-obmm" = "obmm",
	"sobdetector" = "SOB",
	"mobsnvf" = "FOBP",
	"vafsnvf" = "VAFF"
)

## Initialize data.frame to collect each models scores and ground truth label data from each sample into one
all_scores_truths <- data.frame()

## Evaluate the combined results for each models across all samples and tissue type 
message("Combining scores across all samples")
for (path in scores_truths_datasets_paths){

	## Retrieve sample name from the path
	sample_name <- unlist(strsplit(path, "/"))[4]
	base <- basename(path)
	model <- gsub("-scores", "", unlist(strsplit(base, "_"))[3],)
	
	## Read in the table
	scores_truths <- qread(path)

	
	scores_truths$model <- model
	scores_truths$sample_name <- sample_name

	## create common score column name, select necessary columns for evaluations
	colnames(scores_truths)[colnames(scores_truths) == model_score_columns[[model]]] <- "score"
	scores_truths <- scores_truths[, c('chrom','pos','ref','alt','score','snv','truth','model','sample_name')]
	

	## append the data for this sample to the bigger dataframe
	all_scores_truths <- rbind(all_scores_truths, scores_truths)

}



## Stratify across liver and colon samples
message("Stratifying combined scores based on tissue type")
all_colon_samples_scores_truths <- all_scores_truths[grepl("Colon", all_scores_truths$sample_name), ]
all_liver_samples_scores_truths <- all_scores_truths[grepl("Liver", all_scores_truths$sample_name), ]



## Save the combined score_truth set
message("Saving combined and stratified scores")
out_dir <- file.path(main.outdir, "model-scores_truths")
qwrite(all_scores_truths, file.path(out_dir, "all_samples_all-models-scores_truths.tsv"))
qwrite(all_colon_samples_scores_truths, file.path(out_dir, "all_colon_samples_all-models-scores_truths.tsv"))
qwrite(all_liver_samples_scores_truths, file.path(out_dir, "all_liver_samples_all-models-scores_truths.tsv"))


## Function to obtain ROC and PRC coordinates and AUC table
## @ Param model_scores_truths_df data.frame of variants with model scores, ground truth annotation, and model names
## @ Param model_name string with the name of model to be evaluated
## @ Param sample_name string identifier for the sample(s) to be evaluatied
## @ Return list consisting of auc, roc, and prc table for each model
get_eval_metrics <- function(model_scores_truths_df, model_name, sample_id){

	per_model_scores_truths <- all_scores_truths[all_scores_truths$model == model_name, ]

	per_model_eval <- with(per_model_scores_truths, evalmod(scores = score, labels = truth, modnames = model_name))

	per_model_eval_df <- as.data.frame(per_model_eval)
	per_model_eval_df$dsid <- NULL

	# Stratify ROC PRC coordinates
	per_model_roc <- per_model_eval_df[per_model_eval_df$type == "ROC", ]
	per_model_prc <- per_model_eval_df[per_model_eval_df$type == "PRC", ]

	# Get AUC
	per_model_auc <- auc(per_model_eval)

	## Format table
	per_model_auc$dsids <- NULL
	per_model_auc <- reshape(per_model_auc, idvar = c("modnames"), timevar = "curvetypes", direction = "wide")
	colnames(per_model_auc) <- gsub("aucs\\.ROC", "auroc", colnames(per_model_auc))
	colnames(per_model_auc) <- gsub("aucs\\.PRC", "auprc", colnames(per_model_auc))
	colnames(per_model_auc) <- gsub("modnames", "model", colnames(per_model_auc))
	per_model_auc$sample_id <- sample_id
	per_model_auc <- per_model_auc[, c("sample_id", "model", "auroc", "auprc")]

	# Return eval data
	list(
		precrec_eval_object = per_model_eval,
		auc = per_model_auc,
		roc = per_model_roc,
		prc = per_model_prc
	)
}

tissues <- c("Colon", "Liver")
models <- unique(all_scores_truths$model)

## Initialize data structures to compile evaluation results from each samples
all_model_overall_auc <- data.frame()
all_model_overall_roc <- data.frame()
all_model_overall_prc <- data.frame()
all_model_overall_precrec_eval <- list()

all_model_overall_colon_auc <- data.frame()
all_model_overall_colon_roc <- data.frame()
all_model_overall_colon_prc <- data.frame()
all_model_overall_colon_precrec_eval <- list()

all_model_overall_liver_auc <- data.frame()
all_model_overall_liver_roc <- data.frame()
all_model_overall_liver_prc <- data.frame()
all_model_overall_liver_precrec_eval <- list()


# Evaluate each tissue type and model
for (model_name in models){

	message(sprintf("Evaluating overall performance of %s:", model_name))
	
	message("	Across All Samples")
	all_sample_eval <- get_eval_metrics(all_scores_truths, model_name, "all-samples_aggregated")
	
	message("	Across All Colon Samples")
	all_colon_samples_eval <- get_eval_metrics(all_colon_samples_scores_truths, model_name, "all-colon-samples_aggregated")
	
	message("	Across All Liver Samples")
	all_liver_samples_eval <- get_eval_metrics(all_liver_samples_scores_truths, model_name, "all-liver-samples_aggregated")
	
	
	## Compile ROC, PRC, and AUC
	message("Compiling ROC, PRC, and AUC")
	### Across all samples
	all_model_overall_auc <- rbind(all_model_overall_auc, all_sample_eval$auc)
	all_model_overall_roc <- rbind(all_model_overall_roc, all_sample_eval$roc)
	all_model_overall_prc <- rbind(all_model_overall_prc, all_sample_eval$prc)
	all_model_overall_precrec_eval[[model_name]] <- all_sample_eval$precrec_eval_object

	### Across all colon samples
	all_model_overall_colon_auc <- rbind(all_model_overall_colon_auc, all_colon_samples_eval$auc)
	all_model_overall_colon_roc <- rbind(all_model_overall_colon_roc, all_colon_samples_eval$roc)
	all_model_overall_colon_prc <- rbind(all_model_overall_colon_prc, all_colon_samples_eval$prc)
	all_model_overall_colon_precrec_eval[[model_name]] <- all_sample_eval$precrec_eval_object

	### Across all liver samples
	all_model_overall_liver_auc <- rbind(all_model_overall_liver_auc, all_liver_samples_eval$auc)
	all_model_overall_liver_roc <- rbind(all_model_overall_liver_roc, all_liver_samples_eval$roc)
	all_model_overall_liver_prc <- rbind(all_model_overall_liver_prc, all_liver_samples_eval$prc)
	all_model_overall_liver_precrec_eval[[model_name]] <- all_sample_eval$precrec_eval_object

}


## Set output directory and save roc prc and auc evaluation for each model across:
message("Saving compiled evaluation metrics")
### all samples
out_dir <- file.path(main.outdir, "roc-prc-auc", "precrec")
qwrite(all_model_overall_precrec_eval, file.path(out_dir, "all_samples_all_models_precrec_eval.rds"))
qwrite(all_model_overall_auc, file.path(out_dir, "all_samples_all_models_auc_table.tsv"))
qwrite(all_model_overall_roc, file.path(out_dir, "all_samples_all_models_roc_coordinates.tsv"))
qwrite(all_model_overall_prc, file.path(out_dir, "all_samples_all_models_prc_coordinates.tsv"))

### all colon samples
qwrite(all_model_overall_colon_precrec_eval, file.path(out_dir, "all_colon_samples_all_models_precrec_eval.rds"))
qwrite(all_model_overall_colon_auc, file.path(out_dir, "all_colon_samples_all_models_auc_table.tsv"))
qwrite(all_model_overall_colon_roc, file.path(out_dir, "all_colon_samples_all_models_roc_coordinates.tsv"))
qwrite(all_model_overall_colon_prc, file.path(out_dir, "all_colon_samples_all_models_prc_coordinates.tsv"))

### all liver samples
qwrite(all_model_overall_liver_precrec_eval, file.path(out_dir, "all_liver_samples_all_models_precrec_eval.rds"))
qwrite(all_model_overall_liver_auc, file.path(out_dir, "all_liver_samples_all_models_auc_table.tsv"))
qwrite(all_model_overall_liver_roc, file.path(out_dir, "all_liver_samples_all_models_roc_coordinates.tsv"))
qwrite(all_model_overall_liver_prc, file.path(out_dir, "all_liver_samples_all_models_prc_coordinates.tsv"))

