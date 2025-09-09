#!/usr/bin/env Rscript

# Import necessary libraries and modules
library(io)
library(precrec)
library(argparser)

# Importing module containing common functions for evaluations
source("../R/eval.R")


# Setup Directories

# Directory for FFPE SNVF outputs
dataset_id <- "EGAD00001004066"
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

	# Read in the FFPE-Filtered output for current sample

	mobsnvf_output <- read.delim(file.path(ffpe_snvf.dir, "mobsnvf", sample_name, sprintf("%s.mobsnvf.snv", sample_name)))
	# Retain only C>T mutations in mobsnvf output
	mobsnvf_output <- mobsnvf_output[complete.cases(mobsnvf_output$FOBP), ]


	vafsnvf_output <- read.delim(file.path(ffpe_snvf.dir, "vafsnvf", sample_name, sprintf("%s.vafsnvf.snv", sample_name)))
	vafsnvf_output <- vafsnvf_output[complete.cases(vafsnvf_output$VAFF), ]

	# SOBDetector output column explanations:
	# 		artiStatus: Binary classification made by SOBDetector. Values are "snv" or "artifact"
	# 		SOB: This is the strand oreintation bias score column which ranges from 0 and 1. Exception values: "." or NaN. 
	sobdetector_output <- read.delim(file.path(ffpe_snvf.dir, "sobdetector", sample_name, sprintf("%s.sobdetector.snv", sample_name)))
	# variants ignored by SOBdetector have score of "."
	sobdetector_output <- sobdetector_output[!(sobdetector_output$SOB == "."), ]
	# now, it's safe to convert to numeric
	sobdetector_output$SOB <- as.numeric(sobdetector_output$SOB)
	# SOBdetector score = 0 indicates that it is not artifact
	# variants classified by SOBdetect as "snv" have score of NaN
	sobdetector_output$SOB <- ifelse(is.nan(sobdetector_output$SOB), 0, sobdetector_output$SOB)
	### Precaution to remove the any NA Scores if present
	sobdetector_output <- sobdetector_output[complete.cases(sobdetector_output$SOB), ]


	## GATK Orientation Bias Mixture Model outputs
	gatk_obmm_output <- read_vcf(file.path(vcf.dir, sample_name, sprintf("%s.vcf.gz", sample_name)), columns = c("chrom", "pos", "ref", "alt", "filter"))
	### GATK Orientation Bias mixture model makes binary classification. This is casted into scores 0 and 1
	gatk_obmm_output$obmm <- ifelse(grepl("orientation", gatk_obmm_output$filter), 0, 1)
	gatk_obmm_output <- gatk_obmm_output[complete.cases(gatk_obmm_output$obmm), ]

	# Add common IDs to the outputs
	mobsnvf_output <- add_id(mobsnvf_output)
	vafsnvf_output <- add_id(vafsnvf_output)
	sobdetector_output <- add_id(sobdetector_output)
	gatk_obmm_output <- add_id(gatk_obmm_output)


	# Adjust scores
	## In this set the TRUE label indicates real mutation and FALSE indicates artifacts
	## Hence, scores needs to be adjusted so that higher score represents a real mutation.
	## For our VAFSNVF higher score is signifies real mutation : 
	## For MOBSNVF and SOBDetector lower score signifies real mutation:  
	## Hence, we flip scores for MOBSNVF and SOBDetector to make higher score signify a real mutation
	mobsnvf_output$FOBP <- -mobsnvf_output$FOBP
	sobdetector_output$SOB <- -sobdetector_output$SOB
	

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

	## Add truth labels to model names
	mobsnvf_output <- annotate_truth(mobsnvf_output, truths)
	vafsnvf_output <- annotate_truth(vafsnvf_output, truths)
	sobdetector_output <- annotate_truth(sobdetector_output, truths)
	gatk_obmm_output <- annotate_truth(gatk_obmm_output, truths)
	

	## Evaluate using precrec
	mobsnvf_eval <- with(mobsnvf_output, evalmod(scores = FOBP, labels = truth, modnames = "mobsnvf"))
	vafsnvf_eval <- with(vafsnvf_output, evalmod(scores = VAFF, labels = truth, modnames = "vafsnvf"))
	sobdetector_eval <- with(sobdetector_output, evalmod(scores = SOB, labels = truth, modnames = "sobdetector"))
	gatk_obmm_eval <- with(gatk_obmm_output, evalmod(scores = obmm, labels = truth, modnames = "gatk-obmm"))


	# Get evaluation dataframe. 
	## Columns consisits of:
	## x and y coordinates for making ROC and PRC plots,
	## model name, and plot type (ROC or PRC)
	mobsnvf_eval_df <- as.data.frame(mobsnvf_eval)
	vafsnvf_eval_df <- as.data.frame(vafsnvf_eval)
	sobdetector_eval_df <- as.data.frame(sobdetector_eval)
	gatk_obmm_eval_df <- as.data.frame(gatk_obmm_eval)

	# Combine all evaluation dataframes into one
	all_models_eval_df <- rbind(mobsnvf_eval_df, vafsnvf_eval_df, sobdetector_eval_df, gatk_obmm_eval_df)
	all_models_eval_df$dsid <- NULL

	# Stratify ROC PRC coordinates
	all_models_roc <- all_models_eval_df[all_models_eval_df$type == "ROC", ]
	all_models_prc <- all_models_eval_df[all_models_eval_df$type == "PRC", ]

	# Get AUC for each model
	mobsnvf_auc <- auc(mobsnvf_eval)
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
## @ Return list consisting of auc, roc, and prc table for each model
get_eval_metrics <- function(model_scores_truths_df, model_name){

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
	per_model_auc$sample_id <- sample_name
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
	all_sample_eval <- get_eval_metrics(all_scores_truths, model_name)
	
	message("	Across All Colon Samples")
	all_colon_samples_eval <- get_eval_metrics(all_colon_samples_scores_truths, model_name)
	
	message("	Across All Liver Samples")
	all_liver_samples_eval <- get_eval_metrics(all_liver_samples_scores_truths, model_name)
	
	
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

