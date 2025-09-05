#!/usr/bin/env Rscript

# Import necessary libraries and modules
library(io)
library(precrec)
library(jsonlite)
library(argparser)
library(stringr)

## Importing module containing common functions for evaluations
source("../../../R/eval.R")

### Setup Directories

## Directory for FFPE SNVF outputs
ffpe_snvf.dir <- "../../ffpe-snvf/somatic_vcf"
## Directory for the Somatic VCFs
vcf.dir <- "../../data/somatic_vcf"
## Output directory
main.outdir <- "../../results/somatic_vcf"

### Read Annotation Table

lookup_table <- read.delim("../../annot/sample_annotations.tsv")
## create a sample name column
lookup_table$sample_name <- paste0(gsub(" ", "-", lookup_table$title), "_", lookup_table$sample_alias)



## Stratify annotation table based on FFPE and FF Somatic Variants
ffpe_tumoral <- lookup_table[(lookup_table$preservation == "FFPE" & lookup_table$sample_type == "Tumoral"), ]
frozen_tumoral <- lookup_table[(lookup_table$preservation == "Frozen" & lookup_table$sample_type == "Tumoral"), ]


print(glue("Evaluating: "))

## Perform per sample evaluation
for (index in seq_len(dim(ffpe_tumoral)[1])){

	### Prepare data for evaluation

	## obtain metadata for current sample
	sample_name <- ffpe_tumoral[index, "sample_name"]
	tissue <- ffpe_tumoral[index, "tissue_type"]
	variant_caller <- "MuTect2"

	print(sprintf("%d. %s", index, sample_name))

	# Read in the FFPE-Filtered output for current sample

	mobsnvf_output <- read.delim(file.path(ffpe_snvf.dir, "mobsnvf", sample_name, sprintf("%s.mobsnvf.snv", sample_name)))
	### Retain only C>T mutations in mobsnvf output
	mobsnvf_output <- mobsnvf_output[complete.cases(mobsnvf_output), ]


	vafsnvf_output <- read.delim(file.path(ffpe_snvf.dir, "vafsnvf", sample_name, sprintf("%s.vafsnvf.snv", sample_name)))

	## SOBDetector output column explanations:
	## 		artiStatus: Binary classification made by SOBDetector. Values are "snv" or "artifact"
	## 		SOB: This is the strand oreintation bias score column which ranges from 0 and 1. Exception values: "." or NaN. 
	sobdetector_output <- read.delim(file.path(ffpe_snvf.dir, "sobdetector", sample_name, sprintf("%s.sobdetector.snv", sample_name)))
	## Some Variants are not evaluated by SOBdetector.
	## These are denoted with a "." under the SOB (score) column. These variants are removed.
	sobdetector_output <- sobdetector_output[!(sobdetector_output$SOB == "."), ]
	### SOB Column is changed back to numeric as it was read in as character due to presence of "."
	sobdetector_output$SOB <- as.numeric(sobdetector_output$SOB)
	### Change values where SOB score is nan to 0 as artiStatus column for these values are annotated as SNV i.e real mutation
	### By default, higher SOB score indicates likelihood of artifact. This is why real mutations are set to 0
	sobdetector_output$SOB <- ifelse(is.nan(sobdetector_output$SOB), 0, sobdetector_output$SOB)


	## GATK Orientation Bias Mixture Model outputs
	gatk_obmm_output <- read_vcf(file.path(vcf.dir, sample_name, sprintf("%s.vcf.gz", sample_name)), columns = c("chrom", "pos", "ref", "alt", "filter"))
	### GATK Orientation Bias mixture model makes binary classification. This is casted into scores 0 and 1
	gatk_obmm_output$obmm <- ifelse(grepl("orientation", gatk_obmm_output$filter), 0, 1)

	# Add common IDs to the outputs
	mobsnvf_output <- add_id(mobsnvf_output)
	vafsnvf_output <- add_id(vafsnvf_output)
	sobdetector_output <- add_id(sobdetector_output)
	gatk_obmm_output <- add_id(gatk_obmm_output)


	# Now combine all outputs into a single dataframe
	## This removes C>T mutation from the output of other models
	## This also ensures that a uniform variant set is maintained. As some variants may be omitted by some models
	all_model_scores <- mobsnvf_output[, c("chrom", "pos", "ref", "alt", "snv", "FOBP")]
	all_model_scores <- merge(all_model_scores, vafsnvf_output[, c("snv", "VAFF")], by = "snv", all.x = TRUE)
	all_model_scores <- merge(all_model_scores, sobdetector_output[, c("snv", "SOB")], by = "snv", all.x = TRUE)
	all_model_scores <- merge(all_model_scores, gatk_obmm_output[, c("snv", "obmm")], by = "snv", all.x = TRUE)


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


	# Add truth labels variants in the all_model_scores dataframe based on whether the variant is present in the ground truth
	all_model_scores_truths <- annotate_truth(damaged_sample_variants = all_model_scores, ground_truth_variants = truths)


	# Adjust scores
	## In this set the TRUE label indicates real mutation and FALSE indicates artifacts
	## Hence, scores needs to be adjusted so that higher score represents a real mutation.
	## For our VAFSNVF higher score is signifies real mutation : 
	## For MOBSNVF and SOBDetector lower score signifies real mutation:  
	## Hence, we flip scores for MOBSNVF and SOBDetector to make higher score signify a real muatation
	all_model_scores_truths$FOBP <- -all_model_scores_truths$FOBP
	all_model_scores_truths$SOB <- -all_model_scores_truths$SOB
	
	
	## Create output directory for saving the evaluation variant set for each sample
	out_dir <- file.path(main.outdir, "model_scores_labels_truths", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(all_model_scores_truths, file.path(out_dir, sprintf("%s_model_scores_labels_truths.tsv", sample_name)))

	## Create an output directory for saving the ground truth for each sample
	out_dir <- file.path(main.outdir, "ground_truth", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	qwrite(truths, file.path(out_dir, sprintf("%s_ground_truth_variants.tsv", sample_name)))


	# Evaluation using Precrec
	eval_data_df <- all_model_scores_truths
	precrec_eval_metrics <- get.precrec.eval.metrics(eval_data_df, sample_name, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), model_names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))

	## Create output directory for roc prc and auc evaluation
	out_dir <- file.path(main.outdir, "roc-prc-auc", "precrec", sample_name)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

	qwrite(precrec_eval_metrics$auc_table, file.path(out_dir, sprintf("%s_auc_table.tsv", sample_name)))
	qwrite(precrec_eval_metrics$precrec_eval_object, file.path(out_dir, sprintf("%s_precrec_eval.rds", sample_name)))
	qwrite(precrec_eval_metrics$roc, file.path(out_dir, sprintf("%s_roc_coordinates.tsv", sample_name)))
	qwrite(precrec_eval_metrics$prc, file.path(out_dir, sprintf("%s_prc_coordinates.tsv", sample_name)))

}


## Perform overall evaluation

### Read in paths for the evaluation data for each sample samples
scores_labels_truths_datasets_paths <- Sys.glob("../../results/somatic_vcf/model_scores_labels_truths/*/*.tsv")
### Concatenate score_labels_truth dataset for all samples into one
all_scores_labels_truths <- data.frame()

for (path in scores_labels_truths_datasets_paths){

	## Retrieve sample name from the path
	sample_name <- str_split_1(path, "/")[7]

	## Read in the table
	scores_labels_truths <- qread(path)

	## Annotate the sample name
	scores_labels_truths$sample_name <- sample_name

	## append the data for this sample to the bigger dataframe
	all_scores_labels_truths <- rbind(all_scores_labels_truths, scores_labels_truths)

}


## Stratify based on tissue types
liver_samples_scores_labels_truths <- all_scores_labels_truths[str_detect(all_scores_labels_truths$sample_name, "Liver"), ]
colon_samples_scores_labels_truths <- all_scores_labels_truths[str_detect(all_scores_labels_truths$sample_name, "Colon"), ]

## Create overall eval for:
print("Making aggregated evaluations")
### All Samples
precrec_all_samples_eval <- get.precrec.eval.metrics(all_scores_labels_truths, "EGAD00001004066 All FFPE tomoral Samples", score_columns = c("FOBP", "VAFF", "SOB", "obmm"), model_names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))
### Liver Samples
precrec_liver_samples_eval <- get.precrec.eval.metrics(liver_samples_scores_labels_truths, "EGAD00001004066 all FFPE tomoral liver Samples", score_columns = c("FOBP", "VAFF", "SOB", "obmm"), model_names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))
### Colon Samples
precrec_colon_samples_eval <- get.precrec.eval.metrics(colon_samples_scores_labels_truths, "EGAD00001004066 all FFPE tomoral colon Samples", score_columns = c("FOBP", "VAFF", "SOB", "obmm"), model_names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))


## Write aggregated eval input data to file

out_dir <- file.path(main.outdir, "model_scores_labels_truths")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

qwrite(all_scores_labels_truths, file.path(out_dir, "all_samples_scores_labels_truths.tsv"))
qwrite(liver_samples_scores_labels_truths, file.path(out_dir, "liver_samples_scores_labels_truths.tsv"))
qwrite(colon_samples_scores_labels_truths, file.path(out_dir, "colon_samples_scores_labels_truths.tsv"))


## Write aggregated eval data to file
out_dir <- file.path(main.outdir, "roc-prc-auc", "precrec")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### Aggregated all samples evaluation
qwrite(precrec_all_samples_eval$auc_table, file.path(out_dir, "all_samples_auc_table.tsv"))
qwrite(precrec_all_samples_eval$precrec_eval_object, file.path(out_dir, "all_samples_precrec_eval.rds"))
qwrite(precrec_all_samples_eval$roc, file.path(out_dir, "all_samples_roc_coordinates.tsv"))
qwrite(precrec_all_samples_eval$prc, file.path(out_dir, "all_samples_prc_coordinates.tsv"))

### Aggregated liver samples evaluation
qwrite(precrec_liver_samples_eval$auc_table, file.path(out_dir, "liver_samples_auc_table.tsv"))
qwrite(precrec_liver_samples_eval$precrec_eval_object, file.path(out_dir, "liver_samples_precrec_eval.rds"))
qwrite(precrec_liver_samples_eval$roc, file.path(out_dir, "liver_samples_roc_coordinates.tsv"))
qwrite(precrec_liver_samples_eval$prc, file.path(out_dir, "liver_samples_prc_coordinates.tsv"))

### Aggregated colon samples evaluation 
qwrite(precrec_colon_samples_eval$auc_table, file.path(out_dir, "colon_samples_auc_table.tsv"))
qwrite(precrec_colon_samples_eval$precrec_eval_object, file.path(out_dir, "colon_samples_precrec_eval.rds"))
qwrite(precrec_colon_samples_eval$roc, file.path(out_dir, "colon_samples_roc_coordinates.tsv"))
qwrite(precrec_colon_samples_eval$prc, file.path(out_dir, "colon_samples_prc_coordinates.tsv"))

print("Done")



