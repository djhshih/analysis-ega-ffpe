#!/usr/bin/env Rscript

# This script contains functions for evaluation
## This is a work in progress. The functions still needs to be improved for generalizability across multiple projects and datasets

### Load necessary libraries
library(io)
library(precrec)
library(jsonlite)
library(argparser)
library(tidyverse)
library(glue)


### Evaluations

#### A list of Chromosomes 1-22, X and Y
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))

#### Create unique identifier for each variant
add_id <- function(damaged_sample_variants) {
	damaged_sample_variants$snv <- paste(damaged_sample_variants$chrom, damaged_sample_variants$pos, damaged_sample_variants$ref, damaged_sample_variants$alt, sep="_")
	damaged_sample_variants
}

#### Annotate the presence of a variant in the truth set i.e in the Fresh Frozen Sample
annotate_truth <- function(damaged_sample_variants, ground_truth_variants) {
	damaged_sample_variants$truth <- damaged_sample_variants$snv %in% ground_truth_variants$snv;
	damaged_sample_variants
}

#### Function to filter C>T and G>A mutations
ct_only <- function(df) {
  df[(df$ref == "C" & df$alt == "T") | (df$ref == "G" & df$alt == "A"), ]
}

#### Create a set of SNVs across multiple samples
snv_union <- function(matched_undamaged_sample_paths) {

	truths <- data.frame()

	for (truth_sample_path in matched_undamaged_sample_paths){
		truth_sample <- qread(glue(truth_sample_path), type = "vcf") |> select(chrom, pos, ref, alt)
		truths <- rbind(truths, truth_sample)
	}

	truths <- distinct(truths)

	return(truths)
}

#### Function to classify variants based on a set of Q values and a False Positive cutoff threshold
##### Returns a boolean vector
adaptive_fdr_cut <- function(q, fp.cut) {
	n <- length(q);
	# adaptive threshold chosen to expect less than one false positive
	idx <- order(q);
	under <- which((1:n) * q[idx] < fp.cut);
	if (length(under) > 0) {
		top <- under[length(under)];
		pred <- logical(n);
		##### select the top significant results
		pred[idx[1:top]] <- TRUE;
		pred
	} else {
		##### none passes
		rep(FALSE, n)
	}
}

#### Function to label the models prediction based on the adaptive_fdr_cut function
fdr_cut_pred <- function(df, score, fp.cut=0.5) {
	
    ##### Split C>T and non C>T mutations into two dataframes
    df.ct <- df |> filter(complete.cases(.data[[score]]))
    df.nct <- df |> filter(!complete.cases(.data[[score]])) |> mutate(score = NA, q = NA, pred = TRUE)

    ##### Predict artifacts based on q-values and adaptive_fdr_cut function
    df.ct.q <- df.ct |>
        mutate(
            score = ifelse(.data[[score]] == 0, .Machine$double.eps, .data[[score]]),
            q = p.adjust(score, "BH"),
			## TRUE being real muatation and FALSE being artifacts
            pred = ifelse(adaptive_fdr_cut(q, fp.cut), TRUE, FALSE)
        )

    ##### Combine the C>T and non C>T dataframes back into one single dataframe
    rbind(df.ct.q, df.nct) |>
		arrange(chrom, pos)
}

#### Calculate True/False Positives/Negatives from prediction and ground truth
calc.confusion.matrix <- function(df, pred_col = "pred", truth_col = "truth") {
	TP <- df |> filter((.data[[pred_col]]) & (.data[[truth_col]])) |> count() |> pull()
	FP <- df |> filter((.data[[pred_col]]) & !(.data[[truth_col]])) |> count() |> pull()
	FN <- df |> filter(!(.data[[pred_col]]) & (.data[[truth_col]])) |> count() |> pull()
	TN <- df |> filter(!(.data[[pred_col]]) & !(.data[[truth_col]])) |> count() |> pull()
	
	list(
		TP = TP,
		FP = FP,
		FN = FN,
		TN = TN
	)
}

#### Calculate evaluation metrics based on confusion matrix
calc.eval.metrics <- function(df, pred_col = "pred", truth_col = "truth") {

	confusion_matrix <- calc.confusion.matrix(df, pred_col, truth_col)

	precision <- confusion_matrix$TP / (confusion_matrix$TP + confusion_matrix$FP)
	recall <- confusion_matrix$TP / (confusion_matrix$TP + confusion_matrix$FN)
	sensitivity <- recall  # sensitivity is same as recall
	specificity <- confusion_matrix$TN / (confusion_matrix$TN + confusion_matrix$FP)
	
	list(
		precision = precision,
		recall = recall,
		sensitivity = sensitivity,
		specificity = specificity
	)
}


#### Manual evaluation function without using the precrec library
##### This returns a list of dataframes containing coordinates for the ROC and PRC plots
evaluate.roc.prc <- function(model_scores_truth, score_col = "score", truth_col = "truth", model_name = NA) {

	##### Sort the model's scores
	model_scores_truth <- model_scores_truth[order(model_scores_truth[[score_col]]), ]

	##### Checking the distinct score values
	model_unique_scores <- unique(model_scores_truth[[score_col]])
	thresholds <- c(-Inf, unique(model_scores_truth[[score_col]]), Inf)

	##### Calculate Precision, Recall, TPR and FPR based on each cutoff based on the threshold vector
	results_list <- list()

	for (i in seq_along(thresholds)) {
		score <- thresholds[i]
		
		##### Make Predictions based on the cutoff i.e TRUE for real mutation and False for Artifact
		model_scores_truth_pred <- mutate(model_scores_truth, pred = ifelse(.data[[score_col]] >= score, TRUE, FALSE))
		metrics <- calc.eval.metrics(model_scores_truth_pred, pred_col = "pred", truth_col = "truth")

		##### Collect the result in a list
		results_list[[i]] <- data.frame(
			model = model_name,
			threshold = score,
			tpr = metrics$sensitivity,
			fpr = 1 - metrics$specificity,
			precision = metrics$precision,
			recall = metrics$recall
		)
	}

	##### Combine everything at the end
	all_metrics_df <- bind_rows(results_list)
	##### Replace NaN value for recall = 0 with with the previous precision value to that curve is extrapolated horizontally to the y-axis
	all_metrics_df[nrow(all_metrics_df), "precision"] <- all_metrics_df[nrow(all_metrics_df) - 1, "precision"]

	##### Split into ROC and PRC dataframes
	model_roc <- all_metrics_df |> select(model, tpr, fpr) |> arrange(fpr)
	model_prc <- all_metrics_df |> select(model, precision, recall) |> arrange(recall)

	##### Collect results to a list
	list(
		all_metrics = all_metrics_df,
		roc = model_roc,
		prc = model_prc
	)
}

#### Function to calculate AUC using the trapezoidal rule
calculate.auc <- function(df, x_col, y_col) {
  	##### Using dplyr's lead() function to get the next value in the sequence
	df |>
		arrange(.data[[x_col]]) |>
		mutate(
			##### Calculate the width of each trapezoid (change in x)
			dx = lead(.data[[x_col]]) - .data[[x_col]],
			##### Calculate the average height of each trapezoid (average of y)
			avg_y = (lead(.data[[y_col]]) + .data[[y_col]]) / 2,
			##### Calculate the area of each small trapezoid
			area = dx * avg_y
		) |>
		##### Sum up the areas of all trapezoids, removing the NA from the last entry
		summarise(auc = sum(area, na.rm = TRUE)) |>
		pull(auc)
}

## Function to create multimodel AUROC and AUPRC table based on a dataframe containing ROC and PRC
## Column names for ROC df must be fpr and tpr
## Column names for ROC df must be precision and recall
get.multimodel.auroc.auprc <- function(multimodel_roc_df,  multimodel_prc_df, model_names = NULL, model_col = "model") {

	if(is.null(model_names)) {
		model_names <- unique(multimodel_roc_df$model)
	}

	auc_list <- list()
	model_col <- "model"

	for (i in seq_along(model_names)) {
		model_name <- model_names[i]

		model_roc <- multimodel_roc_df |> filter(.data[[model_col]] == model_name)
		model_prc <- multimodel_prc_df |> filter(.data[[model_col]] == model_name)

		model_auroc <- calculate.auc(model_roc, x_col = "fpr", y_col = "tpr")
		model_auprc <- calculate.auc(model_prc, x_col = "recall", y_col = "precision")

		auc_list[[i]] <- data.frame(
			sample_id = sample_name,
			model = model_name,
			auroc = model_auroc,
			auprc = model_auprc
		)

	}

	aucs <- bind_rows(auc_list)

	return(aucs)
}


#### Creates a multimodel evaluation object using the precrec library
make.multimodel.eval.object <- function(scores_labels_df, score_columns = c("FOBP", "VAFF", "SOB", "msec"), names = c("mobsnvf", "vafsnvf", "sobdetector", "microsec")) {
	all.model.scores <- do.call(join_scores, scores_labels_df[score_columns])
	all.model.labels <- do.call(join_labels, replicate(length(score_columns), scores_labels_df$truth, simplify = FALSE))
	all.model.names <- names

	evalmod(mmdata(all.model.scores, all.model.labels, all.model.names))
}


## This function calculate the evaluation metrics from precrec
## The function returns the ROC, PRC coordinates for each model and the AUC table
get.precrec.eval.metrics <- function(models_scores_truths_df, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), model_names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm")) {
	
	## Create PRECREC all model eval object
	precrec_all_model_eval <- make.multimodel.eval.object(scores_labels_df = models_scores_truths_df, score_columns = score_columns, names = model_names)


	## Retrive the AUCs calculated by precrec and append sample name
	precrec_aucs <- auc(precrec_all_model_eval) |> 
		pivot_wider(names_from = curvetypes, values_from = aucs) |> select(-dsids) |> 
		mutate(sample_id = sample_name) |>
		rename(model = modnames, auroc = ROC, auprc = PRC) |>
		select(sample_id, model, auroc, auprc)


	## Obtain the ROC coordinates from the precrec eval object
	precrec_roc <- as.data.frame(precrec_all_model_eval) |>
		filter(type == "ROC") |> select(-dsid) |> 
		mutate(sample_id = sample_name) |>
		select(sample_id, modname, x, y) |>
		rename(fpr = x, tpr = y, model = modname)


	## Obtain the PRC coordinates from the precrec eval object
	precrec_prc <- as.data.frame(precrec_all_model_eval) |>
		filter(type == "PRC") |> select(-dsid) |> 
		mutate(sample_id = sample_name) |>
		select(sample_id, modname, x, y) |>
		rename(recall = x, precision = y, model = modname)

	## Return the table and the roc and prc coordinates as a list
	list(
		auc_table = precrec_aucs,
		roc = precrec_roc,
		prc = precrec_prc
	)

}