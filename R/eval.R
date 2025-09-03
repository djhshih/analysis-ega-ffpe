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


#### A list of Chromosomes 1-22, X and Y
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))

read_vcf <- function(path, columns = NULL) {
	all_lines <- readLines(path)
	filtered_lines <- grep("^##", all_lines, value = TRUE, invert = TRUE)

	vcf <- read.delim(
		text = filtered_lines,
		sep = "\t",
		check.names = FALSE
	)

	names(vcf)[1] <- sub("^#", "", names(vcf)[1])
	names(vcf) <- tolower(names(vcf))

	if (!is.null(columns)) {
		vcf <- vcf[, columns]
	}

	return(vcf)
}


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
		# Read the VCF data
		truth_sample_full <- qread(glue(truth_sample_path), type = "vcf")
		# Select only the required columns using base R
		truth_sample <- truth_sample_full[, c("chrom", "pos", "ref", "alt")]
		truths <- rbind(truths, truth_sample)
	}

	# Use base R's `unique` for data frames, which is equivalent to dplyr's `distinct`
	truths <- unique(truths)

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
	
    ##### Split C>T and non C>T mutations into two dataframes using base R
    df.ct <- df[complete.cases(df[[score]]), ]
    df.nct <- df[!complete.cases(df[[score]]), ]

    ##### For non C>T, mutate columns using direct assignment
    df.nct$score <- NA
    df.nct$q <- NA
    df.nct$pred <- TRUE

    ##### Predict artifacts based on q-values and adaptive_fdr_cut function
    # The original dataframe is modified by direct assignment
    df.ct$score <- ifelse(df.ct[[score]] == 0, .Machine$double.eps, df.ct[[score]])
    df.ct$q <- p.adjust(df.ct$score, "BH")
	## TRUE being real muatation and FALSE being artifacts
    df.ct$pred <- ifelse(adaptive_fdr_cut(df.ct$q, fp.cut), TRUE, FALSE)

    ##### Combine the C>T and non C>T dataframes back into one single dataframe
    df <- rbind(df.ct, df.nct)
	# Arrange the combined dataframe using base R's `order`
	df <- df[order(df$chrom, df$pos), ]
	df
}


#### Calculate True/False Positives/Negatives from prediction and ground truth
calc.confusion.matrix <- function(df, pred_col = "pred", truth_col = "truth") {
	TP <- nrow(df[df[[pred_col]] & df[[truth_col]], ])
	FP <- nrow(df[df[[pred_col]] & !df[[truth_col]], ])
	FN <- nrow(df[!df[[pred_col]] & df[[truth_col]], ])
	TN <- nrow(df[!df[[pred_col]] & !df[[truth_col]], ])
	
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


#### Creates a multimodel evaluation object using the precrec library
make.multimodel.eval.object <- function(scores_labels_df, score_columns = c("FOBP", "VAFF", "SOB", "msec"), names = c("mobsnvf", "vafsnvf", "sobdetector", "microsec")) {
	all.model.scores <- do.call(join_scores, scores_labels_df[score_columns])
	all.model.labels <- do.call(join_labels, replicate(length(score_columns), scores_labels_df$truth, simplify = FALSE))
	all.model.names <- names

	evalmod(mmdata(all.model.scores, all.model.labels, all.model.names))
}


## This function calculate the evaluation metrics from precrec
## The function returns the ROC, PRC coordinates for each model and the AUC table
get.precrec.eval.metrics <- function(models_scores_truths_df, sample_name, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), model_names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm")) {
	
	## Create PRECREC all model eval object
	precrec_all_model_eval <- make.multimodel.eval.object(scores_labels_df = models_scores_truths_df, score_columns = score_columns, names = model_names)


	## Retrive the AUCs calculated by precrec and append sample name

    # 1. Get AUCs from precrec object
	aucs_df <- auc(precrec_all_model_eval)
    # 2. Reshape from long to wide format using base R's reshape
    aucs_wider <- reshape(aucs_df, idvar = c("modnames", "dsids"), timevar = "curvetypes", direction = "wide")
    # 3. Remove dsids column
    aucs_wider$dsids <- NULL
    # 4. Add sample_id column (using the new function argument)
    aucs_wider$sample_id <- sample_name
    # 5. Rename columns to desired format
    names(aucs_wider)[names(aucs_wider) == "modnames"] <- "model"
    names(aucs_wider)[names(aucs_wider) == "aucs.ROC"] <- "auroc"
    names(aucs_wider)[names(aucs_wider) == "aucs.PRC"] <- "auprc"
    # 6. Select and reorder final columns
	precrec_aucs <- aucs_wider[, c("sample_id", "model", "auroc", "auprc")]


	## Obtain the ROC and PRC coordinates from the precrec eval object

    ## Get full data frame from precrec object
	precrec_df <- as.data.frame(precrec_all_model_eval)
    
	## Obtain the ROC coordinates from the precrec eval object
	# 1. Filter for PRC curve data from precrec_df
	roc_df <- precrec_df[precrec_df$type == "ROC", ]
    # 2. Add sample_id and remove dsid
    roc_df$dsid <- NULL
    roc_df$sample_id <- sample_name
    # 3. Select and rename columns
	roc_df <- roc_df[, c("sample_id", "modname", "x", "y")]
    names(roc_df) <- c("sample_id", "model", "fpr", "tpr")
    precrec_roc <- roc_df


	## Obtain the PRC coordinates from the precrec eval object
    # 1. Filter for PRC curve data from the same precrec_df
	prc_df <- precrec_df[precrec_df$type == "PRC", ]
    # 2. Add sample_id and remove dsid
    prc_df$dsid <- NULL
    prc_df$sample_id <- sample_name
    # 3. Select and rename columns
	prc_df <- prc_df[, c("sample_id", "modname", "x", "y")]
    names(prc_df) <- c("sample_id", "model", "recall", "precision")
    precrec_prc <- prc_df

	## Return the table and the roc and prc coordinates as a list
	list(
		precrec_eval_object = precrec_all_model_eval,
		auc_table = precrec_aucs,
		roc = precrec_roc,
		prc = precrec_prc
	)

}

#### Manual evaluation function without using the precrec library
##### This returns a list of dataframes containing coordinates for the ROC and PRC plots
evaluate.roc.prc <- function(model_scores_truth, score_col = "score", truth_col = "truth", model_name = NA) {

	##### Sort the model's scores
	model_scores_truth <- model_scores_truth[order(model_scores_truth[[score_col]]), ]

	##### Checking the distinct score values
	thresholds <- c(-Inf, unique(model_scores_truth[[score_col]]), Inf)

	##### Calculate Precision, Recall, TPR and FPR based on each cutoff based on the threshold vector
	results_list <- list()

	for (i in seq_along(thresholds)) {
		score <- thresholds[i]
		
		##### Make Predictions based on the cutoff i.e TRUE for real mutation and False for Artifact
		model_scores_truth_pred <- model_scores_truth
        model_scores_truth_pred$pred <- ifelse(model_scores_truth_pred[[score_col]] >= score, TRUE, FALSE)
		
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

	##### Combine everything at the end using base R's `do.call` with `rbind`
	all_metrics_df <- do.call(rbind, results_list)
    
	##### Replace NaN value for precision where recall = 0 with with the previous precision value in the table
	##### This is done so that the curve is extrapolated horizontally to be anchored to the y-axis
	##### Allowing for proper calculation of AUPRC
	all_metrics_df[nrow(all_metrics_df), "precision"] <- all_metrics_df[nrow(all_metrics_df) - 1, "precision"]

	##### Split into ROC and PRC dataframes using base R
	model_roc <- all_metrics_df[, c("model", "tpr", "fpr")]
	model_roc <- model_roc[order(model_roc$fpr), ]

	model_prc <- all_metrics_df[, c("model", "precision", "recall")]
	model_prc <- model_prc[order(model_prc$recall), ]

	##### Collect results to a list
	list(
		all_metrics = all_metrics_df,
		roc = model_roc,
		prc = model_prc
	)
}


#### Function to calculate AUC using the trapezoidal rule
calculate.auc <- function(df, x_col, y_col) {
	# Sort data by x column
	df <- df[order(df[[x_col]]), ]
	
	# Calculate changes in x and average y values
	dx <- diff(df[[x_col]])
	y_vals <- df[[y_col]]
	avg_y <- (y_vals[-1] + y_vals[-length(y_vals)]) / 2
	
	# Calculate trapezoid areas and sum
	areas <- dx * avg_y
	auc <- sum(areas)
	
	return(auc)
}


## Function to create multimodel AUROC and AUPRC table based on a dataframe containing ROC and PRC
## Column names for ROC df must be fpr and tpr
## Column names for ROC df must be precision and recall
get.multimodel.auroc.auprc <- function(multimodel_roc_df,  multimodel_prc_df, sample_name, model_names = NULL, model_col = "model") {

	if(is.null(model_names)) {
		model_names <- unique(multimodel_roc_df$model)
	}

	auc_list <- list()

	for (i in seq_along(model_names)) {
		model_name <- model_names[i]

        # Filter using base R subsetting instead of dplyr::filter
		model_roc <- multimodel_roc_df[multimodel_roc_df[[model_col]] == model_name, ]
		model_prc <- multimodel_prc_df[multimodel_prc_df[[model_col]] == model_name, ]

		model_auroc <- calculate.auc(model_roc, x_col = "fpr", y_col = "tpr")
		model_auprc <- calculate.auc(model_prc, x_col = "recall", y_col = "precision")

		auc_list[[i]] <- data.frame(
			sample_id = sample_name,
			model = model_name,
			auroc = model_auroc,
			auprc = model_auprc
		)

	}

    # Combine list of data frames using base R
	aucs <- do.call(rbind, auc_list)

	return(aucs)
}

