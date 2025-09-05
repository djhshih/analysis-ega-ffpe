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


# Create variant ID based on chrom, pos, ref, and alt
# @param d  data.frame of variants
add_id <- function(d) {
	d$snv <- with(d, paste(chrom, pos, ref, alt, sep="_"))
	d
}


# Annotate called variant with whether each variant is in the ground truth
# variant call set
# @param  d  data.frame of called variants to be annotated
# @param  truth  data.frame of ground truth variants
annotate_truth <- function(d, truth) {
	d$truth <- d$snv %in% truth$snv;
	d
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

