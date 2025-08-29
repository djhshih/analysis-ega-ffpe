#!/usr/bin/env Rscript

# This script contains functions for evaluation
## This is a work in progress. The functions still needs to be improved for generalizability across multiple projects and datasets

### Evaluations

#### Create unique identifier for each variant
add_id <- function(ffpe_variants) {
	ffpe_variants$snv <- paste(ffpe_variants$chrom, ffpe_variants$pos, ffpe_variants$ref, ffpe_variants$alt, sep="_")
	ffpe_variants
}

#### Annotate the presence of a variant in the truth set i.e in the Fresh Frozen Sample
annotate_truth <- function(ffpe_variants, frozen_variants) {
	ffpe_variants$truth <- ffpe_variants$snv %in% frozen_variants$snv;
	ffpe_variants
}

#### Function to filter C>T and G>A mutations
ct_only <- function(df) {
  df[(df$ref == "C" & df$alt == "T") | (df$ref == "G" & df$alt == "A"), ]
}

#### Create a set of SNVs across multiple samples
snv_union <- function(truth_samples) {

	truths <- data.frame()

	for (truth_sample_name in truth_samples){
		truth_sample <- qread(glue("{vcf.dir}/{truth_sample_name}/{truth_sample_name}.vcf.gz"), type = "vcf") |> select(chrom, pos, ref, alt)
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
		# select the top significant results
		pred[idx[1:top]] <- TRUE;
		pred
	} else {
		# none passes
		rep(FALSE, n)
	}
}

#### Function to label the models prediction based on the adaptive_fdr_cut function
fdr_cut_pred <- function(df, score, fp.cut=0.5) {
	
    df.ct <- df |> filter(complete.cases(.data[[score]]))
    df.nct <- df |> filter(!complete.cases(.data[[score]])) |> mutate(score = NA, q = NA, pred = TRUE)

    df.ct.q <- df.ct |>
        mutate(
            score = ifelse(.data[[score]] == 0, .Machine$double.eps, .data[[score]]),
            q = p.adjust(score, "BH"),
			## TRUE being real muatation and FALSE being artifacts
            pred = ifelse(adaptive_fdr_cut(q, fp.cut), TRUE, FALSE)
        )

    rbind(df.ct.q, df.nct) |>
		arrange(chrom, pos)
}

#### Calculate Evaluation metrics based on the models prediction and the ground truth label. Both should be boolean
calc.eval.metrics <- function(df, pred_col = "pred", truth_col = "truth"){

	TP <- df |> filter((.data[[pred_col]]) & (.data[[truth_col]])) |> count() |> pull()
	FP <- df |> filter((.data[[pred_col]]) & !(.data[[truth_col]])) |> count() |> pull()
	FN <- df |> filter(!(.data[[pred_col]]) & (.data[[truth_col]])) |> count() |> pull()
	TN <- df |> filter(!(.data[[pred_col]]) & !(.data[[truth_col]])) |> count() |> pull()

	precision <- TP/(TP + FP)
	recall <- TP/(TP + FN)
	tpr <- recall
	fpr <- FP/(FP + TN)

	list(
		precision = precision,
		recall = recall,
		tp_rate = tpr,
		fp_rate = fpr
	)
}

#### Creates a multimodel evaluation object using the precrec library
make.multimodel.eval.object <- function(scores_labels_df, score_columns = c("FOBP", "VAFF", "SOB", "msec"), names = c("mobsnvf", "vafsnvf", "sobdetector", "microsec")) {
	all.model.scores <- do.call(join_scores, scores_labels_df[score_columns])
	all.model.labels <- do.call(join_labels, replicate(length(score_columns), scores_labels_df$truth, simplify = FALSE))
	all.model.names <- names

	evalmod(mmdata(all.model.scores, all.model.labels, all.model.names))
}

#### Creates a text panel containing all the AUC metrics for each model
make.plot.auc.text <- function(multi.model.eval.object, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "microsec")) {
	
	# Get AUCs to include in plot
	all.model.aucs <- auc(multi.model.eval.object)
	
	# Extract AUROC and AUPRC for each model in model.names
	auroc <- sapply(model.names, function(m) {
		all.model.aucs |> filter(modnames == m & curvetypes == "ROC") |> pull(aucs)
	})
	auprc <- sapply(model.names, function(m) {
		all.model.aucs |> filter(modnames == m & curvetypes == "PRC") |> pull(aucs)
	})

	# Dynamically build AUROC and AUPRC text lines for each model
	auroc_lines <- paste0(model.names, "=", round(auroc[model.names], 3))
	auprc_lines <- paste0(model.names, "=", round(auprc[model.names], 3))

	# Make AUCROC and AUPRC texts to include in the plots
	auc_text <- glue(
		"\nAUROC: \n{paste(auroc_lines, collapse = '\n')} \n\nAUPRC: \n{paste(auprc_lines, collapse = '\n')}"
	)

	# Make plot text
	auc_grob <- textGrob(
		auc_text,
		x = 0, y = 1, just = c("left", "top"),
		gp = grid::gpar(fontsize = 8, fontfamily = "mono")
	)

	list(
		text = auc_text,
		text.plot.object = auc_grob
	)
}