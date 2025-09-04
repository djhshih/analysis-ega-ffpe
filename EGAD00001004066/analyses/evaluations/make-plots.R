#!/usr/bin/env Rscript
source("../../../R/plot.R")
source("../../../R/eval.R")

## Set output directory
outdir_root = "../../results/somatic_vcf/plots"

## Make a vector with paths to all the precrec eval objects
precrec_eval_paths <- Sys.glob("../../results/somatic_vcf/roc-prc-auc/precrec/*/*.rds")

## Create plots for each of the evaluated samples
for (path in precrec_eval_paths) {

	sample_name <- str_split_1(path, "/")[7]

	precrec_eval <- readRDS(path)

	## Read in the variant set to count number and annotate the number of SNVs in the plot
	eval_snv_set <- qread(glue("../../results/somatic_vcf/model_scores_labels_truths/{sample_name}/{sample_name}_model_scores_labels_truths.tsv"))
	snv_count <- nrow(eval_snv_set)

	precrec_roc_plot <- autoplot(precrec_eval, "roc")
	precrec_prc_plot <- autoplot(precrec_eval, "prc")
	plot_auc_text <- make.plot.auc.text(precrec_eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object

	# options(repr.plot.width = 7, repr.plot.height = 5)
	precrec_roc_prc_plot <- make.roc.prc.plot(
		precrec_roc_plot, 
		precrec_prc_plot,
		plot_auc_text, 
		title = "FFPE SNVF Evaluation",
		subtitle = glue("EGAD00001004066 Dataset\n{sample_name}"),
		caption = glue("Only C>T SNVs Evaluated\nC>T SNV count: {snv_count}")
	)

	outdir = glue("{outdir_root}/roc_prc_plots/{sample_name}")
	dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

	qdraw(precrec_roc_prc_plot, glue("{outdir}/{sample_name}_roc_prc_plot.pdf"), width = 7, height = 5)

}

## Make an overall evaluation plot for all samples

precrec_all_samples_eval <- readRDS("../../results/somatic_vcf/roc-prc-auc/precrec/all_samples_precrec_eval.rds")

## Read in the variant set to count number and annotate the number of SNVs in the plot
eval_snv_set <- qread(glue("../../results/somatic_vcf/model_scores_labels_truths/all_samples_all_scores_labels_truths.tsv"))
snv_count <- nrow(eval_snv_set)

precrec_all_samples_roc_plot <- autoplot(precrec_all_samples_eval, "roc")
precrec_all_samples_prc_plot <- autoplot(precrec_all_samples_eval, "prc")
all_samples_plot_auc_text <- make.plot.auc.text(precrec_all_samples_eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object

precrec_roc_prc_plot <- make.roc.prc.plot(
	precrec_all_samples_roc_plot, 
	precrec_all_samples_prc_plot,
	all_samples_plot_auc_text, 
	title = "FFPE SNVF Evaluation",
	subtitle = glue("EGAD00001004066 Dataset\nAll 13 Sample (Liver + Colon)"),
	caption = glue("Only C>T SNVs Evaluated\nC>T SNV count: {snv_count}")
)

outdir = glue("{outdir_root}/roc_prc_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

qdraw(precrec_roc_prc_plot, glue("{outdir}/all_samples_roc_prc_plot.pdf"), width = 7, height = 5)

## Make an overall evaluation plot for all colon samples

precrec_all_samples_eval <- readRDS("../../results/somatic_vcf/roc-prc-auc/precrec/colon_samples_precrec_eval.rds")

## Read in the variant set to count number and annotate the number of SNVs in the plot
eval_snv_set <- qread(glue("../../results/somatic_vcf/model_scores_labels_truths/colon_samples_all_scores_labels_truths.tsv"))
snv_count <- nrow(eval_snv_set)

precrec_all_samples_roc_plot <- autoplot(precrec_all_samples_eval, "roc")
precrec_all_samples_prc_plot <- autoplot(precrec_all_samples_eval, "prc")
all_samples_plot_auc_text <- make.plot.auc.text(precrec_all_samples_eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object

precrec_roc_prc_plot <- make.roc.prc.plot(
	precrec_all_samples_roc_plot, 
	precrec_all_samples_prc_plot,
	all_samples_plot_auc_text, 
	title = "FFPE SNVF Evaluation",
	subtitle = glue("EGAD00001004066 Dataset\n6 Colon Samples"),
	caption = glue("Only C>T SNVs Evaluated\nC>T SNV count: {snv_count}")
)

outdir = glue("{outdir_root}/roc_prc_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

qdraw(precrec_roc_prc_plot, glue("{outdir}/colon_samples_roc_prc_plot.pdf"), width = 7, height = 5)

## Make an overall evaluation plot for all liver samples

precrec_all_samples_eval <- readRDS("../../results/somatic_vcf/roc-prc-auc/precrec/liver_samples_precrec_eval.rds")

## Read in the variant set to count number and annotate the number of SNVs in the plot
eval_snv_set <- qread(glue("../../results/somatic_vcf/model_scores_labels_truths/liver_samples_all_scores_labels_truths.tsv"))
snv_count <- nrow(eval_snv_set)

precrec_all_samples_roc_plot <- autoplot(precrec_all_samples_eval, "roc")
precrec_all_samples_prc_plot <- autoplot(precrec_all_samples_eval, "prc")
all_samples_plot_auc_text <- make.plot.auc.text(precrec_all_samples_eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object

precrec_roc_prc_plot <- make.roc.prc.plot(
	precrec_all_samples_roc_plot, 
	precrec_all_samples_prc_plot,
	all_samples_plot_auc_text, 
	title = "FFPE SNVF Evaluation",
	subtitle = glue("EGAD00001004066 Dataset\n7 Liver Samples"),
	caption = glue("Only C>T SNVs Evaluated\nC>T SNV count: {snv_count}")
)

outdir = glue("{outdir_root}/roc_prc_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

qdraw(precrec_roc_prc_plot, glue("{outdir}/liver_samples_roc_prc_plot.pdf"), width = 7, height = 5)


#---------------

# Make box plots

## Collect the AUCs across all the samples

### Create a vector with paths to all the AUC tables
auc_table_paths <- Sys.glob("../../results/somatic_vcf/roc-prc-auc/precrec/*/*_auc_table.tsv")

### Collect AUCs 
all_auc <- data.frame()

for (path in auc_table_paths){
	auc <- read.delim(path)
	all_auc <- rbind(all_auc, auc)
}

## Stratify by tissue type
colon_samples_auc <- all_auc[str_detect(all_auc$sample_id, "Colon"),]
liver_samples_auc <- all_auc[str_detect(all_auc$sample_id, "Liver"),]

## Create Box Plots for all samples
### AUROC
all_samples_auroc_boxplot <- make_auc_boxplot(
	all_auc, 
	"auroc", 
	title = "Area Under ROC Curve", 
	subtitle = "EGAD00001004066 Dataset\nAll 13 FFPE Tumoral Samples (Liver + Colon)\nMutect2 matched normal mode"
)

### AUPRC
all_samples_auprc_boxplot <- make_auc_boxplot(
	all_auc, 
	"auprc", 
	title = "Area Under Precision-Recall Curve", 
	subtitle = "EGAD00001004066 Dataset\nAll 13 FFPE Tumoral Samples (Liver + Colon)\nMutect2 matched normal mode"
)

## Create Box Plots for liver samples
### AUROC
colon_samples_auroc_boxplot <- make_auc_boxplot(
	colon_samples_auc, 
	"auroc", 
	title = "Area Under ROC Curve", 
	subtitle = "EGAD00001004066 Dataset\n6 FFPE Tumoral Colon Samples\nMutect2 matched normal mode"
)

### AUPRC
colon_samples_auprc_boxplot <- make_auc_boxplot(
	colon_samples_auc, 
	"auprc", 
	title = "Area Under Precision-Recall Curve", 
	subtitle = "EGAD00001004066 Dataset\n6 FFPE Tumoral Colon Samples\nMutect2 matched normal mode"
)

## Create Box Plots for liver samples
### AUROC
liver_samples_auroc_boxplot <- make_auc_boxplot(
	liver_samples_auc, 
	"auroc", 
	title = "Area Under ROC Curve", 
	subtitle = "EGAD00001004066 Dataset\n7 FFPE Tumoral Liver Samples\nMutect2 matched normal mode"
)

### AUPRC
liver_samples_auprc_boxplot <- make_auc_boxplot(
	liver_samples_auc, 
	"auprc", 
	title = "Area Under Precision-Recall Curve", 
	subtitle = "EGAD00001004066 Dataset\n7 FFPE Tumoral Liver Samples\nMutect2 matched normal mode"
)



## Assign output directory
out_dir = glue("{outdir_root}/box_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Save Box Plots to file
qdraw(all_samples_auroc_boxplot, glue("{out_dir}/all_samples_auroc_box_plot.pdf"), width = 5, height = 5)
qdraw(colon_samples_auroc_boxplot, glue("{out_dir}/liver_samples_auroc_box_plot.pdf"), width = 5, height = 5)
qdraw(liver_samples_auroc_boxplot, glue("{out_dir}/colon_samples_auroc_box_plot.pdf"), width = 5, height = 5)

qdraw(all_samples_auprc_boxplot, glue("{out_dir}/all_samples_auprc_box_plot.pdf"), width = 5, height = 5)
qdraw(colon_samples_auprc_boxplot, glue("{out_dir}/liver_samples_auprc_box_plot.pdf"), width = 5, height = 5)
qdraw(liver_samples_auprc_boxplot, glue("{out_dir}/colon_samples_auprc_box_plot.pdf"), width = 5, height = 5)


