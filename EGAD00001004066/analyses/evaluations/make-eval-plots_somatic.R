#!/usr/bin/env Rscript

# Make evaluation plots for FFPE SNVF results

library(io)
library(precrec)
library(jsonlite)
library(argparser)
library(tidyverse)
library(glue)
library(patchwork)
library(grid)
library(hrbrthemes)
library(viridis)


main.outdir <- "../../evaluations/somatic_vcf"
ffpe_snvf.dir <- "../../ffpe-snvf/somatic_vcf"
vcf.dir <- "../../data/somatic_vcf"

source("../../../R/eval.R")
source("../../../R/plot.R")

## A lookup table to provide information about each sample and allow identification of matched fresh frozen samples
lookup_table <- read.delim("../../annot/sample_annotations.tsv") |> mutate(sample_name = glue('{gsub(" ", "-", title)}_{sample_alias}'))

## This is a vector of Standard Chromosomes to be used for filtering if the samples contains decoy variants
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))

## Initialize these dataframes
all.samples.all.models.scores.labels <- data.frame()
all.sample.summary <- data.frame()

## Obtain the FFPE tumoral and FF tumoral samples from the sample annotations
ffpe_tumoral <- lookup_table |> filter(preservation == "FFPE", sample_type == "Tumoral")
frozen_tumoral <- lookup_table |> filter(preservation == "Frozen", sample_type == "Tumoral")

## Evaluate each FFPE samples individually
for (i in seq_len(dim(ffpe_tumoral)[1])) {

	sample_name <- ffpe_tumoral[i, "sample_name"]
	tissue <- ffpe_tumoral[i, "tissue_type"]
	variant_caller <- "MuTect2"

	out_dir <- glue("{main.outdir}/roc-prc-plots/{sample_name}/")
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

	# This dataset has matched FFPE and FF from from the same patient.
	# Ground truth is constructed from the union of all the FF samples within the same tissue
	frozen_tumoral |> filter(tissue_type == tissue)
	truth_samples <- frozen_tumoral |> filter(tissue_type == tissue) |> select(sample_name) |> pull()
    truths <- snv_union(truth_samples) |> add_id()
	
	
	## Read in each models outputs and add id (CHROM_POS_REF_ALT)
	mobsnvf <- read.delim(glue("{ffpe_snvf.dir}/mobsnvf/{sample_name}/{sample_name}.mobsnvf.snv")) |> add_id()
	
	vafsnvf <- read.delim(glue("{ffpe_snvf.dir}/vafsnvf/{sample_name}/{sample_name}.vafsnvf.snv")) |> add_id()

	sobdetector <- read.delim(glue("{ffpe_snvf.dir}/sobdetector/{sample_name}/{sample_name}.sobdetector.snv")) |> 
		add_id() |>
		drop_na(artiStatus) |> 
		mutate(SOB = as.numeric(SOB)) |> 
		mutate(SOB = if_else(is.nan(SOB), 0, SOB))

	gatk.obmm <- read.delim(glue("{vcf.dir}/{sample_name}/{sample_name}.vcf.gz"), comment.char = "#", header = FALSE) |> 
		select(V1, V2, V4, V5, V7) |> 
		rename(chrom = V1, pos = V2, ref = V4, alt = V5, filter = V7) |>
		add_id() |>
		mutate(obmm = if_else(str_detect(filter, "orientation"), 0, 1))

	
	## MOBSNVF was run with damage type being FFPE. So it only evaluates C:G>T:A muattions.
	## The other mutations are given a score of NA.
	## Therefore the drop_na() functions drops all non C:G>T:A mutations.
	## The variant set is kept consistent by joining scores of other models to the mobsnvf variant set.
	all.models.scores.labels <- mobsnvf |>
		filter(chrom %in% standard_chromosomes) |>
		select(chrom, pos, ref, alt, snv, FOBP) |>
		drop_na() |>
		left_join(select(vafsnvf, snv, VAFF), by="snv") |>
		left_join(select(sobdetector, snv, SOB), by="snv") |>
		left_join(select(gatk.obmm, snv, obmm), by="snv") |>
		add_id()


	# Higher score is signifies real mutation : VAFSNVF
    # Lower score signifies real mutation:  MOBSNVF, SOBDetector
	# We flip scores for MOBSNVF and SOBDetector to make higher score signify a real muatation
	all.models.scores.labels <- all.models.scores.labels |>
			mutate(FOBP = -FOBP) |>
			mutate(SOB = -abs(SOB)) |>
			annotate_truth(truths)

	
	
	if (dim(filter(all.models.scores.labels, truth))[1] == 0) {
		print(glue("\n{sample_name} could not be evaluated as it contains 0 positive entries for truth\n"))
		next
	}

	# Create the precrec evaluation object
	all.model.eval <- make.multimodel.eval.object(all.models.scores.labels, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))

	# Get AUCs to include in plot
	all.model.aucs <- auc(all.model.eval)

	# Make AUCROC and AUPRC texts to include in the plots
	plot.auc.text <- make.plot.auc.text(all.model.eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object

	## Make ROC and PRC plot objects
	all.model.roc.plot <- autoplot(all.model.eval, "ROC")
	all.model.prc.plot <- autoplot(all.model.eval, "PRC")

	title <- "FFPE SNVF Evaluation"
	caption <- glue("Only C>T SNVs Evaluated \nC>T SNV count: {dim(all.models.scores.labels)[1]} \nVariants Called using {variant_caller} with matched normals")
	subtitle <-  glue("EGAD00001004066 Dataset\nSample Name: {sample_name}")

	## Combine the plots and text
	all.model.roc.prc.plot <- make.roc.prc.plot(
		all.model.roc.plot, 
		all.model.prc.plot, 
		plot.auc.text, 
		title = title,
		subtitle = subtitle,
		caption = caption
	)

	## Print log
	stdout <- glue("\n\n{i}. Making ROC and PRC plots for:
		FFPE Sample: {sample_name}
		Frozen/Truth Sample:  paste(truth_samples, collapse = ', ')
		SNVs Evaluated: {dim(all.models.scores.labels)[1]}
		truth_n_real : {dim(filter(all.models.scores.labels, truth))[1]}
		truth_n_artifacts : {dim(filter(all.models.scores.labels, !truth))[1]}
		\n")
		
	print(stdout)

	# Mark the cutoff position for different values of fp-cut
	mark_mobsnvf_cutoff(all.models.scores.labels, out_dir)


	## Collect Individual Eval Metrics
	## This is used later for making box plots
    p.vafsnvf <- with(all.models.scores.labels, evalmod(scores = VAFF, labels = truth))
    p.mobsnvf <- with(all.models.scores.labels, evalmod(scores = FOBP, labels = truth))
	p.sobdetector <- with(all.models.scores.labels, evalmod(scores = SOB, labels = truth))
	p.gatk.obmm <- with(all.models.scores.labels, evalmod(scores = obmm, labels = truth))

    res <- list(
			mobsnvf = p.mobsnvf,
			vafsnvf = p.vafsnvf,
			sobdetector = p.sobdetector,
			gatk = p.gatk.obmm,
			all_models = all.model.eval
		)
	

	## Save plot and res
    qdraw(all.model.roc.prc.plot, glue("{out_dir}/{sample_name}_roc_prc.pdf"), width = 7, height = 5)
	qwrite(res, glue("{out_dir}/{sample_name}_eval_res.rds"))
	qwrite(as.data.frame(all.model.eval), glue("{out_dir}/{sample_name}_eval_res.tsv"))
	qwrite(all.models.scores.labels, glue("{out_dir}/{sample_name}_dataset_labeled.tsv"))

	## Create a dataframe summarizing the summary metrics
	summary <- data.frame(
		ffpe_sample_name = sample_name,
		frozen_sample_names = paste(truth_samples, collapse = ', '),
		tissue_type = tissue,
		snvs_evaluated = dim(all.models.scores.labels)[1],
		n_real = dim(filter(all.models.scores.labels, truth))[1],
		n_artifacts = dim(filter(all.models.scores.labels, !truth))[1]
	)

	## Collect the per sample summary into a larger dataframe
	all.sample.summary <- rbind(all.sample.summary, summary)

	## Collect the scores and truth labels for each sample into a larger dataframe, annotating the sample name and workflow type
	all.models.scores.labels <- all.models.scores.labels |> mutate(workflow_type = "mutect2_matched_normal", sample_name = sample_name)
	all.samples.all.models.scores.labels <- rbind(all.samples.all.models.scores.labels, all.models.scores.labels)

}

print(glue("Done\n"))

print(glue("Finished making per sample ROC/PRC plots. Directory: {main.outdir}/roc-prc-plots/\n"))

qwrite(all.sample.summary, glue("{main.outdir}/all_samples-summary.tsv"))
print(glue("Saved summary to: {main.outdir}/all_samples-summary.tsv"))

qwrite(all.samples.all.models.scores.labels, glue("{main.outdir}/all_samples_all_models_scores_labels.tsv"))
print(glue("Saved scores and truth lables of all models to: {main.outdir}/all_samples_all_models_scores_labels.tsv"))

# For debugging in notebook
# options(repr.plot.width = 7, repr.plot.height = 5)
# all.model.roc.prc.plot

# options(repr.plot.width = 7, repr.plot.height = 5)

eval_df <- all.samples.all.models.scores.labels
all.samples.eval <- make.multimodel.eval.object(eval_df, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))
plot.auc.text <- make.plot.auc.text(all.samples.eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object
all.model.roc.plot <- autoplot(all.samples.eval, "ROC")
all.model.prc.plot <- autoplot(all.samples.eval, "PRC")


all.samples.roc.prc.plot <- make.roc.prc.plot(
	all.model.roc.plot, 
	all.model.prc.plot, 
	plot.auc.text, 
	subtitle = glue("EGAD00001004066 Dataset, 13 Samples (Liver + Colon)"),
	caption = glue("Only C>T SNVs Evaluated \nC>T SNV count: {dim(eval_df)[1]}\nMicroSEC is not evaluated")
)

# all.samples.roc.prc.plot
outpath <- glue("{main.outdir}/roc-prc-plots/all_samples_roc_prc_plot.pdf")
print(glue("Saved ROC and PRC plots across all samples to {outpath}\n"))
qdraw(all.samples.roc.prc.plot, outpath, width = 7, height = 5)


## Colon Tissues only

eval_df <- all.samples.all.models.scores.labels |> filter(str_detect(sample_name, "Colon"))
all.samples.eval <- make.multimodel.eval.object(eval_df, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))
plot.auc.text <- make.plot.auc.text(all.samples.eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object
all.model.roc.plot <- autoplot(all.samples.eval, "ROC")
all.model.prc.plot <- autoplot(all.samples.eval, "PRC")

all.samples.roc.prc.plot <- make.roc.prc.plot(
	all.model.roc.plot, 
	all.model.prc.plot, 
	plot.auc.text, 
	subtitle = glue("EGAD00001004066 Dataset, 6 colon Samples"),
	caption = glue("Only C>T SNVs Evaluated \nC>T SNV count: {dim(eval_df)[1]}\nMicroSEC is not evaluated")
)

# all.samples.roc.prc.plot
outpath <- glue("{main.outdir}/roc-prc-plots/colon_samples_roc_prc_plot.pdf")
print(glue("Saved ROC and PRC plots across all samples to {outpath}\n"))
qdraw(all.samples.roc.prc.plot, outpath, width = 7, height = 5)


## Liver Tissues only
eval_df <- all.samples.all.models.scores.labels |> filter(str_detect(sample_name, "Liver"))
all.samples.eval <- make.multimodel.eval.object(eval_df, score_columns = c("FOBP", "VAFF", "SOB", "obmm"), names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))
plot.auc.text <- make.plot.auc.text(all.samples.eval, model.names = c("mobsnvf", "vafsnvf", "sobdetector", "gatk-obmm"))$text.plot.object
all.model.roc.plot <- autoplot(all.samples.eval, "ROC")
all.model.prc.plot <- autoplot(all.samples.eval, "PRC")

all.samples.roc.prc.plot <- make.roc.prc.plot(
	all.model.roc.plot, 
	all.model.prc.plot, 
	plot.auc.text, 
	subtitle = glue("EGAD00001004066 Dataset, 7 Liver Samples"),
	caption = glue("Only C>T SNVs Evaluated \nC>T SNV count: {dim(eval_df)[1]}\nMicroSEC is not evaluated")
)

# all.samples.roc.prc.plot
outpath <- glue("{main.outdir}/roc-prc-plots/liver_samples_roc_prc_plot.pdf")
print(glue("Saved ROC and PRC plots across all samples to {outpath}\n"))
qdraw(all.samples.roc.prc.plot, outpath, width = 7, height = 5)


## Obtain AUROC and AUPRC from the evaluation results saved previously
eval_res_path <- list.files(glue("{main.outdir}/roc-prc-plots"), pattern = "_eval_res\\.rds", recursive = TRUE, full.names = TRUE)

## Append AUROC and AUPRC for each model to a dataframe
model.aucs <- data.frame()

for (path in eval_res_path) {
	parts <- stringr::str_split_1(path, "/")
	sample_id <- sub("_eval_res.rds", "", parts[6])
	variant_caller <- "MuTect2"

	eval_res <- readRDS(path)
	models <- names(eval_res)
	
	for (model in models[models != "all_models"]) {
		aucs <- auc(eval_res[[model]])
		auroc <- aucs |> filter(curvetypes == "ROC") |> pull(aucs)
		auprc <- aucs |> filter(curvetypes == "PRC") |> pull(aucs)
		model.aucs <- bind_rows(
			model.aucs,
			data.frame(
				sample_id = sample_id,
				variant_caller = variant_caller,
				model = model,
				auroc = auroc,
				auprc = auprc,
				stringsAsFactors = FALSE
			)
		)
	}
}

## Save the AUCs
write.table(model.aucs, glue("{main.outdir}/aucs.tsv"), sep = "\t", row.names = FALSE)
print(glue("Saved AUCs of each sample at {main.outdir}/aucs.tsv"))


## Find VCFs where mobsnvf performs worse than other models
print(glue("Determining where MOBSNVF performs worse than other models"))
vcfs.mobsnvf.lower.auroc <- model.aucs |>
  group_by(sample_id) |>
  # Summarize mobsnvf AUROC and the max AUROC of the others
  summarize(
    mob_auroc   = auroc[model == "mobsnvf"],
    others_max  = max(auroc[model != "mobsnvf"])
  ) |> filter(mob_auroc < others_max) |>
  pull(sample_id)

lower.mobsnvf.auroc <- model.aucs |> filter(sample_id %in% vcfs.mobsnvf.lower.auroc)
# lower.mobsnvf.auroc

vcfs.mobsnvf.lower.auprc <- model.aucs |>
  group_by(sample_id) |>
  # compute mobsnvf auprc and the min auprc of the others
  summarize(
    mob_auprc   = auprc[model == "mobsnvf"],
    others_max  = max(auprc[model != "mobsnvf"])
  ) |> filter(mob_auprc < others_max) |>
  pull(sample_id)

lower.mobsnvf.auprc <- model.aucs |> filter(sample_id %in% vcfs.mobsnvf.lower.auprc)
# lower.mobsnvf.auprc

print(glue("\tSaving the AUROC and AUPRC for these samples to a table"))


## Save the AUCs of these specific samples to examine later
write.table(lower.mobsnvf.auroc, glue("{main.outdir}/mobsnvf_lower_auroc.tsv"), sep = "\t", row.names = FALSE)
write.table(lower.mobsnvf.auprc, glue("{main.outdir}/mobsnvf_lower_auprc.tsv"), sep = "\t", row.names = FALSE)
print(glue("\t\t- {main.outdir}/mobsnvf_lower_auroc.tsv"))
print(glue("\t\t- {main.outdir}/mobsnvf_lower_auprc.tsv\n"))



## Make boxplots
scale = 0.4; w = 14*scale; h = 12*scale
options(repr.plot.width = w, repr.plot.height = h)

print(glue("Making AUROC and AUPRC boxplots"))
# Call the plotting function
model.aucs.mutated <- model.aucs |> filter(model != "microsec") |> mutate(model = ifelse(model == "gatk", "gatk-obmm", model))

all_auc_boxplots <- make_all_auc_boxplots(
		model.aucs.mutated, 
		variant_callers = c("MuTect2"), 
		subtitle = glue("EGAD00001004066 Dataset\nMuTect2 with Matched Normal\n13 Sample (Liver + Colon)"),
		scale = 0.4, 
		text_scale = 0.5, 
		w = w, 
		h = h
	)

# Set outdir for saving plots
outdir = glue("{main.outdir}/box-plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#### Save AUROC plots
qdraw(all_auc_boxplots$all$auroc$MuTect2, file = glue("{outdir}/auroc_mutect2_boxplot.pdf"), width = w, height = h)
qdraw(all_auc_boxplots$Liver$auroc$MuTect2, file = glue("{outdir}/auroc_liver_mutect2_boxplot.pdf"), width = w, height = h)
qdraw(all_auc_boxplots$Colon$auroc$MuTect2, file = glue("{outdir}/auroc_colon_mutect2_boxplot.pdf"), width = w, height = h)


#### Save AUPRC plots
qdraw(all_auc_boxplots$all$auprc$MuTect2, file = glue("{outdir}/auprc_mutect2_boxplot.pdf"), width = w, height = h)
qdraw(all_auc_boxplots$Liver$auprc$MuTect2, file = glue("{outdir}/auprc_liver_mutect2_boxplot.pdf"), width = w, height = h)
qdraw(all_auc_boxplots$Colon$auprc$MuTect2, file = glue("{outdir}/auprc_colon_mutect2_boxplot.pdf"), width = w, height = h)

print(glue("Done.\nBox-Plots saved to {outdir}/"))



