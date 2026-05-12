## This script is used to make ROC and PRC plots from the evaluation results aggregated using combine_results.R
source("../common-ffpe-snvf/R/plot.R")

dset_dirs <- c(
	"EGAD00001004066/somatic_filtered-exome-blacklist-micr1234",
	"EGAD00001004066/somatic_filtered-exome-blacklist",
	"EGAD00001004066/somatic_filtered-exome",
	"EGAD00001004066/somatic_filtered"
)

models <- c(
	"mobsnvf" = "MOBSNVF",
	"vafsnvf" = "VAFSNVF",
	"gatk-obmm" = "GATK-OBMM",
	"sobdetector" = "SOBDetector",
	"microsec" = "MicroSEC",
	"ideafix" = "Ideafix",
	"ffpolish" = "FFPolish"
)

for (dset_dir in dset_dirs){

	message(sprintf("Processing: %s", dset_dir))

	## Set output directory
	outdir_root <- file.path(dset_dir, "plots")

	## Make a vector with paths to all the precrec eval objects
	roc_coord_paths <- sort(list.files(file.path(dset_dir, "roc-prc-auc/precrec"), pattern = "all-models_roc_coordinates.tsv", recursive = TRUE, full.names = TRUE))
	prc_coord_paths <- sort(list.files(file.path(dset_dir, "roc-prc-auc/precrec"), pattern = "all-models_prc_coordinates.tsv", recursive = TRUE, full.names = TRUE))

	stopifnot(length(roc_coord_paths) == length(prc_coord_paths))
	if (length(roc_coord_paths) == 0) {
		stop("No eval outputs found for: ", dset_dir)
	}

	## Create plots for each of the evaluated samples
	message("Creating ROC PRC plot for:")

	outdir <- file.path(outdir_root, "roc_prc_plots")
	dir.create(file.path(outdir, "roc"), recursive = TRUE, showWarnings = FALSE)
	dir.create(file.path(outdir, "prc"), recursive = TRUE, showWarnings = FALSE)

	# dset_key   <- sub("/.*", "", dset_dir)
	# plot_title <- dset_author[[dset_key]]

	for (i in seq_along(roc_coord_paths)) {
		sample_name <- match_return_sample_name(roc_coord_paths[i], prc_coord_paths[i])
		baseline_precision <- get_baseline_precision(dset_dir, sample_name)
		message(sprintf("\t %d. %s", i, sample_name))

		roc_coord <- read_delim_arrow(roc_coord_paths[i], delim = "\t")
		prc_coord <- read_delim_arrow(prc_coord_paths[i], delim = "\t")
		roc_coord$model <- unname(ifelse(roc_coord$model %in% names(models), models[roc_coord$model], roc_coord$model))
		prc_coord$model <- unname(ifelse(prc_coord$model %in% names(models), models[prc_coord$model], prc_coord$model))

		plots <- make_roc_prc_plot(
			roc_coord, 
			prc_coord,
			baseline_precision = baseline_precision
			# title = plot_title,
			# subtitle = sample_name
		)

		qdraw(plots$roc_prc, file.path(outdir, paste0(sample_name, "_roc_prc_plot.pdf")), width = 13, height = 7)
		qdraw(plots$roc, file.path(outdir, "roc", paste0(sample_name, "_roc_plot.pdf")), width = 7, height = 7)
		qdraw(plots$prc, file.path(outdir, "prc", paste0(sample_name, "_prc_plot.pdf")), width = 7, height = 7)
	}
}
