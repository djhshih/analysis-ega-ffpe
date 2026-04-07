#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm

repo_root = ".."

## Functions
def get_filtering_summary(
    sample_name: str, 
    model: str, 
    snvf_var_set: str,
	msec_var_set: str,
	msec_filter_col: str,
    filtered_var_set: str, 
    source_variants: pl.DataFrame, 
    msec_artifacts: pl.DataFrame, 
    selected_snvf: pl.DataFrame
) -> dict:

    n_orig = source_variants.height
    n_msec_arti = msec_artifacts.height
    n_selected = selected_snvf.height

    summary = {
        "source_var_set" : snvf_var_set,
		"microsec_var_set": msec_var_set,
        "filtered_var_set" : filtered_var_set,
		"microsec_filter_col" : msec_filter_col,
        "sample_name" : sample_name,
        "model" : model,
        "n_var_source_snvf" : n_orig,
        "n_microsec_artifacts" : n_msec_arti,
        "n_var_selected": n_selected,
        "n_var_removed" : n_orig - n_selected,
        "pct_removed" : ((n_orig - n_selected) / n_orig) * 100
    }

    return summary


def get_ffpe_snvf_paths(dataset: str, variant_set: str) -> list:
	"""
	Returns the path for each FFPE SNVF model's
	results for a specified variant set and dataset
	"""

	paths = sorted(
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/*/*/*.mobsnvf.snv") +
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/*/*/*.vafsnvf.snv") +
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/*/*/*.sobdetector.snv") +
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/*/*/*.ideafix.tsv") +
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/*/*/*.gatk-obmm.tsv") +
		glob.glob(f"{repo_root}/ffpe-snvf/{dataset}/{variant_set}/*/*/*.ffpolish.tsv")
	)

	return paths


def filter_dataset(
		dataset: str, 
		source_variant_set: str,
		msec_variant_set: str,
		new_variant_set: str,
		mesc_filter_col: str = "msec_filter_1234"
	) -> None:

	print(f"Processing Dataset: {dataset} | Variant Set: {source_variant_set} | MicroSEC variant Set: {msec_variant_set} | Output Set: {new_variant_set}")

	vcf_dir = f"{repo_root}/vcf/{dataset}/{new_variant_set}"
	msec_dir = f"{repo_root}/microsec/{dataset}/{msec_variant_set}"

	ffpe_snvf_paths = get_ffpe_snvf_paths(dataset, source_variant_set)
	filtering_summary = []

	for path in tqdm(ffpe_snvf_paths):
		model = path.split("/")[-3]
		sample_name = path.split("/")[-2]
		fname = os.path.basename(path)

		## Read FFPE filtering model's results
		snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)
		
		## Read MicroSEC results and select artifacts
		msec_res_path = f"{msec_dir}/{sample_name}/{sample_name}.microsec.tsv"
		msec_res = pl.read_csv(msec_res_path, separator="\t", infer_schema_length=10000)
		msec_artifacts = msec_res.filter(pl.col(mesc_filter_col).is_not_null())
		
		## Remove variants that is seen in MicroSEC artifacts table
		filtered_snvf = snvf.join(msec_artifacts, left_on = ["chrom", "pos", "ref", "alt"], right_on=["Chr", "Pos", "Ref", "Alt"], how="anti")

		## Write filterd variants to disk
		filtered_snvf_outdir = f"{repo_root}/ffpe-snvf/{dataset}/{new_variant_set}/{model}/{sample_name}"
		os.makedirs(filtered_snvf_outdir, exist_ok=True)

		filtered_snvf.write_csv(f"{filtered_snvf_outdir}/{fname}", separator="\t")

		## Generate Filtering Summary
		sample_filtering_summary = get_filtering_summary(
			sample_name=sample_name,
			model=model,
			snvf_var_set=source_variant_set,
			msec_var_set=msec_variant_set,
			msec_filter_col=mesc_filter_col,
			filtered_var_set=new_variant_set,
			source_variants=snvf,
			msec_artifacts=msec_artifacts,
			selected_snvf=filtered_snvf
		)

		filtering_summary.append(sample_filtering_summary)

	pl.DataFrame(filtering_summary).write_csv(f"{dataset}/{new_variant_set}/{os.path.basename(__file__).split(".")[0]}_filtering-summary.tsv", separator="\t")


# ## SNVF MICR Filtering
filter_dataset(
		dataset = "EGAD00001004066",
		source_variant_set = "somatic_filtered-exome-blacklist",
		msec_variant_set = "somatic_filtered-exome",
		new_variant_set = "somatic_filtered-exome-blacklist-micr1234",
		mesc_filter_col = "msec_filter_1234"
	)

