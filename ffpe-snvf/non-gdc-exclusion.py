#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm
import sys
from typing import Optional

# Local Dependencies
repo_root = ".."
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from common import read_variants

def get_filtering_summary(
    sample_name: str, 
    model: str, 
    orig_var_set: str, 
    new_var_set: str, 
    snvf: pl.DataFrame, 
    target_vars: pl.DataFrame, 
    filtered_snvf: pl.DataFrame
) -> dict:

    n_orig = snvf.height
    n_target = target_vars.height
    n_filtered = filtered_snvf.height

    summary = {
        "original_var_set" : orig_var_set,
        "new_var_set" : new_var_set,
        "sample_name" : sample_name,
        "model" : model,
        "n_var_original" : n_orig,
        "n_var_target_set" : n_target,
        "n_var_filtered_snvf": n_filtered,
        "n_var_removed" : n_orig - n_filtered,
        "pct_removed" : ((n_orig - n_filtered) / n_orig) * 100
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


## SNVF GDC Filtering
def filter_dataset(dataset: str, source_variant_set: str, new_variant_set: Optional[str] = None) -> None:

	if not new_variant_set:
		new_variant_set = f"{source_variant_set}-gdc"

	print(f"Processing Dataset: {dataset} | Variant Set: {source_variant_set}")

	vcf_dir = f"{repo_root}/vcf/{dataset}/{new_variant_set}"

	ffpe_snvf = get_ffpe_snvf_paths(dataset, source_variant_set)

	filtering_summary = []

	for path in tqdm(ffpe_snvf):
		model = path.split("/")[-3]
		sample_name = path.split("/")[-2]
		fname = os.path.basename(path)

		snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)
		
		target_vars_path = f"{vcf_dir}/{sample_name}/{sample_name}.vcf.gz"
		target_vars = read_variants(target_vars_path)

		filtered_snvf = snvf.join(target_vars, on = ["chrom", "pos", "ref", "alt"], how="semi")

		filtered_snvf_outdir = f"{dataset}/{new_variant_set}/{model}/{sample_name}"
		os.makedirs(filtered_snvf_outdir, exist_ok=True)

		sample_filtering_summary = get_filtering_summary(
			sample_name,
			model,
			source_variant_set,
			new_variant_set,
			snvf,
			target_vars,
			filtered_snvf
		)

		filtering_summary.append(sample_filtering_summary)
		filtered_snvf.write_csv(f"{filtered_snvf_outdir}/{fname}", separator="\t")

	pl.DataFrame(filtering_summary).write_csv(f"{dataset}/{new_variant_set}/gdc-exclusion_filtering-summary.tsv", separator="\t")


### Filter specified variant set from each dataset
filter_dataset(
	dataset = "EGAD00001004066", 
    source_variant_set = "somatic_filtered-dp20",
    new_variant_set = "somatic_filtered-dp20-blacklist-gdc"
)
