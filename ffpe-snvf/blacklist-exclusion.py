#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm 

# %% [markdown]
# ## Functions

def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT"]) -> pl.DataFrame:
	"""
	Reads Variants from a VCF file into a Polars DataFrame.
	By default, only reads essential columns.
	Splits multi-allelic variants into separate rows.
	Can be extended to read additional columns as needed or
	read custom variant files with similar structure.
	
	:param path: Path to VCF file
	:type path: str
	:param columns: List of column names to read from the VCF file. Defaults to ["#CHROM", "POS", "REF", "ALT"].
	:type columns: list
	:return: DataFrame containing the variants
	:rtype: DataFrame
	"""
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants

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
        "n_var_target_vcf" : n_target,
        "n_var_filtered_snvf": n_filtered,
        "n_var_removed" : n_orig - n_filtered,
        "pct_removed" : ((n_orig - n_filtered) / n_orig) * 100
    }

    return summary

# %% [markdown]
# ## Setup

repo_root = ".."


# %% [markdown]
# ## SNVF Blacklist Filtering

## Wrapper Function
def filter_dataset(dataset: str, source_var_set: str, new_var_set: str) -> None:

	vcf_dir = f"{repo_root}/vcf/{dataset}/{new_var_set}"


	ffpe_snvf = (
		glob.glob(f"{dataset}/{source_var_set}/*/*/*.ffpolish.tsv") +
		glob.glob(f"{dataset}/{source_var_set}/*/*/*.gatk-obmm.tsv") +
		glob.glob(f"{dataset}/{source_var_set}/*/*/*.ideafix.tsv") +
		glob.glob(f"{dataset}/{source_var_set}/*/*/*.mobsnvf.snv") +
		glob.glob(f"{dataset}/{source_var_set}/*/*/*.sobdetector.snv") +
		glob.glob(f"{dataset}/{source_var_set}/*/*/*.vafsnvf.snv")
	)

	filtering_summary = []

	for path in tqdm(ffpe_snvf):
		model = path.split("/")[-3]
		sample_name = path.split("/")[-2]
		fname = os.path.basename(path)

		snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)
		
		target_vars_path = f"{vcf_dir}/{sample_name}/{sample_name}.vcf.gz"
		target_vars = read_variants(target_vars_path)
		

		filtered_snvf = snvf.join(target_vars, on = ["chrom", "pos", "ref", "alt"], how="semi")

		filtered_snvf_outdir = f"{dataset}/{new_var_set}/{model}/{sample_name}"
		os.makedirs(filtered_snvf_outdir, exist_ok=True)


		sample_filtering_summary = get_filtering_summary(
			sample_name,
			model,
			source_var_set,
			new_var_set,
			snvf,
			target_vars,
			filtered_snvf
		)

		filtering_summary.append(sample_filtering_summary)

		filtered_snvf.write_csv(f"{filtered_snvf_outdir}/{fname}", separator="\t")

	pl.DataFrame(filtering_summary).write_csv(f"{dataset}/{new_var_set}/blacklist-exclusion_filtering-summary.tsv", separator="\t")



# %% [markdown]
# ### Filter specified variant set from each dataset

filter_dataset("EGAD00001004066", "somatic_vcf-dp10", "somatic_filtered-dp10-blacklist")


