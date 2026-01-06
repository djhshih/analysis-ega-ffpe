#!/usr/bin/env python
import polars as pl
import os
import glob
from tqdm import tqdm

def read_variants(path:str, columns: list = ["#CHROM", "POS", "REF", "ALT", "FILTER"]) -> pl.DataFrame:
	variants = (
		pl.read_csv(path, separator="\t", comment_prefix="##", infer_schema_length=1000, columns=columns)
		.rename(lambda x: x.lstrip("#").lower())
		.with_columns(pl.col("alt").str.split(","))
		.explode("alt")
	)
	return variants

dataset = "EGAD00001004066"
variant_type = "somatic_vcf"

vcf_paths = glob.glob(f"../vcf/{dataset}/{variant_type}/*/*.vcf.gz")

for path in tqdm(vcf_paths):
	sample_name = path.split("/")[-2]

	if "Frozen" in sample_name:
		continue

	variants = read_variants(path)

	outdir = f"{dataset}/{variant_type}/gatk-obmm/{sample_name}"
	os.makedirs(outdir, exist_ok=True)

	variants.write_csv(f"{outdir}/{sample_name}.gatk-obmm.tsv", separator="\t")


