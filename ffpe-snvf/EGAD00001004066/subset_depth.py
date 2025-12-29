#!/usr/bin/env python
import polars as pl
import glob
import os


# %%
paths = sorted(
	glob.glob("somatic_vcf/*/*/*.mobsnvf.snv") +
	glob.glob("somatic_vcf/*/*/*.vafsnvf.snv") +
	glob.glob("somatic_vcf/*/*/*.sobdetector.snv") +
	glob.glob("somatic_vcf/*/*/*.microsec.tsv") +
	glob.glob("somatic_vcf/*/*/*.ideafix.tsv") +
	glob.glob("somatic_vcf/*/*/*.gatk-obmm.tsv")
)

paths = [path for path in paths if "Frozen" not in path]


# %%
## Ensure VCF with this filter exists
filter_type = "somatic_vcf-dp10"
dataset = "EGAD00001004066"

for i, path in enumerate(paths):
	
	print(f"{i+1}. Processing file: {path}")
	

	model = os.path.basename(path).split(".")[-2]
  
	sample_name = os.path.basename(path).split(".")[0]

	## Read in the filtered VCF
	vcf_path = glob.glob(f"../../vcf/{dataset}/{filter_type}/{sample_name}/{sample_name}*.vcf.gz")[0]

	try:
		vcf_df = pl.read_csv(vcf_path, separator="\t", comment_prefix="##", null_values=".", columns=["#CHROM", "POS", "REF", "ALT"], infer_schema_length=1000).rename({"#CHROM": "CHROM"}).rename(lambda x: x.lower())
		vcf_df = vcf_df.with_columns(pl.col("alt").str.split(",")).explode("alt")
	except Exception as e:
		raise RuntimeError(f"Failed to read VCF file at {vcf_path}: {e}")

	snvf = pl.read_csv(path, separator="\t", infer_schema_length=1000)

	if model == "microsec":
		subset = snvf.join(vcf_df, left_on=["Chr", "Pos", "Ref", "Alt"], right_on=["chrom", "pos", "ref", "alt"], how="semi")
	else:
		subset = snvf.join(vcf_df, on=["chrom", "pos", "ref", "alt"], how="semi")

	
	# Write output
	output_dir = f"{filter_type}/{model}/{sample_name}"
	os.makedirs(output_dir, exist_ok=True)
	if model in ["microsec", "ideafix", "gatk-obmm"]:
		output_path = f"{output_dir}/{sample_name}.{model}.tsv"
	else:
		output_path = f"{output_dir}/{sample_name}.{model}.snv"
	subset.write_csv(output_path, separator="\t")
  
	print(f"\tWritten subsetted data to {output_path}\n")

print("All files processed.")




