import polars as pl
import glob
import os
import logging
import sys

repo_root = ".."

## Local Dependencies
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from common import read_variants, snv_filter


## Set up logging configuration
logging.basicConfig(
	stream=sys.stdout, 
	level=logging.INFO, 
	format='%(message)s', # %(asctime)s %(levelname)s:
	force=True
)

## Maps dataset accession to dataset names
dataset_dict = {
    "PRJEB8754" : "ENA_Betge15",
    "SRP044740" : "ENA_SRP044740",
    "PRJEB44073" : "ENA_Chong21",
    "SRP044740" : "ENA_Oh15",
	"EGAD00001004066" : "EGA_Bonnet18"
}

def process_dataset(dataset: str, variant_set: str, annot_path: str, vcf_suffix: str = ".vcf.gz") -> None:
	"""
	Process a dataset to generate ground-truth variant tables for FFPE samples.

	For each FFPE sample, compares variants to matched frozen samples from the same case,
	marks variants present in any frozen sample as 'truth', and writes the result to disk.

	Args:
		dataset: Dataset identifier.
		variant_set: Variant set identifier.
		annot_path: Path to the sample annotation table for the specified dataset (TSV).
		vcf_suffix: Suffix for VCF files (default: '.vcf').
	"""
	logging.info(f"Processing dataset: {dataset}, variant set: {variant_set}")
	
	# Read the sample annotation table
	annot = pl.read_csv(annot_path, separator="\t")
	# Remove normal samples from table
	if "sample_type" in annot.columns:
		annot = annot.filter(pl.col("sample_type").str.to_lowercase() != "normal")
	# Get a list of only the FFPE Tumor Samples
	ffpe_samples = annot.filter(pl.col("preservation") == "FFPE").get_column("sample_name").to_list()

	# Find all VCF file paths
	vcf_paths = sorted(glob.glob(f"{repo_root}/vcf/{dataset}/{variant_set}/*/*{vcf_suffix}"))
	# Filter the paths only keeping the FFPE samples
	ffpe_paths = [
		path for path in vcf_paths
		if any(sample in path for sample in ffpe_samples)
	]
	
	logging.info(f"Found {len(ffpe_paths)} FFPE VCF files")


	for i, path in enumerate(ffpe_paths, start=1):
		# Extract FFPE sample name and case ID
		ffpe_sample_name = path.split("/")[-2]
		sample_annot = annot.filter(pl.col("sample_name") == ffpe_sample_name)
		case_id = sample_annot[0, "case_id"]
		preservation = sample_annot[0, "preservation"]

		logging.info(f"{i}. Processing FFPE sample: {ffpe_sample_name}")

		# Read and filter FFPE variants
		ffpe = (
			read_variants(path)
			.pipe(snv_filter)
			.with_columns(
				pl.lit(dataset_dict.get(dataset, dataset)).alias("dataset"),
				pl.lit(ffpe_sample_name).alias("sample_name")
			)
			.select(["dataset", "sample_name", "chrom", "pos", "ref", "alt"])
		)

		# Get matching frozen sample names for the same case
		ff_sample_names = annot.filter(pl.col("preservation") == "Frozen", pl.col("case_id") == case_id).get_column("sample_name")
		logging.info(f"\t{len(ff_sample_names)} matched fresh frozen samples found")

		for j, ff_sample in enumerate(ff_sample_names, start=1):
			# Read frozen sample, filter variants
			logging.info(f"\t\tComparing with frozen sample {j}: {ff_sample}")
			ff_sample_path = f"{repo_root}/vcf/{dataset}/{variant_set}/{ff_sample}/{ff_sample}{vcf_suffix}"
			ff_col_name = f"in_ff_{j}"

			ff = (
				read_variants(ff_sample_path)
				.pipe(snv_filter)
				.with_columns(pl.lit(True).alias(ff_col_name))
			)

			# Join FFPE and frozen sample on variant coordinates
			ffpe = (
				ffpe
				.join(ff, how="left", on=["chrom", "pos", "ref", "alt"])
				.with_columns(pl.col(ff_col_name).fill_null(False))
			)

		# Mark variants present in any frozen sample as "truth"
		ff_cols = [col for col in ffpe.columns if "in_ff_" in col]
		ffpe = ffpe.with_columns(pl.any_horizontal(ff_cols).alias("truth"))

		# Write ground-truth table to output directory
		outdir = f"{repo_root}/ground-truth/{dataset}/{variant_set}/{ffpe_sample_name}"
		os.makedirs(outdir, exist_ok=True)
		outpath = f"{outdir}/{ffpe_sample_name}.ground-truth.tsv"
		ffpe.write_csv(outpath, separator="\t")
		logging.info(f"\tGround-truth written to: {outpath}")
		


process_dataset(
    dataset = "EGAD00001004066",
    variant_set = "somatic_filtered",
    annot_path = f"{repo_root}/annot/EGAD00001004066/sample_annotations.tsv"
)


process_dataset(
    dataset = "EGAD00001004066",
    variant_set = "somatic_filtered-dp20",
    annot_path = f"{repo_root}/annot/EGAD00001004066/sample_annotations.tsv"
)

process_dataset(
    dataset = "EGAD00001004066",
    variant_set = "somatic_filtered-dp20-blacklist",
    annot_path = f"{repo_root}/annot/EGAD00001004066/sample_annotations.tsv"
)
