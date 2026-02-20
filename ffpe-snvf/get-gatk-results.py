import polars as pl
from pathlib import Path
import sys
from tqdm import tqdm

repo_root = Path("..")

## Local Dependencies
sys.path.append(str(repo_root / "common-ffpe-snvf/python"))
from common import read_variants

## Functions
def process_dataset(dataset: str = "EGAD00001004066", variant_set: str = "somatic_filtered") -> None:
	"""Read the variant and filter columns from VCF and write 
	a table for evaluating the gatk orientation bias mixture model.

	Parameters:
		dataset: Dataset identifier (e.g. 'EGAD00001004066').
		variant_set: Variant set folder name (e.g. 'somatic_filtered').
	"""
	outdir_root = (repo_root / "ffpe-snvf" / dataset / variant_set / "gatk-obmm")
	vcf_paths = list((repo_root / "vcf" / dataset / variant_set).glob("*/*FFPE*.vcf.gz"))

	for path in tqdm(vcf_paths):
		sample_name = path.parent.name
		outdir = (outdir_root / sample_name)
		outdir.mkdir(exist_ok=True, parents=True)

		gatk_res = read_variants(path, columns=["#CHROM", "POS", "REF", "ALT", "FILTER"])
		gatk_res.write_csv((outdir / sample_name).with_suffix(".gatk-obmm.tsv"), separator="\t")

## Get GATK OBMM results
process_dataset(dataset = "EGAD00001004066", variant_set = "somatic_filtered")

