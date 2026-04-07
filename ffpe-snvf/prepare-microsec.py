#!/usr/bin/env python
import polars as pl
import os
import glob
import pysam
import sys

repo_root = ".."
# Local Dependencies
sys.path.append(f"{repo_root}/common-ffpe-snvf/python")
from microsec_preparation import *
from common import return_path_if_exists

### Main Execution
#### Resource path and Directory Setup
outdir_root = f"{repo_root}/ffpe-snvf"

# MicroSEC batch script dir
batch_script_dir = "script_microsec"
os.makedirs(batch_script_dir, exist_ok=True)

simple_repeat_path = return_path_if_exists(f"{repo_root}/data/ref/ucsc_simple-repeat_hg38.bed", abs=True)
ref_path = return_path_if_exists(f"{repo_root}/data/ref/Homo_sapiens_assembly38.fasta", abs=True)
microsec_script_path = return_path_if_exists(f"{repo_root}/common-ffpe-snvf/R/microsec.R", abs=True)

## Wrapper function to prepare inputs and execution script for a dataset
def prepare_dataset_inputs(
		dataset,
		variant_set,
		sample_paths: pl.DataFrame, 
		vcf_path_col: str = "vcf_path", 
		bam_path_col: str = "bam_path",
		sample_name_col: str = "sample_name",
		simple_repeat_path: str = simple_repeat_path,
		microsec_script_path: str = microsec_script_path,
		ref_path: str = ref_path,
		batch_script_dir: str = batch_script_dir,
		ct_only: bool = False,
		snv_only: bool = False,
		lazy: bool = False
	):

	microsec_root = f"{outdir_root}/{dataset}/{variant_set}/microsec"

	for i in range(sample_paths.height):

		## sample specific variable setup
		sample_name = sample_paths[i, sample_name_col]
		vcf_path = sample_paths[i, vcf_path_col]
		bam_path = sample_paths[i, bam_path_col]

		print(f"{i+1}. Creating MicroSEC input for {sample_name} from VCF: {vcf_path}")

		###--------------------------------------------


		# Read in the simple repeats bed file to annotate variants later
		simple_repeats_bed = pl.read_csv(
			simple_repeat_path, 
			separator="\t", 
			has_header=False, 
			new_columns=["chrom", "start", "end", "type", "score"],
		).select(["chrom", "start", "end"])
		
		###--------------------------------------------
		
		## Prepare mut_info for the sample
		output_dir = f"{microsec_root}/{sample_name}/inputs"
		os.makedirs(output_dir, exist_ok=True)
		mut_info_path = os.path.abspath(f"{output_dir}/{sample_name}.microsec.mut-info.tsv")

		if not (lazy and os.path.exists(mut_info_path)):
			print("\tPreparing mutation info")
			mut_info = make_mut_info(
				vcf_path=vcf_path, 
				ref_fasta_path=ref_path, 
				sample_name=sample_name,
				simple_repeats_bed=simple_repeats_bed,
				ct_only=ct_only, 
				snv_only=snv_only
			)
		
			## Write mut_info
			mut_info.write_csv(mut_info_path, separator="\t")

		###----------------------------------------------

		sample_info_path = os.path.abspath(f"{microsec_root}/{sample_name}/inputs/{sample_name}.microsec.sample_info.tsv")
		
		## Prepare sample_info for the sample
		if not (lazy and os.path.exists(sample_info_path)):
			print("\t - Preparing sample info")
			sample_info = make_sample_info(sample_name, mut_info_path, bam_path, ref_path, simple_repeat_path)

			## Write sample_info
			pl.DataFrame(sample_info).write_csv(sample_info_path, separator="\t")

		###-----------------------------------------------

		# Write batch execution scripts
		print("\tPreparing execution script")
		script_path = f"{batch_script_dir}/{sample_name}.microsec.sh"
		microsec_outdir = os.path.abspath(f"{microsec_root}/{sample_name}")
		
		with open(script_path, "w") as file:
			file.writelines("#!/usr/bin/env bash\n")
			file.writelines(f"Rscript {microsec_script_path} --sample_info '{sample_info_path}' --output_dir '{microsec_outdir}'\n")



#### Prepare EGAD00001004066
## Create annotation table linking ffpe samples' vcf path and bam path 
dataset = "EGAD00001004066"
variant_set = "somatic_filtered"
ffpe_vcf_paths = [os.path.abspath(path) for path in sorted(glob.glob(f"{repo_root}/vcf/{dataset}/{variant_set}/*/*.vcf.gz")) if "FFPE" in path]

sample_paths = (
	pl.DataFrame({
		"sample_name" : [path.split("/")[-2] for path in ffpe_vcf_paths],
		"vcf_path" : ffpe_vcf_paths,
		"bam_path" : [return_path_if_exists(f"{repo_root}/data/{dataset}/bam/{path.split("/")[-2]}/{path.split("/")[-2]}.bam", abs=True) for path in ffpe_vcf_paths],
	})
)

## Prepare microsec inputs and batch execution scripts
prepare_dataset_inputs(
	dataset = dataset, 
	variant_set = variant_set, 
	sample_paths = sample_paths
)

