#!/usr/bin/env python
import os
import glob
import polars as pl

# %%
def return_path_if_exists(path: str, abs: bool =False) -> str:
	if abs:
		path = os.path.abspath(path)
	if os.path.exists(path):
		return path
	else:
		raise FileNotFoundError(f"Path {path} does not exist.")

# %%
ref_path = return_path_if_exists("../data/ref/Homo_sapiens_assembly38.fasta", abs=True)
vcf_paths = [os.path.abspath(path) for path in sorted(glob.glob("../vcf/EGAD00001004066/somatic_vcf-dp10/*/*.vcf.gz"))]

filtered_outdir = return_path_if_exists("EGAD00001004066/somatic_vcf-dp10", abs=True)

model = "ffpolish"


# %%

for vcf_path in vcf_paths:

	# Skip Frozen samples
	if "Frozen" in vcf_path:
		continue

	vcf_filename = os.path.basename(vcf_path)
	sample_name = vcf_filename.split(".")[0]
	
	bam_path = return_path_if_exists(glob.glob(f"../data/EGAD00001004066/bam/{sample_name}/*.bam")[0], abs=True)
	
	outdir = "script_ffpolish"
	os.makedirs(outdir, exist_ok=True)

	content = [
		"#!/usr/bin/env bash",
		f"ffpolish filter -o {filtered_outdir}/ffpolish/{sample_name} -p {sample_name} {ref_path} {vcf_path} {bam_path}"
	]
	filename = f"{outdir}/{model}_{sample_name}.sh"
	
	with open(filename, "w") as f:
		f.writelines(content)

	print(f"Created {model} filtering script for {vcf_path}")

