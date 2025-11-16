import os
import glob

def return_path_if_exists(path: str, abs: bool =False) -> str:
	if abs:
		path = os.path.abspath(path)
	if os.path.exists(path):
		return path
	else:
		raise FileNotFoundError(f"Path {path} does not exist.")


dataset = "EGAD00001004066"
ref_path = return_path_if_exists("../data/ref/Homo_sapiens_assembly38.fasta", abs=True)

vcf_paths = [return_path_if_exists(path, abs=True) for path in glob.glob(f"../vcf/{dataset}/somatic_vcf/*/*.vcf.gz")]

filtered_outdir = os.path.abspath(f"{dataset}/somatic_vcf")


models = ["mobsnvf", "vafsnvf", "sobdetector", "ideafix"]

for vcf_path in vcf_paths:
    
    vcf_filename = os.path.basename(vcf_path)
    sample_name = vcf_filename.removesuffix(".vcf.gz")
    bam_path = return_path_if_exists(f'../../bam/{dataset}/{sample_name}/{sample_name}.bam', abs=True)
    
    for model in models:
        
        if model == "ideafix":
            outdir = "script_ideafix"
            os.makedirs(outdir, exist_ok=True)
                    
            ## Template script for each model
            template_dir = return_path_if_exists(f"{model}.R", abs=True)
            content = f"#!/usr/bin/env Rscript\nRscript {template_dir} --vcf '{vcf_path}' --ref '{ref_path}' --outdir '{filtered_outdir}'\n"
            filename = f"{outdir}/{model}_{sample_name}.sh"
            
            with open(filename, "w") as f:
                f.write(content)
                
        else:
            
            outdir = "script"
            os.makedirs(outdir, exist_ok=True)
        
            template_dir = return_path_if_exists(f"{model}.sh.template", abs=True)
            content = f"#!/bin/bash\nbash {template_dir} {bam_path} '../{vcf_path}' '../{filtered_outdir}'  '{ref_path}'\n"
            filename = f"{outdir}/{model}_{sample_name}.sh"
            
            with open(filename, "w") as f:
                f.write(content)
            
        print(f"Created {model} filtering script for {vcf_path}")

