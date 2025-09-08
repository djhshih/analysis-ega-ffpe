import os
import glob

vcf_paths = glob.glob("../vcf/EGAD00001004066/somatic_vcf/*/*.vcf.gz")
filtered_outdir = "somatic_vcf"

models = ["mobsnvf", "vafsnvf", "sobdetector"]

outdir = "batch_scripts"
os.makedirs(outdir, exist_ok=True)

for vcf_path in vcf_paths:
    
    vcf_filename = os.path.basename(vcf_path)
    sample_name = vcf_filename.removesuffix(".vcf.gz")
    
    for model in models:
        content = f"#!/bin/bash\nbash ../{model}.sh.template '../../bam/EGAD00001004066/{sample_name}/{sample_name}.bam' '../{vcf_path}' '../{filtered_outdir}'\n"
        filename = f"{outdir}/{model}_{sample_name}.sh"
        
        with open(filename, "w") as f:
            f.write(content)
            
        print(f"Created {model} filtering script for {vcf_path}")

