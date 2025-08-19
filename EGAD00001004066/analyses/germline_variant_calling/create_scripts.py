#!/usr/bin/env python
import os
import glob

outdir = "per_sample_scripts"
os.makedirs(outdir, exist_ok=True)

normal_bams = [path for path in glob.glob("../../data/bam/*/*.bam") if "Normal" in path]

template_path = "call_germline_variants.sh"

for path in normal_bams:
    
    sample_name = path.split("/")[-2]
    
    with open(template_path, "r") as f:
        template = f.read()
    
    script = template.replace("$1", path)
    
    outpath = f"{outdir}/{sample_name}_call-germline-variants.sh"
    
    with open(outpath, "w") as s:
        s.write(script)


