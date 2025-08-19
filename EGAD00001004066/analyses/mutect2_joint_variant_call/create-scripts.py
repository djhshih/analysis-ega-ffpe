#!/usr/bin/env python
import os
import glob
import polars as pl

outdir = "per_sample_scripts"
os.makedirs(outdir, exist_ok=True)

tumor_bam_paths = [path for path in glob.glob("../../data/bam/*/*.bam") if "Tumoral" in path]

template_script = "joint_call.template.sh"

with open(template_script, "r") as f:
    script = f.read()

for path in tumor_bam_paths:
    
    tissue = path.split("/")[-2].split("-")[1]
    
    normal_bam_paths = [path for path in glob.glob("../../data/bam/*/*.bam") if (("Normal" in path) & (f"{tissue}" in path) & ("Frozen" in path))]
    
    normals = "\n\t".join(normal_bam_paths)
    
    sample_script = script.replace('"<<tumor-bam-path>>"', path).replace('"<<normal-bam-paths>>"', normals)
    sample_name = os.path.basename(os.path.dirname(path))
    
    outpath = f"{outdir}/{sample_name}_multi-normal_joint-call.sh"
    
    with open(outpath, "w") as f:
        f.write(sample_script)
    


