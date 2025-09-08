#!/usr/bin/env python3
import os
import glob
# import shutil

bam_paths = glob.glob("cromwell-executions/fastq_align_paired/*/call-bam_sort_coord/execution/*.bam")
bai_paths = glob.glob("cromwell-executions/fastq_align_paired/*/call-bam_sort_coord/execution/*.bai")
bam_bai_paths = sorted(bam_paths + bai_paths)

target_root = "../../../bam/EGAD00001004066"



for path in bam_bai_paths:
    
    basename = os.path.basename(path)
    name = basename.split(".")[0]
    
    target_dir = f"{target_root}/{name}"
    os.makedirs(target_dir, exist_ok=True)
    
    target_path = f"{target_dir}/{basename}"
    
    if os.path.exists(target_path):
        print(f"'{basename}' already exists at '{target_dir}'\n\tSKIPPING...")
        continue
    
    os.link(path, target_path)
    print(f"Created hard link for: {basename}, \n\tat {target_dir}\n")

print("Done.")



