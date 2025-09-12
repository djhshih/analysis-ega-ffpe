import os
import glob
import pandas as pd

target_dir_root = "fq"

fastq_annotation = pd.read_csv("../../../annot/EGAD00001004066/sample_annotations.tsv", sep="\t")
# fastq_annotation

fastq_paths = glob.glob("fq_old/*/*fastq.gz")
# fastq_paths

for fastq_path in fastq_paths:
    file_accession = fastq_path.split("/")[-2]
    
    annot = fastq_annotation[fastq_annotation.fastq_accession_id == file_accession]
    
    fq_fname = annot["fastq_file_name"].iloc[0]
    description = annot["title"].iloc[0].replace(" ", "-")
    alias = annot["sample_alias"].iloc[0]
    
    target_dir = f"{target_dir_root}/{description}_{alias}"
    os.makedirs(target_dir, exist_ok=True)
    
    target_path = f"{target_dir}/{description}_{alias}_{fq_fname}"
    
    os.link(fastq_path, target_path)
    print(f"Linked: \n\t{fastq_path}\nTo:\n\t{target_path}\n------")
    
print("Done.")


