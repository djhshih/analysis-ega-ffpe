import os
import glob
import pandas as pd

target_dir_root = "fq"

fastq_annotation = pd.read_csv("../annot/sample_file.csv")
fastq_annotation

sample_annotation = pd.read_csv("../annot/samples.csv")
sample_annotation

phenotypes = sample_annotation[["alias", "phenotype"]]

fastq_annotation_new = pd.merge(fastq_annotation, phenotypes, left_on="sample_alias", right_on="alias")
fastq_annotation_new.to_csv("../annot/sample_file_phenotype.tsv", sep="\t", index=False)

fastq_paths = glob.glob("fq_old/*/*fastq.gz")

for fastq_path in fastq_paths:
    file_accession = fastq_path.split("/")[-2]
    
    annot = fastq_annotation_new[fastq_annotation_new.file_accession_id == file_accession]
    
    file_name = annot["file_name"].iloc[0]
    phenotype = annot["phenotype"].iloc[0].replace(" ", "-")
    alias = annot["alias"].iloc[0]
    
    target_dir = f"{target_dir_root}/{phenotype}_{alias}"
    os.makedirs(target_dir, exist_ok=True)
    
    target_path = f"{target_dir}/{phenotype}_{alias}_{file_name}"
    
    os.link(fastq_path, target_path)
    print(f"Linked: \n\t{fastq_path}\nTo:\n\t{target_path}\n")
    
print("Done.")


