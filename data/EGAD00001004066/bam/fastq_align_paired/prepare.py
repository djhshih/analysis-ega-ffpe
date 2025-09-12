#!/usr/bin/env python3

# Prepare json input files for WDL fastq_aligned_paired workflow 

import os
import json
import pandas as pd

# ---

# absolute path is required in the json input files for WDL
ref_root = os.path.abspath('../../ref')
ref_fname = 'Homo_sapiens_assembly38.fasta'

infname = '../../../../annot/EGAD00001004066/fastq_annotations.tsv'
rg_root = os.path.abspath('../..rg/rg')
fq_root = os.path.abspath('../../fq')
out_dir = 'inputs'


if not os.path.exists(out_dir):
    os.makedirs(out_dir)

fastq_pd = pd.read_csv(infname, sep="\t")

# ----

samples = []
for _, x in fastq_pd.groupby('sample_alias'):
    alias = x["sample_alias"].iloc[0]
    title = x["title"].iloc[0].replace(" ", "-")
    
    fq1 = x["fastq_file_name"].iloc[0]
    fq2 = x["fastq_file_name"].iloc[1]
    
    rg = (f"{title}_{alias}")
    
    fastqs_r1 = os.path.join(fq_root, f"{title}_{alias}/{title}_{alias}_{fq1}")
    fastqs_r2 = os.path.join(fq_root, f"{title}_{alias}/{title}_{alias}_{fq2}")
    samples.append(
        {
            'sample_id': f"{title}_{alias}",
            "rg_name" : rg,
            'fastqs_r1': fastqs_r1,
            'fastqs_r2': fastqs_r2,
        }
    )

ref_fpath = os.path.join(ref_root, ref_fname)

# base wdl input
base = {
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta': ref_fpath,
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_amb': ref_fpath + '.64.amb',
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_ann': ref_fpath + '.64.ann',
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_alt': ref_fpath + '.64.alt',
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_pac': ref_fpath + '.64.pac',
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_bwt': ref_fpath + '.64.bwt',
    'fastq_align_paired.fastq_bwa_mem_paired.ref_fasta_sa': ref_fpath + '.64.sa',
    'fastq_align_paired.fastq_bwa_mem_paired.cpu': 4,
    'fastq_align_paired.fastq_bwa_mem_paired.memory_gb': 12,
    'fastq_align_paired.fastq_bwa_mem_paired.preemptible_tries': 1,
    'fastq_align_paired.bam_sort_coord.cpu': 4,
    'fastq_align_paired.bam_sort_coord.memory_gb': 4,
    'fastq_align_paired.bam_sort_coord.preemptible_tries': 1,
}

# write wdl input json file for each sample
for x in samples:
    out = base.copy()
    out['fastq_align_paired.sample_id'] = x['sample_id']
    out['fastq_align_paired.fastq_bwa_mem_paired.fastqs_r1'] = [x['fastqs_r1']]
    out['fastq_align_paired.fastq_bwa_mem_paired.fastqs_r2'] = [x['fastqs_r2']]
    out['fastq_align_paired.fastq_bwa_mem_paired.rg_header'] = os.path.join(rg_root, x["rg_name"] + '.rg.txt')
    out_path = os.path.join(out_dir, x['sample_id'] + '.inputs')
    with open(out_path, 'w') as outf:
        outf.write(json.dumps(out, indent=True, sort_keys=True))

