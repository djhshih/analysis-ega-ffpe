#!/bin/bash

set -euo pipefail

root_outdir="somatic_filtered-dp10-blacklist"
mkdir -p $root_outdir

blacklist_path="../../data/blacklists/master_blacklist.bed.gz"
echo -e "Filtering using blacklist: $blacklist_path"

for vcf in somatic_vcf-dp10/*/*.vcf.gz; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $sample_name
    
    outdir="${root_outdir}/${sample_name}"
    mkdir -p $outdir

    echo -e "\nFiltering $filename"
    bcftools view -T ^"$blacklist_path" $vcf -Oz -o "${outdir}/${sample_name}.vcf.gz"
    echo "Filtered VCF saved to ${outdir}/${sample_name}.vcf.gz"

    echo "Indexing ${outdir}/${sample_name}.vcf.gz"
    bcftools index -t "${outdir}/${sample_name}.vcf.gz"

done

echo -e "\nFinished."