#!/bin/bash

set -euo pipefail

root_outdir="somatic_vcf-dp10"
mkdir -p $root_outdir

filter_expression='MIN(FMT/DP)>=10'
echo -e "Filtering Expression: $filter_expression"

for vcf in somatic_vcf/*/*.vcf.gz; do
    
    filename=$(basename $vcf)
    sample_name=${filename%%.*}

    echo $sample_name
    
    outdir="${root_outdir}/${sample_name}"
    mkdir -p $outdir

    echo -e "\nFiltering $filename"
    bcftools view -i "$filter_expression" $vcf -o "${outdir}/${sample_name}.vcf.gz"
    echo "Filtered VCF saved to ${outdir}/${sample_name}.vcf.gz"

    echo "Indexing ${outdir}/${sample_name}.vcf.gz"
    bcftools index -t "${outdir}/${sample_name}.vcf.gz"

done

echo -e "\nFinished."