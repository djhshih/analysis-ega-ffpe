#!/bin/bash

set -euo pipefail

## This is as simplified Germline variant calling workflow based on GATK Haplotype Caller. 
## Due to it's simplified nature BQSR and VQSR is not used

normal_bam=$1
outdir_root="../germline_vcf"

reference_genome="../../../data/ref/Homo_sapiens_assembly38.fasta"

sample_name=$(basename "$(dirname "$normal_bam")")
outdir="${outdir_root}/${sample_name}"
intermediate_dir="${outdir}/intermediates"

# Create output directories
mkdir -p "$outdir"
mkdir -p "$intermediate_dir"

echo "--- Starting simplified germline variant calling workflow ---"
echo "Sample: ${sample_name}"
echo "Output Directory: ${outdir}"
echo "--------------------------------------"


# --- Run HaplotypeCaller ---
echo "Running HaplotypeCaller for ${sample_name}..."
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R "$reference_genome" \
    -I "$normal_bam" \
    -O "${intermediate_dir}/${sample_name}.raw.vcf.gz"

# --- Run VariantFiltration ---
echo "Applying filters to ${sample_name}..."
gatk VariantFiltration \
    -R "$reference_genome" \
    -V "${intermediate_dir}/${sample_name}.raw.vcf.gz" \
    -O "${intermediate_dir}/${sample_name}.filter.vcf.gz" \
    --filter-name "QD_filter" -filter "QD < 2.0" \
    --filter-name "FS_filter" -filter "FS > 60.0" \
    --filter-name "SOR_filter" -filter "SOR > 3.0" \
    --filter-name "MQ_filter" -filter "MQ < 40.0" \
    --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"


# --- Run BCFtools for filtering ---
echo "Selecting variants based on applied filters"
regions="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

bcftools view \
    -r "$regions" \
    -f "PASS" \
    -o "${outdir}/${sample_name}.vcf.gz" \
    "${intermediate_dir}/${sample_name}.filter.vcf.gz"

bcftools index -t "${outdir}/${sample_name}.vcf.gz"

echo "--- Workflow Complete ---"
echo "Final filtered VCF is located at: ${outdir}/${sample_name}.vcf.gz"

