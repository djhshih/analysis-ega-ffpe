#!/bin/bash

set -euo pipefail

normal_bam=$1
outdir_root="../../data/germline_vcf"
bam_out_root="../../data/bqsr_bam"

reference_genome="../../../common_data/ref/Homo_sapiens_assembly38.fasta"
dbsnp="../../../common_data/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
mills_indels="../../../common_data/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"


sample_name=$(basename "$(dirname "$normal_bam")")
outdir="${outdir_root}/${sample_name}"
intermediate_dir="${outdir}/intermediates"

# Create output directories
mkdir -p "$outdir"
mkdir -p "$intermediate_dir"

echo "--- Starting GATK Germline Short Variant Discovery Workflow ---"
echo "Sample: ${sample_name}"
echo "Output Directory: ${outdir}"
echo "--------------------------------------"


# # 1a. Generate the recalibration table
# echo "Running BaseRecalibrator for ${sample_name}..."
# gatk BaseRecalibrator \
#     -R "$reference_genome" \
#     -I "$normal_bam" \
#     --known-sites "$dbsnp" \
#     --known-sites "$mills_indels" \
#     -O "${output_dir}/${sample_name}.recal_data.table"

# # 1b. Apply the recalibration
# echo "Running ApplyBQSR for ${sample_name}..."
# gatk ApplyBQSR \
#     -R "$reference_genome" \
#     -I "$input_bam" \
#     -bqsr "${output_dir}/${sample_name}.recal_data.table" \
#     -O "$analysis_ready_bam"


# --- Run HaplotypeCaller ---
echo "Running HaplotypeCaller for ${sample_name}..."
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R "$reference_genome" \
    -I "$normal_bam" \
    -O "${intermediate_dir}/${sample_name}.raw.vcf"


# echo "Applying hard filters for ${sample_name}..."
# gatk VariantFiltration \
#     -R "$reference_genome" \
#     -V "${outdir}/${sample_name}.vcf" \
#     -O "${outdir}/${sample_name}.hard-filtered.vcf" \
#     --filter-name "QD_filter" -filter "QD < 2.0" \
#     --filter-name "FS_filter" -filter "FS > 60.0" \
#     --filter-name "SOR_filter" -filter "SOR > 3.0" \
#     --filter-name "MQ_filter" -filter "MQ < 40.0" \
#     --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
#     --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"