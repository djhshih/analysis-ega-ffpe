#!/usr/bin/env bash

set -euo pipefail

# --- 1. Configuration ---

# Input/Output Directories
outdir_root="../somatic_vcf"

# Input BAM Files
tumor_bam="<<tumor-bam-path>>"

# An array of all normal BAM files to be used for joint calling.
normal_bams=("<<normal-bam-paths>>"
)

# Reference Genome and Resource Files
reference_genome="../../../data/ref/Homo_sapiens_assembly38.fasta"
wgs_calling_regions="../../../data/ref/wgs_calling_regions.hg38.interval_list"

# GATK Resource Bundle Files
panel_of_normals="../../../data/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
germline_resource="../../../data/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
# This VCF is used for GetPileupSummaries to check for contamination at common variant sites
contamination_vcf="../../../data/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"

# --- 2. Setup - AUTOMATICALLY CONFIGURED ---

# Derive a unique sample name from the tumor BAM path
sample_name=$(basename "$(dirname "$tumor_bam")")
outdir="${outdir_root}/${sample_name}"
intermediate_dir="${outdir}/intermediates"

# Create output directories
mkdir -p "$outdir"
mkdir -p "$intermediate_dir"

echo "--- Starting GATK Somatic Workflow ---"
echo "Tumor Sample: ${sample_name}"
echo "Output Directory: ${outdir}"
echo "--------------------------------------"


# --- 3. Build Mutect2 Arguments for Joint Calling ---

# Prepare arguments for all normal samples
normal_args=()
for normal_bam_path in "${normal_bams[@]}"; do
    # Get the sample name from the BAM header to use with the -normal flag
    normal_sample_name=$(gatk GetSampleName -I "$normal_bam_path" -O /dev/stdout)
    normal_args+=(-I "$normal_bam_path")
    normal_args+=(-normal "$normal_sample_name")
done

# --- 4. Run Mutect2 for Joint Variant Calling ---

echo -e "\n1: Running Mutect2 with 1 tumor and ${#normal_bams[@]} normals...\n"
gatk --java-options "-Xmx8g" Mutect2 \
    -R "$reference_genome" \
    -L "$wgs_calling_regions" \
    -I "$tumor_bam" \
    "${normal_args[@]}" \
    --germline-resource "$germline_resource" \
    --panel-of-normals "$panel_of_normals" \
    --f1r2-tar-gz "${intermediate_dir}/${sample_name}.f1r2.tar.gz" \
    -O "${intermediate_dir}/${sample_name}.raw.vcf.gz"

# --- 5. Learn Read Orientation Model for Filtering ---

echo -e "\n2: Learning read orientation model...\n"
gatk LearnReadOrientationModel \
    -I "${intermediate_dir}/${sample_name}.f1r2.tar.gz" \
    -O "${intermediate_dir}/${sample_name}.read-orientation-model.tar.gz"

# --- 6. Estimate Contamination ---
# This requires pileup summaries from the tumor and ONE matched normal.
# We will use the first normal in the array as the matched normal for this estimation.

echo "\n3: Estimating contamination...\n"
matched_normal_bam="${normal_bams[0]}"
matched_normal_name=$(basename "$(dirname "$matched_normal_bam")")

echo -e "\n\\t- Generating pileup summary for tumor: ${sample_name}\n"
gatk GetPileupSummaries \
    -I "$tumor_bam" \
    -V "$contamination_vcf" \
    -L "$contamination_vcf" \
    -O "${intermediate_dir}/${sample_name}.pileups.table"

echo "\n\t3.1 Generating pileup summary for ONE of the matched normal: ${matched_normal_name}\n"
gatk GetPileupSummaries \
    -I "$matched_normal_bam" \
    -V "$contamination_vcf" \
    -L "$contamination_vcf" \
    -O "${intermediate_dir}/${matched_normal_name}.pileups.table"

echo -e "\n\t3.2Calculating contamination for the tumor sample...\n"
gatk CalculateContamination \
    -I "${intermediate_dir}/${sample_name}.pileups.table" \
    -matched "${intermediate_dir}/${matched_normal_name}.pileups.table" \
    -O "${intermediate_dir}/${sample_name}.contamination.table" \
    --tumor-segmentation "${intermediate_dir}/${sample_name}.segments.table"

# --- 7. Filter Mutect2 Calls ---

echo -e "\nStep 4: Filtering Mutect2 calls...\n"
gatk FilterMutectCalls \
    -V "${intermediate_dir}/${sample_name}.raw.vcf.gz" \
    -R "$reference_genome" \
    --contamination-table "${intermediate_dir}/${sample_name}.contamination.table" \
    --tumor-segmentation "${intermediate_dir}/${sample_name}.segments.table" \
    --ob-priors "${intermediate_dir}/${sample_name}.read-orientation-model.tar.gz" \
    -O "${intermediate_dir}/${sample_name}.filtered.vcf.gz"


echo -e "\nStep 5: Selecting filtered variants with BCFtools...\n"
filter_expression='(FILTER="PASS" | FILTER="orientation")'

bcftools view \
    -i "$filter_expression" \
    -o "${outdir}/${sample_name}.vcf.gz" \
    $vcf

bcftools index -t "${outdir}/${sample_name}.vcf.gz"


echo "--- Workflow Complete ---"
echo "Final filtered VCF is located at: ${outdir}/${sample_name}.filtered.vcf.gz"
