#!/bin/bash

set -euo pipefail

# Function to filter VCF based on exome target regions in provided bed
filter_vcf() {
    # Check if the required arguments are provided
    if [[ $# -lt 2 ]]; then
        echo "Usage: filter_vcf <indir_root> <outdir_root> [exome_bed_path]"
        return 1
    fi

    local indir_root="$1"
    local outdir_root="$2"
    # Use the 3rd argument if provided, otherwise default to the hardcoded path
    local exome_bed_path="${3:-../../data/regions/sureselect-all-exon-v5_hg38_regions_200bp-pad.bed}"

    echo -e "Filtering using exome target bed: $(realpath "$exome_bed_path")"
    echo -e "Input directory: $indir_root"
    echo -e "Output directory: $outdir_root\n"

    local i=1
    for vcf in "$indir_root"/*/*.vcf.gz; do
        
        # Skip if no files match the glob
        [[ -e "$vcf" ]] || continue
        
        local filename="${vcf##*/}"
        local sample_name="${filename%%.*}"

        local outdir="${outdir_root}/${sample_name}"
        mkdir -p "$outdir"

        echo -e "\t$i. Filtering: $filename"
        
        bcftools view -T "$exome_bed_path" -O z -o "${outdir}/${sample_name}.vcf.gz" "$vcf"
        bcftools index -t "${outdir}/${sample_name}.vcf.gz"

        ((i++))

    done

    echo -e "Finished processing: $indir_root.\n"
}

# Filter dataset
filter_vcf "somatic_filtered" "somatic_filtered-exome"
