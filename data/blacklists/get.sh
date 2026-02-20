#!/usr/bin/env bash

set -euo pipefail

wget -O encode_hg38-blacklist.v2.bed.gz -c "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz?raw=true"

wget -O ucsc_unusual_regions.bb -c https://hgdownload.soe.ucsc.edu/gbdb/hg38/problematic/comments.bb
wget -O grcExclusions.bb -c https://hgdownload.soe.ucsc.edu/gbdb/hg38/problematic/grcExclusions.bb
wget -O giab_alldifficultregions.bb -c https://hgdownload.soe.ucsc.edu/gbdb/hg38/problematic/GIAB/alldifficultregions.bb
# wget -O GCA_000001405.15_GRCh38_GRC_exclusions.bed -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_GRC_exclusions.bed

bigBedToBed ucsc_unusual_regions.bb ucsc_unusual_regions.bed
bigBedToBed grcExclusions.bb grcExclusions.bed
bigBedToBed giab_alldifficultregions.bb giab_alldifficultregions.bed

rm ucsc_unusual_regions.bb grcExclusions.bb giab_alldifficultregions.bb

zcat -f encode_hg38-blacklist.v2.bed.gz ucsc_unusual_regions.bed grcExclusions.bed giab_alldifficultregions.bed \
    | cut -f 1-3 \
    | sort -k1,1V -k2,2n \
    | bedtools merge -i - \
    | bgzip -c > master_blacklist.bed.gz

tabix -p bed master_blacklist.bed.gz

