#!/bin/bash

# Extract the first entry from the fastq files

mkdir -p fq_head

for x in ../fq/*/*.fastq.gz; do
	dirname=$(dirname $x)
	sample_name=$(basename $dirname)
	echo $sample_name
	zcat $x | head -4 > fq_head/${sample_name}.fq
done
