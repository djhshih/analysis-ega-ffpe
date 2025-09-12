#!/bin/bash

mkdir -p rg

for f in fq_head/*.fq; do
        fname=${f##*/}
        rgname=${fname%.*}
        sample=${rgname}
        library=$(head -n 1 $f | sed 's/\:/\n/g' | tail -n 1)
        rgsam collect --input=$f --format fastq --qnformat illumina-1.8 --sample $sample --library $library | head -1 | sed 's/\t/\\t/g' > rg/$sample.rg.txt
done
