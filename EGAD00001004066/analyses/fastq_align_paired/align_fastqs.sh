#!/bin/bash
ls inputs/* > dlazy_samples.txt
djobs dlazy_samples.txt cromwell run ../../wdl/fastq_align_paired.wdl -i

pdlazy job -j 4