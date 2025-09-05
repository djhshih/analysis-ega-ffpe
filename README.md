# Analysis EGA FFPE

This repository aims to analyze several datasets from EGA with the goal to evaluate the several FFPE artifact filtering methods in the context of somatic and germline variant calling.

## Datasets
- EGAD00001004066:

    This contains data from 1 patients with match FFPE and Frozen samples from the same tumor in two different tissue type (Colon and Liver). The difference between these sample being the DNA extraction kit used. The dataset provides both tumoral and normal samples from each tissue type.


## Dependencies

- R
- python
- cromwell
- gatk
- dlazy
- bcftools
- samtools
- htspan
- SOBDetector
- samblaster
- sambamba

R libraries:

- io
- precrec
- jsonlite
- argparser
- tidyverse
- glue
- patchwork
- grid
- hrbrthemes
- viridis

Python Libraries:
- polars
- pandas
- numpy
- matplotlib
- seaborn
- pysam



## Issues

Precrec does not handle scores with NA values correctly i.e it does not ignore them during evaluation as it should. Due to this NA values should be removed from the score set before being passed to precrec.





