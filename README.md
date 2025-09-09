# Analysis EGA FFPE

This repository aims to analyze several datasets from EGA with the goal to evaluate the several FFPE artifact filtering tools. Here, FFPE artifact filtering is tested in the context of somatic and germline variant calling. Somatic variants were called using Mutect2 with and without matched normal.

## Studies

- **EGAS00001002631:**
    - **Dataset:** EGAD00001004066
    - **Status:** sent request on EGA on 25/07/07; Material Transfer Agreement (MTA) on 2025-07-22; Access granted on 2025-07-22 
    - **Time required for access:** 15 days
    - **Sequencing Method:** WXS
    - **Publication:** https://doi.org/10.1371/journal.pone.0195471
    - **Description:** This contains data from 1 patients with match FFPE and Frozen samples from the same tumor in two different tissue type (Colon and Liver). The difference between these sample being the DNA extraction kit used. The dataset provides both tumoral and normal samples from each tissue type.


- **EGAS00001007521:**
    - **Dataset:** EGAD00001011331
    - **Status:** sent email to TDOadmin@phsa.ca on 25/07/07; MTA signed on 2025-07-11; waiting
    - **Time required for access:** Access Pending
    - **Sequencing Method:** WGS
    - **Publication:** None
    - **Description:** HGSC cases were selected for which matched fresh frozen and FFPE samples were available from the same tissue specimens.  Fresh frozen, FFPE and normal blood were subject to WGS for the purpose of assessing the possibility of using FFPE WGS in place of fresh frozen for somatic mutation calling. 6 BAM available in the repository. No other information is provided publicly.


- **EGAS00001000737:**
    - **Datasets:** EGAD00001000830, EGAD00001000831, EGAD00001000832, EGAD00001000834 
    - **Status:** sent request on EGA; sent email to Paul Flicek (flicek@ebi.ac.uk), Henrik Hornshøj (hhj@ki.au.dk) and Søren Vang (vang@ki.au.dk) on 25/07/07; waiting
    - **Time required for access:** Access Pending
    - **Sequencing Method:** WXS, RNA-Seq
    - **Publication:** https://doi.org/10.1371/journal.pone.0098187
    - **Description:** Next-Generation Sequencing of RNA and DNA Isolated from Paired Fresh-Frozen and Formalin-Fixed Paraffin-Embedded Samples of Human Cancer and Normal Tissue.

        - EGAD00001000830: The paired FF/FFPE prostate set. Consists of 7 paired FF/FFPE prostate carcinoma samples, stored for 2 – 11 years.

        - EGAD00001000831: The paired FF/FFPE colon set, RNA-Seq. Consists of 19 paired FF/FFPE colorectal carcinoma samples; in 13 of these, matching normal FF colon samples. This set had been stored for 2 – 13 years.
        
        - EGAD00001000832: The paired FF/FFPE bladder set. Consists of 8 paired FF/FFPE bladder carcinoma samples, stored for 5 – 9 years.

        - EGAD00001000834: The paired FF/FFPE colon set, Exome-Seq. Consists of 19 paired FF/FFPE colorectal carcinoma samples; in 13 of these, matching normal FF colon samples. This set had been stored for 2 – 13 years.

        The publication also mentions a test set with a total of 16 samples from normal liver, normal (reactive) tonsil, normal skin (from mastectomy specimen), liver, lung, breast, and bladder carcinomas. Matching FF samples available from the liver and tonsil samples. However, this is not seen in the repository.


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

**R libraries:**

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

**Python Libraries:**
- polars
- pandas
- numpy
- matplotlib
- seaborn
- pysam


## FFPE filtering tools status

- **mobsnvf:** Executed
- **vafsnvf:** Executed
- **SOBdetector:** Executed
- **GATK Orientation Bias Mixture Model**: Executed
- **MicroSEC:** Execution Pending
- **IdeaFix:** Installation Pending
- **FFPolish:** Installation Failed


## To do
- Complete writing Plotting Script
- Run MicroSEC
- Try installing IdeaFix
- Make Mutation Signature plots
- Make summary statistics tables
- Try out Mutect2 tumor only mode on the EGAD00001004066 dataset
- Create correlation table between SNV summary statistics and AUROC & AUPRC for each tool.


## Issues

Precrec does not handle scores with NA values correctly i.e it does not ignore them during evaluation as it should. Due to this NA values should be removed from the score set before being passed to precrec.


