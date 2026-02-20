#!/usr/bin/env bash

set -euo pipefail

Rscript eval-ffpolish_precrec.R
Rscript eval-gatk-obmm_precrec.R
Rscript eval-ideafix_precrec.R
Rscript eval-mobsnvf_precrec.R
Rscript eval-sobdetector_precrec.R
Rscript eval-vafsnvf_precrec.R
Rscript combine_results.R
Rscript make-plots.R
