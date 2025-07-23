#!/bin/bash
# Download data

set -euo pipefail
IFS=$'\n\t'

# NB  creditial file must be at ~/.config/ega/login
# NB  if an error regarding tqdm occurs, upgrade tqdm by
#     pip install tqdm --upgrade

pyega3 -d -cf ~/.config/ega/login fetch EGAD00001004066
