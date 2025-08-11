#!/bin/bash
# Download metadata

set -euo pipefail
IFS=$'\n\t'

# sign-in required

wget -O EGAD00001004066-metadata.zip https://profile.ega-archive.org/datasets/EGAD00001004066/metadata/download?format=csv
