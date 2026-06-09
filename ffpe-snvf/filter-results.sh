#!/usr/bin/env bash

set -euox pipefail

python non-exome-exclusion.py
python blacklist-exclusion.py
python micr-exclusion.py
