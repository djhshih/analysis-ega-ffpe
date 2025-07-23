#!/bin/bash
# Download data

set -euo pipefail
IFS=$'\n\t'

pyega3 -d -cf ~/.config/ega/login fetch EGAD00001004066
