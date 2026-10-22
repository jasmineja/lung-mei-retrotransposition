#!/usr/bin/env bash
# Runs steps 2b-5 on the input tables in ./data (or $MEI_DATA_DIR). Results go to ./results.
set -euo pipefail
python 02_filtering/filter_calls.py
python 03_detection_power/detection_power.py
python 04_statistics/burden_statistics.py
python 05_figures/make_figures.py
