#!/usr/bin/env bash
set -e
python -u bench/calibration.py
python -u src/experiments/phase_scan.py
python -u figs/make_phase_diagram.py
python -u figs/make_dimension_plateaus.py
python -u figs/make_operator_accuracy.py
python -u figs/make_interval_distance.py
