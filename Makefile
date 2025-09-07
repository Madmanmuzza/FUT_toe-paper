.PHONY: bench phase figs all

bench:
\t@echo "Stub bench: generate calibration placeholders"
\tpython scripts/stub_pipeline.py bench

phase:
\t@echo "Stub phase: generate results placeholders"
\tpython scripts/stub_pipeline.py phase

figs:
\t@echo "Stub figs: write run_report.json for I–III demo"
\tpython scripts/stub_pipeline.py figs
