CONFIG?=config/params.yaml
PY=python -u

all: bench phase figs paper

bench:
	$(PY) bench/calibration.py --config $(CONFIG)

phase:
	$(PY) src/experiments/phase_scan.py --config $(CONFIG)

figs:
	$(PY) figs/make_phase_diagram.py
	$(PY) figs/make_dimension_plateaus.py
	$(PY) figs/make_operator_accuracy.py
	$(PY) figs/make_interval_distance.py

paper:
	cd paper && pdflatex -interaction=nonstopmode paper_I.tex && bibtex paper_I && pdflatex -interaction=nonstopmode paper_I.tex && pdflatex -interaction=nonstopmode paper_I.tex

clean:
	rm -rf figs/out/* bench/out/* paper/*.aux paper/*.bbl paper/*.blg paper/*.log paper/*.out paper/*.pdf
