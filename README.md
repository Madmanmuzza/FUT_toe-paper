# FUT — Papers I–X (auto-figures bootstrap)

This is a clean scaffold of your repo with:
- Papers **I–X** as standalone LaTeX (each has a self-contained `.bib`, robust figure handling, references last).
- A workflow that **runs your pipeline** (`make bench`, `make phase`, `make figs`) *and* generates the exact **paper figures** the TeX expects, committing them back into `paperN/figs/`.
- A starter `scripts/generate_paper_figs.py` that reads `run_report.json` and writes the expected filenames for **Papers I–III** now (extend similarly for IV–X as their data exist).

## Quick start (local)
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Your pipeline (replace with your real steps)
make bench
make phase
make figs

# Generate paper figures from run_report.json
python scripts/generate_paper_figs.py

# Compile a paper (needs LaTeX on your system)
cd paper1 && pdflatex paper_I.tex && bibtex paper_I && pdflatex paper_I.tex && pdflatex paper_I.tex
```

## GitHub Actions (automatic figures)
On any push, the workflow `.github/workflows/fig-autobuild.yml` will:
1. Set up Python and install deps from `requirements.txt`.
2. Run `make bench && make phase && make figs` (replace with your real pipeline).
3. Run `python scripts/generate_paper_figs.py` to overwrite placeholder images with **real** ones.
4. Commit updated `paper1/figs/*.png`, `paper2/figs/*.png`, `paper3/figs/*.png` back to the repo.

> Extend `scripts/generate_paper_figs.py` with your outputs for Papers IV–X, then add those paths to the workflow's commit step.

## Where to drop your real code
- Put your actual simulation/analysis code under `bench/`, `src/`, `figs/`, and write numeric summaries to **`run_report.json`**.
- The figure generator expects keys as shown in the example `run_report.json` in repo root (edit to match your outputs).

## Important
- Figures currently include **placeholders generated from `run_report.json`**. Once your real pipeline writes that file (and/or additional CSVs), rerun the workflow and the figures will be replaced automatically.
