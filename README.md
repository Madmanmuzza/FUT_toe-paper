[![CI (smoke)](https://github.com/Madmanmuzza/FUT_toe-paper/actions/workflows/ci.yml/badge.svg)](https://github.com/USER/REPO/actions/workflows/ci.yml)
**Latest release:** see the Releases page (v0.1.0).  
**Paper PDF:** see paper/Paper_I_full_with_figures_v3.pdf.

# FUT — Paper I (Negentropic Birth–Death on Causal Sets)

This repo contains:
- `paper/` — LaTeX + BibTeX (appendix with BDG closed-form)
- `bench/` — Minkowski calibration to build r_d lookup for Myrheim–Meyer
- `src/` — operators & dimension estimators; phase-scan harness
- `figs/` — scripts to generate figures (saved in `figs/out/`)
- `config/` — parameter grids and seeds
- `Makefile` — `make all` builds figures and compiles the paper
**Paper II (camera-ready PDF):** [paper2/Paper_II_camera_ready.pdf](paper2/Paper_II_camera_ready.pdf)
**Paper III (WIP, TeX):** [paper3/paper_III.tex](paper3/paper_III.tex)  
**Paper III (WIP, PDF):** [paper3/Paper_III_WIP.pdf](paper3/Paper_III_WIP.pdf)  
**Latest release:** https://github.com/Madmanmuzza/FUT_toe-paper/releases/latest

## Quick start
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
make bench     # build r_d calibration tables
make phase     # run 5x5 grid phase scan (set N and seeds in config)
make figs      # generate all figures into figs/out/
make paper     # compile LaTeX (needs system LaTeX)
make all       # bench + phase + figs + paper
```

## Notes
- BDG layers are classified by interval size |I(y,x)|, not graph distance.
- Spectral dimension uses **lazy** random walks; window chosen by AIC.
- Myrheim–Meyer uses a **monotone** interpolant from `bench/` calibration.


## Continuous Integration
This repository includes a lightweight GitHub Actions workflow that runs a **smoke test** on each push/PR:
- builds the **Myrheim–Meyer calibration** (CI-scale),
- runs a **tiny phase scan** (CI-scale),
- and generates the **figures** into `figs/out/`.

## Paper II — Continuum Emergence & Invariance Tests
- TeX: `paper2/paper_II.tex`
- Figures: `paper2/figs/`
- PDF: `paper2/Paper_II_full_draft.pdf` (if present)
- Notes: BDG error curve is a placeholder until the curved BDG script is run.
