[![CI (smoke)](https://github.com/USER/REPO/actions/workflows/ci.yml/badge.svg)](https://github.com/USER/REPO/actions/workflows/ci.yml)

# FUT — Paper I (Negentropic Birth–Death on Causal Sets)

This repo contains:
- `paper/` — LaTeX + BibTeX (appendix with BDG closed-form)
- `bench/` — Minkowski calibration to build r_d lookup for Myrheim–Meyer
- `src/` — operators & dimension estimators; phase-scan harness
- `figs/` — scripts to generate figures (saved in `figs/out/`)
- `config/` — parameter grids and seeds
- `Makefile` — `make all` builds figures and compiles the paper

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
