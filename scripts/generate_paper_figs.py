#!/usr/bin/env python3
import json, os, glob
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/.."
OUT = os.path.join(ROOT, "figs", "out")          # your existing pipeline output dir

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def save_fig(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()

def load_run_report():
    rr_path = os.path.join(ROOT, "run_report.json")
    if not os.path.exists(rr_path):
        return {}
    with open(rr_path, "r") as f:
        return json.load(f)

# ---------- Paper I ----------
def paper1(rr):
    tgt = os.path.join(ROOT, "paper1", "figs")
    ensure_dir(tgt)

    # 1) Interval abundance L_k and S^(4)
    L = rr.get("paper1", {}).get("interval_L", None)  # expected e.g. [L1,L2,L3,L4]
    if L and len(L) >= 4:
        k = np.arange(1,5)
        plt.figure(figsize=(6,3.6))
        plt.bar(k, L[:4])
        plt.xticks(k, [f"{ki}" for ki in k])
        plt.xlabel("Interval size k")
        plt.ylabel("Average L_k")
        S4 = 1 - L[0] + 9*L[1] - 16*L[2] + 8*L[3]
        plt.title("Interval Abundance (from run)")
        plt.text(0.55, max(L[:4])*1.02, f"S$^{{(4)}}$≈{S4:.3f}")
        save_fig(os.path.join(tgt, "interval_hist.png"))

    # 2) Myrheim–Meyer dimension vs N
    mm = rr.get("paper1", {}).get("mm_dimension", None)  # {"N":[...], "median":[...], "iqr":[...]}
    if mm and all(k in mm for k in ("N","median","iqr")):
        N = np.array(mm["N"], dtype=float)
        med = np.array(mm["median"], dtype=float)
        iqr = np.array(mm["iqr"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.errorbar(N, med, yerr=iqr/2, fmt='o-', capsize=4)
        plt.xscale("log", base=2)
        plt.xticks(N, [fr"$2^{{{int(np.log2(n))}}}$" for n in N])
        plt.axhline(4.0, linestyle='--')
        plt.xlabel("N (elements)")
        plt.ylabel("MM dimension (median ± IQR/2)")
        plt.title("Myrheim–Meyer vs N (from run)")
        save_fig(os.path.join(tgt, "mm_dimension.png"))

    # 3) Rate stability summaries
    rs = rr.get("paper1", {}).get("rate_stability", None)  # dict with N, link_med, ...
    if rs:
        N = np.array(rs["N"], dtype=float)
        def plot_band(x, y, iqr, fmt, label, alpha=0.2):
            plt.plot(x, y, fmt, label=label)
            y, iqr = np.array(y), np.array(iqr)
            plt.fill_between(x, y - iqr/2, y + iqr/2, alpha=alpha)

        plt.figure(figsize=(6,3.6))
        plot_band(N, rs["link_med"], rs["link_iqr"], 'o-', "Link density (median)")
        plot_band(N, rs["deg_in_med"], rs["deg_in_iqr"], 's-', "Mean in-degree (median)", alpha=0.15)
        plot_band(N, rs["deg_out_med"], rs["deg_out_iqr"], '^-', "Mean out-degree (median)", alpha=0.15)
        plt.xscale("log", base=2)
        plt.xticks(N, [fr"$2^{{{int(np.log2(n))}}}$" for n in N])
        plt.xlabel("N (elements)")
        plt.ylabel("Summary statistic")
        plt.title("Rate stability across size ladder (from run)")
        plt.legend()
        save_fig(os.path.join(tgt, "bd_rate_stability.png"))

# ---------- Paper II ----------
def paper2(rr):
    tgt = os.path.join(ROOT, "paper2", "figs")
    ensure_dir(tgt)

    # 1) BDG L2 error vs ell (log-log)
    e = rr.get("paper2", {}).get("bdg_error", None)  # {"ell":[...], "L2":[...], "slope": p}
    if e and all(k in e for k in ("ell","L2")):
        ell = np.array(e["ell"], dtype=float)
        L2 = np.array(e["L2"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.loglog(ell, L2, 'o-')
        if "slope" in e:
            plt.text(min(ell)*1.05, max(L2)/1.2, f"slope≈{e['slope']:.2f}")
        plt.xlabel("coarse length ℓ")
        plt.ylabel(r"$\|B_{\rm BDG}f - \Box f\|_{L^2}$")
        plt.title("BDG vs continuum error (from run)")
        save_fig(os.path.join(tgt, "bdg_error_curve.png"))

    # 2) Spectral dimension flow
    sd = rr.get("paper2", {}).get("spectral_dim", None)  # {"tau":[...], "ds":[...]}
    if sd and all(k in sd for k in ("tau","ds")):
        tau = np.array(sd["tau"], dtype=float)
        ds  = np.array(sd["ds"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.semilogx(tau, ds, 'o-')
        plt.axhline(4.0, ls='--')
        plt.xlabel("walk steps τ")
        plt.ylabel(r"$d_s(\tau)$")
        plt.title("Spectral dimension (from run)")
        save_fig(os.path.join(tgt, "spectral_flow.png"))

    # 3) Orientation invariance CV
    lor = rr.get("paper2", {}).get("lorentz_CV", None)  # {"angles":[...], "CV":[...]}
    if lor and all(k in lor for k in ("angles","CV")):
        a = np.array(lor["angles"], dtype=float)
        cv = np.array(lor["CV"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.plot(a, cv, 'o-')
        plt.xlabel("orientation index")
        plt.ylabel("CV (cross-edges)")
        plt.title("Orientation invariance (from run)")
        save_fig(os.path.join(tgt, "lorentz_bootstrap.png"))

    # 4) Exact BD curvature layers
    bd = rr.get("paper2", {}).get("bd_layers", None)
    if bd and "L" in bd and len(bd["L"])>=4:
        L = bd["L"]
        k = np.arange(1,5)
        plt.figure(figsize=(6,3.6))
        plt.bar(k, L[:4])
        plt.xticks(k, [f"{ki}" for ki in k])
        plt.xlabel("Layer k")
        plt.ylabel("Average L_k")
        plt.title("Exact BD curvature layers (from run)")
        save_fig(os.path.join(tgt, "bd_curvature_exact.png"))

# ---------- Paper III ----------
def paper3(rr):
    tgt = os.path.join(ROOT, "paper3", "figs")
    ensure_dir(tgt)

    # Friedmann residuals vs ell
    fr = rr.get("paper3", {}).get("friedmann_resid", None)  # {"ell":[...], "abs_resid":[...], "slope":...}
    if fr and all(k in fr for k in ("ell","abs_resid")):
        ell = np.array(fr["ell"], dtype=float)
        res = np.array(fr["abs_resid"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.loglog(ell, res, 'o-')
        if "slope" in fr:
            plt.text(min(ell)*1.05, max(res)/1.2, f"slope≈{fr['slope']:.2f}")
        plt.xlabel("coarse length ℓ")
        plt.ylabel("Friedmann abs residual")
        plt.title("FLRW residuals (from run)")
        save_fig(os.path.join(tgt, "friedmann_check.png"))

    # Poisson residuals vs ell
    pr = rr.get("paper3", {}).get("poisson_resid", None)  # {"ell":[...], "L2":[...]}
    if pr and all(k in pr for k in ("ell","L2")):
        ell = np.array(pr["ell"], dtype=float)
        L2  = np.array(pr["L2"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.loglog(ell, L2, 'o-')
        plt.xlabel("coarse length ℓ")
        plt.ylabel(r"$L^2$ error (Poisson)")
        plt.title("Poisson residuals (from run)")
        save_fig(os.path.join(tgt, "poisson_residuals.png"))

    # Wave dispersion
    wd = rr.get("paper3", {}).get("wave_dispersion", None)  # {"k":[...], "omega_disc":[...], "omega_cont":[...]}
    if wd and all(k in wd for k in ("k","omega_disc","omega_cont")):
        k = np.array(wd["k"], dtype=float)
        od = np.array(wd["omega_disc"], dtype=float)
        oc = np.array(wd["omega_cont"], dtype=float)
        plt.figure(figsize=(6,3.6))
        plt.plot(k, od, 'o', label="discrete")
        plt.plot(k, oc, '-', label="continuum")
        plt.xlabel("k")
        plt.ylabel("ω(k)")
        plt.legend()
        plt.title("Wave dispersion (from run)")
        save_fig(os.path.join(tgt, "dispersion_waves.png"))

def main():
    rr = load_run_report()
    paper1(rr)
    paper2(rr)
    paper3(rr)

if __name__ == "__main__":
    main()
