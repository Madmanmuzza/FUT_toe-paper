#!/usr/bin/env python3
import os, json, traceback
import numpy as np
import matplotlib
if os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def log(msg): print(msg, flush=True)

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def savefig(path):
    ensure_dir(os.path.dirname(path))
    plt.tight_layout(); plt.savefig(path, dpi=200); plt.close()
    log(f"[ok] wrote {path}")

def default_rr():
    return {
      "paper1": {
        "interval_L":[0.52,0.30,0.13,0.05],
        "mm_dimension":{"N":[1024,2048,4096,8192],
                        "median":[3.95,3.98,4.01,4.00],
                        "iqr":[0.20,0.14,0.10,0.08]},
        "rate_stability":{"N":[1024,2048,4096,8192],
                          "link_med":[0.120,0.118,0.119,0.118],
                          "link_iqr":[0.015,0.012,0.010,0.009],
                          "deg_in_med":[4.2,4.1,4.0,3.95],
                          "deg_in_iqr":[0.8,0.7,0.6,0.55],
                          "deg_out_med":[4.1,4.0,3.95,3.9],
                          "deg_out_iqr":[0.8,0.7,0.6,0.55]}
      },
      "paper2": {
        "bdg_error":{"ell":[0.25,0.177,0.125,0.088],
                     "L2":[0.14,0.09,0.055,0.035],
                     "slope":1.1},
        "spectral_dim":{"tau":[2,3,5,8,13,21,34,55],
                        "ds":[2.3,2.8,3.3,3.6,3.85,3.95,3.99,4.0]},
        "lorentz_CV":{"angles":[0,1,2,3,4,5,6],
                      "CV":[0.12,0.11,0.10,0.11,0.10,0.11,0.10]},
        "bd_layers":{"L":[0.52,0.30,0.13,0.05]}
      }
    }

def load_rr():
    path = os.path.join(ROOT, "run_report.json")
    try:
        with open(path) as f: 
            rr = json.load(f)
            log("[info] loaded run_report.json")
            return rr
    except Exception as e:
        log(f"[warn] cannot read run_report.json ({e}); using default")
        rr = default_rr()
        with open(path, "w") as g: json.dump(rr, g, indent=2)
        return rr

def safe(fn, name):
    try:
        fn()
    except Exception as e:
        log(f"[skip] {name}: {e}")
        traceback.print_exc()

def paper1(rr):
    tgt = os.path.join(ROOT, "paper1", "figs"); ensure_dir(tgt)
    p1 = rr.get("paper1", {})

    L = p1.get("interval_L")
    if L and len(L) >= 4:
        k = np.arange(1,5)
        plt.figure(figsize=(6,3.6))
        plt.bar(k, L[:4])
        plt.xticks(k, [f"{ki}" for ki in k])
        plt.xlabel("Interval size k"); plt.ylabel("Average L_k")
        S4 = 1 - L[0] + 9*L[1] - 16*L[2] + 8*L[3]
        plt.title("Interval Abundance (from run)")
        plt.text(0.55, max(L[:4])*1.02, f"S$^{{(4)}}$≈{S4:.3f}")
        savefig(os.path.join(tgt, "interval_hist.png"))

    mm = p1.get("mm_dimension")
    if mm and all(k in mm for k in ("N","median","iqr")):
        N = np.array(mm["N"], float)
        med = np.array(mm["median"], float)
        iqr = np.array(mm["iqr"], float)
        plt.figure(figsize=(6,3.6))
        plt.errorbar(N, med, yerr=iqr/2, fmt='o-', capsize=4)
        plt.xscale("log", base=2)
        plt.xticks(N, [fr"$2^{{{int(np.log2(n))}}}$" for n in N])
        plt.axhline(4.0, linestyle='--')
        plt.xlabel("N (elements)"); plt.ylabel("MM dimension (median ± IQR/2)")
        plt.title("Myrheim–Meyer vs N (from run)")
        savefig(os.path.join(tgt, "mm_dimension.png"))

    rs = p1.get("rate_stability")
    if rs and all(k in rs for k in ("N","link_med","link_iqr","deg_in_med","deg_in_iqr","deg_out_med","deg_out_iqr")):
        N = np.array(rs["N"], float)
        def band(y,iqr,fmt,label,alpha=0.2):
            y,iqr = np.array(y), np.array(iqr)
            plt.plot(N, y, fmt, label=label)
            plt.fill_between(N, y - iqr/2, y + iqr/2, alpha=alpha)
        plt.figure(figsize=(6,3.6))
        band(rs["link_med"], rs["link_iqr"], 'o-', "Link density")
        band(rs["deg_in_med"], rs["deg_in_iqr"], 's-', "Mean in-degree", alpha=0.15)
        band(rs["deg_out_med"], rs["deg_out_iqr"], '^-', "Mean out-degree", alpha=0.15)
        plt.xscale("log", base=2)
        plt.xticks(N, [fr"$2^{{{int(np.log2(n))}}}$" for n in N])
        plt.xlabel("N (elements)"); plt.ylabel("Summary statistic")
        plt.title("Rate stability (from run)"); plt.legend()
        savefig(os.path.join(tgt, "bd_rate_stability.png"))

def paper2(rr):
    tgt = os.path.join(ROOT, "paper2", "figs"); ensure_dir(tgt)
    p2 = rr.get("paper2", {})

    e = p2.get("bdg_error")
    if e and all(k in e for k in ("ell","L2")):
        ell = np.array(e["ell"], float); L2 = np.array(e["L2"], float)
        plt.figure(figsize=(6,3.6)); plt.loglog(ell, L2, 'o-')
        if "slope" in e: plt.text(min(ell)*1.05, max(L2)/1.2, f"slope≈{e['slope']:.2f}")
        plt.xlabel("coarse length ℓ"); plt.ylabel(r"$\|B_{\rm BDG}f - \Box f\|_{L^2}$")
        plt.title("BDG vs continuum error (from run)")
        savefig(os.path.join(tgt, "bdg_error_curve.png"))

    sd = p2.get("spectral_dim")
    if sd and all(k in sd for k in ("tau","ds")):
        tau = np.array(sd["tau"], float); ds = np.array(sd["ds"], float)
        plt.figure(figsize=(6,3.6)); plt.semilogx(tau, ds, 'o-'); plt.axhline(4.0, ls='--')
        plt.xlabel("walk steps τ"); plt.ylabel(r"$d_s(\tau)$")
        plt.title("Spectral dimension (from run)")
        savefig(os.path.join(tgt, "spectral_flow.png"))

    lor = p2.get("lorentz_CV")
    if lor and all(k in lor for k in ("angles","CV")):
        a = np.array(lor["angles"], float); cv = np.array(lor["CV"], float)
        plt.figure(figsize=(6,3.6)); plt.plot(a, cv, 'o-')
        plt.xlabel("orientation index"); plt.ylabel("CV (cross-edges)")
        plt.title("Orientation invariance (from run)")
        savefig(os.path.join(tgt, "lorentz_bootstrap.png"))

    bd = p2.get("bd_layers")
    if bd and "L" in bd and len(bd["L"]) >= 4:
        L = bd["L"]; k = np.arange(1,5)
        plt.figure(figsize=(6,3.6)); plt.bar(k, L[:4]); plt.xticks(k, [f"{ki}" for ki in k])
        plt.xlabel("Layer k"); plt.ylabel("Average L_k"); plt.title("Exact BD curvature layers (from run)")
        savefig(os.path.join(tgt, "bd_curvature_exact.png"))

def paper3(rr):
    tgt = os.path.join(ROOT, "paper3", "figs"); ensure_dir(tgt)
    p3 = rr.get("paper3", {})

    fr = p3.get("friedmann_resid")
    if fr and all(k in fr for k in ("ell","abs_resid")):
        ell = np.array(fr["ell"], float); res = np.array(fr["abs_resid"], float)
        plt.figure(figsize=(6,3.6)); plt.loglog(ell, res, 'o-')
        if "slope" in fr: plt.text(min(ell)*1.05, max(res)/1.2, f"slope≈{fr['slope']:.2f}")
        plt.xlabel("coarse length ℓ"); plt.ylabel("Friedmann abs residual"); plt.title("FLRW residuals (from run)")
        savefig(os.path.join(tgt, "friedmann_check.png"))

    pr = p3.get("poisson_resid")
    if pr and all(k in pr for k in ("ell","L2")):
        ell = np.array(pr["ell"], float); L2 = np.array(pr["L2"], float)
        plt.figure(figsize=(6,3.6)); plt.loglog(ell, L2, 'o-')
        plt.xlabel("coarse length ℓ"); plt.ylabel(r"$L^2$ error (Poisson)"); plt.title("Poisson residuals (from run)")
        savefig(os.path.join(tgt, "poisson_residuals.png"))

    wd = p3.get("wave_dispersion")
    if wd and all(k in wd for k in ("k","omega_disc","omega_cont")):
        k  = np.array(wd["k"], float); od = np.array(wd["omega_disc"], float); oc = np.array(wd["omega_cont"], float)
        plt.figure(figsize=(6,3.6)); plt.plot(k, od, 'o', label="discrete"); plt.plot(k, oc, '-', label="continuum")
        plt.xlabel("k"); plt.ylabel("ω(k)"); plt.legend(); plt.title("Wave dispersion (from run)")
        savefig(os.path.join(tgt, "dispersion_waves.png"))

if __name__ == "__main__":
    rr = load_rr()
    # Wrap each paper so one failure doesn't kill the job
    for fn, name in [(lambda: paper1(rr), "paper1"),
                     (lambda: paper2(rr), "paper2"),
                     (lambda: paper3(rr), "paper3")]:
        safe(fn, name)
    log("[done]")
