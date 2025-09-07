#!/usr/bin/env python3
import os, json, numpy as np
import matplotlib
# will use Agg if set by workflow; otherwise enforce
if os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/.."

def savefig(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"[ok] wrote {path}")

def section(rr, *keys):
    d = rr
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            print(f"[skip] missing key {'/'.join(keys)}")
            return None
        d = d[k]
    return d

def main():
    rr_path = os.path.join(ROOT, "run_report.json")
    if not os.path.exists(rr_path):
        print("[warn] run_report.json not found — nothing to do")
        return 0
    with open(rr_path) as f:
        rr = json.load(f)

    # ---------- Paper I ----------
    tgt = os.path.join(ROOT, "paper1", "figs"); os.makedirs(tgt, exist_ok=True)

    L = section(rr, "paper1", "interval_L")
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

    mm = section(rr, "paper1", "mm_dimension")
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

    rs = section(rr, "paper1", "rate_stability")
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

    # ---------- Paper II ----------
    tgt = os.path.join(ROOT, "paper2", "figs"); os.makedirs(tgt, exist_ok=True)

    e = section(rr, "paper2", "bdg_error")
    if e and all(k in e for k in ("ell","L2")):
        ell = np.array(e["ell"], float)
        L2 = np.array(e["L2"], float)
        plt.figure(figsize=(6,3.6))
        plt.loglog(ell, L2, 'o-')
        if "slope" in e:
            plt.text(min(ell)*1.05, max(L2)/1.2, f"slope≈{e['slope']:.2f}")
        plt.xlabel("coarse length ℓ"); plt.ylabel(r"$\|B_{\rm BDG}f - \Box f\|_{L^2}$")
        plt.title("BDG vs continuum error (from run)")
        savefig(os.path.join(tgt, "bdg_error_curve.png"))

    sd = section(rr, "paper2", "spectral_dim")
    if sd and all(k in sd for k in ("tau","ds")):
        tau = np.array(sd["tau"], float)
        ds  = np.array(sd["ds"], float)
        plt.figure(figsize=(6,3.6))
        plt.semilogx(tau, ds, 'o-')
        plt.axhline(4.0, ls='--')
        plt.xlabel("walk steps τ"); plt.ylabel(r"$d_s(\tau)$")
        plt.title("Spectral dimension (from run)")
        savefig(os.path.join(tgt, "spectral_flow.png"))

    lor = section(rr, "paper2", "lorentz_CV")
    if lor and all(k in lor for k in ("angles","CV")):
        a  = np.array(lor["angles"], float)
        cv = np.array(lor["CV"], float)
        plt.figure(figsize=(6,3.6))
        plt.plot(a, cv, 'o-')
        plt.xlabel("orientation index"); plt.ylabel("CV (cross-edges)")
        plt.title("Orientation invariance (from run)")
        savefig(os.path.join(tgt, "lorentz_bootstrap.png"))

    bd = section(rr, "paper2", "bd_layers")
    if bd and "L" in bd and len(bd["L"]) >= 4:
        L = bd["L"]; k = np.arange(1,5)
        plt.figure(figsize=(6,3.6))
        plt.bar(k, L[:4]); plt.xticks(k, [f"{ki}" for ki in k])
        plt.xlabel("Layer k"); plt.ylabel("Average L_k")
        plt.title("Exact BD curvature layers (from run)")
        savefig(os.path.join(tgt, "bd_curvature_exact.png"))

    # ---------- Paper III ----------
    tgt = os.path.join(ROOT, "paper3", "figs"); os.makedirs(tgt, exist_ok=True)

    fr = section(rr, "paper3", "friedmann_resid")
    if fr and all(k in fr for k in ("ell","abs_resid")):
        ell = np.array(fr["ell"], float)
        res = np.array(fr["abs_resid"], float)
        plt.figure(figsize=(6,3.6))
        plt.loglog(ell, res, 'o-')
        if "slope" in fr:
            plt.text(min(ell)*1.05, max(res)/1.2, f"slope≈{fr['slope']:.2f}")
        plt.xlabel("coarse length ℓ"); plt.ylabel("Friedmann abs residual")
        plt.title("FLRW residuals (from run)")
        savefig(os.path.join(tgt, "friedmann_check.png"))

    pr = section(rr, "paper3", "poisson_resid")
    if pr and all(k in pr for k in ("ell","L2")):
        ell = np.array(pr["ell"], float)
        L2  = np.array(pr["L2"], float)
        plt.figure(figsize=(6,3.6))
        plt.loglog(ell, L2, 'o-')
        plt.xlabel("coarse length ℓ"); plt.ylabel(r"$L^2$ error (Poisson)")
        plt.title("Poisson residuals (from run)")
        savefig(os.path.join(tgt, "poisson_residuals.png"))

    wd = section(rr, "paper3", "wave_dispersion")
    if wd and all(k in wd for k in ("k","omega_disc","omega_cont")):
        k  = np.array(wd["k"], float)
        od = np.array(wd["omega_disc"], float)
        oc = np.array(wd["omega_cont"], float)
        plt.figure(figsize=(6,3.6))
        plt.plot(k, od, 'o', label="discrete")
        plt.plot(k, oc, '-', label="continuum")
        plt.xlabel("k"); plt.ylabel("ω(k)"); plt.legend()
        plt.title("Wave dispersion (from run)")
        savefig(os.path.join(tgt, "dispersion_waves.png"))

    print("[done] figure generation completed")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
