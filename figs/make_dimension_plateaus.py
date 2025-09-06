import json, os
import matplotlib.pyplot as plt
import numpy as np

INP = "figs/out/phase_scan.json"
OUTP = "figs/out/dimension_plateaus.png"

def main():
    with open(INP, "r") as f:
        data = json.load(f)
    # pick the best grid point (closest med_d_MM to 4)
    best = min(data, key=lambda rec: abs(rec.get("med_d_MM", 1e9) - 4.0))
    # Mock-up series: we don't store scale dependence in JSON; plot summary bars
    labels = ["d_MM (median)", "d_s (median)"]
    values = [best.get("med_d_MM", np.nan), best.get("med_d_s", np.nan)]
    errs = [0.5*best.get("iqr_d_MM", 0.0), 0.5*best.get("iqr_d_s", 0.0)]
    x = np.arange(len(labels))
    plt.figure()
    plt.bar(x, values, yerr=errs)
    plt.xticks(x, labels)
    plt.axhline(4.0, linestyle="--")
    plt.ylim(0, max(6, np.nanmax(values)+1))
    plt.title(f"Dimension plateaus (best grid: alpha={best['alpha']}, gamma={best['gamma']})")
    plt.tight_layout()
    plt.savefig(OUTP, dpi=160)
    print("Saved:", OUTP)

if __name__ == "__main__":
    main()
