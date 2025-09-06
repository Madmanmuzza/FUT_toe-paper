import json, os
import matplotlib.pyplot as plt
import numpy as np

INP = "figs/out/phase_scan.json"
OUTP = "figs/out/phase_diagram.png"

def main():
    with open(INP, "r") as f:
        data = json.load(f)
    alphas = sorted(set(d["alpha"] for d in data))
    gammas = sorted(set(d["gamma"] for d in data))
    A = {a:i for i,a in enumerate(alphas)}
    G = {g:j for j,g in enumerate(gammas)}
    Z = np.full((len(alphas), len(gammas)), np.nan)
    for rec in data:
        i, j = A[rec["alpha"]], G[rec["gamma"]]
        # encode survival/plateau proxy as med_d_MM closeness to 4.0
        Z[i,j] = -abs(rec.get("med_d_MM", np.nan) - 4.0)
    plt.figure()
    plt.imshow(Z.T, origin="lower", extent=[min(alphas), max(alphas), min(gammas), max(gammas)], aspect="auto")
    plt.colorbar(label="-(|median d_MM - 4|)")
    plt.xlabel("alpha")
    plt.ylabel("gamma")
    plt.title("Phase diagram (proxy via d_MM)")
    plt.tight_layout()
    plt.savefig(OUTP, dpi=160)
    print("Saved:", OUTP)

if __name__ == "__main__":
    main()
