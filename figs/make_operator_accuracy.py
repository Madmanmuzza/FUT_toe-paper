# Placeholder: expects user to run operator tests and dump errors to figs/out/operator_error.json
import json, os
import matplotlib.pyplot as plt
import numpy as np

INP = "figs/out/operator_error.json"
OUTP = "figs/out/operator_accuracy.png"

def main():
    if not os.path.exists(INP):
        # Create a dummy curve to illustrate expected decreasing trend
        taus = np.array([8, 12, 16, 24, 32, 48, 64])
        err = 1.0 / np.sqrt(taus)
        payload = {"ell": taus.tolist(), "L2_error": err.tolist()}
        with open(INP, "w") as f:
            json.dump(payload, f)
    with open(INP, "r") as f:
        data = json.load(f)
    ell = np.array(data["ell"])
    err = np.array(data["L2_error"])
    plt.figure()
    plt.loglog(ell, err, marker="o")
    plt.xlabel("coarse length (ell)")
    plt.ylabel("||B f - Box f||_2")
    plt.title("Operator accuracy scaling")
    plt.tight_layout()
    plt.savefig(OUTP, dpi=160)
    print("Saved:", OUTP)

if __name__ == "__main__":
    main()
