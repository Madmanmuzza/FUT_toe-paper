# Placeholder: expects chisq across d saved to figs/out/interval_chisq.json
import json, os
import matplotlib.pyplot as plt
import numpy as np

INP = "figs/out/interval_chisq.json"
OUTP = "figs/out/interval_distance.png"

def main():
    if not os.path.exists(INP):
        d = np.array([2,3,4,5,6])
        chi = (d-4.0)**2 + 0.5
        with open(INP, "w") as f:
            json.dump({"d": d.tolist(), "chisq": chi.tolist()}, f)
    with open(INP, "r") as f:
        data = json.load(f)
    d = np.array(data["d"])
    chi = np.array(data["chisq"])
    plt.figure()
    plt.plot(d, chi, marker="o")
    plt.xlabel("dimension d")
    plt.ylabel("chi^2 distance")
    plt.title("Interval abundance distance vs d")
    plt.tight_layout()
    plt.savefig(OUTP, dpi=160)
    print("Saved:", OUTP)

if __name__ == "__main__":
    main()
