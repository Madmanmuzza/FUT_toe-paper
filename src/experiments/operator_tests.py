import sys, os; sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import os, json, math
import numpy as np
import networkx as nx

from typing import Tuple, Dict, List
from utils.causet import Causet
from operators import bdg_operator

OUT = "figs/out"
os.makedirs(OUT, exist_ok=True)

def unit_ball_volume(n: int) -> float:
    # Volume of n-dim unit ball
    from math import pi, gamma
    return (pi ** (n / 2.0)) / gamma(n/2.0 + 1.0)

def diamond_volume(d: int, T: float) -> float:
    # Alexandrov diamond in (1+n)-dim Minkowski, with n=d-1
    n = d - 1
    Vn = unit_ball_volume(n)
    return 2.0 * Vn * (T/2.0) ** (n + 1) / (n + 1)

def sample_in_diamond(T: float, n_space: int, N: int, rng: np.random.Generator) -> np.ndarray:
    pts_t = []
    while len(pts_t) < N:
        t = rng.uniform(-T/2, T/2)
        R = (T/2) - abs(t)
        if R <= 0:
            continue
        if rng.random() < (R / (T/2)) ** n_space:
            pts_t.append((t, R))
    t_list = np.array([t for t, _ in pts_t])
    R_list = np.array([R for _, R in pts_t])
    X = []
    for t, R in zip(t_list, R_list):
        if n_space == 0:
            X.append(np.array([t], dtype=float))
            continue
        # radius r = R * U^{1/n}
        u = rng.uniform()
        r = R * (u ** (1.0 / n_space))
        g = rng.normal(size=n_space)
        g /= np.linalg.norm(g)
        vec = r * g
        X.append(np.concatenate(([t], vec)))
    return np.asarray(X)

def partial_order_from_points(X: np.ndarray) -> nx.DiGraph:
    N, d = X.shape
    G = nx.DiGraph()
    G.add_nodes_from(range(N))
    ts = X[:, 0]
    xs = X[:, 1:]
    for i in range(N):
        for j in range(N):
            if ts[j] > ts[i]:
                if ts[j] - ts[i] >= np.linalg.norm(xs[j] - xs[i]):
                    G.add_edge(i, j)
    G = nx.algorithms.dag.transitive_reduction(G)
    return G

def gaussian_probe_and_box(X: np.ndarray, sigma: float, d: int) -> Tuple[np.ndarray, np.ndarray]:
    """f = exp(-sigma*(t^2 + |x|^2)); Minkowski box = d/dt^2 - sum d/dx_i^2"""
    t = X[:, 0]
    x = X[:, 1:]
    r2_t = t**2
    r2_x = np.sum(x*x, axis=1)
    f = np.exp(-sigma*(r2_t + r2_x))
    n = d - 1
    box_f = (4*sigma*sigma*(r2_t - r2_x) + 2*sigma*(n - 1)) * f
    return f, box_f

def l2_error(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y_true - y_pred)**2)))

def main():
    d = 4
    n_space = d - 1
    T = 1.0
    N_list = [250, 350, 500]
    rng = np.random.default_rng(12345)
    sigma = 0.8
    ells = []
    errs = []
    for N in N_list:
        X = sample_in_diamond(T, n_space, N, rng)
        G = partial_order_from_points(X)
        C = Causet(G, field={i: 0.0 for i in range(G.number_of_nodes())})
        f_vals, box_f_vals = gaussian_probe_and_box(X, sigma=sigma, d=d)
        f_map = {i: float(f_vals[i]) for i in range(len(f_vals))}
        # Use ell ~ (Vol/N)^{1/d}
        Vol = diamond_volume(d, T)
        ell = (Vol / N) ** (1.0 / d)
        Bf_map = bdg_operator(C, f_map, d=d, ell=ell)  # dict node -> value
        Bf = np.array([Bf_map[i] for i in range(len(f_vals))], dtype=float)
        err = l2_error(box_f_vals, Bf)
        ells.append(float(ell)); errs.append(float(err))
    out = {"ell": ells, "L2_error": errs}
    with open(os.path.join(OUT, "operator_error.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("Saved figs/out/operator_error.json with", len(ells), "points.")

if __name__ == "__main__":
    main()
