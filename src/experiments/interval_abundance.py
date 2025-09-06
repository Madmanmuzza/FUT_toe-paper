import sys, os; sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import os, json, numpy as np
import networkx as nx
from typing import Dict, Tuple, List
from utils.causet import Causet

OUT = "figs/out"
os.makedirs(OUT, exist_ok=True)

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

def interval_abundances(G: nx.DiGraph, k_max: int = 6) -> List[int]:
    nodes = list(G.nodes())
    anc = {u: nx.ancestors(G, u) for u in nodes}
    des = {u: nx.descendants(G, u) for u in nodes}
    counts = [0]*(k_max+1)  # counts[k] for k=0..k_max (k=|I(y,x)|)
    for x in nodes:
        Ax = anc[x]
        for y in Ax:
            k = len(des[y].intersection(Ax))
            if 0 <= k <= k_max:
                counts[k] += 1
    return counts

def chisq_distance(a: np.ndarray, b: np.ndarray) -> float:
    # Normalize to probabilities; use small epsilon to avoid div by 0
    eps = 1e-12
    pa = a / max(a.sum(), eps)
    pb = b / max(b.sum(), eps)
    return float(np.sum((pa - pb)**2 / (pb + eps)))

def main():
    rng = np.random.default_rng(2025)
    T = 1.0
    N = 500
    k_max = 6
    # Target "observed" causet at d=4
    X = sample_in_diamond(T, n_space=3, N=N, rng=rng)
    G_target = partial_order_from_points(X)
    obs = np.array(interval_abundances(G_target, k_max=k_max), dtype=float)

    chis = []
    dims = [2,3,4,5,6]
    for d in dims:
        n_space = d - 1
        # Benchmark causet for this d (same N)
        Xb = sample_in_diamond(T, n_space=n_space, N=N, rng=rng)
        Gb = partial_order_from_points(Xb)
        ref = np.array(interval_abundances(Gb, k_max=k_max), dtype=float)
        chis.append(chisq_distance(obs, ref))
    out = {"d": dims, "chisq": [float(c) for c in chis]}
    with open(os.path.join(OUT, "interval_chisq.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("Saved figs/out/interval_chisq.json")

if __name__ == "__main__":
    main()
