import numpy as np
import yaml
import os
from typing import Tuple, Dict, List
import networkx as nx

OUT = "bench/out"
os.makedirs(OUT, exist_ok=True)

def sample_in_diamond(T: float, n_space: int, N: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform sample in an Alexandrov diamond in (1+n_space)D Minkowski.
    Rejection sample t ~ Uniform(-T/2, T/2) with acceptance prob ((T/2 - |t|)/(T/2))^n,
    then sample spatial vector uniformly in n-ball of radius R = T/2 - |t|.
    """
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
    # Sample spatial points in n-ball
    X = []
    for t, R in zip(t_list, R_list):
        # Sample radius with PDF ∝ r^{n-1} on [0,R]: r = R * U^{1/n}
        u = rng.uniform()
        r = R * (u ** (1.0 / n_space)) if n_space > 0 else 0.0
        # Random direction on S^{n-1}
        if n_space == 0:
            vec = np.zeros(0, dtype=float)
        else:
            g = rng.normal(size=n_space)
            g /= np.linalg.norm(g)
            vec = r * g
        X.append(np.concatenate(([t], vec)))
    return np.asarray(X)

def partial_order_from_points(X: np.ndarray) -> nx.DiGraph:
    """Create DAG: p ≺ q iff t_q - t_p >= ||x_q - x_p||_2."""
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
    # Transitive reduction (Hasse diagram)
    G = nx.algorithms.dag.transitive_reduction(G)
    return G

def ordering_fraction_from_graph(G: nx.DiGraph) -> float:
    nodes = list(G.nodes())
    N = len(nodes)
    if N < 2:
        return 0.0
    R = 0
    anc = {u: nx.ancestors(G, u) for u in nodes}
    des = {u: nx.descendants(G, u) for u in nodes}
    for i, u in enumerate(nodes):
        reach = anc[u] | des[u]
        for v in nodes[i+1:]:
            if v in reach:
                R += 1
    return 2.0 * R / (N * (N - 1))

def calibrate_r_table(dims: List[int], Ns: List[int], seeds: List[int], T: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """Return monotone arrays (r_vals, d_vals) for d ∈ dims by extrapolating r(N) in 1/N."""
    r_points = []
    for d in dims:
        n_space = d - 1
        r_by_N = []
        invN = []
        for N in Ns:
            r_runs = []
            for s in seeds:
                rng = np.random.default_rng(s + 1000*d + N)
                X = sample_in_diamond(T=T, n_space=n_space, N=N, rng=rng)
                G = partial_order_from_points(X)
                r_runs.append(ordering_fraction_from_graph(G))
            r_by_N.append(np.mean(r_runs))
            invN.append(1.0 / N)
        # Linear fit r(N) ≈ r_d + a*(1/N)
        A = np.vstack([np.ones_like(invN), invN]).T
        r_d, a = np.linalg.lstsq(A, np.array(r_by_N), rcond=None)[0]
        r_points.append((float(r_d), float(d)))
    # Sort by r to ensure monotone mapping
    r_points.sort(key=lambda t: t[0])
    r_vals = np.array([t[0] for t in r_points], dtype=float)
    d_vals = np.array([t[1] for t in r_points], dtype=float)
    return r_vals, d_vals

def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--config", default="config/params.yaml")
    args = ap.parse_args()
    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
    dims = cfg["sprinkling"]["dims"]
    Ns = cfg["sprinkling"]["Ns"]
    seeds = cfg["sprinkling"]["seeds"]
    T = cfg["sprinkling"]["T"]
    r_vals, d_vals = calibrate_r_table(dims, Ns, seeds, T=T)
    out_npz = os.path.join(OUT, "r_to_d.npz")
    np.savez(out_npz, r_vals=r_vals, d_vals=d_vals)
    print("Saved:", out_npz)

if __name__ == "__main__":
    main()
