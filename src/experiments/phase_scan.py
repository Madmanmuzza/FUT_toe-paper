import sys, os; sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import yaml, os, json
from typing import Dict, List, Tuple
import networkx as nx
from utils.causet import Causet
from dimension import ordering_fraction, invert_r_to_d, spectral_dimension, fit_spectral_dimension

OUT = "figs/out"
os.makedirs(OUT, exist_ok=True)

def negentropic_fold_local(G: nx.DiGraph, psi: Dict[int, float], x: int, h: int = 2) -> float:
    """Local fold: alpha * mean(psi^2) - gamma * mean(log(1+psi^2)) over height-h neighborhood."""
    # BFS up to h steps in undirected view to approximate height-h Hasse neighborhood
    U = G.to_undirected()
    import collections
    q = collections.deque([(x, 0)])
    seen = {x}
    vals = []
    while q:
        u, d = q.popleft()
        vals.append(psi[u])
        if d == h:
            continue
        for v in U.neighbors(u):
            if v not in seen:
                seen.add(v)
                q.append((v, d+1))
    import math
    arr = np.array(vals, dtype=float)
    return float(arr.dot(arr) / len(arr) - np.mean(np.log1p(arr*arr)))

def simulate_trajectory(alpha: float, gamma: float, sigma_eta: float, target_N: int, seed: int, h: int = 2):
    """Very simple placeholder simulator to reach target_N using normalized/saturated rates.

    NOTE: Replace with your high-performance implementation. This version grows from a single element.
    """
    rng = np.random.default_rng(seed)
    # start with a single node
    G = nx.DiGraph()
    G.add_node(0)
    psi = {0: float(rng.normal())}
    def boundary_nodes():
        # treat maximal nodes as boundary
        max_nodes = [u for u in G.nodes if G.out_degree(u) == 0]
        return max_nodes
    def extremal_nodes():
        mins = [u for u in G.nodes if G.in_degree(u) == 0]
        maxs = [u for u in G.nodes if G.out_degree(u) == 0]
        return mins + maxs

    def F_loc(u):
        # scaled by alpha,gamma at call site
        base = negentropic_fold_local(G, psi, u, h=h)
        return alpha * (base + gamma) - gamma * np.log1p((base + 1.0)**2)  # placeholder shaping

    Lb, Ld = 1.0, 1.0
    k_b, k_d = 1.0, 1.0

    while G.number_of_nodes() < target_N:
        bnodes = boundary_nodes()
        enodes = extremal_nodes()
        lam_b = []
        for u in bnodes:
            eta = rng.normal(0.0, sigma_eta)
            lam_b.append((u, (k_b / (1+len(bnodes))) * min(max(F_loc(u)+eta, 0.0), Lb)))
        lam_d = []
        for v in enodes:
            eta = rng.normal(0.0, sigma_eta)
            lam_d.append((v, (k_d / (1+len(enodes))) * min(max(-F_loc(v)+eta, 0.0), Ld)))

        total = sum(x for _, x in lam_b) + sum(x for _, x in lam_d)
        if total <= 1e-12:
            # fallback: force a birth on a random boundary node
            u = bnodes[rng.integers(0, len(bnodes))]
            new = max(G.nodes) + 1
            G.add_edge(u, new)
            psi[new] = float(rng.normal(np.mean([psi.get(w, 0.0) for w in G.predecessors(new)] or [0.0]), 1.0))
            continue

        # choose event
        r = rng.random() * total
        acc = 0.0
        chosen = None
        mode = None
        for u, w in lam_b:
            acc += w
            if r <= acc:
                chosen = u
                mode = "birth"
                break
        if chosen is None:
            for v, w in lam_d:
                acc += w
                if r <= acc:
                    chosen = v
                    mode = "death"
                    break
        if mode == "birth":
            # choose ancestors: take down-set closure of predecessors within h steps
            preds = set()
            frontier = {chosen}
            for _ in range(h):
                new_frontier = set()
                for z in frontier:
                    new_frontier.update(G.predecessors(z))
                preds.update(frontier)
                frontier = new_frontier
            preds = set(preds)
            new = max(G.nodes) + 1
            for a in preds:
                G.add_edge(a, new)
            # transitive reduction to keep Hasse diagram minimal
            # (skipped) transitive reduction for speed
            psi[new] = float(rng.normal(np.mean([psi.get(a, 0.0) for a in preds] or [0.0]), 1.0))
        elif mode == "death":
            # remove only if extremal (safe)
            if G.in_degree(chosen) == 0 or G.out_degree(chosen) == 0:
                G.remove_node(chosen)
                if len(G) == 0:
                    G.add_node(0)
                # renumber? keep sparse ids
        else:
            pass

    return Causet(G, psi)

def diagnostics(causet: Causet, r_to_d_npz: str) -> Dict[str, float]:
    data = np.load(r_to_d_npz)
    r_vals = data["r_vals"]; d_vals = data["d_vals"]
    r = ordering_fraction(causet)
    d_mm = invert_r_to_d(r, (r_vals, d_vals))
    taus = np.unique((np.geomspace(4, 512, num=12)).astype(int))
    t, P = spectral_dimension(causet.undirected(), taus, n_walkers=1000)
    windows = [(i, i + 5) for i in range(0, len(t) - 5)]
    d_s, (i0, i1) = fit_spectral_dimension(t, P, windows)
    return {"r": float(r), "d_MM": float(d_mm), "d_s": float(d_s), "tau_window": float(t[i1] - t[i0])}

def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--config", default="config/params.yaml")
    args = ap.parse_args()
    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
    alpha_list = cfg["grid"]["alpha"]
    gamma_list = cfg["grid"]["gamma"]
    sigma_eta = cfg["grid"]["sigma_eta"]
    target_N = cfg["grid"]["target_N"]
    n_seeds = cfg["grid"]["n_seeds"]
    # calibration table
    r_to_d_npz = "bench/out/r_to_d.npz"
    if not os.path.exists(r_to_d_npz):
        raise SystemExit("Run `make bench` first to create bench/out/r_to_d.npz")
    results = []
    seeds = list(range(n_seeds))
    for a in alpha_list:
        for g in gamma_list:
            vals = []
            for s in seeds:
                C = simulate_trajectory(a, g, sigma_eta, target_N, seed=s)
                vals.append(diagnostics(C, r_to_d_npz))
            # aggregate
            keys = vals[0].keys()
            arr = {k: np.array([v[k] for v in vals], dtype=float) for k in keys}
            med = {f"med_{k}": float(np.median(arr[k])) for k in keys}
            iqr = {f"iqr_{k}": float(np.subtract(*np.percentile(arr[k], [75, 25]))) for k in keys}
            rec = {"alpha": a, "gamma": g}
            rec.update(med); rec.update(iqr)
            results.append(rec)
            print("grid", a, g, "=>", rec)
    out_json = os.path.join(OUT, "phase_scan.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved:", out_json)

if __name__ == "__main__":
    main()
