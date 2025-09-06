from __future__ import annotations
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import numpy as np
import networkx as nx
from utils.causet import Causet

def ordering_fraction(causet: Causet) -> float:
    """Compute ordering fraction r = 2R/[N(N-1)] using reachability."""
    G = causet.dag
    nodes = list(G.nodes())
    N = len(nodes)
    if N < 2:
        return 0.0
    R = 0
    # Precompute ancestors/descendants
    anc = {u: nx.ancestors(G, u) for u in nodes}
    des = {u: nx.descendants(G, u) for u in nodes}
    # Count comparable unordered pairs once
    for i, u in enumerate(nodes):
        reach = anc[u] | des[u]
        for v in nodes[i+1:]:
            if v in reach:
                R += 1
    return 2.0 * R / (N * (N - 1))

def invert_r_to_d(r: float, table_r_to_d: Tuple[np.ndarray, np.ndarray]) -> float:
    """Monotone interpolation mapping ordering fraction -> Myrheim–Meyer dimension.
    Clamps r into the calibrated range to avoid failures when the process yields
    ordering fractions outside the sprinkling band."""
    r_vals, d_vals = table_r_to_d
    r_clamped = float(np.clip(r, float(r_vals.min()), float(r_vals.max())))
    return float(np.interp(r_clamped, r_vals, d_vals))

def myrheim_meyer_dimension(causet: Causet, table_r_to_d: Tuple[np.ndarray, np.ndarray]) -> float:
    r = ordering_fraction(causet)
    return invert_r_to_d(r, table_r_to_d)

def spectral_dimension(
    G: nx.Graph,
    taus: Sequence[int],
    n_walkers: int = 2000,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Estimate spectral dimension via lazy random-walk returns."""
    if rng is None:
        rng = np.random.default_rng()
    nodes = list(G.nodes())
    if not nodes:
        raise ValueError("empty graph")
    taus_arr = np.array(sorted(set(int(t) for t in taus if t > 0)))
    P = np.zeros_like(taus_arr, dtype=float)
    for idx, tau in enumerate(taus_arr):
        returns = 0
        for _ in range(n_walkers):
            s = nodes[rng.integers(0, len(nodes))]
            x = s
            for _ in range(tau):
                if rng.random() < 0.5:
                    pass
                else:
                    nbrs = list(G.neighbors(x))
                    if nbrs:
                        x = nbrs[rng.integers(0, len(nbrs))]
            if x == s:
                returns += 1
        P[idx] = returns / float(n_walkers)
    return taus_arr, P

def fit_spectral_dimension(taus: np.ndarray, P_tau: np.ndarray, candidate_windows: Sequence[Tuple[int, int]]) -> Tuple[float, Tuple[int, int]]:
    """Fit log P ~ -(d_s/2) log τ + c across candidate windows; choose by AIC."""
    if len(taus) != len(P_tau):
        raise ValueError("taus and P_tau length mismatch")
    log_tau = np.log(taus)
    log_P = np.log(P_tau + 1e-300)
    best_aic = np.inf
    best = (np.nan, (0, len(taus) - 1))
    for i0, i1 in candidate_windows:
        x = log_tau[i0:i1 + 1]
        y = log_P[i0:i1 + 1]
        if len(x) < 3 or not np.isfinite(y).all():
            continue
        A = np.vstack([x, np.ones_like(x)]).T
        m, b = np.linalg.lstsq(A, y, rcond=None)[0]
        residuals = y - (m * x + b)
        k = 2
        n = len(x)
        rss = float(np.dot(residuals, residuals))
        aic = n * np.log(rss / n + 1e-30) + 2 * k
        if aic < best_aic:
            best_aic = aic
            d_s = -2.0 * m
            best = (float(d_s), (i0, i1))
    return best
