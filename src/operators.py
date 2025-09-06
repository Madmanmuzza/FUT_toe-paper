from __future__ import annotations
from typing import Dict, List, Tuple, Optional
import math
import networkx as nx
import numpy as np
from utils.causet import Causet

def _S_sphere(n: int) -> float:
    return 2.0 * (math.pi ** ((n + 1) / 2.0)) / math.gamma((n + 1) / 2.0)

def _c_d(d: int) -> float:
    return _S_sphere(d - 2) / (d * (d - 1) * (2.0 ** (d / 2.0 - 1.0)))

def _alpha_beta(d: int) -> Tuple[float, float]:
    cd_pow = _c_d(d) ** (2.0 / d)
    if d % 2 == 0:
        alpha = -2.0 * cd_pow / math.gamma(1.0 + 2.0 / d)
        beta  = (2.0 * math.gamma(d/2.0 + 2.0) * math.gamma(d/2.0 + 1.0) /
                 (math.gamma(2.0 / d) * math.gamma(d))) * cd_pow
    else:
        alpha = -1.0 * cd_pow / math.gamma(1.0 + 2.0 / d)
        beta  = (math.gamma((d + 3.0)/2.0) * math.gamma((d + 1.0)/2.0) /
                 (math.gamma(2.0 / d) * math.gamma(d))) * cd_pow
    return alpha, beta

def _C_coeffs(d: int) -> List[float]:
    n_d = (d // 2) + 2
    C: List[float] = []
    if d % 2 == 0:
        denom = math.gamma(d/2.0 + 2.0)
        for i in range(1, n_d + 1):
            s = 0.0
            for k in range(0, i):
                s += (math.comb(i-1, k) * ((-1.0) ** k) *
                      math.gamma((d/2.0)*(k+1) + 2.0) /
                      (denom * math.gamma(1.0 + (d/2.0)*k)))
            C.append(s)
    else:
        denom = math.gamma((d + 3.0)/2.0)
        for i in range(1, n_d + 1):
            s = 0.0
            for k in range(0, i):
                s += (math.comb(i-1, k) * ((-1.0) ** k) *
                      math.gamma((d/2.0)*(k+1) + 1.5) /
                      (denom * math.gamma(1.0 + (d/2.0)*k)))
            C.append(s)
    return C

def get_bdg_coefficients(d: int) -> Tuple[float, float, List[float]]:
    if d < 2:
        raise ValueError("d must be >= 2")
    a, b = _alpha_beta(d)
    C = _C_coeffs(d)
    return a, b, C

def bdg_layers_by_interval_size(causet: Causet, x: int, max_depth: int) -> List[List[int]]:
    """Compute past layers L_i(x) = { y ≺ x : |I(y,x)| = i-1 } up to max_depth.

    Uses set intersections of descendants/ancestors to count |I(y,x)|.
    """
    G = causet.dag
    layers: List[List[int]] = [[] for _ in range(max_depth + 1)]
    anc_x = nx.ancestors(G, x)
    # Precompute ancestors/descendants sets
    anc = {u: nx.ancestors(G, u) for u in G.nodes}
    des = {u: nx.descendants(G, u) for u in G.nodes}
    for y in anc_x:
        interval_nodes = des[y].intersection(anc_x)
        k = len(interval_nodes)  # |I(y,x)|
        if 0 <= k <= max_depth:
            layers[k].append(y)
    return layers

def bdg_operator(
    causet: Causet,
    f: Dict[int, float],
    d: int,
    ell: float,
    x: Optional[int] = None,
) -> Dict[int, float] | float:
    """Apply BDG operator at node x or all nodes."""
    alpha_d, beta_d, C = get_bdg_coefficients(d)
    max_i = (d // 2) + 2
    if len(C) < max_i:
        raise ValueError("insufficient BDG coefficient list")

    def apply_at(node: int) -> float:
        layers = bdg_layers_by_interval_size(causet, node, max_depth=max_i)
        acc = alpha_d * f[node]
        for i in range(1, max_i + 1):  # i=1..max_i corresponds to |I(y,x)|=i-1
            acc += beta_d * C[i - 1] * sum(f[y] for y in layers[i - 1])
        return acc / (ell ** 2)

    if x is not None:
        return apply_at(x)
    out: Dict[int, float] = {}
    for node in causet.dag.nodes():
        out[node] = apply_at(node)
    return out
