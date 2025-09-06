from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional
import networkx as nx

@dataclass
class Causet:
    """Minimal causet wrapper.

    Attributes
    ----------
    dag : nx.DiGraph
        Transitive reduction of the causal set (Hasse diagram). Must be acyclic.
    field : Dict[int, float]
        Map from node -> ψ value.
    """
    dag: nx.DiGraph
    field: Dict[int, float]

    @property
    def N(self) -> int:
        return self.dag.number_of_nodes()

    def undirected(self) -> nx.Graph:
        """Undirected skeleton used for diffusion/spectral dimension."""
        return self.dag.to_undirected(as_view=False)
