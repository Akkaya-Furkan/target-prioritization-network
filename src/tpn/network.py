from __future__ import annotations

import networkx as nx
import pandas as pd

from .scoring import DATATYPES


def build_graph(disease: str, df: pd.DataFrame, dt_scores: dict[str, dict[str, float]]) -> nx.Graph:
    """Tripartite graph: target -- evidence type -- disease. Edge weight = datatype score."""
    g = nx.Graph()
    g.add_node(disease, kind="disease")
    for d in DATATYPES:
        g.add_node(d, kind="evidence")
    for _, r in df.iterrows():
        g.add_node(r.symbol, kind="target", score=float(r.score_main))
        for d, s in dt_scores[r.symbol].items():
            if d in DATATYPES and s > 0:
                g.add_edge(r.symbol, d, weight=float(s))
    for d in DATATYPES:  # evidence type -> disease, weight = mean datatype score of linked targets
        ws = [g[d][t]["weight"] for t in g[d]]
        if ws:
            g.add_edge(d, disease, weight=sum(ws) / len(ws))
    return g


def centrality_table(g: nx.Graph) -> pd.DataFrame:
    deg = nx.degree_centrality(g)
    btw = nx.betweenness_centrality(g, weight=None)
    clo = nx.closeness_centrality(g)
    rows = [dict(symbol=n, degree_c=deg[n], betweenness=btw[n], closeness=clo[n])
            for n, k in g.nodes(data="kind") if k == "target"]
    return pd.DataFrame(rows)
