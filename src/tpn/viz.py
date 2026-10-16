from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from pyvis.network import Network

COLORS = {"target": "#2a7ab0", "evidence": "#e08a1e", "disease": "#b03a48"}


def make_all(g: nx.Graph, df: pd.DataFrame, fig_dir: Path, meta: dict, top: int = 25):
    fig_dir.mkdir(exist_ok=True, parents=True)
    d = df.sort_values("rank_main")
    keep = set(d.head(top)["symbol"])
    sub = g.subgraph([n for n, k in g.nodes(data="kind") if k != "target" or n in keep]).copy()
    sub.remove_nodes_from([n for n in list(sub) if sub.degree(n) == 0])
    _network_png(sub, fig_dir / "network_top.png", meta)
    _score_png(d.head(15), fig_dir / "score_vs_publications.png", meta)
    _pyvis(sub, fig_dir.parent / "network.html", meta)


def _network_png(g, path, meta):
    fig, ax = plt.subplots(figsize=(11, 8))
    pos = nx.spring_layout(g, seed=7, k=1.4, weight=None, iterations=200)
    for kind, c in COLORS.items():
        ns = [n for n, k in g.nodes(data="kind") if k == kind]
        sizes = [60 + 500 * g.nodes[n].get("score", 0) if kind == "target" else 700 for n in ns]
        nx.draw_networkx_nodes(g, pos, nodelist=ns, node_color=c, node_size=sizes, alpha=.9, ax=ax, label=kind)
    nx.draw_networkx_edges(g, pos, alpha=.25, ax=ax)
    nx.draw_networkx_labels(g, pos, font_size=6.5, ax=ax)
    ax.set_title(f"{meta['disease']}: top targets - evidence type - disease (Open Targets {_ver(meta)})")
    ax.legend(scatterpoints=1, loc="lower left", markerscale=.4); ax.axis("off")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def _score_png(d, path, meta):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    d = d.dropna(subset=["pubmed_disease"])
    sc = ax.scatter(d["pubmed_disease"].astype(float) + 1, d["score_main"], s=40 + 12 * d["diversity"] * 6,
                    c=d["n_ligands"].apply(lambda x: x), cmap="viridis", alpha=.85)
    for _, r in d.iterrows():
        ax.annotate(r.symbol, (float(r.pubmed_disease) + 1, r.score_main), fontsize=8, xytext=(4, 3), textcoords="offset points")
    ax.set_xscale("log"); ax.set_xlabel("PubMed records: gene AND disease (+1, log)")
    ax.set_ylabel("main score (equal weights)")
    fig.colorbar(sc, label="ChEMBL ligands with pChEMBL")
    ax.set_title("Score vs. publication volume (top 15)"); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def _ver(meta):
    o = meta["open_targets"]["dataVersion"]
    return f"{o['year']}.{o['month']}"


def _pyvis(g, path, meta):
    net = Network(height="850px", width="100%", bgcolor="#ffffff", cdn_resources="in_line")
    for n, a in g.nodes(data=True):
        net.add_node(n, label=n, color=COLORS[a["kind"]], title=f"{a['kind']}" + (f" | score {a['score']:.3f}" if "score" in a else ""),
                     size=12 + (60 * a.get("score", 0) if a["kind"] == "target" else 18))
    for u, v, a in g.edges(data=True):
        net.add_edge(u, v, value=a["weight"], title=f"{a['weight']:.2f}")
    net.write_html(str(path))
