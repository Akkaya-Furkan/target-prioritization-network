"""Offline tests: recorded-style fake responses, no network."""
import json
from pathlib import Path

import pandas as pd

from tpn.clients import ChEMBL, OpenTargets, PubMed
from tpn.network import build_graph, centrality_table
from tpn.pipeline import collect, pubmed_counts
from tpn.scoring import add_scores, diversity, sensitivity, under_studied

FIX = json.loads((Path(__file__).parent / "fixtures" / "responses.json").read_text())


class Resp:
    def __init__(self, data, status=200):
        self._d, self.status_code = data, status
    def json(self): return self._d
    def raise_for_status(self): pass


class FakeHttp:
    def __init__(self): self.calls = []
    def post(self, url, json=None, **kw):
        self.calls.append(url); return Resp({"data": FIX["ot_disease"]})
    def get(self, url, params=None, **kw):
        self.calls.append(url)
        if "esearch" in url:
            return Resp({"esearchresult": {"count": FIX["pubmed"][params["term"].split('"')[1]] if 'AND' not in params["term"] else 5}})
        path = url.rsplit("/", 1)[1]
        return Resp(FIX["chembl"][path.replace(".json", "")])


def test_diversity_excludes_known_drug():
    assert diversity({"known_drug": 1.0}) == 0
    assert diversity({"genetic_association": .5, "literature": .1, "known_drug": .9}) == 2 / 6


def test_collect_and_score(monkeypatch):
    monkeypatch.setattr("time.sleep", lambda s: None)
    http = FakeHttp()
    df, dts = collect(OpenTargets(http), ChEMBL(http), "EFO_X", 2, log=lambda *_: None)
    assert list(df.symbol) == ["AAA1", "BBB2"]
    assert df.loc[0, "n_ligands"] == 2 and df.loc[0, "max_phase"] == 4.0
    assert df.loc[1, "n_ligands"] == 0 and df.loc[1, "max_phase"] == 0.0   # no uniprot -> no chembl target
    s = add_scores(df)
    assert s.loc[0, "rank_main"] == 1
    assert (s[["score_main", "score_assoc_heavy", "score_ligand_heavy"]] > 0).all().all()


def test_sensitivity_and_pubmed(monkeypatch):
    monkeypatch.setattr("time.sleep", lambda s: None)
    n = 14
    df = pd.DataFrame({"symbol": [f"G{i}" for i in range(n)], "assoc": [1 - i / 20 for i in range(n)],
                       "diversity": [(i % 6 + 1) / 6 for i in range(n)], "n_ligands": [i * 10 for i in range(n)][::-1]})
    s = add_scores(df)
    sen = sensitivity(s, 5)
    assert set(sen) == {"assoc_heavy", "ligand_heavy"} and sen["assoc_heavy"]["k"] == 5
    assert all(-1 <= v["spearman"] <= 1 for v in sen.values())
    http = FakeHttp(); http.get = lambda url, params=None, **kw: Resp({"esearchresult": {"count": "7"}})
    p = pubmed_counts(PubMed(http), s, "x", 3)
    assert p["pubmed_disease"].notna().sum() == 3
    assert isinstance(under_studied(p.assign(pubmed_disease=p.pubmed_disease.fillna(10**6)), 5), pd.DataFrame)


def test_graph_centrality():
    df = pd.DataFrame({"symbol": ["A", "B"], "score_main": [.5, .4]})
    g = build_graph("D", df, {"A": {"genetic_association": .8, "known_drug": .9}, "B": {"literature": .3}})
    assert g.has_edge("A", "genetic_association") and not g.has_edge("A", "known_drug")
    assert g.has_edge("genetic_association", "D")
    c = centrality_table(g)
    assert set(c.symbol) == {"A", "B"}
