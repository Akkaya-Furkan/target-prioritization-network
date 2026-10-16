"""End-to-end run: disease choice -> Open Targets -> ChEMBL -> score -> PubMed -> network/figures."""
from __future__ import annotations

import datetime as dt
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from .clients import ChEMBL, OpenTargets, PubMed
from .network import build_graph, centrality_table
from .scoring import DATATYPES, WEIGHTS, add_scores, diversity, sensitivity, under_studied

CANDIDATES = ["Parkinson disease", "rheumatoid arthritis", "ulcerative colitis",
              "amyotrophic lateral sclerosis", "psoriasis", "asthma"]
POOL = 200       # targets taken from Open Targets (by overall association score)
PUBMED_N = 30    # top targets (main score) queried in PubMed


def select_disease(ot: OpenTargets, names=CANDIDATES, n=200) -> tuple[dict, pd.DataFrame]:
    """Pick the disease with the densest multi-type evidence among its top-n associated targets."""
    rows = []
    for name in names:
        hit = ot.find_disease(name)
        d = ot.disease_targets(hit["id"], n)
        r = d["associatedTargets"]["rows"]
        nd = [diversity({x["id"]: x["score"] for x in t["datatypeScores"]}) * len(DATATYPES) for t in r]
        rows.append(dict(disease=hit["name"], efo=hit["id"], n_assoc_targets=d["associatedTargets"]["count"],
                         top_n=len(r), mean_datatypes=sum(nd) / len(nd), n_ge3_datatypes=sum(1 for x in nd if x >= 3),
                         mean_assoc=sum(t["score"] for t in r) / len(r)))
    df = pd.DataFrame(rows).sort_values(["n_ge3_datatypes", "n_assoc_targets"], ascending=False).reset_index(drop=True)
    return df.iloc[0].to_dict(), df


def collect(ot, chembl, efo: str, pool: int, log=print, workers: int = 8) -> tuple[pd.DataFrame, dict]:
    d = ot.disease_targets(efo, pool)
    trows = d["associatedTargets"]["rows"]

    def one(t):
        tg = t["target"]
        sc = {x["id"]: x["score"] for x in t["datatypeScores"]}
        up = next((p["id"] for p in tg["proteinIds"] if p["source"] == "uniprot_swissprot"), None)
        ct = chembl.single_protein_target(up) if up else None
        lig = chembl.ligands(ct) if ct else {"n_ligands": 0, "n_activities": 0, "truncated": False}
        ph = chembl.max_phase(ct) if ct else {"max_phase": 0.0, "n_mechanism_molecules": 0, "top_molecule": None}
        tr = sorted({f'{x["modality"]}:{x["label"]}' for x in tg["tractability"] if x["value"]})
        log(f"  chembl {tg['approvedSymbol']}")
        return (dict(symbol=tg["approvedSymbol"], name=tg["approvedName"], ensembl=tg["id"], uniprot=up,
                     chembl_target=ct, assoc=t["score"], diversity=diversity(sc),
                     n_datatypes_all=sum(1 for v in sc.values() if v > 0),
                     tractability="; ".join(tr), **lig, **ph), sc)

    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(one, trows))
    return pd.DataFrame([r for r, _ in res]), {r["symbol"]: sc for r, sc in res}


def pubmed_counts(pm, df: pd.DataFrame, disease_term: str, n: int) -> pd.DataFrame:
    df = df.copy()
    df["pubmed_total"] = pd.NA
    df["pubmed_disease"] = pd.NA
    for idx in df.nsmallest(n, "rank_main").index:
        s = df.at[idx, "symbol"]
        df.at[idx, "pubmed_total"] = pm.count(f'"{s}"[Title/Abstract]')
        df.at[idx, "pubmed_disease"] = pm.count(f'"{s}"[Title/Abstract] AND {disease_term}')
    return df


def run(out: Path, disease_name: str | None = None, ot=None, chembl=None, pm=None, log=print) -> dict:
    ot, chembl, pm = ot or OpenTargets(), chembl or ChEMBL(), pm or PubMed()
    out.mkdir(exist_ok=True, parents=True)
    log("selecting disease")
    best, sel = select_disease(ot)
    sel.to_csv(out / "disease_selection.csv", index=False)
    if disease_name:
        best = sel[sel.disease.str.lower() == disease_name.lower()].iloc[0].to_dict()
    log(f"disease: {best['disease']} ({best['efo']})")
    df, dts = collect(ot, chembl, best["efo"], POOL, log)
    df = add_scores(df)
    df = pubmed_counts(pm, df, f'"{best["disease"]}"[Title/Abstract]', PUBMED_N)
    g = build_graph(best["disease"], df, dts)
    df = df.merge(centrality_table(g), on="symbol")
    df.sort_values("rank_main").to_csv(out / "targets_scored.csv", index=False)
    meta = dict(date=dt.date.today().isoformat(), disease=best["disease"], efo=best["efo"],
                open_targets=ot.meta(), chembl=chembl.version(), pool=POOL, pubmed_n=PUBMED_N,
                weights=WEIGHTS, sensitivity=sensitivity(df),
                under_studied=under_studied(df)["symbol"].tolist(),
                graph=dict(nodes=g.number_of_nodes(), edges=g.number_of_edges()))
    (out / "run_meta.json").write_text(json.dumps(meta, indent=2, default=str))
    from .viz import make_all
    make_all(g, df, Path(out).parent / "figures", meta)
    return meta


if __name__ == "__main__":
    run(Path(sys.argv[1] if len(sys.argv) > 1 else "data"))
