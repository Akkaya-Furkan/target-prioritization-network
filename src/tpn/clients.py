"""Thin HTTP clients for Open Targets, ChEMBL and PubMed (NCBI E-utilities).

Every client takes an injectable `http` object (anything with .get/.post returning
objects with .json()), so tests run with a fake client and make no network calls.
"""
from __future__ import annotations

import time

import httpx

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def default_http() -> httpx.Client:
    return httpx.Client(timeout=90, follow_redirects=True)


def _retry(fn, tries=7):
    for i in range(tries):
        try:
            r = fn()
            r.raise_for_status()
            return r
        except (httpx.HTTPError, ValueError):
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


class OpenTargets:
    def __init__(self, http=None):
        self.http = http or default_http()

    def gql(self, query: str, variables: dict | None = None) -> dict:
        r = _retry(lambda: self.http.post(OT_URL, json={"query": query, "variables": variables or {}}))
        body = r.json()
        if body.get("errors"):
            raise RuntimeError(body["errors"])
        return body["data"]

    def meta(self) -> dict:
        return self.gql("{meta{apiVersion{x y z} dataVersion{year month}}}")["meta"]

    def find_disease(self, name: str) -> dict:
        d = self.gql('query($n:String!){search(queryString:$n,entityNames:["disease"],page:{index:0,size:1}){hits{id name}}}', {"n": name})
        return d["search"]["hits"][0]

    def disease_targets(self, efo: str, size: int) -> dict:
        q = """query($i:String!,$s:Int!){disease(efoId:$i){id name
          associatedTargets(page:{index:0,size:$s}){count rows{score
            target{id approvedSymbol approvedName proteinIds{id source}
              tractability{label modality value}}
            datatypeScores{id score}}}}}"""
        return self.gql(q, {"i": efo, "s": size})["disease"]


class ChEMBL:
    def __init__(self, http=None):
        self.http = http or default_http()

    def get(self, path: str, **params) -> dict:
        return _retry(lambda: self.http.get(f"{CHEMBL}/{path}.json", params=params)).json()

    def version(self) -> dict:
        return self.get("status")

    def single_protein_target(self, uniprot: str) -> str | None:
        d = self.get("target", target_components__accession=uniprot, limit=20)
        for t in d["targets"]:
            if t["target_type"] == "SINGLE PROTEIN" and t["organism"] == "Homo sapiens":
                return t["target_chembl_id"]
        return None

    def ligands(self, target: str, page_cap: int = 12) -> dict:
        """Distinct molecules with a pChEMBL-bearing activity, plus total such activities."""
        mols, total, offset = set(), 0, 0
        for _ in range(page_cap):
            d = self.get("activity", target_chembl_id=target, pchembl_value__isnull="false",
                         only="molecule_chembl_id", limit=1000, offset=offset)
            total = d["page_meta"]["total_count"]
            mols.update(a["molecule_chembl_id"] for a in d["activities"] if a.get("molecule_chembl_id"))
            if not d["page_meta"]["next"]:
                break
            offset += 1000
        return {"n_ligands": len(mols), "n_activities": total, "truncated": total > 1000 * page_cap}

    def max_phase(self, target: str) -> dict:
        """Highest max_phase among molecules with a curated mechanism on this target."""
        mech = self.get("mechanism", target_chembl_id=target, limit=1000)["mechanisms"]
        ids = sorted({m["molecule_chembl_id"] for m in mech if m.get("molecule_chembl_id")})
        if not ids:
            return {"max_phase": 0.0, "n_mechanism_molecules": 0, "top_molecule": None}
        best, best_id = 0.0, None
        for i in range(0, len(ids), 50):
            mols = self.get("molecule", molecule_chembl_id__in=",".join(ids[i:i + 50]), limit=100)["molecules"]
            for m in mols:
                p = float(m.get("max_phase") or 0)
                if p >= best:
                    best, best_id = p, m.get("pref_name") or m["molecule_chembl_id"]
        return {"max_phase": best, "n_mechanism_molecules": len(ids), "top_molecule": best_id}


class PubMed:
    def __init__(self, http=None, api_key=None):
        self.http = http or default_http()
        self.key = api_key

    def count(self, term: str) -> int:
        p = {"db": "pubmed", "term": term, "rettype": "count", "retmode": "json"}
        if self.key:
            p["api_key"] = self.key
        for attempt in range(5):
            r = self.http.get(f"{EUTILS}/esearch.fcgi", params=p)
            if r.status_code == 429:
                time.sleep(1.5 * (attempt + 1))
                continue
            r.raise_for_status()
            time.sleep(0.4)
            return int(r.json()["esearchresult"]["count"])
        raise RuntimeError("PubMed rate limit")
