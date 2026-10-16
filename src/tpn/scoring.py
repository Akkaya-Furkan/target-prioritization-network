"""Transparent target score.

score = assoc^a * diversity^b * ligand^c   (weighted geometric mean, weights sum to 1)
  assoc     Open Targets overall association score (0-1)
  diversity (# datatypes with score>0, excluding known_drug) / (# datatypes, excl. known_drug)
  ligand    log10(1+n_ligands) / log10(1+max n_ligands among candidates)   (0 if no ligands)
known_drug is excluded from diversity to avoid counting ChEMBL-derived drug evidence twice.
"""
from __future__ import annotations

import math

import pandas as pd

DATATYPES = ["genetic_association", "somatic_mutation", "affected_pathway",
             "rna_expression", "literature", "animal_model"]  # known_drug excluded
WEIGHTS = {"main": (1 / 3, 1 / 3, 1 / 3), "assoc_heavy": (0.6, 0.2, 0.2), "ligand_heavy": (0.2, 0.2, 0.6)}


def diversity(datatype_scores: dict[str, float]) -> float:
    return sum(1 for d in DATATYPES if datatype_scores.get(d, 0) > 0) / len(DATATYPES)


def add_scores(df: pd.DataFrame, weights: dict = WEIGHTS) -> pd.DataFrame:
    df = df.copy()
    top = max(df["n_ligands"].max(), 1)
    df["ligand_norm"] = df["n_ligands"].apply(lambda n: math.log10(1 + n) / math.log10(1 + top))
    for name, (a, b, c) in weights.items():
        # small epsilon keeps a zero component from erasing the whole score (reported in README)
        e = 0.02
        df[f"score_{name}"] = ((df["assoc"] + e) ** a) * ((df["diversity"] + e) ** b) * ((df["ligand_norm"] + e) ** c)
        df[f"rank_{name}"] = df[f"score_{name}"].rank(ascending=False, method="min").astype(int)
    return df


def sensitivity(df: pd.DataFrame, k: int = 10) -> dict:
    """Top-k overlap and Spearman rank correlation between the main and alternative weightings."""
    out = {}
    base = set(df.nsmallest(k, "rank_main")["symbol"])
    for name in WEIGHTS:
        if name == "main":
            continue
        alt = set(df.nsmallest(k, f"rank_{name}")["symbol"])
        out[name] = {"top_k_overlap": len(base & alt), "k": k,
                     "spearman": float(df["rank_main"].corr(df[f"rank_{name}"], method="pearson")),  # Pearson on ranks = Spearman
                     "entered": sorted(alt - base), "left": sorted(base - alt)}
    return out


def under_studied(df: pd.DataFrame, pool: int = 30) -> pd.DataFrame:
    """Top-15 targets whose disease-specific PubMed volume is below the median of the top-`pool` targets."""
    ref = df.nsmallest(pool, "rank_main")["pubmed_disease"].median()
    top = df.nsmallest(15, "rank_main")
    return top[top["pubmed_disease"] < ref]
