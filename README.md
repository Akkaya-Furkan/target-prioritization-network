# Target prioritisation network (Open Targets + ChEMBL)

**Question.** In a chosen disease area, which targets stand out on evidence and tractability?
**Answer for ulcerative colitis (UC):** under an equal-weight score, JAK2 and JAK3 rank first; IL23R, TYK2, PPARG, PTGS2, NR3C1 and S1PR1 follow, and the *order below the top two is weight-dependent* (see sensitivity). This is a database/literature exercise, **not** a clinical or treatment recommendation.

All numbers come from live API calls on **2026-10-16**: Open Targets Platform **26.09** (API 26.9.0), **ChEMBL_37** (released 2026-05-01), PubMed via NCBI E-utilities. Raw outputs: [`data/`](data), run metadata: [`data/run_meta.json`](data/run_meta.json).

## Deliverables
- Interactive network: [`network.html`](network.html) (pyvis; top-25 targets, evidence types, disease; hover for scores)
- Figures: [`figures/network_top.png`](figures/network_top.png), [`figures/score_vs_publications.png`](figures/score_vs_publications.png)
- Tables: [`data/targets_scored.csv`](data/targets_scored.csv) (200 targets), [`data/summary_top10.md`](data/summary_top10.md), [`data/disease_selection.csv`](data/disease_selection.csv)

## Method
1. **Disease choice (data density).** Six candidates (UC, rheumatoid arthritis, psoriasis, Parkinson disease, ALS, asthma) were compared on their top-200 associated targets; the rule was *most targets with >=3 evidence types* (known_drug excluded). UC won (67 targets, vs 61 RA, 56 psoriasis, 51 Parkinson, 31 ALS, 22 asthma). Open Targets lists 9,046 UC-associated targets; the top 200 by overall association score form the candidate pool.
2. **Open Targets (GraphQL):** overall association score, per-datatype scores (evidence-type breakdown), tractability labels.
3. **ChEMBL REST:** human single-protein target via UniProt; *ligands* = distinct molecules with a pChEMBL-bearing activity (capped at 12,000 activity records; 4 targets hit the cap); *max phase* = highest `max_phase` among molecules with a curated mechanism on the target.
4. **Graph (NetworkX):** tripartite target - evidence type - disease (207 nodes, 437 edges); degree, betweenness and closeness per target are in the CSV. Centrality in this graph is largely determined by how many evidence types a target has, so it is descriptive, not a second ranking.
5. **Score (explicit, weighted geometric mean):**
   `score = (assoc + e)^a * (diversity + e)^b * (ligand + e)^c`, e = 0.02 so a zero component does not erase the product.
   - `assoc` = Open Targets overall score; `diversity` = evidence types with score>0 among 6 (known_drug excluded to avoid double-counting drug evidence); `ligand` = log10(1+n)/log10(1+max n).
   - Weights (a, b, c): **main 1/3, 1/3, 1/3**; alternatives **assoc-heavy 0.6/0.2/0.2**, **ligand-heavy 0.2/0.2/0.6**.
6. **PubMed volume:** for the top-30 targets, counts of `"SYMBOL"[Title/Abstract]` and of the same AND `"ulcerative colitis"[Title/Abstract]`.

## Top 10 (main score)
*Max phase = highest phase of any molecule with a mechanism on the target, for any indication (not UC-specific).*

| # | Target | Assoc. | Evidence types (of 6) | ChEMBL ligands | Max phase* | Score | Rank (assoc-heavy / ligand-heavy) | PubMed gene | PubMed gene+UC |
|---|---|---|---|---|---|---|---|---|---|
| 1 | JAK2 | 0.73 | 4 | 8,942 | 4 | 0.806 | 1 / 1 | 16,907 | 153 |
| 2 | JAK3 | 0.65 | 3 | 7,463 | 4 | 0.702 | 2 / 2 | 2,434 | 40 |
| 3 | ESR1 | 0.36 | 4 | 4,581 | 4 | 0.627 | 19 / 3 | 6,528 | 27 |
| 4 | CXCR2 | 0.53 | 3 | 1,092 | 3 | 0.609 | 12 / 12 | 3,837 | 35 |
| 5 | IL23R | 0.60 | 3 | 489 | 0 | 0.607 | 6 / 21 | 1,380 | 168 |
| 6 | TYK2 | 0.62 | 2 | 6,175 | 4 | 0.605 | 4 / 4 | 1,848 | 42 |
| 7 | PPARG | 0.63 | 2 | 4,357 | 4 | 0.600 | 3 / 6 | 3,792 | 36 |
| 8 | PTGS2 | 0.63 | 2 | 4,522 | 4 | 0.599 | 5 / 5 | 4,583 | 84 |
| 9 | NR3C1 | 0.62 | 2 | 3,133 | 4 | 0.588 | 7 / 7 | 1,734 | 7 |
| 10 | S1PR1 | 0.59 | 2 | 3,424 | 4 | 0.582 | 11 / 8 | 1,169 | 26 |

## Sensitivity to weights
- Assoc-heavy: top-10 overlap with main = **7/10** (entered: ALOX5, PTGS1, TNF; left: CXCR2, ESR1, S1PR1); Spearman over all 200 targets = 0.978.
- Ligand-heavy: overlap **8/10** (entered: ALOX5, RORC; left: CXCR2, IL23R); Spearman = 0.986.
- So the overall ordering is stable but top-10 *membership* is not. ESR1 is #3 only because of ligand count and evidence breadth (association score 0.36); under assoc-heavy weights it falls to #19. IL23R drops from #5 to #21 under ligand-heavy weights (489 ligands, max phase 0 by the mechanism criterion).

## Under-studied but well-supported
Among the top 15, targets with fewer gene+UC PubMed records than the median of the top-30 pool: **NR3C1 (7), ALOX5 (7), PTGS1 (5), FGFR3 (1), PRKCB (3)**. These are candidates for a closer look, not findings. Counts are text-mining proxies: TNF has 5,429 gene+UC records, JAK2 153. Aliases (e.g. "glucocorticoid receptor" for NR3C1) are not searched, so counts for such genes are undercounted.

## Limitations
- Database/literature analysis only; **no clinical or treatment recommendation**.
- **Publication bias:** Open Targets literature and association scores favour well-studied genes; PubMed counts measure attention, not truth.
- **Weights decide the ranking** (above); the equal weights are a choice, not a calibration.
- Open Targets *known_drug* evidence is excluded from diversity but still contributes to the overall association score.
- ChEMBL ligand count mixes potent and weak compounds and is capped for 4 targets; 77 of 200 pool targets have no human single-protein ChEMBL target and 107 have zero ligands, so they score low on that term by construction. Max phase ignores indication.
- Gene-symbol PubMed queries can match unrelated uses of a symbol.
- Only the top 200 of 9,046 associated targets were scored; only the top 30 were queried in PubMed.

## Reproduce
```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
PYTHONPATH=src .venv/bin/python -m pytest -q tests          # offline, fake clients
PYTHONPATH=src .venv/bin/python -m tpn.pipeline data         # live run, ~8 min
```
Tests use recorded-style fake responses ([`tests/fixtures/responses.json`](tests/fixtures/responses.json)) and make no network calls.
