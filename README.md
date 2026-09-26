# Target prioritisation network (Open Targets + ChEMBL)

**Question.** In a chosen disease area, which targets stand out on evidence and tractability?

## Data
Open Targets Platform (GraphQL), ChEMBL, PubMed publication counts.

## Approach
1. Disease-area selection by data density
2. Associated targets, evidence-type breakdown and association scores
3. Tractability from ChEMBL: known ligands, highest development phase
4. Target–evidence–disease graph with centrality measures (NetworkX)
5. Transparent weighted score; top targets vs publication volume (under-studied but well-supported)
6. Interactive network (pyvis)

## Deliverables
Interactive HTML network, summary table, 2 figures, write-up.

## Limitations
- A database/literature analysis — no clinical or treatment recommendation.
- Association scores carry publication bias; weights change the ranking and are stated explicitly.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
