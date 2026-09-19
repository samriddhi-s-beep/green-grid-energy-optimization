# Green Grid Energy Optimization

A Mixed-Integer Linear Programming (MILP) model that determines the optimal electricity procurement plan for a renewable energy provider, minimizing total cost while satisfying supplier, alliance, and contractual constraints.

## Problem Statement

GreenGrid Energy procures electricity from multiple renewable suppliers across several delivery periods. The model decides **which suppliers to select** and **how much to procure from each**, subject to:

- **Demand satisfaction** — total procurement must meet demand each period
- **Supplier concentration limits** — no single supplier can exceed a set share of demand per period
- **Alliance concentration limits** — supplier groups (alliances) are capped to avoid over-reliance on related suppliers
- **Minimum commitment constraints** — once a supplier is selected, a minimum volume must be procured from them

## Approach

- Formulated as a MILP with continuous procurement variables and binary supplier-selection variables
- Implemented in **R** using the `ompr` modeling framework
- Solved and benchmarked using two solvers: **HiGHS** and **GLPK**
- Tested across three problem scales to evaluate solver performance as complexity grows

## Results

| Dataset | Periods | Suppliers | HiGHS Objective | HiGHS Runtime | GLPK Objective | GLPK Runtime | LP Relaxation |
|---|---|---|---|---|---|---|---|
| Small | 3 | 8 | 26,418.61 | 0.05s | 26,418.61 | 0.12s | 26,418.61 |
| Medium | 300 | 30 | 8,077,539 | ~1,214s (20 min) | 8,101,966 | ~3,627s (60 min) | 7,682,880 |
| Large | 400 | 50 | 10,611,363 | 3,600s (time limit reached) | 10,679,796 | 3,632s (time limit reached) | 10,564,340 |

## Key Findings

- **Small scale:** the MILP solution matches its LP relaxation exactly — integrality constraints don't bind, and both solvers solve in under a second.
- **Medium scale:** HiGHS finds a better solution than GLPK (8.08M vs 8.10M) in roughly a third of the time (20 min vs 60 min).
- **Large scale:** neither solver proves optimality within a 1-hour limit, but HiGHS still returns a better feasible solution than GLPK (10.61M vs 10.68M) in the same time budget.
- **Recommendation:** HiGHS is the stronger solver for this problem class, delivering equal-or-better solution quality with consistently lower runtime as problem size increases.

## Files

- `GGDS_Model.R` — model formulation, solver setup, and benchmarking script
- `GGDS_small.csv`, `GGDS_medium.csv`, `GGDS_large.csv` — input datasets at three problem scales

## How to Run

1. Install R and the required packages: `ompr`, `ompr.roi`, `ROI.plugin.glpk`, `ROI.plugin.highs`, `magrittr`
2. Place the CSV files in the same directory as `GGDS_Model.R`
3. Run the script — it solves all three datasets and prints the results table
