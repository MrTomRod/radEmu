# Without the Firth penalty, a taxon missing from a group breaks the Wald test and nudges score tests

**When it happens.** With `penalize = FALSE`, a taxon that is never observed in one group (or only in a single sample) has no finite estimate: its coefficient runs off to the `max_abs_B` limit (about −50 here). The data say nothing about that coefficient, so the information matrix has a direction with (numerically) zero information, i.e. it is singular. With the Firth penalty on, this does not happen.

**What radEmu does now.**
- **Wald:** `micro_wald` cannot invert the matrix and `emuFit(compute_cis = TRUE)` stops with *"computationally singular"*. No standard errors for any taxon, not just the missing one.
- **Score test:** `get_score_stat` adds a small ridge, `lambda · max(diag(I))` with `lambda = 1e-4, 1e-3, …`, to the whole diagonal until the matrix inverts. That ridge is arbitrary and changes the statistic of every test it touches, by about 2% in the example below.

**Example** ([`repro/singular_information.R`](repro/singular_information.R): 12 samples, 160 taxa, two groups, taxon t3 absent from group 2):

| | released radEmu | with the patch |
|---|---|---|
| Wald | error: computationally singular | runs; se for t3 is `NA`, others normal (t5: 0.315) |
| score test of t1 | 0.002605, ridge `lambda = 1e-4` used | 0.002560, no ridge |

**What the patch does** ([fork branch](https://github.com/MrTomRod/radEmu/tree/singular-information), [diff](https://github.com/MrTomRod/radEmu/compare/339e3a774d8047a78b0e722137b8c4164ffae728...singular-information)): instead of a ridge, it finds the directions the data do not inform and inverts the matrix only on the rest (a pseudo-inverse). Quantities that depend on the uninformed directions get `NA` (here: t3's standard error); everything else is computed normally. Problems without singularity take the old code path and give identical results. `control_fn(information_treatment = "r_ladder")` restores the old behaviour. Five commits, one new file `R/information_pinverse.R`, 13 tests; the package's own suite passes.

**Links**
- radEmu at the pinned commit: ridge in [`get_score_stat.R#L76-L102`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/get_score_stat.R#L76-L102), inversion in [`micro_wald.R#L85-L86`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/micro_wald.R#L85-L86)
