# Subset constraint on a two-group design tests the all-taxa median instead

**What happens.** A user asks to compare each taxon against the smoothed median of a *subset* of taxa (e.g. `function(x) pseudohuber_median(x[1:12])`). On a two-group (discrete) design, `emuFit`'s score test silently compares against the median of *all* taxa instead. The estimates and Wald intervals use the requested subset, so the estimate and the p-value next to it answer different questions.

**Why.** For discrete designs, `emuFit` uses a fast solver that only knows the all-taxa median. It decides whether a constraint is "the all-taxa median" by searching the function's source code for the text `pseudohuber_median`, which a subset constraint also contains. The non-discrete branch just below checks the `reference_set` attribute; the discrete branch does not.

**Example** ([`repro/subset_discrete.R`](repro/subset_discrete.R), 30 samples, 25 taxa, two groups, subset = taxa 1–12):

| null fit | gap to all-taxa median | gap to requested subset median | score statistic | p |
|---|--:|--:|--:|--:|
| default (discrete solver) | −5.6e-10 | **−0.023** | 0.0329 | 0.856 |
| `null_fit_alg = "augmented_lagrangian"` | +0.022 | −4.1e-6 | 0.0021 | 0.963 |

The default run satisfies the wrong constraint exactly and misses the requested one.

**Fix** ([fork branch](https://github.com/MrTomRod/radEmu/tree/fix-subset-constraint-discrete), [diff](https://github.com/MrTomRod/radEmu/compare/339e3a774d8047a78b0e722137b8c4164ffae728...fix-subset-constraint-discrete)): the discrete solver is used only if the constraint carries no `reference_set` attribute and actually *behaves* like the all-taxa median on random input, the way the mean is already detected. Everything else falls back to the augmented Lagrangian, as other constraints already do. Two tests added; existing two-group results are bit-identical and the package's own suite passes.

**Links**
- radEmu at the pinned commit: discrete branch [`emuFit.R#L462-L486`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/emuFit.R#L462-L486), non-discrete branch [`emuFit.R#L488-L545`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/emuFit.R#L488-L545)
