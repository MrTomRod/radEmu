# `fit_null_symmetric_subset` drops an offset, so subset null fits never converge tightly

**What happens.** When the tested taxon is outside the reference subset, `fit_null_symmetric_subset()` updates its other coefficients with `micro_fisher()` but passes `z` without the fixed term `X[, k_constr] * B[k_constr, j]`. The step therefore goes the wrong way. Because the backtracking accepts losses up to `1e-2`, a small downhill step is accepted on every iteration, so the likelihood and statistic oscillate and tight tolerances never converge (the fit then falls back to the augmented Lagrangian).

**Why it matters.** fastEmu's default score tests on non-discrete designs use this function, and almost every tested taxon lies outside the reference set; `emuFit` with any subset constraint is affected too.

**Example** ([`repro/subset_null_offset.R`](repro/subset_null_offset.R), simulated, `n = 40`, `J = 150`, continuous covariate, tight tolerances):

| | released radEmu | patched |
|---|---|---|
| `fit_null_symmetric_subset()` | not converged after 120 iterations, 41 downhill steps, statistic oscillates 5.056–5.074 | converged in 29 iterations, never downhill, 5.049 |
| `fastEmuFit()`, one taxon | 3.22 (augmented Lagrangian: 3.35) | 3.33 |

**Fix.** Pass the offset and keep matrix dimensions (`drop = FALSE`): three changed lines plus a test; the package's own suite passes.

**Links**
- radEmu: [`fit_null_symmetric_subset.R#L546-L594`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/fit_null_symmetric_subset.R#L546-L594)
- Fork branch: https://github.com/MrTomRod/radEmu/tree/fix-subset-null-offset · diff: [compare](https://github.com/MrTomRod/radEmu/compare/339e3a774d8047a78b0e722137b8c4164ffae728...fix-subset-null-offset)
