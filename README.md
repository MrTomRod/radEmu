# Notes for the radEmu / fastEmu authors

I am Thomas Roder; I work 30% as a bioinformatics analyst at the University of Bern and the rest as CEO of my startup, Abrinca Genomics. While using and studying radEmu (2.4.0.0, `339e3a7`) and fastEmu (2.0.3, `4695b88`) closely across both roles and in my own free time, I noticed the points below. Thank you for the methods and the packages. Each file is short, self-contained, and has a minimal R reproduction runnable with the released packages where one was cheap (scripts in [`repro/`](repro/), run with [`repro/run.sh`](repro/run.sh)). Severity tags: **bug**, **minor bug**, **improvement**, **docs**, **question**.

| Note | Summary | Tag |
|---|---|---|
| [subset-constraint-discrete-design](subset-constraint-discrete-design.md) | On discrete (e.g. two-group) designs, a subset pseudo-Huber constraint silently tests the all-taxa median instead. | bug |
| [subset-null-offset](subset-null-offset.md) | `fit_null_symmetric_subset` drops an offset, so subset null fits (incl. fastEmu score tests) step downhill and never converge tightly. Fix on a fork branch. | bug |
| [wald-without-cis-crash](wald-without-cis-crash.md) | `emuFit(compute_cis = FALSE, return_wald_p = TRUE)` crashes because `X_cup` is `NULL`; always crashes in `fastEmuFit`. | bug |
| [fastemu-single-ref-bug](fastemu-single-ref-bug.md) | `fastEmuFit` rejects valid single-reference calls due to `length(ref_set > 1)` (plus two adjacent validation typos). | bug |
| [singular-information](singular-information.md) | Without the Firth penalty, a taxon absent from a group makes the information singular: the score test's ridge shifts the statistic, Wald stops. Patch on a fork branch. | bug |
| [null-fit-local-optima](null-fit-local-optima.md) | A multi-start null fit: score tests 7–10× faster at `J ≥ 500` (≈ 4.4 days → 15 h for all taxa at `J = 2000`) with equal or better null fits in every test. Patch on a fork branch. | improvement |
| [fastemu-reduced-refit](fastemu-reduced-refit.md) | `fastEmuFit` silently overrides a user's `tolerance` in the score-test refits (floor `1e-3`); plus a dead `Y_augmented` subsetting. Patch on a fork branch. | minor bug |
| [penalty-scale](penalty-scale.md) | The Firth penalty adds a fixed amount per cell: fine for read counts, but coverage-scale data inflate false positives (20–48%) and inflated values (×1e3–1e6) break taxa absent from a group. Suggest counts or a 1st-percentile rescaling. | question |
| [small-things](small-things.md) | Minor radEmu items: `tolerance` doc default, `pseudohuber_median` slack, `dplyr` in Suggests, fallback message wording, negative `Y`, penalty scale warning, `tol_discrete` floor, `max_abs_B` forwarding; plus fastEmu documentation (when to prefer exact `emuFit`, `constraint_diff`, reference-set default). | docs + bug |
