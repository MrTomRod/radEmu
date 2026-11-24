# Score tests up to 10× faster, more reliable null fit

**Headline.** A drop-in null-fit algorithm, `null_fit_alg = "multistart"`, made radEmu's score tests **7–10× faster** per test at `J ≥ 500` in our R benchmark. In every test it found an equal or better null fit than the current default: never worse, and better in 33 of 45. Projected to testing all taxa at `J = 2000`, that is **about 15 hours instead of about 4.4 days** on one core. The null fit itself went from minutes to under a second; what remains is almost all `get_score_stat()`.

| Design, `n = 40` | `J` | default (s/test) | multistart (s/test) | speed-up | likelihood better / equal / worse |
|---|--:|--:|--:|--:|:-:|
| two-group | 500 | 8.7 | 0.9 | 10× | 3 / 2 / 0 |
| two-group | 1000 | 34.4 | 3.9 | 9× | 2 / 3 / 0 |
| two-group | 2000 | 193 | 26 | 7.5× | 2 / 3 / 0 |
| continuous | 1000 | 29.6 | 4.2 | 7× | 5 / 0 / 0 |
| continuous | 2000 | 187 | 28 | 7× | 5 / 0 / 0 |

R 4.5.2, one thread, simulated data, median of 5 tested taxa per setting (same taxa for both methods). The default dispatches the constraint-sandwich algorithm here.

**Why the fits get better.** Under the null `B[k, j] = g(B[k, ])` with the pseudo-Huber median, the feasible set is not convex: there is one "valley" for every way the other taxa can split above and below `B[k, j]`. A solver started in one place settles in the nearest valley, so the null fit, and with it the p-value, depends on starting values and tolerances. On `wirbel_otu_small` (taxon 24), tightening the tolerances lands 45 log-likelihood units *lower* and moves `p` from 0.048 to 0.055. Across all 47 wirbel taxa, multistart finds a better null fit for 17–19 of them (up to +120 log-likelihood units), changing 2–3 decisions at 0.05.

**What it does.** Same output structure as the existing null fits, so `score_test()` and `get_score_stat()` are unchanged; opt-in, defaults untouched.
* Starts from several deliberately different points: the full-model estimate, all zeros, and a flat tested row.
* Starts from a very smooth version of the median (almost the mean, which has a single optimum) and sharpens it step by step while following the solution.
* Tries a few arrangements of the other taxa above and below the tested one, then small swaps across the median.
* Polishes every candidate with an exact Newton solver on the reduced problem and keeps the best fit that satisfies the constraint (residual ≤ 5e-16).

**Where it is not faster.** On small two-group designs (`J < 150`) R's discrete algorithm with its loose default tolerance is about 2× faster (wirbel: 0.18 vs 0.37 s/test). At tight tolerances it is about 90× slower.

**Next step, if you are interested.** With the null fit this cheap, `get_score_stat()` is now 97–99% of the time per test at `J = 2000`, because it forms the full `pJ × pJ` information matrix. The information is block-diagonal minus a low-rank term, so the statistic can be computed without forming it, which should give another large speed-up. We would be happy to help with that.

**Links**
* Fork branch: https://github.com/MrTomRod/radEmu/tree/multistart-null-fit · diff: [compare](https://github.com/MrTomRod/radEmu/compare/339e3a774d8047a78b0e722137b8c4164ffae728...multistart-null-fit) (4 commits, 112 new testthat assertions; the full radEmu suite passes, 366 / 0)
* Minimal R repro of the local optimum (wirbel taxon 24, default vs tight): [`repro/wirbel_local_optimum.R`](repro/wirbel_local_optimum.R)
