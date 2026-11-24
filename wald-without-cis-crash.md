# `emuFit(compute_cis = FALSE, return_wald_p = TRUE)` crashes (and always crashes in `fastEmuFit`)

Asking for Wald p-values without confidence intervals fails with `requires numeric/complex matrix/vector arguments` whenever score tests are turned off. An internal design matrix (`X_cup`) is only built when `compute_cis` or `run_score_tests` is true, but the Wald code runs whenever `compute_cis` or `return_wald_p` is true and tries to index a `NULL` matrix. Because `fastEmuFit` first calls `emuFit(run_score_tests = FALSE)`, passing `compute_cis = FALSE, return_wald_p = TRUE` to `fastEmuFit` always crashes even when score tests are on. Including `return_wald_p` in the `if (compute_cis | run_score_tests)` guard at `emuFit.R#L248` fixes both packages.

**Links**
- Repro: [`repro/nocis_crash.R`](repro/nocis_crash.R)
- radEmu: [`emuFit.R#L247-L250`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/emuFit.R#L247-L250) (`X_cup` guard), [`emuFit.R#L368-L380`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/emuFit.R#L368-L380) (Wald call), [`micro_wald.R#L55-L58`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/micro_wald.R#L55-L58)
- fastEmu: [`fastEmuFit.R#L247-L262`](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L247-L262)
- Fork branch: https://github.com/MrTomRod/radEmu/tree/fix-wald-without-cis · diff: [compare](https://github.com/MrTomRod/radEmu/compare/339e3a774d8047a78b0e722137b8c4164ffae728...fix-wald-without-cis) (2 commits, adds 1 test in `tests/testthat/test-emuFit.R`; the package's own suite passes)
