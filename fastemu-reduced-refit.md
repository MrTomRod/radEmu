# `fastEmuFit` ignores a user's `tolerance` in the score-test refits

`fastEmuFit` refits each reduced model with `tolerance = max(tolerance, 0.001)`, so a user who asks for `tolerance = 1e-9` silently gets `1e-3` (looser even than the package default `1e-4`; radEmu's docs say `1e-3`, its code `1e-4`). The floor arrived with commit `6981fbb` (April 2025), which switched these refits to `refit = TRUE`; the `Y_augmented` subsetting just before the call is a leftover from `refit = FALSE` and is now unused.
The floor is a sensible speed default: on a continuous design (`n = 40, J = 80, p = 3`) `1e-4` doubled the time per test (1.11 s vs 0.56 s) for a statistic change of only ~1e-4 relative; on two-group designs it makes no difference. But results can move more on other data (up to 2% on a design with 20 reference taxa), so users should be able to tighten it.
The patch keeps `1e-3` as the refit default, uses an explicitly supplied `tolerance` as given, drops the unused subsetting, and documents this in `?fastEmuFit`.

**Links**
- fastEmu: [`fastEmuFit.R#L98-L102`](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L98-L102), [`fastEmuFit.R#L549-L568`](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L549-L568), commit [`6981fbb`](https://github.com/statdivlab/fastEmu/commit/6981fbbbc0a68eacade4dbd404be25a8b33c4ffb)
- radEmu: [`emuFit.R#L258-L282`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/emuFit.R#L258-L282) (`refit = TRUE` path)
- Fork branch: https://github.com/MrTomRod/fastEmu/tree/respect-user-tolerance · diff: [compare](https://github.com/MrTomRod/fastEmu/compare/4695b88b0a621cff2e2eb94ed7ca0fab223203af...respect-user-tolerance) (1 commit; the package's own suite passes)
