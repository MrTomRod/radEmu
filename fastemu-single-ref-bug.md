# `fastEmuFit`: a single reference taxon does not work

Using one reference taxon (e.g. a spike-in: `constraint_fn = 3, reference_set = 3`) stops with an error saying the two arguments do not match ([`repro/fast_single_ref.R`](repro/fast_single_ref.R)). A few small slips in the argument checks cause this:

- **[L224](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L224):** `length(ref_set > 1)` should be `length(ref_set) > 1`. As written it is `1` (true) for a single taxon, so every valid call errors.
- **[L231](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L231):** checks `"data-driven"` (hyphen) instead of `"data_driven"`, so it is always false.
- **[L234](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L234):** compares `ref_set` with the whole list `constraint_fn` instead of `constraint_fn[[k]]`.
- **Indexing:** `emuFit_check` has already turned `constraint_fn = 3` into `function(x) x[3]` before this loop. So the `is.numeric(constraint_fn[[k]])` checks ([L202](https://github.com/statdivlab/fastEmu/blob/4695b88b0a621cff2e2eb94ed7ca0fab223203af/R/fastEmuFit.R#L202), L231, L233) are never true, and later the function is applied to the reference-set subset (`x[reference_set]`, length 1), where `x[3]` is `NA`.

**Fix** ([fork branch](https://github.com/MrTomRod/fastEmu/tree/fix-single-reference), [diff](https://github.com/MrTomRod/fastEmu/compare/4695b88b0a621cff2e2eb94ed7ca0fab223203af...fix-single-reference)): fixes the three typos, checks the user's original `constraint_fn`, and uses `x[1]` on the one-taxon subset. About ten changed lines plus a test; the package's own suite passes.
