# The Firth penalty is fine for read counts but breaks on coverage-scale or inflated data

**Counts (the usual case).** The Firth penalty adds half a hat value (a number between 0 and 0.5) to every cell of `Y`, mostly where the data are thin. Against read counts that is the intended amount: it keeps estimates finite when a taxon is absent from a group and otherwise stays out of the way.

**The failure modes.** The amount is fixed in `Y`'s own units, while everything else in radEmu is scale-free. So the result depends on the arbitrary units of `Y`:
* **Too small** (coverage, e.g. the trimmed-mean ribosomal-protein coverage of the fastEmu marine analysis; relative abundances): the penalty swamps the data.
* **Too large** (bases instead of reads, absolute abundances such as cells per gram): the penalty vanishes exactly where it is needed.

**How big** (released radEmu, default settings, simulated two-group data with known truth, 200–250 replicates per setting; the coverage row is a Poisson setting with n = 30, the inflated rows a sparse ZINB setting with n = 10, where absent taxa are common):

| Same counts, expressed as | false positives, null taxon (5% nominal) | estimates of taxa absent from one group |
|---|---|---|
| counts | 5% | median \|estimate\| ≈ 3 |
| coverage, median ≈ 0.09 | 20% (48% for rare taxa); estimates shrunk ~22% | ≈ 0.1–0.3 |
| counts × 1e3 | 2–3% | ≈ 10 (+log 1000); Wald calls 100% of such null taxa significant |
| counts × 1e6 | 2–3%; power on these taxa drops from 0.14 to 0.02 | ≈ 17 (+log 1e6) |

**Suggestion.** A sentence in `?emuFit`: supply read counts, or convert coverage back to counts (`coverage × genome or marker length / read length`). When that factor is unknown, rescaling `Y` so that the 1st percentile of its nonzero values is 1 recovered count-scale behaviour in both directions in our simulations (rescaling by the smallest nonzero value also works, but a single tiny outlier cell hijacks it). A longer-term idea for you to consider: scale the augmentation by an estimated dispersion, so that the penalty no longer depends on the units of `Y` (but would also change results for overdispersed counts).

**Links**
- radEmu: [`get_augmentations.R#L25-L33`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/get_augmentations.R#L25-L33), [`get_augmentations_discrete.R`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/get_augmentations_discrete.R), [`emuFit.R#L3`](https://github.com/statdivlab/radEmu/blob/339e3a774d8047a78b0e722137b8c4164ffae728/R/emuFit.R#L3) (accepts any nonnegative `Y`)
- Papers: radEmu §3.1 (Firth augmentation), §7 (coverage data); fastEmu preprint p. 2 (`Y` = trimmed-mean coverage of ribosomal protein L14)
