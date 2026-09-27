# nlmixr2scm 0.4.1

* `runSCM()` results now have class `"nlmixr2scm"`, with new `summary()` and
  `print()` methods. `summary()` reports the options used for the search, a
  comparison of the base, final forward and final backward models (OFV, dOFV,
  AIC, BIC, parameter count and SCM covariates), the step and all-candidate
  tables, and the location of saved models and report files on disk.
  `print()` gives a short overview. The result gains `baseFit`, `options`
  and `outputDir` elements; the existing `summaryTable`, `resFwd` and
  `resBck` elements are unchanged.

* `runSCM()` gained a `retryOnUnderflow` argument (default `TRUE`, appended
  as the last parameter to preserve positional calls). The retry mechanism
  treats a forward candidate whose p-value underflows to zero as an
  unrealistic OFV; for `df = 1` that is any OFV drop above about 70, so
  genuinely strong covariate effects were always retried. Set
  `retryOnUnderflow = FALSE` to skip that criterion and the extra fits; the
  OFV-increase and `maxDeltaOFV` criteria are unaffected.

* Fixed a bug in the forward search mixing two `dOFV` sign conventions. Since 0.3,
  `dOFV` is `candidate OFV - reference OFV` (negative = improvement), but the
  profile-on-stall rescue, retry best-attempt tracking and p-value tie-break
  still assumed the pre-0.3 "positive = improvement" sign. As a result:
  * every improving forward candidate was treated as stalled and refitted
    (roughly tripling forward-search run time), and the rescue refit was
    accepted even when it reached a worse OFV than the original fit;
  * rescued candidates were reported with a flipped `deltObjf` sign, so
    `summaryTable` and the printed "Ref OFV" column were inconsistent;
  * when all retries were exhausted, the attempt with the highest OFV was
    kept instead of the lowest;
  * p-value underflow ties were broken in favour of the smallest OFV drop
    instead of the largest.

# nlmixr2scm 0.4

* Fixed: `runSCM()` wrote `scm_log.txt`, `scm_step_summary.csv` and
  `scm_all_candidates.csv` into the current working directory even when
  called with `saveModels = FALSE`. Those report files are now written only
  when `saveModels = TRUE`, so `saveModels = FALSE` leaves the working
  directory untouched.
* `runSCM()` gained an `rxThreads` argument (appended as the last
  parameter, to preserve positional-call compatibility) giving explicit,
  independent control over the number of rxode2 OpenMP threads used per
  parallel `workers` process. Whenever more than one worker is involved,
  `workers * rxThreads` is checked against the machine's physical core
  count before any candidate fitting starts, and the call aborts with an
  explanatory error if it would be exceeded -- each parallel worker runs
  its own independent rxode2 thread pool, so an unconstrained combination
  could otherwise request far more OS threads than the machine has. A
  single worker is never subject to this check. Requires `nlmixr2utils (>=
  0.3)`.
* **Breaking change:** an existing `runSCM(workers = N, ...)` call with `N
  > 1` (or `workers = "auto"`) that does not also set `rxThreads` may now
  abort on most multi-core machines, where it previously ran silently
  oversubscribed. Set `rxThreads` explicitly (e.g. `rxThreads = 1`) to
  restore the prior behavior.

# nlmixr2scm 0.3

* Changed `dOFV` sign convention: `dOFV` is now reported as `candidate OFV − reference OFV` uniformly for both forward and backward steps. A negative value indicates an OFV improvement (forward addition), a positive value indicates an OFV increase (backward removal). Previously the sign was flipped so that `dOFV` was always positive for "meaningful" changes, which was internally convenient but inconsistent with the standard pharmacometric interpretation. The `maxDeltaOFV` retry criterion now compares against `|dOFV|`.

* Fixed the retry mechanism for unrealistic OFV values (added in 0.2) to accept the best attempt seen across all retries when every attempt is exhausted, rather than whichever attempt happened to run last. Under perturbed-init retries, the last attempt was the best only by coincidence.

* `runSCM()` gained a `centers` argument to fix the centering (reference) value for one or more continuous covariates instead of using the per-dataset median, e.g. `centers = c(BW = 70, CrCL = 95)`. This keeps estimated covariate coefficients on the same reference across datasets and matches data-generating models that used fixed references; covariates not named in `centers` continue to use the median.

* `runSCM()` gained `profileInit` and `profileInitOnStall` (default `TRUE`) arguments to warm-start a forward candidate's covariate coefficient with a cheap frozen 1-D profile (Brent method) before the real estimator runs, with structural thetas and between-subject variability fixed at their parent values. This gives gradient optimisers a nonzero, gradient-informative starting value so ODE covariate candidates no longer stall at the flat zero-effect initial estimate under solver-noise-flattened objectives. `profileInitOnStall` fires the rescue only when a candidate's ordinary fit stalls (OFV improvement over its parent `<= stallTol`, default `0`); `profileInit` forces the profile for every forward candidate.

* Fixed forward and backward winner selection in `runSCM()` to break ties on `deltObjf` instead of row order. For `df = 1`, `dOFV >= ~70.5` underflows `1 - pchisq()` to exactly `0` (forward) or `1` (backward), so several strong candidates could tie on `pchisqr`; `which.min()`/`which.max()` then silently picked whichever candidate happened to sort first rather than the one with the largest (forward) or smallest (backward) OFV change.

* Fixed `runSCM()` erroring on every candidate ("numbers of columns of arguments do not match") against current `rxode2`/`nlmixr2est`, which added a `prior` column to `iniDf`. The internal `.rebuildUiFromPairs()` and `.builduiCovariate()` helpers hardcoded the `iniDf` column set when adding a covariate theta row; they now pad the new row to match `iniDf`'s actual columns instead.

# nlmixr2scm 0.2

* Yaping Liu added as package co-author.
* Fixed `scmAddCatCovariates()` selecting the categorical reference level as the first-appearing level (or, for factors, the first `levels()` entry) rather than the most frequent level. The documented rule, and the one `.makeSCMData()` applies on the `runSCM()` search path, is most-frequent; the helper's `uniqueVals[-1]` is not a frequency rule at all.

* `scmAddCatCovariates()` gained a `freqBy` argument controlling how level frequencies are counted. The default `"id"` counts one vote per subject, matching `.makeSCMData()`, so a few densely sampled subjects cannot outvote the majority. `"observation"` counts one vote per row, which is the only sensible rule for a column that is not unique to a subject (e.g. `CMT`), since per-subject counting would reduce such a column to whichever level landed in each subject's first record. `"auto"` chooses per column by testing whether any subject carries more than one non-missing level.

* `scmAddCatCovariates()` also now takes a `catCutoff` argument (default `0`) to lump rare levels with the reference, sets indicators to `0` rather than `NA` for missing covariate values, errors instead of silently overwriting when an indicator name already exists, deduplicates `catcovarsVec`, no longer creates all-zero indicator columns for unused factor levels, and drops expanded categorical names from the returned covariate vector (previously `updatedData[, updcovarsVec]` could error with "undefined columns selected" when a caller passed a categorical column in both `covarsVec` and `catcovarsVec`). Since the function is internal and unused by the `runSCM()` path, no user-facing search result changes.

* Auto-scaled default theta bounds for `lin`, `log`, and `identity` covariate shapes. Previously the shape-agnostic defaults of `(-5, 5)` could produce numerically catastrophic parameter spaces for physiologically-scaled covariates (e.g. weight in kg), causing both diverged optimizations and spurious minima. Bounds are now scaled by `1/|center|` for `lin` and `identity`, and by `1/|log(center)|` for `log`, giving results comparable to the dimensionless `power` shape.

* New retry mechanism for unrealistic OFV values in `runSCM()`. Six new arguments control the behaviour: `maxRetries` (default `3L`), `maxDeltaOFV` (default `Inf`), `retryPerturbSD` (default `0.5`), `retrySmallInit` (default `0.01`), `retryOFVTolerance` (default `NULL`, auto-detected), and `retryFailOnExhaustion` (default `FALSE`). When a candidate fit produces an OFV that is implausibly high or low, the fit is retried with perturbed or reduced initial estimates. Stochastic estimators (SAEM) get a 10-unit tolerance margin automatically to avoid spurious retries due to Monte Carlo noise.

* Fixed `.idColumn()` returning the wrong column when `id`/`ID` was not the first column of the dataset. `which("id" %in% colNamesLower)` always evaluated to `1` (or `integer(0)`) regardless of the matching column's position, causing `.enrichPairs()` to compute per-subject covariate medians from the wrong column. Replaced with `match("id", colNamesLower)`. Datasets where `id` happened to be the first column (e.g. `nlmixr2data::warfarin`) were unaffected.

* Fixed a bug in `runSCM(searchType = "backward", includedRelations = ...)` where every relation in `includedRelations` was silently appended to the candidate `pairs` table at every backward step, even when already present. The dedupe key for `pairs` lacked the `cov_` prefix carried by theta names, so no row was ever recognised as a duplicate; backward steps iterated over twice the intended candidates and wall-clock time roughly doubled. Both keys are now prefixed with `cov_`.

# nlmixr2scm 0.1

* `runSCM()` provides a standalone home for nlmixr2 stepwise covariate modeling, including SCM model-building helpers, categorical covariate expansion, resume support, and the migrated SCM vignette/tests.
