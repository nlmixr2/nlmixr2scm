## Update to 0.4.1

This is a patch release shortly after 0.4 (published 2026-09-24). It fixes a
correctness bug in 0.4: the forward search mixed two sign conventions for
the change in objective function value. As a result, every forward candidate
was needlessly refitted (roughly tripling run time), a worse refit could
replace a better fit, and the reported summary tables were internally
inconsistent. We would rather not leave users on a release that can report
misleading results, hence the quick update. The release also adds a
`retryOnUnderflow` argument and `summary()`/`print()` methods for results.

## Test environments

* Local: Windows 11 x64, R 4.6.1
* devtools::check_win_devel()
* rhub::rhub_check() (r-devel on ubuntu-latest)

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Justin Wilkins <justin.wilkins@occams.com>'
Days since last update: <n>
```

The short interval since 0.4 (published 2026-09-24) is explained above.

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Example timing

With `--run-donttest`, the examples for `runSCM()` and `summary.nlmixr2scm()`
(which shares the same example) may exceed 5 seconds. `runSCM()` performs a
stepwise covariate search, so its example has to fit a base population model
and then fit one candidate model per step. There is no way to exercise the
function meaningfully in under 5 seconds. The example has been kept as small
as is still representative (one parameter-covariate pair, forward direction
only, two threads) and is wrapped in `\donttest{}`.

## Method references

The methods implemented here are described in the two references cited in the
`Description` field:

* Jonsson EN, Karlsson MO (1998). Automated covariate model building within
  NONMEM. *Pharmaceutical Research* 15, 1463-1468.
  <doi:10.1023/A:1011970125687>
* Lindbom L, Ribbing J, Jonsson EN (2004). Perl-speaks-NONMEM (PsN) - a Perl
  module for NONMEM related programming. *Computer Methods and Programs in
  Biomedicine* 75, 85-94. <doi:10.1016/j.cmpb.2003.11.003>

## Notes for the reviewer

* **Examples.** `runSCM()` requires a fitted population model, so its example
  first fits a base model and then runs a one-pair covariate search. That takes
  longer than a few seconds, so the example is wrapped in `\donttest{}` rather
  than `\dontrun{}`; it is fully executable. It is additionally guarded with
  `@examplesIf requireNamespace("nlmixr2data")` because the example data come
  from a suggested package.

* **Vignette.** The vignette demonstrates a complete forward/backward search,
  which involves fitting many models. It is therefore pre-computed: the
  executable source is kept outside the build as `vignettes/runSCM.Rmd.orig`
  and knitted to a static `vignettes/runSCM.Rmd` by `vignettes/precompute.R`,
  so no models are re-fitted during `R CMD check`.

* **Parallelism.** No checked code uses more than two cores. The shipped
  vignette is static, so it fits nothing; its executable source caps threads
  with `rxode2::setRxThreads(2L)` and `rxThreads = 2L` for when it is
  re-generated. The example passes `rxThreads = 2L` and `workers = 1L`, and
  the tests that fit models are skipped on CRAN.

* **Writing to disk.** `runSCM()` can cache fitted models, but does so only
  when the user opts in via `saveModels = TRUE` (and, in interactive sessions,
  confirms the run). Nothing in the examples, tests, or vignette writes outside
  `tempdir()`.
