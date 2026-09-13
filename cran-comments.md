## Test environments

* Local: Windows 11 x64, R 4.6.1

<!-- Before submitting, add results from at least:
     devtools::check_win_devel(), devtools::check_win_release(),
     and a macOS / Linux run (e.g. rhub::rhub_check() or GitHub Actions). -->

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Justin Wilkins <justin.wilkins@occams.com>'
New submission
```

This is a new submission.

With `--run-donttest` there is a second note:

```
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
       user system elapsed
runSCM  8.3   1.09   21.17
```

`runSCM()` performs a stepwise covariate search, so its example has to fit a
base population model and then fit one candidate model per step. There is no
way to exercise the function meaningfully in under 5 seconds. The example has
been kept as small as is still representative (one parameter-covariate pair,
forward direction only, two threads), and the complete check still runs in
about one minute.

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
