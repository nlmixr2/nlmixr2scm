# nlmixr2scm

`nlmixr2scm` provides stepwise covariate modeling (SCM) for `nlmixr2`
population PK/PD models.

SCM is a likelihood-ratio-based covariate selection workflow that proceeds in
two phases:

1. _Forward inclusion:_ starting from a base model, add the candidate
   covariate-parameter relationship that gives the largest statistically
   significant OFV improvement.
2. _Backward elimination:_ starting from the forward-final model, test each
   included relationship for removal using a stricter threshold.

The result is a reproducible covariate-search workflow that can summarize every
candidate tested, keep accepted and dropped relationships explicit, and return
the final forward and backward models for follow-up review. It's based heavily on the approach of [Perl-speaks-NONMEM](https://github.com/UUPharmacometrics/PsN). 

## Details

`runSCM()` implements the SCM workflow for `nlmixr2` fits. It automatically
builds the covariate terms inside the model body, so continuous covariates can
be centered at their observed medians and categorical covariates can be turned
into indicator columns without having to edit the model by hand.

The package supports:

* full covariate searches (via `varsVec` and `covarsVec`)
* exact candidate relationship specifications (via `pairsVec`)
* multiple built-in continuous-covariate shapes (`"power"`, `"lin"`, `"log"`,
  and `"identity"`)
* fixed (rather than per-dataset-median) covariate centering via `centers`
* automatic categorical covariate preprocessing (through `catvarsVec`)
* optional pre-included relationships (through `includedRelations`)
* saved step tables and fitted candidate models for resuming/reviewing workflows
* optional parallel fitting of candidates through `workers` (based on `future.apply`)
* automatic retries and profile warm-starting for candidates that produce an
  unrealistic OFV or stall at their initial estimate (see "Search robustness" below)

`nlmixr2scm` relies on `nlmixr2utils`, which provides 
shared worker-plan helpers and supporting infrastructure.

## Search robustness

Covariate searches routinely run into candidate relationships whose estimation misbehaves - 
solver noise, stalling at initial
estimates, or a fit landing on an implausible OFV. `runSCM()` has three features
that keep these numerical hiccups from corrupting the search's accept/reject
decisions, enabled by sensible defaults, so (hopefully) little to no tuning is required to benefit
from them.

**Profile warm-starting.** With ODE models in particular, gradient optimizers
(`nlminb`, `lbfgsb3c`) can stall at a covariate's zero-effect initial
estimate when solver noise flattens the outer objective: the fit "converges"
having never actually explored the covariate effect, which would otherwise show
up as a false negative in the search. By default (`profileInitOnStall = TRUE`),
`runSCM()` detects a stalled candidate (its OFV improvement over the parent
model is `<= stallTol`, which cannot happen at a true optimum) and rescues it
with a cheap, frozen 1-D profile (Brent method); the structural parameters and
between-subject variability are fixed at the parent's estimates while only
the covariate coefficient is optimized, giving the estimator a non-zero,
gradient-informative starting value for the refit. Set `profileInit = TRUE` to
warm-start every forward candidate this way, not just stalled ones.

**Automatic retries.** A candidate fit can occasionally land on an implausible
OFV, creating too large a jump to be a genuine covariate effect. `runSCM()` retries
such candidates (up to `maxRetries`, default 3) using a perturbed or near-zero
covariate initial estimate on each attempt, then keeps the best attempt seen
across all retries rather than whichever happened to run last. Stochastic
estimators (SAEM) get a wider OFV tolerance automatically, since Monte Carlo
noise alone would otherwise trigger spurious retries. A p-value that underflows
to zero is treated as implausible by default, so very strong effects are retried
too; set `retryOnUnderflow = FALSE` to skip those extra fits.

**Fixed covariate centers.** By default, the `"power"` and `"lin"` shapes
center on each dataset's observed median, which can drift between datasets fit
with the same model (e.g. simulation studies, or re-fitting as new data
arrives). Pass `centers` (e.g. `centers = c(wt = 70)`) to fix the reference
value for one or more continuous covariates instead, so the estimated
coefficients -- and the structural intercept -- stay on a consistent reference
across runs.

See the vignette (`vignette("runSCM", package = "nlmixr2scm")`) for the full
argument reference and further tips for a robust search.

## Installation

Install the released version from CRAN:

```r
install.packages("nlmixr2scm")
```

Or the development version from GitHub, together with `nlmixr2utils`:

```r
remotes::install_github("nlmixr2/nlmixr2utils")
remotes::install_github("nlmixr2/nlmixr2scm")
```

## Basic Use

Here's a quick example:

```r
library(nlmixr2)
library(nlmixr2data)
library(nlmixr2scm)

pkdata <- nlmixr2data::warfarin[nlmixr2data::warfarin$dvid == "cp", ]

warf_pk <- function() {
  ini({
    tka <- log(1.15)
    tcl <- log(0.135)
    tv <- log(7.0)
    eta.ka ~ 0.40
    eta.cl ~ 0.25
    eta.v ~ 0.10
    prop.err <- 0.10
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ prop(prop.err)
  })
}

fit_base <- nlmixr2(
  warf_pk,
  data = pkdata,
  est = "focei",
  control = nlmixr2est::foceiControl(print = 0),
  table = nlmixr2est::tableControl(cwres = TRUE)
)

scm <- runSCM(
  fit = fit_base,
  data = pkdata,
  varsVec = c("cl", "v"),
  covarsVec = "wt",
  catvarsVec = "sex",
  shapes = c("power","lin"),
  pVal = list(fwd = 0.05, bck = 0.01),
  searchType = "scm",
  saveModels = FALSE,
  print = 0,
  workers = 1L,
  confirm = FALSE
)

scm                  # short overview
summary(scm)         # options, model comparison, step tables, files on disk
scm$summaryTable
scm$resBck[[1]]$parFixedDf
```

## Requirements and Practical Notes

Pass the original, untransformed data to `runSCM()`. Centering, indicator
generation, and any other necessary covariate-shape transformations are generated inside the
SCM tool.

For practical use:

* Pass `data` explicitly whenever the base model does not already reference a
  covariate column that SCM needs to test.
* Use `catvarsVec` for categorical predictors instead of creating dummy columns
  manually.
* Use `pairsVec` when you want exact control over which relationships are
  tested, rather than the full parameter x covariate search space.
* Use `includedRelations` when backward elimination should start from a model
  that contains required relationships beyond what forward inclusion selected.
* Use `saveModels = FALSE` for quick exploratory runs; leave it at the default
  when you want saved step tables, fitted candidates, and resume support.
* Use `confirm = FALSE` in scripts, tests, or non-interactive workflows.

## References

* Jonsson, E.N. & Karlsson, M.O. (1998). Automated covariate model building
  within NONMEM. *Pharmaceutical Research*, 15(9), 1463-1468.
* Lindbom, L., Ribbing, J., & Jonsson, E.N. (2004). Perl-speaks-NONMEM (PsN):
  a Perl module for NONMEM related programming. *Computer Methods and Programs
  in Biomedicine*, 75, 85-94.
* Ribbing, J. & Jonsson, E.N. (2004). Power, selection bias and predictive
  performance of the population pharmacokinetic covariate model. *Journal of
  Pharmacokinetics and Pharmacodynamics*, 31(2), 109-134.

For a fuller worked example, see the package vignette:
`vignette("runSCM", package = "nlmixr2scm")`.

## Credit where it's due

`nlmixr2scm` is based on the [PsN implementation](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/scm_userguide.pdf) of SCM for NONMEM written by Lars Lindbom, 
Niclas Jonsson, Pontus Pihlgren, Mats Karlsson, Andrew Hooker, Kajsa Harling, 
Rikard Nordgren and Svetlana Freiberga.