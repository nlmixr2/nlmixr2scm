skip_if_not(exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE),
            "rxode2 has no event bus")

.evRec <- new.env()
.evListen <- function(env = parent.frame()) {
  .evRec$ev <- list()
  rxode2::rxEventListen("nlmixr2scm-test", function(event, ...) {
    .evRec$ev[[length(.evRec$ev) + 1L]] <- list(event = event, p = list(...))
  })
  withr::defer(rxode2::rxEventUnlisten("nlmixr2scm-test"), envir = env)
}
.evNames <- function() vapply(.evRec$ev, function(e) e$event, character(1))

## the base model of test-scm.R
.evBaseFit <- function() {
  .model <- function() {
    ini({
      tka <- log(1.5)
      tcl <- log(0.04)
      tv <- log(0.5)
      eta.ka ~ 0.5
      eta.cl ~ 0.2
      eta.v ~ 0.1
      prop.err <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ prop(prop.err)
    })
  }
  rxode2::rxEventScope(suppressMessages(nlmixr2est::nlmixr2(
    .model, nlmixr2data::theo_sd, est = "focei",
    control = nlmixr2est::foceiControl(print = 0, calcTables = TRUE)
  )))
}

test_that("runSCM emits one summary, then the selected model linked to the base fit", {
  skip_on_cran()
  withr::local_dir(withr::local_tempdir())
  base_fit <- .evBaseFit()
  .evListen()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "backward",
    saveModels = FALSE,
    includedRelations = list(list(var = "cl", covar = "WT", shapes = "power")),
    workers = 1L
  )
  expect_identical(.evNames(), c("fitResult", "fitComplete"))
  .sum <- .evRec$ev[[1]]$p
  expect_identical(.sum$kind, "scm")
  expect_identical(.sum$fit, base_fit)
  expect_s3_class(.sum$result, "nlmixr2scmSummary")
  expect_false(any(vapply(.sum$result, function(x) inherits(x, "nlmixr2FitCore"), TRUE)))
  expect_true(object.size(.sum$result) < 1e6)
  .fin <- .evRec$ev[[2]]$p
  expect_identical(.fin$source, "scm")
  expect_identical(.fin$object, base_fit)
  expect_identical(.fin$fit, res$resBck[[1]])
  expect_identical(rxode2::rxEventDepth(), 0L)
})

test_that("no fitComplete when SCM keeps the base model; nothing on error", {
  .evListen()
  fit <- structure(list(env = new.env()), class = c("nlmixr2FitCore", "list"))
  .scmEventEnter()
  .scmEventExit(list(summaryTable = data.frame(step = 1), resFwd = list(fit, NULL), resBck = NULL),
                fit, quote(runSCM(fit)))
  expect_identical(.evNames(), "fitResult")
  .scmEventEnter()
  .scmEventExit(NULL, fit, NULL)
  expect_length(.evRec$ev, 1L)
  expect_identical(rxode2::rxEventDepth(), 0L)
})
