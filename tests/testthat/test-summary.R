# summary() / print() methods for runSCM() results.  These use lightweight
# mock fits (only the fields the methods read) so they run without fitting.

# Captures stdout and cli/messages together: the step tables reuse the
# end-of-run printer, which emits its rules through cli (message stream).
.capture_all <- function(expr) {
  msgs <- character()
  out <- utils::capture.output(withCallingHandlers(
    expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  paste(c(out, msgs), collapse = "\n")
}

.mock_fit <- function(objf, covs = character()) {
  thetas <- c("tka", "tcl", "tv", "prop.err", covs)
  list(
    objf = objf,
    AIC = objf + 10,
    BIC = objf + 20,
    iniDf = data.frame(name = thetas),
    finalUiEnv = list(ini = list(est = seq_along(thetas)))
  )
}

.mock_table <- function() {
  data.frame(
    step       = c(1L, 1L, 2L, 1L, 2L),
    covar      = c("wt_power", "sex_female", "wt_power", "wt_power", "wt_power"),
    var        = c("v", "cl", "cl", "cl", "v"),
    shape      = c("power", "cat", "power", "power", "power"),
    objf       = c(449.113, 474.173, 442.679, 448.298, 474.521),
    deltObjf   = c(-25.976, -0.917, -6.434, 5.619, 26.223),
    pchisqr    = c(3.5e-7, 0.338, 0.0112, 0.0178, 3.0e-7),
    included   = c("yes", "no", "yes", "dropped", "retained"),
    searchType = c("forward", "forward", "forward", "backward", "backward"),
    covNames   = c("cov_wt_power_v", "cov_sex_female_cl", "cov_wt_power_cl",
                   "cov_wt_power_cl", "cov_wt_power_v"),
    stringsAsFactors = FALSE
  )
}

.mock_options <- function(searchType = "scm", saveModels = FALSE) {
  list(
    searchType = searchType,
    pVal = list(fwd = 0.05, bck = 0.01),
    estimation = "focei",
    control = "foceiControl (inherited from fit)",
    candidates = data.frame(
      var = c("v", "cl"), covar = c("wt_power", "sex_female"),
      shape = c("power", "cat"), stringsAsFactors = FALSE
    ),
    includedRelations = NULL,
    shapes = "power",
    centers = c(wt = 70),
    catCutoff = 0.05,
    profileInit = FALSE,
    profileInitOnStall = TRUE,
    stallTol = 0,
    maxRetries = 3L,
    maxDeltaOFV = Inf,
    retryPerturbSD = 0.5,
    retrySmallInit = 0.01,
    retryOFVTolerance = 0,
    retryFailOnExhaustion = FALSE,
    retryOnUnderflow = TRUE,
    workers = 1L,
    rxThreads = 2L,
    saveModels = saveModels,
    restart = FALSE
  )
}

.mock_scm <- function(searchType = "scm", outputDir = NULL) {
  st <- .mock_table()
  if (searchType == "forward") st <- st[st$searchType == "forward", ]
  structure(
    list(
      summaryTable = st,
      resFwd = list(.mock_fit(442.679, c("cov_wt_power_v", "cov_wt_power_cl")),
                    st[st$searchType == "forward", ], NULL),
      resBck = if (searchType == "scm") {
        list(.mock_fit(448.298, "cov_wt_power_v"), st[st$searchType == "backward", ])
      },
      baseFit = .mock_fit(475.089),
      options = .mock_options(searchType, saveModels = !is.null(outputDir)),
      outputDir = outputDir
    ),
    class = c("nlmixr2scm", "list")
  )
}

test_that("summary() compares base, final forward and final backward models", {
  s <- summary(.mock_scm())
  expect_s3_class(s, "summary.nlmixr2scm")
  cmp <- s$comparison
  expect_equal(cmp$model, c("Base", "Forward final", "Backward final"))
  expect_equal(cmp$OFV, c(475.089, 442.679, 448.298))
  expect_equal(cmp$dOFV, c(0, 442.679 - 475.089, 448.298 - 475.089))
  expect_equal(cmp$nPar, c(4L, 6L, 5L))
  expect_equal(
    cmp$covariates,
    c("(none)", "wt_power~v [power], wt_power~cl [power]", "wt_power~v [power]")
  )
})

test_that("summary() omits the backward row for a forward-only search", {
  s <- summary(.mock_scm("forward"))
  expect_equal(s$comparison$model, c("Base", "Forward final"))
})

test_that("print(summary()) shows options, comparison, steps and no-files note", {
  out <- .capture_all(print(summary(.mock_scm())))
  expect_match(out, "Options")
  expect_match(out, "p-value fwd / bck")
  expect_match(out, "wt = 70", fixed = TRUE)
  expect_match(out, "Model comparison")
  expect_match(out, "Backward final")
  expect_match(out, "SCM Step Summary", fixed = TRUE)
  expect_match(out, "not saved", fixed = TRUE)
})

test_that("print(summary()) lists output directory, report files and saved models", {
  dir <- normalizePath(withr::local_tempdir()) # as runSCM() stores it
  writeLines("log", file.path(dir, "scm_log.txt"))
  saveRDS(1, file.path(dir, "fwd_1.rds"))
  saveRDS(2, file.path(dir, "fwd_2.rds"))
  out <- .capture_all(print(summary(.mock_scm(outputDir = dir))))
  expect_match(out, dir, fixed = TRUE)
  expect_match(out, "scm_log.txt", fixed = TRUE)
  expect_false(grepl("scm_step_summary.csv", out, fixed = TRUE))
  expect_match(out, "2 files (.rds)", fixed = TRUE)
})

test_that("print() gives a short overview pointing at summary()", {
  out <- .capture_all(print(.mock_scm()))
  expect_match(out, "475.089", fixed = TRUE)
  expect_match(out, "448.298", fixed = TRUE)
  expect_match(out, "wt_power~v [power]", fixed = TRUE)
  expect_match(out, "summary()", fixed = TRUE)
  expect_false(grepl("SCM All Candidates", out, fixed = TRUE))
  expect_invisible(print(.mock_scm()))
})
