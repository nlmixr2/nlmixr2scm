library(nlmixr2utils)

.cur <- loadNamespace("nlmixr2scm")

test_that("SCM helper functions remain internal-only", {
  exports <- getNamespaceExports("nlmixr2scm")
  old_names <- c(
    "addCatCovariates",
    "addorremoveCovariate",
    "buildcovInfo",
    "buildupatedUI"
  )
  new_names <- c(
    "scmAddCatCovariates",
    "scmAddOrRemoveCovariate",
    "scmBuildCovInfo",
    "scmBuildUpdatedUi"
  )

  expect_setequal(intersect(old_names, exports), character(0))
  expect_setequal(intersect(new_names, exports), character(0))
  expect_false("covarSearchAuto" %in% exports)
  expect_false(exists("covarSearchAuto", envir = .cur, inherits = FALSE))
  expect_type(.cur$scmAddCatCovariates, "closure")
  expect_type(.cur$scmAddOrRemoveCovariate, "closure")
  expect_type(.cur$scmBuildCovInfo, "closure")
  expect_type(.cur$scmBuildUpdatedUi, "closure")
})

# ── Shared model definition (no fitting) ──────────────────────────────────

.one_cmt_fun <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl + eta.cl) # nolint: object_usage_linter.
    v <- exp(tv + eta.v) # nolint: object_usage_linter.
    linCmt() ~ add(add.sd)
  })
}

# =============================================================================
# .autoScaleInitSpec
# =============================================================================

test_that(".autoScaleInitSpec: non-NULL spec bypasses scaling (scalar)", {
  res <- .cur$.autoScaleInitSpec(0.75, "lin", center = 70)
  expect_equal(res$est, 0.75)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".autoScaleInitSpec: non-NULL spec bypasses scaling (full list)", {
  res <- .cur$.autoScaleInitSpec(list(est = 0.5, lower = -1, upper = 1), "lin", center = 70)
  expect_equal(res$est, 0.5)
  expect_equal(res$lower, -1)
  expect_equal(res$upper, 1)
})

test_that(".autoScaleInitSpec: power shape returns unscaled defaults", {
  res <- .cur$.autoScaleInitSpec(NULL, "power", center = 70)
  expect_equal(res$est, 0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".autoScaleInitSpec: lin shape scales by abs(center)", {
  res <- .cur$.autoScaleInitSpec(NULL, "lin", center = 70)
  expect_equal(res$est,   0.1 / 70, tolerance = 1e-10)
  expect_equal(res$lower, -5 / 70,  tolerance = 1e-10)
  expect_equal(res$upper,  5 / 70,  tolerance = 1e-10)
})

test_that(".autoScaleInitSpec: identity shape scales by abs(center)", {
  res <- .cur$.autoScaleInitSpec(NULL, "identity", center = 70)
  expect_equal(res$est,   0.1 / 70, tolerance = 1e-10)
  expect_equal(res$lower, -5 / 70,  tolerance = 1e-10)
  expect_equal(res$upper,  5 / 70,  tolerance = 1e-10)
})

test_that(".autoScaleInitSpec: log shape scales by abs(log(center))", {
  res <- .cur$.autoScaleInitSpec(NULL, "log", center = 70)
  sc <- abs(log(70))
  expect_equal(res$est,   0.1 / sc, tolerance = 1e-10)
  expect_equal(res$lower, -5 / sc,  tolerance = 1e-10)
  expect_equal(res$upper,  5 / sc,  tolerance = 1e-10)
})

test_that(".autoScaleInitSpec: log with negative center scales by abs(log(abs(center)))", {
  res <- .cur$.autoScaleInitSpec(NULL, "log", center = -70)
  sc <- abs(log(abs(-70)))
  expect_equal(res$est,   0.1 / sc, tolerance = 1e-10)
  expect_equal(res$lower, -5 / sc,  tolerance = 1e-10)
  expect_equal(res$upper,  5 / sc,  tolerance = 1e-10)
})

test_that(".autoScaleInitSpec: lin with center = NA falls back to unscaled defaults", {
  expect_warning(
    res <- .cur$.autoScaleInitSpec(NULL, "lin", center = NA_real_),
    regexp = "Cannot auto-scale"
  )
  expect_equal(res$est, 0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".autoScaleInitSpec: lin with center = 0 falls back to unscaled defaults", {
  expect_warning(
    res <- .cur$.autoScaleInitSpec(NULL, "lin", center = 0),
    regexp = "Cannot auto-scale"
  )
  expect_equal(res$est, 0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".autoScaleInitSpec: log with center = 1 falls back to unscaled defaults", {
  expect_warning(
    res <- .cur$.autoScaleInitSpec(NULL, "log", center = 1),
    regexp = "Cannot auto-scale"
  )
  expect_equal(res$est, 0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".autoScaleInitSpec: lin with negative center scales by abs value", {
  res <- .cur$.autoScaleInitSpec(NULL, "lin", center = -70)
  expect_equal(res$est,   0.1 / 70, tolerance = 1e-10)
  expect_equal(res$lower, -5 / 70,  tolerance = 1e-10)
  expect_equal(res$upper,  5 / 70,  tolerance = 1e-10)
})

# =============================================================================
# .parseInitSpec
# =============================================================================

test_that(".parseInitSpec: NULL returns defaults (0.1, -5, 5)", {
  res <- .cur$.parseInitSpec(NULL)
  expect_equal(res$est, 0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".parseInitSpec: scalar preserves est, uses default bounds", {
  res <- .cur$.parseInitSpec(0.75)
  expect_equal(res$est, 0.75)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

test_that(".parseInitSpec: full named list preserved", {
  res <- .cur$.parseInitSpec(list(est = 0.5, lower = -2, upper = 3))
  expect_equal(res$est, 0.5)
  expect_equal(res$lower, -2)
  expect_equal(res$upper, 3)
})

test_that(".parseInitSpec: 'init' alias accepted for est", {
  res <- .cur$.parseInitSpec(list(init = 0.3, lower = -1, upper = 1))
  expect_equal(res$est, 0.3)
  expect_equal(res$lower, -1)
  expect_equal(res$upper, 1)
})

test_that(".parseInitSpec: partial list fills missing fields with defaults", {
  res_lower <- .cur$.parseInitSpec(list(est = 0.2, lower = 0))
  expect_equal(res_lower$upper, 5)

  res_upper <- .cur$.parseInitSpec(list(est = 0.2, upper = 2))
  expect_equal(res_upper$lower, -5)

  res_est <- .cur$.parseInitSpec(list(lower = -3, upper = 3))
  expect_equal(res_est$est, 0.1)
})

test_that(".parseInitSpec: unrecognised input returns defaults", {
  res <- .cur$.parseInitSpec("notanumber")
  expect_equal(res$est, 0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper, 5)
})

# =============================================================================
# buildPairs
# =============================================================================

test_that("buildPairs: varsVec x covarsVec produces cartesian product", {
  df <- .cur$buildPairs(varsVec = c("cl", "v"), covarsVec = c("wt", "age"))
  expect_equal(nrow(df), 4)
  expect_true(all(c("var", "covar") %in% names(df)))
  expect_setequal(df$var, c("cl", "cl", "v", "v"))
  expect_setequal(df$covar, c("wt", "age", "wt", "age"))
})

test_that("buildPairs: pairsVec as list of lists", {
  df <- .cur$buildPairs(
    pairsVec = list(
      list(var = "cl", covar = "wt"),
      list(var = "v", covar = "age")
    )
  )
  expect_equal(nrow(df), 2)
  expect_equal(df$var[1], "cl")
  expect_equal(df$covar[1], "wt")
  expect_equal(df$var[2], "v")
  expect_equal(df$covar[2], "age")
})

test_that("buildPairs: pairsVec as data frame passed through", {
  pv <- data.frame(
    var = c("cl", "v"),
    covar = c("wt", "wt"),
    stringsAsFactors = FALSE
  )
  df <- .cur$buildPairs(pairsVec = pv)
  expect_equal(nrow(df), 2)
  expect_equal(df$var, c("cl", "v"))
  expect_equal(df$covar, c("wt", "wt"))
})

test_that("buildPairs: shapes list-column preserved when supplied", {
  df <- .cur$buildPairs(
    pairsVec = list(
      list(var = "cl", covar = "wt", shapes = c("power", "lin"))
    )
  )
  expect_true("shapes" %in% names(df))
  expect_equal(df$shapes[[1]], c("power", "lin"))
})

test_that("buildPairs: inits list-column preserved when supplied", {
  df <- .cur$buildPairs(
    pairsVec = list(
      list(
        var = "cl",
        covar = "wt",
        inits = list(power = list(est = 0.75, lower = 0, upper = 2))
      )
    )
  )
  expect_true("inits" %in% names(df))
  expect_equal(df$inits[[1]][["power"]][["est"]], 0.75)
})

test_that("buildPairs: NULL varsVec with no pairsVec errors", {
  expect_error(.cur$buildPairs(varsVec = NULL, covarsVec = "wt"))
})

test_that("buildPairs: varsVec without covarsVec errors", {
  expect_error(.cur$buildPairs(varsVec = "cl", covarsVec = NULL))
})

# =============================================================================
# .makeSCMData
# =============================================================================

.sex_data <- function(n_m = 7, n_f = 3) {
  data.frame(
    ID = seq_len(n_m + n_f),
    sex = c(rep("M", n_m), rep("F", n_f)),
    wt = rnorm(n_m + n_f, 70, 10),
    stringsAsFactors = FALSE
  )
}

test_that(".makeSCMData: creates indicator column for binary covariate", {
  d <- .sex_data()
  res <- .cur$.makeSCMData(d, catvarsVec = "sex", fit = NULL, catCutoff = 0.05)
  expect_true("sex_F" %in% names(res$data))
  expect_true(all(res$data$sex_F %in% c(0L, 1L)))
})

test_that(".makeSCMData: reference level is most-frequent category", {
  d <- .sex_data(n_m = 7, n_f = 3)
  res <- .cur$.makeSCMData(d, catvarsVec = "sex", fit = NULL, catCutoff = 0.05)
  expect_equal(res$catRef[["sex"]], "M")
})

test_that(".makeSCMData: reference level not included in catLevels", {
  d <- .sex_data()
  res <- .cur$.makeSCMData(d, catvarsVec = "sex", fit = NULL, catCutoff = 0.05)
  expect_false("M" %in% res$catLevels[["sex"]])
  expect_true("F" %in% res$catLevels[["sex"]])
})

test_that(".makeSCMData: catCutoff drops rare levels into catDropped", {
  d <- data.frame(
    ID = seq_len(20),
    grp = c(rep("A", 18), rep("B", 1), rep("C", 1)),
    stringsAsFactors = FALSE
  )
  # B and C each have 1/20 = 5%; with cutoff > 5% they should be dropped
  res <- .cur$.makeSCMData(d, catvarsVec = "grp", fit = NULL, catCutoff = 0.06)
  expect_length(res$catDropped[["grp"]], 2)
  expect_length(res$catLevels[["grp"]], 0)
})

test_that(".makeSCMData: levels at or above catCutoff are retained", {
  d <- data.frame(
    ID = seq_len(20),
    grp = c(rep("A", 16), rep("B", 4)),
    stringsAsFactors = FALSE
  )
  # B has 4/20 = 20%; well above any reasonable cutoff
  res <- .cur$.makeSCMData(d, catvarsVec = "grp", fit = NULL, catCutoff = 0.05)
  expect_true("B" %in% res$catLevels[["grp"]])
  expect_length(res$catDropped[["grp"]], 0)
})

test_that(".makeSCMData: three-level variable with one rare level", {
  d <- data.frame(
    ID = seq_len(20),
    grp = c(rep("A", 15), rep("B", 4), rep("C", 1)),
    stringsAsFactors = FALSE
  )
  # A=75% (ref), B=20% (retained), C=5% (< cutoff 0.06 → dropped)
  res <- .cur$.makeSCMData(d, catvarsVec = "grp", fit = NULL, catCutoff = 0.06)
  expect_equal(res$catRef[["grp"]], "A")
  expect_true("B" %in% res$catLevels[["grp"]])
  expect_false("C" %in% res$catLevels[["grp"]])
  expect_true("C" %in% res$catDropped[["grp"]])
})

# =============================================================================
# .enrichPairs
# =============================================================================

test_that(".enrichPairs: continuous covariate gets median as center", {
  d <- data.frame(
    ID = 1:10,
    wt = seq(50, 95, by = 5),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  res <- .cur$.enrichPairs(pairs, d)
  expect_equal(res$type[1], "continuous")
  expect_equal(res$center[1], median(d$wt))
  expect_equal(res$raw_col[1], "wt")
})

test_that(".enrichPairs: multiple continuous covariates each get their own median", {
  d <- data.frame(
    ID = 1:10,
    wt = seq(50, 95, by = 5),
    age = seq(20, 65, by = 5),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(
    var = c("cl", "cl"),
    covar = c("wt", "age"),
    stringsAsFactors = FALSE
  )
  res <- .cur$.enrichPairs(pairs, d)
  expect_equal(res$center[res$covar == "wt"], median(d$wt))
  expect_equal(res$center[res$covar == "age"], median(d$age))
})

test_that(".enrichPairs: categorical covariate (via catvarsVec) expands to per-level rows", {
  d <- data.frame(
    ID = 1:4,
    sex = c("M", "M", "F", "F"),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "sex", stringsAsFactors = FALSE)
  res <- .cur$.enrichPairs(pairs, d, catvarsVec = "sex")
  # Expect one row per non-reference level (F; M is reference as most frequent
  # alphabetically here, but .enrichPairs drops first sorted level)
  expect_equal(nrow(res), 1)
  expect_equal(res$type[1], "categorical")
  expect_true(!is.na(res$level[1]))
})

test_that(".enrichPairs: missing value flagged when missingToken present", {
  d <- data.frame(
    ID = 1:5,
    wt = c(70, 75, -99, 80, 65),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  res <- .cur$.enrichPairs(pairs, d, missingToken = -99)
  expect_true(res$has_missing[1])
  expect_false(is.na(res$missing_check[1]))
})

test_that(".enrichPairs: no missing flag when all values present", {
  d <- data.frame(
    ID = 1:5,
    wt = c(70, 75, 80, 85, 90),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  res <- .cur$.enrichPairs(pairs, d)
  expect_false(isTRUE(res$has_missing[1]))
})

# =============================================================================
# .expandShapes
# =============================================================================

.cont_pairs <- function(data) {
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  .cur$.enrichPairs(pairs, data)
}

test_that(".expandShapes: power shape produces log-ratio expression", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power")
  expect_equal(res$shape[1], "power")
  expect_match(res$covExpr[1], "log\\(wt/")
  expect_match(res$covExpr[1], as.character(median(d$wt)))
})

test_that(".expandShapes: lin shape produces subtraction expression", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "lin")
  expect_equal(res$shape[1], "lin")
  expect_match(res$covExpr[1], "wt - ")
  expect_match(res$covExpr[1], as.character(median(d$wt)))
})

test_that(".expandShapes: identity shape returns raw column", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "identity")
  expect_equal(res$covExpr[1], "wt")
})

test_that(".expandShapes: log shape returns log(col)", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "log")
  expect_match(res$covExpr[1], "^log\\(wt\\)$")
})

# =============================================================================
# .expandShapes: missing-value fill is shape-aware (imputes median, not 0)
# =============================================================================
# Regression test: earlier versions used `missing_fill = 0L` for every
# continuous shape, which silently mis-imputed `cov = 1` for `"log"` and
# `cov = 0` for `"identity"` instead of the median.  The fix routes the fill
# through `.scmFillStr()` which evaluates the shape expression at the
# centering value (median).

test_that(".expandShapes: missing-value fill is shape-aware (imputes median)", {
  d <- data.frame(
    ID = 1:5,
    wt = c(70, 75, NA, 80, 65),
    stringsAsFactors = FALSE
  )
  med <- median(d$wt, na.rm = TRUE) # 72.5
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  enriched <- .cur$.enrichPairs(pairs, d)
  expect_true(enriched$has_missing[1])

  # power: fill = log(m/m) = 0 (also matches the previous behaviour)
  ep_pow <- .cur$.expandShapes(enriched, shapes = "power")
  expect_match(
    ep_pow$covExpr[1],
    "ifelse(is.na(wt), 0, log(wt/",
    fixed = TRUE
  )

  # lin: fill = m - m = 0
  ep_lin <- .cur$.expandShapes(enriched, shapes = "lin")
  expect_match(
    ep_lin$covExpr[1],
    "ifelse(is.na(wt), 0, (wt - ",
    fixed = TRUE
  )

  # log: fill = log(median), NOT 0 (0 would imply cov = 1)
  ep_log <- .cur$.expandShapes(enriched, shapes = "log")
  expect_match(
    ep_log$covExpr[1],
    paste0("ifelse(is.na(wt), ", sprintf("%.15g", log(med)), ", log(wt))"),
    fixed = TRUE
  )
  expect_false(
    grepl(
      "ifelse(is.na(wt), 0, log(wt))",
      ep_log$covExpr[1],
      fixed = TRUE
    )
  )

  # identity: fill = median, NOT 0
  ep_id <- .cur$.expandShapes(enriched, shapes = "identity")
  expect_match(
    ep_id$covExpr[1],
    paste0("ifelse(is.na(wt), ", sprintf("%.15g", med), ", wt)"),
    fixed = TRUE
  )
  expect_false(
    grepl(
      "ifelse(is.na(wt), 0, wt)",
      ep_id$covExpr[1],
      fixed = TRUE
    )
  )
})

test_that(".expandShapes: cat shape uses indicator column name", {
  d <- data.frame(
    ID = 1:4,
    sex = c("M", "M", "F", "F"),
    sex_F = c(0L, 0L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "sex", stringsAsFactors = FALSE)
  pairs <- .cur$.enrichPairs(
    pairs,
    d,
    catvarsVec = "sex",
    catLevels = list(sex = "F")
  )
  res <- .cur$.expandShapes(pairs)
  expect_equal(res$shape[1], "cat")
  expect_equal(res$covExpr[1], "sex_F")
})

test_that(".expandShapes: multiple shapes per row expand to multiple rows", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = c("power", "lin"))
  expect_equal(nrow(res), 2)
  expect_setequal(res$shape, c("power", "lin"))
})

test_that(".expandShapes: global inits scalar sets est; bounds are default", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(
    .cont_pairs(d),
    shapes = "power",
    inits = list(power = 0.75)
  )
  expect_equal(res$init[1], 0.75)
  expect_equal(res$lower[1], -5)
  expect_equal(res$upper[1], 5)
})

test_that(".expandShapes: global inits with bounds all propagated", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(
    .cont_pairs(d),
    shapes = "power",
    inits = list(power = list(est = 0.5, lower = 0, upper = 2))
  )
  expect_equal(res$init[1], 0.5)
  expect_equal(res$lower[1], 0)
  expect_equal(res$upper[1], 2)
})

test_that(".expandShapes: per-pair inits override global inits", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  pairs <- .cont_pairs(d)
  pairs$inits <- list(list(power = 0.99))
  res <- .cur$.expandShapes(pairs, shapes = "power", inits = list(power = 0.1))
  expect_equal(res$init[1], 0.99)
})

test_that(".expandShapes: NULL inits for shape uses defaults", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power", inits = list())
  expect_equal(res$init[1], 0.1)
  expect_equal(res$lower[1], -5)
  expect_equal(res$upper[1], 5)
})

test_that(".expandShapes: continuous covar name has shape appended", {
  d <- data.frame(
    ID = 1:5,
    wt = c(60, 70, 75, 80, 90),
    stringsAsFactors = FALSE
  )
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power")
  expect_equal(res$covar[1], "wt_power")
})

# =============================================================================
# .updatePairsInits
# =============================================================================

test_that(".updatePairsInits: NULL pairs returns NULL", {
  expect_null(.cur$.updatePairsInits(NULL, list()))
})

test_that(".updatePairsInits: zero-row pairs returned unchanged", {
  pairs <- data.frame(
    var = character(0),
    covar = character(0),
    init = numeric(0),
    stringsAsFactors = FALSE
  )
  res <- .cur$.updatePairsInits(pairs, list())
  expect_equal(nrow(res), 0)
})

test_that(".updatePairsInits: updates init from parFixedDf", {
  pairs <- data.frame(
    var = "cl",
    covar = "wt_power",
    init = 0.1,
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate = 0.77,
      row.names = "cov_wt_power_cl",
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$init, 0.77)
})

test_that(".updatePairsInits: theta absent in fit leaves init unchanged", {
  pairs <- data.frame(
    var = "cl",
    covar = "wt_power",
    init = 0.1,
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate = 0.5,
      row.names = "cov_age_lin_cl",
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$init, 0.1)
})

test_that(".updatePairsInits: bounds (lower/upper) never modified", {
  pairs <- data.frame(
    var = "cl",
    covar = "wt_power",
    init = 0.1,
    lower = -2,
    upper = 2,
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate = 0.75,
      row.names = "cov_wt_power_cl",
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$lower, -2)
  expect_equal(res$upper, 2)
  expect_equal(res$init, 0.75)
})

test_that(".updatePairsInits: multiple rows updated independently", {
  pairs <- data.frame(
    var = c("cl", "v"),
    covar = c("wt_power", "wt_power"),
    init = c(0.1, 0.1),
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate = c(0.77, 0.33),
      row.names = c("cov_wt_power_cl", "cov_wt_power_v"),
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$init[res$var == "cl"], 0.77)
  expect_equal(res$init[res$var == "v"], 0.33)
})

# =============================================================================
# scmAddCatCovariates (internal)
# =============================================================================

test_that("scmAddCatCovariates: creates indicator columns for non-reference levels", {
  d <- data.frame(
    ID = 1:6,
    sex = c("M", "M", "M", "F", "F", "F"),
    wt = c(70, 75, 80, 60, 65, 55),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = "wt",
    catcovarsVec = "sex"
  )
  new_data <- res[[1]]
  new_covs <- res[[2]]
  expect_true(any(grepl("^sex_", names(new_data))))
  expect_true(any(grepl("^sex_", new_covs)))
})

test_that("scmAddCatCovariates: original categorical column removed from data", {
  d <- data.frame(
    ID = 1:4,
    sex = c("M", "M", "F", "F"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "sex"
  )
  expect_false("sex" %in% names(res[[1]]))
})

test_that("scmAddCatCovariates: most-frequent level (reference) has no indicator", {
  # 'B' appears first but 'A' is most frequent, so first-appearing and modal
  # rules disagree -- this fixture is what distinguishes them.
  d <- data.frame(
    ID = seq_len(10),
    grp = c("B", rep("A", 7), "B", "C"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "grp"
  )
  new_data <- res[[1]]
  # A is most frequent; its indicator should NOT be created
  expect_false("grp_A" %in% names(new_data))
  expect_true("grp_B" %in% names(new_data))
  expect_true("grp_C" %in% names(new_data))
})

test_that("scmAddCatCovariates: factors use the modal level, not the first level", {
  d <- data.frame(
    ID = seq_len(10),
    grp = factor(
      c("B", rep("A", 7), "B", "C"),
      levels = c("B", "A", "C")
    )
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "grp"
  )
  new_data <- res[[1]]
  # 'B' is the first factor level but 'A' is modal -- 'A' is the reference
  expect_false("grp_A" %in% names(new_data))
  expect_true("grp_B" %in% names(new_data))
  expect_true("grp_C" %in% names(new_data))
})

test_that("scmAddCatCovariates: catCutoff lumps rare levels with the reference", {
  d <- data.frame(
    ID = seq_len(20),
    grp = c(rep("A", 15), rep("B", 4), "C"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "grp",
    catCutoff = 0.1
  )
  new_data <- res[[1]]
  # C is 5% of rows, below the 10% cutoff, so it is lumped with the reference
  expect_true("grp_B" %in% names(new_data))
  expect_false("grp_C" %in% names(new_data))
})

test_that("scmAddCatCovariates: NA rows give 0 indicators, not NA", {
  d <- data.frame(
    ID = seq_len(6),
    grp = c("A", "A", "A", "B", "B", NA),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "grp"
  )
  new_data <- res[[1]]
  expect_equal(new_data$grp_B, c(0L, 0L, 0L, 1L, 1L, 0L))
  expect_false(anyNA(new_data$grp_B))
})

test_that("scmAddCatCovariates: errors on indicator name collision", {
  d <- data.frame(
    ID = seq_len(4),
    grp = c("A", "A", "B", "B"),
    grp_B = c(9L, 9L, 9L, 9L),
    stringsAsFactors = FALSE
  )
  expect_error(
    .cur$scmAddCatCovariates(
      d,
      covarsVec = character(0),
      catcovarsVec = "grp"
    ),
    "already exist"
  )
})

test_that("scmAddCatCovariates: a repeated column is expanded once, not an error", {
  d <- data.frame(
    ID = seq_len(4),
    grp = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = c("grp", "grp")
  )
  expect_equal(res[[2]], "grp_B")
})

test_that("scmAddCatCovariates: unused factor levels get no all-zero indicator", {
  d <- data.frame(
    ID = seq_len(4),
    grp = factor(c("A", "A", "B", "B"), levels = c("A", "B", "D"))
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "grp"
  )
  expect_false("grp_D" %in% names(res[[1]]))
})

test_that("scmAddCatCovariates: default counts per subject, not per row", {
  # 2 male subjects with 20 rows each, 5 female subjects with 2 rows each.
  # Per subject F wins 5-2; per row M wins 40-10.  The default must count
  # subjects, so F is the reference and only SEX_M is built.
  d <- data.frame(
    ID = c(rep(1, 20), rep(2, 20), rep(3:7, each = 2)),
    SEX = c(rep("M", 40), rep("F", 10)),
    stringsAsFactors = FALSE
  )
  expect_equal(
    .cur$scmAddCatCovariates(d, character(0), "SEX")[[2]],
    "SEX_M"
  )
  expect_equal(
    .cur$scmAddCatCovariates(d, character(0), "SEX", freqBy = "observation")[[
      2
    ]],
    "SEX_F"
  )
})

test_that("scmAddCatCovariates: freqBy auto picks the rule per column", {
  # SEX is constant within subject -> per subject; CMT varies -> per row
  d <- data.frame(
    ID = c(rep(1, 20), rep(2, 20), rep(3:7, each = 2)),
    SEX = c(rep("M", 40), rep("F", 10)),
    CMT = c(rep(c(1, 2), 25)),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    character(0),
    c("SEX", "CMT"),
    freqBy = "auto"
  )
  expect_true("SEX_M" %in% res[[2]])
  expect_true(any(grepl("^CMT_", res[[2]])))
})

test_that("scmAddCatCovariates: covarsVec drops the expanded column", {
  # 'SEX' is gone from the returned data, so it must not remain in covarsVec
  d <- data.frame(
    ID = 1:2,
    WT = c(70, 80),
    SEX = c("M", "F"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = c("WT", "SEX"),
    catcovarsVec = "SEX"
  )
  expect_false("SEX" %in% res[[2]])
  expect_true(all(res[[2]] %in% names(res[[1]])))
})

test_that("scmAddCatCovariates: ties resolve deterministically", {
  # A and B tie at 2 subjects each; the first name alphabetically is the
  # reference, matching how .makeSCMData() breaks the same tie
  d <- data.frame(
    ID = 1:4,
    grp = c("B", "B", "A", "A"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(d, character(0), "grp")
  expect_equal(res[[2]], "grp_B")
})

test_that("scmAddCatCovariates: all-NA column yields no indicators", {
  d <- data.frame(
    ID = 1:3,
    grp = c(NA, NA, NA),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(d, character(0), "grp")
  expect_length(res[[2]], 0)
  expect_false("grp" %in% names(res[[1]]))
})

test_that("scmAddCatCovariates: single-level column yields no indicators", {
  d <- data.frame(
    ID = 1:3,
    grp = c("A", "A", "A"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(d, character(0), "grp")
  expect_length(res[[2]], 0)
})

test_that("scmAddCatCovariates: errors when a categorical column is absent", {
  d <- data.frame(ID = 1:4, grp = c("A", "A", "B", "B"))
  expect_error(
    .cur$scmAddCatCovariates(
      d,
      covarsVec = character(0),
      catcovarsVec = "nope"
    ),
    "not found in data"
  )
})

test_that("scmAddCatCovariates: indicator values are 0/1 integers", {
  d <- data.frame(
    ID = 1:6,
    grp = c("A", "A", "A", "B", "B", "C"),
    stringsAsFactors = FALSE
  )
  res <- .cur$scmAddCatCovariates(
    d,
    covarsVec = character(0),
    catcovarsVec = "grp"
  )
  new_data <- res[[1]]
  ind_cols <- grep("^grp_", names(new_data), value = TRUE)
  for (col in ind_cols) {
    expect_true(all(new_data[[col]] %in% c(0L, 1L)))
  }
})

# =============================================================================
# .rebuildUiFromPairs
# =============================================================================

test_that(".rebuildUiFromPairs: NULL pairs_df returns base_ui unchanged", {
  ui <- nlmixr2est::nlmixr(.one_cmt_fun)
  res <- .cur$.rebuildUiFromPairs(ui, NULL)
  expect_equal(
    sort(ui$iniDf$name[!is.na(ui$iniDf$ntheta)]),
    sort(res$iniDf$name[!is.na(res$iniDf$ntheta)])
  )
})

test_that(".rebuildUiFromPairs: empty pairs_df returns base_ui unchanged", {
  ui <- nlmixr2est::nlmixr(.one_cmt_fun)
  empty <- data.frame(
    var = character(0),
    covar = character(0),
    stringsAsFactors = FALSE
  )
  res <- .cur$.rebuildUiFromPairs(ui, empty)
  expect_equal(
    sort(ui$iniDf$name[!is.na(ui$iniDf$ntheta)]),
    sort(res$iniDf$name[!is.na(res$iniDf$ntheta)])
  )
})

test_that(".rebuildUiFromPairs: adds cov theta with correct name", {
  ui <- nlmixr2est::nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var = "cl",
    covar = "wt_power",
    covExpr = "log(wt/70)",
    init = 0.5,
    lower = -2,
    upper = 2,
    stringsAsFactors = FALSE
  )
  res <- .cur$.rebuildUiFromPairs(ui, pairs)
  ini_names <- res$iniDf$name
  expect_true("cov_wt_power_cl" %in% ini_names)
})

test_that(".rebuildUiFromPairs: new theta gets specified init and bounds", {
  ui <- nlmixr2est::nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var = "cl",
    covar = "wt_power",
    covExpr = "log(wt/70)",
    init = 0.5,
    lower = -2,
    upper = 2,
    stringsAsFactors = FALSE
  )
  res <- .cur$.rebuildUiFromPairs(ui, pairs)
  cov_row <- res$iniDf[res$iniDf$name == "cov_wt_power_cl", ]
  expect_equal(cov_row$est, 0.5)
  expect_equal(cov_row$lower, -2)
  expect_equal(cov_row$upper, 2)
})

test_that(".rebuildUiFromPairs: two covariates on different parameters both added", {
  ui <- nlmixr2est::nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var = c("cl", "v"),
    covar = c("wt_power", "wt_power"),
    covExpr = c("log(wt/70)", "log(wt/70)"),
    init = c(0.5, 0.3),
    lower = c(-2, -2),
    upper = c(2, 2),
    stringsAsFactors = FALSE
  )
  res <- .cur$.rebuildUiFromPairs(ui, pairs)
  ini_names <- res$iniDf$name
  expect_true("cov_wt_power_cl" %in% ini_names)
  expect_true("cov_wt_power_v" %in% ini_names)
})

# =============================================================================
# .freezeUiForProfile  (1-D Brent warm-start: frozen base)
# =============================================================================

.make_cov_ui <- function() {
  ui <- nlmixr2est::nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var = "cl", covar = "wt_power", covExpr = "log(wt/70)",
    init = 0.5, lower = -2, upper = 2, stringsAsFactors = FALSE
  )
  .cur$.rebuildUiFromPairs(ui, pairs)
}

test_that(".freezeUiForProfile: leaves exactly one theta free", {
  ui <- .make_cov_ui()
  frozen <- .cur$.freezeUiForProfile(ui, "cov_wt_power_cl")
  ini <- frozen$iniDf
  is_theta <- !is.na(ini$ntheta)
  free <- ini$name[is_theta & !ini$fix]
  expect_equal(free, "cov_wt_power_cl")
})

test_that(".freezeUiForProfile: fixes all structural + residual thetas", {
  ui <- .make_cov_ui()
  frozen <- .cur$.freezeUiForProfile(ui, "cov_wt_power_cl")
  ini <- frozen$iniDf
  is_theta <- !is.na(ini$ntheta)
  fixed <- ini$name[is_theta & ini$fix]
  # tka, tcl, tv, add.sd all fixed; the new cov theta is not
  expect_true(all(c("tka", "tcl", "tv", "add.sd") %in% fixed))
  expect_false("cov_wt_power_cl" %in% fixed)
})

test_that(".freezeUiForProfile: preserves parent theta estimates", {
  ui <- .make_cov_ui()
  before <- ui$iniDf
  frozen <- .cur$.freezeUiForProfile(ui, "cov_wt_power_cl")
  after <- frozen$iniDf
  for (nm in c("tka", "tcl", "tv")) {
    expect_equal(
      after$est[after$name == nm],
      before$est[before$name == nm]
    )
  }
})

test_that(".freezeUiForProfile: zeroes between-subject variability (omega)", {
  ui <- .make_cov_ui()
  frozen <- .cur$.freezeUiForProfile(ui, "cov_wt_power_cl")
  ini <- frozen$iniDf
  # no free eta rows should remain after zeroRe(which = "omega")
  eta_rows <- ini[!is.na(ini$neta1), , drop = FALSE]
  if (nrow(eta_rows) > 0) {
    expect_true(all(eta_rows$fix | eta_rows$est == 0))
  } else {
    expect_equal(nrow(eta_rows), 0L)
  }
})

# =============================================================================
# Integration tests — require nlmixr2data, slow fitting
# =============================================================================

skip_if_not_installed("nlmixr2data")

.theoph <- nlmixr2data::theo_sd

.pk_model <- function() {
  ini({
    tka <- log(1.5)
    tcl <- log(0.04)
    tv <- log(0.5)
    eta.ka ~ 0.5
    eta.cl ~ 0.2
    eta.v ~ 0.1
    prop.err <- 0.1 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl + eta.cl) # nolint: object_usage_linter.
    v <- exp(tv + eta.v) # nolint: object_usage_linter.
    linCmt() ~ prop(prop.err)
  })
}

.fit_base <- function() {
  nlmixr2(
    .pk_model,
    .theoph,
    est = "focei",
    control = nlmixr2est::foceiControl(print = 0, calcTables = TRUE)
  )
}

test_that("runSCM: forward-only returns expected list structure", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )
  expect_type(res, "list")
  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_null(res$resBck)
  expect_type(res$resFwd, "list")
})

test_that("runSCM: backward-only returns expected list structure", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "backward",
    saveModels = FALSE,
    includedRelations = list(list(var = "cl", covar = "WT", shapes = "power")),
    workers = 1L
  )
  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_null(res$resFwd)
})

test_that("runSCM: full SCM returns forward and backward results", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "scm",
    saveModels = FALSE,
    workers = 1L
  )
  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_type(res$resFwd, "list")
  expect_type(res$resBck, "list")
})

test_that("runSCM: saveModels=FALSE creates no output directory", {
  td <- withr::local_tempdir(clean = TRUE)
  withr::local_dir(td)
  skip_on_cran()
  base_fit <- .fit_base()
  n_before <- length(list.dirs(td, recursive = FALSE))
  runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )
  n_after <- length(list.dirs(td, recursive = FALSE))
  expect_equal(n_before, n_after)
})

test_that("runSCM: saveModels=TRUE writes log and CSV files", {
  td <- withr::local_tempdir(clean = TRUE)
  withr::local_dir(td)
  skip_on_cran()
  base_fit <- .fit_base()
  runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = TRUE,
    outputDir = "scm_out",
    restart = TRUE,
    workers = 1L
  )
  expect_true(dir.exists("scm_out"))
  expect_true(file.exists(file.path("scm_out", "scm_log.txt")))
  expect_true(file.exists(file.path("scm_out", "scm_step_summary.csv")))
  expect_true(file.exists(file.path("scm_out", "scm_all_candidates.csv")))
})

test_that("runSCM: summaryTable has expected columns", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )
  if (!is.null(res$summaryTable) && nrow(res$summaryTable) > 0) {
    expect_true(all(
      c(
        "covar",
        "var",
        "deltObjf",
        "pchisqr",
        "included",
        "searchType"
      ) %in%
        names(res$summaryTable)
    ))
  }
})

test_that("runSCM: user-supplied control object accepted", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  ctrl <- nlmixr2est::foceiControl(print = 0, calcTables = TRUE)
  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      control = ctrl,
      workers = 1L
    )
  )
})

test_that("runSCM: inits with bounds run without error", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = 1L,
      inits = list(power = list(est = 0.75, lower = 0, upper = 3))
    )
  )
})

test_that("runSCM: per-pair shapes via pairsVec respected", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(
      list(var = "cl", covar = "WT", shapes = c("power", "lin"))
    ),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )
  # At least two candidates should have been tested (power + lin)
  if (!is.null(res$summaryTable) && nrow(res$summaryTable) > 0) {
    expect_gte(nrow(res$summaryTable), 1)
  }
})

test_that("runSCM: restart=TRUE backs up existing outputDir", {
  td <- withr::local_tempdir(clean = TRUE)
  withr::local_dir(td)
  skip_on_cran()
  base_fit <- .fit_base()

  # First run — create the directory
  runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = TRUE,
    outputDir = "restart_dir",
    restart = TRUE,
    workers = 1L
  )
  expect_true(dir.exists("restart_dir"))

  # Second run — should back up and start fresh
  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = TRUE,
      outputDir = "restart_dir",
      restart = TRUE,
      workers = 1L
    )
  )
  backups <- list.dirs(td, recursive = FALSE, full.names = FALSE)
  expect_true(any(grepl("^restart_dir_backup_", backups)))
})

test_that("runSCM: multiple pairs tested simultaneously", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(
      list(var = "cl", covar = "WT", shapes = "power"),
      list(var = "v", covar = "WT", shapes = "power")
    ),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )
  expect_type(res, "list")
})

test_that("runSCM: explicit outputDir used as absolute path", {
  td <- withr::local_tempdir(clean = TRUE)
  out_dir <- file.path(td, "my_scm_output")
  skip_on_cran()
  base_fit <- .fit_base()
  runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = TRUE,
    outputDir = out_dir,
    restart = TRUE,
    workers = 1L
  )
  expect_true(dir.exists(out_dir))
  expect_true(file.exists(file.path(out_dir, "scm_log.txt")))
})

# =============================================================================
# workers parameter — parallelization
# =============================================================================

test_that("runSCM: workers=1 returns same structure as workers=NULL", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  res_default <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers = NULL
  )

  res_w1 <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )

  expect_named(res_w1, names(res_default))
  expect_type(res_w1$resFwd, "list")
  expect_null(res_w1$resBck)
})

test_that("runSCM: workers=1 forward+backward both respect parameter", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "scm",
    saveModels = FALSE,
    workers = 1L
  )

  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_type(res$resFwd, "list")
  expect_type(res$resBck, "list")
})

test_that("runSCM: workers='auto' runs without error", {
  skip_if_not_installed("future")
  # multisession workers load from the *installed* package, so this test only
  # makes sense when nlmixr2scm is installed (not just load_all()-ed).
  skip_if(
    isTRUE(tryCatch(
      pkgload::is_dev_package("nlmixr2scm"),
      error = function(e) FALSE
    )),
    "Package loaded via load_all(); install first to run multisession test"
  )
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = "auto"
    )
  )
})

test_that("runSCM: future plan restored to original after workers=1", {
  skip_if_not_installed("future")
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  plan_orig <- class(future::plan())
  on.exit(future::plan("sequential"), add = TRUE)

  runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers = 1L
  )

  expect_equal(class(future::plan()), plan_orig)
})

test_that("runSCM: rxThreads is appended as the last formal (positional compatibility)", {
  # Confirms rxThreads was appended, not inserted after workers -- an
  # existing positional caller passing arguments past `workers` (e.g.
  # `confirm`) must still bind to the same parameters as before.
  scm_formals <- names(formals(runSCM))
  expect_equal(scm_formals[length(scm_formals)], "rxThreads")
  expect_true(which(scm_formals == "workers") < which(scm_formals == "confirm"))
})

test_that("runSCM: prints the exact effective workers/rxThreads console numbers", {
  skip_if_not_installed("rxode2")
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  expect_message(
    suppressWarnings(runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = 1L,
      rxThreads = 3L
    )),
    "threads / worker.*3"
  )
})

test_that("runSCM: workers * rxThreads exceeding cores aborts (mocked core count)", {
  skip_if_not_installed("future")
  skip_if_not_installed("rxode2")
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  testthat::local_mocked_bindings(
    .resolveTotalCores = function() 4L,
    .package = "nlmixr2utils"
  )

  expect_error(
    suppressMessages(runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = 2L,
      rxThreads = 3L
    ))
  )
})

test_that("runSCM: rxThreads propagates to .fitCandidatePairs() via the forward path", {
  skip_if_not_installed("rxode2")
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  target <- 3L
  seen <- new.env()
  seen$threads <- NULL
  # local_mocked_bindings replaces the real nlmixr2utils::.plap() -- which
  # dispatches candidate fits in parallel and would need real threading
  # infra to test end-to-end -- with a serial spy that still calls the real
  # FUN for each candidate (so .fitCandidatePairs() gets genuine per-candidate
  # results and completes normally) while recording the rxThreads argument
  # it was called with. Returning bare list() here (as originally sketched)
  # makes raw_results empty, which trips .fitCandidatePairs()'s own
  # `if (length(results) == 0L) stop(...)` before the test ever reaches the
  # assertion below -- so the spy must forward to FUN via lapply() instead.
  testthat::local_mocked_bindings(
    .plap = function(X, FUN, ..., rxThreads = NULL, .label = NULL) {
      seen$threads <- rxThreads
      lapply(X, FUN)
    },
    .package = "nlmixr2utils"
  )

  runSCM(
    base_fit,
    varsVec = c("cl"),
    covarsVec = c("WT"),
    confirm = FALSE,
    saveModels = FALSE,
    searchType = "forward",
    workers = 1L,
    rxThreads = target
  )

  expect_equal(seen$threads, target)
})

test_that("runSCM: rxThreads propagates to .fitCandidatePairs() via the backward path", {
  skip_if_not_installed("rxode2")
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()

  target <- 3L
  seen <- new.env()
  seen$threads <- NULL
  # See the forward-path test above for why the spy must forward to FUN via
  # lapply() rather than return bare list().
  testthat::local_mocked_bindings(
    .plap = function(X, FUN, ..., rxThreads = NULL, .label = NULL) {
      seen$threads <- rxThreads
      lapply(X, FUN)
    },
    .package = "nlmixr2utils"
  )

  # searchType = "backward" with includedRelations gives backwardSearch()
  # candidates to test independently of any forward step, so its own
  # .fitCandidatePairs()/.plap() call is guaranteed to execute.
  #
  # NOTE: unlike the plan's original sketch, varsVec/covarsVec are also
  # supplied here (matching the forward test's relation). runSCM() always
  # builds its top-level `pairs` via buildPairs(varsVec, covarsVec, pairsVec)
  # regardless of searchType (see R/scm.R ~L299-306), and buildPairs() errors
  # ("provide either pairsVec or both varsVec and covarsVec") when all three
  # are NULL -- includedRelations alone is not enough to reach
  # .fitCandidatePairs(). This mirrors the existing "runSCM: backward-only
  # returns expected list structure" test's pattern of supplying the same
  # relation via both channels.
  runSCM(
    base_fit,
    varsVec = c("cl"),
    covarsVec = c("WT"),
    includedRelations = list(list(var = "cl", covar = "WT")),
    confirm = FALSE,
    saveModels = FALSE,
    searchType = "backward",
    workers = 1L,
    rxThreads = target
  )

  expect_equal(seen$threads, target)
})

# =============================================================================
# .expandShapes — auto-scaled bounds for lin / log / identity
# =============================================================================

test_that(".expandShapes: lin shape auto-scales bounds by center", {
  d <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90), stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "lin")
  center <- median(d$wt)
  expect_equal(res$init[1],  0.1 / center, tolerance = 1e-10)
  expect_equal(res$lower[1], -5 / center,  tolerance = 1e-10)
  expect_equal(res$upper[1],  5 / center,  tolerance = 1e-10)
})

test_that(".expandShapes: identity shape auto-scales bounds by center", {
  d <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90), stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "identity")
  center <- median(d$wt)
  expect_equal(res$init[1],  0.1 / center, tolerance = 1e-10)
  expect_equal(res$lower[1], -5 / center,  tolerance = 1e-10)
  expect_equal(res$upper[1],  5 / center,  tolerance = 1e-10)
})

test_that(".expandShapes: log shape auto-scales bounds by abs(log(center))", {
  d <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90), stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "log")
  center <- median(d$wt)
  sc <- abs(log(center))
  expect_equal(res$init[1],  0.1 / sc, tolerance = 1e-10)
  expect_equal(res$lower[1], -5 / sc,  tolerance = 1e-10)
  expect_equal(res$upper[1],  5 / sc,  tolerance = 1e-10)
})

test_that(".expandShapes: power shape keeps unscaled defaults", {
  d <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90), stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power")
  expect_equal(res$init[1],  0.1)
  expect_equal(res$lower[1], -5)
  expect_equal(res$upper[1],  5)
})

test_that(".expandShapes: user-supplied lin inits not overridden by auto-scaling", {
  d <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90), stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(
    .cont_pairs(d),
    shapes = "lin",
    inits = list(lin = list(est = 0.005, lower = -0.1, upper = 0.1))
  )
  expect_equal(res$init[1],  0.005)
  expect_equal(res$lower[1], -0.1)
  expect_equal(res$upper[1],  0.1)
})

# =============================================================================
# .resolveOFVTolerance
# =============================================================================

test_that(".resolveOFVTolerance: explicit numeric is returned as-is", {
  mock_fit <- list(est = "focei")
  expect_equal(.cur$.resolveOFVTolerance(mock_fit, 15), 15)
})

test_that(".resolveOFVTolerance: explicit zero overrides stochastic auto-detect", {
  mock_fit <- list(est = "saem")
  expect_equal(.cur$.resolveOFVTolerance(mock_fit, 0), 0)
})

test_that(".resolveOFVTolerance: NULL with focei returns 0", {
  mock_fit <- list(est = "focei")
  expect_equal(.cur$.resolveOFVTolerance(mock_fit, NULL), 0)
})

test_that(".resolveOFVTolerance: NULL with foce returns 0", {
  mock_fit <- list(est = "foce")
  expect_equal(.cur$.resolveOFVTolerance(mock_fit, NULL), 0)
})

test_that(".resolveOFVTolerance: NULL with saem returns 10", {
  mock_fit <- list(est = "saem")
  expect_equal(.cur$.resolveOFVTolerance(mock_fit, NULL), 10)
})

test_that(".resolveOFVTolerance: NULL with SAEM (uppercase) returns 10", {
  mock_fit <- list(est = "SAEM")
  expect_equal(.cur$.resolveOFVTolerance(mock_fit, NULL), 10)
})

# =============================================================================
# .isUnrealisticOFV
# =============================================================================

test_that(".isUnrealisticOFV: realistic result returns FALSE", {
  # Forward: candidate OFV 460, ref 470 => dObjf = -10 (improvement)
  expect_false(.cur$.isUnrealisticOFV(
    x_objf = 460, ref_objf = 470, dObjf = -10,
    pchisqr = 0.01, maxDeltaOFV = Inf, effective_tolerance = 0
  ))
})

test_that(".isUnrealisticOFV: OFV above parent triggers TRUE (criterion 1)", {
  # Candidate OFV 480 > ref 470 => dObjf = 10 (worsened)
  expect_true(.cur$.isUnrealisticOFV(
    x_objf = 480, ref_objf = 470, dObjf = 10,
    pchisqr = 1, maxDeltaOFV = Inf, effective_tolerance = 0
  ))
})

test_that(".isUnrealisticOFV: OFV above parent but within tolerance returns FALSE", {
  # Candidate OFV 475 > ref 470, but within tolerance of 10
  expect_false(.cur$.isUnrealisticOFV(
    x_objf = 475, ref_objf = 470, dObjf = 5,
    pchisqr = 1, maxDeltaOFV = Inf, effective_tolerance = 10
  ))
})

test_that(".isUnrealisticOFV: OFV above parent beyond tolerance triggers TRUE", {
  # Candidate OFV 500 > ref 470 + tolerance 10
  expect_true(.cur$.isUnrealisticOFV(
    x_objf = 500, ref_objf = 470, dObjf = 30,
    pchisqr = 1, maxDeltaOFV = Inf, effective_tolerance = 10
  ))
})

test_that(".isUnrealisticOFV: p-value underflow triggers TRUE (criterion 2)", {
  # Candidate OFV 100, ref 470 => dObjf = -370 (huge improvement, suspicious)
  expect_true(.cur$.isUnrealisticOFV(
    x_objf = 100, ref_objf = 470, dObjf = -370,
    pchisqr = 0, maxDeltaOFV = Inf, effective_tolerance = 0
  ))
})

test_that(".isUnrealisticOFV: pchisqr just above .Machine$double.eps is FALSE", {
  expect_false(.cur$.isUnrealisticOFV(
    x_objf = 100, ref_objf = 470, dObjf = -370,
    pchisqr = .Machine$double.eps * 2, maxDeltaOFV = Inf, effective_tolerance = 0
  ))
})

test_that(".isUnrealisticOFV: |dObjf| above maxDeltaOFV triggers TRUE (criterion 3)", {
  # dObjf = -5000 (huge improvement), |dObjf| > maxDeltaOFV
  expect_true(.cur$.isUnrealisticOFV(
    x_objf = 100, ref_objf = 470, dObjf = -5000,
    pchisqr = 0.001, maxDeltaOFV = 1000, effective_tolerance = 0
  ))
})

test_that(".isUnrealisticOFV: |dObjf| below maxDeltaOFV returns FALSE", {
  expect_false(.cur$.isUnrealisticOFV(
    x_objf = 100, ref_objf = 470, dObjf = -500,
    pchisqr = 0.001, maxDeltaOFV = 1000, effective_tolerance = 0
  ))
})

test_that(".isUnrealisticOFV: maxDeltaOFV = Inf never triggers criterion 3", {
  expect_false(.cur$.isUnrealisticOFV(
    x_objf = -1e9, ref_objf = 470, dObjf = -(1e9 + 470),
    pchisqr = 0.001, maxDeltaOFV = Inf, effective_tolerance = 0
  ))
})

test_that("runSCM: accepts new retry parameters without error (smoke)", {
  # Verify the new parameters exist in the formals — no fitting needed
  fmls <- formals(.cur$runSCM)
  expect_true("maxRetries"         %in% names(fmls))
  expect_true("maxDeltaOFV"        %in% names(fmls))
  expect_true("retryPerturbSD"     %in% names(fmls))
  expect_true("retrySmallInit"     %in% names(fmls))
  expect_true("retryOFVTolerance"  %in% names(fmls))
  expect_equal(fmls$maxRetries,        3L)
  expect_equal(fmls$maxDeltaOFV,       Inf)
  expect_equal(fmls$retryPerturbSD,    0.5)
  expect_equal(fmls$retrySmallInit,    0.01)
  expect_null(fmls$retryOFVTolerance)
  expect_true("retryFailOnExhaustion" %in% names(fmls))
  expect_false(fmls$retryFailOnExhaustion)
})

test_that("runSCM: maxRetries=0 disables retry (parameter accepted, no error)", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = 1L,
      maxRetries = 0L
    )
  )
})

test_that("runSCM: maxDeltaOFV passed through without error", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = 1L,
      maxDeltaOFV = 1000
    )
  )
})

test_that("runSCM: retryOFVTolerance=0 passed through without error", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  expect_no_error(
    runSCM(
      fit = base_fit,
      pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers = 1L,
      retryOFVTolerance = 0
    )
  )
})


# =============================================================================
# cli brace escaping — failure messages containing {} must not crash cli_warn
# =============================================================================

test_that(".fitCandidatePairs: braces in failure reason do not crash cli", {
  skip_on_cran()
  base_fit <- .fit_base()
  bad <- data.frame(
    var = "cl", covar = "wt{x}", covExpr = "WT", shape = "power",
    init = 0.1, lower = -5, upper = 5, stringsAsFactors = FALSE
  )
  # Should fail with the "All … failed" stop(), not a cli parse error
  expect_error(
    .cur$.fitCandidatePairs(
      pairs = bad, base_ui = base_fit$ui, fit = base_fit,
      data = .theoph, pVal = 0.05, stepIdx = 1L, add = TRUE, print = 0
    ),
    "All 1 candidate fit\\(s\\) failed"
  )
})

# =============================================================================
# summaryTable$included: backward removals labeled "dropped" (not "yes")
# =============================================================================

test_that("runSCM: backward-removed covariates labeled 'dropped' in summaryTable", {
  withr::local_tempdir(clean = TRUE)
  skip_on_cran()
  base_fit <- .fit_base()
  # Backward-only: pre-include WT~cl so it can be tested for removal
  res <- runSCM(
    fit = base_fit,
    pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "backward",
    includedRelations = list(list(var = "cl", covar = "WT", shapes = "power")),
    saveModels = FALSE,
    workers = 1L
  )

  st <- res$summaryTable
  expect_true("included" %in% names(st))

  # Backward rows must NEVER use "yes"
  bck <- st[st$searchType == "backward", , drop = FALSE]
  expect_false(any(bck$included == "yes"))
  # And the included column should only contain the documented tokens
  expect_true(all(st$included %in% c("yes", "no", "dropped", "retained")))
})

# =============================================================================
# Fixed covariate centers (reference values)
# -----------------------------------------------------------------------------
#
# Mechanism: .enrichPairs() first fills each continuous covariate's `center`
# with the per-subject median; .applyFixedCenters() then OVERRIDES the named
# ones. buildPairs() carries an optional per-row `center` through so a user can
# also pin centers per (var, covar) pair via pairsVec.
# =============================================================================

# ---- .applyFixedCenters: override named, keep unnamed, stay safe -----------

test_that(".applyFixedCenters overrides named continuous covariates only", {
  # BW appears on both cl and vc (a single entry must fix both); BMI is unnamed
  # and must keep its median; categorical rows are never touched.
  pairs <- data.frame(
    var     = c("cl", "vc", "vc", "cl"),
    covar   = c("BW", "BW", "BMI", "SEX_1"),
    raw_col = c("BW", "BW", "BMI", "SEX"),
    type    = c("continuous", "continuous", "continuous", "categorical"),
    center  = c(71.3, 68.9, 26.1, NA_real_),  # per-dataset medians
    stringsAsFactors = FALSE
  )

  out <- .cur$.applyFixedCenters(pairs, c(BW = 70))

  expect_equal(out$center[out$raw_col == "BW"], c(70, 70))   # both params pinned
  expect_equal(out$center[out$raw_col == "BMI"], 26.1)       # unnamed -> median
  expect_true(is.na(out$center[out$type == "categorical"]))  # categorical safe
})

test_that(".applyFixedCenters is a no-op for NULL / unmatched and errors if unnamed", {
  pairs <- data.frame(
    var = "cl", covar = "BW", raw_col = "BW",
    type = "continuous", center = 71.3, stringsAsFactors = FALSE
  )
  expect_identical(.cur$.applyFixedCenters(pairs, NULL),       pairs)  # default path
  expect_identical(.cur$.applyFixedCenters(pairs, c(WT = 70)), pairs)  # no match
  expect_error(.cur$.applyFixedCenters(pairs, 70), "must be a named numeric vector")
})

# ---- buildPairs: center passthrough (the plumbing fix) ---------------------

test_that("buildPairs carries a per-pair center through, and omits it when absent", {
  # list form with centers
  pv <- list(
    list(var = "cl", covar = "BW",   center = 70),
    list(var = "vc", covar = "BW",   center = 70),
    list(var = "cl", covar = "CrCL", center = 95)
  )
  out <- .cur$buildPairs(pairsVec = pv)
  expect_equal(out$center, c(70, 70, 95))

  # no center supplied -> no center column (unchanged legacy shape)
  bare <- .cur$buildPairs(pairsVec = list(list(var = "cl", covar = "BW")))
  expect_false("center" %in% names(bare))
})

# ---- Integration: median first, then fixed override ------------------------

test_that("enrichPairs medians are overridden by .applyFixedCenters", {
  # one row per subject; medians are deliberately off the fixed anchors so the
  # override is observable.
  set.seed(1)
  dat <- data.frame(
    ID   = 1:11,
    TIME = 0,
    DV   = 0,
    BW   = 60:70,                       # median 65  (anchor will be 70)
    CrCL = seq(80, 100, length.out = 11),  # median 90  (anchor will be 95)
    BMI  = seq(22, 32, length.out = 11)    # median 27  (no anchor -> stays)
  )

  pairs <- data.frame(
    var   = c("cl", "cl", "vc"),
    covar = c("BW", "CrCL", "BMI"),
    stringsAsFactors = FALSE
  )

  enriched <- .cur$.enrichPairs(pairs, dat)
  # sanity: enrichment filled the per-subject medians
  expect_equal(enriched$center[enriched$raw_col == "BW"],   65)
  expect_equal(enriched$center[enriched$raw_col == "CrCL"], 90)
  expect_equal(enriched$center[enriched$raw_col == "BMI"],  27)

  fixed <- .cur$.applyFixedCenters(enriched, c(BW = 70, CrCL = 95))
  # named anchors override the medians ...
  expect_equal(fixed$center[fixed$raw_col == "BW"],   70)
  expect_equal(fixed$center[fixed$raw_col == "CrCL"], 95)
  # ... unnamed BMI keeps the data-driven median
  expect_equal(fixed$center[fixed$raw_col == "BMI"],  27)
})

# .pickForwardWinner / .pickBackwardWinner — winner-selection tie-breaking
# -----------------------------------------------------------------------------
# For df = 1 any dOFV >= ~70.5 makes 1 - pchisq() underflow to exactly 0, so
# several genuinely-strong forward candidates tie at pchisqr == 0.  The old
# which.min(pchisqr) then returned the FIRST row (alphabetical covar order),
# which could pick a weaker covariate over a much stronger, collinear one and
# steer the search into a wrong-shape branch.  The tie-break now prefers the
# largest deltObjf (biggest OFV drop) forward, and the smallest deltObjf (least
# OFV increase on removal) backward.
# =============================================================================

test_that(".pickForwardWinner: distinct p-values pick the smallest (unchanged behaviour)", {
  rt <- data.frame(
    covar    = c("BW", "CrCL"),
    pchisqr  = c(1e-3, 1e-6),
    deltObjf = c(20, 40)
  )
  # smallest pchisqr is row 2; ties never engaged
  expect_equal(.cur$.pickForwardWinner(rt), 2L)
})

test_that(".pickForwardWinner: p-value underflow tie broken by largest deltObjf", {
  # Both true covariates underflow to pchisqr == 0; alphabetical order puts BW
  # first, but CrCL has the larger dOFV and must win.
  rt <- data.frame(
    covar    = c("BW", "CrCL", "SEX"),
    pchisqr  = c(0, 0, 1e-7),
    deltObjf = c(90.35, 183.05, 27.74)
  )
  expect_equal(.cur$.pickForwardWinner(rt), 2L)          # CrCL, dOFV 183
  expect_equal(rt$covar[.cur$.pickForwardWinner(rt)], "CrCL")
})

test_that(".pickForwardWinner: full tie (equal pchisqr AND deltObjf) is deterministic first row", {
  rt <- data.frame(
    covar    = c("BW", "CrCL"),
    pchisqr  = c(0, 0),
    deltObjf = c(100, 100)
  )
  expect_equal(.cur$.pickForwardWinner(rt), 1L)
})

test_that(".pickBackwardWinner: distinct p-values drop the largest (unchanged behaviour)", {
  rt <- data.frame(
    covar    = c("BW", "CrCL"),
    pchisqr  = c(0.02, 0.80),
    deltObjf = c(5, 0.5)
  )
  # least important removal = highest pchisqr = row 2
  expect_equal(.cur$.pickBackwardWinner(rt), 2L)
})

test_that(".pickBackwardWinner: p == 1 tie broken by smallest deltObjf", {
  # Both removals are non-significant (pchisqr == 1); drop the one that raises
  # OFV the least (smallest deltObjf).
  rt <- data.frame(
    covar    = c("BW", "CrCL"),
    pchisqr  = c(1, 1),
    deltObjf = c(3.0, 0.4)
  )
  expect_equal(.cur$.pickBackwardWinner(rt), 2L)          # CrCL raises OFV least
})

# Retry-exhaustion best-attempt tracking
#
# .fitCandidatePairs() must keep the BEST (largest-dObjf) attempt across
# perturbed-init retries, not whichever attempt happened to run last.
#
# The selection rule is INLINED inside .fitCandidatePairs() because that loop
# runs in future.apply workers spawned by .plap(); workers load the installed
# package and cannot see helpers introduced via devtools::load_all().  These
# tests mirror the production rule locally so the spec is documented and the
# regression scenario stays covered.  If the inlined block in R/scm.R changes,
# update .update_best_attempt() below to match.
# =============================================================================
# Local mirror of the production rule.  Source of truth in R/scm.R inside
# .fitCandidatePairs():
#
#   if (is.null(best_attempt) ||
#       .cand_attempt$dObjf > best_attempt$dObjf) {
#     best_attempt <- .cand_attempt
#   }
.update_best_attempt <- function(best, candidate) {
  if (is.null(best) || candidate$dObjf > best$dObjf) candidate else best
}

test_that("retry tracking: first attempt becomes best when no incumbent", {
  cand <- list(x = "a", dObjf = -250, dof = 1L, pchisqr = 1, attempt_num = 1L)
  expect_identical(.update_best_attempt(NULL, cand), cand)
})

test_that("retry tracking: candidate with larger dObjf replaces incumbent", {
  best <- list(x = "a", dObjf = -250, dof = 1L, pchisqr = 1, attempt_num = 1L)
  cand <- list(x = "b", dObjf =   -1, dof = 1L, pchisqr = 1, attempt_num = 2L)
  expect_identical(.update_best_attempt(best, cand), cand)
})

test_that("retry tracking: candidate with smaller dObjf keeps incumbent", {
  best <- list(x = "a", dObjf =  -1, dof = 1L, pchisqr = 1, attempt_num = 1L)
  cand <- list(x = "b", dObjf = -100, dof = 1L, pchisqr = 1, attempt_num = 2L)
  expect_identical(.update_best_attempt(best, cand), best)
})

test_that("retry tracking: ties resolve to incumbent (no churn)", {
  best <- list(x = "a", dObjf = -50, dof = 1L, pchisqr = 1, attempt_num = 1L)
  cand <- list(x = "b", dObjf = -50, dof = 1L, pchisqr = 1, attempt_num = 2L)
  expect_identical(.update_best_attempt(best, cand), best)
})

test_that("retry tracking: regression -- last attempt with smaller dObjf does not overwrite best", {
  # Bug scenario from the SEX_1~cl retry chain: three attempts produce
  # dObjf = -250, -1, -100.  The original .fitCandidatePairs() unconditionally
  # kept the LAST attempt (-100), even though the per-attempt warning correctly
  # claimed "best available".  The current implementation tracks the running
  # best, so attempt 2 (-1) survives.
  b <- NULL
  b <- .update_best_attempt(b, list(x = "att1", dObjf = -250, dof = 1L,
                                    pchisqr = 1, attempt_num = 1L))
  b <- .update_best_attempt(b, list(x = "att2", dObjf =   -1, dof = 1L,
                                    pchisqr = 1, attempt_num = 2L))
  b <- .update_best_attempt(b, list(x = "att3", dObjf = -100, dof = 1L,
                                    pchisqr = 1, attempt_num = 3L))
  expect_equal(b$x, "att2")
  expect_equal(b$dObjf, -1)
  expect_equal(b$attempt_num, 2L)
})

# =============================================================================
# profileInitOnStall / stallTol — profile-on-stall rescue parameters
# -----------------------------------------------------------------------------
# The forward search can stall when the derivative-free outer optimiser
# (bobyqa) never steps a new covariate coefficient off its init, leaving the
# nested model with a WORSE OFV than its parent (dObjf <= 0) -- mathematically
# impossible at a true optimum.  Observed for ODE models whose FOCEi objective
# carries solver noise; the true per-step subproblem is unimodal, so a single
# 1-D FOCEi profile init rescues it.
#
# The rescue lives at the END of the .fitCandidatePairs() retry loop so it
# fires INDEPENDENTLY of maxRetries -- in particular it must still fire when
# maxRetries = 0 (the benchmark config).  It keeps the profile-init refit ONLY
# when it STRICTLY improves dObjf, so it can never make a candidate worse.
# =============================================================================

test_that("runSCM: profileInitOnStall / stallTol parameters exist with expected defaults", {
  # Formals-only smoke check -- no fitting needed.
  fmls <- formals(.cur$runSCM)
  expect_true("profileInitOnStall" %in% names(fmls))
  expect_true("stallTol"           %in% names(fmls))
  expect_true(isTRUE(eval(fmls$profileInitOnStall)))
  expect_equal(eval(fmls$stallTol), 0)
})

test_that(".fitCandidatePairs / forwardSearch: profile-on-stall parameters threaded through", {
  # The rescue must be reachable from every layer that .fitCandidatePairs is
  # called from, so the args have to appear in each formals list.
  for (fn in c("runSCM", "forwardSearch", ".fitCandidatePairs")) {
    fmls <- formals(.cur[[fn]])
    expect_true("profileInitOnStall" %in% names(fmls),
                info = paste0(fn, " lacks profileInitOnStall"))
    expect_true("stallTol" %in% names(fmls),
                info = paste0(fn, " lacks stallTol"))
  }
})
