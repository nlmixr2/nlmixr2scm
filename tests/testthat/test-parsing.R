library(nlmixr2utils)

# ==== getThetaName

.cur <- loadNamespace("nlmixr2scm")

test_that("get the population parameter from variable name", {
  ## Compartment specifications to test
  # simple one compartment with ka,cl,v
  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45
      label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui1 <- nlmixr2est::nlmixr(one.cmt)
  ui <- ui1
  varName <- "ka"

  funstring1 <- .cur$.getThetaName(ui, varName)
  funstring2 <- "tka"

  expect_equal(funstring1, funstring2)
})


test_that("get the population parameter from variable name", {
  two.compartment <- function() {
    ini({
      tcl <- log(53.4) # Log Cl
      tv1 <- log(73.6)
      tv2 <- log(320) # Log V
      tQ <- log(191)
      eta.cl ~ 0.43^2
      eta.v1 ~ 0.48^2
      eta.v2 ~ 0.49^2
      eta.Q ~ 0.36^2
      prop.sd <- 0.44^2
    })
    # and a model block with the error specification and model specification
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      v2 <- exp(tv2 + eta.v2)
      Q <- exp(tQ + eta.Q)
      linCmt() ~ prop(prop.sd)
    })
  }

  ui2 <- nlmixr2est::nlmixr(two.compartment)
  ui <- ui2
  varName <- "ka"
  expect_error(.cur$.getThetaName(ui, varName))
})

test_that("get the population parameter from variable name", {
  # tainted model
  tainted <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45 # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    # and a model block with the error specification and model specification
    model({
      ka <- exp(tka + eta.ka)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - exp(tcl + eta.cl) / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  tui <- nlmixr2est::nlmixr(tainted)
  ui <- tui
  varName <- "cl"
  expect_error(.cur$.getThetaName(ui, varName))
})


# ==== addCovariate

test_that("Add covariate to the ui", {
  ## Compartment specifications to test
  # simple one compartment with ka,cl,v
  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45
      label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui1 <- nlmixr2est::nlmixr(one.cmt)
  ui <- ui1
  varName <- "ka"
  covariate <- "WT"

  res <- .cur$scmAddOrRemoveCovariate(ui, varName, covariate, add = TRUE)
  res_strs <- vapply(
    res,
    function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "),
    character(1L)
  )

  # New multi-line format: intermediate variable, aggregate, then parameter line
  expect_true(any(grepl("WT_ka <- WT \\* cov_WT_ka", res_strs)))
  expect_true(any(grepl("cov_ka <- WT_ka", res_strs)))
  expect_true(any(grepl(
    "ka <- exp\\(tka \\+ eta\\.ka \\+ cov_ka\\)",
    res_strs
  )))
})


test_that("Add covariate to the ui", {
  two.compartment <- function() {
    ini({
      tcl <- log(53.4) # Log Cl
      tv1 <- log(73.6)
      tv2 <- log(320) # Log V
      tQ <- log(191)
      eta.cl ~ 0.43^2
      eta.v1 ~ 0.48^2
      eta.v2 ~ 0.49^2
      eta.Q ~ 0.36^2
      prop.sd <- 0.44^2
    })
    # and a model block with the error specification and model specification
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      v2 <- exp(tv2 + eta.v2)
      Q <- exp(tQ + eta.Q)
      linCmt() ~ prop(prop.sd)
    })
  }

  ui2 <- nlmixr2est::nlmixr(two.compartment)
  ui <- ui2
  varName <- "cl"
  covariate <- "WT"

  res2 <- .cur$scmAddOrRemoveCovariate(ui, varName, covariate, add = TRUE)
  res_strs2 <- vapply(
    res2,
    function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "),
    character(1L)
  )

  # New multi-line format: intermediate variable, aggregate, then parameter line
  expect_true(any(grepl("WT_cl <- WT \\* cov_WT_cl", res_strs2)))
  expect_true(any(grepl("cov_cl <- WT_cl", res_strs2)))
  expect_true(any(grepl(
    "cl <- exp\\(tcl \\+ eta\\.cl \\+ cov_cl\\)",
    res_strs2
  )))
})

test_that("Add covariate to the ui", {
  two.compartment <- function() {
    ini({
      tcl <- log(53.4) # Log Cl
      tv1 <- log(73.6)
      tv2 <- log(320) # Log V
      tQ <- log(191)
      eta.cl ~ 0.43^2
      eta.v1 ~ 0.48^2
      eta.v2 ~ 0.49^2
      eta.Q ~ 0.36^2
      prop.sd <- 0.44^2
    })
    # and a model block with the error specification and model specification
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      v2 <- exp(tv2 + eta.v2)
      Q <- exp(tQ + eta.Q)
      linCmt() ~ prop(prop.sd)
    })
  }

  ui2 <- nlmixr2est::nlmixr(two.compartment)
  ui <- ui2
  varName <- "v"
  covariate <- "WT"

  expect_error(.cur$scmAddOrRemoveCovariate(ui, varName, covariate, add = TRUE))
})


test_that("Add covariate to the ui", {
  # tainted model
  tainted <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45 # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    # and a model block with the error specification and model specification
    model({
      ka <- exp(tka + eta.ka)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - exp(tcl + eta.cl) / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  tui <- nlmixr2est::nlmixr(tainted)
  ui <- tui
  varName <- "v1"
  covariate <- "WT"
  expect_error(.cur$scmAddOrRemoveCovariate(ui, varName, covariate, add = TRUE))
})


# ==== idColumn

test_that("Extract column corresponding to  Individual", {
  funstring1 <- .idColumn(Theoph)
  funstring2 <- "ID"

  expect_equal(funstring1, funstring2)
})

test_that(".idColumn() finds id regardless of column position", {
  # Regression test for Bug A:
  #   uidCol <- colNames[which("id" %in% colNamesLower)]
  # `"id" %in% colNamesLower` is a single TRUE/FALSE, so `which()` returns 1
  # (or integer(0)) instead of the *index* of the matching column. The function
  # therefore returned the first column whenever an id-like column existed but
  # was not in position 1.

  # id IS the first column -- passes even with the buggy code (by accident)
  d_first <- data.frame(ID = 1:3, TIME = 0:2, WT = c(70, 80, 90))
  expect_equal(.cur$.idColumn(d_first), "ID")

  # id is NOT the first column -- buggy code returns "TIME"
  d_second <- data.frame(TIME = 0:2, ID = 1:3, WT = c(70, 80, 90))
  expect_equal(.cur$.idColumn(d_second), "ID")

  # id is in the middle, lowercase -- buggy code returns "TIME"
  d_middle <- data.frame(TIME = 0:2, AMT = c(100, 0, 0),
                         id = 1:3, WT = c(70, 80, 90))
  expect_equal(.cur$.idColumn(d_middle), "id")

  # mixed case is preserved
  d_mixed <- data.frame(WT = c(70, 80, 90), Id = 1:3)
  expect_equal(.cur$.idColumn(d_mixed), "Id")
})

test_that(".idColumn() falls back to 'ID' when no id-like column exists", {
  d <- data.frame(TIME = 0:2, WT = c(70, 80, 90))
  expect_equal(.cur$.idColumn(d), "ID")
})


# ==== Build ui from the covariate

test_that("Build ui from the covariate", {
  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45
      label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui1 <- nlmixr2est::nlmixr(one.cmt)
  ui <- ui1
  varName <- "ka"
  covariate <- "WT"

  funstring1 <- intersect(
    (.builduiCovariate(ui, varName, covariate, add = TRUE))$iniDf$name,
    "cov_WT_ka"
  )
  funstring2 <- "cov_WT_ka"
  funstring3 <- (.builduiCovariate(
    ui,
    varName,
    covariate,
    add = TRUE
  ))$covariates
  funstring4 <- "WT"
  expect_equal(funstring1, funstring2)
  expect_equal(funstring3, funstring4)
})

test_that("Build ui from the covariate sequentially on a tainted ui", {
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  ui <- nlmixr2est::nlmixr(one.cmt)
  ui <- .builduiCovariate(ui, "ka", "WT", add = TRUE)
  ui <- .builduiCovariate(ui, "cl", "WT", add = TRUE)
  ui_strs <- vapply(
    ui$lstExpr,
    function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "),
    character(1L)
  )

  expect_true(all(c("cov_WT_ka", "cov_WT_cl") %in% ui$iniDf$name))
  expect_true(all(
    c("cov_WT_ka", "cov_WT_cl") %in%
      ui$muRefCovariateDataFrame$covariateParameter
  ))
  expect_true(any(grepl("cov_ka <- WT_ka", ui_strs, fixed = TRUE)))
  expect_true(any(grepl("cov_cl <- WT_cl", ui_strs, fixed = TRUE)))
})

test_that("Build ui from the covariate", {
  two.compartment <- function() {
    ini({
      tcl <- log(53.4) # Log Cl
      tv1 <- log(73.6)
      tv2 <- log(320) # Log V
      tQ <- log(191)
      eta.cl ~ 0.43^2
      eta.v1 ~ 0.48^2
      eta.v2 ~ 0.49^2
      eta.Q ~ 0.36^2
      prop.sd <- 0.44^2
    })
    # and a model block with the error specification and model specification
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      v2 <- exp(tv2 + eta.v2)
      Q <- exp(tQ + eta.Q)
      linCmt() ~ prop(prop.sd)
    })
  }

  ui2 <- nlmixr2est::nlmixr(two.compartment)
  ui <- ui2
  varName <- "ka"
  expect_error(.builduiCovariate(ui, varName, covariate, add = TRUE))
})


# ==== Build covInfo list from varsVec and covarsVec

test_that("Build covInfo list from varsVec and covarsVec", {
  varsVec <- c("ka", "cl", "v")
  covarsVec <- c("WT", "BMI")

  funstring1 <- .cur$scmBuildCovInfo(varsVec, covarsVec)
  funstring2 <- .cur$scmBuildCovInfo(varsVec, covarsVec)[[1]]
  expect_length(funstring1, 6)
  expect_length(funstring2, 2)
})

test_that("Build covInfo list from varsVec and covarsVec", {
  varsVec <- c("cl", "v1", "v2", "Q")
  covarsVec <- c("WT", "BMI")

  funstring1 <- .cur$scmBuildCovInfo(varsVec, covarsVec)
  expect_error(expect_length(.cur$scmBuildCovInfo(varsVec, covarsVec), 6))
})


# ==== Build updated from the covariate and variable vector list

test_that("Build updated from the covariate and variable vector list", {
  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45
      label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui1 <- nlmixr2est::nlmixr(one.cmt)
  ui <- ui1
  varsVec <- c("ka", "cl", "v")
  covarsVec <- c("WT", "BMI")

  funstring1 <- intersect(
    (.cur$scmBuildUpdatedUi(
      ui1,
      varsVec,
      covarsVec,
      add = TRUE,
      indep = FALSE
    ))$iniDf$name,
    c(
      "cov_WT_ka",
      "cov_WT_cl",
      "cov_WT_v",
      "cov_BMI_ka",
      "cov_BMI_cl",
      "cov_BMI_v"
    )
  )
  funstring2 <- c(
    "cov_WT_ka",
    "cov_WT_cl",
    "cov_WT_v",
    "cov_BMI_ka",
    "cov_BMI_cl",
    "cov_BMI_v"
  )
  ui2 <- .cur$scmBuildUpdatedUi(
    ui1,
    varsVec,
    covarsVec,
    add = TRUE,
    indep = FALSE
  )
  funstring3 <- vapply(
    ui2$lstExpr,
    function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "),
    character(1L)
  )
  expect_equal(funstring1, funstring2)
  expect_true(any(grepl("cov_ka <- WT_ka \\+ BMI_ka", funstring3)))
  expect_true(any(grepl("cov_cl <- WT_cl \\+ BMI_cl", funstring3)))
  expect_true(all(
    c("cov_WT_ka", "cov_BMI_ka", "cov_WT_cl", "cov_BMI_cl") %in%
      ui2$muRefCovariateDataFrame$covariateParameter
  ))
})

test_that("Build ui from the covariate", {
  two.compartment <- function() {
    ini({
      tcl <- log(53.4) # Log Cl
      tv1 <- log(73.6)
      tv2 <- log(320) # Log V
      tQ <- log(191)
      eta.cl ~ 0.43^2
      eta.v1 ~ 0.48^2
      eta.v2 ~ 0.49^2
      eta.Q ~ 0.36^2
      prop.sd <- 0.44^2
    })
    # and a model block with the error specification and model specification
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      v2 <- exp(tv2 + eta.v2)
      Q <- exp(tQ + eta.Q)
      linCmt() ~ prop(prop.sd)
    })
  }

  ui2 <- nlmixr2est::nlmixr(two.compartment)
  ui <- ui2
  varsVec <- "ka"
  covarsVec <- c("WT", "BMI")
  expect_error(
    .cur$scmBuildUpdatedUi(ui, varsVec, covarsVec, add = TRUE, indep = FALSE)
  )
})


# ==== Make dummy variable cols and updated covarsVec

test_that("Make dummy variable cols and updated covarsVec", {
  # CMT is record-level, not subject-level: every subject's first record is a
  # dose (CMT 1), so under the default per-ID counting it collapses to a single
  # level and no indicator column is built.
  covarsVec <- c("WT", "BMI")
  catcovarsVec <- "CMT"
  funstring1 <- .cur$scmAddCatCovariates(
    nlmixr2data::theo_sd,
    covarsVec,
    catcovarsVec
  )[[2]]
  expect_equal(funstring1, c("WT", "BMI"))
})

test_that("Make dummy variable cols and updated covarsVec, freqBy observation", {
  # theo_sd$CMT is 12 dose rows at level 1 against 132 observation rows at
  # level 2, so counting by observation makes level 2 the reference and builds
  # only CMT_1.
  covarsVec <- c("WT", "BMI")
  catcovarsVec <- "CMT"
  funstring1 <- .cur$scmAddCatCovariates(
    nlmixr2data::theo_sd,
    covarsVec,
    catcovarsVec,
    freqBy = "observation"
  )[[2]]
  expect_equal(intersect(funstring1, "CMT_1"), "CMT_1")
  expect_length(intersect(funstring1, "CMT_2"), 0)
})

test_that("Make dummy variable cols and updated covarsVec, freqBy auto", {
  # CMT varies within subject, so "auto" selects per-observation counting
  covarsVec <- c("WT", "BMI")
  catcovarsVec <- "CMT"
  funstring1 <- .cur$scmAddCatCovariates(
    nlmixr2data::theo_sd,
    covarsVec,
    catcovarsVec,
    freqBy = "auto"
  )[[2]]
  expect_equal(intersect(funstring1, "CMT_1"), "CMT_1")
})
