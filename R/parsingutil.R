#' Built-in covariate shape expression builders for SCM
#'
#' A named list of functions \code{function(col, center, level)} that return the
#' covariate expression string embedded inside the model body.  Each function is
#' multiplied by the covariate theta: \code{cov_theta * covExpr}.
#'
#' \describe{
#'   \item{power}{  \code{log(col/center)} -- power/allometric scaling;
#'     on the exponentiated scale this gives
#'     \code{baseline * (col/center)^theta}.}
#'   \item{lin}{    \code{(col - center)} -- linear on the log scale
#'     (exponential in natural scale);
#'     \code{baseline * exp(theta * (col - center))}.}
#'   \item{log}{    \code{log(col)} -- uncentered log transform.}
#'   \item{identity}{\code{col} -- raw covariate on the log scale.}
#'   \item{cat}{    \code{col_level} -- pre-computed 0/1 indicator column
#'     (e.g. \code{sex_male}).  The column must already exist in the dataset,
#'     typically created during SCM categorical preprocessing.  The theta multiplies
#'     this indicator directly on the log scale, giving mu-referencing
#'     compatible code: \code{log(par) = theta_pop + eta + theta_cat * indicator}.}
#' }
#'
#' Users may pass additional shapes via the \code{customShapes} argument of
#' \code{runSCM()}.  Custom shape functions must accept named arguments
#' \code{col}, \code{center}, and \code{level}.
#'
#' @noRd
.SCM_SHAPES <- list(
  power = function(col, center, level = NULL) {
    paste0("log(", col, "/", center, ")")
  },
  lin = function(col, center, level = NULL) {
    paste0("(", col, " - ", center, ")")
  },
  log = function(col, center, level = NULL) paste0("log(", col, ")"),
  identity = function(col, center, level = NULL) col,
  cat = function(col, center = NULL, level) {
    # Return the name of the pre-computed 0/1 indicator column (e.g. "sex_male").
    # The column must exist in the dataset before model fitting, typically
    # created during SCM categorical preprocessing.  This keeps the model body mu-referencing
    # compatible with log-linear parameterisation.
    paste0(col, "_", level)
  }
)

#' Rebuild a UI function from metadata, ini, and model expressions
#'
#' @param ui rxode2 UI object
#' @param ini ini expression
#' @param model model expression
#' @return function that rebuilds the UI when called
#' @noRd
.getUiFunFromIniAndModel <- function(ui, ini, model) {
  .ls <- ls(ui$meta, all.names = TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(
      text = paste(
        "quote(",
        .ls[.i],
        "<-",
        deparse1(ui$meta[[.ls[.i]]]),
        ")"
      )
    ))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- ini
  .ret[[.len + 3]] <- model
  .retf <- function() {}
  body(.retf) <- as.call(.ret)
  .retf
}

#' Parse an init specification into a canonical list(est, lower, upper)
#'
#' Accepts either a scalar numeric (lower/upper default to -5/5) or a named
#' list with elements \code{est} (or \code{init}), \code{lower}, and
#' \code{upper}.  \code{NULL} is treated as the null-effect default
#' (est=0.1, lower=-5, upper=5).
#' @param spec scalar, named list, or NULL
#' @return list(est, lower, upper)
#' @noRd
.parseInitSpec <- function(spec) {
  if (is.null(spec)) {
    return(list(est = 0.1, lower = -5, upper = 5))
  }
  if (is.numeric(spec) && length(spec) == 1L) {
    return(list(est = spec, lower = -5, upper = 5))
  }
  if (is.list(spec)) {
    est <- if (!is.null(spec[["est"]])) {
      spec[["est"]]
    } else if (!is.null(spec[["init"]])) {
      spec[["init"]]
    } else {
      0.1
    }
    lower <- if (!is.null(spec[["lower"]])) spec[["lower"]] else -5
    upper <- if (!is.null(spec[["upper"]])) spec[["upper"]] else 5
    return(list(est = est, lower = lower, upper = upper))
  }
  list(est = 0.1, lower = -5, upper = 5)
}

#' Auto-scale covariate theta init spec for non-dimensionless shapes
#'
#' For \code{lin}, \code{log}, and \code{identity} shapes the default theta
#' bounds \code{(-5, 5)} are wildly inappropriate when the covariate has
#' physiological units (e.g. weight in kg).  This function scales the defaults
#' to be numerically equivalent to the \code{power} shape's defaults.
#'
#' The scaling is derived from the first-order equivalence
#' \code{log(x/c) \approx (x-c)/c}, so
#' \code{theta_lin = theta_power / |center|}.
#'
#' @param spec scalar, named list, or \code{NULL}. When non-\code{NULL},
#'   passed directly to \code{.parseInitSpec()} with no scaling.
#' @param shape character; one of the built-in shape names
#' @param center numeric or \code{NULL}; centering value (covariate median)
#' @param label character; used in warning messages to identify the pair
#' @return list(est, lower, upper)
#' @noRd
.autoScaleInitSpec <- function(spec, shape, center, label = "") {
  if (!is.null(spec)) {
    return(.parseInitSpec(spec))
  }
  if (!shape %in% c("lin", "log", "identity")) {
    return(.parseInitSpec(NULL))
  }
  if (is.null(center) || is.na(center) || center == 0) {
    cli::cli_warn(c(
      "!" = "Cannot auto-scale '{shape}' bounds{if (nzchar(label)) paste0(' for ', label) else ''}: center is NA, zero, or missing.",
      "i" = "Using unscaled defaults. Supply custom 'inits' for this shape to suppress."
    ))
    return(.parseInitSpec(NULL))
  }
  scale <- switch(shape,
    lin      = abs(center),
    identity = abs(center),
    log      = {
      lc <- abs(log(abs(center)))
      if (lc < .Machine$double.eps) {
        cli::cli_warn(c(
          "!" = "Cannot auto-scale 'log' bounds{if (nzchar(label)) paste0(' for ', label) else ''}: |log(center)| is approximately 0 (center is approximately 1).",
          "i" = "Using unscaled defaults."
        ))
        return(.parseInitSpec(NULL))
      }
      lc
    }
  )
  list(est = 0.1 / scale, lower = -5 / scale, upper = 5 / scale)
}

#' Build the covariate expression string for a given shape name
#'
#' @param shape  character; key into \code{.SCM_SHAPES} or \code{customShapes}
#' @param col    character; data column name (e.g. \code{"wt"}, \code{"sex"})
#' @param center numeric or NULL; centering value for continuous shapes
#' @param level  character or NULL; level for categorical shapes
#' @param customShapes named list of additional shape builder functions; merged
#'   with \code{.SCM_SHAPES} (custom entries take precedence)
#' @return character expression string
#' @noRd
.scmCovExpr <- function(
  shape,
  col,
  center = NULL,
  level = NULL,
  customShapes = NULL
) {
  all_shapes <- c(customShapes, .SCM_SHAPES) # custom takes precedence
  fn <- all_shapes[[shape]]
  if (is.null(fn)) {
    stop(
      "Unknown SCM shape '",
      shape,
      "'. Available: ",
      paste(names(all_shapes), collapse = ", "),
      call. = FALSE
    )
  }
  fn(col = col, center = center, level = level)
}

#' Population-typical missing-value fill string for a covariate shape
#'
#' Returns the model-body literal that the covariate expression takes when
#' the covariate equals its centering value (the observed median for the
#' built-in continuous shapes).  This is used by \code{.expandShapes()} to
#' wrap continuous covariate expressions in an \code{ifelse()} guard so that
#' missing observations contribute the population-typical effect (i.e. as if
#' the subject had \code{cov = median}).
#'
#' Built-in shapes are hard-coded for clarity and numerical efficiency:
#' \describe{
#'   \item{power}{   \code{"0"}  (\code{log(m/m) = 0})}
#'   \item{lin}{     \code{"0"}  (\code{m - m   = 0})}
#'   \item{log}{     numeric literal for \code{log(m)}}
#'   \item{identity}{numeric literal for \code{m}}
#' }
#' Custom shapes fall back to invoking the shape function with \code{col} set
#' to the formatted center value (so the returned string is the shape
#' expression evaluated symbolically at \code{cov = center}).
#'
#' Using \code{"0"} unconditionally (as a previous implementation did) is
#' only correct for the centered shapes (\code{power}, \code{lin}); for
#' \code{log} it would correspond to imputing \code{cov = 1} (since
#' \code{log(1) = 0}) and for \code{identity} it would impute \code{cov = 0}
#' -- neither of which is the population median.
#'
#' @param shape character; shape name (one of the built-ins or a custom key)
#' @param center numeric or \code{NULL}; centering value (covariate median)
#' @param level character or \code{NULL}; passed through to custom shape functions
#' @param customShapes optional named list of additional shape builder functions
#' @return character string suitable for embedding in an \code{ifelse()} body
#' @noRd
.scmFillStr <- function(shape, center, level = NULL, customShapes = NULL) {
  if (is.null(center) || !is.finite(center)) {
    return("0")
  }
  switch(
    shape,
    power = "0",
    lin = "0",
    log = if (center > 0) sprintf("%.15g", log(center)) else "0",
    identity = sprintf("%.15g", center),
    {
      shape_fn <- if (!is.null(customShapes) && shape %in% names(customShapes)) {
        customShapes[[shape]]
      } else {
        .SCM_SHAPES[[shape]]
      }
      if (is.null(shape_fn)) {
        return("0")
      }
      shape_fn(
        col = sprintf("%.15g", center),
        center = center,
        level = level
      )
    }
  )
}

#' Rebuild a UI by applying a set of covariate relationships to a clean base UI
#'
#' Starting from \code{base_ui} (always the original model without any
#' covariates), applies each row of \code{pairs_df} in sequence via
#' \code{.builduiCovariate()}.  This avoids the tainted-UI problem that arises
#' when trying to add a second covariate to a model that already has one.
#'
#' @param base_ui  the clean base rxode2 UI (no covariates added)
#' @param pairs_df data frame with columns \code{var}, \code{covar}, and
#'   optionally \code{covExpr} and \code{init}
#' @return updated UI with all covariate relationships applied, or
#'   \code{base_ui} unchanged when \code{pairs_df} is \code{NULL} or empty
#' @noRd
.rebuildUiFromPairs <- function(base_ui, pairs_df) {
  if (is.null(pairs_df) || nrow(pairs_df) == 0L) {
    return(base_ui)
  }

  base_ui <- rxode2::rxUiDecompress(base_ui)

  # Single-pass expansion from the clean base model.  Calling .builduiCovariate()
  # in a loop fails because the first call produces a tainted UI (tv moves from
  # pureMuRef to taintMuRef), so the second call can no longer find tv in
  # muRefDataFrame and produces a malformed "+".  Here we build the full covDf
  # for ALL pairs at once and process the base model body in one pass.
  cov_rows <- lapply(seq_len(nrow(pairs_df)), function(i) {
    var_i <- pairs_df$var[i]
    covar_i <- pairs_df$covar[i]
    expr_i <- {
      cv <- if ("covExpr" %in% names(pairs_df)) {
        pairs_df$covExpr[i]
      } else {
        NA_character_
      }
      if (is.null(cv) || is.na(cv)) covar_i else cv
    }
    theta_i <- .getThetaName(base_ui, var_i)
    data.frame(
      theta = theta_i,
      covariate = expr_i,
      covariateParameter = paste0("cov_", covar_i, "_", var_i),
      stringsAsFactors = FALSE
    )
  })
  .covDf <- do.call(rbind, cov_rows)

  .murefDf <- base_ui$muRefDataFrame
  .split <- base_ui$getSplitMuModel
  .pars <- intersect(
    c(names(.split$pureMuRef), names(.split$taintMuRef)),
    .murefDf$theta
  )
  .model <- nlmixr2est::.saemDropMuRefFromModel(base_ui)

  .model_list <- lapply(
    .model,
    .expandRefMu,
    murefDf = .murefDf,
    covDf = .covDf,
    pars = .pars
  )
  .model_list <- .injectCovBlocks(.model_list, .covDf, .murefDf)

  .newModel <- eval(parse(
    text = paste0(
      "quote(model({",
      paste0(as.character(.model_list), collapse = "\n"),
      "}))"
    )
  ))

  .iniDf <- base_ui$iniDf
  nthetaLength <- length(which(!is.na(.iniDf$ntheta)))
  new_ini_rows <- lapply(seq_len(nrow(.covDf)), function(i) {
    init_i <- if ("init" %in% names(pairs_df)) pairs_df$init[i] else NA_real_
    if (is.null(init_i) || is.na(init_i)) {
      init_i <- 0.1
    }
    lower_i <- if ("lower" %in% names(pairs_df)) pairs_df$lower[i] else NA_real_
    if (is.null(lower_i) || is.na(lower_i)) {
      lower_i <- -5
    }
    upper_i <- if ("upper" %in% names(pairs_df)) pairs_df$upper[i] else NA_real_
    if (is.null(upper_i) || is.na(upper_i)) {
      upper_i <- 5
    }
    data.frame(
      ntheta = as.integer(nthetaLength + i),
      neta1 = NA_character_,
      neta2 = NA_character_,
      name = .covDf$covariateParameter[i],
      lower = lower_i,
      est = init_i,
      upper = upper_i,
      fix = FALSE,
      label = NA_character_,
      backTransform = NA_character_,
      condition = NA_character_,
      err = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  new_rows <- do.call(rbind, new_ini_rows)
  # rxode2/nlmixr2est occasionally add columns to iniDf (e.g. "prior"); pad
  # rather than hardcode the full column set so new columns don't break rbind.
  for (col in setdiff(names(.iniDf), names(new_rows))) new_rows[[col]] <- NA
  .ini <- rbind(.iniDf, new_rows[, names(.iniDf), drop = FALSE])
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)

  .ui <- rxode2::rxUiDecompress(.getUiFunFromIniAndModel(
    base_ui,
    .ini,
    .newModel
  )())
  .annotateCovariateUi(.ui, .covDf)
}

# Add or remove one covariate relationship from an SCM UI object.
scmAddOrRemoveCovariate <- function(
  ui,
  varName,
  covariate,
  add = TRUE,
  covExpr = NULL
) {
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)

  checkmate::assertCharacter(varName, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(covariate, len = 1, any.missing = FALSE)

  if (inherits(try(str2lang(varName)), "try-error")) {
    stop("`varName` must be a valid R expression", call. = FALSE)
  }
  if (inherits(try(str2lang(covariate)), "try-error")) {
    stop("`varName` must be a valid R expression", call. = FALSE)
  }
  .pop <- .getThetaName(ui, varName = varName)
  .cov <- paste0("cov_", covariate, "_", varName)
  # covExpr is the expression embedded in the model (may differ from the
  # covariate name used for theta naming, e.g. "log(wt/70)" or
  # 'ifelse(sex=="male",1,0)').  Defaults to the raw covariate name.
  .expr <- if (!is.null(covExpr)) covExpr else covariate

  if (add) {
    .covdf <- rbind(
      ui$muRefCovariateDataFrame,
      data.frame(theta = .pop, covariate = .expr, covariateParameter = .cov)
    )
  } else {
    .covdf <- ui$muRefCovariateDataFrame[
      ui$muRefCovariateDataFrame$covariateParameter != .cov,
    ]
  }
  .split <- ui$getSplitMuModel
  .pars <- intersect(
    c(names(.split$pureMuRef), names(.split$taintMuRef)),
    ui$muRefDataFrame$theta
  )
  .model <- nlmixr2est::.saemDropMuRefFromModel(ui)
  .model_list <- lapply(
    .model,
    .expandRefMu,
    murefDf = ui$muRefDataFrame,
    covDf = .covdf,
    pars = .pars
  )
  if (nrow(.covdf) > 0L) {
    .model_list <- .injectCovBlocks(.model_list, .covdf, ui$muRefDataFrame)
  }
  .model_list
}


#' Given model expression, Expand population expression
#'
#' @param x model expression
#' @param murefDf MU referencing data frame
#' @param covDf covariate referencing data frame
#' @param pars theta parameters
#'
#' @noRd
#' @author Matthew Fidler, Vishal Sarsani
#' @noRd
.expandRefMu <- function(x, murefDf, covDf, pars) {
  if (is.name(x)) {
    currparam <- as.character(x)
    if (currparam %in% pars) {
      return(str2lang(.expandPopExpr(currparam, murefDf, covDf)))
    }
  } else if (is.call(x)) {
    return(as.call(c(
      list(x[[1]]),
      lapply(x[-1], .expandRefMu, murefDf = murefDf, covDf = covDf, pars = pars)
    )))
  }
  x
}


#' Derive the PK parameter name from the population theta name
#'
#' @param popParam character; population theta name (e.g. \code{"tv"})
#' @param murefDf  the muRefDataFrame from the UI (columns \code{theta}, \code{eta})
#' @return character; PK parameter name (e.g. \code{"v"})
#' @noRd
.covVarNameFromTheta <- function(popParam, murefDf) {
  eta <- murefDf$eta[murefDf$theta == popParam]
  if (length(eta) == 0L) {
    return(popParam)
  }
  sub("^eta\\.", "", eta[[1L]])
}

#' Derive the intermediate variable name from a covariateParameter string
#'
#' Strips the \code{"cov_"} prefix so that e.g. \code{"cov_wt_power_v"} becomes
#' \code{"wt_power_v"} and \code{"cov_sex_male_v"} becomes \code{"sex_male_v"}.
#'
#' @param covariateParameter character; theta name from \code{muRefCovariateDataFrame}
#' @return character; intermediate variable name
#' @noRd
.covIntermediateName <- function(covariateParameter) {
  sub("^cov_", "", covariateParameter)
}

#' Build the covariate block lines for a single parameter
#'
#' Generates one intermediate assignment per covariate and a final
#' \code{cov_<var>} aggregate line.  Example output for weight (power) and
#' sex on \code{v}:
#' \preformatted{
#'   wt_power_v = log(wt/70.5) * cov_wt_power_v
#'   sex_male_v = sex_male * cov_sex_male_v
#'   cov_v = wt_power_v + sex_male_v
#' }
#'
#' @param popParam character; population theta name (e.g. \code{"tv"})
#' @param murefDf  muRefDataFrame
#' @param covDf    full covariate reference data frame
#' @param factor   optional scalar multiplied into each covariate theta (LASSO)
#' @return character vector of lines, or \code{character(0)} when no covariates
#' @noRd
.buildCovBlock <- function(popParam, murefDf, covDf, factor = NULL) {
  varName <- .covVarNameFromTheta(popParam, murefDf)
  rows <- covDf[covDf$theta == popParam, , drop = FALSE]
  if (nrow(rows) == 0L) {
    return(character(0L))
  }

  int_names <- .covIntermediateName(rows$covariateParameter)

  int_lines <- if (!is.null(factor)) {
    paste0(
      int_names,
      " <- ",
      rows$covariate,
      " * ",
      rows$covariateParameter,
      " * ",
      factor
    )
  } else {
    paste0(int_names, " <- ", rows$covariate, " * ", rows$covariateParameter)
  }

  agg_line <- paste0(
    "cov_",
    varName,
    " <- ",
    paste(int_names, collapse = " + ")
  )
  c(int_lines, agg_line)
}

#' Carry covariate metadata onto a rebuilt UI
#'
#' rxode2 preserves the compiled model text but does not reconstruct
#' \code{muRefCovariateDataFrame} from the generated covariate helper lines, so
#' we attach the explicit covariate mapping we used to build the model.
#'
#' @param ui rebuilt UI
#' @param covDf covariate reference data frame
#' @return rebuilt UI with covariate metadata attached
#' @noRd
.annotateCovariateUi <- function(ui, covDf) {
  covDf <- covDf[, c("theta", "covariate", "covariateParameter"), drop = FALSE]
  rownames(covDf) <- NULL

  ui$muRefCovariateDataFrame <- covDf
  ui$muRefTable <- cbind(
    ui$muRefDataFrame,
    covariates = vapply(
      ui$muRefDataFrame$theta,
      function(theta) {
        rows <- covDf[covDf$theta == theta, , drop = FALSE]
        if (nrow(rows) == 0L) {
          return(NA_character_)
        }
        paste(
          paste0(rows$covariate, "*", rows$covariateParameter),
          collapse = " + "
        )
      },
      character(1L)
    ),
    stringsAsFactors = FALSE
  )
  ui
}

#' Inject covariate block lines into the model expression list
#'
#' For each PK parameter that has entries in \code{covDf}, inserts the
#' intermediate variable assignments and the \code{cov_<var>} aggregate line
#' immediately before that parameter's assignment line.
#'
#' @param model_list list of R language objects (model body statements)
#' @param covDf      covariate reference data frame
#' @param murefDf    muRefDataFrame
#' @param factor     optional LASSO constraint factor
#' @return modified \code{model_list}
#' @noRd
.injectCovBlocks <- function(model_list, covDf, murefDf, factor = NULL) {
  params_with_covs <- unique(covDf$theta)
  strs <- vapply(
    model_list,
    function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "),
    character(1L)
  )

  # Collect (insertion_index, block_exprs) pairs; insert in reverse order so
  # earlier insertions do not shift the indices of later ones.
  insertions <- vector("list", length(params_with_covs))

  for (k in seq_along(params_with_covs)) {
    popParam <- params_with_covs[[k]]
    varName <- .covVarNameFromTheta(popParam, murefDf)
    cov_var <- paste0("cov_", varName)

    block <- .buildCovBlock(popParam, murefDf, covDf, factor)
    if (length(block) == 0L) {
      next
    }

    idx <- which(grepl(cov_var, strs, fixed = TRUE))
    if (length(idx) == 0L) {
      next
    }
    idx <- idx[[1L]]

    block_exprs <- lapply(block, function(ln) {
      tryCatch(str2lang(ln), error = function(e) parse(text = ln)[[1L]])
    })

    insertions[[k]] <- list(idx = idx, exprs = block_exprs)
  }

  insertions <- Filter(Negate(is.null), insertions)
  if (length(insertions) == 0L) {
    return(model_list)
  }

  # Sort descending so we insert from the bottom up
  ord <- order(vapply(insertions, `[[`, integer(1L), "idx"), decreasing = TRUE)
  insertions <- insertions[ord]

  for (ins in insertions) {
    idx <- ins$idx
    exprs <- ins$exprs
    model_list <- c(
      model_list[seq_len(idx - 1L)],
      exprs,
      model_list[idx:length(model_list)]
    )
  }

  model_list
}

#' Expand population expression given the mu reference and covariate reference Data frames
#'
#' When covariates are present the reference is replaced with
#' \code{theta+eta+cov_<var>}; the detailed covariate expressions are
#' generated as separate lines by \code{.buildCovBlock()} /
#' \code{.injectCovBlocks()} in \code{scmAddOrRemoveCovariate()}.
#'
#' @param popParam  population parameter variable
#' @param murefDf MU referencing data frame
#' @param covDf covariate referencing data frame
#' @param factor unused; kept for signature compatibility
#'
#' @return expanded expression string
#' @author Matthew Fidler,Vishal Sarsani
#' @noRd
.expandPopExpr <- function(popParam, murefDf, covDf, factor = NULL) {
  .par1 <- murefDf[murefDf$theta == popParam, ]
  .res <- paste0(.par1$theta, "+", .par1$eta)
  .w <- which(covDf$theta == popParam)
  if (length(.w) > 0L) {
    varName <- .covVarNameFromTheta(popParam, murefDf)
    .res <- paste0(.res, "+cov_", varName)
  }
  .res
}

#' get the population parameter from variable name
#'
#' @param ui compiled rxode2 nlmir2 model or fi
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return population parameter variable
#'
#' @author Matthew Fidler
#' @noRd
.getThetaName <- function(ui, varName) {
  .split <- ui$getSplitMuModel
  if (varName %in% names(.split$pureMuRef)) {
    return(varName)
  }
  .w <- which(.split$pureMuRef == varName)
  if (length(.w) == 1) {
    return(names(.split$pureMuRef)[.w])
  }
  if (varName %in% names(.split$taintMuRef)) {
    return(varName)
  }
  # taintMuRef values are PK param names (e.g. "v") whose theta has a covariate;
  # the original check only looked at names (theta names), so we also check values
  .w2 <- which(.split$taintMuRef == varName)
  if (length(.w2) == 1) {
    return(names(.split$taintMuRef)[.w2])
  }
  # Fallback: use muRefDataFrame which maps theta names to eta names (e.g. tv -> eta.v).
  # This works even after a covariate has been added (tv moves from pureMuRef to
  # taintMuRef and taintMuRef values may not be plain PK param name strings).
  # Guard: only proceed if varName is an actual named LHS assignment in the model body
  # (e.g. `ka <- exp(tka + eta.ka)`). Without this check, models where eta.cl is
  # declared in the ini block but cl is used inline in an ODE (not as a named variable)
  # would incorrectly match and return a theta name.
  .lhs_vars <- vapply(
    ui$lstExpr,
    function(e) {
      if (
        is.call(e) && identical(e[[1L]], as.symbol("<-")) && is.name(e[[2L]])
      ) {
        as.character(e[[2L]])
      } else {
        NA_character_
      }
    },
    character(1L)
  )
  .lhs_vars <- .lhs_vars[!is.na(.lhs_vars)]
  if (!(varName %in% .lhs_vars)) {
    stop("'", varName, "'", "has not been found in the model ui", call. = FALSE)
  }
  .mdf <- ui$muRefDataFrame
  .eta_name <- paste0("eta.", varName)
  .w3 <- which(.mdf$eta == .eta_name)
  if (length(.w3) >= 1) {
    return(.mdf$theta[.w3[1]])
  }
  stop("'", varName, "'", "has not been found in the model ui", call. = FALSE)
}


#' Given a data frame extract column corresponding to  Individual
#'
#' @param data given data frame
#'
#' @return column name of individual
#' @noRd
#' @author  Vishal Sarsani
.idColumn <- function(data) {
  # check if it is a dataframe
  checkmate::assertDataFrame(data, col.names = "named")
  # Extract individual ID from column names
  colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uidCol <- colNames[match("id", colNamesLower)]
  } else {
    uidCol <- "ID"
  }
  uidCol
}

#' Build ui from the covariate
#'
#' @import lotri
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate name used for theta naming (e.g. \code{"wt"},
#'   \code{"sex_male"})
#' @param add boolean indicating if the covariate needs to be added or removed
#' @param type \code{"continuous"} (default) or \code{"categorical"}.  Controls
#'   how the model expression is generated.
#' @param center for continuous covariates, the centering value (typically the
#'   population median).  When supplied the model expression becomes
#'   \code{cov_theta * log(raw_col / center)}.
#' @param raw_col the data column name to use in the model expression.  Defaults
#'   to \code{covariate}.  For categorical covariates this is the original
#'   factor/character column (e.g. \code{"sex"}), while \code{covariate} carries
#'   the level-suffixed name used for the theta (e.g. \code{"sex_male"}).
#' @param level for categorical covariates, the level being tested (e.g.
#'   \code{"male"}).  The model expression becomes
#'   \code{cov_theta * ifelse(raw_col == "level", 1, 0)}.
#' @param covExpr optional pre-computed covariate expression string.  When
#'   supplied, the \code{type}/\code{center}/\code{level} shape logic is
#'   skipped and this string is used verbatim as the expression multiplied by
#'   the covariate theta.  Intended for use by \code{.fitCandidatePairs()} after
#'   \code{.expandShapes()} has already determined the expression.
#' @param init numeric starting estimate for the covariate theta; default
#'   \code{0} (no covariate effect).  Set via \code{.expandShapes()} from the
#'   \code{inits} argument of \code{runSCM()}.
#' @return ui with added or removed covariate
#' @noRd
#' @author  Vishal Sarsani
.builduiCovariate <- function(
  ui,
  varName,
  covariate,
  add = TRUE,
  type = "continuous",
  center = NULL,
  raw_col = NULL,
  level = NULL,
  covExpr = NULL,
  init = 0.1,
  lower = -5,
  upper = 5
) {
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$finalUiEnv
  }
  ui <- rxode2::assertRxUi(ui)
  ui <- rxode2::rxUiDecompress(ui)

  checkmate::assertCharacter(varName, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(covariate, len = 1, any.missing = FALSE)

  if (inherits(try(str2lang(varName)), "try-error")) {
    stop("`varName` must be a valid R expression", call. = FALSE)
  }
  if (inherits(try(str2lang(covariate)), "try-error")) {
    stop("`covariate` must be a valid R expression", call. = FALSE)
  }

  # The data column used in the model expression (may differ from covariate for
  # categorical: covariate = "sex_male", raw_col = "sex", level = "male")
  rc <- if (!is.null(raw_col)) raw_col else covariate

  # Build the covariate expression that appears in the model body.
  # For continuous (centred) covariates the expression is log(col / center).
  # For categorical covariates the expression is the pre-computed 0/1 indicator
  # column name, which must already exist in the data.
  # For plain (uncentred) the raw column name is used as the expression.
  # When covExpr is supplied by the caller (e.g. from .expandShapes()), use it
  # directly rather than recomputing from type/center/level.
  if (is.null(covExpr)) {
    covExpr <- if (type == "categorical" && !is.null(level)) {
      paste0(rc, "_", level)
    } else if (type == "continuous" && !is.null(center)) {
      paste0("log(", rc, "/", center, ")")
    } else {
      rc
    }
  }

  covName <- paste0("cov_", covariate, "_", varName)
  .pop <- .getThetaName(ui, varName = varName)
  .expr <- if (!is.null(covExpr)) covExpr else covariate
  .covdf <- if (add) {
    rbind(
      ui$muRefCovariateDataFrame,
      data.frame(theta = .pop, covariate = .expr, covariateParameter = covName)
    )
  } else {
    ui$muRefCovariateDataFrame[
      ui$muRefCovariateDataFrame$covariateParameter != covName,
    ]
  }

  # add covariate
  if (add) {
    lst <- scmAddOrRemoveCovariate(
      ui,
      varName,
      covariate,
      add = TRUE,
      covExpr = covExpr
    )
    .newModel <- eval(parse(
      text = paste0(
        "quote(model({",
        paste0(as.character(lst), collapse = "\n"),
        "}))"
      )
    ))
    nthetaLength <- length(which(!is.na(ui$iniDf$ntheta)))
    .ini <- ui$iniDf
    new_row <- data.frame(
      ntheta = as.integer(nthetaLength + 1),
      neta1 = NA_character_,
      neta2 = NA_character_,
      name = covName,
      lower = lower,
      est = init,
      upper = upper,
      fix = FALSE,
      label = NA_character_,
      backTransform = NA_character_,
      condition = NA_character_,
      err = NA_character_
    )
    # rxode2/nlmixr2est occasionally add columns to iniDf (e.g. "prior"); pad
    # rather than hardcode the full column set so new columns don't break rbind.
    for (col in setdiff(names(.ini), names(new_row))) new_row[[col]] <- NA
    .ini <- rbind(.ini, new_row[, names(.ini), drop = FALSE])
  } else {
    # remove covariate
    lst <- scmAddOrRemoveCovariate(ui, varName, covariate, add = FALSE)
    .newModel <- eval(parse(
      text = paste0(
        "quote(model({",
        paste0(as.character(lst), collapse = "\n"),
        "}))"
      )
    ))
    .ini <- ui$iniDf[ui$iniDf$name != covName, ]
  }

  # build ui
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)
  .ui <- rxode2::rxUiDecompress(.getUiFunFromIniAndModel(ui, .ini, .newModel)())
  .annotateCovariateUi(.ui, .covdf)
}

# Build covInfo list from varsVec and covarsVec.
scmBuildCovInfo <- function(varsVec, covarsVec) {
  checkmate::assert_character(varsVec, min.len = 1)
  checkmate::assert_character(covarsVec, min.len = 1)
  possiblePerms <- expand.grid(varsVec, covarsVec)
  possiblePerms <-
    list(
      as.character(possiblePerms[[1]]),
      as.character(possiblePerms[[2]])
    )
  names(possiblePerms) <- c("vars", "covars")
  covInfo <- list() # reversivle listVarName!!
  for (item in Map(list, possiblePerms$vars, possiblePerms$covars)) {
    listVarName <- paste0(item[[2]], item[[1]])
    covInfo[[listVarName]] <- list(varName = item[[1]], covariate = item[[2]])
  }
  covInfo
}


# Build an updated UI from variable and covariate vectors.
scmBuildUpdatedUi <- function(
  ui,
  varsVec,
  covarsVec,
  add = TRUE,
  indep = FALSE
) {
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$finalUiEnv
  }
  ui <- rxode2::assertRxUi(ui)
  ui <- rxode2::rxUiDecompress(ui)
  # construct covInfo
  covInfo <- scmBuildCovInfo(varsVec, covarsVec)
  # check if the covInfo is a list
  checkmate::assert_list(covInfo)
  covSearchRes <- list()
  covsAdded <- list() # to keep track of covariates added and store in a file
  if (add) {
    # Add covariates one after other
    if (indep) {
      ui_base <- ui # snapshot before the loop so each candidate starts fresh
      for (i in seq_along(covInfo)) {
        x <- covInfo[[i]]
        covName <- paste0("cov_", x$covariate, "_", x$varName)
        ui_candidate <- tryCatch(
          {
            res <- do.call(.builduiCovariate, c(ui_base, x))
            res # to return 'ui'
          },
          error = function(error_message) {
            message("error  while  ADDING covariate")
            message(error_message)
            message("skipping this covariate")
            res # return NA otherwise (instead of NULL)
          }
        )
        covSearchRes[[i]] <- list(
          ui_candidate,
          c(x$covariate, x$varName),
          covName
        )[[1]]
      }
      covSearchRes
    } else {
      ## Add all at once
      covsAddedIdx <- 1
      for (x in covInfo) {
        covName <- paste0("cov_", x$covariate, "_", x$varName)
        if (length(covsAdded) == 0) {
          covsAdded[[covsAddedIdx]] <- c(x$covariate, x$varName)
          ui <- tryCatch(
            {
              res <- do.call(.builduiCovariate, c(ui, x))
              res # to return 'ui'
            },
            error = function(error_message) {
              message("error  while SIMULTANEOUSLY ADDING covariates")
              message(error_message)
              message("skipping this covariate")
              res # return NA otherwise (instead of NULL)
            }
          )
          covSearchRes[[covsAddedIdx]] <- list(
            ui,
            covsAdded[[covsAddedIdx]],
            covName
          )
          covsAddedIdx <- covsAddedIdx + 1
        } else {
          covsAdded[[covsAddedIdx]] <-
            c(covsAdded[[covsAddedIdx - 1]], x$covariate, x$varName)
          ui <- tryCatch(
            {
              res <- do.call(.builduiCovariate, c(ui, x))
              res # to return 'ui'
            },
            error = function(error_message) {
              message("error  while SIMULTANEOUSLY ADDING covariates")
              message(error_message)
              message("skipping this covariate")
              res # return NA otherwise (instead of NULL)
            }
          )
          covSearchRes[[covsAddedIdx]] <- list(
            ui,
            covsAdded[[covsAddedIdx]],
            covName
          )
          covsAddedIdx <- covsAddedIdx + 1
        }
        covSearchRes
      }
      covSearchRes[length(covSearchRes)][[1]][[1]]
    }
  } else {
    # remove covariate one after other
    for (i in seq_along(covInfo)) {
      x <- covInfo[[i]]
      covName <- paste0("cov_", x$covariate, "_", x$varName)
      ui <- tryCatch(
        {
          res <- do.call(.builduiCovariate, c(ui, x, add = FALSE))
          res # to return 'ui'
        },
        error = function(error_message) {
          message("error  while  ADDING covariate")
          message(error_message)
          message("skipping this covariate")
          res # return NA otherwise (instead of NULL)
        }
      )
      covSearchRes[[i]] <- list(ui, c(x$covariate, x$varName), covName)[[1]]
    }
    covSearchRes
  }
}


#' Create dummy variable columns and return updated covariate metadata
#'
#' The reference level is the most frequent level, and non-reference levels
#' whose proportion falls below \code{catCutoff} are lumped with the reference
#' and receive no indicator column.  Previously the reference was whichever
#' level happened to appear first, which is not a frequency rule at all.
#'
#' How frequency is counted depends on the column, because the two available
#' rules are each right for a different kind of covariate:
#' \describe{
#'   \item{per subject}{One vote per unique ID, taken from each subject's first
#'     record.  Used when the column is constant within every subject, i.e. a
#'     genuine subject-level covariate such as sex or race.  This is the rule
#'     \code{.makeSCMData()} applies on the \code{runSCM()} search path, so
#'     both routes agree on the reference level even when subjects contribute
#'     unequal numbers of rows.  Counting such a column by row would let a few
#'     densely sampled subjects outvote the majority.}
#'   \item{per observation}{One vote per data row.  Used when the column varies
#'     within a subject, e.g. \code{CMT}.  Per-subject counting is meaningless
#'     for such a column -- every subject would contribute only whichever level
#'     happened to land in its first record.}
#' }
#' Both rules are available via \code{freqBy}, which defaults to \code{"id"} --
#' the rule \code{.makeSCMData()} uses, and the right one for the subject-level
#' covariates SCM actually tests.  Use \code{"observation"} for a column that is
#' not unique to a subject, or \code{"auto"} to choose per column by testing
#' whether any subject carries more than one non-missing level.  When
#' \code{data} has no recognisable ID column each row is treated as its own
#' subject, so all three settings reduce to per-row counting.
#'
#' @param data data frame containing the categorical columns
#' @param covarsVec character vector of covariate names to extend.  Any entry
#'   naming a column in \code{catcovarsVec} is dropped, since that column no
#'   longer exists in the returned data.
#' @param catcovarsVec character vector of categorical columns to expand
#' @param catCutoff minimum proportion of subjects (or of rows, when counting
#'   per observation) in a non-reference level for that level to get an
#'   indicator column.  Default \code{0}, which keeps every observed
#'   non-reference level.
#' @param freqBy how to count level frequencies: \code{"id"} (default) counts
#'   one vote per subject, \code{"observation"} one vote per row, and
#'   \code{"auto"} picks per column -- per subject when the column is constant
#'   within every subject, per observation when it varies within a subject.
#' @return list of the updated data (original categorical columns removed) and
#'   the updated covariate vector
#' @noRd
scmAddCatCovariates <- function(
  data,
  covarsVec,
  catcovarsVec,
  catCutoff = 0,
  freqBy = c("id", "observation", "auto")
) {
  # check for valid inputs
  checkmate::assert_data_frame(data)
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(catcovarsVec)
  checkmate::assert_number(catCutoff, lower = 0, upper = 1)
  freqBy <- match.arg(freqBy)

  idCol <- .idColumn(data)
  if (!idCol %in% names(data)) {
    idCol <- NULL
  }

  newcatvars <- character(0)
  # Deduplicate so a repeated column does not trip the collision check below
  for (col in unique(catcovarsVec)) {
    if (!col %in% names(data)) {
      stop(
        "Categorical covariate '",
        col,
        "' not found in data.",
        call. = FALSE
      )
    }

    keep <- !is.na(data[[col]])
    vals <- as.character(data[[col]])[keep]
    if (length(vals) == 0L) {
      next
    }

    # Choose the counting rule for this column.  Per-subject counting is right
    # for subject-level covariates; per-observation is the only sensible rule
    # for a column that is not unique to a subject, since otherwise a subject
    # contributes only whichever level landed in its first record.
    byId <- !is.null(idCol) && freqBy != "observation"
    if (byId && freqBy == "auto") {
      ids <- data[[idCol]][keep]
      distinctPairs <- !duplicated(data.frame(ids, vals))
      byId <- anyDuplicated(ids[distinctPairs]) == 0L
    }
    if (byId) {
      # One vote per subject, taken from that subject's first record
      vals <- vals[!duplicated(data[[idCol]][keep])]
    }

    tbl <- sort(table(vals), decreasing = TRUE)
    props <- tbl / sum(tbl)
    ref <- names(tbl)[[1L]] # most frequent level = reference
    lvls <- names(props)[props >= catCutoff & names(props) != ref]
    if (length(lvls) == 0L) {
      next
    }

    # Collision check -- stop rather than silently overwrite
    collisions <- paste0(col, "_", lvls)[
      paste0(col, "_", lvls) %in% names(data)
    ]
    if (length(collisions) > 0L) {
      stop(
        "Cannot create indicator column(s): ",
        paste(collisions, collapse = ", "),
        " -- name(s) already exist in data.",
        call. = FALSE
      )
    }

    for (lev in lvls) {
      colname <- paste0(col, "_", lev)
      data[[colname]] <- as.integer(
        !is.na(data[[col]]) & as.character(data[[col]]) == lev
      )
      newcatvars <- c(newcatvars, colname)
    }
  }
  # Remove original categorical variables
  updatedData <- data[, !(names(data) %in% catcovarsVec), drop = FALSE]
  # Update entire covarsvec with added categorical variables.  Drop the
  # expanded columns themselves, which no longer exist in updatedData --
  # keeping them would make updatedData[, updcovarsVec] an error.
  updcovarsVec <- c(setdiff(covarsVec, catcovarsVec), newcatvars)
  list(updatedData, updcovarsVec)
}
