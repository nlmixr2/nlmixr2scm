#' Summarise a stepwise covariate modeling run
#'
#' \code{summary()} collects a human-readable report of a \code{\link{runSCM}}
#' result: the options used for the search, a comparison of the base model with
#' the final forward and final backward models, the step and all-candidate
#' tables, the covariates in the final model, and where models and report files
#' were saved on disk (if anywhere).  \code{print()} on the \code{runSCM()}
#' result itself gives a short overview.
#'
#' @param object,x an object of class \code{"nlmixr2scm"} returned by
#'   \code{\link{runSCM}} (or, for the \code{print()} method of the summary, an
#'   object of class \code{"summary.nlmixr2scm"}).
#' @param ... unused; for compatibility with the generics.
#' @return \code{summary()} returns an object of class
#'   \code{"summary.nlmixr2scm"}, a list with elements \code{comparison} (a
#'   data frame with one row per model: \code{model}, \code{OFV}, \code{dOFV}
#'   relative to the base model, \code{AIC}, \code{BIC}, \code{nPar} and the
#'   SCM \code{covariates} it contains), \code{summaryTable}, \code{options},
#'   \code{outputDir} and \code{finalFit}.  The \code{print()} methods return
#'   their input invisibly.
#' @examples
#' \dontrun{
#' scm <- runSCM(fit, varsVec = c("cl", "v"), covarsVec = "wt")
#' scm           # short overview
#' summary(scm)  # full report
#' }
#' @export
summary.nlmixr2scm <- function(object, ...) {
  st <- object$summaryTable
  fits <- list(
    "Base" = object$baseFit,
    "Forward final" = if (!is.null(object$resFwd)) object$resFwd[[1]],
    "Backward final" = if (!is.null(object$resBck)) object$resBck[[1]]
  )
  fits <- fits[!vapply(fits, is.null, logical(1))]
  base_objf <- object$baseFit$objf

  comparison <- data.frame(
    model = names(fits),
    OFV = vapply(fits, function(f) f$objf, numeric(1)),
    dOFV = vapply(fits, function(f) f$objf - base_objf, numeric(1)),
    AIC = vapply(fits, function(f) f$AIC, numeric(1)),
    BIC = vapply(fits, function(f) f$BIC, numeric(1)),
    nPar = vapply(fits, function(f) length(f$finalUiEnv$ini$est), integer(1)),
    covariates = vapply(fits, .scmFitCovariates, character(1), st = st),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      comparison = comparison,
      summaryTable = st,
      options = object$options,
      outputDir = object$outputDir,
      finalFit = fits[[length(fits)]]
    ),
    class = "summary.nlmixr2scm"
  )
}

#' @rdname summary.nlmixr2scm
#' @export
print.summary.nlmixr2scm <- function(x, ...) {
  cat(cli::rule(left = "nlmixr2scm summary"), "\n", sep = "")
  opts <- .scmOptionLines(x$options)
  cat(cli::rule(left = "Options"), "\n", sep = "")
  .catAligned(opts)

  cat(cli::rule(left = "Model comparison"), "\n", sep = "")
  cmp <- x$comparison
  disp <- data.frame(
    Model = cmp$model,
    OFV = formatC(cmp$OFV, format = "f", digits = 3),
    dOFV = formatC(cmp$dOFV, format = "f", digits = 3),
    AIC = formatC(cmp$AIC, format = "f", digits = 3),
    BIC = formatC(cmp$BIC, format = "f", digits = 3),
    Params = cmp$nPar,
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE, right = FALSE)
  cat("\nSCM covariates:\n")
  .catAligned(stats::setNames(cmp$covariates, paste0("  ", cmp$model)))

  .printFinalSCMSummary(
    x$summaryTable,
    finalFit = x$finalFit,
    searchType = x$options$searchType,
    outputDir = NULL
  )

  cat(cli::rule(left = "Files"), "\n", sep = "")
  cat(.scmFileLines(x$outputDir), sep = "\n")
  invisible(x)
}

#' @rdname summary.nlmixr2scm
#' @export
print.nlmixr2scm <- function(x, ...) {
  final <- if (!is.null(x$resBck)) x$resBck[[1]] else x$resFwd[[1]]
  base_objf <- x$baseFit$objf
  o <- x$options
  cat(cli::rule(left = "nlmixr2scm result"), "\n", sep = "")
  cat(sprintf(
    "Search      : %s (%s), %d candidate relation%s\n",
    o$searchType, o$estimation, nrow(o$candidates),
    if (nrow(o$candidates) == 1L) "" else "s"
  ))
  cat(sprintf(
    "OFV         : base %.3f -> final %.3f (dOFV %.3f)\n",
    base_objf, final$objf, final$objf - base_objf
  ))
  cat("Covariates  : ", .scmFitCovariates(final, x$summaryTable), "\n", sep = "")
  if (!is.null(x$outputDir)) cat("Saved to    : ", x$outputDir, "\n", sep = "")
  cat("Use summary() for options, model comparison and step tables.\n")
  invisible(x)
}

# -- Helpers -------------------------------------------------------------------

#' Print a named character vector as aligned "name : value" lines, wrapping
#' long values under the value column
#' @noRd
.catAligned <- function(x) {
  w <- max(nchar(names(x)))
  width <- max(getOption("width", 80L) - w - 3L, 20L)
  for (i in seq_along(x)) {
    lines <- strwrap(x[[i]], width = width)
    if (length(lines) == 0L) lines <- ""
    cat(sprintf("%-*s : %s", w, names(x)[i], lines[1]), sep = "\n")
    if (length(lines) > 1L) {
      cat(paste0(strrep(" ", w + 3L), lines[-1]), sep = "\n")
    }
  }
}

#' Relation label, e.g. "wt_power~cl [power]"
#' @noRd
.scmRelLabel <- function(covar, var, shape) {
  lbl <- paste0(covar, "~", var)
  has_shape <- !is.na(shape) & nzchar(shape)
  lbl[has_shape] <- paste0(lbl[has_shape], " [", shape[has_shape], "]")
  lbl
}

#' Comma-separated SCM relations contained in a fit
#'
#' Covariate thetas are matched to their relation through the summary table's
#' \code{covNames}; any \code{cov_} theta not in the table is shown by name.
#' @noRd
.scmFitCovariates <- function(fit, st) {
  thetas <- fit$iniDf$name
  lbl <- character(0)
  if (!is.null(st) && nrow(st) > 0L) {
    st <- st[!duplicated(st$covNames), , drop = FALSE]
    lbl <- stats::setNames(.scmRelLabel(st$covar, st$var, st$shape), st$covNames)
  }
  cov_thetas <- thetas[thetas %in% names(lbl) | startsWith(thetas, "cov_")]
  if (length(cov_thetas) == 0L) return("(none)")
  shown <- ifelse(cov_thetas %in% names(lbl), lbl[cov_thetas], cov_thetas)
  paste(shown, collapse = ", ")
}

#' Named character vector of the options used for a search
#' @noRd
.scmOptionLines <- function(o) {
  rel_list <- function(pairs) {
    if (is.null(pairs) || nrow(pairs) == 0L) return("none")
    shape <- if ("shape" %in% names(pairs)) pairs$shape else NA_character_
    paste0(nrow(pairs), ": ",
           paste(.scmRelLabel(pairs$covar, pairs$var, shape), collapse = ", "))
  }
  centers <- if (length(o$centers) == 0L) {
    "median of each continuous covariate"
  } else {
    paste0(paste(names(o$centers), "=", o$centers, collapse = ", "),
           " (others: median)")
  }
  warm_start <- if (isTRUE(o$profileInit)) {
    "all forward candidates"
  } else if (isTRUE(o$profileInitOnStall)) {
    paste0("stalled candidates only (stallTol = ", o$stallTol, ")")
  } else {
    "off"
  }
  retries <- paste0(
    "up to ", o$maxRetries, "; maxDeltaOFV = ", o$maxDeltaOFV,
    "; OFV tolerance = ", o$retryOFVTolerance,
    "; p-value underflow ", if (isTRUE(o$retryOnUnderflow)) "retried" else "not retried",
    "; on exhaustion ", if (isTRUE(o$retryFailOnExhaustion)) "fail" else "keep best attempt"
  )
  c(
    "Search type" = o$searchType,
    "p-value fwd / bck" = paste(o$pVal$fwd, "/", o$pVal$bck),
    "Estimation" = paste0(o$estimation, "; control: ", o$control),
    "Candidates" = rel_list(o$candidates),
    "Included relations" = rel_list(o$includedRelations),
    "Shapes" = paste(o$shapes, collapse = ", "),
    "Centers" = centers,
    "Categorical cutoff" = as.character(o$catCutoff),
    "Profile warm-start" = warm_start,
    "Retries" = retries,
    "Parallel" = paste0(
      "workers = ", if (is.null(o$workers)) "sequential" else o$workers,
      "; rxThreads per worker = ", o$rxThreads
    ),
    "Save models" = paste0(
      if (isTRUE(o$saveModels)) "yes" else "no",
      if (isTRUE(o$restart)) " (restarted)" else ""
    )
  )
}

#' Lines describing where the run's files are on disk
#' @noRd
.scmFileLines <- function(outputDir) {
  if (is.null(outputDir)) {
    return("None: models and reports were not saved (saveModels = FALSE).")
  }
  if (!dir.exists(outputDir)) {
    return(paste0("Output directory: ", outputDir, " (no longer exists)"))
  }
  reports <- c("scm_log.txt", "scm_step_summary.csv", "scm_all_candidates.csv")
  reports <- reports[file.exists(file.path(outputDir, reports))]
  n_rds <- length(list.files(outputDir, pattern = "\\.rds$",
                             recursive = TRUE, ignore.case = TRUE))
  c(
    paste0("Output directory: ", outputDir),
    paste0("Report files    : ",
           if (length(reports)) paste(reports, collapse = ", ") else "none"),
    paste0("Saved models    : ", n_rds, " file", if (n_rds == 1L) "" else "s",
           " (.rds)")
  )
}
