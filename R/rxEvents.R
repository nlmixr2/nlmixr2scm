## Emitting events on the rxode2 event bus (see rxode2::rxEventListen()), so
## loggers such as nlmixr2log store one SCM result instead of one run per
## candidate model.  runSCM() enters the bus scope (candidate fits are silent;
## the workers of nlmixr2utils::.plap() enter it themselves) and on exit
## emits, in this order:
##
## * one `fitResult` (kind "scm") for the base fit with a small summary (the
##   step table, the selected model's text and the output directory), and
## * when SCM selected a different model, one `fitComplete` for it with
##   `object` = the base fit and `source = "scm"`, so it becomes its own run
##   linked to the base fit.
##
## Nothing happens when rxode2 has no event bus.

#' @noRd
.scmEventBus <- function() {
  exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE)
}

#' @noRd
.scmEventEnter <- function() {
  if (.scmEventBus()) getExportedValue("rxode2", ".rxEventEnter")()
  invisible()
}

#' Leave runSCM()'s scope and emit its result and selected model
#' @noRd
.scmEventExit <- function(result, fit, call) {
  if (!.scmEventBus()) {
    return(invisible())
  }
  .exit <- getExportedValue("rxode2", ".rxEventExit")
  if (is.null(result) || !inherits(fit, "nlmixr2FitCore")) {
    return(.exit())
  }
  .final <- .scmEventFinal(result)
  .summary <- tryCatch(.scmEventSummary(result, .final), error = function(e) NULL)
  if (is.null(.summary)) {
    return(.exit())
  }
  .exit("fitResult", fit = fit, result = .summary, kind = "scm", call = call, fun = "runSCM")
  if (inherits(.final, "nlmixr2FitCore") && !identical(.final$env, fit$env)) {
    getExportedValue("rxode2", "rxEventEmit")(
      "fitComplete", fit = .final, object = fit, call = call, objName = NULL,
      source = "scm", fun = "runSCM"
    )
  }
  invisible()
}

#' The model SCM selected (backward step, else forward step), or NULL
#' @noRd
.scmEventFinal <- function(result) {
  .pick <- function(step) if (is.list(step) && length(step)) step[[1]] else NULL
  .final <- .pick(result$resBck)
  if (is.null(.final)) .final <- .pick(result$resFwd)
  .final
}

#' A fit-free summary of an SCM run
#' @noRd
.scmEventSummary <- function(result, final) {
  .model <- if (inherits(final, "nlmixr2FitCore")) {
    tryCatch(
      paste(deparse(as.function(final$finalUi)), collapse = "\n"),
      error = function(e) NA_character_
    )
  } else {
    NA_character_
  }
  structure(
    list(
      summaryTable = as.data.frame(result$summaryTable),
      finalModel = .model,
      outputDir = result$outputDir
    ),
    class = "nlmixr2scmSummary"
  )
}
