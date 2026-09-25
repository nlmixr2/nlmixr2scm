#' Stepwise Covariate Model-selection (SCM) method
#'
#' @param fit an nlmixr2 'fit' object
#' @param varsVec character vector of PK parameters to which covariates may be
#'   added (e.g. \code{c("ka", "cl", "v")}).  Ignored when \code{pairsVec} is
#'   provided.
#' @param covarsVec character vector of continuous covariate names to test.
#'   Ignored when \code{pairsVec} is provided.
#' @param catvarsVec character vector of categorical covariate names.  Dummy
#'   columns are created internally and appended to \code{covarsVec} before
#'   building the candidate set.
#' @param pairsVec optional data frame with columns \code{var} and \code{covar}
#'   (or a list of \code{list(var=, covar=)} items) that explicitly enumerates
#'   the parameter-covariate pairs to test.  When supplied, the cartesian
#'   product of \code{varsVec} x \code{covarsVec} is bypassed.  Categorical
#'   dummy columns must already be present in the dataset; use the expanded
#'   covariate name (e.g. \code{covar = "sex_male"}).
#' @param pVal a named list with names \code{fwd} and \code{bck} for the
#'   forward and backward p-value thresholds; default
#'   \code{list(fwd = 0.05, bck = 0.01)}
#' @param inits named list mapping shape names to initial theta estimates (and
#'   optionally bounds) for the covariate parameter.  Each entry may be either
#'   a scalar estimate (e.g. \code{list(power = 0.1, lin = 0.01, cat = 0.01)})
#'   or a named list with elements \code{est}, \code{lower}, and \code{upper}
#'   (e.g. \code{list(power = list(est = 0.1, lower = -5, upper = 5),
#'   lin = 0.01)}).  Applied globally to all pairs; a per-pair \code{inits}
#'   element in a \code{pairsVec} list item takes precedence.  Shapes not
#'   listed default to \code{est = 0, lower = -Inf, upper = Inf}.  When a
#'   covariate is accepted, the estimated value from that step is automatically
#'   used as the starting estimate for all subsequent steps (warm-starting);
#'   bounds are preserved across steps.
#' @param missingToken sentinel value used to identify missing covariate
#'   observations in addition to \code{NA}.  Default \code{NA} (only true
#'   \code{NA} values are treated as missing).  Set to a numeric value (e.g.
#'   \code{-99}) or a character string (e.g. \code{"."}) to also flag that
#'   sentinel as missing.  When missing observations are detected for a
#'   covariate, the generated model expression is wrapped in an
#'   \code{ifelse()} guard so that the covariate contribution corresponds to
#'   the population-typical value (the observed median for continuous
#'   covariates, the mode for categorical).  For continuous covariates the
#'   fill is shape-aware: it is the shape expression evaluated at
#'   \code{cov = median}, which is \code{0} for \code{"power"} and
#'   \code{"lin"} (centered shapes), \code{log(median)} for \code{"log"},
#'   and \code{median} itself for \code{"identity"}.  For categorical
#'   covariates the indicator is set to the mode of the observed levels.
#' @param data data frame containing all subjects and columns required by the
#'   search (covariates, dose/observation records, etc.).  When \code{NULL}
#'   (default), \code{nlme::getData(fit)} is used, which may omit covariate
#'   columns that are not referenced by the base model.  Pass the data set
#'   explicitly whenever the base model does not reference a covariate column
#'   that will be used in the search (e.g. a categorical covariate such as
#'   \code{sex}).  No column transformations should be applied before passing
#'   the data: centering and other shape transformations are generated inside
#'   the model body by the SCM machinery.
#' @param outputDir character; path of the subdirectory used to store all
#'   fitted model \code{.rds} files.  When \code{NULL} (default), the name is
#'   derived automatically as \code{<fit_name>_scm_<N>} where \code{N} is one
#'   greater than the number of existing \code{<fit_name>_scm_*} directories in
#'   the current working directory (so repeated runs produce
#'   \code{fit_scm_1}, \code{fit_scm_2}, ...).  If the resolved directory
#'   already exists, the search stops with an error; set
#'   \code{restart = TRUE} to back it up and start fresh instead.  Ignored
#'   when \code{saveModels = FALSE}.
#' @param saveModels logical; if \code{TRUE} (default) every fitted candidate
#'   model accepted during forward or backward search is written to
#'   \code{outputDir} as an \code{.rds} file.  Set to \code{FALSE} to run the
#'   search without writing any files.
#' @param verbose logical; if \code{TRUE}, print additional detail at each
#'   search step: the full candidate table (with covariate expression and
#'   starting estimate) before fitting, the complete results table for every
#'   candidate after fitting (sorted by p-value), and--when a covariate is
#'   accepted--the accepted model's fixed-effect estimates and diagonal omega
#'   variances.  Default \code{FALSE} (only the summary line and accepted-
#'   model best row are printed, as at present).
#' @param includedRelations optional specification of covariate-parameter
#'   relationships that must be present in the model at the start of backward
#'   elimination (analogous to PsN's \code{[included_relations]}).  Can be a
#'   data frame with columns \code{var} and \code{covar} (and optionally
#'   \code{shape} or \code{shapes}), or a list of
#'   \code{list(var=, covar=, shapes=)} items in the same format as
#'   \code{pairsVec}.  Relations already in the forward-search result are left
#'   untouched; missing ones are added to the model and a single re-fit is
#'   performed before backward elimination begins.  Included relations are
#'   eligible for removal during backward elimination like any other covariate.
#'   Ignored when \code{searchType = "forward"}.
#' @param shapes character vector of shape names to test for each continuous
#'   covariate-parameter pair.  Built-in shapes are \code{"power"} (default;
#'   \code{log(cov/median)}), \code{"lin"} (\code{cov - median}), \code{"log"}
#'   (\code{log(cov)}), and \code{"identity"} (raw covariate).  Per-pair
#'   overrides can be specified via the \code{shapes} element of individual
#'   \code{pairsVec} list items.  Categorical covariates always use the
#'   \code{"cat"} shape regardless of this setting.
#' @param customShapes named list of additional shape builder functions of the
#'   form \code{function(col, center, level) -> character}.  These are merged
#'   with the built-in shapes and take precedence over them.
#' @param centers optional named numeric vector of fixed reference (centering)
#'   values for continuous covariates, e.g. \code{c(BW = 70, CrCL = 95)}.  When
#'   supplied, the named covariates are centered on these fixed values in every
#'   fit instead of the per-dataset median.  This keeps the estimated covariate
#'   coefficients (and the structural intercept) on the SAME reference across
#'   datasets and matches a data-generating model that used the same fixed
#'   references.  Names are matched against the raw covariate column; unmatched
#'   names are ignored and covariates not named here fall back to the median.
#' @param control optional control object (e.g. \code{saemControl()},
#'   \code{foceiControl()}) passed to every \code{nlmixr2()} call during the
#'   search.  When \code{NULL} (default), the control settings are inherited
#'   from the base \code{fit} object.  The \code{print} argument always
#'   overrides \code{control$print} regardless of which source is used.
#' @param print integer; how often (in iterations) the optimiser prints
#'   progress for each candidate model fit.  Default \code{100}.  Set to
#'   \code{0} to suppress all iteration output, or \code{1} to print every
#'   iteration.  Always applied as \code{control$print}, overriding any value
#'   in a user-supplied \code{control} object.
#' @param catCutoff minimum proportion of subjects that must belong to a
#'   non-reference level for it to be tested.  Levels below this threshold are
#'   lumped with the reference and excluded from the candidate set.  Applied
#'   per unique subject (not per observation row).  Default \code{0.05} (5\%).
#'   Set to \code{0} to test all non-reference levels regardless of frequency.
#' @param searchType one of \code{'scm'}, \code{'forward'}, or
#'   \code{'backward'}; default \code{'scm'}
#' @param restart logical; if \code{TRUE} the existing cache is backed up and
#'   the search restarts from scratch; default \code{FALSE}
#' @param workers integer or \code{NULL}; number of parallel workers to use
#'   when fitting candidate models.  When \code{NULL} (default), the current
#'   \code{future} plan is used unchanged.  A positive integer temporarily
#'   switches to a \code{multisession} plan with that many workers for the
#'   duration of the search.
#' @param confirm logical; if \code{TRUE} (default) and the session is
#'   interactive, the user is prompted to confirm before the search begins.
#'   Set to \code{FALSE} to skip the prompt (useful in scripts or tests).
#' @param profileInit logical; if \code{TRUE}, each forward candidate's new
#'   covariate coefficient is warm-started with a cheap 1-D FOCEi profile on a
#'   frozen base (all structural thetas fixed at their parent estimates and the
#'   between-subject variability FIXED -- not zeroed -- at its parent values)
#'   before the real estimator runs.  Keeping the random effects intact matters:
#'   profiling a covariate on a fixed-effect-only model is misspecified and can
#'   return the wrong sign.  This supplies gradient optimisers (\code{nlminb},
#'   \code{lbfgsb3c}) with a nonzero, gradient-informative starting value so
#'   they do not stall at the flat zero-effect point.  bobyqa is never involved
#'   in the caller; the profiled value is handed back to the fit's own
#'   estimator.  Default \code{FALSE}.
#' @param profileInitOnStall logical; if \code{TRUE} (default), a forward
#'   candidate whose ordinary fit \emph{stalls} -- i.e. the nested model's OFV
#'   improvement over its parent is \code{<= stallTol} (a nested model can never
#'   be genuinely worse than its parent at a true optimum, so this signals the
#'   outer optimiser never left the covariate init) -- triggers a one-shot
#'   \code{.profileCovInit()} frozen 1-D profile.  The profiled coefficient is
#'   then used as the init for a rescue refit.  Unlike \code{profileInit}, this
#'   fires only when a stall is detected, so healthy candidates (e.g. analytic
#'   \code{linCmt()} fits) pay no extra cost.  It is the fix for ODE covariate
#'   candidates that stall at their init under solver-noise-flattened outer
#'   objectives.  Default \code{TRUE}.
#' @param stallTol numeric; OFV-improvement threshold below which a forward
#'   candidate is considered stalled and eligible for the profile rescue.
#'   Default \code{0} (any nested model no better than its parent).
#' @param maxRetries integer; maximum number of retry attempts per candidate
#'   when the OFV is deemed unrealistic.  Default \code{3L}.  Set to \code{0}
#'   to disable the retry mechanism entirely.
#' @param maxDeltaOFV numeric; absolute ceiling on plausible OFV change.  A
#'   candidate whose \code{|dOFV|} exceeds this value is flagged as unrealistic
#'   and retried.  Default \code{Inf} (disabled; only criterion 1 and 2 apply).
#' @param retryPerturbSD numeric; standard deviation of the log-normal
#'   perturbation applied to \code{cov_init} on odd-numbered retries.  Default
#'   \code{0.5}.
#' @param retrySmallInit numeric; covariate theta init used on even-numbered
#'   retries as a near-zero safe-harbour.  Default \code{0.01}.
#' @param retryOFVTolerance numeric or \code{NULL}; the margin (in OFV units)
#'   by which the candidate OFV must exceed the parent before a retry is
#'   triggered (criterion 1).  \code{NULL} (default) auto-detects: \code{10}
#'   for SAEM (stochastic noise), \code{0} for all other estimators.
#' @param retryFailOnExhaustion logical; when \code{TRUE} and all retry
#'   attempts are exhausted with an unrealistic OFV, the candidate is marked
#'   as failed and excluded from the search.  When \code{FALSE} (default) a
#'   warning is emitted and the best available result is accepted.
#' @param rxThreads integer, \code{"auto"}, or \code{NULL}; number of rxode2
#'   OpenMP threads to use per worker while fitting candidate models.  When
#'   \code{NULL} (default), the current \code{rxode2::getRxThreads()} value
#'   is used and applied consistently to every worker.  \code{"auto"}
#'   divides the total core count evenly across the effective number of
#'   workers.  Whenever more than one worker is involved, \code{workers *
#'   rxThreads} must not exceed the machine's core count -- \code{runSCM()}
#'   aborts with an explanatory error before any fitting starts if it
#'   would, since each parallel worker is a separate process that runs its
#'   own independent rxode2 thread pool.  A single worker is never subject
#'   to this check.  Note this means an existing \code{workers > 1} or
#'   \code{workers = "auto"} call that does not also set \code{rxThreads}
#'   may now abort where it previously ran silently oversubscribed; set
#'   \code{rxThreads} explicitly (e.g. \code{rxThreads = 1}) to restore the
#'   prior behavior.
#'
#' @return A list with elements \code{summaryTable} (combined forward and
#'   backward results), \code{resFwd} (list of final fit and step table from
#'   forward search), and \code{resBck} (same for backward search).
#'
#' @export
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani, Justin Wilkins
#'
#' @examplesIf requireNamespace("nlmixr2data", quietly = TRUE)
#' \donttest{
#' # Fitting the base model and running even a one-pair search takes longer
#' # than a few seconds, so this example is wrapped in \donttest{}.
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     label("Ka")
#'     tcl <- log(c(0, 2.7, 100))
#'     label("Cl")
#'     tv <- 3.45
#'     label("V")
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' fit <- nlmixr2est::nlmixr2(
#'   one.cmt, nlmixr2data::theo_sd,
#'   est = "focei", control = nlmixr2est::foceiControl(print = 0)
#' )
#'
#' # Test one parameter-covariate pair.  saveModels = FALSE keeps the search
#' # from writing anything to disk.
#' res <- runSCM(fit,
#'   pairsVec = list(list(var = "cl", covar = "WT", shapes = "power")),
#'   searchType = "forward",
#'   saveModels = FALSE,
#'   confirm = FALSE,
#'   print = 0,
#'   workers = 1L,
#'   rxThreads = 2L
#' )
#' res$summaryTable
#' }
runSCM <- function(
  fit,
  varsVec = NULL,
  covarsVec = NULL,
  pVal = list(fwd = 0.05, bck = 0.01),
  catvarsVec = NULL,
  pairsVec = NULL,
  shapes = "power",
  customShapes = NULL,
  centers = NULL,
  inits = list(),
  missingToken = NA,
  includedRelations = NULL,
  data = NULL,
  catCutoff = 0.05,
  outputDir = NULL,
  saveModels = TRUE,
  verbose = FALSE,
  control = NULL,
  print = 100,
  searchType = c("scm", "forward", "backward"),
  restart = FALSE,
  workers = NULL,
  confirm = TRUE,
  profileInit = FALSE,
  profileInitOnStall = TRUE,
  stallTol = 0,
  maxRetries = 3L,
  maxDeltaOFV = Inf,
  retryPerturbSD = 0.5,
  retrySmallInit = 0.01,
  retryOFVTolerance = NULL,
  retryFailOnExhaustion = FALSE,
  rxThreads = NULL
) {
  if (!is.numeric(stats::AIC(fit))) {
    cli::cli_alert_danger(
      "the 'fit' object needs to have an objective function value associated with it"
    )
    cli::cli_alert_info(
      "try computing 'AIC(fitobject)' in console to compute and store the corresponding OBJF value"
    )
    stop(
      "aborting...objf value not associated with the current 'fit' object",
      call. = FALSE
    )
  }

  if (is.null(data)) {
    data <- nlme::getData(fit)
  }

  # Build categorical indicator columns and trim dataset to needed columns.
  # Must happen before buildPairs/enrichPairs so the indicator columns exist
  # in data when the model is fitted.
  cat_levels <- NULL
  cat_ref <- NULL
  cat_dropped <- NULL
  if (!is.null(catvarsVec) && length(catvarsVec) > 0) {
    scm_prep <- .makeSCMData(
      data,
      catvarsVec,
      fit = fit,
      catCutoff = catCutoff,
      covarsVec = covarsVec
    )
    data <- scm_prep$data
    cat_levels <- scm_prep$catLevels
    cat_ref <- scm_prep$catRef
    cat_dropped <- scm_prep$catDropped
  }

  # When using varsVec/covarsVec, merge catvarsVec into covarsVec so that
  # buildPairs creates (var, catcovar) base rows that .enrichPairs then
  # expands into per-level rows.  pairsVec users must include cat covariates
  # explicitly (either pre-expanded with level info, or as raw column names
  # paired with catvarsVec for auto-expansion).
  if (!is.null(varsVec) && is.null(pairsVec)) {
    covarsVec <- c(covarsVec, catvarsVec)
  }
  pairs <- buildPairs(
    varsVec = varsVec,
    covarsVec = covarsVec,
    pairsVec = pairsVec
  )
  pairs <- .enrichPairs(
    pairs,
    data,
    catvarsVec = catvarsVec,
    missingToken = missingToken,
    catLevels = cat_levels
  )
  pairs <- .applyFixedCenters(pairs, centers)
  pairs <- .expandShapes(
    pairs,
    shapes = shapes,
    customShapes = customShapes,
    inits = inits
  )

  # Process includedRelations through the same pipeline so it carries covExpr
  included_pairs <- NULL
  if (!is.null(includedRelations)) {
    included_pairs <- buildPairs(pairsVec = includedRelations)
    included_pairs <- .enrichPairs(
      included_pairs,
      data,
      catvarsVec = catvarsVec,
      missingToken = missingToken,
      catLevels = cat_levels
    )
    included_pairs <- .applyFixedCenters(included_pairs, centers)
    included_pairs <- .expandShapes(
      included_pairs,
      shapes = shapes,
      customShapes = customShapes,
      inits = inits
    )
  }

  searchType <- match.arg(searchType)

  if (!all(names(pVal) %in% c("fwd", "bck"))) {
    stop("pVal must be a list with names 'fwd' and 'bck'")
  }

  # -- Resolve output dir name before summary (no dir creation yet) ----------
  if (saveModels && is.null(outputDir)) {
    .fit_nm <- gsub("[^A-Za-z0-9_.]", "_", as.character(substitute(fit)))
    .n_exist <- sum(grepl(
      paste0("^", .fit_nm, "_scm_[0-9]+$"),
      list.dirs(".", full.names = FALSE, recursive = FALSE)
    ))
    outputDir <- paste0(.fit_nm, "_scm_", .n_exist + 1L)
  }
  if (!is.null(outputDir)) {
    outputDir <- normalizePath(outputDir, mustWork = FALSE)
  }

  # -- Pre-SCM summary and confirmation --------------------------------------
  n_pairs <- nrow(pairs)
  base_objf <- tryCatch(round(fit$objf, 3), error = function(e) NA) # nolint: object_usage_linter.
  base_n_par <- tryCatch(
    # nolint: object_usage_linter.
    length(fit$finalUiEnv$ini$est),
    error = function(e) NA
  )
  out_lbl <- if (saveModels) outputDir else "(none \u2014 saveModels = FALSE)" # nolint: object_usage_linter.

  ctrl_lbl <- if (!is.null(control)) {
    # nolint: object_usage_linter.
    paste0(class(control)[1], " (user-supplied)")
  } else {
    paste0(class(fit$control)[1], " (inherited from fit)")
  }

  cli::cli_rule(left = "SCM Summary")
  cli::cli_inform(c(
    "i" = "Estimation method : {fit$est}",
    "i" = "Control           : {ctrl_lbl}",
    "i" = "Base model OFV    : {base_objf}",
    "i" = "Base model params : {base_n_par}",
    "i" = "Search type       : {searchType}",
    "i" = "p-value fwd / bck : {pVal$fwd} / {pVal$bck}",
    "i" = "Output folder     : {out_lbl}",
    "i" = "Parameters        : {length(unique(pairs$var))}",
    "  " = "({paste(unique(pairs$var), collapse = ', ')})",
    "i" = "Covariates        : {length(unique(pairs$covar))}",
    "  " = "({paste(unique(pairs$covar), collapse = ', ')})",
    "i" = "Total candidates  : {n_pairs}"
  ))

  # Categorical covariate breakdown
  if (!is.null(cat_levels) && length(cat_levels) > 0) {
    cli::cli_rule(left = "Categorical covariates")
    for (cat_col in names(cat_levels)) {
      cat_ref_lbl <- if (!is.null(cat_ref[[cat_col]])) {
        cat_ref[[cat_col]]
      } else {
        "?"
      } # nolint: object_usage_linter.
      cat_lvls <- cat_levels[[cat_col]]
      cat_drop <- if (!is.null(cat_dropped[[cat_col]])) {
        cat_dropped[[cat_col]]
      } else {
        character(0)
      }
      cat_ind <- paste0(cat_col, "_", cat_lvls) # nolint: object_usage_linter.
      cli::cli_inform(c(
        "i" = "{cat_col}: reference = '{cat_ref_lbl}'"
      ))
      if (length(cat_lvls) > 0) {
        cli::cli_inform(
          "  Indicators ({length(cat_lvls)}): {paste(cat_ind, collapse = ', ')}"
        )
      } else {
        cli::cli_inform("  No qualifying levels to test")
      }
      if (length(cat_drop) > 0) {
        cli::cli_inform(
          "  Dropped (<{catCutoff*100}%): {paste(cat_drop, collapse = ', ')}"
        )
      }
    }
  }

  cli::cli_rule(left = "Relationships to test")
  for (pair_j in seq_len(n_pairs)) {
    pair_sl <- if (
      # nolint: object_usage_linter.
      "shape" %in%
        names(pairs) &&
        !is.na(pairs$shape[pair_j]) &&
        nzchar(pairs$shape[pair_j])
    ) {
      paste0(" [", pairs$shape[pair_j], "]")
    } else {
      ""
    }
    cli::cli_inform(
      "  {pair_j}. {pairs$covar[pair_j]} ~ {pairs$var[pair_j]}{pair_sl}"
    )
  }
  cli::cli_rule()

  if (isTRUE(confirm) && interactive()) {
    .ans <- readline("Proceed with SCM? [y/N]: ")
    if (!tolower(trimws(.ans)) %in% c("y", "yes")) {
      cli::cli_alert_warning("SCM cancelled.")
      return(invisible(NULL))
    }
  }

  # -- Create output directory (name already resolved above) -----------------
  if (saveModels) {
    if (!dir.exists(outputDir)) {
      dir.create(outputDir)
    } else if (isTRUE(restart)) {
      backupDir <- paste0(
        outputDir,
        "_backup_",
        format(Sys.time(), "%Y%m%d_%H%M%S")
      )
      file.rename(outputDir, backupDir)
      cli::cli_alert_info("Prior cache moved to {backupDir}")
      dir.create(outputDir)
    } else {
      stop(
        "Output directory '",
        outputDir,
        "' already exists. ",
        "Use restart = TRUE to back it up and start fresh, or set ",
        "saveModels = FALSE to run without saving.",
        call. = FALSE
      )
    }
    cli::cli_alert_info("Saving models to: {outputDir}")
  } else {
    # saveModels = FALSE: no directory needed; use cwd as placeholder for forwardSearch/backwardSearch
    if (is.null(outputDir)) outputDir <- getwd()
  }

  # The placeholder above is the user's working directory, which must never be
  # written to.  Report files are emitted only when saveModels = TRUE.
  reportDir <- if (saveModels) outputDir else NULL

  effective_rx_threads <- nlmixr2utils::resolveRxThreads(workers, rxThreads)
  nlmixr2utils::.withWorkerPlan(workers, rxThreads = effective_rx_threads, {
    # nolint: object_usage_linter.
    if (searchType == "scm") {
      resFwd <- forwardSearch(
        pairs,
        fit,
        data,
        pVal = pVal$fwd,
        outputDir = outputDir,
        restart = restart,
        saveModels = saveModels,
        verbose = verbose,
        control = control,
        print = print,
        profileInit = profileInit,
        profileInitOnStall = profileInitOnStall,
        stallTol = stallTol,
        maxRetries = maxRetries,
        maxDeltaOFV = maxDeltaOFV,
        retryPerturbSD = retryPerturbSD,
        retrySmallInit = retrySmallInit,
        retryOFVTolerance = retryOFVTolerance,
        retryFailOnExhaustion = retryFailOnExhaustion,
        rxThreads = effective_rx_threads
      )
      resBck <- backwardSearch(
        pairs,
        fitorig = fit,
        fitupdated = resFwd[[1]],
        data = data,
        pVal = pVal$bck,
        reFitCovars = FALSE,
        context_pairs = resFwd[[3]],
        includedRelations = included_pairs,
        outputDir = outputDir,
        restart = restart,
        saveModels = saveModels,
        verbose = verbose,
        control = control,
        print = print,
        maxRetries = maxRetries,
        maxDeltaOFV = maxDeltaOFV,
        retryPerturbSD = retryPerturbSD,
        retrySmallInit = retrySmallInit,
        retryOFVTolerance = retryOFVTolerance,
        retryFailOnExhaustion = retryFailOnExhaustion,
        rxThreads = effective_rx_threads
      )
      data <- NULL # release temporary SCM dataset
      summaryTable <- Reduce(rbind, list(resFwd[[2]], resBck[[2]]))
      .printFinalSCMSummary(
        summaryTable,
        finalFit = resBck[[1]],
        searchType = searchType,
        outputDir = reportDir,
        metadata = list(
          baseFit = fit,
          pVal = pVal,
          catCutoff = catCutoff
        )
      )
      list(summaryTable = summaryTable, resFwd = resFwd, resBck = resBck)
    } else if (searchType == "forward") {
      resFwd <- forwardSearch(
        pairs,
        fit,
        data,
        pVal = pVal$fwd,
        outputDir = outputDir,
        restart = restart,
        saveModels = saveModels,
        verbose = verbose,
        control = control,
        print = print,
        profileInit = profileInit,
        profileInitOnStall = profileInitOnStall,
        stallTol = stallTol,
        maxRetries = maxRetries,
        maxDeltaOFV = maxDeltaOFV,
        retryPerturbSD = retryPerturbSD,
        retrySmallInit = retrySmallInit,
        retryOFVTolerance = retryOFVTolerance,
        retryFailOnExhaustion = retryFailOnExhaustion,
        rxThreads = effective_rx_threads
      )
      data <- NULL # release temporary SCM dataset
      .printFinalSCMSummary(
        resFwd[[2]],
        finalFit = resFwd[[1]],
        searchType = searchType,
        outputDir = reportDir,
        metadata = list(
          baseFit = fit,
          pVal = pVal,
          catCutoff = catCutoff
        )
      )
      list(summaryTable = resFwd[[2]], resFwd = resFwd, resBck = NULL)
    } else {
      resBck <- backwardSearch(
        pairs,
        fitorig = fit,
        data = data,
        pVal = pVal$bck,
        reFitCovars = TRUE,
        includedRelations = included_pairs,
        outputDir = outputDir,
        restart = restart,
        saveModels = saveModels,
        verbose = verbose,
        control = control,
        print = print,
        maxRetries = maxRetries,
        maxDeltaOFV = maxDeltaOFV,
        retryPerturbSD = retryPerturbSD,
        retrySmallInit = retrySmallInit,
        retryOFVTolerance = retryOFVTolerance,
        retryFailOnExhaustion = retryFailOnExhaustion,
        rxThreads = effective_rx_threads
      )
      data <- NULL # release temporary SCM dataset
      .printFinalSCMSummary(
        resBck[[2]],
        finalFit = resBck[[1]],
        searchType = searchType,
        outputDir = reportDir,
        metadata = list(
          baseFit = fit,
          pVal = pVal,
          catCutoff = catCutoff
        )
      )
      list(summaryTable = resBck[[2]], resFwd = NULL, resBck = resBck)
    }
  })
}

# -- Helpers -------------------------------------------------------------------

#' Print a formatted end-of-SCM summary table and write log files
#'
#' Prints to the console:
#' \enumerate{
#'   \item A step table (one row per step -- best candidate selected).
#'   \item A full candidate table (every model tested).
#'   \item A final-model conclusion.
#' }
#' When \code{outputDir} is provided, writes three files:
#' \describe{
#'   \item{\code{scm_step_summary.csv}}{Step table as CSV.}
#'   \item{\code{scm_all_candidates.csv}}{All tested models as CSV.}
#'   \item{\code{scm_log.txt}}{Human-readable log of the full run.}
#' }
#'
#' @param summaryTable  combined resTableComplete from forward + backward search
#' @param finalFit      the nlmixr2 fit object at the end of the search
#' @param searchType    one of \code{"scm"}, \code{"forward"}, \code{"backward"}
#' @param outputDir     path to output directory; \code{NULL} skips file writing
#' @param metadata      optional list with elements \code{baseFit}, \code{pVal},
#'   and \code{catCutoff} used to enrich the log header
#' @noRd
.printFinalSCMSummary <- function(
  summaryTable,
  finalFit,
  searchType,
  outputDir = NULL,
  metadata = NULL
) {
  if (is.null(summaryTable) || nrow(summaryTable) == 0L) {
    cli::cli_inform("No SCM steps were recorded.")
    return(invisible(NULL))
  }

  # -- Helper: build relation label ------------------------------------------
  .rel_lbl <- function(covar, var, shape) {
    lbl <- paste0(covar, "~", var)
    if (!is.na(shape) && nzchar(shape)) {
      lbl <- paste0(lbl, " [", shape, "]")
    }
    lbl
  }

  .ref_objf <- function(st_vec, objf_vec, delta_vec) {
    objf_vec - delta_vec
  }

  # -- Step summary (best per step) ------------------------------------------
  step_keys <- paste0(summaryTable$searchType, "_", summaryTable$step)
  best_rows <- do.call(
    rbind,
    lapply(unique(step_keys), function(k) {
      grp <- summaryTable[step_keys == k, , drop = FALSE]
      acc <- grp[grp$included %in% c("yes", "dropped"), , drop = FALSE]
      if (nrow(acc) > 0L) {
        return(acc[1L, , drop = FALSE])
      }
      # Exit step: for forward show most-significant candidate (lowest p);
      # for backward show least-significant (highest p -- last one that couldn't be dropped).
      if (grp$searchType[1L] == "backward") {
        grp[which.max(grp$pchisqr), , drop = FALSE]
      } else {
        grp[which.min(grp$pchisqr), , drop = FALSE]
      }
    })
  )
  fwd_rows <- best_rows[best_rows$searchType == "forward", , drop = FALSE]
  bck_rows <- best_rows[best_rows$searchType == "backward", , drop = FALSE]
  if (nrow(fwd_rows) > 0L) {
    fwd_rows <- fwd_rows[order(fwd_rows$step), , drop = FALSE]
  }
  if (nrow(bck_rows) > 0L) {
    bck_rows <- bck_rows[order(bck_rows$step), , drop = FALSE]
  }
  best_rows <- rbind(fwd_rows, bck_rows)

  step_ref_objf <- .ref_objf(
    best_rows$searchType,
    best_rows$objf,
    best_rows$deltObjf
  )
  step_rel <- mapply(
    .rel_lbl,
    best_rows$covar,
    best_rows$var,
    best_rows$shape
  )
  step_dir <- ifelse(best_rows$searchType == "forward", "Forward", "Backward")
  step_decision <- ifelse(
    best_rows$included %in% c("yes", "dropped"),
    ifelse(best_rows$searchType == "forward", "Added", "Removed"),
    ifelse(best_rows$searchType == "backward", "Retained", "Not selected")
  )

  # -- All-candidates table ---------------------------------------------------
  all_sorted <- summaryTable[
    order(
      match(summaryTable$searchType, c("forward", "backward")),
      summaryTable$step,
      summaryTable$pchisqr
    ),
  ]
  all_ref_objf <- .ref_objf(
    all_sorted$searchType,
    all_sorted$objf,
    all_sorted$deltObjf
  )
  all_rel <- mapply(
    .rel_lbl,
    all_sorted$covar,
    all_sorted$var,
    all_sorted$shape
  )
  all_dir <- ifelse(all_sorted$searchType == "forward", "Forward", "Backward")
  all_decision <- ifelse(
    all_sorted$included %in% c("yes", "dropped"),
    ifelse(all_sorted$searchType == "forward", "Added", "Removed"),
    ifelse(all_sorted$searchType == "backward", "Retained", "Not selected")
  )

  # -- Final model sets ------------------------------------------------------
  fwd_acc <- summaryTable[
    summaryTable$searchType == "forward" & summaryTable$included == "yes",
    ,
    drop = FALSE
  ]
  bck_rem <- summaryTable[
    summaryTable$searchType == "backward" & summaryTable$included == "dropped",
    ,
    drop = FALSE
  ]
  bck_keys <- paste0(bck_rem$covar, "_", bck_rem$var)
  fwd_keys <- paste0(fwd_acc$covar, "_", fwd_acc$var)

  if (searchType == "backward") {
    final_acc <- NULL
  } else {
    final_acc <- fwd_acc[!fwd_keys %in% bck_keys, , drop = FALSE]
  }

  final_objf <- tryCatch(round(finalFit$objf, 3), error = function(e) NA)

  # -- Console: step summary -------------------------------------------------
  rel_w <- max(nchar(c("Relation", step_rel)))
  hdr <- sprintf(
    "%-9s  %-4s  %-*s  %10s  %10s  %8s  %9s  %s",
    "Direction",
    "Step",
    rel_w,
    "Relation",
    "Ref OFV",
    "OFV",
    "dOFV",
    "p-value",
    "Decision"
  )
  cli::cli_rule(left = "SCM Step Summary")
  cat(hdr, "\n")
  cat(strrep("-", nchar(hdr)), "\n")
  for (i in seq_len(nrow(best_rows))) {
    cat(sprintf(
      "%-9s  %-4d  %-*s  %10.3f  %10.3f  %8.3f  %9.4f  %s\n",
      step_dir[i],
      best_rows$step[i],
      rel_w,
      step_rel[i],
      step_ref_objf[i],
      best_rows$objf[i],
      best_rows$deltObjf[i],
      best_rows$pchisqr[i],
      step_decision[i]
    ))
  }

  # -- Console: all candidates -----------------------------------------------
  rel_w2 <- max(nchar(c("Relation", all_rel)))
  hdr2 <- sprintf(
    "%-9s  %-4s  %-*s  %10s  %10s  %8s  %9s  %s",
    "Direction",
    "Step",
    rel_w2,
    "Relation",
    "Ref OFV",
    "OFV",
    "dOFV",
    "p-value",
    "Decision"
  )
  cli::cli_rule(left = "SCM All Candidates")
  cat(hdr2, "\n")
  cat(strrep("-", nchar(hdr2)), "\n")
  prev_key <- ""
  for (i in seq_len(nrow(all_sorted))) {
    cur_key <- paste0(all_dir[i], "_", all_sorted$step[i])
    if (cur_key != prev_key && i > 1L) {
      cat("\n")
    }
    prev_key <- cur_key
    cat(sprintf(
      "%-9s  %-4d  %-*s  %10.3f  %10.3f  %8.3f  %9.4f  %s\n",
      all_dir[i],
      all_sorted$step[i],
      rel_w2,
      all_rel[i],
      all_ref_objf[i],
      all_sorted$objf[i],
      all_sorted$deltObjf[i],
      all_sorted$pchisqr[i],
      all_decision[i]
    ))
  }

  # -- Console: final model --------------------------------------------------
  cli::cli_rule(left = "Final model")
  if (searchType == "backward") {
    if (nrow(bck_rem) == 0L) {
      cli::cli_inform("  No covariates were removed.")
    } else {
      cli::cli_inform("  Removed:")
      for (j in seq_len(nrow(bck_rem))) {
        rel_lbl_j <- .rel_lbl(
          bck_rem$covar[j],
          bck_rem$var[j],
          bck_rem$shape[j]
        ) # nolint: object_usage_linter.
        cli::cli_inform("    {rel_lbl_j}")
      }
    }
  } else {
    if (is.null(final_acc) || nrow(final_acc) == 0L) {
      cli::cli_inform("  No covariates were retained in the final model.")
    } else {
      cli::cli_inform("  Retained:")
      for (j in seq_len(nrow(final_acc))) {
        rel_lbl_j <- .rel_lbl(
          final_acc$covar[j],
          final_acc$var[j],
          final_acc$shape[j]
        )
        cli::cli_inform("    {rel_lbl_j}")
      }
    }
  }
  cli::cli_inform(c("i" = "Final model OFV: {final_objf}"))
  cli::cli_rule()

  # -- File output -----------------------------------------------------------
  if (!is.null(outputDir) && dir.exists(outputDir)) {
    # -- scm_step_summary.csv ----------------------------------------------
    step_csv <- data.frame(
      Direction = step_dir,
      Step = best_rows$step,
      Relation = step_rel,
      Ref_OFV = round(step_ref_objf, 3),
      OFV = round(best_rows$objf, 3),
      dOFV = round(best_rows$deltObjf, 3),
      p_value = round(best_rows$pchisqr, 4),
      Decision = step_decision,
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      step_csv,
      file.path(outputDir, "scm_step_summary.csv"),
      row.names = FALSE
    )

    # -- scm_all_candidates.csv --------------------------------------------
    all_csv <- data.frame(
      Direction = all_dir,
      Step = all_sorted$step,
      Relation = all_rel,
      Ref_OFV = round(all_ref_objf, 3),
      OFV = round(all_sorted$objf, 3),
      dOFV = round(all_sorted$deltObjf, 3),
      p_value = round(all_sorted$pchisqr, 4),
      AIC = round(all_sorted$AIC, 3),
      BIC = round(all_sorted$BIC, 3),
      numParams = all_sorted$numParams,
      covarEffect = round(all_sorted$covarEffect, 4),
      bsvReduction = round(all_sorted$bsvReduction, 4),
      Decision = all_decision,
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      all_csv,
      file.path(outputDir, "scm_all_candidates.csv"),
      row.names = FALSE
    )

    # -- scm_log.txt -------------------------------------------------------
    log <- character(0)
    .ln <- function(...) {
      log <<- c(log, paste0(...))
    }
    .sep <- function(ch = "-", w = 70) .ln(strrep(ch, w))

    .sep("=")
    .ln("SCM LOG")
    .sep("=")
    .ln("Generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    .ln("Search type: ", searchType)

    if (!is.null(metadata)) {
      bf <- metadata$baseFit
      .ln(
        "Estimation : ",
        tryCatch(bf$est, error = function(e) "unknown")
      )
      .ln(
        "Base OFV   : ",
        tryCatch(round(bf$objf, 3), error = function(e) "unknown")
      )
      .ln(
        "Base params: ",
        tryCatch(length(bf$finalUiEnv$ini$est), error = function(e) "unknown")
      )
      pv <- metadata$pVal
      if (!is.null(pv)) {
        if (!is.null(pv$fwd)) {
          .ln("p-val fwd  : ", pv$fwd)
        }
        if (!is.null(pv$bck)) .ln("p-val bck  : ", pv$bck)
      }
    }
    .ln("Output dir : ", outputDir)
    .ln("")

    # Unique relationships actually tested
    tested_rels <- unique(paste0(
      summaryTable$covar,
      "_",
      summaryTable$var,
      "_",
      summaryTable$shape
    ))
    .ln("Relationships tested (", length(tested_rels), "):")
    uniq_rows <- summaryTable[
      !duplicated(
        paste0(summaryTable$covar, summaryTable$var, summaryTable$shape)
      ),
    ]
    for (k in seq_len(nrow(uniq_rows))) {
      .ln(
        "  ",
        k,
        ". ",
        .rel_lbl(uniq_rows$covar[k], uniq_rows$var[k], uniq_rows$shape[k])
      )
    }

    # Step-by-step detail
    for (phase in c("forward", "backward")) {
      phase_rows <- summaryTable[
        summaryTable$searchType == phase,
        ,
        drop = FALSE
      ]
      if (nrow(phase_rows) == 0L) {
        next
      }
      .ln("")
      .sep()
      .ln(toupper(phase), " SEARCH")
      .sep()
      pv_thresh <- if (!is.null(metadata$pVal)) {
        if (phase == "forward") metadata$pVal$fwd else metadata$pVal$bck
      } else {
        NA
      }

      for (s in sort(unique(phase_rows$step))) {
        grp <- phase_rows[phase_rows$step == s, , drop = FALSE]
        grp <- grp[order(grp$pchisqr), , drop = FALSE]
        grp_ref <- .ref_objf(grp$searchType, grp$objf, grp$deltObjf)
        accepted <- any(grp$included %in% c("yes", "dropped"))
        .ln("")
        .ln(
          "Step ",
          s,
          " (",
          nrow(grp),
          " candidate",
          if (nrow(grp) != 1L) "s" else "",
          "):"
        )
        for (r in seq_len(nrow(grp))) {
          tag <- if (grp$included[r] %in% c("yes", "dropped")) {
            if (phase == "forward") " [ADDED]" else " [REMOVED]"
          } else {
            ""
          }
          .ln(sprintf(
            "  %-30s  Ref OFV=%10.3f  OFV=%10.3f  dOFV=%8.3f  p=%8.4f%s",
            .rel_lbl(grp$covar[r], grp$var[r], grp$shape[r]),
            grp_ref[r],
            grp$objf[r],
            grp$deltObjf[r],
            grp$pchisqr[r],
            tag
          ))
        }
        if (!accepted) {
          thr_str <- if (!is.na(pv_thresh)) {
            paste0(" (p \u2264 ", pv_thresh, ")")
          } else {
            ""
          }
          .ln(
            "  \u2192 No candidate met threshold",
            thr_str,
            ". ",
            tools::toTitleCase(phase),
            " search complete."
          )
        }
      }
    }

    # Final conclusion
    .ln("")
    .sep("=")
    .ln("FINAL MODEL")
    .sep("=")
    if (searchType == "backward") {
      if (nrow(bck_rem) == 0L) {
        .ln("No covariates were removed.")
      } else {
        .ln("Removed:")
        for (j in seq_len(nrow(bck_rem))) {
          .ln(
            "  ",
            .rel_lbl(bck_rem$covar[j], bck_rem$var[j], bck_rem$shape[j])
          )
        }
      }
    } else {
      if (is.null(final_acc) || nrow(final_acc) == 0L) {
        .ln("No covariates were retained.")
      } else {
        .ln("Retained:")
        for (j in seq_len(nrow(final_acc))) {
          .ln(
            "  ",
            .rel_lbl(final_acc$covar[j], final_acc$var[j], final_acc$shape[j])
          )
        }
      }
    }
    .ln("Final OFV: ", final_objf)

    writeLines(log, file.path(outputDir, "scm_log.txt"))
    cli::cli_alert_success(
      "Log files written to {outputDir}: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv"
    )
  }

  invisible(NULL)
}

#' Build a pairs data frame from varsVec/covarsVec or an explicit pairsVec
#'
#' @param varsVec character vector of PK parameter names
#' @param covarsVec character vector of covariate names
#' @param pairsVec optional data frame with columns \code{var} and \code{covar},
#'   or a list of \code{list(var=, covar=)} items.  May additionally carry a
#'   \code{center} column / element to fix the reference (centering) value for a
#'   continuous covariate; when absent the per-dataset median is used by
#'   \code{.enrichPairs()}.
#' @return data frame with columns \code{var} and \code{covar}, one row per pair
#' @noRd
buildPairs <- function(varsVec = NULL, covarsVec = NULL, pairsVec = NULL) {
  if (!is.null(pairsVec)) {
    if (is.data.frame(pairsVec)) {
      stopifnot(all(c("var", "covar") %in% names(pairsVec)))
      keep <- intersect(c("var", "covar", "shapes", "inits", "center"),
                        names(pairsVec))
      df <- pairsVec[, keep, drop = FALSE]
      # If caller used a 'shape' character column instead of a 'shapes' list-col,
      # convert it so .expandShapes() can dispatch per-row shapes correctly.
      if ("shape" %in% names(pairsVec) && !"shapes" %in% names(df)) {
        df$shapes <- as.list(pairsVec$shape)
      }
      return(df)
    }
    # Build column vectors manually to support optional list-columns
    n <- length(pairsVec)
    var_vec <- character(n)
    covar_vec <- character(n)
    shapes_list <- vector("list", n)
    inits_list <- vector("list", n)
    center_vec <- rep(NA_real_, n)
    has_shapes <- FALSE
    has_inits <- FALSE
    has_center <- FALSE
    for (k in seq_len(n)) {
      p <- pairsVec[[k]]
      var_vec[k] <- if (is.list(p)) p$var else p[[1]]
      covar_vec[k] <- if (is.list(p)) p$covar else p[[2]]
      if (is.list(p) && !is.null(p$shapes)) {
        shapes_list[[k]] <- p$shapes
        has_shapes <- TRUE
      }
      if (is.list(p) && !is.null(p$inits)) {
        inits_list[[k]] <- p$inits # named list, e.g. list(power=0.75, lin=0.01)
        has_inits <- TRUE
      }
      if (is.list(p) && !is.null(p$center)) {
        center_vec[k] <- as.numeric(p$center)
        has_center <- TRUE
      }
    }
    df <- data.frame(var = var_vec, covar = covar_vec, stringsAsFactors = FALSE)
    if (has_shapes) {
      df$shapes <- shapes_list
    }
    if (has_inits) {
      df$inits <- inits_list
    }
    if (has_center) {
      df$center <- center_vec
    }
    return(df)
  }
  if (is.null(varsVec) || is.null(covarsVec)) {
    stop("provide either pairsVec or both varsVec and covarsVec")
  }
  expand.grid(var = varsVec, covar = covarsVec, stringsAsFactors = FALSE)
}

#' Build a minimal SCM dataset with categorical indicator columns
#'
#' For each column in \code{catvarsVec}:
#' \enumerate{
#'   \item The reference level is the most frequent level assessed
#'     \emph{per unique subject} (not per observation row).
#'   \item Non-reference levels whose subject-proportion is below
#'     \code{catCutoff} are dropped and lumped with the reference.
#'   \item A 0/1 indicator column is created for each qualifying level,
#'     named \code{<col>_<level>} (e.g. \code{race_black}).
#' }
#' The returned dataset is trimmed to columns needed for fitting: those present
#' in the base model's data, plus the covariate columns under test, plus the
#' new indicator columns.  The original categorical columns are retained so
#' that \code{.enrichPairs()} can inspect them if needed.
#'
#' @param data       data frame containing all subjects and covariates
#' @param catvarsVec character vector of categorical covariate column names
#' @param fit        nlmixr2 fit object (used to identify columns needed by
#'   the base model via \code{nlme::getData(fit)})
#' @param catCutoff  minimum subject proportion for a non-reference level to be
#'   tested; levels below this are lumped with the reference.  Default
#'   \code{0.05}.
#' @param covarsVec  optional character vector of continuous covariate column
#'   names to retain in the trimmed dataset
#' @return a list with elements:
#'   \describe{
#'     \item{\code{data}}{trimmed dataset with indicator columns added}
#'     \item{\code{catLevels}}{named list -- one entry per \code{catvarsVec}
#'       column, containing the qualifying non-reference level names in
#'       frequency-descending order}
#'   }
#' @noRd
.makeSCMData <- function(
  data,
  catvarsVec,
  fit,
  catCutoff = 0.05,
  covarsVec = NULL
) {
  id_col <- .idColumn(data)
  cat_levels <- list()
  cat_ref <- list()
  cat_dropped <- list()

  for (col in catvarsVec) {
    if (!col %in% names(data)) {
      stop(
        "Categorical covariate '",
        col,
        "' not found in data.",
        call. = FALSE
      )
    }

    # Frequency computed per subject, not per observation row
    subj_vals <- as.character(
      data[[col]][!duplicated(data[[id_col]])]
    )
    subj_vals <- subj_vals[!is.na(subj_vals)]

    if (length(subj_vals) == 0L) {
      cat_levels[[col]] <- character(0L)
      cat_ref[[col]] <- NA_character_
      cat_dropped[[col]] <- character(0L)
      next
    }

    tbl <- sort(table(subj_vals), decreasing = TRUE)
    props <- tbl / sum(tbl)
    ref <- names(tbl)[[1L]] # most frequent level = reference

    lvls <- names(props)[props >= catCutoff & names(props) != ref]
    dropped <- names(props)[props < catCutoff & names(props) != ref]

    cat_ref[[col]] <- ref
    cat_dropped[[col]] <- dropped

    if (length(lvls) == 0L) {
      cat_levels[[col]] <- character(0L)
      next
    }

    # Collision check -- stop rather than silently overwrite
    ind_names <- paste0(col, "_", lvls)
    collisions <- ind_names[ind_names %in% names(data)]
    if (length(collisions) > 0L) {
      stop(
        "Cannot create indicator column(s): ",
        paste(collisions, collapse = ", "),
        " -- name(s) already exist in data. ",
        "Rename the conflicting column(s) before calling runSCM().",
        call. = FALSE
      )
    }

    for (lev in lvls) {
      data[[paste0(col, "_", lev)]] <- as.integer(
        !is.na(data[[col]]) & as.character(data[[col]]) == lev
      )
    }

    cat_levels[[col]] <- lvls
  }

  # Trim to columns needed for fitting
  ind_cols <- unlist(
    lapply(catvarsVec, function(col) {
      lvls <- cat_levels[[col]]
      if (length(lvls) == 0L) character(0L) else paste0(col, "_", lvls)
    }),
    use.names = FALSE
  )
  fit_cols <- tryCatch(names(nlme::getData(fit)), error = function(e) {
    character(0L)
  })
  keep_cols <- unique(c(
    intersect(fit_cols, names(data)),
    intersect(c(catvarsVec, covarsVec), names(data)),
    ind_cols
  ))

  list(
    data = data[, keep_cols, drop = FALSE],
    catLevels = cat_levels,
    catRef = cat_ref,
    catDropped = cat_dropped
  )
}

#' Enrich a pairs data frame with covariate type, centering values, and levels
#'
#' For continuous covariates the median of the covariate column is computed from
#' \code{data} and stored in \code{center}; the model expression will be
#' \code{cov_theta * log(covariate / center)}.  For categorical covariates each
#' row is expanded into one row per non-reference level and the model expression
#' will be \code{cov_theta * raw_col_level} where \code{raw_col_level} is a
#' pre-computed 0/1 indicator column (e.g. \code{sex_male}) that must already
#' exist in the dataset (created during SCM categorical preprocessing).
#'
#' When \code{missingToken} is supplied (or the covariate column contains
#' \code{NA}s), three extra columns are added to support wrapping the
#' covariate expression in an \code{ifelse()} guard inside
#' \code{.expandShapes()}:
#' \describe{
#'   \item{has_missing}{logical -- whether any observations are missing}
#'   \item{missing_check}{character -- the \code{ifelse} condition string,
#'     e.g. \code{"is.na(wt)"} or \code{"is.na(wt) | wt == -99"}}
#'   \item{missing_fill}{integer placeholder (always \code{0L} for
#'     continuous rows; \code{0} or \code{1} for categorical rows based on
#'     whether the observed mode equals the indicator level).  For
#'     continuous covariates the actual shape-aware fill literal is built
#'     later in \code{.expandShapes()} via \code{.scmFillStr()}, because
#'     the appropriate value depends on the shape (e.g. \code{log(median)}
#'     for \code{"log"}, \code{median} for \code{"identity"}) and the
#'     shape is not yet known at \code{.enrichPairs()} time.}
#' }
#'
#' @param pairs data frame with at least \code{var} and \code{covar} columns
#' @param data  data frame used for fitting (source of medians and factor levels)
#' @param catvarsVec character vector of categorical covariate column names;
#'   rows whose \code{covar} matches an entry are expanded per level
#' @param missingToken sentinel value treated as missing in addition to
#'   \code{NA}; default \code{NA} (only true \code{NA}s detected)
#' @param catLevels optional named list produced by \code{.makeSCMData()},
#'   mapping each categorical column name to the qualifying non-reference
#'   levels.  When supplied, level selection uses this list directly instead
#'   of deriving levels from \code{data} -- ensuring only levels that passed
#'   the subject-frequency cutoff are tested.
#' Override median centers with fixed reference values
#'
#' Given the enriched \code{pairs} (which already carries \code{type},
#' \code{center}, and \code{raw_col}), replace the \code{center} of any
#' continuous row whose raw covariate column is named in \code{centers} with
#' the supplied fixed value.  A no-op when \code{centers} is \code{NULL}.
#' @noRd
.applyFixedCenters <- function(pairs, centers = NULL) {
  if (is.null(centers) || length(centers) == 0 || nrow(pairs) == 0) {
    return(pairs)
  }
  if (is.null(names(centers)) || any(!nzchar(names(centers)))) {
    stop("'centers' must be a named numeric vector, e.g. c(BW = 70, CrCL = 95)")
  }
  for (nm in names(centers)) {
    idx <- which(pairs$type == "continuous" & pairs$raw_col == nm)
    if (length(idx) > 0) {
      pairs$center[idx] <- as.numeric(centers[[nm]])
    }
  }
  pairs
}

#' @return enriched data frame with additional columns \code{type},
#'   \code{center}, \code{raw_col}, \code{level}, and (when missing values
#'   are present) \code{has_missing}, \code{missing_check}, \code{missing_fill}
#' @noRd
.enrichPairs <- function(
  pairs,
  data,
  catvarsVec = NULL,
  missingToken = NA,
  catLevels = NULL
) {
  if (!"type" %in% names(pairs)) {
    pairs$type <- "continuous"
  }
  if (!"center" %in% names(pairs)) {
    pairs$center <- NA_real_
  }
  if (!"raw_col" %in% names(pairs)) {
    pairs$raw_col <- pairs$covar
  }
  if (!"level" %in% names(pairs)) {
    pairs$level <- NA_character_
  }

  # Expand categorical covariate rows into one row per non-reference level
  if (!is.null(catvarsVec) && length(catvarsVec) > 0) {
    keep <- list()
    extra <- list()
    for (i in seq_len(nrow(pairs))) {
      rc <- pairs$covar[i]
      if (rc %in% catvarsVec) {
        if (!is.null(catLevels) && rc %in% names(catLevels)) {
          lvls <- catLevels[[rc]]
        } else {
          col_data <- data[[rc]]
          lvls <- if (is.factor(col_data)) {
            levels(col_data)
          } else {
            sort(unique(as.character(col_data[!is.na(col_data)])))
          }
          lvls <- lvls[-1] # drop first (alphabetical) level as reference
        }
        for (lev in lvls) {
          row <- pairs[i, , drop = FALSE]
          row$covar <- paste0(rc, "_", lev)
          row$type <- "categorical"
          row$center <- NA_real_
          row$raw_col <- rc
          row$level <- lev
          extra[[length(extra) + 1]] <- row
        }
      } else {
        keep[[length(keep) + 1]] <- pairs[i, , drop = FALSE]
      }
    }
    pairs <- do.call(rbind, c(keep, extra))
    rownames(pairs) <- NULL
  }

  # Compute medians for continuous covariates that don't already have a center.
  # Use one row per subject so that subjects with many observations do not
  # disproportionately influence the centering value.
  subj_data <- data[!duplicated(data[[.idColumn(data)]]), , drop = FALSE]
  for (i in seq_len(nrow(pairs))) {
    if (pairs$type[i] == "continuous" && is.na(pairs$center[i])) {
      col <- pairs$raw_col[i]
      if (col %in% names(subj_data)) {
        pairs$center[i] <- median(subj_data[[col]], na.rm = TRUE)
      }
    }
  }

  # -- Missing value detection ------------------------------------------------
  # For each row, check whether the covariate has NA or missingToken values.
  # When found, record three columns consumed by .expandShapes() to wrap the
  # covariate expression in an ifelse() guard:
  #
  #   has_missing   - logical flag
  #   missing_check - condition string for the ifelse() (e.g. "is.na(wt)")
  #   missing_fill  - replacement value when condition is TRUE:
  #                     continuous  -> placeholder 0L; the actual shape-aware
  #                                    fill literal (e.g. log(median) for
  #                                    the log shape, median for identity)
  #                                    is built later by .scmFillStr() in
  #                                    .expandShapes() once the shape is known.
  #                     categorical -> 1 if mode == level, else 0 (used as-is
  #                                    by .expandShapes()).
  pairs$has_missing <- FALSE
  pairs$missing_fill <- 0L
  pairs$missing_check <- NA_character_

  for (i in seq_len(nrow(pairs))) {
    rc <- pairs$raw_col[i]
    if (!rc %in% names(data)) {
      next
    }
    col_data <- data[[rc]]

    # Identify missing observations (NA and/or sentinel token)
    is_missing <- is.na(col_data)
    if (!is.na(missingToken)) {
      is_missing <- is_missing | (!is.na(col_data) & col_data == missingToken)
    }
    if (!any(is_missing)) {
      next
    }

    pairs$has_missing[i] <- TRUE

    # Build the ifelse() condition string
    check <- paste0("is.na(", rc, ")")
    if (!is.na(missingToken)) {
      # Quote the token if the data column is character/factor or the token
      # itself is a string; leave unquoted for numeric tokens on numeric columns
      quote_tok <- is.character(missingToken) ||
        inherits(col_data, c("character", "factor"))
      if (quote_tok) {
        check <- paste0(check, " | ", rc, ' == "', missingToken, '"')
      } else {
        check <- paste0(check, " | ", rc, " == ", missingToken)
      }
    }
    pairs$missing_check[i] <- check

    # Fill value: continuous -> 0 (median imputation on log scale)
    #             categorical -> 1 if mode of valid values equals this level
    if (pairs$type[i] == "continuous") {
      pairs$missing_fill[i] <- 0L
    } else {
      valid_vals <- as.character(col_data[!is_missing])
      if (length(valid_vals) == 0L) {
        pairs$missing_fill[i] <- 0L
      } else {
        mode_val <- names(sort(table(valid_vals), decreasing = TRUE))[[1L]]
        pairs$missing_fill[i] <- as.integer(mode_val == pairs$level[i])
      }
    }
  }

  pairs
}

#' Expand a pairs data frame into one row per (pair, shape)
#'
#' For continuous covariates each shape in \code{shapes} (or the per-row
#' \code{shapes} list-column when present) produces one output row.  The
#' \code{covar} column is suffixed with the shape name so that theta identifiers
#' remain unique (e.g. \code{"wt"} with shapes \code{c("power","lin")} yields
#' rows with \code{covar = "wt_power"} and \code{covar = "wt_lin"}).  Categorical
#' rows are not expanded; they always use the \code{"cat"} shape and their
#' \code{covar} name is unchanged.
#'
#' The \code{covExpr} column is pre-computed by \code{.scmCovExpr()} so that
#' downstream code can pass it directly to \code{.builduiCovariate()} without
#' repeating shape logic.
#'
#' @param pairs   enriched pairs data frame (output of \code{.enrichPairs()})
#' @param shapes  character vector of shape names applied globally to continuous
#'   pairs with no per-row override; default \code{"power"}
#' @param customShapes named list of additional shape builder functions; merged
#'   with \code{.SCM_SHAPES} (custom entries take precedence)
#' @param inits   named list mapping shape names to initial theta estimates,
#'   e.g. \code{list(power = 0.1, lin = 0.01, cat = 0.01)}.  Used as the
#'   global default; a per-pair \code{inits} list-column (set in
#'   \code{buildPairs()} from \code{pairsVec} items) takes precedence on a
#'   shape-by-shape basis, allowing different starting values per shape even
#'   within the same relationship.  Shapes not covered by either source
#'   default to \code{0}.
#' @return expanded pairs data frame with additional columns \code{shape},
#'   \code{covExpr}, \code{init}, \code{lower}, and \code{upper}; the
#'   \code{shapes} and \code{inits} list-columns (if present) are removed
#' @noRd
.expandShapes <- function(
  pairs,
  shapes = "power",
  customShapes = NULL,
  inits = list()
) {
  has_per_row <- "shapes" %in% names(pairs)
  has_per_inits <- "inits" %in% names(pairs)
  has_miss_col <- "has_missing" %in% names(pairs)
  result_rows <- vector("list", nrow(pairs) * length(shapes))
  out_idx <- 0L

  for (i in seq_len(nrow(pairs))) {
    row <- pairs[i, , drop = FALSE]
    cov_type <- row$type
    raw_col <- row$raw_col
    center <- if (is.na(row$center)) NULL else row$center
    level <- if (is.na(row$level)) NULL else row$level
    base_covar <- row$covar

    # Per-pair inits: named list keyed by shape (e.g. list(power=0.75, lin=0.01))
    pair_inits <- if (has_per_inits) row$inits[[1]] else list()

    # Missing-value guard (populated by .enrichPairs())
    has_miss <- has_miss_col && isTRUE(row$has_missing)
    miss_check <- if (has_miss) row$missing_check else NA_character_
    miss_fill <- if (has_miss) row$missing_fill else 0L

    if (cov_type == "categorical") {
      expr <- .scmCovExpr(
        "cat",
        raw_col, # nolint: object_usage_linter.
        level = level,
        customShapes = customShapes
      )
      if (has_miss && !is.na(miss_check)) {
        expr <- paste0("ifelse(", miss_check, ", ", miss_fill, ", ", expr, ")")
      }
      spec_cat <- if (!is.null(pair_inits[["cat"]])) {
        pair_inits[["cat"]]
      } else if (!is.null(inits[["cat"]])) {
        inits[["cat"]]
      } else {
        NULL
      }
      parsed_cat <- .parseInitSpec(spec_cat) # nolint: object_usage_linter.
      row$shape <- "cat"
      row$covExpr <- expr
      row$init <- parsed_cat$est
      row$lower <- parsed_cat$lower
      row$upper <- parsed_cat$upper
      if (has_per_row) {
        row$shapes <- NULL
      }
      if (has_per_inits) {
        row$inits <- NULL
      }
      out_idx <- out_idx + 1L
      result_rows[[out_idx]] <- row
    } else {
      # Determine which shapes to test for this continuous row
      row_shapes <- if (has_per_row && length(row$shapes[[1]]) > 0) {
        row$shapes[[1]]
      } else {
        shapes
      }
      for (sh in row_shapes) {
        r <- row
        expr <- .scmCovExpr(
          sh,
          raw_col, # nolint: object_usage_linter.
          center = center,
          customShapes = customShapes
        )
        if (has_miss && !is.na(miss_check)) {
          # Shape-aware population-typical fill: imputing the median means
          # evaluating the shape expression at cov = center.  For built-in
          # shapes this gives 0 for power/lin (centered), log(center) for
          # log, and center itself for identity.  Using 0 unconditionally
          # (as the earlier implementation did) would correspond to
          # imputing cov = 1 for the log shape and cov = 0 for identity --
          # neither is the population median.
          fill_str <- .scmFillStr(
            shape = sh,
            center = center,
            level = level,
            customShapes = customShapes
          )
          expr <- paste0(
            "ifelse(",
            miss_check,
            ", ",
            fill_str,
            ", ",
            expr,
            ")"
          )
        }
        spec_sh <- if (!is.null(pair_inits[[sh]])) {
          pair_inits[[sh]]
        } else if (!is.null(inits[[sh]])) {
          inits[[sh]]
        } else {
          NULL
        }
        parsed_sh <- .autoScaleInitSpec( # nolint: object_usage_linter.
          spec_sh, sh, center,
          label = paste0(sh, "~", base_covar)
        )
        r$shape <- sh
        r$covExpr <- expr
        r$covar <- paste0(base_covar, "_", sh)
        r$init <- parsed_sh$est
        r$lower <- parsed_sh$lower
        r$upper <- parsed_sh$upper
        if (has_per_row) {
          r$shapes <- NULL
        }
        if (has_per_inits) {
          r$inits <- NULL
        }
        out_idx <- out_idx + 1L
        result_rows[[out_idx]] <- r
      }
    }
  }

  result <- do.call(rbind, result_rows[seq_len(out_idx)])
  rownames(result) <- NULL
  result
}

#' Update the \code{init} column of a pairs data frame from a fitted model's estimates
#'
#' After a covariate is accepted or a step completes, update the \code{init}
#' values for all rows whose \code{cov_<covar>_<var>} theta is present in
#' \code{fit}.  Rows not found in the fit are left unchanged.  Bounds
#' (\code{lower}, \code{upper}) are never modified -- only the starting estimate.
#'
#' @param pairs data frame with at minimum \code{var} and \code{covar} columns
#' @param fit   nlmixr2 fit object
#' @return \code{pairs} with \code{init} updated from \code{fit$parFixedDf}
#' @noRd
.updatePairsInits <- function(pairs, fit) {
  if (is.null(pairs) || nrow(pairs) == 0L) {
    return(pairs)
  }
  pf <- tryCatch(fit$parFixedDf, error = function(e) NULL)
  if (is.null(pf)) {
    return(pairs)
  }
  if (!"init" %in% names(pairs)) {
    pairs$init <- 0
  }
  for (i in seq_len(nrow(pairs))) {
    theta_nm <- paste0("cov_", pairs$covar[i], "_", pairs$var[i])
    est_val <- tryCatch(pf[theta_nm, "Estimate"], error = function(e) NA_real_)
    if (!is.null(est_val) && length(est_val) == 1L && !is.na(est_val)) {
      pairs$init[i] <- est_val
    }
  }
  pairs
}

#' Freeze a candidate UI for a cheap 1-D profile of one covariate theta
#'
#' Fixes every population theta except \code{freeTheta} at its current
#' estimate (inherited from the parent fit via \code{base_ui$iniDf}) and
#' zeroes all between-subject variability, yielding a fast fixed-effect model
#' whose only free parameter is the new covariate coefficient.  This is the
#' "frozen base" on which the 1-D Brent profile runs.
#'
#' @param ui        candidate rxode2 UI (base + context + new covariate)
#' @param freeTheta name of the single theta left free (e.g. "cov_BW_cl")
#' @return a decompressed rxode2 UI with etas removed and all but one theta
#'   fixed
#' @noRd
.freezeUiForProfile <- function(ui, freeTheta) {
  ui <- rxode2::rxUiDecompress(ui)
  ini <- ui$iniDf
  isTheta <- !is.na(ini$ntheta)
  ini$fix[isTheta] <- TRUE
  ini$fix[isTheta & ini$name == freeTheta] <- FALSE
  ui$iniDf <- ini
  # Drop omega (between-subject variability) so the profile is a fast
  # fixed-effect fit; residual error stays fixed at the parent value.
  ui <- rxode2::zeroRe(ui, which = "omega")
  ui
}

#' Profile a single new covariate theta with a 1-D Brent search
#'
#' Builds the candidate UI (base + accepted context + the one new covariate),
#' freezes it via \code{.freezeUiForProfile()}, and runs a derivative-free 1-D
#' Brent optimisation (\code{est = "optim"}) over the single free covariate
#' coefficient.  Used to seed gradient optimisers (nlminb / lbfgsb3c) with a
#' nonzero, gradient-informative starting value so they do not stall at the
#' flat zero-effect point.  bobyqa is never involved: the profiled value is
#' handed back to the fit's own estimator by the caller.
#'
#' @param base_ui   clean base UI (\code{fit$finalUiEnv})
#' @param ctx_df    accepted-covariate context pairs (or \code{NULL})
#' @param nam_var   parameter name (e.g. "cl")
#' @param nam_covar covariate name (e.g. "BW")
#' @param cov_expr  covariate model expression string (or \code{NULL})
#' @param cov_init  fallback initial value for the new theta (Brent start)
#' @param cov_lower lower bound (Brent bracket)
#' @param cov_upper upper bound (Brent bracket)
#' @param covNames  theta name of the new covariate (e.g. "cov_BW_cl")
#' @param data      modelling data frame
#' @return profiled estimate (finite numeric strictly inside the bracket) or
#'   \code{NA_real_} when profiling fails or lands on a bound
#' @noRd
.profileCovInit <- function(base_ui, ctx_df, nam_var, nam_covar, cov_expr,
                            cov_init, cov_lower, cov_upper, covNames, data) {
  cand_df <- data.frame(
    var     = nam_var,
    covar   = nam_covar,
    covExpr = if (!is.null(cov_expr)) cov_expr else nam_covar,
    init    = cov_init,
    lower   = cov_lower,
    upper   = cov_upper,
    stringsAsFactors = FALSE
  )
  all_pairs <- if (!is.null(ctx_df)) rbind(ctx_df, cand_df) else cand_df
  ui        <- .rebuildUiFromPairs(base_ui, all_pairs)
  frozen    <- .freezeUiForProfile(ui, covNames)

  # est = "optim" with a single free parameter dispatches to stats::optim
  # method = "Brent", a safe derivative-free 1-D search over [lower, upper].
  # calcTables = FALSE skips the (unneeded) output tables; the profiled theta
  # is read from the fit core's fixed-effect vector regardless.  Fall back to
  # a bare Brent control if the installed nlmixr2est rejects an argument.
  ctrl <- tryCatch(
    nlmixr2est::optimControl(
      method     = "Brent",
      calcTables = FALSE,
      print      = 0L
    ),
    error = function(e) nlmixr2est::optimControl(method = "Brent")
  )
  pf <- suppressWarnings(
    nlmixr2est::nlmixr2(frozen, data, "optim", control = ctrl)
  )
  val <- tryCatch(unname(pf$theta[[covNames]]), error = function(e) NA_real_)
  if (is.null(val) || length(val) != 1L || !is.finite(val)) {
    return(NA_real_)
  }
  # Reject boundary solutions: a theta pinned to the bracket edge is not a
  # trustworthy warm start (caller falls back to the default init).
  if (val <= cov_lower || val >= cov_upper) {
    return(NA_real_)
  }
  val
}

#' Filter a pairs data frame to only those pairs present in a fitted model
#'
#' @param pairs data frame with columns \code{var} and \code{covar}
#' @param fit nlmixr2 fit object
#' @return subset of \code{pairs} whose \code{cov_<covar>_<var>} theta exists in \code{fit}
#' @noRd
.pairsInModel <- function(pairs, fit) {
  mrd <- tryCatch(fit$finalUiEnv$muRefCovariateDataFrame, error = function(e) {
    NULL
  })
  if (!is.null(mrd) && nrow(mrd) > 0 && "covariateParameter" %in% names(mrd)) {
    present <- mrd$covariateParameter
  } else {
    ini_names <- tryCatch(fit$finalUiEnv$iniDf$name, error = function(e) {
      character(0)
    })
    present <- ini_names[grepl("^cov_", ini_names)]
  }
  expected <- paste0("cov_", pairs$covar, "_", pairs$var)
  pairs[expected %in% present, , drop = FALSE]
}

#' Compute % reduction in omega variance for the eta associated with a variable
#'
#' @param fit_without fit without the covariate (higher omega expected)
#' @param fit_with    fit with the covariate
#' @param nam_var     variable name (e.g. \code{"cl"})
#' @return numeric % reduction, or \code{NA_real_} if the eta cannot be resolved
#' @noRd
.bsvReduction <- function(fit_without, fit_with, nam_var) {
  tryCatch(
    {
      theta_name <- .getThetaName(fit_without$finalUiEnv, nam_var)
      mrd <- fit_without$muRefDataFrame
      eta <- mrd$eta[mrd$theta == theta_name]
      if (length(eta) == 0) {
        return(NA_real_)
      }
      eta <- eta[1]
      om_without <- diag(fit_without$omega)[eta]
      om_with <- diag(fit_with$omega)[eta]
      if (is.na(om_without) || om_without == 0) {
        return(NA_real_)
      }
      (om_without - om_with) / om_without * 100
    },
    error = function(e) NA_real_
  )
}

#' Resolve the effective OFV tolerance for criterion 1 of .isUnrealisticOFV()
#'
#' For deterministic estimators any OFV worsening is suspect; for stochastic
#' estimators (SAEM) Monte Carlo noise can legitimately land a few OFV units
#' above the parent even for a well-converged, uninformative candidate.  A
#' tolerance of 10 is safely above typical SAEM noise (1-5 units) while
#' catching genuine numerical failures (hundreds to millions of units).
#'
#' @param fit   reference nlmixr2 fit object (used for \code{fit$est})
#' @param retryOFVTolerance numeric or \code{NULL}; when non-\code{NULL}
#'   returned as-is; when \code{NULL} auto-detected from \code{fit$est}
#' @return numeric scalar
#' @noRd
.resolveOFVTolerance <- function(fit, retryOFVTolerance) {
  if (!is.null(retryOFVTolerance)) return(as.numeric(retryOFVTolerance))
  if (tolower(fit$est) %in% c("saem")) 10 else 0
}

#' Check whether a candidate OFV result is unrealistic and should trigger a retry
#'
#' @param x_objf           OFV of the candidate fit
#' @param ref_objf         OFV of the reference (parent) fit
#' @param dObjf            OFV change: candidate minus reference (negative = improvement)
#' @param pchisqr          chi-squared p-value for the improvement
#' @param maxDeltaOFV      user-specified ceiling on plausible absolute OFV change
#' @param effective_tolerance margin added to ref_objf before criterion 1 fires
#' @return logical scalar
#' @noRd
.isUnrealisticOFV <- function(x_objf, ref_objf, dObjf, pchisqr,
                               maxDeltaOFV, effective_tolerance) {
  x_objf > ref_objf + effective_tolerance ||
    pchisqr < .Machine$double.eps ||
    abs(dObjf) > maxDeltaOFV
}

#' Fit all candidate models for one SCM step
#'
#' Builds a candidate UI for every row of \code{pairs} (adding or removing the
#' covariate-parameter relationship), fits each one, and assembles a per-row
#' result list.  Rows whose candidate cannot be built or fitted are silently
#' dropped.
#'
#' @param pairs    data frame with columns \code{var} and \code{covar}
#' @param ui       current model UI (rxode2 environment)
#' @param fit      current reference fit (nlmixr2FitCore)
#' @param data     data frame passed to \code{nlmixr2}
#' @param pVal     p-value threshold (used only to compute the \code{qchisqr} label)
#' @param stepIdx  integer step counter stored in the stats table
#' @param add      logical; \code{TRUE} = forward (add covariate),
#'                 \code{FALSE} = backward (remove covariate)
#' @return list of \code{list(var, covar, fit, stats)}, one per successful candidate
#' @noRd
.fitCandidatePairs <- function(
  pairs,
  base_ui,
  context_pairs = NULL,
  fit,
  data,
  pVal,
  stepIdx,
  add,
  control = NULL,
  print = 100,
  profileInit = FALSE,
  profileInitOnStall = TRUE,
  stallTol = 0,
  maxRetries = 3L,
  maxDeltaOFV = Inf,
  retryPerturbSD = 0.5,
  retrySmallInit = 0.01,
  effective_tolerance = 0,
  retryFailOnExhaustion = FALSE,
  rxThreads = NULL
) {
  raw_results <- nlmixr2utils::.plap(
    # nolint: object_usage_linter.
    seq_len(nrow(pairs)),
    function(i) {
      nam_var <- pairs$var[i]
      nam_covar <- pairs$covar[i]
      covNames <- paste0("cov_", nam_covar, "_", nam_var)

      # Extract pre-computed covExpr, shape, init, lower, upper (set by .expandShapes())
      cov_expr <- if ("covExpr" %in% names(pairs)) pairs$covExpr[i] else NULL
      cov_shape <- if ("shape" %in% names(pairs)) {
        pairs$shape[i]
      } else {
        NA_character_
      }
      cov_init <- if ("init" %in% names(pairs) && !is.na(pairs$init[i])) {
        pairs$init[i]
      } else {
        0.1
      }
      cov_lower <- if ("lower" %in% names(pairs) && !is.na(pairs$lower[i])) {
        pairs$lower[i]
      } else {
        -5
      }
      cov_upper <- if ("upper" %in% names(pairs) && !is.na(pairs$upper[i])) {
        pairs$upper[i]
      } else {
        5
      }

      cand_control <- tryCatch(
        {
          ctrl <- if (!is.null(control)) control else fit$control
          ctrl$print <- print
          ctrl
        },
        error = function(e) NULL
      )

      direction <- if (add) "Forward" else "Backward"
      shape_lbl <-
        if (!is.na(cov_shape) && nzchar(cov_shape)) {
          paste0(" [", cov_shape, "]")
        } else {
          ""
        }
      cli::cli_rule()
      cli::cli_inform(c(
        ">" = "{direction} step {stepIdx}, candidate {i}/{nrow(pairs)}: {nam_covar} ~ {nam_var}{shape_lbl}"
      ))
      cli::cli_rule()

      orig_cov_init <- cov_init

      # Build ctx_df once -- context pairs do not change between retry attempts
      ctx_df <- if (!is.null(context_pairs) && nrow(context_pairs) > 0) {
        data.frame(
          var    = context_pairs$var,
          covar  = context_pairs$covar,
          covExpr = if ("covExpr" %in% names(context_pairs)) {
            context_pairs$covExpr
          } else {
            context_pairs$covar
          },
          init  = if ("init"  %in% names(context_pairs)) context_pairs$init  else rep(0.1, nrow(context_pairs)),
          lower = if ("lower" %in% names(context_pairs)) context_pairs$lower else rep(-5,  nrow(context_pairs)),
          upper = if ("upper" %in% names(context_pairs)) context_pairs$upper else rep(5,   nrow(context_pairs)),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }

      # Optional 1-D Brent warm-start: profile the single new covariate theta
      # on a frozen base (all other thetas fixed at their parent estimates,
      # between-subject variability zeroed) to obtain a gradient-informative
      # starting value before the real estimator runs.  Forward path only.
      if (isTRUE(profileInit) && add) {
        prof_init <- tryCatch(
          .profileCovInit(
            base_ui   = base_ui,
            ctx_df    = ctx_df,
            nam_var   = nam_var,
            nam_covar = nam_covar,
            cov_expr  = cov_expr,
            cov_init  = orig_cov_init,
            cov_lower = cov_lower,
            cov_upper = cov_upper,
            covNames  = covNames,
            data      = data
          ),
          error = function(e) NA_real_
        )
        if (is.finite(prof_init)) {
          cli::cli_inform(c(
            "i" = "Brent warm-start {covNames}: {round(orig_cov_init, 4)} -> {round(prof_init, 4)}"
          ))
          orig_cov_init <- prof_init
        }
      }

      # For backward elimination: build the candidate UI once (no init to vary)
      ui_cand <- if (!add) {
        cp_minus <- if (!is.null(context_pairs) && nrow(context_pairs) > 0) {
          context_pairs[
            !(context_pairs$var == nam_var & context_pairs$covar == nam_covar),
            , drop = FALSE
          ]
        } else {
          NULL
        }
        tryCatch(
          .rebuildUiFromPairs(base_ui, cp_minus),
          error = function(e) list(
            .failed = TRUE, .reason = conditionMessage(e),
            .pair = paste0(nam_covar, " ~ ", nam_var)
          )
        )
      } else {
        NULL  # forward path builds ui inside the loop per attempt
      }
      if (isTRUE(ui_cand$.failed)) return(ui_cand)

      loop_result <- NULL
      # Best attempt seen across retries (forward path only).  
      # "Best" = the candidate with the largest dObjf. 
      best_attempt <- NULL

      for (attempt in 0:maxRetries) {
        cov_init_attempt <- if (attempt == 0L) {
          orig_cov_init
        } else if (attempt %% 2L == 1L) {
          orig_cov_init * exp(stats::rnorm(1L, 0, retryPerturbSD))
        } else {
          retrySmallInit
        }

        ui_attempt <- tryCatch(
          {
            if (add) {
              cand_df_a <- data.frame(
                var     = nam_var,
                covar   = nam_covar,
                covExpr = if (!is.null(cov_expr)) cov_expr else nam_covar,
                init    = cov_init_attempt,
                lower   = cov_lower,
                upper   = cov_upper,
                stringsAsFactors = FALSE
              )
              all_pairs_a <- if (!is.null(ctx_df)) rbind(ctx_df, cand_df_a) else cand_df_a
              .rebuildUiFromPairs(base_ui, all_pairs_a)
            } else {
              ui_cand
            }
          },
          error = function(e) {
            list(
              .failed = TRUE,
              .reason = conditionMessage(e),
              .pair   = paste0(nam_covar, " ~ ", nam_var)
            )
          }
        )
        if (isTRUE(ui_attempt$.failed)) return(ui_attempt)

        x <- tryCatch(
          suppressWarnings(nlmixr2est::nlmixr2(
            ui_attempt,
            data,
            fit$est,
            control = cand_control
          )),
          error = function(e) {
            list(
              .failed = TRUE,
              .reason = conditionMessage(e),
              .pair   = paste0(nam_covar, " ~ ", nam_var)
            )
          }
        )

        if (isTRUE(x$.failed)) {
          if (attempt < maxRetries) next else return(x)
        }

        dObjf <- x$objf - fit$objf
        if (is.na(dObjf)) {
          if (attempt < maxRetries) {
            next
          } else {
            return(list(
              .failed = TRUE,
              .reason = "NA objective function value",
              .pair   = paste0(nam_covar, " ~ ", nam_var)
            ))
          }
        }

        dof <- if (add) {
          length(x$finalUiEnv$ini$est) - length(fit$finalUiEnv$ini$est)
        } else {
          length(fit$finalUiEnv$ini$est) - length(x$finalUiEnv$ini$est)
        }
        # The LRT is one-sided: only the direction indicating a meaningful change
        # (OFV drop for forward, OFV increase for backward) yields a test
        # statistic.  A change in the wrong direction gets p = 1.
        chisq_stat <- if (add) -dObjf else dObjf
        pchisqr <- if (chisq_stat > 0) 1 - stats::pchisq(chisq_stat, df = dof) else 1

        # For backward elimination the unrealistic-OFV criteria are not applied:
        # removing a significant covariate legitimately produces a large OFV
        # increase (pchisqr approximately 0), which would incorrectly trigger criterion 2.
        # Backward steps always accept the first converged result.
        if (!add) {
          loop_result <- list(x = x, dObjf = dObjf, dof = dof, pchisqr = pchisqr)
          break
        }

        # Selection rule: candidate wins on the first attempt (best is NULL)
        # or when its dObjf (candidate - parent OFV) is strictly below the
        # incumbent's, i.e. it reached a lower OFV; ties keep the
        # incumbent.  Without this tracker, the exhaustion branch silently
        # kept whichever attempt happened to be last -- under perturbed-init
        # retries that was the best only by coincidence.  The corresponding
        # regression tests live in tests/testthat/test-scm.R.
        .cand_attempt <- list(
          x = x, dObjf = dObjf, dof = dof, pchisqr = pchisqr,
          attempt_num = attempt + 1L
        )
        if (is.null(best_attempt) ||
            .cand_attempt$dObjf < best_attempt$dObjf) {
          best_attempt <- .cand_attempt
        }

        if (!.isUnrealisticOFV(
          x$objf, fit$objf, dObjf, pchisqr, maxDeltaOFV, effective_tolerance
        )) {
          loop_result <- list(x = x, dObjf = dObjf, dof = dof, pchisqr = pchisqr)
          break
        }

        trigger <- if (x$objf > fit$objf + effective_tolerance) {
          paste0(
            "OFV increased vs parent (",
            round(x$objf, 3), " > ", round(fit$objf, 3),
            if (effective_tolerance > 0) paste0(" + ", effective_tolerance) else "",
            ")"
          )
        } else if (pchisqr < .Machine$double.eps) {
          paste0("p-value underflow (dOFV = ", round(dObjf, 3), ")")
        } else {
          paste0("|dOFV| (", round(abs(dObjf), 3), ") exceeds maxDeltaOFV (", maxDeltaOFV, ")")
        }

        if (attempt < maxRetries) {
          next_strategy <- if ((attempt + 1L) %% 2L == 1L) "perturbed" else "small"
 
          cli::cli_warn(c(
            "!" = paste0(
              "{nam_covar} ~ {nam_var}: unrealistic OFV on attempt ",
              attempt + 1L, "/", maxRetries + 1L, ": ", trigger, "."
            ),
            "i" = "Retrying with {next_strategy} init."
          ))
        } else {
          if (retryFailOnExhaustion) {
            return(list(
              .failed = TRUE,
              .reason = paste0(
                "OFV unrealistic after ", maxRetries,
                if (maxRetries == 1L) " retry" else " retries",
                ": ", trigger
              ),
              .pair = paste0(nam_covar, " ~ ", nam_var)
            ))
          } else {
 
            cli::cli_warn(c("!" = paste0(
              "{nam_covar} ~ {nam_var}: unrealistic OFV after all ",
              maxRetries + 1L, " attempt",
              if (maxRetries + 1L == 1L) "" else "s",
              ": ", trigger, ". Accepting best available result (attempt ",
              best_attempt$attempt_num, "/", maxRetries + 1L,
              ", dOFV = ", round(best_attempt$dObjf, 3), ")."
            )))
            loop_result <- list(
              x = best_attempt$x, dObjf = best_attempt$dObjf,
              dof = best_attempt$dof, pchisqr = best_attempt$pchisqr
            )
            break
          }
        }
      }

      x       <- loop_result$x
      dObjf   <- loop_result$dObjf
      dof     <- loop_result$dof
      pchisqr <- loop_result$pchisqr

      # -- Profile-on-stall rescue (forward path only) ----------------------
      # Independent of the retry budget: runs AFTER the retry loop so it fires
      # even when maxRetries = 0 (the fast benchmark setting).  A forward
      # candidate whose accepted fit still stalled -- OFV improvement
      # (-dObjf) <= stallTol, i.e. the nested model is no better than its
      # parent, which is impossible at a
      # true optimum -- means the outer optimiser never moved the covariate off
      # its init.  We run ONE frozen 1-D .profileCovInit() (all other thetas
      # fixed, BSV fixed-but-present) to obtain a basin-correct coefficient,
      # then refit the full candidate from that init.  The rescue result is
      # kept ONLY if it strictly lowers dObjf, so it can never make a
      # candidate worse.  Healthy candidates (-dObjf > stallTol, e.g. analytic
      # linCmt fits) skip this entirely and pay no cost.
      if (add && isTRUE(profileInitOnStall) &&
          is.finite(dObjf) && -dObjf <= stallTol) {
        prof_val <- tryCatch(
          .profileCovInit(
            base_ui   = base_ui,
            ctx_df    = ctx_df,
            nam_var   = nam_var,
            nam_covar = nam_covar,
            cov_expr  = cov_expr,
            cov_init  = orig_cov_init,
            cov_lower = cov_lower,
            cov_upper = cov_upper,
            covNames  = covNames,
            data      = data
          ),
          error = function(e) NA_real_
        )
        if (is.finite(prof_val)) {
          cand_df_p <- data.frame(
            var     = nam_var,
            covar   = nam_covar,
            covExpr = if (!is.null(cov_expr)) cov_expr else nam_covar,
            init    = prof_val,
            lower   = cov_lower,
            upper   = cov_upper,
            stringsAsFactors = FALSE
          )
          all_pairs_p <- if (!is.null(ctx_df)) rbind(ctx_df, cand_df_p) else cand_df_p
          x_p <- tryCatch(
            {
              ui_p <- .rebuildUiFromPairs(base_ui, all_pairs_p)
              suppressWarnings(nlmixr2est::nlmixr2(
                ui_p, data, fit$est, control = cand_control
              ))
            },
            error = function(e) NULL
          )
          if (!is.null(x_p) && !isTRUE(x_p$.failed) &&
              !is.null(x_p$objf) && is.finite(x_p$objf)) {
            dObjf_p <- x_p$objf - fit$objf
            if (is.finite(dObjf_p) && dObjf_p < dObjf) {
              dof_p <- length(x_p$finalUiEnv$ini$est) -
                length(fit$finalUiEnv$ini$est)
              pchisqr_p <- if (dObjf_p < 0) {
                1 - stats::pchisq(-dObjf_p, df = dof_p)
              } else {
                1
              }
              cli::cli_inform(c(
                "v" = paste0(
                  "{nam_covar} ~ {nam_var}: profile-on-stall rescue ",
                  "dOFV ", round(dObjf, 3), " -> ", round(dObjf_p, 3),
                  " (profile init ", round(prof_val, 4), ")."
                )
              ))
              x       <- x_p
              dObjf   <- dObjf_p
              dof     <- dof_p
              pchisqr <- pchisqr_p
            }
          }
        }
      }

      covarEffect <- if (add) {
        x$parFixedDf[covNames, "Estimate"]
      } else {
        fit$parFixedDf[covNames, "Estimate"]
      }
      bsvReduction <- if (add) {
        .bsvReduction(fit_without = fit, fit_with = x, nam_var = nam_var)
      } else {
        .bsvReduction(fit_without = x, fit_with = fit, nam_var = nam_var)
      }

      stats <- data.frame(
        step        = stepIdx,
        covar       = nam_covar,
        var         = nam_var,
        shape       = cov_shape,
        objf        = x$objf,
        deltObjf    = dObjf,
        AIC         = x$AIC,
        BIC         = x$BIC,
        numParams   = length(x$finalUiEnv$ini$est),
        qchisqr     = stats::qchisq(1 - pVal, dof),
        pchisqr     = pchisqr,
        # `included` token semantics:
        #   forward  : "yes"      = added to model      | "no"       = tested, not added
        #   backward : "dropped"  = removed from model  | "retained" = tested, kept in model
        included    = if (add) "no" else "retained",
        searchType  = if (add) "forward" else "backward",
        covNames    = covNames,
        covarEffect = covarEffect,
        bsvReduction = bsvReduction,
        stringsAsFactors = FALSE
      )

      list(var = nam_var, covar = nam_covar, fit = x, stats = stats)
    },
    rxThreads = rxThreads,
    .label = function(i) {
      sh <- if (
        "shape" %in%
          names(pairs) &&
          !is.na(pairs$shape[i]) &&
          nzchar(pairs$shape[i])
      ) {
        paste0(" [", pairs$shape[i], "]")
      } else {
        ""
      }
      paste0(pairs$covar[i], " ~ ", pairs$var[i], sh)
    }
  )

  # Separate failures from successes and report
  is_fail <- vapply(raw_results, function(r) isTRUE(r$.failed), logical(1))
  failures <- raw_results[is_fail]
  results <- raw_results[!is_fail]

  if (length(failures) > 0L) {
    fail_msgs <- vapply(
      failures,
      function(f) paste0(f$.pair, ": ", f$.reason),
      character(1)
    )
    
    # Escape any curly braces in the reason strings before passing to cli,
    # because error messages from nlmixr2/cli themselves may contain `{`/`}`
    # which cli would try to parse as glue expressions.
    fail_msgs_safe <- gsub("\\{", "{{", gsub("\\}", "}}", fail_msgs))
    cli::cli_warn(c(
      "!" = "{length(failures)} candidate fit(s) failed at step {stepIdx}:",
      setNames(fail_msgs_safe, rep("x", length(fail_msgs_safe)))
    ))
  }

  if (length(results) == 0L) {
    stop(
      "All ",
      nrow(pairs),
      " candidate fit(s) failed at step ",
      stepIdx,
      ". Cannot continue.",
      call. = FALSE
    )
  }

  results
}

# -- Search functions -----------------------------------------------------------

#' Forward covariate search
#'
#' @param pairs     data frame with columns \code{var} and \code{covar}
#' @param fit       base nlmixr2 fit object
#' @param data      data frame for fitting
#' @param pVal      p-value threshold for inclusion
#' @param outputDir cache directory path
#' @param restart   logical; \code{FALSE} = resume from cache (only when
#'   \code{saveModels = TRUE})
#' @param saveModels logical; if \code{FALSE} all \code{saveRDS}/\code{readRDS}
#'   operations are skipped
#' @param verbose logical; if \code{TRUE} print candidate table before each
#'   step, full results table after fitting, and parameter/omega detail when a
#'   covariate is accepted
#' Pick the winning forward-search candidate
#'
#' Winner = most significant candidate (smallest p-value). For \code{df = 1}
#' any \code{dOFV >= ~70.5} makes \code{1 - pchisq()} underflow to exactly 0,
#' so several strong candidates tie at \code{pchisqr == 0}. Such ties are
#' broken by the most negative \code{deltObjf} (candidate minus reference OFV,
#' i.e. the biggest OFV drop) rather than by row
#' order, which would otherwise pick the first candidate alphabetically.
#'
#' @param resTable candidate results table with \code{pchisqr} and
#'   \code{deltObjf} columns
#' @return integer row index of the winning candidate
#' @noRd
.pickForwardWinner <- function(resTable) {
  order(resTable$pchisqr, resTable$deltObjf)[1]
}

#' Pick the backward-search candidate to drop
#'
#' The candidate to drop is the one whose removal causes the LEAST significant
#' OFV increase (highest p-value = least important). Removals with
#' \code{dOFV <= 0} all report \code{pchisqr == 1}, so such ties are broken by
#' the SMALLEST \code{deltObjf} (least OFV increase on removal) rather than by
#' row order.
#'
#' @param resTable candidate results table with \code{pchisqr} and
#'   \code{deltObjf} columns
#' @return integer row index of the candidate to drop
#' @noRd
.pickBackwardWinner <- function(resTable) {
  order(-resTable$pchisqr, resTable$deltObjf)[1]
}

#' @return \code{list(final_fit, resTableComplete)}
#' @noRd
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
forwardSearch <- function(
  pairs,
  fit,
  data,
  pVal = 0.05,
  outputDir,
  restart = FALSE,
  saveModels = TRUE,
  verbose = FALSE,
  control = NULL,
  print = 100,
  profileInit = FALSE,
  profileInitOnStall = TRUE,
  stallTol = 0,
  maxRetries = 3L,
  maxDeltaOFV = Inf,
  retryPerturbSD = 0.5,
  retrySmallInit = 0.01,
  retryOFVTolerance = NULL,
  retryFailOnExhaustion = FALSE,
  rxThreads = NULL
) {
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  if (saveModels && missing(outputDir)) {
    stop("please specify outputDir for forward search")
  }

  effective_tolerance <- .resolveOFVTolerance(fit, retryOFVTolerance)
  base_ui <- fit$finalUiEnv # fixed -- never updated
  orig_pairs <- pairs # preserved for resume reconstruction
  accepted_pairs <- NULL
  resTableComplete <- NULL
  stepIdx <- 1

  # When not saving there is no cache to read; treat as a clean start
  if (!saveModels) {
    restart <- TRUE
  }

  if (!restart) {
    fileExistsTab <- list.files(
      paste0("./", outputDir),
      pattern = "forward_step_[0-9]+_table_[a-zA-Z0-9_]+\\.rds"
    )
    fileExistsFit <- list.files(
      paste0("./", outputDir),
      pattern = "forward_step_[0-9]+_fit_[a-zA-Z0-9_]+\\.rds"
    )
    fileExistsCT <- list.files(
      paste0("./", outputDir),
      pattern = "forward_step_[0-9]+_completetable_[a-zA-Z0-9_]+\\.rds"
    )

    if (length(fileExistsTab) == 0) {
      restart <- TRUE
    } else {
      resumeRows <- lapply(fileExistsTab, function(f) {
        readRDS(paste0(outputDir, "/", f))
      })
      resumeTable <- data.table::rbindlist(resumeRows)
      fit <- readRDS(paste0(
        outputDir,
        "/",
        fileExistsFit[[length(fileExistsFit)]]
      ))
      accepted_pairs <- .pairsInModel(orig_pairs, fit)
      resTableComplete <- readRDS(paste0(
        outputDir,
        "/",
        fileExistsCT[[length(fileExistsCT)]]
      ))
      stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1
      for (i in seq_len(nrow(resumeTable))) {
        pairs <- pairs[
          !(pairs$var == resumeTable$var[i] &
            pairs$covar == resumeTable$covar[i]),
        ]
      }
      cli::cli_alert_success(
        "loaded forward search data from disk, resuming ..."
      )
    }
  }

  cli::cli_h1("starting forward search...")

  while (nrow(pairs) > 0) {
    # -- verbose: current model code + candidates before fitting -----------
    if (verbose) {
      cli::cli_h2("Forward step {stepIdx}: current model")
      tryCatch(
        print(fit$finalUiEnv),
        error = function(e) {
          cli::cli_alert_warning(
            "Model code unavailable: {conditionMessage(e)}"
          )
        }
      )
      cli::cli_h2("Forward step {stepIdx}: {nrow(pairs)} candidate(s) to test")
      v_cols <- intersect(
        c("var", "covar", "shape", "init", "covExpr"),
        names(pairs)
      )
      print(pairs[, v_cols, drop = FALSE])
    }

    results <- .fitCandidatePairs(
      pairs,
      base_ui,
      accepted_pairs,
      fit,
      data,
      pVal,
      stepIdx,
      add = TRUE,
      control = control,
      print = print,
      profileInit = profileInit,
      profileInitOnStall = profileInitOnStall,
      stallTol = stallTol,
      maxRetries = maxRetries,
      maxDeltaOFV = maxDeltaOFV,
      retryPerturbSD = retryPerturbSD,
      retrySmallInit = retrySmallInit,
      effective_tolerance = effective_tolerance,
      retryFailOnExhaustion = retryFailOnExhaustion,
      rxThreads = rxThreads
    )
    if (length(results) == 0) {
      break
    }

    resTable <- do.call(rbind, lapply(results, `[[`, "stats"))
    # Winner = most significant candidate (smallest p-value). Ties at
    # pchisqr == 0 (p-value underflow for strong candidates) are broken by the
    # most negative deltObjf (biggest OFV drop); see .pickForwardWinner().
    bestIdx <- .pickForwardWinner(resTable)
    bestRow <- resTable[bestIdx, , drop = FALSE]

    # -- verbose: show all candidate results -------------------------------
    if (verbose) {
      cli::cli_h2("Forward step {stepIdx}: all candidate results")
      v_cols <- intersect(
        c(
          "covar",
          "var",
          "shape",
          "objf",
          "deltObjf",
          "pchisqr",
          "AIC",
          "BIC",
          "bsvReduction",
          "covarEffect"
        ),
        names(resTable)
      )
      print(resTable[order(resTable$pchisqr), v_cols, drop = FALSE])
    }

    if (bestRow$pchisqr <= pVal) {
      resTable[bestIdx, "included"] <- "yes"
      bestRow[, "included"] <- "yes"

      cli::cli_h1("best model at step {stepIdx}:")
      print(bestRow)

      fit <- results[[bestIdx]]$fit

      # -- verbose: updated model code + parameter estimates + omega --------
      if (verbose) {
        cli::cli_h2("Accepted model -- updated model code:")
        tryCatch(
          print(fit$finalUiEnv),
          error = function(e) {
            cli::cli_alert_warning(
              "Model code unavailable: {conditionMessage(e)}"
            )
          }
        )
        cli::cli_h2("Accepted model -- fixed effects:")
        print(fit$parFixedDf)
        cli::cli_h2("Accepted model -- omega diagonal (BSV variances):")
        print(round(diag(fit$omega), 4))
      }

      acc_var <- results[[bestIdx]]$var
      acc_covar <- results[[bestIdx]]$covar

      # Track accepted pair for rebuild-from-base at next forward step
      acc_row <- pairs[
        pairs$var == acc_var & pairs$covar == acc_covar,
        ,
        drop = FALSE
      ]
      if (nrow(acc_row) > 0L) {
        accepted_pairs <- rbind(accepted_pairs, acc_row[1L, , drop = FALSE])
      }
      # Warm-start: update all accepted pairs' inits from the current fit so
      # that subsequent candidates use estimated rather than initial values.
      accepted_pairs <- .updatePairsInits(accepted_pairs, fit)

      # Remove the accepted pair AND all other shapes of the same raw covariate
      # on the same parameter -- only one parameterisation of a given
      # covariate-parameter relationship may enter the model.
      if ("raw_col" %in% names(pairs)) {
        acc_raw <- pairs$raw_col[
          pairs$var == acc_var & pairs$covar == acc_covar
        ]
        acc_raw <- if (length(acc_raw) > 0L) acc_raw[1L] else NA_character_
        if (!is.na(acc_raw)) {
          dropped <- pairs$covar[
            pairs$var == acc_var &
              pairs$raw_col == acc_raw &
              pairs$covar != acc_covar
          ]
          if (length(dropped) > 0L) {
            cli::cli_alert_info(
              "Dropping alternative shape(s) for {acc_raw}~{acc_var}: {paste(dropped, collapse=', ')}"
            )
          }
          pairs <- pairs[!(pairs$var == acc_var & pairs$raw_col == acc_raw), ]
        } else {
          pairs <- pairs[!(pairs$var == acc_var & pairs$covar == acc_covar), ]
        }
      } else {
        pairs <- pairs[!(pairs$var == acc_var & pairs$covar == acc_covar), ]
      }

      resTableComplete <- rbind(resTableComplete, resTable)
      if (saveModels) {
        key <- paste0(acc_covar, "_", acc_var)
        saveRDS(
          fit,
          paste0(outputDir, "/forward_step_", stepIdx, "_fit_", key, ".rds")
        )
        saveRDS(
          bestRow,
          paste0(outputDir, "/forward_step_", stepIdx, "_table_", key, ".rds")
        )
        saveRDS(
          resTableComplete,
          paste0(
            outputDir,
            "/forward_step_",
            stepIdx,
            "_completetable_",
            key,
            ".rds"
          )
        )
      }

      cli::cli_h2("accepted {acc_covar}~{acc_var}")
      stepIdx <- stepIdx + 1
    } else {
      cli::cli_h1("OFV did not improve, exiting forward search ...")
      resTableComplete <- rbind(resTableComplete, resTable)
      break
    }
  }

  cli::cli_h2(cli::col_red("forward search complete"))
  list(fit, resTableComplete, accepted_pairs)
}

#' Backward covariate search
#'
#' @param pairs       data frame with columns \code{var} and \code{covar}
#' @param fitorig     the original base fit object
#' @param fitupdated  the fit after forward search (if any)
#' @param data        data frame for fitting
#' @param pVal        p-value threshold for retaining a covariate
#' @param reFitCovars logical; if \code{TRUE} all pairs are added simultaneously
#'   before backward elimination (used for standalone backward search)
#' @param includedRelations enriched pairs data frame (output of the pipeline
#'   in \code{runSCM}) of relationships that must be present in the
#'   model before backward elimination begins.  Any that are missing are added
#'   and the model is re-fit once.  \code{NULL} = no forced inclusions.
#' @param outputDir   cache directory path
#' @param restart     logical; \code{FALSE} = resume from cache (only when
#'   \code{saveModels = TRUE})
#' @param saveModels  logical; if \code{FALSE} all \code{saveRDS}/\code{readRDS}
#'   operations are skipped
#' @param verbose     logical; if \code{TRUE} print candidate table before each
#'   step, full results table after fitting, and parameter/omega detail when a
#'   covariate is retained
#' @return \code{list(final_fit, resTableComplete)}
#' @noRd
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
backwardSearch <- function(
  pairs,
  fitorig,
  fitupdated,
  data,
  pVal = 0.01,
  reFitCovars = FALSE,
  context_pairs = NULL,
  includedRelations = NULL,
  outputDir,
  restart = FALSE,
  saveModels = TRUE,
  verbose = FALSE,
  control = NULL,
  print = 100,
  maxRetries = 3L,
  maxDeltaOFV = Inf,
  retryPerturbSD = 0.5,
  retrySmallInit = 0.01,
  retryOFVTolerance = NULL,
  retryFailOnExhaustion = FALSE,
  rxThreads = NULL
) {
  if (!inherits(fitorig, "nlmixr2FitCore")) {
    stop("'fitorig' needs to be a nlmixr2 fit")
  }
  if (saveModels && missing(outputDir)) {
    stop("please specify outputDir for backward search")
  }

  effective_tolerance <- .resolveOFVTolerance(fitorig, retryOFVTolerance)
  base_ui <- fitorig$finalUiEnv # fixed -- never updated
  resTableComplete <- NULL
  stepIdx <- 1

  # Determine starting model
  if (!missing(fitupdated)) {
    fwd_thetas <- fitupdated$iniDf$name[!is.na(fitupdated$iniDf$ntheta)]
    orig_thetas <- fitorig$iniDf$name[!is.na(fitorig$iniDf$ntheta)]
    if (all(fwd_thetas %in% orig_thetas)) {
      # No forward covariates accepted. Skip backward unless includedRelations
      # forces something into the model first.
      if (is.null(includedRelations) || nrow(includedRelations) == 0) {
        cli::cli_alert_warning(
          "no covariates added in the forward search, skipping backward search"
        )
        return(list(fitorig, NULL))
      }
      cli::cli_alert_info(
        "no forward covariates accepted; proceeding for includedRelations"
      )
    }
    fit <- fitupdated
  } else {
    fit <- fitorig
  }

  bck_control <- tryCatch(
    {
      ctrl <- if (!is.null(control)) control else fit$control
      ctrl$print <- print
      ctrl
    },
    error = function(e) NULL
  )

  if (reFitCovars) {
    xmod <- .rebuildUiFromPairs(base_ui, pairs) # nolint: object_usage_linter.
    fit <- suppressWarnings(nlmixr2est::nlmixr2(
      xmod,
      data,
      fit$est,
      control = bck_control
    ))
    context_pairs <- .pairsInModel(pairs, fit)
  }

  # -- includedRelations --------------------------------------------------------
  # Ensure all specified relations are in the model.  Relations already present
  # (from forward search) are left untouched.  Any that are missing are added
  # to context_pairs and the model is rebuilt from base and re-fit once before
  # the backward loop begins.  Included relations are then unioned into `pairs`
  # so that .pairsInModel() can find them and they are eligible for removal.
  if (!is.null(includedRelations) && nrow(includedRelations) > 0) {
    # Derive context_pairs now if not yet set (needed to rebuild below)
    if (is.null(context_pairs)) {
      context_pairs <- .pairsInModel(pairs, fit)
    }

    cur_thetas <- fit$iniDf$name[!is.na(fit$iniDf$ntheta)]
    inc_keys <- paste0(
      "cov_",
      includedRelations$covar,
      "_",
      includedRelations$var
    )
    missing_rel <- includedRelations[!inc_keys %in% cur_thetas, , drop = FALSE]

    if (nrow(missing_rel) > 0) {
      cli::cli_alert_info(
        "Adding {nrow(missing_rel)} includedRelation(s) not in current model: {paste(missing_rel$covar, missing_rel$var, sep='~', collapse=', ')}" # nolint: line_length_linter.
      )

      # Merge missing relations into context_pairs and rebuild from base
      new_context <- rbind(context_pairs, missing_rel)
      cp_key <- paste0(new_context$var, "_", new_context$covar)
      context_pairs <- new_context[!duplicated(cp_key), , drop = FALSE]

      xmod <- .rebuildUiFromPairs(base_ui, context_pairs) # nolint: object_usage_linter.
      fit <- suppressWarnings(nlmixr2est::nlmixr2(
        xmod,
        data,
        fit$est,
        control = bck_control
      ))
    }

    # Union includedRelations into the candidate pairs so .pairsInModel()
    # picks them up.  De-duplicate on (covar, var) key.
    #
    # NOTE: `pairs$covar` and `includedRelations$covar` have BOTH already
    # been rewritten by .expandShapes() / categorical expansion to include
    # the shape / level suffix (e.g. "BW_power", "BW_log", "SEX_1").
    # `inc_keys` is built with the "cov_" prefix to compare against
    # `cur_thetas` (which are theta names like "cov_BW_power_cl").  The
    # `pair_keys` MUST therefore also carry the "cov_" prefix or no row
    # will ever be considered a duplicate, and every includedRelation
    # gets appended on top of the matching pair already in `pairs`,
    # silently doubling every candidate at every backward step.
    pair_keys <- paste0("cov_", pairs$covar, "_", pairs$var)
    new_rows <- includedRelations[!inc_keys %in% pair_keys, , drop = FALSE]
    if (nrow(new_rows) > 0) {
      pairs <- rbind(pairs, new_rows)
    }
  }

  # Restrict to pairs that are actually in the current model
  pairs <- .pairsInModel(pairs, fit)
  if (nrow(pairs) == 0) {
    cli::cli_alert_warning(
      "no candidate pairs found in the model, skipping backward search"
    )
    return(list(fit, NULL))
  }

  # Derive context_pairs from the model if not yet set
  if (is.null(context_pairs)) {
    context_pairs <- pairs
  }

  # When not saving there is no cache to read; treat as a clean start
  if (!saveModels) {
    restart <- TRUE
  }

  if (!restart) {
    fileExistsTab <- list.files(
      paste0("./", outputDir),
      pattern = "backward_step_[0-9]+_table_[a-zA-Z0-9_]+\\.rds"
    )
    fileExistsFit <- list.files(
      paste0("./", outputDir),
      pattern = "backward_step_[0-9]+_fit_[a-zA-Z0-9_]+\\.rds"
    )
    fileExistsCT <- list.files(
      paste0("./", outputDir),
      pattern = "backward_step_[0-9]+_completetable_[a-zA-Z0-9_]+\\.rds"
    )

    if (length(fileExistsTab) == 0) {
      restart <- TRUE
    } else {
      resumeRows <- lapply(fileExistsTab, function(f) {
        readRDS(paste0(outputDir, "/", f))
      })
      resumeTable <- data.table::rbindlist(resumeRows)
      fit <- readRDS(paste0(
        outputDir,
        "/",
        fileExistsFit[[length(fileExistsFit)]]
      ))
      context_pairs <- .pairsInModel(pairs, fit)
      resTableComplete <- readRDS(paste0(
        outputDir,
        "/",
        fileExistsCT[[length(fileExistsCT)]]
      ))
      stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1
      for (i in seq_len(nrow(resumeTable))) {
        pairs <- pairs[
          !(pairs$var == resumeTable$var[i] &
            pairs$covar == resumeTable$covar[i]),
        ]
      }
      cli::cli_alert_success(
        "loaded backward search data from disk, resuming ..."
      )
    }
  }

  cli::cli_h1("starting backward search...")

  while (nrow(pairs) > 0) {
    # -- verbose: current model code + candidates before fitting -----------
    if (verbose) {
      cli::cli_h2("Backward step {stepIdx}: current model")
      tryCatch(
        print(fit$finalUiEnv),
        error = function(e) {
          cli::cli_alert_warning(
            "Model code unavailable: {conditionMessage(e)}"
          )
        }
      )
      cli::cli_h2("Backward step {stepIdx}: {nrow(pairs)} candidate(s) to test")
      v_cols <- intersect(
        c("var", "covar", "shape", "init", "covExpr"),
        names(pairs)
      )
      print(pairs[, v_cols, drop = FALSE])
    }

    results <- .fitCandidatePairs(
      pairs,
      base_ui,
      context_pairs,
      fit,
      data,
      pVal,
      stepIdx,
      add = FALSE,
      control = control,
      print = print,
      maxRetries = maxRetries,
      maxDeltaOFV = maxDeltaOFV,
      retryPerturbSD = retryPerturbSD,
      retrySmallInit = retrySmallInit,
      effective_tolerance = effective_tolerance,
      retryFailOnExhaustion = retryFailOnExhaustion,
      rxThreads = rxThreads
    )
    if (length(results) == 0) {
      break
    }

    resTable <- do.call(rbind, lapply(results, `[[`, "stats"))
    # Backward: drop the candidate whose removal is LEAST significant (highest
    # p-value). Ties at pchisqr == 1 (dOFV <= 0) are broken by the smallest
    # deltObjf; see .pickBackwardWinner().
    bestIdx <- .pickBackwardWinner(resTable)
    bestRow <- resTable[bestIdx, , drop = FALSE]

    # -- verbose: show all candidate results -------------------------------
    if (verbose) {
      cli::cli_h2("Backward step {stepIdx}: all candidate results")
      v_cols <- intersect(
        c(
          "covar",
          "var",
          "shape",
          "objf",
          "deltObjf",
          "pchisqr",
          "AIC",
          "BIC",
          "bsvReduction",
          "covarEffect"
        ),
        names(resTable)
      )
      print(resTable[
        order(resTable$pchisqr, decreasing = TRUE),
        v_cols,
        drop = FALSE
      ])
    }

    # Remove the covariate only when its OFV impact is non-significant
    # (p > pVal means the OFV increase from removal does not exceed threshold).
    if (bestRow$pchisqr > pVal) {
      # Use "dropped" (not "yes") so the raw summaryTable is self-documenting:
      # forward+"yes" = added; backward+"dropped" = removed.
      resTable[bestIdx, "included"] <- "dropped"
      bestRow[, "included"] <- "dropped"

      cli::cli_h1("removing covariate at step {stepIdx}:")
      print(bestRow)

      # fit is now the reduced model (with the removed covariate dropped)
      fit <- results[[bestIdx]]$fit

      # -- verbose: updated model code + parameter estimates + omega --------
      if (verbose) {
        cli::cli_h2("Reduced model -- updated model code:")
        tryCatch(
          print(fit$finalUiEnv),
          error = function(e) {
            cli::cli_alert_warning(
              "Model code unavailable: {conditionMessage(e)}"
            )
          }
        )
        cli::cli_h2("Reduced model -- fixed effects:")
        print(fit$parFixedDf)
        cli::cli_h2("Reduced model -- omega diagonal (BSV variances):")
        print(round(diag(fit$omega), 4))
      }

      ret_var <- results[[bestIdx]]$var
      ret_covar <- results[[bestIdx]]$covar
      pairs <- pairs[!(pairs$var == ret_var & pairs$covar == ret_covar), ]

      # Update context_pairs to reflect the covariate that was removed
      context_pairs <- context_pairs[
        !(context_pairs$var == ret_var & context_pairs$covar == ret_covar),
        ,
        drop = FALSE
      ]
      # Warm-start remaining pairs from the reduced model's estimates
      context_pairs <- .updatePairsInits(context_pairs, fit)

      resTableComplete <- rbind(resTableComplete, resTable)
      if (saveModels) {
        key <- paste0(ret_covar, "_", ret_var)
        saveRDS(
          fit,
          paste0(outputDir, "/backward_step_", stepIdx, "_fit_", key, ".rds")
        )
        saveRDS(
          bestRow,
          paste0(outputDir, "/backward_step_", stepIdx, "_table_", key, ".rds")
        )
        saveRDS(
          resTableComplete,
          paste0(
            outputDir,
            "/backward_step_",
            stepIdx,
            "_completetable_",
            key,
            ".rds"
          )
        )
      }

      cli::cli_h2("removed {ret_covar}~{ret_var}")
      stepIdx <- stepIdx + 1
    } else {
      cli::cli_h1(
        "All remaining covariates significant (p <= {pVal}), backward search complete."
      )
      resTableComplete <- rbind(resTableComplete, resTable)
      break
    }
  }

  cli::cli_h2(cli::col_red("backward search complete"))
  list(fit, resTableComplete)
}
