# Precompute the runSCM vignette.
#
# runSCM.Rmd.orig fits real models and runs a full forward/backward search,
# which takes several minutes -- far too long to re-run inside `R CMD check`.
# Knitting it here bakes the results into runSCM.Rmd as static output, so the
# shipped vignette is plain markdown that pandoc renders in seconds.
#
# Run this from the package root whenever runSCM.Rmd.orig changes:
#
#   Rscript vignettes/precompute.R
#
# Both this script and runSCM.Rmd.orig are excluded from the build via
# .Rbuildignore; only the generated runSCM.Rmd ships.

# Make sure we are knitting against the version being developed, not whatever
# older nlmixr2scm happens to be installed -- otherwise chunk errors get baked
# into the shipped vignette.
if (utils::packageVersion("nlmixr2scm") != read.dcf("DESCRIPTION", "Version")[1]) {
  stop(
    "Installed nlmixr2scm (", utils::packageVersion("nlmixr2scm"),
    ") differs from the source tree (", read.dcf("DESCRIPTION", "Version")[1],
    "). Install the current source before precomputing.",
    call. = FALSE
  )
}

old <- setwd("vignettes")

# Abort on a chunk error rather than embedding it in the output.
knitr::opts_chunk$set(error = FALSE)

knitr::knit("runSCM.Rmd.orig", output = "runSCM.Rmd")

# print()/summary() report the absolute output directory; replace the local
# working-directory prefix so no machine-specific path ships in the vignette.
rmd <- readLines("runSCM.Rmd", encoding = "UTF-8")
for (prefix in paste0(normalizePath(getwd()), c("\\", "/"))) {
  rmd <- gsub(prefix, "/path/to/project/", rmd, fixed = TRUE)
}
writeLines(rmd, "runSCM.Rmd", useBytes = TRUE)

# The vignette demonstrates saveModels = TRUE, which writes a fitted-model
# cache into the working directory.  Those artifacts are a by-product of
# precomputation, not vignette source, so drop them.
unlink(list.files(".", pattern = "_scm_[0-9]+(_backup_.*)?$"), recursive = TRUE)
unlink(c("scm_log.txt", "scm_step_summary.csv", "scm_all_candidates.csv"))

setwd(old)
