# OSMSS package installer
# Installs the R packages OSMSS needs. Packages that are already installed are
# left as they are (no updates). Writes a short report for the installer wizard.
#
# Usage:  Rscript install_packages.R <report_file>
# The same package list is given as a one-line command in README.md (Installation, step 2).

args   <- commandArgs(trailingOnly = TRUE)
report <- if (length(args) >= 1) args[1] else file.path(tempdir(), "osmss_install_report.txt")

say <- function(...) { cat(..., "\n", sep = ""); flush.console() }

pkgs <- c("shiny", "shinyjs", "bslib", "DT", "plotly", "sortable", "htmltools", "jsonlite", "xml2",
          "httr", "base64enc", "dplyr", "tidyr", "readr", "stringr", "glue", "scales", "data.table",
          "future", "promises", "ggplot2", "ggvenn", "pheatmap", "e1071", "MASS", "statmod", "ruv",
          "bio3d", "NGLVieweR", "msigdbr", "rstudioapi", "limma", "proDA", "QFeatures", "MsCoreUtils",
          "SummarizedExperiment", "S4Vectors", "impute", "PhosR", "Biostrings", "msa", "drawProteins",
          "clusterProfiler")
# Packages from other repositories, as name = repository URL (none needed at present;
# msigdbr downloads the MSigDB gene sets itself, so msigdbdf is not required).
extra <- character(0)

# ---- 1. A writable library ------------------------------------------------------------
# Rscript cannot answer R's "use a personal library?" question, so create it up front.
writable <- function(d) dir.exists(d) && file.access(d, 2) == 0
lib <- Filter(writable, .libPaths())[1]
if (is.na(lib) || is.null(lib)) {
  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
}
say("R ", as.character(getRversion()), " | installing into: ", lib)

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = max(600, getOption("timeout")),
        Ncpus = max(1L, parallel::detectCores(logical = FALSE) - 1L, na.rm = TRUE))   # parallel source builds (Linux)
is_installed <- function(p) nzchar(system.file(package = p))

all_pkgs <- c(pkgs, names(extra))
pre      <- all_pkgs[vapply(all_pkgs, is_installed, logical(1))]
missing  <- setdiff(pkgs, pre)
say(length(pre), " of ", length(all_pkgs), " packages already installed; ",
    length(missing) + sum(!names(extra) %in% pre), " to install.")

# ---- 2. Install what is missing -------------------------------------------------------
errors <- character(0)
if (!is_installed("BiocManager")) {
  say("Installing BiocManager ...")
  tryCatch(install.packages("BiocManager", lib = lib),
           error = function(e) errors <<- c(errors, paste("BiocManager:", conditionMessage(e))))
}
if (length(missing) > 0 && is_installed("BiocManager")) {
  say("Installing ", length(missing), " packages from CRAN and Bioconductor (this can take 10-30 minutes) ...")
  tryCatch(BiocManager::install(missing, lib = lib, update = FALSE, ask = FALSE),
           error = function(e) errors <<- c(errors, conditionMessage(e)))
}
for (p in names(extra)) if (!is_installed(p)) {
  say("Installing ", p, " ...")
  tryCatch(install.packages(p, lib = lib, repos = c(extra[[p]], "https://cloud.r-project.org")),
           error = function(e) errors <<- c(errors, paste0(p, ": ", conditionMessage(e))))
}

# ---- 3. Check every package loads ---------------------------------------------------
say("Checking that each package loads ...")
not_installed <- all_pkgs[!vapply(all_pkgs, is_installed, logical(1))]
no_load <- character(0)
for (p in setdiff(all_pkgs, not_installed)) {
  ok <- tryCatch(suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE)), error = function(e) FALSE)
  if (!isTRUE(ok)) no_load <- c(no_load, p)
}
new_inst <- setdiff(setdiff(all_pkgs, pre), not_installed)
failed   <- union(not_installed, no_load)

out <- c(paste0("R_VERSION\t", getRversion()),
         paste0("LIBRARY\t", lib),
         paste0("PREINSTALLED\t", paste(pre, collapse = ",")),
         paste0("INSTALLED\t", paste(new_inst, collapse = ",")),
         paste0("FAILED\t", paste(failed, collapse = ",")),
         paste0("ERRORS\t", gsub("[\r\n\t]+", " ", paste(errors, collapse = " | "))),
         paste0("STATUS\t", if (length(failed) == 0) "OK" else "FAILED"))
writeLines(out, report)
say(if (length(failed) == 0) "All packages are installed and load." else
      paste0("Packages that could not be installed or loaded: ", paste(failed, collapse = ", ")))
say("DONE")
