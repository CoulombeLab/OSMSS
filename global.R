# global.R =====================================================================
# One Stop MS Shop (OSMSS) - Residue PTM Prediction + Local MS Integration
# Shared libraries, constants, and helper functions used across all modules.
# ============================================================================

# OSMSS release version: shown on the welcome screen and printed in the methods text.
OSMSS_VERSION <- "3.1"

# === PERSISTENT FILE LOGGING ==================================================
# Writes structured log entries to logs/ directory. Survives crashes.
# Each day gets its own file: logs/osmss_YYYY-MM-DD.log
# Format: [YYYY-MM-DD HH:MM:SS] [LEVEL] [MODULE] [SESSION] message
# ==============================================================================
LOG_DIR <- normalizePath(file.path(getwd(), "logs"), winslash = "/", mustWork = FALSE)

# Ensure logs directory exists at startup
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)

#' Write a single structured entry to the daily log file
#' @param msg     Character message to log
#' @param level   Severity: "INFO", "WARN", "ERROR"
#' @param module  Module name (e.g., "DE", "PCA", "TMT_Volcano")
#' @param session_id  Short session token (first 8 chars of session$token)
osmss_file_log <- function(msg, level = "INFO", module = "System",
                           session_id = "--------") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_file  <- file.path(LOG_DIR, paste0("osmss_", Sys.Date(), ".log"))
  entry     <- paste0("[", timestamp, "] [", level, "] [", module, "] [",
                       session_id, "] ", msg)
  tryCatch({
    # Purge .log and .txt files in LOG_DIR that are >= 2 days old
    cutoff <- Sys.time() - 2 * 86400
    old_files <- list.files(LOG_DIR, pattern = "\\.(log|txt)$",
                            full.names = TRUE, recursive = FALSE)
    old_files <- old_files[file.mtime(old_files) < cutoff]
    if (length(old_files) > 0) file.remove(old_files)

    if (!file.exists(log_file)) {
      file.create(log_file, showWarnings = FALSE)
      Sys.chmod(log_file, mode = "0664")
    }
    cat(entry, "\n", file = log_file, append = TRUE, sep = "")
  }, error = function(e) message("LOG WRITE FAILED: ", e$message))
}

#' Auto-detect severity & module from legacy logger messages.
#' Many existing calls use prefixes like "!!! [TKL] ..." or ">>> [EVO] ...".
#' This parses those into structured (level, module, clean_msg) fields.
osmss_parse_log_msg <- function(txt) {
  level  <- "INFO"
  module <- "System"
  msg    <- txt

  # Detect severity from prefix

  if (grepl("^!!!", txt))           level <- "ERROR"
  else if (grepl("^\\?\\?\\?", txt)) level <- "WARN"
  else if (grepl("^>>>", txt))      level <- "INFO"
  else if (grepl("^\\.\\.\\.\\.", txt)) level <- "DEBUG"

  # Extract [MODULE] tag if present, e.g., ">>> [TKL] Prediction..."
  m <- regmatches(txt, regexpr("\\[([A-Za-z0-9_-]+)\\]", txt))
  if (length(m) == 1 && nchar(m) > 0) {
    module <- gsub("\\[|\\]", "", m)
  }

  # Strip leading markers for a cleaner file message

  msg <- sub("^(!{3}|>{3}|\\?{3}|\\.{3})\\s*", "", msg)
  msg <- sub("^\\[[A-Za-z0-9_-]+\\]\\s*", "", msg)

  list(level = level, module = module, msg = msg)
}

# Log app startup
osmss_file_log("=== OSMSS Application Starting ===", level = "INFO", module = "Startup")

# === LIBRARIES ================================================================
options(shiny.maxRequestSize = 1024 * 1024^2)  # 1 GB upload limit (MaxQuant evidence.txt); keep in sync with app.R
library(shiny)
library(shinyjs)
library(bslib)
# library(ggplot2)    # Removed: plotly imports ggplot2 as a dependency
library(jsonlite)
# library(ggrepel)    # Removed: no longer used
library(dplyr)
library(readr)
library(stringr)
# library(ggnewscale) # Removed: no longer used
# library(ggtext)      # Removed: no longer used
library(scales)
library(httr)
library(tidyr)
library(plotly)
library(DT)
library(xml2)
library(sortable)
library(glue)
library(data.table)  # fread for high-performance file I/O
library(future)      # async background workers (ExtendedTask)
library(promises)    # promise-based async for Shiny ExtendedTask
# Bioconductor  -" ProDA / QFeatures (Phase 2 DE engine)
# Install via: BiocManager::install(c("proDA", "QFeatures", "MsCoreUtils"))
#
# NOTE: Use system.file() instead of requireNamespace() for availability checks.
# requireNamespace() loads the full package namespace (~4-5 s for proDA,
# ~3.5 s for MSstats) even though we only need to know if it is installed.
# system.file() is a pure filesystem check (~0 ms) and avoids that cost.
.pkg_installed <- function(pkg) nchar(system.file(package = pkg)) > 0L
.PRODA_AVAILABLE <- .pkg_installed("proDA")       &&
                    .pkg_installed("QFeatures")    &&
                    .pkg_installed("MsCoreUtils")  &&
                    .pkg_installed("SummarizedExperiment")
.IMPUTE_AVAILABLE  <- .pkg_installed("impute")

# PhosR  -" reference data for Kinase-Site Analysis (motifs, PSP substrates, kinase families),
# extracted on first use by ensure_phosr_reference(). Install via: BiocManager::install("PhosR")
.PHOSR_AVAILABLE <- .pkg_installed("PhosR")
if (!.PHOSR_AVAILABLE)
  message(">>> SYSTEM: PhosR not installed. Kinase-Site Analysis will not be available. ",
          "Install: BiocManager::install('PhosR')")

# SAINTexpress  -" AP-MS interaction scoring (intensity mode)
.SAINT_AVAILABLE <- {
  saint_dir <- file.path("SAINT", "SAINTexpress_v3.6.3__2018-03-09", "bin")
  if (.Platform$OS.type == "windows") {
    file.exists(file.path(saint_dir, "SAINTexpress-int.exe"))
  } else {
    file.exists(file.path(saint_dir, "SAINTexpress-int"))
  }
}

.saint_int_binary <- function() {
  saint_dir <- file.path("SAINT", "SAINTexpress_v3.6.3__2018-03-09", "bin")
  if (.Platform$OS.type == "windows") {
    file.path(saint_dir, "SAINTexpress-int.exe")
  } else {
    file.path(saint_dir, "SAINTexpress-int")
  }
}

# Invariant-protein median normalization (ip_sub = "invariant"): the ProDA
# normalization of every MNAR-aware route (AP-MS fractionated / cross-condition /
# trap control with a factorial design, and non-IP "proteins entirely absent").
# .INVARIANT_MIN_REFS        -" hard minimum panel size; below this the pipeline
#                             falls back to global median normalization.
# .INVARIANT_MIN_REFS_IDEAL  -" preferred panel size; between MIN and IDEAL the
#                             panel is used with a low-confidence warning.
# .INVARIANT_PEAK_RATIO  -" a second density peak at least this tall relative
#                         to the main one (and > 2 bandwidths away) makes the
#                         mode ambiguous: the centre is the candidate peak
#                         nearest zero.
# .INVARIANT_CENTRE_GAP  -" mode and median centres further apart than this
#                         (log2, between-condition scale) are reported.
.INVARIANT_MIN_REFS       <- 10L
.INVARIANT_MIN_REFS_IDEAL <- 20L
.INVARIANT_PEAK_RATIO <- 0.5
.INVARIANT_CENTRE_GAP <- 0.1

# === ASYNC WORKER BACKEND =====================================================
# ExtendedTask workers: 3 background R sessions for concurrent heavy tasks.
# Each worker is a full R process (~200-400 MB). Increase if server RAM allows.
tryCatch(
  {
    future::plan(future::multisession, workers = 3L)
    message(">>> SYSTEM: future::multisession backend initialized (workers=3)")
  },
  error = function(e) {
    message(">>> SYSTEM WARNING: future::multisession failed (", e$message,
            ") - falling back to sequential. Async tasks will block the session.")
    future::plan(future::sequential)
  }
)

.PRODA_LOADED <- FALSE
if (.PRODA_AVAILABLE) {
  # Defer library() calls until the user actually selects ProDA.
  # This shaves ~5-8 s off app startup for the majority of sessions
  # that never touch ProDA.
  ensure_proda_loaded <- function() {
    if (!.PRODA_LOADED) {
      suppressPackageStartupMessages({
        library(proDA)
        library(QFeatures)
        library(MsCoreUtils)
        library(SummarizedExperiment)
      })
      .PRODA_LOADED <<- TRUE
      message(">>> SYSTEM: ProDA + QFeatures loaded (on first use).")
    }
    invisible(TRUE)
  }
  message(">>> SYSTEM: ProDA + QFeatures available (will load on first use).")
} else {
  ensure_proda_loaded <- function() {
    stop("proDA / QFeatures not installed. ",
         "Install: BiocManager::install(c('proDA','QFeatures','MsCoreUtils'))",
         call. = FALSE)
  }
  message(">>> SYSTEM: proDA / QFeatures not installed. ",
          "Label-free DE analysis will not be available. ",
          "Install: BiocManager::install(c('proDA','QFeatures','MsCoreUtils'))")
}

# Null-coalescing operator (returns rhs when lhs is NULL)
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

# -----------------------------------------------------------------------------
# ptm_parent_acc_of()
# Parent accession from a PTM site key.  The three summarizers build site IDs
# differently, so the pattern must match the one that built the key:
#   pd_tmt    ACC_S32_<anything>
#   maxquant  ACC_S32_M1[_SharedInference]
#   default   ACC_S32[+S36]_Phospho...   (Spectronaut; '+'-aware, multi-site safe)
# Single source of truth for this parse - callers must not re-derive it.
# -----------------------------------------------------------------------------
ptm_parent_acc_of <- function(keys, software = NULL) {
  sw <- as.character(software %||% "spectronaut")
  if (identical(sw, "pd_tmt"))
    sub("_[A-Za-z][0-9]+_.*$", "", keys)
  else if (identical(sw, "maxquant"))
    sub("_[A-Za-z][0-9]+_M[0-9]+(_SharedInference)?$", "", keys)
  else
    sub("_(?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+_Phospho.*$", "", keys, perl = TRUE)
}

# -----------------------------------------------------------------------------
# ptm_has_global() / ptm_layers() / ptm_primary_layer()
# A DIA PTM cache is built in one of two modes, chosen on the first screen of the
# DIA PTM upload wizard and stored once as cache$ptm_mode:
#   "with_global" - phospho + global proteome: stoichiometry-adjusted, absolute PTM
#                   and protein layers; adjusted is the primary view
#   "ptm_only"    - phospho alone: the absolute PTM layer only
# A cache without ptm_mode is "with_global": caches saved before the mode existed,
# and TMT PTM caches, always carry a global proteome.  Consumers ask these helpers
# which layers exist rather than testing matrices for NULL.  Layer names are the
# dia_ptm_limma_de() result names.
# -----------------------------------------------------------------------------
ptm_has_global <- function(cache) !identical(cache$ptm_mode, "ptm_only")

ptm_layers <- function(cache)
  if (ptm_has_global(cache)) c("adjusted", "unadjusted", "protein") else "unadjusted"

ptm_primary_layer <- function(cache) ptm_layers(cache)[[1L]]

# Display label of a layer, as the PTM DE view selector shows it.
PTM_LAYER_LABELS <- c(adjusted   = "Adjusted to Total Protein",
                      unadjusted = "Absolute PTM",
                      protein    = "Protein View")

# Residue set from a PTM site key, e.g. "S44" or "S44+S48".
ptm_site_part_of <- function(key) {
  m <- regmatches(key, regexpr("(?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+(?=_Phospho)",
                               key, perl = TRUE))
  if (length(m) == 1L && nzchar(m)) m else ""
}

# -----------------------------------------------------------------------------
# ptm_cache_fingerprint()
# Cheap, deterministic identity for a DIA/TMT PTM summarization cache.
#
# Two uses, both about not showing one dataset's statistics against another's:
#   1. server.R watches this value and clears rv$dia_ptm_de_results whenever it
#      changes, so a stored DE result can never outlive the cache it was computed
#      on.  DE results are keyed by contrast label ("<test> vs <ref>"), which is
#      NOT unique across datasets.
#   2. The PTM DE .rds handlers stamp it on save and check it on load, so a bundle
#      saved against a different dataset is not published to the shared store.
#
# Deliberately built ONLY from fields that describe the dataset.  It must NOT
# change when server_ptm_de.R writes the k-NN memo (ptm_abs_matrix_imputed,
# global_matrix_imputed, impute_stats, ptm_impute_cap) back into the cache -
# that mutation fires the same reactive, and a fingerprint sensitive to it would
# wipe the store immediately after the run that filled it.
#
# Returns NA_character_ for a NULL cache, which callers must treat as
# "cannot verify".
# -----------------------------------------------------------------------------
ptm_cache_fingerprint <- function(cache) {
  if (is.null(cache)) return(NA_character_)
  .n <- function(x) if (is.null(x) || length(x) == 0L) "NA" else as.character(x[1])
  samples <- cache$matched_samples
  if (is.null(samples)) samples <- colnames(cache$ptm_abs_matrix)
  samples <- if (is.null(samples)) "NA" else paste(sort(as.character(samples)), collapse = ",")
  paste(
    .n(cache$software),
    .n(cache$n_sites_input),
    .n(cache$n_sites_output),
    .n(cache$n_sites_adjusted),
    .n(cache$prob_threshold),
    .n(nrow(cache$ptm_abs_matrix)),
    .n(nrow(cache$global_matrix)),
    samples,
    sep = "|"
  )
}

# Heavy packages  -" loaded on first use via requireNamespace() + pkg::fn()
# drawProteins, bio3d, msa, Biostrings, NGLVieweR, limma, statmod, PhosR (data only)

# Helper: demand-load a package on first use (deferred from startup)
# Uses library() to fully attach  -" needed for packages with .onAttach hooks (e.g., msa)
require_pkg <- function(pkg) {
  if (!pkg %in% .packages()) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(paste0("Package '", pkg, "' is required but not installed."), call. = FALSE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# === RELATIVE PATH CONFIGURATION =============================================
DATA_DIR  <- "data"
CACHE_DIR <- "cache"

# Residue Prediction paths
DEFAULT_ST_MATRIX   <- file.path(DATA_DIR, "matrices", "ser_thr", "ser_thr_all_log2_matrices.tsv")
DEFAULT_Y_MATRIX    <- file.path(DATA_DIR, "matrices", "tyrosine", "tyrosine_all_log2_matrices.tsv")
DEFAULT_ST_FAV      <- file.path(DATA_DIR, "matrices", "ser_thr", "st_favorability.tsv")
DEFAULT_BG_DIR      <- file.path(DATA_DIR, "Background")
PATH_KINASE_FAMILIES <- file.path(DATA_DIR, "Kinase_Families.tsv")  # built by data/Kinase_Families_build.R (CC BY 4.0 sources)
PATH_KINASE_SUB     <- file.path(DATA_DIR, "Kinase_Substrate_Dataset", "Kinase_Substrate_Dataset")
PATH_SITES          <- file.path(DATA_DIR, "Phosphorylation_site_dataset", "Phosphorylation_site_dataset")

# Attribution for external web services (shown above the plots and in plot titles)
ATTR_SCANSITE <- "Data generated using the Scansite 4 public API (PMID: 12824383)"
ATTR_1433     <- "Data generated using the 14-3-3-Pred public API (PMID: 25735772)"

# === USER-DOWNLOADED DATA (not redistributed with OSMSS; CC BY-NC-SA) =========
# The Kinase Library matrices and the PhosphoSitePlus datasets are downloaded by
# each user from the original source. Modules that need them check first and
# show missing_data_modal() if they are absent.
TKL_REPO_URL     <- "https://github.com/TheKinaseLibrary/kinase-library"
TKL_REPO_DIR     <- "src/kinase_library/databases/matrices"
PSP_DOWNLOAD_URL <- "https://www.phosphosite.org/staticDownloads"

# PSP files may be kept unzipped (as downloaded, no extension), as .txt, or as .gz
resolve_psp_file <- function(path) {
  cand <- c(path, paste0(path, ".txt"), paste0(path, ".gz"))
  hit  <- cand[file.exists(cand)]
  if (length(hit) == 0) NULL else hit[1]
}

# Reads a PSP table whether or not it still has PSP's notice lines above the header
read_psp_table <- function(path) {
  f <- resolve_psp_file(path)
  if (is.null(f)) return(NULL)
  hdr <- grep("^GENE\t", readLines(f, n = 20, warn = FALSE))[1]
  if (is.na(hdr)) hdr <- 1L
  read.delim(f, skip = hdr - 1L, header = TRUE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

# === PhosR REFERENCE DATA (KS scoring) =========================================
# Not shipped with OSMSS: built on first use from the user's installed PhosR
# package (Bioconductor; BiocManager::install("PhosR")). Replaces the former
# shipped data/KSEA/phosr_*.rds files.
PHOSR_REF_DIR <- file.path(DATA_DIR, "KSEA")
phosr_ref_files <- function() {
  stats::setNames(file.path(PHOSR_REF_DIR, c("phosr_kinase_family.rds", "phosr_motif_human.rds",
                                             "phosr_motif_mouse.rds", "phosr_substrate_human.rds",
                                             "phosr_substrate_mouse.rds", "phosr_motif_weights_mouse_native.rds")),
                  c("kinase_family", "motif_human", "motif_mouse", "substrate_human", "substrate_mouse",
                    "weights_mouse_native"))
}

ensure_phosr_reference <- function(logger = message) {
  f <- phosr_ref_files()
  if (all(file.exists(f))) return(TRUE)
  if (!requireNamespace("PhosR", quietly = TRUE)) {
    logger("!!! [PhosR] Package not installed. Install with BiocManager::install('PhosR').")
    return(FALSE)
  }
  logger(">>> [PhosR] Building reference files from the installed PhosR package (first use)...")
  tryCatch({
    e <- new.env(parent = emptyenv())
    for (d in c("KinaseFamily", "KinaseMotifs", "PhosphoSitePlus"))
      utils::data(list = d, package = "PhosR", envir = e)
    # Mouse PWMs carry the human motif weights (NumInputSeq) by default, as in earlier
    # OSMSS versions; PhosR's own mouse weights are kept separately as a user option.
    motif_mouse <- e$motif.mouse.list
    weights_mouse_native <- motif_mouse$NumInputSeq
    motif_mouse$NumInputSeq <- e$motif.human.list$NumInputSeq
    dir.create(PHOSR_REF_DIR, showWarnings = FALSE, recursive = TRUE)
    saveRDS(e$KinaseFamily,      f[["kinase_family"]])
    saveRDS(e$motif.human.list,  f[["motif_human"]])
    saveRDS(motif_mouse,         f[["motif_mouse"]])
    saveRDS(e$PhosphoSite.human, f[["substrate_human"]])
    saveRDS(e$PhosphoSite.mouse, f[["substrate_mouse"]])
    saveRDS(weights_mouse_native, f[["weights_mouse_native"]])
    logger(">>> [PhosR] Reference files ready.")
    all(file.exists(f))
  }, error = function(err) {
    logger(paste0("!!! [PhosR] Could not build reference files: ", conditionMessage(err)))
    FALSE
  })
}

tkl_data_available <- function() file.exists(DEFAULT_ST_MATRIX) && file.exists(DEFAULT_Y_MATRIX)
psp_data_available <- function() !is.null(resolve_psp_file(PATH_SITES)) && !is.null(resolve_psp_file(PATH_KINASE_SUB))

# Modal telling the user where to download the missing data and where to put it
missing_data_modal <- function(need_tkl = FALSE, need_psp = FALSE) {
  data_abs <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
  # Create the destination folders so users only need to drop the files in
  dest <- c(if (need_tkl) dirname(c(DEFAULT_ST_MATRIX, DEFAULT_Y_MATRIX)),
            if (need_psp) dirname(c(PATH_SITES, PATH_KINASE_SUB)))
  for (d in unique(dest)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  tkl_part <- if (need_tkl) tagList(
    h5(icon("book"), " The Kinase Library (TKL)"),
    p("Download these three files from the ",
      tags$a(href = TKL_REPO_URL, target = "_blank", "Kinase Library GitHub repository"),
      " (folder ", tags$code(TKL_REPO_DIR), "):"),
    tags$ul(
      tags$li(tags$code("ser_thr/ser_thr_all_log2_matrices.tsv"), " and ",
              tags$code("ser_thr/st_favorability.tsv"), " → ",
              tags$code(file.path(data_abs, "matrices", "ser_thr"))),
      tags$li(tags$code("tyrosine/tyrosine_all_log2_matrices.tsv"), " → ",
              tags$code(file.path(data_abs, "matrices", "tyrosine")))
    ),
    p(class = "text-muted", style = "font-size: 0.85em;",
      "The Kinase Library data are licensed CC BY-NC-SA 3.0 (non-commercial use).")
  )
  psp_part <- if (need_psp) tagList(
    h5(icon("atom"), " PhosphoSitePlus (PSP)"),
    p("Download these two datasets from the ",
      tags$a(href = PSP_DOWNLOAD_URL, target = "_blank", "PhosphoSitePlus downloads page"),
      " (free registration required), and unzip them:"),
    tags$ul(
      tags$li(tags$code("Phosphorylation_site_dataset"), " → ",
              tags$code(file.path(data_abs, "Phosphorylation_site_dataset"))),
      tags$li(tags$code("Kinase_Substrate_Dataset"), " → ",
              tags$code(file.path(data_abs, "Kinase_Substrate_Dataset")))
    ),
    p(class = "text-muted", style = "font-size: 0.85em;",
      "PhosphoSitePlus data are licensed CC BY-NC-SA 3.0 (non-commercial use).")
  )
  modalDialog(
    title = tagList(icon("triangle-exclamation"), " Reference data needed"),
    p("This analysis uses reference data that OSMSS cannot redistribute. ",
      "Download the files below and place them in your OSMSS data folder for full functionality."),
    tkl_part, psp_part,
    if (need_tkl) p(class = "text-muted", style = "font-size: 0.85em;",
      "Restart OSMSS after adding the Kinase Library files."),
    easyClose = TRUE, footer = modalButton("Close")
  )
}

# Local MS / FragPipe paths
FRAGPIPE_DB_FILE <- file.path(DATA_DIR, "FragpipeDB", "FragpipeDB.csv")

# === GLOBAL CONSTANTS =========================================================
BINDER_TERMS <- c("SH2", "PTB", "14-3-3", "WW", "WD40", "FHA", "BRCT",
                  "Polo-box", "MH2", "PDZ", "C2", "Chromodomain", "Bromodomain",
                  "Grb2", "Crk", "Nck", "Shc", "Plcg", "Pik3", "Stat", "SOCS", "Cbl")

GLOBAL_RESOURCES <- list(ST_MAT=NULL, Y_MAT=NULL, BIAS_S=NULL, BIAS_T=NULL, BG_DIR=NULL, LOADED=FALSE)

# Lazy accessor: loads TKL matrices on first call, then returns cached result
ensure_global_resources <- function() {
  if (!GLOBAL_RESOURCES$LOADED) {
    message(">>> SYSTEM: Lazy-loading TKL Matrices into Memory...")
    res <- load_global_resources(DEFAULT_ST_MATRIX, DEFAULT_Y_MATRIX, DEFAULT_ST_FAV, DEFAULT_BG_DIR)
    if (!is.null(res) && res$LOADED) {
      GLOBAL_RESOURCES <<- res
      message(">>> SYSTEM: TKL Matrices Ready.")
    } else {
      message("!!! SYSTEM: Failed to load default matrices. Check 'data' folder structure.")
    }
  }
  invisible(GLOBAL_RESOURCES$LOADED)
}

# === SHARED PD / TMT CSV HELPERS =============================================
# Centralised parsing functions used by TMT Volcano, TMT Trends, and any
# future module that reads Proteome Discoverer multi-plex CSV exports.
# Defined once here to avoid dangerous duplication.
# =============================================================================

# Collapse multi-line CSV fields  -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -"
# PD exports can wrap GO/Pfam annotations across multiple lines inside quoted
# fields. This merges those continuation lines back into the opening line by
# tracking open-quote state via quote-parity.
collapse_multiline_csv <- function(all_lines) {
  # Fast path: if no line has an odd number of quotes, no multi-line fields
  n_quotes <- nchar(gsub('[^"]', '', all_lines))
  if (!any(n_quotes %% 2L == 1L)) return(all_lines)

  # Slow path: pre-allocate list to avoid O(n^2) vector growth
  result   <- vector("list", length(all_lines))
  j        <- 0L
  buffer   <- ""
  in_quote <- FALSE
  for (ln in all_lines) {
    if (in_quote) {
      buffer <- paste0(buffer, " ", ln)
    } else {
      if (nchar(buffer) > 0L) { j <- j + 1L; result[[j]] <- buffer }
      buffer <- ln
    }
    n_q <- nchar(gsub('[^"]', '', ln))
    if (n_q %% 2L == 1L) in_quote <- !in_quote
  }
  if (nchar(buffer) > 0L) { j <- j + 1L; result[[j]] <- buffer }
  unlist(result[seq_len(j)])
}

# Protein-row detection in a Proteome Discoverer nested CSV  -" -" -" -" -" -" -" -" -" -" -"
# A line is a protein row when the field in the "FDR Confidence" column holds
# High, Medium or Low. Fields are counted the CSV way: a quoted field may contain
# commas ("DUR1,2") and doubled quotes, and counts as one field. The previous rule
# counted raw commas, so a quoted comma in any field before the FDR Confidence
# column shifted the count and the protein was silently dropped.
# For lines without double-quote characters the result is identical to the old
# rule. The column index is taken from the full header (blank header cells
# included), because it is applied to raw lines; the old rule took it after
# blank header names were removed, which mis-counted when a blank preceded it.
# Without an FDR Confidence column the old fallback (second field) is kept.
# Returns the line numbers (> 1) of protein rows.
pd_protein_row_idx <- function(all_lines, header_line = all_lines[1]) {
  hdr_all <- as.character(
    read.csv(text = header_line, header = FALSE,
             stringsAsFactors = FALSE, check.names = FALSE)
  )
  conf_col_idx <- which(grepl("FDR Confidence", hdr_all))[1]
  if (is.na(conf_col_idx)) conf_col_idx <- 2L
  .field     <- '(?:"(?:[^"]|"")*"|[^,"]*)'   # quoted field, or unquoted field
  conf_regex <- paste0("^(?:", .field, ",){", conf_col_idx - 1L, "}",
                       '"?(?:High|Medium|Low)"?,')
  idx <- grep(conf_regex, all_lines, perl = TRUE)
  idx[idx > 1]
}

# Parse the structure of a single PD TMT CSV  -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -"
# Returns a list with: name, all_lines, header_line, headers, protein_idx,
# pep_hdr_idx, pep_data_idx, abund_cols, ratio_cols, pval_cols,
# log2_ratio_cols, channel_info, plex_prefix, has_peptides.
parse_tmt_file_structure <- function(all_lines, file_name = "unknown") {
  header_line <- all_lines[1]
  protein_headers <- as.character(
    read.csv(text = header_line, header = FALSE,
             stringsAsFactors = FALSE, check.names = FALSE)
  )
  protein_headers <- protein_headers[!is.na(protein_headers) &
                                       nchar(trimws(protein_headers)) > 0]

  # Detect protein rows via FDR Confidence column (quote-aware; see pd_protein_row_idx)
  protein_idx <- pd_protein_row_idx(all_lines, header_line)

  # Nested peptide data (lines starting with comma)
  pep_hdr_idx  <- grep("^,.*Annotated Sequence", all_lines)
  all_nested   <- grep("^,", all_lines)
  pep_data_idx <- setdiff(all_nested, pep_hdr_idx)

  # TMT channel abundance columns
  # Raw abundance columns ONLY. The "Abundances (Scaled)" fallback was removed
  # deliberately: those values are already normalised by Proteome Discoverer, and
  # log2 + median centring on top of them is not a valid quantification. Scaled
  # columns are still detected so the failure can be reported precisely.
  abund_cols  <- grep("^Abundance: F\\d+:", protein_headers, value = TRUE)
  scaled_cols <- grep("^Abundances? \\(Scaled\\): F\\d+:", protein_headers, value = TRUE)
  abundance_source <- if (length(abund_cols) > 0) "raw"
                      else if (length(scaled_cols) > 0) "scaled_only"
                      else "none"

  ratio_cols      <- grep("^Abundance Ratio:", protein_headers, value = TRUE)
  ratio_cols      <- ratio_cols[!grepl("Variability|P-Value", ratio_cols)]
  pval_cols       <- grep("Abundance Ratio Adj\\. P-Value:", protein_headers, value = TRUE)
  log2_ratio_cols <- grep("^Abundance Ratio \\(log2\\):", protein_headers, value = TRUE)

  plex_prefix <- unique(sub("^Abundances?.*?(F\\d+):.*", "\\1", abund_cols))

  channel_info <- data.frame(
    ColumnName = abund_cols,
    Tag        = sub("^.*F\\d+: ([^,]+),.*", "\\1", abund_cols),
    Condition  = sub("^.*F\\d+: [^,]+, [^,]+, ([^,]+),.*", "\\1", abund_cols),
    BioRep     = sub("^.*F\\d+: [^,]+, [^,]+, [^,]+, (\\d+).*", "\\1", abund_cols),
    stringsAsFactors = FALSE
  )

  list(
    name         = file_name,
    all_lines    = all_lines,
    header_line  = header_line,
    headers      = protein_headers,
    protein_idx  = protein_idx,
    pep_hdr_idx  = pep_hdr_idx,
    pep_data_idx = pep_data_idx,
    abund_cols   = abund_cols,
    abundance_source = abundance_source,
    n_scaled_cols    = length(scaled_cols),
    ratio_cols   = ratio_cols,
    pval_cols    = pval_cols,
    log2_ratio_cols = log2_ratio_cols,
    channel_info = channel_info,
    plex_prefix  = plex_prefix,
    has_peptides = length(pep_data_idx) > 10
  )
}

# Multi-plex PD input comes either as one file holding every plex (file IDs
# F1, F2, ... in the column names) or as one file per plex. Separate PD studies
# usually all write "F1", so two per-plex files can carry identical channel
# names ("Abundance: F1: 126, Sample, A, 1"); channels are matched to the
# annotation by name, so the second file's channels would land on the first
# file's. When names collide, every file ID of the later file is renumbered to
# an unused one ("F1" -> "F3") in its header lines and parsed fields, so all
# channels stay distinct. Files without collisions are returned unchanged.
tmt_disambiguate_plex_files <- function(files, logger = message) {
  if (length(files) < 2L) return(files)
  .fid_num <- function(x) suppressWarnings(as.integer(sub("^F", "", x)))
  used_ids <- unlist(lapply(files, function(f) .fid_num(f$plex_prefix)))
  next_id  <- max(c(0L, used_ids), na.rm = TRUE)
  seen     <- character(0)
  for (fi in seq_along(files)) {
    f <- files[[fi]]
    if (fi > 1L && any(f$abund_cols %in% seen)) {
      new_ids <- paste0("F", next_id + seq_along(f$plex_prefix))
      next_id <- next_id + length(f$plex_prefix)
      id_map  <- setNames(new_ids, f$plex_prefix)
      # Rewrites ": Fn:" tokens in one pass, so an ID already renamed is never
      # renamed again (F1 -> F2 while the file also has an F2).
      .remap <- function(x) {
        m <- gregexpr(": F\\d+:", x)
        regmatches(x, m) <- lapply(regmatches(x, m), function(v) {
          if (length(v) == 0L) return(v)
          id  <- sub("^: (F\\d+):$", "\\1", v)
          hit <- id %in% names(id_map)
          id[hit] <- unname(id_map[id[hit]])
          paste0(": ", id, ":")
        })
        x
      }
      hdr_lines <- unique(c(1L, f$pep_hdr_idx))
      f$all_lines[hdr_lines]   <- .remap(f$all_lines[hdr_lines])
      f$header_line            <- f$all_lines[1L]
      f$headers                <- .remap(f$headers)
      f$abund_cols             <- .remap(f$abund_cols)
      f$channel_info$ColumnName <- .remap(f$channel_info$ColumnName)
      f$plex_prefix            <- unname(id_map[f$plex_prefix])
      f$plex_id_renamed        <- id_map
      logger(paste0("TMT Upload: channel names in ", f$name, " repeat those of an earlier file; ",
                    "file ID(s) ", paste(names(id_map), collapse = ","), " renamed to ",
                    paste(id_map, collapse = ","), "."))
      files[[fi]] <- f
    }
    seen <- c(seen, f$abund_cols)
  }
  files
}

# Default plex for each channel in the design grid: the file it came from when
# several files were uploaded (one file per plex), or its PD file ID when one
# file holds several (F1 -> plex 1, F2 -> plex 2, in order of first appearance).
# Only a default: the user can change every channel's plex in the grid.
tmt_default_channel_plex <- function(files, channels) {
  out <- setNames(rep(1L, length(channels)), channels)
  if (length(files) > 1L) {
    for (fi in seq_along(files))
      out[intersect(files[[fi]]$abund_cols, channels)] <- fi
  } else {
    fid <- ifelse(grepl(": F\\d+:", channels), sub("^.*?: (F\\d+):.*$", "\\1", channels), NA)
    if (length(unique(stats::na.omit(fid))) > 1L)
      out[!is.na(fid)] <- match(fid[!is.na(fid)], unique(stats::na.omit(fid)))
  }
  out
}

# Parse peptide-level data from a file struct  -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -"
# Key: we strip only the leading comma from the raw CSV lines and let
# read.csv() handle quoted fields (column names with commas) correctly.
parse_pd_peptide_data <- function(f) {
  # Fast path: return cached data.frame if available (flat peptide build)
  if (!is.null(f$.peptide_df_cache)) return(f$.peptide_df_cache)

  all_lines    <- f$all_lines
  pep_hdr_idx  <- f$pep_hdr_idx
  pep_data_idx <- f$pep_data_idx
  if (length(pep_hdr_idx) == 0 || length(pep_data_idx) < 2) return(NULL)

  pep_header <- sub("^,+", "", all_lines[pep_hdr_idx[1]])
  pep_lines  <- sub("^,+", "", all_lines[pep_data_idx])
  pep_block  <- paste(c(pep_header, pep_lines), collapse = "\n")

  tryCatch(
    read.csv(text = pep_block, header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, na.strings = c("", "NA"), fill = TRUE),
    error = function(e) NULL
  )
}

# Parse protein-level data from a file struct  -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -" -"
# Returns list(protein_df, gene_names, accessions) or NULL.
# Uses the 3-strategy gene name resolution (user column → GN= → common cols).
parse_pd_protein_data <- function(f, rv = NULL) {
  # Fast path: return cached protein data if available (flat peptide build)
  if (!is.null(f$.protein_df_cache)) {
    protein_df <- f$.protein_df_cache
    accessions <- as.character(protein_df[["Accession"]])
    n_prot     <- length(accessions)
    gene_names <- rep("", n_prot)

    # Still run gene resolution strategies on cached data
    if (!is.null(rv) && !isTRUE(rv$pd_has_gn) && !is.null(rv$pd_gene_col) &&
        rv$pd_gene_col %in% colnames(protein_df)) {
      candidate <- as.character(protein_df[[rv$pd_gene_col]])
      candidate[is.na(candidate)] <- ""
      if (sum(candidate != "") / max(n_prot, 1) >= 0.3) gene_names <- candidate
    }
    if (sum(gene_names != "") / max(n_prot, 1) < 0.3) {
      descriptions <- protein_df[["Description"]]
      if (!is.null(descriptions)) {
        gn <- stringr::str_extract(descriptions, "GN=([^ ]+)")
        gn <- stringr::str_replace(gn, "^GN=", "")
        gn[is.na(gn)] <- ""
        if (sum(gn != "") / max(n_prot, 1) >= 0.3) gene_names <- gn
      }
    }
    if (sum(gene_names != "") / max(n_prot, 1) < 0.3) {
      for (gc in c("Gene Symbol", "Gene.Symbol", "Gene", "GeneName", "Gene.Name")) {
        gc_match <- grep(paste0("^", gc, "$"), colnames(protein_df), value = TRUE)[1]
        if (!is.na(gc_match)) {
          candidate <- as.character(protein_df[[gc_match]])
          candidate[is.na(candidate)] <- ""
          if (sum(candidate != "") / max(n_prot, 1) >= 0.3) { gene_names <- candidate; break }
        }
      }
    }
    empty_idx <- which(gene_names == "")
    if (length(empty_idx) > 0) gene_names[empty_idx] <- accessions[empty_idx]

    return(list(protein_df = protein_df, gene_names = gene_names, accessions = accessions))
  }

  if (length(f$protein_idx) == 0) return(NULL)

  header_line   <- f$all_lines[1]
  protein_lines <- f$all_lines[f$protein_idx]
  text_block    <- paste(c(header_line, protein_lines), collapse = "\n")

  protein_df <- tryCatch(
    read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, na.strings = c("", "NA"), fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(protein_df)) return(NULL)

  accessions <- as.character(protein_df[["Accession"]])
  n_prot     <- length(accessions)
  gene_names <- rep("", n_prot)
  tried_col  <- NULL

  # Strategy 1: User-selected gene column
  if (!is.null(rv) && !isTRUE(rv$pd_has_gn) && !is.null(rv$pd_gene_col) &&
      rv$pd_gene_col %in% colnames(protein_df)) {
    tried_col <- rv$pd_gene_col
    candidate <- as.character(protein_df[[rv$pd_gene_col]])
    candidate[is.na(candidate)] <- ""
    if (sum(candidate != "") / max(n_prot, 1) >= 0.3) {
      gene_names <- candidate
    }
  }

  # Strategy 2: GN= from Description
  if (sum(gene_names != "") / max(n_prot, 1) < 0.3) {
    descriptions <- protein_df[["Description"]]
    if (!is.null(descriptions)) {
      gn <- stringr::str_extract(descriptions, "GN=([^ ]+)")
      gn <- stringr::str_replace(gn, "^GN=", "")
      gn[is.na(gn)] <- ""
      if (sum(gn != "") / max(n_prot, 1) >= 0.3) {
        gene_names <- gn
      }
    }
  }

  # Strategy 3: Common gene-name columns
  if (sum(gene_names != "") / max(n_prot, 1) < 0.3) {
    gene_col_candidates <- c("Gene Symbol", "Gene.Symbol", "Gene",
                             "GeneName", "Gene.Name")
    for (gc in gene_col_candidates) {
      if (!is.null(tried_col) && gc == tried_col) next
      gc_match <- grep(paste0("^", gc, "$"), colnames(protein_df), value = TRUE)[1]
      if (!is.na(gc_match)) {
        candidate <- as.character(protein_df[[gc_match]])
        candidate[is.na(candidate)] <- ""
        if (sum(candidate != "") / max(n_prot, 1) >= 0.3) {
          gene_names <- candidate
          break
        }
      }
    }
  }

  empty_idx <- which(gene_names == "")
  if (length(empty_idx) > 0) {
    gene_names[empty_idx] <- accessions[empty_idx]
  }

  list(protein_df = protein_df, gene_names = gene_names, accessions = accessions)
}

# ===========================================================================
# Gene column auto-detection (shared by all modules)
# ===========================================================================
# Priority: 1) "Gene Symbol" / "Gene.Symbol" column, 2) GN= tags in
#           Description, 3) fall back to asking the user.
#
# Args:
#   headers       -" character vector of column names from the CSV
#   parsed_file   -" optional; a parsed file list (from parse_tmt_file_structure)
#                  used to sample Description rows for GN= check
#   sample_df     -" optional; a data.frame already parsed from the CSV (server.R
#                  builds this from raw lines). If provided AND parsed_file is
#                  NULL, GN= check uses this data.frame's Description column.
#   rv            -" reactive values list (optional; used for context)
#   logger        -" optional logging function
#
# Returns a list:
#   $source   -" "column" | "gn" | "ask"
#   $column   -" the auto-detected column name (for "column"), NULL otherwise
detect_gene_column <- function(headers, parsed_file = NULL, sample_df = NULL,
                               rv = NULL, logger = NULL) {
  log <- function(msg) { if (is.function(logger)) logger(msg) }

  # 1) Check for a Gene Symbol column in the headers
  gene_sym_col <- intersect(c("Gene Symbol", "Gene.Symbol"), headers)
  if (length(gene_sym_col) > 0) {
    log(paste0("Gene auto-detect: found column '", gene_sym_col[1], "'"))
    return(list(source = "column", column = gene_sym_col[1]))
  }

  # 2) Check for GN= tags in Description
  desc_col <- NULL
  if (!is.null(parsed_file)) {
    pdata <- parse_pd_protein_data(parsed_file, rv)
    if (!is.null(pdata) && "Description" %in% colnames(pdata$protein_df)) {
      desc_col <- pdata$protein_df[["Description"]]
    }
  } else if (!is.null(sample_df) && "Description" %in% colnames(sample_df)) {
    desc_col <- sample_df[["Description"]]
  }

  if (!is.null(desc_col)) {
    gn_hits <- sum(grepl("GN=", desc_col, ignore.case = FALSE), na.rm = TRUE)
    gn_pct  <- if (length(desc_col) > 0) gn_hits / length(desc_col) else 0
    if (!is.na(gn_pct) && gn_pct >= 0.5) {
      log(paste0("Gene auto-detect: GN= tags found in ",
                 round(gn_pct * 100), "% of descriptions"))
      return(list(source = "gn", column = NULL))
    } else {
      log(paste0("Gene auto-detect: GN= only in ",
                 round(gn_pct * 100), "%  -  insufficient"))
    }
  }

  # 3) Neither found  -" user must choose
  log("Gene auto-detect: no Gene Symbol column or GN= tags found")
  return(list(source = "ask", column = NULL))
}

# Search box + "Select All" for a checkboxGroupInput (wizard group pickers).
# The search keeps items whose visible label contains the term (case-
# insensitive); "Select All" checks only the items currently shown. Both run
# client-side (the .cbx-filter handlers in ui.R) so that checked items hidden
# by the filter keep their selection, which a server re-render would drop.
# `search_id`  -" inputId for the search textInput (reset it with updateTextInput)
# `target_id`  -" inputId of the checkboxGroupInput being filtered
checkbox_filter_controls <- function(search_id, target_id,
                                     placeholder = "Search samples (contains)...") {
  tags$div(class = "cbx-filter d-flex gap-1 align-items-start",
           `data-target` = target_id,
    tags$div(class = "flex-grow-1",
      textInput(search_id, NULL, width = "100%", placeholder = placeholder)),
    tags$button(type = "button",
                class = "btn btn-outline-primary btn-sm cbx-filter-select-all",
                icon("check-double"), " Select All")
  )
}

# Build a modal dialog asking the user to pick the Gene ID column.
# `radio_id`     -" inputId for the radioButtons (must be unique per module)
# `confirm_id`   -" inputId for the Confirm actionButton
# `headers`      -" character vector of available columns
build_gene_col_modal <- function(radio_id, confirm_id, headers) {
  # Best-guess default (Gene Symbol already ruled out by detect_gene_column)
  default <- if ("Gene" %in% headers) "Gene"
             else if ("GeneName" %in% headers) "GeneName"
             else if ("Gene.Name" %in% headers) "Gene.Name"
             else headers[1]
  modalDialog(
    title = tags$span(icon("triangle-exclamation"), " Gene ID Column Required"),
    tagList(
      tags$div(class = "alert alert-warning",
        icon("circle-info"),
        " No ", tags$b("Gene Symbol"), " column or ", tags$code("GN="),
        " tags were found. Please select the column containing Gene IDs."
      ),
      tags$div(
        style = "max-height: 350px; overflow-y: auto; border: 1px solid #eee; padding: 8px; border-radius: 4px;",
        radioButtons(radio_id, label = "Select Gene ID Column:",
                     choices  = setNames(headers, headers),
                     selected = default)
      )
    ),
    size = "m", easyClose = FALSE,
    footer = tagList(
      actionButton(confirm_id, "Confirm",
                   class = "btn-success", icon = icon("check"))
    )
  )
}

# Surface a warning to the user when possible, without requiring a reactive
# context. Inside a Shiny session this shows a notification; outside one (tests,
# scripts, batch use) getDefaultReactiveDomain() is NULL and the showNotification
# call is skipped. The log line and the base warning() always fire, so the
# message is never lost. Purely additive: no caller depends on its return value.
.osmss_user_warn <- function(msg, duration = 15) {
  try(warning(msg, call. = FALSE), silent = TRUE)
  try({
    if (exists("logger", mode = "function")) logger(paste0("[WARN] ", msg))
  }, silent = TRUE)
  try({
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(msg, type = "warning", duration = duration)
    }
  }, silent = TRUE)
  invisible(NULL)
}

# Whether plex offsets can be removed from a set of TMT channels without erasing
# the biology. limma::removeBatchEffect estimates each plex's offset from that
# plex's own channels, so a plex left holding a single channel (typical once the
# Norm/bridge channels are dropped for display) has its offset equal to that
# channel's whole deviation from its condition mean: the correction sets the
# channel to the condition mean. With one channel per plex every sample in a
# condition collapsed onto identical values in the heatmap and PCA. limma only
# warns ("Coefficients not estimable"), so the callers' tryCatch never fired.
# `batch_vec` is the plex of each channel to be corrected. TRUE only when there
# are >= 2 plexes and every plex holds >= 2 of those channels.
tmt_plex_correctable <- function(batch_vec) {
  n <- table(batch_vec[!is.na(batch_vec)])
  length(n) >= 2L && all(n >= 2L)
}

# Gene labels for TMT display modules, from the summarization cache's gene map
# (cache$gene_names, keyed by accession). This is the map DE uses; it is built
# at summarization time and restored with saved experiments, so Trend Clusters
# and Heatmap no longer re-parse the raw files with the live gene-column state.
# Returns an unnamed character vector parallel to `accessions`; accessions with
# no gene name are labelled with the accession itself.
tmt_cache_gene_labels <- function(cache, accessions) {
  accessions <- as.character(accessions)
  gl <- cache$gene_names
  # Unnamed maps aligned with the matrix rows (older caches) are keyed here.
  if (!is.null(gl) && is.null(names(gl)) && !is.null(cache$protein_matrix) &&
      length(gl) == nrow(cache$protein_matrix))
    names(gl) <- rownames(cache$protein_matrix)
  if (is.null(gl) || is.null(names(gl))) return(accessions)
  mapped <- unname(as.character(gl[accessions]))
  ifelse(!is.na(mapped) & nzchar(mapped), mapped, accessions)
}

# ===========================================================================
# EXPERIMENT RESET (see rv$reset_experiment_state in server.R)
# ===========================================================================
# rv keys a wizard upload leaves alone: the single-shot tools (Rank FC, Trends
# Analysis) and the PTM workflows have their own uploads and are not part of
# the global-proteome dataset being replaced. A saved-experiment import still
# resets everything.
RESET_KEEP_ON_UPLOAD <- c("^ss_pd_", "^dia_ptm_", "^tmt_ptm_")

# Register a module's own reactiveValues with the experiment reset. On every
# reset, each field is restored to the value it had when the module started
# (fields added later are set to NULL), except the fields named in `keep`.
# Call once, at module start-up, after the reactiveValues are created.
register_experiment_reset <- function(rv, store, keep = character(0)) {
  init <- isolate(reactiveValuesToList(store))
  isolate(rv$on_experiment_reset)(function() {
    for (k in isolate(names(store))) {
      if (k %in% keep) next
      store[[k]] <- if (k %in% names(init)) init[[k]] else NULL
    }
  })
}

# ===========================================================================
# SHARED DATA-FILTER HELPERS  -" used by DE and Rank FC modules
# ===========================================================================

# Build a standardised filter-editor modal dialog.
# prefix    -" unique string for input IDs (e.g. "de_f", "rfcf")
# headers   -" character vector of column names to offer for filtering
# existing_conditions  -" current list of condition lists (for group detection)
build_filter_modal <- function(prefix, headers, existing_conditions = list()) {
  existing_groups <- unique(vapply(
    existing_conditions, function(f) f$group, character(1)
  ))
  if (length(existing_groups) == 0) existing_groups <- character(0)
  next_group <- paste0("Group ", length(existing_groups) + 1)
  group_choices <- c(existing_groups, next_group)

  # Determine default within-group logic for the last group

  default_logic <- "OR"
  if (length(existing_groups) > 0) {
    last_grp <- tail(existing_groups, 1)
    grp_conds <- Filter(function(f) f$group == last_grp, existing_conditions)
    if (length(grp_conds) > 0 && !is.null(grp_conds[[1]]$logic))
      default_logic <- grp_conds[[1]]$logic
  }

  grp_id     <- paste0(prefix, "_group")
  logic_id   <- paste0(prefix, "_group_logic")
  col_id     <- paste0(prefix, "_column")
  op_id      <- paste0(prefix, "_operator")
  val_id     <- paste0(prefix, "_value")
  confirm_id <- paste0(prefix, "_confirm")

  modalDialog(
    title = tags$span(icon("filter"), " Data Filters"),
    tagList(
      helpText(class = "small text-muted mb-2",
        "Filter proteins by any column. ",
        "Conditions within a group use the selected logic (AND/OR); groups are AND'd."),
      fluidRow(
        column(8,
          tags$label(class = "form-label fw-bold", "Filter Group"),
          selectInput(grp_id, label = NULL,
                      choices  = group_choices,
                      selected = if (length(existing_groups) > 0) tail(existing_groups, 1)
                                 else next_group)
        ),
        column(4,
          tags$label(class = "form-label fw-bold", "Within-Group Logic"),
          radioButtons(logic_id, label = NULL,
                       choices = c("OR" = "OR", "AND" = "AND"),
                       selected = default_logic, inline = TRUE)
        )
      ),
      tags$div(class = "mb-3",
        tags$label(class = "form-label fw-bold", "Column"),
        tags$div(style = "max-height: 250px; overflow-y: auto; border: 1px solid #eee; padding: 6px; border-radius: 4px;",
          radioButtons(col_id, label = NULL,
                       choices  = setNames(headers, headers),
                       selected = headers[1])
        )
      ),
      fluidRow(
        column(6, selectInput(op_id, "Operator",
                  choices = c(">" = ">", ">=" = ">=", "<" = "<", "<=" = "<=",
                              "==" = "==", "!=" = "!=",
                              "contains" = "contains",
                              "not contains" = "not_contains",
                              "is NA (empty)" = "isNA",
                              "is not NA" = "notNA"))),
        column(6, textInput(val_id, "Value",
                            placeholder = "e.g. 5 (leave blank for isNA/notNA)"))
      ),
      # JS: disable value input when isNA/notNA selected
      tags$script(HTML(sprintf("
        $(document).on('change', '#%s', function() {
          var op = $(this).val();
          if (op === 'isNA' || op === 'notNA') {
            $('#%s').val('').prop('disabled', true)
              .attr('placeholder', 'Not needed for this operator');
          } else {
            $('#%s').prop('disabled', false)
              .attr('placeholder', 'e.g. 5');
          }
        });
      ", op_id, val_id, val_id)))
    ),
    size = "m", easyClose = TRUE,
    footer = tagList(
      modalButton("Cancel"),
      actionButton(confirm_id, "Add Condition", class = "btn-success",
                   icon = icon("check"))
    )
  )
}

# Render grouped condition badges with AND/OR indicators and remove buttons.
# conds   -" list of condition lists (each with id, group, column, op, value, logic)
# prefix  -" string prepended to remove-button IDs (e.g. "de_frm_", "rfcf_rm_")
render_filter_badges <- function(conds, prefix) {
  groups <- unique(vapply(conds, function(f) f$group, character(1)))
  group_tags <- lapply(groups, function(grp) {
    grp_conds <- Filter(function(f) f$group == grp, conds)
    # Determine within-group logic from stored conditions
    grp_logic <- "OR"
    for (f in grp_conds) {
      if (!is.null(f$logic) && f$logic %in% c("OR", "AND")) {
        grp_logic <- f$logic; break
      }
    }
    cond_tags <- lapply(grp_conds, function(f) {
      op_label <- switch(f$op,
        "contains"     = "contains",
        "not_contains" = "not contains",
        "isNA"         = "is NA",
        "notNA"        = "is not NA",
        f$op)
      badge_text <- if (f$op %in% c("isNA", "notNA"))
        paste0(f$column, " ", op_label)
      else
        paste0(f$column, " ", op_label, " ", f$value)
      tags$div(class = "d-flex align-items-center gap-1 mb-1",
        tags$span(class = "badge bg-light text-dark",
                  style = "font-size: 0.8em; word-break: break-all;",
                  badge_text),
        actionButton(paste0(prefix, f$id), label = NULL,
                     icon = icon("xmark"),
                     class = "btn-outline-danger btn-sm",
                     style = "padding: 1px 5px; font-size: 0.7em;")
      )
    })
    separated <- list()
    logic_badge_class <- if (grp_logic == "AND") "badge bg-primary" else "badge bg-warning text-dark"
    for (j in seq_along(cond_tags)) {
      separated <- c(separated, list(cond_tags[[j]]))
      if (j < length(cond_tags))
        separated <- c(separated, list(
          tags$span(class = logic_badge_class,
                    style = "font-size: 0.75em;", grp_logic)
        ))
    }
    tags$div(class = "mb-2",
      style = "background: #f8f9fa; border-left: 3px solid #18bc9c; padding: 6px 8px; border-radius: 4px;",
      tags$strong(style = "font-size: 0.8em; color: #6c757d;", grp),
      tagList(separated)
    )
  })
  final_tags <- list()
  for (k in seq_along(group_tags)) {
    final_tags <- c(final_tags, list(group_tags[[k]]))
    if (k < length(group_tags))
      final_tags <- c(final_tags, list(
        tags$div(class = "text-center mb-1",
          tags$span(class = "badge bg-primary", style = "font-size: 0.8em;", "AND"))
      ))
  }
  tagList(final_tags)
}

# Evaluate filter/highlight conditions against a data.frame.
# Returns a logical vector (TRUE = keep / matches).
# logic = "and_or"  → group-aware (AND across groups, per-group AND/OR)
# logic = "or"      → simple OR across all conditions (for highlights)
apply_filter_conditions <- function(df, conditions, logic = "and_or") {
  if (length(conditions) == 0) return(rep(TRUE, nrow(df)))

  eval_mask <- function(cond, df) {
    col <- cond$column
    op  <- cond$op
    val <- cond$value
    if (!(col %in% colnames(df))) return(rep(FALSE, nrow(df)))
    col_data <- df[[col]]

    # isNA / notNA operators
    if (op == "isNA")  return(is.na(col_data) | as.character(col_data) == "")
    if (op == "notNA") return(!is.na(col_data) & as.character(col_data) != "")

    if (op %in% c("contains", "not_contains")) {
      col_text <- as.character(col_data)
      col_text[is.na(col_text)] <- ""
      m <- grepl(val, col_text, ignore.case = TRUE)
      if (op == "not_contains") m <- !m
      return(m)
    }

    num_data <- suppressWarnings(as.numeric(col_data))
    num_val  <- suppressWarnings(as.numeric(val))
    if (is.na(num_val)) {
      col_text <- as.character(col_data)
      col_text[is.na(col_text)] <- ""
      return(switch(op,
        "=="  = col_text == val,
        "!="  = col_text != val,
        rep(FALSE, length(col_text))
      ))
    }
    switch(op,
      ">"  = !is.na(num_data) & num_data >  num_val,
      ">=" = !is.na(num_data) & num_data >= num_val,
      "<"  = !is.na(num_data) & num_data <  num_val,
      "<=" = !is.na(num_data) & num_data <= num_val,
      "==" = !is.na(num_data) & num_data == num_val,
      "!=" = !is.na(num_data) & num_data != num_val,
      rep(FALSE, length(num_data))
    )
  }

  if (logic == "and_or") {
    groups <- unique(vapply(conditions, function(f) f$group, character(1)))
    group_masks <- lapply(groups, function(grp) {
      grp_conds <- Filter(function(f) f$group == grp, conditions)
      grp_logic <- "OR"
      for (f in grp_conds) {
        if (!is.null(f$logic) && f$logic %in% c("OR", "AND")) {
          grp_logic <- f$logic; break
        }
      }
      masks <- lapply(grp_conds, function(cond) eval_mask(cond, df))
      if (grp_logic == "AND") Reduce(`&`, masks) else Reduce(`|`, masks)
    })
    return(Reduce(`&`, group_masks))
  } else {
    or_mask <- rep(FALSE, nrow(df))
    for (cond in conditions) or_mask <- or_mask | eval_mask(cond, df)
    return(or_mask)
  }
}

# ===========================================================================
# SHARED TMT HELPERS  -" used by Volcano, Trends, Heatmap, PCA modules
# ===========================================================================

# Short display label for TMT abundance column names
tmt_short_label <- function(col_name) {
  lbl <- sub("^Abundances?[^:]*:\\s*", "", col_name)
  lbl
}

# Extract TMT tag from column name (e.g. "126", "127N")
tmt_extract_tag <- function(col_name) {
  m <- regmatches(col_name, regexec("F\\d+:\\s*([^,]+)", col_name))[[1]]
  if (length(m) >= 2) trimws(m[2]) else col_name
}

# =============================================================================
# PTM MODIFICATION-SITE PARSING
# =============================================================================

# Parses a "Modifications in Master Proteins" (MMP) string into per-site rows
# with accession, site, modification and localization probability.
parse_ptm_mmp_sites <- function(mmp_txt, tgt_pattern, absent_prob = 0,
                                prob_scale = NULL) {
  if (is.null(mmp_txt) || is.na(mmp_txt) || trimws(mmp_txt) == "")
    return(data.frame(Accession = character(0), Site = character(0),
                      Modification = character(0),
                      Localization_Prob = numeric(0),
                      stringsAsFactors = FALSE))
  txt    <- trimws(as.character(mmp_txt))
  chunks <- unlist(strsplit(txt, ";\\s*(?=[A-Z0-9][A-Z0-9._-]*\\s+\\d+x)",
                            perl = TRUE))
  chunks <- trimws(chunks[nchar(trimws(chunks)) > 0])
  if (length(chunks) == 0) chunks <- txt

  acc_v <- character(0); site_v <- character(0)
  mod_v <- character(0); prob_v <- numeric(0)

  for (ch in chunks) {
    m <- regmatches(ch, regexec(
      "^([A-Z0-9][A-Z0-9._-]*)\\s+\\d+x([^\\[]+)\\s*\\[([^\\]]+)\\]",
      ch, perl = TRUE))[[1]]
    if (length(m) < 4) next
    acc        <- trimws(m[2])
    mod        <- trimws(m[3])
    site_block <- trimws(m[4])
    if (!grepl(tgt_pattern, mod, ignore.case = TRUE)) next
    for (se in strsplit(site_block, ";\\s*")[[1]]) {
      se <- trimws(se)
      sm <- regmatches(se, regexec(
        "([A-Za-z])(\\d+)(?:\\(([0-9.]+)\\))?", se, perl = TRUE))[[1]]
      if (length(sm) >= 3 && nchar(sm[2]) > 0 && nchar(sm[3]) > 0) {
        site_label <- paste0(sm[2], sm[3])
        prob_val   <- if (length(sm) >= 4 && nchar(sm[4]) > 0)
          suppressWarnings(as.numeric(sm[4])) else absent_prob
        if (is.na(prob_val) && !is.na(absent_prob)) prob_val <- absent_prob
      } else {
        site_label <- "Ambiguous"; prob_val <- absent_prob
      }
      acc_v  <- c(acc_v,  acc)
      site_v <- c(site_v, site_label)
      mod_v  <- c(mod_v,  mod)
      prob_v <- c(prob_v, prob_val)
    }
  }
  if (is.null(prob_scale)) {
    nonzero <- prob_v[!is.na(prob_v) & prob_v > 0]
    if (length(nonzero) > 0 && all(nonzero <= 1)) prob_v <- prob_v * 100
  } else if (identical(prob_scale, "fraction")) {
    prob_v <- prob_v * 100
  }
  data.frame(Accession = acc_v, Site = site_v, Modification = mod_v,
             Localization_Prob = prob_v, stringsAsFactors = FALSE)
}

# =============================================================================
# TMT FLAT PEPTIDE HELPERS (proteome  -" no PTM fields)
# =============================================================================

# Show a column-mapping modal for TMT flat peptide input.
# `col_names`   -" character vector of column names from the uploaded file.
# `ns_prefix`   -" module-specific prefix (e.g. "tmt_flat_").
# `file_name`   -" original filename for display.
# `file_idx`    -" 1-based index of the current file.
# `n_files`     -" total number of files being uploaded.
flat_tmt_show_mapping_modal <- function(session, col_names, ns_prefix,
                                        file_name = "File 1",
                                        file_idx  = 1L,
                                        n_files   = 1L) {
  col_choices <- setNames(col_names, col_names)
  none_choice <- c("Unknown/None" = "")
  contained_choice <- c("Use Peptide Sequence" = "contained")

  # Auto-detect likely columns
  auto_acc    <- grep("(Master Protein Accession|Protein.Accession|ProteinAccession|Accession)",
                      col_names, value = TRUE, ignore.case = TRUE)[1]
  auto_seq    <- grep("(Annotated Sequence|Sequence|Peptide)", col_names, value = TRUE,
                      ignore.case = TRUE)[1]
  auto_gene   <- grep("(Gene Symbol|Gene\\.Symbol|Gene$|GeneName|GeneID)",
                      col_names, value = TRUE, ignore.case = TRUE)[1]
  auto_charge <- grep("(Charge|Precursor Charge)", col_names, value = TRUE,
                      ignore.case = TRUE)[1]
  auto_abund  <- grep("(?i)(abundance|intensity|signal)", col_names, value = TRUE, perl = TRUE)

  id_suffix <- paste0("f", file_idx, "_")
  mid       <- paste0(ns_prefix, id_suffix)
  title_text <- paste0("Map Columns \u2014 File ", file_idx, " of ", n_files, ": ", file_name)

  showModal(modalDialog(
    title = tags$span(icon("table-columns"), paste0(" ", title_text)),
    size = "l", easyClose = FALSE,
    tagList(
      tags$div(class = "alert alert-info py-2 px-3 small mb-3",
        icon("circle-info"),
        " Map columns to the required fields. Each file is mapped individually.",
        if (n_files > 1 && file_idx == 1)
          tags$span(class = "ms-2 text-primary",
            icon("copy"), " An 'Apply to all' option is available below.")
        else NULL
      ),

      # ── Core Identifiers ──
      tags$h6(icon("key"), " Core Identifiers"),
      selectInput(paste0(mid, "map_accession"),
                  "Protein Accession Column:",
                  choices = col_choices,
                  selected = if (!is.na(auto_acc)) auto_acc else col_names[1]),
      selectInput(paste0(mid, "map_sequence"),
                  "Peptide Sequence Column:",
                  choices = col_choices,
                  selected = if (!is.na(auto_seq)) auto_seq else col_names[1]),
      selectInput(paste0(mid, "map_gene"),
                  "Gene Symbol / Description Column:",
                  choices = c(none_choice, col_choices),
                  selected = if (!is.na(auto_gene)) auto_gene else ""),
      selectInput(paste0(mid, "map_charge"),
                  "Precursor Charge Column:",
                  choices = c(contained_choice, col_choices),
                  selected = if (!is.na(auto_charge)) auto_charge else "contained"),

      tags$hr(),

      # ── Quantitative Data ──
      tags$h6(icon("chart-bar"), " Quantitative Data"),
      helpText(class = "small text-muted",
        "Select ONLY raw, un-normalized abundance/intensity columns. ",
        "Each selected column becomes one TMT channel."),
      tags$div(
        style = "max-height: 250px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 4px; padding: 8px;",
        checkboxGroupInput(paste0(mid, "map_abundance"),
                           label = NULL,
                           choices = col_choices,
                           selected = auto_abund)
      ),

      # ── Apply-all checkbox (only file 1 in multi-file upload) ──
      if (n_files > 1 && file_idx == 1) {
        tagList(
          tags$hr(),
          checkboxInput(paste0(ns_prefix, "apply_all"),
                        "Apply this mapping to all remaining files",
                        value = FALSE)
        )
      } else NULL
    ),
    footer = tagList(
      # An observed button, not modalButton(): the caller decides what Cancel
      # returns to. The TMT wizard reopens itself; a bare modalButton closed
      # this dialog and left no wizard on screen.
      actionButton(paste0(ns_prefix, "flat_cancel"), "Cancel", class = "btn-secondary"),
      actionButton(paste0(ns_prefix, "flat_confirm"),
                   paste0("Confirm File ", file_idx),
                   class = "btn-success", icon = icon("check"))
    )
  ))
}

# Collect mapping values set in the flat_tmt_show_mapping_modal.
# Returns list(accession_col, sequence_col, gene_col, charge_col, abundance_cols).
flat_tmt_collect_mapping <- function(input, ns_prefix, file_idx = 1L) {
  id_suffix <- paste0("f", file_idx, "_")
  mid       <- paste0(ns_prefix, id_suffix)
  list(
    accession_col  = input[[paste0(mid, "map_accession")]],
    sequence_col   = input[[paste0(mid, "map_sequence")]],
    gene_col       = input[[paste0(mid, "map_gene")]],
    charge_col     = input[[paste0(mid, "map_charge")]],
    abundance_cols = input[[paste0(mid, "map_abundance")]]
  )
}

# Build a mock PD-nested CSV character vector from a flat TMT peptide table.
# `pep_df`     -" data.frame of raw peptide-level data (all columns character).
# `mapping`    -" list from flat_tmt_collect_mapping.
# `plex_idx`   -" 1-based plex index; controls F-number in abundance column names.
# Returns a character vector of CSV lines parseable by parse_tmt_file_structure.
flat_tmt_build_mock_pd <- function(pep_df, mapping, plex_idx = 1L) {
  work <- pep_df

  # ── Rename mapped columns to canonical names ──
  rename_col <- function(df, from, to) {
    if (!is.null(from) && nchar(from) > 0 && from %in% colnames(df) && from != to) {
      colnames(df)[colnames(df) == from] <- to
    }
    df
  }
  work <- rename_col(work, mapping$accession_col, "Master Protein Accessions")
  work <- rename_col(work, mapping$sequence_col,  "Annotated Sequence")
  if (!is.null(mapping$gene_col) && nchar(mapping$gene_col) > 0)
    work <- rename_col(work, mapping$gene_col, "Gene Symbol")

  # Handle charge: use mapped column or default to "2"
  if (!is.null(mapping$charge_col) && mapping$charge_col != "contained" &&
      nchar(mapping$charge_col) > 0 && mapping$charge_col %in% colnames(work)) {
    colnames(work)[colnames(work) == mapping$charge_col] <- "Charge"
  } else {
    if (!"Charge" %in% colnames(work)) work$Charge <- "2"
  }

  # ── Resolve user-selected abundance columns ──
  raw_abund <- mapping$abundance_cols
  raw_abund <- raw_abund[raw_abund %in% colnames(work)]
  if (length(raw_abund) == 0) stop("No abundance columns found in flat TMT file.")

  # Rename abundance columns to PD-style: Abundance: F<plex>: <label>, Sample, Cond, 1
  # The label is the original column name with commas replaced, so tmt_extract_tag()
  # returns the readable original label in the design grid.
  clean_labels <- gsub(",", " ", raw_abund)
  pd_style_names <- paste0("Abundance: F", plex_idx, ": ", clean_labels, ", Sample, Cond, 1")

  # Guard against accidental duplicates
  if (anyDuplicated(pd_style_names)) {
    for (j in which(duplicated(pd_style_names)))
      pd_style_names[j] <- paste0(pd_style_names[j], "_dup", j)
  }

  # Rename in work frame
  colnames(work)[match(raw_abund, colnames(work))] <- pd_style_names
  abund_cols <- pd_style_names

  # Ensure abundance columns are numeric
  for (ac in abund_cols) work[[ac]] <- suppressWarnings(as.numeric(work[[ac]]))

  # Defaults and validation
  if (!"Charge" %in% colnames(work)) work$Charge <- "2"
  acc_col <- "Master Protein Accessions"
  if (!acc_col %in% colnames(work)) stop("Accession column not found after renaming.")
  work[[acc_col]][is.na(work[[acc_col]]) | nchar(trimws(work[[acc_col]])) == 0] <- "Unknown"

  # Drop multi-protein peptides
  multi <- grepl(";", work[[acc_col]], fixed = TRUE)
  if (any(multi)) work <- work[!multi, , drop = FALSE]
  if (nrow(work) == 0) stop("No peptides remaining after removing multi-protein entries.")

  # ── Gene lookup ──
  gene_vals <- rep("", nrow(work))
  if ("Gene Symbol" %in% colnames(work)) {
    gene_vals <- as.character(work[["Gene Symbol"]])
    gene_vals[is.na(gene_vals)] <- ""
  }
  unique_accessions <- unique(work[[acc_col]])
  gene_lookup <- data.frame(
    Accession = as.character(work[[acc_col]]),
    Gene      = gene_vals,
    stringsAsFactors = FALSE
  )
  gene_lookup <- gene_lookup[!duplicated(gene_lookup$Accession), , drop = FALSE]
  gene_map <- setNames(gene_lookup$Gene, gene_lookup$Accession)

  # ── Summarize: sum peptide abundances per protein ──
  prot_rows <- list()
  for (acc in unique_accessions) {
    mask <- work[[acc_col]] == acc
    sub  <- work[mask, , drop = FALSE]
    sums <- vapply(abund_cols, function(ac) sum(sub[[ac]], na.rm = TRUE), numeric(1))
    gene <- gene_map[acc]
    if (is.na(gene) || gene == "") gene <- acc
    desc <- paste0("GN=", gene)
    prot_rows[[length(prot_rows) + 1L]] <- c(acc, "High", gene, desc, as.character(sums))
  }

  # ── Build PD-nested CSV text ──
  .quote_if_comma <- function(x) {
    vapply(x, function(v) if (grepl(",", v, fixed = TRUE)) paste0('"', v, '"') else v,
           character(1))
  }

  prot_header_cols <- c("Accession", "FDR Confidence", "Gene Symbol", "Description", abund_cols)
  prot_header_line <- paste(.quote_if_comma(prot_header_cols), collapse = ",")

  # Ensure required peptide columns exist
  if (!"Confidence" %in% colnames(work)) work$Confidence <- "High"

  pep_keep_cols <- c("Annotated Sequence", "Master Protein Accessions",
                     "Confidence", "Charge", abund_cols)
  pep_keep_cols <- pep_keep_cols[pep_keep_cols %in% colnames(work)]

  pep_header_line <- paste0(",", paste(.quote_if_comma(pep_keep_cols), collapse = ","))

  all_lines <- prot_header_line
  for (acc in unique_accessions) {
    idx      <- match(acc, unique_accessions)
    prot_row <- prot_rows[[idx]]
    all_lines <- c(all_lines, paste(.quote_if_comma(prot_row), collapse = ","))
    all_lines <- c(all_lines, pep_header_line)
    mask    <- work[[acc_col]] == acc
    sub_pep <- work[mask, pep_keep_cols, drop = FALSE]
    for (r in seq_len(nrow(sub_pep))) {
      vals <- vapply(pep_keep_cols, function(cn) {
        v <- as.character(sub_pep[[cn]][r])
        if (is.na(v)) "" else if (grepl(",", v, fixed = TRUE)) paste0('"', v, '"') else v
      }, character(1))
      all_lines <- c(all_lines, paste0(",", paste(vals, collapse = ",")))
    }
  }

  all_lines
}

# === MODULE 0: SHARED HELPERS ================================================
# These functions are used by BOTH Residue Prediction and Local MS modules.
# They are defined once here to avoid duplication.

sanitize_uniprot_id <- function(uniprot_id) {
  raw_id <- as.character(uniprot_id)
  cleaned <- gsub("[^A-Za-z0-9_-]", "", raw_id)
  return(cleaned)
}

ensure_trailing_newline <- function(file_path) {
  if (!file.exists(file_path)) return(FALSE)
  file_info <- file.info(file_path)
  if (is.na(file_info$size) || file_info$size == 0) return(FALSE)
  con <- file(file_path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = -1, origin = "end")
  last_char <- readChar(con, 1, useBytes = TRUE)
  if (last_char != "\n") {
    cat("\n", file = file_path, append = TRUE)
  }
  return(TRUE)
}

fetch_uniprot_seq <- function(uniprot_id) {
  safe_id <- sanitize_uniprot_id(uniprot_id)
  if (nchar(safe_id) == 0) return(NULL)
  url <- paste0("https://rest.uniprot.org/uniprotkb/", safe_id, ".json")
  tryCatch({
    resp <- jsonlite::fromJSON(httr::content(httr::GET(url), "text", encoding="UTF-8"))
    return(resp$sequence$value)
  }, error = function(e) return(NULL))
}

# Convenience wrapper  -" used by Local MS module
get_uniprot_length <- function(uniprot_id) {
  seq_str <- fetch_uniprot_seq(uniprot_id)
  if (is.null(seq_str)) return(NA_integer_)
  return(nchar(seq_str))
}

get_seq_context <- function(sequence, position, window=4) {
  if(position < 1 || position > nchar(sequence)) return("Invalid Index")
  # Pad with X to avoid truncated context at termini
  padded <- paste0(paste(rep("X", window), collapse = ""),
                   sequence,
                   paste(rep("X", window), collapse = ""))
  pos_pad <- position + window
  prefix  <- substr(padded, pos_pad - window, pos_pad - 1)
  residue <- substr(padded, pos_pad, pos_pad)
  suffix  <- substr(padded, pos_pad + 1, pos_pad + window)
  return(paste0(prefix, "(", residue, ")", suffix))
}

load_global_resources <- function(path_st, path_y, path_fav, dir_bg) {
  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] DEBUG: Starting load_global_resources..."))
  if(!file.exists(path_st) || !file.exists(path_y)) {
    message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] DEBUG: !!! Files not found."))
    return(NULL)
  }
  st_raw <- suppressMessages(read_tsv(path_st, show_col_types = FALSE))
  st_mat <- as.matrix(st_raw[, -1]); rownames(st_mat) <- st_raw[[1]]
  y_raw <- suppressMessages(read_tsv(path_y, show_col_types = FALSE))
  y_mat <- as.matrix(y_raw[, -1]); rownames(y_mat) <- y_raw[[1]]
  bias_vec_s <- setNames(rep(0, nrow(st_mat)), rownames(st_mat))
  bias_vec_t <- setNames(rep(0, nrow(st_mat)), rownames(st_mat))
  if (file.exists(path_fav)) {
    fav_df <- suppressMessages(read_tsv(path_fav, show_col_types = FALSE))
    common_kinases <- intersect(fav_df[[1]], rownames(st_mat))
    if("s" %in% tolower(colnames(fav_df))) {
      s_col <- grep("s", colnames(fav_df), ignore.case=TRUE, value=TRUE)[1]
      for(k in common_kinases) { val <- as.numeric(fav_df[fav_df[[1]] == k, s_col]); if(!is.na(val) && val > 0) bias_vec_s[k] <- log2(val) }
    }
    if("t" %in% tolower(colnames(fav_df))) {
      t_col <- grep("t", colnames(fav_df), ignore.case=TRUE, value=TRUE)[1]
      for(k in common_kinases) { val <- as.numeric(fav_df[fav_df[[1]] == k, t_col]); if(!is.na(val) && val > 0) bias_vec_t[k] <- log2(val) }
    }
  }
  return(list(ST_MAT = st_mat, Y_MAT = y_mat, BIAS_S = bias_vec_s, BIAS_T = bias_vec_t, BG_DIR = dir_bg, LOADED=TRUE))
}

fetch_tkl_kinase_list <- function() {
  ensure_global_resources()
  if(!GLOBAL_RESOURCES$LOADED) return(NULL)
  st_kins <- if(!is.null(GLOBAL_RESOURCES$ST_MAT)) rownames(GLOBAL_RESOURCES$ST_MAT) else c()
  y_kins  <- if(!is.null(GLOBAL_RESOURCES$Y_MAT)) rownames(GLOBAL_RESOURCES$Y_MAT) else c()
  return(sort(unique(c(st_kins, y_kins))))
}

# Kinase group table (Johnson 2023 + Yaron-Barir 2024 Supplementary Table 1, UniProt pkinfam; CC BY 4.0)
load_kinase_families <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)
  df <- read.delim(file_path, stringsAsFactors = FALSE, na.strings = "", quote = "")
  df$Kinase <- toupper(df$Kinase)
  df[!duplicated(df$Kinase), ]
}

# Named vector of upper-case group labels (AGC, CAMK, ..., TK); unmatched kinases -> "OTHER"
assign_kinase_family <- function(kinase_vec, fam_df) {
  if (is.null(fam_df)) return(rep("Unknown", length(kinase_vec)))
  g <- toupper(fam_df$Group[match(toupper(kinase_vec), fam_df$Kinase)])
  g[is.na(g) | g == ""] <- "OTHER"
  stats::setNames(g, kinase_vec)
}

get_3d_sort_order <- function(pdb_path, residue_positions) {
  if (is.null(pdb_path) || !file.exists(pdb_path)) return(NULL)
  require_pkg("bio3d")
  pdb <- suppressWarnings(bio3d::read.pdb(pdb_path))
  sele <- bio3d::atom.select(pdb, resno = residue_positions, elety = "CA")
  if (length(sele$atom) < 3) return(NULL)
  xyz <- pdb$xyz[sele$xyz]; coords_mat <- matrix(xyz, ncol = 3, byrow = TRUE); rownames(coords_mat) <- residue_positions
  d_mat <- dist(coords_mat); hc <- hclust(d_mat)
  return(as.numeric(hc$labels[hc$order]))
}

download_pdb_r_global <- function(uniprot_id) {
  safe_id <- sanitize_uniprot_id(uniprot_id)
  if (nchar(safe_id) == 0) return(NULL)
  dest <- file.path(tempdir(), paste0(safe_id, ".pdb"))
  if (file.exists(dest) && file.info(dest)$size > 1000) {
    ensure_trailing_newline(dest)
    return(dest)
  }
  versions <- c('v6', 'v4', 'v3', 'v2', 'v1')
  for (v in versions) {
    url <- sprintf('https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_%s.pdb', safe_id, v)
    tryCatch({
      download.file(url, dest, quiet = TRUE, mode = "wb")
      if (file.exists(dest) && file.info(dest)$size > 1000) {
        ensure_trailing_newline(dest)
        return(dest)
      }
    }, error = function(e) {})
  }
  return(NULL)
}

# === LOCAL MS: FRAGPIPE HELPERS ==============================================
# Modification parsing, classification, and confidence scoring for FragPipe data.
# These are specific to the Local MS module but defined globally for consistency.

# Normalize modification names for consistency
normalize_mod_name <- function(mod_name) {
  # Collapse all "Unannotated mass-shift <number>" variants into one
  if (grepl("^Unannotated mass-shift", mod_name)) return("Unannotated mass-shift")
  # Rename bare "Hex" (hexose sugar +162 Da) to explicit name
  if (mod_name == "Hex") return("Hexose Sugar")
  return(mod_name)
}

parse_modifications <- function(mod_string) {
  if (is.na(mod_string) || mod_string == "") return(list())
  mods <- list()
  matches <- str_extract_all(mod_string, "Mod\\d+:\\s*([^,(]+)")[[1]]
  if (length(matches) > 0) {
    for (match in matches) {
      mod_name <- str_extract(match, "Mod\\d+:\\s*(.+)$")
      mod_name <- str_remove(mod_name, "^Mod\\d+:\\s*")
      mod_name <- str_trim(mod_name)
      if (mod_name != "") mods[[length(mods) + 1]] <- normalize_mod_name(mod_name)
    }
  }
  return(unique(mods))
}

parse_mass_error_for_mod <- function(obs_mod_string, mod_name) {
  if (is.na(obs_mod_string) || obs_mod_string == "") return(NA_real_)
  escaped_mod <- str_replace_all(mod_name, "([\\[\\](){}.*+?^$|\\\\])", "\\\\\\1")
  pattern <- paste0(escaped_mod, "\\s*\\(PeakApex:\\s*([0-9.eE+-]+),\\s*Theoretical:\\s*([0-9.eE+-]+)\\)")
  match <- str_match(obs_mod_string, pattern)
  if (!is.na(match[1, 1])) {
    peak_apex <- as.numeric(match[1, 2])
    theoretical <- as.numeric(match[1, 3])
    if (!is.na(peak_apex) && !is.na(theoretical)) return(abs(peak_apex - theoretical))
  }
  return(NA_real_)
}

extract_residue_positions <- function(localization_string) {
  if (is.na(localization_string) || localization_string == "") return(integer(0))
  chars <- strsplit(localization_string, "")[[1]]
  positions <- which(grepl("[a-z]", chars))
  return(positions)
}

split_mapped_genes <- function(mapped_genes_string) {
  if (is.na(mapped_genes_string) || mapped_genes_string == "") return(NULL)
  genes <- str_split(mapped_genes_string, ",")[[1]] %>% str_trim() %>% unique()
  return(genes)
}

expand_position_string <- function(position_string) {
  if (is.na(position_string) || position_string == "") return(integer(0))
  parts <- str_split(position_string, ",")[[1]]
  positions <- suppressWarnings(as.integer(str_trim(parts)))
  positions <- positions[!is.na(positions)]
  return(positions)
}

# Modification classification lists
KNOWN_BIOLOGICAL <- c(
  "Phospho", "Phosphorylation", "Acetyl", "Acetylation",
  "Ubiquitinylation", "Ubiquitin", "SUMOylation", "SUMO",
  "Citrullination", "ADP-Ribosylation",
  "Glycosylation", "HexNAc", "Hexose Sugar",
  "Nitrosylation", "Hydroxylation",
  "Palmitoylation", "Myristoylation",
  "Crotonylation", "Succinylation", "Propionylation",
  "Butyrylation", "Malonylation", "Glutarylation",
  "TriMethylation", "Quinone"
)

KNOWN_ARTIFACTS <- c(
  "First isotopic peak", "Isotopic peak error",
  "Second isotopic peak", "Third isotopic peak",
  "Carbamidomethyl/Addition of G",
  "Carbamylation", "Carbamyl", "Formylation",
  "Sodium adduct", "Sodium", "Potassium adduct",
  "Carbon Adduct", "Sulfonation",
  "Aminoethylbenzenesulfonylation", "Persulfide",
  "Half of a disulfide bridge", "Propionamide",
  "Met loss", "Acetald",
  "Replacement of proton", "Replacement of 2 protons",
  "Zinc", "Addition of"
)

CONTEXT_DEPENDENT <- c(
  "Deamidation", "Oxidation", "Dioxidation", "Trioxidation",
  "Methyl", "DiMethyl/Ethyl", "DiMethyl", "Carboxymethylation",
  "Ammonia loss", "Didehydrobutyrine/Water loss",
  "Amidation", "Water loss",
  "Pyrrolidinone", "Dehydroalanine", "Didehydro",
  "H-2 O1", "Kynurenine",
  "Arginine not cleaved", "Lysine not cleaved"
)

classify_mod_type <- function(mod_name) {
  mod_lower <- tolower(mod_name)
  for (bio in KNOWN_BIOLOGICAL) {
    if (grepl(tolower(bio), mod_lower, fixed = TRUE)) return("biological")
  }
  for (art in KNOWN_ARTIFACTS) {
    if (grepl(tolower(art), mod_lower, fixed = TRUE)) return("artifact")
  }
  for (ctx in CONTEXT_DEPENDENT) {
    if (grepl(tolower(ctx), mod_lower, fixed = TRUE)) return("context")
  }
  if (grepl("unannotated mass-shift", mod_lower, fixed = TRUE)) return("artifact")
  return("context")
}

# === MODULE 1: EVOLUTION ======================================================
get_evolutionary_profile_shiny <- function(uniprot_id, gene_name, logger, progress_update) {
  require_pkg("Biostrings")
  require_pkg("msa")
  logger(paste0(">>> [EVO] Analysis for ", gene_name))
  progress_update("Conservation: Fetching human sequence...")
  human_seq <- fetch_uniprot_seq(uniprot_id)
  if(is.null(human_seq)) return(NULL)
  
  target_map <- list("gorilla_gorilla"="Gorilla", "marmota_marmota_marmota"="Marmot", "panthera_pardus"="Leopard", "ursus_maritimus"="Polar Bear", "sus_scrofa"="Pig", "equus_caballus"="Horse", "tursiops_truncatus"="Dolphin", "balaenoptera_musculus"="Blue Whale", "monodelphis_domestica"="Opossum", "meleagris_gallopavo"="Turkey", "chelonoidis_abingdonii"="Tortoise")
  seq_list <- list(Human = human_seq)
  ens_url <- paste0("https://rest.ensembl.org/homology/symbol/human/", gene_name)
  resp_ens <- httr::GET(ens_url, query=list(content_type="application/json", type="orthologues"))
  
  if(httr::status_code(resp_ens) == 200) {
    json_dat <- jsonlite::fromJSON(httr::content(resp_ens, "text", encoding="UTF-8"))
    if(length(json_dat$data) > 0 && "homologies" %in% names(json_dat$data)) {
      homologies <- json_dat$data$homologies[[1]]
      if("target" %in% names(homologies)) homologies <- homologies$target
      cols <- colnames(homologies); species_col <- cols[grep("species$", cols)][1]; id_col <- cols[grep("protein_id$", cols)][1]
      found_ids <- list()
      for(i in 1:nrow(homologies)) {
        sp <- homologies[i, species_col]; id <- homologies[i, id_col]
        if(sp %in% names(target_map) && !is.na(id) && id != "") found_ids[[target_map[[sp]]]] <- id
      }
      logger(paste0("... Found ", length(found_ids), " orthologs."))
      
      if(length(found_ids) > 0) {
        progress_update(paste0("Conservation: Downloading ", length(found_ids), " ortholog sequences..."))
        logger(paste0("... Batch downloading ", length(found_ids), " sequences..."))
        ids_to_fetch <- as.character(unlist(found_ids))
        id_to_species <- setNames(names(found_ids), ids_to_fetch)
        
        max_retries <- 3
        for(attempt in 1:max_retries) {
          tryCatch({
            resp_bulk <- httr::POST(
              url = "https://rest.ensembl.org/sequence/id",
              body = list(ids = ids_to_fetch),
              encode = "json",
              httr::add_headers("Content-Type" = "application/json", "Accept" = "application/json")
            )
            
            if(httr::status_code(resp_bulk) == 429 && attempt < max_retries) {
              wait_time <- 2^attempt
              logger(paste0("... Ensembl rate limited (HTTP 429). Waiting ", wait_time, "s before retry (attempt ", attempt, "/", max_retries, ")..."))
              Sys.sleep(wait_time)
              next
            }
            
            if(httr::status_code(resp_bulk) == 200) {
              bulk_data <- httr::content(resp_bulk, as = "parsed")
              for(item in bulk_data) {
                if(!is.null(item$id) && !is.null(item$seq)) {
                  sp_name <- id_to_species[[item$id]]
                  if(!is.null(sp_name) && nchar(item$seq) > 50) {
                    seq_list[[sp_name]] <- item$seq
                  }
                }
              }
              break
            } else {
              logger(paste0("!!! Ensembl Bulk Error: ", httr::status_code(resp_bulk)))
              break
            }
          }, error = function(e) logger(paste0("!!! Bulk Download Failed: ", e$message)))
        }
      }
    }
  }
  
  if(length(seq_list) < 3) return(NULL)
  
  tryCatch({
    progress_update(paste0("Conservation: Running ClustalOmega MSA (", length(seq_list), " sequences)..."))
    aa_set <- Biostrings::AAStringSet(unlist(seq_list))
    my_aln <- msa::msa(aa_set, method = "ClustalOmega")
    mat <- as.matrix(my_aln); mode(mat) <- "character"
    if(is.null(rownames(mat)) && nrow(mat)==length(seq_list)) rownames(mat) <- names(seq_list)
    if(! "Human" %in% rownames(mat)) return(NULL)
    human_row <- mat["Human", ]; valid_indices <- which(human_row != "-")
    final_scores <- numeric(length(valid_indices))
    other <- setdiff(rownames(mat), "Human")
    for(i in seq_along(valid_indices)) {
      col_idx <- valid_indices[i]
      matches <- sum(mat[other, col_idx] == human_row[col_idx])
      final_scores[i] <- if(length(other)>0) matches/length(other) else 0
    }
    return(data.frame(Position = 1:length(final_scores), Conservation = round(final_scores, 3)))
  }, error = function(e) return(NULL))
}

# === MODULE 2: ACCESSIBILITY (PURE R / SASA / SHRAKE-RUPLEY) ==================
calculate_sasa_r_only <- function(pdb_path, probe_radius = 1.4, n_points = 300) {
  if(!file.exists(pdb_path)) stop("PDB file not found")
  require_pkg("bio3d")
  
  pdb <- suppressWarnings(bio3d::read.pdb(pdb_path))
  
  vdw_radii <- c("C" = 1.70, "N" = 1.55, "O" = 1.52, "S" = 1.80, "H" = 1.20, "P" = 1.80, "FE" = 1.80)
  get_radius <- function(ele) {
    ele <- gsub("[0-9 ]", "", ele)
    r <- vdw_radii[ele]
    if (is.na(r)) r <- 1.80
    return(r)
  }
  
  indices <- seq(0, n_points - 1) + 0.5
  phi <- acos(1 - 2 * indices / n_points)
  theta <- pi * (1 + 5^0.5) * indices
  sphere_points <- cbind(x = cos(theta) * sin(phi), y = sin(theta) * sin(phi), z = cos(phi))
  
  atoms <- pdb$atom[pdb$atom$elesy != "H", ] 
  coords <- as.matrix(atoms[, c("x", "y", "z")])
  n_atoms <- nrow(coords)
  
  elements <- atoms$elesy
  effective_radii <- sapply(elements, get_radius) + probe_radius
  atom_areas <- numeric(n_atoms)
  
  for (i in 1:n_atoms) {
    current_center <- coords[i, ]
    current_radius <- effective_radii[i]
    dists_sq <- rowSums((t(t(coords) - current_center))^2)
    neighbor_indices <- which(dists_sq > 0 & dists_sq < (current_radius + effective_radii)^2)
    
    if (length(neighbor_indices) == 0) {
      atom_areas[i] <- 4 * pi * (current_radius^2)
      next
    }
    
    test_points <- t(t(sphere_points) * current_radius + current_center)
    is_accessible <- rep(TRUE, n_points)
    
    for (j in neighbor_indices) {
      r_j_sq <- effective_radii[j]^2
      d_test_j <- rowSums((t(t(test_points) - coords[j,]))^2)
      is_accessible[d_test_j < r_j_sq] <- FALSE
      if (!any(is_accessible)) break
    }
    atom_areas[i] <- (sum(is_accessible) / n_points) * (4 * pi * current_radius^2)
  }
  
  atoms$sasa <- atom_areas
  res_sasa <- aggregate(sasa ~ resno + resid + chain, data = atoms, sum)
  res_sasa <- res_sasa[order(res_sasa$chain, res_sasa$resno), ]
  colnames(res_sasa) <- c("Position", "Residue", "Chain", "SASA_Abs")
  res_sasa$SASA_Abs <- round(res_sasa$SASA_Abs, 2)
  return(res_sasa)
}

get_accessibility_shiny <- function(uniprot_id, logger, progress_update) {
  require_pkg("bio3d")
  logger(paste0(">>> [SASA] Structural Analysis (Pure R) for ", uniprot_id))
  
  progress_update("SASA: Downloading PDB structure...")
  dest <- download_pdb_r_global(uniprot_id)
  if(is.null(dest)) { logger("!!! [SASA] No PDB found."); return(NULL) }
  
  tryCatch({
    progress_update("SASA: Parsing PDB coordinates...")
    pdb <- suppressWarnings(bio3d::read.pdb(dest))
    pdb_resno <- pdb$atom[pdb$calpha, "resno"]
    xyz <- pdb$xyz[bio3d::atom.select(pdb, "calpha")$xyz]
    d_mat <- as.matrix(dist(matrix(xyz, ncol=3, byrow=TRUE)))
    contact_metric <- 20 - (rowSums(d_mat < 10) - 1)
    
    logger("... Calculating Solvent Accessibility (Shrake-Rupley)...")
    progress_update("SASA: Shrake-Rupley calculation (this may take a moment)...")
    sasa_df <- calculate_sasa_r_only(dest)
    
    final_df <- data.frame(Position = pdb_resno, Contact_Pred = contact_metric)
    final_df$DSSP <- 0
    match_idx <- match(final_df$Position, sasa_df$Position)
    final_df$DSSP[!is.na(match_idx)] <- sasa_df$SASA_Abs[na.omit(match_idx)]
    
    return(final_df)
    
  }, error = function(e) {
    logger(paste("!!! [SASA] Calculation Error:", e$message))
    return(NULL)
  })
}

# === MODULE 3: TKL PREDICTOR (NATIVE R) =======================================
get_tkl_shiny <- function(uniprot_id, gene_name, acc_df, cons_df, logger, progress_update, specific_positions = NULL) {
  logger(">>> [TKL] Prediction (Native R)")
  
  ensure_global_resources()
  if (!GLOBAL_RESOURCES$LOADED) {
    logger("!!! [TKL] Matrices not loaded. Check Settings.")
    return(NULL)
  }
  
  if (!dir.exists(CACHE_DIR)) {
    dir.create(CACHE_DIR, recursive = TRUE)
  }
  cache_file <- file.path(CACHE_DIR, paste0(gene_name, "_TKL_Full_Preds.csv"))
  clean_data <- NULL
  
  if(is.null(specific_positions)) {
    if(file.exists(cache_file)) {
      if (file.info(cache_file)$size > 100) {
        logger("... Loading from cache...")
        clean_data <- read.csv(cache_file, stringsAsFactors = FALSE)
      } else { unlink(cache_file) }
    }
  }
  
  if (is.null(clean_data)) {
    progress_update("TKL: Fetching Sequence...")
    seq_str <- fetch_uniprot_seq(uniprot_id)
    if(is.null(seq_str)) { logger("!!! [TKL] Seq Error"); return(NULL) }
    
    progress_update("TKL: Calculating Scores...")
    seq_chars <- unlist(strsplit(seq_str, ""))
    seq_len <- length(seq_chars)

    # Pad sequence with dummy residues to prevent terminal window truncation
    PAD <- 5L  # max window radius across all PSSM modes
    seq_chars_padded <- c(rep("X", PAD), seq_chars, rep("X", PAD))
    seq_len_padded   <- length(seq_chars_padded)
    
    all_st_indices <- which(seq_chars %in% c("S", "T"))
    all_y_indices  <- which(seq_chars == "Y")
    
    if (!is.null(specific_positions)) {
      all_st_indices <- intersect(all_st_indices, specific_positions)
      all_y_indices  <- intersect(all_y_indices, specific_positions)
    }
    
    total_sites <- length(all_st_indices) + length(all_y_indices)
    progress_update(paste0("TKL: Scoring ", total_sites, " phosphosites across all kinases..."))
    
    final_results <- list()
    
    run_vectorized <- function(indices, mode) {
      if (length(indices) == 0) return(NULL)
      MAT <- if(mode == "ST") GLOBAL_RESOURCES$ST_MAT else GLOBAL_RESOURCES$Y_MAT
      if(mode == "ST") { left <- 5; right <- 4 } else { left <- 5; right <- 5 }
      out_list <- list()
      
      for (i in indices) {
        residue <- seq_chars[i]
        i_pad <- i + PAD  # shift into padded coordinate space
        cols_to_fetch <- c()
        offsets <- c(seq(-left, -1), seq(1, right))
        for (off in offsets) {
          p <- i_pad + off
          aa <- seq_chars_padded[p]
          if (aa != "X") {
            col_name <- paste0(off, aa)
            if (col_name %in% colnames(MAT)) cols_to_fetch <- c(cols_to_fetch, col_name)
          }
        }
        if(length(cols_to_fetch) > 0) {
          sub_mat <- MAT[, cols_to_fetch, drop=FALSE]
          raw_scores <- rowSums(sub_mat, na.rm=TRUE)
        } else { raw_scores <- rep(0, nrow(MAT)) }
        
        if (mode == "ST") {
          bias_vec <- if(residue == "S") GLOBAL_RESOURCES$BIAS_S else GLOBAL_RESOURCES$BIAS_T
          final_scores <- raw_scores + bias_vec
        } else { final_scores <- raw_scores }
        
        df_site <- data.frame(
          Protein = uniprot_id, Gene = gene_name,
          Position = i, Residue = residue, Kinase = names(final_scores),
          Score = round(unname(final_scores), 4), Percentile = NA, 
          stringsAsFactors = FALSE
        )
        out_list[[length(out_list)+1]] <- df_site
      }
      return(bind_rows(out_list))
    }
    
    if(length(all_st_indices) > 0) {
      progress_update(paste0("TKL: Scoring ", length(all_st_indices), " Ser/Thr sites..."))
      final_results[[1]] <- run_vectorized(all_st_indices, "ST")
    }
    if(length(all_y_indices) > 0) {
      progress_update(paste0("TKL: Scoring ", length(all_y_indices), " Tyr sites..."))
      final_results[[2]] <- run_vectorized(all_y_indices, "Y")
    }
    
    clean_data <- bind_rows(final_results) %>% filter(Score > 0)
    
    if(is.null(specific_positions)) write.csv(clean_data, cache_file, row.names=FALSE)
  }
  
  if(is.null(clean_data) || nrow(clean_data)==0) return(NULL)
  
  clean_data$Position <- as.integer(clean_data$Position)
  if(!is.null(acc_df)) clean_data$DSSP <- acc_df$DSSP[match(clean_data$Position, acc_df$Position)]
  if(!is.null(cons_df)) clean_data$Conservation <- cons_df$Conservation[match(clean_data$Position, cons_df$Position)]
  
  return(clean_data)
}

# === MODULE 5: PSP & 14-3-3 ==================================================
get_psp_shiny <- function(uniprot_id, gene_symbol, acc_df, cons_df, logger, progress_update, selected_evidence=c("LT_LIT", "MS_LIT", "MS_CST")) {
  logger(">>> [PSP] Querying Known Sites")
  if (!psp_data_available()) { logger("!!! [PSP] PhosphoSitePlus files not found in data/."); return(NULL) }
  raw_sites <- read_psp_table(PATH_SITES)
  my_sites <- raw_sites %>% filter(GENE == gene_symbol, ORGANISM == "human")
  if(nrow(my_sites) == 0) return(NULL)
  
  clean_sites <- my_sites %>% mutate(Residue = substr(MOD_RSD, 1, 1), Position = as.numeric(str_extract(MOD_RSD, "[0-9]+")))
  cols <- c("LT_LIT", "MS_LIT", "MS_CST")
  for(c in cols) if(!(c %in% names(clean_sites))) clean_sites[[c]] <- 0
  clean_sites$Score <- rowSums(clean_sites[, cols], na.rm=TRUE)
  clean_sites <- clean_sites %>% filter(Score > 0) %>% select(Position, Residue, Score) %>% distinct()
  
  raw_kinases <- read_psp_table(PATH_KINASE_SUB)
  colnames(raw_kinases) <- c("GENE", "KINASE", "KIN_ACC", "KIN_ORG", "SUB_NAME", "SUB_GENE_ID", "SUB_ACC", "SUB_GENE", "SUB_ORG", "SUB_MOD", "GRP", "SEQ", "DOM", "INVIVO", "INVITRO", "CST")
  
  clean_kinases <- raw_kinases %>% filter(SUB_GENE == gene_symbol, SUB_ORG == "human") %>% mutate(Position = as.numeric(str_extract(SUB_MOD, "[0-9]+")), Residue = substr(SUB_MOD, 1, 1)) %>% select(KINASE, Position, Residue) %>% distinct()
  
  merged <- clean_sites %>% left_join(clean_kinases, by=c("Position", "Residue")) %>% group_by(Position, Residue, Score) %>% summarise(Kinase = paste(unique(na.omit(KINASE)), collapse=", "), .groups="drop") %>% mutate(Protein = uniprot_id)
  
  if(!is.null(acc_df)) merged$DSSP <- acc_df$DSSP[match(merged$Position, acc_df$Position)]
  if(!is.null(cons_df)) merged$Conservation <- cons_df$Conservation[match(merged$Position, cons_df$Position)]
  return(merged)
}

get_1433_data <- function(uniprot_id, logger) {
  safe_id <- sanitize_uniprot_id(uniprot_id)
  if (nchar(safe_id) == 0) return(NULL)
  logger(paste0(">>> [14-3-3] Fetching predictions for ", safe_id))
  url_1433 <- paste0("http://www.compbio.dundee.ac.uk/1433pred/pid=", safe_id, "&out=json")
  tryCatch({
    seq_str <- fetch_uniprot_seq(safe_id)
    if(is.null(seq_str)) { logger("!!! [14-3-3] Sequence fetch failed."); return(NULL) }
    json_out <- fromJSON(url_1433)
    df_pred <- as.data.frame(json_out)
    if (nrow(df_pred) == 0) return(NULL)
    positions <- as.numeric(gsub("[^0-9]", "", df_pred$Site))
    residues <- substring(seq_str, positions, positions)
    out <- data.frame(
      Position = positions, Residue = residues, Score = as.numeric(df_pred$Consensus),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$Score), ]
    return(out)
  }, error = function(e) { logger(paste("!!! [14-3-3] Error:", e$message)); return(NULL) })
}

# === MODULE 6: MASTER PLOTTER ================================================
# Returns a plotly object with interactive hover on lollipop points.
# Visual layout:
#   - Protein chain backbone (grey bar)
#   - Colored domain rectangles
#   - Lollipop stems (segments) colored by Source
#   - Lollipop heads (markers) with:
#       border color = Source (blue/red)
#       fill color   = Conservation gradient (magenta→white→green)
#       size         = SASA accessibility
#   - Top-site text annotations
#   - Known-kinase italic annotations

# Helper: map conservation (0-1) to a magenta→white→green color
.cons_color <- function(cons) {
  # NA conservation → neutral light gray (data unavailable)
  is_missing <- is.na(cons)
  cons <- pmin(1, pmax(0, ifelse(is_missing, 0, cons)))
  r <- ifelse(cons < 0.5, 255, as.integer(255 * (1 - cons) * 2))
  g <- ifelse(cons < 0.5, as.integer(255 * cons * 2), as.integer(128 + 127 * (cons - 0.5) * 2))
  b <- ifelse(cons < 0.5, as.integer(255 * cons * 2), 0)
  ifelse(is_missing, "rgb(200,200,200)", paste0("rgb(", r, ",", g, ",", b, ")"))
}

# Helper: map DSSP (already transformed) to marker size in px
.dssp_to_size <- function(dssp, size_lim, vis_range) {
  lo <- size_lim[1]; hi <- size_lim[2]
  # NOTE: use if/else (not ifelse)  -" ifelse(scalar, vector, ...) truncates to length 1
  frac <- if (hi > lo) (dssp - lo) / (hi - lo) else rep(0, length(dssp))
  frac <- pmin(1, pmax(0, frac))
  vis_range[1] + frac * (vis_range[2] - vis_range[1])
}

plot_master_shiny <- function(pred_data, known_data, struct_df, uniprot_id, query_val, mode_type, selected_domains=NULL,
                              filter_score_pct=0, filter_cons_pct=0, filter_acc_pct=0) {
  
  # ---------- DATA PREPARATION (unchanged logic) ----------
  scale_ratio <- 1
  
  if(mode_type == "1433") {
    combined_data <- pred_data %>% 
      mutate(Source = "Prediction", Score_Scaled = Score, Score_Raw = Score, Plot_Pos = Position)
    subtitle_text <- paste0("<span style='color:blue'>14-3-3 Predictions (Dundee)</span> | ", ATTR_1433)
    y_label <- "Consensus Score"
    
  } else {
    if(!is.null(pred_data)) pred_data <- pred_data %>% filter(Protein == uniprot_id)
    if(!is.null(known_data)) known_data <- known_data %>% filter(Protein == uniprot_id)
    
    combined_data <- NULL
    
    if (grepl("Combined", mode_type)) {
      pred_clean <- pred_data %>% 
        filter(grepl(toupper(query_val), toupper(Kinase), fixed=TRUE), Score > 0) %>% 
        mutate(Source="Prediction")
      known_clean <- known_data %>% filter(Score > 0) %>% mutate(Source="Known")
      
      if (grepl("Log2", mode_type)) {
        pred_clean$Score <- pmax(0, pred_clean$Score) 
        known_clean$Score <- pmax(0, log2(known_clean$Score))
        y_label <- paste0("Log2(", query_val, " Prediction)")
      } else {
        y_label <- paste0(query_val, " Prediction")
      }
      
      if(nrow(pred_clean) == 0 && nrow(known_clean) == 0) return(NULL)
      
      max_pred <- max(pred_clean$Score, na.rm=TRUE); if(is.infinite(max_pred)) max_pred <- 5
      max_known <- max(known_clean$Score, na.rm=TRUE); if(is.infinite(max_known)) max_known <- 1
      scale_ratio <- (max_pred * 0.75) / max_known
      known_clean$Score_Scaled <- known_clean$Score * scale_ratio
      pred_clean$Score_Scaled <- pred_clean$Score
      
      combined_data <- bind_rows(pred_clean, known_clean) %>% mutate(Plot_Pos = Position)
      subtitle_text <- paste0("<span style='color:blue'>Predicted (", query_val, ")</span> &amp; <span style='color:darkred'>Known Sites</span>")
      
    } else {
      plot_data <- if(!is.null(pred_data)) pred_data else known_data
      if(grepl("TKL", mode_type)) {
        plot_data <- plot_data %>% filter(grepl(toupper(query_val), toupper(Kinase), fixed=TRUE), Score > 0)
        y_label <- paste0(query_val, " Score")
        combined_data <- plot_data %>% mutate(Source="Prediction", Score_Scaled=Score, Score_Raw=Score, Plot_Pos=Position)
      } else {
        combined_data <- plot_data %>% mutate(Source="Known", Score_Scaled=Score, Score_Raw=Score, Plot_Pos=Position)
        y_label <- "Literature Evidence"
      }
      subtitle_text <- paste0("Source: ", mode_type)
    }
  }
  
  if(is.null(combined_data) || nrow(combined_data) == 0) return(NULL)
  
  # Ensure Conservation and DSSP columns exist (may be absent if cons_df/acc_df were NULL)
  if (!"Conservation" %in% names(combined_data)) combined_data$Conservation <- NA_real_
  if (!"DSSP"         %in% names(combined_data)) combined_data$DSSP         <- NA_real_
  
  combined_data <- combined_data %>% 
    group_by(Plot_Pos, Residue, Source) %>% 
    summarise(
      Score_Scaled = Score_Scaled[which.max(abs(Score_Scaled))],
      Score = Score[which.max(abs(Score))],
      Kinase = dplyr::first(Kinase),
      Conservation = dplyr::first(Conservation),
      DSSP = dplyr::first(DSSP),
      .groups = "drop"
    )
  
  # ---------- PERCENTILE FILTERING ----------
  if (filter_score_pct > 0 && nrow(combined_data) > 0) {
    score_vals <- combined_data$Score_Scaled
    threshold <- quantile(score_vals, probs = filter_score_pct / 100, na.rm = TRUE)
    combined_data <- combined_data %>% filter(Score_Scaled >= threshold)
  }
  if (filter_cons_pct > 0 && nrow(combined_data) > 0) {
    cons_vals <- combined_data$Conservation
    if (!all(is.na(cons_vals))) {
      threshold <- quantile(cons_vals, probs = filter_cons_pct / 100, na.rm = TRUE)
      combined_data <- combined_data %>% filter(!is.na(Conservation) & Conservation >= threshold)
    }
  }
  if (filter_acc_pct > 0 && nrow(combined_data) > 0) {
    acc_vals <- combined_data$DSSP
    if (!all(is.na(acc_vals))) {
      threshold <- quantile(acc_vals, probs = filter_acc_pct / 100, na.rm = TRUE)
      combined_data <- combined_data %>% filter(!is.na(DSSP) & DSSP >= threshold)
    }
  }
  if (nrow(combined_data) == 0) return(NULL)

  max_abs_score <- max(abs(combined_data$Score_Scaled), na.rm=TRUE)
  domain_h <- max_abs_score * 0.15
  chain_h  <- domain_h * 0.25
  
  # DSSP transformation
  dssp_raw_present <- "DSSP" %in% names(combined_data) && !all(is.na(combined_data$DSSP))
  if(!dssp_raw_present) {
    combined_data$DSSP <- 0
    combined_data$DSSP_Original <- 0
    size_lim <- c(0, 1); vis_range <- c(4, 8)
  } else {
    combined_data$DSSP_Original <- pmin(combined_data$DSSP, 150)
    combined_data$DSSP <- pmin(combined_data$DSSP, 150) ^ 1.5
    size_lim <- c(0, 150^1.5); vis_range <- c(4, 22)
  }
  
  # Domain filtering
  domains <- subset(struct_df, type %in% c("DOMAIN", "REGION", "MOTIF", "REPEAT"))
  if(nrow(domains) > 0) domains <- domains[order(domains$begin), ]
  
  if(!is.null(selected_domains) && nrow(domains) > 0) {
    domains$selection_label <- paste0(domains$description, " [", domains$begin, "-", domains$end, "]")
    domains <- domains[domains$selection_label %in% selected_domains, ]
  } else { domains <- domains[FALSE, ] }
  
  if(nrow(domains) > 0) {
    ordered_levels <- unique(domains$description)
    domains$description <- factor(domains$description, levels = ordered_levels)
  }
  
  # Top sites for text annotations
  top_sites <- combined_data %>% mutate(label_txt = paste0(Residue, Plot_Pos))
  if(mode_type == "1433") {
    top_sites <- top_sites %>% arrange(desc(Score)) %>% slice_head(n=5)
  } else if (grepl("Combined", mode_type)) {
    preds <- top_sites %>% filter(Source=="Prediction") %>% arrange(desc(Score_Scaled)) %>% slice_head(n=5)
    known <- top_sites %>% filter(Source=="Known")
    top_sites <- bind_rows(preds, known)
  } else {
    top_sites <- top_sites %>% arrange(desc(abs(Score))) %>% slice_head(n=5)
  }
  
  kinase_labels <- combined_data %>% filter(Source == "Known", !is.na(Kinase), Kinase != "")
  
  # ---------- PLOTLY CONSTRUCTION ----------
  source_colors <- c("Prediction" = "blue", "Known" = "red")
  
  # Compute marker aesthetics
  combined_data$marker_border <- ifelse(combined_data$Source == "Prediction", "blue", "red")
  combined_data$marker_fill   <- .cons_color(combined_data$Conservation)
  combined_data$marker_size   <- .dssp_to_size(combined_data$DSSP, size_lim, vis_range)
  
  # Build hover text
  combined_data$hover_text <- paste0(
    "<b>", combined_data$Residue, combined_data$Plot_Pos, "</b><br>",
    "Source: ", combined_data$Source, "<br>",
    "Score: ", round(combined_data$Score, 3), "<br>",
    "Conservation: ", ifelse(is.na(combined_data$Conservation), "N/A",
                              round(combined_data$Conservation, 3)), "<br>",
    "SASA: ", if (dssp_raw_present) round(combined_data$DSSP_Original, 1) else "N/A", "<br>",
    ifelse(!is.na(combined_data$Kinase) & combined_data$Kinase != "",
           paste0("Kinase: ", combined_data$Kinase), "")
  )
  
  # Protein length from CHAIN
  chain_data <- subset(struct_df, type == "CHAIN")
  protein_length <- if (nrow(chain_data) > 0) max(chain_data$end) else max(combined_data$Plot_Pos, na.rm = TRUE)
  
  p <- plot_ly()
  
  # --- Layer 1: Protein chain backbone (grey bar) ---
  if (nrow(chain_data) > 0) {
    for (i in seq_len(nrow(chain_data))) {
      p <- p %>% add_trace(
        type = "scatter", mode = "lines", x = c(chain_data$begin[i], chain_data$end[i]),
        y = c(0, 0), line = list(width = 10, color = "rgba(200,200,200,0.7)"),
        showlegend = FALSE, hoverinfo = "text",
        text = paste0("Chain: ", chain_data$description[i],
                       " (", chain_data$begin[i], "\u2013", chain_data$end[i], ")")
      )
    }
  }
  
  # --- Layer 2: Domain rectangles ---
  domain_palette <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                       "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
                       "#ccebc5", "#ffed6f")
  if (nrow(domains) > 0) {
    dom_names <- unique(as.character(domains$description))
    dom_colors <- setNames(
      domain_palette[((seq_along(dom_names) - 1) %% length(domain_palette)) + 1],
      dom_names
    )
    shapes <- list()
    annotations <- list()
    for (i in seq_len(nrow(domains))) {
      d <- domains[i, ]
      shapes[[length(shapes) + 1]] <- list(
        type = "rect", x0 = d$begin, x1 = d$end,
        y0 = -domain_h, y1 = domain_h,
        fillcolor = dom_colors[as.character(d$description)],
        opacity = 0.8, line = list(color = "black", width = 1),
        layer = "below"
      )
    }
    # Add domain legend entries (invisible markers with matching fill)
    for (dname in dom_names) {
      p <- p %>% add_trace(
        type = "scatter", mode = "markers",
        x = c(NA), y = c(NA),
        marker = list(size = 12, color = dom_colors[dname],
                      symbol = "square",
                      line = list(color = "black", width = 1)),
        name = dname, legendgroup = paste0("dom_", dname),
        showlegend = TRUE, hoverinfo = "none"
      )
    }
  } else {
    shapes <- list()
    annotations <- list()
  }
  
  # --- Layer 3 & 4: Lollipop stems + heads per Source ---
  for (src in unique(combined_data$Source)) {
    src_data <- combined_data %>% filter(Source == src)
    src_color <- source_colors[src]
    
    # Stems: one invisible segment trace per point (for the line)
    stem_x <- as.numeric(rbind(src_data$Plot_Pos, src_data$Plot_Pos, NA))
    stem_y <- as.numeric(rbind(rep(0, nrow(src_data)), src_data$Score_Scaled, NA))
    p <- p %>% add_trace(
      type = "scatter", mode = "lines",
      x = stem_x, y = stem_y,
      line = list(color = src_color, width = 1),
      showlegend = FALSE, hoverinfo = "skip"
    )
    
    # Heads: markers with conservation fill, source border, DSSP size
    p <- p %>% add_trace(
      type = "scatter", mode = "markers",
      x = src_data$Plot_Pos, y = src_data$Score_Scaled,
      marker = list(
        size = src_data$marker_size,
        sizemode = "diameter",
        color = src_data$marker_fill,
        line = list(color = src_color, width = 1.5),
        opacity = 0.9
      ),
      name = src,
      legendgroup = src,
      text = src_data$hover_text,
      hoverinfo = "text"
    )
  }
  
  # --- Layer 5: Top-site text annotations ---
  for (i in seq_len(nrow(top_sites))) {
    ts <- top_sites[i, ]
    label_color <- source_colors[ts$Source]
    nudge <- if (ts$Score_Scaled > 0) max_abs_score * 0.08 else -max_abs_score * 0.08
    annotations[[length(annotations) + 1]] <- list(
      x = ts$Plot_Pos, y = ts$Score_Scaled + nudge,
      text = paste0(ts$Residue, ts$Plot_Pos),
      showarrow = TRUE, arrowhead = 0, arrowwidth = 0.5, arrowcolor = "#999",
      ax = 0, ay = if (ts$Score_Scaled > 0) -20 else 20,
      font = list(size = 11, color = label_color),
      xanchor = "center"
    )
  }
  
  # --- Layer 6: Known-kinase italic labels ---
  if (nrow(kinase_labels) > 0) {
    for (i in seq_len(nrow(kinase_labels))) {
      kl <- kinase_labels[i, ]
      annotations[[length(annotations) + 1]] <- list(
        x = kl$Plot_Pos, y = kl$Score_Scaled,
        text = paste0("<i>", kl$Kinase, "</i>"),
        showarrow = TRUE, arrowhead = 0, arrowwidth = 0.5, arrowcolor = "#cc5555",
        ax = 35, ay = -25,
        font = list(size = 10, color = "red"),
        xanchor = "left"
      )
    }
  }
  
  # --- Conservation colour-bar legend ---
  # Compute y_range early (needed by SASA legend positioning below)
  y_range <- c(min(0, min(combined_data$Score_Scaled, na.rm = TRUE)) * 1.15,
                max_abs_score * 1.25)

  # Use the actual lollipop marker data with a colorbar
  # We add one combined marker trace that carries the colorbar
  p <- p %>% add_trace(
    type = "scatter", mode = "markers",
    x = combined_data$Plot_Pos,
    y = combined_data$Score_Scaled,
    marker = list(
      size = 0.1, opacity = 0,
      color = ifelse(is.na(combined_data$Conservation), 0.5, combined_data$Conservation),
      colorscale = list(c(0, "magenta"), c(0.5, "white"), c(1, "darkgreen")),
      cmin = 0, cmax = 1,
      colorbar = list(title = list(text = "Conservation", font = list(size = 11)),
                       x = 1.02, y = 0.5, len = 0.4, thickness = 12,
                       tickvals = c(0, 0.5, 1), ticktext = c("Low", "0.5", "High"))
    ),
    showlegend = FALSE, hoverinfo = "skip"
  )

  # --- SASA size legend (annotation-based with example circles) ---
  sasa_breaks   <- c(10, 50, 100, 150)
  sasa_labels_t <- c("10", "50", "100", "150+")
  sasa_sizes    <- .dssp_to_size(sasa_breaks ^ 1.5, size_lim, vis_range)
  # Place legend circles off-screen right, stacked vertically
  legend_x <- protein_length * 1.08
  legend_y_start <- y_range[2] * 0.95
  legend_y_step  <- max_abs_score * 0.12

  # Title annotation for SASA legend
  annotations[[length(annotations) + 1]] <- list(
    x = legend_x, y = legend_y_start + legend_y_step,
    text = "<b>SASA (\u00c5\u00b2)</b>", showarrow = FALSE,
    font = list(size = 10, color = "#333"), xanchor = "left"
  )
  for (i in seq_along(sasa_breaks)) {
    lx <- legend_x
    ly <- legend_y_start - (i - 1) * legend_y_step
    p <- p %>% add_trace(
      type = "scatter", mode = "markers+text",
      x = c(lx), y = c(ly),
      marker = list(size = sasa_sizes[i], sizemode = "diameter",
                     color = "rgba(180,180,180,0.5)",
                     line = list(color = "grey", width = 1)),
      text = sasa_labels_t[i], textposition = "middle right",
      textfont = list(size = 9, color = "#555"),
      showlegend = FALSE, hoverinfo = "skip"
    )
  }
  
  # --- Layout ---
  
  # Title
  plot_title <- paste0("Profile: ", uniprot_id)
  
  # Dual y-axis for Combined mode
  if (grepl("Combined", mode_type)) {
    p <- p %>% layout(
      title = list(text = paste0(plot_title, "<br><sup>", subtitle_text, "</sup>"),
                    font = list(size = 15)),
      xaxis = list(title = "Residue Position",
                    autorange = TRUE, showgrid = TRUE, zeroline = FALSE),
      yaxis = list(title = list(text = y_label, font = list(color = "blue")),
                    tickfont = list(color = "blue"),
                    autorange = TRUE, showgrid = TRUE, zeroline = FALSE),
      yaxis2 = list(title = list(text = "Known Evidence", font = list(color = "red")),
                     tickfont = list(color = "red"),
                     overlaying = "y", side = "right",
                     autorange = TRUE, showgrid = FALSE),
      shapes = shapes, annotations = annotations,
      legend = list(orientation = "h", x = 0, y = -0.15, xanchor = "left"),
      hovermode = "closest",
      margin = list(l = 70, r = 100, t = 80, b = 60)
    )
  } else {
    p <- p %>% layout(
      title = list(text = paste0(plot_title, "<br><sup>", subtitle_text, "</sup>"),
                    font = list(size = 15)),
      xaxis = list(title = "Residue Position",
                    autorange = TRUE, showgrid = TRUE, zeroline = FALSE),
      yaxis = list(title = y_label,
                    autorange = TRUE, showgrid = TRUE, zeroline = FALSE),
      shapes = shapes, annotations = annotations,
      legend = list(orientation = "h", x = 0, y = -0.15, xanchor = "left"),
      hovermode = "closest",
      margin = list(l = 70, r = 100, t = 80, b = 60)
    )
  }
  
  p
}

# === MODULE 7: MULTI-SITE ANALYZER ===========================================
parse_sites <- function(input_str) {
  if (is.null(input_str) || input_str == "") return(NULL)
  parts <- unlist(strsplit(input_str, ","))
  cleaned <- as.numeric(trimws(parts))
  cleaned <- sort(unique(cleaned[!is.na(cleaned)]))
  if (length(cleaned) == 0) return(NULL)
  return(cleaned)
}

# === MODULE 9: SCANSITE =======================================================
get_scansite_predictions <- function(uniprot_id, sequence, logger) {
  safe_id <- sanitize_uniprot_id(uniprot_id)
  if (nchar(safe_id) == 0) return(NULL)
  logger(paste0(">>> [Scansite] Querying MIT Server for ", safe_id, "..."))
  base_url <- "https://scansite4.mit.edu/webservice/proteinscan"
  clean_seq <- gsub("[^A-Z]", "", sequence)
  request_url <- paste0(base_url, "/identifier=", safe_id, "/sequence=", clean_seq, "/motifclass=MAMMALIAN", "/stringency=Low")
  
  tryCatch({
    resp <- GET(request_url, timeout(30))
    if (status_code(resp) != 200) { logger(paste("!!! Scansite Error:", status_code(resp))); return(NULL) }
    
    xml_content <- content(resp, as = "parsed", encoding = "UTF-8")
    sites <- xml_find_all(xml_content, ".//predictedSite")
    
    if (length(sites) == 0) { logger("!!! No motifs found."); return(NULL) }
    
    results <- lapply(sites, function(site) {
      pos_num  <- as.numeric(gsub("[^0-9]", "", xml_text(xml_find_first(site, ".//site"))))
      perc_num <- as.numeric(xml_text(xml_find_first(site, ".//percentile")))
      name_txt <- xml_text(xml_find_first(site, ".//motifName"))
      nick_txt <- xml_text(xml_find_first(site, ".//motifNickname"))
      grp_txt  <- xml_text(xml_find_first(site, ".//motifGroup"))
      
      best_label <- if(!is.na(name_txt) && name_txt != "") name_txt else 
        if(!is.na(nick_txt) && nick_txt != "") nick_txt else 
          if(!is.na(grp_txt)  && grp_txt  != "") grp_txt  else "Unknown"
      
      search_str <- paste(best_label, grp_txt, nick_txt)
      type <- "Upstream Kinase"
      if (any(sapply(BINDER_TERMS, function(x) grepl(x, search_str, ignore.case = TRUE)))) {
        type <- "Downstream Binder"
      }
      
      data.frame(
        Position = pos_num,
        Percentile = perc_num,
        Motif_Name = best_label,
        Category = type,
        stringsAsFactors = FALSE
      )
    })
    
    final_df <- bind_rows(results) %>% filter(!is.na(Position), !is.na(Percentile))
    return(final_df)
    
  }, error = function(e) { logger(paste("!!! Scansite Connection Error:", e$message)); return(NULL) })
}

report_site_status <- function(target_sites, seq_str, found_data) {
  seq_len <- nchar(seq_str)
  status_list <- list()
  for(site in target_sites) {
    if(site > seq_len || site < 1) {
      status_list[[as.character(site)]] <- "Invalid Index"
    } else {
      context <- get_seq_context(seq_str, site)
      hits <- found_data %>% filter(Position == site)
      if(nrow(hits) == 0) {
        status_list[[as.character(site)]] <- paste0("No Predictions. Context: ", context)
      } else {
        binders <- hits %>% filter(Category == "Downstream Binder")
        kinases <- hits %>% filter(Category == "Upstream Kinase")
        desc <- c()
        if(nrow(binders) > 0) desc <- c(desc, paste(nrow(binders), "Binders"))
        if(nrow(kinases) > 0) desc <- c(desc, paste(nrow(kinases), "Kinases"))
        if(length(desc) == 0) desc <- c("Low score hits only")
        status_list[[as.character(site)]] <- paste0("Found: ", paste(desc, collapse=", "), " | Context: ", context)
      }
    }
  }
  return(paste(sapply(names(status_list), function(x) paste0("Site ", x, ": ", status_list[[x]])), collapse = "\n"))
}

# === LOCAL MS: FRAGPIPE UPDATER ===============================================
extract_experiment_id <- function(spectrum_value) {
  sub("\\..*$", "", as.character(spectrum_value))
}

update_fragpipe_db <- function(psm_path, peptide_path, user_id,
                                existing_db = NULL, progress_fn = NULL) {
  
  log_msg <- function(msg) {
    if (!is.null(progress_fn)) progress_fn(msg)
    cat(msg, "\n")
  }
  
  if (!file.exists(psm_path)) {
    return(list(success = FALSE, message = paste0("PSM file not found: ", psm_path), n_entries = 0))
  }
  if (!file.exists(peptide_path)) {
    return(list(success = FALSE, message = paste0("Peptide file not found: ", peptide_path), n_entries = 0))
  }
  if (is.null(user_id) || nchar(trimws(user_id)) == 0) {
    return(list(success = FALSE, message = "UserID cannot be empty.", n_entries = 0))
  }
  user_id <- trimws(user_id)
  
  log_msg("Reading PSM file...")
  psm_df <- read_tsv(psm_path, show_col_types = FALSE)
  log_msg(paste0("  PSM rows: ", nrow(psm_df)))
  
  if (!"Spectrum" %in% colnames(psm_df)) {
    return(list(success = FALSE, message = "PSM file missing 'Spectrum' column.", n_entries = 0))
  }
  
  experiment_id <- extract_experiment_id(psm_df$Spectrum[1])
  log_msg(paste0("  ExperimentID detected: ", experiment_id))
  
  if (!is.null(existing_db) && "ExperimentID" %in% colnames(existing_db)) {
    existing_ids <- unique(existing_db$ExperimentID)
    if (experiment_id %in% existing_ids) {
      return(list(
        success = FALSE,
        message = paste0("ExperimentID '", experiment_id,
                         "' already exists in the database. Import aborted to prevent duplicates."),
        n_entries = 0
      ))
    }
  }
  
  log_msg("Reading peptide file...")
  peptide_df <- read_tsv(peptide_path, show_col_types = FALSE)
  log_msg(paste0("  Peptide rows: ", nrow(peptide_df)))
  
  log_msg("Building peptide modification lookup...")
  peptide_mods_list <- vector("list", nrow(peptide_df))
  for (i in seq_len(nrow(peptide_df))) {
    row <- peptide_df[i, ]
    obs_mods_str <- as.character(row$`Observed Modifications`)
    if (!is.na(obs_mods_str) && obs_mods_str != "") {
      mods <- parse_modifications(obs_mods_str)
      if (length(mods) > 0) {
        peptide_mods_list[[i]] <- tibble(
          Peptide = rep(as.character(row$Peptide), length(mods)),
          Modification = mods,
          Gene = rep(as.character(row$Gene), length(mods))
        )
      }
    }
  }
  peptide_mods <- bind_rows(peptide_mods_list)
  peptide_mods_unique <- peptide_mods %>% distinct(Peptide, Modification, Gene)
  peptide_mod_keys <- paste(peptide_mods_unique$Peptide,
                            peptide_mods_unique$Modification,
                            peptide_mods_unique$Gene, sep = "|")
  peptide_mod_set <- new.env(hash = TRUE, parent = emptyenv())
  for (k in peptide_mod_keys) peptide_mod_set[[k]] <- TRUE
  log_msg(paste0("  Unique peptide-mod-gene combos: ", nrow(peptide_mods_unique)))
  
  log_msg("Processing PSM modifications...")
  results_list <- vector("list", nrow(psm_df))
  entry_count <- 0L
  
  for (i in seq_len(nrow(psm_df))) {
    row <- psm_df[i, ]
    obs_mods_str <- as.character(row$`Observed Modifications`)
    if (is.na(obs_mods_str) || obs_mods_str == "") next
    mods <- parse_modifications(obs_mods_str)
    if (length(mods) == 0) next
    
    peptide <- as.character(row$Peptide)
    primary_gene <- as.character(row$Gene)
    delta_mass <- as.numeric(row$`Delta Mass`)
    protein_start <- as.integer(row$`Protein Start`)
    protein_end <- as.integer(row$`Protein End`)
    peptideprophet_prob <- as.numeric(row$`PeptideProphet Probability`)
    expectation <- as.numeric(row$Expectation)
    hyperscore <- as.numeric(row$Hyperscore)
    best_score_with <- as.numeric(row$`Best Score with Delta Mass`)
    best_score_without <- as.numeric(row$`Best Score without Delta Mass`)
    is_unique_raw <- tolower(as.character(row$`Is Unique`))
    is_unique <- is_unique_raw == "true"
    missed_cleavages <- as.integer(row$`Number of Missed Cleavages`)
    
    localization_str <- as.character(row$`MSFragger Localization`)
    has_localization <- !is.na(localization_str) && nchar(localization_str) > 0
    residue_positions <- extract_residue_positions(localization_str)
    residue_location <- if (length(residue_positions) > 0) {
      paste(residue_positions, collapse = ", ")
    } else { NA_character_ }
    
    modification_location_in_protein <- if (!is.na(protein_start) && length(residue_positions) > 0) {
      paste(protein_start + residue_positions - 1, collapse = ", ")
    } else { NA_character_ }
    
    mapped_genes_str <- as.character(row$`Mapped Genes`)
    mapped_genes <- split_mapped_genes(mapped_genes_str)
    all_genes <- unique(c(primary_gene, if (!is.null(mapped_genes)) mapped_genes))
    
    for (mod in mods) {
      if (mod == "First isotopic peak") next
      mass_error <- parse_mass_error_for_mod(obs_mods_str, mod)
      for (gene in all_genes) {
        lookup_key <- paste(peptide, mod, gene, sep = "|")
        in_peptide_file <- exists(lookup_key, envir = peptide_mod_set, inherits = FALSE)
        other_genes <- all_genes[all_genes != gene]
        additional_hits <- if (length(other_genes) > 0) paste(other_genes, collapse = ", ") else ""
        
        entry_count <- entry_count + 1L
        results_list[[entry_count]] <- data.frame(
          Gene = gene, Peptide = peptide, Modification = mod,
          Residue_Location = residue_location,
          Protein_Start = protein_start, Protein_End = protein_end,
          Modification_Location_in_Protein = modification_location_in_protein,
          Delta_Mass = delta_mass, Additional_Gene_Hits = additional_hits,
          PeptideProphet_Prob = peptideprophet_prob,
          Expectation = expectation, Hyperscore = hyperscore,
          Best_Score_With_DM = best_score_with,
          Best_Score_Without_DM = best_score_without,
          Is_Unique = is_unique, Missed_Cleavages = missed_cleavages,
          Has_Localization = has_localization, Mass_Error = mass_error,
          In_Peptide_File = in_peptide_file,
          stringsAsFactors = FALSE
        )
      }
    }
    if (i %% 2000 == 0) {
      log_msg(paste0("  Processed ", i, " PSM rows (", entry_count, " entries)..."))
    }
  }
  
  fragpipe_db <- as_tibble(bind_rows(results_list[seq_len(entry_count)]))
  log_msg(paste0("  Total raw entries: ", entry_count))
  
  log_msg("Aggregating diagnostics...")
  peptide_totals <- fragpipe_db %>%
    group_by(Gene, Peptide) %>%
    summarise(Total_Peptide_Identified = n(), .groups = "drop")
  mod_totals <- fragpipe_db %>%
    group_by(Gene, Peptide, Modification) %>%
    summarise(Times_Mod_Identified = n(), .groups = "drop")
  
  diagnostics_agg <- fragpipe_db %>%
    group_by(Gene, Peptide, Modification, Residue_Location, Protein_Start, Protein_End,
             Modification_Location_in_Protein) %>%
    summarise(
      Delta_Mass = dplyr::first(Delta_Mass),
      Additional_Gene_Hits = dplyr::first(Additional_Gene_Hits),
      Best_PeptideProphet = max(PeptideProphet_Prob, na.rm = TRUE),
      Best_Expectation = min(Expectation, na.rm = TRUE),
      Best_Hyperscore = max(Hyperscore, na.rm = TRUE),
      Best_Score_With_DM = max(Best_Score_With_DM, na.rm = TRUE),
      Best_Score_Without_DM = min(Best_Score_Without_DM, na.rm = TRUE),
      Is_Unique = any(Is_Unique),
      Min_Missed_Cleavages = min(Missed_Cleavages, na.rm = TRUE),
      Has_Localization = any(Has_Localization),
      Mass_Error = suppressWarnings(min(Mass_Error, na.rm = TRUE)),
      In_Peptide_File = any(In_Peptide_File),
      .groups = "drop"
    ) %>%
    mutate(
      Best_PeptideProphet   = ifelse(is.infinite(Best_PeptideProphet), NA_real_, Best_PeptideProphet),
      Best_Expectation      = ifelse(is.infinite(Best_Expectation), NA_real_, Best_Expectation),
      Best_Hyperscore       = ifelse(is.infinite(Best_Hyperscore), NA_real_, Best_Hyperscore),
      Best_Score_With_DM    = ifelse(is.infinite(Best_Score_With_DM), NA_real_, Best_Score_With_DM),
      Best_Score_Without_DM = ifelse(is.infinite(Best_Score_Without_DM), NA_real_, Best_Score_Without_DM),
      Min_Missed_Cleavages  = ifelse(is.infinite(Min_Missed_Cleavages), NA_integer_, Min_Missed_Cleavages),
      Mass_Error            = ifelse(is.infinite(Mass_Error), NA_real_, Mass_Error)
    )
  
  new_data <- diagnostics_agg %>%
    left_join(peptide_totals, by = c("Gene", "Peptide")) %>%
    left_join(mod_totals, by = c("Gene", "Peptide", "Modification"))
  
  log_msg("Calculating confidence scores...")
  new_data <- new_data %>%
    mutate(
      Mod_Type_Class = sapply(Modification, classify_mod_type),
      C1 = case_when(Mod_Type_Class == "biological" ~ 25, Mod_Type_Class == "context" ~ 12, TRUE ~ 0),
      C2 = ifelse(Has_Localization, 20, 0),
      sr = ifelse(!is.na(Best_Score_With_DM) & !is.na(Best_Score_Without_DM) & Best_Score_Without_DM > 0,
                  Best_Score_With_DM / Best_Score_Without_DM, NA_real_),
      C3 = case_when(
        !is.na(sr) ~ pmin(20, pmax(0, (sr - 1) * 20)),
        !is.na(Best_Score_With_DM) & Best_Score_With_DM > 0 &
          (is.na(Best_Score_Without_DM) | Best_Score_Without_DM == 0) ~ 20,
        TRUE ~ 0
      ),
      C4 = ifelse(!is.na(Best_PeptideProphet), Best_PeptideProphet * 15, 0),
      C5 = ifelse(Is_Unique, 10, 0),
      C6a = case_when(is.na(Mass_Error) ~ 0, Mass_Error < 0.001 ~ 4,
                      Mass_Error < 0.005 ~ 3, Mass_Error < 0.01 ~ 2, Mass_Error < 0.05 ~ 1, TRUE ~ 0),
      C6b = case_when(is.na(Min_Missed_Cleavages) ~ 0, Min_Missed_Cleavages == 0 ~ 3,
                      Min_Missed_Cleavages == 1 ~ 1, TRUE ~ 0),
      C6c = ifelse(In_Peptide_File, 3, 0),
      Confidence_Score = round(C1 + C2 + C3 + C4 + C5 + C6a + C6b + C6c, 1),
      Biological_Confidence = case_when(
        Confidence_Score >= 70 ~ "High", Confidence_Score >= 45 ~ "Mixed",
        Confidence_Score >= 25 ~ "Low", TRUE ~ "Very Low"
      ),
      UserID = user_id,
      ExperimentID = experiment_id
    ) %>%
    select(-C1, -C2, -sr, -C3, -C4, -C5, -C6a, -C6b, -C6c) %>%
    arrange(Gene, Peptide, Modification)
  
  log_msg(paste0("Processed ", nrow(new_data), " entries for ExperimentID: ", experiment_id))
  
  return(list(
    success = TRUE,
    message = paste0("Successfully imported ", nrow(new_data), " entries.\n",
                     "ExperimentID: ", experiment_id, "\n",
                     "UserID: ", user_id, "\n",
                     "Unique genes: ", n_distinct(new_data$Gene), "\n",
                     "Unique modifications: ", n_distinct(new_data$Modification)),
    n_entries = nrow(new_data),
    data      = new_data
  ))
}

# === LOCAL MS: PLOTLY MODIFICATION VIEWER =====================================

MOD_COLORS <- c(
  "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
  "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
  "#dcbeff", "#9A6324", "#fffac8", "#800000", "#aaffc3",
  "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#e6beff"
)

.get_mod_colors <- function(mod_names) {
  n <- length(mod_names)
  colors <- MOD_COLORS[((seq_len(n) - 1) %% length(MOD_COLORS)) + 1]
  setNames(colors, mod_names)
}

.prepare_pos_df <- function(gene_df) {
  gene_df %>%
    rowwise() %>%
    mutate(Position = list(expand_position_string(Modification_Location_in_Protein))) %>%
    unnest(cols = c(Position)) %>%
    ungroup() %>%
    filter(!is.na(Position)) %>%
    mutate(
      .idx = Position - Protein_Start + 1,
      Residue = ifelse(.idx >= 1 & .idx <= nchar(Peptide),
                       substr(Peptide, .idx, .idx), "")
    ) %>%
    select(-.idx)
}

.add_backbone <- function(p, protein_length, show_thresholds = TRUE) {
  p <- p %>% add_segments(
    x = 1, xend = protein_length, y = 0, yend = 0,
    line = list(color = "#333333", width = 2),
    showlegend = FALSE, hoverinfo = "none"
  )
  if (show_thresholds) {
    for (thresh in list(
      list(y = 70, label = "High (70)", color = "rgba(0,128,0,0.3)"),
      list(y = 45, label = "Mixed (45)", color = "rgba(255,165,0,0.3)"),
      list(y = 25, label = "Low (25)",   color = "rgba(255,0,0,0.3)")
    )) {
      p <- p %>% add_segments(
        x = 1, xend = protein_length,
        y = thresh$y, yend = thresh$y,
        line = list(color = thresh$color, width = 1, dash = "dash"),
        showlegend = FALSE, hoverinfo = "text", text = thresh$label
      )
    }
  }
  p
}

.make_annotations <- function(protein_length) {
  list(
    list(x = protein_length + 1, y = 70, text = "High", showarrow = FALSE,
         xanchor = "left", font = list(size = 9, color = "green")),
    list(x = protein_length + 1, y = 45, text = "Mixed", showarrow = FALSE,
         xanchor = "left", font = list(size = 9, color = "orange")),
    list(x = protein_length + 1, y = 25, text = "Low", showarrow = FALSE,
         xanchor = "left", font = list(size = 9, color = "red"))
  )
}

# Fetch UniProt domain/feature data for a given ID
fetch_uniprot_domains <- function(uniprot_id) {
  require_pkg("drawProteins")
  tryCatch({
    feat <- drawProteins::get_features(uniprot_id)
    df <- drawProteins::feature_to_dataframe(feat)
    domains <- subset(df, type %in% c("DOMAIN", "REGION", "MOTIF", "REPEAT"))
    if(nrow(domains) > 0) {
      domains <- domains[order(domains$begin), ]
      domains$label <- paste0(domains$description, " [", domains$begin, "-", domains$end, "]")
    }
    domains
  }, error = function(e) {
    data.frame(type = character(0), begin = integer(0), end = integer(0),
               description = character(0), label = character(0))
  })
}

# Pre-defined domain colors (ColorBrewer Set3 style)
DOMAIN_COLORS <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                   "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
                   "#ccebc5", "#ffed6f")

.get_domain_colors <- function(domain_names) {
  n <- length(domain_names)
  colors <- DOMAIN_COLORS[((seq_len(n) - 1) %% length(DOMAIN_COLORS)) + 1]
  setNames(colors, domain_names)
}

.prepare_coverage_data <- function(gene_df) {
  gene_df %>%
    filter(!is.na(Protein_Start), !is.na(Protein_End)) %>%
    group_by(Modification, Mod_Type_Class, Peptide, Protein_Start, Protein_End) %>%
    summarise(
      Best_Score = max(Confidence_Score, na.rm = TRUE),
      Best_Confidence = Biological_Confidence[which.max(Confidence_Score)],
      .groups = "drop"
    )
}

build_aggregated_plot <- function(fragpipe_df, gene, uniprot_id, protein_length, domains_df = NULL, y_mode = "count") {
  gene_df     <- fragpipe_df %>% filter(Gene == gene)
  pos_df      <- .prepare_pos_df(gene_df)
  coverage_df <- .prepare_coverage_data(gene_df)

  if (nrow(pos_df) == 0 && nrow(coverage_df) == 0) {
    return(list(plot = NULL, error = paste0("No modification data for gene '", gene, "'.")))
  }

  has_experiment <- "ExperimentID" %in% colnames(gene_df)

  agg <- NULL
  if (nrow(pos_df) > 0) {
    if (has_experiment) {
      agg <- pos_df %>%
        group_by(Position, Modification, Mod_Type_Class) %>%
        summarise(
          Best_Score = max(Confidence_Score, na.rm = TRUE),
          N_Experiments = n_distinct(ExperimentID),
          Experiments = paste(unique(ExperimentID), collapse = ", "),
          Best_Confidence = Biological_Confidence[which.max(Confidence_Score)],
          Residue = dplyr::first(Residue),
          .groups = "drop"
        )
    } else {
      agg <- pos_df %>%
        group_by(Position, Modification, Mod_Type_Class) %>%
        summarise(
          Best_Score = max(Confidence_Score, na.rm = TRUE),
          N_Experiments = 1L,
          Experiments = "single",
          Best_Confidence = Biological_Confidence[which.max(Confidence_Score)],
          Residue = dplyr::first(Residue),
          .groups = "drop"
        )
    }

    # Y value depends on mode
    if (y_mode == "count") {
      agg <- agg %>% mutate(Y_Val = N_Experiments)
    } else {
      agg <- agg %>% mutate(Y_Val = Best_Score)
    }

    agg <- agg %>%
      mutate(
        HoverText = paste0(
          "<b>Position: ", Residue, Position, "</b><br>",
          "Modification: ", Modification, "<br>",
          "Best Confidence: ", Best_Confidence, " (", Best_Score, "/100)<br>",
          "Experiments: ", N_Experiments
        )
      )
  }

  top_mods <- if (!is.null(agg)) unique(agg$Modification) else character(0)
  cov_mods <- if (nrow(coverage_df) > 0) unique(coverage_df$Modification) else character(0)
  all_mods <- unique(c(top_mods, cov_mods))
  color_map <- .get_mod_colors(all_mods)

  all_class_map <- character(length(all_mods))
  for (i in seq_along(all_mods)) {
    mod <- all_mods[i]
    cls <- NULL
    if (!is.null(agg)) {
      rows <- agg %>% filter(Modification == mod)
      if (nrow(rows) > 0) cls <- rows$Mod_Type_Class[1]
    }
    if (is.null(cls) && nrow(coverage_df) > 0) {
      rows <- coverage_df %>% filter(Modification == mod)
      if (nrow(rows) > 0) cls <- rows$Mod_Type_Class[1]
    }
    all_class_map[i] <- ifelse(is.null(cls), "context", cls)
  }

  ordered_cov_mods <- all_mods[all_mods %in% cov_mods]
  n_cov_mods <- length(ordered_cov_mods)

  # Threshold lines only in confidence mode
  show_thresh <- (y_mode == "confidence")
  p_top <- plot_ly()
  p_top <- .add_backbone(p_top, protein_length, show_thresholds = show_thresh)
  n_fixed_traces <- if (show_thresh) 4L else 1L

  # Helper: initial visibility based on class (default = biological only)
  .init_vis <- function(mod) {
    idx <- which(all_mods == mod)
    if (length(idx) > 0 && all_class_map[idx] == "biological") TRUE else "legendonly"
  }

  if (!is.null(agg)) {
    for (mod in top_mods) {
      mod_data <- agg %>% filter(Modification == mod)
      p_top <- p_top %>% add_trace(
        data = mod_data, type = "bar",
        x = ~Position, y = ~Y_Val, name = mod,
        customdata = ~HoverText,
        hovertemplate = "%{customdata}<extra></extra>",
        marker = list(
          color = color_map[mod],
          opacity = 0.85,
          line = list(width = 0.5, color = "rgba(0,0,0,0.3)")
        ),
        legendgroup = mod,
        visible = .init_vis(mod),
        width = 1
      )
    }
  }

  cov_only_mods <- setdiff(cov_mods, top_mods)
  for (mod in cov_only_mods) {
    p_top <- p_top %>% add_trace(
      type = "scatter", mode = "markers",
      x = c(NA_real_), y = c(NA_real_), name = mod,
      marker = list(color = color_map[mod], size = 10, symbol = "square"),
      legendgroup = mod,
      showlegend = TRUE,
      visible = .init_vis(mod)
    )
  }

  # Y-axis config depends on mode  -" use biological-only range since default dropdown is "Likely Biological"
  if (y_mode == "count") {
    bio_agg <- if (!is.null(agg)) agg %>% filter(Mod_Type_Class == "biological") else NULL
    max_y <- if (!is.null(bio_agg) && nrow(bio_agg) > 0) max(bio_agg$Y_Val, na.rm = TRUE) + 2 else 5
    y_axis <- list(title = "Number of Experiments",
                   showticklabels = TRUE, showgrid = TRUE,
                   range = c(0, max_y), dtick = 1)
  } else {
    y_axis <- list(title = "Biological Confidence Score",
                   showticklabels = TRUE, showgrid = TRUE,
                   range = c(0, 105))
  }

  p_top <- p_top %>% layout(
    yaxis = y_axis,
    xaxis = list(showticklabels = FALSE, showgrid = TRUE,
                 range = c(0, protein_length + 1)),
    barmode = "overlay", hovermode = "closest"
  )

  has_coverage <- nrow(coverage_df) > 0
  has_domains <- !is.null(domains_df) && nrow(domains_df) > 0

  # --- DOMAIN TRACK (middle subplot) ---
  if (has_domains) {
    p_domain <- plot_ly()
    domain_names <- unique(domains_df$description)
    dom_colors <- .get_domain_colors(domain_names)
    dom_annotations <- list()
    min_label_width <- protein_length * 0.08

    for (d in domain_names) {
      d_rows <- domains_df %>% filter(description == d)
      for (i in seq_len(nrow(d_rows))) {
        row <- d_rows[i, ]
        hover <- paste0("<b>", d, "</b><br>Range: ", row$begin, "\u2013", row$end)
        p_domain <- p_domain %>% add_trace(
          type = "scatter", mode = "lines",
          x = c(row$begin, row$end), y = c(1, 1),
          line = list(width = 24, color = dom_colors[d]),
          name = d, legendgroup = paste0("dom_", d),
          showlegend = FALSE,
          hovertext = c(hover, hover), hoverinfo = "text"
        )
        dom_width <- row$end - row$begin
        mid_x <- (row$begin + row$end) / 2
        if (dom_width >= min_label_width) {
          dom_annotations <- c(dom_annotations, list(list(
            x = mid_x, y = 1, text = d,
            showarrow = FALSE, xanchor = "center", yanchor = "middle",
            font = list(size = 10, color = "black")
          )))
        } else {
          dom_annotations <- c(dom_annotations, list(list(
            x = mid_x, y = 1, ax = -40, ay = -25,
            text = d, showarrow = TRUE,
            arrowhead = 2, arrowsize = 0.8, arrowwidth = 1,
            xanchor = "right", yanchor = "bottom",
            font = list(size = 9, color = "#333")
          )))
        }
      }
    }
    p_domain <- p_domain %>% layout(
      yaxis = list(showticklabels = FALSE, showgrid = FALSE,
                   range = c(0.4, 1.6), fixedrange = TRUE, title = ""),
      xaxis = list(showticklabels = FALSE, showgrid = FALSE,
                   range = c(0, protein_length + 1)),
      annotations = dom_annotations
    )
  }

  # --- COVERAGE TRACK (bottom subplot) ---
  if (has_coverage) {
    p_bottom <- plot_ly()
    mod_y_map <- setNames(seq_len(n_cov_mods), ordered_cov_mods)

    for (mod in ordered_cov_mods) {
      mod_data <- coverage_df %>% filter(Modification == mod)
      y_level <- mod_y_map[mod]

      x_vals <- numeric(0)
      y_vals <- numeric(0)
      hover_vals <- character(0)

      for (i in seq_len(nrow(mod_data))) {
        row <- mod_data[i, ]
        hover <- paste0(
          "<b>", mod, "</b><br>",
          "Peptide: ", row$Peptide, "<br>",
          "Range: ", row$Protein_Start, "\u2013", row$Protein_End, "<br>",
          "Confidence: ", row$Best_Confidence, " (", row$Best_Score, "/100)"
        )
        x_vals <- c(x_vals, row$Protein_Start, row$Protein_End, NA)
        y_vals <- c(y_vals, y_level, y_level, NA)
        hover_vals <- c(hover_vals, hover, hover, "")
      }

      p_bottom <- p_bottom %>% add_trace(
        type = "scatter", mode = "lines",
        x = x_vals, y = y_vals,
        line = list(width = 12, color = color_map[mod]),
        name = mod, legendgroup = mod,
        showlegend = FALSE,
        visible = .init_vis(mod),
        hovertext = hover_vals, hoverinfo = "text"
      )
    }

    # Initial coverage labels/range match default "Likely Biological" filter
    init_cov_keep <- sapply(ordered_cov_mods, function(mod) {
      cls <- unique(coverage_df$Mod_Type_Class[coverage_df$Modification == mod])
      any(cls == "biological")
    })
    init_cov_labels <- ifelse(init_cov_keep, ordered_cov_mods, "")
    init_cov_visible_y <- which(init_cov_keep)
    if (length(init_cov_visible_y) > 0) {
      init_cov_lo <- min(init_cov_visible_y) - 0.7
      init_cov_hi <- max(init_cov_visible_y) + 0.7
    } else {
      init_cov_lo <- 0.3
      init_cov_hi <- 1.7
    }

    p_bottom <- p_bottom %>% layout(
      yaxis = list(
        tickvals = seq_len(n_cov_mods),
        ticktext = init_cov_labels,
        title = "",
        showgrid = FALSE,
        range = c(init_cov_lo, init_cov_hi),
        tickfont = list(size = 9)
      ),
      xaxis = list(title = "Protein position",
                    range = c(0, protein_length + 1), showgrid = TRUE)
    )
  }

  # --- ASSEMBLE SUBPLOTS ---
  sub_plots <- list(p_top)
  sub_heights <- c(0.55)
  if (has_domains) {
    sub_plots <- c(sub_plots, list(p_domain))
    sub_heights <- c(sub_heights, 0.10)
  }
  if (has_coverage) {
    # Use biological-only mod count for initial coverage fraction (matches default dropdown)
    n_init_bio_cov <- sum(init_cov_keep)
    cov_h <- min(0.35, max(0.08, n_init_bio_cov * 0.04))
    sub_plots <- c(sub_plots, list(p_bottom))
    sub_heights <- c(sub_heights, cov_h)
  }
  sub_heights <- sub_heights / sum(sub_heights)

  if (length(sub_plots) > 1) {
    p <- subplot(sub_plots, nrows = length(sub_plots), shareX = TRUE,
                 heights = sub_heights)
  } else {
    p <- p_top
  }

  n_top_trace <- length(top_mods) + length(cov_only_mods)
  n_cov_trace <- if (has_coverage) n_cov_mods else 0L
  n_dom_trace <- if (has_domains) length(unique(domains_df$description)) else 0L

  top_class_map <- character(n_top_trace)
  for (i in seq_along(top_mods)) {
    top_class_map[i] <- all_class_map[which(all_mods == top_mods[i])]
  }
  for (i in seq_along(cov_only_mods)) {
    top_class_map[length(top_mods) + i] <- all_class_map[which(all_mods == cov_only_mods[i])]
  }

  cov_class_map <- character(n_cov_trace)
  if (has_coverage) {
    for (i in seq_along(ordered_cov_mods)) {
      cov_class_map[i] <- all_class_map[which(all_mods == ordered_cov_mods[i])]
    }
  }

  fixed_vis <- rep(TRUE, n_fixed_traces)
  dom_vis <- rep(TRUE, n_dom_trace)
  make_vis <- function(classes) {
    top_vis <- ifelse(top_class_map %in% classes, TRUE, "legendonly")
    dom_part <- as.list(dom_vis)
    cov_vis <- if (n_cov_trace > 0) {
      ifelse(cov_class_map %in% classes, TRUE, "legendonly")
    } else {
      logical(0)
    }
    c(as.list(fixed_vis), as.list(top_vis), dom_part, as.list(cov_vis))
  }

  # Determine which yaxis the coverage track is on for tick label updates
  cov_yaxis_key <- NULL
  if (has_coverage) {
    subplot_idx <- 1L  # p_top is always first
    if (has_domains) subplot_idx <- subplot_idx + 1L  # p_domain second
    subplot_idx <- subplot_idx + 1L  # p_bottom
    cov_yaxis_key <- if (subplot_idx == 1) "yaxis" else paste0("yaxis", subplot_idx)
  }

  # Pre-compute per-class max Y for bar chart auto-ranging
  .class_max_y <- function(classes) {
    if (is.null(agg) || nrow(agg) == 0) return(5)
    filtered <- agg %>% filter(Mod_Type_Class %in% classes)
    if (nrow(filtered) == 0) return(5)
    max(filtered$Y_Val, na.rm = TRUE) + 2
  }

  # Store the initial subplot domains for dynamic resizing
  # plotly subplot uses: bottom-up ordering with a 2% gap
  # yaxis2 (coverage) at bottom, yaxis (bar chart) at top
  gap <- 0.04  # standard plotly subplot gap
  initial_bar_domain <- NULL
  initial_cov_domain <- NULL
  initial_dom_domain <- NULL
  if (has_coverage) {
    norm_h <- sub_heights  # already normalized
    if (has_domains) {
      # 3 subplots: bar (top), domain (middle), coverage (bottom)
      cov_top <- norm_h[3]
      dom_bot <- cov_top + gap
      dom_top <- dom_bot + norm_h[2]
      bar_bot <- dom_top + gap
      bar_top <- 1.0
      initial_cov_domain <- c(0, cov_top)
      initial_dom_domain <- c(dom_bot, dom_top)
      initial_bar_domain <- c(bar_bot, bar_top)
    } else {
      # 2 subplots: bar (top), coverage (bottom)
      cov_top <- norm_h[2]
      bar_bot <- cov_top + gap
      bar_top <- 1.0
      initial_cov_domain <- c(0, cov_top)
      initial_bar_domain <- c(bar_bot, bar_top)
    }
  }

  # Build per-class layout updates: coverage tick labels + range, bar Y range, and domain resizing
  make_layout_update <- function(classes) {
    updates <- list()

    # Auto-range the bar chart Y axis
    class_max <- .class_max_y(classes)
    if (y_mode == "count") {
      updates[["yaxis"]] <- list(
        title = "Number of Experiments",
        showticklabels = TRUE, showgrid = TRUE,
        range = c(0, class_max), dtick = 1
      )
    } else {
      updates[["yaxis"]] <- list(
        title = "Biological Confidence Score",
        showticklabels = TRUE, showgrid = TRUE,
        range = c(0, 105)
      )
    }

    # Update coverage track axis: filter labels + compress range + resize domains
    if (has_coverage && !is.null(cov_yaxis_key)) {
      keep <- cov_class_map %in% classes
      filtered_labels <- ifelse(keep, ordered_cov_mods, "")
      n_visible <- sum(keep)
      visible_y <- which(keep)

      if (n_visible > 0) {
        y_lo <- min(visible_y) - 0.7
        y_hi <- max(visible_y) + 0.7
      } else {
        y_lo <- 0.3
        y_hi <- 1.7
        n_visible <- 1
      }

      # Dynamically resize subplot proportions: give coverage less space
      # when fewer mods are visible
      new_cov_frac <- min(0.35, max(0.08, n_visible * 0.04))
      dom_frac <- if (has_domains) 0.08 else 0
      bar_frac <- 1.0 - new_cov_frac - dom_frac - (if (has_domains) 2 * gap else gap)

      new_cov_domain <- c(0, new_cov_frac)
      if (has_domains) {
        new_dom_bot <- new_cov_frac + gap
        new_dom_top <- new_dom_bot + dom_frac
        new_bar_bot <- new_dom_top + gap
        # Update domain axis too
        dom_yaxis_key <- "yaxis2"
        cov_actual_key <- "yaxis3"
        updates[[dom_yaxis_key]] <- list(
          domain = c(new_dom_bot, new_dom_top),
          showticklabels = FALSE, showgrid = FALSE,
          range = c(0.4, 1.6), fixedrange = TRUE, title = ""
        )
      } else {
        new_bar_bot <- new_cov_frac + gap
        cov_actual_key <- cov_yaxis_key
      }

      updates[["yaxis"]]$domain <- c(new_bar_bot, 1.0)

      updates[[cov_actual_key]] <- list(
        tickvals = seq_len(n_cov_mods),
        ticktext = as.list(filtered_labels),
        title = "", showgrid = FALSE,
        range = c(y_lo, y_hi),
        domain = new_cov_domain,
        tickfont = list(size = 9)
      )
    }

    updates
  }

  # Title and annotations depend on mode
  if (y_mode == "count") {
    plot_title <- paste0(gene, " (", uniprot_id, ") \u2013 Experiment Frequency")
    base_annots <- list()
  } else {
    plot_title <- paste0(gene, " (", uniprot_id, ") \u2013 Confidence Map")
    base_annots <- .make_annotations(protein_length)
  }

  p <- p %>% layout(
    title = list(text = plot_title, font = list(size = 16)),
    legend = list(title = list(text = "Modification<br>(Click to filter)"),
                  x = 1.02, y = 1, xanchor = "left"),
    annotations = base_annots,
    margin = list(l = 150, b = 60, t = 80),
    hovermode = "closest",
    barmode = "overlay",
    updatemenus = list(
      list(
        type = "dropdown", direction = "down",
        xanchor = "left", yanchor = "top", x = 0, y = 1.15,
        active = 1,
        buttons = list(
          list(label = "Show All",
               method = "update",
               args = list(
                 list(visible = make_vis(c("biological", "context", "artifact"))),
                 c(list(annotations = base_annots),
                   make_layout_update(c("biological", "context", "artifact")))
               )),
          list(label = "Likely Biological",
               method = "update",
               args = list(
                 list(visible = make_vis("biological")),
                 c(list(annotations = base_annots),
                   make_layout_update("biological"))
               )),
          list(label = "Unlikely Biological",
               method = "update",
               args = list(
                 list(visible = make_vis("context")),
                 c(list(annotations = base_annots),
                   make_layout_update("context"))
               )),
          list(label = "Likely Artifact",
               method = "update",
               args = list(
                 list(visible = make_vis("artifact")),
                 c(list(annotations = base_annots),
                   make_layout_update("artifact"))
               ))
        )
      )
    )
  )

  return(list(plot = p, error = NULL))
}

build_modification_plot <- function(fragpipe_df, gene, uniprot_id,
                                     y_mode = "count",
                                     domains_df = NULL) {
  protein_length <- get_uniprot_length(uniprot_id)
  if (is.na(protein_length)) {
    return(list(plot = NULL,
                error = paste0("Unable to fetch UniProt sequence for ", uniprot_id,
                               ". Check the ID and your internet connection.")))
  }

  gene_df <- fragpipe_df %>% filter(Gene == gene)
  if (nrow(gene_df) == 0) {
    return(list(plot = NULL,
                error = paste0("No entries found for gene '", gene, "' in the database.")))
  }

  build_aggregated_plot(fragpipe_df, gene, uniprot_id, protein_length, domains_df, y_mode)
}

# === SHARED HEATMAP HELPERS ===================================================
# Used by both server_heatmap.R and server_tmt_heatmap.R
# ==============================================================================

#' Build a plotly heatmap from a numeric matrix.
#' Row z-score normalization + hierarchical clustering on rows.
#' Column ordering: hierarchical clustering (default) or group order.
#'
#' @param mat Numeric matrix (rows=genes, cols=samples), already subsetted
#' @param title_txt Plot title
#' @param mode_label Mode description for metadata
#' @param col_order_mode "cluster" (default) or "group"
#' @param col_groups Named character vector (column_name -> group_name), optional
#' @param log_fn Logger function (single string arg), defaults to message()
#' @return list(plot, height, n_genes, n_samples, mode)
build_heatmap_plotly <- function(mat, title_txt = "Heatmap", mode_label = "",
                                 col_order_mode = "cluster", col_groups = NULL,
                                 scale_mode = "zscore", na_display = "zero",
                                 log_fn = message) {
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

  is_zscore_display <- scale_mode %in% c("zscore", "zscore_pre")

  if (identical(scale_mode, "zscore")) {
    # Row z-score normalization (with zero-variance purge)
    row_means <- rowMeans(mat, na.rm = TRUE)
    row_sds   <- apply(mat, 1, sd, na.rm = TRUE)

    # Variance Gate: purge zero-variance rows before dividing
    row_var <- row_sds^2
    row_var[is.na(row_var)] <- 0
    keep <- row_var > 0
    mat      <- mat[keep, , drop = FALSE]
    row_means <- row_means[keep]
    row_sds   <- row_sds[keep]

    if (nrow(mat) < 2) return(NULL)

    z_mat <- (mat - row_means) / row_sds

    # Impute residual NAs to 0 (z-score mean = neutral) instead of dropping rows.
    # With na_display = "gap" they stay missing so plotly renders them as holes
    # instead of at the exact centre of the diverging ramp. Non-finite values
    # are normalised to NA in that mode so the clustering zero-fill below still
    # sees the same cells it always did.
    if (identical(na_display, "zero")) {
      z_mat[is.na(z_mat) | !is.finite(z_mat)] <- 0
    } else {
      z_mat[!is.finite(z_mat)] <- NA_real_
    }
  } else if (identical(scale_mode, "zscore_pre")) {
    # Pre-computed z-score: 0 is the row mean, so a neutral fill is valid here.
    z_mat <- mat
    if (identical(na_display, "zero")) z_mat[is.na(z_mat)] <- 0
  } else {
    # Normalized abundance: leave NAs as NAs. Filling with 0 on a log2 scale
    # plants a false floor far below any real value and rescales the colour
    # ramp. Plotly renders NA cells as gaps.
    z_mat <- mat
  }

  if (nrow(z_mat) < 2) return(NULL)

  # Whether clustering actually ran. A row count is not proof of it: the
  # tryCatch blocks below swallow errors and leave the matrix unordered, and the
  # caller had no way to tell a successful ordering from a caught failure. Set
  # TRUE on the success path only, and returned so provenance can report it.
  row_clustering_applied <- FALSE
  col_clustering_applied <- FALSE

  # Hierarchical clustering on rows. The zero-fill applies ONLY in z-score
  # display modes, where 0 is the row mean and therefore a neutral stand-in. In
  # abundance mode the values are log2 intensities around 15-25, so a zero-fill
  # plants a false floor about 20 units below any real value and dominates the
  # distance - the same argument the display path above already makes. Left as
  # NA, stats::dist excludes NA pairs and rescales by the number of
  # contributing coordinates, which is the principled behaviour.
  z_for_dist <- z_mat
  if (is_zscore_display) z_for_dist[is.na(z_for_dist)] <- 0
  if (nrow(z_for_dist) > 2) {
    tryCatch({
      d <- dist(z_for_dist, method = "euclidean")
      hc <- hclust(d, method = "ward.D2")
      z_mat <- z_mat[hc$order, , drop = FALSE]
      row_clustering_applied <- TRUE
    }, error = function(e) {
      log_fn(paste0("Heatmap: Row clustering failed: ", e$message))
    })
  }

  # Column ordering
  if (col_order_mode == "manual") {
    # Keep columns in their current order (user-specified)
  } else if (col_order_mode == "group" && !is.null(col_groups) && length(col_groups) > 0) {
    mat_cols <- colnames(z_mat)
    matched_groups <- col_groups[mat_cols]
    matched_groups[is.na(matched_groups)] <- "Ungrouped"
    group_levels <- unique(col_groups)
    if ("Ungrouped" %in% matched_groups && !"Ungrouped" %in% group_levels) {
      group_levels <- c(group_levels, "Ungrouped")
    }
    col_order <- order(match(matched_groups, group_levels))
    z_mat <- z_mat[, col_order, drop = FALSE]
  } else {
    z_col <- t(z_mat)
    if (is_zscore_display) z_col[is.na(z_col)] <- 0   # see the row block above
    if (nrow(z_col) > 2) {
      tryCatch({
        d_col <- dist(z_col, method = "euclidean")
        hc_col <- hclust(d_col, method = "ward.D2")
        z_mat <- z_mat[, hc_col$order, drop = FALSE]
        col_clustering_applied <- TRUE
      }, error = function(e) {
        log_fn(paste0("Heatmap: Column clustering failed: ", e$message))
      })
    }
  }

  # Hover text
  hover_mat <- matrix("", nrow = nrow(z_mat), ncol = ncol(z_mat))
  for (r in seq_len(nrow(z_mat))) {
    for (cc in seq_len(ncol(z_mat))) {
      val_text  <- if (is.na(z_mat[r, cc])) "NA" else round(z_mat[r, cc], 3)
      val_label <- if (is_zscore_display) "Z-score: " else "log2 Abundance: "
      hover_mat[r, cc] <- paste0(
        "<b>", rownames(z_mat)[r], "</b><br>",
        "Sample: ", colnames(z_mat)[cc], "<br>",
        val_label, val_text
      )
    }
  }

  # Color scale
  finite_vals <- z_mat[is.finite(z_mat)]
  if (length(finite_vals) == 0) return(NULL)
  max_abs <- max(abs(finite_vals), na.rm = TRUE)

  if (is_zscore_display) {
    colorscale <- list(
      list(0,    "#2166ac"),
      list(0.25, "#67a9cf"),
      list(0.5,  "#f7f7f7"),
      list(0.75, "#ef8a62"),
      list(1,    "#b2182b")
    )
    colorbar_title <- "Z-score"
    z_min <- -max_abs
    z_max <-  max_abs
  } else {
    # Viridis-inspired sequential scale for normalized abundance
    colorscale <- list(
      list(0,    "#440154"),
      list(0.25, "#3b528b"),
      list(0.5,  "#21918c"),
      list(0.75, "#5ec962"),
      list(1,    "#fde725")
    )
    colorbar_title <- "log2 Abund."
    z_min <- min(finite_vals, na.rm = TRUE)
    z_max <- max(finite_vals, na.rm = TRUE)
  }

  calc_height <- max(500, min(8000, 80 + nrow(z_mat) * 16))

  p <- plot_ly(
    z = z_mat,
    x = colnames(z_mat),
    y = rownames(z_mat),
    type = "heatmap",
    colorscale = colorscale,
    zmin = z_min, zmax = z_max,
    hoverinfo = "text", text = hover_mat,
    showscale = TRUE,
    colorbar = list(title = colorbar_title, len = 0.5)
  ) %>% layout(
    title = list(text = title_txt, font = list(size = 15)),
    xaxis = list(title = "", tickangle = -45, tickfont = list(size = 10)),
    yaxis = list(title = "", tickfont = list(size = 9), autorange = "reversed"),
    margin = list(l = 160, b = 130, t = 60)
  )

  list(plot = p, height = calc_height,
       n_genes = nrow(z_mat), n_samples = ncol(z_mat),
       mode = mode_label,
       row_clustering_applied = row_clustering_applied,
       col_clustering_applied = col_clustering_applied)
}

#' Collapse replicate columns by group mean.
#' @param mat Numeric matrix (rows=genes, cols=samples)
#' @param col_groups Named character vector (column_name -> group_name)
#' @return Averaged matrix with one column per group, or NULL if < 2 groups
heatmap_group_average <- function(mat, col_groups) {
  mat_cols <- colnames(mat)
  matched_groups <- col_groups[mat_cols]
  matched_groups[is.na(matched_groups)] <- "Ungrouped"

  group_levels <- unique(col_groups)
  if ("Ungrouped" %in% matched_groups && !"Ungrouped" %in% group_levels) {
    group_levels <- c(group_levels, "Ungrouped")
  }
  unique_groups <- group_levels[group_levels %in% matched_groups]

  if (length(unique_groups) < 2) return(NULL)

  avg_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = length(unique_groups))
  rownames(avg_mat) <- rownames(mat)
  colnames(avg_mat) <- unique_groups
  for (g in seq_along(unique_groups)) {
    grp_cols <- mat_cols[matched_groups == unique_groups[g]]
    if (length(grp_cols) == 1) {
      avg_mat[, g] <- mat[, grp_cols]
    } else if (length(grp_cols) > 1) {
      avg_mat[, g] <- rowMeans(mat[, grp_cols, drop = FALSE], na.rm = TRUE)
    }
  }
  avg_mat
}


# === SHARED TREND CLUSTER HELPERS =============================================
# Used by both server_trends.R and server_trend_clusters.R
# ==============================================================================

#' Row-wise Z-score normalization of a matrix.
#' Zero-variance rows are purged before scaling to prevent NaN from sd()==0.
#' A final complete.cases() filter catches any rogue NaN survivors.
#' @param mat Numeric matrix (rows=proteins, cols=conditions), log2 scale
#' @param col_labels Character vector of condition labels for column names
#' @return List with z (Z-scored matrix) and kept_idx (integer indices of rows retained)
rowwise_zscore <- function(mat, col_labels) {
  n_input <- nrow(mat)
  n_na    <- sum(is.na(mat))

  # Variance Gate: identify zero-variance rows
  row_var <- apply(mat, 1, var, na.rm = TRUE)
  row_var[is.na(row_var)] <- 0
  keep <- row_var > 0
  n_zero_var <- sum(!keep)
  mat <- mat[keep, , drop = FALSE]

  # Execute Scaling (NA-safe: proteins with partial NAs keep their valid values)
  z <- t(apply(mat, 1, function(row) {
    m <- mean(row, na.rm = TRUE); s <- sd(row, na.rm = TRUE)
    (row - m) / s
  }))
  colnames(z) <- col_labels

  # Keep rows that have at least 2 non-NA z-scores
  has_data <- rowSums(!is.na(z)) >= 2
  n_sparse <- sum(!has_data)
  z <- z[has_data, , drop = FALSE]

  message(paste0("[rowwise_zscore] ", n_input, " input rows, ",
                 n_na, " NAs in matrix, ",
                 n_zero_var, " zero-var dropped, ",
                 n_sparse, " sparse (<2 values) dropped, ",
                 nrow(z), " kept"))

  # Return both the matrix and the logical mask for callers to subset metadata
  kept_idx <- which(keep)[has_data]
  list(z = z, kept_idx = kept_idx)
}

#' Build long-form data.frame from Z-score matrix + cluster assignments,
#' compute cluster means, sort clusters by last-condition direction.
#'
#' @param zscore_mat Z-scored numeric matrix (rows=proteins, cols=conditions)
#' @param accessions Character vector of accession IDs
#' @param gene_names Character vector of gene names
#' @param col_labels Character vector of condition labels
#' @param clusters Integer vector of cluster assignments (length = nrow)
#' @param memberships Matrix of fuzzy c-means membership scores (nrow x k), or NULL
#' @return list(long_df, cluster_means, cluster_info, cluster_map)
build_cluster_long_df <- function(zscore_mat, accessions, gene_names,
                                   col_labels, clusters, memberships = NULL) {
  n_proteins <- nrow(zscore_mat)
  n_conds    <- ncol(zscore_mat)

  long_df <- data.frame(
    Accession = rep(accessions, each = n_conds),
    GeneName  = rep(gene_names, each = n_conds),
    Condition = rep(col_labels, times = n_proteins),
    CondIdx   = rep(seq_len(n_conds), times = n_proteins),
    ZScore    = as.vector(t(zscore_mat)),
    Cluster   = rep(clusters, each = n_conds),
    stringsAsFactors = FALSE
  )

  # Add membership scores (fuzzy c-means) or default 1.0

  if (!is.null(memberships)) {
    max_membership <- apply(memberships, 1, max)
    long_df$Membership <- rep(max_membership, each = n_conds)
  } else {
    long_df$Membership <- 1.0
  }

  # Cluster means
  cluster_means <- aggregate(ZScore ~ Cluster + Condition + CondIdx,
                              data = long_df, FUN = mean)
  colnames(cluster_means)[colnames(cluster_means) == "ZScore"] <- "MeanZ"

  # Cluster info
  cluster_info <- as.data.frame(table(Cluster = clusters))
  colnames(cluster_info) <- c("Cluster", "n_proteins")
  cluster_info$Cluster <- as.integer(as.character(cluster_info$Cluster))

  # Sort clusters by mean direction at last condition (descending).
  # Cluster numbering therefore depends on the user-chosen condition order and
  # is DISPLAY ORDER ONLY. Cluster IDs are not stable between runs, nor between
  # two orderings of the same conditions, and must not be cited as identifiers.
  last_cond_means <- cluster_means[cluster_means$CondIdx == n_conds, ]
  cluster_order <- last_cond_means$Cluster[order(-last_cond_means$MeanZ)]
  cluster_map <- setNames(seq_along(cluster_order), cluster_order)
  long_df$Cluster       <- cluster_map[as.character(long_df$Cluster)]
  cluster_means$Cluster <- cluster_map[as.character(cluster_means$Cluster)]
  cluster_info$Cluster  <- cluster_map[as.character(cluster_info$Cluster)]

  # Remap membership matrix columns if fuzzy c-means
  if (!is.null(memberships)) {
    orig_ids <- as.integer(names(cluster_map))
    orig_ids <- orig_ids[orig_ids <= ncol(memberships)]
    memberships <- memberships[, orig_ids, drop = FALSE]
    colnames(memberships) <- paste0("Cluster_", seq_len(ncol(memberships)))
  }

  # Labeled factors
  cluster_labels <- paste0(
    "Cluster ", cluster_info$Cluster[order(cluster_info$Cluster)],
    " (n=", cluster_info$n_proteins[order(cluster_info$Cluster)], ")")

  long_df$ClusterLabel <- factor(
    paste0("Cluster ", long_df$Cluster,
           " (n=", cluster_info$n_proteins[match(long_df$Cluster, cluster_info$Cluster)], ")"),
    levels = cluster_labels)
  cluster_means$ClusterLabel <- factor(
    paste0("Cluster ", cluster_means$Cluster,
           " (n=", cluster_info$n_proteins[match(cluster_means$Cluster, cluster_info$Cluster)], ")"),
    levels = cluster_labels)

  # Condition as ordered factor
  long_df$Condition       <- factor(long_df$Condition, levels = col_labels)
  cluster_means$Condition <- factor(cluster_means$Condition, levels = col_labels)

  list(long_df = long_df, cluster_means = cluster_means,
       cluster_info = cluster_info, cluster_map = cluster_map,
       memberships = memberships)
}

#' Render plotly trend profile subplots.
#'
#' @param long_df From build_cluster_long_df()
#' @param cluster_means From build_cluster_long_df()
#' @param cluster_info From build_cluster_long_df()
#' @param col_labels Condition label vector
#' @param k Number of clusters
#' @param highlight Gene name to highlight (NULL or "")
#' @param exact_match Logical, exact vs regex match for highlight
#' @param title Plot title string
#' @param method "fcm" or "mfuzz" (controls membership display)
#' @param mem_thresh Membership threshold for filtering (Mfuzz only, 0 = off)
#' @return plotly figure
render_trend_subplots <- function(long_df, cluster_means, cluster_info,
                                   col_labels, k,
                                   highlight = NULL, exact_match = FALSE,
                                   title = NULL, method = "fcm",
                                   mem_thresh = 0) {
  clusters_sorted <- sort(unique(long_df$Cluster))

  plot_list <- lapply(clusters_sorted, function(cl) {
    cl_data  <- long_df[long_df$Cluster == cl, ]
    cl_means <- cluster_means[cluster_means$Cluster == cl, ]
    cl_means <- cl_means[order(cl_means$CondIdx), ]
    n_prot   <- cluster_info$n_proteins[cluster_info$Cluster == cl]

    # Membership threshold filter (fuzzy c-means only)
    if (method %in% c("fcm", "mfuzz") && mem_thresh > 0) {
      cl_data <- cl_data[cl_data$Membership >= mem_thresh, ]
    }
    n_displayed <- length(unique(cl_data$Accession))

    p <- plot_ly()

    # Consolidated protein lines (single trace with NA breaks)
    unique_acc <- unique(cl_data$Accession)
    n_conds    <- length(col_labels)
    all_x    <- vector("numeric",   length(unique_acc) * (n_conds + 1))
    all_y    <- vector("numeric",   length(unique_acc) * (n_conds + 1))
    all_text <- vector("character", length(unique_acc) * (n_conds + 1))
    idx <- 1L
    for (acc in unique_acc) {
      acc_df <- cl_data[cl_data$Accession == acc, ]
      acc_df <- acc_df[order(acc_df$CondIdx), ]
      gene   <- acc_df$GeneName[1]
      n_r    <- nrow(acc_df)
      rng    <- idx:(idx + n_r - 1L)
      all_x[rng] <- acc_df$CondIdx
      all_y[rng] <- acc_df$ZScore
      all_text[rng] <- paste0(
        "<b>Gene:</b> ", gene,
        "<br><b>Accession:</b> ", acc,
        "<br><b>Condition:</b> ", acc_df$Condition,
        "<br><b>Z-Score:</b> ", round(acc_df$ZScore, 3),
        if (method %in% c("fcm", "mfuzz") && "Membership" %in% colnames(acc_df))
          paste0("<br><b>Membership:</b> ", round(acc_df$Membership, 3))
        else ""
      )
      idx <- idx + n_r
      all_x[idx]    <- NA_real_
      all_y[idx]    <- NA_real_
      all_text[idx] <- ""
      idx <- idx + 1L
    }
    all_x    <- all_x[seq_len(idx - 1L)]
    all_y    <- all_y[seq_len(idx - 1L)]
    all_text <- all_text[seq_len(idx - 1L)]

    p <- add_trace(p,
      x = all_x, y = all_y,
      type = "scatter", mode = "lines",
      line = list(color = "rgba(170,170,170,0.15)", width = 0.8),
      hoverinfo = "text", text = all_text,
      showlegend = FALSE, connectgaps = FALSE
    )

    # Cluster average (thick red)
    mean_hover <- paste0(
      "<b>Cluster ", cl, " Average</b> (all ", n_prot, " members)",
      "<br><b>Condition:</b> ", cl_means$Condition,
      "<br><b>Mean Z-Score:</b> ", round(cl_means$MeanZ, 3)
    )
    p <- add_trace(p,
      x = cl_means$CondIdx, y = cl_means$MeanZ,
      type = "scatter", mode = "lines+markers",
      line = list(color = "#e74c3c", width = 3),
      marker = list(color = "#e74c3c", size = 7),
      hoverinfo = "text", text = mean_hover,
      name = paste0("Cluster ", cl, " mean"),
      showlegend = FALSE
    )

    # Highlight gene traces
    if (!is.null(highlight) && nchar(highlight) > 0) {
      if (exact_match) {
        hl_df <- cl_data[toupper(cl_data$GeneName) == toupper(highlight), ]
      } else {
        hl_df <- cl_data[grepl(highlight, cl_data$GeneName, ignore.case = TRUE), ]
      }
      if (nrow(hl_df) > 0) {
        for (hacc in unique(hl_df$Accession)) {
          h_df <- hl_df[hl_df$Accession == hacc, ]
          h_df <- h_df[order(h_df$CondIdx), ]
          hl_hover <- paste0(
            "<b>\u2b50 ", h_df$GeneName[1], "</b>",
            "<br><b>Accession:</b> ", hacc,
            "<br><b>Condition:</b> ", h_df$Condition,
            "<br><b>Z-Score:</b> ", round(h_df$ZScore, 3)
          )
          p <- add_trace(p,
            x = h_df$CondIdx, y = h_df$ZScore,
            type = "scatter", mode = "lines+markers",
            line = list(color = "#f39c12", width = 2.5),
            marker = list(color = "#f39c12", size = 6),
            hoverinfo = "text", text = hl_hover,
            showlegend = FALSE
          )
        }
      }
    }

    # Annotation for subplot title
    anno_text <- if (method %in% c("fcm", "mfuzz") && n_displayed < n_prot) {
      paste0("<b>Cluster ", cl, "</b> (", n_displayed, "/", n_prot, ")")
    } else {
      paste0("<b>Cluster ", cl, "</b> (n=", n_prot, ")")
    }

    p <- layout(p,
      xaxis = list(tickvals = seq_along(col_labels), ticktext = col_labels,
                    tickangle = 45, title = ""),
      yaxis = list(title = if (cl == clusters_sorted[1]) "Z-Score" else ""),
      annotations = list(list(
        text = anno_text, xref = "paper", yref = "paper",
        x = 0.5, y = 1.08, showarrow = FALSE,
        font = list(size = 12, color = "#2c3e50")
      ))
    )
    p
  })

  n_cols <- min(3, length(plot_list))
  n_rows <- ceiling(length(plot_list) / n_cols)

  fig <- subplot(plot_list, nrows = n_rows, shareX = TRUE, shareY = TRUE,
                  titleY = TRUE, margin = 0.06)

  if (is.null(title)) {
    title <- paste0("Protein Trend Profiles (k = ", k, ")")
  }

  fig <- layout(fig,
    title = list(text = title, font = list(size = 16, color = "#2c3e50")),
    showlegend = FALSE, hovermode = "closest",
    margin = list(t = 80, b = 80)
  )
  fig
}

# ==============================================================================
# SHARED TMT DESIGN GRID  -" used by Heatmap, PCA, Trends, Volcano
# Registers arch-fork, design-grid modal, and all grid event handlers.
#
# Arguments:
#   input, session         -" Shiny input/session
#   prefix                 -" Input ID prefix, e.g. "thm_", "tpca_", "ttrnd_", "tmt_"
#   store                  -" Module-local reactiveValues (must have: all_channels,
#                           n_plexes, has_bridge, design_conds, design_option,
#                           annotation_map, conditions, status)
#   rv                     -" Shared cross-module reactiveValues
#   logger                 -" Logger function
#   module_name            -" Display name for logs, e.g. "TMT Heatmap"
#   design_source_tag      -" Tag for rv$tmt_design_source, e.g. "tmt_heatmap"
#   on_design_finished     -" Optional callback(ann, NL, real_conds) for module-
#                           specific post-design actions (e.g. volcano comparisons)
#   cond_placeholder       -" Placeholder text for "Add Condition" input
#
# Returns: list with $show_arch_fork and $show_design_grid functions.
# ==============================================================================
# Annotation map from the design grid (one row per channel): the plex number and
# condition chosen for each channel. Used by the grid's Finish Design handler.
tmt_grid_annotation <- function(channels, grid_plex, grid_cond) {
  ann_rows <- list(); biorep_counter <- list()
  # Norm channels per plex: a plex with one keeps the label "MixN_Norm"; with
  # several they are numbered so that every channel's BioReplicate stays unique.
  n_norm_plex <- tapply(grid_cond == "Norm", grid_plex, sum)
  for (i in seq_along(channels)) {
    ch   <- channels[i]; tag <- tmt_extract_tag(ch); cond <- grid_cond[i]
    pid  <- grid_plex[i]; mix <- paste0("Mix", pid)
    key  <- paste0(mix, "_", cond)
    biorep_counter[[key]] <- (biorep_counter[[key]] %||% 0L) + 1L
    bnum <- biorep_counter[[key]]
    biorep_label <- if (cond != "Norm") paste0(mix, "_", cond, "_", bnum)
                    else if (n_norm_plex[[as.character(pid)]] > 1L) paste0(mix, "_Norm_", bnum)
                    else paste0(mix, "_Norm")
    ann_rows[[i]] <- data.frame(ColumnName = ch, Channel = tag, Condition = cond,
                                 BioReplicate = biorep_label, Mixture = mix,
                                 Run = paste0(mix, "_run"), stringsAsFactors = FALSE)
  }
  do.call(rbind, ann_rows)
}

# Methods-text sentence for the channels marked "Exclude" in the design grid
# (character(0) when none). Shared by the TMT protein and TMT PTM methods text.
tmt_excluded_sentence <- function(excluded) {
  n <- length(excluded)
  if (n == 0L) return(character(0))
  paste0(n, " TMT channel", if (n != 1L) "s" else "", " (",
         paste(vapply(excluded, tmt_extract_tag, character(1), USE.NAMES = FALSE),
               collapse = ", "), ") ",
         if (n != 1L) "were" else "was",
         " excluded from the analysis at the experimental-design step and ",
         "not used in normalization, summarization, imputation or statistical testing.")
}

tmt_register_design_grid <- function(input, session, prefix, store, rv, logger,
                                      module_name, design_source_tag,
                                      on_design_finished = NULL,
                                      cond_placeholder = "e.g. Treated") {

  gp <- function(suffix) paste0(prefix, suffix)

  # ── Arch Fork Modal ──────────────────────────────────────────────────────
  .show_arch_fork <- function() {
    showModal(modalDialog(
      title = tags$span(icon("sitemap"), " Experimental Design"),
      tagList(
        tags$div(class = "alert alert-info py-2 px-3",
          tags$p(class = "mb-1", tags$b("How many TMT plexes does your experiment have?")),
          tags$p(class = "mb-0 small text-muted",
            "Single plex = one TMT set. Multiple plexes = samples split across runs.")
        ),
        tags$hr(),
        tags$div(class = "d-flex gap-3 justify-content-center",
          actionButton(gp("arch_single"), tags$span(icon("1"), " Single Plex"),
                       class = "btn-primary btn-lg px-4"),
          actionButton(gp("arch_multi"), tags$span(icon("layer-group"), " Multiple Plexes"),
                       class = "btn-warning btn-lg px-4")
        )
      ),
      size = "m", easyClose = FALSE, footer = NULL
    ))
  }

  observeEvent(input[[gp("arch_single")]], {
    removeModal()
    store$n_plexes <- 1L; store$has_bridge <- FALSE; store$design_conds <- character(0)
    store$grid_assignments <- NULL; store$paint_cond <- NULL
    Sys.sleep(0.2)
    if (length(store$all_channels) == 0) {
      showNotification("No channels detected.", type = "error"); return()
    }
    .show_design_grid()
  })

  observeEvent(input[[gp("arch_multi")]], {
    removeModal()
    n_f <- max(2L, tmt_default_channel_plex(store$files, store$all_channels))
    store$n_plexes <- as.integer(n_f)
    store$has_bridge <- FALSE; store$design_conds <- character(0)
    store$grid_assignments <- NULL; store$paint_cond <- NULL
    Sys.sleep(0.2)
    if (length(store$all_channels) == 0) {
      showNotification("No channels detected.", type = "error"); return()
    }
    .show_design_grid()
  })

  # ── Design Grid Modal ───────────────────────────────────────────────────
  # `keep_search`: TRUE when the grid is rebuilt from inside itself (add
  # condition, bridge, plexes), so the channel search term survives; a fresh
  # open starts with an empty search.
  .show_design_grid <- function(keep_search = FALSE) {
    channels <- store$all_channels
    np       <- store$n_plexes
    conds    <- store$design_conds
    is_multi <- np > 1
    has_br   <- store$has_bridge
    search_term <- if (keep_search) input[[gp("grid_search")]] %||% "" else ""

    # "Exclude" is always offered: an excluded channel is left out of the
    # annotation map, so it is never normalized, summarized, modeled or shown.
    cond_choices <- c(conds, if (has_br) "Norm", "Exclude")

    # ── Initialise / maintain persistent assignment vector ──────────────
    if (is.null(store$grid_assignments) || length(store$grid_assignments) != length(channels)) {
      store$grid_assignments <- rep("", length(channels))
    }
    asgn <- store$grid_assignments
    asgn[!(asgn %in% c(cond_choices, ""))] <- ""
    store$grid_assignments <- asgn
    assignments <- store$grid_assignments
    paint_cond  <- store$paint_cond

    # ── Channel → file index map (for auto-plex defaults) ───────────────
    ch_file_idx <- tmt_default_channel_plex(store$files, channels)

    # ── Colour palette (12 distinct colours rotating by condition index) ─
    .COND_COLORS <- c("#2196F3","#E91E63","#FF9800","#4CAF50","#9C27B0",
                      "#00BCD4","#795548","#F44336","#3F51B5","#CDDC39",
                      "#607D8B","#FF5722")
    .cond_bg <- function(cn) {
      if (cn == "Norm") return("#607D8B")
      if (cn == "Exclude") return("#343a40")
      idx <- match(cn, conds)
      if (is.na(idx)) return("#aaaaaa")
      .COND_COLORS[((idx - 1L) %% length(.COND_COLORS)) + 1L]
    }

    # ── JS state: bake current state into client-side handler ─────────
    cond_map  <- lapply(cond_choices, function(cn) list(name = cn, color = .cond_bg(cn)))
    asgn_json <- jsonlite::toJSON(as.list(assignments), auto_unbox = TRUE)
    cmap_json <- jsonlite::toJSON(cond_map,             auto_unbox = TRUE)
    pfx_json  <- jsonlite::toJSON(gp(""),               auto_unbox = TRUE)
    pc_json   <- if (is.null(paint_cond)) "null" else
                   jsonlite::toJSON(paint_cond, auto_unbox = TRUE)

    # ── Paint-mode palette buttons (one per condition) ───────────────────
    palette_buttons <- tags$div(
      class = "d-flex flex-wrap gap-2 mb-1 p-2 rounded",
      style = "background:#f8f9fa;",
      tags$span(class = "small text-muted me-1 align-self-center",
                icon("paintbrush"), " Paint:"),
      lapply(seq_along(cond_choices), function(j) {
        cn        <- cond_choices[j]
        is_active <- !is.null(paint_cond) && paint_cond == cn
        bg        <- .cond_bg(cn)
        tags$button(
          type              = "button",
          class             = "btn btn-sm grid-paint-btn",
          `data-cond-name`  = cn,
          `data-cond-color` = bg,
          style             = paste0(
            "background-color:", bg, "; color:white; border:2px solid ",
            if (is_active) "rgba(0,0,0,0.5);" else "transparent;",
            if (is_active) " box-shadow:0 0 0 2px white, 0 0 0 4px rgba(0,0,0,0.4);" else ""
          ),
          cn
        )
      })
    )

    # ── Grid rows ────────────────────────────────────────────────────────
    grid_rows <- lapply(seq_along(channels), function(i) {
      ch  <- channels[i]
      tag <- tmt_extract_tag(ch)
      lbl <- tmt_short_label(ch)

      # Auto-plex default from file index (sticky user override preserved)
      plex_choices <- as.character(seq_len(np))
      file_idx     <- ch_file_idx[[ch]]
      default_plex <- if (!is.na(file_idx) && as.character(file_idx) %in% plex_choices)
                        as.character(file_idx) else plex_choices[1]
      prev_plex    <- input[[paste0(gp("grid_plex_"), i)]]
      plex_sel     <- if (!is.null(prev_plex) && prev_plex %in% plex_choices) prev_plex
                      else default_plex

      # Condition cell  -" initial colour set by R; JS takes over for live updates
      assigned <- if (i <= length(assignments)) assignments[[i]] else ""
      has_cond <- !is.na(assigned) && nchar(assigned) > 0
      cond_bg  <- if (has_cond) .cond_bg(assigned) else "#e9ecef"
      cond_txt <- if (has_cond) assigned else "\u2014 unassigned \u2014"
      cond_col <- if (has_cond) "white" else "#6c757d"

      cond_cell <- tags$button(
        type          = "button",
        class         = "btn btn-sm w-100 grid-ch-btn",
        `data-ch-idx` = as.character(i - 1L),
        style         = paste0(
          "background-color:", cond_bg, "; color:", cond_col,
          "; border:1px solid rgba(0,0,0,0.12); text-align:left; cursor:pointer;"
        ),
        cond_txt
      )

      tags$tr(class = "grid-ch-row",
        tags$td(class = "grid-ch-label",
                style = "vertical-align:middle; font-size:0.85em; white-space:normal; word-break:break-word;",
                tags$span(class = "badge bg-light text-dark", tag),
                " ", tags$span(class = "text-muted", lbl)),
        tags$td(style = "vertical-align:middle; width:90px;",
          if (is_multi) selectInput(paste0(gp("grid_plex_"), i), label = NULL,
                                    choices = plex_choices, selected = plex_sel, width = "80px")
          else tags$span(class = "badge bg-secondary", "1")
        ),
        tags$td(style = "vertical-align:middle; min-width:160px;",
          cond_cell
        )
      )
    })

    plex_section <- NULL
    if (is_multi) {
      plex_section <- tags$div(class = "d-flex align-items-end gap-2 mb-3",
        tags$div(style = "flex:0 0 120px;",
          numericInput(gp("grid_n_plexes"), "Number of Plexes:",
                       value = np, min = 2, max = 20, step = 1, width = "100%")),
        actionButton(gp("grid_update_plexes"), "Update",
                     class = "btn-outline-secondary btn-sm mb-3", icon = icon("refresh"))
      )
    }

    # Condition badges with delete (X) button, colored by palette
    cond_badges <- if (length(conds) > 0)
      tags$div(class = "d-flex flex-wrap gap-1 mb-2",
        lapply(seq_along(conds), function(j) {
          cn <- conds[j]
          tags$span(
            class = "badge d-inline-flex align-items-center gap-1 pe-1",
            style = paste0("background-color:", .cond_bg(cn), "; color:white;"),
            cn,
            tags$button(
              type    = "button",
              class   = "btn-close btn-close-white",
              style   = "font-size:0.55em; opacity:0.85;",
              title   = paste0("Remove \u2018", cn, "\u2019"),
              onclick = paste0("Shiny.setInputValue('", gp("grid_del_cond"),
                               "', ", j, ", {priority:'event'});")
            )
          )
        }))
    else NULL

    cond_section <- tags$div(class = "mb-3",
      tags$label(class = "form-label fw-bold", icon("tags"), " Conditions"),
      cond_badges,
      palette_buttons,
      tags$small(id = gp("grid_paint_status"), class = "text-muted d-block mb-2",
        if (!is.null(paint_cond) && nchar(paint_cond) > 0)
          tags$span(icon("paintbrush"), paste0(" Painting \u2018", paint_cond,
            "\u2019 \u2014 click rows to assign. Click the button again to exit paint mode."))
        else
          "Click a condition above to activate paint mode, then click channel rows to assign."
      ),
      tags$small(class = "text-muted d-block mb-2",
        tags$b("Exclude"), " leaves a channel out of the analysis entirely (e.g. an IgG ",
        "control or an unused channel): it is not normalized, summarized, tested or shown, ",
        "and is listed in the methods text. Decide exclusions before looking at results."),
      tags$div(class = "d-flex align-items-end gap-2",
        tags$div(style = "flex:0 0 200px;",
          textInput(gp("grid_new_cond"), label = NULL,
                    placeholder = cond_placeholder, width = "100%")),
        actionButton(gp("grid_add_cond"), "Add Condition",
                     class = "btn-outline-primary btn-sm", icon = icon("plus"))
      )
    )

    bridge_section <- NULL
    if (is_multi) {
      bridge_section <- tags$div(class = "mb-3",
        tags$label(class = "form-label fw-bold", icon("bridge"), " Pooled Bridge Channel?"),
        tags$div(class = "d-flex gap-2",
          actionButton(gp("grid_bridge_yes"),
                       if (has_br) tags$b("Yes (active)") else "Yes",
                       class = if (has_br) "btn-success btn-sm" else "btn-outline-secondary btn-sm",
                       icon = icon("check")),
          actionButton(gp("grid_bridge_no"),
                       if (!has_br) tags$b("No (active)") else "No",
                       class = if (!has_br) "btn-success btn-sm" else "btn-outline-secondary btn-sm",
                       icon = icon("xmark"))
        ),
        if (has_br) tags$p(class = "small text-muted mt-1",
          "Assign the pooled channel(s) to ", tags$b("Norm"),
          " (at least one per plex; several are averaged).")
      )
    }

    # ── Client-side JS: instant paint/assign without server round-trips ─
    grid_js <- paste0(
      "(function() {\n",
      "  var prefix = ", pfx_json, ";\n",
      "  var assignments = ", asgn_json, ";\n",
      "  var condMap = ", cmap_json, ";\n",
      "  var paintCond = ", pc_json, ";\n",
      "\n",
      "  function getCondColor(name) {\n",
      "    for (var i = 0; i < condMap.length; i++) {\n",
      "      if (condMap[i].name === name) return condMap[i].color;\n",
      "    }\n",
      "    return '#aaaaaa';\n",
      "  }\n",
      "\n",
      "  function updatePaletteUI() {\n",
      "    document.querySelectorAll('.grid-paint-btn').forEach(function(btn) {\n",
      "      var cond   = btn.getAttribute('data-cond-name');\n",
      "      var active = (paintCond !== null && paintCond === cond);\n",
      "      btn.style.border    = active ? '2px solid rgba(0,0,0,0.5)' : '2px solid transparent';\n",
      "      btn.style.boxShadow = active ? '0 0 0 2px white, 0 0 0 4px rgba(0,0,0,0.4)' : 'none';\n",
      "    });\n",
      "    var statusEl = document.getElementById(prefix + 'grid_paint_status');\n",
      "    if (statusEl) {\n",
      "      statusEl.textContent = paintCond\n",
      "        ? 'Painting \\u2018' + paintCond + '\\u2019 \\u2014 click rows to assign. Click the button again to exit paint mode.'\n",
      "        : 'Click a condition above to activate paint mode, then click channel rows to assign.';\n",
      "    }\n",
      "    document.querySelectorAll('.grid-ch-btn').forEach(function(btn) {\n",
      "      btn.style.cursor = paintCond ? 'pointer' : 'default';\n",
      "    });\n",
      "    var shownBtn = document.querySelector('.grid-paint-shown-btn');\n",
      "    if (shownBtn) shownBtn.disabled = !paintCond;\n",
      "  }\n",
      "\n",
      "  // Channel search: show rows whose visible channel text contains the term\n",
      "  function applyFilter() {\n",
      "    var box  = document.getElementById(prefix + 'grid_search');\n",
      "    var term = box ? box.value.trim().toLowerCase() : '';\n",
      "    var rows = document.querySelectorAll('.grid-ch-row');\n",
      "    var nShown = 0;\n",
      "    rows.forEach(function(row) {\n",
      "      var txt  = row.querySelector('.grid-ch-label').textContent.toLowerCase();\n",
      "      var show = (term === '' || txt.indexOf(term) !== -1);\n",
      "      row.style.display = show ? '' : 'none';\n",
      "      if (show) nShown++;\n",
      "    });\n",
      "    var countEl = document.getElementById(prefix + 'grid_search_count');\n",
      "    if (countEl) countEl.textContent = (term === '') ? rows.length + ' channels'\n",
      "      : nShown + ' of ' + rows.length + ' channels shown';\n",
      "  }\n",
      "\n",
      "  function updateChannelBtn(btn, idx) {\n",
      "    var asgn = assignments[idx] || '';\n",
      "    if (asgn.length > 0) {\n",
      "      btn.style.backgroundColor = getCondColor(asgn);\n",
      "      btn.style.color           = 'white';\n",
      "      btn.textContent           = asgn;\n",
      "    } else {\n",
      "      btn.style.backgroundColor = '#e9ecef';\n",
      "      btn.style.color           = '#6c757d';\n",
      "      btn.textContent           = '\\u2014 unassigned \\u2014';\n",
      "    }\n",
      "  }\n",
      "\n",
      "  var modal = document.querySelector('.modal-body');\n",
      "  if (!modal) return;\n",
      "  modal.addEventListener('click', function(e) {\n",
      "    var paintBtn = e.target.closest('.grid-paint-btn');\n",
      "    if (paintBtn) {\n",
      "      var cond = paintBtn.getAttribute('data-cond-name');\n",
      "      paintCond = (paintCond === cond) ? null : cond;\n",
      "      Shiny.setInputValue(prefix + 'grid_paint_sync', paintCond === null ? '' : paintCond);\n",
      "      updatePaletteUI();\n",
      "      return;\n",
      "    }\n",
      "    var chBtn = e.target.closest('.grid-ch-btn');\n",
      "    if (chBtn) {\n",
      "      if (!paintCond) return;\n",
      "      var idx = parseInt(chBtn.getAttribute('data-ch-idx'));\n",
      "      assignments[idx] = (assignments[idx] === paintCond) ? '' : paintCond;\n",
      "      updateChannelBtn(chBtn, idx);\n",
      "      Shiny.setInputValue(prefix + 'grid_asgn_sync', JSON.stringify(assignments));\n",
      "      return;\n",
      "    }\n",
      "    // Paint every channel the search currently shows (assign, never toggle)\n",
      "    if (e.target.closest('.grid-paint-shown-btn')) {\n",
      "      if (!paintCond) return;\n",
      "      document.querySelectorAll('.grid-ch-row').forEach(function(row) {\n",
      "        if (row.style.display === 'none') return;\n",
      "        var btn = row.querySelector('.grid-ch-btn');\n",
      "        var idx = parseInt(btn.getAttribute('data-ch-idx'));\n",
      "        assignments[idx] = paintCond;\n",
      "        updateChannelBtn(btn, idx);\n",
      "      });\n",
      "      Shiny.setInputValue(prefix + 'grid_asgn_sync', JSON.stringify(assignments));\n",
      "    }\n",
      "  });\n",
      "  var searchBox = document.getElementById(prefix + 'grid_search');\n",
      "  if (searchBox) searchBox.addEventListener('input', applyFilter);\n",
      "\n",
      "  updatePaletteUI();\n",
      "  applyFilter();\n",
      "})();\n"
    )

    title_text <- if (is_multi) "Multi-Plex Design" else "Single Plex Design"
    showModal(modalDialog(
      title = tags$span(icon("table-cells-large"), " ", title_text),
      tagList(plex_section, cond_section, bridge_section, tags$hr(),
        tags$label(class = "form-label fw-bold", icon("table"), " Assign each channel:"),
        tags$div(class = "d-flex gap-2 align-items-start",
          tags$div(class = "flex-grow-1",
            textInput(gp("grid_search"), NULL, value = search_term, width = "100%",
                      placeholder = "Search channels (contains)...")),
          tags$button(type = "button",
                      class = "btn btn-outline-primary btn-sm grid-paint-shown-btn",
                      title = "Assign the active paint condition to every channel shown",
                      icon("fill-drip"), " Paint shown")
        ),
        tags$small(id = gp("grid_search_count"), class = "text-muted d-block mb-1"),
        tags$div(style = "max-height:400px; overflow-y:auto;",
          tags$table(class = "table table-sm table-hover mb-0",
            tags$thead(tags$tr(tags$th("Channel"), tags$th("Plex"), tags$th("Condition"))),
            tags$tbody(grid_rows)
          )
        )
      ),
      size = "l", easyClose = FALSE,
      footer = tagList(modalButton("Cancel"),
        actionButton(gp("grid_finish"), "Finish Design",
                     class = "btn-success", icon = icon("check")))
    ))
    shinyjs::runjs(grid_js)
  }

  # ── Grid action handlers ────────────────────────────────────────────────
  observeEvent(input[[gp("grid_update_plexes")]], {
    new_np <- input[[gp("grid_n_plexes")]]
    if (is.null(new_np) || is.na(new_np) || new_np < 2) {
      showNotification("Enter at least 2.", type = "warning"); return()
    }
    store$n_plexes <- as.integer(new_np)
    removeModal(); Sys.sleep(0.15); .show_design_grid(keep_search = TRUE)
  })

  observeEvent(input[[gp("grid_add_cond")]], {
    nm <- trimws(input[[gp("grid_new_cond")]] %||% "")
    if (nchar(nm) == 0) { showNotification("Enter a condition name.", type = "warning"); return() }
    if (tolower(nm) %in% c("norm", "exclude")) {
      showNotification(paste0("'", nm, "' is reserved."), type = "warning"); return()
    }
    if (nm %in% store$design_conds) { showNotification("Already exists.", type = "warning"); return() }
    store$design_conds <- c(store$design_conds, nm)
    removeModal(); Sys.sleep(0.15); .show_design_grid(keep_search = TRUE)
  })

  observeEvent(input[[gp("grid_bridge_yes")]], {
    store$has_bridge <- TRUE; removeModal(); Sys.sleep(0.15); .show_design_grid(keep_search = TRUE)
  })
  observeEvent(input[[gp("grid_bridge_no")]], {
    store$has_bridge <- FALSE; removeModal(); Sys.sleep(0.15); .show_design_grid(keep_search = TRUE)
  })

  # ── Delete condition badge ──────────────────────────────────────────────
  observeEvent(input[[gp("grid_del_cond")]], {
    idx   <- as.integer(input[[gp("grid_del_cond")]])
    conds <- store$design_conds
    if (idx < 1L || idx > length(conds)) return()
    cn_to_remove       <- conds[idx]
    store$design_conds <- conds[-idx]
    asgn <- store$grid_assignments
    if (!is.null(asgn)) {
      asgn[!is.na(asgn) & asgn == cn_to_remove] <- ""
      store$grid_assignments <- asgn
    }
    if (!is.null(store$paint_cond) && store$paint_cond == cn_to_remove)
      store$paint_cond <- NULL
    removeModal(); Sys.sleep(0.15); .show_design_grid(keep_search = TRUE)
  }, ignoreInit = TRUE)

  # ── JS sync: paint state (no modal rebuild) ────────────────────────────
  observeEvent(input[[gp("grid_paint_sync")]], {
    val <- input[[gp("grid_paint_sync")]]
    store$paint_cond <- if (is.null(val) || val == "") NULL else as.character(val)
  }, ignoreInit = TRUE)

  # ── JS sync: assignment vector (no modal rebuild) ──────────────────────
  observeEvent(input[[gp("grid_asgn_sync")]], {
    json_str <- input[[gp("grid_asgn_sync")]]
    if (!is.null(json_str) && nchar(json_str) > 0) {
      asgn <- jsonlite::fromJSON(json_str)
      store$grid_assignments <- as.character(asgn)
    }
  }, ignoreInit = TRUE)

  # ── Finish Design ───────────────────────────────────────────────────────
  observeEvent(input[[gp("grid_finish")]], {
    channels <- store$all_channels
    n        <- length(channels)
    np       <- store$n_plexes
    conds    <- store$design_conds
    has_br   <- store$has_bridge
    is_multi <- np > 1

    if (length(conds) == 0 && !has_br) {
      showNotification("Add at least one condition.", type = "warning"); return()
    }

    grid_plex <- integer(n); grid_cond <- character(n)
    asgn_vec  <- store$grid_assignments
    for (i in seq_len(n)) {
      gp_val <- if (is_multi) input[[paste0(gp("grid_plex_"), i)]] else 1L
      grid_plex[i] <- as.integer(gp_val %||% 1L)
      grid_cond[i] <- if (!is.null(asgn_vec) && i <= length(asgn_vec) && !is.na(asgn_vec[i]))
                        asgn_vec[i]
                      else ""
    }

    empty <- which(grid_cond == "" | grid_cond == "-- select --")
    if (length(empty) > 0) {
      showNotification(paste0(length(empty), " channel(s) unassigned. Assign a condition, ",
                              "or 'Exclude' to leave a channel out."), type = "warning"); return()
    }

    # Excluded channels are dropped here, before the annotation map is built, so
    # every check below and every downstream step sees only the kept channels.
    excl          <- grid_cond == "Exclude"
    excluded_chs  <- channels[excl]
    channels      <- channels[!excl]
    grid_plex     <- grid_plex[!excl]
    grid_cond     <- grid_cond[!excl]

    real_conds <- unique(grid_cond[grid_cond != "Norm"])
    if (length(real_conds) < 2) {
      showNotification("Need at least 2 experimental conditions.", type = "warning"); return()
    }

    # One or more pooled reference channels per plex. apply_irs_normalization()
    # takes their per-protein mean (log2), so a plex with two pools uses both.
    if (has_br) {
      for (pid in seq_len(np)) {
        nn <- sum(grid_cond == "Norm" & grid_plex == pid)
        if (nn < 1) {
          showNotification(paste0("Plex ", pid, ": needs at least 1 Norm (has 0)."),
                           type = "error"); return()
        }
      }
    }

    if (is_multi && !has_br) {
      ref <- NULL; ok <- TRUE
      for (pid in seq_len(np)) {
        ct <- sort(table(grid_cond[grid_plex == pid]))
        if (is.null(ref)) ref <- ct else if (!identical(ct, ref)) { ok <- FALSE; break }
      }
      if (!ok) {
        showModal(modalDialog(
          title = tags$span(icon("triangle-exclamation"), " Unbalanced Design"),
          tags$div(class = "alert alert-danger",
            "Multi-plex without bridge requires balanced conditions across plexes."),
          size = "m", easyClose = TRUE, footer = modalButton("OK")))
        return()
      }
    }

    ann <- tmt_grid_annotation(channels, grid_plex, grid_cond)

    has_bridge_in_ann <- any(ann$Condition == "Norm")
    NL <- if (!is_multi) 1L else if (has_bridge_in_ann) 2L else 3L

    store$design_option     <- NL
    store$annotation_map    <- ann
    store$excluded_channels <- excluded_chs
    conds_final <- sort(real_conds)
    store$conditions     <- conds_final
    removeModal()

    # Share with other TMT modules
    rv$tmt_annotation_map     <- ann
    rv$tmt_excluded_channels  <- excluded_chs
    rv$tmt_design_option      <- NL
    rv$tmt_conditions         <- conds_final
    rv$tmt_design_source      <- design_source_tag

    # Invalidate TMT master cache (design changed)
    if (!is.null(rv$tmt_master_cache)) {
      rv$tmt_master_cache        <- NULL
      rv$tmt_master_cache_design <- NULL
      rv$tmt_master_cache_ann    <- NULL
      showModal(modalDialog(
        title = "Preprocessing Parameter Changed",
        "The experimental design was modified. The cached TMT summarization has been ",
        "invalidated. Please re-run your analysis.",
        easyClose = FALSE,
        footer = modalButton("OK")
      ))
      logger(paste0(module_name, ": Design changed  -  TMT master cache invalidated"))
    }

    opt_label <- switch(as.character(NL),
      "1" = "Single Plex", "2" = "Multi-Plex + Bridge", "3" = "Multi-Plex Balanced")
    logger(paste0(module_name, ": Design  -  ", opt_label, ", ",
                  nrow(ann), " channels, ", length(real_conds), " conditions",
                  if (length(excluded_chs) > 0)
                    paste0(", ", length(excluded_chs), " channel(s) excluded: ",
                           paste(tmt_short_label(excluded_chs), collapse = "; "))
                  else "", "."))
    store$status <- paste0("Design: ", opt_label, ". ", length(real_conds), " conditions.")
    showNotification("Experimental design saved!", type = "message")

    # Module-specific callback (e.g. volcano comparison choosers, trends cond_order)
    if (is.function(on_design_finished)) {
      on_design_finished(ann, NL, conds_final)
    }
  })

  # Return callable handles
  list(show_arch_fork = .show_arch_fork, show_design_grid = .show_design_grid)
}

# =============================================================================

# =============================================================================
# LOOKUP_GENE_NAMES_UNIPROT
# Retrieves gene names from UniProt using the official ID Mapping API.
# Submits all accessions in a single multipart POST to /idmapping/run,
# polls for completion, then downloads results as TSV.
# Supports up to 100,000 accessions in a single job (UniProt limit).
#
# Flow:
#   1) POST /idmapping/run  (multipart, as per official R example)
#   2) Poll /idmapping/status/{jobId} until FINISHED
#   3) GET  /idmapping/uniprotkb/results/stream/{jobId}?fields=accession,gene_primary&format=tsv
#
# Args:
#   accessions  -" character vector of UniProt accession numbers
#   logger      -" logging function
#
# Returns: named character vector  (names = accessions, values = gene symbols)
#          Accessions for which no gene was found keep themselves as value.
# =============================================================================
lookup_gene_names_uniprot <- function(accessions, logger = message) {
  accessions <- unique(accessions[!is.na(accessions) & nzchar(accessions)])
  if (length(accessions) == 0L)
    return(setNames(character(0), character(0)))

  if (!requireNamespace("httr", quietly = TRUE)) {
    logger("httr not available \u2014 cannot query UniProt REST API")
    return(setNames(accessions, accessions))
  }

  orig_accessions <- accessions
  clean_acc <- trimws(accessions)
  valid <- grepl("^[A-Za-z0-9][-A-Za-z0-9]{4,14}$", clean_acc)
  if (sum(valid) == 0L) {
    logger("UniProt gene lookup: no valid UniProt accessions found after cleaning")
    return(setNames(orig_accessions, orig_accessions))
  }
  if (sum(!valid) > 0L)
    logger(paste0("UniProt gene lookup: skipping ", sum(!valid),
                  " non-UniProt-format accessions (e.g. '",
                  clean_acc[!valid][1], "')"))

  query_acc <- clean_acc[valid]
  logger(paste0("UniProt gene lookup: submitting ", length(query_acc),
                " accessions to ID Mapping API (single job)..."))

  # ── Step 1: Submit ID mapping job (multipart, as per official UniProt R docs) ──
  submit_resp <- tryCatch(
    httr::POST(
      url    = "https://rest.uniprot.org/idmapping/run",
      body   = list(ids  = paste(query_acc, collapse = ","),
                    from = "UniProtKB_AC-ID",
                    to   = "UniProtKB"),
      encode = "multipart",
      httr::accept_json(),
      httr::user_agent("OSMSS/1.0")
    ),
    error = function(e) { logger(paste0("UniProt API submit error: ", e$message)); NULL }
  )
  if (is.null(submit_resp) || httr::status_code(submit_resp) != 200L) {
    sc <- if (!is.null(submit_resp)) httr::status_code(submit_resp) else "NULL"
    body_msg <- tryCatch(httr::content(submit_resp, "text", encoding = "UTF-8"), error = function(e) "")
    logger(paste0("UniProt API submit failed (HTTP ", sc, "): ", substr(body_msg, 1, 300)))
    return(setNames(orig_accessions, orig_accessions))
  }

  job_id <- tryCatch(
    httr::content(submit_resp, as = "parsed")[["jobId"]],
    error = function(e) NULL
  )
  if (is.null(job_id) || !nzchar(job_id)) {
    logger("UniProt API: no job ID returned")
    return(setNames(orig_accessions, orig_accessions))
  }
  logger(paste0("UniProt API: job ", job_id, " submitted \u2014 polling for results..."))

  # ── Step 2: Poll until results appear (max ~90 s) ──
  status_url <- paste0("https://rest.uniprot.org/idmapping/status/", job_id)
  complete   <- FALSE
  for (attempt in seq_len(30L)) {
    Sys.sleep(3)
    sr <- tryCatch(
      httr::GET(status_url, httr::accept_json(), httr::user_agent("OSMSS/1.0")),
      error = function(e) NULL
    )
    if (is.null(sr)) next
    st <- tryCatch(httr::content(sr, as = "parsed"), error = function(e) list())
    # results or failedIds present → job done (may also be a 303 redirect followed by httr)
    if (!is.null(st[["results"]]) || !is.null(st[["failedIds"]])) {
      complete <- TRUE; break
    }
    js <- st[["jobStatus"]]
    if (identical(js, "FINISHED")) { complete <- TRUE; break }
    if (identical(js, "ERROR"))    {
      logger("UniProt API: job failed on server side"); break
    }
  }
  if (!complete) {
    logger("UniProt API: timeout waiting for results")
    return(setNames(orig_accessions, orig_accessions))
  }

  # ── Step 3: Fetch results as TSV ──
  results_url <- paste0(
    "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/",
    job_id,
    "?fields=accession,gene_primary&format=tsv"
  )
  tsv_resp <- tryCatch(
    httr::GET(results_url, httr::user_agent("OSMSS/1.0")),
    error = function(e) { logger(paste0("UniProt API fetch error: ", e$message)); NULL }
  )
  if (is.null(tsv_resp) || httr::status_code(tsv_resp) != 200L) {
    sc2 <- if (!is.null(tsv_resp)) httr::status_code(tsv_resp) else "NULL"
    logger(paste0("UniProt API: failed to fetch results TSV (HTTP ", sc2, ")"))
    return(setNames(orig_accessions, orig_accessions))
  }

  tsv_text <- httr::content(tsv_resp, "text", encoding = "UTF-8")
  result_df <- tryCatch(
    read.delim(text = tsv_text, stringsAsFactors = FALSE, header = TRUE),
    error = function(e) { logger(paste0("UniProt TSV parse error: ", e$message)); NULL }
  )
  if (is.null(result_df) || nrow(result_df) == 0L) {
    logger("UniProt API: empty result set")
    return(setNames(orig_accessions, orig_accessions))
  }

  # ── Step 4: Parse TSV  ──
  # idmapping stream returns columns: From, Entry, Gene Names (primary)
  # read.delim converts "Gene Names (primary)" -> "Gene.Names..primary."
  acc_col  <- if ("From" %in% names(result_df)) "From" else names(result_df)[1]
  gene_col <- NULL
  for (cn in c("Gene.Names..primary.", "Gene.Names.primary.",
               "Gene.Names", "gene_primary")) {
    if (cn %in% names(result_df)) { gene_col <- cn; break }
  }
  if (is.null(gene_col)) {
    non_id <- setdiff(names(result_df), c("From", "Entry"))
    if (length(non_id) > 0L) gene_col <- non_id[1L]
  }
  if (is.null(gene_col)) {
    logger(paste0("UniProt API: unexpected columns: ", paste(names(result_df), collapse = ", ")))
    return(setNames(orig_accessions, orig_accessions))
  }

  out       <- setNames(orig_accessions, orig_accessions)
  acc_vals  <- trimws(as.character(result_df[[acc_col]]))
  gene_vals <- trimws(as.character(result_df[[gene_col]]))
  for (i in seq_along(acc_vals)) {
    a <- acc_vals[i]; g <- gene_vals[i]
    if (a %in% names(out) && !is.na(g) && nzchar(g)) out[a] <- g
  }

  n_found <- sum(out != names(out))
  logger(paste0("UniProt gene lookup: resolved ", n_found, " / ",
                length(orig_accessions), " accessions to gene names"))
  out
}

# =============================================================================
# QF_SUMMARIZE_TO_PROTEIN_MATRIX
# Aggregates peptide-level data to protein-level using QFeatures /
# MsCoreUtils robustSummary (MASS::rlm). Preserves genuine NAs (proteins with
# all-NA peptides in a sample remain NA in the output matrix).
#
# Args:
#   input_df          -" data.frame of peptide data
#   input_type        -" "flat_peptide" (the only supported value)
#   protein_col       -" column name containing protein accessions
#   peptide_col       -" column name containing peptide sequences
#   intensity_cols    -" character vector of intensity column names
#   sample_names        -" optional names for output matrix columns
#   extra_rowdata_cols  -" optional character vector of additional columns from
#                        input_df to carry through as rowData (e.g. metadata)
#   qvalue_col        -" optional column name for Q-value filtering
#   qvalue_cutoff     -" Q-value threshold (default 0.01)
#   logger            -" logging function
#
# Returns: list(protein_matrix, acc_to_gene, n_peptides_input, n_proteins_output)
#          OR list(error = "message") on failure.
# =============================================================================
qf_summarize_to_protein_matrix <- function(input_df,
                                            input_type,
                                            protein_col,
                                            peptide_col,
                                            intensity_cols,
                                            sample_names            = NULL,
                                            extra_rowdata_cols      = NULL,
                                            # Unreachable from the TMT route: no
                                            # TMT caller supplies a qvalue_col, so
                                            # the filter below never fires there and
                                            # peptide-level confidence is whatever
                                            # the search engine export already
                                            # applied. Disclosure, not wiring -
                                            # activating it would change results.
                                            qvalue_col              = NULL,
                                            qvalue_cutoff           = 0.01,
                                            normalize_peptide_level = FALSE,
                                            logger                  = message) {
  if (!isTRUE(.PRODA_AVAILABLE))
    return(list(error = paste0(
      "QFeatures / MsCoreUtils not installed. ",
      "Install via: BiocManager::install(c('proDA','QFeatures','MsCoreUtils'))")))
  ensure_proda_loaded()

  # ── Validate inputs ──
  if (!protein_col %in% colnames(input_df))
    return(list(error = paste0("qf_summarize: protein column '", protein_col, "' not found.")))
  if (!peptide_col %in% colnames(input_df))
    return(list(error = paste0("qf_summarize: peptide column '", peptide_col, "' not found.")))
  missing_int <- setdiff(intensity_cols, colnames(input_df))
  if (length(missing_int) > 0)
    return(list(error = paste0("qf_summarize: intensity columns missing: ",
                               paste(missing_int, collapse = ", "))))

  # ── Q-value filter ──
  if (!is.null(qvalue_col) && qvalue_col %in% colnames(input_df)) {
    n_before <- nrow(input_df)
    input_df <- input_df[
      is.na(input_df[[qvalue_col]]) | input_df[[qvalue_col]] <= qvalue_cutoff, ]
    logger(paste0("qf_summarize: Q-value filter (\u2264", qvalue_cutoff, "): ",
                  n_before, " \u2192 ", nrow(input_df), " rows"))
    if (nrow(input_df) < 10)
      return(list(error = paste0("qf_summarize: fewer than 10 rows after Q-value filtering.")))
  }

  # ── Build intensity matrix (rows = features, cols = samples) ──
  int_mat <- as.matrix(input_df[, intensity_cols, drop = FALSE])
  storage.mode(int_mat) <- "numeric"
  # Feature-level record of cells that were present but non-positive. These are
  # about to become NA and would otherwise be indistinguishable from genuine
  # gaps. Carried through aggregation and returned as zero_origin_mask so
  # downstream imputation can decline to invent a value for them.
  .zero_feat <- !is.na(int_mat) & int_mat <= 0
  int_mat[int_mat <= 0] <- NA
  int_mat <- log2(int_mat)

  # Column-wise median centering at the feature (peptide/PSM) level, BEFORE
  # aggregation to proteins.  Only applied when normalize_peptide_level = TRUE.
  # This is the QFeatures-vignette-correct order: normalise features first,
  # then aggregate via robustSummary, so aggregation averages already-
  # normalised values rather than raw ones.
  # Default is FALSE to preserve existing behaviour for callers (e.g. TMT)
  # that apply their own normalisation after aggregation.
  if (isTRUE(normalize_peptide_level)) {
    col_medians   <- apply(int_mat, 2, median, na.rm = TRUE)
    global_median <- median(col_medians, na.rm = TRUE)
    int_mat       <- sweep(int_mat, 2, col_medians - global_median)
    logger("qf_summarize: column-wise median centering applied at feature level (before aggregation)")
  }

  if (!is.null(sample_names) && length(sample_names) == ncol(int_mat))
    colnames(int_mat) <- sample_names

  n_pep_input <- nrow(int_mat)
  logger(paste0("qf_summarize: ", n_pep_input, " feature rows, ",
                ncol(int_mat), " samples"))

  # ── Build rowData for QFeatures ──
  rd_cols <- c(protein_col, peptide_col)
  extra_meta <- if (!is.null(extra_rowdata_cols))
    intersect(extra_rowdata_cols, colnames(input_df)) else character(0)
  if (length(extra_meta) > 0) rd_cols <- c(rd_cols, extra_meta)
  row_data <- as.data.frame(input_df[, rd_cols, drop = FALSE], stringsAsFactors = FALSE)
  colnames(row_data)[1:2] <- c("Protein", "Peptide")

  # Strip semicolon-separated protein group suffixes (PD writes "P06310; A0A075B6S2");
  # keep only the first (canonical master) accession.
  row_data$Protein <- trimws(sub(";.*", "", row_data$Protein))

  # Drop rows where Protein or Peptide is NA/empty  -" aggregateFeatures cannot
  # group by NA fcol values and will throw "missing value where TRUE/FALSE needed"
  valid_rows <- !is.na(row_data$Protein) & nzchar(row_data$Protein) &
                !is.na(row_data$Peptide) & nzchar(row_data$Peptide)
  n_dropped <- 0L
  if (!all(valid_rows)) {
    n_dropped <- sum(!valid_rows)
    logger(paste0("qf_summarize: dropping ", n_dropped,
                  " rows with NA/empty Protein or Peptide"))
    row_data <- row_data[valid_rows, , drop = FALSE]
    int_mat  <- int_mat[valid_rows, , drop = FALSE]
    .zero_feat <- .zero_feat[valid_rows, , drop = FALSE]
  }
  if (nrow(row_data) < 1)
    return(list(error = "qf_summarize: no rows remain after removing NA protein/peptide entries."))

  rownames(int_mat)  <- make.unique(row_data$Peptide)
  rownames(row_data) <- rownames(int_mat)

  se_feat <- tryCatch(
    SummarizedExperiment::SummarizedExperiment(
      assays   = list(intensity = int_mat),
      rowData  = row_data
    ),
    error = function(e) { logger(paste0("qf_summarize: SummarizedExperiment error: ", e$message)); NULL }
  )
  if (is.null(se_feat))
    return(list(error = "qf_summarize: failed to build SummarizedExperiment."))

  # ── Build QFeatures object ──
  named_list <- list(se_feat); names(named_list) <- "peptides"
  # colData rownames must match SE colnames for QFeatures to map them correctly
  col_df <- S4Vectors::DataFrame(sample = colnames(int_mat))
  rownames(col_df) <- colnames(int_mat)
  qf <- tryCatch(
    QFeatures::QFeatures(named_list, colData = col_df),
    error = function(e) { logger(paste0("qf_summarize: QFeatures build error: ", e$message)); NULL }
  )
  if (is.null(qf))
    return(list(error = "qf_summarize: failed to build QFeatures object."))

  # ── Aggregation using robustSummary with colMeans fallback ──
  # robustSummary calls MASS::rlm which requires ≥2 rows (features) per group.
  # For singleton groups (1 peptide per protein) we fall back
  # to column means, which is equivalent and avoids the contrast-level error.
  # Three estimators can coexist in one output matrix - robustSummary, the
  # column-mean fallback, and the singleton pass-through - and until now there
  # was no record of how many proteins used which. Counted in an environment so
  # the closures below are untouched.
  .agg_n <- new.env(parent = emptyenv())
  .agg_n$colmeans  <- 0L
  .agg_n$singleton <- 0L

  .safe_robust <- function(x, ...) {
    if (nrow(x) <= 1L) {
      .agg_n$singleton <- .agg_n$singleton + 1L
      v        <- x[1L, , drop = TRUE]
      v[is.nan(v)] <- NA_real_
      return(v)
    }
    tryCatch(
      MsCoreUtils::robustSummary(x, ...),
      error = function(e) {
        .agg_n$colmeans <- .agg_n$colmeans + 1L
        apply(x, 2L, mean, na.rm = TRUE)
      }
    )
  }

  aggregate_safe <- function(qf_obj, from, to, by_col) {
    tryCatch(
      QFeatures::aggregateFeatures(
        object = qf_obj,
        i      = from,
        fcol   = by_col,
        name   = to,
        fun    = .safe_robust,
        na.rm  = TRUE
      ),
      error = function(e) { logger(paste0("qf_summarize: aggregateFeatures error: ", e$message)); NULL }
    )
  }

  # One-step: peptide → protein
  qf <- aggregate_safe(qf, from = "peptides", to = "proteins", by_col = "Protein")
  if (is.null(qf))
    return(list(error = "qf_summarize: peptide \u2192 protein aggregation failed."))

  # ── Extract protein matrix ──
  prot_se <- qf[["proteins"]]
  if (is.null(prot_se))
    return(list(error = "qf_summarize: protein-level assay missing after aggregation."))

  # aggregateFeatures names the aggregated assay slot "assay" (not "intensity"),
  # so use the first assay slot by index rather than by name.
  protein_matrix <- SummarizedExperiment::assay(prot_se, 1L)
  if (!is.null(sample_names) && length(sample_names) == ncol(protein_matrix))
    colnames(protein_matrix) <- sample_names

  logger(paste0("qf_summarize: protein matrix ", nrow(protein_matrix), " \u00d7 ",
                ncol(protein_matrix), " (proteins \u00d7 samples)"))

  # ── Build acc_to_gene mapping ──
  # Try to extract GN= from description columns if present
  gene_col <- NULL
  for (cn in c("Gene Symbol", "Gene.Symbol", "Description", "Protein Description",
               "Protein.Description")) {
    if (cn %in% colnames(input_df)) { gene_col <- cn; break }
  }
  all_accs <- rownames(protein_matrix)
  if (!is.null(gene_col)) {
    gene_lookup <- unique(input_df[, c(protein_col, gene_col), drop = FALSE])
    gene_lookup <- gene_lookup[!is.na(gene_lookup[[protein_col]]), ]
    gene_vals <- setNames(gene_lookup[[gene_col]], gene_lookup[[protein_col]])
    if (grepl("GN=", paste(gene_vals[1:min(3, length(gene_vals))], collapse = " "))) {
      # Description-style: extract GN= tag
      gene_vals <- stringr::str_replace(
        stringr::str_extract(gene_vals, "GN=([^ ]+)"), "^GN=", "")
      gene_vals[is.na(gene_vals)] <- ""
    }
    acc_to_gene <- gene_vals[all_accs]
    acc_to_gene[is.na(acc_to_gene)] <- all_accs[is.na(acc_to_gene)]
  } else {
    acc_to_gene <- setNames(all_accs, all_accs)
  }

  # Flag: TRUE when gene names came from the data; FALSE when all gene names
  # are just the accession number (i.e. no gene column was found in data).
  genes_from_data <- !is.null(gene_col) &&
    any(acc_to_gene != names(acc_to_gene), na.rm = TRUE)

  # ── Per-protein peptide count (for single-peptide flagging) ──
  pep_count <- tryCatch({
    rc <- SummarizedExperiment::rowData(prot_se)$.n
    if (!is.null(rc)) setNames(as.integer(rc), rownames(protein_matrix)) else NULL
  }, error = function(e) NULL)

  # Aggregate the feature-level zero record to protein x channel. A protein cell
  # is zero-origin only when NONE of its contributing features had an observed
  # value in that channel AND at least one was non-positive. Uses the same
  # grouping key as the peptide -> protein aggregation step.
  zero_origin_mask <- tryCatch({
    .grp  <- as.character(row_data$Protein)
    .zc   <- rowsum(.zero_feat * 1L,        group = .grp, reorder = FALSE)
    .oc   <- rowsum((!is.na(int_mat)) * 1L, group = .grp, reorder = FALSE)
    .zm   <- (.oc == 0L) & (.zc > 0L)
    .idx  <- match(rownames(protein_matrix), rownames(.zm))
    .out  <- matrix(FALSE, nrow = nrow(protein_matrix), ncol = ncol(protein_matrix),
                    dimnames = dimnames(protein_matrix))
    .ok   <- !is.na(.idx)
    if (any(.ok)) .out[.ok, ] <- .zm[.idx[.ok], , drop = FALSE]
    .out
  }, error = function(e) {
    logger(paste0("qf_summarize: zero-origin mask unavailable: ", e$message))
    NULL
  })

  list(
    protein_matrix        = protein_matrix,
    zero_origin_mask      = zero_origin_mask,
    acc_to_gene           = acc_to_gene,
    genes_from_data       = genes_from_data,
    n_peptides_input      = n_pep_input,
    n_proteins_output     = nrow(protein_matrix),
    n_rows_dropped        = n_dropped,
    pep_count_per_protein = pep_count,
    # Summarization estimator provenance. Three estimators coexist; these say
    # how many groups used each fallback rather than robustSummary.
    n_agg_fallback_colmeans = .agg_n$colmeans,
    n_agg_singleton         = .agg_n$singleton,
    # Where peptide-level FDR control came from. "inherited_from_search_engine"
    # means no q-value filter ran here: confidence is whatever the upstream
    # export (e.g. Proteome Discoverer) already applied.
    peptide_fdr_source      = if (!is.null(qvalue_col) &&
                                  qvalue_col %in% colnames(input_df))
                                "qvalue_filtered_here"
                              else "inherited_from_search_engine"
  )
}

# =============================================================================
# LF DDA UPLOAD AUTO-DETECTION
# The label-free wizard takes one file and decides what it is from its header
# row: a Proteome Discoverer protein export (nested CSV) or MaxQuant
# proteinGroups.txt. evidence.txt is asked for only after proteinGroups.txt is
# accepted. Anything else is turned away with a message saying what it looks
# like and what is accepted.
# =============================================================================

# Header row of a delimited text file (comma or tab, whichever is more common
# on the first line). character(0) when the file cannot be read.
lf_read_header <- function(path) {
  first <- tryCatch(readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) character(0))
  if (length(first) == 0L) return(character(0))
  first <- sub("^\ufeff", "", first)   # byte-order mark
  n_tab <- lengths(regmatches(first, gregexpr("\t", first, fixed = TRUE)))
  n_com <- lengths(regmatches(first, gregexpr(",",  first, fixed = TRUE)))
  hdr <- if (n_tab > n_com) {
    strsplit(first, "\t", fixed = TRUE)[[1]]
  } else {
    tryCatch(as.character(utils::read.csv(text = first, header = FALSE,
                                          stringsAsFactors = FALSE, check.names = FALSE)),
             error = function(e) character(0))
  }
  hdr <- trimws(gsub('^"|"$', "", hdr))
  hdr[!is.na(hdr) & nzchar(hdr)]
}

# Returns list(kind, message, headers).
#   kind: "pd_nested" | "mq_proteingroups" | "mq_evidence" | "fp_combined_protein" |
#         "fp_workflow" | "unsupported"
#   message: user-facing explanation; NULL for the two accepted kinds.
detect_lf_upload_type <- function(path) {
  hdr  <- lf_read_header(path)
  .out <- function(kind, message = NULL) list(kind = kind, message = message, headers = hdr)
  if (length(hdr) == 0L)
    return(.out("unsupported", "The file is empty or its header row could not be read."))
  has  <- function(x)  all(x %in% hdr)
  hasr <- function(re) any(grepl(re, hdr))
  accepted <- paste0("The label-free DDA upload accepts a Proteome Discoverer protein ",
                     "export (.csv, with peptides nested under each protein), MaxQuant ",
                     "proteinGroups.txt or FragPipe combined_protein.tsv.")

  # FragPipe: combined_protein.tsv (accepted); fragpipe.workflow is asked for after it
  if (grepl("^#\\s*FragPipe", hdr[1]))
    return(.out("fp_workflow", paste0(
      "This looks like fragpipe.workflow. Upload combined_protein.tsv first; the workflow file is ",
      "asked for after it.")))
  if (has(c("Protein", "Protein ID")) && hasr(" Intensity$") &&
      hasr("^Combined (Total Peptides|Spectral Count)$"))
    return(.out("fp_combined_protein"))
  if (has(c("Spectrum", "Peptide")) && hasr("^Spectrum File$|^Protein ID$"))
    return(.out("unsupported", paste(
      "This looks like a FragPipe PSM file (psm.tsv); OSMSS does not work with PSM files.",
      "Upload FragPipe\u2019s combined_protein.tsv instead.")))
  if (hasr("^(Peptide Sequence|Modified Sequence)$") && hasr(" Intensity$"))
    return(.out("unsupported", paste(
      "This looks like FragPipe\u2019s peptide- or ion-level table.",
      "Upload FragPipe\u2019s combined_protein.tsv instead.")))

  # MaxQuant
  if (has(c("Protein IDs", "Majority protein IDs")))
    return(.out("mq_proteingroups"))
  if (has(c("Raw file", "Sequence", "Proteins")) && hasr("^Intensity$|^MS/MS count$|^Type$"))
    return(.out("mq_evidence", paste0(
      "This looks like MaxQuant evidence.txt. Upload proteinGroups.txt first; ",
      "evidence.txt is asked for after it.")))
  if (has(c("Raw file", "Scan number")))
    return(.out("unsupported", paste("This looks like MaxQuant msms.txt.", accepted)))
  if (has(c("Sequence", "Proteins", "Leading razor protein")))
    return(.out("unsupported", paste("This looks like MaxQuant peptides.txt.", accepted)))

  # Proteome Discoverer
  if ("Accession" %in% hdr && hasr("FDR Confidence|^Abundances?[ :(]")) {
    # TMT reporter channels are labelled by tag: "Abundance: F1: 126, ..."
    if (hasr("^Abundances?( \\([^)]*\\))?: F[0-9]+: (1[23][0-9][NC]?|TMT[^,]*),"))
      return(.out("unsupported", paste0(
        "This looks like a TMT export from Proteome Discoverer. ",
        "Use the TMT upload instead.")))
    return(.out("pd_nested"))
  }
  if ("Annotated Sequence" %in% hdr) {
    if (hasr("^Spectrum File$|^First Scan$"))
      return(.out("unsupported", paste(
        "This looks like a Proteome Discoverer PSM export; OSMSS does not work with PSM files.",
        accepted)))
    return(.out("unsupported", paste(
      "This looks like a Proteome Discoverer peptide-groups export. Export the Proteins",
      "table with the peptide level nested under it instead.", accepted)))
  }

  # Other software
  if (hasr("^PG\\.ProteinGroups$|^Protein\\.Group$|^Precursor\\.Id$|^R\\.FileName$"))
    return(.out("unsupported", paste0(
      "This looks like DIA output (DIA-NN or Spectronaut). Use the DIA upload instead.")))
  if (hasr("MaxLFQ Intensity$"))
    return(.out("unsupported", paste(
      "This looks like a FragPipe export, which the label-free DDA upload does not read.",
      accepted)))

  .out("unsupported", paste("The file type could not be recognised from its column headers.",
                            accepted))
}

# Experiment names carried by proteinGroups.txt quantity columns
# ("Intensity <experiment>", or "LFQ intensity <experiment>" when a file has no
# raw Intensity columns), from the raw header text.
mq_pg_experiments <- function(hdr) {
  ex <- sub("^Intensity ", "", grep("^Intensity .", hdr, value = TRUE))
  if (length(ex) == 0L)
    ex <- sub("^LFQ intensity ", "", grep("^LFQ intensity .", hdr, value = TRUE))
  ex
}

# proteinGroups.txt quantity column (make.names form) for one evidence.txt
# Raw.file: Raw.file -> Experiment (when evidence.txt has that column), then
# raw "Intensity.<exp>" preferred over "LFQ.intensity.<exp>". NA when none.
# Shared by the evidence.txt upload check and the MaxQuant matrix build, so
# the upload rejects exactly the files the build could not use.
mq_resolve_quant_col <- function(raw_file_name, raw_to_exp, cn_norm) {
  exp_name <- if (!is.null(raw_to_exp) && raw_file_name %in% names(raw_to_exp))
    raw_to_exp[[raw_file_name]] else raw_file_name
  exp_safe <- make.names(exp_name)
  candidates <- c(
    paste0("Intensity.", exp_safe),
    paste0("LFQ.intensity.", exp_safe)
  )
  for (cand in candidates) {
    idx <- which(cn_norm == cand)
    if (length(idx) > 0) return(cn_norm[idx[1]])
  }
  # Partial suffix match ("Intensity" is case-sensitive, so it never matches "LFQ.intensity.")
  matches <- grep(paste0("Intensity\\.", exp_safe, "$"), cn_norm, value = TRUE)
  if (length(matches) > 0) return(matches[1])
  matches <- grep(paste0("LFQ.intensity\\.", exp_safe, "$"), cn_norm, value = TRUE)
  if (length(matches) > 0) return(matches[1])
  NA_character_
}

# =============================================================================
# BUILD_PRODA_PROTEIN_MATRIX
# Routes to the correct matrix-building strategy based on input mode / format.
# Returns: list(protein_matrix, acc_to_gene, n_peptides_input, n_proteins_output)
#          OR list(error = "message") on failure.
#
# Args:
#   rv            -" Shiny reactiveValues (or plain list)  -" read-only
#   cond_groups   -" named list: group_name → c(col1, col2, ...)
#   logger        -" logging function
# =============================================================================
build_proda_protein_matrix <- function(rv, cond_groups,
                                       shared_peptides = "unique_only",
                                       logger = message) {

  # Hard error, not a fallback: an unrecognised value would otherwise be treated
  # silently as one of the accepted settings and misreported in provenance.
  # Consumed ONLY by the peptide-aggregating route (B). The protein-direct
  # routes (A, MQ) ignore it and record NA_character_ instead.
  if (!(length(shared_peptides) == 1L && !is.na(shared_peptides) &&
        shared_peptides %in% c("unique_only", "keep_all", "razor")))
    stop("build_proda_protein_matrix: shared_peptides must be one of ",
         "\"unique_only\", \"keep_all\" or \"razor\".")

  input_mode   <- rv$lf_input_mode   %||% "peptide"
  input_format <- rv$lf_input_format %||% "pd_nested"
  data_types   <- rv$wizard_data_types %||% character(0)

  # Flatten all column names in the right order and build sample/condition mapping
  group_names  <- names(cond_groups)
  all_cols     <- unlist(cond_groups, use.names = FALSE)
  sample_names <- character(0)
  design_rows  <- list()
  idx <- 0L
  for (gi in seq_along(cond_groups)) {
    grp_name <- group_names[gi]
    for (ri in seq_along(cond_groups[[gi]])) {
      idx <- idx + 1L
      col_nm <- cond_groups[[gi]][ri]
      sample_nm <- paste0(grp_name, "_Rep", ri)
      sample_names <- c(sample_names, sample_nm)
      design_rows[[idx]] <- data.frame(
        sample       = col_nm,
        sample_name  = sample_nm,
        condition    = grp_name,
        bioreplicate = paste0(grp_name, "_", ri),
        stringsAsFactors = FALSE
      )
    }
  }
  design_df <- do.call(rbind, design_rows)

  # ── Route A: PD nested + protein data type available ──
  # Extract protein-level abundances directly (no peptide aggregation needed)
  if (identical(input_mode, "peptide") && identical(input_format, "pd_nested") &&
      "protein" %in% data_types) {

    logger("ProDA matrix: PD nested  -  extracting protein-level abundances")
    pd_raw   <- rv$pd_raw_lines
    pd_idx   <- rv$pd_protein_line_idx
    pd_hdr   <- if (length(pd_idx) > 0) pd_raw[1] else NULL

    if (is.null(pd_raw) || length(pd_idx) < 1)
      return(list(error = "ProDA matrix: PD protein data not available."))

    block   <- paste(c(pd_raw[1], pd_raw[pd_idx]), collapse = "\n")
    prot_df <- tryCatch(
      read.csv(text = block, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, na.strings = c("", "NA")),
      error = function(e) NULL
    )
    if (is.null(prot_df) || nrow(prot_df) < 1)
      return(list(error = "ProDA matrix: PD protein CSV parse failed."))

    # Gene names
    accessions <- as.character(prot_df[["Accession"]])
    if (!isTRUE(rv$pd_has_gn) && !is.null(rv$pd_gene_col) &&
        rv$pd_gene_col %in% colnames(prot_df)) {
      gene_names <- as.character(prot_df[[rv$pd_gene_col]])
    } else {
      descs      <- prot_df[["Description"]]
      gene_names <- stringr::str_replace(stringr::str_extract(descs, "GN=([^ ]+)"), "^GN=", "")
    }
    gene_names[is.na(gene_names)] <- ""
    acc_to_gene_full <- setNames(gene_names, accessions)
    # Flag: TRUE when gene names came from the data (not just accessions)
    .genes_from_data <- any(nzchar(gene_names))

    # Find abundance columns matching cond_groups
    prot_cn      <- colnames(prot_df)
    abund_cols   <- grep("^Abundance", prot_cn, ignore.case = TRUE, value = TRUE)
    abund_cols   <- abund_cols[!grepl(
      "ratio|p-value|variability|normalized|normalised|count|grouped|scaled",
      abund_cols, ignore.case = TRUE)]

    # Map requested columns to actual protein-level columns
    mat_cols  <- character(length(all_cols))
    for (ci in seq_along(all_cols)) {
      col_nm <- all_cols[ci]
      if (col_nm %in% prot_cn) {
        mat_cols[ci] <- col_nm
      } else {
        f_id      <- stringr::str_extract(col_nm, "F\\d+")
        if (!is.na(f_id)) {
          match_c <- grep(paste0("\\b", f_id, "\\b"), abund_cols, value = TRUE)
          if (length(match_c) > 1)
            return(list(error = paste0(
              "ProDA matrix: ambiguous abundance column for '", col_nm,
              "' - ", f_id, " matches ", length(match_c), " columns: ",
              paste(match_c, collapse = "; "))))
          if (length(match_c) > 0) mat_cols[ci] <- match_c[1]
        }
      }
    }
    unmapped <- which(mat_cols == "")
    if (length(unmapped) > 0)
      return(list(error = paste0("ProDA matrix: could not map ",
                                 length(unmapped), " PD column(s).")))

    int_mat <- prot_df[, mat_cols, drop = FALSE]
    int_mat <- apply(int_mat, 2, function(x) suppressWarnings(as.numeric(x)))
    int_mat[int_mat <= 0] <- NA
    int_mat <- log2(int_mat)
    rownames(int_mat) <- accessions
    colnames(int_mat) <- design_df$sample_name

    # Filter: require at least 1 non-NA per protein
    keep <- apply(int_mat, 1, function(r) sum(!is.na(r)) >= 1)
    int_mat <- int_mat[keep, , drop = FALSE]

    acc_to_gene <- acc_to_gene_full[rownames(int_mat)]
    logger(paste0("ProDA matrix (PD prot): ", nrow(int_mat), " proteins \u00d7 ",
                  ncol(int_mat), " samples"))
    return(list(
      protein_matrix    = int_mat,
      acc_to_gene       = acc_to_gene,
      genes_from_data   = .genes_from_data,
      n_peptides_input  = NA_integer_,
      n_proteins_output = nrow(int_mat),
      design_df         = design_df,
      input_route       = "pd_protein_direct",  # no peptide aggregation performed
      # No peptides are aggregated here, so no shared-peptide filter can run.
      shared_peptide_handling = NA_character_
    ))
  }

  # ── Route B: PD nested + peptide data only (no protein type) ──
  if (identical(input_mode, "peptide") && identical(input_format, "pd_nested")) {
    logger("ProDA matrix: PD nested + peptide  -  parsing peptide-level data for QFeatures")
    pd_raw <- rv$pd_raw_lines
    if (is.null(pd_raw))
      return(list(error = "ProDA matrix: pd_raw_lines not available."))

    pep_hdr_idx  <- grep("^,.*Annotated Sequence", pd_raw)
    all_nested   <- grep("^,", pd_raw)
    pep_data_idx <- setdiff(all_nested, pep_hdr_idx)

    if (length(pep_data_idx) < 10)
      return(list(error = "ProDA matrix: no peptide-level data found in CSV."))

    pep_block <- paste(c(pd_raw[pep_hdr_idx[1]], pd_raw[pep_data_idx]), collapse = "\n")
    pep_df <- tryCatch(
      read.csv(text = pep_block, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, na.strings = c("", "NA")),
      error = function(e) NULL
    )
    if (is.null(pep_df) || nrow(pep_df) < 10)
      return(list(error = "ProDA matrix: peptide CSV parse failed."))

    pep_cn     <- colnames(pep_df)
    master_col <- grep("Master Protein Accessions", pep_cn, value = TRUE)[1]
    seq_col    <- grep("Annotated Sequence",        pep_cn, value = TRUE)[1]
    if (is.na(master_col) || is.na(seq_col))
      return(list(error = "ProDA matrix: required peptide columns missing."))

    # Find matching abundance columns
    pep_abund  <- grep("^Abundance", pep_cn, ignore.case = TRUE, value = TRUE)
    pep_abund  <- pep_abund[!grepl(
      "ratio|p-value|variability|normalized|normalised|count|grouped|scaled",
      pep_abund, ignore.case = TRUE)]

    col_map <- list()
    for (pcol in all_cols) {
      f_id <- stringr::str_extract(pcol, "F\\d+")
      if (!is.na(f_id)) {
        m <- grep(f_id, pep_abund, value = TRUE)
        if (length(m) >= 1) col_map[[pcol]] <- m[1]
      }
      if (is.null(col_map[[pcol]]) && pcol %in% pep_cn) col_map[[pcol]] <- pcol
    }
    int_cols <- unlist(lapply(all_cols, function(x) col_map[[x]]))
    if (any(is.na(int_cols)))
      return(list(error = "ProDA matrix: could not map all peptide abundance columns."))

    # -- Shared-peptide handling (Route B) --
    # A PD peptide export carries no unique/razor verdict of its own, so the
    # semicolon-delimited master accession list is the evidence used. Same
    # helper as the TMT peptide route, deliberately not duplicated.
    .shared_b   <- tmt_flag_shared_peptides(pep_df, master_col)
    n_pep_total_b  <- nrow(pep_df)
    n_pep_shared_b <- sum(.shared_b$shared)
    n_pep_removed_b <- 0L
    n_razor_assigned_b <- 0L
    n_razor_dropped_b  <- 0L
    .rz_prot_b <- character(0)
    .rz_asg_b  <- logical(0)
    if (identical(shared_peptides, "unique_only")) {
      pep_df <- pep_df[!.shared_b$shared, , drop = FALSE]
      n_pep_removed_b <- n_pep_shared_b
      logger(paste0("ProDA matrix: shared-peptide filter (PD nested peptide): ",
                    n_pep_total_b, " -> ", nrow(pep_df),
                    " rows (peptides matching >1 master protein removed)"))
      if (nrow(pep_df) == 0L)
        return(list(error = paste0(
          "ProDA matrix: no peptides left after the shared-peptide filter - every ",
          "peptide matched more than one master protein. Re-run with ",
          "'Keep all peptides' to include them.")))
    } else if (identical(shared_peptides, "razor")) {
      # Parsimony assignment, MaxQuant's rule: each shared peptide goes to the
      # candidate protein supported by the most unique peptide SEQUENCES. Same
      # helper as the TMT peptide route, deliberately not duplicated. This route
      # builds one complete table, so this already operates on the pooled data.
      .rz_b <- tmt_razor_assign_shared(as.character(pep_df[[seq_col]]),
                                       .shared_b$candidates)
      n_razor_assigned_b <- sum(.rz_b$assigned)
      n_razor_dropped_b  <- sum(.rz_b$dropped)
      .keep_b    <- !.rz_b$dropped
      pep_df     <- pep_df[.keep_b, , drop = FALSE]
      .rz_prot_b <- .rz_b$protein[.keep_b]
      .rz_asg_b  <- .rz_b$assigned[.keep_b]
      logger(paste0("ProDA matrix: parsimony (razor) assignment (PD nested ",
                    "peptide): ", n_razor_assigned_b, " shared row(s) assigned by ",
                    "unique-peptide support, ", n_razor_dropped_b,
                    " dropped for lack of unique support; ", nrow(pep_df),
                    " rows remain"))
      if (nrow(pep_df) == 0L)
        return(list(error = paste0(
          "ProDA matrix: no peptides left after parsimony (razor) assignment - no ",
          "candidate protein had any unique-peptide support. Re-run with ",
          "'Keep all peptides' to include them.")))
    }

    pep_df$CleanProtein <- sapply(strsplit(pep_df[[master_col]], ";"), function(x) trimws(x[1]))

    # Razor overwrites the reassigned rows ONLY; unique rows keep the value the
    # first-accession expression above gave them.
    if (identical(shared_peptides, "razor") && any(.rz_asg_b))
      pep_df$CleanProtein[.rz_asg_b] <- .rz_prot_b[.rz_asg_b]

    qf_result <- qf_summarize_to_protein_matrix(
      input_df                = pep_df,
      input_type              = "flat_peptide",
      protein_col             = "CleanProtein",
      peptide_col             = seq_col,
      intensity_cols          = int_cols,
      sample_names            = design_df$sample_name,
      # Peptide-level centering (QFeatures vignette order). On this label-free
      # route it does not decide the final per-sample alignment: robustSummary
      # carries a per-sample shift through to the protein value unchanged, and
      # lf_core_summarization_proda() then normalizes the protein matrix (global
      # median or invariant-protein panel), which replaces it. The data are
      # therefore normalized once, at the protein level; the methods text
      # describes only that step. Checked on synthetic data 2026-10-05: final
      # protein values with and without this centering differ by one common
      # constant (~0.002 log2 protein-to-protein, from proDA's NA handling).
      # When the input was normalized upstream the protein step is skipped, so this
      # centering is skipped too: otherwise it would be the one normalization left in
      # effect. (PD input is never normalized upstream, so in practice this stays on.)
      normalize_peptide_level = lf_apply_normalization_effective(rv),
      logger                  = logger
    )
    if (!is.null(qf_result$error)) return(qf_result)
    logger(paste0("ProDA matrix: peptide-level centering above is provisional; the protein-level ",
                  "normalization applied at summarization sets the final per-sample alignment ",
                  "(one normalization in effect)."))
    qf_result$design_df <- design_df
    qf_result$shared_peptide_handling     <- shared_peptides
    qf_result$n_peptides_total            <- n_pep_total_b
    qf_result$n_peptides_shared           <- n_pep_shared_b
    qf_result$n_peptides_removed_shared   <- n_pep_removed_b
    qf_result$shared_detection_method     <- "master_accession_count"
    qf_result$shared_detection_applicable <- TRUE
    qf_result$n_peptides_razor_assigned             <- n_razor_assigned_b
    qf_result$n_peptides_dropped_no_unique_evidence <- n_razor_dropped_b
    qf_result$razor_basis <- if (identical(shared_peptides, "razor"))
                               "unique_peptide_count" else NA_character_
    return(qf_result)
  }

  # ── Route MQ: MaxQuant proteinGroups.txt raw Intensity direct (LFQ intensity only as fallback) ──────────────
  # When input is MaxQuant, read LFQ Intensity (or Intensity) columns directly
  # from proteinGroups.txt. cond_groups maps Raw.file names to groups; we use
  # rv$mq_rawfile_to_experiment to resolve the LFQ column name per run, then
  # build a protein x sample log2 matrix without going through MSstats.
  if (identical(rv$wizard_search_program %||% "", "maxquant")) {
    logger("ProDA matrix: MaxQuant -- reading proteinGroups.txt LFQ columns")

    pg_path <- rv$mq_pg_path %||% rv$mq_proteingroups_path
    if (is.null(pg_path) || !file.exists(pg_path))
      return(list(error = "ProDA matrix (MQ): proteinGroups.txt not available."))

    pg_raw <- tryCatch(
      data.table::fread(pg_path, data.table = FALSE),
      error = function(e) NULL
    )
    if (is.null(pg_raw) || nrow(pg_raw) < 1)
      return(list(error = "ProDA matrix (MQ): failed to read proteinGroups.txt."))

    cn_orig <- colnames(pg_raw)
    cn_norm <- make.names(cn_orig, unique = TRUE)

    # Identify protein ID column. proteinGroups.txt intensities are already
    # quantified by MaxQuant under its own razor/unique setting, so there is no
    # protein-ID choice to make here: the leading ID of each group is used.
    # (rv$mq_protein_id_type from older sessions/bundles is ignored.)
    prot_id_col <- intersect(c("Protein.IDs", "Majority.protein.IDs"), cn_norm)[1]
    if (is.na(prot_id_col))
      return(list(error = "ProDA matrix (MQ): proteinGroups.txt has no 'Protein IDs' column."))

    # Resolve Raw.file -> LFQ column name via Experiment mapping
    raw_to_exp <- rv$mq_rawfile_to_experiment
    # Raw "Intensity" columns are preferred (2026-10-05, Erez): they are not
    # normalized by MaxQuant, so OSMSS normalizes them once (global median or
    # invariant-protein panel). "LFQ intensity" (MaxLFQ) values are already
    # normalized by MaxQuant; normalizing them again in OSMSS would be a second
    # normalization, so they are used only when a file has no raw Intensity
    # columns. Benchmark (12 OpDEA DDA datasets): raw Intensity + OSMSS median
    # beat MaxLFQ with or without OSMSS normalization (dda_norm_check/README.md).
    .resolve_lfq_col <- function(raw_file_name)
      mq_resolve_quant_col(raw_file_name, raw_to_exp, cn_norm)

    # Build column mapping: design_df$sample (Raw.file) -> LFQ column name
    lfq_cols <- vapply(design_df$sample, .resolve_lfq_col,
                       character(1L), USE.NAMES = TRUE)
    failed <- which(is.na(lfq_cols) | !nzchar(lfq_cols))
    if (length(failed) > 0)
      return(list(error = paste0(
        "ProDA matrix (MQ): could not find LFQ column(s) for: ",
        paste(design_df$sample[failed], collapse = ", "))))

    # Record which value type each sample was read from. The resolver prefers
    # raw Intensity and falls back to LFQ intensity only when no raw column
    # exists; a mix across samples is not comparable.
    quant_types <- ifelse(grepl("^LFQ\\.intensity\\.", lfq_cols), "LFQ intensity", "Intensity")
    quant_type  <- if (length(unique(quant_types)) == 1L) quant_types[1] else "mixed"
    if (identical(quant_type, "mixed"))
      return(list(error = paste0(
        "ProDA matrix (MQ): some samples have raw Intensity columns and others only ",
        "LFQ intensity (", paste(design_df$sample[quant_types == "LFQ intensity"], collapse = ", "),
        "). Re-export proteinGroups.txt from MaxQuant with Intensity columns for all experiments.")))
    if (identical(quant_type, "LFQ intensity"))
      logger(paste0("!!! ProDA matrix (MQ): no raw 'Intensity' columns found - using 'LFQ intensity' ",
                    "(already normalized by MaxQuant); OSMSS normalization is applied on top"))
    else
      logger("ProDA matrix (MQ): using raw 'Intensity' columns (normalized once, by OSMSS)")

    # Filter: remove reverse hits and only-identified-by-site groups; remove
    # potential contaminants only when the user opted in (rv$mq_remove_contaminants,
    # default FALSE). MaxQuant flags a group '+' when it contains ANY contaminant
    # entry, so the flag also catches e.g. keratins or UPS1 albumin led by a
    # sample protein. Counts are sequential (a row flagged twice is counted under
    # the first flag applied).
    remove_contam <- isTRUE(rv$mq_remove_contaminants)
    contam_idx    <- which(cn_norm == "Potential.contaminant")
    n_flagged_contam <- if (length(contam_idx) > 0)
      sum(pg_raw[[contam_idx[1]]] %in% "+") else NA_integer_
    keep_rows <- rep(TRUE, nrow(pg_raw))
    all_flags     <- c("Potential.contaminant", "Reverse", "Only.identified.by.site")
    flag_patterns <- if (remove_contam) all_flags else setdiff(all_flags, "Potential.contaminant")
    n_removed_flags <- setNames(integer(length(all_flags)), all_flags)
    for (fp in flag_patterns) {
      idx <- which(cn_norm == fp)
      if (length(idx) > 0) {
        flag_ok <- is.na(pg_raw[[idx[1]]]) | pg_raw[[idx[1]]] != "+"
        n_removed_flags[[fp]] <- sum(keep_rows & !flag_ok)
        keep_rows <- keep_rows & flag_ok
      }
    }
    pg_filt <- pg_raw[keep_rows, , drop = FALSE]
    colnames(pg_filt) <- cn_norm   # normalize so character-name lookups work
    logger(paste0("ProDA matrix (MQ): ", nrow(pg_filt), " proteins after filters (",
                  if (remove_contam) paste0(n_removed_flags[["Potential.contaminant"]], " contaminant, ")
                  else paste0("contaminants kept [", n_flagged_contam, " flagged], "),
                  n_removed_flags[["Reverse"]], " reverse, ",
                  n_removed_flags[["Only.identified.by.site"]], " only-by-site removed)"))

    # Build log2 matrix
    int_mat <- matrix(NA_real_, nrow = nrow(pg_filt), ncol = nrow(design_df))
    colnames(int_mat) <- design_df$sample_name
    for (ci in seq_len(nrow(design_df))) {
      col_nm <- lfq_cols[[ci]]
      if (!is.na(col_nm) && col_nm %in% cn_norm) {
        vals <- suppressWarnings(as.numeric(pg_filt[[col_nm]]))
        vals[vals <= 0] <- NA_real_
        int_mat[, ci] <- log2(vals)
      }
    }

    # Protein accession: leading ID of the group. When MaxQuant ran without an
    # identifier parse rule the IDs are full UniProt headers
    # ("sp|P02533|K1C14_HUMAN"); reduce those to the accession ("P02533").
    accessions <- sub(";.*", "", trimws(as.character(pg_filt[[prot_id_col]])))
    accessions <- sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", accessions)

    # Gene names: 'Gene names' column when MaxQuant parsed it, otherwise the GN=
    # tag in 'Fasta headers' (same extraction as the PD Description route).
    gene_names  <- rep(NA_character_, length(accessions))
    gene_source <- "none"
    gene_col_idx <- which(cn_norm == "Gene.names")
    if (length(gene_col_idx) > 0) {
      gene_names  <- sub(";.*", "", trimws(as.character(pg_filt[[gene_col_idx[1]]])))
      gene_source <- "gene_names_column"
    }
    fasta_idx <- which(cn_norm == "Fasta.headers")
    if (length(fasta_idx) > 0) {
      need <- is.na(gene_names) | !nzchar(gene_names)
      if (any(need)) {
        gn_fasta <- stringr::str_replace(
          stringr::str_extract(as.character(pg_filt[[fasta_idx[1]]][need]), "GN=([^ ;]+)"),
          "^GN=", "")
        gene_names[need] <- gn_fasta
        if (identical(gene_source, "none") && any(!is.na(gn_fasta)))
          gene_source <- "fasta_headers_gn"
      }
    }
    gene_names[is.na(gene_names) | !nzchar(gene_names)] <-
      accessions[is.na(gene_names) | !nzchar(gene_names)]

    if (anyDuplicated(accessions))
      return(list(error = paste0(
        "ProDA matrix (MQ): leading protein ID is not unique across protein groups (e.g. '",
        accessions[anyDuplicated(accessions)], "').")))

    rownames(int_mat) <- accessions
    acc_to_gene_mq <- setNames(gene_names, accessions)

    # Filter: require at least 1 non-NA observation
    keep2 <- apply(int_mat, 1, function(r) sum(!is.na(r)) >= 1L)
    int_mat <- int_mat[keep2, , drop = FALSE]
    acc_to_gene_mq <- acc_to_gene_mq[rownames(int_mat)]

    logger(paste0("ProDA matrix (MQ): ", nrow(int_mat), " proteins x ",
                  ncol(int_mat), " samples"))
    genes_from_data_mq <- any(acc_to_gene_mq != names(acc_to_gene_mq), na.rm = TRUE)
    logger(paste0("ProDA matrix (MQ): gene names from ", gene_source, " for ",
                  sum(acc_to_gene_mq != names(acc_to_gene_mq), na.rm = TRUE), " / ",
                  length(acc_to_gene_mq), " proteins"))
    return(list(
      protein_matrix    = int_mat,
      acc_to_gene       = acc_to_gene_mq,
      genes_from_data   = genes_from_data_mq,
      n_peptides_input  = NA_integer_,
      n_proteins_output = nrow(int_mat),
      design_df         = design_df,
      input_route       = "maxquant_lfq_direct",
      # Carried onto the cache for provenance / methods text.
      mq_stats          = list(
        quant_type             = quant_type,
        n_groups_input         = nrow(pg_raw),
        n_removed_contaminant  = n_removed_flags[["Potential.contaminant"]],
        contaminants_removed   = remove_contam,
        n_flagged_contaminant  = n_flagged_contam,
        n_removed_reverse      = n_removed_flags[["Reverse"]],
        n_removed_only_by_site = n_removed_flags[["Only.identified.by.site"]],
        n_removed_no_values    = sum(!keep2),
        n_proteins_output      = nrow(int_mat),
        gene_source            = gene_source,
        # match between runs, detected at evidence.txt upload (server.R); NA when unknown
        mbr_on                 = rv$mq_mbr_on %||% NA,
        n_match_features       = rv$mq_n_match_features %||% NA_integer_,
        n_evidence_features    = rv$mq_n_evidence_features %||% NA_integer_),
      # No peptides are aggregated here, so no shared-peptide filter can run.
      shared_peptide_handling = NA_character_
    ))
  }

  list(error = paste0("ProDA matrix: unsupported input mode/format combination: ",
                      input_mode, "/", input_format))
}

# === INVARIANT-PROTEIN MEDIAN NORMALIZATION ===================================
# Reference-panel selection for ip_sub = "invariant". Candidate reference
# proteins must be detected in every run (all conditions, including negative
# controls), be highly abundant, have a low within-condition SD, and a low
# between-condition amplitude. The amplitude is evaluated on condition means
# centred on the bulk of completely detected proteins, not on the raw values.
#
# Why centre: when a large share of the proteins shifts in one direction
# (knockouts, strong perturbations, bait enrichment), loading and MaxLFQ's own
# run normalization move the truly unchanged proteins away from a zero
# between-condition difference. A raw amplitude gate then rejects every
# unchanged protein and can admit only changed proteins whose shift happens to
# cancel the offset (benchmark: Benchmark/semisynth_KO, R/08_ebn_relative_test.R).
# Centring removes any common offset before invariance is judged.
#
# Centre: per condition, the deviation of each complete protein's condition
# mean from its mean over conditions; its density mode (stats::density
# defaults, n = 4096). Unchanged proteins form the tallest, sharpest peak even
# when they are a minority, so the mode stays on them where the median does not.
# Fallback: if a second density peak at least .INVARIANT_PEAK_RATIO as tall as
# the main one lies more than two bandwidths away (the mode is ambiguous), the
# candidate peak nearest zero is used, i.e. the one implying the smallest
# loading correction (in AP-MS, the background; in a bimodal knockout, the
# unchanged peak unless the loading offset exceeds half the knockout effect).
# The median is not used: with two comparable peaks it falls between them.
# Mode and median centres further apart than .INVARIANT_CENTRE_GAP on the
# between-condition scale are reported (the mode is kept). A panel below
# .INVARIANT_MIN_REFS falls back to global median normalization.
# ==============================================================================
.invariant_centre <- function(dev, cond, peak_ratio_max = .INVARIANT_PEAK_RATIO) {
  conds  <- unique(as.character(cond))
  mode_c <- setNames(rep(NA_real_, length(conds)), conds)
  med_c  <- mode_c
  ratio  <- mode_c
  near_c <- mode_c    # candidate peak nearest zero (main or a competing one)
  for (k in conds) {
    x <- dev[cond == k]
    x <- x[is.finite(x)]
    if (length(x) == 0L) next
    med_c[k] <- stats::median(x)
    if (length(x) >= 10L && stats::sd(x) > 0) {
      d      <- stats::density(x, n = 4096)
      i_main <- which.max(d$y)
      mode_c[k] <- d$x[i_main]
      pk  <- which(diff(sign(diff(d$y))) == -2L) + 1L     # local maxima
      far <- pk[abs(d$x[pk] - d$x[i_main]) > 2 * d$bw]
      ratio[k] <- if (length(far)) max(d$y[far]) / d$y[i_main] else 0
      cand <- c(i_main, far[d$y[far] >= peak_ratio_max * d$y[i_main]])
      near_c[k] <- d$x[cand[which.min(abs(d$x[cand]))]]
    }
  }
  # A condition without a usable density (too few values) keeps its median.
  near_c[is.na(near_c)] <- med_c[is.na(near_c)]
  ambiguous <- any(is.na(mode_c)) || any(ratio >= peak_ratio_max, na.rm = TRUE)
  # Gap between the two centres on the between-condition scale: the largest
  # disagreement over condition pairs of (centre_a - centre_b).
  .pair <- function(v) outer(v, v, "-")
  gap <- if (any(is.na(mode_c)) || any(is.na(med_c))) NA_real_ else
    max(abs(.pair(mode_c) - .pair(med_c)))
  list(centre        = if (ambiguous) near_c else mode_c,
       method        = if (ambiguous) "peak_nearest_zero" else "mode",
       mode_centre   = mode_c,
       median_centre = med_c,
       peak_ratio    = if (all(is.na(ratio))) NA_real_ else max(ratio, na.rm = TRUE),
       gap           = gap,
       ambiguous     = ambiguous)
}

select_invariant_protein_panel <- function(msstats_input,
                                           intensity_is_log2 = FALSE,
                                           logger = message) {
  valid <- if (isTRUE(intensity_is_log2)) {
    msstats_input[is.finite(msstats_input$Intensity), ]
  } else {
    msstats_input[!is.na(msstats_input$Intensity) &
                  msstats_input$Intensity > 0, ]
  }
  valid$log2Int <- if (isTRUE(intensity_is_log2)) {
    valid$Intensity
  } else log2(valid$Intensity)

  # Protein x Run medians
  prot_run_med <- aggregate(log2Int ~ ProteinName + Run + Condition,
                            data = valid, FUN = median, na.rm = TRUE)

  all_runs <- unique(msstats_input$Run)
  n_runs   <- length(all_runs)

  # --- Gate 1: Completeness -- protein must appear in every run ---
  prot_run_count <- tapply(prot_run_med$Run, prot_run_med$ProteinName,
                           function(x) length(unique(x)))
  complete_prots <- names(prot_run_count)[prot_run_count == n_runs]

  if (length(complete_prots) < .INVARIANT_MIN_REFS) {
    logger(paste0("Invariant-protein normalization: only ", length(complete_prots),
                  " proteins detected in all ", n_runs,
                  " runs (need ≥", .INVARIANT_MIN_REFS, "). Cannot build reference panel."))
    return(NULL)
  }
  if (length(complete_prots) < .INVARIANT_MIN_REFS_IDEAL) {
    logger(paste0("Invariant-protein normalization: warning — only ",
                  length(complete_prots), " complete proteins (ideal ≥",
                  .INVARIANT_MIN_REFS_IDEAL, "); reference panel may have low confidence"))
  }

  prot_run_med <- prot_run_med[prot_run_med$ProteinName %in% complete_prots, ]

  # --- Gate 2: High Abundance -- global median intensity per protein ---
  prot_global_med <- aggregate(log2Int ~ ProteinName, data = prot_run_med,
                               FUN = median, na.rm = TRUE)
  names(prot_global_med)[2] <- "GlobalMedian"

  # --- Gate 3: Relative amplitude (condition means centred on the bulk) ---
  prot_cond_mean <- aggregate(log2Int ~ ProteinName + Condition,
                              data = prot_run_med, FUN = mean, na.rm = TRUE)
  .pm       <- ave(prot_cond_mean$log2Int, prot_cond_mean$ProteinName)
  .cond_chr <- as.character(prot_cond_mean$Condition)
  .ctr      <- .invariant_centre(prot_cond_mean$log2Int - .pm, .cond_chr)
  prot_cond_mean$log2Int <- prot_cond_mean$log2Int - .ctr$centre[.cond_chr]
  .gap_flag <- isTRUE(.ctr$gap > .INVARIANT_CENTRE_GAP)
  logger(paste0("Invariant-protein normalization: condition means centred on the ",
                if (identical(.ctr$method, "mode")) "density mode" else "density peak nearest zero",
                " of ", length(complete_prots), " complete proteins (",
                paste0(names(.ctr$centre), " ", sprintf("%+.3f", .ctr$centre),
                       collapse = ", "), " log2)",
                if (isTRUE(.ctr$ambiguous))
                  "; density mode ambiguous (competing peak), peak nearest zero used" else "",
                if (!is.na(.ctr$gap))
                  paste0("; mode-median gap ", sprintf("%.3f", .ctr$gap), " log2") else "",
                if (.gap_flag) " — exceeds the reporting threshold" else ""))

  # Amplitude = max(centred condition means) - min(centred condition means)
  prot_amplitude <- aggregate(log2Int ~ ProteinName, data = prot_cond_mean,
                              FUN = function(x) max(x) - min(x))
  names(prot_amplitude)[2] <- "Amplitude"

  # --- Gate 4: Universal Within-Group Variance Filter ---
  conditions <- unique(prot_run_med$Condition)
  sd_list <- list()
  for (cond in conditions) {
    cond_data <- prot_run_med[prot_run_med$Condition == cond, ]
    cond_sd   <- aggregate(log2Int ~ ProteinName, data = cond_data,
                           FUN = sd, na.rm = TRUE)
    names(cond_sd)[2] <- paste0("SDlog2_", cond)
    sd_list[[cond]] <- cond_sd
  }
  sd_merged <- sd_list[[1]]
  if (length(sd_list) > 1) {
    for (i in 2:length(sd_list)) {
      sd_merged <- merge(sd_merged, sd_list[[i]], by = "ProteinName", all = TRUE)
    }
  }
  sd_cols <- grep("^SDlog2_", colnames(sd_merged), value = TRUE)
  sd_merged$MaxSDlog2 <- apply(sd_merged[, sd_cols, drop = FALSE], 1,
                               function(x) max(x, na.rm = TRUE))
  .n_nonfinite_sd <- sum(!is.finite(sd_merged$MaxSDlog2))
  if (.n_nonfinite_sd > 0L) {
    sd_merged$MaxSDlog2[!is.finite(sd_merged$MaxSDlog2)] <- Inf
    logger(paste0("Invariant-protein normalization: ", .n_nonfinite_sd,
                  " protein(s) had no usable within-condition SD and were",
                  " excluded from the reference panel"))
  }

  ref_stats <- merge(prot_global_med, prot_amplitude, by = "ProteinName")
  ref_stats <- merge(ref_stats, sd_merged[, c("ProteinName", "MaxSDlog2")],
                     by = "ProteinName")

  abundance_cutoff_top25 <- quantile(ref_stats$GlobalMedian, 0.75, na.rm = TRUE)
  abundance_cutoff_top50 <- quantile(ref_stats$GlobalMedian, 0.50, na.rm = TRUE)

  # --- Threshold relaxation schedule ---
  relaxation_steps <- list(
    list(sd = 0.15, amp = 0.1, abund = abundance_cutoff_top25),
    list(sd = 0.20, amp = 0.1, abund = abundance_cutoff_top25),
    list(sd = 0.25, amp = 0.1, abund = abundance_cutoff_top25),
    list(sd = 0.30, amp = 0.1, abund = abundance_cutoff_top25),
    list(sd = 0.30, amp = 0.2, abund = abundance_cutoff_top25),
    list(sd = 0.30, amp = 0.3, abund = abundance_cutoff_top25),
    list(sd = 0.30, amp = 0.5, abund = abundance_cutoff_top25),
    list(sd = 0.30, amp = 0.5, abund = abundance_cutoff_top50)
  )

  selected  <- character(0)
  step_used <- 0
  best_n    <- 0L
  for (i in seq_along(relaxation_steps)) {
    s <- relaxation_steps[[i]]
    mask <- ref_stats$MaxSDlog2 <= s$sd &
            ref_stats$Amplitude <= s$amp &
            ref_stats$GlobalMedian >= s$abund
    candidates <- ref_stats$ProteinName[mask]
    if (length(candidates) > best_n) best_n <- length(candidates)
    if (length(candidates) >= .INVARIANT_MIN_REFS) {
      selected  <- as.character(candidates)
      step_used <- i
      break
    }
  }

  if (length(selected) < .INVARIANT_MIN_REFS) {
    logger(paste0("Invariant-protein normalization: after all relaxation steps, best",
                  " candidate count was ", best_n, " (need ≥", .INVARIANT_MIN_REFS,
                  "). No invariant reference panel meets the criteria on this",
                  " dataset; the pipeline will use global median normalization",
                  " instead. This is an expected and acceptable outcome, not an error."))
    return(NULL)
  }
  if (length(selected) < .INVARIANT_MIN_REFS_IDEAL) {
    logger(paste0("Invariant-protein normalization: warning — only ",
                  length(selected), " reference proteins selected (ideal ≥",
                  .INVARIANT_MIN_REFS_IDEAL, "); normalization confidence is reduced"))
  }
  if (step_used > 1) {
    s <- relaxation_steps[[step_used]]
    logger(paste0("Invariant-protein normalization: thresholds relaxed (step ", step_used,
                  "): SD≤", s$sd, " log2, amplitude≤", s$amp,
                  ", abundance≥", round(s$abund, 2)))
  }
  logger(paste0("Invariant-protein normalization: selected ", length(selected),
                " reference proteins (step ", step_used, " of ",
                length(relaxation_steps), ") across ",
                length(conditions), " conditions"))

  .amp_by_step <- vapply(relaxation_steps,
                         function(s) as.numeric(s$amp %||% NA_real_), numeric(1))
  .strict_max_step <- if (length(.amp_by_step) == 0L ||
                          all(is.na(.amp_by_step))) {
    1L
  } else {
    max(which(.amp_by_step <= .amp_by_step[[1]]))
  }

  # Ladder attributes keep their historical ebn_* names (read by provenance,
  # methods text and saved caches), plus the centring diagnostics.
  attr(selected, "ebn_step")            <- step_used
  attr(selected, "ebn_n_steps")         <- length(relaxation_steps)
  attr(selected, "ebn_thresholds")      <- relaxation_steps[[step_used]]
  attr(selected, "ebn_n_selected")      <- length(selected)
  attr(selected, "ebn_strict")          <- step_used <= .strict_max_step
  attr(selected, "ebn_strict_max_step") <- .strict_max_step
  attr(selected, "invariant_centre_method")    <- .ctr$method
  attr(selected, "invariant_centre")           <- .ctr$centre
  attr(selected, "invariant_mode_centre")      <- .ctr$mode_centre
  attr(selected, "invariant_median_centre")    <- .ctr$median_centre
  attr(selected, "invariant_peak_ratio")       <- .ctr$peak_ratio
  attr(selected, "invariant_centre_gap")       <- .ctr$gap
  attr(selected, "invariant_centre_disagree")  <- .gap_flag
  attr(selected, "invariant_mode_ambiguous")   <- isTRUE(.ctr$ambiguous)

  return(selected)
}

# ProDA-based summarization for the invariant-protein and global-median routes.
# Takes a pre-built protein matrix (proteins × samples, log2, NAs OK) and
# applies invariant-protein (ip_sub = "invariant") or global (ip_sub = "no_norm")
# median normalisation via proDA::median_normalization.
#
# Returns: list(protein_matrix, design_df, acc_to_gene, ebn_ref_vector,
#               pipeline_label, n_prots_input, n_prots_processed,
#               normalization_method, processed)
#          OR list(error = "message") on failure.
# =============================================================================
lf_core_summarization_proda <- function(protein_matrix,
                                         design_df,
                                         acc_to_gene,
                                         ip_sub,
                                         igg_conditions = character(0),
                                         ebn_ref_vector = NULL,
                                         # FALSE: the input was normalized upstream (FragPipe
                                         # with IonQuant normalization, per fragpipe.workflow);
                                         # the global median step is skipped. Ignored on the
                                         # invariant-protein route, whose panel normalization
                                         # is always applied.
                                         apply_normalization = TRUE,
                                         logger         = message) {

  if (!isTRUE(.PRODA_AVAILABLE))
    return(list(error = paste0(
      "proDA not installed. ",
      "Install via: BiocManager::install(c('proDA','QFeatures','MsCoreUtils'))")))
  ensure_proda_loaded()

  # ── Step 1: IgG exclusion removed ──
  # All conditions (including IgG/controls) are included in reference-panel
  # selection. Truly invariant proteins should be invariant across ALL
  # conditions. If fewer than .INVARIANT_MIN_REFS references are found, the
  # pipeline automatically falls back to global median normalization.
  igg_conditions <- character(0)  # kept for API compat, always empty

  # ── Step 2: Filter proteins with zero non-NA values ──
  n_prots_input <- nrow(protein_matrix)
  keep_prots    <- apply(protein_matrix, 1, function(r) sum(!is.na(r)) >= 1)
  protein_matrix <- protein_matrix[keep_prots, , drop = FALSE]
  n_prots_processed <- nrow(protein_matrix)
  logger(paste0("ProDA core: ", n_prots_input, " \u2192 ", n_prots_processed, " proteins after filter"))

  if (n_prots_processed < 10)
    return(list(error = paste0("ProDA core: too few proteins after filtering (",
                               n_prots_processed, ").")))

  # ── Step 3: reference-panel selection and normalisation ──
  ebn_ref_vector_out    <- NULL
  normalization_method  <- "none"
  ebn_fallback          <- FALSE   # set TRUE if the panel fails and falls back to global median
  ebn_step              <- NA_integer_   # relaxation step that produced the panel
  ebn_thresholds        <- NULL          # thresholds in force at that step
  ebn_n_steps           <- NA_integer_   # length of the relaxation ladder
  n_ebn_singleton       <- NA_integer_   # candidates with <2 observed values
  n_ebn_zero_var        <- NA_integer_   # candidates with 2+ values, var == 0
  ebn_strict            <- NA             # panel drawn before amplitude relaxed
  ebn_strict_max_step   <- NA_integer_
  ebn_n_selected_prematch <- NA_integer_  # panel size as the selector returned it
  ebn_n_spike_in_rows     <- NA_integer_  # members surviving the matrix match
  ebn_n_selected        <- NA_integer_
  # EBN (ip_sub = "ebn") was replaced by invariant-protein median normalization.
  # A stale ip_sub (e.g. restored from an old session) must not run unnormalized.
  if (identical(ip_sub, "ebn"))
    return(list(error = paste0(
      "ProDA core: EBN normalization has been replaced by invariant-protein ",
      "median normalization. Please re-run the analysis wizard.")))
  # .inv: route requested (fixed); .inv_route(): still on it (FALSE after fallback)
  .inv       <- identical(ip_sub, "invariant")
  .inv_route <- function() identical(ip_sub, "invariant")
  inv_diag   <- NULL

  if (.inv_route()) {
    if (is.null(ebn_ref_vector)) {
      logger("ProDA core: selecting invariant-protein reference panel from protein matrix...")
      # Build MSstats-compatible long-format adapter for select_invariant_protein_panel().
      # Include IgG/control samples  -" they ARE the source of background reference
      # proteins. The panel is invariant across ALL conditions (bait + IgG),
      # which is exactly the definition of non-specific background.
      # Zero-variance proteins used to be excluded here on the grounds that they
      # are imputed-minimum artifacts. This route imputes nothing, so that
      # hazard cannot arise, and excluding invariant proteins from an
      # invariance-based reference panel has no methodological basis. They are
      # now counted and reported instead of dropped; Gate 1 in
      # select_invariant_protein_panel still enforces presence in every run.
      n_prots  <- nrow(protein_matrix)
      prot_ids <- rownames(protein_matrix)
      samp_ids <- colnames(protein_matrix)
      cond_map <- setNames(design_df$condition, design_df$sample_name)

      # Count, but no longer exclude, the two degenerate populations.
      n_obs_ebn   <- apply(protein_matrix, 1, function(r) sum(!is.na(r)))
      obs_var_ebn <- apply(protein_matrix, 1, function(r) {
        x <- r[!is.na(r)]; if (length(x) < 2) 0 else var(x)
      })
      n_ebn_singleton <- sum(n_obs_ebn < 2L)
      n_ebn_zero_var  <- sum(n_obs_ebn >= 2L &
                             !is.na(obs_var_ebn) & obs_var_ebn == 0)
      logger(paste0("ProDA core: invariant-protein panel candidate pool = all ", length(prot_ids),
                    " proteins (", n_ebn_singleton,
                    " with <2 observed values, ", n_ebn_zero_var,
                    " with 2+ observed values and zero variance; both retained)"))

      long_rows <- list()
      idx <- 0L
      for (pi in seq_along(prot_ids)) {
        for (si in seq_along(samp_ids)) {
          val <- protein_matrix[pi, si]
          if (!is.na(val)) {
            idx <- idx + 1L
            long_rows[[idx]] <- list(
              ProteinName = prot_ids[pi],
              Run         = samp_ids[si],
              Condition   = cond_map[[samp_ids[si]]] %||% samp_ids[si],
              Intensity   = val
            )
          }
        }
      }
      ebn_long <- do.call(rbind.data.frame, long_rows)

      ebn_ref_vector <- select_invariant_protein_panel(ebn_long,
                                                   intensity_is_log2 = TRUE,
                                                   logger = logger)
      # Read the panel criteria off the returned vector before the fallback
      # branch below can set it to NULL. attr() on NULL is NULL, not an error.
      ebn_step       <- attr(ebn_ref_vector, "ebn_step")       %||% NA_integer_
      ebn_thresholds <- attr(ebn_ref_vector, "ebn_thresholds")
      ebn_n_selected <- attr(ebn_ref_vector, "ebn_n_selected") %||% NA_integer_
      ebn_n_steps    <- attr(ebn_ref_vector, "ebn_n_steps")    %||% NA_integer_
      ebn_strict          <- attr(ebn_ref_vector, "ebn_strict")          %||% NA
      ebn_strict_max_step <- attr(ebn_ref_vector, "ebn_strict_max_step") %||% NA_integer_
      ebn_n_selected_prematch <- length(ebn_ref_vector %||% character(0))
      inv_diag <- list(
        invariant_centre_method   = attr(ebn_ref_vector, "invariant_centre_method")   %||% NA_character_,
        invariant_centre          = attr(ebn_ref_vector, "invariant_centre"),
        invariant_mode_centre     = attr(ebn_ref_vector, "invariant_mode_centre"),
        invariant_median_centre   = attr(ebn_ref_vector, "invariant_median_centre"),
        invariant_peak_ratio      = attr(ebn_ref_vector, "invariant_peak_ratio")      %||% NA_real_,
        invariant_centre_gap      = attr(ebn_ref_vector, "invariant_centre_gap")      %||% NA_real_,
        invariant_centre_disagree = isTRUE(attr(ebn_ref_vector, "invariant_centre_disagree")),
        invariant_mode_ambiguous  = isTRUE(attr(ebn_ref_vector, "invariant_mode_ambiguous")))
      if (is.null(ebn_ref_vector) || length(ebn_ref_vector) < .INVARIANT_MIN_REFS) {
        logger(paste0("ProDA core: invariant-protein panel could not identify \u2265", .INVARIANT_MIN_REFS, " invariant reference proteins \u2014 auto-fallback to global median"))
        ip_sub       <- "no_norm"
        ebn_fallback <- TRUE
        ebn_ref_vector <- NULL
      }
      logger(paste0("ProDA core: invariant-protein panel selected ", length(ebn_ref_vector),
                    " reference proteins"))
    }

    # Guard: only attempt invariant-protein normalization if not already fallen back
    if (.inv_route()) {
      spike_in_rows <- which(rownames(protein_matrix) %in% ebn_ref_vector)
      # The panel the ladder accepted is not necessarily the panel that reaches
      # normalization: members absent from the filtered matrix drop out here, so
      # a panel that passed the floor can arrive below it. Report both sizes so
      # the reason a fallback fired is visible.
      .n_panel_pre   <- length(ebn_ref_vector %||% character(0))
      ebn_n_spike_in_rows <- length(spike_in_rows)
      logger(paste0("ProDA core: invariant-protein panel ", .n_panel_pre,
                    " selected \u2192 ", ebn_n_spike_in_rows,
                    " present in the filtered matrix"))
      if (ebn_n_spike_in_rows < .n_panel_pre)
        logger(paste0("ProDA core: invariant-protein panel shrink at matrix intersection: ",
                      .n_panel_pre - ebn_n_spike_in_rows,
                      " reference protein(s) selected by the ladder are not in",
                      " the filtered matrix"))
      if (length(spike_in_rows) < .INVARIANT_MIN_REFS) {
        logger(paste0("ProDA core: invariant-protein panel reference proteins not in filtered matrix \u2014 auto-fallback to global median"))
        ip_sub       <- "no_norm"
        ebn_fallback <- TRUE
      }
    }

    if (.inv_route()) {
    protein_matrix <- tryCatch(
      proDA::median_normalization(protein_matrix, spike_in_rows = spike_in_rows),
      error = function(e) {
        logger(paste0("ProDA: median_normalization error: ", e$message)); NULL
      }
    )
    if (is.null(protein_matrix))
      return(list(error = "ProDA core: median_normalization failed."))

    ebn_ref_vector_out   <- ebn_ref_vector
    # Two independent reasons to call the panel low-confidence: too few members,
    # or drawn after the between-condition amplitude gate had already relaxed.
    # A higher floor buys panel size by spending invariance, so the step reached
    # is a signal in its own right and not a proxy for the count.
    .ebn_few    <- length(ebn_ref_vector) < .INVARIANT_MIN_REFS_IDEAL
    .ebn_loose  <- identical(ebn_strict, FALSE)
    if (.ebn_few || .ebn_loose) {
      normalization_method <- "invariant_low_confidence"
      pipeline_label       <- paste0("ProDA + invariant-protein median normalization, ",
                                     length(ebn_ref_vector), " refs \u2014 low confidence")
      logger(paste0("ProDA core: invariant-protein median normalisation applied (",
                    length(ebn_ref_vector), " refs; low confidence \u2014 ",
                    paste(c(if (.ebn_few) paste0("below ideal ",
                                                 .INVARIANT_MIN_REFS_IDEAL),
                            if (.ebn_loose) paste0("drawn at ladder step ",
                                                   ebn_step,
                                                   ", past the strict-amplitude",
                                                   " boundary of step ",
                                                   ebn_strict_max_step)),
                          collapse = "; "), ")"))
    } else {
      normalization_method <- "invariant_median"
      pipeline_label       <- "ProDA + invariant-protein median normalization"
      logger(paste0("ProDA core: invariant-protein median normalisation applied"))
    }
    }  # end inner if (.inv_route())  -" invariant-protein normalization

  }  # end outer if (.inv_route())

  if (!.inv && !isTRUE(apply_normalization)) {
    # Global route, input normalized upstream (FragPipe/IonQuant): no OSMSS normalization.
    normalization_method <- "none_input_normalized"
    pipeline_label       <- "Global (ProDA; no OSMSS normalization: input normalized upstream)"
    logger("ProDA core: global median normalization NOT applied (input normalized upstream)")
  } else if (!.inv_route()) {
    # Global median  -" applies to: no_norm path AND invariant-panel fallback path.
    # Apply a global column-median shift so that ProDA's EM starts from
    # comparable per-sample baselines.
    .norm_failed <- FALSE
    protein_matrix <- tryCatch(
      proDA::median_normalization(protein_matrix),
      error = function(e) {
        .norm_failed <<- TRUE
        logger(paste0("ProDA: global median normalization failed: ", e$message,
                      " \u2014 proceeding without normalization"))
        protein_matrix
      }
    )
    if (isTRUE(ebn_fallback)) {
      # Same tryCatch as the branch below: a caught error means the shift was
      # never applied, so the label must not claim it was.
      normalization_method <- if (isTRUE(.norm_failed))
        "invariant_fallback_global_median_FAILED_not_applied"
      else "invariant_fallback_global_median"
      pipeline_label       <- "Global (ProDA + invariant-protein panel not found \u2192 Global Median fallback)"
      logger("ProDA core: global median normalization applied (invariant-protein panel fallback)")
    } else {
      # A caught normalization error leaves the matrix unnormalized; say so
      # rather than claiming the shift was applied.
      normalization_method <- if (isTRUE(.norm_failed))
        "global_median_FAILED_not_applied" else "global_median"
      pipeline_label       <- "Global (ProDA + Global Median)"
      logger("ProDA core: global column-median normalisation applied (no_norm path)")
    }
  }

  # ── Step 4: Build MSstats-compatible ProteinLevelData shim ──
  # This allows downstream modules (heatmap, PCA, trends, rank_fc) to work
  # without modification, consuming processed$ProteinLevelData as normal.
  prot_ids <- rownames(protein_matrix)
  samp_ids <- colnames(protein_matrix)
  cond_map <- setNames(design_df$condition, design_df$sample_name)
  subj_map <- setNames(design_df$bioreplicate, design_df$sample_name)

  # Build ProteinLevelData via vectorized construction.
  # NOTE: do.call(rbind, list_of_lists) + as.data.frame() produces list columns,
  # not atomic vectors. Downstream modules (heatmap, PCA) subscript with character
  # vectors (e.g. acc_to_gene[proteins_v]) which fails with "invalid subscript
  # type 'list'" when the column is a list. Use direct vector construction instead.
  n_prots_out <- length(prot_ids)
  n_samps_out <- length(samp_ids)
  run_vec     <- rep(samp_ids, each = n_prots_out)
  group_vec   <- cond_map[run_vec]
  group_vec[is.na(group_vec)] <- run_vec[is.na(group_vec)]
  subj_vec    <- subj_map[run_vec]
  subj_vec[is.na(subj_vec)] <- run_vec[is.na(subj_vec)]
  pld <- data.frame(
    Protein        = rep(prot_ids, times = n_samps_out),
    LogIntensities = as.numeric(protein_matrix),   # column-major: all prots of samp1, then samp2, ...
    originalRUN    = run_vec,
    RUN            = run_vec,
    GROUP          = group_vec,
    SUBJECT        = subj_vec,
    censored       = FALSE,
    predicted      = FALSE,
    stringsAsFactors = FALSE
  )

  processed <- list(ProteinLevelData = pld)

  acc_to_gene_filt <- acc_to_gene[prot_ids]
  acc_to_gene_filt[is.na(acc_to_gene_filt)] <- prot_ids[is.na(acc_to_gene_filt)]

  out <- list(
    protein_matrix       = protein_matrix,
    design_df            = design_df,
    acc_to_gene          = acc_to_gene_filt,
    ebn_ref_vector       = ebn_ref_vector_out,
    pipeline_label       = pipeline_label,
    n_prots_input        = n_prots_input,
    n_prots_processed    = n_prots_processed,
    normalization_method = normalization_method,
    ebn_fallback         = isTRUE(ebn_fallback),
    ebn_step             = ebn_step,
    ebn_thresholds       = ebn_thresholds,
    ebn_n_selected       = ebn_n_selected,
    ebn_n_steps          = ebn_n_steps,
    ebn_candidates_singleton = n_ebn_singleton,
    ebn_candidates_zero_var  = n_ebn_zero_var,
    ebn_strict               = ebn_strict,
    ebn_strict_max_step      = ebn_strict_max_step,
    ebn_n_selected_prematch  = ebn_n_selected_prematch,
    ebn_n_spike_in_rows      = ebn_n_spike_in_rows,
    processed            = processed
  )
  if (.inv) out <- c(out, inv_diag %||% list(invariant_centre_method = NA_character_))
  out
}

# =============================================================================
# LF_CORE_SUMMARIZATION_SAINT
# Prepares a SAINTexpress-ready protein matrix cache from a build_proda_protein_matrix()
# result.  Unlike the ProDA path, NO normalization is applied  -" SAINTexpress has
# its own internal normalization model and applying a panel normalization beforehand would distort it.
# The function builds the same MSstats-compatible ProteinLevelData shim so PCA,
# heatmap, and trends modules work unchanged.
#
# Arguments:
#   protein_matrix   -" log2 abundance matrix (rows=proteins, cols=samples)
#   design_df        -" data.frame with columns Run, Condition, BioReplicate
#   acc_to_gene      -" named character vector Accession → GeneName
#   logger           -" logging function (default: message)
#
# Returns:
#   list with keys: protein_matrix, design_df, acc_to_gene, normalization_method,
#                   pipeline_label, n_prots_input, n_prots_processed,
#                   ebn_ref_vector, ebn_fallback, processed, aggregation_method
#   OR list(error = "message") on failure.
# =============================================================================
lf_core_summarization_saint <- function(protein_matrix,
                                         design_df,
                                         acc_to_gene,
                                         aggregation_method = "pd_protein_direct",
                                         logger             = message) {

  # ── Step 1: Filter proteins with zero non-NA values ──
  n_prots_input <- nrow(protein_matrix)
  keep_prots    <- apply(protein_matrix, 1, function(r) sum(!is.na(r)) >= 1)
  protein_matrix <- protein_matrix[keep_prots, , drop = FALSE]

  n_prots_processed <- nrow(protein_matrix)
  logger(paste0("SAINT core: ", n_prots_input, " \u2192 ", n_prots_processed,
                " proteins after filtering (no normalization applied)"))

  if (n_prots_processed < 5)
    return(list(error = paste0("SAINT core: too few proteins after filtering (",
                               n_prots_processed, ").")))

  # ── Step 3: Build MSstats-compatible ProteinLevelData shim ──
  # Same format as ProDA shim so PCA/heatmap/trends work unchanged.
  prot_ids   <- rownames(protein_matrix)
  sample_ids <- colnames(protein_matrix)

  # Map sample columns to conditions via design_df
  # build_proda_protein_matrix() returns lowercase: sample_name, condition, bioreplicate
  run_col  <- if ("sample_name" %in% colnames(design_df)) "sample_name" else "Run"
  cond_col <- if ("condition"   %in% colnames(design_df)) "condition"   else "Condition"
  subj_col <- if ("bioreplicate" %in% colnames(design_df)) "bioreplicate" else
              if ("BioReplicate" %in% colnames(design_df)) "BioReplicate" else NA_character_
  run_to_cond <- setNames(design_df[[cond_col]], design_df[[run_col]])
  run_to_subj <- if (!is.na(subj_col))
    setNames(as.character(design_df[[subj_col]]), design_df[[run_col]])
  else
    setNames(seq_len(nrow(design_df)), design_df[[run_col]])

  long_rows <- lapply(prot_ids, function(prot) {
    vals <- protein_matrix[prot, sample_ids]
    data.frame(
      Protein       = prot,
      LogIntensities = as.numeric(vals),
      originalRUN   = sample_ids,
      RUN           = sample_ids,
      GROUP         = run_to_cond[sample_ids],
      SUBJECT       = run_to_subj[sample_ids],
      censored      = FALSE,
      predicted     = FALSE,
      stringsAsFactors = FALSE
    )
  })
  pld <- do.call(rbind, long_rows)
  processed <- list(ProteinLevelData = pld)

  acc_to_gene_filt <- acc_to_gene[prot_ids]
  acc_to_gene_filt[is.na(acc_to_gene_filt)] <- prot_ids[is.na(acc_to_gene_filt)]

  list(
    protein_matrix       = protein_matrix,
    design_df            = design_df,
    acc_to_gene          = acc_to_gene_filt,
    ebn_ref_vector       = character(0),
    pipeline_label       = "IP (SAINTexpress)",
    n_prots_input        = n_prots_input,
    n_prots_processed    = n_prots_processed,
    normalization_method = "none",
    ebn_fallback         = FALSE,
    aggregation_method   = aggregation_method,
    processed            = processed
  )
}

# =============================================================================
# APPLY_IRS_NORMALIZATION
# Performs Internal Reference Scaling (IRS; Plubell et al. MCP 2017) in
# log2 space. Bridge/reference channels (Condition == "Norm" in annotation_map)
# are used as per-protein, per-plex scale references. For each protein, the
# plex-specific bridge mean is subtracted and the global bridge mean is added,
# aligning intensities across plexes multiplicatively (additively in log2 space).
# Falls back to column-wise median centering if no bridge channels are found.
#
# Args:
#   protein_matrix  -" log2-transformed proteins × channels matrix
#   annotation_map  -" data.frame with ColumnName, Condition, Mixture columns
#   logger          -" logging function
# Returns: IRS-corrected protein_matrix (same dimensions)
# =============================================================================
apply_irs_normalization <- function(protein_matrix, annotation_map, logger = message,
                                    uncorrectable = c("keep", "na")) {
  bridge_idx  <- !is.na(annotation_map$Condition) & annotation_map$Condition == "Norm"
  bridge_cols <- intersect(annotation_map$ColumnName[bridge_idx], colnames(protein_matrix))

  if (length(bridge_cols) == 0) {
    logger("[IRS] No bridge channels (Condition='Norm') found  -  falling back to median centering.")
    smed <- apply(protein_matrix, 2, median, na.rm = TRUE)
    gmed <- median(smed, na.rm = TRUE)
    for (j in seq_len(ncol(protein_matrix)))
      if (!is.na(smed[j]))
        protein_matrix[, j] <- protein_matrix[, j] + (gmed - smed[j])
    return(protein_matrix)
  }

  plexes <- unique(annotation_map$Mixture)
  n_prot <- nrow(protein_matrix)

  # ── Per-plex bridge mean per protein ──
  plex_bridge_mean <- matrix(NA_real_, nrow = n_prot, ncol = length(plexes),
                              dimnames = list(rownames(protein_matrix), plexes))
  for (pm in plexes) {
    plex_ann <- annotation_map[!is.na(annotation_map$Mixture) &
                                  annotation_map$Mixture == pm, , drop = FALSE]
    pb_cols  <- intersect(plex_ann$ColumnName[plex_ann$Condition == "Norm"],
                          colnames(protein_matrix))
    if (length(pb_cols) == 0) next
    plex_bridge_mean[, pm] <- rowMeans(protein_matrix[, pb_cols, drop = FALSE], na.rm = TRUE)
  }

  # ── Global bridge mean per protein (arithmetic mean in log2 space) ──
  global_bridge_mean <- rowMeans(plex_bridge_mean, na.rm = TRUE)

  # ── Apply IRS offset to every channel in each plex ──
  # uncorrectable = "keep" is the historical behaviour: a feature with no bridge
  # value in a plex keeps its UNCORRECTED values there, mixed into an otherwise
  # IRS-aligned row. That is a silent scale inconsistency, harmless only while the
  # caller pre-filters to complete cases. "na" instead marks the feature-by-plex
  # block missing, which is what MSstatsTMT's reference-channel normalization does:
  # the run is dropped for that feature, the feature survives on its other runs.
  .na_uncorrectable <- identical(match.arg(uncorrectable), "na")
  n_corr           <- 0L
  n_blanked        <- 0L
  n_plex_no_bridge <- 0L
  for (pm in plexes) {
    # plexes comes from unique(annotation_map$Mixture), which retains NA. Under
    # "keep" an NA entry was inert (it yields no columns and no valid rows); it is
    # skipped explicitly here so the blanking branch below cannot act on it.
    if (is.na(pm)) next
    plex_ann  <- annotation_map[!is.na(annotation_map$Mixture) &
                                   annotation_map$Mixture == pm, , drop = FALSE]
    plex_cols <- intersect(plex_ann$ColumnName, colnames(protein_matrix))
    pbm       <- plex_bridge_mean[, pm]
    valid     <- !is.na(pbm) & !is.na(global_bridge_mean)
    if (.na_uncorrectable && length(plex_cols) > 0L && any(!valid)) {
      # A plex carrying no bridge channel at all cannot be aligned for ANY feature.
      # That is a design problem, not per-feature missingness, so it is counted and
      # warned about separately rather than disappearing into the per-block total.
      if (!any(plex_ann$Condition == "Norm", na.rm = TRUE)) {
        n_plex_no_bridge <- n_plex_no_bridge + 1L
        logger(paste0("[IRS] WARNING: plex ", pm, " carries no bridge (Norm) channel. ",
                      "Its values cannot be aligned to the other plexes and were set ",
                      "missing for all ", sum(!valid), " feature(s)."))
      }
      for (col in plex_cols) protein_matrix[!valid, col] <- NA_real_
      n_blanked <- n_blanked + sum(!valid)
    }
    if (!any(valid)) next
    for (col in plex_cols)
      protein_matrix[valid, col] <- protein_matrix[valid, col] -
                                    pbm[valid] + global_bridge_mean[valid]
    n_corr <- n_corr + sum(valid)
  }

  logger(paste0("[IRS] Normalization complete: ", length(plexes), " plex(es), ",
                length(bridge_cols), " bridge channel(s), ",
                n_corr, " protein-plex correction(s) applied",
                if (.na_uncorrectable)
                  paste0("; ", n_blanked, " feature-plex block(s) had no bridge value ",
                         "and were set missing")
                else "", "."))
  # Attached rather than returned in a list: every caller assigns the result
  # straight back to a matrix, and changing the return type would break all three.
  # NOTE the attribute does NOT survive subsetting - read it immediately.
  #
  # Attached ONLY under the "na" policy. Under "keep" the returned object must be
  # identical() to what this function returned before the policy existed,
  # attributes included, so that the Global TMT routes - which never pass a policy
  # - are provably untouched. Every consumer already treats a NULL irs_stats as
  # "no statistics available", so nothing downstream needs a matching change.
  if (.na_uncorrectable)
    attr(protein_matrix, "irs_stats") <- list(
      n_corrected      = n_corr,
      n_blanked        = n_blanked,
      n_plex_no_bridge = n_plex_no_bridge,
      policy           = "na"
    )
  protein_matrix
}

# =============================================================================
# TMT_LIMMA_PROTEIN_SUMMARIZE
# Summarizes protein-level TMT data for the limma pipeline (all design options).
# Takes protein abundance columns from parsed TMT files, log2-transforms,
# and applies IRS normalization (design_opt=2) or column-wise median centering
# (design_opt=1 or 3). No imputation  -" deferred to the comparison step.
#
# Args:
#   design_opt  -" integer: 1 = single-plex, 2 = multi-plex IRS bridge,
#                3 = multi-plex run-covariate balanced
#
# Returns: list(protein_matrix, gene_names, conditions, annotation_map,
#               design_option, n_proteins, n_channels, route)
#          OR list(error = "message") on failure.
# =============================================================================
# Proteins in the uploaded file without a reporter value in any channel (identified but not quantified, e.g. PSMs
# excluded from quantification by Proteome Discoverer, or zeros only): listed apart from quantified-but-sparse proteins.
TMT_STAGE_NO_VALUES  <- "Input (no abundance values)"
TMT_REASON_NO_VALUES <- "Identified but not quantified in the uploaded file (no reporter value in any channel)"
TMT_FILTER_STAGES    <- c("Summarization (plex coverage)", "Before DE (missingness filter)", TMT_STAGE_NO_VALUES)

# Coverage of proteins removed by a TMT coverage filter, for the DE module's "Not Tested & Low Coverage" tab: the
# plexes a protein was quantified in, and its observed / annotated channels per condition (pooled channels included)
# and per plex. `mat` holds the values the filter saw (NA = missing), columns = annotation_map$ColumnName.
# Reporting only: nothing here enters a model.
tmt_coverage_profile <- function(mat, annotation_map, proteins, gene_lookup = NULL) {
  proteins <- intersect(proteins, rownames(mat))
  if (!length(proteins)) return(NULL)
  ann <- annotation_map[annotation_map$ColumnName %in% colnames(mat) & !is.na(annotation_map$Condition) &
                          nzchar(annotation_map$Condition), , drop = FALSE]
  obs  <- !is.na(mat[proteins, ann$ColumnName, drop = FALSE])
  plex <- if ("Mixture" %in% names(ann)) as.character(ann$Mixture) else rep("Mix1", nrow(ann))
  pl_ids <- unique(plex)
  pl_q <- vapply(pl_ids, function(p) rowSums(obs[, plex == p, drop = FALSE]) > 0, logical(length(proteins)))
  if (is.null(dim(pl_q))) pl_q <- matrix(pl_q, nrow = length(proteins))
  gn <- if (!is.null(gene_lookup)) unname(gene_lookup[proteins]) else rep(NA_character_, length(proteins))
  gn[is.na(gn) | !nzchar(gn)] <- proteins[is.na(gn) | !nzchar(gn)]
  out <- data.frame(Protein = proteins, Gene = gn, Plexes_quantified = paste0(rowSums(pl_q), " of ", length(pl_ids)),
                    stringsAsFactors = FALSE, check.names = FALSE)
  for (cn in unique(ann$Condition)) { k <- ann$Condition == cn
    out[[paste0(cn, " (channels observed)")]] <- paste0(rowSums(obs[, k, drop = FALSE]), "/", sum(k)) }
  if (length(pl_ids) > 1L) for (p in pl_ids) { k <- plex == p
    out[[paste0(p, " (channels observed)")]] <- paste0(rowSums(obs[, k, drop = FALSE]), "/", sum(k)) }
  rownames(out) <- NULL
  out
}

# The low-coverage table (compute_tmt_de()$one_arm_df) as the DE module shows and exports it: readable column names,
# and whether each protein is in the DE results / volcano. A protein with at least one value in both compared
# conditions was tested (real p-value, in the FDR correction and the volcano); one with no value at all in a compared
# condition has no estimable contrast, no p-value, and is not plotted. `result_df` = compute_tmt_de()$result_df.
tmt_low_coverage_display <- function(one_arm_df, result_df) {
  if (is.null(one_arm_df) || !nrow(one_arm_df)) return(one_arm_df)
  .cmp  <- strsplit(as.character(one_arm_df$Comparison[1]), " vs ", fixed = TRUE)[[1]]
  .test <- .cmp[1]; .ref <- if (length(.cmp) > 1) .cmp[2] else "reference"
  m <- match(one_arm_df$Protein, as.character(result_df$Protein))
  adjp <- if (!is.null(result_df$adj.pvalue)) result_df$adj.pvalue[m] else rep(NA_real_, nrow(one_arm_df))
  est  <- if (!is.null(result_df$Diag_estimable)) result_df$Diag_estimable[m] else !is.na(adjp)
  est[is.na(est)] <- FALSE
  empty <- ifelse(one_arm_df$n_valid_test == 0, .test, ifelse(one_arm_df$n_valid_ref == 0, .ref, NA_character_))
  out <- data.frame(Protein = one_arm_df$Protein, Gene = one_arm_df$Gene, Comparison = one_arm_df$Comparison,
                    check.names = FALSE, stringsAsFactors = FALSE)
  out[[paste0(.test, " (test): channels observed")]]      <- paste0(one_arm_df$n_valid_test, "/", one_arm_df$n_test)
  out[[paste0(.ref, " (reference): channels observed")]]  <- paste0(one_arm_df$n_valid_ref, "/", one_arm_df$n_ref)
  out[["Adequately covered in (>= 50% of channels)"]]     <- one_arm_df$DetectedIn
  out[["In DE table & volcano"]] <- ifelse(est, "Yes (tested)",
                                           paste0("No: no value in ", ifelse(is.na(empty), "one condition", empty),
                                                  " (no estimable contrast, no p-value)"))
  out[["adj. p"]] <- ifelse(est, adjp, NA_real_)
  out[["Zero values in input (compared channels)"]] <- one_arm_df$n_zero_origin
  out[["Peptides"]] <- one_arm_df$n_peptides
  out
}

tmt_limma_protein_summarize <- function(files, annotation_map, rv,
                                        design_opt = 1L, logger = message,
                                        phospho_files = NULL,
                                        min_plex_frac = 0.5,
                                        prot_qvalue_max = NULL,
                                        norm_channels = NULL,
                                        irs_uncorrectable = c("keep", "na"),
                                        apply_centering = TRUE) {

  if (is.null(annotation_map) || nrow(annotation_map) == 0)
    return(list(error = "tmt_limma_protein_summarize: annotation_map is empty."))
  if (length(files) == 0)
    return(list(error = "tmt_limma_protein_summarize: no files provided."))

  target_cols <- annotation_map$ColumnName
  gene_lookup <- character(0)
  combined    <- list()

  for (fi in seq_along(files)) {
    f     <- files[[fi]]
    pdata <- parse_pd_protein_data(f, rv)
    if (is.null(pdata) || nrow(pdata$protein_df) < 1) {
      logger(paste0("tmt_limma_protein_summarize: skipping file ", f$name,
                    " (no protein data)"))
      next
    }
    prot_df    <- pdata$protein_df
    gene_names <- pdata$gene_names

    if (!is.null(prot_qvalue_max)) {
      .pq_col <- grep("^Exp\\. q-value: Combined$", colnames(prot_df),
                      value = TRUE)[1]
      if (is.na(.pq_col))
        .pq_col <- grep("^Exp\\. q-value", colnames(prot_df), value = TRUE)[1]
      if (!is.na(.pq_col)) {
        .n_before <- nrow(prot_df)
        .pqv      <- suppressWarnings(as.numeric(prot_df[[.pq_col]]))
        .keep_q   <- !is.na(.pqv) & .pqv < prot_qvalue_max
        prot_df    <- prot_df[.keep_q, , drop = FALSE]
        gene_names <- gene_names[.keep_q]
        logger(paste0("tmt_limma_protein_summarize: protein q-value < ",
                      prot_qvalue_max, " retained ", nrow(prot_df), " of ",
                      .n_before, " proteins in ", f$name))
      } else {
        logger(paste0("tmt_limma_protein_summarize: no Exp. q-value column in ",
                      f$name, "; no protein-level confidence filter applied."))
      }
    }

    # Accession column
    acc_col <- grep("^Accession$", colnames(prot_df), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(acc_col)) acc_col <- colnames(prot_df)[2]
    prot_ids <- as.character(prot_df[[acc_col]])

    # Strip semicolon-separated accession groups (keep first)
    prot_ids <- trimws(sub(";.*", "", prot_ids))

    # Build gene lookup
    for (k in seq_along(prot_ids)) {
      gn <- if (!is.null(gene_names) && length(gene_names) >= k) gene_names[k] else prot_ids[k]
      if (!is.na(gn) && nzchar(gn)) gene_lookup[prot_ids[k]] <- gn
    }

    # Extract abundance columns present in this file
    file_cols    <- colnames(prot_df)
    file_targets <- intersect(target_cols, file_cols)

    # Fallback: PD assigns a unique FN: file-tracking ID per workflow, so the
    # phospho enrichment files (used to build annotation_map) often have a
    # different fraction number than the global files (e.g. F4 vs F2, F1 vs F3).
    # Normalise by stripping ": FN:" and re-match.
    #
    # IMPORTANT  -" per-plex safety: when phospho_files is supplied, restrict
    # the match to channels from the SAME plex (phospho_files[[fi]]$abund_cols)
    # so a global file for plex 2 cannot accidentally match plex 1 channels.
    # (After stripping FN: all plexes share identical channel names when the
    # TMT design and sample labels are the same across plexes.)
    if (length(file_targets) == 0) {
      plex_tgt <- if (!is.null(phospho_files) &&
                      fi <= length(phospho_files) &&
                      !is.null(phospho_files[[fi]]$abund_cols)) {
        phospho_files[[fi]]$abund_cols   # exact channels for this plex
      } else {
        target_cols
      }
      strip_fid <- function(x) sub(": F\\d+:", ":", x)
      norm_tgt  <- strip_fid(plex_tgt)
      norm_file <- strip_fid(file_cols)
      m_idx     <- match(norm_tgt, norm_file)
      has_m     <- !is.na(m_idx)
      if (any(has_m)) {
        colnames(prot_df)[m_idx[has_m]] <- plex_tgt[has_m]
        file_targets <- plex_tgt[has_m]
        logger(paste0("tmt_limma_protein_summarize: FN-normalised ",
                      sum(has_m), " column match(es) in ", f$name))
      }
    }

    if (length(file_targets) == 0) {
      logger(paste0("tmt_limma_protein_summarize: no matching abundance columns in ", f$name))
      next
    }

    sub_mat <- matrix(NA_real_, nrow = nrow(prot_df), ncol = length(file_targets),
                      dimnames = list(prot_ids, file_targets))
    for (j in seq_along(file_targets)) {
      sub_mat[, j] <- suppressWarnings(as.numeric(prot_df[[file_targets[j]]]))
    }

    combined[[fi]] <- sub_mat
  }

  if (length(combined) == 0)
    return(list(error = "tmt_limma_protein_summarize: no valid protein data extracted from files."))

  # ── Merge matrices across files by accession union ──
  all_accs <- unique(unlist(lapply(combined, rownames)))
  protein_matrix <- matrix(NA_real_, nrow = length(all_accs), ncol = length(target_cols),
                           dimnames = list(all_accs, target_cols))

  for (mat in combined) {
    row_idx <- match(rownames(mat), all_accs)
    col_idx <- match(colnames(mat), target_cols)
    valid_r <- !is.na(row_idx)
    valid_c <- !is.na(col_idx)
    if (any(valid_r) && any(valid_c)) {
      protein_matrix[row_idx[valid_r], col_idx[valid_c]] <-
        mat[valid_r, valid_c, drop = FALSE]
    }
  }

  # ── Sanitize zeros BEFORE log2 (zeros → NA to avoid log2(0) = -Inf) ──
  protein_matrix[protein_matrix <= 0] <- NA

  # ── Log2 transform ──
  protein_matrix <- log2(protein_matrix)

  # ── Defense-in-depth: remove any residual non-finite values ──
  protein_matrix[!is.finite(protein_matrix)] <- NA

  # ── Plex coverage filter (multi-plex only) ──
  # TMT missingness is block-structured: a protein is quantified in a plex or it
  # is not, and if not, ALL that plex's channels are missing together. Files are
  # merged by accession union, so a protein seen in 1 of 4 plexes would otherwise
  # enter with 75% of its channels empty. Require coverage in enough plexes
  # instead of filling those blocks in later.
  n_proteins_input   <- nrow(protein_matrix)
  n_dropped_coverage <- 0L
  coverage_dropped   <- NULL          # proteins removed here, with their coverage (reporting only)
  min_plexes_required <- NA_integer_
  plex_ids <- unique(annotation_map$Mixture[!is.na(annotation_map$Mixture)])

  if (length(plex_ids) >= 2L && nrow(protein_matrix) > 0L) {
    plex_present <- matrix(FALSE, nrow = nrow(protein_matrix), ncol = length(plex_ids),
                           dimnames = list(rownames(protein_matrix), as.character(plex_ids)))
    for (pi in seq_along(plex_ids)) {
      pcols <- intersect(
        annotation_map$ColumnName[!is.na(annotation_map$Mixture) &
                                    annotation_map$Mixture == plex_ids[pi]],
        colnames(protein_matrix))
      if (length(pcols) == 0L) next
      plex_present[, pi] <- rowSums(!is.na(protein_matrix[, pcols, drop = FALSE])) > 0L
    }
    min_plexes <- max(1L, ceiling(min_plex_frac * length(plex_ids)))
    keep_cov   <- rowSums(plex_present) >= min_plexes
    n_dropped_coverage <- sum(!keep_cov)
    min_plexes_required <- min_plexes
    if (n_dropped_coverage > 0L) {
      # Listed, not dropped silently: the DE module's "Not Tested & Low Coverage" tab shows these proteins with the
      # plexes and channels they were seen in. Computed before normalization, on the values the filter saw.
      coverage_dropped <- tmt_coverage_profile(protein_matrix, annotation_map, rownames(protein_matrix)[!keep_cov],
                                               gene_lookup)
      if (!is.null(coverage_dropped)) {
        .npl <- rowSums(plex_present[!keep_cov, , drop = FALSE])
        # quantified in no plex at all: in the file, identified, but without any reporter value
        coverage_dropped <- data.frame(coverage_dropped[1:2],
          Stage  = ifelse(.npl == 0, TMT_STAGE_NO_VALUES, "Summarization (plex coverage)"),
          Reason = ifelse(.npl == 0, TMT_REASON_NO_VALUES,
                          paste0("Quantified in ", .npl, " of ", length(plex_ids), " plexes; at least ", min_plexes,
                                 " (", round(100 * min_plex_frac), "% of plexes) required")),
          coverage_dropped[-(1:2)], check.names = FALSE, stringsAsFactors = FALSE)
      }
      protein_matrix <- protein_matrix[keep_cov, , drop = FALSE]
      all_accs       <- rownames(protein_matrix)
      logger(paste0("tmt_limma_protein_summarize: ", n_dropped_coverage,
                    " proteins removed (quantified in fewer than ", min_plexes,
                    " of ", length(plex_ids), " plexes)."))
    }
  }

  # ── Normalization ──
  # Two distinct jobs, applied in sequence (not either/or):
  #  1. Column median centering equalises per-channel loading differences.
  #  2. IRS (opt=2) aligns plexes to each other via the bridge channels.
  # IRS applies the SAME per-protein offset to every channel in a plex, so it
  # cannot correct within-plex loading. The published workflow is sample-loading
  # normalization followed by IRS; previously step 1 was skipped whenever IRS ran.
  # apply_centering = FALSE (single plex only): the user reported in the wizard
  # that the abundances were already normalized across channels before upload
  # (see tmt_normalization_check()); the values are then used as they are.
  # Multi-plex designs always center: IRS cannot correct within-plex loading.
  .center <- isTRUE(apply_centering) || !identical(as.integer(design_opt), 1L)
  sample_medians <- apply(protein_matrix, 2, median, na.rm = TRUE)
  .const_cols    <- if (is.null(norm_channels)) colnames(protein_matrix) else
                    intersect(norm_channels, colnames(protein_matrix))
  if (length(.const_cols) == 0L) .const_cols <- colnames(protein_matrix)
  global_median  <- median(sample_medians[.const_cols], na.rm = TRUE)
  if (.center) {
    for (j in seq_len(ncol(protein_matrix))) {
      if (!is.na(sample_medians[j]))
        protein_matrix[, j] <- protein_matrix[, j] + (global_median - sample_medians[j])
    }
  }
  norm_applied <- if (.center) "column-wise median centering (log2)"
                  else "none in OSMSS (input reported as already normalized across channels)"
  if (!.center)
    logger("tmt_limma_protein_summarize: column-wise median centering NOT applied (input reported as already normalized across channels).")

  prot_irs_stats <- NULL
  if (identical(as.integer(design_opt), 2L)) {
    protein_matrix <- apply_irs_normalization(protein_matrix, annotation_map, logger,
                                              uncorrectable = match.arg(irs_uncorrectable))
    # Read immediately: the attribute does not survive the subsetting below.
    prot_irs_stats <- attr(protein_matrix, "irs_stats")
    norm_applied   <- paste0(norm_applied, " followed by IRS bridge-channel alignment")
  }

  # ── Fill in gene names for any accessions not found ──
  missing_gn <- setdiff(all_accs, names(gene_lookup))
  for (acc in missing_gn) gene_lookup[acc] <- acc

  conds <- unique(annotation_map$Condition[annotation_map$Condition != "Norm" &
                                           !is.na(annotation_map$Condition)])

  logger(paste0("tmt_limma_protein_summarize: ", nrow(protein_matrix), " proteins x ",
                ncol(protein_matrix), " channels"))

  list(
    protein_matrix = protein_matrix,
    gene_names     = gene_lookup,
    conditions     = sort(conds),
    annotation_map = annotation_map,
    design_option  = as.integer(design_opt),
    n_proteins           = nrow(protein_matrix),
    n_proteins_input     = n_proteins_input,
    n_dropped_coverage   = n_dropped_coverage,
    coverage_dropped     = coverage_dropped,
    min_plexes_required  = min_plexes_required,
    irs_stats            = prot_irs_stats,
    n_plexes             = length(plex_ids),
    min_plex_frac        = min_plex_frac,
    normalization        = norm_applied,
    normalization_applied = .center,
    n_channels     = ncol(protein_matrix),
    route          = "limma_protein"
  )
}

# tmt_flag_shared_peptides() --------------------------------------------------
# Identify peptides that map to more than one master protein, from the
# semicolon-delimited master accession column of a Proteome Discoverer peptide
# export. This is the only shared-peptide evidence a peptide-level export
# carries: PD reports unique/razor verdicts on PSM exports only.
#
# Returns list(shared = logical vector, n_candidates = integer vector,
#              candidates = list of character vectors).
tmt_flag_shared_peptides <- function(pep_df, prot_col) {
  acc_raw <- if (!is.null(prot_col) && prot_col %in% colnames(pep_df))
    as.character(pep_df[[prot_col]]) else rep(NA_character_, nrow(pep_df))
  acc_raw[is.na(acc_raw)] <- ""

  candidates <- lapply(strsplit(acc_raw, ";", fixed = TRUE), function(x) {
    x <- trimws(x)
    x[nzchar(x)]
  })
  n_candidates <- vapply(candidates, length, integer(1))

  # n_candidates == 0 (NA or empty accession) is flagged shared, not unique:
  # such a row cannot be attributed to any protein, and the first-accession
  # expression downstream would silently make it the protein "".
  list(shared       = n_candidates != 1L,
       n_candidates = n_candidates,
       candidates   = candidates)
}

# tmt_razor_assign_shared() ---------------------------------------------------
# Parsimony (razor) assignment for peptides matching more than one master
# protein, as MaxQuant applies Occam's razor: each shared peptide goes to the
# candidate protein supported by the most UNIQUE peptide sequences in the
# pooled dataset.
#
# This is deliberately NOT the rule in the PSM razor block of
# qf_summarize_to_protein_matrix. That one counts, per peptide sequence, how
# many rows name each protein - a quantity that exists only when each row
# carries a single protein, i.e. PSM input. In a peptide export one row carries
# the whole candidate list, so every candidate would tie at one and which.max()
# would return the first listed accession: exactly the positional behaviour
# razor assignment exists to replace.
#
# Deterministic: unique support counts DISTINCT peptide sequences, so a peptide
# observed in several files is not double-weighted, and candidates are sorted
# lexicographically before the maximum is taken, so equal support resolves
# alphabetically rather than by position in the export.
#
# Args:
#   peptides   - character vector, one peptide sequence per row
#   candidates - list of character vectors, parsed master accession candidates
#                per row (tmt_flag_shared_peptides()$candidates)
#
# Returns: list(protein = character, one assignment per row (NA where the row
#                         cannot be razored),
#               assigned = logical, TRUE for rows reassigned by parsimony,
#               dropped  = logical, TRUE where no candidate has unique support
#                          (or there is no candidate at all),
#               unique_support = named integer vector of unique-peptide counts).
tmt_razor_assign_shared <- function(peptides, candidates) {
  n      <- length(candidates)
  n_cand <- vapply(candidates, length, integer(1))
  peptides <- as.character(peptides)

  # Unique support: the number of distinct peptide sequences for which an
  # accession is the SOLE candidate.
  sole <- rep(NA_character_, n)
  if (any(n_cand == 1L))
    sole[n_cand == 1L] <- vapply(candidates[n_cand == 1L], `[`, character(1), 1L)
  has_sole <- !is.na(sole)
  unique_support <- integer(0)
  if (any(has_sole)) {
    pairs <- unique(cbind(peptides[has_sole], sole[has_sole]))
    tab   <- table(pairs[, 2L])
    unique_support <- stats::setNames(as.integer(tab), names(tab))
  }

  protein  <- rep(NA_character_, n)
  assigned <- rep(FALSE, n)
  dropped  <- rep(FALSE, n)
  protein[has_sole] <- sole[has_sole]

  for (i in which(n_cand > 1L)) {
    cand <- sort(candidates[[i]])                 # lexicographic tie-break
    sup  <- if (length(unique_support))
              unname(unique_support[match(cand, names(unique_support))]) else
              rep(NA_integer_, length(cand))
    sup[is.na(sup)] <- 0L
    if (all(sup == 0L)) {
      # No candidate has any unique support: the protein this peptide would
      # create is identified by shared peptides only. Not razorable - drop it.
      dropped[i] <- TRUE
    } else {
      protein[i]  <- cand[which.max(sup)]
      assigned[i] <- TRUE
    }
  }
  # No candidate at all (NA or empty accession): not attributable to anything.
  dropped[n_cand == 0L] <- TRUE

  list(protein        = protein,
       assigned       = assigned,
       dropped        = dropped,
       unique_support = unique_support)
}

# =============================================================================
# TMT_LIMMA_QF_SUMMARIZE
# Summarizes PSM/peptide-level TMT data for the limma pipeline (all design options).
# Applies co-isolation interference filter (<=interference_threshold%),
# delegates to qf_summarize_to_protein_matrix() with peptide-level normalization,
# then applies IRS (design_opt=2) or column-wise median centering (opt=1/3).
#
# Args:
#   design_opt             -" integer: 1 = single-plex, 2 = multi-plex IRS bridge,
#                           3 = multi-plex run-covariate balanced
#   interference_threshold  -" numeric: max allowed isolation interference % (default 50)
#
# Returns: same structure as tmt_limma_protein_summarize, plus:
#   n_peptides_input, interference_filtered, interference_threshold,
#   pep_count_per_protein
#          OR list(error = "message") on failure.
# =============================================================================
tmt_limma_qf_summarize <- function(files, annotation_map, data_level,
                                   rv, design_opt = 1L,
                                   interference_threshold = 50,
                                   shared_peptides = "unique_only",
                                   logger = message,
                                   apply_centering = TRUE) {
  # apply_centering = FALSE (single plex only): no peptide-level median centering,
  # as in tmt_limma_protein_summarize().
  .center <- isTRUE(apply_centering) || !identical(as.integer(design_opt), 1L)

  # Hard error, not a fallback: an unrecognised value here would otherwise be
  # silently treated as one of the accepted settings and misreported in
  # provenance as a filter that did not run.
  if (!(length(shared_peptides) == 1L && !is.na(shared_peptides) &&
        shared_peptides %in% c("unique_only", "keep_all", "razor")))
    stop("tmt_limma_qf_summarize: shared_peptides must be one of ",
         "\"unique_only\", \"keep_all\" or \"razor\".")

  if (!isTRUE(.PRODA_AVAILABLE))
    return(list(error = paste0(
      "QFeatures / MsCoreUtils not installed. ",
      "Install via: BiocManager::install(c('proDA','QFeatures','MsCoreUtils'))")))

  if (is.null(annotation_map) || nrow(annotation_map) == 0)
    return(list(error = "tmt_limma_qf_summarize: annotation_map is empty."))
  if (length(files) == 0)
    return(list(error = "tmt_limma_qf_summarize: no files provided."))

  .scaled_only <- vapply(files, function(f)
    identical(f$abundance_source, "scaled_only"), logical(1))
  if (any(.scaled_only))
    return(list(error = paste0(
      "This PD export contains only 'Abundances (Scaled)' columns for ",
      sum(.scaled_only), " file(s). OSMSS requires raw 'Abundance: F#:' columns: ",
      "scaled abundances are already normalised by Proteome Discoverer, and ",
      "re-normalising them produces invalid quantification. Re-export from PD ",
      "with raw abundances included.")))

  target_cols <- annotation_map$ColumnName
  gene_lookup <- character(0)

  # ── Collect peptide/PSM data from all files ──
  all_pep_rows <- list()
  interference_col_found <- FALSE
  # Rows kept BECAUSE the interference value was missing, counted separately
  # from rows kept because they passed the threshold. The two are different
  # claims and the methods text should be able to state both numerically.
  n_interf_kept_na   <- 0L
  n_interf_kept_pass <- 0L
  # Shared-peptide accounting, accumulated across files. Independent of the
  # interference counters above: the two filters run in sequence and both are
  # reportable.
  n_pep_total  <- 0L
  n_pep_shared <- 0L
  shared_detection_method <- NULL

  for (fi in seq_along(files)) {
    f    <- files[[fi]]
    pdata <- tryCatch(
      parse_pd_peptide_data(f),
      error = function(e) {
        logger(paste0("tmt_limma_qf_summarize: parse error in ", f$name, ": ", e$message))
        NULL
      }
    )
    if (is.null(pdata) || nrow(pdata) < 1) {
      logger(paste0("tmt_limma_qf_summarize: no peptide data in ", f$name))
      next
    }

    pep_df <- pdata

    # Protein column (Master Protein Accessions or Accession)
    prot_col <- grep("Master.Protein.Accessions|Protein.Accessions|^Accession$",
                     colnames(pep_df), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(prot_col)) prot_col <- colnames(pep_df)[1]

    # Peptide column (Annotated Sequence or Sequence)
    pep_col <- grep("Annotated.Sequence|Sequence$", colnames(pep_df),
                    value = TRUE, ignore.case = TRUE)[1]
    if (is.na(pep_col)) pep_col <- colnames(pep_df)[2]

    # Interference column
    int_col <- grep("Isolation.Interference|Co.isolation", colnames(pep_df),
                    value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(int_col)) interference_col_found <- TRUE

    # ── Apply co-isolation interference filter ──
    if (!is.na(int_col)) {
      interf_vals  <- suppressWarnings(as.numeric(pep_df[[int_col]]))
      keep_interf  <- is.na(interf_vals) | interf_vals <= interference_threshold
      n_interf_kept_na   <- n_interf_kept_na   + sum(is.na(interf_vals))
      n_interf_kept_pass <- n_interf_kept_pass +
                            sum(!is.na(interf_vals) &
                                interf_vals <= interference_threshold)
      n_before     <- nrow(pep_df)
      pep_df       <- pep_df[keep_interf, , drop = FALSE]
      logger(paste0("tmt_limma_qf_summarize: interference filter (",
                    f$name, "): ", n_before, " -> ", nrow(pep_df),
                    " rows (\u2264", interference_threshold, "% isolation interference)"))
    }

    # Find abundance columns present in this file
    file_cols    <- colnames(pep_df)
    file_targets <- intersect(target_cols, file_cols)
    if (length(file_targets) == 0) {
      logger(paste0("tmt_limma_qf_summarize: no matching intensity cols in ", f$name))
      next
    }

    # -- Shared-peptide handling --
    # Runs after the interference filter, so the two row counts are independent
    # and both reportable. A peptide-level export carries no unique/razor
    # verdict of its own, so the master-accession list is the evidence used.
    .shared_flags <- tmt_flag_shared_peptides(pep_df, prot_col)
    n_pep_total   <- n_pep_total  + nrow(pep_df)
    n_pep_shared  <- n_pep_shared + sum(.shared_flags$shared)
    .this_method  <- "master_accession_count"
    if (is.null(shared_detection_method)) {
      shared_detection_method <- .this_method
    } else if (!identical(shared_detection_method, .this_method)) {
      shared_detection_method <- "mixed"
      logger(paste0("tmt_limma_qf_summarize: shared-peptide detection method ",
                    "differs between files - recorded as 'mixed'."))
    }

    if (identical(shared_peptides, "unique_only")) {
      n_before_shared <- nrow(pep_df)
      pep_df <- pep_df[!.shared_flags$shared, , drop = FALSE]
      logger(paste0("tmt_limma_qf_summarize: shared-peptide filter (",
                    f$name, "): ", n_before_shared, " -> ", nrow(pep_df),
                    " rows (peptides matching >1 master protein removed; ",
                    "detection: ", .this_method, ")"))
      if (nrow(pep_df) == 0L) {
        logger(paste0("tmt_limma_qf_summarize: no unique peptides left in ",
                      f$name, " after the shared-peptide filter"))
        next
      }
    }

    # Build standardized data.frame for aggregation
    sub_df <- data.frame(
      Protein = trimws(sub(";.*", "", as.character(pep_df[[prot_col]]))),
      Peptide = as.character(pep_df[[pep_col]]),
      # Full, untruncated master accession list. Travels through rbind into
      # combined_pep_df so parsimony (razor) assignment can work on the pooled
      # table. NOT passed to qf_summarize_to_protein_matrix, whose protein_col
      # stays "Protein"; that function selects its columns by name, so the extra
      # column is inert there.
      MasterAccessions = as.character(pep_df[[prot_col]]),
      stringsAsFactors = FALSE
    )
    for (col in file_targets) {
      sub_df[[col]] <- suppressWarnings(as.numeric(pep_df[[col]]))
    }
    # Fill missing target columns with NA
    for (col in setdiff(target_cols, file_targets)) sub_df[[col]] <- NA_real_

    all_pep_rows[[fi]] <- sub_df
  }

  if (!interference_col_found) {
    logger(paste0("TMT Limma QF: Isolation Interference column not found \u2014 ",
                  "proceeding without filter. WARNING: co-isolation interference ",
                  "may distort fold changes."))
  }

  all_pep_rows <- Filter(Negate(is.null), all_pep_rows)
  if (length(all_pep_rows) == 0) {
    # Name the shared-peptide filter when it is the plausible cause, so the user
    # is not left guessing which of the two row filters emptied the table.
    if (identical(shared_peptides, "unique_only") && n_pep_shared >= n_pep_total &&
        n_pep_total > 0L)
      return(list(error = paste0(
        "tmt_limma_qf_summarize: no peptide data available after filtering - ",
        "every peptide matched more than one master protein and was removed by ",
        "the shared-peptide filter ('Unique peptides only'). Re-run with ",
        "'Keep all peptides' to include them.")))
    return(list(error = "tmt_limma_qf_summarize: no peptide data available after filtering."))
  }

  combined_pep_df <- do.call(rbind, all_pep_rows)
  logger(paste0("tmt_limma_qf_summarize: combined peptide table: ",
                nrow(combined_pep_df), " rows"))

  # One file per plex: each file holds the same peptide again, with values only in
  # its own plex's channels (tmt_disambiguate_plex_files() keeps the channel sets
  # disjoint). Rows of one peptide from different files are combined into a single
  # feature, as in one export that holds every plex, so both upload layouts give
  # the same peptide table. The k-th row of a Protein + Peptide pair in one file is
  # combined with the k-th row of that pair in the others; repeats within a file
  # stay separate features, as before. A single file is left as it is.
  n_pep_merged_across_files <- 0L
  if (length(all_pep_rows) > 1L) {
    .key <- unlist(lapply(all_pep_rows, function(d) {
      k <- paste(d$Protein, d$Peptide, sep = "\r")
      paste(k, stats::ave(seq_along(k), k, FUN = seq_along), sep = "\r")
    }), use.names = FALSE)
    if (anyDuplicated(.key)) {
      .grp    <- match(.key, unique(.key))
      .merged <- combined_pep_df[!duplicated(.grp), , drop = FALSE]
      for (col in target_cols) {
        v   <- combined_pep_df[[col]]
        ok  <- which(!is.na(v))
        out <- rep(NA_real_, nrow(.merged))
        out[.grp[ok]] <- v[ok]
        .merged[[col]] <- out
      }
      n_pep_merged_across_files <- nrow(combined_pep_df) - nrow(.merged)
      combined_pep_df <- .merged
      rownames(combined_pep_df) <- NULL
      logger(paste0("tmt_limma_qf_summarize: ", n_pep_merged_across_files,
                    " row(s) of peptides repeated across files combined; ",
                    nrow(combined_pep_df), " peptide features"))
    }
  }

  # -- Parsimony (razor) assignment --
  # Runs on the POOLED table, not per file, so a peptide cannot be razored one
  # way in plex 1 and another in plex 2 and leave the protein spuriously
  # single-plex. Uses the untruncated MasterAccessions column carried through
  # the rbind above.
  n_razor_assigned <- 0L
  n_razor_dropped  <- 0L
  if (identical(shared_peptides, "razor")) {
    .pool_flags <- tmt_flag_shared_peptides(combined_pep_df, "MasterAccessions")
    .razor      <- tmt_razor_assign_shared(combined_pep_df$Peptide,
                                           .pool_flags$candidates)
    n_razor_assigned <- sum(.razor$assigned)
    n_razor_dropped  <- sum(.razor$dropped)
    # Reassigned rows only: unique rows keep the Protein the loop gave them.
    .reassign <- which(.razor$assigned)
    if (length(.reassign))
      combined_pep_df$Protein[.reassign] <- .razor$protein[.reassign]
    if (any(.razor$dropped))
      combined_pep_df <- combined_pep_df[!.razor$dropped, , drop = FALSE]
    logger(paste0("tmt_limma_qf_summarize: parsimony (razor) assignment on the ",
                  "pooled table: ", n_razor_assigned, " shared row(s) assigned by ",
                  "unique-peptide support, ", n_razor_dropped,
                  " dropped for lack of unique support; ",
                  nrow(combined_pep_df), " rows remain"))
    if (nrow(combined_pep_df) == 0L)
      return(list(error = paste0(
        "tmt_limma_qf_summarize: no peptide data available after parsimony ",
        "(razor) assignment - no candidate protein had any unique-peptide ",
        "support. Re-run with 'Keep all peptides' to include them.")))
  }

  # ── Delegate to qf_summarize_to_protein_matrix ──
  qf_result  <- qf_summarize_to_protein_matrix(
    input_df                = combined_pep_df,
    input_type              = "flat_peptide",
    protein_col             = "Protein",
    peptide_col             = "Peptide",
    intensity_cols          = target_cols,
    sample_names            = target_cols,
    normalize_peptide_level = .center,   # normalize before aggregation (QFeatures vignette order)
    logger                  = logger
  )

  if (!is.null(qf_result$error))
    return(list(error = paste0("tmt_limma_qf_summarize: ", qf_result$error)))

  protein_matrix   <- qf_result$protein_matrix
  zero_origin_mask <- qf_result$zero_origin_mask
  # qf_summarize already applied log2 and set non-finite to NA; defense-in-depth:
  protein_matrix[!is.finite(protein_matrix)] <- NA

  # ── Post-aggregation normalization ──
  # Peptide-level centering was already applied in qf_summarize_to_protein_matrix
  # (normalize_peptide_level=TRUE), which is the QFeatures-vignette-correct
  # position: normalise features, THEN aggregate. The pipeline used to centre a
  # second time here, at protein level, which is the same correction applied
  # twice. That second pass has been removed.
  #
  # The shift it WOULD have applied is still computed and logged. It is the
  # evidence that removal was safe and the number to check against real data: a
  # large value means aggregation is introducing channel offsets that
  # peptide-level centering did not remove, which is worth investigating rather
  # than silently correcting. On mock data (1000 proteins, 3 peptides each,
  # 10-20 channels) it sits at 0.03-0.09 log2 units.
  sample_medians <- apply(protein_matrix, 2, median, na.rm = TRUE)
  global_median  <- median(sample_medians, na.rm = TRUE)
  residual_column_shift_max <- suppressWarnings(
    max(abs(global_median - sample_medians), na.rm = TRUE))
  if (!is.finite(residual_column_shift_max)) residual_column_shift_max <- NA_real_
  logger(paste0("tmt_limma_qf_summarize: residual protein-level column shift ",
                "(computed, NOT applied): max |global median - column median| = ",
                signif(residual_column_shift_max, 4), " log2 units."))
  norm_applied <- if (.center) "peptide-level column median centering (log2), applied before aggregation"
                  else "none in OSMSS (input reported as already normalized across channels)"
  if (!.center)
    logger("tmt_limma_qf_summarize: peptide-level median centering NOT applied (input reported as already normalized across channels).")

  if (identical(as.integer(design_opt), 2L)) {
    protein_matrix <- apply_irs_normalization(protein_matrix, annotation_map, logger)
    norm_applied   <- paste0(norm_applied, " followed by IRS bridge-channel alignment")
  }

  # ── Cross-plex correction on opt 3: deliberately NONE at the matrix level ──
  # Do not "fix" the apparent gap here. Sample-loading (median) centering plus a
  # plex term in the design matrix IS the field-standard no-bridge treatment,
  # and the pipeline already has both. Adding a pseudo-reference IRS on top
  # would apply the same per-protein additive plex correction twice - once in
  # normalization, once via the covariate - leaving the covariate to estimate
  # approximately zero while the pre-adjustment silently consumes degrees of
  # freedom the model never accounts for.
  #
  # The real gap is that non-model consumers (heatmaps, trend clusters, exports)
  # see uncorrected plex offsets. That is what display_matrix below is for:
  # batch term in the model for inference, removeBatchEffect for display only.

  # ── Display-only batch-corrected matrix (>=2 plexes) ──
  # Mirrors make_pca_thumbnail(): same batch factor, same condition design, same
  # tryCatch-to-uncorrected fallback. NEVER pass this to lmFit - the plex term
  # in the model is what removes the offset for inference, and correcting twice
  # would consume degrees of freedom the model does not know about.
  display_matrix <- protein_matrix
  display_matrix_batch_corrected <- FALSE
  .mix_vec <- annotation_map$Mixture[match(colnames(protein_matrix),
                                           annotation_map$ColumnName)]
  .n_plex_disp <- length(unique(.mix_vec[!is.na(.mix_vec)]))
  if (.n_plex_disp >= 2L && requireNamespace("limma", quietly = TRUE)) {
    .cond_vec <- annotation_map$Condition[match(colnames(protein_matrix),
                                                annotation_map$ColumnName)]
    .cond_vec[is.na(.cond_vec) | !nzchar(.cond_vec)] <- "Unknown"
    display_matrix <- tryCatch({
      .dm <- limma::removeBatchEffect(protein_matrix,
                                      batch  = factor(.mix_vec),
                                      design = model.matrix(~ factor(.cond_vec)))
      display_matrix_batch_corrected <- TRUE
      .dm
    }, error = function(e) {
      logger(paste0("tmt_limma_qf_summarize: display removeBatchEffect failed (",
                    e$message, ") - display matrix left uncorrected."))
      protein_matrix
    })
    logger(paste0("tmt_limma_qf_summarize: display-only batch-corrected matrix ",
                  if (isTRUE(display_matrix_batch_corrected)) "built" else "NOT built",
                  " (", .n_plex_disp, " plexes). Analysis matrix is unchanged."))
  }

  pep_count <- qf_result$pep_count_per_protein

  # ── Build gene name lookup from protein-level PD data ──
  # qf_result$acc_to_gene is derived from peptide rows which have no Gene Symbol
  # or GN= information. Pull gene names from the protein-level table in each PD
  # file instead, using the same strategy as tmt_limma_protein_summarize.
  gene_lookup <- character(0)
  for (fi in seq_along(files)) {
    f       <- files[[fi]]
    pdata_p <- tryCatch(parse_pd_protein_data(f, rv), error = function(e) NULL)
    if (!is.null(pdata_p)) {
      for (k in seq_along(pdata_p$accessions)) {
        acc <- pdata_p$accessions[k]
        gn  <- if (length(pdata_p$gene_names) >= k) pdata_p$gene_names[k] else acc
        if (!is.na(acc) && nzchar(acc) && !(acc %in% names(gene_lookup)) &&
            !is.na(gn) && nzchar(gn)) {
          gene_lookup[acc] <- gn
        }
      }
    }
  }
  if (length(gene_lookup) == 0) {
    # Fallback: QF-derived acc_to_gene (accessions only) or protein matrix rownames
    gene_lookup <- qf_result$acc_to_gene %||%
                   setNames(rownames(protein_matrix), rownames(protein_matrix))
  }
  # Fill in any proteins from QF aggregation not present in protein-level table
  missing_gn <- setdiff(rownames(protein_matrix), names(gene_lookup))
  for (acc in missing_gn) gene_lookup[acc] <- acc

  conds <- unique(annotation_map$Condition[annotation_map$Condition != "Norm" &
                                           !is.na(annotation_map$Condition)])

  logger(paste0("tmt_limma_qf_summarize: ", nrow(protein_matrix), " proteins x ",
                ncol(protein_matrix), " channels"))

  list(
    protein_matrix         = protein_matrix,
    zero_origin_mask       = zero_origin_mask,
    gene_names             = gene_lookup,
    conditions             = sort(conds),
    annotation_map         = annotation_map,
    design_option          = as.integer(design_opt),
    n_proteins             = nrow(protein_matrix),
    n_channels             = ncol(protein_matrix),
    n_peptides_input       = qf_result$n_peptides_input,
    n_proteins_output      = qf_result$n_proteins_output,
    interference_filtered  = interference_col_found,
    interference_threshold = interference_threshold,
    pep_count_per_protein  = pep_count,
    # Same key name as tmt_limma_protein_summarize, so the two summarization
    # routes are symmetric. Provenance only - no consumer reads it yet.
    normalization          = norm_applied,
    normalization_applied  = .center,
    residual_column_shift_max = residual_column_shift_max,
    # Summarization estimator provenance, carried up from qf_summarize so the
    # DE-side provenance can state how many groups used which estimator.
    n_agg_fallback_colmeans = qf_result$n_agg_fallback_colmeans,
    n_agg_singleton         = qf_result$n_agg_singleton,
    peptide_fdr_source      = qf_result$peptide_fdr_source,
    n_interf_kept_na        = n_interf_kept_na,
    n_interf_kept_pass      = n_interf_kept_pass,
    # Rows of peptides repeated across per-plex files, combined into one feature.
    n_peptides_merged_across_files = n_pep_merged_across_files,
    # Shared-peptide provenance. The methods text reads these; a run whose
    # setting is absent predates this feature and did what "keep_all" does.
    shared_peptide_handling   = shared_peptides,
    n_peptides_total          = n_pep_total,
    n_peptides_shared         = n_pep_shared,
    n_peptides_removed_shared = if (identical(shared_peptides, "unique_only"))
                                  n_pep_shared else 0L,
    shared_detection_method   = shared_detection_method %||% NA_character_,
    # Proteins present under keep_all and absent under unique_only. The exact
    # quantity is a difference of two protein-matrix row sets, so it needs a
    # second summarization pass; that is not worth doubling the runtime, and a
    # peptide-table proxy would be a different number under this key's name.
    n_proteins_lost_to_shared_filter = if (identical(shared_peptides, "keep_all"))
                                         0L else NA_integer_,
    # Parsimony (razor) provenance. Zero unless shared_peptides == "razor".
    n_peptides_razor_assigned             = n_razor_assigned,
    n_peptides_dropped_no_unique_evidence = n_razor_dropped,
    razor_basis = if (identical(shared_peptides, "razor"))
                    "unique_peptide_count" else NA_character_,
    # DISPLAY ONLY. Plex offsets removed with limma::removeBatchEffect so that
    # heatmaps, trend clusters and exports do not show uncorrected between-plex
    # steps. Must NEVER be passed to lmFit or any other model fit: inference
    # uses protein_matrix with the plex term in the design.
    display_matrix                 = display_matrix,
    display_matrix_batch_corrected = display_matrix_batch_corrected,
    cross_plex_normalization       = switch(as.character(as.integer(design_opt)),
                                            "1" = "none_single_plex",
                                            "2" = "irs_bridge",
                                            "3" = "model_covariate",
                                            "unknown"),
    route                  = "limma_qf"
  )
}

# === INITIAL LOAD =============================================================
# TKL matrices are now deferred  -" loaded on first use via ensure_global_resources()
message(">>> SYSTEM: TKL Matrix loading deferred to first use (lazy-load).")

# === ASYNC OVERLAY HELPERS ====================================================
# Called from module servers to show/hide the full-page blocking overlay while
# an ExtendedTask runs in a background worker. Requires shinyjs (loaded above).

show_async_overlay <- function(session, label = "Processing...", detail = "") {
  shinyjs::html("async_task_label", label)
  shinyjs::html("async_task_detail", detail)
  shinyjs::show("async_task_overlay")
}

hide_async_overlay <- function(session) {
  shinyjs::hide("async_task_overlay")
  shinyjs::html("async_task_label", "Processing...")
  shinyjs::html("async_task_detail", "")
}

# === ASYNC COMPUTE FUNCTIONS ==================================================
# Pure (side-effect-free) compute functions called from ExtendedTask workers.
# They accept plain-R arguments (no reactiveValues) and return named lists.

# compute_proda_de -------------------------------------------------------------
# Runs the full ProDA probabilistic-dropout DE pipeline on a pre-built protein
# matrix cache. Returns list(success, error, result_df, plot_data, acc_to_gene,
# provenance, messages). Called from a future_promise() inside proda_de_task.
compute_proda_de <- function(
    cache,
    ctrl_label,
    cond_label,
    cond_groups,
    factorial_design,
    batch_for_model,
    batch_groups,
    batch_filter_enabled,
    subject_pairing,
    expect_mnar,
    ip_sub,
    igg_conditions   = character(0),
    logger_fn        = message
) {
  library(proDA)
  library(QFeatures)

  msgs <- character(0)   # warning/info messages to replay in main thread

  # proDA initialises missing cells with rnorm() and sets no seed of its own
  # (proDA.R: Y_compl[is.na(Y)] <- rnorm(...)), so the same input returns
  # different p-values on every run. Declared as a local, not a package-level
  # constant, so it is captured when this function is sent to a future worker.
  .proda_seed <- 42L

  design_df      <- cache$design_df
  protein_matrix <- cache$protein_matrix
  acc_to_gene    <- cache$acc_to_gene

  sanitize_cond <- function(x) gsub("[^A-Za-z0-9_.]", "_", x)

  ctrl_safe <- sanitize_cond(ctrl_label)
  cond_safe <- sanitize_cond(cond_label)

  design_df$condition <- sanitize_cond(design_df$condition)
  all_conds_safe      <- unique(design_df$condition)

  other_conds <- setdiff(all_conds_safe, ctrl_safe)
  design_df$condition <- factor(design_df$condition,
                                levels = c(ctrl_safe, other_conds))
  all_samps  <- design_df$sample_name[design_df$sample_name %in% colnames(protein_matrix)]
  design_sub <- design_df[design_df$sample_name %in% all_samps, , drop = FALSE]
  rownames(design_sub) <- design_sub$sample_name
  mat_sub    <- protein_matrix[, all_samps, drop = FALSE]

  if (!ctrl_safe %in% levels(design_sub$condition)) {
    return(list(success = FALSE,
                error = paste0("Control '", ctrl_label, "' not found after sanitization."),
                messages = msgs))
  }
  if (!cond_safe %in% levels(design_sub$condition)) {
    return(list(success = FALSE,
                error = paste0("Condition '", cond_label, "' not found after sanitization."),
                messages = msgs))
  }
  if (ncol(mat_sub) < 2) {
    return(list(success = FALSE,
                error = "Need at least 2 samples for ProDA DE.",
                messages = msgs))
  }

  # ── Pairwise subset (AP/co-IP only) ───────────────────────────────────────
  # For AP/co-IP experiments (IgG conditions are defined), restrict ProDA to
  # only the two compared conditions before fitting. In these experiments the
  # IgG and bait groups are so biologically heterogeneous that fitting on all
  # conditions simultaneously inflates the EB variance prior, crushing every
  # t-statistic and pushing all adj. p-values toward 1.
  # For global proteomics (igg_conditions is empty) we keep all conditions in
  # the model so the pooled EB variance prior is estimated correctly.
  pairwise_only <- length(igg_conditions) > 0 &&
                   (is.null(factorial_design) || !isTRUE(factorial_design$enabled) ||
                    length(factorial_design$assignments) == 0 ||
                    !nzchar(factorial_design$f1_name %||% "") ||
                    !nzchar(factorial_design$f2_name %||% ""))
  .variance_prior_scope  <- "all_conditions"
  .n_conds_dropped_model <- 0L
  if (pairwise_only) {
    keep_samps <- design_sub$sample_name[
      as.character(design_sub$condition) %in% c(ctrl_safe, cond_safe)]
    if (length(keep_samps) >= 2 && length(keep_samps) < ncol(mat_sub)) {
      n_dropped_conds <- length(unique(as.character(design_sub$condition))) - 2L
      logger_fn(paste0("[ProDA] Pairwise subset: restricting to ",
                       ctrl_label, " + ", cond_label, " (",
                       length(keep_samps), " samples); dropped ",
                       n_dropped_conds, " other condition(s) from model"))
      design_sub           <- design_sub[design_sub$sample_name %in% keep_samps, , drop = FALSE]
      design_sub$condition <- droplevels(design_sub$condition)
      mat_sub              <- mat_sub[, keep_samps, drop = FALSE]
      .variance_prior_scope  <- "pairwise_subset"
      .n_conds_dropped_model <- as.integer(n_dropped_conds)

    }
  }

  # \u2500\u2500 Contrast-level minimum evidence \u2500\u2500
  # Summarization keeps anything observed in >=1 sample because the contrast is
  # unknown at that stage. After the pairwise column subset, proteins observed
  # only in a dropped condition are entirely-NA rows carrying no information
  # about this contrast. proDA is built to model incomplete rows, so only rows
  # with nothing at all in either compared arm are removed; sparse-but-fittable
  # rows are kept and flagged instead.
  .ev_ctrl_samps <- design_sub$sample_name[
    as.character(design_sub$condition) == ctrl_safe]
  .ev_cond_samps <- design_sub$sample_name[
    as.character(design_sub$condition) == cond_safe]
  .n_obs_ctrl <- if (length(.ev_ctrl_samps) > 0) {
    rowSums(!is.na(mat_sub[, .ev_ctrl_samps, drop = FALSE]))
  } else rep(0L, nrow(mat_sub))
  .n_obs_cond <- if (length(.ev_cond_samps) > 0) {
    rowSums(!is.na(mat_sub[, .ev_cond_samps, drop = FALSE]))
  } else rep(0L, nrow(mat_sub))
  # Keyed by protein name, not by position: the drop below changes both the
  # length and the order of mat_sub relative to proda_res.
  # LowEvidence = fewer than two observations in AT LEAST ONE compared arm
  # (pmin). The earlier rule (pmax) required fewer than two in BOTH arms and
  # missed proteins absent from one arm, whose fold changes rest almost
  # entirely on the dropout model. Display/export flag only: it changes no
  # fit, p-value, fold change or score.
  .low_evidence_by_prot <- setNames(pmin(.n_obs_ctrl, .n_obs_cond) < 2L,
                                    rownames(mat_sub))
  .contrast_empty <- (.n_obs_ctrl + .n_obs_cond) == 0
  .n_contrast_empty_removed <- as.integer(sum(.contrast_empty))
  if (.n_contrast_empty_removed > 0L) {
    mat_sub <- mat_sub[!.contrast_empty, , drop = FALSE]
    logger_fn(paste0("[ProDA] Contrast minimum-evidence filter: removed ",
                     .n_contrast_empty_removed,
                     " protein(s) with no observation in either ", ctrl_label,
                     " or ", cond_label, "; ", nrow(mat_sub), " retained"))
  }

  logger_fn(paste0("[ProDA] Fitting model: ", nrow(mat_sub), " proteins \u00d7 ",
                   ncol(mat_sub), " samples (max_iter=20, verbose=TRUE)"))

  .proda_last_iter     <- 0L
  .proda_fit_strategy  <- "full"
  .proda_max_iter_used <- 20L
  .proda_df_corrected  <- FALSE
  .proda_hessian_events <- 0L

  # Local fit helper  -" no setProgress in worker context
  .proda_fit <- function(mat, col_data, design = ~ condition,
                         moderate_loc = TRUE, moderate_var = TRUE,
                         max_iter = 20L) {
    set.seed(.proda_seed)
    .proda_hessian_events <<- 0L
    withCallingHandlers(
      proDA::proDA(
        data                    = mat,
        design                  = design,
        col_data                = col_data,
        data_is_log_transformed = TRUE,
        moderate_location       = moderate_loc,
        moderate_variance       = moderate_var,
        max_iter                = max_iter,
        verbose                 = TRUE
      ),
      message = function(m) {
        msg <- trimws(conditionMessage(m))
        logger_fn(paste0("[ProDA] ", msg))
        if (grepl("^Starting iter:", msg)) {
          iter_num <- suppressWarnings(
            as.integer(sub("^Starting iter:\\s*", "", msg)))
          if (!is.na(iter_num)) .proda_last_iter <<- iter_num
        }
        invokeRestart("muffleMessage")
      },
      # proDA's invert_hessian_matrix() catches a failed solve() of one protein's
      # coefficient Hessian, gives that protein infinite coefficient variance for
      # the step and carries on. With verbose = TRUE it then reports the caught
      # error with warning(err); because err is an error object, the condition is
      # of class "error" and the tryCatch(error = ...) around every .proda_fit()
      # call used to treat it as a failed fit, dropping to the early-stop fallback
      # (max_iter = 1 on HYEtims735). Only a condition raised through warning()
      # offers a "muffleWarning" restart, so this handler absorbs exactly those
      # reports; a genuine stop() inside proDA has no such restart and still
      # propagates to the caller's tryCatch.
      error = function(e) {
        if (!is.null(findRestart("muffleWarning"))) {
          .proda_hessian_events <<- .proda_hessian_events + 1L
          if (.proda_hessian_events <= 3L)
            logger_fn(paste0("[ProDA] Singular coefficient Hessian at iter ",
                             .proda_last_iter, " handled inside proDA (",
                             conditionMessage(e), "); fit continues"))
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  fit <- NULL

  # ── Factorial design setup ─────────────────────────────────────────────────
  fd           <- factorial_design
  is_factorial <- !is.null(fd) && isTRUE(fd$enabled) &&
                  length(fd$assignments) > 0 &&
                  nzchar(fd$f1_name %||% "") && nzchar(fd$f2_name %||% "") &&
                  length(fd$f1_levels) >= 2 && length(fd$f2_levels) >= 2

  if (is_factorial) {
    f1_name   <- fd$f1_name
    f2_name   <- fd$f2_name
    f1_levels <- fd$f1_levels
    f2_levels <- fd$f2_levels
    f1_ref    <- sanitize_cond(f1_levels[1])
    f1_treat  <- sanitize_cond(f1_levels[length(f1_levels)])
    assigns_safe <- setNames(fd$assignments, sanitize_cond(names(fd$assignments)))
    design_sub[[f1_name]] <- factor(
      sapply(as.character(design_sub$condition), function(g) {
        a <- assigns_safe[[g]]
        if (is.null(a) || !nzchar(a$f1)) NA_character_ else sanitize_cond(a$f1)
      }),
      levels = c(f1_ref, setdiff(sanitize_cond(f1_levels), f1_ref))
    )
    design_sub[[f2_name]] <- factor(
      sapply(as.character(design_sub$condition), function(g) {
        a <- assigns_safe[[g]]
        if (is.null(a) || !nzchar(a$f2)) NA_character_ else sanitize_cond(a$f2)
      }),
      levels = sanitize_cond(f2_levels)
    )
    n_missing_f1 <- sum(is.na(design_sub[[f1_name]]))
    n_missing_f2 <- sum(is.na(design_sub[[f2_name]]))
    if (n_missing_f1 > 0 || n_missing_f2 > 0) {
      logger_fn(paste0("[ProDA] Factorial: ", n_missing_f1 + n_missing_f2,
                       " samples have incomplete factor assignments \u2014 falling back to one-way model"))
      msgs <- c(msgs, paste0("Factorial design: ", n_missing_f1 + n_missing_f2,
                              " group(s) have incomplete factor assignments. Falling back to one-way model."))
      is_factorial <- FALSE
    } else {
      logger_fn(paste0("[ProDA] Factorial model: ", f1_name, " [",
                       paste(levels(design_sub[[f1_name]]), collapse = "/"), "] ~ ",
                       f2_name, " [", paste(levels(design_sub[[f2_name]]), collapse = "/"), "]"))
    }
  }

  # ── Batch covariate ────────────────────────────────────────────────────────
  batch_in_model <- FALSE
  if (isTRUE(batch_for_model) && !is.null(batch_groups)) {
    bg               <- batch_groups
    batch_vec        <- rep(NA_character_, nrow(design_sub))
    names(batch_vec) <- design_sub$sample_name
    for (bl in names(bg)) {
      batch_vec[names(batch_vec) %in% bg[[bl]]] <- bl
    }
    if (!any(is.na(batch_vec))) {
      batch_vec_safe   <- gsub("[^A-Za-z0-9_.]", "_", batch_vec)
      design_sub$batch <- factor(batch_vec_safe)
      bc_tbl     <- table(design_sub$batch, design_sub$condition)
      confounded <- any(colSums(bc_tbl > 0L) < 2L)
      if (confounded) {
        logger_fn("[ProDA] Batch confounded with condition \u2014 removing from model")
        msgs <- c(msgs, "Batch is confounded with condition; excluded from ProDA model.")
        design_sub$batch <- NULL
      } else {
        batch_in_model <- TRUE
        logger_fn(paste0("[ProDA] Batch added to model: ",
                         paste(levels(design_sub$batch), collapse = " / ")))
      }
    } else {
      logger_fn("[ProDA] Batch vector has NAs for some samples \u2014 batch not added to model")
    }
  }

  use_interaction <- is_factorial && isTRUE(fd$interaction)

  # ── Batch-spanning protein filter ──────────────────────────────────────────
  if (batch_in_model && isTRUE(batch_filter_enabled)) {
    cell_labels <- paste0(as.character(design_sub$batch), ":",
                          as.character(design_sub$condition))
    names(cell_labels) <- design_sub$sample_name
    unique_cells <- unique(cell_labels)
    n_cells      <- length(unique_cells)
    keep_prot <- apply(mat_sub, 1L, function(row_vals) {
      valid_cells <- cell_labels[!is.na(row_vals[names(cell_labels)])]
      length(unique(valid_cells)) >= n_cells
    })
    n_dropped <- sum(!keep_prot)
    if (n_dropped > 0L) {
      logger_fn(paste0("[ProDA] Batch-spanning filter: removed ", n_dropped,
                       " proteins missing in \u22651 batch\u00d7condition cell (",
                       sum(keep_prot), " remain, ", n_cells, " cells)"))
      mat_sub <- mat_sub[keep_prot, , drop = FALSE]
    } else {
      logger_fn("[ProDA] Batch-spanning filter: all proteins pass (no removals)")
    }
  }

  # ── Interaction rank check ─────────────────────────────────────────────────
  if (use_interaction && batch_in_model) {
    logger_fn("[ProDA] Interaction+batch design check \u2014 sample table:")
    for (i in seq_len(nrow(design_sub))) {
      r <- design_sub[i, ]
      logger_fn(paste0("  ", r$sample_name,
                       "  cond=", r$condition,
                       "  ", fd$f1_name, "=", r[[fd$f1_name]],
                       "  ", fd$f2_name, "=", r[[fd$f2_name]],
                       "  batch=", r$batch))
    }
    mm_test <- tryCatch(
      model.matrix(
        as.formula(paste("~", fd$f1_name, "*", fd$f2_name, "+ batch")),
        data = design_sub
      ),
      error = function(e) {
        logger_fn(paste0("[ProDA] model.matrix error: ", e$message)); NULL
      }
    )
    if (!is.null(mm_test)) {
      logger_fn(paste0("[ProDA] Model matrix columns: ",
                       paste(colnames(mm_test), collapse = ", ")))
      cor_mat <- tryCatch(cor(mm_test[, -1L, drop = FALSE]), error = function(e) NULL)
      if (!is.null(cor_mat)) {
        pairs <- which(abs(cor_mat) > 0.999 & upper.tri(cor_mat), arr.ind = TRUE)
        if (nrow(pairs) > 0) {
          for (p in seq_len(nrow(pairs))) {
            logger_fn(paste0("[ProDA] Near-perfect collinearity: '",
                             colnames(mm_test[-1L][pairs[p, 1]]), "' ~ '",
                             colnames(mm_test[-1L][pairs[p, 2]]), "'"))
          }
        } else {
          logger_fn("[ProDA] No pairwise collinearity detected (>0.999)")
        }
      }
      rcond <- tryCatch(1 / kappa(mm_test, exact = TRUE), error = function(e) 0)
      logger_fn(paste0("[ProDA] Design matrix rcond=", signif(rcond, 4),
                       "  (threshold=", signif(.Machine$double.eps * 10, 3), ")"))
      if (rcond < .Machine$double.eps * 10) {
        return(list(
          success  = FALSE,
          error    = paste0(
            "Cannot fit interaction model with current batch assignment. ",
            "The batch column is linearly dependent on the ",
            fd$f1_name, ":", fd$f2_name,
            " interaction term. Try the additive model or reassign batches."),
          messages = msgs))
      }
    }
  }

  # ── Build design formula ───────────────────────────────────────────────────
  subject_in_model <- FALSE   # initialise before conditional (may stay FALSE for factorial path)
  design_formula <- if (is_factorial) {
    op <- if (use_interaction) "*" else "+"
    as.formula(paste("~",
                     if (batch_in_model) "batch + " else "",
                     fd$f1_name, op, fd$f2_name))
  } else {
    if (!is.null(subject_pairing) && length(subject_pairing) > 0) {
      subj_vec <- setNames(character(nrow(design_sub)), design_sub$sample_name)
      for (sn in names(subj_vec)) {
        if (sn %in% names(subject_pairing) && nzchar(subject_pairing[[sn]]))
          subj_vec[[sn]] <- subject_pairing[[sn]]
      }
      if (!any(!nzchar(subj_vec))) {
        design_sub$subject <- factor(gsub("[^A-Za-z0-9_.]", "_", subj_vec))
        sc_tbl <- table(design_sub$subject, design_sub$condition)
        if (any(rowSums(sc_tbl > 0L) < 2L)) {
          logger_fn("[ProDA] Subject confounded with condition \u2014 removing from model")
          msgs <- c(msgs, "Subject pairing is confounded with condition; excluded from ProDA model.")
          design_sub$subject <- NULL
        } else {
          subject_in_model <- TRUE
          logger_fn(paste0("[ProDA] Subject added to model: ",
                           paste(levels(design_sub$subject), collapse = " / ")))
        }
      } else {
        logger_fn("[ProDA] Subject pairing has blanks for some samples \u2014 not added to model")
      }
    }
    if (batch_in_model && subject_in_model) ~ batch + subject + condition
    else if (batch_in_model)               ~ batch + condition
    else if (subject_in_model)             ~ subject + condition
    else                                   ~ condition
  }
  logger_fn(paste0("[ProDA] Formula: ", deparse(design_formula),
                   if (is_factorial) paste0(" (",
                     if (use_interaction) "interaction" else "additive", ")") else ""))

  # ── Fit (phase 1: full moderated) ─────────────────────────────────────────
  fit <- tryCatch(
    .proda_fit(mat_sub, design_sub, design = design_formula,
               moderate_loc = TRUE, moderate_var = TRUE),
    error = function(e) {
      logger_fn(paste0("[ProDA] Full model failed at iter ", .proda_last_iter,
                       " (", e$message, ")")); NULL
    }
  )

  # Phase 2: descending early-stop
  if (is.null(fit)) {
    crash_iter <- .proda_last_iter
    for (try_iter in seq(max(crash_iter - 2L, 1L), 1L, by = -1L)) {
      logger_fn(paste0("[ProDA] Retrying moderated model with max_iter=", try_iter))
      .proda_last_iter     <- 0L
      .proda_fit_strategy  <- "early_stop"
      .proda_max_iter_used <- try_iter
      fit <- tryCatch(
        .proda_fit(mat_sub, design_sub, design = design_formula,
                   moderate_loc = TRUE, moderate_var = TRUE, max_iter = try_iter),
        error = function(e) {
          logger_fn(paste0("[ProDA] max_iter=", try_iter, " failed (", e$message, ")")); NULL
        }
      )
      if (!is.null(fit)) {
        msgs <- c(msgs, paste0(
          "The proDA model did not converge at the requested 20 iterations. ",
          "Results come from an early-stopped fit at max_iter = ", try_iter, "."))
        break
      }
    }
  }

  # Phase 3: unmoderated
  if (is.null(fit)) {
    logger_fn("[ProDA] All moderated attempts failed \u2014 retrying with moderation disabled")
    .proda_last_iter     <- 0L
    .proda_fit_strategy  <- "unmoderated"
    .proda_max_iter_used <- 20L
    fit <- tryCatch(
      .proda_fit(mat_sub, design_sub, design = design_formula,
                 moderate_loc = FALSE, moderate_var = FALSE),
      error = function(e) {
        logger_fn(paste0("[ProDA] Unmoderated model also failed: ", e$message)); NULL
      }
    )
    if (!is.null(fit))
      msgs <- c(msgs, paste0(
        "Variance and location moderation were disabled for this proDA fit ",
        "because every moderated attempt failed."))
  }

  if (is.null(fit)) {
    err_msg <- if (isTRUE(batch_in_model)) {
      "Statistical model fitting failed (singular matrix). The batch column may be linearly dependent on condition. Try reassigning batches or removing batch correction."
    } else if (isTRUE(use_interaction)) {
      paste0("The ", fd$f1_name, " \u00d7 ", fd$f2_name,
             " interaction model could not be fit. The design matrix is rank-deficient. Try re-running with simple conditions.")
    } else {
      "ProDA could not converge. Check that each condition has \u22652 replicates with valid values."
    }
    return(list(success = FALSE, error = err_msg, messages = msgs))
  }

  # proDA::convergence() returns list(successful, iterations, error); there are
  # no 'converged' or 'delta' fields, so reading those always gave FALSE / empty.
  conv <- tryCatch(proDA::convergence(fit), error = function(e) NULL)
  if (!is.null(conv)) {
    conv_delta <- tryCatch(signif(as.numeric(conv$error), 4), error = function(e) NA_real_)
    logger_fn(paste0("[ProDA] Converged=", isTRUE(conv$successful),
                     "  iterations=", conv$iterations %||% NA,
                     "  final_delta=", conv_delta,
                     "  singular_hessian_events_handled=", .proda_hessian_events))
    if (!isTRUE(conv$successful))
      msgs <- c(msgs, paste0(
        "The proDA model reported that it had not converged after ",
        conv$iterations %||% NA, " iterations. Moderated variance estimates ",
        "and the p-values derived from them should be treated with caution."))
  }

  # ── Contrast string ────────────────────────────────────────────────────────
  contrast_str <- if (is_factorial) {
    ctrl_f1     <- sanitize_cond(assigns_safe[[ctrl_safe]]$f1 %||% "")
    ctrl_f2     <- sanitize_cond(assigns_safe[[ctrl_safe]]$f2 %||% "")
    cond_f1     <- sanitize_cond(assigns_safe[[cond_safe]]$f1 %||% "")
    cond_f2     <- sanitize_cond(assigns_safe[[cond_safe]]$f2 %||% "")
    f2_ref_safe <- sanitize_cond(f2_levels[1])
    cell_terms <- function(f1v, f2v, sign_char) {
      terms <- character(0)
      f1c <- if (!identical(f1v, f1_ref))      paste0(f1_name, f1v)  else NULL
      f2c <- if (!identical(f2v, f2_ref_safe)) paste0(f2_name, f2v)  else NULL
      ic  <- if (use_interaction && !is.null(f1c) && !is.null(f2c))
               paste0("`", f1_name, f1v, ":", f2_name, f2v, "`") else NULL
      for (t in c(f1c, f2c, ic)) {
        if (!is.null(t)) terms <- c(terms, paste0(sign_char, t))
      }
      terms
    }
    all_terms <- c(cell_terms(cond_f1, cond_f2, ""),
                   cell_terms(ctrl_f1, ctrl_f2, "-"))
    coef_map <- list()
    for (t in all_terms) {
      sign <- if (startsWith(t, "-")) -1L else 1L
      nm   <- sub("^-", "", t)
      coef_map[[nm]] <- (coef_map[[nm]] %||% 0L) + sign
    }
    final_coef <- Filter(function(v) v != 0L, coef_map)
    if (length(final_coef) == 0) {
      return(list(
        success  = FALSE,
        error    = paste0("'", ctrl_label, "' and '", cond_label,
                          "' map to identical model coefficients \u2014 no contrast to test."),
        messages = msgs))
    }
    parts <- mapply(function(nm, v) {
      if (v == 1L) nm
      else if (v == -1L) paste0("-", nm)
      else paste0(v, "*", nm)
    }, names(final_coef), final_coef, SIMPLIFY = TRUE)
    paste(parts, collapse = " + ")
  } else {
    paste0("condition", cond_safe)
  }
  logger_fn(paste0("[ProDA] contrast: ", contrast_str, " (",
                   if (is_factorial) paste0(
                     if (use_interaction) "interaction" else "additive",
                     " factorial \u2014 ", cond_label, " vs ", ctrl_label)
                   else "one-way", ")"))

  proda_res <- tryCatch(
    proDA::test_diff(fit, contrast = contrast_str),
    error = function(e) {
      logger_fn(paste0("[ProDA] test_diff error: ", e$message)); NULL
    }
  )
  if (is.null(proda_res) || nrow(proda_res) == 0) {
    return(list(success = FALSE,
                error   = "ProDA test_diff failed \u2014 see log.",
                messages = msgs))
  }
  if (all(is.na(proda_res$pval))) {
    logger_fn("[ProDA] test_diff returned all-NA p-values \u2014 model fit is unusable")
    return(list(success = FALSE,
                error   = "ProDA model produced no usable p-values \u2014 see log.",
                messages = msgs))
  }

  # ── Rescue NA proteins ─────────────────────────────────────────────────────
  # Which test produced each p-value. Starts as the proDA model Wald test and
  # is overwritten wherever a rescue or a default value took its place, so the
  # p-value column stops mixing procedures without saying so.
  test_source <- rep("proda_wald", nrow(proda_res))
  na_mask <- is.na(proda_res$diff)
  if (any(na_mask)) {
    rescue_ctrl <- design_sub$sample_name[as.character(design_sub$condition) == ctrl_safe]
    rescue_cond <- design_sub$sample_name[as.character(design_sub$condition) == cond_safe]
    n_rescued <- 0L
    for (i in which(na_mask)) {
      prot <- proda_res$name[i]
      if (prot %in% rownames(mat_sub)) {
        cv <- mat_sub[prot, rescue_ctrl]; cv <- cv[!is.na(cv)]
        dv <- mat_sub[prot, rescue_cond]; dv <- dv[!is.na(dv)]
        if (length(cv) >= 1 && length(dv) >= 1) {
          proda_res$diff[i]          <- mean(dv) - mean(cv)
          proda_res$avg_abundance[i] <- mean(c(cv, dv))
          test_source[i]             <- "fc_only"
          n_rescued <- n_rescued + 1L
        }
        if (length(cv) >= 2 && length(dv) >= 2) {
          tt <- tryCatch(t.test(dv, cv), error = function(e) NULL)
          if (!is.null(tt)) {
            proda_res$pval[i]        <- tt$p.value
            proda_res$t_statistic[i] <- unname(tt$statistic)
            test_source[i]           <- "welch_t"
          }
        }
      }
    }
    proda_res$adj_pval <- p.adjust(proda_res$pval, method = "BH")
    if (n_rescued > 0L)
      logger_fn(paste0("[ProDA] Rescued ", n_rescued, "/", sum(na_mask),
                       " NA-coefficient proteins via direct fold-change computation"))
  }

  # ── Empirical Bayes df correction ──────────────────────────────────────────
  # proDA::test_diff already returns p-values that account for the moderated
  # degrees of freedom. The moderated recomputation below is retained as a
  # DIAGNOSTIC only: it is no longer written back over proda_res$pval /
  # $adj_pval, because doing so applied a second, independent moderation on
  # top of the package's own, and did so conditionally, so two datasets with
  # identical settings could be analysed by different procedures.
  .proda_med_p_native     <- NA_real_
  .proda_med_p_moderated  <- NA_real_
  tryCatch({
    feat_params  <- proDA::feature_parameters(fit)
    mod_df       <- feat_params$df
    raw_df       <- nrow(proDA::design(fit)) - ncol(proDA::design(fit))
    t_stats      <- proda_res$t_statistic
    vp_df <- tryCatch(fit@variance_prior_df, error = function(e) NULL)
    # No arbitrary constant: if the variance prior df is unavailable the
    # moderated df is simply unknown, and is reported as such.
    df_cap <- if (!is.null(vp_df) && is.finite(vp_df) && vp_df > 0) {
      vp_df + raw_df
    } else {
      NA_real_
    }
    mod_df[is.infinite(mod_df)] <- df_cap
    mod_df[is.na(mod_df)]       <- raw_df
    .proda_med_p_native <- median(proda_res$pval, na.rm = TRUE)
    if (!is.null(mod_df) && !all(is.na(mod_df)) &&
        median(mod_df, na.rm = TRUE) > raw_df) {
      pval_mod <- 2 * pt(-abs(t_stats), df = mod_df)
      pval_mod[is.na(pval_mod)] <- 1
      adj_pval_mod <- p.adjust(pval_mod, method = "BH")
      .proda_med_p_moderated <- median(pval_mod, na.rm = TRUE)
      logger_fn(paste0("[ProDA] Wald df correction computed for diagnostics ",
                       "only; proDA p-values retained. raw df=", raw_df,
                       ", median moderated df=",
                       round(median(mod_df, na.rm = TRUE), 1),
                       " \u2192 median raw p ",
                       signif(median(proda_res$pval, na.rm = TRUE), 3),
                       " \u2192 ", signif(median(pval_mod, na.rm = TRUE), 3)))
      # proda_res$pval / $adj_pval are deliberately NOT overwritten. The
      # moderated df is legitimate metadata and is kept.
      proda_res$df       <- mod_df
      .proda_df_corrected <- FALSE
    } else {
      logger_fn(paste0("[ProDA] Wald df diagnostic not computed: moderated df not larger than raw df=",
                       raw_df))
    }
  }, error = function(e) {
    logger_fn(paste0("[ProDA] Wald df correction failed: ", e$message,
                     " \u2014 using native p-values"))
  })

  # ── Build result table ─────────────────────────────────────────────────────
  build_result <- tryCatch({
    n_res     <- nrow(proda_res)
    all_prots <- as.character(proda_res$name)

    # Flag, do not drop: proDA can legitimately fit a single-observation arm.
    # Indexed by name, never by position, per the alignment note above.
    low_evidence_vec <- unname(.low_evidence_by_prot[all_prots])
    low_evidence_vec[is.na(low_evidence_vec)] <- FALSE

    # Proteins the model could not test at all. Capture them before the
    # defaults below (p = 1, log2FC = 0) make them indistinguishable from
    # proteins that were tested and came out null.
    untested <- is.na(proda_res$pval) | is.na(proda_res$diff)
    test_source[untested] <- "untested"

    adj_p   <- as.numeric(proda_res$adj_pval);    adj_p[is.na(adj_p)]   <- 1
    raw_p   <- as.numeric(proda_res$pval);        raw_p[is.na(raw_p)]   <- 1
    log2fc  <- as.numeric(proda_res$diff);        log2fc[is.na(log2fc)] <- 0
    avg_int <- as.numeric(proda_res$avg_abundance); avg_int[is.na(avg_int)] <- 0

    ctrl_samps <- design_sub$sample_name[as.character(design_sub$condition) == ctrl_safe]
    cond_samps <- design_sub$sample_name[as.character(design_sub$condition) == cond_safe]
    n_ctrl     <- length(ctrl_samps)
    n_cond     <- length(cond_samps)

    ctrl_mat <- mat_sub[all_prots[all_prots %in% rownames(mat_sub)],
                        ctrl_samps, drop = FALSE]
    cond_mat <- mat_sub[all_prots[all_prots %in% rownames(mat_sub)],
                        cond_samps, drop = FALSE]

    ctrl_detect <- apply(ctrl_mat, 1, function(r) sum(!is.na(r))) / n_ctrl
    cond_detect <- apply(cond_mat, 1, function(r) sum(!is.na(r))) / n_cond
    ctrl_detect <- ctrl_detect[all_prots]; cond_detect <- cond_detect[all_prots]
    ctrl_detect[is.na(ctrl_detect)] <- 0;  cond_detect[is.na(cond_detect)] <- 0

    # Reported fold change is proDA's own test_diff()$diff on every row,
    # including rows where one group has no observations: proDA estimates that
    # group's mean from its dropout curve and returns a finite diff by
    # construction, so no floor substitution is applied (Pass F, F1/F2).
    display_fc <- log2fc
    finite_fc_vals <- display_fc[is.finite(display_fc)]
    fc_cap <- if (length(finite_fc_vals) > 0) max(abs(finite_fc_vals), na.rm = TRUE) + 1 else 15
    display_fc[display_fc ==  Inf] <-  fc_cap
    display_fc[display_fc == -Inf] <- -fc_cap
    display_fc[is.na(display_fc)]  <- 0

    # Per-replicate columns
    ctrl_rep_df <- as.data.frame(
      t(ctrl_mat[all_prots[all_prots %in% rownames(ctrl_mat)], , drop = FALSE]))
    cond_rep_df <- as.data.frame(
      t(cond_mat[all_prots[all_prots %in% rownames(cond_mat)], , drop = FALSE]))
    ctrl_rep_df <- as.data.frame(t(ctrl_rep_df))
    cond_rep_df <- as.data.frame(t(cond_rep_df))
    if (nrow(ctrl_rep_df) > 0) {
      ctrl_rep_df <- ctrl_rep_df[match(all_prots, rownames(ctrl_rep_df)), , drop = FALSE]
      for (si in seq_len(n_ctrl)) names(ctrl_rep_df)[si] <- paste0("Ctrl_Rep", si)
      ctrl_rep_df$Accession <- all_prots
      ctrl_rep_df[, grep("^Ctrl_", names(ctrl_rep_df))] <-
        round(ctrl_rep_df[, grep("^Ctrl_", names(ctrl_rep_df))], 2)
    }
    if (nrow(cond_rep_df) > 0) {
      cond_rep_df <- cond_rep_df[match(all_prots, rownames(cond_rep_df)), , drop = FALSE]
      for (si in seq_len(n_cond)) names(cond_rep_df)[si] <- paste0("Cond_Rep", si)
      cond_rep_df$Accession <- all_prots
      cond_rep_df[, grep("^Cond_", names(cond_rep_df))] <-
        round(cond_rep_df[, grep("^Cond_", names(cond_rep_df))], 2)
    }

    gene_names <- acc_to_gene[all_prots]
    gene_names[is.na(gene_names) | gene_names == ""] <-
      all_prots[is.na(gene_names) | gene_names == ""]
    pi_adj <- pmax(adj_p, .Machine$double.xmin); pi_adj[is.na(pi_adj)] <- 1

    result_df <- data.frame(
      GeneName        = gene_names,
      Accession       = all_prots,
      AvgIntensity    = round(avg_int, 3),
      Log2FC          = round(display_fc, 3),
      CtrlDetected    = paste0(round(ctrl_detect * n_ctrl), "/", n_ctrl),
      CondDetected    = paste0(round(cond_detect * n_cond), "/", n_cond),
      PiScore         = round(abs(display_fc) * (-log10(pi_adj)), 3),
      AdjPvalue       = signif(adj_p, 4),
      RawPvalue       = signif(raw_p, 4),
      TestMethod      = test_source,
      LowEvidence     = low_evidence_vec,
      # proDA test_diff() per-protein totals over ALL samples in the fit (not
      # per group). Display/export only; nothing filters on them.
      NObs            = as.integer(proda_res$n_obs),
      NApprox         = round(as.numeric(proda_res$n_approx), 3),
      stringsAsFactors = FALSE
    )
    if (exists("ctrl_rep_df") && nrow(ctrl_rep_df) > 0)
      result_df <- merge(result_df, ctrl_rep_df, by = "Accession", all.x = TRUE)
    if (exists("cond_rep_df") && nrow(cond_rep_df) > 0)
      result_df <- merge(result_df, cond_rep_df, by = "Accession", all.x = TRUE)
    col_order <- c("GeneName", "Accession", "AvgIntensity", "Log2FC",
                   sort(grep("^Ctrl_Rep|^Cond_Rep", colnames(result_df), value = TRUE)),
                   "CtrlDetected", "CondDetected", "PiScore", "AdjPvalue",
                   "RawPvalue", "TestMethod", "LowEvidence",
                   "NObs", "NApprox")
    col_order <- intersect(col_order, colnames(result_df))
    result_df <- result_df[, col_order, drop = FALSE]
    # Most significant first, larger effect breaking ties
    result_df <- result_df[order(result_df$AdjPvalue, -abs(result_df$Log2FC)), ]
    rownames(result_df) <- NULL

    plot_data <- data.frame(
      Accession    = result_df$Accession,
      GeneName     = result_df$GeneName,
      FC_linear    = 2 ^ result_df$Log2FC,
      Log2FC       = result_df$Log2FC,
      AdjPValue    = result_df$AdjPvalue,
      NegLog10P    = -log10(pmax(result_df$AdjPvalue, .Machine$double.xmin)),
      IsCapped     = FALSE,
      # Nothing is imputed on this route. The volcano's IsImputed slot carries the
      # zero-detection marker: one group has no detections, so the estimate
      # leans on proDA's dropout model rather than on measured values.
      IsImputed    = grepl("^0/", result_df$CtrlDetected) |
                     grepl("^0/", result_df$CondDetected),
      AvgIntensity = result_df$AvgIntensity,
      stringsAsFactors = FALSE
    )

    list(result_df        = result_df,
         plot_data        = plot_data,
         n_ctrl           = n_ctrl,
         n_cond           = n_cond,
         batch_in_model   = batch_in_model,
         subject_in_model = subject_in_model,
         n_zero_obs_ctrl    = sum(round(ctrl_detect * n_ctrl) == 0L),
         n_zero_obs_cond    = sum(round(cond_detect * n_cond) == 0L),
         n_test_proda_wald  = sum(test_source == "proda_wald"),
         n_test_welch_t     = sum(test_source == "welch_t"),
         n_test_fc_only     = sum(test_source == "fc_only"),
         n_test_untested    = sum(test_source == "untested"),
         n_low_evidence     = sum(low_evidence_vec))
  }, error = function(e) {
    logger_fn(paste0("[ProDA] Result processing error: ", e$message))
    list(error = e$message)
  })

  if (!is.null(build_result$error)) {
    return(list(success = FALSE,
                error   = paste0("ProDA post-processing failed: ", build_result$error),
                messages = msgs))
  }

  # ── Provenance ─────────────────────────────────────────────────────────────
  proda_v <- tryCatch(as.character(packageVersion("proDA")),    error = function(e) NA_character_)
  qf_v    <- tryCatch(as.character(packageVersion("QFeatures")), error = function(e) NA_character_)

  # Aggregation actually performed on this route. build_proda_protein_matrix
  # records the ingestion route and the summarize step carries it on the cache
  # as aggregation_method; the protein-level routes aggregate nothing at all,
  # so "robustSummary" must not be claimed for them.
  .agg_from_route <- function(route) {
    if (is.null(route) || length(route) != 1L || is.na(route) || !nzchar(route))
      return(NA_character_)
    switch(as.character(route),
           "pd_protein_direct"   = "none (protein-level input)",
           "maxquant_lfq_direct" = "none (protein-level input)",
           "robustSummary"       = "robustSummary",
           NA_character_)
  }

  # Reference-panel criteria that actually applied. The variability field
  # is named by whatever select_invariant_protein_panel used on this run, so read
  # the name off the thresholds list rather than assuming one.
  .ebn_thr      <- cache$ebn_thresholds
  .ebn_var_name <- if (is.list(.ebn_thr))
    setdiff(names(.ebn_thr), c("amp", "abund")) else character(0)
  .ebn_var_val  <- if (length(.ebn_var_name) == 1L)
    as.numeric(.ebn_thr[[.ebn_var_name[1]]]) else NA_real_

  provenance <- list(
    pipeline             = "ip",
    ip_sub               = ip_sub,
    skip_impute          = FALSE,
    pipeline_label       = cache$pipeline_label %||% "IP (ProDA)",
    lod_n_imputed        = 0L,
    n_prots_input        = cache$n_prots_input     %||% NA_integer_,
    n_prots_processed    = cache$n_prots_processed %||% NA_integer_,
    n_prots_output       = nrow(build_result$result_df),
    from_cache           = TRUE,
    n_conditions_total   = length(cond_groups),
    igg_conditions       = igg_conditions,
    ctrl_group           = ctrl_label,
    cond_group           = cond_label,
    analysis_mode        = "proda",
    de_engine            = "proda",
    normalization_method = cache$normalization_method %||% "none",
    # Search software of the label-free input, and FragPipe's workflow settings
    norm_software         = cache$norm_software         %||% NA_character_,
    fragpipe_settings     = cache$fragpipe_settings,
    software_used        = "proDA",
    pvalue_adjustment    = "BH",
    ebn_ref_proteins     = cache$ebn_ref_vector     %||% character(0),
    ebn_step             = cache$ebn_step           %||% NA_integer_,
    ebn_n_selected       = cache$ebn_n_selected     %||% NA_integer_,
    ebn_n_steps          = cache$ebn_n_steps        %||% NA_integer_,
    ebn_candidates_singleton = as.integer(
      cache$ebn_candidates_singleton %||% NA_integer_),
    ebn_candidates_zero_var  = as.integer(
      cache$ebn_candidates_zero_var  %||% NA_integer_),
    ebn_strict           = cache$ebn_strict          %||% NA,
    ebn_strict_max_step  = as.integer(
      cache$ebn_strict_max_step %||% NA_integer_),
    ebn_n_selected_prematch = as.integer(
      cache$ebn_n_selected_prematch %||% NA_integer_),
    ebn_n_spike_in_rows  = as.integer(
      cache$ebn_n_spike_in_rows %||% NA_integer_),
    ebn_threshold_amplitude = if (is.list(.ebn_thr))
      as.numeric(.ebn_thr$amp %||% NA_real_) else NA_real_,
    ebn_threshold_abundance = if (is.list(.ebn_thr))
      as.numeric(.ebn_thr$abund %||% NA_real_) else NA_real_,
    variance_prior_scope            = .variance_prior_scope,
    n_conditions_dropped_from_model = .n_conds_dropped_model,
    # Log2FC is proDA::test_diff()$diff on every row (Pass F); rows with no
    # observations in one arm of the contrast are counted, not floored.
    fc_source            = "proDA::test_diff()$diff",
    n_zero_obs_ctrl      = as.integer(build_result$n_zero_obs_ctrl %||% NA_integer_),
    n_zero_obs_cond      = as.integer(build_result$n_zero_obs_cond %||% NA_integer_),
    n_test_proda_wald    = as.integer(build_result$n_test_proda_wald  %||% NA_integer_),
    n_test_welch_t       = as.integer(build_result$n_test_welch_t     %||% NA_integer_),
    n_test_fc_only       = as.integer(build_result$n_test_fc_only     %||% NA_integer_),
    n_test_untested      = as.integer(build_result$n_test_untested    %||% NA_integer_),
    n_contrast_empty_removed = .n_contrast_empty_removed,
    n_low_evidence       = as.integer(build_result$n_low_evidence     %||% NA_integer_),
    aggregation_method   = cache$aggregation_method %||% NA_character_,
    proda_version        = proda_v,
    qfeatures_version    = qf_v,
    mq_stats             = cache$mq_stats,
    shared_peptide_stats = cache$shared_peptide_stats,
    proda_fit_strategy         = .proda_fit_strategy,
    proda_seed                 = .proda_seed,
    proda_max_iter_used        = .proda_max_iter_used,
    proda_df_corrected         = .proda_df_corrected,
    proda_median_pval_native   = .proda_med_p_native,
    proda_median_pval_moderated_diagnostic = .proda_med_p_moderated,
    proda_converged            = if (!is.null(conv)) isTRUE(conv$successful) else NA,
    proda_iterations           = if (!is.null(conv)) as.integer(conv$iterations %||% NA) else NA_integer_,
    proda_hessian_singular_events = .proda_hessian_events,
    normalization              = cache$normalization_method %||% "none",
    experiment_type            = if (identical(ip_sub, "invariant")) "ip_invariant" else "ip",
    mb_impute                  = FALSE,
    lod_imputation             = FALSE,
    starting_proteins          = as.numeric(cache$n_prots_input %||% NA),
    dropped_proteins           = as.numeric(
      (cache$n_prots_input %||% nrow(build_result$result_df)) -
        nrow(build_result$result_df)),
    drop_reason          = character(),
    summary_method       = .agg_from_route(cache$aggregation_method),
    censored_int         = "NA (genuine missing, not imputed)",
    bait_gene            = character(),
    bait_accs            = character(),
    norm_strategy        = character(),
    imputation_method    = "none",
    batch_corrected      = isTRUE(build_result$batch_in_model),
    # Batch label -> ORIGINAL sample names (design_df$sample), for the methods
    # text only; the model used the internal "<Group>_Rep<i>" names above.
    batch_groups         = if (isTRUE(build_result$batch_in_model)) {
      .orig_of <- if ("sample" %in% colnames(design_df))
                    setNames(as.character(design_df$sample), design_df$sample_name) else NULL
      lapply(batch_groups, function(s) {
        o <- if (is.null(.orig_of)) s else unname(.orig_of[s])
        ifelse(is.na(o) | !nzchar(o), s, o)
      })
    } else NULL,
    factorial_enabled    = isTRUE(is_factorial),
    interaction_model    = isTRUE(use_interaction),
    factorial_f1         = if (isTRUE(is_factorial)) fd$f1_name %||% "" else character(),
    factorial_f2         = if (isTRUE(is_factorial)) fd$f2_name %||% "" else character(),
    expect_mnar          = isTRUE(expect_mnar),
    ebn_fallback         = isTRUE(cache$ebn_fallback),
    subject_pairing_used = isTRUE(build_result$subject_in_model),
    subject_n            = if (isTRUE(build_result$subject_in_model) &&
                                "subject" %in% colnames(design_sub))
                             nlevels(design_sub$subject) else NA_integer_
  )

  provenance[[paste0("ebn_threshold_",
                     if (length(.ebn_var_name) == 1L) .ebn_var_name[1] else "cv")]] <-
    .ebn_var_val

  # Invariant-protein median normalization: how the invariance gate was centred.
  if (identical(ip_sub, "invariant")) {
    for (.k in c("invariant_centre_method", "invariant_centre", "invariant_mode_centre",
                 "invariant_median_centre", "invariant_peak_ratio", "invariant_centre_gap",
                 "invariant_centre_disagree", "invariant_mode_ambiguous"))
      provenance[[.k]] <- cache[[.k]]
  }

  .rdf <- build_result$result_df
  logger_fn(paste0("DE [ProDA]: ", nrow(.rdf), " proteins. adj.p<0.05: ",
                   sum(.rdf$AdjPvalue < 0.05, na.rm = TRUE),
                   " LowEvidence=", sum(.rdf$LowEvidence, na.rm = TRUE)))

  list(
    success     = TRUE,
    error       = NULL,
    result_df   = build_result$result_df,
    plot_data   = build_result$plot_data,
    acc_to_gene = acc_to_gene,
    provenance  = provenance,
    messages    = msgs
  )
}

# compute_lf_summarization -----------------------------------------------------
# Runs LF protein summarization for the Wizard's "Summarize" step.
# All heavy work (lf_core_summarization_proda/saint, maxquant/psm/lf_run_summarization)
# is pure computation  -" no reactive values.
# Returns list(success, error, result, genes_from_data,
#              rep_quality_flag, messages).
compute_lf_summarization <- function(
    mode,                    # "proda" | "saint" | "maxquant"
    # ProDA / SAINT / MaxQuant: pre-built matrix (from build_proda_protein_matrix in main thread)
    protein_matrix     = NULL,
    design_df          = NULL,
    acc_to_gene        = NULL,
    aggregation_method = NULL,
    genes_from_data    = FALSE,
    # Common IP/IgG
    ip_sub             = NULL,
    igg_conditions     = character(0),
    skip_impute        = FALSE,
    apply_normalization = TRUE,   # FALSE: skip OSMSS global median (input normalized upstream)
    logger_fn          = message
) {
  msgs <- character(0)

  result <- tryCatch({
    switch(mode,
      proda = {
        logger_fn("[Wizard] Summarization engine: ProDA (lf_core_summarization_proda)")
        res <- lf_core_summarization_proda(
          protein_matrix = protein_matrix,
          design_df      = design_df,
          acc_to_gene    = acc_to_gene,
          ip_sub         = ip_sub,
          igg_conditions = igg_conditions,
          apply_normalization = apply_normalization,
          logger         = logger_fn
        )
        res$genes_from_data <- isTRUE(genes_from_data)
        res
      },
      saint = {
        logger_fn("[Wizard] Summarization engine: SAINTexpress (lf_core_summarization_saint)")
        res <- lf_core_summarization_saint(
          protein_matrix     = protein_matrix,
          design_df          = design_df,
          acc_to_gene        = acc_to_gene,
          aggregation_method = aggregation_method %||% "pd_protein_direct",
          logger             = logger_fn
        )
        res$genes_from_data <- isTRUE(genes_from_data)
        res
      },
      maxquant = {
        logger_fn("[Wizard] Summarization engine: ProDA (MaxQuant Route MQ)")
        res <- lf_core_summarization_proda(
          protein_matrix = protein_matrix,
          design_df      = design_df,
          acc_to_gene    = acc_to_gene,
          ip_sub         = ip_sub,
          igg_conditions = igg_conditions,
          apply_normalization = apply_normalization,
          logger         = logger_fn
        )
        res$genes_from_data <- isTRUE(genes_from_data)
        res
      },
      stop(paste0("Unknown summarization mode: ", mode))
    )
  }, error = function(e) list(error = e$message))

  if (!is.null(result$error)) {
    return(list(success = FALSE, error = result$error, messages = msgs))
  }

  # ── Replicate quality flag ─────────────────────────────────────────────────
  rep_quality_flag <- local({
    pm        <- result$protein_matrix
    design_df <- result$design_df
    if (is.null(pm) || !is.matrix(pm) || ncol(pm) < 2) return(NULL)
    per_rep_na <- colMeans(is.na(pm))
    if (!is.null(design_df) &&
        "condition"   %in% colnames(design_df) &&
        "sample_name" %in% colnames(design_df)) {
      n_per_grp   <- as.integer(round(median(table(design_df$condition))))
      by_group    <- split(design_df$sample_name, design_df$condition)
      best_excess <- 0
      worst_samp  <- colnames(pm)[which.max(per_rep_na)]
      worst_rate  <- max(per_rep_na)
      for (grp_samps in by_group) {
        grp_samps <- intersect(grp_samps, colnames(pm))
        if (length(grp_samps) < 2) next
        grp_rates <- per_rep_na[grp_samps]
        med_rate  <- median(grp_rates)
        max_rate  <- max(grp_rates)
        excess    <- max_rate - med_rate
        if (excess > best_excess && max_rate > 0.30) {
          best_excess <- excess
          worst_samp  <- names(grp_rates)[which.max(grp_rates)]
          worst_rate  <- max_rate
        }
      }
      eff_n <- if (best_excess >= 0.20) n_per_grp - 1L else n_per_grp
      severity <-
        if      (eff_n <= 1L && best_excess >= 0.20)                                          "block"
        else if (eff_n <= 2L && n_per_grp <= 3L && best_excess >= 0.25)                       "block"
        else if (eff_n <= 2L && n_per_grp <= 3L && worst_rate > 0.30 && best_excess > 0.05)  "block"
        else if (eff_n <= 2L && n_per_grp <= 3L && best_excess >= 0.15)                       "warn"
        else                                                                                    "ok"
      list(severity          = severity,
           worst_rate        = worst_rate,
           worst_samp        = worst_samp,
           eff_n             = eff_n,
           n_per_grp         = n_per_grp,
           na_outlier_excess = round(best_excess, 3))
    } else {
      worst_rate <- max(per_rep_na)
      worst_samp <- colnames(pm)[which.max(per_rep_na)]
      n_per_grp  <- max(1L, as.integer(floor(ncol(pm) / 2L)))
      eff_n      <- floor(n_per_grp * (1 - worst_rate))
      severity   <- if (eff_n <= 1L || (worst_rate > 0.30 && n_per_grp <= 3L)) "block"
                    else if (eff_n <= 2L && n_per_grp <= 3L)                    "warn"
                    else                                                          "ok"
      list(severity          = severity,
           worst_rate        = worst_rate,
           worst_samp        = worst_samp,
           eff_n             = eff_n,
           n_per_grp         = n_per_grp,
           na_outlier_excess = worst_rate)
    }
  })

  list(
    success          = TRUE,
    error            = NULL,
    result           = result,
    genes_from_data  = isTRUE(result$genes_from_data),
    rep_quality_flag = rep_quality_flag,
    messages         = msgs
  )
}

# compute_saint_de -------------------------------------------------------------
# Runs the SAINTexpress DE pipeline in a background worker.
# All file I/O uses temp dirs (no reactive state needed).
# Returns list(success, error, result_df, plot_data, acc_to_gene,
#              provenance, messages).
compute_saint_de <- function(
    cache,
    ctrl_label,
    cond_label,
    ctrl_cols,
    cond_cols,
    bait_gene,
    cond_groups,
    saint_binary_path,
    logs_dir       = NULL,
    logger_fn      = message
) {
  msgs <- character(0)

  protein_matrix <- cache$protein_matrix
  acc_to_gene    <- cache$acc_to_gene

  # ── Generate SAINTexpress input files in temp dir ─────────────────────────
  work_dir <- file.path(tempdir(),
                        paste0("saint_", format(Sys.time(), "%Y%m%d%H%M%S"),
                               "_", as.integer(Sys.getpid())))
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

  # Inner helper: writes SAINTexpress input files (self-contained for worker isolation)
  .gen_saint_files <- function(pm, a2g, ctrl_c, cond_c, bait, wdir) {
    preys      <- rownames(pm)
    gene_names <- a2g[preys]
    gene_names[is.na(gene_names) | gene_names == ""] <-
      preys[is.na(gene_names) | gene_names == ""]
    bait_prey_idx <- which(gene_names == bait)
    if (length(bait_prey_idx) > 0) {
      pm         <- pm[-bait_prey_idx, , drop = FALSE]
      preys      <- preys[-bait_prey_idx]
      gene_names <- gene_names[-bait_prey_idx]
    }
    bait_rows <- data.frame(
      IP   = c(cond_c, ctrl_c),
      Bait = bait,
      Type = c(rep("T", length(cond_c)), rep("C", length(ctrl_c))),
      stringsAsFactors = FALSE
    )
    write.table(bait_rows, file.path(wdir, "bait.dat"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    prey_rows <- data.frame(
      Accession = preys, Length = 1000L, GeneName = gene_names,
      stringsAsFactors = FALSE
    )
    write.table(prey_rows, file.path(wdir, "prey.dat"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    all_cols     <- c(cond_c, ctrl_c)
    mat_sub      <- pm[, all_cols, drop = FALSE]
    was_detected <- is.finite(mat_sub)          # TRUE where protein was observed
    mat_lin      <- 2^mat_sub
    mat_lin[!was_detected] <- NA                # keep true absences as NA
    global_max <- max(mat_lin, na.rm = TRUE)
    if (is.finite(global_max) && global_max > 1e8)
      mat_lin <- mat_lin / global_max * 1e8
    mat_lin[was_detected & !is.na(mat_lin) & mat_lin < 1] <- 1   # floor only detected values
    int_rows <- do.call(rbind, lapply(all_cols, function(col) {
      vals  <- mat_lin[, col]
      keep  <- !is.na(vals)                     # omit undetected proteins (NA)
      data.frame(IP = col, Bait = bait, Prey = preys[keep],
                 Intensity = vals[keep], stringsAsFactors = FALSE)
    }))
    write.table(int_rows, file.path(wdir, "interaction.dat"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    invisible(pm)    # return trimmed matrix (bait row removed)
  }

  trimmed_matrix <- tryCatch(
    .gen_saint_files(protein_matrix, acc_to_gene, ctrl_cols, cond_cols,
                     bait_gene, work_dir),
    error = function(e) {
      return(list(error = paste0("SAINT input file generation failed: ", e$message)))
    }
  )
  if (is.list(trimmed_matrix) && !is.null(trimmed_matrix$error)) {
    return(list(success = FALSE, error = trimmed_matrix$error, messages = msgs))
  }

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  # Copy input files to logs_dir if provided
  if (!is.null(logs_dir) && dir.exists(logs_dir)) {
    for (ff in c("bait.dat", "prey.dat", "interaction.dat")) {
      src <- file.path(work_dir, ff)
      if (file.exists(src))
        file.copy(src, file.path(logs_dir,
                  paste0("saint_", sub("\\.", "_", ff), "_", ts, ".txt")),
                  overwrite = TRUE)
    }
  }

  # ── Execute SAINTexpress binary ───────────────────────────────────────────
  binary <- tryCatch(normalizePath(saint_binary_path, mustWork = TRUE),
                     error = function(e) {
                       return(list(error = paste0("SAINT binary not found: ", e$message)))
                     })
  if (is.list(binary) && !is.null(binary$error)) {
    return(list(success = FALSE, error = binary$error, messages = msgs))
  }

  logger_fn(paste0("[SAINT] Executing: ", binary))
  old_wd <- getwd()
  setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)

  result_lines <- system2(
    command = binary,
    args    = c("-L", as.character(length(ctrl_cols)),
                "interaction.dat", "prey.dat", "bait.dat"),
    stdout  = TRUE,
    stderr  = TRUE,
    wait    = TRUE
  )
  exit_code <- attr(result_lines, "status") %||% 0L

  # Save stderr to logs_dir
  if (!is.null(logs_dir) && dir.exists(logs_dir)) {
    log_path <- file.path(logs_dir, paste0("saint_stderr_", ts, ".txt"))
    writeLines(result_lines, log_path)
  }

  list_file <- file.path(work_dir, "list.txt")
  if (exit_code != 0 || !file.exists(list_file)) {
    err_msg <- paste(result_lines, collapse = "\n")
    return(list(success = FALSE,
                error   = paste0("SAINTexpress failed (exit ", exit_code, "): ",
                                 substr(err_msg, 1, 500)),
                messages = msgs))
  }

  # Copy SAINTexpress output (list.txt) to logs_dir if provided
  if (!is.null(logs_dir) && dir.exists(logs_dir)) {
    file.copy(list_file,
              file.path(logs_dir, paste0("saint_list_txt_", ts, ".txt")),
              overwrite = TRUE)
  }

  saint_out <- tryCatch(
    read.delim(list_file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(saint_out) || nrow(saint_out) == 0) {
    return(list(success = FALSE,
                error   = "SAINTexpress produced an empty output file (list.txt).",
                messages = msgs))
  }
  saint_fc <- if ("FoldChange" %in% colnames(saint_out))
    suppressWarnings(as.numeric(saint_out$FoldChange)) else rep(NA_real_, nrow(saint_out))

  logger_fn(paste0("[SAINT] Done: ", nrow(saint_out), " interactions scored."))

  # ── Build result table ───────────────────
  n_cond    <- length(cond_cols)
  n_ctrl    <- length(ctrl_cols)
  bfdr_raw  <- as.numeric(saint_out$BFDR); bfdr_raw[is.na(bfdr_raw)] <- 1
  nonzero   <- bfdr_raw[bfdr_raw > 0]
  pseudo_min <- if (length(nonzero) > 0) min(nonzero) / 10 else 1e-10
  bfdr_safe <- ifelse(bfdr_raw == 0, pseudo_min, bfdr_raw)

  pm_idx       <- match(saint_out$Prey, rownames(protein_matrix))
  is_unmatched <- is.na(pm_idx)
  if (any(is_unmatched)) {
    return(list(success = FALSE,
                error   = paste0("SAINT returned ", sum(is_unmatched),
                                 " prey ID(s) absent from the protein matrix; ",
                                 "accession strings may have been altered by the binary."),
                messages = msgs))
  }
  ctrl_mat_raw <- protein_matrix[pm_idx, ctrl_cols, drop = FALSE]
  cond_mat_raw <- protein_matrix[pm_idx, cond_cols, drop = FALSE]
  ctrl_detect  <- rowSums(!is.na(ctrl_mat_raw))
  cond_detect  <- rowSums(!is.na(cond_mat_raw))
  single_replicate <- cond_detect == 1L
  is_exclusive <- (ctrl_detect == 0) & (cond_detect > 0)
  is_ctrl_only <- (cond_detect == 0) & (ctrl_detect > 0)
  fc_origin    <- rep("measured", length(is_exclusive))
  fc_origin[is_exclusive] <- "control-absent"
  fc_origin[is_ctrl_only] <- "bait-absent"
  fc_origin[is_unmatched] <- "unmatched"

  mean_cond <- rowMeans(cond_mat_raw, na.rm = TRUE)
  mean_ctrl <- rowMeans(ctrl_mat_raw, na.rm = TRUE)
  log2fc    <- mean_cond - mean_ctrl
  excl_floor   <- if (sum(!is.na(ctrl_mat_raw)) > 10)
    quantile(ctrl_mat_raw, 0.05, na.rm = TRUE)
  else
    min(ctrl_mat_raw, na.rm = TRUE) - 1
  excl_floor_src <- if (sum(!is.na(ctrl_mat_raw)) > 10) "ctrl_q05" else "ctrl_min_minus1"
  excl_floor_src <- if (is.finite(excl_floor)) excl_floor_src else "cond_min_minus3"
  excl_floor <- if (is.finite(excl_floor)) excl_floor else min(mean_cond, na.rm = TRUE) - 3
  log2fc[is_exclusive] <- mean_cond[is_exclusive] - excl_floor
  finite_fc <- log2fc[is.finite(log2fc)]
  fc_cap    <- if (length(finite_fc) > 0) max(abs(finite_fc), na.rm = TRUE) + 1 else 15
  log2fc[log2fc ==  Inf] <-  fc_cap
  log2fc[log2fc == -Inf] <- -fc_cap
  log2fc[is_ctrl_only]   <- NA_real_
  log2fc[is.nan(log2fc)] <- NA_real_

  avg_spec <- as.numeric(saint_out$AvgIntensity)
  avg_spec[is.na(avg_spec) | avg_spec <= 0] <- .Machine$double.xmin
  avg_int  <- log2(avg_spec)

  n_res <- nrow(saint_out)

  # Confidence tiers come from SAINTexpress's own SaintScore
  sp          <- as.numeric(saint_out$SaintScore); sp[is.na(sp)] <- 0
  confidence  <- rep("Low", n_res)
  confidence[sp >= 0.5]                            <- "Medium"
  confidence[sp >= 0.8]                            <- "High"
  confidence[is_exclusive & sp >= 0.8]             <- "Exclusive"
  confidence[is_exclusive & sp >= 0.5 & sp < 0.8] <- "High"
  confidence[is_exclusive & sp <  0.5]             <- "Medium"

  accessions  <- saint_out$Prey
  gene_names  <- saint_out$PreyGene
  na_gene     <- is.na(gene_names) | gene_names == ""
  gene_names[na_gene] <- accessions[na_gene]
  pi_adj   <- pmax(bfdr_safe, 1e-300)
  pi_score <- round(abs(log2fc) * (-log10(pi_adj)), 3)

  ctrl_rep_df <- data.frame(Accession = accessions, stringsAsFactors = FALSE)
  cond_rep_df <- data.frame(Accession = accessions, stringsAsFactors = FALSE)
  ctrl_mat_sub <- protein_matrix[match(accessions, rownames(protein_matrix)),
                                  ctrl_cols, drop = FALSE]
  cond_mat_sub <- protein_matrix[match(accessions, rownames(protein_matrix)),
                                  cond_cols, drop = FALSE]
  # Use actual sample column names so the output is self-documenting
  for (si in seq_len(n_ctrl))
    ctrl_rep_df[[ctrl_cols[si]]] <- round(ctrl_mat_sub[, si], 2)
  for (si in seq_len(n_cond))
    cond_rep_df[[cond_cols[si]]] <- round(cond_mat_sub[, si], 2)

  result_df <- data.frame(
    GeneName         = gene_names,
    Accession        = accessions,
    Log2AvgIntensity = round(avg_int, 3),   # log2 of SAINTexpress AvgSpec
    Log2FC           = round(log2fc, 3),
    FC_Origin        = fc_origin,
    CtrlDetected     = paste0(ctrl_detect, "/", n_ctrl),
    CondDetected     = paste0(cond_detect, "/", n_cond),
    SingleReplicate  = single_replicate,
    PiScore          = pi_score,
    AdjPvalue        = signif(bfdr_safe, 4),  # floor-corrected BFDR (internal use)
    Confidence       = confidence,            # tier label (internal use for colors)
    SaintScore       = round(sp, 4),
    SAINT_FoldChange = signif(saint_fc, 4),
    BFDR             = signif(bfdr_raw, 4),   # raw Bayesian FDR from SAINTexpress
    stringsAsFactors = FALSE
  )
  result_df <- merge(result_df, ctrl_rep_df, by = "Accession", all.x = TRUE)
  result_df <- merge(result_df, cond_rep_df, by = "Accession", all.x = TRUE)
  col_order <- c("GeneName", "Accession", "Log2AvgIntensity", "Log2FC", "FC_Origin",
                 ctrl_cols, cond_cols,
                 "CtrlDetected", "CondDetected", "SingleReplicate",
                 "PiScore", "AdjPvalue", "Confidence",
                 "SaintScore", "SAINT_FoldChange", "BFDR")
  col_order <- intersect(col_order, colnames(result_df))
  result_df <- result_df[, col_order, drop = FALSE]
  # Highest SaintScore first, lower BFDR breaking ties
  result_df <- result_df[order(-result_df$SaintScore, result_df$BFDR), ]
  rownames(result_df) <- NULL

  plot_data <- data.frame(
    Accession    = result_df$Accession,
    GeneName     = result_df$GeneName,
    FC_linear    = 2 ^ result_df$Log2FC,
    Log2FC       = result_df$Log2FC,
    AdjPValue    = result_df$AdjPvalue,
    NegLog10P    = -log10(pmax(result_df$AdjPvalue, .Machine$double.xmin)),
    IsCapped     = FALSE,
    IsImputed    = result_df$FC_Origin == "control-absent",
    FC_Origin    = result_df$FC_Origin,
    SingleReplicate = result_df$SingleReplicate,
    Confidence   = result_df$Confidence,
    AvgIntensity = result_df$Log2AvgIntensity,
    stringsAsFactors = FALSE
  )

  # ── Raw abundance export df ──────────────────────────────────────────────
  raw_wide <- as.data.frame(protein_matrix)
  raw_wide$Accession <- rownames(raw_wide)

  provenance <- list(
    pipeline             = "ip",
    ip_sub               = "saint",
    de_engine            = "saint",
    analysis_mode        = "saint",
    pipeline_label       = "IP (SAINTexpress)",
    software_used        = "SAINTexpress-int v3.6.3",
    pvalue_adjustment    = "BFDR (Bayesian FDR)",
    normalization_method = "none",
    ctrl_group           = ctrl_label,
    cond_group           = cond_label,
    bait_protein         = bait_gene,
    n_prots_input        = nrow(cache$protein_matrix),
    n_prots_processed    = nrow(protein_matrix),
    n_prots_output       = nrow(result_df),
    n_conditions_total   = 2L,
    expect_mnar          = TRUE,
    from_cache           = TRUE,
    saint_binary         = saint_binary_path,
    aggregation_method   = cache$aggregation_method %||% "pd_protein_direct",
    mq_stats             = cache$mq_stats,
    shared_peptide_stats = cache$shared_peptide_stats,
    norm_software        = cache$norm_software %||% NA_character_,
    fragpipe_settings    = cache$fragpipe_settings,
    normalization        = "none (OSMSS); SAINTexpress internal normalization applied during scoring",
    experiment_type      = "ip_saint",
    result_schema        = "saint_v2",
    n_pre_saint_input    = cache$pre_saint_input %||% NA_integer_,
    n_dropped_pre_saint  = cache$pre_saint_dropped %||% NA_integer_,
    skip_impute          = TRUE,
    lod_n_imputed        = sum(fc_origin == "control-absent"),
    lod_imputation       = any(fc_origin == "control-absent"),
    fc_floor_value       = unname(excl_floor),
    fc_floor_source      = excl_floor_src,
    n_single_replicate   = sum(single_replicate),
    n_bait_absent        = sum(is_ctrl_only),
    n_unmatched          = sum(is_unmatched),
    mb_impute            = FALSE,
    batch_corrected      = FALSE,
    batch_groups         = NULL,
    factorial_enabled    = FALSE,
    subject_pairing_used = FALSE,
    starting_proteins    = as.numeric(nrow(cache$protein_matrix)),
    dropped_proteins     = as.numeric(nrow(cache$protein_matrix) - nrow(protein_matrix))
  )

  list(
    success      = TRUE,
    error        = NULL,
    result_df    = result_df,
    plot_data    = plot_data,
    acc_to_gene  = acc_to_gene,
    raw_abund_df = raw_wide,
    ctrl_cols    = ctrl_cols,
    cond_cols    = cond_cols,
    bait_gene    = bait_gene,
    provenance   = provenance,
    messages     = msgs
  )
}

# compute_tmt_summarization ----------------------------------------------------
# Runs TMT protein summarization in a background worker.
# Supports two routes: "limma_protein", "limma_qf".
# Returns list(success, error, route, result, messages).
# -----------------------------------------------------------------------------
# tmt_normalization_check()
# Single-plex TMT counterpart of dia_normalization_check(): do the uploaded
# reporter abundances look normalized across channels already? Raw reporter
# abundances (PD "Abundance: F#:" columns, other PSM software) normally are not,
# and OSMSS median-centers them. They can still have been normalized before
# upload (another tool, a script); centering such data a second time can shift
# the unchanged background when a large share of the proteome differs between
# conditions. Runs on the protein rows of the upload, log2 and before centering
# (tmt_limma_protein_summarize(apply_centering = FALSE)), and compares channel
# medians between replicate channels of each condition. Thresholds are
# dia_normalization_check()'s (calibrated on Spectronaut DIA exports). Advice
# only: OSMSS centers unless the user reports the input as normalized.
# Returns dia_normalization_check()'s list.
# -----------------------------------------------------------------------------
tmt_normalization_check <- function(files, annotation_map, rv = list(),
                                    logger = function(...) invisible(NULL)) {
  out <- list(verdict = "inconclusive", within_sd = NA_real_, n_complete = 0L,
              n_needed = NA_integer_, n_runs_used = 0L, between_offsets = NULL, message = "")
  if (is.null(annotation_map) || nrow(annotation_map) == 0L || length(files) == 0L) {
    out$message <- "No data or no channel assignments yet."
    return(out)
  }
  bio <- annotation_map[!is.na(annotation_map$Condition) & nzchar(annotation_map$Condition) &
                          !annotation_map$Condition %in% "Norm", , drop = FALSE]
  pm <- tmt_limma_protein_summarize(files = files, annotation_map = annotation_map, rv = rv,
                                    design_opt = 1L, logger = logger, apply_centering = FALSE)
  if (!is.null(pm$error)) {
    out$message <- paste0("The check could not read protein abundances: ", pm$error)
    return(out)
  }
  dia_normalization_check(2^pm$protein_matrix, split(bio$ColumnName, bio$Condition),
                          feature_label = "proteins", unit = "channel",
                          normalized_label = "cross-channel normalized")
}

compute_tmt_summarization <- function(
    route,
    files,
    annotation_map,
    design_opt      = 1L,
    wizard_data_level = NULL,
    shared_peptides = "unique_only",
    rv_snapshot     = list(),   # plain list with pd_has_gn, pd_gene_col etc.
    logger_fn       = message,
    apply_centering = TRUE      # FALSE: single plex reported as already normalized
) {
  msgs <- character(0)

  result <- tryCatch({
    switch(route,
      limma_protein = {
        logger_fn("[TMT] Summarization route: limma_protein")
        tmt_limma_protein_summarize(
          files          = files,
          annotation_map = annotation_map,
          rv             = rv_snapshot,
          design_opt     = design_opt,
          logger         = logger_fn,
          apply_centering = apply_centering
        )
      },
      limma_qf = {
        logger_fn("[TMT] Summarization route: limma_qf")
        tmt_limma_qf_summarize(
          files           = files,
          annotation_map  = annotation_map,
          data_level      = wizard_data_level,
          rv              = rv_snapshot,
          design_opt      = design_opt,
          shared_peptides = shared_peptides,
          logger          = logger_fn,
          apply_centering = apply_centering
        )
      },
      {
        stop(paste0("Unknown TMT summarization route: '", route,
                    "'. Expected 'limma_protein' or 'limma_qf'."))
      }
    )
  }, error = function(e) list(error = e$message))

  if (!is.null(result$error)) {
    return(list(success = FALSE, route = route, error = result$error, messages = msgs))
  }

  list(success = TRUE, route = route, result = result, error = NULL, messages = msgs)
}
# =============================================================================
# TMT_PREPARE_MATRIX_FOR_DE
# Shared pre-model spine for every TMT design option. Filtering, imputation,
# zero-origin handling and count provenance used to be written once per design
# option, which is how the three routes drifted apart. They are written once
# here instead.
#
# The design option selects ONE thing: the imputation neighbourhood.
#   1 -> global k-NN over every annotated channel (single plex)
#   2 -> none (bridge-normalised multi-plex)
#   3 -> within-plex k-NN (run-covariate multi-plex)
# The coverage filter, the pre-imputation snapshot, the zero-origin hold-back
# and the counters are unconditional.
#
# Returns: list(matrix, observed_matrix, n_before, n_kept, n_filled_global,
#               n_filled_within_plex, n_zero_held, impute_pkg_ok,
#               n_channels_over_colmax, route_imputation)
# =============================================================================
tmt_prepare_matrix_for_de <- function(protein_matrix,
                                      design_df,
                                      design_opt,
                                      batch_groups         = NULL,
                                      zero_origin_mask     = NULL,
                                      min_valid_frac       = 0.5,
                                      max_missing_per_cond = 1L,
                                      plex_min_obs_frac    = 0.5,
                                      zero_is_absence      = TRUE,
                                      knn_k                = 10L,
                                      logger               = message) {

  n_before        <- nrow(protein_matrix)
  n_global_filled <- 0L
  n_within_plex   <- 0L
  n_zero_held     <- 0L
  impute_pkg_ok   <- NA
  n_channels_over_colmax <- 0L
  route_imputation <- "none"
  # Silent-degradation counters: both of these were logged but never reached
  # provenance, so a run could quietly lose zero-origin protection or a plex's
  # imputation with no record in the reported numbers.
  n_zero_mask_dropped <- 0L
  n_plex_knn_failed   <- 0L

  pm_g  <- protein_matrix
  gv_g  <- design_df$condition[match(colnames(pm_g), design_df$sample_name)]
  cds_g <- unique(gv_g[!is.na(gv_g) & nzchar(gv_g)])

  .frac_valid <- function(m, cc) {
    idx <- which(gv_g == cc)
    rowSums(!is.na(m[, idx, drop = FALSE])) / length(idx)
  }
  .n_missing <- function(m, cc) {
    idx <- which(gv_g == cc)
    rowSums(is.na(m[, idx, drop = FALSE]))
  }

  # 1. Global filter: quantifiable in at least one condition. Deliberately
  #    permissive - it defines the protein universe, not the tested set.
  #    Runs BEFORE any imputation, on every route: a filter applied to imputed
  #    data judges invented values and rejects almost nothing.
  vf_g    <- vapply(cds_g, function(cc) .frac_valid(pm_g, cc),
                    numeric(nrow(pm_g)))
  if (is.null(dim(vf_g))) vf_g <- matrix(vf_g, nrow = nrow(pm_g))
  keep_g  <- apply(vf_g, 1L, function(x) any(x >= min_valid_frac, na.rm = TRUE))
  keep_g[is.na(keep_g)] <- FALSE
  pm_g    <- pm_g[keep_g, , drop = FALSE]
  n_kept  <- nrow(pm_g)
  logger(paste0("[DE TMT] Global filter (>=", round(min_valid_frac * 100),
                "% valid in >=1 condition): ", length(keep_g), " -> ",
                n_kept, " proteins."))
  if (n_kept < 2L)
    stop("Fewer than 2 proteins pass the global missingness filter.")

  zmask_g <- zero_origin_mask
  if (!is.null(zmask_g)) {
    # Degrade rather than fail. A mask that does not cover the annotated
    # channel set (stale or partially migrated bundle) is dropped, which
    # reverts to the same behaviour as no mask at all.
    zmask_g <- tryCatch({
      if (!all(rownames(pm_g) %in% rownames(zmask_g)) ||
          !all(colnames(pm_g) %in% colnames(zmask_g)))
        stop("mask does not cover the annotated channel set")
      zmask_g[rownames(pm_g), colnames(pm_g), drop = FALSE]
    }, error = function(e) {
      logger(paste0("[DE TMT] zero-origin mask unusable (", e$message,
                    ") - proceeding without zero-origin exclusion."))
      NULL
    })
    if (is.null(zmask_g)) n_zero_mask_dropped <- 1L
  }
  # The pre-imputation snapshot is captured immediately after the filter and
  # before any imputation, on every route. Everything downstream that claims to
  # count observations reads this, not the model matrix.
  pm_pre_impute <- pm_g

  # ---- Global k-NN imputation (single plex only) ----------------------------
  # Neighbourhood = every annotated channel. Only meaningful when there is one
  # plex; with two or more, cross-plex neighbours carry the plex offset into
  # the fill. The multi-plex route uses the within-plex block below.
  if (identical(as.integer(design_opt), 1L)) {

    # Eligible only if lightly missing in EVERY condition, not just the two
    # being compared - which is what makes the imputed matrix independent of
    # the contrast the user selected.
    if (requireNamespace("impute", quietly = TRUE)) {
      impute_pkg_ok    <- TRUE
      route_imputation <- "global_knn"
      mn_g <- vapply(cds_g, function(cc) .n_missing(pm_g, cc),
                     numeric(nrow(pm_g)))
      if (is.null(dim(mn_g))) mn_g <- matrix(mn_g, nrow = nrow(pm_g))
      elig_g <- apply(mn_g, 1L,
                      function(x) all(x <= max_missing_per_cond, na.rm = TRUE)) &
                (rowSums(is.na(pm_g)) > 0L)
      elig_g[is.na(elig_g)] <- FALSE

      if (sum(elig_g) >= 2L && sum(is.na(pm_g[elig_g, , drop = FALSE])) > 0L) {
        # colmax is evaluated over the FULL filtered matrix rather than the
        # dense eligible submatrix, so a badly-covered channel that the
        # submatrix used to hide can now abort the whole call. Say which
        # channel before trying, otherwise the failure reads as a generic
        # "k-NN failed" with nothing imputed anywhere.
        .col_na_frac <- colMeans(is.na(pm_g))
        n_channels_over_colmax <- sum(.col_na_frac > 0.8, na.rm = TRUE)
        if (n_channels_over_colmax > 0L)
          logger(paste0(
            "[DE TMT] WARNING: ", n_channels_over_colmax,
            " channel(s) exceed 80% missing across the filtered matrix (",
            paste(colnames(pm_g)[.col_na_frac > 0.8], collapse = ", "),
            "). impute.knn's colmax = 0.8 will reject the entire matrix and ",
            "no values will be imputed on any protein."))
        imp_g <- tryCatch(
          # Donor pool is the full filtered matrix, not the eligible rows.
          # Restricting donors to near-complete proteins skews them toward the
          # high-abundance, better-measured end of the range and biases the
          # fills upward. Canonical k-NN searches the whole feature space and
          # writes back selectively; the write-back below is what stays
          # restricted to eligible rows.
          impute::impute.knn(pm_g, k = knn_k,
                             rowmax = 1, colmax = 0.8,
                             maxp = Inf, rng.seed = 42L)$data,
          error = function(e) {
            logger("[DE TMT] Global k-NN failed - NAs retained.")
            NULL
          }
        )
        if (!is.null(imp_g)) {
          # imp_g is full-size; everything below operates on the eligible slice.
          imp_e <- imp_g[elig_g, , drop = FALSE]
          if (isTRUE(zero_is_absence) && !is.null(zmask_g)) {
            zsub <- zmask_g[elig_g, , drop = FALSE]
            n_zero_held <- n_zero_held +
                           sum(zsub & is.na(pm_g[elig_g, , drop = FALSE]))
            imp_e[zsub] <- NA_real_
          }
          # Count fills within eligible rows only; values computed for
          # ineligible rows are donors, not results, and are discarded.
          n_global_filled <- sum(is.na(pm_g[elig_g, , drop = FALSE])) -
                             sum(is.na(imp_e))
          pm_g[elig_g, ] <- imp_e
        }
      }
      logger(paste0("[DE TMT] Global k-NN (k=", knn_k, "): ", n_global_filled,
                    " values filled across all conditions; ", n_zero_held,
                    " zero-origin cells left as NA; ", sum(is.na(pm_g)),
                    " NAs retained for lmFit."))
    } else {
      impute_pkg_ok <- FALSE
      .osmss_user_warn(paste0(
        "The 'impute' package is not available, so no missing values were ",
        "imputed. This analysis ran complete-case only: proteins with any ",
        "missing value contribute fewer observations to the model, and results ",
        "will differ from a run with imputation enabled. Install the 'impute' ",
        "Bioconductor package and re-run."))
      logger("[DE TMT] 'impute' unavailable - no imputation performed.")
    }
  }

  # ---- Within-plex k-NN (run-covariate multi-plex) --------------------------
  # Impute independently within each plex so that k-NN borrows from proteins
  # measured in the same run rather than across runs. Proteins too sparse
  # within a plex retain their NAs; lmFit fits them on available observations.
  if (identical(as.integer(design_opt), 3L) && !is.null(batch_groups) &&
      requireNamespace("impute", quietly = TRUE)) {
    impute_pkg_ok    <- TRUE
    route_imputation <- "within_plex_knn"
    pm_wk            <- pm_g
    n_plexes_imputed <- 0L
    # Same wording as the global block, so a dropped mask reads identically on
    # every route rather than being silently skipped here.
    if (isTRUE(zero_is_absence) && is.null(zmask_g))
      logger(paste0("[DE TMT] Opt=3: no usable zero-origin mask ",
                    "- proceeding without zero-origin exclusion."))
    for (plex_chans in batch_groups) {
      plex_cols <- intersect(plex_chans, colnames(pm_wk))
      if (length(plex_cols) < 2L) next
      sub_plex  <- pm_wk[, plex_cols, drop = FALSE]
      if (sum(is.na(sub_plex)) == 0L) next
      # Row eligibility, per plex. Without this a protein seen in 1 of 10
      # channels receives 9 fabricated values.
      elig_plex <- (rowSums(!is.na(sub_plex)) >=
                      ceiling(plex_min_obs_frac * length(plex_cols))) &
                   (rowSums(is.na(sub_plex)) > 0L)
      elig_plex[is.na(elig_plex)] <- FALSE
      if (!any(elig_plex)) next
      imp_plex <- tryCatch(
        # rowmax = 1: never swap k-NN for impute.knn's sample-mean fallback.
        # Rows too sparse to impute keep their NAs and are handled downstream.
        # Imputed over the FULL plex submatrix so the donor pool is every
        # protein in the plex, not just the eligible ones; the write-back
        # below is what stays restricted.
        impute::impute.knn(sub_plex, k = knn_k,
                           rowmax = 1, colmax = 0.8,
                           maxp = Inf, rng.seed = 42L)$data,
        error = function(e) {
          logger("[DE TMT] within-plex k-NN failed - NAs retained.")
          n_plex_knn_failed <<- n_plex_knn_failed + 1L
          sub_plex
        }
      )
      .imp_e <- imp_plex[elig_plex, , drop = FALSE]
      # Zero-origin hold-back, identical rule to the global block: a cell whose
      # only evidence was a non-positive reporter ion is not filled by a
      # nearest-neighbour (MAR) method.
      if (isTRUE(zero_is_absence) && !is.null(zmask_g)) {
        .z_hold <- zmask_g[elig_plex, plex_cols, drop = FALSE] &
                   is.na(sub_plex[elig_plex, , drop = FALSE])
        n_zero_held    <- n_zero_held + sum(.z_hold)
        .imp_e[.z_hold] <- NA_real_
      }
      .n_fill_plex <- sum(is.na(sub_plex[elig_plex, , drop = FALSE])) -
                      sum(is.na(.imp_e))
      pm_wk[elig_plex, plex_cols] <- .imp_e
      n_within_plex    <- n_within_plex + .n_fill_plex
      n_plexes_imputed <- n_plexes_imputed + 1L
      logger(paste0("[DE TMT] Opt=3 plex ", n_plexes_imputed, ": ",
                    sum(elig_plex), " of ", nrow(sub_plex),
                    " proteins eligible (>=",
                    ceiling(plex_min_obs_frac * length(plex_cols)), " of ",
                    length(plex_cols), " channels observed); ",
                    .n_fill_plex, " values filled."))
    }
    if (n_plexes_imputed > 0L)
      logger(paste0("[DE TMT] Opt=3: within-plex kNN imputation applied to ",
                    n_plexes_imputed, " plex(es)."))
    pm_g <- pm_wk
  }

  list(
    matrix                 = pm_g,
    observed_matrix        = pm_pre_impute,
    n_before               = n_before,
    n_kept                 = n_kept,
    n_filled_global        = n_global_filled,
    n_filled_within_plex   = n_within_plex,
    n_zero_held            = n_zero_held,
    impute_pkg_ok          = impute_pkg_ok,
    n_channels_over_colmax = n_channels_over_colmax,
    route_imputation       = route_imputation,
    n_zero_mask_dropped    = n_zero_mask_dropped,
    n_plex_knn_failed      = n_plex_knn_failed,
    zero_origin_mask       = zmask_g
  )
}


# compute_tmt_de ---------------------------------------------------------------
# Runs TMT differential expression in a background worker via limma.
# All design options handled: 1=single-plex, 2=multi-plex IRS bridge,
# 3=multi-plex run-covariate balanced.
# Returns list(success, error, route, result_df, plot_data, analysis_mode,
#              cond_test, cond_ref, provenance_tmt_volcano, provenance_de,
#              messages).
compute_tmt_de <- function(
    tmt_master_cache,
    design_opt,
    cond_test,
    cond_ref,
    tmt_analysis_route = NULL,
    acc_to_gene        = NULL,
    logger_fn          = message
) {
  msgs <- character(0)
  design_opt <- as.integer(design_opt %||% 1L)

  result <- tryCatch({
    # All design options (1 = single-plex, 2 = multi-plex IRS, 3 = run-covariate)
    # are handled by limma via dia_limma_de.
    logger_fn(paste0("[DE TMT] Design option ", design_opt, " (",
                     switch(as.character(design_opt),
                       "1" = "single-plex",
                       "2" = "multi-plex IRS bridge",
                       "3" = "multi-plex run-covariate",
                       paste0("opt", design_opt)), "): using limma"))

    batch_groups    <- NULL
    batch_for_model <- FALSE

    # limma_protein / limma_qf cache: protein_matrix already normalised
    if (!is.null(tmt_master_cache$protein_matrix)) {
      pm  <- tmt_master_cache$protein_matrix
      ann <- tmt_master_cache$annotation_map
      # Pass ALL real sample channels, not just the two compared conditions.
      # dia_limma_de(fit_all_groups = TRUE) fits one model over all of them and
      # extracts the requested contrast, which substantially increases residual
      # df. Bridge ("Norm") channels are excluded from the model.
      use_chans <- ann$ColumnName[!is.na(ann$Condition) &
                                    !ann$Condition %in% "Norm" &
                                    nzchar(ann$Condition)]
      use_chans <- intersect(use_chans, colnames(pm))
      if (!all(c(cond_test, cond_ref) %in% ann$Condition[ann$ColumnName %in% use_chans]))
        stop("Selected conditions are not present among the annotated channels.")
      if (length(use_chans) < 3)
        stop("Not enough channels for the selected conditions.")
      design_df <- data.frame(
        sample_name = use_chans,
        condition   = ann$Condition[match(use_chans, ann$ColumnName)],
        stringsAsFactors = FALSE
      )
      tmt_cache <- list(
      protein_matrix = pm[, use_chans, drop = FALSE],
      design_df      = design_df,
      acc_to_gene    = acc_to_gene,
      # Peptides aggregated per protein. Reporting only - the
      # one-condition-missing table flags single-peptide proteins, which appear
      # there disproportionately for mechanical reasons. Never used in the model.
      pep_count      = tmt_master_cache$pep_count_per_protein
    )
      # For design_opt=3: build run-covariate batch_groups from Mixture column
      if (identical(design_opt, 3L)) {
        mix_ann  <- ann[ann$ColumnName %in% use_chans, , drop = FALSE]
        mix_list <- split(mix_ann$ColumnName, mix_ann$Mixture)
        if (length(mix_list) >= 2) {
          batch_groups    <- mix_list
          batch_for_model <- TRUE
          logger_fn(paste0("[DE TMT] Run-covariate: ", length(mix_list),
                           " plexes added as batch covariate."))
        } else {
          logger_fn("[DE TMT] Run-covariate: only 1 plex in selected channels, ignoring batch covariate.")
        }
      }
    }

    # ── Track protein count before missingness filter ─────────────────────────
    n_proteins_before <- nrow(tmt_cache$protein_matrix)

    # ---- Shared pre-model spine (all design options) ------------------------
    # Filter, imputation, zero-origin hold-back and count provenance live in
    # tmt_prepare_matrix_for_de(). The design option selects only the imputation
    # neighbourhood; everything else is route-independent.
    #
    # Pipeline parameters, held here so they reach both the spine and prov_tv.
    # Minimum observed fraction within a plex for a protein to be eligible for
    # within-plex k-NN. 0.5 matches min_valid_frac used by the global filter and
    # by dia_limma_de, so one coverage rule governs the whole pipeline. Passing
    # the testing filter earns a model fit, not a fabricated value.
    plex_min_obs_frac <- 0.5
    # Do not invent a value where the instrument reported zero. Set to FALSE to
    # treat non-positive intensities as ordinary missingness.
    .zero_is_absence  <- TRUE
    .prep <- tmt_prepare_matrix_for_de(
      protein_matrix       = tmt_cache$protein_matrix,
      design_df            = design_df,
      design_opt           = design_opt,
      batch_groups         = if (isTRUE(batch_for_model)) batch_groups else NULL,
      zero_origin_mask     = tmt_master_cache$zero_origin_mask,
      min_valid_frac       = 0.5,
      max_missing_per_cond = 1L,
      plex_min_obs_frac    = plex_min_obs_frac,
      zero_is_absence      = .zero_is_absence,
      knn_k                = 10L,
      logger               = logger_fn
    )
    tmt_cache$protein_matrix <- .prep$matrix
    pm_pre_impute            <- .prep$observed_matrix
    zmask_g                  <- .prep$zero_origin_mask
    n_global_kept            <- .prep$n_kept
    n_global_filled          <- .prep$n_filled_global
    n_values_imputed_within_plex <- .prep$n_filled_within_plex
    n_zero_held              <- .prep$n_zero_held
    impute_pkg_ok            <- .prep$impute_pkg_ok
    n_channels_over_colmax   <- .prep$n_channels_over_colmax
    route_imputation         <- .prep$route_imputation
    # The spine always filters, on every route.
    global_filter_applied    <- TRUE

    # ---- Proteins not tested (reporting only) --------------------------------
    # Every protein a coverage filter removed, listed with the coverage it had so
    # the user can judge it: the missingness filter just applied (all routes) and,
    # on the protein route of a multi-plex design, the plex-coverage filter applied
    # at summarization (tmt_master_cache$coverage_dropped). Shown in the DE
    # module's "Not Tested & Low Coverage" tab. Nothing here enters the model.
    filtered_df <- NULL
    .pm_all <- tmt_master_cache$protein_matrix
    if (!is.null(.pm_all) && !is.null(tmt_master_cache$annotation_map)) {
      .gf_drop <- setdiff(rownames(.pm_all), rownames(.prep$matrix))
      .gf <- tmt_coverage_profile(.pm_all, tmt_master_cache$annotation_map, .gf_drop,
                                  acc_to_gene %||% tmt_master_cache$gene_names)
      if (!is.null(.gf)) {
        # best-covered condition, as the filter saw it (sample channels only)
        .conds <- unique(design_df$condition)
        .ch_of <- lapply(stats::setNames(.conds, .conds),
                         function(cn) design_df$sample_name[design_df$condition == cn])
        .nob <- vapply(.ch_of, function(ch) rowSums(!is.na(.pm_all[.gf$Protein, ch, drop = FALSE])),
                       numeric(nrow(.gf)))
        if (is.null(dim(.nob))) .nob <- matrix(.nob, nrow = nrow(.gf), dimnames = list(NULL, .conds))
        .fr <- sweep(.nob, 2, lengths(.ch_of), "/")
        .bi <- max.col(.fr, ties.method = "first")
        .ix <- cbind(seq_len(nrow(.gf)), .bi)
        .best <- paste0(.conds[.bi], " ", .nob[.ix], "/", lengths(.ch_of)[.bi], " (", round(100 * .fr[.ix]), "%)")
        # no value in ANY annotated channel (pooled channels included): never quantified
        .ann_cols <- intersect(tmt_master_cache$annotation_map$ColumnName, colnames(.pm_all))
        .none <- rowSums(!is.na(.pm_all[.gf$Protein, .ann_cols, drop = FALSE])) == 0
        .gf <- data.frame(.gf[1:2], Stage = ifelse(.none, TMT_STAGE_NO_VALUES, "Before DE (missingness filter)"),
                          Reason = ifelse(.none, TMT_REASON_NO_VALUES,
                                          paste0("Fewer than 50% of channels observed in every condition (best: ", .best,
                                                 "); at least 50% in one condition required")),
                          .gf[-(1:2)], check.names = FALSE, stringsAsFactors = FALSE)
      }
      .cd <- tmt_master_cache$coverage_dropped
      filtered_df <- if (!is.null(.cd) && !is.null(.gf) && identical(names(.cd), names(.gf))) rbind(.cd, .gf)
                     else if (!is.null(.gf)) .gf else .cd
      if (!is.null(filtered_df)) {
        # quantified-but-sparse proteins first, never-quantified ones last
        filtered_df <- filtered_df[order(match(filtered_df$Stage, TMT_FILTER_STAGES), filtered_df$Protein), , drop = FALSE]
        rownames(filtered_df) <- NULL
      }
    }


    # ---- Per-pair coverage flag (reporting only) ---------------------------
    # A protein solid in one arm and sparse in the other still has an estimable
    # contrast and a genuine p-value; its weakness is a large standard error
    # from few observations, which is exactly what limma's per-protein
    # df.residual accounting handles. Nothing is removed from the fit here.
    # This block only measures per-pair coverage so low-coverage proteins can be
    # flagged in the results table and tabulated for the user.
    one_arm_df <- NULL
    if (!is.null(pm_pre_impute)) {
      gv_p  <- design_df$condition[match(colnames(pm_pre_impute),
                                         design_df$sample_name)]
      i_t   <- which(gv_p == cond_test)
      i_r   <- which(gv_p == cond_ref)
      f_t   <- rowSums(!is.na(pm_pre_impute[, i_t, drop = FALSE])) / length(i_t)
      f_r   <- rowSums(!is.na(pm_pre_impute[, i_r, drop = FALSE])) / length(i_r)
      ok_t  <- f_t >= 0.5
      ok_r  <- f_r >= 0.5
      testable <- ok_t & ok_r
      one_arm  <- xor(ok_t, ok_r)

      if (any(one_arm)) {
        .oa <- rownames(pm_pre_impute)[one_arm]
        # Zero-origin gaps within the contrasted channels: the reporter ion was
        # present but non-positive, i.e. at or below the noise floor of an
        # otherwise-quantified spectrum. Distinguishes an intensity-driven gap
        # from a gap where no contributing peptide was measured at all. NA when
        # no mask is available (limma_protein route, older bundles).
        .zm_oa <- if (!is.null(zmask_g)) {
          rowSums(zmask_g[one_arm, c(i_t, i_r), drop = FALSE])
        } else {
          rep(NA_integer_, sum(one_arm))
        }
        # Peptide count per protein, when the cache carries it. Single-peptide
        # proteins are mechanically over-represented here and should be read
        # with that in mind.
        .pc_oa <- if (!is.null(tmt_cache$pep_count)) {
          unname(tmt_cache$pep_count[.oa])
        } else {
          rep(NA_integer_, length(.oa))
        }
        one_arm_df <- data.frame(
          Protein       = .oa,
          Gene          = unname((tmt_cache$acc_to_gene %||%
                                    setNames(.oa, .oa))[.oa]),
          n_valid_test  = rowSums(!is.na(pm_pre_impute[one_arm, i_t, drop = FALSE])),
          n_valid_ref   = rowSums(!is.na(pm_pre_impute[one_arm, i_r, drop = FALSE])),
          n_test        = length(i_t),
          n_ref         = length(i_r),
          DetectedIn    = ifelse(ok_t[one_arm], cond_test, cond_ref),
          n_zero_origin = .zm_oa,
          n_peptides    = .pc_oa,
          Comparison    = paste0(cond_test, " vs ", cond_ref),
          stringsAsFactors = FALSE
        )
        rownames(one_arm_df) <- NULL
      }

      logger_fn(paste0("[DE TMT] Pair coverage (", cond_test, " vs ",
                       cond_ref, "): ", sum(testable),
                       " with >=50% valid in both arms, ",
                       sum(one_arm), " flagged low coverage (adequate in one arm ",
                       "only - still fitted and tested), ",
                       sum(!testable & !one_arm), " sparse in both."))
      if (sum(testable) < 2L)
        stop(paste0("Fewer than 2 proteins are testable for ", cond_test,
                    " vs ", cond_ref, "."))
    }

    # ---- Per-condition variance diagnostic ---------------------------------
    # fit_all_groups = TRUE pools residual variance across conditions. That is
    # the standard limma workflow and is what runs, but it assumes a common
    # residual variance across groups. Measure the assumption and say so;
    # arrayWeights (below) is what compensates when it is violated.
    # Warn, never stop: nothing in the code branches on this value.
    residual_var_ratio <- NA_real_
    if (!is.null(pm_pre_impute)) {
      .gvv <- design_df$condition[match(colnames(pm_pre_impute),
                                        design_df$sample_name)]
      .cc_rows <- stats::complete.cases(pm_pre_impute)
      if (sum(.cc_rows) >= 20L) {
        .cm <- pm_pre_impute[.cc_rows, , drop = FALSE]
        .sd_by_cond <- vapply(
          unique(.gvv[!is.na(.gvv) & nzchar(.gvv)]),
          function(cc) {
            .ix <- which(.gvv == cc)
            if (length(.ix) < 2L) return(NA_real_)
            stats::median(apply(.cm[, .ix, drop = FALSE], 1L, stats::sd,
                                na.rm = TRUE), na.rm = TRUE)
          }, numeric(1))
        .sd_ok <- .sd_by_cond[is.finite(.sd_by_cond) & .sd_by_cond > 0]
        if (length(.sd_ok) >= 2L) {
          residual_var_ratio <- max(.sd_ok) / min(.sd_ok)
          logger_fn(paste0(
            "[DE TMT] Per-condition residual SD (median over ", sum(.cc_rows),
            " complete proteins): ",
            paste0(names(.sd_ok), "=", signif(.sd_ok, 3), collapse = ", "),
            " | max/min = ", signif(residual_var_ratio, 3), ". ",
            if (residual_var_ratio < 2)
              paste0("Under 2x is unremarkable; pooling residual variance across ",
                     "conditions is well justified here.")
            else if (residual_var_ratio <= 3)
              paste0("Between 2x and 3x - worth a sentence in the methods ",
                     "section; sample weights are absorbing the difference.")
            else
              paste0("Above 3x - the array weights are doing substantial work ",
                     "and that should be stated explicitly in the methods.")))
        }
      }
    }

    # ---- Pre-model design assertions ---------------------------------------
    # WARNINGS, NEVER STOPS. Real datasets in use include an n = 1 condition and
    # must still run. The point is that the consequence is named in plain terms
    # at the moment it applies, not that the analysis is blocked.
    design_warn_low_replication  <- FALSE
    design_warn_batch_confounded <- FALSE
    design_warn_low_residual_df  <- FALSE
    {
      .dcond <- design_df$condition[!is.na(design_df$condition) &
                                      nzchar(design_df$condition)]
      .ncond <- table(.dcond)

      # 1. Replicate count.
      .thin <- names(.ncond)[.ncond < 2L]
      if (length(.thin) > 0L) {
        design_warn_low_replication <- TRUE
        .osmss_user_warn(paste0(
          "Condition", if (length(.thin) > 1L) "s" else "", " ",
          paste(.thin, collapse = ", "), " ",
          if (length(.thin) > 1L) "have" else "has",
          " fewer than 2 samples. Such a condition contributes a coefficient to ",
          "the model but zero residual degrees of freedom, so its moderated ",
          "variance leans entirely on the empirical Bayes prior borrowed from ",
          "the other conditions. Fold changes involving it are estimable but ",
          "their standard errors are not supported by replication."))
      }

      # 2. Batch estimability: is each condition seen in >=2 distinct plexes?
      if (isTRUE(batch_for_model) && !is.null(batch_groups)) {
        .bmap <- setNames(
          rep(names(batch_groups), times = vapply(batch_groups, length, integer(1))),
          unlist(batch_groups))
        .pl_per_cond <- tapply(
          .bmap[design_df$sample_name], design_df$condition,
          function(z) length(unique(z[!is.na(z)])))
        .conf <- names(.pl_per_cond)[!is.na(.pl_per_cond) & .pl_per_cond < 2L]
        if (length(.conf) > 0L) {
          design_warn_batch_confounded <- TRUE
          .osmss_user_warn(paste0(
            "Condition", if (length(.conf) > 1L) "s" else "", " ",
            paste(.conf, collapse = ", "),
            " appear", if (length(.conf) == 1L) "s" else "",
            " in only one plex. The plex covariate is therefore partially ",
            "confounded with condition, and contrasts involving ",
            if (length(.conf) > 1L) "these conditions" else "this condition",
            " may be biased: the model cannot separate the plex effect from the ",
            "biological effect. The covariate is still fitted."))
        }
      }

      # 3. Parameter count vs samples.
      .n_par <- length(unique(.dcond)) +
                (if (isTRUE(batch_for_model) && !is.null(batch_groups))
                   max(0L, length(batch_groups) - 1L) else 0L)
      .n_smp <- nrow(design_df)
      if (.n_par > (.n_smp - 2L)) {
        design_warn_low_residual_df <- TRUE
        .osmss_user_warn(paste0(
          "The design estimates ", .n_par, " parameters from ", .n_smp,
          " samples, leaving ", .n_smp - .n_par,
          " residual degrees of freedom. Variance estimates at this point rest ",
          "almost entirely on the empirical Bayes prior rather than on the data ",
          "for each protein."))
      }
    }

    # dia_limma_de() is shared with the DIA route. No TMT route imputes inside
    # it any more: opt 1 imputes in the contrast-independent global block above,
    # opt 3 in the within-plex block, and opt 2 does not impute at all.
# Imputation is now performed once, contrast-independently, in the global
    # block above (single-plex) or the within-plex block (opt 3). Multi-plex
    # never imputed here. max_missing_per_cond = 0L makes dia_limma_de's
    # eligibility test unsatisfiable, so its internal k-NN cannot fire on any
    # TMT route. Where the global pre-pass has run, dia_limma_de's own filter is
    # switched off entirely (see the argument comments below): the pre-pass owns
    # the protein universe and the fit spans every condition. Where it has not
    # run, require_both = TRUE remains the operative per-contrast filter.
    limma_res <- dia_limma_de(
      cache                    = tmt_cache,
      cond_test                = cond_test,
      cond_ref                 = cond_ref,
      # The spine's filter has already defined the protein universe on EVERY
      # route, so dia_limma_de's own filter is deliberately made a no-op (a
      # fraction >= 0 is always TRUE). That keeps the fitted set identical
      # across contrasts, which is the whole point of the pre-pass. The spine
      # guarantees every surviving row has >=50% valid in at least one
      # condition, so no degenerate all-NA row reaches lmFit.
      min_valid_frac           = 0,
      require_both             = FALSE,
      fit_all_groups           = TRUE,
      max_missing_per_cond     = 0L,
      knn_k                    = 10L,
      batch_groups             = batch_groups,
      batch_for_model          = batch_for_model,
      observed_mat             = pm_pre_impute,
      use_array_weights        = TRUE,
      logger                   = logger_fn
    )
    if (!is.null(limma_res$error)) stop(limma_res$error)

    rdf <- limma_res$result_df

    # ---- Low-coverage flag (reporting) -------------------------------------
    # one_arm_df is a reporting table, not an exclusion list. The proteins in it
    # are in the results table with real p-values; the flag travels with the
    # results so a reader can see which hits rest on one well-covered arm.
    rdf$LowCoverage <- if (!is.null(one_arm_df))
                         as.character(rdf$Protein) %in% one_arm_df$Protein
                       else rep(FALSE, nrow(rdf))

    # n_tested must count rows that actually produced a test, not rows handed to
    # lmFit. A protein with no observations in an arm returns NA and leaves the
    # BH universe on its own; it should not be counted as tested. This is the
    # same quantity as Diag_estimable below - computed once, used by both.
    limma_res$n_tested <- sum(!is.na(rdf$pvalue))
    logger_fn(paste0("[DE TMT] Fitting set: ", nrow(rdf),
                     " rows entered the fit, ", limma_res$n_tested,
                     " returned an estimable contrast, ",
                     sum(rdf$LowCoverage), " flagged low coverage."))

    # ── Variance-moderation diagnostics ──────────────────────────────────────
    # eBayes shares variance across proteins, so changing WHICH proteins enter
    # the model changes every protein's p-value even when its own data is
    # identical. These columns expose that so pipeline versions can be compared.
    .fit <- limma_res$fit_object
    if (!is.null(.fit) && !is.null(.fit$coefficients)) {
      .ids <- rownames(.fit$coefficients)
      .idcol <- if ("Protein" %in% names(rdf)) "Protein" else
                if ("Accession" %in% names(rdf)) "Accession" else names(rdf)[1]
      .i <- match(as.character(rdf[[.idcol]]), .ids)

      # trend=TRUE makes s2.prior a vector; robust=TRUE makes df.prior a vector.
      .vec <- function(x, n) if (length(x) == 1L) rep(x, n) else x
      .n   <- length(.ids)

      rdf$Diag_sigma       <- .fit$sigma[.i]
      rdf$Diag_df_residual <- .fit$df.residual[.i]
      rdf$Diag_df_prior    <- .vec(.fit$df.prior,  .n)[.i]
      rdf$Diag_s2_prior    <- .vec(.fit$s2.prior,  .n)[.i]
      rdf$Diag_s2_post     <- if (!is.null(.fit$s2.post)) .fit$s2.post[.i] else NA_real_
      rdf$Diag_t           <- if (!is.null(.fit$t)) .fit$t[.i, 1] else NA_real_
      rdf$Diag_n_obs       <- rowSums(!is.na(limma_res$observed_mat))[
                                match(as.character(rdf[[.idcol]]),
                                      rownames(limma_res$observed_mat))]
      # TRUE where the contrast was estimable and a p-value was returned. Rows
      # with no observations at all in one compared arm return NA and are absent
      # from the BH universe; this column makes that visible in the results.
      rdf$Diag_estimable   <- !is.na(rdf$pvalue)

      # Summary to the UI log panel
      .dfp <- .vec(.fit$df.prior, .n)
      .s2p <- .vec(.fit$s2.prior, .n)
      logger_fn(paste0(
        "[DE TMT DIAG] proteins in model: ", .n,
        " | df.prior median: ", signif(median(.dfp, na.rm = TRUE), 4),
        " | s2.prior median: ", signif(median(.s2p, na.rm = TRUE), 4),
        " | df.residual range: ", paste(range(.fit$df.residual, na.rm = TRUE), collapse = "-"),
        " | df.residual <= 2: ", sum(.fit$df.residual <= 2, na.rm = TRUE),
        " | NA adj.p: ", sum(is.na(rdf$adj.pvalue))))
    }
    # ---- Restore pre-imputation observation counts (ALL routes) ------------
    # dia_limma_de builds raw_mat, Diag_n_obs and n_valid_<condition> from the
    # matrix it was handed. On every TMT route that matrix may already have been
    # imputed upstream - globally on opt 1, within-plex on opt 3 - so those
    # counts would report imputed cells as observations.
    # Recompute them from pm_pre_impute, which is the post-filter,
    # pre-imputation matrix. Must run AFTER the diagnostics block that assigns
    # Diag_n_obs and BEFORE plot_df is built.
    if (!is.null(pm_pre_impute)) {
      .idc <- if ("Protein" %in% names(rdf)) "Protein" else names(rdf)[1]
      .ids <- as.character(rdf[[.idc]])
      .m   <- match(.ids, rownames(pm_pre_impute))
      .gvq <- design_df$condition[match(colnames(pm_pre_impute),
                                        design_df$sample_name)]
      for (grp in unique(.gvq[!is.na(.gvq) & nzchar(.gvq)])) {
        .cols <- which(.gvq == grp)
        rdf[[paste0("n_valid_", make.names(grp))]] <-
          rowSums(!is.na(pm_pre_impute[, .cols, drop = FALSE]))[.m]
      }
      rdf$Diag_n_obs <- rowSums(!is.na(pm_pre_impute))[.m]
      logger_fn("[DE TMT] n_valid_* and Diag_n_obs recomputed pre-imputation.")
    }
    rdf$GeneName  <- rdf$Gene
    rdf$Accession <- rdf$Protein
    rdf$Log2FC    <- rdf$log2FC
    rdf$AdjPvalue <- rdf$adj.pvalue
    # PiScore = |log2FC| x -log10(adjusted p). This is an OSMSS ranking score, a
    # MODIFIED pi-value. It is NOT the Xiao et al. definition, which uses the
    # unadjusted p-value and a different weighting. Do not describe it as a
    # "pi-value" without qualification. Kept as computed: recomputing would
    # change every ranking and plot ordering across the app.
    rdf$PiScore   <- abs(rdf$log2FC) * (-log10(pmax(rdf$adj.pvalue, 1e-300)))

    # ── Single-peptide protein flagging (▼ in volcano) ──
    pep_counts <- tmt_master_cache$pep_count_per_protein
    if (!is.null(pep_counts) && length(pep_counts) > 0) {
      rdf$SinglePeptide <- pep_counts[rdf$Protein] == 1L
      rdf$SinglePeptide[is.na(rdf$SinglePeptide)] <- FALSE
    } else {
      rdf$SinglePeptide <- FALSE
    }

    n_sig <- sum(rdf$adj.pvalue < 0.05, na.rm = TRUE)

    # ---- Single-peptide sensitivity (reporting only) -----------------------
    # Single-peptide quantification is noisier, and the volcano already flags
    # it. This says how much of the significant set rests on it, so the reader
    # can judge without the tested set being changed. Not an exclusion rule:
    # excluding them would change which proteins are tested and is a separate
    # decision. The field-standard treatment for peptide-count-dependent
    # variance is DEqMS, which would be a new dependency.
    n_tested_single_peptide      <- sum(rdf$SinglePeptide, na.rm = TRUE)
    n_significant_single_peptide <- sum(rdf$SinglePeptide &
                                        rdf$adj.pvalue < 0.05, na.rm = TRUE)
    logger_fn(paste0(
      "[DE TMT] Single-peptide sensitivity: ", n_significant_single_peptide,
      " of ", n_sig, " significant proteins are single-peptide (",
      n_tested_single_peptide, " single-peptide proteins tested)."))
    plot_df <- data.frame(
      Accession     = rdf$Protein,
      GeneName      = rdf$Gene,
      Log2FC        = rdf$log2FC,
      AdjPValue     = rdf$adj.pvalue,
      NegLog10P     = -log10(pmax(rdf$adj.pvalue, 1e-300)),
      IsCapped      = FALSE,
      SinglePeptide = rdf$SinglePeptide,
      stringsAsFactors = FALSE
    )

    # ── Attach diagnostics safely ────────────────────────────────────────────
    # Any of these can be NULL or zero-length depending on what limma returned,
    # which is fatal inside data.frame(). Normalise length before assigning.
    .fit_col <- function(x, n) {
      if (is.null(x) || length(x) == 0L) return(rep(NA_real_, n))
      if (length(x) == 1L)               return(rep(as.numeric(x), n))
      if (length(x) != n)                return(rep(NA_real_, n))
      as.numeric(x)
    }
    .n_row <- nrow(plot_df)
    for (.cn in c("PValue", "AveExpr", "Diag_sigma", "Diag_df_residual",
                  "Diag_df_prior", "Diag_s2_prior", "Diag_s2_post",
                  "Diag_t", "Diag_n_obs")) {
      .src <- switch(.cn,
                     # dia_limma_de() names the raw p-value column "pvalue"
                     # (limma's topTable name "P.Value" is not carried over).
                     "PValue"  = rdf[["pvalue"]],
                     "AveExpr" = rdf[["AveExpr"]],
                     rdf[[.cn]])
      plot_df[[.cn]] <- .fit_col(.src, .n_row)
    }
    # Assigned outside the loop above on purpose: .fit_col() coerces with
    # as.numeric(), which would turn this logical into 0/1.
    plot_df$Diag_estimable <- if (!is.null(rdf$Diag_estimable) &&
                                  length(rdf$Diag_estimable) == .n_row)
                                rdf$Diag_estimable else rep(NA, .n_row)
    logger_fn(paste0("[DE TMT DIAG] column lengths - ",
                     paste(vapply(c("pvalue","AveExpr","Diag_sigma","Diag_df_residual",
                                    "Diag_df_prior","Diag_s2_prior","Diag_s2_post",
                                    "Diag_t","Diag_n_obs"),
                                  function(k) paste0(k, "=", length(rdf[[k]])),
                                  character(1)), collapse = ", "),
                     " | nrow=", .n_row))
    logger_fn(paste0(
      "[DE TMT] Protein reconciliation: ", n_proteins_before, " summarized -> ",
      n_global_kept %||% NA, " after global filter -> ", nrow(rdf),
      " fitted -> ", limma_res$n_tested %||% NA,
      " with an estimable contrast (the remainder have no observations in one ",
      "of the compared arms and return NA). Of the fitted rows, ",
      if (!is.null(one_arm_df)) nrow(one_arm_df) else 0L,
      " are flagged low coverage - tested, not excluded."))
    prov_tv <- list(
      analysis_mode          = tmt_analysis_route %||% "limma",
      design_option          = design_opt,
      comparison             = limma_res$comparison,
      n_proteins                 = limma_res$n_tested,
      n_significant              = n_sig,
      n_proteins_dropped_missing = n_proteins_before - (limma_res$n_tested %||% n_proteins_before),
      software_used              = "limma",
      # Summarization-level metadata (pulled from cache for methods text)
      data_level             = tmt_master_cache$route %||% tmt_analysis_route,
      interference_filtered  = tmt_master_cache$interference_filtered  %||% FALSE,
      interference_threshold = tmt_master_cache$interference_threshold %||% 50L,
      n_peptides_input       = tmt_master_cache$n_peptides_input,
      within_plex_knn        = identical(design_opt, 3L) && isTRUE(batch_for_model),
      knn_k                      = 10L,
      cross_plex_impute_disabled = !identical(design_opt, 1L),
      n_plexes                   = if (!is.null(batch_groups)) length(batch_groups) else 1L,
      # Pooled reference (Norm) channels per plex, range over plexes (design 2).
      n_bridge_per_plex          = local({
        .a <- tmt_master_cache$annotation_map
        if (is.null(.a) || !all(c("Condition", "Mixture") %in% colnames(.a)) ||
            !any(.a$Condition == "Norm", na.rm = TRUE)) NA_integer_
        else range(tapply(.a$Condition == "Norm", .a$Mixture, sum, na.rm = TRUE))
      }),
      # Plex-coverage filter at summarization (protein route, multi-plex) and the
      # number of proteins listed as not tested (both filters).
      n_proteins_plex_coverage_dropped = tmt_master_cache$n_dropped_coverage %||% 0L,
      plex_coverage_min_plexes   = tmt_master_cache$min_plexes_required %||% NA_integer_,
      plex_coverage_n_plexes     = tmt_master_cache$n_plexes %||% NA_integer_,
      n_proteins_not_tested_listed = if (is.null(filtered_df)) 0L else nrow(filtered_df),
      n_conditions_model         = length(unique(design_df$condition)),
      df_residual                = if (!is.null(limma_res$fit_object))
                                     stats::median(limma_res$fit_object$df.residual, na.rm = TRUE)
                                   else NA_integer_,
      imputation_k               = 10L,
      # Contrast-independent pre-pass counts. These reconcile as:
      #   n_proteins_before -> n_proteins_global_filter -> n_proteins (estimable),
      # the second step accounting only for rows with no observations at all in
      # one compared arm. n_proteins_one_arm_excluded is a legacy key name kept
      # for existing consumers; nothing is excluded by it any more - it counts
      # the low-coverage flag, which is reported alongside the results.
      global_filter_applied      = isTRUE(global_filter_applied),
      n_proteins_global_filter   = n_global_kept,
      n_proteins_one_arm_excluded = if (!is.null(one_arm_df)) nrow(one_arm_df) else 0L,
      n_proteins_unestimable     = sum(!rdf$Diag_estimable),
      n_proteins_low_coverage    = sum(rdf$LowCoverage),
      n_values_imputed_global    = n_global_filled,
      n_values_imputed_within_plex = n_values_imputed_within_plex,
      plex_min_obs_frac          = plex_min_obs_frac,
      route_imputation           = route_imputation,
      n_zero_origin_held         = n_zero_held,
      n_channels_over_colmax     = n_channels_over_colmax,
      array_weights_applied      = isTRUE(limma_res$array_weights_applied),
      array_weights              = limma_res$array_weights %||% NA_real_,
      residual_var_ratio         = residual_var_ratio,
      # Silent-degradation counters (T12). Each records a fallback that was
      # previously logged but absent from the reported provenance.
      n_array_weights_failed     = limma_res$n_array_weights_failed %||% 0L,
      n_zero_mask_dropped        = .prep$n_zero_mask_dropped,
      n_plex_knn_failed          = .prep$n_plex_knn_failed,
      n_agg_fallback_colmeans    = tmt_master_cache$n_agg_fallback_colmeans %||% NA_integer_,
      n_agg_singleton            = tmt_master_cache$n_agg_singleton %||% NA_integer_,
      peptide_fdr_source         = tmt_master_cache$peptide_fdr_source %||%
                                     "inherited_from_search_engine",
      n_tested_single_peptide      = n_tested_single_peptide,
      n_significant_single_peptide = n_significant_single_peptide,
      piscore_definition           = "abs_log2FC_x_neglog10_adjP",
      design_warn_low_replication  = design_warn_low_replication,
      design_warn_batch_confounded = design_warn_batch_confounded,
      design_warn_low_residual_df  = design_warn_low_residual_df,
      n_interf_kept_na           = tmt_master_cache$n_interf_kept_na   %||% NA_integer_,
      n_interf_kept_pass         = tmt_master_cache$n_interf_kept_pass %||% NA_integer_,
      # Shared-peptide provenance, carried from the summarization cache
      # (tmt_limma_qf_summarize). render_tmt_volcano_methods() reads these; without
      # them it fell back to "keep_all" and reported shared peptides as retained
      # even when they had been discarded. A peptide-route cache without the key
      # predates the setting and kept every peptide, so it still reads "keep_all";
      # the protein route has no peptides and gets NA (sentence omitted).
      shared_peptide_handling    = tmt_master_cache$shared_peptide_handling %||%
                                     (if (identical(tmt_master_cache$route, "limma_qf"))
                                        "keep_all" else NA_character_),
      n_peptides_total           = tmt_master_cache$n_peptides_total          %||% NA_integer_,
      n_peptides_shared          = tmt_master_cache$n_peptides_shared         %||% NA_integer_,
      n_peptides_removed_shared  = tmt_master_cache$n_peptides_removed_shared %||% NA_integer_,
      shared_detection_method    = tmt_master_cache$shared_detection_method   %||% NA_character_,
      n_peptides_razor_assigned  = tmt_master_cache$n_peptides_razor_assigned %||% NA_integer_,
      n_peptides_dropped_no_unique_evidence =
                                   tmt_master_cache$n_peptides_dropped_no_unique_evidence %||% NA_integer_,
      cross_plex_normalization   = tmt_master_cache$cross_plex_normalization,
      # Channel normalization (tmt_normalization_check(), upload wizard step 6).
      # A cache without the key predates the option and was centered.
      normalization_applied      = tmt_master_cache$normalization_applied %||% TRUE,
      input_normalized_upstream  = tmt_master_cache$input_normalized_upstream %||% NA,
      norm_check_verdict         = tmt_master_cache$norm_check_verdict    %||% NA_character_,
      norm_check_within_sd       = tmt_master_cache$norm_check_within_sd  %||% NA_real_,
      norm_check_n_complete      = tmt_master_cache$norm_check_n_complete %||% NA_integer_,
      # Channels marked "Exclude" in the design grid (never in the matrix).
      excluded_channels          = tmt_master_cache$excluded_channels %||% character(0),
      impute_pkg_available       = impute_pkg_ok,
      timestamp              = as.character(Sys.time())
    )
    prov_de <- list(
      source        = "tmt",
      engine        = "limma",
      analysis_mode = tmt_analysis_route %||% "limma",
      comparison    = limma_res$comparison,
      n_tested      = limma_res$n_tested,
      n_sig         = n_sig,
      timestamp     = as.character(Sys.time())
    )
    logger_fn(paste0("[DE TMT] limma done: ", limma_res$n_tested, " proteins, ",
                     n_sig, " significant (design_option=", design_opt, ")."))
    list(
      route                  = "limma",
      result_df              = rdf,
      one_arm_df             = one_arm_df,
      filtered_df            = filtered_df,
      plot_data              = plot_df,
      analysis_mode          = "tmt",
      cond_test              = cond_test,
      cond_ref               = cond_ref,
      provenance_tmt_volcano = prov_tv,
      provenance_de          = prov_de
    )
  }, error = function(e) list(error = e$message))

  if (!is.null(result$error)) {
    return(list(success = FALSE, error = result$error, messages = msgs))
  }
  c(list(success = TRUE, error = NULL, messages = msgs), result)
}


# === REUSABLE UI/SERVER MODULES ===============================================
source("module_methods_text_ui.R",     local = FALSE)
source("module_methods_text_server.R", local = FALSE)

# =============================================================================
# DIA LABEL-FREE FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# detect_dia_upload_type()
# Which DIA protein-level export a file is, from its header row alone. The
# rules use exact column names and are checked in this order:
#  1. Spectronaut: one of Spectronaut's protein-group ID columns,
#     PG.ProteinGroups or PG.ProteinAccessions, bare or with Perseus' "T: "
#     prefix. DIA-NN never writes these (it uses Protein.Group, Protein.Ids),
#     although its report.tsv has PG.Quantity and PG.Q.Value, so a "PG." prefix
#     alone would not tell the two apart. Within Spectronaut:
#     a. T: PG.ProteinGroups / T: PG.ProteinAccessions -> Pivot Report exported
#        through Perseus
#     b. any EG.PTM* or EG.TotalQuantity column -> PTM report (DIA PTM upload);
#        a protein-group report has no precursor-level (EG.) columns
#     c. R.FileName or a bare PG.Quantity column -> long ("normal") report,
#        rejected. A Pivot Report names the run inside each quantity column
#        ("[1] <run>.PG.Quantity"), so these columns occur only in long reports.
#     d. "[N] <run>.PG.Quantity" columns -> Pivot Report
#  2. DIA-NN: Protein.Group (Spectronaut never writes it). Within DIA-NN:
#     a. Run or File.Name -> main report (one row per precursor and run), rejected
#     b. Precursor.Id -> precursor matrix (report.pr_matrix.tsv), rejected
#     c. otherwise -> protein-group matrix (report.pg_matrix.tsv)
#  3. A label-free DDA file (detect_lf_upload_type) -> pointed to that upload
# Before these, files that have no header row: DIA-NN's parquet main report
# (starts with the bytes "PAR1"), and the two optional second files asked for
# after a DIA-NN matrix: report.log.txt (first line "DIA-NN <version>
# (Data-Independent Acquisition by Neural Networks)") and FragPipe's
# fragpipe.workflow (first line "# FragPipe (<version>) ..." or "# Workflow: ...").
# FragPipe's DIA output is DIA-NN's: dia-quant-output/report.pg_matrix.tsv.
# Returns list(kind, message, headers).
#   kind: "spectronaut_pivot" | "spectronaut_perseus" | "diann_pg_matrix" | "unsupported"
#   message: user-facing explanation; NULL for the accepted kinds.
# -----------------------------------------------------------------------------
detect_dia_upload_type <- function(path) {
  diann_main <- paste0(
    "This is DIA-NN’s main report, with one row per precursor in each run. OSMSS needs ",
    "DIA-NN’s protein-group matrix, report.pg_matrix.tsv, which DIA-NN saves next to the main ",
    "report (DIA-NN 1.x: tick ‘Generate matrices’ / --matrices). Upload report.pg_matrix.tsv, ",
    "and optionally report.log.txt from the same run.")
  magic <- tryCatch(readBin(path, "raw", 4L), error = function(e) raw(0))
  if (identical(magic, charToRaw("PAR1")))
    return(list(kind = "unsupported", message = diann_main, headers = character(0)))
  hdr  <- lf_read_header(path)
  .out <- function(kind, message = NULL) list(kind = kind, message = message, headers = hdr)
  if (length(hdr) == 0L)
    return(.out("unsupported", "The file is empty or its header row could not be read."))
  has  <- function(x)  any(x %in% hdr)
  hasr <- function(re) any(grepl(re, hdr))
  accepted <- paste0("The DIA upload accepts a Spectronaut protein-group Pivot Report ",
                     "(PG.Quantity as the cell value), directly or exported through Perseus, ",
                     "or DIA-NN’s report.pg_matrix.tsv.")
  if (grepl("^DIA-NN .*Data-Independent Acquisition by Neural Networks", hdr[1]))
    return(.out("unsupported", paste0(
      "This is DIA-NN’s log (report.log.txt). Upload report.pg_matrix.tsv first; ",
      "the log is asked for after it.")))
  if (grepl("^#\\s*(FragPipe|Workflow:)", hdr[1]))
    return(.out("unsupported", paste0(
      "This is FragPipe’s fragpipe.workflow. For a FragPipe DIA run, upload ",
      "dia-quant-output/report.pg_matrix.tsv first; the workflow is asked for after it.")))

  # 1. Spectronaut
  if (has(c("PG.ProteinGroups", "PG.ProteinAccessions",
            "T: PG.ProteinGroups", "T: PG.ProteinAccessions"))) {
    if (has(c("T: PG.ProteinGroups", "T: PG.ProteinAccessions")))
      return(.out("spectronaut_perseus"))
    if (hasr("EG\\.PTM|EG\\.TotalQuantity"))
      return(.out("unsupported", paste0(
        "This looks like a Spectronaut PTM report. Use the Label Free PTM upload for it; ",
        "this upload takes protein-group quantities. ", accepted)))
    if (has(c("R.FileName", "PG.Quantity")))
      return(.out("unsupported", paste0(
        "This is a Spectronaut long-format (normal) report, with one row per precursor or ",
        "fragment in each run. OSMSS needs the protein-group Pivot Report instead. Please open ",
        "the .sne file in Spectronaut (the Spectronaut Viewer is free) and export the Pivot Report ",
        "with PG.Quantity as the cell value, or export with the OSMSS_Global.rs report scheme ",
        "provided by the OSMSS authors. Alternatively, ask your proteomics core for the Pivot ",
        "Report or for the original Spectronaut (.sne) file.")))
    if (hasr("^(\\[\\d+\\] )?.+\\.PG\\.Quantity$"))
      return(.out("spectronaut_pivot"))
    return(.out("unsupported", paste0(
      "This is a Spectronaut export, but it has no per-run PG.Quantity columns ",
      "(“[1] <run>.PG.Quantity”). Export the Pivot Report with PG.Quantity as the cell value.")))
  }

  # 2. DIA-NN
  if (has("Protein.Group")) {
    if (has(c("Run", "File.Name")))
      return(.out("unsupported", diann_main))
    if (has("Precursor.Id"))
      return(.out("unsupported", paste0(
        "This is DIA-NN’s precursor matrix (report.pr_matrix.tsv). Upload the protein-group ",
        "matrix, report.pg_matrix.tsv, from the same folder instead.")))
    return(.out("diann_pg_matrix"))
  }

  # 3. Label-free DDA
  if (detect_lf_upload_type(path)$kind %in%
      c("pd_nested", "mq_proteingroups", "mq_evidence", "fp_combined_protein", "fp_workflow"))
    return(.out("unsupported", paste(
      "This looks like a label-free DDA file. Use the Label Free DDA upload instead.", accepted)))

  .out("unsupported", paste("The file type could not be recognised from its column headers.", accepted))
}

# -----------------------------------------------------------------------------
# parse_spectronaut_pivot()
# Parse a Spectronaut 2D pivot report (TSV) into a list suitable for DIA
# summarization. Expects columns PG.Genes, PG.ProteinAccessions, and
# sample-level columns named [N] SampleName.PG.Quantity.
# Returns: list(protein_matrix, gene_names, accessions, sample_names,
#               n_proteins, n_samples, multi_acc_map)
# -----------------------------------------------------------------------------
parse_spectronaut_pivot <- function(file_path, software = "spectronaut",
                                    dedup_order = c("conditional", "always"),
                                    logger = message) {
  dedup_order <- match.arg(dedup_order)
  logger(paste0("[DIA] Parsing Spectronaut pivot file: ", basename(file_path)))

  # 1. Read header to identify columns
  header_raw <- readLines(file_path, n = 1, warn = FALSE)
  all_cols   <- strsplit(header_raw, "\t")[[1]]

  gene_col_idx  <- which(all_cols == "PG.Genes")[1]
  acc_col_idx   <- which(all_cols == "PG.ProteinAccessions")[1]

  if (is.na(gene_col_idx))
    return(list(error = "Column 'PG.Genes' not found in file header."))
  if (is.na(acc_col_idx))
    return(list(error = "Column 'PG.ProteinAccessions' not found in file header."))

  # 2. Find quantity columns: ends with .PG.Quantity
  quant_idx  <- grep("\\.PG\\.Quantity$", all_cols)
  if (length(quant_idx) < 2)
    return(list(error = paste0("Fewer than 2 PG.Quantity columns found. ",
                               "Ensure the Spectronaut export uses Pivot Report format ",
                               "with PG.Quantity as the cell value.")))

  # 3. Extract sample names: strip leading [N] prefix and .PG.Quantity suffix
  raw_col_names <- all_cols[quant_idx]
  sample_names  <- sub("^\\.PG\\.Quantity$", "",
                   sub("^\\[\\d+\\]\\s*", "", raw_col_names))
  sample_names  <- sub("\\.PG\\.Quantity$", "", sample_names)

  # 4. Read only needed columns via data.table::fread
  cols_to_read <- c(gene_col_idx, acc_col_idx, quant_idx)
  dt <- tryCatch(
    data.table::fread(file_path, select = cols_to_read,
                      na.strings = c("NA", "NaN", "", "Filtered"),
                      data.table = FALSE),
    error = function(e) list(error = paste0("File read error: ", conditionMessage(e)))
  )
  if (!is.null(dt$error)) return(dt)

  # 5. Rename columns for clarity
  colnames(dt)[1] <- "Gene"
  colnames(dt)[2] <- "Accession"
  colnames(dt)[seq(3, ncol(dt))] <- sample_names

  # 6. Primary accession/gene: first entry before ";"
  dt$Accession_primary <- sub(";.*", "", dt$Accession)
  dt$Gene_primary      <- sub(";.*", "", dt$Gene)

  # Handle duplicate primary accessions: keep first occurrence
  # Order by completeness before deduplicating so the retained row is the
  # better-quantified one rather than whichever appeared first in the file.
  .quant_cols <- setdiff(names(dt), c("Accession_primary", "Gene_primary",
                                      "Gene", "Accession"))
  .quant_cols <- .quant_cols[vapply(dt[.quant_cols], is.numeric, logical(1))]
  # Only reorder when there is actually a duplicate to resolve. Row order
  # propagates downstream (it is an input to impute.knn neighbour selection), and
  # this parser is shared with the DIA PTM route, so an unconditional reorder is a
  # side effect well beyond the scope of duplicate resolution.
  # dedup_order = "conditional" (default) is that behaviour, byte-identical for
  # every existing caller.  "always" sorts unconditionally, so which row
  # represents a protein group no longer depends on whether the file happens to
  # contain a duplicate; the tiebreak on Accession_primary makes it reproducible.
  if (length(.quant_cols) > 0 &&
      (identical(dedup_order, "always") ||
       anyDuplicated(dt$Accession_primary) > 0L)) {
    .completeness <- rowSums(!is.na(dt[, .quant_cols, drop = FALSE]))
    dt <- if (identical(dedup_order, "always"))
            dt[order(-.completeness, dt$Accession_primary), , drop = FALSE]
          else
            dt[order(-.completeness), , drop = FALSE]
  }
  dup_mask <- duplicated(dt$Accession_primary)
  n_dup_removed <- sum(dup_mask)
  if (n_dup_removed > 0) {
    logger(paste0("[DIA] Warning: ", n_dup_removed,
                  " protein groups sharing a leading accession were removed (keeping the most completely quantified)."))
    dt <- dt[!dup_mask, , drop = FALSE]
  }

  # 7. Build outputs
  acc_to_gene   <- setNames(dt$Gene_primary, dt$Accession_primary)
  multi_acc_map <- setNames(dt$Accession, dt$Accession_primary)

  # 8. Build numeric matrix (proteins x samples)
  mat <- as.matrix(dt[, sample_names, drop = FALSE])
  rownames(mat) <- dt$Accession_primary
  class(mat)    <- "numeric"

  # 9. Validate
  n_prots   <- nrow(mat)
  n_samples <- ncol(mat)
  if (n_prots == 0)
    return(list(error = "No protein rows found after parsing."))
  if (n_samples < 2)
    return(list(error = "Fewer than 2 samples detected."))

  pct_na <- mean(is.na(mat))
  if (pct_na > 0.80)
    logger(paste0("[DIA] Warning: ", round(pct_na * 100, 1),
                  "% of intensity values are missing. Check Q-value filtering."))

  logger(paste0("[DIA] Parse complete: ", n_prots, " proteins x ",
                n_samples, " samples."))

  list(
    protein_matrix = mat,
    gene_names     = acc_to_gene,
    accessions     = rownames(mat),
    sample_names   = sample_names,
    n_proteins     = n_prots,
    n_samples      = n_samples,
    n_dup_removed  = n_dup_removed,
    multi_acc_map  = multi_acc_map
  )
}

# -----------------------------------------------------------------------------
# parse_maxquant_global()
# Parse a Spectronaut protein-group pivot report that was exported through
# Perseus (TSV/TXT). The columns are Spectronaut's, with Perseus-style prefixes
#   (T: = text, N: = numeric, C: = categorical). Quantity columns are named
#   "[N] SampleName.PG.Quantity". Accessions in "T: PG.ProteinGroups".
# Perseus exports are usually log2; such values are returned on the linear
# scale (input_log2 = TRUE), so that downstream steps treat both exports alike.
# Returns same structure as parse_spectronaut_pivot(), plus secondary_to_primary
# and input_log2.
# -----------------------------------------------------------------------------
parse_maxquant_global <- function(file_path, logger = message) {
  logger(paste0("[MQ] Parsing MaxQuant global file: ", basename(file_path)))

  # Read header via readLines to avoid fread separator mis-detection on
  # column names that contain brackets (e.g. "[1] SampleName.PG.Quantity").
  header_raw <- readLines(file_path, n = 1L, warn = FALSE)
  all_cols   <- strsplit(header_raw, "\t")[[1]]

  # ── Locate ID columns (with or without Perseus T: prefix) ──
  acc_col_idx  <- which(all_cols %in% c("T: PG.ProteinGroups", "PG.ProteinGroups"))[1]
  gene_col_idx <- which(all_cols %in% c("T: PG.Genes", "PG.Genes"))[1]

  if (is.na(acc_col_idx))
    return(list(error = "Column 'PG.ProteinGroups' (or 'T: PG.ProteinGroups') not found in file header."))
  if (is.na(gene_col_idx))
    return(list(error = "Column 'PG.Genes' (or 'T: PG.Genes') not found in file header."))

  # ── Locate quantity columns: [N] SampleName  OR  [N] SampleName.PG.Quantity ──
  # Quantity columns are distinguished from all other [N]-numbered columns by
  # the absence of a leading prefix (N:, T:, C:). Both export formats are
  # handled: old format appends .PG.Quantity; new format uses bare names.
  # The leading [ may also be absent in some exports (1] vs [1]).
  quant_idx <- grep("^\\[?\\d+\\]", all_cols)
  if (length(quant_idx) < 2)
    return(list(error = paste0("Fewer than 2 quantity columns found. ",
                               "Expected columns starting with '[N] SampleName' ",
                               "(optionally followed by .PG.Quantity).")))

  # Exclude any columns that fread may have grouped in but have annotation prefixes
  quant_idx <- quant_idx[!grepl("^[TNC]:\\s", all_cols[quant_idx])]
  if (length(quant_idx) < 2)
    return(list(error = paste0("Fewer than 2 quantity columns found after excluding ",
                               "annotation (T:/N:/C:) columns.")))

  # ── Extract sample names ──
  raw_col_names <- all_cols[quant_idx]
  sample_names  <- sub("^\\[?\\d+\\]\\s*", "", raw_col_names)  # strip [N] or N] prefix
  sample_names  <- sub("\\.PG\\.Quantity$", "", sample_names)    # strip suffix if present
  sample_names  <- trimws(sample_names)

  # ── Read file: select only the columns we need ──
  cols_to_read <- c(acc_col_idx, gene_col_idx, quant_idx)
  dt <- tryCatch(
    data.table::fread(file_path, select = cols_to_read,
                      na.strings = c("NA", "NaN", "", "Filtered"),
                      data.table = FALSE),
    error = function(e) list(error = paste0("File read error: ", conditionMessage(e)))
  )
  if (!is.null(dt$error)) return(dt)

  colnames(dt)[1] <- "Accession"
  colnames(dt)[2] <- "Gene"
  colnames(dt)[seq(3, ncol(dt))] <- sample_names

  # ── Primary accession = first before ";" ──
  dt$Accession_primary <- sub(";.*", "", dt$Accession)
  dt$Gene_primary      <- sub(";.*", "", dt$Gene)

  # ── Remove duplicate primary accessions ──
  dup_mask <- duplicated(dt$Accession_primary)
  if (any(dup_mask)) {
    logger(paste0("[MQ] Warning: ", sum(dup_mask),
                  " duplicate primary accessions removed (keeping first)."))
    dt <- dt[!dup_mask, , drop = FALSE]
  }

  acc_to_gene   <- setNames(dt$Gene_primary, dt$Accession_primary)
  multi_acc_map <- setNames(dt$Accession, dt$Accession_primary)

  # ── Build secondary -> primary accession lookup ──
  secondary_to_primary <- character(0)
  for (i in seq_len(nrow(dt))) {
    all_accs <- trimws(strsplit(dt$Accession[i], ";")[[1]])
    primary  <- dt$Accession_primary[i]
    for (acc in all_accs) {
      secondary_to_primary[acc] <- primary
    }
  }

  # ── Build numeric matrix ──
  mat <- as.matrix(dt[, sample_names, drop = FALSE])
  rownames(mat) <- dt$Accession_primary
  class(mat)    <- "numeric"

  n_prots   <- nrow(mat)
  n_samples <- ncol(mat)
  if (n_prots == 0)
    return(list(error = "No protein rows found after parsing."))
  if (n_samples < 2)
    return(list(error = "Fewer than 2 samples detected."))

  # log2 or linear? The two scales do not overlap: log2 protein intensities stay below ~40,
  # linear ones run into the thousands and above. log2 input is returned as
  # 2^x, so dia_limma_summarize() and dia_normalization_check() log2 it back.
  .finite <- mat[is.finite(mat)]
  if (length(.finite) == 0L)
    return(list(error = "No numeric quantity values found."))
  input_log2 <- max(.finite) < 100
  if (input_log2) {
    mat <- 2^mat
    logger("[MQ] Quantities are log2 (Perseus); converted to linear scale for OSMSS's own log2 step.")
  } else if (any(.finite < 0)) {
    return(list(error = paste0("The quantity columns hold negative values on a linear scale ",
                               "(maximum ", signif(max(.finite), 3), "). Expected Perseus log2 ",
                               "values or linear Spectronaut PG.Quantity values.")))
  }

  pct_na <- mean(is.na(mat))
  if (pct_na > 0.80)
    logger(paste0("[MQ] Warning: ", round(pct_na * 100, 1),
                  "% of intensity values are missing. Check filtering settings."))

  logger(paste0("[MQ] Global parse complete: ", n_prots, " proteins x ", n_samples, " samples."))

  list(
    protein_matrix       = mat,
    gene_names           = acc_to_gene,
    accessions           = rownames(mat),
    sample_names         = sample_names,
    n_proteins           = n_prots,
    n_samples            = n_samples,
    multi_acc_map        = multi_acc_map,
    secondary_to_primary = secondary_to_primary,
    input_log2           = input_log2
  )
}

# -----------------------------------------------------------------------------
# parse_diann_pg_matrix()
# Parse DIA-NN's protein-group matrix (report.pg_matrix.tsv): one row per
# protein group, one column per run, holding DIA-NN's protein quantities
# (MaxLFQ in 1.x, QuantUMS in 2.x; cross-run normalized unless --no-norm),
# which DIA-NN filters at 1% protein-group FDR. Only Protein.Group is required:
# deposited matrices are often trimmed to Protein.Group plus the runs (e.g.
# PXD071300). The annotation columns DIA-NN writes are named in
# DIANN_PG_ANNOTATION_COLS; every other column is a run and must be numeric.
# Run columns are named as DIA-NN was given the runs, often full paths; sample
# names drop the folder and the raw-file extension. Protein groups without a
# gene name are labelled by their leading accession (as on the TMT route).
# Returns the same structure as parse_spectronaut_pivot(), plus n_no_gene.
# -----------------------------------------------------------------------------
DIANN_PG_ANNOTATION_COLS <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                              "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")

# Run file name as a sample name: folder and raw-file extension removed. FragPipe
# hands DIA-NN "<run>_uncalibrated.mzML" for Thermo .raw input (FragPipe-GUI
# CmdDiann.java), so that suffix is removed as well.
diann_run_name <- function(x) {
  x <- sub("_uncalibrated\\.mzml$", "", sub("^.*[\\\\/]", "", x), ignore.case = TRUE)
  sub("\\.(raw|mzml|d|wiff2?|dia)$", "", x, ignore.case = TRUE)
}

parse_diann_pg_matrix <- function(file_path, logger = message) {
  logger(paste0("[DIA-NN] Parsing protein-group matrix: ", basename(file_path)))
  all_cols  <- strsplit(sub("\r$", "", readLines(file_path, n = 1L, warn = FALSE)), "\t")[[1]]
  run_files <- setdiff(all_cols, DIANN_PG_ANNOTATION_COLS)
  if (length(run_files) < 2L)
    return(list(error = "Fewer than 2 run columns found next to Protein.Group."))
  sample_names <- diann_run_name(run_files)
  if (anyDuplicated(sample_names))
    return(list(error = paste0("Two run columns give the same sample name after removing folder and ",
                               "extension: ", paste(unique(sample_names[duplicated(sample_names)]),
                                                    collapse = ", "), ".")))

  dt <- tryCatch(
    data.table::fread(file = file_path, sep = "\t", select = intersect(c("Protein.Group", "Genes", run_files), all_cols),
                      na.strings = c("NA", "NaN", ""), data.table = FALSE, showProgress = FALSE),
    error = function(e) list(error = paste0("File read error: ", conditionMessage(e))))
  if (!is.null(dt$error)) return(dt)
  not_num <- run_files[!vapply(dt[run_files], is.numeric, logical(1))]
  if (length(not_num) > 0L)
    return(list(error = paste0("These columns are neither DIA-NN annotation columns nor numeric run ",
                               "quantities: ", paste(utils::head(not_num, 5L), collapse = ", "), ".")))

  dt$Accession_primary <- sub(";.*", "", dt$Protein.Group)
  genes <- if ("Genes" %in% names(dt)) sub(";.*", "", dt$Genes) else rep(NA_character_, nrow(dt))
  no_gene <- is.na(genes) | !nzchar(trimws(genes))
  genes[no_gene] <- dt$Accession_primary[no_gene]
  if (any(no_gene))
    logger(paste0("[DIA-NN] ", sum(no_gene), " protein groups have no gene name; labelled by accession."))
  dt$Gene_primary <- genes
  dt$no_gene      <- no_gene

  # Groups sharing a leading accession: keep the most completely quantified
  # (as parse_spectronaut_pivot() does)
  if (anyDuplicated(dt$Accession_primary) > 0L)
    dt <- dt[order(-rowSums(!is.na(dt[, run_files, drop = FALSE]))), , drop = FALSE]
  dup_mask      <- duplicated(dt$Accession_primary)
  n_dup_removed <- sum(dup_mask)
  if (n_dup_removed > 0L) {
    logger(paste0("[DIA-NN] Warning: ", n_dup_removed,
                  " protein groups sharing a leading accession were removed (keeping the most completely quantified)."))
    dt <- dt[!dup_mask, , drop = FALSE]
  }

  mat <- as.matrix(dt[, run_files, drop = FALSE])
  dimnames(mat) <- list(dt$Accession_primary, sample_names)
  if (nrow(mat) == 0L)
    return(list(error = "No protein rows found after parsing."))
  .finite <- mat[is.finite(mat)]
  if (length(.finite) > 0L && max(.finite) < 100)
    return(list(error = paste0("The quantities look log2-transformed (maximum ", signif(max(.finite), 3),
                               "); DIA-NN writes linear quantities. Upload the matrix as DIA-NN saved it.")))

  logger(paste0("[DIA-NN] Parse complete: ", nrow(mat), " proteins x ", ncol(mat), " samples."))
  list(
    protein_matrix = mat,
    gene_names     = setNames(dt$Gene_primary, dt$Accession_primary),
    accessions     = rownames(mat),
    sample_names   = sample_names,
    n_proteins     = nrow(mat),
    n_samples      = ncol(mat),
    n_dup_removed  = n_dup_removed,
    multi_acc_map  = setNames(dt$Protein.Group, dt$Accession_primary),
    n_no_gene      = sum(dt$no_gene)
  )
}

# -----------------------------------------------------------------------------
# diann_read_log()
# DIA-NN's report.log.txt. Line 1 is "DIA-NN <version> (Data-Independent
# Acquisition by Neural Networks)"; the command line follows ("diann.exe --f
# <run> --f <run> ... --qvalue 0.01 ..."; DIA-NN 1.8.1 logs in Benchmark/opDEA-
# DDA_DIA_TMT Single/Raw ProteomeXchange/10482353/DIA). Cross-run normalization
# is DIA-NN's default; --no-norm switches it off, --global-norm / --no-rt-norm
# and --sig-norm change its kind. run_files are the runs DIA-NN loaded ("Loading
# run <path>" lines; the --f arguments when there are none), matched against
# the matrix's run columns by the caller. FragPipe passes the runs in a file
# list (--cfg filelist_diann.txt) and writes to dia-quant-output, which marks a
# DIA-NN run within FragPipe.
# Returns list(version, normalization, norm_mode, run_files, run_within, source)
# or list(error); diann_read_fp_workflow() returns the same fields.
# -----------------------------------------------------------------------------
diann_norm_mode <- function(flags) {
  if ("no-norm" %in% flags) "none"
  else if (any(c("global-norm", "no-rt-norm") %in% flags)) "global"
  else if ("sig-norm" %in% flags) "signal-dependent"
  else "RT-dependent"
}

diann_read_log <- function(path) {
  x <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  x <- sub("\r$", "", x)
  if (length(x) == 0L || !grepl("^DIA-NN .*Data-Independent Acquisition by Neural Networks", x[1]))
    return(list(error = paste0("This is neither DIA-NN\u2019s report.log.txt (first line \u2018DIA-NN <version> ",
                               "(Data-Independent Acquisition by Neural Networks)\u2019) nor FragPipe\u2019s ",
                               "fragpipe.workflow.")))
  cmd <- x[grepl("\\s--(f|dir|cfg|fasta|lib|out)\\s", x)][1]
  if (is.na(cmd))
    return(list(error = "The DIA-NN command line was not found in this log."))
  args   <- trimws(strsplit(paste0(" ", cmd), "\\s+--")[[1]][-1])
  flags  <- sub("\\s.*$", "", args)
  loaded <- unique(trimws(sub("^.*Loading run ", "", grep("Loading run ", x, value = TRUE, fixed = TRUE))))
  list(
    version       = sub("^DIA-NN ([^ ]+).*$", "\\1", x[1]),
    normalization = !"no-norm" %in% flags,
    norm_mode     = diann_norm_mode(flags),
    run_files     = if (length(loaded)) loaded else trimws(sub("^f\\s+", "", args[flags == "f"])),
    run_within    = if (grepl("filelist_diann\\.txt|dia-quant-output", cmd)) "FragPipe" else NA_character_,
    source        = "report.log.txt"
  )
}

# fragpipe.workflow of a FragPipe DIA run (DIA-NN quantification, Quant (DIA)
# tab). FragPipe adds no normalization flag of its own; the user's additional
# DIA-NN options (diann.cmd-opts) are passed through (FragPipe-GUI
# CmdDiann.java), so the quantities are normalized unless those hold --no-norm.
# The workflow does not list the runs (run_files is empty).
diann_read_fp_workflow <- function(path) {
  x    <- sub("\r$", "", tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0)))
  vals <- fp_workflow_values(x)
  if (!any(startsWith(names(vals), "diann.")))
    return(list(error = "This FragPipe workflow has no DIA-NN settings."))
  if (identical(unname(tolower(vals["diann.run-dia-nn"])), "false"))
    return(list(error = "DIA-NN quantification was not run in this FragPipe workflow."))
  opts  <- if ("diann.cmd-opts" %in% names(vals)) vals[["diann.cmd-opts"]] else ""
  flags <- sub("\\s.*$", "", trimws(strsplit(paste0(" ", opts), "\\s+--")[[1]][-1]))
  ver   <- fp_workflow_version(x)
  list(
    version       = NA_character_,
    normalization = !"no-norm" %in% flags,
    norm_mode     = diann_norm_mode(flags),
    run_files     = character(0),
    run_within    = if (is.na(ver)) "FragPipe" else paste("FragPipe", ver),
    source        = "fragpipe.workflow"
  )
}

# -----------------------------------------------------------------------------
# dia_normalization_check()
# Does an uploaded protein matrix look cross-run normalized already?
# Spectronaut applies cross-run normalization by default (Quantification >
# Cross Run Normalization; strategy "Automatic" = local regression for < 500
# runs), so its PG.Quantity exports are normally normalized. A second,
# column-wise median centering in OSMSS is then redundant, and in comparisons
# where a large share of the proteome differs between conditions it shifts the
# unchanged background (DIA benchmark, 2026-10-04: offsets up to 0.15 log2,
# FDP 0.56 vs 0.39 without it).
#
# The check: per-run median of log2 intensities over proteins quantified in
# every assigned run, centred within each condition. Normalization removes
# run-to-run (loading / spray) differences between replicates of a condition,
# so the spread of these within-condition offsets is small. It does not judge
# between-condition offsets, which can be real composition differences.
# Calibration (6 Spectronaut exports, normalization on: SD 0.011-0.019; the
# same runs re-quantified without normalization: SD 0.026-0.267). Checked on
# 7 DIA-NN 1.8.1 reports (2026-10-14): normalized PG.MaxLFQ, as in
# report.pg_matrix.tsv, SD 0.004-0.022 (2 of 7 borderline); not normalized
# PG.Quantity SD 0.030-0.207 (Benchmark/opDEA- DDA_DIA_TMT Single/DIA/diann_norm_check).
#   SD <= 0.02          consistent with normalized input
#   0.02 < SD <= 0.03    borderline (well-loaded runs can look normalized)
#   SD > 0.03            suggests the input was not normalized
# Returns list(verdict, within_sd, n_complete, n_needed, n_runs_used,
#              between_offsets, message). verdict is one of "normalized",
#              "borderline", "not_normalized", "inconclusive".
# unit / normalized_label only change the wording of the message (TMT passes
# "channel" / "cross-channel normalized"); the defaults give the DIA wording.
# -----------------------------------------------------------------------------
dia_normalization_check <- function(protein_matrix, cond_groups,
                                    sd_normalized = 0.02, sd_not_normalized = 0.03,
                                    min_complete_n = 100L, min_complete_frac = 0.05,
                                    feature_label = "proteins",
                                    unit = "run", normalized_label = "cross-run normalized") {
  .unit_cap <- paste0(toupper(substr(unit, 1L, 1L)), substr(unit, 2L, nchar(unit)))
  out <- list(verdict = "inconclusive", within_sd = NA_real_, n_complete = 0L,
              n_needed = NA_integer_, n_runs_used = 0L, between_offsets = NULL,
              message = "")
  if (is.null(protein_matrix) || length(cond_groups) < 1L) {
    out$message <- "No data or no condition groups yet."
    return(out)
  }
  runs <- intersect(unname(unlist(cond_groups)), colnames(protein_matrix))
  cond <- setNames(rep(names(cond_groups), vapply(cond_groups, length, integer(1))),
                   unlist(cond_groups))[runs]
  rep_groups <- names(which(table(cond) >= 2L))
  runs <- runs[cond %in% rep_groups]; cond <- cond[runs]
  out$n_runs_used <- length(runs)
  if (length(runs) < 2L) {
    out$message <- paste0("No condition has at least 2 ", unit, "s, so replicate ", unit, "s cannot be compared.")
    return(out)
  }
  x <- protein_matrix[, runs, drop = FALSE]
  x[!is.finite(x) | x <= 0] <- NA
  x <- log2(x)
  cc <- stats::complete.cases(x)
  out$n_complete <- sum(cc)
  out$n_needed   <- max(as.integer(min_complete_n), ceiling(min_complete_frac * nrow(x)))
  if (out$n_complete < out$n_needed) {
    out$message <- paste0("Only ", out$n_complete, " ", feature_label, " are quantified in every ", unit, " (",
                          out$n_needed, " needed), so the check is inconclusive.")
    return(out)
  }
  med <- apply(x[cc, , drop = FALSE], 2, median)
  off <- med - median(med)
  resid <- unlist(tapply(off, cond, function(v) v - mean(v)))
  out$within_sd <- stats::sd(resid)
  out$between_offsets <- tapply(off, cond, mean)
  out$verdict <- if (out$within_sd <= sd_normalized) "normalized"
                 else if (out$within_sd <= sd_not_normalized) "borderline"
                 else "not_normalized"
  out$message <- paste0(
    .unit_cap, " medians of the ", out$n_complete, " ", feature_label, " quantified in every ", unit, " differ between ",
    "replicate ", unit, "s of the same condition by SD ", formatC(out$within_sd, format = "f", digits = 3),
    " log2. ",
    switch(out$verdict,
      normalized     = paste0("This is consistent with ", normalized_label, " input."),
      borderline     = paste0("This is borderline: well-loaded ", unit, "s can look normalized without being normalized."),
      not_normalized = paste0("This suggests the input was NOT ", normalized_label, ".")))
  out
}

# -----------------------------------------------------------------------------
# Label-free normalization: OSMSS normalizes (global median) unless the input
# was normalized upstream. PD and MaxQuant: never - OSMSS reads their raw
# columns ("Abundance", "Intensity"), not the normalized ones they also export.
# FragPipe: as recorded in the required fragpipe.workflow (ionquant.normalization);
# IonQuant normalizes by default, and normalizing its intensities a second time
# raised false positives in the OpDEA DDA benchmark (Benchmark/opDEA- DDA_DIA_TMT
# Single/DDA Global Median/dda_norm_check/fragpipe_norm, 2026-10-13).
# Skipping applies on the ProDA global route only: the invariant-protein panel
# (IP route) is always applied, and SAINT normalizes internally.
lf_normalization_skippable <- function(de_engine, ip_sub)
  identical(de_engine %||% "proda", "proda") && !identical(ip_sub, "invariant")

# Effective setting for a summarization: FALSE only when the input was
# normalized upstream (rv$lf_apply_normalization FALSE) and skipping applies.
lf_apply_normalization_effective <- function(rv)
  !(lf_normalization_skippable(rv$lf_de_engine, rv$lf_ip_sub) && isFALSE(rv$lf_apply_normalization))

# -----------------------------------------------------------------------------
# FragPipe label-free input
# combined_protein.tsv is converted to Proteome Discoverer protein layout and
# taken through the PD protein-direct route (build_proda_protein_matrix,
# Route A). Checked on the 12 OpDEA DDA datasets (FragPipe 20.0): the values
# read this way equal OpDEA's FragPipe protein tables exactly
# (dda_norm_check/fragpipe_norm/R/04_real_fragpipe_files.R).

# fragpipe.workflow: "key=value" lines; the first line names the version,
# e.g. "# FragPipe (20.0) runtime properties". fp_workflow_values() returns the
# settings as a named character vector, fp_workflow_version() the version.
fp_workflow_values <- function(x) {
  kv <- x[grepl("=", x, fixed = TRUE) & !startsWith(trimws(x), "#")]
  setNames(trimws(sub("^[^=]*=", "", kv)), trimws(sub("=.*$", "", kv)))
}
fp_workflow_version <- function(x) {
  ver <- regmatches(x[1], regexpr("(?<=FragPipe \\()[^)]+", x[1], perl = TRUE))
  if (length(ver)) ver else NA_character_
}

# The IonQuant settings OSMSS uses, or list(error = ...) when the file cannot be used.
fp_read_workflow <- function(path) {
  x <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  if (length(x) == 0L) return(list(error = "The workflow file could not be read."))
  vals <- fp_workflow_values(x)
  g    <- function(k) if (k %in% names(vals)) vals[[k]] else NA_character_
  if (!any(startsWith(names(vals), "ionquant.")))
    return(list(error = "This is not a FragPipe workflow file: it has no IonQuant settings."))
  if (identical(tolower(g("ionquant.run-ionquant")), "false"))
    return(list(error = paste0("IonQuant was not run in this FragPipe workflow, so combined_protein.tsv has no ",
                               "MS1 intensities to analyse.")))
  if (identical(tolower(g("ionquant.use-labeling")), "true"))
    return(list(error = "This FragPipe workflow used isotope labeling; the label-free upload needs a label-free run."))
  norm <- g("ionquant.normalization")
  if (!norm %in% c("0", "1"))
    return(list(error = paste0("The workflow file does not state whether IonQuant normalized the intensities ",
                               "(no ionquant.normalization setting).")))
  list(version       = fp_workflow_version(x),
       normalization = identical(norm, "1"),
       mbr           = identical(g("ionquant.mbr"), "1"),
       top_n         = suppressWarnings(as.integer(g("ionquant.tp"))),
       uniqueness    = switch(g("ionquant.uniqueness") %||% "", "0" = "unique + razor", "1" = "unique only",
                              "2" = "all", NA_character_))
}

# combined_protein.tsv -> PD-layout CSV lines. Uses the per-sample
# "<sample> Intensity" columns (IonQuant), not "<sample> MaxLFQ Intensity".
# Columns become "Abundance: <sample>" (step 4 lists them); Accession is
# "Protein ID", Gene Symbol is "Gene". Decoy rows (rev_) are dropped.
fp_combined_protein_to_pd <- function(path) {
  cp <- tryCatch(data.table::fread(path, sep = "\t", quote = "", data.table = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(cp) || !all(c("Protein", "Protein ID") %in% colnames(cp)))
    return(list(error = "combined_protein.tsv could not be read (no 'Protein ID' column)."))
  ic <- grep(" Intensity$", colnames(cp), value = TRUE)
  ic <- ic[!grepl(" MaxLFQ Intensity$", ic)]
  if (length(ic) < 2L)
    return(list(error = "combined_protein.tsv has fewer than 2 '<sample> Intensity' columns. Was IonQuant run?"))
  dec <- grepl("^rev_", cp$Protein)
  cp  <- cp[!dec, , drop = FALSE]
  acc <- as.character(cp[["Protein ID"]])
  if (anyDuplicated(acc)) acc <- as.character(cp$Protein)
  samples <- sub(" Intensity$", "", ic)
  ab <- paste0("Abundance: ", samples)
  pd <- data.frame(Accession = acc,
                   Description = if ("Description" %in% colnames(cp)) as.character(cp$Description) else "",
                   `Gene Symbol` = if ("Gene" %in% colnames(cp)) as.character(cp$Gene) else "",
                   check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_along(ic)) pd[[ab[i]]] <- suppressWarnings(as.numeric(cp[[ic[i]]]))
  tc <- textConnection("pd_lines", "w", local = TRUE)
  utils::write.csv(pd, tc, row.names = FALSE, na = ""); close(tc)
  list(lines = pd_lines, protein_idx = seq_len(nrow(pd)) + 1L, headers = colnames(pd),
       samples = samples, n_proteins = nrow(pd), n_decoys_removed = sum(dec),
       has_gene = "Gene" %in% colnames(cp))
}

# Methods sentence for FragPipe input (provenance$fragpipe_settings = fp_read_workflow()).
fp_methods_sentence <- function(fs) {
  if (is.null(fs)) return(paste0("Protein intensities were taken from the ‘Intensity’ columns of FragPipe’s ",
                                 "combined_protein.tsv (IonQuant); the MaxLFQ intensity columns were not used."))
  ions <- if (isTRUE(is.finite(fs$top_n)) && fs$top_n > 0L) paste0("the top ", fs$top_n, " ions") else "all ions"
  paste0("Protein intensities were taken from the ‘Intensity’ columns of FragPipe’s combined_protein.tsv",
         if (!is.na(fs$version %||% NA)) paste0(" (FragPipe ", fs$version, ")") else "",
         ": IonQuant MS1 quantification summed over ", ions, " of each protein",
         if (!is.na(fs$uniqueness %||% NA)) paste0(", using ", fs$uniqueness, " peptides") else "",
         ", with match-between-runs ", if (isTRUE(fs$mbr)) "on" else "off",
         ". The MaxLFQ intensity columns were not used. Zero values were set to missing and intensities were ",
         "log2-transformed.")
}

# -----------------------------------------------------------------------------
# dia_limma_summarize()
# Build a ProteinLevelData (PLD) shim from DIA parsed data.
# Applies: replace 0->NA, log2 transform, missingness filter (>=2 non-NA),
# column-wise median centering (unless apply_centering = FALSE: the input was
# already cross-run normalized upstream, e.g. by Spectronaut; the centering
# constants are still computed and reported, but not applied). The default
# (TRUE) keeps every existing caller, including the DIA PTM route, unchanged.
# Returns: list with processed PLD, protein_matrix, design_df, acc_to_gene,
#          pipeline_label, n_prots_input, n_prots_processed, software, normalization
# -----------------------------------------------------------------------------
dia_limma_summarize <- function(parsed_data, cond_groups,
                                 subject_pairing = NULL,
                                 factorial_design = NULL,
                                 norm_min_complete_n    = 100L,
                                 norm_min_complete_frac = 0.05,
                                 apply_centering        = TRUE,
                                 logger = message) {
  logger(paste0("[DIA] Starting summarization (log2",
                if (isTRUE(apply_centering)) " + median centering" else
                  "; no median centering: input already cross-run normalized", ")..."))

  mat           <- parsed_data$protein_matrix
  acc_to_gene   <- parsed_data$gene_names
  n_prots_input <- nrow(mat)

  # 1. Replace exact zeros with NA
  mat[mat == 0] <- NA

  # 2. log2 transform
  mat <- log2(mat)
  mat[is.infinite(mat) | is.nan(mat)] <- NA

  # 3. Missingness filter: keep proteins with ≥50% valid values in at least one
  #    condition group (field standard; matches Perseus/MaxQuant default behaviour).
  #    This is group-aware  -" a protein absent from one condition but fully detected
  #    in another is retained, as it may be a genuine exclusive hit.
  keep_global <- apply(mat, 1, function(row) {
    any(vapply(cond_groups, function(g_samples) {
      s  <- intersect(g_samples, colnames(mat))
      ns <- length(s)
      ns > 0 && (sum(!is.na(row[s])) / ns) >= 0.5
    }, logical(1)))
  })
  mat <- mat[keep_global, , drop = FALSE]
  logger(paste0("[DIA] Removed ", n_prots_input - nrow(mat),
                " proteins failing the \u226550% valid-in-\u22651-group filter (",
                n_prots_input, " \u2192 ", nrow(mat), ")."))

  # 4. Column-wise median centering
  # Constants are derived from COMPLETE CASES where enough exist. Computing them
  # with na.rm = TRUE over each sample's own detected set makes the correction
  # biased in the direction of the missingness: a high-dropout sample takes its
  # median from a higher-abundance subset and is over-corrected.
  #
  # The sufficiency rule is BOTH absolute and relative, because this pipeline
  # must serve targeted panels (few hundred proteins) and full proteomes alike:
  # a fixed count of 200 is unreachable for the former and potentially
  # unrepresentative for the latter.
  #
  # Both bases are always computed and their discrepancy recorded, so that every
  # run self-reports whether the choice of basis mattered on THAT dataset rather
  # than relying on it having been anticipated during development.
  .cc     <- stats::complete.cases(mat)
  .n_cc   <- sum(.cc)
  .n_need <- max(as.integer(norm_min_complete_n),
                 ceiling(norm_min_complete_frac * nrow(mat)))

  .med_all <- apply(mat, 2, median, na.rm = TRUE)
  .med_cc  <- if (.n_cc > 0L) apply(mat[.cc, , drop = FALSE], 2, median)
              else rep(NA_real_, ncol(mat))

  .use_cc <- .n_cc >= .n_need && all(is.finite(.med_cc))

  if (.use_cc) {
    col_medians <- .med_cc
    .norm_basis <- paste0("complete cases (n=", .n_cc, " of ", nrow(mat),
                          "; threshold ", .n_need, ")")
  } else {
    col_medians <- .med_all
    .norm_basis <- paste0("all observed values (", .n_cc, " complete cases of ",
                          nrow(mat), "; threshold ", .n_need, " not met)")
  }

  # Discrepancy between the two bases, expressed on the constants themselves
  # (log2 units). This is the number that tells a user whether complete-case
  # centering changed their result. Because centering constants enter a contrast
  # as a per-column offset, the effect on any log2FC is bounded by the spread of
  # this difference across the contrasted columns.
  .norm_delta <- if (all(is.finite(.med_cc)))
                   max(abs((.med_cc - median(.med_cc, na.rm = TRUE)) -
                           (.med_all - median(.med_all, na.rm = TRUE))))
                 else NA_real_

  global_median <- median(col_medians, na.rm = TRUE)
  # Largest per-sample shift the centering applies (or would apply), log2 units
  .center_max_shift <- max(abs(col_medians - global_median), na.rm = TRUE)
  if (isTRUE(apply_centering)) {
    mat <- sweep(mat, 2, col_medians - global_median)
    logger(paste0("[DIA] Median centering constants computed on ", .norm_basis, "."))
  } else {
    logger(paste0("[DIA] Median centering NOT applied (input already cross-run normalized). ",
                  "Constants on ", .norm_basis, " would have shifted samples by up to ",
                  signif(.center_max_shift, 3), " log2."))
  }
  if (isTRUE(apply_centering) && is.finite(.norm_delta))
    logger(paste0("[DIA] Centering-basis discrepancy: max |complete-case - ",
                  "all-observed| constant = ", signif(.norm_delta, 3),
                  " log2 units. Fold changes differ from the alternative basis ",
                  "by at most this amount."))

  n_prots_processed <- nrow(mat)

  # 5. Build design_df
  all_samples <- unname(unlist(cond_groups))
  cond_map    <- setNames(
    rep(names(cond_groups), times = vapply(cond_groups, length, integer(1))),
    unlist(cond_groups)
  )

  design_df <- data.frame(
    sample_name  = all_samples,
    condition    = cond_map[all_samples],
    stringsAsFactors = FALSE
  )

  if (!is.null(subject_pairing)) {
    design_df$bioreplicate <- subject_pairing[all_samples]
  } else {
    design_df$bioreplicate <- ave(
      seq_along(all_samples),
      cond_map[all_samples],
      FUN = function(i) paste0(cond_map[all_samples][i[1]], "_", seq_along(i))
    )
  }

  # 6. Build ProteinLevelData (PLD) shim compatible with LF DDA downstream patterns
  pld_rows <- lapply(seq_len(nrow(mat)), function(i) {
    prot  <- rownames(mat)[i]
    vals  <- mat[i, ]
    runs  <- names(vals)
    valid <- !is.na(vals)
    data.frame(
      Protein        = prot,
      LogIntensities = vals[valid],
      RUN            = runs[valid],
      GROUP          = cond_map[runs[valid]],
      SUBJECT        = if (!is.null(subject_pairing)) {
        subject_pairing[runs[valid]]
      } else {
        design_df$bioreplicate[match(runs[valid], design_df$sample_name)]
      },
      stringsAsFactors = FALSE
    )
  })
  pld <- do.call(rbind, pld_rows)

  logger(paste0("[DIA] Summarization complete: ", n_prots_processed,
                " proteins retained (from ", n_prots_input, " input)."))

  # Subset matrix to only assigned samples
  valid_samples <- all_samples[all_samples %in% colnames(mat)]
  mat_out       <- mat[, valid_samples, drop = FALSE]

  list(
    processed          = list(ProteinLevelData = pld),
    protein_matrix     = mat_out,
    design_df          = design_df,
    acc_to_gene        = acc_to_gene,
    pipeline_label     = if (isTRUE(apply_centering)) "DIA (limma, median-centered)"
                         else "DIA (limma, input cross-run normalized upstream)",
    n_prots_input      = n_prots_input,
    n_prots_processed  = n_prots_processed,
    software           = "spectronaut",
    normalization      = if (isTRUE(apply_centering)) "median_centering" else "none_upstream_normalized",
    normalization_applied      = isTRUE(apply_centering),
    normalization_basis        = .norm_basis,
    normalization_used_cc      = .use_cc,
    normalization_n_complete   = .n_cc,
    normalization_delta_log2   = .norm_delta,
    normalization_max_shift_log2 = .center_max_shift
  )
}

# -----------------------------------------------------------------------------
# dia_fit_scope_guard()
# Decides the model scope for dia_limma_de(fit_all_groups = "auto"): fit all
# conditions in one model (pooled residual variance, more degrees of freedom)
# unless the conditions differ clearly in typical within-condition variability.
# Per condition with >= min_n samples: median, over proteins with >= 2 values in
# that condition, of the within-condition SD (log2). ratio = max / min over those
# conditions. ratio <= threshold (1.5) -> all conditions; otherwise only the two
# compared conditions are fitted. Fewer than 2 conditions with >= min_n samples:
# the ratio cannot be assessed and the two compared conditions are fitted.
# Two-condition designs give the same fit either way.
# Returns list(all_groups, ratio, threshold, n_conditions, n_assessed,
#              median_sd (named), reason, message).
# -----------------------------------------------------------------------------
dia_fit_scope_guard <- function(mat, design_df, threshold = 1.5, min_n = 3L) {
  ok    <- !is.na(design_df$condition) & nzchar(design_df$condition) & design_df$sample_name %in% colnames(mat)
  dd    <- design_df[ok, , drop = FALSE]
  conds <- unique(dd$condition)
  sizes <- table(dd$condition)[conds]
  big   <- conds[sizes >= min_n]
  med   <- vapply(big, function(cn) {
    x  <- mat[, dd$sample_name[dd$condition == cn], drop = FALSE]
    nv <- rowSums(!is.na(x))
    s  <- apply(x[nv >= 2L, , drop = FALSE], 1, stats::sd, na.rm = TRUE)
    stats::median(s, na.rm = TRUE)
  }, numeric(1))
  ratio <- if (length(med) >= 2L && all(is.finite(med)) && min(med) > 0) max(med) / min(med) else NA_real_
  out <- list(threshold = threshold, n_conditions = length(conds), n_assessed = length(big),
              median_sd = med, ratio = ratio)
  if (length(conds) <= 2L) {
    out$all_groups <- TRUE; out$reason <- "two_conditions"
    out$message <- "two conditions in the design, so all-conditions and two-condition fits coincide."
  } else if (is.na(ratio)) {
    out$all_groups <- FALSE; out$reason <- "not_assessable"
    out$message <- paste0("within-condition variability could not be compared (fewer than 2 conditions with >= ",
                          min_n, " samples); only the two compared conditions are fitted.")
  } else if (ratio <= threshold) {
    out$all_groups <- TRUE; out$reason <- "homogeneous"
    out$message <- paste0("all ", length(conds), " conditions fitted in one model (typical within-condition SD ",
                          "differs between conditions by a factor of ", formatC(ratio, format = "f", digits = 2),
                          ", within ", threshold, ").")
  } else {
    out$all_groups <- FALSE; out$reason <- "heterogeneous"
    out$message <- paste0("only the two compared conditions fitted (typical within-condition SD differs between ",
                          "conditions by a factor of ", formatC(ratio, format = "f", digits = 2), ", above ", threshold, ").")
  }
  out
}

# -----------------------------------------------------------------------------
# dia_limma_de()
# Differential expression for DIA protein-level data using limma.
# eBayes(trend=TRUE, robust=TRUE). Supports pairing (duplicateCorrelation),
# batch in model, and factorial design.
# Returns: list(result_df, fit_object, contrast, n_tested, n_imputed, comparison)
# -----------------------------------------------------------------------------
dia_limma_de <- function(cache, cond_test, cond_ref,
                          min_valid_frac   = 0.5,
                          require_both     = TRUE,
                          fit_all_groups   = FALSE,
                          impute_max_miss_per_cond = NULL,   # DEPRECATED - see max_missing_per_cond
                          max_missing_per_cond     = 1L,
                          knn_k                    = 10L,
                          batch_groups     = NULL,
                          batch_for_model  = FALSE,
                          subject_pairing  = NULL,
                          factorial_design = NULL,
                          perseus_lod      = FALSE,
                          observed_mat     = NULL,
                          use_array_weights = FALSE,
                          logger           = message) {
  if (is.null(cache) || is.null(cache$protein_matrix))
    return(list(error = "No DIA summarization cache found. Please run the wizard first."))

  mat         <- cache$protein_matrix
  design_df   <- cache$design_df
  acc_to_gene <- cache$acc_to_gene
  # Degradation guard. A NULL or empty map made `genes <- acc_to_gene[prots]`
  # return a zero-length vector, and the data.frame() built from it failed with
  # "arguments imply differing number of rows" - an opaque error in place of the
  # intended accession fallback. An empty NAMED character vector subsets to all
  # NA instead, which the existing genes[is.na(genes)] <- prots[is.na(genes)]
  # line then fills with accessions. One guard here covers both exits, since
  # both read this same variable. No effect when the map is populated.
  if (is.null(acc_to_gene) || length(acc_to_gene) == 0L) {
    logger(paste0("[DIA] No accession-to-gene map supplied - gene names will ",
                  "fall back to accessions for every protein."))
    acc_to_gene <- setNames(character(0), character(0))
  }

  # Determine which samples to include.
  #
  # fit_all_groups = TRUE fits ONE model across every condition in the design and
  # extracts only the requested contrast. This substantially increases residual df
  # and gives every comparison one consistent empirical Bayes prior.
  #
  # CURRENT STATE - read before changing the default:
  #   - TMT (compute_tmt_de) passes TRUE. Enabled and validated.
  #   - DIA passes FALSE (the default). This is DELIBERATE, not an oversight.
  #     DIA routes proteins failing the comparison filter to dia_nmar_fisher()
  #     for presence/absence analysis. Fitting across all conditions changes
  #     which contrasts are estimable and therefore which proteins are routed
  #     where, so enabling it for DIA requires re-validating that interaction.
  #     It also strengthens the common-residual-variance assumption across a
  #     larger and more heterogeneous set of conditions.
  #
  # Do not flip the default without scoping that work separately.
  #
  # fit_all_groups = "auto" (passed by the DIA DE screen since 2026-10-04): fit all
  # conditions unless dia_fit_scope_guard() finds that the conditions differ clearly
  # in typical within-condition variability (max / min ratio > 1.5), in which case
  # only the two compared conditions are fitted. The scoped interactions above were
  # handled as follows: the test filter (and therefore NMAR routing) is evaluated on
  # the two compared conditions in either scope, and k-NN imputation fills only the
  # compared conditions' columns, so other conditions enter the fit as observed.
  fit_scope_guard <- NULL
  if (identical(fit_all_groups, "auto")) {
    fit_scope_guard <- dia_fit_scope_guard(mat, design_df)
    fit_all_groups  <- isTRUE(fit_scope_guard$all_groups)
    logger(paste0("[DIA] Model scope: ", fit_scope_guard$message))
  }
  .is_factorial <- !is.null(factorial_design) && isTRUE(factorial_design$enabled)
  if (!is.null(factorial_design) && isTRUE(factorial_design$enabled)) {
    use_samples <- colnames(mat)
  } else if (isTRUE(fit_all_groups)) {
    use_samples <- design_df$sample_name[!is.na(design_df$condition) &
                                           nzchar(design_df$condition)]
    use_samples <- use_samples[use_samples %in% colnames(mat)]
    if (!all(c(cond_test, cond_ref) %in% design_df$condition[
          design_df$sample_name %in% use_samples]))
      return(list(error = paste0(
        "fit_all_groups: '", cond_test, "' or '", cond_ref,
        "' not present among the design conditions.")))
    logger(paste0("[DIA] Fitting across all ", length(unique(
      design_df$condition[design_df$sample_name %in% use_samples])),
      " conditions (", length(use_samples), " samples); extracting ",
      cond_test, " vs ", cond_ref, "."))
  } else {
    groups_needed <- c(cond_test, cond_ref)
    use_samples   <- design_df$sample_name[design_df$condition %in% groups_needed]
    use_samples   <- use_samples[use_samples %in% colnames(mat)]
  }

  sub_mat    <- mat[, use_samples, drop = FALSE]
  sub_design <- design_df[design_df$sample_name %in% use_samples, , drop = FALSE]

# 1. Missingness filter
  #    Perseus LOD mode : ≥1 valid value anywhere (filter is applied AFTER imputation
  #                       in Perseus / paper workflow, so effectively no filter)
  #    Standard mode    : ≥50% valid in BOTH compared conditions.
  #                       Proteins solid in one condition but absent/sparse in the
  #                       other are deliberately NOT tested here — they fall through
  #                       to dia_nmar_fisher() for presence/absence analysis rather
  #                       than having an entire arm imputed.
  group_vec <- sub_design$condition[match(colnames(sub_mat), sub_design$sample_name)]

  keep_rows <- if (isTRUE(perseus_lod)) {
    # Require at least 1 real detection across all samples;
    # all-NA rows would produce pure LOD noise in both groups → no signal.
    apply(sub_mat, 1, function(row) any(!is.na(row)))
  } else {
    # Evaluate coverage only on the two conditions being contrasted. In factorial
    # mode sub_mat spans every group in the experiment, so this must be explicit.
    i_test <- which(group_vec == cond_test)
    i_ref  <- which(group_vec == cond_ref)
    if (length(i_test) == 0L || length(i_ref) == 0L)
      return(list(error = paste0(
        "Could not find samples for both conditions (", cond_test,
        ", ", cond_ref, ") in the summarization cache.")))
    # min_valid_frac is a fraction, but what is actually applied is the integer
    # it resolves to at this replicate count (>=2 of 3 = 67%, >=2 of 4 = 50%).
    # Record both so the methods text can report the rule that was applied.
    .min_n_test <- ceiling(min_valid_frac * length(i_test))
    .min_n_ref  <- ceiling(min_valid_frac * length(i_ref))
    apply(sub_mat, 1, function(row) {
      ok_test <- sum(!is.na(row[i_test])) / length(i_test) >= min_valid_frac
      ok_ref  <- sum(!is.na(row[i_ref]))  / length(i_ref)  >= min_valid_frac
      if (isTRUE(require_both)) ok_test && ok_ref else ok_test || ok_ref
    })
  }

  orig_sub_mat <- sub_mat[keep_rows, , drop = FALSE]
  sub_mat      <- orig_sub_mat
  n_tested     <- nrow(sub_mat)

  if (n_tested == 0)
    return(list(error = paste0(
      "No proteins pass the missingness filter. ",
      if (isTRUE(perseus_lod))
        "(Perseus LOD mode: all proteins are 100% missing?)"
      else
        "(Standard mode: needs \u226550% valid in BOTH compared conditions)")))

  # 2. Imputation for remaining NAs
  #    Perseus LOD mode : per-column downshifted Gaussian (mean - 1.8*SD, width 0.3*SD)
  #    Standard mode    : k-NN (impute package)
  n_before_imp <- sum(is.na(sub_mat))
  n_knn_rows   <- 0L
  n_knn_values <- 0L
  n_left_na    <- 0L
  # Resolved eligibility count for provenance. Overridden below when the
  # deprecated fraction is supplied; initialised here so the return value is
  # defined even when the k-NN branch never runs.
  .max_miss    <- max_missing_per_cond
  if (isTRUE(perseus_lod)) {
    if (n_before_imp > 0) {
      set.seed(42L)
      mat_lod <- sub_mat
      for (j in seq_len(ncol(mat_lod))) {
        col    <- mat_lod[, j]
        obs    <- col[!is.na(col)]
        if (length(obs) < 2L) next
        col_sd  <- sd(obs)
        col_lod <- mean(obs) - 1.8 * col_sd
        na_idx  <- which(is.na(col))
        mat_lod[na_idx, j] <- rnorm(length(na_idx),
                                    mean = col_lod,
                                    sd   = 0.3 * col_sd)
      }
      sub_mat <- mat_lod
      logger(paste0("[DIA] Perseus LOD: filled ", n_before_imp - sum(is.na(sub_mat)),
                  " of ", n_before_imp, " NA values."))
      # n_imputed is derived from n_knn_values, which is only set in the standard
      # branch. Without this, LOD mode - the mode that imputes everything - reports
      # zero imputed values into the methods text.
      #
      # Count actual fills rather than assuming the loop filled everything: columns
      # with fewer than two observed values are skipped by the downshift loop above,
      # so their NAs survive. Assuming a complete fill puts a false number into the
      # generated methods text.
      n_knn_values <- n_before_imp - sum(is.na(sub_mat))
      n_left_na    <- sum(is.na(sub_mat))
      if (n_left_na > 0L)
        logger(paste0("[DIA] Perseus LOD: ", n_left_na, " missing values could not be ",
                      "filled (columns with fewer than two observed values) and were ",
                      "retained."))
    }
} else if (n_before_imp > 0) {

    if (!requireNamespace("impute", quietly = TRUE)) {
      logger("[DIA] WARNING: 'impute' package unavailable - no imputation performed. NAs passed to lmFit.")
      n_left_na <- n_before_imp
    } else {

      # ── Fix 2: per-condition eligibility ──────────────────────────────────
      # Only rows with LIGHT missingness inside BOTH compared conditions are
      # imputed. Everything heavier keeps its NAs; lmFit fits those rows on the
      # observations that actually exist. This is stricter than the testing
      # filter on purpose: passing the test filter earns you a model fit, not a
      # fabricated value.
      i_test <- which(group_vec == cond_test)
      i_ref  <- which(group_vec == cond_ref)

      # Eligibility is a COUNT of missing values per condition, not a fraction.
      # A fraction resolves to a different integer at every replicate count and
      # can collapse onto the testing filter (e.g. at n=3 both mean "<=1
      # missing"), silently removing the intended margin. A count behaves
      # identically at every n. Set max_missing_per_cond = 0L to disable
      # imputation entirely.
      n_miss_test <- rowSums(is.na(sub_mat[, i_test, drop = FALSE]))
      n_miss_ref  <- rowSums(is.na(sub_mat[, i_ref,  drop = FALSE]))

      # Back-compatibility: if a caller still passes the deprecated fraction,
      # resolve it to a count against the smaller condition and warn once.
      .max_miss <- max_missing_per_cond
      if (!is.null(impute_max_miss_per_cond)) {
        .n_min    <- min(length(i_test), length(i_ref))
        .max_miss <- floor(impute_max_miss_per_cond * .n_min)
        logger(paste0("[DIA] NOTE: impute_max_miss_per_cond is deprecated. ",
                      "Resolved ", impute_max_miss_per_cond, " to a count of ",
                      .max_miss, " against the smaller condition (n=", .n_min,
                      "). Use max_missing_per_cond instead."))
      }

      elig <- (n_miss_test <= .max_miss) &
              (n_miss_ref  <= .max_miss) &
              (n_miss_test + n_miss_ref > 0)      # skip complete rows
      elig[is.na(elig)] <- FALSE

      # All-conditions fit (not factorial): impute only the compared conditions'
      # columns. Eligibility is judged on those two conditions alone, so filling
      # other conditions' gaps would invent values the filter never vetted.
      # Two-condition fits are unchanged (the columns ARE the compared conditions).
      .imp_cols <- if (isTRUE(fit_all_groups) && !.is_factorial) sort(c(i_test, i_ref))
                   else seq_len(ncol(sub_mat))

      n_knn_rows   <- sum(elig)
      n_knn_values <- sum(is.na(sub_mat[elig, .imp_cols, drop = FALSE]))
      n_left_na    <- n_before_imp - n_knn_values

      if (n_knn_rows >= 2L && n_knn_values > 0L) {
        elig_mat <- sub_mat[elig, .imp_cols, drop = FALSE]
        imp_res  <- tryCatch(
          # ── Fix 3: rowmax can no longer fire ────────────────────────────
          # Eligibility is decided above, per condition. rowmax is set to 1 so
          # impute.knn never silently swaps k-NN for its sample-mean fallback.
          # maxp = Inf keeps it doing true global k-NN instead of the two-means
          # clustering approximation it defaults to above 1500 rows.
          impute::impute.knn(elig_mat, k = knn_k,
                             rowmax = 1, colmax = 0.8,
                             maxp = Inf, rng.seed = 42L)$data,
          error = function(e) { logger("[DIA] k-NN imputation failed - NAs retained."); NULL }
        )
        if (!is.null(imp_res)) {
          sub_mat[elig, .imp_cols] <- imp_res
        } else {
          n_knn_values <- 0L
          n_left_na    <- n_before_imp
        }
      } else {
        n_knn_values <- 0L
        n_left_na    <- n_before_imp
      }

      logger(paste0("[DIA] Imputation: ", n_knn_values, " values in ",
                    n_knn_rows, " proteins filled by k-NN (k=", knn_k,
                    "); ", n_left_na, " NAs retained for lmFit."))
    }
  }

  n_imputed <- n_knn_values   # honest count: what was actually filled

  # ── Perseus LOD mode: simple s0-adjusted Welch t-test, no covariates ────────
  # Matches the paper's methodology: plain 2-group comparison, s0=0.1 denominator
  # fudge (SAM-style), BH FDR.  No limma eBayes / no batch/covariate terms.
  if (isTRUE(perseus_lod)) {
    comp_label <- paste0(cond_test, " vs ", cond_ref)
    s0_val     <- 0.1   # Perseus default used in paper

    test_cols <- use_samples[group_vec == cond_test]
    ref_cols  <- use_samples[group_vec == cond_ref]

    prot_stats <- t(apply(sub_mat, 1, function(row) {
      v_t <- row[test_cols]
      v_r <- row[ref_cols]
      n_t <- length(v_t);  n_r <- length(v_r)
      m_t <- mean(v_t, na.rm = TRUE)
      m_r <- mean(v_r, na.rm = TRUE)
      s_t <- if (n_t > 1) var(v_t, na.rm = TRUE) else 0
      s_r <- if (n_r > 1) var(v_r, na.rm = TRUE) else 0
      se  <- sqrt(s_t / n_t + s_r / n_r)
      df  <- if (se > 0) {
        num <- (s_t / n_t + s_r / n_r)^2
        den <- (s_t / n_t)^2 / (n_t - 1) + (s_r / n_r)^2 / (n_r - 1)
        if (den > 0) num / den else n_t + n_r - 2
      } else n_t + n_r - 2
      t_mod <- (m_t - m_r) / (se + s0_val)
      pval  <- if (se + s0_val > 0) 2 * pt(-abs(t_mod), df = max(df, 1)) else 1
      c(log2fc   = m_t - m_r,
        ave      = mean(c(v_t, v_r), na.rm = TRUE),
        t_stat   = t_mod,
        se       = se,
        pvalue   = pval,
        df_welch = df)    # per-protein Welch-Satterthwaite df (for boundary curve)
    }))

    # ── Permutation FDR (Variant E: median null count, closest to Perseus) ──────
    # Uses median(colSums(|t_null| >= delta)) as the expected false-positive count.
    # Median is more robust than mean against rare extreme permutations, giving
    # ~1,363 sig at FDR<=0.05 vs the paper's 1,426  -" the closest reproducible match.
    logger("[DIA] Perseus LOD: computing permutation FDR (250 perms, median null)...")
    .lod_n_perm   <- 250L
    .lod_n_t      <- length(test_cols)
    .lod_n_r      <- length(ref_cols)
    .lod_all_cols <- c(test_cols, ref_cols)
    .lod_n_all    <- .lod_n_t + .lod_n_r
    .lod_abs_t    <- abs(prot_stats[, "t_stat"])
    set.seed(42L)
    .lod_null_mat <- vapply(seq_len(.lod_n_perm), function(b) {
      perm_idx <- sample(.lod_n_all)
      p_t      <- .lod_all_cols[perm_idx[seq_len(.lod_n_t)]]
      p_r      <- .lod_all_cols[perm_idx[(.lod_n_t + 1L):.lod_n_all]]
      mat_t    <- sub_mat[, p_t, drop = FALSE]
      mat_r    <- sub_mat[, p_r, drop = FALSE]
      m_t      <- rowMeans(mat_t, na.rm = TRUE)
      m_r      <- rowMeans(mat_r, na.rm = TRUE)
      s_t      <- rowSums(sweep(mat_t, 1L, m_t)^2, na.rm = TRUE) / max(.lod_n_t - 1L, 1L)
      s_r      <- rowSums(sweep(mat_r, 1L, m_r)^2, na.rm = TRUE) / max(.lod_n_r - 1L, 1L)
      se_p     <- sqrt(s_t / .lod_n_t + s_r / .lod_n_r)
      (m_t - m_r) / (se_p + s0_val)
    }, numeric(nrow(sub_mat)))                          # nrow x n_perm
    .lod_abs_null <- abs(.lod_null_mat)
    .lod_deltas   <- sort(unique(.lod_abs_t), decreasing = TRUE)
    .lod_exp_fp   <- vapply(.lod_deltas,
                            function(d) median(colSums(.lod_abs_null >= d, na.rm = TRUE)),
                            numeric(1))
    .lod_obs_cnt  <- vapply(.lod_deltas, function(d) sum(.lod_abs_t >= d), numeric(1))
    .lod_fdr_cv   <- ifelse(.lod_obs_cnt == 0L, 1.0,
                            pmin(.lod_exp_fp / .lod_obs_cnt, 1.0))
    .lod_fdr_mono <- rev(cummin(rev(.lod_fdr_cv)))
    .lod_rank_ord <- rank(-.lod_abs_t, ties.method = "first")
    perm_fdr_lod  <- .lod_fdr_mono[.lod_rank_ord]
    rm(.lod_null_mat, .lod_abs_null)   # free memory

    prots  <- rownames(sub_mat)
    genes  <- acc_to_gene[prots]
    genes[is.na(genes)] <- prots[is.na(genes)]

    # adj.pvalue now holds the permutation FDR - the value that actually
    # determines significance. The raw (s0-modified Welch) p-value is kept in
    # BOTH `pvalue` and `raw_pvalue`; the volcano plots raw_pvalue on the Y-axis
    # because the permutation FDR is step-like and unsuitable for a continuous
    # axis. Previously adj.pvalue held the raw p, so any consumer thresholding
    # adj.pvalue silently treated a raw p-value as an FDR.
    result_df <- data.frame(
      Protein    = prots,
      Gene       = genes,
      log2FC     = prot_stats[, "log2fc"],
      AveExpr    = prot_stats[, "ave"],
      t_stat     = prot_stats[, "t_stat"],
      pvalue     = prot_stats[, "pvalue"],
      raw_pvalue = prot_stats[, "pvalue"],   # raw p, for the volcano Y-axis only
      adj.pvalue = perm_fdr_lod,             # permutation FDR (significance)
      perm_fdr   = perm_fdr_lod,             # retained for existing consumers
      SE         = prot_stats[, "se"],
      df_welch   = prot_stats[, "df_welch"], # per-protein Welch df (for boundary curve)
      Comparison = comp_label,
      stringsAsFactors = FALSE
    )

    # Detection counts (pre-imputation)
    for (grp in c(cond_test, cond_ref)) {
      grp_cols <- use_samples[group_vec == grp]
      grp_cols <- grp_cols[grp_cols %in% colnames(orig_sub_mat)]
      col_name <- paste0("n_valid_", make.names(grp))
      result_df[[col_name]] <- apply(
        orig_sub_mat[, grp_cols, drop = FALSE], 1,
        function(x) sum(!is.na(x))
      )
    }

    result_df$PiScore <- abs(result_df$log2FC) *
      (-log10(pmax(perm_fdr_lod, 1e-300)))   # PiScore uses perm_fdr magnitude

    logger(paste0("[DIA] Perseus LOD t-test complete: ", nrow(result_df),
                  " proteins tested (s0=", s0_val,
                  ", permutation FDR \u2014 median null, 250 perms), ",
                  sum(perm_fdr_lod <= 0.05), " sig at FDR\u22640.05, ",
                  n_imputed, " NA values imputed."))

    return(list(
      result_df   = result_df,
      fit_object  = NULL,
      contrast    = comp_label,
      n_tested    = n_tested,
      n_imputed   = n_imputed,
      comparison  = comp_label,
      raw_mat     = orig_sub_mat,   # legacy key, retained for existing consumers
      imputed_mat = sub_mat,        # post-imputation (NAs filled)
      observed_mat = orig_sub_mat,  # nothing was imputed upstream on this exit
      model_mat   = sub_mat         # exactly what the LOD t-test consumed
    ))
  }

  # ── Standard mode: limma with covariates / eBayes ─────────────────────────

  # 3. Build design matrix
  group_f <- factor(group_vec)

  if (!is.null(factorial_design) && isTRUE(factorial_design$enabled)) {
    f1_map <- setNames(factorial_design$assignments$f1,
                       factorial_design$assignments$condition)
    f2_map <- setNames(factorial_design$assignments$f2,
                       factorial_design$assignments$condition)
    fused  <- paste(f1_map[group_vec], f2_map[group_vec], sep = "_")
    group_f <- factor(fused)
  }

  batch_in_model <- FALSE   # what actually happened, for provenance

  if (!is.null(batch_groups) && isTRUE(batch_for_model) &&
      nlevels(factor(setNames(
        rep(names(batch_groups), times = vapply(batch_groups, length, integer(1))),
        unlist(batch_groups))[use_samples])) >= 2L) {
    batch_map <- setNames(
      rep(names(batch_groups), times = vapply(batch_groups, length, integer(1))),
      unlist(batch_groups)
    )
    batch_vec  <- factor(batch_map[use_samples])
    batch_in_model <- TRUE
    design_mat <- model.matrix(~ 0 + group_f + batch_vec)
    # Rename group columns to match the non-batch convention so contrast strings
    # built with make.names(cond_test) / make.names(cond_ref) resolve correctly.
    # model.matrix names them "group_f<level>" but the contrast builder uses
    # make.names("<level>")  -" strip the "group_f" prefix before sanitising.
    n_group <- nlevels(group_f)
    colnames(design_mat)[seq_len(n_group)] <- make.names(levels(group_f))
    if (ncol(design_mat) > n_group)
      colnames(design_mat)[(n_group + 1L):ncol(design_mat)] <-
        make.names(colnames(design_mat)[(n_group + 1L):ncol(design_mat)])
  } else {
    if (!is.null(batch_groups) && isTRUE(batch_for_model))
      logger(paste0("[DIA] Batch covariate requested but only one batch is represented ",
                    "among the compared samples - batch term omitted from the model."))
    design_mat <- model.matrix(~ 0 + group_f)
    colnames(design_mat) <- make.names(levels(group_f))
  }

  # 3b. Sample/group variance weights (opt-in; DIA leaves this FALSE)
  # eBayes with a single fit assumes a common residual variance across groups.
  # If one condition is genuinely noisier, the pooled prior is contaminated in
  # both directions. The field's answer is weighting (Ritchie et al. 2006), not
  # refitting on subsets - selecting between pooled and pairwise after looking
  # at the data is model selection on the outcome and invalidates the p-values.
  # Note this is a different problem from eBayes(robust = TRUE), which
  # down-weights individual PROTEINS with aberrant variance and does nothing
  # about a condition that is uniformly noisier across the proteome.
  aw <- NULL
  n_array_weights_failed <- 0L
  if (isTRUE(use_array_weights)) {
    # arrayWeights needs rows without missing values. Estimate on the complete-case
    # subset (substantial after the global filter), then apply to the full fit.
    .cc <- stats::complete.cases(sub_mat)
    if (sum(.cc) >= 50L) {
      aw <- tryCatch(
        limma::arrayWeights(sub_mat[.cc, , drop = FALSE], design_mat,
                            var.group = group_vec),
        error = function(e) {
          logger("[DIA] arrayWeights failed - unweighted fit.")
          n_array_weights_failed <<- 1L
          NULL
        }
      )
      if (!is.null(aw) && (!all(is.finite(aw)) || any(aw <= 0)))
        { logger("[DIA] arrayWeights returned non-finite weights - unweighted fit.")
          n_array_weights_failed <- 1L; aw <- NULL }
    } else {
      logger("[DIA] Too few complete rows for arrayWeights - unweighted fit.")
      n_array_weights_failed <- 1L
    }
    if (!is.null(aw))
      logger(paste0("[DIA] arrayWeights (var.group): ",
                    paste(round(aw, 3), collapse = ", ")))
  }

  # 4. Fit model (with optional duplicateCorrelation for pairing)
  if (!is.null(subject_pairing)) {
    block_vec <- subject_pairing[use_samples]
    logger(paste0("[DIA] Paired design: running duplicateCorrelation with ",
                  length(unique(block_vec)), " subjects/blocks."))
    corfit    <- tryCatch(
      # weights is part of duplicateCorrelation's signature (limma >= 3.x); when
      # aw is NULL this is the documented default, so unweighted behaviour is
      # byte-identical to before.
      limma::duplicateCorrelation(sub_mat, design_mat, block = block_vec,
                                  weights = aw),
      error = function(e) NULL
    )
    # Retained NAs can leave duplicateCorrelation with too few complete blocks,
    # in which case it returns NaN or a degenerate consensus WITHOUT erroring.
    # lmFit(correlation = NaN) then yields silent garbage, so validate first.
    cc <- if (!is.null(corfit)) corfit$consensus.correlation else NA_real_
    if (!is.null(corfit) && is.finite(cc) && abs(cc) < 0.99) {
      logger(paste0("[DIA] duplicateCorrelation consensus: ", round(cc, 4)))
      fit <- limma::lmFit(sub_mat, design_mat,
                          block = block_vec,
                          correlation = cc,
                          weights = aw)
    } else if (!is.null(corfit)) {
      logger(paste0("[DIA] duplicateCorrelation returned an unusable consensus (", cc,
                    ") - falling back to unpaired fit. This can occur when retained ",
                    "NAs leave too few complete blocks."))
      fit <- limma::lmFit(sub_mat, design_mat, weights = aw)
    } else {
      logger("[DIA] duplicateCorrelation failed  -  falling back to unpaired fit.")
      fit <- limma::lmFit(sub_mat, design_mat, weights = aw)
    }
  } else {
    fit <- limma::lmFit(sub_mat, design_mat, weights = aw)
  }

  # 5. Build contrasts
  #    Factorial mode fits across ALL fused factor-level groups (better variance
  #    estimation), but the contrast tested is always the one the user selected.
  #    Previously this took combn()[[1]], i.e. the first pair alphabetically,
  #    which was unrelated to cond_test / cond_ref.
  group_levels <- levels(group_f)
  fac_lvl_test <- NULL
  fac_lvl_ref  <- NULL

  if (!is.null(factorial_design) && isTRUE(factorial_design$enabled)) {

    if (!(cond_test %in% names(f1_map)) || !(cond_ref %in% names(f1_map)))
      return(list(error = paste0(
        "Factorial design: condition(s) '", cond_test, "' / '", cond_ref,
        "' have no factor-level assignment. Re-open the wizard and assign ",
        "factor levels to every condition.")))

    fac_lvl_test <- paste(f1_map[[cond_test]], f2_map[[cond_test]], sep = "_")
    fac_lvl_ref  <- paste(f1_map[[cond_ref]],  f2_map[[cond_ref]],  sep = "_")

    if (identical(fac_lvl_test, fac_lvl_ref))
      return(list(error = paste0(
        "Factorial design: '", cond_test, "' and '", cond_ref,
        "' map to the same factor-level combination (", fac_lvl_test,
        "), so they cannot be contrasted. Check the factor assignments.")))

    missing_lvl <- setdiff(c(fac_lvl_test, fac_lvl_ref), group_levels)
    if (length(missing_lvl) > 0)
      return(list(error = paste0(
        "Factorial design: factor level(s) '", paste(missing_lvl, collapse = "', '"),
        "' not present in the fitted model.")))

    contrast_string <- paste0(make.names(fac_lvl_test), " - ", make.names(fac_lvl_ref))
    contrast_mat    <- limma::makeContrasts(
      contrasts = contrast_string,
      levels    = design_mat
    )
    coef_use <- 1L
    logger(paste0("[DIA] Factorial contrast: ", contrast_string,
                  "  (", cond_test, " vs ", cond_ref, ")"))

  } else {
    contrast_string <- paste0(make.names(cond_test), " - ", make.names(cond_ref))
    contrast_mat    <- limma::makeContrasts(
      contrasts = contrast_string,
      levels    = design_mat
    )
    coef_use <- 1L
  }

  # 6. Apply contrasts, eBayes
  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit3 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

  # 7. topTable
  top <- limma::topTable(fit3, coef = coef_use,
                          number = Inf, sort.by = "none",
                          adjust.method = "BH")

  # 8. Format result_df
  prots <- rownames(top)
  genes <- acc_to_gene[prots]
  genes[is.na(genes)] <- prots[is.na(genes)]

  comp_label <- if (!is.null(factorial_design) && isTRUE(factorial_design$enabled)) {
    paste0(cond_test, " vs ", cond_ref,
           " [factor levels: ", fac_lvl_test, " - ", fac_lvl_ref, "]")
  } else {
    paste0(cond_test, " vs ", cond_ref)
  }

  result_df <- data.frame(
    Protein    = prots,
    Gene       = genes,
    log2FC     = top$logFC,
    AveExpr    = top$AveExpr,
    t_stat     = top$t,
    pvalue     = top$P.Value,
    adj.pvalue = top$adj.P.Val,
    SE         = sqrt(fit3$s2.post) * fit3$stdev.unscaled[, coef_use],
    Comparison = comp_label,
    stringsAsFactors = FALSE
  )

  # 9. Detection counts (pre-imputation non-NA per group)
  for (grp in unique(group_vec)) {
    grp_cols <- use_samples[group_vec == grp]
    grp_cols <- grp_cols[grp_cols %in% colnames(orig_sub_mat)]
    col_name <- paste0("n_valid_", make.names(grp))
    result_df[[col_name]] <- apply(
      orig_sub_mat[, grp_cols, drop = FALSE], 1,
      function(x) sum(!is.na(x))
    )
  }

  # 10. PiScore
  result_df$PiScore <- abs(result_df$log2FC) *
    (-log10(pmax(result_df$adj.pvalue, 1e-300)))

  logger(paste0("[DIA] DE complete: ", nrow(result_df),
                " proteins tested, ", n_imputed, " NA values imputed."))

  list(
    result_df    = result_df,
    fit_object   = fit3,
    contrast     = comp_label,
    n_tested     = n_tested,
    n_imputed    = n_imputed,
    batch_in_model = batch_in_model,
    n_knn_rows   = n_knn_rows,
    n_left_na    = n_left_na,
    knn_k        = knn_k,
    max_missing_per_cond = .max_miss,
    min_valid_frac       = min_valid_frac,
    require_both         = require_both,
    array_weights_applied = !is.null(aw),
    array_weights         = if (!is.null(aw)) aw else NA_real_,
    n_array_weights_failed = n_array_weights_failed,
    min_n_test     = .min_n_test,
    min_n_ref      = .min_n_ref,
    n_test_samples = length(i_test),
    n_ref_samples  = length(i_ref),
    comparison   = comp_label,
    # Model scope ("auto" guard): what was fitted and why
    fit_all_groups          = isTRUE(fit_all_groups) || .is_factorial,
    fit_scope_guard         = fit_scope_guard,
    n_samples_fitted        = length(use_samples),
    n_conditions_fitted     = length(unique(group_vec)),
    # Exported per-sample abundances stay on the two compared conditions in an
    # all-conditions fit (the model matrix below keeps every fitted sample).
    raw_mat      = if (isTRUE(fit_all_groups) && !.is_factorial)
                     orig_sub_mat[, sort(c(i_test, i_ref)), drop = FALSE] else orig_sub_mat,
    imputed_mat  = if (isTRUE(fit_all_groups) && !.is_factorial)
                     sub_mat[, sort(c(i_test, i_ref)), drop = FALSE] else sub_mat,
    observed_mat = if (!is.null(observed_mat))
                     observed_mat[intersect(rownames(orig_sub_mat),
                                            rownames(observed_mat)), , drop = FALSE]
                   else orig_sub_mat,
    model_mat    = sub_mat         # exactly what was passed to lmFit
  )
}

# =============================================================================
# DIA PTM HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# ptm_quant_matrix()
# Per-sample linear quantity matrix of a parsed Spectronaut PTM file, for
# dia_normalization_check() in the PTM upload wizard: rows are precursors
# (parse_spectronaut_phospho: EG.TotalQuantity (Settings)) or PTM sites
# (parse_maxquant_ptm on a Spectronaut PTM pivot: PTM.Quantity), columns are
# the parser's sample names. Values are returned as read (no transformation).
# -----------------------------------------------------------------------------
ptm_quant_matrix <- function(parsed) {
  cols <- parsed$all_cols[parsed$quant_idx]
  if (is.null(parsed$raw_df) || length(cols) < 2L || !all(cols %in% colnames(parsed$raw_df)))
    return(NULL)
  m <- vapply(cols, function(cl) suppressWarnings(as.numeric(parsed$raw_df[[cl]])),
              numeric(nrow(parsed$raw_df)))
  m <- matrix(m, nrow = nrow(parsed$raw_df), dimnames = list(NULL, parsed$sample_names))
  m
}

# -----------------------------------------------------------------------------
# parse_spectronaut_phospho()
# Parse a Spectronaut phospho DIA report (precursor-level, EG.TotalQuantity).
# Returns: list(raw_dt, sample_names, n_rows, n_samples,
#               accessions_col, genes_col, pep_pos_col)
# -----------------------------------------------------------------------------
parse_spectronaut_phospho <- function(filepath, logger = message) {
  all_cols <- colnames(data.table::fread(filepath, nrows = 0L))

  # ── Validate required ID columns ──
  missing_id <- setdiff(c("PG.ProteinAccessions", "PG.Genes", "PEP.PeptidePosition",
                           "EG.ModifiedSequence"), all_cols)
  if (length(missing_id) > 0) {
    stop(paste0("Phospho DIA file is missing required column(s): ",
                paste(missing_id, collapse = ", ")))
  }

  # ── Detect per-sample metric columns ──
  quant_idx <- grep("\\.EG\\.TotalQuantity\\s*\\(Settings\\)$", all_cols)
  pos_idx   <- grep("\\.EG\\.PTMPositions\\s*\\[Phospho\\s*\\(STY\\)\\]$", all_cols)
  prob_idx  <- grep("\\.EG\\.PTMProbabilities\\s*\\[Phospho\\s*\\(STY\\)\\]$", all_cols)
  site_idx  <- grep("\\.EG\\.PTMSites\\s*\\[Phospho\\s*\\(STY\\)\\]$", all_cols)

  for (nm in c("EG.TotalQuantity", "EG.PTMPositions", "EG.PTMProbabilities", "EG.PTMSites")) {
    idx <- switch(nm,
      "EG.TotalQuantity"    = quant_idx,
      "EG.PTMPositions"     = pos_idx,
      "EG.PTMProbabilities" = prob_idx,
      "EG.PTMSites"         = site_idx
    )
    if (length(idx) < 2)
      stop(paste0("Phospho DIA file needs >=2 per-sample columns for ", nm,
                  "; found ", length(idx), "."))
  }

  # ── Extract sample names from EG.TotalQuantity columns ──
  raw_quant_cols <- all_cols[quant_idx]
  sample_names   <- sub("^\\[\\d+\\]\\s*", "", raw_quant_cols)
  sample_names   <- sub("\\.EG\\.TotalQuantity\\s*\\(Settings\\)$", "", sample_names)
  sample_names   <- trimws(sample_names)

  # ── Read file ──
  dt <- data.table::fread(filepath,
                           na.strings = c("NA", "NaN", "", "Filtered"),
                           header = TRUE)
  dt <- as.data.frame(dt, stringsAsFactors = FALSE)

  n_rows    <- nrow(dt)
  n_samples <- length(sample_names)
  logger(paste0("[DIA PTM] Phospho file parsed: ", n_rows, " rows x ", n_samples, " samples."))

  list(
    raw_df       = dt,
    sample_names = sample_names,
    n_rows       = n_rows,
    n_samples    = n_samples,
    all_cols     = all_cols,
    quant_idx    = quant_idx,
    pos_idx      = pos_idx,
    prob_idx     = prob_idx,
    site_idx     = site_idx
  )
}

# -----------------------------------------------------------------------------
# parse_maxquant_ptm()
# Parse a MaxQuant/Perseus PTM-level report (TSV/TXT).
# Expected format: site-level export with Perseus-style column prefixes.
#   Quantity columns: "[N] SampleName.PTM.Quantity"
#   Site ID:          "C: PTM.CollapseKey"  (e.g., ACC_S941_M1)
#   Mod type:         "C: PTM.ModificationTitle"
#   Protein (resolved): "T: PTM.ProteinId"
#   Flanking region:  "T: PTM.FlankingRegion"
#   Site position:    "T: PTM.SiteLocation"
#   Protein groups:   "T: PG.ProteinGroups", "T: PG.Genes"
#
# Probability column: In MaxQuant/Perseus exports this may appear as
#   "C: PTM.LocalizationProbability", "C: PTM.SiteProbability",
#   "C: PTM.LocProbabilityScore", or "PTM.Probability".
#   If absent, has_probability = FALSE and the caller should warn the user.
#
# Auto-filters to Phospho (STY) rows. Reports n_dropped_nonphospho.
#
# Returns: list(raw_df, sample_names, n_rows, n_samples,
#               collapse_key_col, protein_id_col, gene_col,
#               site_aa_col, site_loc_col, flanking_col,
#               prob_col, has_probability,
#               n_dropped_nonphospho, n_shared_rows)
# -----------------------------------------------------------------------------
parse_maxquant_ptm <- function(filepath, logger = message) {
  logger(paste0("[MQ PTM] Parsing MaxQuant PTM file: ", basename(filepath)))

  # Read header via readLines to avoid fread separator mis-detection on
  # column names that contain brackets (e.g. "[1] SampleName.PTM.Quantity").
  header_raw <- readLines(filepath, n = 1L, warn = FALSE)
  all_cols   <- strsplit(header_raw, "\t")[[1]]

  # ── Helper: find column by candidate names ──
  .find_col <- function(candidates) {
    idx <- which(all_cols %in% candidates)[1]
    if (is.na(idx)) NA_character_ else all_cols[idx]
  }

  # ── Locate quantity columns: [N] SampleName  OR  [N] SampleName.PTM.Quantity ──
  # Quantity columns are distinguished from all other [N]-numbered columns by
  # the absence of a leading prefix (N:, T:, C:). Both export formats handled.
  # The leading [ may also be absent in some exports (1] vs [1]).
  quant_idx <- grep("^\\[?\\d+\\]", all_cols)
  quant_idx <- quant_idx[!grepl("^[TNC]:\\s", all_cols[quant_idx])]
  if (length(quant_idx) < 2)
    stop(paste0("Fewer than 2 quantity columns found. ",
                "Expected columns starting with '[N] SampleName' ",
                "(optionally followed by .PTM.Quantity)."))

  # ── Locate required ID columns ──
  collapse_key_col <- .find_col(c("T: PTM.CollapseKey", "C: PTM.CollapseKey", "PTM.CollapseKey"))
  mod_title_col    <- .find_col(c("C: PTM.ModificationTitle", "T: PTM.ModificationTitle", "PTM.ModificationTitle"))
  protein_id_col   <- .find_col(c("T: PTM.ProteinId", "C: PTM.ProteinId", "PTM.ProteinId"))
  gene_col         <- .find_col(c("T: PG.Genes", "PG.Genes"))
  prot_groups_col  <- .find_col(c("T: PG.ProteinGroups", "PG.ProteinGroups"))
  site_aa_col      <- .find_col(c("T: PTM.SiteAA", "C: PTM.SiteAA", "PTM.SiteAA"))
  site_loc_col     <- .find_col(c("T: PTM.SiteLocation", "C: PTM.SiteLocation", "PTM.SiteLocation"))
  flanking_col     <- .find_col(c("T: PTM.FlankingRegion", "C: PTM.FlankingRegion", "PTM.FlankingRegion"))

  missing_req <- character(0)
  if (is.na(collapse_key_col)) missing_req <- c(missing_req, "PTM.CollapseKey")
  if (is.na(protein_id_col))   missing_req <- c(missing_req, "PTM.ProteinId")
  if (is.na(gene_col))         missing_req <- c(missing_req, "PG.Genes")
  if (length(missing_req) > 0)
    stop(paste0("MaxQuant PTM file is missing required column(s): ",
                paste(missing_req, collapse = ", ")))

  # ── Locate probability column (may be absent in pre-filtered exports) ──
  prob_candidates <- c(
    "C: PTM.LocalizationProbability", "PTM.LocalizationProbability",
    "C: PTM.SiteProbability",         "PTM.SiteProbability",
    "C: PTM.LocProbabilityScore",     "PTM.LocProbabilityScore",
    "C: PTM.Probability",             "PTM.Probability",
    "Localization prob",              "C: Localization prob"
  )
  prob_col        <- .find_col(prob_candidates)
  has_probability <- !is.na(prob_col)

  # ── Extract sample names ──
  raw_quant_cols <- all_cols[quant_idx]
  sample_names   <- sub("^\\[?\\d+\\]\\s*", "", raw_quant_cols)  # strip [N] or N] prefix
  sample_names   <- sub("\\.PTM\\.Quantity$", "", sample_names)   # strip suffix if present
  sample_names   <- trimws(sample_names)

  # ── Read file ──
  dt <- tryCatch(
    data.table::fread(filepath,
                      na.strings = c("NA", "NaN", "", "Filtered"),
                      data.table = FALSE, header = TRUE),
    error = function(e) stop(paste0("File read error: ", conditionMessage(e)))
  )
  # Ensure column names match our readLines-parsed all_cols (fread may mis-detect
  # the separator on files where column names contain '[' brackets).
  if (ncol(dt) == length(all_cols)) colnames(dt) <- all_cols

  n_rows_total <- nrow(dt)

  # ── Filter to Phospho (STY) only ──
  n_dropped_nonphospho <- 0L
  if (!is.na(mod_title_col) && mod_title_col %in% colnames(dt)) {
    phospho_mask       <- grepl("Phospho", dt[[mod_title_col]], ignore.case = TRUE)
    n_dropped_nonphospho <- sum(!phospho_mask)
    if (n_dropped_nonphospho > 0) {
      logger(paste0("[MQ PTM] Dropped ", n_dropped_nonphospho, " non-Phospho row(s) ",
                    "(e.g., Carbamidomethyl). Retaining only Phospho (STY)."))
    }
    dt <- dt[phospho_mask, , drop = FALSE]
  }

  n_rows <- nrow(dt)
  if (n_rows == 0)
    stop("No Phospho (STY) rows found after filtering. Check PTM.ModificationTitle column.")

  # ── Count shared-peptide rows ──
  n_shared_rows <- 0L
  if (!is.na(prot_groups_col) && prot_groups_col %in% colnames(dt)) {
    n_shared_rows <- sum(grepl(";", dt[[prot_groups_col]], fixed = TRUE))
  }

  n_samples <- length(sample_names)
  logger(paste0("[MQ PTM] Parse complete: ", n_rows, " phospho rows x ", n_samples,
                " samples. Has probability column: ", has_probability, "."))

  list(
    raw_df               = dt,
    sample_names         = sample_names,
    n_rows               = n_rows,
    n_samples            = n_samples,
    n_dropped_nonphospho = n_dropped_nonphospho,
    n_shared_rows        = n_shared_rows,
    collapse_key_col     = collapse_key_col,
    protein_id_col       = protein_id_col,
    gene_col             = gene_col,
    prot_groups_col      = prot_groups_col,
    site_aa_col          = site_aa_col,
    site_loc_col         = site_loc_col,
    flanking_col         = flanking_col,
    prob_col             = if (has_probability) prob_col else NULL,
    has_probability      = has_probability,
    quant_idx            = quant_idx,
    all_cols             = all_cols
  )
}

# -----------------------------------------------------------------------------
# build_pld_from_matrix()
# Convert a log2-abundance matrix + design_df into a ProteinLevelData data.frame.
# Columns match MSstatsPTM ProteinLevelData: Protein, Abundance, originalRUN,
# Condition, BioReplicate.
# -----------------------------------------------------------------------------
build_pld_from_matrix <- function(mat, design_df) {
  rows <- lapply(seq_len(nrow(mat)), function(i) {
    prot  <- rownames(mat)[i]
    vals  <- mat[i, ]
    runs  <- names(vals)
    valid <- !is.na(vals)
    if (!any(valid)) return(NULL)
    cond_vec <- design_df$condition[match(runs[valid], design_df$sample_name)]
    bio_vec  <- design_df$bioreplicate[match(runs[valid], design_df$sample_name)]
    data.frame(
      Protein      = prot,
      Abundance    = as.numeric(vals[valid]),
      originalRUN  = runs[valid],
      Condition    = cond_vec,
      BioReplicate = bio_vec,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# dia_ptm_summarize()
# Core processing pipeline for DIA phospho data.
# Steps: parse raw per-sample PTM cols -> consensus site ID (Phase 1) ->
#        per-sample quantity collection (Phase 2) -> prob filter ->
#        raw intensity summation per (site, sample) -> log2 transform ->
#        stoichiometric adjustment vs global protein matrix (skipped in
#        PTM-only mode, global_parsed = NULL).
#
# NOTE: Multiple precursors (charge states / modified variants) mapping to the
# same consensus site ID are aggregated by RAW INTENSITY SUMMATION before
# log2 transformation. TMP medpolish is NOT used for this roll-up because the
# median of log2-transformed values suppresses a dominant charge state when
# minor charge states are present (standard in field: MaxQuant, MSstats,
# Spectronaut all use summation for precursor -> peptide/site aggregation).
#
# Parameters:
#   phospho_parsed  : output of parse_spectronaut_phospho()
#   global_parsed   : output of parse_spectronaut_pivot(), or NULL for PTM-only
#                     mode (the cache then has ptm_mode = "ptm_only" and NULL
#                     adjusted / global layers; see ptm_layers())
#   cond_groups     : named list condition -> c(sample_names)
#   prob_threshold  : numeric 0-100 (default 75)
#   subject_pairing : named char sample->subject, or NULL
#   logger          : logging function
#
# Returns: list (see "Return structure" at end)
# -----------------------------------------------------------------------------
dia_ptm_summarize <- function(phospho_parsed, global_parsed, cond_groups,
                               prob_threshold  = 75,
                               # FLOOR, not the threshold: the effective per-condition
                               # requirement is max(this, ceiling(0.5 * n_samples_in_cond)).
                               min_valid_per_cond = 2L,
                               # "medpolish": Tukey median polish over the
                               # peptidoform x sample matrix per site (default).
                               # "sum": legacy linear-scale summation keyed on the
                               # site identifier - reproduces pre-Pass-4 behaviour.
                               site_summarization = c("medpolish", "sum"),
                               subject_pairing = NULL,
                               # Were the Spectronaut exports already cross-run
                               # normalized (Spectronaut default)? FALSE (default) keeps
                               # the earlier behaviour: all three centerings below.
                               #   ptm_input_normalized    -> skip centering of the site matrix
                               #   global_input_normalized -> skip centering of the global matrix
                               #   both                    -> also skip centering of the adjusted
                               #                              (site - protein) matrix: the per-run
                               #                              offsets of two normalized experiments
                               #                              cancel in the subtraction up to a constant
                               ptm_input_normalized    = FALSE,
                               global_input_normalized = FALSE,
                               logger          = message) {
  site_summarization <- match.arg(site_summarization)

  logger(paste0("[DIA PTM] Starting summarization (prob_threshold=", prob_threshold, "%)..."))

  dt           <- phospho_parsed$raw_df
  sample_names <- phospho_parsed$sample_names
  all_cols     <- phospho_parsed$all_cols

  # PTM-only mode (no global file): every phospho sample is used and steps 5-6
  # (protein summary, stoichiometric adjustment) are skipped.
  has_global <- !is.null(global_parsed)

  # Matched samples (intersection with global)
  matched_samples <- if (has_global) intersect(sample_names, global_parsed$sample_names)
                     else sample_names
  if (length(matched_samples) < 2)
    stop(if (has_global) "Fewer than 2 matched samples between global and phospho files."
         else "Fewer than 2 samples in the phospho file.")

  # Build cond_map restricted to matched samples
  cond_map <- setNames(
    rep(names(cond_groups), times = vapply(cond_groups, length, integer(1))),
    unlist(cond_groups)
  )
  cond_map <- cond_map[intersect(names(cond_map), matched_samples)]

  # ---------------------------------------------------------------------------
  # STEP 1: Pre-compute consensus site assignment per row (Phase 1)
  # ---------------------------------------------------------------------------
  # For each precursor row, gather localization probabilities from ALL matched
  # samples and take the per-candidate-position maximum across samples.
  # The top-N candidates (by multiplicity) chosen from these consensus
  # probabilities are used as the site ID for EVERY sample in that row,
  # guaranteeing all samples in the row map to the same site.
  # Probability ties are broken in N-to-C terminal order (R's which.max).
  # ---------------------------------------------------------------------------

  # Build column name lookup for each sample (once, before row loop)
  .sample_cols <- list()
  for (s in matched_samples) {
    find_col <- function(pattern) {
      hits <- grep(paste0(gsub(".", "\\.", s, fixed = TRUE), pattern), all_cols, value = TRUE)
      if (length(hits) == 0) return(NA_character_)
      hits[1]
    }
    quant_col <- find_col("\\.EG\\.TotalQuantity\\s*\\(Settings\\)$")
    pos_col   <- find_col("\\.EG\\.PTMPositions\\s*\\[Phospho\\s*\\(STY\\)\\]$")
    prob_col  <- find_col("\\.EG\\.PTMProbabilities\\s*\\[Phospho\\s*\\(STY\\)\\]$")
    site_col  <- find_col("\\.EG\\.PTMSites\\s*\\[Phospho\\s*\\(STY\\)\\]$")
    if (anyNA(c(quant_col, pos_col, prob_col, site_col))) {
      logger(paste0("[DIA PTM] Warning: could not find all 4 metric columns for sample '", s, "'  -  skipping."))
      next
    }
    .sample_cols[[s]] <- list(quant = quant_col, pos = pos_col,
                               prob = prob_col,  site = site_col)
  }
  matched_samples <- names(.sample_cols)  # restrict to samples with all 4 cols found

  # Phase 1: consensus site ID per precursor row
  consensus_cache <- vector("list", nrow(dt))
  for (r in seq_len(nrow(dt))) {
    # Obtain candidate positions and residues from the first non-missing sample
    ref_pos   <- NULL
    ref_sites <- NULL
    n_cand    <- 0L
    for (s in matched_samples) {
      pos_raw <- dt[[.sample_cols[[s]]$pos]][r]
      if (is.na(pos_raw)) next
      pos_str   <- as.character(pos_raw)
      pos_clean <- gsub("[^0-9;]", "", pos_str)
      pos_clean <- gsub(";+", ";", pos_clean)
      pos_clean <- gsub("^;|;$", "", pos_clean)
      if (!nzchar(pos_clean)) next
      pi <- suppressWarnings(as.integer(strsplit(pos_clean, ";")[[1]]))
      pi <- pi[!is.na(pi) & pi > 0]
      if (length(pi) == 0) next
      sc <- trimws(strsplit(as.character(dt[[.sample_cols[[s]]$site]][r]), ";")[[1]])
      n_cand <- min(length(pi), length(sc))
      if (n_cand == 0) next
      ref_pos   <- pi[seq_len(n_cand)]
      ref_sites <- sc[seq_len(n_cand)]
      break
    }
    if (is.null(ref_pos)) next

    # Accumulate max probability per candidate across all samples
    consensus_probs <- rep(0, n_cand)
    for (s in matched_samples) {
      prob_raw <- dt[[.sample_cols[[s]]$prob]][r]
      if (is.na(prob_raw)) next
      prob_str <- as.character(prob_raw)
      if (!nzchar(prob_str)) next
      pv   <- suppressWarnings(as.numeric(strsplit(prob_str, ";")[[1]]))
      n_pv <- min(length(pv), n_cand)
      if (n_pv == 0) next
      consensus_probs[seq_len(n_pv)] <- pmax(consensus_probs[seq_len(n_pv)],
                                              pv[seq_len(n_pv)], na.rm = TRUE)
    }

    # Phospho multiplicity from EG.ModifiedSequence
    mod_seq_raw    <- dt[["EG.ModifiedSequence"]][r]
    mod_seq_str    <- if (!is.na(mod_seq_raw)) as.character(mod_seq_raw) else ""
    phospho_charge <- {
      m <- gregexpr("[Phospho (STY)]", mod_seq_str, fixed = TRUE)[[1]]
      if (identical(m, -1L)) 1L else max(1L, length(m))
    }

    # Select top-N candidates by consensus probability; re-sort by position (N→C)
    n_select       <- min(phospho_charge, n_cand)
    sorted_by_prob <- order(consensus_probs, decreasing = TRUE)
    selected_idx   <- sort(sorted_by_prob[seq_len(n_select)])
    consensus_prob <- min(consensus_probs[selected_idx])  # weakest-link for threshold
    if (is.na(consensus_prob)) next

    # Protein accessions and peptide start positions
    acc_str       <- as.character(dt[["PG.ProteinAccessions"]][r])
    gen_str       <- as.character(dt[["PG.Genes"]][r])
    pp_str        <- as.character(dt[["PEP.PeptidePosition"]][r])
    accessions    <- trimws(strsplit(acc_str, ";")[[1]])
    genes         <- trimws(strsplit(gen_str, ";")[[1]])
    pep_positions <- suppressWarnings(as.integer(trimws(strsplit(pp_str, ";")[[1]])))
    if (length(pep_positions) < length(accessions)) {
      pep_positions <- c(pep_positions,
                         rep(pep_positions[length(pep_positions)],
                             length(accessions) - length(pep_positions)))
    }
    pep_positions <- pep_positions[seq_len(length(accessions))]
    pep_positions[is.na(pep_positions)] <- pep_positions[1]

    # Build consensus site ID
    is_shared    <- length(accessions) > 1
    sel_rel_pos  <- ref_pos[selected_idx]
    sel_residues <- ref_sites[selected_idx]

    if (!is_shared) {
      abs_pos_sel <- pep_positions[1] + sel_rel_pos - 1L
      sites_part  <- paste(paste0(sel_residues, abs_pos_sel), collapse = "+")
      site_id     <- paste0(accessions[1], "_", sites_part, "_Phospho_M", phospho_charge)
      gene_lbl    <- genes[1]
    } else {
      per_prot_parts <- vapply(seq_along(accessions), function(pi) {
        abs_pos_i <- pep_positions[pi] + sel_rel_pos - 1L
        sites_i   <- paste(paste0(sel_residues, abs_pos_i), collapse = "+")
        paste0(accessions[pi], "_", sites_i, "_Phospho_M", phospho_charge)
      }, character(1))
      site_id  <- paste(per_prot_parts, collapse = ";")
      gene_lbl <- paste(genes, collapse = ";")
    }

    consensus_cache[[r]] <- list(
      site_id        = site_id,
      gene_lbl       = gene_lbl,
      acc_str        = acc_str,
      consensus_prob = consensus_prob,
      is_shared      = is_shared,
      mod_seq        = mod_seq_str,   # diagnostic: peptidoform multiplicity
      selected_idx   = selected_idx   # diagnostic: per-sample localization check
    )
  }

  n_consensus_rows <- sum(!vapply(consensus_cache, is.null, logical(1)))
  logger(paste0("[DIA PTM] Phase 1 consensus: site IDs computed for ",
                n_consensus_rows, " / ", nrow(dt), " precursor rows."))

  # ---------------------------------------------------------------------------
  # STEP 1b: Collect per-sample quantities using consensus site IDs (Phase 2)
  # ---------------------------------------------------------------------------
  collector <- vector("list", nrow(dt) * length(matched_samples))
  coll_idx  <- 0L

  for (s in matched_samples) {
    quant_col <- .sample_cols[[s]]$quant
    for (r in seq_len(nrow(dt))) {
      cc <- consensus_cache[[r]]
      if (is.null(cc)) next
      quant_val <- dt[[quant_col]][r]
      if (is.na(quant_val) || quant_val <= 0) next
      # Diagnostic: this sample's OWN localization probability for this row,
      # taken at the selected candidate positions (weakest link, as the consensus
      # does).  NA when the sample reported no probability string for the row.
      own_prob <- {
        .p_raw <- dt[[.sample_cols[[s]]$prob]][r]
        if (is.na(.p_raw)) NA_real_ else {
          .p_str <- as.character(.p_raw)
          if (!nzchar(.p_str)) NA_real_ else {
            .pv   <- suppressWarnings(as.numeric(strsplit(.p_str, ";")[[1]]))
            .sel  <- cc$selected_idx[cc$selected_idx <= length(.pv)]
            .vsel <- .pv[.sel]
            .vsel <- .vsel[!is.na(.vsel)]
            if (length(.vsel) == 0L) NA_real_ else min(.vsel)
          }
        }
      }
      coll_idx <- coll_idx + 1L
      collector[[coll_idx]] <- list(
        row_idx   = r,
        sample    = s,
        site_id   = cc$site_id,
        gene_lbl  = cc$gene_lbl,
        acc_str   = cc$acc_str,
        quant_val = quant_val,
        best_prob = cc$consensus_prob,
        is_shared = cc$is_shared,
        mod_seq   = cc$mod_seq,
        own_prob  = own_prob
      )
    }
  }

  if (coll_idx == 0) stop("No valid PTM observations found after initial parsing.")
  collector <- collector[seq_len(coll_idx)]

  # Convert to data.frame
  coll_df <- data.frame(
    row_idx   = vapply(collector, `[[`, integer(1),  "row_idx"),
    sample    = vapply(collector, `[[`, character(1), "sample"),
    site_id   = vapply(collector, `[[`, character(1), "site_id"),
    gene_lbl  = vapply(collector, `[[`, character(1), "gene_lbl"),
    acc_str   = vapply(collector, `[[`, character(1), "acc_str"),
    quant_val = vapply(collector, `[[`, numeric(1),   "quant_val"),
    best_prob = vapply(collector, `[[`, numeric(1),   "best_prob"),
    is_shared = vapply(collector, `[[`, logical(1),   "is_shared"),
    mod_seq   = vapply(collector, `[[`, character(1), "mod_seq"),
    own_prob  = vapply(collector, `[[`, numeric(1),   "own_prob"),
    stringsAsFactors = FALSE
  )

  n_sites_input <- length(unique(coll_df$site_id))
  logger(paste0("[DIA PTM] Raw collection: ", nrow(coll_df), " observations, ",
                n_sites_input, " unique sites."))

  # Diagnostic: peptidoform multiplicity per site.
  # STEP 3 aggregates precursor rows by SUMMING on site_id, not on peptidoform,
  # so any distinct EG.ModifiedSequence mapping to the same site set is summed in.
  # Measured here only; the aggregation itself is unchanged.
  # Computed twice: over ALL collected sites here, and again below over only the
  # sites that survive the localization filter, since the retained figure is the
  # one that describes the analysed data.
  .pf_count <- function(df) {
    if (nrow(df) == 0L) return(list(n_multi = 0L, max_pf = 0L))
    .per_site <- tapply(df$mod_seq, df$site_id,
                        function(x) length(unique(x[nzchar(x)])))
    list(n_multi = sum(.per_site > 1L, na.rm = TRUE),
         max_pf  = if (length(.per_site) > 0L) max(.per_site, na.rm = TRUE) else 0L)
  }
  .pf_all <- .pf_count(coll_df)
  n_sites_multi_peptidoform_all <- .pf_all$n_multi
  max_peptidoforms_per_site_all <- .pf_all$max_pf
  logger(paste0("[DIA PTM] Peptidoform multiplicity (all collected sites): ",
                n_sites_multi_peptidoform_all,
                " site(s) carry more than one distinct EG.ModifiedSequence; ",
                "maximum per site = ", max_peptidoforms_per_site_all, "."))

  # ---------------------------------------------------------------------------
  # STEP 2: Probability threshold  -" max across samples per site
  # ---------------------------------------------------------------------------
  max_prob <- tapply(coll_df$best_prob, coll_df$site_id, max, na.rm = TRUE)
  thresh   <- prob_threshold / 100
  passing  <- names(max_prob)[max_prob >= thresh]
  coll_df  <- coll_df[coll_df$site_id %in% passing, ]

  n_sites_passing_prob <- length(unique(coll_df$site_id))
  logger(paste0("[DIA PTM] Probability filter (>=", prob_threshold, "%): ",
                n_sites_input, " -> ", n_sites_passing_prob, " sites."))

  # Peptidoform multiplicity again, now over the retained sites only.
  .pf_ret <- .pf_count(coll_df)
  n_sites_multi_peptidoform_retained <- .pf_ret$n_multi
  max_peptidoforms_per_site_retained <- .pf_ret$max_pf
  logger(paste0("[DIA PTM] Peptidoform multiplicity (retained sites): ",
                n_sites_multi_peptidoform_retained,
                " of ", n_sites_passing_prob,
                " site(s) carry more than one distinct EG.ModifiedSequence; ",
                "maximum per site = ", max_peptidoforms_per_site_retained, "."))

  # Diagnostic: per-sample localization evidence.
  # The consensus rule is per site across all samples, so a passing site's
  # intensity is collected from samples whose own probability for that row was
  # missing or below threshold.  Counted here only; no observation is skipped.
  n_obs_collected           <- nrow(coll_df)
  n_obs_no_own_prob         <- sum(is.na(coll_df$own_prob))
  n_obs_below_own_threshold <- sum(!is.na(coll_df$own_prob) & coll_df$own_prob < thresh)
  logger(paste0("[DIA PTM] Per-sample localization: of ", n_obs_collected,
                " retained (site, sample) observations, ", n_obs_no_own_prob,
                " (", if (n_obs_collected > 0L)
                  round(100 * n_obs_no_own_prob / n_obs_collected, 1) else 0,
                "%) had no probability recorded for that sample; ",
                n_obs_below_own_threshold,
                " (", if (n_obs_collected > 0L)
                  round(100 * n_obs_below_own_threshold / n_obs_collected, 1) else 0,
                "%) had a recorded probability below the threshold."))

  # Store max probability per site (0-1 scale)
  site_to_prob <- max_prob[passing]

  # ---------------------------------------------------------------------------
  # STEP 3: Build site × sample quantity matrix (raw intensity summation)
  # ---------------------------------------------------------------------------
  # Multiple precursor rows per (site, sample)  -" typically different charge
  # states of the same phosphopeptide  -" are aggregated by SUMMING the raw
  # intensities before log2 transformation.  This preserves the linear
  # additive physics of the mass spectrometer: a dominant +2 charge state
  # at 14 530 counts is not suppressed by co-detected +3/+4 traces at ~20
  # counts.  Summing first, then log2-transforming, is the approach used by
  # MaxQuant, Spectronaut (precursor→peptide), MSstats, and DIA-NN.
  # ---------------------------------------------------------------------------
  # NOTE (Pass 4): the header above describes site_summarization == "sum", which
  # is no longer the default. Charge-state summation within a peptidoform still
  # happens exactly as described; what changed is that the per-peptidoform sums
  # are combined by median polish instead of being summed across peptidoforms.
  all_sites   <- sort(unique(coll_df$site_id))
  ptm_abs_mat <- matrix(NA_real_,
                        nrow = length(all_sites),
                        ncol = length(matched_samples),
                        dimnames = list(all_sites, matched_samples))

  # Under site_summarization == "medpolish" the site estimate is instead obtained
  # by Tukey median polish over the peptidoform x sample matrix for that site:
  # charge states of one peptidoform are still summed on the linear scale (same
  # molecule), the per-(peptidoform, sample) sums are log2-transformed, and the
  # site value for sample j is overall + col[j].  A peptidoform that appears or
  # disappears between conditions then loads onto the ROW (peptidoform) term
  # rather than moving the sample estimate, which a plain sum over the site key
  # cannot do.  "sum" reproduces the pre-Pass-4 behaviour exactly.
  n_medpolish_fallback <- 0L

  for (site in all_sites) {
    rows_s <- coll_df[coll_df$site_id == site, ]

    # Number of contributing observations per sample, computed ONCE from the raw
    # rows.  This is the NA mask: a (site, sample) cell with zero contributing
    # observations must end as NA under BOTH methods.  medpolish returns
    # overall + col[j] for every column it was given, including empty ones, so
    # without this mask missingness would be turned into fabricated data.
    .ok      <- !is.na(rows_s$quant_val) & rows_s$quant_val > 0
    .n_obs   <- vapply(matched_samples,
                       function(s) sum(.ok & rows_s$sample == s), integer(1))
    .present <- matched_samples[.n_obs > 0L]
    if (length(.present) == 0L) next

    # Linear-scale sums per (peptidoform, sample); with "sum" the peptidoform
    # dimension is collapsed, which is the legacy behaviour.
    .grp  <- if (identical(site_summarization, "medpolish")) rows_s$mod_seq
             else rep("", nrow(rows_s))
    .grp[is.na(.grp) | !nzchar(.grp)] <- "(unknown_peptidoform)"
    .sums <- tapply(rows_s$quant_val[.ok],
                    list(.grp[.ok], rows_s$sample[.ok]), sum)

    .site_vals <- setNames(rep(NA_real_, length(.present)), .present)

    if (identical(site_summarization, "medpolish") && nrow(.sums) > 1L) {
      .lm <- log2(.sums[, intersect(colnames(.sums), .present), drop = FALSE])
      .lm[is.infinite(.lm) | is.nan(.lm)] <- NA
      # medpolish is unstable on all-NA rows/columns; drop them and map back by name.
      .lm <- .lm[rowSums(!is.na(.lm)) > 0L, colSums(!is.na(.lm)) > 0L, drop = FALSE]
      .mp <- if (nrow(.lm) > 1L && ncol(.lm) > 0L)
               tryCatch(stats::medpolish(.lm, na.rm = TRUE, trace.iter = FALSE),
                        error = function(e) NULL)
             else NULL
      if (!is.null(.mp)) {
        .est <- .mp$overall + .mp$col
        names(.est) <- colnames(.lm)
        if (all(is.finite(.est))) {
          .site_vals[names(.est)] <- .est
        } else {
          .mp <- NULL
        }
      }
      if (is.null(.mp)) {
        n_medpolish_fallback <- n_medpolish_fallback + 1L
        .tot <- colSums(.sums[, intersect(colnames(.sums), .present), drop = FALSE],
                        na.rm = TRUE)
        .site_vals[names(.tot)] <- ifelse(.tot > 0, log2(.tot), NA_real_)
      }
    } else {
      # "sum", or a single peptidoform under "medpolish".  A one-row matrix has no
      # row effect to estimate and medpolish degenerates on it, so the single-
      # peptidoform half of the data takes this same path: log2 of the summed
      # linear intensities, identical to the legacy result.
      .tot <- colSums(.sums[, intersect(colnames(.sums), .present), drop = FALSE],
                      na.rm = TRUE)
      .site_vals[names(.tot)] <- ifelse(.tot > 0, log2(.tot), NA_real_)
    }

    ptm_abs_mat[site, .present] <- .site_vals[.present]
  }

  if (identical(site_summarization, "medpolish"))
    logger(paste0("[DIA PTM] Site summarization: Tukey median polish over the ",
                  "peptidoform x sample matrix (stats::medpolish, defaults ",
                  "eps = 0.01, maxiter = 10); ", n_medpolish_fallback,
                  " site(s) fell back to linear-scale summation."))
  else
    logger("[DIA PTM] Site summarization: linear-scale summation keyed on site identifier.")

  # Replace Inf / NaN with NA
  ptm_abs_mat[is.infinite(ptm_abs_mat) | is.nan(ptm_abs_mat)] <- NA

  # Remove sites where no condition reaches its effective valid-value threshold.
  # Rule as implemented: a site is kept if AT LEAST ONE condition carries valid
  # values in at least half that condition's samples, subject to a floor of
  # min_valid_per_cond.  The threshold therefore SCALES with replicate count:
  #   effective[g] = max(min_valid_per_cond, ceiling(0.5 * n_samples_in_g))
  # At n = 3 that is 2 (the floor, unchanged); at n = 10 it is 5.  This is the
  # same >=50%-valid-in-at-least-one-group rule dia_limma_summarize() applies on
  # the global DIA route, so the two routes no longer disagree.
  # min_valid_per_cond is a FLOOR, not the threshold.
  .min_valid <- as.integer(min_valid_per_cond)
  .cond_names       <- names(cond_groups)
  .cond_sample_sets <- lapply(.cond_names, function(cond)
    intersect(cond_groups[[cond]], matched_samples))
  names(.cond_sample_sets) <- .cond_names
  min_valid_per_cond_effective <- vapply(.cond_sample_sets, function(samps)
    max(.min_valid, as.integer(ceiling(0.5 * length(samps)))), integer(1))
  names(min_valid_per_cond_effective) <- .cond_names
  .cond_valid <- vapply(.cond_names, function(cond) {
    samps <- .cond_sample_sets[[cond]]
    if (length(samps) == 0L) return(rep(FALSE, nrow(ptm_abs_mat)))
    rowSums(!is.na(ptm_abs_mat[, samps, drop = FALSE])) >=
      min_valid_per_cond_effective[[cond]]
  }, logical(nrow(ptm_abs_mat)))
  .keep          <- if (is.matrix(.cond_valid)) apply(.cond_valid, 1, any) else .cond_valid
  ptm_abs_mat    <- ptm_abs_mat[.keep, , drop = FALSE]
  all_sites      <- rownames(ptm_abs_mat)
  n_sites_output <- nrow(ptm_abs_mat)
  logger(paste0("[DIA PTM] Missingness filter: floor = ", .min_valid,
                ", effective per-condition thresholds (>=50% of samples) = ",
                paste0(names(min_valid_per_cond_effective), ":",
                       as.integer(min_valid_per_cond_effective), collapse = ", "), "."))
  logger(paste0("[DIA PTM] After missingness filter: ", n_sites_output, " sites retained."))

  # Normalize PTM matrix with column-wise median centering.
  # dia_limma_summarize() applies median centering to the global protein matrix;
  # without this step ptm_abs_mat is on raw log2 scale while global_mat is
  # centered, making the per-sample stoichiometric subtraction scale-asymmetric.
  # Skipped when the phospho export was already cross-run normalized in Spectronaut
  # (the constants are still computed and the shift they would apply is reported).
  # In PTM-only mode this is the only cross-sample normalization of the sites.
  .ptm_col_med <- apply(ptm_abs_mat, 2, median, na.rm = TRUE)
  .ptm_grd_med <- median(.ptm_col_med, na.rm = TRUE)
  ptm_centering_max_shift <- max(abs(.ptm_grd_med - .ptm_col_med), na.rm = TRUE)
  if (!isTRUE(ptm_input_normalized)) {
    for (.j in seq_len(ncol(ptm_abs_mat)))
      if (!is.na(.ptm_col_med[.j]))
        ptm_abs_mat[, .j] <- ptm_abs_mat[, .j] + (.ptm_grd_med - .ptm_col_med[.j])
  } else {
    logger(paste0("[DIA PTM] Site-matrix median centering NOT applied (phospho export cross-run ",
                  "normalized); it would have shifted samples by up to ",
                  signif(ptm_centering_max_shift, 3), " log2."))
  }

  # Gene label per site, from the collector (both modes)
  site_to_gene <- setNames(coll_df$gene_lbl[match(all_sites, coll_df$site_id)], all_sites)
  .no_lbl <- is.na(site_to_gene)
  site_to_gene[.no_lbl] <- sub("_Phospho.*$", "", all_sites[.no_lbl])

  # PTM-only mode: no protein layer, no adjusted matrix and no shared-peptide
  # expansion (a shared site stays one absolute row; ptm_de_annotate_views()
  # flags it Is_Shared).  Steps 5-6 fill these when a global file was given.
  global_cache          <- NULL
  global_mat            <- NULL
  ptm_adj_mat           <- NULL
  adj_to_abs_map        <- NULL
  n_shared_expanded     <- 0L
  n_no_parent_protein   <- 0L
  adj_centering_delta   <- NULL
  adj_centering_applied <- NULL

  if (has_global) {
    # ---------------------------------------------------------------------------
    # STEP 5: Global protein matrix via dia_limma_summarize
    # ---------------------------------------------------------------------------
    global_cache <- tryCatch(
      dia_limma_summarize(global_parsed, cond_groups,
                          subject_pairing = subject_pairing,
                          apply_centering = !isTRUE(global_input_normalized),
                          logger = logger),
      error = function(e) stop(paste0("Global DIA summarization failed: ", e$message))
    )
    global_mat <- global_cache$protein_matrix  # log2, median-centered, proteins × samples
    # Restrict to matched_samples
    global_mat <- global_mat[, intersect(colnames(global_mat), matched_samples), drop = FALSE]

    # ---------------------------------------------------------------------------
    # STEP 6: Stoichiometric adjustment + shared peptide expansion
    # ---------------------------------------------------------------------------
    # For each site_id in ptm_abs_mat:
    #   non-shared -> single adjusted row keyed as site_id
    #   shared     -> expand into one row per protein (unique _SharedInference keys)

    adj_rows       <- list()
    adj_to_abs_map <- character(0)   # named: adjusted row name -> absolute site_id
    # n_shared_expanded and n_no_parent_protein (adjusted rows that are all-NA
    # because the parent protein was never quantified - they fail the DE
    # missingness filter and would otherwise disappear silently) start at 0 above.

    for (site in all_sites) {
      is_shared_site <- grepl(";", site, fixed = TRUE)

      if (!is_shared_site) {
        # ── Non-shared site ──
        # Parse parent accession: everything before the first _<Residue><Pos>[+...] _Phospho
        parent_acc <- sub("_(?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+_Phospho.*$", "", site, perl = TRUE)

        ptm_row  <- ptm_abs_mat[site, ]
        if (parent_acc %in% rownames(global_mat)) {
          prot_row <- global_mat[parent_acc, ]
          adj_row  <- ptm_row - prot_row
        } else {
          adj_row  <- rep(NA_real_, length(ptm_row))
          names(adj_row) <- names(ptm_row)
          n_no_parent_protein <- n_no_parent_protein + 1L
        }
        adj_rows[[site]]       <- adj_row
        adj_to_abs_map[site]   <- site

      } else {
        # ── Shared peptide: expand into per-protein hypothesis rows ──
        # site is e.g. "KRT6A_S44_Phospho;KRT6B_S44_Phospho"
        per_prot_sites <- trimws(strsplit(site, ";")[[1]])
        ptm_row        <- ptm_abs_mat[site, ]

        for (pps in per_prot_sites) {
          # ── COLLISION PREVENTION: unique row name ──
          adj_row_name <- paste0(pps, "_SharedInference")
          parent_acc   <- sub("_(?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+_Phospho.*$", "", pps, perl = TRUE)

          if (parent_acc %in% rownames(global_mat)) {
            prot_row <- global_mat[parent_acc, ]
            adj_row  <- ptm_row - prot_row
          } else {
            adj_row  <- rep(NA_real_, length(ptm_row))
            names(adj_row) <- names(ptm_row)
            n_no_parent_protein <- n_no_parent_protein + 1L
          }

          adj_rows[[adj_row_name]]     <- adj_row
          adj_to_abs_map[adj_row_name] <- site
          n_shared_expanded            <- n_shared_expanded + 1L
        }
      }
    }

    # Build adjusted matrix
    adj_names    <- names(adj_rows)
    ptm_adj_mat  <- matrix(NA_real_,
                            nrow = length(adj_names),
                            ncol = length(matched_samples),
                            dimnames = list(adj_names, matched_samples))
    for (nm in adj_names) {
      v <- adj_rows[[nm]]
      ptm_adj_mat[nm, names(v)] <- v
    }

    # Column-wise median centering of the adjusted matrix.
    # ptm_abs_mat is centered on an all-observed basis; global_mat is centered by
    # dia_limma_summarize() on a complete-case basis when enough complete rows
    # exist.  The two offsets do not cancel under subtraction, so the residue
    # enters every adjusted contrast as a per-column constant.  Same construction
    # as the PTM matrix centering above.
    # When both exports were cross-run normalized upstream (and neither matrix was
    # centered here), there is no residue of mismatched centering bases to remove,
    # so this third centering is skipped as well.
    .adj_col_med <- apply(ptm_adj_mat, 2, median, na.rm = TRUE)
    .adj_grd_med <- median(.adj_col_med, na.rm = TRUE)
    adj_centering_delta <- .adj_grd_med - .adj_col_med
    adj_centering_applied <- !(isTRUE(ptm_input_normalized) && isTRUE(global_input_normalized))
    if (adj_centering_applied) {
      for (.j in seq_len(ncol(ptm_adj_mat)))
        if (!is.na(.adj_col_med[.j]))
          ptm_adj_mat[, .j] <- ptm_adj_mat[, .j] + (.adj_grd_med - .adj_col_med[.j])
    }
    logger(paste0("[DIA PTM] Adjusted matrix centering", if (!adj_centering_applied) " NOT applied (both exports normalized upstream; would be)" else "", ": median shift range ",
                  round(min(adj_centering_delta, na.rm = TRUE), 3), " to ",
                  round(max(adj_centering_delta, na.rm = TRUE), 3),
                  " log2 (spread ",
                  round(diff(range(adj_centering_delta, na.rm = TRUE)), 3), ")."))

    logger(paste0("[DIA PTM] Stoichiometry: ", nrow(ptm_adj_mat), " adjusted rows (",
                  n_shared_expanded, " from shared peptide expansion)."))
    logger(paste0("[DIA PTM] Parent protein not quantified for ", n_no_parent_protein,
                  " adjusted row(s) - these are all-NA and drop out at the DE missingness filter."))
  } else {
    logger("[DIA PTM] PTM-only mode: no global file - protein summary and stoichiometric adjustment skipped.")
  }

  # ---------------------------------------------------------------------------
  # STEP 7: Build design_df + missingness map
  # ---------------------------------------------------------------------------
  all_samples <- intersect(unlist(cond_groups), matched_samples)
  bio_rep_vec <- ave(seq_along(all_samples),
                     cond_map[all_samples],
                     FUN = function(i) paste0(cond_map[all_samples][i[1]], "_", seq_along(i)))

  if (!is.null(subject_pairing)) {
    bio_rep_vec <- subject_pairing[all_samples]
  }

  design_df <- data.frame(
    sample_name  = all_samples,
    condition    = cond_map[all_samples],
    bioreplicate = bio_rep_vec,
    stringsAsFactors = FALSE
  )

  # Missingness map: for each site × condition, fraction of non-NA samples
  conds <- unique(design_df$condition)
  miss_rows <- lapply(all_sites, function(site) {
    row_mat <- ptm_abs_mat[site, all_samples, drop = FALSE]
    lapply(conds, function(cond) {
      samps_in_cond <- design_df$sample_name[design_df$condition == cond]
      vals_in_cond  <- row_mat[, samps_in_cond, drop = FALSE]
      n_total   <- ncol(vals_in_cond)
      n_detected <- sum(!is.na(vals_in_cond))
      data.frame(Site_ID = site, Condition = cond,
                 Detected = n_detected, Total = n_total,
                 Fraction = if (n_total > 0) n_detected / n_total else 0,
                 stringsAsFactors = FALSE)
    })
  })
  missingness_map <- do.call(rbind, unlist(miss_rows, recursive = FALSE))

  logger(paste0("[DIA PTM] Summarization complete: ", n_sites_output, " sites, ",
                length(matched_samples), " samples, ",
                length(conds), " conditions."))

  # ---------------------------------------------------------------------------
  # STEP 8: Return
  # ---------------------------------------------------------------------------
  list(
    ptm_mode          = if (has_global) "with_global" else "ptm_only",
    ptm_abs_matrix    = ptm_abs_mat,
    ptm_adj_matrix    = ptm_adj_mat,
    global_matrix     = global_mat,
    design_df         = design_df,
    cond_groups       = cond_groups,
    acc_to_gene       = global_cache$acc_to_gene,
    site_to_gene      = site_to_gene,
    site_to_prob      = site_to_prob,
    adj_to_abs_map    = adj_to_abs_map,
    missingness_map   = missingness_map,
    n_sites_input     = n_sites_input,
    n_sites_passing_prob = n_sites_passing_prob,
    n_sites_output    = n_sites_output,
    n_sites_adjusted  = if (has_global) nrow(ptm_adj_mat) else NULL,
    n_shared_expanded = n_shared_expanded,
    n_no_parent_protein = n_no_parent_protein,
    n_sites_multi_peptidoform_all      = n_sites_multi_peptidoform_all,
    max_peptidoforms_per_site_all      = max_peptidoforms_per_site_all,
    n_sites_multi_peptidoform_retained = n_sites_multi_peptidoform_retained,
    max_peptidoforms_per_site_retained = max_peptidoforms_per_site_retained,
    n_obs_collected           = n_obs_collected,
    n_obs_no_own_prob         = n_obs_no_own_prob,
    n_obs_below_own_threshold = n_obs_below_own_threshold,
    # the shift actually applied to the adjusted matrix (NULL when skipped, and
    # both NULL in PTM-only mode, which has no adjusted matrix)
    adj_centering_delta = if (isTRUE(adj_centering_applied)) adj_centering_delta else NULL,
    adj_centering_delta_not_applied = if (isFALSE(adj_centering_applied)) adj_centering_delta else NULL,
    ptm_centering_applied    = !isTRUE(ptm_input_normalized),
    global_centering_applied = if (has_global) !isTRUE(global_input_normalized) else NULL,
    adj_centering_applied    = adj_centering_applied,
    ptm_centering_max_shift_log2    = ptm_centering_max_shift,
    global_centering_max_shift_log2 = global_cache$normalization_max_shift_log2,
    min_valid_per_cond  = .min_valid,
    min_valid_per_cond_effective = min_valid_per_cond_effective,
    site_summarization  = site_summarization,
    n_medpolish_fallback = n_medpolish_fallback,
    normalization_basis      = global_cache$normalization_basis,
    normalization_delta_log2 = global_cache$normalization_delta_log2,
    matched_samples   = matched_samples,
    prob_threshold    = prob_threshold,
    pipeline_label    = if (has_global) "DIA PTM (limma, intensity summation + stoichiometric adjustment)"
                        else "DIA PTM, PTM only (limma, intensity summation; no stoichiometric adjustment)",
    software          = "spectronaut",
    normalization     = if (has_global) "raw_intensity_summation_log2_stoichiometric_adjustment"
                        else "raw_intensity_summation_log2"
  )
}

# -----------------------------------------------------------------------------
# maxquant_ptm_summarize()
# Processing pipeline for MaxQuant/Perseus PTM data.
# Unlike dia_ptm_summarize(), site IDs are taken directly from PTM.CollapseKey  -"
# no consensus site localization is needed. Stoichiometric adjustment is applied
# identically to the Spectronaut pipeline.
#
# If the PTM file has no probability column, all sites are treated as passing;
# the caller (wizard) must inform the user before invoking this function.
#
# Returns: same cache structure as dia_ptm_summarize(), with software = "maxquant"
# -----------------------------------------------------------------------------
maxquant_ptm_summarize <- function(ptm_parsed, global_parsed, cond_groups,
                                    prob_threshold  = 75,
                                    subject_pairing = NULL,
                                    acquisition_mode = "dda",   # "dda" | "dia"
                                    # Spectronaut PTM pivot path only (acquisition_mode "dia"):
                                    # were the PTM / global exports already cross-run normalized?
                                    # Same rules as dia_ptm_summarize(). Ignored on the DDA path,
                                    # which stays bit-identical.
                                    ptm_input_normalized    = FALSE,
                                    global_input_normalized = FALSE,
                                    logger          = message) {
  .is_dia     <- identical(acquisition_mode, "dia")
  .ptm_norm   <- .is_dia && isTRUE(ptm_input_normalized)
  .glob_norm  <- .is_dia && isTRUE(global_input_normalized)
  logger(paste0("[MQ PTM] Starting summarization (prob_threshold=",
                prob_threshold, "%, acquisition=", acquisition_mode, ")..."))

  dt           <- ptm_parsed$raw_df
  sample_names <- ptm_parsed$sample_names
  all_cols     <- ptm_parsed$all_cols
  key_col      <- ptm_parsed$collapse_key_col
  prot_id_col  <- ptm_parsed$protein_id_col
  gene_col     <- ptm_parsed$gene_col
  has_prob     <- ptm_parsed$has_probability
  prob_col     <- ptm_parsed$prob_col

  # ── Matched samples (intersection with global) ──
  global_samples  <- global_parsed$sample_names
  matched_samples <- intersect(sample_names, global_samples)

  # Positional fallback: when PTM and global files have different sample name
  # conventions (e.g. pp_01 vs wp_01) but the same number and order of samples.
  ptm_to_global_map <- NULL
  if (length(matched_samples) < 2 &&
      length(sample_names) == length(global_samples) &&
      length(sample_names) >= 2) {
    ptm_to_global_map <- setNames(global_samples, sample_names)
    matched_samples   <- sample_names
    logger(paste0("[MQ PTM] Sample names differ between PTM and global files. ",
                  "Using positional matching (", length(sample_names), " samples). ",
                  "PTM: ", paste(head(sample_names, 3), collapse = ", "),
                  " -> Global: ", paste(head(global_samples, 3), collapse = ", "), " ..."))
  } else if (length(matched_samples) < 2) {
    stop(paste0("Fewer than 2 matched samples between global and PTM files. ",
                "PTM samples: ", paste(head(sample_names, 4), collapse = ", "),
                " | Global samples: ", paste(head(global_samples, 4), collapse = ", ")))
  }

  cond_map <- setNames(
    rep(names(cond_groups), times = vapply(cond_groups, length, integer(1))),
    unlist(cond_groups)
  )
  cond_map <- cond_map[intersect(names(cond_map), matched_samples)]

  # ── Quantity column map for matched samples ──
  # Handles both: '[N] SampleName.PTM.Quantity' and bare '[N] SampleName' / 'N] SampleName'.
  .quant_col_for_sample <- function(sname) {
    escaped <- gsub("([+.])", "\\\\\\1", sname)
    # Try suffixed format first (old: [N] SampleName.PTM.Quantity)
    hits <- grep(paste0("^\\[?\\d+\\]\\s*", escaped, "\\.PTM\\.Quantity$"),
                 all_cols, value = TRUE)
    # Fall back to bare format (new: [N] SampleName  or  N] SampleName)
    if (length(hits) == 0)
      hits <- grep(paste0("^\\[?\\d+\\]\\s*", escaped, "$"),
                   all_cols, value = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }
  quant_col_map <- vapply(matched_samples, .quant_col_for_sample, character(1))
  matched_samples <- matched_samples[!is.na(quant_col_map)]
  quant_col_map   <- quant_col_map[!is.na(quant_col_map)]

  if (length(matched_samples) < 2)
    stop("Could not find PTM.Quantity columns for >=2 matched samples.")

  # ── Gene/protein lookup from data ──
  site_keys    <- as.character(dt[[key_col]])
  protein_ids  <- as.character(dt[[prot_id_col]])
  gene_labels  <- if (!is.na(gene_col) && gene_col %in% colnames(dt))
    as.character(dt[[gene_col]]) else rep("", nrow(dt))
  gene_labels  <- sub(";.*", "", gene_labels)  # primary gene only

  n_sites_input <- length(unique(site_keys))

  # ── Build per-site probability (max across available column) ──
  if (has_prob && !is.null(prob_col) && prob_col %in% colnames(dt)) {
    raw_probs    <- suppressWarnings(as.numeric(dt[[prob_col]]))
    raw_probs[is.na(raw_probs)] <- 0
    # Group by site key, take max
    site_prob_max <- tapply(raw_probs, site_keys, max, na.rm = TRUE)
  } else {
    # No probability column: treat all sites as passing (prob = 1.0)
    unique_keys   <- unique(site_keys)
    site_prob_max <- setNames(rep(1.0, length(unique_keys)), unique_keys)
  }

  # ── Build coll_df (site x sample observations) ──
  rows_list <- vector("list", length(matched_samples))
  for (si in seq_along(matched_samples)) {
    s      <- matched_samples[si]
    q_col  <- quant_col_map[si]
    vals   <- if (q_col %in% colnames(dt)) dt[[q_col]] else rep(NA_real_, nrow(dt))
    valid  <- !is.na(vals) & vals > 0
    if (!any(valid)) next
    rows_list[[si]] <- data.frame(
      sample   = s,
      site_id  = site_keys[valid],
      prot_id  = protein_ids[valid],
      gene_lbl = gene_labels[valid],
      quant    = as.numeric(vals[valid]),
      stringsAsFactors = FALSE
    )
  }
  coll_df <- do.call(rbind, rows_list[!vapply(rows_list, is.null, logical(1))])
  if (is.null(coll_df) || nrow(coll_df) == 0)
    stop("No valid PTM observations found after initial parsing.")

  # ── Probability filter ──
  thresh   <- prob_threshold / 100
  passing  <- names(site_prob_max)[site_prob_max >= thresh]
  coll_df  <- coll_df[coll_df$site_id %in% passing, ]
  n_sites_output <- length(unique(coll_df$site_id))
  logger(paste0("[MQ PTM] Probability filter (>=", prob_threshold, "%): ",
                n_sites_input, " -> ", n_sites_output, " sites."))

  site_to_prob <- site_prob_max[passing]

  # ── Build site × sample raw-sum matrix, then log2 ──
  all_sites   <- sort(unique(coll_df$site_id))
  ptm_abs_mat <- matrix(NA_real_,
                        nrow = length(all_sites),
                        ncol = length(matched_samples),
                        dimnames = list(all_sites, matched_samples))

  for (site in all_sites) {
    rows_s <- coll_df[coll_df$site_id == site, ]
    for (s in matched_samples) {
      vals <- rows_s$quant[rows_s$sample == s]
      vals <- vals[!is.na(vals) & vals > 0]
      if (length(vals) == 0) next
      total_raw <- sum(vals)
      if (total_raw > 0) ptm_abs_mat[site, s] <- log2(total_raw)
    }
  }
  ptm_abs_mat[is.infinite(ptm_abs_mat) | is.nan(ptm_abs_mat)] <- NA

  # Remove sites where no condition reaches its valid-value threshold.
  # ROUTE GATE: this function serves both maxquant wizard paths. On the DIA
  # (Spectronaut PTM pivot) path the threshold SCALES with replicate count -
  # max(2, ceiling(0.5 * n_samples_in_cond)) - matching dia_ptm_summarize() and
  # the global DIA route. The DDA MaxQuant path keeps the plain floor of 2 and
  # is bit-identical to before; it is excluded under rule 4.
  .mq_floor         <- 2L
  .cond_names       <- names(cond_groups)
  .cond_sample_sets <- lapply(.cond_names, function(cond)
    intersect(cond_groups[[cond]], matched_samples))
  names(.cond_sample_sets) <- .cond_names
  min_valid_per_cond_effective <- vapply(.cond_sample_sets, function(samps)
    if (identical(acquisition_mode, "dia"))
      max(.mq_floor, as.integer(ceiling(0.5 * length(samps))))
    else
      .mq_floor,
    integer(1))
  names(min_valid_per_cond_effective) <- .cond_names
  .cond_valid <- vapply(.cond_names, function(cond) {
    samps <- .cond_sample_sets[[cond]]
    if (length(samps) == 0L) return(rep(FALSE, nrow(ptm_abs_mat)))
    rowSums(!is.na(ptm_abs_mat[, samps, drop = FALSE])) >=
      min_valid_per_cond_effective[[cond]]
  }, logical(nrow(ptm_abs_mat)))
  .keep          <- if (is.matrix(.cond_valid)) apply(.cond_valid, 1, any) else .cond_valid
  ptm_abs_mat    <- ptm_abs_mat[.keep, , drop = FALSE]
  all_sites      <- rownames(ptm_abs_mat)
  n_sites_output <- nrow(ptm_abs_mat)
  logger(paste0("[MQ PTM] Missingness filter: floor = ", .mq_floor,
                if (identical(acquisition_mode, "dia"))
                  ", effective per-condition thresholds (>=50% of samples) = "
                else
                  " (DDA path: plain floor, no scaling) = ",
                paste0(names(min_valid_per_cond_effective), ":",
                       as.integer(min_valid_per_cond_effective), collapse = ", "), "."))
  logger(paste0("[MQ PTM] After missingness filter: ", n_sites_output, " sites retained."))

  # Normalize PTM matrix with column-wise median centering.
  # For the DDA path this must precede MinProb imputation so that LOD values
  # are computed relative to already-normalized column means and SDs.
  # Spectronaut pivot path: skipped when the PTM export was cross-run normalized.
  .ptm_col_med <- apply(ptm_abs_mat, 2, median, na.rm = TRUE)
  .ptm_grd_med <- median(.ptm_col_med, na.rm = TRUE)
  ptm_centering_max_shift <- max(abs(.ptm_grd_med - .ptm_col_med), na.rm = TRUE)
  if (!.ptm_norm) {
    for (.j in seq_len(ncol(ptm_abs_mat)))
      if (!is.na(.ptm_col_med[.j]))
        ptm_abs_mat[, .j] <- ptm_abs_mat[, .j] + (.ptm_grd_med - .ptm_col_med[.j])
  } else {
    logger(paste0("[MQ PTM] Site-matrix median centering NOT applied (Spectronaut PTM export cross-run ",
                  "normalized); it would have shifted samples by up to ", signif(ptm_centering_max_shift, 3), " log2."))
  }

  # ── Build site_to_gene lookup ──
  site_to_gene <- character(0)
  for (site in all_sites) {
    idx <- match(site, coll_df$site_id)
    g   <- if (!is.na(idx)) coll_df$gene_lbl[idx] else ""
    if (is.na(g) || g == "") g <- sub("_[A-Za-z][0-9]+_M[0-9]+$", "", site)
    site_to_gene[site] <- g
  }

  # ── Global protein matrix ──
  # When positional matching is active, cond_groups uses PTM sample names but
  # dia_limma_summarize needs global sample names.  Remap for that call only.
  global_cond_groups <- if (!is.null(ptm_to_global_map)) {
    lapply(cond_groups, function(samps) {
      mapped <- ptm_to_global_map[samps]
      unname(mapped[!is.na(mapped)])
    })
  } else {
    cond_groups
  }
  # Also remap subject_pairing keys if positional matching is active
  global_subject_pairing <- if (!is.null(ptm_to_global_map) && !is.null(subject_pairing)) {
    new_keys  <- ptm_to_global_map[names(subject_pairing)]
    valid     <- !is.na(new_keys)
    setNames(subject_pairing[valid], new_keys[valid])
  } else {
    subject_pairing
  }
  global_cache <- tryCatch(
    dia_limma_summarize(global_parsed, global_cond_groups,
                        subject_pairing = global_subject_pairing,
                        apply_centering = !.glob_norm,
                        logger = logger),
    error = function(e) stop(paste0("Global summarization failed: ", e$message))
  )
  global_mat <- global_cache$protein_matrix
  # Select and align global matrix columns to matched_samples.
  # In the positional-match case, global columns have different names (e.g. wp_01
  # vs pp_01); rename them to PTM sample names so downstream indexing is consistent.
  if (!is.null(ptm_to_global_map)) {
    global_cols_ordered <- ptm_to_global_map[matched_samples]  # PTM->global name
    global_cols_ordered <- global_cols_ordered[global_cols_ordered %in% colnames(global_mat)]
    global_mat <- global_mat[, global_cols_ordered, drop = FALSE]
    colnames(global_mat) <- names(global_cols_ordered)         # rename to PTM names
  } else {
    global_mat <- global_mat[, intersect(colnames(global_mat), matched_samples), drop = FALSE]
  }

  # Build secondary->primary lookup (for non-primary PTM.ProteinId resolution)
  sec2pri <- global_parsed$secondary_to_primary %||% character(0)

  # ── DDA: MinProb imputation on both matrices before stoichiometric subtraction ──
  # Applied independently per matrix using the Perseus-style per-column downshifted
  # Gaussian: values drawn from N(mean(col) - 1.8*sd(col), (0.3*sd(col))^2).
  # This handles MNAR missingness characteristic of DDA data.
  if (identical(acquisition_mode, "dda")) {
    .minprob_impute <- function(mat, label) {
      n_na_before <- sum(is.na(mat))
      if (n_na_before > 0L) {
        set.seed(42L)
        for (j in seq_len(ncol(mat))) {
          col <- mat[, j]
          obs <- col[!is.na(col)]
          if (length(obs) < 2L) next
          col_sd  <- sd(obs)
          col_lod <- mean(obs) - 1.8 * col_sd
          na_idx  <- which(is.na(col))
          mat[na_idx, j] <- rnorm(length(na_idx), mean = col_lod, sd = 0.3 * col_sd)
        }
        logger(paste0("[MQ PTM DDA] MinProb imputation on ", label,
                      ": filled ", n_na_before, " NA values across ",
                      ncol(mat), " samples."))
      } else {
        logger(paste0("[MQ PTM DDA] MinProb imputation on ", label,
                      ": no NAs found, skipped."))
      }
      mat
    }
    ptm_abs_mat <- .minprob_impute(ptm_abs_mat, "PTM matrix")
    global_mat  <- .minprob_impute(global_mat,  "global matrix")
  }

  # ── Stoichiometric adjustment ──
  # Each MaxQuant row is already per-protein (PTM.ProteinId resolved).
  # Extract parent accession from CollapseKey: strip "_AA999_M1" suffix.
  adj_rows         <- list()
  adj_to_abs_map   <- character(0)
  n_shared_expanded <- 0L

  for (site in all_sites) {
    ptm_row    <- ptm_abs_mat[site, ]

    # Extract parent accession from site ID (format: ACC_S941_M1)
    parent_acc <- sub("_[A-Za-z][0-9]+_M[0-9]+$", "", site)

    # Resolve through secondary->primary if needed
    if (!parent_acc %in% rownames(global_mat) && parent_acc %in% names(sec2pri)) {
      parent_acc <- sec2pri[[parent_acc]]
    }

    if (parent_acc %in% rownames(global_mat)) {
      prot_row <- global_mat[parent_acc, ]
      adj_row  <- ptm_row - prot_row
    } else {
      adj_row  <- rep(NA_real_, length(ptm_row))
      names(adj_row) <- names(ptm_row)
    }

    adj_rows[[site]]     <- adj_row
    adj_to_abs_map[site] <- site

    # Mark as shared if PTM.ProteinId in a multi-protein group (detected from CollapseKey
    # matching multiple entries in coll_df with same site but different prot_ids)
    prot_ids_for_site <- unique(coll_df$prot_id[coll_df$site_id == site])
    if (length(prot_ids_for_site) > 1)
      n_shared_expanded <- n_shared_expanded + 1L
  }

  # Build adjusted matrix
  adj_names   <- names(adj_rows)
  ptm_adj_mat <- matrix(NA_real_,
                         nrow = length(adj_names),
                         ncol = length(matched_samples),
                         dimnames = list(adj_names, matched_samples))
  for (nm in adj_names) {
    v <- adj_rows[[nm]]
    ptm_adj_mat[nm, names(v)] <- v
  }

  # Column-wise median centering of the adjusted matrix - DIA (pivot) path only.
  # ptm_abs_mat is centered above on an all-observed basis; global_mat is centered
  # by dia_limma_summarize() on a complete-case basis when enough complete rows
  # exist.  The two offsets do not cancel under subtraction, so the residue enters
  # every adjusted contrast as a per-column constant.  Same construction as the
  # PTM matrix centering above, and as the DIA route's equivalent step.
  #
  # ROUTE GATE: this function serves both maxquant wizard paths.  The DDA path
  # MinProb-imputes both matrices before the subtraction, so it has no NAs and is
  # out of scope; it must remain bit-identical.  Only acquisition_mode == "dia"
  # (the Spectronaut PTM pivot path) is centered here.
  # Also skipped on that path when BOTH exports were cross-run normalized upstream
  # (no mismatched centering bases to reconcile; see dia_ptm_summarize()).
  adj_centering_delta <- NULL
  adj_centering_delta_not_applied <- NULL
  if (identical(acquisition_mode, "dia") && .ptm_norm && .glob_norm) {
    .adj_col_med <- apply(ptm_adj_mat, 2, median, na.rm = TRUE)
    adj_centering_delta_not_applied <- median(.adj_col_med, na.rm = TRUE) - .adj_col_med
    logger(paste0("[MQ PTM] Adjusted matrix centering NOT applied (both Spectronaut exports normalized ",
                  "upstream); it would have shifted samples by up to ",
                  round(max(abs(adj_centering_delta_not_applied), na.rm = TRUE), 3), " log2."))
  } else if (identical(acquisition_mode, "dia")) {
    .adj_col_med <- apply(ptm_adj_mat, 2, median, na.rm = TRUE)
    .adj_grd_med <- median(.adj_col_med, na.rm = TRUE)
    adj_centering_delta <- .adj_grd_med - .adj_col_med
    for (.j in seq_len(ncol(ptm_adj_mat)))
      if (!is.na(.adj_col_med[.j]))
        ptm_adj_mat[, .j] <- ptm_adj_mat[, .j] + (.adj_grd_med - .adj_col_med[.j])
    logger(paste0("[MQ PTM] Adjusted matrix centering: median shift range ",
                  round(min(adj_centering_delta, na.rm = TRUE), 3), " to ",
                  round(max(adj_centering_delta, na.rm = TRUE), 3),
                  " log2 (spread ",
                  round(diff(range(adj_centering_delta, na.rm = TRUE)), 3), ")."))
  }

  logger(paste0("[MQ PTM] Stoichiometry: ", nrow(ptm_adj_mat), " adjusted rows."))

  # ── Design data.frame ──
  all_samples <- intersect(unlist(cond_groups), matched_samples)
  if (!is.null(subject_pairing)) {
    bio_rep_vec <- subject_pairing[all_samples]
  } else {
    bio_rep_vec <- ave(seq_along(all_samples),
                       cond_map[all_samples],
                       FUN = function(i) paste0(cond_map[all_samples][i[1]], "_", seq_along(i)))
  }
  design_df <- data.frame(
    sample_name  = all_samples,
    condition    = cond_map[all_samples],
    bioreplicate = bio_rep_vec,
    stringsAsFactors = FALSE
  )

  # ── Missingness map ──
  conds     <- unique(design_df$condition)
  miss_rows <- lapply(all_sites, function(site) {
    row_mat <- ptm_abs_mat[site, all_samples, drop = FALSE]
    lapply(conds, function(cond) {
      samps_in_cond <- design_df$sample_name[design_df$condition == cond]
      vals_in_cond  <- row_mat[, samps_in_cond, drop = FALSE]
      n_total   <- ncol(vals_in_cond)
      n_detected <- sum(!is.na(vals_in_cond))
      data.frame(Site_ID = site, Condition = cond,
                 Detected = n_detected, Total = n_total,
                 Fraction = if (n_total > 0) n_detected / n_total else 0,
                 stringsAsFactors = FALSE)
    })
  })
  missingness_map <- do.call(rbind, unlist(miss_rows, recursive = FALSE))

  logger(paste0("[MQ PTM] Summarization complete: ", n_sites_output, " sites, ",
                length(matched_samples), " samples, ", length(conds), " conditions."))

  list(
    ptm_abs_matrix    = ptm_abs_mat,
    ptm_adj_matrix    = ptm_adj_mat,
    global_matrix     = global_mat,
    design_df         = design_df,
    cond_groups       = cond_groups,
    acc_to_gene       = global_cache$acc_to_gene,
    site_to_gene      = site_to_gene,
    site_to_prob      = site_to_prob,
    adj_to_abs_map    = adj_to_abs_map,
    missingness_map   = missingness_map,
    n_sites_input     = n_sites_input,
    n_sites_output    = n_sites_output,
    n_sites_adjusted  = nrow(ptm_adj_mat),
    n_shared_expanded = n_shared_expanded,
    adj_centering_delta = adj_centering_delta,   # NULL on the DDA path, and when skipped
    adj_centering_delta_not_applied = adj_centering_delta_not_applied,
    ptm_centering_applied    = !.ptm_norm,
    global_centering_applied = !.glob_norm,
    adj_centering_applied    = identical(acquisition_mode, "dia") && !(.ptm_norm && .glob_norm),
    ptm_centering_max_shift_log2    = ptm_centering_max_shift,
    global_centering_max_shift_log2 = global_cache$normalization_max_shift_log2,
    min_valid_per_cond  = .mq_floor,
    min_valid_per_cond_effective = min_valid_per_cond_effective,
    matched_samples   = matched_samples,
    prob_threshold    = prob_threshold,
    has_probability   = has_prob,
    pipeline_label    = if (identical(acquisition_mode, "dda"))
                          "MaxQuant PTM DDA (MinProb imputation + limma, stoichiometric adjustment)"
                        else
                          "MaxQuant PTM DIA (limma, intensity summation + stoichiometric adjustment)",
    software          = "maxquant",
    acquisition_mode  = acquisition_mode,
    normalization     = "raw_intensity_summation_log2_stoichiometric_adjustment"
  )
}

# =============================================================================
# tmt_ptm_limma_summarize()
# Summarizes paired TMT phosphoproteomics (Proteome Discoverer nested CSV)
# into the dia_ptm_master_cache format so all downstream PTM modules can
# consume the result without modification.
#
# Inputs:
#   ptm_parsed_files   -" list of parsed phospho-enriched PD file structs
#   prot_parsed_files  -" list of parsed global profiling PD file structs
#   annotation_map     -" data.frame: ColumnName, Channel, Condition,
#                         BioReplicate, Mixture, Run
#                       Norm rows = bridge channels (excluded from
#                       DE but used for IRS normalization when design_opt=2)
#   target_ptms        -" character vector of PTM types to include (e.g. "Phospho")
#   prob_threshold     -" ptmRS/phosphoRS localization probability cutoff (%)
#   design_opt         -" integer: 1 = single-plex (median centering),
#                         2 = multi-plex IRS bridge, 3 = multi-plex run-covariate
#   subject_pairing    -" NULL or named character: sample_name -> subject_id
#   rv                 -" list(pd_has_gn, pd_gene_col) for gene name resolution
#   logger             -" logging function
#
# Returns: list in the same format as dia_ptm_summarize() / dia_ptm_master_cache
# =============================================================================
tmt_ptm_limma_summarize <- function(ptm_parsed_files,
                                     prot_parsed_files,
                                     annotation_map,
                                     target_ptms    = "Phospho",
                                     prob_threshold = 75,
                                     pep_qvalue_max = 0.01,
                                     prot_qvalue_max = 0.01,
                                     interference_threshold = 50,
                                     min_valid_per_cond = 2L,
                                     # Plex-coverage fractions. Default 0, which disables the
                                     # filter: max(1L, ceiling(0 * n_plexes)) is 1, i.e. "present
                                     # in at least one plex", which every row already satisfies.
                                     # This matches MSstatsTMT, which applies no run-coverage
                                     # filter and instead reduces the model for a feature measured
                                     # in a single run. Retention is governed by the per-condition
                                     # valid-value rule and, at testing, by the both-arms rule.
                                     # A site living entirely within one plex has both compared
                                     # conditions inside the same batch and is the LEAST exposed to
                                     # between-plex variation, so deleting it inverted the intent.
                                     # Design 2 no longer forces 1.0: IRS marks the plexes it
                                     # cannot correct instead of the caller deleting the site.
                                     min_plex_frac  = 0,
                                     prot_min_plex_frac = 0,
                                     design_opt     = 1L,
                                     subject_pairing = NULL,
                                     rv             = list(),
                                     logger         = message) {
  logger(paste0("[TMT PTM] Starting summarization (", length(ptm_parsed_files),
                " phospho plex(es), prob_threshold=", prob_threshold, "%)..."))

  # ── A. Identify biological channels ──────────────────────────────────────
  bio_ann <- annotation_map[
    !is.na(annotation_map$Condition) &
    !annotation_map$Condition %in% c("Norm", ""),
    , drop = FALSE
  ]
  target_cols <- bio_ann$ColumnName    # biological channel column names
  all_cols    <- annotation_map$ColumnName  # all columns incl. bridge

  if (length(target_cols) == 0)
    return(list(error = "No biological channels found in annotation_map."))

  tgt_pattern <- paste(target_ptms, collapse = "|")

  # ── B. Collect site-level observations from phospho (peptide-level) files ──
  # For each peptide row: parse site IDs from "Modifications in Master Proteins",
  # collect raw abundances per channel, accumulate into a flat data.frame.
  acc_parts <- list()

  interference_col_found <- FALSE
  n_interf_kept_na       <- 0L
  n_interf_kept_pass     <- 0L
  n_interf_removed       <- 0L
  n_razor_reassigned     <- 0L
  # Peptide-group identification-confidence filter, accumulated across files so
  # the methods text can state the basis used and the retention, rather than
  # leaving both only in the log. Recording only; the filter itself is unchanged.
  n_pep_groups_input     <- 0L
  n_pep_groups_conf_kept <- 0L
  pep_confidence_basis   <- "none"

  # Razor support map: accession -> total PSMs, summed across EVERY phospho file.
  # Built once so the razor winner for a given set of master proteins is the same
  # in every plex. Built per file, a shared site can be assigned to one protein in
  # one plex and another in the next, splitting one site into two half-populated
  # rows that each look single-plex.
  .razor_psm <- local({
    .acc_tot <- c()
    for (.f in ptm_parsed_files) {
      .pd <- tryCatch(parse_pd_protein_data(.f, rv), error = function(e) NULL)
      if (is.null(.pd) || is.null(.pd$protein_df)) next
      .acc_c <- grep("^Accession$", colnames(.pd$protein_df),
                     value = TRUE, ignore.case = TRUE)[1]
      .psm_c <- grep("^# PSMs$", colnames(.pd$protein_df), value = TRUE)[1]
      if (is.na(.acc_c) || is.na(.psm_c)) next
      .a <- trimws(sub(";.*", "", as.character(.pd$protein_df[[.acc_c]])))
      .n <- suppressWarnings(as.numeric(.pd$protein_df[[.psm_c]]))
      .ok <- !is.na(.a) & nzchar(.a) & !is.na(.n)
      if (!any(.ok)) next
      .agg <- tapply(.n[.ok], .a[.ok], sum, na.rm = TRUE)
      for (.k in names(.agg))
        .acc_tot[.k] <- (if (.k %in% names(.acc_tot)) .acc_tot[[.k]] else 0) + .agg[[.k]]
    }
    if (length(.acc_tot) == 0L) NULL else .acc_tot
  })
  if (is.null(.razor_psm))
    logger(paste0("[TMT PTM] No protein-level PSM counts available in any phospho ",
                  "file; shared peptides will be assigned to the first listed ",
                  "master protein."))
  else
    logger(paste0("[TMT PTM] Razor support map built across ",
                  length(ptm_parsed_files), " phospho file(s): ",
                  length(.razor_psm), " accessions."))

  for (fi in seq_along(ptm_parsed_files)) {
    f <- ptm_parsed_files[[fi]]
    pep_df <- tryCatch(parse_pd_peptide_data(f), error = function(e) NULL)
    if (is.null(pep_df) || nrow(pep_df) == 0) {
      logger(paste0("[TMT PTM] Phospho file ", fi, ": no peptide data, skipping."))
      next
    }

    # Peptide-level identification confidence. Qvality q-value where present,
    # else the Confidence tier. Absent columns leave the file unfiltered and are
    # reported, never silently treated as passing.
    .n_pep_before <- nrow(pep_df)
    n_pep_groups_input <- n_pep_groups_input + .n_pep_before
    .q_col <- grep("^Qvality q-value$", colnames(pep_df), value = TRUE)[1]
    .c_col <- grep("^Confidence$",      colnames(pep_df), value = TRUE)[1]
    if (!is.na(.q_col)) {
      .qv    <- suppressWarnings(as.numeric(pep_df[[.q_col]]))
      pep_df <- pep_df[!is.na(.qv) & .qv < pep_qvalue_max, , drop = FALSE]
      pep_confidence_basis <- "qvalue"
      logger(paste0("[TMT PTM] Phospho file ", fi, ": peptide q-value < ",
                    pep_qvalue_max, " retained ", nrow(pep_df), " of ",
                    .n_pep_before, " peptide groups."))
    } else if (!is.na(.c_col)) {
      pep_df <- pep_df[toupper(trimws(as.character(pep_df[[.c_col]]))) == "HIGH", ,
                       drop = FALSE]
      # A file with no q-value column does not downgrade a q-value basis already
      # established by an earlier file; "mixed" is recorded so the methods text
      # never claims a single basis that was not applied throughout.
      pep_confidence_basis <- if (identical(pep_confidence_basis, "qvalue"))
        "mixed" else "confidence"
      logger(paste0("[TMT PTM] Phospho file ", fi, ": no Qvality q-value column; ",
                    "Confidence = High retained ", nrow(pep_df), " of ",
                    .n_pep_before, " peptide groups."))
    } else {
      if (!identical(pep_confidence_basis, "none")) pep_confidence_basis <- "mixed"
      logger(paste0("[TMT PTM] Phospho file ", fi, ": no peptide confidence column ",
                    "found; no peptide-level confidence filter applied."))
    }
    n_pep_groups_conf_kept <- n_pep_groups_conf_kept + nrow(pep_df)
    if (nrow(pep_df) == 0L) {
      logger(paste0("[TMT PTM] Phospho file ", fi,
                    ": no peptide groups passed the confidence filter, skipping."))
      next
    }

    # Co-isolation interference filter, mirroring tmt_limma_qf_summarize().
    # Missing values are kept and counted separately from passing values: the
    # two are different claims and the methods text must be able to state both.
    .int_col <- grep("Isolation.Interference|Co.isolation", colnames(pep_df),
                     value = TRUE, ignore.case = TRUE)[1]
    if (!is.na(.int_col)) {
      interference_col_found <- TRUE
      .iv        <- suppressWarnings(as.numeric(pep_df[[.int_col]]))
      .keep_int  <- is.na(.iv) | .iv <= interference_threshold
      n_interf_kept_na   <- n_interf_kept_na +
                            sum(is.na(.iv))
      n_interf_kept_pass <- n_interf_kept_pass +
                            sum(!is.na(.iv) & .iv <= interference_threshold)
      n_interf_removed   <- n_interf_removed + sum(!.keep_int)
      .n_before_int <- nrow(pep_df)
      pep_df        <- pep_df[.keep_int, , drop = FALSE]
      logger(paste0("[TMT PTM] Phospho file ", fi, ": interference filter (<= ",
                    interference_threshold, "%): ", .n_before_int, " -> ",
                    nrow(pep_df), " peptide groups."))
      if (nrow(pep_df) == 0L) {
        logger(paste0("[TMT PTM] Phospho file ", fi,
                      ": no peptide groups passed the interference filter, skipping."))
        next
      }
    }

    mmp_col <- grep("^Modifications in Master Proteins$",
                    colnames(pep_df), value = TRUE)[1]
    if (is.na(mmp_col)) {
      logger(paste0("[TMT PTM] Phospho file ", fi,
                    ": 'Modifications in Master Proteins' column not found, skipping."))
      next
    }

    # Columns in this file that match annotation_map
    file_abund_cols <- intersect(all_cols, colnames(pep_df))
    if (length(file_abund_cols) == 0) {
      logger(paste0("[TMT PTM] Phospho file ", fi, ": no matching abundance columns."))
      next
    }

    # Cache parsed MMP strings to avoid re-parsing duplicates
    mmp_vals  <- pep_df[[mmp_col]]
    uniq_mmps <- unique(mmp_vals[!is.na(mmp_vals) & nchar(trimws(mmp_vals)) > 0])
    # Localization-probability scale is a property of the file, decided once over
    # every parenthesised value in the column, not per modification string.
    .all_probs <- suppressWarnings(as.numeric(unlist(regmatches(
      uniq_mmps, gregexpr("(?<=\\()[0-9.]+(?=\\))", uniq_mmps, perl = TRUE)))))
    .all_probs <- .all_probs[!is.na(.all_probs) & .all_probs > 0]
    file_prob_scale <- if (length(.all_probs) == 0L) "percent"
                       else if (all(.all_probs <= 1)) "fraction" else "percent"
    logger(paste0("[TMT PTM] Phospho file ", fi, ": localization probability scale ",
                  "detected as ", file_prob_scale, " (", length(.all_probs),
                  " values inspected)."))
    mmp_cache <- lapply(uniq_mmps, function(v)
      parse_ptm_mmp_sites(v, tgt_pattern, absent_prob = NA_real_,
                          prob_scale = file_prob_scale))
    names(mmp_cache) <- uniq_mmps

    # Peptidoform identity columns, resolved once per file rather than per row.
    .seq_col <- grep("^Annotated Sequence$", colnames(pep_df), value = TRUE)[1]
    .mod_col <- grep("^Modifications$",      colnames(pep_df), value = TRUE)[1]

    for (r in seq_len(nrow(pep_df))) {
      mv <- mmp_vals[r]
      if (is.na(mv) || nchar(trimws(mv)) == 0) next
      sites_df <- mmp_cache[[mv]]
      if (is.null(sites_df) || nrow(sites_df) == 0) next

      # Collect intensities for this peptide row across all matching channels
      row_vals <- suppressWarnings(
        as.numeric(unlist(pep_df[r, file_abund_cols, drop = FALSE])))
      valid_mask <- !is.na(row_vals) & row_vals > 0
      if (!any(valid_mask)) next
      valid_cols <- file_abund_cols[valid_mask]
      valid_vals <- row_vals[valid_mask]

      # Peptidoform identity for this peptide-group row: annotated sequence plus
      # its full modification string. Used to count distinct peptidoforms per
      # site; falls back to the row index when the columns are absent.
      # .seq_col / .mod_col are resolved once per file, above the row loop.
      .row_peptidoform <- paste0(
        if (!is.na(.seq_col)) as.character(pep_df[[.seq_col]][r]) else "",
        "|",
        if (!is.na(.mod_col)) as.character(pep_df[[.mod_col]][r]) else "",
        "|", fi, ":", r)

      # Razor: a peptide mapping to several master proteins is assigned to the
      # one with the most supporting PSMs. Ties and unknown counts fall back to
      # the first listed accession, matching PD's own ordering.
      .accs_here <- unique(sites_df$Accession[!is.na(sites_df$Accession) &
                                                nchar(sites_df$Accession) > 0])
      .is_shared_pep <- length(.accs_here) > 1L
      .alt_accs      <- ""
      .all_accs      <- paste(.accs_here, collapse = ";")
      if (.is_shared_pep) {
        .sup <- if (!is.null(.razor_psm)) .razor_psm[.accs_here] else
                rep(NA_real_, length(.accs_here))
        .win <- if (all(is.na(.sup))) .accs_here[1] else
                .accs_here[which.max(replace(.sup, is.na(.sup), -Inf))]
        .alt_accs  <- paste(setdiff(.accs_here, .win), collapse = ";")
        # Razor winner first, then the alternates. The winner is the key and the
        # protein the stoichiometric adjustment used; the rest are equally valid
        # mappings of the same measurement and are kept for labelling and search.
        .all_accs  <- paste(c(.win, setdiff(.accs_here, .win)), collapse = ";")
        sites_df   <- sites_df[sites_df$Accession == .win, , drop = FALSE]
        n_razor_reassigned <- n_razor_reassigned + 1L
      }
      if (nrow(sites_df) == 0L) next

      for (si in seq_len(nrow(sites_df))) {
        acc  <- sites_df$Accession[si]
        site <- sites_df$Site[si]
        mod  <- sites_df$Modification[si]
        prob <- sites_df$Localization_Prob[si]
        if (site == "Ambiguous" || is.na(acc) || nchar(acc) == 0) next

        site_id <- paste0(acc, "_", site, "_", mod)
        acc_parts[[length(acc_parts) + 1L]] <- data.frame(
          site_id     = site_id,
          col_name    = valid_cols,
          intensity   = valid_vals,
          prob        = prob,
          peptidoform = .row_peptidoform,
          shared_peptide = .is_shared_pep,
          alt_accessions = .alt_accs,
          all_accessions = .all_accs,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!interference_col_found)
    logger(paste0("[TMT PTM] Isolation Interference column not found in any ",
                  "phospho file - proceeding without the co-isolation filter. ",
                  "WARNING: co-isolation interference may compress fold changes. ",
                  "The thresholds applied inside Proteome Discoverer are not ",
                  "recorded in this export."))

  if (length(acc_parts) == 0)
    return(list(error = paste0(
      "No PTM sites detected in phospho files. ",
      "Check that 'Modifications in Master Proteins' column is present and contains ",
      "the selected PTM type (", paste(target_ptms, collapse = "/"), ").")))

  coll_df <- do.call(rbind, acc_parts)
  n_sites_input <- length(unique(coll_df$site_id))
  n_obs_no_prob <- sum(is.na(coll_df$prob))
  if (n_obs_no_prob == nrow(coll_df))
    return(list(error = paste0(
      "No ptmRS/phosphoRS localization probabilities were found in the ",
      "'Modifications in Master Proteins' column. Enable the ptmRS (or IMP-ptmRS) ",
      "node in the Proteome Discoverer consensus workflow and re-export. ",
      "Lowering the probability threshold will not help.")))
  logger(paste0("[TMT PTM] Raw collection: ", nrow(coll_df),
                " observations, ", n_sites_input, " unique sites."))

  # ---- C. Localization threshold - site-level filter --------------------------
  # Site-level localization decision (class-I convention): a site qualifies on its
  # maximum probability across all observations, and once it qualifies its
  # intensities are taken from every channel in which it was measured.
  .site_max_prob <- tapply(coll_df$prob, coll_df$site_id, function(p) {
    p <- p[!is.na(p)]
    if (length(p) == 0L) NA_real_ else max(p)
  })
  passing_sites <- names(.site_max_prob)[
    !is.na(.site_max_prob) & .site_max_prob >= prob_threshold]

  n_obs_before_loc <- nrow(coll_df)
  coll_df          <- coll_df[coll_df$site_id %in% passing_sites, , drop = FALSE]
  n_obs_below_own_threshold <- sum(!is.na(coll_df$prob) &
                                     coll_df$prob < prob_threshold)
  n_obs_no_own_prob         <- sum(is.na(coll_df$prob))

  n_sites_pass <- length(passing_sites)
  site_to_prob <- .site_max_prob[passing_sites]

  # Computed after the localization filter so it cannot name sites absent from
  # the output matrix.
  shared_pep_sites <- unique(coll_df$site_id[which(as.logical(coll_df$shared_peptide))])

  # Site -> every master protein the peptide mapped to, razor winner first.
  # One row per site is kept (the razor key); the full list is retained for
  # labelling and search, mirroring MaxQuant's Phospho (STY)Sites convention of
  # one site row carrying a semicolon-separated Proteins list plus a designated
  # leading protein.
  # Union across every observation of the site, not the first one. Peptidoforms
  # of the same residue can carry different accession sets while resolving to the
  # same razor winner, in which case taking v[1] makes the alternates depend on
  # row order. Order within the union follows first appearance, so the razor
  # winner stays first in all_accessions.
  .aa_union <- tapply(coll_df$all_accessions, coll_df$site_id,
                      function(v) { v <- unique(unlist(strsplit(v[nzchar(v)], ";",
                                                                fixed = TRUE)))
                                    if (length(v) == 0L) NA_character_
                                    else paste(v, collapse = ";") })
  site_all_accessions <- setNames(as.character(.aa_union), names(.aa_union))
  .alt_union <- tapply(coll_df$alt_accessions, coll_df$site_id,
                       function(v) { v <- unique(unlist(strsplit(v[nzchar(v)], ";",
                                                                 fixed = TRUE)))
                                     if (length(v) == 0L) ""
                                     else paste(v, collapse = ";") })
  site_alt_accessions <- setNames(as.character(.alt_union), names(.alt_union))

  logger(paste0("[TMT PTM] Localization filter (site max >=", prob_threshold, "%): ",
                n_sites_input, " -> ", n_sites_pass, " sites; ",
                nrow(coll_df), " of ", n_obs_before_loc, " observations retained (",
                n_obs_below_own_threshold, " below the threshold individually, ",
                n_obs_no_own_prob, " with no recorded probability)."))

  if (n_sites_pass == 0)
    return(list(error = paste0(
      "No sites passed probability threshold (", prob_threshold, "%). ",
      "Try lowering the threshold or verify ptmRS probabilities in your data.")))

  # Loading constants are derived from the peptide-level observations that passed
  # the confidence and localization filters, before rollup and before the
  # valid-value filter thins the matrix. Biological channels only.
  .obs_bio <- coll_df[coll_df$col_name %in% target_cols, , drop = FALSE]
  pre_rollup_constants <- if (nrow(.obs_bio) > 0L) {
    .cm <- tapply(log2(.obs_bio$intensity), .obs_bio$col_name,
                  function(v) median(v[is.finite(v)], na.rm = TRUE))
    .gm <- median(.cm, na.rm = TRUE)
    .gm - .cm
  } else NULL

  # ---- D. Rollup: peptidoform summation per (site, channel), then log2 --------
  # Every row in coll_df belongs to a site whose maximum localization probability
  # met the threshold; individual rows below it are retained.  Among those
  # candidates, every peptidoform intensity mapping to the same (site, channel)
  # pair is summed on the linear scale and the sum is log2-transformed.  No
  # peptidoform is selected over another and none is discarded.  When only one
  # peptidoform exists for a site-channel pair the sum is that single value.
  logger("[TMT PTM DEBUG] Phase D: rollup by linear-scale summation across peptidoforms")
  # Strip embedded null bytes that can appear in PD CSV string columns.
  # This MUST run before all_sites is derived. It used to run after, so the matrix
  # row names were the UNSTRIPPED site ids while coll_key was built from stripped
  # ones - match(site_id, all_sites) would then miss and those intensities would be
  # dropped silently. n_peptidoforms_by_site, computed further down from the same
  # (mutated) column, had the same mismatch. Latent: it fires only on an export
  # that actually carries embedded nulls.
  .nul_chr <- suppressWarnings(rawToChar(as.raw(0)))
  if (nchar(.nul_chr, type = "bytes") > 0L) {
    coll_df$site_id  <- gsub(.nul_chr, "", coll_df$site_id,  fixed = TRUE)
    coll_df$col_name <- gsub(.nul_chr, "", coll_df$col_name, fixed = TRUE)
  }
  all_sites   <- sort(unique(coll_df$site_id))
  ptm_abs_mat <- matrix(NA_real_,
    nrow = length(all_sites), ncol = length(all_cols),
    dimnames = list(all_sites, all_cols))

  coll_key <- paste(coll_df$site_id, coll_df$col_name, sep = "||")

  # Sum all peptidoform intensities mapping to a site within a channel, on the
  # linear scale, then log2. Distinct peptidoforms are summed together
  # (multiplicity classes are collapsed, not split).
  sum_per_key <- tapply(coll_df$intensity, coll_key, function(v) sum(v, na.rm = TRUE))

  for (nm in names(sum_per_key)) {
    int_val <- sum_per_key[[nm]]
    parts   <- strsplit(nm, "||", fixed = TRUE)[[1]]
    site_id <- parts[1]; col_nm <- parts[2]
    ri <- match(site_id, all_sites)
    ci <- match(col_nm, all_cols)
    if (!is.na(ri) && !is.na(ci) && is.finite(int_val) && int_val > 0)
      ptm_abs_mat[ri, ci] <- log2(int_val)
  }

  # Peptidoform support per site, recorded for the results table.
  n_peptidoforms_by_site <- tapply(coll_df$peptidoform, coll_df$site_id,
                                   function(v) length(unique(v)))
  # Scalar for the methods text: how many sites rest on more than one distinct
  # peptidoform, and are therefore summed across peptidoforms rather than being a
  # single measurement. Computed before the missingness filters, so A3 recomputes
  # it against the surviving rows at the point of return.
  n_sites_multi_peptidoform_all <- sum(n_peptidoforms_by_site > 1L, na.rm = TRUE)
  ptm_abs_mat[!is.finite(ptm_abs_mat)] <- NA

  # Retain a site quantified in at least min_valid_per_cond channels of at least
  # one condition. Deliberately looser than the per-contrast DE filter.
  ptm_abs_bio <- ptm_abs_mat[, target_cols, drop = FALSE]
  .cond_of_col <- setNames(bio_ann$Condition, bio_ann$ColumnName)[colnames(ptm_abs_bio)]
  .keep_site <- apply(ptm_abs_bio, 1, function(x) {
    .n <- tapply(!is.na(x), .cond_of_col, sum)
    any(.n >= min_valid_per_cond, na.rm = TRUE)
  })
  ptm_abs_bio <- ptm_abs_bio[.keep_site, , drop = FALSE]
  # Recorded before the plex-coverage filter runs, so the two filters can be
  # reported as separate claims instead of one combined drop.
  n_sites_after_valid_filter <- nrow(ptm_abs_bio)

  # Plex coverage: TMT missingness is block-structured, so a site quantified in
  # too few plexes enters with whole plexes empty. Same rule as the protein matrix.
  .plex_ids <- unique(annotation_map$Mixture[!is.na(annotation_map$Mixture)])
  n_sites_dropped_plex <- 0L
  # Realised coverage threshold, for the methods text. NA when the filter does
  # not run (fewer than two plexes, or an empty matrix).
  n_plexes_detected <- length(.plex_ids)
  site_min_plex_n   <- NA_integer_
  # One rule for every design. IRS no longer needs the caller to pre-filter to
  # complete cases: under the "na" policy it marks the site-by-plex blocks it
  # cannot align, so the site survives on the plexes where it can be corrected.
  # At the default fraction of 0 this filter is inert.
  site_min_plex_frac <- min_plex_frac
  if (length(.plex_ids) >= 2L && nrow(ptm_abs_bio) > 0L) {
    .present <- vapply(.plex_ids, function(p) {
      .pc <- intersect(annotation_map$ColumnName[!is.na(annotation_map$Mixture) &
                                                   annotation_map$Mixture == p],
                       colnames(ptm_abs_bio))
      if (length(.pc) == 0L) rep(FALSE, nrow(ptm_abs_bio)) else
        rowSums(!is.na(ptm_abs_bio[, .pc, drop = FALSE])) > 0L
    }, logical(nrow(ptm_abs_bio)))
    if (is.null(dim(.present))) .present <- matrix(.present, nrow = nrow(ptm_abs_bio))
    .min_plex <- max(1L, ceiling(site_min_plex_frac * length(.plex_ids)))
    site_min_plex_n <- as.integer(.min_plex)
    .keep_cov <- rowSums(.present) >= .min_plex
    n_sites_dropped_plex <- sum(!.keep_cov)
    if (n_sites_dropped_plex > 0L) {
      ptm_abs_bio <- ptm_abs_bio[.keep_cov, , drop = FALSE]
      logger(paste0("[TMT PTM] Plex coverage: ", n_sites_dropped_plex,
                    " sites removed (quantified in fewer than ", .min_plex,
                    " of ", length(.plex_ids), " plexes; completeness fraction ",
                    site_min_plex_frac, "). A fraction of 0 disables this filter."))
    }
  }
  ptm_abs_mat <- ptm_abs_mat[rownames(ptm_abs_bio), , drop = FALSE]  # keep all cols for IRS
  all_sites   <- rownames(ptm_abs_mat)
  n_sites_output <- nrow(ptm_abs_mat)
  logger(paste0("[TMT PTM] After missingness filter (>=", min_valid_per_cond,
                " valid values in at least one condition): ", n_sites_output,
                " sites retained."))

  if (n_sites_output == 0)
    return(list(error = paste0(
      "No sites retained after missingness filter (require >=", min_valid_per_cond,
      " non-NA values in at least one condition).")))

  # ── E. Global protein matrix from protein-level data ─────────────────────
  logger("[TMT PTM DEBUG] Phase E: protein matrix from global files")
  glob_result <- tmt_limma_protein_summarize(
    files           = prot_parsed_files,
    annotation_map  = annotation_map,
    rv              = rv,
    design_opt      = design_opt,
    logger          = logger,
    phospho_files   = ptm_parsed_files,  # anchors per-plex FN-normalisation
    min_plex_frac   = prot_min_plex_frac,
    # The protein matrix must follow the same IRS policy as the site matrix, or
    # the stoichiometric subtraction mixes corrected site values against
    # uncorrected protein values.
    irs_uncorrectable = if (identical(as.integer(design_opt), 2L)) "na" else "keep",
    prot_qvalue_max = prot_qvalue_max,
    norm_channels   = target_cols
  )
  if (!is.null(glob_result$error))
    return(list(error = paste0("Global protein summarization failed: ", glob_result$error)))

  global_mat  <- glob_result$protein_matrix   # log2, normalized, proteins x all_cols
  acc_to_gene <- glob_result$gene_names        # named char: accession -> gene

  # ── F. Normalize PTM abs matrix (IRS or median centering) ─────────────────
  logger("[TMT PTM DEBUG] Phase F: normalize PTM matrix")
  # Loading correction first, then IRS where a bridge design is in use. IRS
  # applies one offset per feature per plex and cannot correct within-plex
  # channel loading, so it never replaces the centering step. Constants are
  # derived from biological channels only; bridge channels stay in the
  # matrix for IRS but do not set the reference level.
  .norm_cols <- intersect(target_cols, colnames(ptm_abs_mat))
  if (length(.norm_cols) == 0L) .norm_cols <- colnames(ptm_abs_mat)
  sm      <- apply(ptm_abs_mat[, .norm_cols, drop = FALSE], 2, median, na.rm = TRUE)
  gm      <- median(sm, na.rm = TRUE)
  .all_sm <- apply(ptm_abs_mat, 2, median, na.rm = TRUE)
  # T9: prefer the pre-rollup peptide-observation constants where the channel has
  # an entry; fall back to the matrix-derived shift otherwise.
  logger(paste0("[TMT PTM] Centering basis: ",
                if (!is.null(pre_rollup_constants)) "pre-rollup peptide observations"
                else "site matrix", "."))
  normalization_basis <- if (!is.null(pre_rollup_constants))
    "pre-rollup peptide observations" else "site matrix"

  # Per-channel shift, resolved once so that every column ends on ONE reference
  # level. The pre-rollup constants are derived from individual peptide log2
  # intensities and the matrix medians from log2 summed-site intensities; the two
  # populations sit on different scales, so a channel taking one constant and a
  # channel taking the other would end up on different levels. gm_target is the
  # level the biological channels actually reach after their own shift, and
  # channels with no pre-rollup entry (bridge) are brought to that level.
  .shift_of <- setNames(rep(NA_real_, ncol(ptm_abs_mat)), colnames(ptm_abs_mat))
  for (.cn in .norm_cols) {
    .pre <- if (!is.null(pre_rollup_constants) &&
                .cn %in% names(pre_rollup_constants))
              pre_rollup_constants[[.cn]] else NA_real_
    .shift_of[[.cn]] <- if (!is.na(.pre)) .pre
                        else if (.cn %in% names(sm)) gm - sm[[.cn]]
                        else NA_real_
  }
  .bio_post <- vapply(.norm_cols, function(.cn) {
    if (is.na(.shift_of[[.cn]]) || is.na(.all_sm[[.cn]])) NA_real_
    else unname(.all_sm[[.cn]] + .shift_of[[.cn]])
  }, numeric(1))
  gm_target <- median(.bio_post, na.rm = TRUE)
  if (!is.finite(gm_target)) gm_target <- gm
  for (.cn in setdiff(colnames(ptm_abs_mat), .norm_cols))
    if (!is.na(.all_sm[[.cn]])) .shift_of[[.cn]] <- gm_target - .all_sm[[.cn]]

  for (j in seq_len(ncol(ptm_abs_mat))) {
    .shift <- .shift_of[[colnames(ptm_abs_mat)[j]]]
    if (!is.na(.shift)) ptm_abs_mat[, j] <- ptm_abs_mat[, j] + .shift
  }
  site_irs_stats <- NULL
  if (identical(as.integer(design_opt), 2L)) {
    ptm_abs_mat <- apply_irs_normalization(ptm_abs_mat, annotation_map, logger,
                                           uncorrectable = "na")
    # Read here, before the column subsetting below: attributes do not survive it.
    site_irs_stats <- attr(ptm_abs_mat, "irs_stats")
  }

  # Restrict both matrices to biological channels
  shared_bio_cols  <- intersect(target_cols, colnames(ptm_abs_mat))
  shared_glob_cols <- intersect(shared_bio_cols, colnames(global_mat))
  matched_samples  <- shared_glob_cols

  ptm_abs_mat <- ptm_abs_mat[, matched_samples, drop = FALSE]
  global_mat  <- global_mat[, matched_samples, drop = FALSE]

  if (length(matched_samples) == 0)
    return(list(error = "No overlapping biological channels between PTM and global matrices."))

  # Post-IRS coverage re-check (design 2 only; a no-op elsewhere, and gated so the
  # other designs are provably untouched by this edit). IRS marks a site's values
  # missing in any plex where it had no bridge value, and that happens after the
  # summarization missingness filter has already run. Re-apply the same rule here
  # so n_sites_output remains a count of sites that are actually usable, and record
  # the additional drop as its own number rather than folding it into the first.
  # ptm_abs_mat's columns are exactly matched_samples at this point.
  n_sites_dropped_post_irs <- 0L
  if (identical(as.integer(design_opt), 2L) && nrow(ptm_abs_mat) > 0L) {
    .cond_post <- setNames(bio_ann$Condition, bio_ann$ColumnName)[matched_samples]
    .keep_post <- apply(ptm_abs_mat, 1, function(x) {
      .n <- tapply(!is.na(x), .cond_post, sum)
      any(.n >= min_valid_per_cond, na.rm = TRUE)
    })
    n_sites_dropped_post_irs <- sum(!.keep_post)
    if (n_sites_dropped_post_irs > 0L) {
      ptm_abs_mat    <- ptm_abs_mat[.keep_post, , drop = FALSE]
      all_sites      <- rownames(ptm_abs_mat)
      n_sites_output <- nrow(ptm_abs_mat)
      logger(paste0("[TMT PTM] Post-IRS coverage: ", n_sites_dropped_post_irs,
                    " site(s) fell below >=", min_valid_per_cond,
                    " valid values in every condition once uncorrectable plex blocks ",
                    "were set missing; ", n_sites_output, " sites retained."))
    }
    if (n_sites_output == 0)
      return(list(error = paste0(
        "No sites retained after the post-IRS coverage re-check: every site lost too ",
        "many channels when plexes with no bridge value were set missing. Check that ",
        "each plex carries a bridge (Norm) channel.")))
  }

  # ── G. Site-to-gene lookup ────────────────────────────────────────────────
  logger("[TMT PTM DEBUG] Phase G: site-to-gene lookup")
  site_to_gene <- setNames(character(length(all_sites)), all_sites)
  .gene_of_acc <- function(a) {
    g <- if (!is.null(acc_to_gene) && a %in% names(acc_to_gene))
      acc_to_gene[[a]] else a
    if (!is.na(g) && nzchar(g)) g else a
  }
  for (site in all_sites) {
    # Parse parent accession: everything before _<Residue><Pos>_<Mod>
    acc_parent <- sub("_[A-Za-z][0-9]+_.*$", "", site)
    # Every master protein this site's peptide mapped to, razor winner first.
    # Semicolon-joined: consumers take the first element for display and the
    # whole string is searchable. Falls back to the parent accession alone when
    # the site carries no recorded accession list.
    .accs <- if (!is.null(site_all_accessions) &&
                 site %in% names(site_all_accessions) &&
                 !is.na(site_all_accessions[[site]]) &&
                 nzchar(site_all_accessions[[site]]))
      strsplit(site_all_accessions[[site]], ";", fixed = TRUE)[[1]]
    else acc_parent
    .accs <- unique(c(acc_parent, .accs))
    site_to_gene[[site]] <- paste(unique(vapply(.accs, .gene_of_acc, character(1))),
                                  collapse = ";")
  }

  # ── H. Stoichiometric adjustment + shared peptide expansion ──────────────
  logger("[TMT PTM DEBUG] Phase H: stoichiometric adjustment")
  adj_rows       <- list()
  adj_to_abs_map <- character(0)
  n_shared_expanded <- 0L
  n_no_parent_protein <- 0L

  for (site in all_sites) {
    ptm_row    <- ptm_abs_mat[site, matched_samples]
    parent_acc <- sub("_[A-Za-z][0-9]+_.*$", "", site)
    if (parent_acc %in% rownames(global_mat)) {
      prot_row <- global_mat[parent_acc, matched_samples]
      adj_row  <- ptm_row - prot_row
    } else {
      adj_row  <- setNames(rep(NA_real_, length(matched_samples)), matched_samples)
      n_no_parent_protein <- n_no_parent_protein + 1L
    }
    adj_rows[[site]]     <- adj_row
    adj_to_abs_map[site] <- site
  }

  adj_names   <- names(adj_rows)
  ptm_adj_mat <- matrix(NA_real_,
    nrow = length(adj_names), ncol = length(matched_samples),
    dimnames = list(adj_names, matched_samples))
  for (nm in adj_names) {
    v <- adj_rows[[nm]]
    cm <- intersect(names(v), matched_samples)
    ptm_adj_mat[nm, cm] <- v[cm]
  }

  # The PTM and protein matrices were centred independently, so their per-column
  # offsets do not cancel under subtraction. Remove the residual offset.
  .adj_med   <- apply(ptm_adj_mat, 2, median, na.rm = TRUE)
  .adj_grand <- median(.adj_med, na.rm = TRUE)
  adj_centering_delta <- .adj_grand - .adj_med
  for (j in seq_len(ncol(ptm_adj_mat)))
    if (!is.na(adj_centering_delta[j]))
      ptm_adj_mat[, j] <- ptm_adj_mat[, j] + adj_centering_delta[j]
  logger(paste0("[TMT PTM] Adjusted matrix re-centred; shift range ",
                signif(diff(range(adj_centering_delta, na.rm = TRUE)), 3),
                " log2 units."))

  logger(paste0("[TMT PTM] Stoichiometry: ", nrow(ptm_adj_mat), " adjusted rows (",
                n_shared_expanded, " from shared peptide expansion)."))

  # ── I. Design df (biological samples only) ───────────────────────────────
  logger("[TMT PTM DEBUG] Phase I: build design df")
  bio_cols_in_map <- bio_ann[bio_ann$ColumnName %in% matched_samples, , drop = FALSE]
  cond_map        <- setNames(bio_cols_in_map$Condition, bio_cols_in_map$ColumnName)

  cond_groups <- list()
  for (cond in unique(bio_cols_in_map$Condition)) {
    samps <- bio_cols_in_map$ColumnName[bio_cols_in_map$Condition == cond]
    if (length(samps) > 0) cond_groups[[cond]] <- samps
  }

  bio_rep_vec <- ave(seq_along(matched_samples),
                     cond_map[matched_samples],
                     FUN = function(i) paste0(cond_map[matched_samples][i[1]], "_", seq_along(i)))
  if (!is.null(subject_pairing))
    bio_rep_vec <- subject_pairing[matched_samples]

  mixture_col <- bio_cols_in_map$Mixture[match(matched_samples, bio_cols_in_map$ColumnName)]
  design_df <- data.frame(
    sample_name  = matched_samples,
    condition    = cond_map[matched_samples],
    bioreplicate = bio_rep_vec,
    mixture      = mixture_col,
    stringsAsFactors = FALSE
  )

  # ── J. Missingness map ────────────────────────────────────────────────────
  logger("[TMT PTM DEBUG] Phase J: missingness map")
  conds <- unique(design_df$condition)
  miss_rows_list <- lapply(all_sites, function(site) {
    if (!site %in% rownames(ptm_abs_mat)) return(NULL)
    row_mat <- ptm_abs_mat[site, matched_samples, drop = FALSE]
    lapply(conds, function(cond) {
      samps_cond <- design_df$sample_name[design_df$condition == cond]
      vals       <- row_mat[, samps_cond, drop = FALSE]
      n_total    <- ncol(vals)
      n_detected <- sum(!is.na(vals))
      data.frame(Site_ID = site, Condition = cond,
                 Detected = n_detected, Total = n_total,
                 Fraction = if (n_total > 0) n_detected / n_total else 0,
                 stringsAsFactors = FALSE)
    })
  })
  missingness_map <- do.call(rbind, unlist(miss_rows_list, recursive = FALSE))

  logger(paste0("[TMT PTM] Summarization complete: ", n_sites_output, " sites, ",
                length(matched_samples), " biological channels, ",
                length(conds), " conditions."))

  # ── K. Return cache matching dia_ptm_master_cache contract ───────────────
  norm_label <- if (identical(as.integer(design_opt), 2L))
    "log2_IRS_stoichiometric_adjustment"
  else
    "log2_median_centering_stoichiometric_adjustment"

  # Site count resting on more than one peptidoform, restricted to the sites that
  # actually survived every filter. n_sites_multi_peptidoform_all was computed on
  # the pre-filter set; the methods text describes the output matrix, so it is the
  # surviving rows that must be counted.
  n_sites_multi_peptidoform <- if (length(all_sites) > 0L)
    sum(n_peptidoforms_by_site[intersect(names(n_peptidoforms_by_site),
                                         all_sites)] > 1L, na.rm = TRUE)
  else 0L
  # Which abundance column family the phospho exports actually carried. "raw" is
  # the required case; "scaled_only" means PD had already normalised the values,
  # which the methods text must disclose.
  .abund_src <- unique(vapply(ptm_parsed_files, function(.f) {
    .v <- .f$abundance_source
    if (is.null(.v) || length(.v) == 0L) "unknown" else as.character(.v)[1]
  }, character(1)))

  list(
    ptm_abs_matrix    = ptm_abs_mat,
    ptm_adj_matrix    = ptm_adj_mat,
    global_matrix     = global_mat,
    design_df         = design_df,
    cond_groups       = cond_groups,
    acc_to_gene       = acc_to_gene,
    site_to_gene      = site_to_gene,
    site_to_prob      = site_to_prob,
    n_peptidoforms_by_site = n_peptidoforms_by_site,
    adj_to_abs_map    = adj_to_abs_map,
    missingness_map   = missingness_map,
    n_sites_input     = n_sites_input,
    n_obs_collected   = n_obs_before_loc,
    n_obs_below_own_threshold = n_obs_below_own_threshold,
    n_obs_no_own_prob = n_obs_no_own_prob,
    # --- WP2: recorded for the methods text; none of these affect the analysis ---
    n_pep_groups_input     = n_pep_groups_input,
    n_pep_groups_conf_kept = n_pep_groups_conf_kept,
    pep_confidence_basis   = pep_confidence_basis,
    pep_qvalue_max         = pep_qvalue_max,
    prot_qvalue_max        = prot_qvalue_max,
    abundance_source       = paste(.abund_src, collapse = ";"),
    n_sites_passing_prob   = n_sites_pass,
    n_sites_after_valid_filter = n_sites_after_valid_filter,
    n_sites_multi_peptidoform  = n_sites_multi_peptidoform,
    n_shared_peptide_sites     = length(shared_pep_sites),
    n_plexes_detected      = n_plexes_detected,
    site_min_plex_n        = site_min_plex_n,
    # --- WP3: plex-coverage policy and IRS partial-correction accounting ---
    n_sites_dropped_post_irs = n_sites_dropped_post_irs,
    irs_policy             = if (!is.null(site_irs_stats)) site_irs_stats$policy else NA_character_,
    irs_n_blanked          = if (!is.null(site_irs_stats)) site_irs_stats$n_blanked else NA_integer_,
    irs_n_plex_no_bridge   = if (!is.null(site_irs_stats)) site_irs_stats$n_plex_no_bridge else NA_integer_,
    prot_irs_n_blanked     = if (!is.null(glob_result$irs_stats)) glob_result$irs_stats$n_blanked else NA_integer_,
    # --- end WP2 additions ---
    n_sites_output    = n_sites_output,
    n_sites_adjusted  = nrow(ptm_adj_mat),
    n_shared_expanded = n_shared_expanded,
    n_razor_reassigned = n_razor_reassigned,
    shared_peptide_sites = shared_pep_sites,
    site_all_accessions  = site_all_accessions,
    site_alt_accessions  = site_alt_accessions,
    n_no_parent_protein = n_no_parent_protein,
    matched_samples   = matched_samples,
    prob_threshold    = prob_threshold,
    min_valid_per_cond = as.integer(min_valid_per_cond),
    min_plex_frac        = min_plex_frac,
    prot_min_plex_frac = prot_min_plex_frac,
    site_min_plex_frac = site_min_plex_frac,
    n_sites_dropped_plex = n_sites_dropped_plex,
    interference_filtered   = interference_col_found,
    interference_threshold  = interference_threshold,
    n_interf_kept_na        = n_interf_kept_na,
    n_interf_kept_pass      = n_interf_kept_pass,
    n_interf_removed        = n_interf_removed,
    design_opt        = as.integer(design_opt),
    normalization_basis = normalization_basis,
    adj_centering_delta = adj_centering_delta,
    pipeline_label    = "TMT PTM (limma, stoichiometric adjustment)",
    software          = "pd_tmt",
    normalization     = norm_label
  )
}

# -----------------------------------------------------------------------------
# dia_ptm_impute_cache()
# k-NN imputation of the DIA PTM cache matrices, computed ONCE over every sample
# in the cache instead of once per contrast.  Imputing inside the DE step makes
# the donor pool "whichever two conditions are being compared", so the same site
# receives different filled values - and a different p-value - in A-vs-B than in
# A-vs-C.  Computing it here removes that dependence.
#
# Eligibility is the same rule as dia_ptm_limma_de()'s knn_max_missing_per_cond -
# a COUNT of missing values allowed per condition, matching dia_limma_de()'s
# max_missing_per_cond - but evaluated in EVERY condition of the cache so that
# it, too, is independent of the contrast.  Ineligible rows keep their NAs and
# are fitted by lmFit on the observations that exist.
#
# The DE-stage missingness filter is deliberately NOT applied here: it stays
# contrast-scoped, since a site observed only in a third condition should not be
# tested in A-vs-B.  Only the donor pool moves.
#
# Returns list(ptm_abs_matrix_imputed, global_matrix_imputed, impute_stats), or
# NULL for a TMT PTM cache or when the impute package is unavailable.
# -----------------------------------------------------------------------------
dia_ptm_impute_cache <- function(cache, k = 10L,
                                 knn_max_missing_per_cond = 1L,
                                 logger = message) {
  if (is.null(cache) || isTRUE(cache$software == "pd_tmt")) return(NULL)
  if (!requireNamespace("impute", quietly = TRUE)) {
    logger("[DIA PTM] 'impute' package unavailable - cache-level k-NN skipped.")
    return(NULL)
  }
  design_df <- cache$design_df
  if (is.null(design_df) || is.null(design_df$sample_name)) return(NULL)

  .impute_one <- function(mat, mat_label) {
    if (is.null(mat) || nrow(mat) < 2L)
      return(list(mat = mat, n_skipped = 0L, n_filled = 0L,
                  n_rows = NA_integer_, n_eligible = NA_integer_))
    samps <- intersect(design_df$sample_name, colnames(mat))
    if (length(samps) < 2L)
      return(list(mat = mat, n_skipped = 0L, n_filled = 0L,
                  n_rows = NA_integer_, n_eligible = NA_integer_))
    sub      <- mat[, samps, drop = FALSE]
    n_before <- sum(is.na(sub))
    if (n_before == 0L)
      return(list(mat = mat, n_skipped = 0L, n_filled = 0L,
                  n_rows = NA_integer_, n_eligible = NA_integer_))

    conds <- design_df$condition[match(samps, design_df$sample_name)]
    elig  <- rep(TRUE, nrow(sub))
    for (.g in unique(conds[!is.na(conds)])) {
      .idx <- which(conds == .g)
      if (length(.idx) == 0L) next
      elig <- elig &
        (rowSums(is.na(sub[, .idx, drop = FALSE])) <= knn_max_missing_per_cond)
    }
    elig[is.na(elig)] <- FALSE
    n_skipped <- sum(!elig)

    if (sum(elig) >= 2L && any(is.na(sub[elig, , drop = FALSE]))) {
      .imp <- tryCatch(
        impute::impute.knn(sub[elig, , drop = FALSE], k = k,
                           rowmax = 1, colmax = 0.8,
                           maxp = Inf, rng.seed = 42L)$data,
        error = function(e) {
          logger(paste0("[DIA PTM] Cache k-NN failed for ", mat_label,
                        " - NAs retained."))
          NULL
        })
      if (!is.null(.imp)) sub[elig, ] <- .imp
    }
    n_filled <- n_before - sum(is.na(sub))
    out <- mat
    out[, samps] <- sub
    logger(paste0("[DIA PTM] Cache k-NN (", mat_label, "): filled ", n_filled,
                  " of ", n_before, " missing values across ", length(samps),
                  " samples; ", n_skipped,
                  " row(s) above the per-condition cap kept their NAs."))
    list(mat = out, n_skipped = n_skipped, n_filled = n_filled,
         n_rows = nrow(sub), n_eligible = sum(elig))
  }

  .a <- .impute_one(cache$ptm_abs_matrix, "absolute PTM")
  .g <- .impute_one(cache$global_matrix,  "protein")

  list(
    ptm_abs_matrix_imputed = .a$mat,
    global_matrix_imputed  = .g$mat,
    ptm_impute_cap         = knn_max_missing_per_cond,
    impute_stats = list(
      k                        = k,
      knn_max_missing_per_cond = knn_max_missing_per_cond,
      scope                    = "cache (all samples in design_df)",
      n_rows_skipped_ptm       = .a$n_skipped,
      n_values_filled_ptm      = .a$n_filled,
      n_rows_ptm               = .a$n_rows,
      n_eligible_ptm           = .a$n_eligible,
      n_rows_skipped_protein   = .g$n_skipped,
      n_values_filled_protein  = .g$n_filled,
      n_rows_protein           = .g$n_rows,
      n_eligible_protein       = .g$n_eligible
    )
  )
}

# -----------------------------------------------------------------------------
# dia_ptm_limma_de()
# Run limma-trend DE on DIA PTM data (adjusted, absolute, and protein matrices).
# Returns a list:
#   list(adjusted, unadjusted, protein, comparison,
#        n_tested_adjusted, n_tested_unadjusted, n_tested_protein)
# Each sub-element is a data.frame with columns:
#   Protein, log2FC, SE, DF, pvalue, adj.pvalue, Label
# -----------------------------------------------------------------------------
dia_ptm_limma_de <- function(cache, cond_test, cond_ref,
                              knn_max_missing_per_cond = 1L,
                              # FLOOR, not the threshold: on the spectronaut and
                              # pivot routes the effective per-group requirement is
                              # max(this, ceiling(0.5 * n_samples_in_group)).
                              de_min_valid_per_cond    = 2L,
                              batch_groups     = NULL,
                              batch_for_model  = FALSE,
                              subject_pairing  = NULL,
                              logger           = message) {

  if (is.null(cache) || is.null(cache$ptm_abs_matrix))
    return(list(error = "No DIA PTM cache found. Please run the upload wizard first."))

  design_df   <- cache$design_df
  comp_label  <- paste0(cond_test, " vs ", cond_ref)

  # TMT data has MNAR missingness (below-detection dropout); lmFit handles NAs
  # natively on available data, which is statistically more honest than filling
  # in with k-NN neighbours.  DIA data tends toward MCAR/MAR (stochastic
  # peptide sampling), so k-NN imputation is defensible there.
  use_knn_imputation <- !isTRUE(cache$software == "pd_tmt")

  # The Adjusted table is produced by fitting limma on cache$ptm_adj_matrix, the
  # per-sample difference (PTM log2 minus protein log2).
  #
  # FINDING, recorded so it is not re-derived: contrast-level variance
  # propagation in the manner of MSstatsPTM - deriving the adjusted result from
  # the Absolute PTM and Protein fits with SE = sqrt(SE_ptm^2 + SE_prot^2) and a
  # Welch-Satterthwaite df - was implemented and tested here, and is NOT
  # applicable to this pipeline.  MSstatsPTM propagates two variances because its
  # protein term comes from a separate group-level model.  Here the adjustment is
  # a PER-SAMPLE subtraction, so the adjusted matrix's residual variance already
  # contains the protein-level variance and limma estimates it from the data;
  # there is nothing to add back.  Measured on complete mock data: log2FC
  # identical to the matrix fit, median SE ratio 0.997 scattered either side of
  # one.  Under missingness the propagated form is actively anti-conservative
  # (median SE ratio 0.85, hit count 10 -> 17) because the absolute and protein
  # fits are k-NN imputed while the adjusted fit is not.  Do not reintroduce it.

  # Collapse guard for the k-NN eligibility cap.  When a condition has few enough
  # samples that n - cap is at or below the testing filter's own requirement, the
  # cap admits everything that passed that filter and does nothing.  Warn only -
  # behaviour is deliberately unchanged.
  knn_cap_degenerate <- FALSE
  if (use_knn_imputation && !is.null(design_df$condition)) {
    .n_per_cond <- table(design_df$condition[
      design_df$sample_name %in% colnames(cache$ptm_abs_matrix)])
    for (.g in names(.n_per_cond)) {
      .n_g <- as.integer(.n_per_cond[[.g]])
      if (.n_g - knn_max_missing_per_cond <= de_min_valid_per_cond) {
        knn_cap_degenerate <- TRUE
        logger(paste0("[DIA PTM DE] WARNING: k-NN eligibility cap is degenerate for ",
                      "condition '", .g, "': n = ", .n_g, ", cap = ",
                      knn_max_missing_per_cond, ", so the cap requires >= ",
                      .n_g - knn_max_missing_per_cond,
                      " valid values while de_min_valid_per_cond is ",
                      de_min_valid_per_cond,
                      ". The cap admits every row that passed the testing filter."))
      }
    }
  }

  # Cache-scoped imputation.  Reuse the memoised matrices when the caller has
  # already stored them and they were built with the same eligibility cap;
  # otherwise compute them now (lazily, on the first DE run - not at
  # summarization time, so existing prerun bundles keep loading).
  # Never runs for pd_tmt: use_knn_imputation is FALSE there and
  # dia_ptm_impute_cache() returns NULL for that cache anyway.
  # The memo is invalidated whenever the STORED cap (cache$ptm_impute_cap, written
  # alongside the matrices) differs from the cap in force now.  A bundle restored
  # from disk can carry matrices imputed under a different cap, so this must read
  # a stored value and never assume one.
  imputed_cache <- NULL
  if (use_knn_imputation) {
    .stored_cap <- cache$ptm_impute_cap %||%
                   cache$impute_stats$knn_max_missing_per_cond
    if (!is.null(cache$ptm_abs_matrix_imputed) &&
        (!ptm_has_global(cache) || !is.null(cache$global_matrix_imputed)) &&
        !is.null(.stored_cap) &&
        isTRUE(all.equal(.stored_cap, knn_max_missing_per_cond))) {
      imputed_cache <- list(
        ptm_abs_matrix_imputed = cache$ptm_abs_matrix_imputed,
        global_matrix_imputed  = cache$global_matrix_imputed,
        impute_stats           = cache$impute_stats,
        ptm_impute_cap         = .stored_cap
      )
    } else {
      imputed_cache <- dia_ptm_impute_cache(
        cache, k = 10L,
        knn_max_missing_per_cond = knn_max_missing_per_cond,
        logger = logger)
    }
  }

  # Mirror of the knn_cap_degenerate guard above: that one catches a cap too LOOSE
  # to do any work, this one catches a cap so TIGHT that essentially nothing is
  # imputed.  n_rows / n_eligible are NA when a matrix had no missing values to
  # begin with (the DDA MaxQuant path, already MinProb-imputed at summarization),
  # and the guard is skipped in that case - there was nothing to exclude.
  # Warning only: the cap, its default and every eligibility decision are unchanged.
  knn_cap_exhaustive <- FALSE
  if (use_knn_imputation && !is.null(imputed_cache$impute_stats)) {
    .st <- imputed_cache$impute_stats
    for (.m in c("ptm", "protein")) {
      .n_rows <- .st[[paste0("n_rows_", .m)]]
      .n_elig <- .st[[paste0("n_eligible_", .m)]]
      if (is.null(.n_rows) || is.null(.n_elig) ||
          is.na(.n_rows) || is.na(.n_elig) || .n_rows == 0L) next
      if (.n_elig == 0L || .n_elig < 0.10 * .n_rows) {
        knn_cap_exhaustive <- TRUE
        .n_per_cond2 <- table(design_df$condition[
          design_df$sample_name %in% colnames(cache$ptm_abs_matrix)])
        logger(paste0("[DIA PTM DE] WARNING: k-NN eligibility cap excluded ",
                      .n_rows - .n_elig, " of ", .n_rows, " row(s) (",
                      round(100 * (.n_rows - .n_elig) / .n_rows, 1),
                      "%) from the ", .m, " matrix. Cap = ",
                      knn_max_missing_per_cond, " missing value(s) per condition; ",
                      "replicates per condition: ",
                      paste0(names(.n_per_cond2), "=", as.integer(.n_per_cond2),
                             collapse = ", "),
                      ". Imputation filled little or nothing."))
      }
    }
  }

  # Parent accession from an adjusted-matrix row key.  Delegates to the shared
  # ptm_parent_acc_of() so the DE route, PTM Single and any future consumer parse
  # site keys with one implementation.  Signature and behaviour are unchanged.
  .parent_acc_of <- function(keys) ptm_parent_acc_of(keys, cache$software)

  .run_limma_on_mat <- function(mat, mat_label, use_knn = use_knn_imputation,
                                detect_mat = NULL, use_trend = TRUE,
                                imputed_mat = NULL) {
    groups_needed <- c(cond_test, cond_ref)
    use_samples   <- design_df$sample_name[design_df$condition %in% groups_needed]
    use_samples   <- use_samples[use_samples %in% colnames(mat)]
    if (length(use_samples) < 2)
      return(list(error = paste0("Not enough samples for ", mat_label)))

    sub_mat    <- mat[, use_samples, drop = FALSE]
    sub_design <- design_df[design_df$sample_name %in% use_samples, , drop = FALSE]
    group_vec  <- sub_design$condition[match(colnames(sub_mat), sub_design$sample_name)]

    # Missingness filter: keep a row if AT LEAST ONE of the two contrasted groups
    # carries valid values in at least half that group's samples, subject to a
    # floor of de_min_valid_per_cond.  Each group therefore gets its own
    # threshold from its OWN size, which matters when the groups are unbalanced.
    #
    # The thresholds are precomputed ONCE here, outside the row loop, and indexed
    # BY NAME against tapply()'s output inside the loop.  tapply() returns values
    # in factor-level order; relying on that order matching a separately computed
    # vector is exactly how this goes silently wrong, so it is never done.
    #
    # ROUTE GATE: scaling applies to the spectronaut and Spectronaut-pivot routes
    # only.  pd_tmt (TMT PTM) and the DDA MaxQuant path keep the plain floor and
    # are bit-identical to before.  cache$acquisition_mode exists only on the
    # maxquant caches; it is NULL elsewhere, which identical() handles correctly.
    .de_min_valid <- as.integer(de_min_valid_per_cond)
    .scale_thr    <- !isTRUE(cache$software == "pd_tmt") &&
                     !identical(cache$acquisition_mode, "dda")
    .grp_sizes    <- table(group_vec)
    de_min_valid_effective <- vapply(names(.grp_sizes), function(.g)
      if (.scale_thr)
        max(.de_min_valid, as.integer(ceiling(0.5 * .grp_sizes[[.g]])))
      else
        .de_min_valid,
      integer(1))
    names(de_min_valid_effective) <- names(.grp_sizes)

    # At n = 2 per condition a both-arms requirement of 2 is a complete-case
    # analysis, and the presence/absence branch has no attainable significance
    # (smallest two-sided Fisher p at 2 vs 2 is 0.333). Warn; do not alter it.
    if (isTRUE(cache$software == "pd_tmt") && any(.grp_sizes <= 2L))
      logger(paste0("[DIA PTM DE] WARNING: ", mat_label, " - condition(s) ",
                    paste(names(.grp_sizes)[.grp_sizes <= 2L], collapse = ", "),
                    " have <= 2 replicates. With a both-arms requirement of ",
                    .de_min_valid, " this is a complete-case analysis (no missing ",
                    "values tolerated), and the presence/absence set cannot reach ",
                    "significance at this replicate count."))

    # TMT PTM requires the minimum in BOTH arms: a fold change needs both arms,
    # and single-arm sites belong in the presence/absence analysis, not here.
    # Other routes keep the any-arm rule and are unchanged.
    .require_both <- isTRUE(cache$software == "pd_tmt")
    keep_rows  <- apply(sub_mat, 1, function(row) {
      .n_valid <- tapply(!is.na(row), group_vec, sum)
      if (.require_both)
        all(.n_valid >= de_min_valid_effective[names(.n_valid)])
      else
        any(.n_valid >= de_min_valid_effective[names(.n_valid)])
    })
    if (!any(keep_rows)) return(list(error = paste0("No rows pass missingness filter for ", mat_label)))
    orig_sub_mat <- sub_mat[keep_rows, , drop = FALSE]
    sub_mat      <- orig_sub_mat
    n_tested     <- nrow(sub_mat)

    # Detection fractions per group  -" computed BEFORE k-NN imputation.
    # Produces Found_<cond_test> and Found_<cond_ref> (e.g. "7/10") plus
    # high_missingness = TRUE when >75% of samples in any group are missing.
    # These rows sit above the k-NN eligibility cap and are never imputed, so
    # the flag records missingness, not imputation.
    # When detect_mat is supplied, detection is read from that matrix instead of
    # the matrix being fitted, so the adjusted table reports site missingness and
    # not "site OR parent protein missing".  Row keys are resolved through
    # cache$adj_to_abs_map; rows that do not resolve fall back to the fitted
    # matrix and are counted.  detect_mat = NULL preserves the original behaviour.
    detect_src        <- orig_sub_mat
    n_detect_fallback <- 0L
    if (!is.null(detect_mat) &&
        all(colnames(orig_sub_mat) %in% colnames(detect_mat))) {
      .src_keys <- rownames(orig_sub_mat)
      .abs_map  <- cache$adj_to_abs_map
      .res_keys <- if (!is.null(.abs_map)) unname(.abs_map[.src_keys])
                   else rep(NA_character_, length(.src_keys))
      .res_keys[is.na(.res_keys)] <- .src_keys[is.na(.res_keys)]
      .hit              <- .res_keys %in% rownames(detect_mat)
      n_detect_fallback <- sum(!.hit)
      if (any(.hit))
        detect_src[.hit, ] <- detect_mat[.res_keys[.hit],
                                         colnames(orig_sub_mat), drop = FALSE]
      if (n_detect_fallback > 0L)
        logger(paste0("[DIA PTM DE] ", mat_label, ": ", n_detect_fallback,
                      " row(s) unresolved against the detection matrix - detection taken from the fitted matrix."))
    }

    detect_df <- {
      d_res <- data.frame(row.names = rownames(orig_sub_mat))
      for (.g in c(cond_test, cond_ref)) {
        .g_cols <- which(group_vec == .g)
        if (length(.g_cols) == 0) next
        .g_mat  <- detect_src[, .g_cols, drop = FALSE]
        .n_obs  <- rowSums(!is.na(.g_mat))
        d_res[[paste0("Found_", .g)]] <- paste0(.n_obs, "/", ncol(.g_mat))
      }
      d_res$high_missingness <- apply(detect_src, 1, function(.row)
        any(tapply(is.na(.row), group_vec, mean) > 0.75))
      d_res
    }

    # k-NN imputation  -" DIA only (MCAR/MAR missingness assumption).
    # For TMT data (MNAR), NAs are left as-is and lmFit fits each row on
    # its available observations, preserving the true missingness structure.
    n_before_imp       <- sum(is.na(sub_mat))
    n_knn_rows_skipped <- 0L
    n_knn_values       <- 0L
    knn_scope          <- "none"
    if (use_knn && !is.null(imputed_mat) &&
        all(rownames(orig_sub_mat) %in% rownames(imputed_mat)) &&
        all(colnames(orig_sub_mat) %in% colnames(imputed_mat))) {
      # Cache-scoped imputation already done over every sample: just subset it.
      sub_mat      <- imputed_mat[rownames(orig_sub_mat),
                                  colnames(orig_sub_mat), drop = FALSE]
      n_knn_values <- n_before_imp - sum(is.na(sub_mat))
      knn_scope    <- "cache"
      logger(paste0("[DIA PTM DE] ", mat_label,
                    ": using cache-scoped k-NN matrix (", n_knn_values, " of ",
                    n_before_imp, " missing values already filled)."))
    } else if (use_knn && n_before_imp > 0 &&
        requireNamespace("impute", quietly = TRUE)) {
      knn_scope <- "contrast"
      k_val   <- 10L
      # Per-condition eligibility.  rowmax stays at 1 deliberately: a lower rowmax
      # makes impute.knn swap k-NN for its column-mean fallback, which deflates
      # variance.  But at rowmax = 1 no row is ever excluded, so a site with 2
      # observed values out of 10 - the minimum that passes the missingness filter
      # above - would be filled in the other 8 cells and still receive a p-value.
      # Rows with more than knn_max_missing_per_cond missing values in either
      # contrasted condition are therefore routed out of imputation and keep their
      # NAs for lmFit to handle on the observations that exist.  The cap is a
      # COUNT, not a fraction, matching dia_limma_de()'s max_missing_per_cond: a
      # fraction resolves to a different integer at every replicate count.
      .i_test <- which(group_vec == cond_test)
      .i_ref  <- which(group_vec == cond_ref)
      .f_test <- if (length(.i_test) > 0L)
                   rowSums(is.na(sub_mat[, .i_test, drop = FALSE]))
                 else rep(0L, nrow(sub_mat))
      .f_ref  <- if (length(.i_ref) > 0L)
                   rowSums(is.na(sub_mat[, .i_ref, drop = FALSE]))
                 else rep(0L, nrow(sub_mat))
      elig <- (.f_test <= knn_max_missing_per_cond) &
              (.f_ref  <= knn_max_missing_per_cond)
      elig[is.na(elig)] <- FALSE
      n_knn_rows_skipped <- sum(!elig)

      if (sum(elig) >= 2L && any(is.na(sub_mat[elig, , drop = FALSE]))) {
        imp_res <- tryCatch(
          impute::impute.knn(sub_mat[elig, , drop = FALSE], k = k_val,
                             rowmax = 1, colmax = 0.8,
                             maxp = Inf, rng.seed = 42L)$data,
          error = function(e) { logger("[DIA PTM] k-NN imputation failed - NAs retained."); NULL }
        )
        if (!is.null(imp_res)) sub_mat[elig, ] <- imp_res
      }
      n_knn_values <- n_before_imp - sum(is.na(sub_mat))
      logger(paste0("[DIA PTM DE] ", mat_label, ": k-NN filled ", n_knn_values,
                    " of ", n_before_imp, " missing values; ", n_knn_rows_skipped,
                    " row(s) exceeded the per-condition missingness cap (",
                    round(knn_max_missing_per_cond, 3),
                    ") and kept their NAs."))
    }

    # Build design matrix
    group_f    <- factor(group_vec)

    if (!is.null(batch_groups) && isTRUE(batch_for_model)) {
      batch_map  <- setNames(
        rep(names(batch_groups), times = vapply(batch_groups, length, integer(1))),
        unlist(batch_groups)
      )
      batch_vec  <- factor(batch_map[use_samples])
      design_mat <- model.matrix(~ 0 + group_f + batch_vec)
      colnames(design_mat) <- make.names(sub("^group_f", "", colnames(design_mat)))
    } else {
      design_mat <- model.matrix(~ 0 + group_f)
      colnames(design_mat) <- make.names(levels(group_f))
    }

    # Fit model (with optional duplicateCorrelation for pairing)
    if (!is.null(subject_pairing)) {
      block_vec <- subject_pairing[use_samples]
      corfit    <- tryCatch(
        limma::duplicateCorrelation(sub_mat, design_mat, block = block_vec),
        error = function(e) NULL
      )
      cc <- if (!is.null(corfit)) corfit$consensus.correlation else NA_real_
      if (!is.null(corfit) && is.finite(cc) && abs(cc) < 0.99) {
        fit <- limma::lmFit(sub_mat, design_mat,
                            block = block_vec,
                            correlation = cc)
      } else {
        if (!is.null(corfit))
          logger(paste0("[DIA PTM] Unusable duplicateCorrelation consensus (", cc,
                        ") - falling back to unpaired fit."))
        fit <- limma::lmFit(sub_mat, design_mat)
      }
    } else {
      fit <- limma::lmFit(sub_mat, design_mat)
    }

    # Contrasts
    test_name <- make.names(cond_test)
    ref_name  <- make.names(cond_ref)
    contrast_string <- paste0(test_name, " - ", ref_name)
    contrast_mat    <- tryCatch(
      limma::makeContrasts(contrasts = contrast_string, levels = design_mat),
      error = function(e) NULL
    )
    if (is.null(contrast_mat))
      return(list(error = paste0("Could not build contrast for ", mat_label, ": ", contrast_string)))

    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    # trend = TRUE regresses the prior variance on Amean.  That is correct for the
    # absolute and protein matrices, where Amean is a mean log-abundance.  On the
    # adjusted matrix Amean is a mean log-RATIO centred near zero and carries no
    # abundance information, so the caller passes use_trend = FALSE there.
    fit3 <- limma::eBayes(fit2, trend = isTRUE(use_trend), robust = TRUE)

    top  <- limma::topTable(fit3, coef = 1, number = Inf,
                            sort.by = "none", adjust.method = "BH")

    prots <- rownames(top)
    df_vec <- fit3$df.residual + fit3$df.prior
    se_vec <- sqrt(fit3$s2.post) * fit3$stdev.unscaled[prots, 1]

    result_df <- data.frame(
      Protein    = prots,
      log2FC     = top$logFC,
      SE         = se_vec,
      DF         = df_vec,
      pvalue     = top$P.Value,
      adj.pvalue = top$adj.P.Val,
      Label      = comp_label,
      stringsAsFactors = FALSE
    )
    result_df <- cbind(result_df, detect_df[result_df$Protein, , drop = FALSE])

    logger(paste0("[DIA PTM DE] ", mat_label, ": ", n_tested, " sites tested."))
    list(result_df = result_df, n_tested = n_tested,
         n_estimable = sum(!is.na(result_df$log2FC)),
         n_detect_fallback = n_detect_fallback,
         n_knn_rows_skipped = n_knn_rows_skipped,
         n_knn_values = n_knn_values,
         knn_scope = knn_scope,
         de_min_valid_effective = de_min_valid_effective)
  }

  # ── Tier C: sites that cannot be tested by limma ──────────────────────────
  # A site failing the valid-value filter is dropped from the fit below and would
  # otherwise vanish without a result.  Test instead whether its DETECTION pattern
  # differs between the two conditions (Fisher's exact on detected / not-detected),
  # BH-corrected within this excluded set only.  This is the same treatment the
  # Global DIA route gives proteins failing the same filter.
  #
  # The correction is independent of the one applied to the limma-tested set, so no
  # joint false discovery rate is claimed across the two.  Run on the ABSOLUTE
  # matrix only: NAs in the adjusted matrix conflate "site not detected" with
  # "parent protein not quantified".
  #
  # The threshold derivation below duplicates .run_limma_on_mat()'s deliberately,
  # so that this function cannot perturb the fitted path.  Do not factor it out.
  .presence_absence <- function(mat) {
    if (is.null(mat) || nrow(mat) == 0L) return(NULL)
    use_samples <- design_df$sample_name[design_df$condition %in% c(cond_test, cond_ref)]
    use_samples <- use_samples[use_samples %in% colnames(mat)]
    if (length(use_samples) < 2) return(NULL)
    sub_mat   <- mat[, use_samples, drop = FALSE]
    sub_des   <- design_df[design_df$sample_name %in% use_samples, , drop = FALSE]
    group_vec <- sub_des$condition[match(colnames(sub_mat), sub_des$sample_name)]

    .de_min_valid <- as.integer(de_min_valid_per_cond)
    .scale_thr    <- !isTRUE(cache$software == "pd_tmt") &&
                     !identical(cache$acquisition_mode, "dda")
    .grp_sizes    <- table(group_vec)
    .thr <- vapply(names(.grp_sizes), function(.g)
      if (.scale_thr) max(.de_min_valid, as.integer(ceiling(0.5 * .grp_sizes[[.g]])))
      else .de_min_valid, integer(1))
    names(.thr) <- names(.grp_sizes)

    .require_both <- isTRUE(cache$software == "pd_tmt")
    keep <- apply(sub_mat, 1, function(row) {
      .n_valid <- tapply(!is.na(row), group_vec, sum)
      if (.require_both) all(.n_valid >= .thr[names(.n_valid)])
      else               any(.n_valid >= .thr[names(.n_valid)])
    })
    excl <- sub_mat[!keep, , drop = FALSE]
    if (nrow(excl) == 0L)
      return(list(table = NULL, n = 0L, has_power = NA, min_achievable_p = NA_real_))

    i_t <- which(group_vec == cond_test)
    i_r <- which(group_vec == cond_ref)
    n_t <- length(i_t); n_r <- length(i_r)
    if (n_t == 0L || n_r == 0L) return(NULL)

    d_t <- rowSums(!is.na(excl[, i_t, drop = FALSE]))
    d_r <- rowSums(!is.na(excl[, i_r, drop = FALSE]))

    # Smallest two-sided p attainable at this replicate count, taken from the most
    # extreme table.  When it exceeds 0.05 no row in this set can reach
    # significance, so p-values are reported as NA and the set is a ranked
    # descriptive list rather than a test result.
    min_p <- tryCatch(
      stats::fisher.test(matrix(c(n_t, 0L, 0L, n_r), nrow = 2))$p.value,
      error = function(e) NA_real_)
    has_power <- isTRUE(is.finite(min_p) && min_p <= 0.05)

    pv <- if (!has_power) rep(NA_real_, nrow(excl)) else
      vapply(seq_len(nrow(excl)), function(i) {
        tryCatch(stats::fisher.test(
          matrix(c(d_t[i], n_t - d_t[i], d_r[i], n_r - d_r[i]), nrow = 2))$p.value,
          error = function(e) NA_real_)
      }, numeric(1))

    out <- data.frame(
      Protein        = rownames(excl),
      Detection_diff = (d_t / n_t) - (d_r / n_r),
      pvalue         = pv,
      adj.pvalue     = if (has_power) stats::p.adjust(pv, method = "BH")
                       else rep(NA_real_, nrow(excl)),
      Label          = comp_label,
      stringsAsFactors = FALSE
    )
    out[[paste0("Found_", cond_test)]] <- paste0(d_t, "/", n_t)
    out[[paste0("Found_", cond_ref)]]  <- paste0(d_r, "/", n_r)
    out <- out[order(-abs(out$Detection_diff)), , drop = FALSE]
    list(table = out, n = nrow(out), has_power = has_power, min_achievable_p = min_p)
  }
  nmar_res <- tryCatch(.presence_absence(cache$ptm_abs_matrix),
                       error = function(e) { logger(paste0(
                         "[DIA PTM DE] presence/absence analysis failed - ",
                         conditionMessage(e))); NULL })
  if (!is.null(nmar_res) && !is.null(nmar_res$n))
    logger(paste0("[DIA PTM DE] Presence/absence: ", nmar_res$n,
                  " site(s) failed the valid-value filter and were analysed by ",
                  "Fisher's exact test (",
                  if (isTRUE(nmar_res$has_power)) "BH-corrected within this set"
                  else "no attainable significance at this replicate count; p reported as NA",
                  ")."))

  # Stoichiometric-adjusted matrix: leave NAs for lmFit to fit on available
  # observations (k-NN on a derived difference matrix is less well-motivated
  # than on raw abundance matrices). Raw and protein matrices still use k-NN.
  # A PTM-only cache (ptm_layers() = "unadjusted") has no adjusted or protein
  # matrix: only the Absolute PTM fit runs, and the other two results are NULL.
  .layers   <- ptm_layers(cache)
  adj_res   <- if ("adjusted" %in% .layers)
                 .run_limma_on_mat(cache$ptm_adj_matrix,  "Adjusted",     use_knn = FALSE,
                                   detect_mat = cache$ptm_abs_matrix,
                                   use_trend  = FALSE)
  unadj_res <- .run_limma_on_mat(cache$ptm_abs_matrix,  "Absolute PTM",
                                 imputed_mat = imputed_cache$ptm_abs_matrix_imputed)
  prot_res  <- if ("protein" %in% .layers)
                 .run_limma_on_mat(cache$global_matrix,   "Protein",
                                   imputed_mat = imputed_cache$global_matrix_imputed)

  if (!is.null(adj_res$error))   return(list(error = adj_res$error))
  if (!is.null(unadj_res$error)) return(list(error = unadj_res$error))
  if (!is.null(prot_res$error))  return(list(error = prot_res$error))

  # Was the parent protein quantified at all?  Distinguishes "site missing" from
  # "protein never measured" on the adjusted table.
  if (!is.null(adj_res))
    adj_res$result_df$Protein_Available <-
      .parent_acc_of(adj_res$result_df$Protein) %in% rownames(cache$global_matrix)

  list(
    adjusted            = adj_res$result_df,
    unadjusted          = unadj_res$result_df,
    protein             = prot_res$result_df,
    comparison          = comp_label,
    n_tested_adjusted   = adj_res$n_tested,
    n_tested_unadjusted = unadj_res$n_tested,
    n_tested_protein    = prot_res$n_tested,
    n_estimable_adjusted   = adj_res$n_estimable,
    n_estimable_unadjusted = unadj_res$n_estimable,
    n_estimable_protein    = prot_res$n_estimable,
    n_detect_fallback_adjusted = adj_res$n_detect_fallback,
    knn_max_missing_per_cond   = knn_max_missing_per_cond,
    n_knn_rows_skipped_unadjusted = unadj_res$n_knn_rows_skipped,
    n_knn_rows_skipped_protein    = prot_res$n_knn_rows_skipped,
    n_knn_values_unadjusted       = unadj_res$n_knn_values,
    n_knn_values_protein          = prot_res$n_knn_values,
    knn_scope                     = unadj_res$knn_scope,
    imputed_cache                 = imputed_cache,
    de_min_valid_per_cond         = as.integer(de_min_valid_per_cond),
    de_require_both_arms          = isTRUE(cache$software == "pd_tmt"),
    # Same per-group thresholds for every fit (they depend on group sizes only);
    # read from the primary fit
    de_min_valid_per_cond_effective = (adj_res %||% unadj_res)$de_min_valid_effective,
    knn_cap_degenerate            = knn_cap_degenerate,
    knn_cap_exhaustive            = knn_cap_exhaustive,
    nmar_table                    = nmar_res$table,
    n_nmar                        = nmar_res$n %||% 0L,
    nmar_has_power                = nmar_res$has_power,
    nmar_min_achievable_p         = nmar_res$min_achievable_p
  )
}

# -----------------------------------------------------------------------------
# ptm_de_batch_groups() / ptm_de_annotate_views() / ptm_de_store_entry() /
# ptm_de_backfill_contrasts()
#
# Shared PTM DE plumbing, lifted out of server_ptm_de.R's run observer.
#
# Why: PTM Single reports the DE tab's moderated statistics and never fits its
# own model, so a contrast the DE tab had not run showed NA fold changes even for
# a site detected in every replicate of every condition. Single could not fit it
# itself because the annotation helpers were observer-locals - nothing outside
# that 400-line observer could call them. They live here now so that a contrast
# fitted from the DE tab and the same contrast fitted on demand from PTM Single
# go through ONE code path and produce identical store entries.
#
# The three annotation closures below are byte-identical to the definitions they
# replace, including their four-space indentation. That is deliberate: it makes
# the move provable by diff. Do not re-indent them.
# -----------------------------------------------------------------------------

# Plex batch covariate: TMT PTM with >=2 plexes and non-IRS normalization
# (design_opt != 2). IRS already equalises plexes, so no batch term there.
ptm_de_batch_groups <- function(cache, design_option) {
  if (isTRUE(design_option == "tmt_ptm") &&
      "mixture" %in% colnames(cache$design_df) &&
      !identical(as.integer(cache$design_opt), 2L)) {
    bg <- split(cache$design_df$sample_name, cache$design_df$mixture)
    if (length(bg) >= 2L) bg else NULL
  } else NULL
}

# Annotate the three dia_ptm_limma_de() result tables with display columns.
ptm_de_annotate_views <- function(res, cache) {
  gene_map <- cache$site_to_gene    # abs_site_id -> gene label(s)
  acc_gene <- cache$acc_to_gene     # accession   -> gene name
  adj_map  <- cache$adj_to_abs_map  # adj_row_name -> abs_site_id
  # A NULL gene map made .to_first_gene() error on `if (!is.na(g) ...)` because
  # NULL[x] is NULL and the condition is length 0. An empty named vector yields
  # NA and falls through to the accession-parsing branch, which is the intended
  # degraded behaviour. No effect when gene_map is present.
  if (is.null(gene_map)) gene_map <- setNames(character(0), character(0))

    .to_first_gene <- function(abs_id) {
      g <- gene_map[abs_id]
      if (!is.na(g) && nzchar(g))
        strsplit(g, ";", fixed = TRUE)[[1]][1L]
      else {
        # Fall back to the accession, using the SAME two patterns the Protein
        # column is parsed with a few lines below, in the same order. The old
        # single-residue-only pattern matched neither ACC_S44+S48_Phospho nor a
        # MaxQuant ACC_S12_M1, so for a site with no gene-map entry Gene became
        # the whole site id and PTM_Label doubled up (ACC_S44+S48_Phospho_S44+S48
        # _Phospho). Only the no-gene-name branch is affected.
        .p <- strsplit(abs_id, ";", fixed = TRUE)[[1]][1L]
        .g <- sub("_(?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+_Phospho.*$", "", .p, perl = TRUE)
        if (identical(.g, .p))
          .g <- sub("_[A-Za-z][0-9]+_M[0-9]+$", "", .p)
        .g
      }
    }

    # For adjusted/unadjusted: Protein is a matrix row name; parse out
    # accession + site from the site_id (resolving through adj_map for adjusted)
    .annotate_ptm_view <- function(df, resolve_map = NULL) {
      if (is.null(df) || nrow(df) == 0L) return(df)
      row_keys <- df$Protein
      abs_ids  <- if (!is.null(resolve_map)) {
        resolved <- resolve_map[row_keys]
        ifelse(is.na(resolved), row_keys, resolved)
      } else {
        row_keys
      }
      primary         <- vapply(abs_ids, function(x)
        strsplit(x, ";", fixed = TRUE)[[1]][1L], character(1))
      df$Gene         <- vapply(abs_ids, .to_first_gene, character(1))
      # ── Spectronaut format: ACC_S32_Phospho[...] ─────────────────────────
      prot_out <- sub("_(?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+_Phospho.*$", "", primary, perl = TRUE)
      site_out <- sub("^.*_((?:[A-Za-z][0-9]+\\+)*[A-Za-z][0-9]+)_Phospho.*$", "\\1", primary, perl = TRUE)
      # ── MaxQuant fallback: ACC_S32_M1 → protein=ACC, site=S32 ───────────
      mq_mask  <- (prot_out == primary)   # regex did not match — still full ID
      if (any(mq_mask)) {
        prot_out[mq_mask] <- sub("_[A-Za-z][0-9]+_M[0-9]+$", "",          primary[mq_mask])
        site_out[mq_mask] <- sub("^.*_([A-Za-z][0-9]+)_M[0-9]+$", "\\1", primary[mq_mask])
      }
      df$Protein      <- prot_out
      df$Site         <- site_out
      df$Modification <- "Phospho"
      df$PTM_Label    <- paste0(df$Gene, "_", df$Site, "_Phospho")
      df$PTM_ID       <- row_keys
      # Every master protein the site's peptide mapped to, razor winner first.
      # One row per site is tested; the alternates are equally valid mappings of
      # the same measurement and are shown so a site can be found under any of
      # its proteins. The column is rendered, so DT's search box matches it.
      .alt_map <- cache$site_alt_accessions
      .all_gene <- gene_map[abs_ids]
      df$Also_Maps_To <- vapply(seq_along(abs_ids), function(i) {
        g <- .all_gene[i]
        if (!is.na(g) && grepl(";", g, fixed = TRUE)) {
          paste(strsplit(g, ";", fixed = TRUE)[[1]][-1L], collapse = ", ")
        } else if (!is.null(.alt_map) && abs_ids[i] %in% names(.alt_map) &&
                   nzchar(.alt_map[[abs_ids[i]]])) {
          gsub(";", ", ", .alt_map[[abs_ids[i]]], fixed = TRUE)
        } else ""
      }, character(1))
      # Shared-peptide sites: one measured intensity vector mapping to more than
      # one master protein. Read from the cache where available, since the site
      # key no longer encodes sharing.
      .shared_set <- cache$shared_peptide_sites
      df$Is_Shared    <- grepl("_SharedInference$", row_keys) |
                         grepl(";", abs_ids, fixed = TRUE) |
                         nzchar(df$Also_Maps_To) |
                         (!is.null(.shared_set) & abs_ids %in% .shared_set)
      df
    }

    # For protein view: Protein is accession; derive Gene from acc_to_gene
    .annotate_prot_view <- function(df) {
      if (is.null(df) || nrow(df) == 0L) return(df)
      df$Gene         <- vapply(df$Protein, function(p)
        if (!is.null(acc_gene) && p %in% names(acc_gene)) acc_gene[[p]] else p,
        character(1))
      df$Site         <- NA_character_
      df$Modification <- NA_character_
      df$PTM_Label    <- df$Gene
      df$PTM_ID       <- df$Protein
      df
    }

  list(
    adjusted   = .annotate_ptm_view(res$adjusted,   adj_map),
    unadjusted = .annotate_ptm_view(res$unadjusted),
    protein    = .annotate_prot_view(res$protein)
  )
}

# The entry published to rv$dia_ptm_de_results. One definition, so a DE-tab run
# and an on-demand run from PTM Single cannot drift apart in shape.
ptm_de_store_entry <- function(views, res, cache, cond_test, cond_ref,
                               engine = "dia_ptm_limma_de", run_at = Sys.time()) {
  list(
    adjusted   = views$adjusted,
    unadjusted = views$unadjusted,
    protein    = views$protein,
    # Sites that fail the per-comparison valid-value filter cannot be fitted and
    # are assessed by Fisher's exact test on detection counts instead. Keyed on
    # ABSOLUTE matrix row names, since .presence_absence() runs on
    # cache$ptm_abs_matrix - consumers must join on abs_id, not on the adjusted id.
    nmar       = list(
      table            = res$nmar_table,
      has_power        = isTRUE(res$nmar_has_power),
      min_achievable_p = as.numeric(res$nmar_min_achievable_p %||% NA_real_),
      key              = "abs_id"
    ),
    comparison = paste0(as.character(cond_test), " vs ", as.character(cond_ref)),
    cond_test  = as.character(cond_test),
    cond_ref   = as.character(cond_ref),
    join_key   = "PTM_ID",
    engine     = engine,
    cache_fingerprint = ptm_cache_fingerprint(cache),
    run_at     = run_at
  )
}

# Fit and publish every contrast against ref_cond that is not already in `store`.
# Contrasts already present are NEVER refitted, so this is idempotent and cheap
# on a second call. Returns the updated store plus what ran, what failed, and the
# k-NN memo the caller should write back into rv$dia_ptm_master_cache.
# Touches no reactive value itself - the caller owns every write.
ptm_de_backfill_contrasts <- function(cache, ref_cond, test_conds, store,
                                      design_option,
                                      subject_pairing = NULL,
                                      logger   = function(...) invisible(NULL),
                                      progress = NULL) {
  .none <- list(store = store, ran = character(0), failed = character(0),
                cache_update = NULL)
  if (is.null(cache) || is.null(ref_cond) || !nzchar(as.character(ref_cond)[1]) ||
      length(test_conds) == 0L) return(.none)
  if (is.null(store)) store <- list()
  labels <- paste0(test_conds, " vs ", ref_cond)
  todo   <- test_conds[!(labels %in% names(store))]
  if (length(todo) == 0L) return(list(store = store, ran = character(0),
                                      failed = character(0), cache_update = NULL))

  .bg  <- ptm_de_batch_groups(cache, design_option)
  .bfm <- !is.null(.bg)
  .dmv <- as.integer(cache$min_valid_per_cond %||%
                       formals(dia_ptm_limma_de)$de_min_valid_per_cond)
  ran <- character(0); failed <- character(0); cache_update <- NULL

  for (tc in todo) {
    .lbl <- paste0(tc, " vs ", ref_cond)
    if (is.function(progress))
      tryCatch(progress(.lbl), error = function(e) NULL)
    res <- tryCatch(
      dia_ptm_limma_de(
        cache           = cache,
        cond_test       = tc,
        cond_ref        = ref_cond,
        de_min_valid_per_cond = .dmv,
        batch_groups    = .bg,
        batch_for_model = .bfm,
        subject_pairing = subject_pairing,
        logger          = logger
      ),
      error = function(e) list(error = conditionMessage(e)))
    if (!is.null(res$error)) {
      failed <- c(failed, .lbl)
      logger(paste0("PTM DE backfill: '", .lbl, "' failed - ", res$error))
      next
    }
    # Same k-NN memoisation the DE tab does, so later contrasts reuse one donor
    # pool. Cache-scoped and seeded, so the filled values do not depend on which
    # tab computed them first.
    if (!is.null(res$imputed_cache) && is.null(cache$ptm_abs_matrix_imputed)) {
      cache$ptm_abs_matrix_imputed <- res$imputed_cache$ptm_abs_matrix_imputed
      cache$global_matrix_imputed  <- res$imputed_cache$global_matrix_imputed
      cache$impute_stats           <- res$imputed_cache$impute_stats
      cache$ptm_impute_cap         <- res$imputed_cache$ptm_impute_cap
      cache_update                 <- res$imputed_cache
    }
    store[[.lbl]] <- ptm_de_store_entry(
      ptm_de_annotate_views(res, cache), res, cache, tc, ref_cond)
    ran <- c(ran, .lbl)
  }
  list(store = store, ran = ran, failed = failed, cache_update = cache_update)
}

# -----------------------------------------------------------------------------
# lod_impute_fc_scatter()
# Perseus-style LOD imputation on a log2 matrix, then compute condition FCs.
#
# Parameters:
#   mat            -" matrix (features × samples), log2 scale, NA = not detected
#   samp_test      -" character: sample names belonging to test condition
#   samp_ref       -" character: sample names belonging to reference condition
#   lod_downshift  -" numeric: LOD placed this many SDs below per-column mean (default 1.8)
#   lod_width      -" numeric: Gaussian noise width as fraction of per-column SD (default 0.3)
#
# Returns data.frame:
#   feature, log2fc, n_obs_test, n_total_test, n_obs_ref, n_total_ref, imp_class
#   imp_class: "none" | "partial" | "dropout_one_cond"
#
# Only features with >=1 observed value in the combined set are returned.
# LOD is per-column; columns with <2 observations fall back to global matrix stats.
# -----------------------------------------------------------------------------
lod_impute_fc_scatter <- function(mat, samp_test, samp_ref,
                                   lod_downshift = 1.8,
                                   lod_width     = 0.3,
                                   stats_mat     = NULL) {
  # stats_mat: matrix from which the per-column LOD (mean, SD) is derived.
  # Defaults to `mat`. Callers that pass a biased subset of rows -- e.g.
  # dia_nmar_fisher(), which passes only the proteins EXCLUDED from limma --
  # must supply the full matrix here. Deriving mu/SD from sparsely-detected
  # low-abundance proteins places the downshifted LOD far below the true
  # sample floor and inflates every resulting fold change. Perseus computes
  # the downshift from the full sample distribution.
  samp_test <- intersect(samp_test, colnames(mat))
  samp_ref  <- intersect(samp_ref,  colnames(mat))
  if (length(samp_test) == 0 || length(samp_ref) == 0) return(NULL)

  both_samps <- c(samp_test, samp_ref)
  sub_mat    <- mat[, both_samps, drop = FALSE]

  # Keep only features with >=1 valid value across both conditions
  has_any <- rowSums(!is.na(sub_mat)) >= 1L
  sub_mat <- sub_mat[has_any, , drop = FALSE]
  if (nrow(sub_mat) == 0L) return(NULL)

  # Global fallback stats (used when a column has <2 observations).
  # Derived from stats_mat when supplied, for the reason given above.
  .stats_src   <- if (!is.null(stats_mat)) stats_mat else sub_mat
  all_obs      <- as.vector(.stats_src[!is.na(.stats_src)])
  global_mu    <- if (length(all_obs) >= 2) mean(all_obs) else 0
  global_sd    <- if (length(all_obs) >= 2) sd(all_obs)   else 1

  # Per-column LOD imputation.
  # Seed once, outside the loop: seeding inside restarts the same stream for
  # every column, so columns with equal NA counts received identical noise
  # vectors, making the added jitter perfectly correlated across samples.
  set.seed(42L)
  imp_mat <- sub_mat
  for (j in seq_len(ncol(imp_mat))) {
    na_idx <- which(is.na(imp_mat[, j]))
    if (length(na_idx) == 0L) next
    col_nm <- colnames(imp_mat)[j]
    obs_j  <- if (!is.null(stats_mat) && col_nm %in% colnames(stats_mat))
                stats_mat[, col_nm] else imp_mat[-na_idx, j]
    obs_j  <- obs_j[!is.na(obs_j)]
    if (length(obs_j) >= 2L) {
      mu_j <- mean(obs_j)
      sd_j <- sd(obs_j)
    } else {
      mu_j <- global_mu
      sd_j <- global_sd
    }
    lod_j <- mu_j - lod_downshift * sd_j
    # Add small Gaussian noise so LOD-filled values don't collapse to a point
    imp_mat[na_idx, j] <- lod_j + rnorm(length(na_idx), 0, lod_width * max(sd_j, 1e-6))
  }

  # Condition means from imputed matrix
  mean_test <- rowMeans(imp_mat[, samp_test, drop = FALSE], na.rm = TRUE)
  mean_ref  <- rowMeans(imp_mat[, samp_ref,  drop = FALSE], na.rm = TRUE)
  log2fc    <- mean_test - mean_ref

  # Count observed (non-imputed) values per condition
  n_obs_test <- rowSums(!is.na(sub_mat[, samp_test, drop = FALSE]))
  n_obs_ref  <- rowSums(!is.na(sub_mat[, samp_ref,  drop = FALSE]))
  n_tot_test <- length(samp_test)
  n_tot_ref  <- length(samp_ref)

  imp_class <- ifelse(
    n_obs_test == 0L | n_obs_ref == 0L, "dropout_one_cond",
    ifelse(n_obs_test < n_tot_test | n_obs_ref < n_tot_ref, "partial", "none")
  )

  data.frame(
    feature      = rownames(sub_mat),
    log2fc       = log2fc,
    n_obs_test   = n_obs_test,
    n_total_test = n_tot_test,
    n_obs_ref    = n_obs_ref,
    n_total_ref  = n_tot_ref,
    imp_class    = imp_class,
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
}

# -----------------------------------------------------------------------------
# dia_nmar_fisher()
# NMAR presence/absence analysis for DIA proteins EXCLUDED from limma DE.
#
# For each protein that did not pass the DE missingness filter, tests whether
# the detection pattern across groups is non-random using Fisher's exact test
# on a 2×2 table: detected/missing × test/ref.  FC is estimated from
# LOD-imputed values (lod_impute_fc_scatter), so even 100%-absent conditions
# receive a finite, directional estimate.
#
# Parameters:
#   cache            -" dia_master_cache (must have $protein_matrix, $design_df,
#                     $acc_to_gene)
#   cond_test        -" character: test condition name
#   cond_ref         -" character: reference condition name
#   tested_proteins  -" character: protein IDs already in the limma result (excluded
#                     from this analysis)
#   cond_groups      -" named list: condition → sample names
#
# Returns data.frame with columns:
#   Protein, Gene, log2fc, n_obs_test, n_total_test, n_obs_ref, n_total_ref,
#   imp_class, fisher_p, fisher_adj_p, pct_test, pct_ref
# -----------------------------------------------------------------------------
dia_nmar_fisher <- function(cache, cond_test, cond_ref,
                             tested_proteins, cond_groups,
                             lod_downshift = 1.8, lod_width = 0.3,
                             logger = message) {
  mat         <- cache$protein_matrix
  acc_to_gene <- cache$acc_to_gene

  samp_test <- intersect(cond_groups[[cond_test]], colnames(mat))
  samp_ref  <- intersect(cond_groups[[cond_ref]],  colnames(mat))
  if (length(samp_test) == 0 || length(samp_ref) == 0) return(NULL)

  # Restrict to proteins NOT tested by limma
  excl_prots <- setdiff(rownames(mat), tested_proteins)
  if (length(excl_prots) == 0) return(NULL)

  sub_mat <- mat[excl_prots, c(samp_test, samp_ref), drop = FALSE]

  # Remove proteins with zero observations in both groups (truly absent everywhere)
  has_any <- rowSums(!is.na(sub_mat)) >= 1L
  sub_mat <- sub_mat[has_any, , drop = FALSE]
  if (nrow(sub_mat) == 0L) return(NULL)

  # LOD-imputed FC estimates
  fc_df <- lod_impute_fc_scatter(mat[rownames(sub_mat), , drop = FALSE],
                                  samp_test, samp_ref,
                                  lod_downshift = lod_downshift,
                                  lod_width     = lod_width,
                                  stats_mat     = mat[, c(samp_test, samp_ref), drop = FALSE])
  if (is.null(fc_df) || nrow(fc_df) == 0) return(NULL)

  # Fisher's exact test per protein: 2x2 table [detected, missing] x [test, ref]
  n_test <- length(samp_test)
  n_ref  <- length(samp_ref)

  # Minimum achievable two-sided Fisher p for this design. When BH-adjusted
  # significance is unreachable, the p-values are not interpretable and are
  # suppressed rather than reported as non-significant.
  .min_achievable_p <- 2 / choose(n_test + n_ref, n_test)
  .fisher_has_power <- .min_achievable_p <= 0.05
  if (!.fisher_has_power)
    logger(paste0("[NMAR] Design has no resolution for Fisher's exact test ",
                  "(minimum achievable p = ", signif(.min_achievable_p, 3),
                  "). P-values suppressed; results presented as a ranked list."))

  fisher_p <- vapply(fc_df$feature, function(prot) {
    obs_t <- sum(!is.na(sub_mat[prot, samp_test]))
    obs_r <- sum(!is.na(sub_mat[prot, samp_ref]))
    tbl <- matrix(c(obs_t,       n_test - obs_t,
                    obs_r,       n_ref  - obs_r),
                  nrow = 2,
                  dimnames = list(c("detected", "missing"),
                                  c("test", "ref")))
    # Suppress warnings for tables where both marginals force outcome
    suppressWarnings(fisher.test(tbl)$p.value)
  }, numeric(1))

  fisher_adj_p <- p.adjust(fisher_p, method = "BH")

  # Honour the stated behaviour: when the design cannot reach significance, the
  # p-values are not interpretable. Blank them rather than returning numbers
  # that will be read as non-significant findings. Column names are preserved so
  # existing consumers keep working; results remain a ranked descriptive list
  # via log2fc, pct_test and pct_ref.
  if (!.fisher_has_power) {
    fisher_p     <- rep(NA_real_, length(fisher_p))
    fisher_adj_p <- rep(NA_real_, length(fisher_adj_p))
  }

  genes <- acc_to_gene[fc_df$feature]
  genes[is.na(genes)] <- fc_df$feature[is.na(genes)]

  data.frame(
    Protein      = fc_df$feature,
    Gene         = genes,
    log2fc       = fc_df$log2fc,
    n_obs_test   = fc_df$n_obs_test,
    n_total_test = fc_df$n_total_test,
    n_obs_ref    = fc_df$n_obs_ref,
    n_total_ref  = fc_df$n_total_ref,
    pct_test     = round(fc_df$n_obs_test / fc_df$n_total_test * 100, 1),
    pct_ref      = round(fc_df$n_obs_ref  / fc_df$n_total_ref  * 100, 1),
    imp_class    = fc_df$imp_class,
    fisher_p     = fisher_p,
    fisher_adj_p = fisher_adj_p,
    min_achievable_p = .min_achievable_p,
    fisher_has_power = .fisher_has_power,
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
}

# ===========================================================================
# PRERUN EXPERIMENTS  -" shared helpers
# ===========================================================================

# PRERUN_CACHE_DIR  -" server-side directory where preview sidecars are written
# at export time. The Prerun Experiments module defaults to scanning this folder.
PRERUN_CACHE_DIR <- normalizePath(
  file.path(getwd(), "prerun_cache"),
  winslash = "/", mustWork = FALSE
)
if (!dir.exists(PRERUN_CACHE_DIR))
  dir.create(PRERUN_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# make_pca_thumbnail()
# Generates a compact PCA scatter from a summarization bundle.
# Returns a base64 "data:image/png;base64,..." string or NULL on failure.
# Called at export time when the cache is already in memory (no extra cost).
# Requires: ggplot2, base64enc, tidyr.
# ---------------------------------------------------------------------------
make_pca_thumbnail <- function(bundle, width_px = 420L, height_px = 310L) {
  tryCatch({
    type      <- bundle$type %||% "unknown"
    mat       <- NULL
    group_map <- NULL   # named list: group_name -> character vector of sample names
    .is_dia_ptm_route <- FALSE   # set below; selects the complete-case PCA path

    # ── 1. Extract sample × feature matrix + group map ──────────────────────
    if (identical(type, "lf_summarization")) {
      # ProteinLevelData lives inside $processed (MSstats/SAINT shim path)
      pld <- bundle$cache$processed$ProteinLevelData %||% bundle$cache$ProteinLevelData
      if (is.null(pld)) return(NULL)
      pld_df <- as.data.frame(pld)
      pld_df <- pld_df[!is.na(pld_df$LogIntensities), ]
      wide <- tryCatch(
        tidyr::pivot_wider(
          pld_df[, c("Protein", "originalRUN", "LogIntensities")],
          names_from  = "originalRUN",
          values_from = "LogIntensities",
          values_fn   = mean
        ), error = function(e) NULL)
      if (is.null(wide) || nrow(wide) == 0L) return(NULL)
      mat <- as.matrix(wide[, -1L, drop = FALSE])
      rownames(mat) <- wide[[1L]]
      # Build group_map directly from ProteinLevelData so run IDs always match.
      # cond_groups stores the original PD column names (e.g. "Abundance: F2: ..."),
      # which differ from the renamed originalRUN values (e.g. "GFP_Rep1").
      if ("GROUP" %in% colnames(pld_df) && "originalRUN" %in% colnames(pld_df)) {
        rg        <- unique(pld_df[!is.na(pld_df$GROUP), c("originalRUN", "GROUP")])
        group_map <- split(rg$originalRUN, rg$GROUP)
      } else {
        group_map <- bundle$cond_groups %||% bundle$meta$groups
      }

    } else if (identical(type, "dia_summarization")) {
      mat <- bundle$cache$protein_matrix
      group_map <- bundle$cond_groups

    } else if (identical(type, "tmt_summarization")) {
      if (!is.null(bundle$limma_cache$protein_matrix)) {
        mat <- bundle$limma_cache$protein_matrix
      } else if (!is.null(bundle$cache$ProteinLevelData)) {
        pld_df <- as.data.frame(bundle$cache$ProteinLevelData)
        pld_df <- pld_df[!is.na(pld_df$Abundance), ]
        wide <- tryCatch(
          tidyr::pivot_wider(
            pld_df[, c("Protein", "Run", "Abundance")],
            names_from  = "Run",
            values_from = "Abundance",
            values_fn   = mean
          ), error = function(e) NULL)
        if (!is.null(wide) && nrow(wide) > 0L) {
          mat <- as.matrix(wide[, -1L, drop = FALSE])
          rownames(mat) <- wide[[1L]]
        }
      }
      ann <- bundle$annotation_map; conds <- bundle$conditions
      if (!is.null(ann) && !is.null(conds)) {
        # TMT annotation_map uses ColumnName (not Run) for the matrix column IDs
        group_map <- lapply(setNames(conds, conds),
                            function(cond) ann$ColumnName[ann$Condition == cond])
        # Apply plex batch correction when ≥2 plexes (mirrors server_pca.R behaviour)
        if (!is.null(mat) && "Mixture" %in% colnames(ann)) {
          mixture_vec <- ann$Mixture[match(colnames(mat), ann$ColumnName)]
          n_plexes    <- length(unique(mixture_vec[!is.na(mixture_vec)]))
          if (n_plexes >= 2L && requireNamespace("limma", quietly = TRUE)) {
            sample_grps <- ann$Condition[match(colnames(mat), ann$ColumnName)]
            mat <- tryCatch(
              limma::removeBatchEffect(
                mat,
                batch  = factor(mixture_vec),
                design = model.matrix(~ factor(sample_grps))
              ),
              error = function(e) mat
            )
          }
        }
      }

    } else if (identical(type, "ptm_summarization")) {
      if ((identical(bundle$route, "dia_ptm") ||
           identical(bundle$route, "tmt_ptm")) &&
          !is.null(bundle$dia_ptm_cache$ptm_abs_matrix)) {
        # LF/TMT PTM path  -" absolute log2 matrix already in cache
        mat       <- bundle$dia_ptm_cache$ptm_abs_matrix
        group_map <- bundle$cond_groups
        # dia_ptm only  -" NOT tmt_ptm, whose thumbnail stays on the shared path.
        .is_dia_ptm_route <- identical(bundle$route, "dia_ptm")
        # Apply plex batch correction for TMT PTM with >=2 plexes
        if (identical(bundle$route, "tmt_ptm")) {
          ddf <- bundle$dia_ptm_cache$design_df
          if (!is.null(ddf) && "mixture" %in% colnames(ddf)) {
            mix_map <- setNames(as.character(ddf$mixture),
                                as.character(ddf$sample_name))
            mix_vec <- mix_map[colnames(mat)]
            n_plx   <- length(unique(mix_vec[!is.na(mix_vec)]))
            if (n_plx >= 2L && requireNamespace("limma", quietly = TRUE)) {
              cond_map2 <- setNames(as.character(ddf$condition),
                                    as.character(ddf$sample_name))
              cond_vec  <- cond_map2[colnames(mat)]
              cond_vec[is.na(cond_vec)] <- "Unknown"
              mat <- tryCatch(
                limma::removeBatchEffect(
                  mat,
                  batch  = factor(mix_vec),
                  design = model.matrix(~ factor(cond_vec))
                ),
                error = function(e) mat
              )
            }
          }
        }
      } else {
        # MSstatsPTM (TMT) path
        ptm_int <- bundle$cache$PTMIntensities
        if (!is.null(ptm_int)) {
          pld_df   <- as.data.frame(ptm_int)
          pld_df   <- pld_df[!is.na(pld_df$Abundance), ]
          feat_col <- if ("PTM" %in% colnames(pld_df)) "PTM" else colnames(pld_df)[1L]
          wide <- tryCatch(
            tidyr::pivot_wider(
              pld_df[, c(feat_col, "Run", "Abundance")],
              names_from  = "Run",
              values_from = "Abundance",
              values_fn   = mean
            ), error = function(e) NULL)
          if (!is.null(wide) && nrow(wide) > 0L) {
            mat <- as.matrix(wide[, -1L, drop = FALSE])
            rownames(mat) <- wide[[1L]]
          }
        }
        ann <- bundle$annotation_map; conds <- bundle$conditions
        if (!is.null(ann) && !is.null(conds))
          group_map <- lapply(setNames(conds, conds),
                              function(cond) ann$Run[ann$Condition == cond])
      }
    }

    if (is.null(mat) || nrow(mat) < 3L || ncol(mat) < 3L) return(NULL)

    # ── 2. SAINT/IP complete cases; else filter (>=50% non-NA) + row-median imputation ───────────────
    # The property that matters is "this is an IP/AP-MS experiment, where
    # absence is signal", not "this is SAINT". meta$pipeline carries that flag
    # ("ip" | "global"); pipeline_label must not be used for it, because two of
    # its values begin with "Global" on runs that are IP experiments. The two
    # original SAINT tests are OR-ed in so that route cannot regress.
    .is_saint_ip <- identical(type, "lf_summarization") &&
      (identical(bundle$meta$de_engine %||% "", "saint") ||
       identical(bundle$cache$pipeline_label %||% "", "IP (SAINTexpress)"))
    .is_ip_route <- identical(type, "lf_summarization") &&
      (identical(bundle$meta$pipeline %||% "", "ip") || .is_saint_ip)
    if (.is_ip_route) {
      mat_cc <- mat[stats::complete.cases(mat), , drop = FALSE]
      message("make_pca_thumbnail: IP complete-case PCA on ",
              nrow(mat_cc), " proteins (no imputation)")
      if (nrow(mat_cc) < 50L) {
        # SAINT keeps its existing "no thumbnail" outcome. For a non-SAINT IP
        # run, no plot at all is worse than the imputed one, so fall through to
        # the >=50% + row-median path. This is the one branch where the QC
        # thumbnail still contains fabricated values, so it is logged
        # explicitly rather than flagged on the return value: the return is a
        # scalar data-URI and every consumer of it lives outside this file.
        if (.is_saint_ip) return(NULL)
        message("make_pca_thumbnail: IMPUTED THUMBNAIL \u2014 only ", nrow(mat_cc),
                " complete case(s), below the 50 required; fell back to the ",
                ">=50% non-NA filter plus per-protein row-median imputation. ",
                "This QC PCA contains fabricated values and must not be read ",
                "as evidence of group separation.")
      } else {
        mat <- mat_cc
      }
    }
    if (.is_dia_ptm_route) {
      # DIA PTM: complete cases only. The >=50% + row-median path imputes
      # group-blind, which pulls every sample toward the row centre and flatters
      # group separation on a matrix that is heavily missing. Mirrors the
      # decision already taken for the SAINT and ProDA routes.
      mat <- mat[stats::complete.cases(mat), , drop = FALSE]
      message("make_pca_thumbnail: DIA PTM complete-case PCA on ", nrow(mat),
              " sites (no imputation)")
      if (nrow(mat) < 3L) return(NULL)
    } else {
      keep <- rowMeans(!is.na(mat)) >= 0.5
      mat  <- mat[keep, , drop = FALSE]
      if (nrow(mat) < 3L) return(NULL)
      for (i in seq_len(nrow(mat))) {
        nas <- is.na(mat[i, ])
        if (any(nas)) mat[i, nas] <- median(mat[i, !nas], na.rm = TRUE)
      }
    }
    if (anyNA(mat)) return(NULL)

    # ── 3. PCA ──────────────────────────────────────────────────────────────
    pca_res <- prcomp(t(mat), center = TRUE, scale. = FALSE)
    scores  <- as.data.frame(pca_res$x[, 1:2])
    scores$Sample <- rownames(scores)
    var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

    scores$Group <- if (!is.null(group_map) && length(group_map) > 0L) {
      g <- rep("Other", nrow(scores))
      for (grp in names(group_map))
        g[scores$Sample %in% group_map[[grp]]] <- grp
      g
    } else "All Samples"

    # ── 4. ggplot2 scatter ──────────────────────────────────────────────────
    p <- ggplot2::ggplot(scores,
           ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::labs(x = paste0("PC1 (", var_exp[1], "%)"),
                    y = paste0("PC2 (", var_exp[2], "%)")) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        legend.position  = "bottom",
        legend.title     = ggplot2::element_blank(),
        legend.text      = ggplot2::element_text(size = 7),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin      = ggplot2::margin(4, 4, 4, 4)
      )

    # ── 5. Base64-encode PNG ────────────────────────────────────────────────
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp), add = TRUE)
    grDevices::png(tmp, width = width_px, height = height_px, res = 96)
    print(p)
    grDevices::dev.off()
    paste0("data:image/png;base64,", base64enc::base64encode(tmp))

  }, error = function(e) {
    message("make_pca_thumbnail: ", e$message)
    NULL
  })
}

# ---------------------------------------------------------------------------
# count_bundle_features()
# Returns the number of proteins (or PTM sites) in a summarization bundle.
# ---------------------------------------------------------------------------
count_bundle_features <- function(bundle) {
  tryCatch({
    type <- bundle$type %||% ""
    if (identical(type, "lf_summarization"))
      # ProteinLevelData lives inside $processed (MSstats/SAINT shim path)
      length(unique((bundle$cache$processed$ProteinLevelData %||%
                     bundle$cache$ProteinLevelData)$Protein))
    else if (identical(type, "dia_summarization"))
      nrow(bundle$cache$protein_matrix)
    else if (identical(type, "tmt_summarization")) {
      if (!is.null(bundle$limma_cache$protein_matrix))
        nrow(bundle$limma_cache$protein_matrix)
      else
        length(unique(bundle$cache$ProteinLevelData$Protein))
    } else if (identical(type, "ptm_summarization")) {
      # n_sites_output is the POST-SUMMARIZATION site count: it is taken after the
      # summarization missingness filter but BEFORE the per-contrast DE missingness
      # filter and BEFORE shared-peptide expansion.  It is a bundle/preview feature
      # count and can never equal the DE panel's n_tested_* figures, which are a
      # different quantity by construction.  Do not "reconcile" the two.
      if (identical(bundle$route, "dia_ptm"))
        bundle$dia_ptm_cache$n_sites_output %||%
          nrow(bundle$dia_ptm_cache$ptm_adj_matrix) %||% NA_integer_
      else if (identical(bundle$route, "tmt_ptm"))
        bundle$dia_ptm_cache$n_sites_output %||%
          nrow(bundle$dia_ptm_cache$ptm_adj_matrix) %||% NA_integer_
      else
        nrow(bundle$cache$PTMIntensities)
    } else
      NA_integer_
  }, error = function(e) NA_integer_)
}

# ---------------------------------------------------------------------------
# save_prerun_preview()
# Writes a lightweight .preview.rds sidecar to PRERUN_CACHE_DIR (server-side
# cache) AND returns the preview list invisibly so callers can also write it
# to a second download (e.g. as a ZIP companion or a separate downloadButton).
# ---------------------------------------------------------------------------
save_prerun_preview <- function(bundle, download_filename, logger_fn = message) {
  # Build the preview object first  -" independent of whether the cache save succeeds.
  preview <- list(
    type             = bundle$type,
    meta             = bundle$meta,
    route            = bundle$route,
    pca_thumbnail    = bundle$pca_thumbnail,
    n_features       = bundle$n_features,
    cond_groups      = bundle$cond_groups,
    conditions       = bundle$conditions,
    annotation_map   = bundle$annotation_map,
    data_format      = bundle$data_format,
    design_option    = bundle$design_option,
    .source_filename = download_filename
  )
  # Attempt to persist to the server-side cache (used when running locally).
  # Failure here must NOT prevent the preview from being returned to the caller
  # (who writes it into the download ZIP).
  tryCatch({
    preview_fname <- sub("\\.rds$", ".preview.rds",
                         download_filename, ignore.case = TRUE)
    saveRDS(preview, file.path(PRERUN_CACHE_DIR, preview_fname))
    logger_fn(paste0("Prerun preview saved: ", preview_fname))
  }, error = function(e) {
    logger_fn(paste0("Warning: could not save prerun preview to cache \u2014 ", e$message))
  })
  invisible(preview)
}
