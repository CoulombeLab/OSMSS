#!/usr/bin/env Rscript
# app.R ========================================================================
# ONE STOP MS SHOP (OSMSS) - Modular Launcher
# Integrates Residue PTM Prediction + Local MS (FragPipe) Analysis

options(keep.source = TRUE)
# Sets working directory to source file location if running in RStudio
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Increase upload size limit for large proteomics files (1 GB, e.g. MaxQuant evidence.txt).
# Shiny's limit is app-wide; global.R sets the same value (it is sourced after this line).
options(shiny.maxRequestSize = 1024 * 1024^2)

library(shiny)

# Source components into the global environment so all modules
# are visible regardless of how the app is launched (runApp vs Rscript)
source("global.R",                    local = .GlobalEnv)

# Log all module source operations so load failures are traceable
.source_module <- function(file) {
  tryCatch({
    source(file, local = .GlobalEnv)
    osmss_file_log(paste0("Loaded module: ", file), level = "INFO", module = "Startup")
  }, error = function(e) {
    osmss_file_log(paste0("FAILED to load module: ", file, " — ", e$message),
                   level = "ERROR", module = "Startup")
    stop(e)
  })
}

.source_module("server_residue_prediction.R")
.source_module("server_local_ms.R")
.source_module("server_de.R")
.source_module("server_pca.R")
.source_module("server_heatmap.R")
.source_module("server_rank_fc.R")
.source_module("server_trends.R")
.source_module("server_trend_clusters.R")
.source_module("server_fragpipe_db.R")
.source_module("server_tmt_upload.R")
.source_module("server_ptm_de.R")
.source_module("server_ptm_trends.R")
.source_module("server_ptm_single.R")
.source_module("server_protein_coverage.R")
.source_module("server_gsea.R")
.source_module("ui.R")
.source_module("server.R")
if (file.exists("Troubleshoot/osmss_tracer.R")) {
  source("Troubleshoot/osmss_tracer.R",       local = .GlobalEnv)
  source("Troubleshoot/osmss_tracer_hooks.R", local = .GlobalEnv)
  osmss_trace_auto()
}

osmss_file_log("All modules loaded — launching app", level = "INFO", module = "Startup")

.server <- if (exists(".osmss") && isTRUE(.osmss$active)) osmss_wrap_server(server) else server
shinyApp(ui = ui, server = .server, options = list(launch.browser = TRUE))
