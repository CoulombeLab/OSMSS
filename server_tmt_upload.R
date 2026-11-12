# server_tmt_upload.R
# ==============================================================================
# TMT Upload Wizard Module
# Extracted from server_tmt_volcano.R (Phase 1 of OSMSS v4→v5 refactor).
#
# Responsibilities:
#   - File upload and parsing (PD nested + flat peptide paths)
#   - Experimental design grid (via tmt_register_design_grid)
#   - TMT Setup Wizard (steps 1-4, 6: program → data level → upload → plex → summary)
#   - Summarization dispatch (MSstatsTMT, limma_qf, limma_protein)
#   - Export / import of summarization bundles (.rds)
#
# On wizard completion, writes:
#   rv$data_source            <- "tmt"
#   rv$tmt_wizard_completed   <- TRUE
#   rv$tmt_analysis_route
#   rv$tmt_wizard_search_program
#   rv$tmt_wizard_data_level
#   rv$tmt_annotation_map
#   rv$tmt_conditions
#   rv$tmt_design_option
#   (and clears LF DDA + DIA caches — last-upload-wins)
# ==============================================================================

server_tmt_upload <- function(input, output, session, rv, logger) {

  # Module-local reactive values (upload / wizard / design only; no DE fields)
  tmt_upl_rv <- reactiveValues(
    # Per-file parsed data
    files          = list(),
    n_files        = 0,
    all_channels   = character(0),
    # Study design
    design_option  = NULL,
    annotation_map = NULL,
    # Channels the user marked "Exclude" in the design grid (absent from
    # annotation_map); carried into the summarization cache for the methods text.
    excluded_channels = character(0),
    conditions     = NULL,
    n_plexes       = 1L,
    design_conds   = character(0),
    has_bridge     = FALSE,
    # Analysis routing
    analysis_route = NULL,          # "limma_protein" | "limma_qf"
    limma_cache    = NULL,
    # Status
    status         = "Upload TMT CSV file(s) to begin.",
    # Wizard state
    wizard_completed       = FALSE,
    wizard_modal_open      = FALSE,
    wizard_step            = 1L,
    wizard_search_program  = "pd",
    wizard_data_level      = "peptide",
    wizard_data_levels_raw = NULL,
    # Shared-peptide handling for the peptide route. NA_character_ on the
    # protein route, where no peptides exist and the setting is meaningless.
    shared_peptides        = "unique_only",
    # Single plex: did the user report the abundances as already normalized across
    # channels (step 6)? FALSE = OSMSS median-centers (default). pending_norm holds
    # the answer and the data check for the summarization in flight.
    input_normalized       = FALSE,
    pending_norm           = NULL,
    wizard_reachable_steps = 1:6,
    # Flat peptide upload state
    flat_dfs        = list(),
    flat_names      = character(0),
    flat_n_files    = 0L,
    flat_mappings   = list(),
    flat_current_idx = 1L,
    # Design grid transient state
    grid_assignments = NULL,
    paint_cond       = NULL
  )
  # A TMT upload resets the experiment itself (.do_finalize_upload), so the
  # wizard's own progress and its pending flat-peptide input survive the reset.
  register_experiment_reset(rv, tmt_upl_rv, keep = c(
    "wizard_modal_open", "wizard_step", "wizard_reachable_steps",
    "wizard_search_program", "wizard_data_level", "wizard_data_levels_raw",
    "shared_peptides",
    "flat_dfs", "flat_names", "flat_n_files", "flat_mappings", "flat_current_idx"))

  # TMT summarization async task -----------------------------------------------
  tmt_wiz_task <- ExtendedTask$new(function(
      route, files, annotation_map, design_opt, wizard_data_level, rv_snapshot,
      shared_peptides, apply_centering = TRUE) {
    promises::future_promise({
      compute_tmt_summarization(
        route             = route,
        files             = files,
        annotation_map    = annotation_map,
        design_opt        = design_opt,
        wizard_data_level = wizard_data_level,
        shared_peptides   = shared_peptides,
        rv_snapshot       = rv_snapshot,
        apply_centering   = apply_centering
      )
    }, seed = TRUE)
  })
  tmt_wiz_cancelled <- reactiveVal(FALSE)

  # TMT wizard summarization result observer -----------------------------------
  observeEvent(tmt_wiz_task$result(), {
    hide_async_overlay(session)
    if (isTRUE(tmt_wiz_cancelled())) {
      tmt_wiz_cancelled(FALSE)
      return()
    }
    res <- tmt_wiz_task$result()
    for (msg in res$messages %||% character(0)) {
      showNotification(msg, type = "warning", duration = 8)
    }
    if (!isTRUE(res$success)) {
      showNotification(paste0("TMT summarization failed: ", res$error),
                       type = "error", duration = 15)
      return()
    }
    result <- res$result
    route  <- res$route
    # Channel-normalization answer and data check (step 6), carried in the cache
    # so the DE provenance and the methods text can report them.
    .pn <- tmt_upl_rv$pending_norm
    if (!is.null(.pn)) for (.k in names(.pn)) result[[.k]] <- .pn[[.k]]
    # Channels left out at the design step (not in annotation_map): reported in
    # the methods text, and kept with the cache so an exported bundle carries them.
    result$excluded_channels <- tmt_upl_rv$excluded_channels %||% character(0)
    tmt_upl_rv$limma_cache   <- result
    rv$tmt_limma_cache       <- result
    tmt_upl_rv$conditions    <- result$conditions
    rv$tmt_conditions        <- result$conditions

    # Populate rv$tmt_master_cache (ProteinLevelData long format) so that the
    # Trend Clusters, Heatmap, and PCA modules — which still use the legacy
    # MSstatsTMT-style structure — work without modification.
    {
      ann_mc     <- result$annotation_map
      pm_mc      <- result$protein_matrix
      bio_ann_mc <- ann_mc[!is.na(ann_mc$Condition) &
                             !ann_mc$Condition %in% "Norm", , drop = FALSE]
      shared_mc  <- intersect(bio_ann_mc$ColumnName, colnames(pm_mc))
      if (length(shared_mc) > 0) {
        pld_mc <- do.call(rbind, lapply(shared_mc, function(ch) {
          row_mc <- bio_ann_mc[bio_ann_mc$ColumnName == ch, , drop = FALSE]
          data.frame(
            Protein      = rownames(pm_mc),
            Condition    = row_mc$Condition,
            BioReplicate = row_mc$BioReplicate,
            Abundance    = pm_mc[, ch],
            Run          = row_mc$Mixture,   # plex ID — used for batch correction
            stringsAsFactors = FALSE
          )
        }))
        rv$tmt_master_cache        <- list(ProteinLevelData = pld_mc)
        rv$tmt_master_cache_design <- result$design_option
        rv$tmt_master_cache_ann    <- ann_mc
      }
    }
    label <- switch(route,
      limma_protein = "Protein-level",
      limma_qf      = "QFeatures",
      route
    )
    showNotification(paste0(label, " summarization complete!"),
                     type = "message", duration = 5)
    logger(paste0("TMT Wizard: ", route, " summarization stored in cache."))
  })

  # TMT upload summarization result observer (removed — MSstatsTMT sidebar path deprecated)
  # This observer was bound to the legacy tmt_upload_task (MSstatsTMT-only).
  # The wizard tmt_wiz_task observer above handles all current routes.

  # TMT cancel observers -------------------------------------------------------
  observeEvent(rv$async_cancel_requested, {
    if (isTRUE(rv$async_cancel_requested)) {
      if (identical(tmt_wiz_task$status(), "running")) {
        tmt_wiz_cancelled(TRUE)
        rv$async_cancel_requested <- FALSE
      }
    }
  })

  # =========================================================================
  # HELPERS
  # =========================================================================
  .short_label <- tmt_short_label
  .extract_tag <- tmt_extract_tag

  # Enable/disable tmt_run_btn (still present in UI during Phase 1).
  # Will be removed when DE is unified into server_de.R (Phase 3d).
  .enable_run_btn <- function() {
    shinyjs::enable("tmt_run_btn")
    shinyjs::removeClass("tmt_run_btn", "btn-secondary")
    shinyjs::addClass("tmt_run_btn", "btn-success")
  }

  .disable_run_btn <- function() {
    shinyjs::disable("tmt_run_btn")
    shinyjs::removeClass("tmt_run_btn", "btn-success")
    shinyjs::addClass("tmt_run_btn", "btn-secondary")
  }

  # =========================================================================
  # FINALIZE UPLOAD (shared by PD and flat-peptide paths)
  # =========================================================================
  .do_finalize_upload <- function(parsed_files) {
    # New data replaces the whole experiment (all data types, every module's
    # results) before any of it is written. Nothing is reset if no file parsed.
    if (length(parsed_files) > 0) {
      rv$reset_experiment_state(keep = RESET_KEEP_ON_UPLOAD)
      logger("TMT Upload: New files uploaded - previous experiment state reset")
    }

    # One file per plex: keep channel names distinct across files (see
    # tmt_disambiguate_plex_files()); a no-op for one file or distinct names.
    parsed_files <- tmt_disambiguate_plex_files(parsed_files, logger)
    .renamed <- Filter(Negate(is.null), lapply(parsed_files, `[[`, "plex_id_renamed"))
    if (length(.renamed) > 0)
      showNotification(paste0("The files use the same channel names, so the file IDs of ",
                              length(.renamed), " file(s) were renumbered to keep every ",
                              "plex's channels separate (e.g. ",
                              names(.renamed[[1]])[1], " -> ", .renamed[[1]][[1]], ")."),
                       type = "message", duration = 8)

    tmt_upl_rv$files   <- parsed_files
    tmt_upl_rv$n_files <- length(parsed_files)

    if (length(parsed_files) > 0) {
      f1 <- parsed_files[[1]]
      rv$pd_protein_headers  <- f1$headers
      rv$pd_raw_lines        <- f1$all_lines
      rv$pd_protein_line_idx <- f1$protein_idx
      rv$pd_ratio_cols       <- f1$ratio_cols
      rv$pd_pval_cols        <- f1$pval_cols
      rv$pd_other_cols       <- setdiff(f1$headers, c(f1$ratio_cols, f1$pval_cols))
      rv$pd_filename         <- f1$name
      rv$pd_gene_col         <- NULL

      rv$tmt_files        <- parsed_files
      rv$tmt_n_files      <- length(parsed_files)
      rv$tmt_files_source <- "volcano"

      all_ch <- unique(unlist(lapply(parsed_files, function(f) f$abund_cols)))
      tmt_upl_rv$all_channels <- all_ch
      tmt_upl_rv$status <- paste0(length(parsed_files), " file(s) loaded, ",
                                  length(all_ch), " channels detected.")

      all_conds <- unique(unlist(lapply(parsed_files, function(f) f$channel_info$Condition)))
      tmt_upl_rv$conditions <- sort(all_conds)
      rv$tmt_conditions <- tmt_upl_rv$conditions

      gd <- detect_gene_column(f1$headers, parsed_file = f1, rv = rv, logger = logger)
      if (isTRUE(tmt_upl_rv$wizard_modal_open)) {
        logger("TMT Upload: Wizard is open; skipping automatic design modal.")
      } else if (gd$source == "column") {
        rv$pd_has_gn   <- FALSE
        rv$pd_gene_col <- gd$column
        .show_arch_fork_modal()
      } else if (gd$source == "gn") {
        rv$pd_has_gn <- TRUE
        .show_arch_fork_modal()
      } else {
        rv$pd_has_gn <- FALSE
        if (is.null(rv$pd_gene_col)) {
          showModal(build_gene_col_modal("pd_gene_col_radio",
                                         "tmt_gene_col_confirm", f1$headers))
        } else {
          .show_arch_fork_modal()
        }
      }
    }
  }

  # =========================================================================
  # FILE UPLOAD HANDLER (multi-file)
  # =========================================================================
  observeEvent(input$tmt_upload_csv_files, {
    req(input$tmt_upload_csv_files)

    file_info <- input$tmt_upload_csv_files
    n <- nrow(file_info)
    logger(paste0("TMT Upload: ", n, " file(s) uploaded"))

    # ── Flat peptide branch ──
    if (isTRUE(input$tmt_upload_input_format == "flat_peptide")) {
      tmt_upl_rv$flat_dfs <- lapply(seq_len(n), function(i)
        readr::read_csv(file_info$datapath[i],
                        col_types = readr::cols(.default = "c"),
                        show_col_types = FALSE)
      )
      tmt_upl_rv$flat_names       <- file_info$name
      tmt_upl_rv$flat_n_files     <- n
      tmt_upl_rv$flat_mappings    <- vector("list", n)
      tmt_upl_rv$flat_current_idx <- 1L
      flat_tmt_show_mapping_modal(session, colnames(tmt_upl_rv$flat_dfs[[1]]),
                                  "tmt_flat_",
                                  file_name = tmt_upl_rv$flat_names[1],
                                  file_idx  = 1L,
                                  n_files   = n)
      return()
    }

    parsed_files <- list()

    withProgress(message = "Parsing TMT CSV file(s)...", value = 0, {
      for (i in seq_len(n)) {
        incProgress(1 / n, detail = paste0("File ", i, " of ", n, ": ", file_info$name[i]))

        path <- file_info$datapath[i]
        name <- file_info$name[i]

        all_lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) {
          logger(paste0("!!! [TMT_Upload] File read error (", name, "): ", e$message))
          NULL
        })
        if (is.null(all_lines) || length(all_lines) < 2) {
          showNotification(paste0("Could not read file: ", name), type = "error")
          next
        }

        all_lines <- collapse_multiline_csv(all_lines)
        parsed_files[[i]] <- parse_tmt_file_structure(all_lines, name)

        pf <- parsed_files[[i]]
        logger(paste0("TMT Upload: '", name, "' - ", length(pf$protein_idx), " proteins, ",
                      length(pf$pep_data_idx), " peptides, ", length(pf$abund_cols),
                      " channels (", paste(pf$plex_prefix, collapse = ","), ")"))
      }
    })

    .do_finalize_upload(parsed_files)
  })

  # ── Flat peptide: return to the wizard ──
  # The column-mapping dialog is a modal, so opening it REPLACES the wizard
  # modal (Shiny shows one modal at a time). The PD route uploads inside the
  # wizard and never leaves it; the flat route did, and nothing brought the
  # wizard back: after "Converting flat peptide files..." the screen was empty
  # (.do_finalize_upload deliberately shows nothing while the wizard is open).
  # Reopen it at step 3, where the upload status shows the loaded file(s).
  .reopen_wizard_step3 <- function() {
    if (isTRUE(tmt_upl_rv$wizard_modal_open)) {
      tmt_upl_rv$wizard_step <- 3L
      showModal(.tmt_wiz_build_modal())
    }
  }

  # Cancel in the mapping dialog: discard the pending flat upload, back to step 3
  observeEvent(input$tmt_flat_flat_cancel, {
    removeModal()
    tmt_upl_rv$flat_dfs         <- list()
    tmt_upl_rv$flat_mappings    <- list()
    tmt_upl_rv$flat_n_files     <- 0L
    tmt_upl_rv$flat_current_idx <- 1L
    logger("TMT Upload (flat): column mapping cancelled")
    .reopen_wizard_step3()
  }, ignoreInit = TRUE)

  # ── Flat peptide: column-mapping confirm observer ──
  observeEvent(input$tmt_flat_flat_confirm, {
    idx <- if (is.null(tmt_upl_rv$flat_current_idx)) 1L else tmt_upl_rv$flat_current_idx
    mapping <- flat_tmt_collect_mapping(input, "tmt_flat_", file_idx = idx)

    if (is.null(mapping$abundance_cols) || length(mapping$abundance_cols) == 0) {
      showNotification("Select at least one abundance column.", type = "warning")
      return()
    }

    tmt_upl_rv$flat_mappings[[idx]] <- mapping
    removeModal()

    n_f <- if (is.null(tmt_upl_rv$flat_n_files)) 1L else tmt_upl_rv$flat_n_files
    if (isTRUE(input$tmt_flat_apply_all) && idx == 1L && n_f > 1L) {
      for (j in seq(2L, n_f)) tmt_upl_rv$flat_mappings[[j]] <- mapping
      idx <- n_f
    }

    if (idx < n_f) {
      next_idx <- idx + 1L
      tmt_upl_rv$flat_current_idx <- next_idx
      flat_tmt_show_mapping_modal(session, colnames(tmt_upl_rv$flat_dfs[[next_idx]]),
                                  "tmt_flat_",
                                  file_name = tmt_upl_rv$flat_names[next_idx],
                                  file_idx  = next_idx,
                                  n_files   = n_f)
    } else {
      parsed_files <- vector("list", n_f)
      ok <- TRUE
      withProgress(message = "Converting flat peptide files...", value = 0, {
        for (i in seq_len(n_f)) {
          incProgress(1 / n_f, detail = paste0("File ", i, " of ", n_f))
          mock_lines <- tryCatch(
            flat_tmt_build_mock_pd(tmt_upl_rv$flat_dfs[[i]], tmt_upl_rv$flat_mappings[[i]],
                                   plex_idx = i),
            error = function(e) {
              logger(paste0("!!! TMT Upload flat build error (file ", i, "): ", e$message))
              NULL
            }
          )
          if (is.null(mock_lines)) {
            showNotification(paste0("Error converting file ", i,
                                    ": check column mappings."), type = "error")
            ok <- FALSE; break
          }
          parsed_files[[i]] <- parse_tmt_file_structure(mock_lines, tmt_upl_rv$flat_names[i])
          pf <- parsed_files[[i]]
          logger(paste0("TMT Upload (flat): '", tmt_upl_rv$flat_names[i], "' — ",
                        length(pf$abund_cols), " channels"))
        }
      })
      if (ok) .do_finalize_upload(parsed_files)
      # Success or failure, the mapping dialog is gone: bring the wizard back
      # (on failure the error notification is shown and the user can re-upload).
      .reopen_wizard_step3()
    }
  }, ignoreInit = TRUE)

  # After gene column confirmed, proceed to architectural fork
  observeEvent(input$tmt_gene_col_confirm, {
    rv$pd_gene_col <- input$pd_gene_col_radio
    removeModal()
    Sys.sleep(0.3)
    .show_arch_fork_modal()
  }, ignoreInit = TRUE)

  # =========================================================================
  # SYNC: When other TMT modules upload files, sync here
  # =========================================================================
  observeEvent(rv$tmt_files_source, {
    if (!isTRUE(rv$tmt_files_source %in% c("trends", "tmt_pca", "tmt_heatmap"))) return()
    if (length(rv$tmt_files) == 0) return()

    tmt_upl_rv$files      <- rv$tmt_files
    tmt_upl_rv$n_files    <- rv$tmt_n_files
    tmt_upl_rv$conditions <- rv$tmt_conditions

    all_ch <- unique(unlist(lapply(rv$tmt_files, function(f) f$abund_cols)))
    tmt_upl_rv$all_channels <- all_ch
    tmt_upl_rv$status <- paste0(rv$tmt_n_files,
                                " file(s) synced. ", length(all_ch), " channels.")

    f1 <- rv$tmt_files[[1]]
    rv$pd_protein_headers  <- f1$headers
    rv$pd_raw_lines        <- f1$all_lines
    rv$pd_protein_line_idx <- f1$protein_idx
    rv$pd_ratio_cols       <- f1$ratio_cols
    rv$pd_pval_cols        <- f1$pval_cols
    rv$pd_other_cols       <- setdiff(f1$headers, c(f1$ratio_cols, f1$pval_cols))
    rv$pd_filename         <- f1$name

    pdata <- parse_pd_protein_data(f1, rv)
    has_gn <- TRUE
    if (!is.null(pdata) && "Description" %in% colnames(pdata$protein_df)) {
      descs <- pdata$protein_df[["Description"]]
      gn_hits <- sum(grepl("GN=", descs, ignore.case = FALSE), na.rm = TRUE)
      gn_pct  <- if (length(descs) > 0) gn_hits / length(descs) else 0
      if (!is.na(gn_pct) && gn_pct < 0.5) has_gn <- FALSE
    }
    rv$pd_has_gn <- has_gn

    logger(paste0("TMT Upload: Synced ", rv$tmt_n_files,
                  " file(s) from another TMT module"))
  }, ignoreInit = TRUE)

  # =========================================================================
  # SYNC: When other TMT modules set design, sync here
  # =========================================================================
  observeEvent(rv$tmt_design_source, {
    src <- rv$tmt_design_source
    if (is.null(src) || src == "volcano") return()
    ann <- rv$tmt_annotation_map
    if (is.null(ann) || nrow(ann) == 0) return()

    tmt_upl_rv$annotation_map    <- ann
    tmt_upl_rv$excluded_channels <- rv$tmt_excluded_channels %||% character(0)
    tmt_upl_rv$design_option     <- rv$tmt_design_option
    conds <- sort(unique(ann$Condition[ann$Condition != "Norm"]))
    tmt_upl_rv$conditions   <- conds
    tmt_upl_rv$design_conds <- conds
    tmt_upl_rv$has_bridge   <- any(ann$Condition == "Norm")
    tmt_upl_rv$n_plexes     <- length(unique(ann$Mixture))

    logger(paste0("TMT Upload: Design synced from ", src))
  }, ignoreInit = TRUE)

  # =========================================================================
  # FILE STATUS OUTPUT
  # =========================================================================
  output$tmt_upload_file_status <- renderUI({
    if (tmt_upl_rv$n_files == 0) {
      tags$p(class = "text-muted small", icon("circle-xmark"), " No files loaded")
    } else {
      file_lines <- lapply(tmt_upl_rv$files, function(f) {
        pep_badge <- if (f$has_peptides) {
          tags$span(class = "badge bg-success", "peptides")
        } else {
          tags$span(class = "badge bg-warning text-dark", "protein-only")
        }
        tags$li(
          tags$span(class = "small", f$name),
          tags$span(class = "text-muted small",
                    paste0(" (", length(f$protein_idx), " proteins, ",
                           paste(f$plex_prefix, collapse = "+"), ")")),
          " ", pep_badge
        )
      })
      tagList(
        tags$p(class = "small text-success mb-1",
               icon("circle-check"),
               paste0(" ", tmt_upl_rv$n_files, " file(s) loaded")),
        tags$ul(class = "list-unstyled mb-0 ps-2", file_lines)
      )
    }
  })

  # =========================================================================
  # SHARED DESIGN GRID (arch fork + grid + finish)
  # =========================================================================
  .dg <- tmt_register_design_grid(
    input = input, session = session,
    prefix = "tmt_", store = tmt_upl_rv, rv = rv, logger = logger,
    module_name = "TMT Upload", design_source_tag = "volcano",
    cond_placeholder = "e.g. Treated",
    on_design_finished = function(ann, NL, conds_final) {
      if (isTRUE(tmt_upl_rv$wizard_modal_open) || tmt_upl_rv$wizard_step == 4L) {
        design_opt <- tmt_upl_rv$design_option
        data_lvl   <- tmt_upl_rv$wizard_data_level
        tmt_upl_rv$analysis_route <- if (identical(data_lvl, "protein")) "limma_protein" else "limma_qf"
        tmt_upl_rv$wizard_step       <- 6L
        tmt_upl_rv$wizard_modal_open <- TRUE
        showModal(.tmt_wiz_build_modal())
        logger(paste0("TMT Wizard: Design complete. Route = ", tmt_upl_rv$analysis_route,
                      ". Reopening wizard at step 6."))
      }
    }
  )
  .show_arch_fork_modal <- .dg$show_arch_fork

  # =========================================================================
  # TMT WIZARD — modal builder and step UIs
  # =========================================================================

  .tmt_wiz_build_modal <- function() {
    modalDialog(
      title     = uiOutput("tmt_wiz_modal_title"),
      uiOutput("tmt_wiz_step_content"),
      footer    = uiOutput("tmt_wiz_footer"),
      size      = "l",
      easyClose = FALSE
    )
  }

  .tmt_wiz_step1_ui <- function() {
    tagList(
      tags$h5("Which software generated your data?"),
      radioButtons("tmt_wiz_search_program", NULL,
        choiceNames = list(
          tags$div(class = "p-2 border rounded mb-1",
            style = "background:#ebf5fb; border-color:#3498db !important;",
            tags$strong(style = "color:#3498db;",
              icon("database"), " Proteome Discoverer"),
            tags$br(),
            tags$small(class = "text-muted",
              "Nested CSV with protein + peptide data from Thermo PD.")
          ),
          tags$div(class = "p-2 border rounded mb-1",
            style = "background:#f9f9f9; border-color:#aaa !important;",
            tags$strong(icon("file-csv"), " Other Software (flat peptide table)"),
            tags$br(),
            tags$small(class = "text-muted",
              "Flat peptide table (e.g., FragPipe, Spectronaut, DIA-NN).")
          )
        ),
        choiceValues = c("pd", "other"),
        selected = isolate(tmt_upl_rv$wizard_search_program) %||% "pd"
      )
    )
  }

  .tmt_wiz_step2_ui <- function() {
    sp   <- isolate(tmt_upl_rv$wizard_search_program) %||% "pd"
    prev <- isolate(tmt_upl_rv$wizard_data_levels_raw)
    if (is.null(prev) || length(prev) == 0) {
      prev <- if (identical(sp, "pd")) c("protein", "peptide") else c("peptide")
    }
    # NA_character_ whenever the last pass through this step chose protein only.
    shared_prev <- isolate(tmt_upl_rv$shared_peptides)
    if (is.null(shared_prev) || is.na(shared_prev)) shared_prev <- "unique_only"
    tagList(
      tags$h5("What data types are in your export?"),
      tags$p(class = "small text-muted",
        "Select all that apply."),
      checkboxInput("tmt_wiz_data_protein",
        tagList(icon("table-cells"), tags$strong(" Protein-Level Data"),
          tags$span(class = "text-muted d-block small",
            "Pre-summarized protein abundances per TMT channel.")),
        value = "protein" %in% prev),
      checkboxInput("tmt_wiz_data_peptide",
        tagList(icon("tablets"), tags$strong(" Peptide-Level Data"),
          tags$span(class = "text-muted d-block small",
            "Peptide intensities (recommended for most workflows).")),
        value = "peptide" %in% prev),
      tags$div(class = "alert alert-info py-2 px-3 small mt-2",
        icon("circle-info"),
        " When both Protein and Peptide are selected, the recommended peptide-level pipeline will be used."),
      # Peptide route only. On the protein route there are no peptides and the
      # setting is meaningless, so the control is hidden and the stored value is
      # forced to NA_character_ (see the step 2 Next handler).
      conditionalPanel(
        condition = "input.tmt_wiz_data_peptide",
        tags$hr(class = "my-2"),
        tags$h6(class = "mt-2", "Peptides matching more than one protein"),
        radioButtons("tmt_wiz_shared_peptides", NULL,
          choiceNames = list(
            tags$div(class = "p-2 border rounded mb-1",
              tags$strong(icon("filter"), " Discard shared peptides"),
              tags$span(class = "badge bg-success ms-2", "Recommended"),
              tags$br(),
              tags$small(class = "text-muted",
                "Peptides whose Master Protein Accessions list more than one ",
                "protein are removed before protein-level summarization. This is ",
                "the MSstatsTMT default (useUniquePeptide = TRUE).")),
            tags$div(class = "p-2 border rounded mb-1",
              tags$strong(icon("triangle-exclamation"), " Keep all peptides"),
              tags$br(),
              tags$small(class = "text-muted",
                "Shared peptides are kept and attributed to the ",
                tags$strong("first accession listed"), " in the export. That is ",
                "position in a string, not evidence: it is not razor assignment, ",
                "and it can flatten the fold change of the protein that truly ",
                "contributed the signal.")),
            tags$div(class = "p-2 border rounded mb-1",
              tags$strong(icon("scale-balanced"), " Assign by parsimony (razor)"),
              tags$br(),
              tags$small(class = "text-muted",
                "Each shared peptide is assigned to the single candidate protein ",
                "with the most unique-peptide evidence across the dataset. ",
                "Peptides whose candidates have no unique support at all are ",
                "removed. This assignment is derived here from the exported ",
                "peptide table; it is not taken from the search engine."))
          ),
          choiceValues = c("unique_only", "keep_all", "razor"),
          selected = shared_prev
        )
      )
    )
  }

  .tmt_wiz_step3_ui <- function() {
    sp <- isolate(tmt_upl_rv$wizard_search_program)
    tagList(
      tags$h5("Upload your TMT file(s)"),
      tags$p(class = "small text-muted",
        "Upload one CSV/TSV per plex. Multi-select is supported."),
      if (identical(sp, "pd")) {
        tagList(
          tags$div(style = "display:none;",
            radioButtons("tmt_upload_input_format", NULL,
              choices = c("Proteome Discoverer (Nested)" = "pd_nested"),
              selected = "pd_nested")),
          tags$div(class = "alert alert-secondary py-2 px-3 small mb-2",
            icon("database"),
            " Files will be parsed as Proteome Discoverer nested CSV format.")
        )
      } else {
        tagList(
          tags$div(style = "display:none;",
            radioButtons("tmt_upload_input_format", NULL,
              choices = c("Flat Peptide Table" = "flat_peptide"),
              selected = "flat_peptide")),
          tags$p(class = "small text-info",
            icon("circle-info"),
            " Other software uses a flat peptide table.")
        )
      },
      fileInput("tmt_upload_csv_files", "Upload CSV(s)",
        accept = c(".csv", ".tsv"), width = "100%", multiple = TRUE),
      uiOutput("tmt_upload_file_status")
    )
  }

  .tmt_wiz_step4_ui <- function() {
    tagList(
      tags$h5("Single Plex or Multiple Plexes?"),
      radioButtons("tmt_wiz_plex_type", NULL,
        choiceNames = list(
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("vial"), " Single Plex"),
            tags$br(),
            tags$small(class = "text-muted",
              "All samples in one TMT experiment. ",
              "Each channel is one biological sample.")),
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("vials"), " Multiple Plexes"),
            tags$br(),
            tags$small(class = "text-muted",
              "Samples split across multiple TMT experiments (plexes). ",
              "Optional bridge/reference channel(s) for cross-plex normalization."))
        ),
        choiceValues = c("single", "multi"),
        selected = "single"
      ),
      # Plain-language guide for users who received a file from a core facility
      # and do not know what a "plex" is.
      tags$div(class = "alert alert-info py-2 px-3 small mt-2",
        tags$strong(icon("circle-question"), " Not sure which to pick?"),
        tags$p(class = "mb-1 mt-1",
          "The core facility measures TMT samples in ", tags$strong("sets"),
          ". Each set holds a limited number of samples, usually up to 16 or 18."),
        tags$ul(class = "mb-1 ps-3",
          tags$li(tags$strong("Single Plex:"),
            " all your samples were measured together in ", tags$strong("one set"),
            ". This is usually the case if you sent ",
            tags$strong("18 or fewer samples, all at the same time"), "."),
          tags$li(tags$strong("Multiple Plexes:"),
            " your samples were measured in ", tags$strong("two or more sets"),
            ". This happens when you had more samples than fit in one set, or ",
            "when you sent samples in separate batches at different times ",
            "(for example, a first experiment and a follow-up).")
        ),
        tags$p(class = "mb-0",
          "Still unsure? Ask your core facility: ",
          tags$em("\"Were my samples run as one TMT set or several?\""))
      )
    )
  }

  # ── Step 6 (single plex): were the abundances already normalized across channels? ──
  # Raw reporter abundances are normally not, so the box starts unticked and OSMSS
  # median-centers. tmt_normalization_check() runs on the uploaded protein rows as
  # soon as this step is shown, and its result is displayed next to the box as
  # advice. The user ticks the box only when the data were normalized before upload.
  .tmt_wiz_norm_question_ui <- function() {
    tags$div(class = "card mb-3",
      tags$div(class = "card-body py-2",
        tags$h6(class = "card-title", "Channel normalization"),
        tags$p(class = "small mb-2",
          "Raw TMT reporter abundances are normally not normalized, and OSMSS centers each ",
          "channel on the median of the channel medians (log2). If the abundances were already ",
          "normalized across channels before upload (by another tool or script), a second ",
          "centering can shift proteins that did not change, especially when a large part of the ",
          "proteome differs between conditions."),
        uiOutput("tmt_wiz_norm_check"),
        checkboxInput("tmt_wiz_input_normalized",
          "The abundances were already normalized across channels (skip OSMSS median centering)",
          value = isTRUE(tmt_upl_rv$input_normalized)),
        tags$p(class = "text-muted small mb-0",
          "Leave unticked unless you know the values were normalized before upload. ",
          "The data check compares channel medians between replicate channels of each condition; ",
          "it cannot tell loading differences from real differences between conditions.")
      )
    )
  }

  tmt_norm_check <- reactive({
    req(tmt_upl_rv$wizard_step == 6L, tmt_upl_rv$n_files > 0,
        identical(as.integer(tmt_upl_rv$design_option), 1L), !is.null(tmt_upl_rv$annotation_map))
    tryCatch(
      tmt_normalization_check(tmt_upl_rv$files, tmt_upl_rv$annotation_map,
                              rv = list(pd_has_gn = isTRUE(rv$pd_has_gn), pd_gene_col = rv$pd_gene_col)),
      error = function(e) list(verdict = "inconclusive", within_sd = NA_real_, n_complete = 0L,
                               message = paste0("The check could not run: ", conditionMessage(e))))
  })

  output$tmt_wiz_norm_check <- renderUI({
    chk <- tmt_norm_check()
    cls <- switch(chk$verdict, normalized = "alert-success", borderline = "alert-warning",
                  not_normalized = "alert-secondary", "alert-secondary")
    ic  <- switch(chk$verdict, normalized = "circle-info", borderline = "triangle-exclamation",
                  "circle-check")
    tags$div(class = paste("alert py-2 px-3 small mb-2", cls),
      icon(ic), tags$strong(" Data check: "), chk$message,
      switch(chk$verdict,
        normalized     = tags$span(" If you know the values were normalized before upload, tick the box below; ",
                                   "otherwise OSMSS centers them as usual."),
        borderline     = tags$span(" Tick the box only if you know the values were normalized before upload."),
        not_normalized = tags$span(" OSMSS median centering is recommended (box unticked)."),
        NULL))
  })

  observeEvent(input$tmt_wiz_input_normalized, {
    tmt_upl_rv$input_normalized <- isTRUE(input$tmt_wiz_input_normalized)
  }, ignoreInit = TRUE)

  output$tmt_wiz_norm_step <- renderUI({
    skip <- isTRUE(input$tmt_wiz_input_normalized %||% tmt_upl_rv$input_normalized)
    if (identical(tmt_upl_rv$analysis_route, "limma_qf"))
      if (skip) "No median centering in OSMSS (input reported as already normalized across channels)"
      else "Peptide-level column-wise median centering (before aggregation)"
    else
      if (skip) "No median centering in OSMSS (input reported as already normalized across channels)"
      else "Column-wise median centering (each channel shifted to the median of the channel medians)"
  })

  .tmt_wiz_step6_summary_ui <- function() {
    data_level  <- tmt_upl_rv$wizard_data_level %||% "?"
    search_prog <- tmt_upl_rv$wizard_search_program %||% "pd"
    design_opt  <- tmt_upl_rv$design_option
    plex_type <- if (identical(as.integer(design_opt), 1L)) "Single Plex"
                 else if (identical(as.integer(design_opt), 2L)) "Multi-Plex + Bridge"
                 else "Multi-Plex (Balanced)"
    route <- tmt_upl_rv$analysis_route %||% "?"

    # NULL on the protein route: the row is omitted rather than reporting a
    # filter that could not have run.
    shared_pep <- tmt_upl_rv$shared_peptides
    shared_pep_label <- if (is.null(shared_pep) || is.na(shared_pep)) NULL else switch(
      shared_pep,
      unique_only = "Discarded (peptides matching >1 protein removed)",
      keep_all    = "Kept, attributed to first listed accession",
      razor       = "Assigned by parsimony (most unique-peptide evidence)",
      as.character(shared_pep))

    route_label <- switch(route,
      limma_protein = "Limma (Protein-Level)",
      limma_qf      = "QFeatures + Limma",
      "Unknown"
    )
    # Steps as the code runs them: summarization in tmt_limma_protein_summarize()
    # / tmt_limma_qf_summarize() + qf_summarize_to_protein_matrix(); filter and
    # imputation in tmt_prepare_matrix_for_de(); the model in compute_tmt_de()
    # -> dia_limma_de(). Keep in step with those functions when they change.
    dopt <- as.integer(design_opt)
    .filter_step <- tags$li(paste0(
      "Global missingness filter: a protein is kept if at least 50% of its channels ",
      "are observed in at least one condition (the same protein set for every contrast)"))
    .impute_step <- tags$li(
      if (identical(dopt, 1L))
        paste0("k-NN imputation (k = 10) over all channels: values are filled only for ",
               "proteins with at most 1 missing channel in every condition; other ",
               "missing values stay missing and limma fits the observed values",
               if (identical(route, "limma_qf"))
                 "; cells whose reporter ions were present but non-positive are not filled"
               else "")
      else if (identical(dopt, 2L))
        "No imputation (missing values stay missing; limma fits the observed values)"
      else
        paste0("k-NN imputation (k = 10) within each plex, for proteins with at least 50% ",
               "of that plex's channels observed; other missing values stay missing"))
    .limma_step <- tags$li(paste0(
      "Limma DE: one model over all conditions",
      if (identical(dopt, 3L)) " plus a plex covariate" else "",
      "; condition-level sample weights (arrayWeights by condition, estimated on ",
      "proteins with no missing values, when at least 50 exist) \u2192 lmFit \u2192 ",
      "contrasts.fit \u2192 eBayes(trend = TRUE, robust = TRUE) \u2192 BH adjustment within ",
      "the selected contrast"))
    route_steps <- switch(route,
      limma_protein = tags$ol(class = "small mb-0",
        tags$li("Zero and negative abundances \u2192 missing; log2 transformation"),
        if (!identical(dopt, 1L))
          tags$li("Plex coverage filter: a protein must be quantified in at least half of the plexes"),
        if (identical(dopt, 1L)) uiOutput("tmt_wiz_norm_step", container = tags$li)
        else tags$li(if (identical(dopt, 2L))
          "Column-wise median centering, then IRS normalization (bridge-channel alignment across plexes)"
        else
          "Column-wise median centering (each channel shifted to the median of the channel medians)"),
        .filter_step, .impute_step, .limma_step
      ),
      limma_qf = tags$ol(class = "small mb-0",
        tags$li("Isolation interference filter (\u226450%; rows with no value are kept), if the column is present"),
        tags$li("Shared peptides handled as set above"),
        tags$li("Zero and negative intensities \u2192 missing; log2 transformation"),
        if (identical(dopt, 1L)) uiOutput("tmt_wiz_norm_step", container = tags$li)
        else tags$li("Peptide-level column-wise median centering (before aggregation)"),
        tags$li(paste0("QFeatures aggregation to proteins via robustSummary (single-peptide proteins ",
                       "taken as is; column means for a protein where robustSummary fails)")),
        if (identical(dopt, 2L))
          tags$li("IRS normalization (bridge-channel alignment across plexes)"),
        .filter_step, .impute_step, .limma_step
      ),
      tags$p("Unknown route")
    )

    conds <- tmt_upl_rv$conditions %||% character(0)
    ann   <- tmt_upl_rv$annotation_map
    excl  <- tmt_upl_rv$excluded_channels %||% character(0)

    tagList(
      tags$h5("Review & Confirm"),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body py-2",
          tags$h6(class = "card-title", "Configuration Summary"),
          tags$table(class = "table table-sm mb-0",
            tags$tbody(
              tags$tr(tags$td(class = "fw-bold", "Source:"),
                      tags$td(if (identical(search_prog, "pd")) "Proteome Discoverer"
                              else "Other (flat peptide table)")),
              tags$tr(tags$td(class = "fw-bold", "Data Level:"),
                      tags$td(data_level)),
              tags$tr(tags$td(class = "fw-bold", "Design:"),
                      tags$td(plex_type)),
              # Omitted on the protein route, where no peptide filter can run.
              if (!is.null(shared_pep_label))
                tags$tr(tags$td(class = "fw-bold", "Shared peptides:"),
                        tags$td(shared_pep_label)),
              tags$tr(tags$td(class = "fw-bold", "Analysis Route:"),
                      tags$td(route_label)),
              tags$tr(tags$td(class = "fw-bold", "Files:"),
                      tags$td(paste0(tmt_upl_rv$n_files, " file(s)"))),
              tags$tr(tags$td(class = "fw-bold", "Conditions:"),
                      tags$td(paste(conds, collapse = ", "))),
              tags$tr(tags$td(class = "fw-bold", "Channels:"),
                      tags$td(paste0(if (!is.null(ann)) nrow(ann) else 0, " analyzed"))),
              if (length(excl) > 0)
                tags$tr(tags$td(class = "fw-bold", "Excluded:"),
                        tags$td(paste0(length(excl), " channel(s): ",
                                       paste(tmt_short_label(excl), collapse = "; "))))
            )
          )
        )
      ),
      if (identical(dopt, 1L)) .tmt_wiz_norm_question_ui(),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body py-2",
          tags$h6(class = "card-title", "Pipeline Steps"),
          route_steps
        )
      ),
      tags$p(class = "small text-muted",
        "Click Summarize to process your data. Once complete, the ",
        tags$strong("Done"), " button will turn green."),
      actionButton("tmt_wiz_summarize_btn",
        tagList(icon("gears"), " Summarize Data"),
        class = "btn-primary w-100 mb-2"),
      uiOutput("tmt_wiz_summarize_status")
    )
  }

  # ── Wizard title bar ──
  output$tmt_wiz_modal_title <- renderUI({
    step      <- tmt_upl_rv$wizard_step
    reachable <- tmt_upl_rv$wizard_reachable_steps
    pos       <- match(step, reachable)
    if (is.na(pos)) pos <- 1L
    n <- length(reachable)
    tagList(
      tags$span(class = "fw-bold", "TMT Setup Wizard"),
      tags$span(class = "text-muted small ms-2",
        paste0(" \u2014 Step ", pos, " of ", n)),
      tags$div(class = "progress mt-1", style = "height: 4px;",
        tags$div(class = "progress-bar bg-success",
          style = paste0("width:", round(100 * pos / n), "%;")))
    )
  })

  # ── Wizard step content ──
  output$tmt_wiz_step_content <- renderUI({
    step <- tmt_upl_rv$wizard_step
    switch(as.character(step),
      "1" = .tmt_wiz_step1_ui(),
      "2" = .tmt_wiz_step2_ui(),
      "3" = .tmt_wiz_step3_ui(),
      "4" = .tmt_wiz_step4_ui(),
      "6" = .tmt_wiz_step6_summary_ui(),
      tags$p(class = "text-muted", paste0("Unknown step: ", step))
    )
  })

  # ── Wizard footer ──
  output$tmt_wiz_footer <- renderUI({
    step      <- tmt_upl_rv$wizard_step
    reachable <- tmt_upl_rv$wizard_reachable_steps
    pos       <- match(step, reachable)
    if (is.na(pos)) pos <- 1L
    n         <- length(reachable)
    is_last   <- (pos == n)

    cancel_btn <- modalButton("Cancel")
    back_btn   <- if (pos > 1 && step != 6L) {
      actionButton("tmt_wiz_back", "Back",
        class = "btn-secondary", icon = icon("arrow-left"))
    } else if (step == 6L) {
      actionButton("tmt_wiz_back", "Back",
        class = "btn-secondary", icon = icon("arrow-left"))
    } else NULL

    if (is_last || step == 6L) {
      has_cache <- !is.null(rv$tmt_master_cache) || !is.null(tmt_upl_rv$limma_cache)
      done_btn  <- actionButton("tmt_wiz_close", "Done",
        class = if (has_cache) "btn-success" else "btn-secondary",
        icon  = icon("check"))
      if (!has_cache) done_btn <- shinyjs::disabled(done_btn)
      tagList(cancel_btn, back_btn, done_btn)
    } else {
      next_btn <- actionButton("tmt_wiz_next", "Next",
        class = "btn-primary", icon = icon("arrow-right"))
      tagList(cancel_btn, back_btn, next_btn)
    }
  })

  # ── Wizard summarize status badge ──
  output$tmt_wiz_summarize_status <- renderUI({
    has_cache <- !is.null(tmt_upl_rv$limma_cache)
    if (has_cache) {
      tags$div(class = "alert alert-success py-1 px-3 small mb-0 mt-1",
        icon("circle-check"), " Summarization complete.")
    } else {
      NULL
    }
  })

  # ── Wizard completed flag ──
  output$tmt_wiz_completed_flag <- reactive({
    isTRUE(tmt_upl_rv$wizard_completed)
  })
  outputOptions(output, "tmt_wiz_completed_flag", suspendWhenHidden = FALSE)

  # ── Sidebar badges (post-wizard) ──
  output$tmt_wiz_sidebar_badges <- renderUI({
    route      <- tmt_upl_rv$analysis_route %||% "?"
    data_lvl   <- tmt_upl_rv$wizard_data_level %||% "?"
    search_prg <- tmt_upl_rv$wizard_search_program %||% "?"
    design_opt <- tmt_upl_rv$design_option
    conds      <- tmt_upl_rv$conditions %||% character(0)
    n_ch       <- length(tmt_upl_rv$all_channels)
    n_excl     <- length(tmt_upl_rv$excluded_channels)

    plex_str <- if (identical(as.integer(design_opt), 1L)) "Single Plex"
                else if (identical(as.integer(design_opt), 2L)) "Multi-Plex + Bridge"
                else "Multi-Plex (Balanced)"
    route_str <- switch(route,
      limma_protein = "Limma (Protein-Level)",
      limma_qf      = "QFeatures + Limma",
      route
    )

    tags$div(class = "card mb-2",
      tags$div(class = "card-body py-2 px-2",
        tags$h6(class = "card-title mb-1", icon("check-circle"), " Setup Complete"),
        tags$table(class = "table table-sm table-borderless mb-0",
          style = "font-size:0.82em;",
          tags$tbody(
            tags$tr(tags$td(class = "text-muted", "Source:"),
                    tags$td(if (identical(search_prg, "pd")) "Proteome Discoverer"
                            else "Other (flat peptide table)")),
            tags$tr(tags$td(class = "text-muted", "Level:"), tags$td(data_lvl)),
            tags$tr(tags$td(class = "text-muted", "Design:"), tags$td(plex_str)),
            tags$tr(tags$td(class = "text-muted", "Route:"), tags$td(route_str)),
            tags$tr(tags$td(class = "text-muted", "Files:"),
                    tags$td(paste0(tmt_upl_rv$n_files, " file(s)"))),
            tags$tr(tags$td(class = "text-muted", "Conditions:"),
                    tags$td(paste(conds, collapse = ", "))),
            tags$tr(tags$td(class = "text-muted", "Channels:"),
                    tags$td(paste0(n_ch, " total",
                                   if (n_excl > 0) paste0(", ", n_excl, " excluded") else "")))
          )
        )
      )
    )
  })

  # =========================================================================
  # TMT WIZARD — navigation handlers
  # =========================================================================

  observeEvent(input$tmt_wiz_open, {
    tmt_upl_rv$wizard_step       <- 1L
    tmt_upl_rv$wizard_modal_open <- TRUE
    showModal(.tmt_wiz_build_modal())
  })

  observeEvent(input$tmt_wiz_reopen, {
    tmt_upl_rv$wizard_step       <- 1L
    tmt_upl_rv$wizard_modal_open <- TRUE
    showModal(.tmt_wiz_build_modal())
  })

  # Next button
  observeEvent(input$tmt_wiz_next, {
    step <- tmt_upl_rv$wizard_step

    if (step == 1L) {
      if (!is.null(input$tmt_wiz_search_program))
        tmt_upl_rv$wizard_search_program <- input$tmt_wiz_search_program
      tmt_upl_rv$wizard_step <- 2L

    } else if (step == 2L) {
      pep_sel  <- isTRUE(input$tmt_wiz_data_peptide)
      prot_sel <- isTRUE(input$tmt_wiz_data_protein)
      if (!pep_sel && !prot_sel) {
        showNotification("Please select at least one data type.", type = "warning")
        return()
      }
      tmt_upl_rv$wizard_data_levels_raw <- c(
        if (prot_sel) "protein",
        if (pep_sel)  "peptide"
      )
      tmt_upl_rv$wizard_data_level <- if (pep_sel) "peptide" else "protein"

      # Shared-peptide handling. Meaningless on the protein route, so the stored
      # value is forced to NA_character_ there and nothing downstream reports a
      # filter that could not have run.
      shared_new <- if (pep_sel)
        (input$tmt_wiz_shared_peptides %||% "unique_only") else NA_character_
      shared_old <- tmt_upl_rv$shared_peptides
      tmt_upl_rv$shared_peptides <- shared_new

      # Invalidate the TMT master cache, exactly as a design change does: this
      # setting changes the protein matrix, so a stale cache would silently mix
      # two settings.
      if (!identical(shared_old, shared_new) && !is.null(rv$tmt_master_cache)) {
        rv$tmt_master_cache        <- NULL
        rv$tmt_master_cache_design <- NULL
        rv$tmt_master_cache_ann    <- NULL
        showModal(modalDialog(
          title = "Preprocessing Parameter Changed",
          "The shared-peptide setting was modified. The cached TMT summarization has been ",
          "invalidated. Please re-run your analysis.",
          easyClose = FALSE,
          footer = modalButton("OK")
        ))
        logger(paste0("TMT Upload: Shared-peptide handling changed  -  ",
                      "TMT master cache invalidated"))
      }

      tmt_upl_rv$wizard_step <- 3L

    } else if (step == 3L) {
      if (tmt_upl_rv$n_files == 0) {
        showNotification("Please upload at least one TMT CSV file.", type = "warning")
        return()
      }
      tmt_upl_rv$wizard_step <- 4L

    } else if (step == 4L) {
      plex_type <- input$tmt_wiz_plex_type %||% "single"
      removeModal()

      if (identical(plex_type, "single")) {
        tmt_upl_rv$n_plexes          <- 1L
        tmt_upl_rv$has_bridge        <- FALSE
        tmt_upl_rv$design_option     <- 1L
        tmt_upl_rv$wizard_modal_open <- TRUE
        if (length(tmt_upl_rv$design_conds) == 0) tmt_upl_rv$design_conds <- character(0)
        .dg$show_design_grid()
      } else {
        n_f <- max(2L, tmt_default_channel_plex(tmt_upl_rv$files, tmt_upl_rv$all_channels))
        tmt_upl_rv$n_plexes         <- as.integer(n_f)
        tmt_upl_rv$has_bridge       <- FALSE
        tmt_upl_rv$design_conds     <- character(0)
        tmt_upl_rv$grid_assignments <- NULL
        tmt_upl_rv$paint_cond       <- NULL
        tmt_upl_rv$wizard_modal_open <- TRUE
        if (length(tmt_upl_rv$all_channels) == 0) {
          showNotification("No channels detected. Please upload TMT files first.",
                           type = "error")
          return()
        }
        Sys.sleep(0.2)
        .dg$show_design_grid()
      }
    }
  })

  # Back button
  observeEvent(input$tmt_wiz_back, {
    step <- tmt_upl_rv$wizard_step
    if (step == 2L)      tmt_upl_rv$wizard_step <- 1L
    else if (step == 3L) tmt_upl_rv$wizard_step <- 2L
    else if (step == 4L) tmt_upl_rv$wizard_step <- 3L
    else if (step == 6L) tmt_upl_rv$wizard_step <- 4L
  })

  # Done — mark wizard complete, set rv$data_source, clear other caches
  observeEvent(input$tmt_wiz_close, {
    removeModal()
    tmt_upl_rv$wizard_modal_open <- FALSE
    tmt_upl_rv$wizard_completed  <- TRUE

    # ── Sync to shared reactive values ──
    rv$tmt_wizard_completed      <- TRUE
    rv$tmt_analysis_route        <- tmt_upl_rv$analysis_route
    rv$tmt_wizard_search_program <- tmt_upl_rv$wizard_search_program
    rv$tmt_wizard_data_level     <- tmt_upl_rv$wizard_data_level
    rv$tmt_annotation_map        <- tmt_upl_rv$annotation_map
    rv$tmt_conditions            <- tmt_upl_rv$conditions
    rv$tmt_design_option         <- tmt_upl_rv$design_option

    # ── Set active data source ──
    # The previous experiment (LF, DIA, older TMT) was already cleared when
    # these files were uploaded (.do_finalize_upload -> reset_experiment_state).
    rv$data_source <- "tmt"

    # ── Enable DE run button if cache is already available ──
    has_cache <- !is.null(rv$tmt_master_cache) || !is.null(tmt_upl_rv$limma_cache)
    if (has_cache) .enable_run_btn()

    logger(paste0("TMT Wizard: Completed. Route = ", tmt_upl_rv$analysis_route,
                  ", Design = ", tmt_upl_rv$design_option,
                  ". rv$data_source set to 'tmt'."))
  })

  # =========================================================================
  # WIZARD SUMMARIZE BUTTON (step 6)
  # =========================================================================
  observeEvent(input$tmt_wiz_summarize_btn, {
    req(tmt_upl_rv$n_files > 0)

    if (is.null(tmt_upl_rv$design_option) || is.null(tmt_upl_rv$annotation_map)) {
      showNotification("Please define an experimental design first.", type = "warning")
      return()
    }

    if (identical(tmt_wiz_task$status(), "running")) {
      showNotification("TMT summarization already running.", type = "warning")
      return()
    }

    route <- tmt_upl_rv$analysis_route %||% "limma_qf"

    rv_snapshot <- list(
      pd_has_gn  = isTRUE(rv$pd_has_gn),
      pd_gene_col = rv$pd_gene_col
    )

    # Single plex only: skip median centering when the user reported the input as
    # already normalized across channels. Multi-plex designs always center.
    .single <- identical(as.integer(tmt_upl_rv$design_option), 1L)
    apply_centering <- !(.single && isTRUE(input$tmt_wiz_input_normalized))
    .chk <- if (.single) tryCatch(tmt_norm_check(), error = function(e) NULL) else NULL
    tmt_upl_rv$pending_norm <- list(
      input_normalized_upstream = if (.single) isTRUE(input$tmt_wiz_input_normalized) else NA,
      norm_check_verdict        = .chk$verdict    %||% NA_character_,
      norm_check_within_sd      = .chk$within_sd  %||% NA_real_,
      norm_check_n_complete     = .chk$n_complete %||% NA_integer_)

    logger(paste0("TMT Wizard: Starting async summarization (route=", route,
                  if (!apply_centering) ", median centering skipped: input reported as normalized" else "", ")"))
    show_async_overlay(session,
                       label  = paste0("TMT summarization (", route, ")"),
                       detail = "Please wait...")
    tmt_wiz_cancelled(FALSE)
    tmt_wiz_task$invoke(
      route             = route,
      files             = tmt_upl_rv$files,
      annotation_map    = tmt_upl_rv$annotation_map,
      design_opt        = tmt_upl_rv$design_option,
      wizard_data_level = tmt_upl_rv$wizard_data_level,
      shared_peptides   = tmt_upl_rv$shared_peptides %||% "unique_only",
      rv_snapshot       = rv_snapshot,
      apply_centering   = apply_centering
    )
  })

  # =========================================================================
  # SYNC Run button state when shared cache changes
  # =========================================================================
  observeEvent(rv$tmt_master_cache, {
    if (!is.null(rv$tmt_master_cache)) {
      .enable_run_btn()
    } else if (is.null(tmt_upl_rv$limma_cache)) {
      .disable_run_btn()
    }
  }, ignoreInit = TRUE)

  # =========================================================================
  # EXPORT SUMMARIZATION (.rds)
  # =========================================================================
  output$tmt_upload_export_summarization <- downloadHandler(
    filename = function() {
      paste0("TMT_summarization_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      has_limma <- !is.null(tmt_upl_rv$limma_cache)
      if (!has_limma) {
        showNotification("No summarization to export. Run Summarize Data first.",
                         type = "warning")
        return()
      }
      lc   <- tmt_upl_rv$limma_cache
      ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
      base <- paste0("TMT_summarization_", ts)
      bundle <- list(
        type           = "tmt_summarization",
        route          = tmt_upl_rv$analysis_route,
        limma_cache    = lc,
        design_option  = lc$design_option %||% tmt_upl_rv$design_option,
        annotation_map = lc$annotation_map %||% tmt_upl_rv$annotation_map,
        conditions     = tmt_upl_rv$conditions,
        data_level     = tmt_upl_rv$wizard_data_level,
        search_program = tmt_upl_rv$wizard_search_program,
        # Recorded so an imported bundle reports the setting it was built with,
        # rather than whatever the current default happens to be.
        shared_peptides = lc$shared_peptide_handling %||% tmt_upl_rv$shared_peptides,
        exported_at    = Sys.time()
      )
      bundle$pca_thumbnail <- make_pca_thumbnail(bundle)
      bundle$n_features    <- count_bundle_features(bundle)
      rds_path     <- file.path(tempdir(), paste0(base, ".rds"))
      preview_path <- file.path(tempdir(), paste0(base, ".preview.rds"))
      saveRDS(bundle, rds_path)
      preview <- save_prerun_preview(bundle, paste0(base, ".rds"), logger)
      saveRDS(preview, preview_path)
      utils::zip(zipfile = file, files = c(rds_path, preview_path), flags = "-j")
      unlink(c(rds_path, preview_path))
      logger("TMT Upload: Summarization exported to ZIP")
    }
  )

# =========================================================================
  # IMPORT SUMMARIZATION (.rds) (From Prerun/Saved Session)
  # =========================================================================
  observeEvent(rv$prerun_load_tmt_bundle, {
    req(rv$prerun_load_tmt_bundle)
    bundle <- rv$prerun_load_tmt_bundle
    
    logger("TMT Upload: Importing prerun TMT bundle...")
    
    # 1. Restore local module state
    tmt_upl_rv$analysis_route        <- bundle$route
    tmt_upl_rv$design_option         <- bundle$design_option
    tmt_upl_rv$annotation_map        <- bundle$annotation_map
    tmt_upl_rv$excluded_channels     <- bundle$limma_cache$excluded_channels %||% character(0)
    tmt_upl_rv$conditions            <- bundle$conditions
    tmt_upl_rv$wizard_data_level     <- bundle$data_level
    tmt_upl_rv$wizard_search_program <- bundle$search_program
    # Bundles predating the shared-peptide setting carry no key. Those runs kept
    # every peptide, so they restore as "keep_all" - defaulting them to the new
    # default would misreport what they actually did.
    tmt_upl_rv$shared_peptides       <- bundle$shared_peptides %||% "keep_all"
    tmt_upl_rv$limma_cache           <- bundle$limma_cache
    tmt_upl_rv$wizard_completed      <- TRUE
    
    # 2. Restore global shared state (rv)
    rv$tmt_wizard_completed      <- TRUE
    rv$tmt_analysis_route        <- bundle$route
    rv$tmt_wizard_search_program <- bundle$search_program
    rv$tmt_wizard_data_level     <- bundle$data_level
    rv$tmt_annotation_map        <- bundle$annotation_map
    rv$tmt_conditions            <- bundle$conditions
    rv$tmt_design_option         <- bundle$design_option
    rv$tmt_limma_cache           <- bundle$limma_cache
    
    # 3. Reconstruct the legacy tmt_master_cache for Trend Clusters
    if (!is.null(bundle$limma_cache)) {
      res        <- bundle$limma_cache
      ann_mc     <- res$annotation_map
      pm_mc      <- res$protein_matrix
      bio_ann_mc <- ann_mc[!is.na(ann_mc$Condition) & 
                             !ann_mc$Condition %in% "Norm", , drop = FALSE]
      shared_mc  <- intersect(bio_ann_mc$ColumnName, colnames(pm_mc))
      
      if (length(shared_mc) > 0) {
        pld_mc <- do.call(rbind, lapply(shared_mc, function(ch) {
          row_mc <- bio_ann_mc[bio_ann_mc$ColumnName == ch, , drop = FALSE]
          data.frame(
            Protein          = rownames(pm_mc),
            Condition        = row_mc$Condition,
            BioReplicate     = row_mc$BioReplicate,
            Abundance        = pm_mc[, ch],
            Run              = row_mc$Mixture,   
            stringsAsFactors = FALSE
          )
        }))
        rv$tmt_master_cache        <- list(ProteinLevelData = pld_mc)
        rv$tmt_master_cache_design <- res$design_option
        rv$tmt_master_cache_ann    <- res$annotation_map
      }
    }
    
    # 4. Update UI helpers
    .enable_run_btn()
    
    # 5. Set the active data source so downstream modules recognize TMT is active
    rv$data_source <- "tmt"
    
    # 6. Clean up the trigger to prevent re-firing
    rv$prerun_load_tmt_bundle <- NULL
    
    logger(paste0("TMT Upload: Prerun import complete. Route = ", bundle$route))
  }, ignoreInit = TRUE)
}
