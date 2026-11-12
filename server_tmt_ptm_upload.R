# server_tmt_ptm_upload.R
# ==============================================================================
# TMT PTM Upload Wizard Module
#
# Handles TMT phosphoproteomics data upload (Proteome Discoverer nested CSV),
# experimental design, and limma-based summarization using
# tmt_ptm_limma_summarize() defined in global.R.
#
# On wizard completion:
#   rv$dia_ptm_master_cache      <- tmt_ptm_limma_summarize() output
#   rv$data_source               <- "tmt_ptm"
#   rv$tmt_ptm_wizard_completed  <- TRUE
#   rv$ptm_data$source           <- "tmt_ptm_upload"  (triggers downstream sync)
#   rv$ptm_data$design_option    <- "tmt_ptm"
#
# Stores result in rv$dia_ptm_master_cache so all existing downstream PTM
# modules (PTM DE, PTM Trends, PTM Single, KS Scoring) work without change.
# ==============================================================================

server_tmt_ptm_upload <- function(input, output, session, rv, logger) {

  N_STEPS <- 7L
  STEP_LABELS <- c(
    "Input Type",
    "Upload Files",
    "PTM Types",
    "Plex Design",
    "Probability Threshold",
    "Summary & Summarize",
    "Review"
  )

  # ── Module-local reactive values ──────────────────────────────────────────
  ptm_upl_rv <- reactiveValues(
    # Wizard navigation
    wizard_step        = 1L,
    wizard_modal_open  = FALSE,
    wizard_completed   = FALSE,
    # File state (per-plex lists)
    n_plexes           = 1L,
    # Raw upload info — stored immediately on file select, cleared after parse
    phospho_raw        = list(),    # per-plex: list(datapath, name)
    global_raw         = list(),    # per-plex: list(datapath, name)
    # Parsed structs — populated after async parse completes
    phospho_files      = list(),    # parsed file structs
    global_files       = list(),    # parsed file structs
    phospho_statuses   = list(),    # per-plex: list(uploaded/n_peptides/n_channels/error)
    global_statuses    = list(),    # per-plex: list(uploaded/n_proteins/n_channels/error)
    parsing_in_progress = FALSE,   # TRUE while parse task is running
    # PTM types
    detected_ptms      = character(0),
    selected_ptms      = "Phospho",
    # Design grid state (required fields for tmt_register_design_grid)
    all_channels       = character(0),
    has_bridge         = FALSE,
    design_conds       = character(0),
    design_option      = NULL,
    annotation_map     = NULL,
    excluded_channels  = character(0),   # "Exclude" in the grid; not in annotation_map
    conditions         = character(0),
    grid_assignments   = NULL,
    paint_cond         = NULL,
    # Probability threshold
    prob_threshold     = 75,
    min_valid_per_cond = 2L,
    interference_threshold = 50,
    # Summarization results
    summarization_done = FALSE,
    n_sites            = 0L,
    n_sites_adjusted   = 0L
  )

  # ── Completed flag output (for conditionalPanel in ui.R) ──────────────────
  output$tmt_ptm_wiz_completed_flag <- reactive({
    isTRUE(rv$tmt_ptm_wizard_completed)
  })
  outputOptions(output, "tmt_ptm_wiz_completed_flag", suspendWhenHidden = FALSE)

  # ── Register design grid (arch fork + channel assignment grid) ────────────
  .dg <- tmt_register_design_grid(
    input             = input,
    session           = session,
    prefix            = "tmt_ptm_",
    store             = ptm_upl_rv,
    rv                = rv,
    logger            = logger,
    module_name       = "TMT PTM Upload",
    design_source_tag = "tmt_ptm_upload",
    cond_placeholder  = "e.g. Treated",
    on_design_finished = function(ann, NL, conds_final) {
      ptm_upl_rv$annotation_map    <- ann
      ptm_upl_rv$conditions        <- conds_final
      # Reopen wizard at Step 5 (Probability Threshold) after design grid closes
      ptm_upl_rv$wizard_step       <- 5L
      ptm_upl_rv$wizard_modal_open <- TRUE
      showModal(.build_modal())
      logger(paste0("[TMT PTM Wizard] Design complete — ",
                    length(conds_final), " condition(s). Reopening at Step 5."))
    }
  )

  # ── Modal builder ──────────────────────────────────────────────────────────
  .build_modal <- function() {
    modalDialog(
      title     = uiOutput("tmt_ptm_wiz_title"),
      uiOutput("tmt_ptm_wiz_content"),
      footer    = uiOutput("tmt_ptm_wiz_footer"),
      size      = "l",
      easyClose = FALSE
    )
  }

  # ============================================================================
  # ASYNC FILE PARSE TASK
  # Triggered when user clicks Next on step 2. Parses all uploaded plexes in a
  # background worker so the UI stays responsive. Returns parsed structs +
  # summary counts for every plex.
  # ============================================================================
  tmt_ptm_parse_task <- ExtendedTask$new(function(phospho_raws, global_raws) {
    promises::future_promise({
      # Parse phospho files
      phospho_out <- lapply(phospho_raws, function(f) {
        tryCatch({
          lines  <- readLines(f$datapath, warn = FALSE)
          lines  <- collapse_multiline_csv(lines)
          struct <- parse_tmt_file_structure(lines, f$name)
          pep_df <- parse_pd_peptide_data(struct)
          struct$.peptide_df_cache <- pep_df   # cache so detect_ptm_types is fast
          list(
            struct     = struct,
            n_peptides = if (!is.null(pep_df)) nrow(pep_df) else 0L,
            n_channels = length(struct$abund_cols %||% character(0))
          )
        }, error = function(e) list(error = conditionMessage(e), name = f$name))
      })
      # Parse global files
      global_out <- lapply(global_raws, function(f) {
        tryCatch({
          lines  <- readLines(f$datapath, warn = FALSE)
          lines  <- collapse_multiline_csv(lines)
          struct <- parse_tmt_file_structure(lines, f$name)
          list(
            struct     = struct,
            n_proteins = length(struct$protein_idx),
            n_channels = length(struct$abund_cols %||% character(0))
          )
        }, error = function(e) list(error = conditionMessage(e), name = f$name))
      })
      list(success = TRUE, phospho = phospho_out, global = global_out)
    }, seed = TRUE)
  })

  # ── Parse task result handler ──────────────────────────────────────────────
  observeEvent(tmt_ptm_parse_task$result(), {
    hide_async_overlay(session)
    ptm_upl_rv$parsing_in_progress <- FALSE

    res <- tmt_ptm_parse_task$result()
    if (!isTRUE(res$success)) {
      showNotification(
        paste0("File parsing failed: ", res$error %||% "unknown error"),
        type = "error", duration = 15)
      ptm_upl_rv$wizard_step <- 2L
      showModal(.build_modal())
      return()
    }

    # Apply results per plex
    parse_errors <- character(0)
    for (i in seq_along(res$phospho)) {
      r <- res$phospho[[i]]
      if (!is.null(r$error)) {
        parse_errors <- c(parse_errors, paste0("Plex ", i, " phospho: ", r$error))
        ptm_upl_rv$phospho_statuses[[i]] <- list(error = r$error)
      } else {
        ptm_upl_rv$phospho_files[[i]]    <- r$struct
        ptm_upl_rv$phospho_raw[[i]]      <- NULL   # clear raw after successful parse
        ptm_upl_rv$phospho_statuses[[i]] <- list(
          n_peptides = r$n_peptides,
          n_channels = r$n_channels
        )
      }
    }
    for (i in seq_along(res$global)) {
      r <- res$global[[i]]
      if (!is.null(r$error)) {
        parse_errors <- c(parse_errors, paste0("Plex ", i, " global: ", r$error))
        ptm_upl_rv$global_statuses[[i]] <- list(error = r$error)
      } else {
        ptm_upl_rv$global_files[[i]]    <- r$struct
        ptm_upl_rv$global_raw[[i]]      <- NULL   # clear raw after successful parse
        ptm_upl_rv$global_statuses[[i]] <- list(
          n_proteins = r$n_proteins,
          n_channels = r$n_channels
        )
      }
    }

    # Expose phospho_files as store$files so tmt_register_design_grid can assign
    # ch_file_idx correctly (maps each channel to its source plex).
    ptm_upl_rv$files <- Filter(Negate(is.null), ptm_upl_rv$phospho_files)

    .refresh_channels()
    .detect_ptm_types()   # fast — uses .peptide_df_cache stored in structs

    if (length(parse_errors) > 0) {
      for (e in parse_errors) showNotification(e, type = "error", duration = 10)
      ptm_upl_rv$wizard_step <- 2L   # stay on step 2 so user can fix
    } else {
      ptm_upl_rv$wizard_step <- 3L
      logger(paste0("[TMT PTM] Parsed ", length(res$phospho), " phospho + ",
                    length(res$global), " global file(s). Advancing to step 3."))
    }
    showModal(.build_modal())
  })

  # ── Open / Reopen wizard ───────────────────────────────────────────────────
  observeEvent(input$tmt_ptm_wiz_open, {
    ptm_upl_rv$wizard_step       <- 1L
    ptm_upl_rv$wizard_modal_open <- TRUE
    showModal(.build_modal())
  })

  observeEvent(input$tmt_ptm_wiz_reopen, {
    ptm_upl_rv$wizard_step       <- 1L
    ptm_upl_rv$wizard_modal_open <- TRUE
    showModal(.build_modal())
  })

  # ── Title (progress bar) ───────────────────────────────────────────────────
  output$tmt_ptm_wiz_title <- renderUI({
    step        <- ptm_upl_rv$wizard_step %||% 1L
    visible_step <- min(step, N_STEPS)
    lbl         <- STEP_LABELS[[visible_step]]
    pct         <- round(visible_step / N_STEPS * 100)
    tagList(
      tags$div(class = "d-flex justify-content-between align-items-center",
        tags$span(icon("flask-vial"), " TMT PTM Upload Wizard"),
        tags$span(class = "text-muted small",
          paste0("Step ", visible_step, " / ", N_STEPS, " \u2014 ", lbl))
      ),
      tags$div(class = "progress mt-2", style = "height:5px;",
        tags$div(class = "progress-bar bg-success", role = "progressbar",
          style = paste0("width:", pct, "%;"),
          `aria-valuenow` = pct, `aria-valuemin` = 0, `aria-valuemax` = 100)
      )
    )
  })

  # ── Content router ─────────────────────────────────────────────────────────
  output$tmt_ptm_wiz_content <- renderUI({
    switch(as.character(ptm_upl_rv$wizard_step %||% 1L),
      "1" = .step1_ui(),
      "2" = .step2_ui(),
      "3" = .step3_ui(),
      "4" = .step4_ui(),
      "5" = .step5_ui(),
      "6" = .step6_ui(),
      "7" = .step7_ui(),
      tags$p("Unknown wizard step.")
    )
  })

  # ── Footer ─────────────────────────────────────────────────────────────────
  output$tmt_ptm_wiz_footer <- renderUI({
    step       <- ptm_upl_rv$wizard_step %||% 1L
    cancel_btn <- modalButton("Cancel")

    back_btn <- if (step > 1L && step < 7L)
      actionButton("tmt_ptm_wiz_back", "\u2190 Back", class = "btn-secondary")
    else NULL

    if (step == 7L) {
      done_cls <- if (ptm_upl_rv$summarization_done) "btn-success" else "btn-secondary"
      done_btn <- actionButton("tmt_ptm_wiz_done", "\u2714 Done",
                               class = done_cls)
      tagList(tags$div(class = "d-flex justify-content-between",
        tags$div(cancel_btn), tags$div(done_btn)))

    } else if (step == 6L) {
      done_cls <- if (ptm_upl_rv$summarization_done) "btn-success" else "btn-secondary"
      done_btn <- actionButton("tmt_ptm_wiz_done", "\u2714 Done",
                               class = done_cls,
                               disabled = if (!ptm_upl_rv$summarization_done) NA else NULL)
      tagList(tags$div(class = "d-flex justify-content-between",
        tags$div(cancel_btn),
        tags$div(back_btn, done_btn)))

    } else {
      next_btn <- actionButton("tmt_ptm_wiz_next", "Next \u2192", class = "btn-primary")
      tagList(tags$div(class = "d-flex justify-content-between",
        tags$div(cancel_btn),
        tags$div(back_btn, next_btn)))
    }
  })

  # ============================================================================
  # STEP 1 — Input Type
  # ============================================================================
  .step1_ui <- function() {
    tagList(
      tags$h5("Select Input Type"),
      tags$p(class = "text-muted small",
        "Choose the software that generated your phosphoproteomics data."),
      fluidRow(
        # PD — enabled and selected
        column(6,
          tags$div(class = "card border-success mb-2",
            style = "cursor:pointer;",
            tags$div(class = "card-body py-3",
              tags$div(class = "d-flex align-items-start",
                tags$div(class = "me-3 mt-1", style = "font-size:1.8rem; color:#27ae60;",
                  icon("database")),
                tags$div(
                  tags$strong("Proteome Discoverer"),
                  tags$span(class = "badge bg-success ms-1 small", "Selected"),
                  tags$p(class = "text-muted small mb-0 mt-1",
                    "Nested CSV export with protein and peptide sections. ",
                    "Upload paired phospho-enriched and global profiling files.")
                )
              )
            )
          )
        )
      ),
      tags$div(class = "alert alert-info mt-1 small py-2 px-3",
        icon("circle-info"),
        " Upload one paired set of files per TMT plex. ",
        "The phospho-enriched file must contain the ",
        tags$strong("'Modifications in Master Proteins'"),
        " column with ptmRS/phosphoRS probabilities.")
    )
  }

  # ============================================================================
  # STEP 2 — Upload Files
  # ============================================================================
  .step2_ui <- function() {
    n <- ptm_upl_rv$n_plexes %||% 1L

    plex_rows <- lapply(seq_len(n), function(i) {
      phos_st <- ptm_upl_rv$phospho_statuses[i][[1L]]
      glob_st <- ptm_upl_rv$global_statuses[i][[1L]]

      phos_badge <- if (!is.null(phos_st)) {
        if (!is.null(phos_st$error))
          tags$div(class = "alert alert-danger py-1 px-2 small mb-0 mt-1",
            icon("circle-xmark"), " ", phos_st$error)
        else if (isTRUE(phos_st$uploaded))
          tags$div(class = "alert alert-info py-1 px-2 small mb-0 mt-1",
            icon("clock"), " ", phos_st$filename, " — will parse on Next")
        else
          tags$div(class = "alert alert-success py-1 px-2 small mb-0 mt-1",
            icon("circle-check"), " ",
            paste0(phos_st$n_peptides, " peptides, ",
                   phos_st$n_channels, " channels"))
      }
      glob_badge <- if (!is.null(glob_st)) {
        if (!is.null(glob_st$error))
          tags$div(class = "alert alert-danger py-1 px-2 small mb-0 mt-1",
            icon("circle-xmark"), " ", glob_st$error)
        else if (isTRUE(glob_st$uploaded))
          tags$div(class = "alert alert-info py-1 px-2 small mb-0 mt-1",
            icon("clock"), " ", glob_st$filename, " — will parse on Next")
        else
          tags$div(class = "alert alert-success py-1 px-2 small mb-0 mt-1",
            icon("circle-check"), " ",
            paste0(glob_st$n_proteins, " proteins, ",
                   glob_st$n_channels, " channels"))
      }

      tags$div(class = "card mb-2",
        tags$div(class = "card-body py-2 px-3",
          tags$strong(paste0("Plex ", i), class = "d-block mb-2"),
          fluidRow(
            column(6,
              tags$small(class = "text-muted d-block mb-1",
                icon("microscope"), " Phospho-Enriched CSV:"),
              fileInput(paste0("tmt_ptm_phos_", i), NULL,
                        accept = c(".csv", "text/csv"), width = "100%"),
              phos_badge
            ),
            column(6,
              tags$small(class = "text-muted d-block mb-1",
                icon("globe"), " Global Profiling CSV:"),
              fileInput(paste0("tmt_ptm_glob_", i), NULL,
                        accept = c(".csv", "text/csv"), width = "100%"),
              glob_badge
            )
          )
        )
      )
    })

    tagList(
      tags$h5("Upload Paired Files"),
      tags$p(class = "text-muted small",
        "For each TMT plex, upload one Phospho-Enriched and one Global Profiling ",
        "PD nested CSV. Both must come from the same plex experiment."),
      tagList(plex_rows),
      if (n < 6L)
        actionButton("tmt_ptm_add_plex", icon("plus"), " Add Plex",
                     class = "btn-outline-secondary btn-sm w-100 mt-1")
    )
  }

  # ── File upload handlers (dynamic for up to 6 plexes) ─────────────────────
  # Files are NOT parsed here — parsing is deferred to when the user clicks
  # Next on step 2. This keeps the UI fully responsive while files are selected.
  for (.i in seq_len(6L)) {
    local({
      pi <- .i

      observeEvent(input[[paste0("tmt_ptm_phos_", pi)]], {
        f <- input[[paste0("tmt_ptm_phos_", pi)]]
        req(f)
        ptm_upl_rv$phospho_raw[[pi]]      <- list(datapath = f$datapath, name = f$name)
        ptm_upl_rv$phospho_statuses[[pi]] <- list(uploaded = TRUE, filename = f$name)
        ptm_upl_rv$phospho_files[[pi]]    <- NULL   # clear any previous parse result
      }, ignoreInit = TRUE)

      observeEvent(input[[paste0("tmt_ptm_glob_", pi)]], {
        f <- input[[paste0("tmt_ptm_glob_", pi)]]
        req(f)
        ptm_upl_rv$global_raw[[pi]]      <- list(datapath = f$datapath, name = f$name)
        ptm_upl_rv$global_statuses[[pi]] <- list(uploaded = TRUE, filename = f$name)
        ptm_upl_rv$global_files[[pi]]    <- NULL   # clear any previous parse result
      }, ignoreInit = TRUE)
    })
  }

  observeEvent(input$tmt_ptm_add_plex, {
    n <- (ptm_upl_rv$n_plexes %||% 1L) + 1L
    ptm_upl_rv$n_plexes <- n
    showModal(.build_modal())
  }, ignoreInit = TRUE)

  .refresh_channels <- function() {
    .per_file <- lapply(Filter(Negate(is.null), ptm_upl_rv$phospho_files),
                        function(f) f$abund_cols %||% character(0))
    all_ch <- unique(unlist(.per_file))
    # Channel identity is the PD header string. If two plexes share a header the
    # columns collide into one and their data compete for the same matrix cell.
    if (length(.per_file) > 1L) {
      .dupes <- unlist(.per_file)[duplicated(unlist(.per_file))]
      if (length(.dupes) > 0) {
        showNotification(paste0(
          "Channel names are shared between plexes (", length(unique(.dupes)),
          " duplicated). Each plex must export distinct channel column names ",
          "(different F numbers or sample labels); otherwise their measurements ",
          "cannot be told apart. First duplicate: ", unique(.dupes)[1]),
          type = "error", duration = 25)
        logger(paste0("[TMT PTM] Duplicate channel names across plexes: ",
                      paste(unique(.dupes), collapse = ", ")))
      }
    }
    ptm_upl_rv$all_channels <- all_ch
  }

  .detect_ptm_types <- function() {
    skip_pat  <- "^(TMT|Carbamidomethyl|Oxidation)"
    all_types <- character(0)
    for (f in Filter(Negate(is.null), ptm_upl_rv$phospho_files)) {
      # Use cached peptide df if available (avoids re-parsing large files)
      pep_df <- f$.peptide_df_cache
      if (is.null(pep_df))
        pep_df <- tryCatch(parse_pd_peptide_data(f), error = function(e) NULL)
      if (is.null(pep_df)) next
      mod_col <- grep("^Modifications$", colnames(pep_df), value = TRUE)[1]
      if (is.na(mod_col)) next
      mod_strings <- pep_df[[mod_col]]
      mod_strings <- mod_strings[!is.na(mod_strings) & nchar(mod_strings) > 0]
      hits      <- unlist(regmatches(mod_strings,
                    gregexpr("\\d+x([A-Za-z][A-Za-z0-9]*)", mod_strings)))
      mod_names <- unique(gsub("^\\d+x", "", hits))
      all_types <- unique(c(all_types, mod_names))
    }
    ptm_types <- sort(all_types[!grepl(skip_pat, all_types, ignore.case = TRUE)])
    if (length(ptm_types) == 0) ptm_types <- "Phospho"
    ptm_upl_rv$detected_ptms <- ptm_types
    # This release analyzes phosphosites only (the site-level pipeline is
    # Phospho-specific downstream), so the selection is always Phospho.
    ptm_upl_rv$selected_ptms <- "Phospho"
  }

  # ============================================================================
  # STEP 3 — PTM Type Selection
  # ============================================================================
  .step3_ui <- function() {
    ptm_types <- ptm_upl_rv$detected_ptms
    if (length(ptm_types) == 0) ptm_types <- "Phospho"
    other_types <- setdiff(ptm_types, "Phospho")

    tagList(
      tags$h5("Target PTM Type"),
      tags$div(class = "alert alert-info py-2 px-3 small mb-3",
        icon("circle-info"),
        " This release analyzes phosphosites. Peptides carrying other ",
        "modifications from the ",
        "phospho-enriched file will not be included in the site-level matrix."),
      radioButtons("tmt_ptm_type_selector",
        label = "PTM type:",
        choices  = "Phospho",
        selected = "Phospho"),
      if (!("Phospho" %in% ptm_types))
        tags$div(class = "alert alert-warning py-2 px-3 small mt-2",
          icon("triangle-exclamation"),
          " No phospho modifications were found in the phospho-enriched file. ",
          "Check that the file is the phosphopeptide export."),
      if (length(other_types) > 0L)
        tags$p(class = "text-muted small mt-2",
          icon("circle-info"),
          " Also detected: ", paste(other_types, collapse = ", "),
          ". Other modification types are not analyzed in this release.")
    )
  }

  observeEvent(input$tmt_ptm_type_selector, {
    sel <- input$tmt_ptm_type_selector
    if (!is.null(sel) && length(sel) > 0) ptm_upl_rv$selected_ptms <- sel
  }, ignoreInit = TRUE)

  # ============================================================================
  # STEP 4 — Plex Design
  # ============================================================================
  .step4_ui <- function() {
    tagList(
      tags$h5("Experimental Design"),
      tags$p(class = "text-muted small mb-3",
        "Select your TMT plex configuration. Clicking Next will open the channel ",
        "assignment grid where you assign conditions to each TMT channel."),
      radioButtons("tmt_ptm_plex_type", NULL,
        choiceNames = list(
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("vial"), " Single Plex"),
            tags$br(),
            tags$small(class = "text-muted",
              "All samples in one TMT experiment. No bridge channel needed.")),
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("vials"), " Multiple Plexes"),
            tags$br(),
            tags$small(class = "text-muted",
              "Samples split across multiple TMT experiments. ",
              "Optionally include a bridge/reference channel for IRS normalization."))
        ),
        choiceValues = c("single", "multi"),
        selected = "single"
      ),
      tags$div(class = "alert alert-warning py-2 px-3 small mt-2",
        icon("triangle-exclamation"),
        " After clicking Next, the wizard will temporarily close while you assign ",
        "conditions to each TMT channel. It will reopen automatically when done.")
    )
  }

  # ============================================================================
  # STEP 5 — Probability Threshold
  # ============================================================================
  .step5_ui <- function() {
    cur_val <- ptm_upl_rv$prob_threshold %||% 75
    tagList(
      tags$h5("PTM Site Localization Probability Threshold"),
      tags$p(class = "text-muted small",
        "Filter phosphopeptides by their ptmRS or phosphoRS localization probability. ",
        "Only sites whose maximum probability across all samples meets this threshold ",
        "will be included in the quantification matrix."),
      sliderInput("tmt_ptm_prob_threshold",
        label = "Minimum localization probability (%):",
        min = 0, max = 100, value = cur_val, step = 5),
      numericInput("tmt_ptm_min_valid",
        label = "Minimum valid values required per condition:",
        value = ptm_upl_rv$min_valid_per_cond %||% 2L,
        min = 1, max = 10, step = 1),
      tags$div(class = "alert alert-secondary py-2 px-3 small mt-2",
        icon("circle-info"),
        " Differential testing requires this many valid values in ",
        tags$strong("both"), " compared conditions. With 2 replicates per ",
        "condition a value of 2 means no missing values are tolerated at all ",
        "(complete-case analysis), and sites failing it cannot reach ",
        "significance in the presence/absence analysis either."),
      numericInput("tmt_ptm_interference",
        label = "Maximum co-isolation interference (%):",
        value = ptm_upl_rv$interference_threshold %||% 50,
        min = 0, max = 100, step = 5),
      tags$div(class = "alert alert-secondary py-2 px-3 small mt-2",
        icon("circle-info"),
        " Applied to the Isolation Interference column when the export contains ",
        "it. Peptide groups with no recorded value are kept and counted ",
        "separately. If the column is absent, no filter is applied and the ",
        "thresholds set inside Proteome Discoverer are not recorded in the ",
        "export - obtain them from the facility that ran the search."),
      tags$div(class = "alert alert-info py-2 px-3 small mt-2",
        icon("circle-info"),
        " Recommended: 75% (default). This is applied to the ",
        tags$strong("'Modifications in Master Proteins'"),
        " column from Proteome Discoverer. Lower values retain more sites but ",
        "may include ambiguously localized phosphorylations.")
    )
  }

  observeEvent(input$tmt_ptm_prob_threshold, {
    val <- input$tmt_ptm_prob_threshold
    if (!is.null(val) && is.numeric(val))
      ptm_upl_rv$prob_threshold <- val
  }, ignoreInit = TRUE)

  observeEvent(input$tmt_ptm_min_valid, {
    val <- input$tmt_ptm_min_valid
    if (!is.null(val) && is.numeric(val) && is.finite(val) && val >= 1)
      ptm_upl_rv$min_valid_per_cond <- as.integer(val)
  }, ignoreInit = TRUE)

  observeEvent(input$tmt_ptm_interference, {
    val <- input$tmt_ptm_interference
    if (!is.null(val) && is.numeric(val) && is.finite(val) && val >= 0)
      ptm_upl_rv$interference_threshold <- as.numeric(val)
  }, ignoreInit = TRUE)

  # ============================================================================
  # STEP 6 — Summary & Summarize
  # ============================================================================
  .step6_ui <- function() {
    ann        <- ptm_upl_rv$annotation_map
    conds      <- ptm_upl_rv$conditions %||% character(0)
    design_opt <- ptm_upl_rv$design_option
    n_ch       <- if (!is.null(ann))
      nrow(ann[!ann$Condition %in% c("Norm", "") & !is.na(ann$Condition), ])
    else 0L
    prob_thresh <- ptm_upl_rv$prob_threshold %||% 75
    ptm_sel     <- ptm_upl_rv$selected_ptms %||% "Phospho"
    n_phos      <- length(Filter(Negate(is.null), ptm_upl_rv$phospho_files))
    n_glob      <- length(Filter(Negate(is.null), ptm_upl_rv$global_files))

    plex_label <- switch(as.character(as.integer(design_opt %||% 1L)),
      "1" = "Single Plex (median centering)",
      "2" = "Multi-Plex + IRS Bridge Normalization",
      "3" = "Multi-Plex (Plex as Run Covariate)",
      "Single Plex (median centering)"
    )

    pipeline_steps <- tags$ol(class = "small mb-0",
      tags$li("Log2 transformation of raw TMT channel abundances (zeros \u2192 NA)"),
      tags$li(
        if (identical(as.integer(design_opt), 2L))
          "Column-wise median centering, then IRS bridge-channel alignment, applied independently to the phospho-enriched and global protein matrices"
        else
          "Column-wise median centering applied independently to the phospho-enriched and global protein matrices"
      ),
      tags$li(paste0("ptmRS localization filter (site maximum \u2265 ", prob_thresh,
        "%; qualifying sites keep all their observations)")),
      tags$li("Site-level rollup: peptidoform intensities summed within each site and channel, then log2"),
      tags$li(paste0("Missingness filter: a site is kept when it has at least ",
                     as.integer(ptm_upl_rv$min_valid_per_cond %||% 2L),
                     " valid values in at least one condition")),
      tags$li("No plex-coverage filter: a site is not required to appear in a minimum number of plexes"),
      tags$li("Stoichiometric adjustment: log\u2082(phospho site) \u2212 log\u2082(global protein)"),
      if (identical(as.integer(design_opt), 3L))
        tags$li("Plex run covariate added to limma design matrix")
      else NULL,
      tags$li("limma eBayes: lmFit \u2192 eBayes(trend=TRUE, robust=TRUE) \u2192 topTable(adj.BH)")
    )

    tagList(
      tags$h5("Summary & Summarize"),
      if (identical(as.integer(design_opt), 2L))
        tags$div(class = "alert alert-info py-2 px-3 small mb-3",
          icon("circle-info"),
          tags$strong(" IRS Bridge design selected: "),
          "IRS corrects between-plex technical variation but does ",
          tags$strong("not"), " correct within-plex channel loading differences, ",
          "so column-median centering is applied first and IRS second. A site with no ",
          "bridge value in a given plex cannot be aligned there; that site-by-plex block ",
          "is set missing and the site is kept on the plexes where it can be corrected, ",
          "rather than the whole site being discarded. Both steps are ",
          "performed here on the raw ",
          tags$em("Abundance:"),
          " columns - leave normalization off in Proteome Discoverer, as any ",
          "normalization applied there is written to separate columns that are ",
          "not read."
        ),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body py-2",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$tbody(
              tags$tr(tags$td(tags$strong("Source")),
                      tags$td("Proteome Discoverer (Nested CSV)")),

              tags$tr(tags$td(tags$strong("Phospho-enriched files")),
                      tags$td(paste0(n_phos, " file(s)"))),
              tags$tr(tags$td(tags$strong("Global profiling files")),
                      tags$td(paste0(n_glob, " file(s)"))),
              tags$tr(tags$td(tags$strong("Target PTMs")),
                      tags$td(paste(ptm_sel, collapse = ", "))),
              tags$tr(tags$td(tags$strong("Design")),
                      tags$td(plex_label)),
              tags$tr(tags$td(tags$strong("Conditions")),
                      tags$td(if (length(conds) > 0) paste(conds, collapse = ", ")
                              else tags$em(class = "text-muted", "not set"))),
              tags$tr(tags$td(tags$strong("Biological channels")),
                      tags$td(n_ch)),
              tags$tr(tags$td(tags$strong("Prob. threshold")),
                      tags$td(paste0(prob_thresh, "%"))),
              tags$tr(tags$td(tags$strong("DE engine")),
                      tags$td("limma eBayes (trend=TRUE, robust=TRUE)"))
            )
          )
        )
      ),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body py-2",
          tags$h6(class = "card-title mb-2", icon("sitemap"), " Pipeline Steps"),
          pipeline_steps
        )
      ),
      actionButton("tmt_ptm_wiz_summarize_btn",
        tagList(icon("gears"), " Summarize"),
        class = "btn-primary w-100 mb-2"),
      uiOutput("tmt_ptm_wiz_summarize_status")
    )
  }

  output$tmt_ptm_wiz_summarize_status <- renderUI({
    if (ptm_upl_rv$summarization_done) {
      tags$div(class = "alert alert-success py-2 px-3 small mt-1",
        icon("circle-check"), " Summarization complete: ",
        tags$strong(ptm_upl_rv$n_sites), " PTM sites, ",
        tags$strong(ptm_upl_rv$n_sites_adjusted), " stoichiometry-adjusted rows."
      )
    } else {
      tags$p(class = "text-muted small mt-1",
        "Click Summarize to build the site-level quantification matrices.")
    }
  })

  # ============================================================================
  # STEP 7 — Review
  # ============================================================================
  .step7_ui <- function() {
    cache <- rv$dia_ptm_master_cache
    if (is.null(cache)) {
      return(tags$div(class = "alert alert-warning",
        icon("triangle-exclamation"),
        " Summarization has not been run. Please go back and click Summarize."))
    }

    tagList(
      tags$h5("Review & Finalize"),
      tags$div(class = "alert alert-success py-2 px-3 mb-3",
        icon("circle-check"),
        " Summarization complete. Click ",
        tags$strong("\u2714 Done"),
        " to send data to all PTM analysis modules."),
      tags$div(class = "card",
        tags$div(class = "card-body py-2",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$tbody(
              tags$tr(tags$td(tags$strong("Pipeline")),
                      tags$td("TMT PTM \u2014 limma eBayes, stoichiometric adjustment")),
              tags$tr(tags$td(tags$strong("PTM sites (output)")),
                      tags$td(cache$n_sites_output)),
              tags$tr(tags$td(tags$strong("Stoichiometry-adjusted rows")),
                      tags$td(cache$n_sites_adjusted)),
              tags$tr(tags$td(tags$strong("Sites dropped, plex coverage")),
                      tags$td(if (identical(as.numeric(cache$site_min_plex_frac %||% 0), 0))
                                "not applied"
                              else cache$n_sites_dropped_plex %||% 0L)),
              tags$tr(tags$td(tags$strong("IRS: site-plex blocks set missing")),
                      tags$td(if (is.na(cache$irs_n_blanked %||% NA_integer_))
                                "n/a" else cache$irs_n_blanked)),
              tags$tr(tags$td(tags$strong("Sites dropped after IRS re-check")),
                      tags$td(cache$n_sites_dropped_post_irs %||% 0L)),
              tags$tr(tags$td(tags$strong("Shared peptides reassigned (razor)")),
                      tags$td(cache$n_razor_reassigned %||% 0L)),
              tags$tr(tags$td(tags$strong("Sites with no parent protein")),
                      tags$td(cache$n_no_parent_protein %||% 0L)),
              tags$tr(tags$td(tags$strong("Biological channels")),
                      tags$td(length(cache$matched_samples))),
              tags$tr(tags$td(tags$strong("Conditions")),
                      tags$td(paste(names(cache$cond_groups), collapse = ", "))),
              tags$tr(tags$td(tags$strong("Prob. threshold applied")),
                      tags$td(paste0(cache$prob_threshold, "%"))),
              tags$tr(tags$td(tags$strong("Normalization")),
                      tags$td(cache$normalization))
            )
          )
        )
      ),
      tags$p(class = "text-muted small mt-3",
        icon("arrow-right"), " Navigate to ",
        tags$strong("PTM DE"), ", ",
        tags$strong("PTM Trends"), ", or ",
        tags$strong("PTM Single"),
        " in the sidebar to begin analysis.")
    )
  }

  # ============================================================================
  # NAVIGATION: Next button
  # ============================================================================
  observeEvent(input$tmt_ptm_wiz_next, {
    step <- ptm_upl_rv$wizard_step %||% 1L

    if (step == 1L) {
      # Step 1: Input type (PD always selected — no validation needed)
      ptm_upl_rv$wizard_step <- 2L
      showModal(.build_modal())

    } else if (step == 2L) {
      # Validate plex 1 files are present (raw-uploaded or already parsed)
      phos_ok <- !is.null(ptm_upl_rv$phospho_raw[1][[1L]]) ||
                 !is.null(ptm_upl_rv$phospho_files[1][[1L]])
      glob_ok <- !is.null(ptm_upl_rv$global_raw[1][[1L]])  ||
                 !is.null(ptm_upl_rv$global_files[1][[1L]])
      if (!phos_ok) {
        showNotification("Upload a Phospho-Enriched CSV for Plex 1 first.",
                         type = "warning")
        return()
      }
      if (!glob_ok) {
        showNotification("Upload a Global Profiling CSV for Plex 1 first.",
                         type = "warning")
        return()
      }

      # Collect files that still need parsing
      phospho_raws <- Filter(Negate(is.null), ptm_upl_rv$phospho_raw)
      global_raws  <- Filter(Negate(is.null), ptm_upl_rv$global_raw)

      if (length(phospho_raws) == 0 && length(global_raws) == 0) {
        # Everything already parsed (e.g. user went Back and returned)
        ptm_upl_rv$wizard_step <- 3L
        showModal(.build_modal())
        return()
      }

      # Launch async parse — closes the modal, shows overlay
      n_pl <- max(length(phospho_raws), length(global_raws))
      removeModal()
      show_async_overlay(
        session,
        label  = paste0("Parsing ", n_pl, " plex", if (n_pl > 1) "es" else "", "\u2026"),
        detail = "Reading PD nested CSV files. This may take a moment for large files."
      )
      ptm_upl_rv$parsing_in_progress <- TRUE
      tmt_ptm_parse_task$invoke(phospho_raws, global_raws)

    } else if (step == 3L) {
      sel <- input$tmt_ptm_type_selector
      if (is.null(sel) || length(sel) == 0) {
        showNotification("Select at least one PTM type to continue.", type = "warning")
        return()
      }
      ptm_upl_rv$selected_ptms <- sel
      ptm_upl_rv$wizard_step   <- 4L
      showModal(.build_modal())

    } else if (step == 4L) {
      # Step 4: Close wizard, open design grid sub-modal
      plex_type <- input$tmt_ptm_plex_type %||% "single"
      if (length(ptm_upl_rv$all_channels) == 0) {
        showNotification("No TMT channels detected. Please upload files in Step 2 first.",
                         type = "error")
        return()
      }
      removeModal()
      ptm_upl_rv$wizard_modal_open <- TRUE

      if (identical(plex_type, "single")) {
        ptm_upl_rv$n_plexes      <- 1L
        ptm_upl_rv$has_bridge    <- FALSE
        ptm_upl_rv$design_option <- 1L
        if (length(ptm_upl_rv$design_conds) == 0)
          ptm_upl_rv$design_conds <- character(0)
        .dg$show_design_grid()
      } else {
        n_f <- length(Filter(Negate(is.null), ptm_upl_rv$phospho_files))
        ptm_upl_rv$n_plexes         <- max(n_f, 2L)
        ptm_upl_rv$has_bridge       <- FALSE
        ptm_upl_rv$design_conds     <- character(0)
        ptm_upl_rv$grid_assignments <- NULL
        ptm_upl_rv$paint_cond       <- NULL
        Sys.sleep(0.15)
        .dg$show_design_grid()
      }
      # Wizard will reopen at Step 5 via on_design_finished callback

    } else if (step == 5L) {
      ptm_upl_rv$prob_threshold <- input$tmt_ptm_prob_threshold %||% 75
      ptm_upl_rv$min_valid_per_cond <-
        as.integer(input$tmt_ptm_min_valid %||% 2L)
      ptm_upl_rv$interference_threshold <-
        as.numeric(input$tmt_ptm_interference %||% 50)
      ptm_upl_rv$wizard_step   <- 6L
      showModal(.build_modal())

    } else if (step == 6L) {
      if (!ptm_upl_rv$summarization_done) {
        showNotification("Please click Summarize before proceeding.", type = "warning")
        return()
      }
      ptm_upl_rv$wizard_step <- 7L
      showModal(.build_modal())
    }
  })

  # ============================================================================
  # NAVIGATION: Back button
  # ============================================================================
  observeEvent(input$tmt_ptm_wiz_back, {
    step <- ptm_upl_rv$wizard_step %||% 1L
    new_step <- switch(as.character(step),
      "2" = 1L, "3" = 2L, "4" = 3L, "5" = 4L, "6" = 5L, 1L)
    ptm_upl_rv$wizard_step <- new_step
    showModal(.build_modal())
  })

  # ============================================================================
  # SUMMARIZE BUTTON
  # ============================================================================
  observeEvent(input$tmt_ptm_wiz_summarize_btn, {
    ann <- ptm_upl_rv$annotation_map
    if (is.null(ann) || nrow(ann) == 0) {
      showNotification(
        "Please complete the experimental design first (Step 4).",
        type = "warning")
      return()
    }

    phospho_files <- Filter(Negate(is.null), ptm_upl_rv$phospho_files)
    global_files  <- Filter(Negate(is.null), ptm_upl_rv$global_files)

    if (length(phospho_files) == 0) {
      showNotification("No phospho-enriched files loaded.", type = "warning")
      return()
    }
    if (length(global_files) == 0) {
      showNotification("No global profiling files loaded.", type = "warning")
      return()
    }

    # ---- Design validation -------------------------------------------------
    .d_opt  <- as.integer(ptm_upl_rv$design_option %||% 1L)
    .plexes <- unique(ann$Mixture[!is.na(ann$Mixture)])
    .bio    <- ann[!is.na(ann$Condition) &
                   !ann$Condition %in% c("Norm", ""), , drop = FALSE]

    if (.d_opt == 3L && length(.plexes) >= 2L && nrow(.bio) > 0L) {
      .tab  <- table(.bio$Condition, .bio$Mixture)
      .conf <- rownames(.tab)[rowSums(.tab > 0) < 2L]
      if (length(.conf) > 0) {
        showNotification(paste0(
          "Plex and condition are confounded: condition(s) ",
          paste(.conf, collapse = ", "),
          " appear in only one plex. A plex covariate cannot separate plex from ",
          "condition in this design, so the comparison cannot be corrected for ",
          "plex. Re-assign channels in Step 4, or use a bridge-channel design."),
          type = "error", duration = 25)
        logger(paste0("[TMT PTM Wizard] Design validation failed - plex/condition ",
                      "confounded for: ", paste(.conf, collapse = ", ")))
        return()
      }
    }

    if (length(.plexes) >= 2L && nrow(.bio) > 0L) {
      .tab2 <- table(.bio$Condition, .bio$Mixture)
      if (length(unique(as.vector(.tab2))) > 1L)
        showNotification(paste0(
          "Conditions are unevenly represented across plexes. The analysis will ",
          "proceed, but plex effects are estimated less precisely."),
          type = "warning", duration = 12)
    }
    # Channel-name collision across plexes. Channel identity is the PD header
    # string, so two plexes exporting the same header collapse into one matrix
    # column and their measurements compete for the same cell. Blocking, not
    # advisory: the result would be silently wrong.
    .ch_per_file <- lapply(phospho_files, function(.f) .f$abund_cols %||% character(0))
    if (length(.ch_per_file) > 1L) {
      .all_ch <- unlist(.ch_per_file)
      .dupes  <- unique(.all_ch[duplicated(.all_ch)])
      if (length(.dupes) > 0) {
        showNotification(paste0(
          "Cannot summarize: ", length(.dupes), " channel name(s) are shared ",
          "between plexes, so their measurements cannot be told apart. Each plex ",
          "must export distinct channel column names (different F numbers or ",
          "sample labels). First duplicate: ", .dupes[1]),
          type = "error", duration = 30)
        logger(paste0("[TMT PTM Wizard] Summarization blocked - duplicate channel ",
                      "names across plexes: ", paste(.dupes, collapse = ", ")))
        return()
      }
    }
    # ---- End design validation ---------------------------------------------

    # ---- Header vs grid assignment cross-check -----------------------------
    # Header sample labels and grid condition names are different vocabularies
    # ("TPA 0 hr" vs "0h"), so string equality is meaningless. What is checkable
    # is the PARTITION: channels sharing a header label should share a grid
    # condition. A header group splitting across two grid conditions means one
    # channel was assigned differently from its siblings. Merging several header
    # groups into one grid condition is legitimate and is not reported.
    .hdr <- do.call(rbind, lapply(phospho_files, function(.f) {
      if (is.null(.f$channel_info)) NULL else .f$channel_info
    }))
    if (!is.null(.hdr) && nrow(.hdr) > 0) {
      .hmap <- setNames(trimws(as.character(.hdr$Condition)),
                        as.character(.hdr$ColumnName))
      .cmp  <- .bio[.bio$ColumnName %in% names(.hmap), , drop = FALSE]
      .hlab <- unname(.hmap[.cmp$ColumnName])
      # When the header carries no parseable condition field, sub() returns the
      # whole column name and every channel forms its own group. No split is
      # possible, so the check is silent rather than noisy.
      .parsed <- nrow(.cmp) > 0 && length(unique(.hlab)) < nrow(.cmp)
      if (!.parsed) {
        logger(paste0("[TMT PTM Wizard] Header conditions not parseable as groups ",
                      "(", length(unique(.hlab)), " distinct labels for ",
                      nrow(.cmp), " channels); partition check skipped."))
      } else {
        .split <- character(0)
        for (.h in unique(.hlab)) {
          .g <- unique(trimws(as.character(.cmp$Condition[.hlab == .h])))
          .g <- .g[!is.na(.g) & nzchar(.g)]
          if (length(.g) > 1L)
            .split <- c(.split, paste0("'", .h, "' -> ",
                                       paste(sprintf("'%s'", .g), collapse = " + ")))
        }
        if (length(.split) > 0) {
          showNotification(paste0(
            length(.split), " sample group(s) from the file header were split ",
            "across more than one condition in Step 4. This usually means a ",
            "channel was assigned to the wrong condition. Check Step 4 before ",
            "continuing: ", paste(utils::head(.split, 5), collapse = "; "),
            if (length(.split) > 5) " ..." else ""),
            type = "warning", duration = 25)
          logger(paste0("[TMT PTM Wizard] Header/grid partition split on ",
                        length(.split), " header group(s): ",
                        paste(.split, collapse = " | ")))
        }
      }
    }
    # ---- End cross-check ---------------------------------------------------

    prob_threshold <- ptm_upl_rv$prob_threshold %||% 75
    design_opt     <- as.integer(ptm_upl_rv$design_option %||% 1L)
    target_ptms    <- ptm_upl_rv$selected_ptms %||% "Phospho"

    withProgress(message = "Summarizing TMT PTM data...", value = 0.1, {
      incProgress(0.1, detail = "Parsing sites and building matrices...")

      .ptm_debug_stack <- NULL
      result <- tryCatch(
        withCallingHandlers(
          tmt_ptm_limma_summarize(
            ptm_parsed_files  = phospho_files,
            prot_parsed_files = global_files,
            annotation_map    = ann,
            target_ptms        = target_ptms,
            prob_threshold     = prob_threshold,
            min_valid_per_cond = as.integer(ptm_upl_rv$min_valid_per_cond %||% 2L),
            interference_threshold =
              as.numeric(ptm_upl_rv$interference_threshold %||% 50),
            design_opt         = design_opt,
            subject_pairing   = NULL,
            rv                = list(
              pd_has_gn  = isTRUE(rv$pd_has_gn),
              pd_gene_col = rv$pd_gene_col
            ),
            logger = logger
          ),
          error = function(.e) {
            .ptm_debug_stack <<- tryCatch({
              calls <- sys.calls()
              fns <- vapply(calls, function(x)
                tryCatch(deparse(x[[1L]], nlines = 1L)[1L], error = function(e2) "?"),
                character(1L))
              paste(fns, collapse = " -> ")
            }, error = function(e2) "<stack unavailable>")
          }
        ),
        error = function(e) list(error = conditionMessage(e))
      )
      incProgress(0.8, detail = "Finalizing...")
    })

    if (!is.null(result$error)) {
      showNotification(paste0("Summarization failed: ", result$error),
                       type = "error", duration = 15)
      logger(paste0("[TMT PTM Wizard] Summarization error: ", result$error))
      if (!is.null(.ptm_debug_stack))
        logger(paste0("[TMT PTM Wizard] Call stack: ", .ptm_debug_stack))
      return()
    }

    rv$dia_ptm_master_cache            <- result
    ptm_upl_rv$summarization_done      <- TRUE
    ptm_upl_rv$n_sites                 <- result$n_sites_output %||% 0L
    ptm_upl_rv$n_sites_adjusted        <- result$n_sites_adjusted %||% 0L

    logger(paste0("[TMT PTM Wizard] Summarization complete — ",
                  result$n_sites_output, " sites, ",
                  result$n_sites_adjusted, " adjusted rows, ",
                  length(result$matched_samples), " channels."))
    showModal(.build_modal())
  })

  # ============================================================================
  # DONE button — finalize and route data to downstream modules
  # ============================================================================
  observeEvent(input$tmt_ptm_wiz_done, {
    cache <- rv$dia_ptm_master_cache
    if (is.null(cache)) {
      showNotification("Please run Summarize first.", type = "warning")
      return()
    }

    removeModal()
    ptm_upl_rv$wizard_modal_open <- FALSE
    ptm_upl_rv$wizard_completed  <- TRUE

    # Set shared reactive values
    rv$tmt_ptm_wizard_completed <- TRUE
    rv$data_source              <- "tmt_ptm"

    # Populate rv$ptm_data for downstream sync (server_ptm_de.R sync observer)
    rv$ptm_data$design_option <- "tmt_ptm"
    rv$ptm_data$conditions    <- sort(unique(cache$design_df$condition))
    rv$ptm_data$data_format   <- "tmt_ptm"
    rv$ptm_data$n_files       <- length(ptm_upl_rv$phospho_files) * 2L
    rv$ptm_data$gene_lookup   <- NULL
    rv$ptm_data$ptm_long      <- NULL
    rv$ptm_data$protein_long  <- NULL
    rv$ptm_data$annotation_map <- ptm_upl_rv$annotation_map
    # Trigger sync — set source last so the observer fires with all fields ready
    rv$ptm_data$source        <- NULL
    rv$ptm_data$source        <- "tmt_ptm_upload"

    # Update provenance — write to $dia_ptm so render_dia_ptm_de_methods() reads it
    # (section = "dia_ptm_de" in server_ptm_de.R reads prov$dia_ptm, not prov$ptm_de)
    rv$provenance$dia_ptm <- list(
      acquisition_mode  = "TMT_PD",
      software          = "Proteome Discoverer (TMT)",
      n_sites_input     = cache$n_sites_input,
      n_sites_output    = cache$n_sites_output,
      n_sites_adjusted  = cache$n_sites_adjusted,
      n_shared_expanded = cache$n_shared_expanded %||% 0L,
      prob_threshold    = as.numeric(ptm_upl_rv$prob_threshold %||% 75),
      interference_filtered  = isTRUE(cache$interference_filtered),
      interference_threshold = cache$interference_threshold %||% NA_real_,
      n_interf_removed       = cache$n_interf_removed %||% 0L,
      n_interf_kept_na       = cache$n_interf_kept_na %||% 0L,
      n_interf_kept_pass     = cache$n_interf_kept_pass %||% 0L,
      min_valid_per_cond     = as.integer(ptm_upl_rv$min_valid_per_cond %||% 2L),
      condition_names   = names(cache$cond_groups),
      n_plexes          = as.integer(ptm_upl_rv$n_plexes %||% 1L),
      design_option     = as.integer(ptm_upl_rv$design_option %||% 1L),
      selected_ptms     = ptm_upl_rv$selected_ptms %||% "Phospho",
      normalization     = cache$normalization,
      pairing           = !is.null(ptm_upl_rv$subject_pairing),
      # ---- WP2: everything the methods text must read rather than assume ----
      # Input provenance
      abundance_source       = cache$abundance_source %||% NA_character_,
      # Peptide-group identification confidence
      n_pep_groups_input     = cache$n_pep_groups_input     %||% NA_integer_,
      n_pep_groups_conf_kept = cache$n_pep_groups_conf_kept %||% NA_integer_,
      pep_confidence_basis   = cache$pep_confidence_basis   %||% NA_character_,
      pep_qvalue_max         = cache$pep_qvalue_max         %||% NA_real_,
      prot_qvalue_max        = cache$prot_qvalue_max        %||% NA_real_,
      # Localization: observation-level detail behind the site-level decision
      n_obs_collected           = cache$n_obs_collected           %||% NA_integer_,
      n_obs_below_own_threshold = cache$n_obs_below_own_threshold %||% NA_integer_,
      n_obs_no_own_prob         = cache$n_obs_no_own_prob         %||% NA_integer_,
      n_sites_passing_prob      = cache$n_sites_passing_prob      %||% NA_integer_,
      # Shared-peptide handling (razor; no expansion on this route)
      n_razor_reassigned     = cache$n_razor_reassigned  %||% 0L,
      n_shared_peptide_sites = cache$n_shared_peptide_sites %||% 0L,
      # Rollup
      n_sites_multi_peptidoform = cache$n_sites_multi_peptidoform %||% NA_integer_,
      # Missingness and plex coverage
      n_sites_after_valid_filter = cache$n_sites_after_valid_filter %||% NA_integer_,
      n_sites_dropped_plex   = cache$n_sites_dropped_plex %||% 0L,
      n_plexes_detected      = cache$n_plexes_detected    %||% NA_integer_,
      site_min_plex_n        = cache$site_min_plex_n      %||% NA_integer_,
      n_sites_dropped_post_irs = cache$n_sites_dropped_post_irs %||% 0L,
      irs_policy               = cache$irs_policy             %||% NA_character_,
      irs_n_blanked            = cache$irs_n_blanked          %||% NA_integer_,
      irs_n_plex_no_bridge     = cache$irs_n_plex_no_bridge   %||% NA_integer_,
      prot_irs_n_blanked       = cache$prot_irs_n_blanked     %||% NA_integer_,
      site_min_plex_frac     = cache$site_min_plex_frac   %||% NA_real_,
      min_plex_frac          = cache$min_plex_frac        %||% NA_real_,
      prot_min_plex_frac     = cache$prot_min_plex_frac   %||% NA_real_,
      # Normalization and adjustment
      normalization_basis    = cache$normalization_basis  %||% NA_character_,
      adj_centering_delta    = cache$adj_centering_delta,
      n_no_parent_protein    = cache$n_no_parent_protein  %||% 0L,
      # ---- end WP2 additions ----
      # Channels marked "Exclude" in the design grid (absent from annotation_map)
      excluded_channels      = ptm_upl_rv$excluded_channels %||% character(0),
      timestamp         = as.character(Sys.time())
    )

    logger(paste0("[TMT PTM Wizard] Completed. Sites=", cache$n_sites_output,
                  ", Channels=", length(cache$matched_samples),
                  ", Design=", ptm_upl_rv$design_option,
                  ". rv$ptm_data$source <- 'tmt_ptm_upload'."))
  })

  # ============================================================================
  # SIDEBAR BADGES (rendered post-wizard completion)
  # ============================================================================
  output$tmt_ptm_wiz_sidebar_badges <- renderUI({
    conds      <- ptm_upl_rv$conditions %||% character(0)
    n_sites    <- ptm_upl_rv$n_sites %||% 0L
    design_opt <- ptm_upl_rv$design_option %||% 1L
    ann        <- ptm_upl_rv$annotation_map

    plex_str <- switch(as.character(as.integer(design_opt)),
      "1" = "Single Plex",
      "2" = "Multi-Plex + IRS",
      "3" = "Multi-Plex (Run Covariate)",
      "Single Plex"
    )
    n_bio <- if (!is.null(ann))
      nrow(ann[!ann$Condition %in% c("Norm","") & !is.na(ann$Condition), ])
    else 0L

    tags$div(class = "card mb-2",
      tags$div(class = "card-body py-2 px-2",
        tags$h6(class = "card-title mb-1",
          icon("circle-check"), " TMT PTM Upload Complete"),
        tags$table(
          class = "table table-sm table-borderless mb-0",
          style = "font-size:0.82em;",
          tags$tbody(
            tags$tr(tags$td(class = "text-muted pe-2", "Design:"),
                    tags$td(plex_str)),
            tags$tr(tags$td(class = "text-muted pe-2", "Route:"),
                    tags$td("limma eBayes")),
            tags$tr(tags$td(class = "text-muted pe-2", "Sites:"),
                    tags$td(n_sites)),
            tags$tr(tags$td(class = "text-muted pe-2", "Channels:"),
                    tags$td(paste0(n_bio, " biological"))),
            tags$tr(tags$td(class = "text-muted pe-2", "Conditions:"),
                    tags$td(paste(conds, collapse = ", ")))
          )
        )
      )
    )
  })

  # ============================================================================
  # MAIN CONTENT STATUS (rendered in the tabPanel area)
  # ============================================================================
  output$tmt_ptm_upload_main_status <- renderUI({
    if (!isTRUE(rv$tmt_ptm_wizard_completed)) {
      tags$div(class = "alert alert-info py-2 px-3 mt-3",
        icon("wand-magic-sparkles"),
        " Open the ",
        tags$strong("TMT PTM Upload Wizard"),
        " from the sidebar to upload and configure your TMT phosphoproteomics data.")
    } else {
      cache <- rv$dia_ptm_master_cache
      if (is.null(cache)) return(NULL)
      tags$div(class = "alert alert-success py-2 px-3 mt-3",
        icon("circle-check"), " Upload complete: ",
        tags$strong(cache$n_sites_output, " phosphosites"),
        " across ",
        tags$strong(length(cache$matched_samples), " channels"),
        " ready for PTM DE, Trends, and Single-site analysis.")
    }
  })

  # ============================================================================
  # EXPORT SUMMARIZATION (.zip with full .rds + .preview.rds)
  # Bundle: type="ptm_summarization", route="tmt_ptm"
  # Mirrors server_dia_ptm_upload.R so server_prerun.R can load it identically.
  # ============================================================================
  output$tmt_ptm_upload_export_btn <- downloadHandler(
    filename = function() {
      paste0("TMT_PTM_summarization_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      cache <- rv$dia_ptm_master_cache
      if (is.null(cache)) {
        showNotification("No summarization to export. Complete the wizard first.",
                         type = "warning")
        return()
      }
      ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
      base <- paste0("TMT_PTM_summarization_", ts)

      bundle <- list(
        type            = "ptm_summarization",
        route           = "tmt_ptm",
        dia_ptm_cache   = cache,
        cond_groups     = cache$cond_groups,
        conditions      = sort(unique(cache$design_df$condition)),
        design_option   = "tmt_ptm",
        data_format     = "tmt_ptm",
        annotation_map  = ptm_upl_rv$annotation_map,
        n_plexes        = as.integer(ptm_upl_rv$n_plexes %||% 1L),
        selected_ptms   = ptm_upl_rv$selected_ptms %||% "Phospho",
        subject_pairing = NULL,
        meta            = list(
          software       = "pd_tmt",
          pipeline_label = "TMT PTM \u2014 Proteome Discoverer (limma)",
          n_sites        = cache$n_sites_output,
          n_samples      = length(cache$matched_samples)
        ),
        provenance      = reactiveValuesToList(rv$provenance),
        exported_at     = Sys.time()
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
      logger("[TMT PTM] Summarization exported to ZIP")
    }
  )

  # ============================================================================
  # PRERUN LOAD — restore a previously exported TMT PTM bundle.
  # Triggered by server_prerun.R setting rv$prerun_load_tmt_ptm_bundle.
  # Replicates the same rv assignments as the wizard-done observer so all
  # downstream PTM modules (PCA, DE, Trends, Single) sync correctly.
  # ============================================================================
  observeEvent(rv$prerun_load_tmt_ptm_bundle, {
    bundle <- rv$prerun_load_tmt_ptm_bundle
    req(!is.null(bundle))
    cache <- bundle$dia_ptm_cache
    if (is.null(cache)) {
      showNotification("Bundle is missing TMT PTM cache data.", type = "error")
      rv$prerun_load_tmt_ptm_bundle <- NULL
      return()
    }

    withProgress(message = "Loading TMT PTM experiment\u2026", value = 0.2, {
      incProgress(0.5, detail = "Restoring session\u2026")

      rv$dia_ptm_master_cache     <- cache
      rv$tmt_ptm_wizard_completed <- TRUE
      rv$data_source              <- "tmt_ptm"

      # Restore provenance if present
      if (!is.null(bundle$provenance$dia_ptm))
        rv$provenance$dia_ptm <- bundle$provenance$dia_ptm

      # Sync local wizard reactiveValues so sidebar badges render correctly
      ptm_upl_rv$wizard_completed  <- TRUE
      ptm_upl_rv$summarization_done <- TRUE
      ptm_upl_rv$n_plexes          <- bundle$n_plexes        %||% 1L
      ptm_upl_rv$prob_threshold    <- cache$prob_threshold   %||% 75
      ptm_upl_rv$selected_ptms     <- bundle$selected_ptms   %||% "Phospho"
      ptm_upl_rv$n_sites           <- cache$n_sites_output   %||% 0L
      ptm_upl_rv$n_sites_adjusted  <- cache$n_sites_adjusted %||% 0L
      ptm_upl_rv$conditions        <- sort(unique(cache$design_df$condition))
      ptm_upl_rv$annotation_map    <- bundle$annotation_map

      # Populate rv$ptm_data for cross-module sync
      rv$ptm_data$design_option  <- "tmt_ptm"
      rv$ptm_data$conditions     <- sort(unique(cache$design_df$condition))
      rv$ptm_data$data_format    <- "tmt_ptm"
      rv$ptm_data$n_files        <- (bundle$n_plexes %||% 1L) * 2L
      rv$ptm_data$gene_lookup    <- NULL
      rv$ptm_data$annotation_map <- bundle$annotation_map
      rv$ptm_data$ptm_long       <- NULL
      rv$ptm_data$protein_long   <- NULL

      incProgress(1, detail = "Done")
    })

    # Set source last — triggers all downstream observeEvent(rv$ptm_data$source)
    rv$ptm_data$source            <- NULL
    rv$ptm_data$source            <- "tmt_ptm_upload"
    rv$prerun_load_tmt_ptm_bundle <- NULL   # clear trigger
    logger(paste0("[TMT PTM] Prerun bundle loaded \u2014 sites=", cache$n_sites_output))
  }, ignoreNULL = TRUE)
}

