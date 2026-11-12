# server_dia_upload.R
# ==============================================================================
# DIA Upload Wizard Module
# Extracted from server.R (Phase 1 of OSMSS v4→v5 refactor).
#
# Responsibilities:
#   - DIA data file parsing (Spectronaut pivot report, directly or exported through Perseus)
#   - Condition group assignment
#   - Paired / repeated-measures design
#   - Factorial design (optional)
#   - Summarization via dia_limma_summarize()
#   - Batch correction setup (step 6)
#   - Export / import of summarization bundles (.rds)
#
# Writes to rv$dia_* on wizard completion.
# On wizard completion (dia_wiz_done), additionally sets:
#   rv$data_source  <- "dia"
#   (and clears LF DDA + TMT caches — last-upload-wins)
# ==============================================================================

server_dia_upload <- function(input, output, session, rv, logger) {

  # ── Reactive outputs for conditionalPanel conditions ──
  output$dia_wiz_completed_flag <- reactive({ isTRUE(rv$dia_wizard_completed) })
  outputOptions(output, "dia_wiz_completed_flag", suspendWhenHidden = FALSE)

  output$dia_batch_active <- reactive({
    isTRUE(rv$dia_wizard_completed) && !is.null(rv$dia_batch_groups)
  })
  outputOptions(output, "dia_batch_active", suspendWhenHidden = FALSE)

  # rv$dia_software, from the detected file type: "spectronaut" (Pivot Report),
  # "perseus" (the same report exported through Perseus) or "diann"
  # (report.pg_matrix.tsv).
  .dia_sw_label <- function(sw) switch(sw %||% "spectronaut",
    perseus = "Spectronaut (via Perseus)", diann = "DIA-NN", "Spectronaut")
  .dia_file_label <- function(sw) switch(sw %||% "spectronaut",
    perseus = "Spectronaut Pivot Report (via Perseus)", diann = "DIA-NN protein-group matrix",
    "Spectronaut Pivot Report")

  # DIA-NN's optional report.log.txt or FragPipe's fragpipe.workflow, read and
  # (when it lists the runs) matched to the matrix
  .dia_log_ok <- function() !is.null(rv$dia_diann_log) && is.null(rv$dia_diann_log$error)
  # e.g. "DIA-NN 1.8.1", "DIA-NN 1.9, run within FragPipe", "FragPipe 22.0"
  .dia_log_who <- function(lg)
    paste(c(if (!is.na(lg$version)) paste("DIA-NN", lg$version), if (!is.na(lg$run_within)) lg$run_within),
          collapse = ", run within ")

  # ── Build wizard modal shell ──
  .dia_wiz_build_modal <- function() {
    modalDialog(
      title     = uiOutput("dia_wiz_modal_title"),
      uiOutput("dia_wiz_step_content"),
      footer    = uiOutput("dia_wiz_footer"),
      size      = "l",
      easyClose = FALSE
    )
  }

  observeEvent(input$dia_wiz_open, {
    rv$dia_wizard_step       <- 1L
    rv$dia_wizard_modal_open <- TRUE
    showModal(.dia_wiz_build_modal())
  })

  observeEvent(input$dia_wiz_reopen, {
    rv$dia_wizard_step       <- 1L
    rv$dia_wizard_modal_open <- TRUE
    showModal(.dia_wiz_build_modal())
  })

  # ── Modal title with step progress indicator ──
  output$dia_wiz_modal_title <- renderUI({
    step    <- rv$dia_wizard_step
    n_steps <- 6L
    step_labels <- c("File Upload", "Condition Groups",
                     "Paired Samples", "Factorial Design",
                     "Summary & Summarize", "Batch Correction")
    lbl <- if (step >= 1L && step <= n_steps) step_labels[[step]] else ""
    tagList(
      tags$div(class = "d-flex justify-content-between align-items-center",
        tags$span(icon("atom"), " DIA Upload Wizard"),
        tags$span(class = "text-muted small",
                  paste0("Step ", step, " / ", n_steps, " — ", lbl))
      ),
      tags$div(class = "progress mt-2", style = "height: 5px;",
        tags$div(class = "progress-bar bg-info", role = "progressbar",
                 style = paste0("width:", round(step / n_steps * 100), "%;"))
      )
    )
  })

  # ── Step content router ──
  output$dia_wiz_step_content <- renderUI({
    switch(as.character(rv$dia_wizard_step),
      "1" = .dia_wiz_step1_ui(),
      "2" = .dia_wiz_step2_ui(),
      "3" = .dia_wiz_step3_ui(),
      "4" = .dia_wiz_step4_ui(),
      "5" = .dia_wiz_step5_ui(),
      "6" = .dia_wiz_step6_ui(),
      tags$p("Unknown step.")
    )
  })

  # ── Step UI helpers ──
  # Step 1: the export type is detected from the file (detect_dia_upload_type, global.R).
  .dia_wiz_step1_ui <- function() {
    # Blue question mark with a hover hint (bslib tooltip, as in the label-free DDA upload)
    .hint <- function(text)
      bslib::tooltip(
        tags$span(class = "ms-1", style = "cursor: help; color: #3498db;",
                  icon("circle-question")),
        text, placement = "right")
    tagList(
      tags$h5("Upload DIA Data File"),
      tags$div(class = "alert alert-info",
        icon("circle-info"), tags$strong(" Protein-Level Analysis Only."),
        tags$p(class = "mb-0 mt-1",
          "Due to the size of DIA data files and the structure of DIA quantification, ",
          "this module is limited to pre-summarized protein-level exports.")
      ),
      tags$div(class = "border rounded p-2 mb-3 small",
        tags$div(class = "fw-bold mb-1", "OSMSS supports the following files:"),
        tags$ul(class = "mb-0 ps-3",
          tags$li("Spectronaut Pivot Report (2D).",
            .hint(paste0("Columns PG.ProteinAccessions, PG.Genes and [N] SampleName.PG.Quantity. ",
                         "No imputation in Spectronaut; Q-value Data Filtering applied ",
                         "(Protein Group Q-value < 0.01)."))),
          tags$li("Spectronaut report exported through Perseus.",
            .hint(paste0("Columns T: PG.ProteinGroups, T: PG.Genes and [N] SampleName(.PG.Quantity). ",
                         "log2 or linear values are recognised from the numbers. ",
                         "No imputation in Perseus."))),
          tags$li("DIA-NN protein-group matrix (report.pg_matrix.tsv), also from FragPipe DIA.",
            .hint(paste0("Columns Protein.Group and one column per run; Genes and DIA-NN\u2019s other ",
                         "annotation columns are optional. FragPipe saves it in dia-quant-output. ",
                         "Optional, after the matrix: report.log.txt from the same DIA-NN run, or ",
                         "fragpipe.workflow from the same FragPipe run; either records whether DIA-NN ",
                         "normalized the quantities.")))
        ),
        tags$div(class = "text-muted mt-1",
          "The file type is detected automatically from its column headers.")
      ),
      fileInput("dia_wiz_file", "Select file:",
                accept = c(".tsv", ".txt", ".csv")),
      uiOutput("dia_wiz_qval_ui"),
      uiOutput("dia_wiz_parse_status"),
      uiOutput("dia_wiz_diann_log_ui")
    )
  }

  # Spectronaut's Q-value filtering is the user's to report; DIA-NN filters its matrix itself
  output$dia_wiz_qval_ui <- renderUI({
    if (identical(rv$dia_software, "diann")) return(NULL)
    checkboxInput("dia_wiz_qval_applied",
                  "Q-value filtering (PG Q-value < 0.01) was applied in Spectronaut",
                  value = TRUE)
  })

  # DIA-NN: optional report.log.txt (DIA-NN) or fragpipe.workflow (FragPipe DIA),
  # asked for once the matrix is accepted (its status in its own output, so this
  # file box is not re-rendered by that upload)
  output$dia_wiz_diann_log_ui <- renderUI({
    req(identical(rv$dia_software, "diann"), !is.null(rv$dia_raw_data), is.null(rv$dia_raw_data$error))
    tagList(
      fileInput("dia_wiz_diann_log",
                "Optional: upload report.log.txt (DIA-NN) or fragpipe.workflow (FragPipe) from the same run:",
                accept = c(".txt", ".log", ".workflow"), width = "100%"),
      helpText(class = "small text-muted mt-0",
        "report.log.txt is saved next to the matrix (in FragPipe: dia-quant-output); fragpipe.workflow ",
        "in FragPipe\u2019s output folder. Either records whether DIA-NN normalized the quantities across ",
        "runs, which decides whether OSMSS normalizes them again. Without one, OSMSS assumes DIA-NN\u2019s ",
        "default (normalized) and checks the data in step 5."),
      uiOutput("dia_wiz_diann_log_status")
    )
  })

  output$dia_wiz_diann_log_status <- renderUI({
    lg <- rv$dia_diann_log
    if (is.null(lg))
      return(tags$div(class = "alert alert-secondary py-2 px-3 small mt-2",
        icon("circle-info"), " No log or workflow uploaded (optional)."))
    if (!is.null(lg$error))
      return(tags$div(class = "alert alert-danger py-2 px-3 small mt-2",
        icon("circle-xmark"), " ", lg$error, " You can continue without it."))
    tags$div(class = "alert alert-success py-2 px-3 small mt-2",
      icon("circle-check"), tags$b(paste0(" ", rv$dia_diann_log_filename, " (", .dia_log_who(lg), "): ")),
      if (isTRUE(lg$normalization))
        paste0("DIA-NN cross-run normalization was ON (", lg$norm_mode, "), so OSMSS will not normalize ",
               "the quantities again.")
      else "DIA-NN cross-run normalization was OFF (--no-norm), so OSMSS will apply median centering.",
      if (length(lg$runs_not_in_matrix) > 0L)
        tags$div(class = "text-muted mt-1", paste0(
          length(lg$runs_not_in_matrix), " run(s) in the log are not in the matrix (DIA-NN could not ",
          "load them or found nothing in them): ", paste(lg$runs_not_in_matrix, collapse = ", "), ".")),
      if (length(lg$run_files) == 0L)
        tags$div(class = "text-muted mt-1",
          "This file does not list the runs, so it could not be matched to the matrix\u2019s runs."))
  })

  observeEvent(input$dia_wiz_diann_log, {
    req(input$dia_wiz_diann_log, identical(rv$dia_software, "diann"))
    f     <- input$dia_wiz_diann_log
    first <- tryCatch(readLines(f$datapath, n = 1L, warn = FALSE), error = function(e) "")
    lg    <- if (grepl("^#\\s*(FragPipe|Workflow:)", first[1] %||% "")) diann_read_fp_workflow(f$datapath)
             else diann_read_log(f$datapath)
    # Every run in the matrix must be one the log's DIA-NN run loaded; runs the log
    # names but the matrix lacks are runs DIA-NN skipped (e.g. "cannot load the file")
    if (is.null(lg$error) && length(lg$run_files) > 0L) {
      log_runs <- diann_run_name(lg$run_files)
      not_in_log <- setdiff(rv$dia_raw_data$sample_names, log_runs)
      if (length(not_in_log) > 0L)
        lg <- list(error = paste0("This log is not from the DIA-NN run that wrote the matrix: ",
                                  length(not_in_log), " of the matrix\u2019s runs are not in it (e.g. ",
                                  paste(utils::head(not_in_log, 3L), collapse = ", "), ")."))
      else lg$runs_not_in_matrix <- setdiff(log_runs, rv$dia_raw_data$sample_names)
    }
    rv$dia_diann_log          <- lg
    rv$dia_diann_log_filename <- f$name
    logger(paste0("[DIA] DIA-NN ", lg$source %||% "log", " ", f$name, ": ",
                  if (!is.null(lg$error)) paste0("rejected - ", lg$error)
                  else paste0(.dia_log_who(lg), ", normalization ", lg$norm_mode)))
  })

  .dia_wiz_step2_ui <- function() {
    tagList(
      tags$h5("Assign Condition Groups"),
      tags$p(class = "text-muted small",
        "Assign each sample to a condition group. Minimum 2 groups, ",
        "each with at least 1 sample. At least one group must have >=2 samples."),
      fluidRow(
        column(5,
          textInput("dia_wiz_new_group_name", "New group name:", placeholder = "e.g. Control"),
          checkbox_filter_controls("dia_wiz_sample_search", "dia_wiz_sample_sel"),
          uiOutput("dia_wiz_sample_checklist"),
          actionButton("dia_wiz_add_group_btn", "Add Group",
                       class = "btn-primary btn-sm mt-1", icon = icon("plus"))
        ),
        column(7,
          tags$strong("Defined groups:"),
          uiOutput("dia_wiz_groups_display")
        )
      ),
      uiOutput("dia_wiz_groups_status")
    )
  }

  .dia_wiz_step3_ui <- function() {
    tagList(
      tags$h5("Paired / Repeated Measures Design (Optional)"),
      checkboxInput("dia_wiz_pairing_enabled", "Enable paired/repeated measures design",
                    value = FALSE),
      conditionalPanel(
        condition = "input.dia_wiz_pairing_enabled",
        tags$div(class = "alert alert-info",
          icon("circle-info"),
          " Samples from the same subject must be measured in >=2 conditions. ",
          "limma uses duplicateCorrelation() to account for within-subject correlation."
        ),
        uiOutput("dia_wiz_pairing_inputs"),
        uiOutput("dia_wiz_pairing_validation")
      )
    )
  }

  .dia_wiz_step4_ui <- function() {
    n_groups <- length(rv$dia_cond_groups)
    if (n_groups <= 2L) {
      return(tags$div(class = "alert alert-secondary",
        icon("circle-info"),
        " Factorial design requires more than 2 groups. This step is skipped ",
        "for your current group configuration (", n_groups, " groups)."))
    }
    tagList(
      tags$h5("Factorial Design (Optional)"),
      checkboxInput("dia_wiz_factorial_enabled",
                    "Enable 2x2 factorial design", value = FALSE),
      conditionalPanel(
        condition = "input.dia_wiz_factorial_enabled",
        fluidRow(
          column(6,
            textInput("dia_wiz_f1_name", "Factor 1 name:", placeholder = "e.g. Treatment"),
            textInput("dia_wiz_f1_levels", "Factor 1 levels (comma-sep):",
                      placeholder = "e.g. Vehicle,Drug")
          ),
          column(6,
            textInput("dia_wiz_f2_name", "Factor 2 name:", placeholder = "e.g. Timepoint"),
            textInput("dia_wiz_f2_levels", "Factor 2 levels (comma-sep):",
                      placeholder = "e.g. 4h,24h")
          )
        ),
        tags$strong("Assign factors to conditions:"),
        uiOutput("dia_wiz_factorial_assignments")
      )
    )
  }

  .dia_wiz_step5_ui <- function() {
    sw_lbl  <- .dia_sw_label(rv$dia_software)
    n_samp  <- if (!is.null(rv$dia_protein_matrix)) ncol(rv$dia_protein_matrix) else 0L
    n_prot  <- if (!is.null(rv$dia_protein_matrix)) nrow(rv$dia_protein_matrix) else 0L
    n_grp   <- length(rv$dia_cond_groups)
    paired  <- isTRUE(input$dia_wiz_pairing_enabled)
    tagList(
      tags$h5("Summary & Summarize"),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$tbody(
              tags$tr(tags$td(tags$strong("Software")),
                      tags$td(sw_lbl)),
              tags$tr(tags$td(tags$strong("File")),
                      tags$td(rv$dia_filename %||% "(none)")),
              tags$tr(tags$td(tags$strong("Proteins")),
                      tags$td(paste0(n_prot, " x ", n_samp, " samples"))),
              tags$tr(tags$td(tags$strong("Groups")),
                      tags$td(paste(names(rv$dia_cond_groups), collapse = ", "))),
              tags$tr(tags$td(tags$strong("Pairing")),
                      tags$td(if (paired) "Yes (duplicateCorrelation)" else "No")),
              tags$tr(tags$td(tags$strong("Engine")),
                      tags$td("limma (trend=TRUE, robust=TRUE)")),
              tags$tr(tags$td(tags$strong("Normalization")),
                      tags$td(uiOutput("dia_wiz_norm_row", inline = TRUE)))
            )
          )
        )
      ),
      .dia_wiz_norm_question_ui(),
      actionButton("dia_wiz_summarize_btn", "Summarize",
                   class = "btn-success w-100", icon = icon("play")),
      uiOutput("dia_wiz_summarize_status")
    )
  }

  # ── Step 5: was the input already cross-run normalized? ──
  # Spectronaut and DIA-NN both normalize across runs by default, so the box
  # starts ticked and OSMSS then skips its own median centering. With DIA-NN's
  # report.log.txt the answer is read from the log instead and there is no box.
  # dia_normalization_check() runs on the uploaded matrix and the assigned groups
  # as soon as this step is shown, and its result is displayed to substantiate
  # (or question) the answer. Asked for exports through Perseus too: their
  # quantities are Spectronaut's.
  .dia_wiz_norm_question_ui <- function() {
    is_diann <- identical(rv$dia_software, "diann")
    sw_name  <- if (is_diann) "DIA-NN" else "Spectronaut"
    why <- paste0("Normalizing such data a second time can shift proteins that did not change, ",
                  "especially when a large part of the proteome differs between conditions.")
    if (.dia_log_ok()) {
      lg <- rv$dia_diann_log
      return(tags$div(class = "card mb-3",
        tags$div(class = "card-body",
          tags$strong("Cross-run normalization in DIA-NN"),
          tags$p(class = "small mb-2 mt-1",
            paste0("From ", rv$dia_diann_log_filename, " (", .dia_log_who(lg), "): "),
            if (isTRUE(lg$normalization))
              paste0("normalization was on (", lg$norm_mode, "). ", why,
                     " OSMSS therefore uses the quantities as they are and does not apply its own median centering.")
            else "DIA-NN was run with --no-norm. OSMSS applies column-wise median centering (log2)."),
          uiOutput("dia_wiz_norm_check"))))
    }
    tags$div(class = "card mb-3",
      tags$div(class = "card-body",
        tags$strong(paste("Cross-run normalization in", sw_name)),
        tags$p(class = "small mb-2 mt-1",
          if (is_diann)
            paste0("DIA-NN normalizes protein quantities across runs by default (RT-dependent ",
                   "normalization; switched off only with --no-norm). No report.log.txt was uploaded, ",
                   "so this cannot be read from the run itself. ")
          else
            paste0("Spectronaut normalizes protein quantities across runs by default ",
                   "(Quantification > Cross Run Normalization; strategy 'Automatic', which uses ",
                   "local regression normalization for fewer than 500 runs). "),
          why, " If the input was normalized, OSMSS therefore uses the quantities as they are and ",
          "does not apply its own median centering.",
          if (identical(rv$dia_software, "perseus"))
            " For an export through Perseus this concerns the Spectronaut quantities that went into Perseus."),
        uiOutput("dia_wiz_norm_check"),
        checkboxInput("dia_wiz_input_normalized",
          paste0("The ", sw_name, " quantities were cross-run normalized (", sw_name, " default)"),
          # default Yes; when returning to this step after summarizing, show the answer used
          value = if (is.null(rv$dia_master_cache)) TRUE
                  else !isTRUE(rv$dia_master_cache$normalization_applied %||% TRUE)),
        tags$p(class = "text-muted small mb-0",
          if (is_diann) "Untick only if you know DIA-NN was run with --no-norm. "
          else "Untick only if you know cross-run normalization was switched off in Spectronaut. ",
          "OSMSS then applies column-wise median centering (log2).")
      )
    )
  }

  dia_norm_check <- reactive({
    req(rv$dia_wizard_step == 5L, !is.null(rv$dia_raw_data), is.null(rv$dia_raw_data$error))
    dia_normalization_check(rv$dia_raw_data$protein_matrix, rv$dia_cond_groups)
  })

  output$dia_wiz_norm_check <- renderUI({
    chk <- dia_norm_check()
    cls <- switch(chk$verdict, normalized = "alert-success", borderline = "alert-warning",
                  not_normalized = "alert-danger", "alert-secondary")
    ic  <- switch(chk$verdict, normalized = "circle-check", inconclusive = "circle-info",
                  "triangle-exclamation")
    tags$div(class = paste("alert py-2 px-3 small", cls),
      icon(ic), tags$strong(" Data check: "), chk$message,
      if (chk$verdict %in% c("borderline", "not_normalized") && !.dia_log_ok())
        tags$span(paste0(" If cross-run normalization was switched off in ",
                         if (identical(rv$dia_software, "diann")) "DIA-NN" else "Spectronaut",
                         ", untick the box below.")))
  })

  # Will OSMSS apply its own median centering? DIA-NN with its log: as the log
  # records; otherwise as the user answered in step 5.
  .dia_apply_centering <- function() {
    if (.dia_log_ok()) return(!isTRUE(rv$dia_diann_log$normalization))
    !isTRUE(input$dia_wiz_input_normalized %||% TRUE)
  }

  output$dia_wiz_norm_row <- renderUI({
    if (.dia_apply_centering()) "Column-wise median centering (log2)"
    else paste0("None in OSMSS (input cross-run normalized in ",
                if (identical(rv$dia_software, "diann")) "DIA-NN" else "Spectronaut", ")")
  })

  .dia_wiz_step6_ui <- function() {
    tagList(
      tags$h5("Batch Correction (Optional)"),
      conditionalPanel(
        condition = "!output.dia_wiz_batch_pca_ready",
        tags$div(class = "text-center py-3",
          tags$span(class = "text-muted", icon("circle-notch"), " Loading PCA..."))
      ),
      conditionalPanel(
        condition = "output.dia_wiz_batch_pca_ready",
        plotlyOutput("dia_wiz_batch_pca_plot", height = "300px")
      ),
      checkboxInput("dia_wiz_batch_enabled", "Enable batch correction", value = FALSE),
      conditionalPanel(
        condition = "input.dia_wiz_batch_enabled",
        numericInput("dia_wiz_batch_n", "Number of batches:", value = 2, min = 2, max = 10),
        uiOutput("dia_wiz_batch_assignments"),
        uiOutput("dia_wiz_batch_confound_status")
      ),
      tags$p(class = "text-muted small mt-2",
        "Click 'Done' to complete the wizard. You can return to change settings at any time.")
    )
  }

  # ── Footer ──
  output$dia_wiz_footer <- renderUI({
    step <- rv$dia_wizard_step
    back_btn <- if (step > 1L)
      actionButton("dia_wiz_back", "\u2190 Back", class = "btn-secondary")
    else NULL
    fwd_btn <- if (step < 6L)
      actionButton("dia_wiz_next", "Next \u2192", class = "btn-primary")
    else
      actionButton("dia_wiz_done", "\u2714 Done", class = "btn-success")
    cancel_btn <- modalButton("Cancel")
    tagList(
      tags$div(class = "d-flex justify-content-between",
        tags$div(cancel_btn),
        tags$div(back_btn, fwd_btn)
      )
    )
  })

  # ── Parse status ──
  output$dia_wiz_parse_status <- renderUI({
    raw <- rv$dia_raw_data
    if (is.null(raw)) {
      tags$p(class = "text-muted small", "No file uploaded yet.")
    } else if (!is.null(raw$error)) {
      tags$div(class = "alert alert-danger", icon("circle-xmark"), " ", raw$error)
    } else {
      tags$div(class = "alert alert-success",
        icon("circle-check"), " Detected: ", tags$strong(.dia_file_label(rv$dia_software)),
        ". Parsed: ",
        tags$strong(raw$n_proteins, " proteins"), " x ",
        tags$strong(raw$n_samples, " samples"),
        if (isTRUE(raw$input_log2))
          tags$div(class = "small mt-1",
            "The quantities are log2 (as Perseus exports them); OSMSS does not log-transform them again."),
        if (isTRUE(raw$n_no_gene > 0L))
          tags$div(class = "small mt-1", paste0(
            raw$n_no_gene, " protein groups have no gene name in the matrix and are labelled by their accession.")))
    }
  })

  # ── File input observer: parse on upload ──
  observeEvent(input$dia_wiz_file, {
    req(input$dia_wiz_file)
    f   <- input$dia_wiz_file
    det <- detect_dia_upload_type(f$datapath)
    logger(paste0("[DIA] Upload detected as: ", det$kind))
    sw  <- switch(det$kind, spectronaut_perseus = "perseus", spectronaut_pivot = "spectronaut",
                  diann_pg_matrix = "diann", NULL)
    result <- if (is.null(sw)) {
      list(error = det$message)
    } else {
      tryCatch(
        switch(sw,
          perseus     = parse_maxquant_global(f$datapath, logger = logger),
          diann       = parse_diann_pg_matrix(f$datapath, logger = logger),
          spectronaut = parse_spectronaut_pivot(f$datapath, logger = logger)),
        error = function(e) list(error = conditionMessage(e))
      )
    }
    # New data replaces the whole experiment (all data types, every module's
    # results). The wizard's own progress is kept.
    if (is.null(result$error)) {
      rv$reset_experiment_state(keep = c(RESET_KEEP_ON_UPLOAD,
        "^dia_wizard_step$", "^dia_wizard_modal_open$"))
      logger("[DIA] New file uploaded - previous experiment state reset")
      rv$dia_software <- sw
    }
    rv$dia_raw_data       <- result
    rv$dia_filename       <- f$name
    if (is.null(result$error)) {
      rv$dia_protein_matrix <- result$protein_matrix
      rv$dia_sample_names   <- result$sample_names
      rv$dia_gene_names     <- result$gene_names
      rv$dia_accessions     <- result$accessions
      rv$dia_cond_groups    <- list()
    }
  })

  # ── Step 2: Group management ──
  output$dia_wiz_sample_checklist <- renderUI({
    all_samples <- rv$dia_sample_names
    if (is.null(all_samples)) return(tags$p(class = "text-muted small", "Upload a file first."))
    assigned <- unlist(rv$dia_cond_groups)
    available <- setdiff(all_samples, assigned)
    if (length(available) == 0) return(tags$p(class = "text-muted small", "All samples assigned."))
    checkboxGroupInput("dia_wiz_sample_sel", "Available samples:",
                       choices = available, selected = NULL)
  })

  output$dia_wiz_groups_display <- renderUI({
    grps <- rv$dia_cond_groups
    if (length(grps) == 0) return(tags$p(class = "text-muted small", "No groups defined."))
    tagList(lapply(names(grps), function(gname) {
      tags$div(class = "card mb-1",
        tags$div(class = "card-body py-1 px-2 d-flex justify-content-between align-items-center",
          tags$span(tags$strong(gname), ": ",
                    paste(grps[[gname]], collapse = ", ")),
          actionButton(paste0("dia_wiz_rm_grp_", make.names(gname)),
                       label = "", icon = icon("trash"),
                       class = "btn-danger btn-sm",
                       `data-grpname` = gname)
        )
      )
    }))
  })

  observeEvent(input$dia_wiz_add_group_btn, {
    gname  <- trimws(input$dia_wiz_new_group_name %||% "")
    samps  <- input$dia_wiz_sample_sel
    if (nchar(gname) == 0 || length(samps) == 0) return()
    grps           <- rv$dia_cond_groups
    grps[[gname]]  <- samps
    rv$dia_cond_groups <- grps
    updateTextInput(session, "dia_wiz_new_group_name", value = "")
    updateTextInput(session, "dia_wiz_sample_search",  value = "")
  })

  # Dynamic remove-group buttons
  observe({
    grps <- rv$dia_cond_groups
    lapply(names(grps), function(gname) {
      btn_id <- paste0("dia_wiz_rm_grp_", make.names(gname))
      observeEvent(input[[btn_id]], {
        grps_now <- rv$dia_cond_groups
        grps_now[[gname]] <- NULL
        rv$dia_cond_groups <- grps_now
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  output$dia_wiz_groups_status <- renderUI({
    grps    <- rv$dia_cond_groups
    n_grps  <- length(grps)
    if (n_grps < 2L)
      return(tags$div(class = "alert alert-warning py-2 mt-2",
                      icon("triangle-exclamation"), " Need at least 2 groups."))
    grp_sizes <- vapply(grps, length, integer(1))
    if (!any(grp_sizes >= 2L))
      return(tags$div(class = "alert alert-warning py-2 mt-2",
                      icon("triangle-exclamation"),
                      " At least one group must have >=2 samples."))
    tags$div(class = "alert alert-success py-2 mt-2",
             icon("circle-check"), " ", n_grps, " groups defined. Ready to proceed.")
  })

  # ── Step 3: Pairing inputs ──
  output$dia_wiz_pairing_inputs <- renderUI({
    req(isTRUE(input$dia_wiz_pairing_enabled))
    samples <- unlist(rv$dia_cond_groups)
    if (length(samples) == 0) return(NULL)
    tagList(
      tags$p(class = "text-muted small",
        "Assign a Subject ID to each sample. Samples from the same subject share the same ID."),
      lapply(samples, function(sn) {
        fluidRow(
          column(6, tags$label(class = "col-form-label col-form-label-sm", sn)),
          column(6, selectizeInput(
            paste0("dia_wiz_pair_", make.names(sn)),
            label = NULL,
            choices = character(0),
            options = list(create = TRUE, placeholder = "Subject ID")
          ))
        )
      })
    )
  })

  output$dia_wiz_pairing_validation <- renderUI({
    req(isTRUE(input$dia_wiz_pairing_enabled))
    samples <- unlist(rv$dia_cond_groups)
    if (length(samples) == 0) return(NULL)
    grp_map <- setNames(
      rep(names(rv$dia_cond_groups),
          vapply(rv$dia_cond_groups, length, integer(1))),
      unlist(rv$dia_cond_groups)
    )
    subject_ids <- vapply(samples, function(sn) {
      input[[paste0("dia_wiz_pair_", make.names(sn))]] %||% ""
    }, character(1))
    bad_subjects <- character(0)
    for (subj in unique(subject_ids[subject_ids != ""])) {
      subj_samps <- samples[subject_ids == subj]
      n_conds    <- length(unique(grp_map[subj_samps]))
      if (n_conds < 2L) bad_subjects <- c(bad_subjects, subj)
    }
    if (length(bad_subjects) > 0) {
      tags$div(class = "alert alert-warning mt-2",
        icon("triangle-exclamation"), " Subject(s) ",
        paste(bad_subjects, collapse = ", "),
        " appear in fewer than 2 conditions.")
    } else if (all(subject_ids == "")) {
      tags$div(class = "alert alert-warning mt-2",
        icon("triangle-exclamation"), " No subject IDs assigned.")
    } else {
      tags$div(class = "alert alert-success mt-2",
               icon("circle-check"), " Pairing looks valid.")
    }
  })

  # ── Step 4: Factorial assignment dropdowns ──
  output$dia_wiz_factorial_assignments <- renderUI({
    req(isTRUE(input$dia_wiz_factorial_enabled))
    f1_levels <- Filter(function(x) nchar(x) > 0, trimws(strsplit(input$dia_wiz_f1_levels %||% "", ",")[[1]]))
    f2_levels <- Filter(function(x) nchar(x) > 0, trimws(strsplit(input$dia_wiz_f2_levels %||% "", ",")[[1]]))
    if (length(f1_levels) == 0 || length(f2_levels) == 0) return(NULL)
    grp_names <- names(rv$dia_cond_groups)
    tagList(lapply(grp_names, function(gn) {
      fluidRow(
        column(4, tags$label(class = "col-form-label col-form-label-sm mt-2",
                             tags$strong(gn))),
        column(4, selectInput(paste0("dia_wiz_fac1_", make.names(gn)),
                              paste0(input$dia_wiz_f1_name %||% "Factor 1", ":"),
                              choices = f1_levels)),
        column(4, selectInput(paste0("dia_wiz_fac2_", make.names(gn)),
                              paste0(input$dia_wiz_f2_name %||% "Factor 2", ":"),
                              choices = f2_levels))
      )
    }))
  })

  # ── Step 5: Summarize button ──
  output$dia_wiz_summarize_status <- renderUI({
    if (is.null(rv$dia_master_cache)) {
      tags$p(class = "text-muted small mt-2",
             "Click 'Summarize' to process the data.")
    } else {
      m <- rv$dia_master_cache
      tags$div(class = "alert alert-success mt-2",
        icon("circle-check"), " Summarization complete: ",
        tags$strong(m$n_prots_processed, " proteins"),
        " (from ", m$n_prots_input, " input). Normalization: ",
        if (isFALSE(m$normalization_applied)) "none in OSMSS (input cross-run normalized)."
        else "column-wise median centering."
      )
    }
  })

  observeEvent(input$dia_wiz_summarize_btn, {
    raw <- rv$dia_raw_data
    if (is.null(raw) || !is.null(raw$error)) {
      showNotification("Please upload a valid file first.", type = "warning")
      return()
    }
    if (length(rv$dia_cond_groups) < 2L) {
      showNotification("Define at least 2 groups first.", type = "warning")
      return()
    }
    subj_pairing <- NULL
    if (isTRUE(input$dia_wiz_pairing_enabled)) {
      samples <- unlist(rv$dia_cond_groups)
      subj_pairing <- setNames(
        vapply(samples, function(sn)
          input[[paste0("dia_wiz_pair_", make.names(sn))]] %||% sn,
          character(1)),
        samples
      )
      rv$dia_subject_pairing <- subj_pairing
    } else {
      rv$dia_subject_pairing <- NULL
    }
    fac_design <- NULL
    if (isTRUE(input$dia_wiz_factorial_enabled)) {
      gnames    <- names(rv$dia_cond_groups)
      f1_assign <- setNames(
        vapply(gnames, function(gn) input[[paste0("dia_wiz_fac1_", make.names(gn))]] %||% "", character(1)),
        gnames)
      f2_assign <- setNames(
        vapply(gnames, function(gn) input[[paste0("dia_wiz_fac2_", make.names(gn))]] %||% "", character(1)),
        gnames)
      fac_design <- list(
        enabled     = TRUE,
        f1_name     = input$dia_wiz_f1_name %||% "Factor1",
        f1_levels   = strsplit(input$dia_wiz_f1_levels %||% "", ",")[[1]] |> trimws(),
        f2_name     = input$dia_wiz_f2_name %||% "Factor2",
        f2_levels   = strsplit(input$dia_wiz_f2_levels %||% "", ",")[[1]] |> trimws(),
        assignments = list(
          condition = gnames,
          f1        = unname(f1_assign),
          f2        = unname(f2_assign)
        )
      )
      rv$dia_factorial_design <- fac_design
    } else {
      rv$dia_factorial_design <- NULL
    }

    apply_centering <- .dia_apply_centering()
    norm_chk <- tryCatch(dia_normalization_check(raw$protein_matrix, rv$dia_cond_groups),
                         error = function(e) list(verdict = "inconclusive", message = conditionMessage(e)))
    withProgress(message = "Summarizing DIA data...", value = 0.3, {
      result <- tryCatch(
        dia_limma_summarize(
          parsed_data      = raw,
          cond_groups      = rv$dia_cond_groups,
          subject_pairing  = subj_pairing,
          factorial_design = fac_design,
          apply_centering  = apply_centering,
          logger           = logger
        ),
        error = function(e) list(error = conditionMessage(e))
      )
      incProgress(0.7)
    })

    if (!is.null(result$error)) {
      showNotification(paste0("Summarization failed: ", result$error), type = "error")
      return()
    }

    rv$dia_master_cache <- result
    rv$dia_master_cache_meta <- list(
      pipeline      = "dia",
      software      = rv$dia_software %||% "spectronaut",
      engine        = "limma",
      n_proteins    = result$n_prots_processed,
      has_pairing   = !is.null(subj_pairing),
      has_factorial = !is.null(fac_design),
      timestamp     = Sys.time()
    )
    rv$dia_cache_stale <- FALSE

    via_perseus <- identical(rv$dia_software, "perseus")
    is_diann    <- identical(rv$dia_software, "diann")
    log_ok      <- .dia_log_ok()
    rv$provenance$dia_global <- list(
      input_file           = rv$dia_filename,
      software             = if (is_diann) "DIA-NN" else "Spectronaut",
      software_version     = if (log_ok) rv$dia_diann_log$version else NA_character_,
      run_within           = if (log_ok) rv$dia_diann_log$run_within else NA_character_,
      exported_via         = if (via_perseus) "Perseus" else NA_character_,
      input_log2           = isTRUE(raw$input_log2),
      data_level           = "protein",
      # DIA-NN filters report.pg_matrix.tsv itself (1% protein-group FDR)
      qvalue_filtered      = if (is_diann) NA else isTRUE(input$dia_wiz_qval_applied),
      normalization_source = if (log_ok) rv$dia_diann_log$source else "user",
      upstream_norm_mode   = if (log_ok) rv$dia_diann_log$norm_mode else NA_character_,
      n_proteins_input     = result$n_prots_input,
      n_dup_removed        = rv$dia_raw_data$n_dup_removed %||% 0L,
      n_proteins_processed = result$n_prots_processed,
      n_samples            = length(rv$dia_sample_names),
      n_conditions         = length(rv$dia_cond_groups),
      condition_names      = names(rv$dia_cond_groups),
      normalization        = if (apply_centering) "column-wise median centering (log2)"
                             else "none in OSMSS (input cross-run normalized upstream)",
      normalization_applied        = apply_centering,
      input_normalized_upstream    = !apply_centering,
      norm_check_verdict           = norm_chk$verdict %||% NA_character_,
      norm_check_within_sd         = norm_chk$within_sd %||% NA_real_,
      norm_check_n_complete        = norm_chk$n_complete %||% NA_integer_,
      normalization_max_shift_log2 = result$normalization_max_shift_log2 %||% NA_real_,
      normalization_basis      = result$normalization_basis      %||% NA_character_,
      normalization_used_cc    = result$normalization_used_cc    %||% NA,
      normalization_n_complete = result$normalization_n_complete %||% NA_integer_,
      normalization_delta_log2 = result$normalization_delta_log2 %||% NA_real_,
      engine               = "limma (eBayes: trend=TRUE, robust=TRUE)",
      pairing              = !is.null(subj_pairing),
      factorial            = !is.null(fac_design),
      timestamp            = as.character(Sys.time())
    )
    logger("[DIA] Wizard step 5: summarization complete.")
  })

  # ── Step 6: Batch PCA ──
  output$dia_wiz_batch_pca_ready <- reactive({
    isTRUE(rv$dia_wizard_step == 6L) && !is.null(rv$dia_master_cache)
  })
  outputOptions(output, "dia_wiz_batch_pca_ready", suspendWhenHidden = FALSE)

  output$dia_wiz_batch_pca_plot <- renderPlotly({
    req(rv$dia_wizard_step == 6L)
    req(!is.null(rv$dia_master_cache))
    mat    <- rv$dia_master_cache$protein_matrix
    grps   <- rv$dia_cond_groups
    cond_map <- setNames(
      rep(names(grps), vapply(grps, length, integer(1))),
      unlist(grps)
    )
    mat <- apply(mat, 1, function(r) { r[is.na(r)] <- median(r, na.rm = TRUE); r }) |> t()
    mat <- mat[complete.cases(mat), , drop = FALSE]
    if (nrow(mat) < 3L) return(plotly::plotly_empty())
    pca_res  <- prcomp(t(mat), center = TRUE, scale. = FALSE)
    pcs      <- as.data.frame(pca_res$x[, 1:min(2, ncol(pca_res$x))])
    pcs$Sample <- rownames(pcs)
    pcs$Group  <- cond_map[pcs$Sample]
    var_exp    <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
    plotly::plot_ly(pcs, x = ~PC1, y = ~PC2, color = ~Group,
                    text = ~Sample, type = "scatter", mode = "markers+text",
                    textposition = "top center") |>
      plotly::layout(
        xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
        yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
        title = "PCA Preview (pre-batch correction)",
        showlegend = TRUE
      )
  })

  output$dia_wiz_batch_assignments <- renderUI({
    req(isTRUE(input$dia_wiz_batch_enabled))
    n_batches <- as.integer(input$dia_wiz_batch_n %||% 2L)
    bat_choices <- paste0("Batch ", seq_len(n_batches))
    samples     <- unlist(rv$dia_cond_groups)
    tagList(
      tags$p(class = "text-muted small mt-2", "Assign each sample to a batch:"),
      lapply(samples, function(sn) {
        fluidRow(
          column(6, tags$label(class = "col-form-label col-form-label-sm", sn)),
          column(6, selectInput(paste0("dia_wiz_batch_sel_", make.names(sn)),
                                label = NULL, choices = bat_choices,
                                selected = bat_choices[1]))
        )
      })
    )
  })

  output$dia_wiz_batch_confound_status <- renderUI({
    req(isTRUE(input$dia_wiz_batch_enabled))
    n_batches   <- as.integer(input$dia_wiz_batch_n %||% 2L)
    bat_choices <- paste0("Batch ", seq_len(n_batches))
    samples     <- unlist(rv$dia_cond_groups)
    grp_map     <- setNames(
      rep(names(rv$dia_cond_groups),
          vapply(rv$dia_cond_groups, length, integer(1))),
      unlist(rv$dia_cond_groups)
    )
    bat_assign <- vapply(samples, function(sn)
      input[[paste0("dia_wiz_batch_sel_", make.names(sn))]] %||% bat_choices[1],
      character(1))
    cond_vec   <- grp_map[samples]
    tbl        <- table(bat_assign, cond_vec)
    confounded <- any(colSums(tbl > 0L) < 2L)
    if (confounded) {
      tags$div(class = "alert alert-warning py-2 px-3 small",
        icon("triangle-exclamation"),
        tags$strong(" Confounding detected."),
        " Some conditions appear in only one batch. Batch will be used for ",
        tags$em("visualization only"), " (not added to the limma model).")
    } else {
      tags$div(class = "alert alert-success py-2 px-3 small",
        icon("circle-check"),
        " No confounding detected. Batch will be added to the limma model.")
    }
  })

  # ── Navigation: Next ──
  observeEvent(input$dia_wiz_next, {
    step <- rv$dia_wizard_step
    if (step == 1L) {
      if (is.null(rv$dia_raw_data) || !is.null(rv$dia_raw_data$error)) {
        showNotification("Please upload a valid file.", type = "warning")
        return()
      }
      rv$dia_wizard_step <- 2L
    } else if (step == 2L) {
      grps <- rv$dia_cond_groups
      if (length(grps) < 2L) {
        showNotification("Define at least 2 condition groups.", type = "warning")
        return()
      }
      if (!any(vapply(grps, length, integer(1)) >= 2L)) {
        showNotification("At least one group must have >=2 samples.", type = "warning")
        return()
      }
      rv$dia_wizard_step <- 3L
    } else if (step == 3L) {
      rv$dia_wizard_step <- 4L
    } else if (step == 4L) {
      rv$dia_wizard_step <- 5L
    } else if (step == 5L) {
      if (is.null(rv$dia_master_cache)) {
        showNotification("Please run Summarize first.", type = "warning")
        return()
      }
      # The normalization answer was changed after summarizing: the cache no longer matches it
      if (!identical(isTRUE(rv$dia_master_cache$normalization_applied %||% TRUE), .dia_apply_centering())) {
        showNotification("The normalization setting changed after summarizing. Please click Summarize again.",
                         type = "warning")
        return()
      }
      rv$dia_wizard_step <- 6L
    }
    showModal(.dia_wiz_build_modal())
  })

  # ── Navigation: Back ──
  observeEvent(input$dia_wiz_back, {
    step <- rv$dia_wizard_step
    rv$dia_wizard_step <- max(1L, step - 1L)
    showModal(.dia_wiz_build_modal())
  })

  # ── Done — set rv$data_source, clear other caches ──
  observeEvent(input$dia_wiz_done, {
    if (isTRUE(input$dia_wiz_batch_enabled)) {
      n_batches   <- as.integer(input$dia_wiz_batch_n %||% 2L)
      bat_choices <- paste0("Batch ", seq_len(n_batches))
      samples     <- unlist(rv$dia_cond_groups)
      bat_assign  <- vapply(samples, function(sn)
        input[[paste0("dia_wiz_batch_sel_", make.names(sn))]] %||% bat_choices[1],
        character(1))
      bat_groups <- lapply(bat_choices, function(b)
        names(bat_assign)[bat_assign == b])
      names(bat_groups) <- bat_choices
      rv$dia_batch_groups <- bat_groups

      grp_map <- setNames(
        rep(names(rv$dia_cond_groups),
            vapply(rv$dia_cond_groups, length, integer(1))),
        unlist(rv$dia_cond_groups)
      )
      cond_vec   <- grp_map[samples]
      tbl        <- table(bat_assign, cond_vec)
      confounded <- any(colSums(tbl > 0L) < 2L)
      rv$dia_batch_corrected  <- TRUE
      rv$dia_batch_for_model  <- !confounded
    } else {
      rv$dia_batch_groups    <- NULL
      rv$dia_batch_corrected <- FALSE
      rv$dia_batch_for_model <- FALSE
    }

    rv$dia_wizard_completed  <- TRUE
    rv$dia_wizard_modal_open <- FALSE
    removeModal()
    logger("[DIA] Wizard completed.")

    # ── Set active data source ──
    # The previous experiment (LF, TMT, older DIA) was already cleared when
    # this file was uploaded (dia_wiz_file -> reset_experiment_state).
    rv$data_source <- "dia"

    logger("[DIA] rv$data_source set to 'dia'.")
  })

  # ── Post-wizard sidebar badges ──
  output$dia_wiz_sidebar_badges <- renderUI({
    req(isTRUE(rv$dia_wizard_completed))
    m <- rv$dia_master_cache_meta
    tagList(
      tags$div(class = "alert alert-success py-2 px-3 small",
        icon("circle-check"), " Wizard complete.",
        tags$br(),
        if (!is.null(m)) {
          tagList(
            tags$span(class = "badge bg-primary me-1", m$n_proteins, " proteins"),
            tags$span(class = "badge bg-secondary me-1",
                      length(rv$dia_cond_groups), " groups"),
            if (isTRUE(m$has_pairing))
              tags$span(class = "badge bg-info me-1", "Paired")
          )
        }
      )
    )
  })

  # ── Upload main status ──
  output$dia_upload_main_status <- renderUI({
    if (!isTRUE(rv$dia_wizard_completed)) {
      return(tags$div(class = "alert alert-secondary",
        icon("circle-info"), " Complete the Upload Wizard to get started."))
    }
    m <- rv$dia_master_cache_meta
    tags$div(class = "alert alert-success",
      icon("circle-check"), " Data ready.",
      tags$br(),
      paste0(m$n_proteins, " proteins | ",
             length(rv$dia_cond_groups), " groups | ",
             if (isTRUE(m$has_pairing)) "Paired" else "Unpaired",
             " | ", .dia_sw_label(m$software)))
  })

  # ── Prerun Experiments: DIA bundle load trigger ──────────────────────────────
  observeEvent(rv$prerun_load_dia_bundle, {
    bundle <- rv$prerun_load_dia_bundle
    req(!is.null(bundle))
    withProgress(message = "Loading DIA experiment...", value = 0.3, {
      rv$dia_master_cache      <- bundle$cache
      rv$dia_master_cache_meta <- bundle$meta
      rv$dia_cond_groups       <- bundle$cond_groups
      rv$dia_subject_pairing   <- bundle$subject_pairing
      rv$dia_factorial_design  <- bundle$factorial_design
      rv$dia_batch_groups      <- bundle$batch_groups
      rv$dia_batch_corrected   <- bundle$batch_corrected
      rv$dia_batch_for_model   <- bundle$batch_for_model
      rv$dia_sample_names      <- colnames(bundle$cache$protein_matrix)
      rv$dia_gene_names        <- bundle$cache$acc_to_gene
      rv$dia_accessions        <- rownames(bundle$cache$protein_matrix)
      rv$dia_filename          <- bundle$meta$input_file %||% "(imported)"
      rv$dia_wizard_completed  <- TRUE
      rv$data_source           <- "dia"
      # Restore summarization provenance so the methods generator can report
      # normalization and pre-filter counts for imported experiments.
      rv$provenance$dia_global <- bundle$provenance_dia_global %||% list(
        input_file           = bundle$meta$input_file %||% "(imported)",
        software             = "Spectronaut",
        data_level           = "protein",
        n_proteins_input     = bundle$cache$n_prots_input,
        n_proteins_processed = bundle$cache$n_prots_processed,
        n_samples            = ncol(bundle$cache$protein_matrix),
        n_conditions         = length(bundle$cond_groups),
        condition_names      = names(bundle$cond_groups),
        normalization        = "column-wise median centering (log2)",
        engine               = "limma (eBayes: trend=TRUE, robust=TRUE)",
        timestamp            = as.character(Sys.time())
      )
      incProgress(1)
    })
    rv$prerun_load_dia_bundle <- NULL
    showNotification("DIA experiment loaded.", type = "message", duration = 4)
    logger("[DIA] Bundle loaded from Prerun Experiments")
  }, ignoreNULL = TRUE)

  # ── RDS Export ──
  output$dia_upload_export_summ_btn <- downloadHandler(
    filename = function() {
      paste0("DIA_summarization_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
      base <- paste0("DIA_summarization_", ts)
      bundle <- list(
        type             = "dia_summarization",
        version          = 1L,
        cache            = rv$dia_master_cache,
        meta             = rv$dia_master_cache_meta,
        cond_groups      = rv$dia_cond_groups,
        subject_pairing  = rv$dia_subject_pairing,
        factorial_design = rv$dia_factorial_design,
        batch_groups     = rv$dia_batch_groups,
        batch_corrected  = rv$dia_batch_corrected,
        batch_for_model  = rv$dia_batch_for_model,
        provenance_dia_global = rv$provenance$dia_global
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
      logger("DIA: Summarization exported")
    }
  )

}
