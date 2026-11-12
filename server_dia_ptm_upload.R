# server_dia_ptm_upload.R
# ==============================================================================
# Label Free PTM Upload Wizard Module
# Handles label-free phosphoproteomics data upload for PTM analysis.
# Supports DIA input from the Spectronaut precursor-level report only. The DDA (MaxQuant / Perseus) and Spectronaut
# PTM pivot options were removed from the wizard 2026-10-10 (untested; to be rebuilt) - see
# Backup_ptmddaremoved_20261010_081944 and Backup_ptmpivotremoved_20261010_083338.
#
# Step 1 also picks the analysis mode: "PTM with global proteomics" (default) or
# "PTM only" (no global file; absolute PTM only - see ptm_layers() in global.R).
#
# Responsibilities:
#   - Phospho file parsing (Spectronaut DIA or MaxQuant DDA)
#   - Global proteomics file parsing (Spectronaut pivot or MaxQuant matrix; not in PTM-only mode)
#   - Condition group assignment
#   - Paired / repeated-measures design
#   - Probability threshold selection
#   - Summarization + stoichiometric adjustment (DDA: MinProb imputation first)
#   - Export of rv$dia_ptm_master_cache
#   - Population of rv$ptm_data for cross-module sync
#
# On wizard completion:
#   rv$data_source          <- "dia_ptm"
#   rv$dia_ptm_master_cache <- <summarization cache list>
#   rv$ptm_data$source      <- "dia_ptm_upload"  (triggers sync in PTM modules)
#   rv$ptm_data$design_option <- "dia_ptm"
#   rv$provenance$dia_ptm   <- <provenance list>
# ==============================================================================

server_dia_ptm_upload <- function(input, output, session, rv, logger) {

  N_STEPS <- 7L
  STEP_LABELS <- c(
    "Input Type",
    "File Upload",
    "Condition Groups",
    "Paired Samples",
    "Probability Threshold",
    "Summary & Summarize",
    "Review"
  )

  # Analysis mode chosen on step 1: "with_global" (default) or "ptm_only"
  .wiz_ptm_only <- function() identical(rv$dia_ptm_mode, "ptm_only")

  # ── Reactive outputs for conditionalPanel conditions ──
  output$dia_ptm_wiz_completed_flag <- reactive({ isTRUE(rv$dia_ptm_wizard_completed) })
  outputOptions(output, "dia_ptm_wiz_completed_flag", suspendWhenHidden = FALSE)

  # ── Main content panel status (shown in the tabPanel body) ──
  output$dia_ptm_upload_main_status <- renderUI({
    if (!isTRUE(rv$dia_ptm_wizard_completed)) {
      tags$div(class = "alert alert-info py-2 px-3 mt-3",
        icon("wand-magic-sparkles"),
        " Open the DIA PTM Upload Wizard in the sidebar to get started.")
    } else {
      cache <- rv$dia_ptm_master_cache
      tags$div(class = "alert alert-success py-2 px-3 mt-3",
        icon("circle-check"), " Upload complete. ",
        tags$strong(cache$n_sites_output, " phosphosites"),
        " across ",
        tags$strong(length(cache$matched_samples), " samples"),
        " are ready for analysis in the PTM modules.")
    }
  })

  # ── Build wizard modal shell ──
  .wiz_build_modal <- function() {
    modalDialog(
      title     = uiOutput("dia_ptm_wiz_modal_title"),
      uiOutput("dia_ptm_wiz_step_content"),
      footer    = uiOutput("dia_ptm_wiz_footer"),
      size      = "l",
      easyClose = FALSE
    )
  }

  observeEvent(input$dia_ptm_wiz_open, {
    rv$dia_ptm_wizard_step       <- 1L
    rv$dia_ptm_wizard_modal_open <- TRUE
    showModal(.wiz_build_modal())
  })

  observeEvent(input$dia_ptm_wiz_reopen, {
    rv$dia_ptm_wizard_step       <- 1L
    rv$dia_ptm_wizard_modal_open <- TRUE
    showModal(.wiz_build_modal())
  })

  # ── Modal title with step progress indicator ──
  output$dia_ptm_wiz_modal_title <- renderUI({
    step <- rv$dia_ptm_wizard_step %||% 1L
    lbl  <- if (step >= 1L && step <= N_STEPS) STEP_LABELS[[step]] else ""
    tagList(
      tags$div(class = "d-flex justify-content-between align-items-center",
        tags$span(icon("flask"), " DIA PTM Upload Wizard"),
        tags$span(class = "text-muted small",
                  paste0("Step ", step, " / ", N_STEPS, " \u2014 ", lbl))
      ),
      tags$div(class = "progress mt-2", style = "height: 5px;",
        tags$div(class = "progress-bar bg-success", role = "progressbar",
                 style = paste0("width:", round(step / N_STEPS * 100), "%;"))
      )
    )
  })

  # ── Step content router ──
  output$dia_ptm_wiz_step_content <- renderUI({
    switch(as.character(rv$dia_ptm_wizard_step %||% 1L),
      "1" = .wiz_step1_ui(),
      "2" = .wiz_step2_ui(),
      "3" = .wiz_step3_ui(),
      "4" = .wiz_step4_ui(),
      "5" = .wiz_step5_ui(),
      "6" = .wiz_step6_ui(),
      "7" = .wiz_step7_ui(),
      tags$p("Unknown step.")
    )
  })

  # ---------------------------------------------------------------------------
  # STEP 1 — Input Type
  # ---------------------------------------------------------------------------
  .wiz_step1_ui <- function() {
    current_sw <- rv$dia_ptm_software %||% "spectronaut"
    mode       <- rv$dia_ptm_mode %||% "with_global"
    tagList(
      # ── Analysis mode: with or without a global proteome ──
      # One pipeline either way: PTM-only runs the same summarization and DE as the
      # with-global mode and leaves out the protein summary and the stoichiometric
      # adjustment (ptm_layers() in global.R).
      tags$h5("Select Analysis Mode"),
      fluidRow(
        column(6,
          actionButton("dia_ptm_wiz_mode_global", "PTM with global proteomics",
                       class = if (mode == "with_global") "btn-success w-100" else "btn-outline-secondary w-100",
                       icon = icon("layer-group")),
          tags$p(class = "text-muted small mt-1",
            "Phospho file + global protein file from the same samples. Results: adjusted to total ",
            "protein (stoichiometry), absolute PTM and protein.")
        ),
        column(6,
          actionButton("dia_ptm_wiz_mode_ptmonly", "PTM only",
                       class = if (mode == "ptm_only") "btn-success w-100" else "btn-outline-secondary w-100",
                       icon = icon("flask")),
          tags$p(class = "text-muted small mt-1",
            "Phospho file alone. Results: absolute PTM only. A change in a site then mixes ",
            "phosphorylation change with change in its protein's abundance.")
        )
      ),
      tags$hr(),
      tags$h5("Select Input Type"),
      tags$p(class = "text-muted small",
        "Choose the acquisition mode and software that generated your files. ",
        "This determines both parsing and the analysis pipeline."),
      fluidRow(
        # ── Option 1: Spectronaut DIA precursor-level (DEFAULT / RECOMMENDED) ──
        # Only the default Spectronaut precursor-level report is offered (Erez, 2026-10-10). Removed from this step:
        #  - "DDA - MaxQuant / Perseus": never tested or troubleshot; read Spectronaut-style PTM.CollapseKey columns,
        #    not MaxQuant's own site table (Backup_ptmddaremoved_20261010_081944);
        #  - "DIA - Spectronaut PTM Pivot": not tested either (Backup_ptmpivotremoved_20261010_083338).
        # Both are to be rebuilt later. Their "maxquant" / "spectronaut_pivot" branches below are unreachable from the
        # UI and are kept for that rebuild; parse_spectronaut_pivot() is still used to read the GLOBAL protein file of
        # the precursor route and must stay.
        column(12,
          tags$div(
            class = if (current_sw == "spectronaut") "card border-success mb-2" else "card mb-2",
            style = "cursor:pointer;",
            id    = "dia_ptm_sw_card_spectronaut",
            tags$div(class = "card-body py-3",
              tags$div(class = "d-flex align-items-center",
                tags$div(class = "me-3", style = "font-size:2rem;", icon("flask")),
                tags$div(
                  tags$strong("DIA \u2014 Spectronaut"),
                  tags$span(class = "badge bg-success ms-1 small", "Recommended"),
                  tags$p(class = "text-muted small mb-0 mt-1",
                    "Precursor-level phospho report. ",
                    "Highest quantitative resolution: TMP normalization, ",
                    "per-precursor localization probabilities, kNN imputation.")
                )
              ),
              if (current_sw == "spectronaut")
                tags$div(class = "mt-2", tags$span(class = "badge bg-success", icon("check"), " Selected"))
            )
          ),
          actionButton("dia_ptm_wiz_sw_spectronaut", "Select DIA Spectronaut",
                       class = if (current_sw == "spectronaut") "btn-success w-100" else "btn-outline-secondary w-100",
                       icon = icon("flask"))
        )
      ),
      tags$div(class = "alert alert-info mt-2",
        icon("circle-info"),
        tags$strong(" DIA Spectronaut (precursor): "), "TMP normalization + kNN imputation + limma eBayes. Highest quality."
      )
    )
  }

  observeEvent(input$dia_ptm_wiz_sw_spectronaut, {
    rv$dia_ptm_software        <- "spectronaut"
    rv$dia_ptm_acquisition     <- "dia"
    rv$dia_ptm_phospho_raw     <- NULL
    rv$dia_ptm_global_raw      <- NULL
    rv$dia_ptm_master_cache    <- NULL
    showModal(.wiz_build_modal())
  })

  # Switching mode drops the global file and any summarization built in the other mode
  .wiz_set_mode <- function(mode) {
    if (!identical(rv$dia_ptm_mode %||% "with_global", mode)) {
      rv$dia_ptm_mode            <- mode
      rv$dia_ptm_global_raw      <- NULL
      rv$dia_ptm_global_filename <- NULL
      rv$dia_ptm_master_cache    <- NULL
    }
    showModal(.wiz_build_modal())
  }
  observeEvent(input$dia_ptm_wiz_mode_global,  .wiz_set_mode("with_global"))
  observeEvent(input$dia_ptm_wiz_mode_ptmonly, .wiz_set_mode("ptm_only"))

  # ---------------------------------------------------------------------------
  # STEP 2 — File Upload
  # ---------------------------------------------------------------------------
  .wiz_step2_ui <- function() {
    sw <- rv$dia_ptm_software %||% "spectronaut"
    if (sw == "spectronaut") {
      tagList(
        tags$h5("Upload DIA Phospho Files (Spectronaut \u2014 Precursor-Level)"),
        tags$div(class = "alert alert-secondary mb-3",
          tags$strong("Spectronaut Phospho Report requirements:"),
          tags$ul(class = "mb-0 mt-1",
            tags$li("Export format: Precursor-level report"),
            tags$li("Required columns: PG.ProteinAccessions, PG.Genes, PEP.PeptidePosition, EG.ModifiedSequence"),
            tags$li("Per-sample columns: EG.TotalQuantity (Settings), ",
                    "EG.PTMPositions [Phospho (STY)], ",
                    "EG.PTMProbabilities [Phospho (STY)], ",
                    "EG.PTMSites [Phospho (STY)]"),
            tags$li("No imputation in Spectronaut"),
            tags$li("Q-value filtering applied (FDR < 0.01)")
          )
        ),
        if (.wiz_ptm_only())
          tags$div(class = "alert alert-info py-2",
            icon("circle-info"), tags$strong(" PTM only: "),
            "no global protein file. Sites are analysed as absolute PTM abundances, without ",
            "stoichiometric adjustment."),
        fluidRow(
          column(6,
            tags$h6("Phospho-Enriched File:"),
            fileInput("dia_ptm_wiz_phospho_file",
                      "Select phospho report (.tsv / .txt / .csv):",
                      accept = c(".tsv", ".txt", ".csv")),
            uiOutput("dia_ptm_wiz_phospho_status")
          ),
          if (!.wiz_ptm_only())
            column(6,
              tags$h6("Global Proteomics File:"),
              fileInput("dia_ptm_wiz_global_file",
                        "Select global protein report (.tsv / .txt / .csv):",
                        accept = c(".tsv", ".txt", ".csv")),
              uiOutput("dia_ptm_wiz_global_status")
            )
        )
      )
    } else if (sw == "spectronaut_pivot") {
      tagList(
        tags$h5("Upload Spectronaut DIA PTM Pivot Files"),
        tags$div(class = "alert alert-warning mb-3",
          icon("triangle-exclamation"),
          tags$strong(" Precursor-level export is strongly recommended for Spectronaut DIA data."),
          tags$br(),
          "The PTM Pivot Table has already summarized intensities across precursors inside Spectronaut, ",
          "which prevents TMP-based normalization and limits quantitative flexibility. ",
          "If your Spectronaut project allows it, go to ",
          tags$strong("Report > Phospho (Site) Pivot Report"),
          " and switch to ",
          tags$strong("Precursor-level"),
          " export, then use the ",
          tags$em("DIA \u2014 Spectronaut"),
          " option instead.",
          tags$br(), tags$br(),
          "The pivot format is supported and uses kNN imputation + limma eBayes, ",
          "but will produce less precise stoichiometric ratios."
        ),
        tags$div(class = "alert alert-secondary mb-3",
          tags$strong("Spectronaut PTM Pivot Report requirements:"),
          tags$ul(class = "mb-0 mt-1",
            tags$li("PTM file: site-level pivot with C: PTM.CollapseKey, T: PTM.ProteinId, T: PG.Genes"),
            tags$li("Per-sample PTM columns: [N] SampleName.PTM.Quantity"),
            tags$li("Global file: protein-group pivot with [N] SampleName.PG.Quantity"),
            tags$li("Global ID columns: T: PG.ProteinGroups, T: PG.Genes"),
            tags$li("Non-Phospho rows (e.g., Carbamidomethyl) will be auto-removed"),
            tags$li(tags$strong("DIA pipeline: "),
                    "kNN imputation on both matrices before stoichiometric adjustment (no MinProb downshift).")
          )
        ),
        fluidRow(
          column(6,
            tags$h6("PTM Pivot File (site-level):"),
            fileInput("dia_ptm_wiz_phospho_file",
                      "Select PTM pivot report (.tsv / .txt / .csv):",
                      accept = c(".tsv", ".txt", ".csv")),
            uiOutput("dia_ptm_wiz_phospho_status")
          ),
          column(6,
            tags$h6("Global Proteomics Pivot File:"),
            fileInput("dia_ptm_wiz_global_file",
                      "Select global protein pivot (.tsv / .txt / .csv):",
                      accept = c(".tsv", ".txt", ".csv")),
            uiOutput("dia_ptm_wiz_global_status")
          )
        )
      )
    } else {
      tagList(
        tags$h5("Upload DDA PTM Files (MaxQuant / Perseus)"),
        tags$div(class = "alert alert-secondary mb-3",
          tags$strong("MaxQuant / Perseus DDA PTM Report requirements:"),
          tags$ul(class = "mb-0 mt-1",
            tags$li("PTM file: site-level export with columns PTM.CollapseKey, PTM.ModificationTitle, PTM.ProteinId"),
            tags$li("Per-sample PTM columns: [N] SampleName.PTM.Quantity  or  N] SampleName (bare format)"),
            tags$li("Global file: protein-group level with [N] SampleName.PG.Quantity or N] SampleName"),
            tags$li("ID columns: T: PG.ProteinGroups, T: PG.Genes"),
            tags$li("Non-Phospho PTM rows (e.g., Carbamidomethyl) will be auto-removed"),
            tags$li("Probability column (PTM.LocalizationProbability) is optional — ",
                    "if absent, all sites are treated as passing the threshold"),
            tags$li(tags$strong("DDA pipeline: "),
                    "MinProb imputation (Perseus-style, per-column downshifted Gaussian) ",
                    "applied to both matrices before stoichiometric adjustment.")
          )
        ),
        fluidRow(
          column(6,
            tags$h6("PTM File (site-level):"),
            fileInput("dia_ptm_wiz_phospho_file",
                      "Select PTM report (.tsv / .txt / .csv):",
                      accept = c(".tsv", ".txt", ".csv")),
            uiOutput("dia_ptm_wiz_phospho_status")
          ),
          column(6,
            tags$h6("Global Proteomics File:"),
            fileInput("dia_ptm_wiz_global_file",
                      "Select global protein report (.tsv / .txt / .csv):",
                      accept = c(".tsv", ".txt", ".csv")),
            uiOutput("dia_ptm_wiz_global_status")
          )
        )
      )
    }
  }

  output$dia_ptm_wiz_phospho_status <- renderUI({
    d  <- rv$dia_ptm_phospho_raw
    sw <- rv$dia_ptm_software %||% "spectronaut"
    if (is.null(d))  return(tags$p(class = "text-muted small", "No file uploaded yet."))
    if (!is.null(d$error))
      return(tags$div(class = "alert alert-danger py-2", icon("circle-xmark"), " ", d$error))
    if (sw %in% c("maxquant", "spectronaut_pivot")) {
      no_prob_msg <- if (sw == "spectronaut_pivot")
        "No per-sample probability column detected. Sites are assumed pre-filtered by Spectronaut during export."
      else
        "No probability column detected. All sites will be treated as above threshold. Multiplicity disambiguation (M1 vs M2) will not be available."
      tagList(
        tags$div(class = "alert alert-success py-2",
          icon("circle-check"), " Parsed: ",
          tags$strong(d$n_rows, " site rows"), " \u00d7 ",
          tags$strong(d$n_samples, " samples.")),
        if (d$n_dropped_nonphospho > 0)
          tags$div(class = "alert alert-info py-2",
            icon("circle-info"), " ",
            tags$strong(d$n_dropped_nonphospho), " non-Phospho rows removed (e.g., Carbamidomethyl)."),
        if (!d$has_probability)
          tags$div(class = "alert alert-warning py-2",
            icon("triangle-exclamation"),
            tags$strong(" No probability column detected."),
            " ", no_prob_msg)
      )
    } else {
      tags$div(class = "alert alert-success py-2",
        icon("circle-check"), " Parsed: ",
        tags$strong(d$n_rows, " precursor rows"), " \u00d7 ",
        tags$strong(d$n_samples, " samples."))
    }
  })

  output$dia_ptm_wiz_global_status <- renderUI({
    d <- rv$dia_ptm_global_raw
    if (is.null(d))  return(tags$p(class = "text-muted small", "No file uploaded yet."))
    if (!is.null(d$error))
      return(tags$div(class = "alert alert-danger py-2", icon("circle-xmark"), " ", d$error))
    tags$div(class = "alert alert-success py-2",
      icon("circle-check"), " Parsed: ",
      tags$strong(d$n_proteins, " proteins"), " \u00d7 ",
      tags$strong(d$n_samples, " samples."))
  })

  observeEvent(input$dia_ptm_wiz_phospho_file, {
    req(input$dia_ptm_wiz_phospho_file)
    f  <- input$dia_ptm_wiz_phospho_file
    sw <- rv$dia_ptm_software %||% "spectronaut"
    result <- if (sw %in% c("maxquant", "spectronaut_pivot")) {
      tryCatch(
        parse_maxquant_ptm(f$datapath, logger = logger),
        error = function(e) list(error = conditionMessage(e))
      )
    } else {
      tryCatch(
        parse_spectronaut_phospho(f$datapath, logger = logger),
        error = function(e) list(error = conditionMessage(e))
      )
    }
    rv$dia_ptm_phospho_raw      <- result
    rv$dia_ptm_phospho_filename <- f$name
  })

  observeEvent(input$dia_ptm_wiz_global_file, {
    req(input$dia_ptm_wiz_global_file)
    f  <- input$dia_ptm_wiz_global_file
    sw <- rv$dia_ptm_software %||% "spectronaut"
    result <- if (sw %in% c("maxquant", "spectronaut_pivot")) {
      tryCatch(
        parse_maxquant_global(f$datapath, logger = logger),
        error = function(e) list(error = conditionMessage(e))
      )
    } else {
      tryCatch(
        parse_spectronaut_pivot(f$datapath, dedup_order = "always", logger = logger),
        error = function(e) list(error = conditionMessage(e))
      )
    }
    rv$dia_ptm_global_raw      <- result
    rv$dia_ptm_global_filename <- f$name
  })

  # ---------------------------------------------------------------------------
  # STEP 3 — Condition Groups
  # ---------------------------------------------------------------------------
  .wiz_step3_ui <- function() {
    tagList(
      tags$h5("Assign Condition Groups"),
      tags$p(class = "text-muted small",
        "Assign each sample to a condition group. Minimum 2 groups, ",
        "at least one must have >=2 samples. Sample names come from the phospho file."),
      fluidRow(
        column(5,
          textInput("dia_ptm_wiz_new_group_name", "New group name:",
                    placeholder = "e.g. Control"),
          checkbox_filter_controls("dia_ptm_wiz_sample_search", "dia_ptm_wiz_sample_sel"),
          uiOutput("dia_ptm_wiz_sample_checklist"),
          actionButton("dia_ptm_wiz_add_group_btn", "Add Group",
                       class = "btn-primary btn-sm mt-1", icon = icon("plus"))
        ),
        column(7,
          tags$strong("Defined groups:"),
          uiOutput("dia_ptm_wiz_groups_display")
        )
      ),
      uiOutput("dia_ptm_wiz_groups_status")
    )
  }

  output$dia_ptm_wiz_sample_checklist <- renderUI({
    all_samples <- rv$dia_ptm_phospho_raw$sample_names
    if (is.null(all_samples))
      return(tags$p(class = "text-muted small", "Upload phospho file first."))
    assigned  <- unlist(rv$dia_ptm_cond_groups)
    available <- setdiff(all_samples, assigned)
    if (length(available) == 0)
      return(tags$p(class = "text-muted small", "All samples assigned."))
    checkboxGroupInput("dia_ptm_wiz_sample_sel", "Available samples:",
                       choices = available, selected = NULL)
  })

  output$dia_ptm_wiz_groups_display <- renderUI({
    grps <- rv$dia_ptm_cond_groups
    if (length(grps) == 0)
      return(tags$p(class = "text-muted small", "No groups defined."))
    tagList(lapply(names(grps), function(gname) {
      tags$div(class = "card mb-1",
        tags$div(class = "card-body py-1 px-2 d-flex justify-content-between align-items-center",
          tags$span(tags$strong(gname), ": ",
                    paste(grps[[gname]], collapse = ", ")),
          actionButton(paste0("dia_ptm_wiz_rm_grp_", make.names(gname)),
                       label = "", icon = icon("trash"),
                       class = "btn-danger btn-sm",
                       `data-grpname` = gname)
        )
      )
    }))
  })

  observeEvent(input$dia_ptm_wiz_add_group_btn, {
    gname <- trimws(input$dia_ptm_wiz_new_group_name %||% "")
    samps <- input$dia_ptm_wiz_sample_sel
    if (nchar(gname) == 0 || length(samps) == 0) return()
    grps          <- rv$dia_ptm_cond_groups
    grps[[gname]] <- samps
    rv$dia_ptm_cond_groups <- grps
    updateTextInput(session, "dia_ptm_wiz_new_group_name", value = "")
    updateTextInput(session, "dia_ptm_wiz_sample_search",  value = "")
  })

  observe({
    grps <- rv$dia_ptm_cond_groups
    lapply(names(grps), function(gname) {
      btn_id <- paste0("dia_ptm_wiz_rm_grp_", make.names(gname))
      observeEvent(input[[btn_id]], {
        grps_now          <- rv$dia_ptm_cond_groups
        grps_now[[gname]] <- NULL
        rv$dia_ptm_cond_groups <- grps_now
      }, ignoreInit = TRUE, once = TRUE)
    })
  })

  output$dia_ptm_wiz_groups_status <- renderUI({
    grps     <- rv$dia_ptm_cond_groups
    n_grps   <- length(grps)
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

  # ---------------------------------------------------------------------------
  # STEP 4 — Paired Samples
  # ---------------------------------------------------------------------------
  .wiz_step4_ui <- function() {
    tagList(
      tags$h5("Paired / Repeated Measures Design (Optional)"),
      checkboxInput("dia_ptm_wiz_pairing_enabled",
                    "Enable paired/repeated measures design", value = FALSE),
      conditionalPanel(
        condition = "input.dia_ptm_wiz_pairing_enabled",
        tags$div(class = "alert alert-info",
          icon("circle-info"),
          " Samples from the same subject must appear in >=2 conditions. ",
          "limma uses duplicateCorrelation() to account for within-subject correlation."
        ),
        uiOutput("dia_ptm_wiz_pairing_inputs"),
        uiOutput("dia_ptm_wiz_pairing_validation")
      )
    )
  }

  output$dia_ptm_wiz_pairing_inputs <- renderUI({
    req(isTRUE(input$dia_ptm_wiz_pairing_enabled))
    samples <- unlist(rv$dia_ptm_cond_groups)
    if (length(samples) == 0) return(NULL)
    tagList(
      tags$p(class = "text-muted small",
             "Assign a Subject ID to each sample (shared ID = same subject):"),
      lapply(samples, function(sn) {
        fluidRow(
          column(6, tags$label(class = "col-form-label col-form-label-sm", sn)),
          column(6, selectizeInput(
            paste0("dia_ptm_wiz_pair_", make.names(sn)),
            label    = NULL,
            choices  = character(0),
            options  = list(create = TRUE, placeholder = "Subject ID")
          ))
        )
      })
    )
  })

  output$dia_ptm_wiz_pairing_validation <- renderUI({
    req(isTRUE(input$dia_ptm_wiz_pairing_enabled))
    samples <- unlist(rv$dia_ptm_cond_groups)
    if (length(samples) == 0) return(NULL)
    grp_map <- setNames(
      rep(names(rv$dia_ptm_cond_groups),
          vapply(rv$dia_ptm_cond_groups, length, integer(1))),
      unlist(rv$dia_ptm_cond_groups)
    )
    subject_ids <- vapply(samples, function(sn)
      input[[paste0("dia_ptm_wiz_pair_", make.names(sn))]] %||% "",
      character(1))
    bad_subjects <- character(0)
    for (subj in unique(subject_ids[subject_ids != ""])) {
      subj_samps <- samples[subject_ids == subj]
      if (length(unique(grp_map[subj_samps])) < 2L)
        bad_subjects <- c(bad_subjects, subj)
    }
    if (length(bad_subjects) > 0)
      return(tags$div(class = "alert alert-warning mt-2",
        icon("triangle-exclamation"), " Subject(s) ",
        paste(bad_subjects, collapse = ", "),
        " appear in fewer than 2 conditions."))
    if (all(subject_ids == ""))
      return(tags$div(class = "alert alert-warning mt-2",
        icon("triangle-exclamation"), " No subject IDs assigned."))
    tags$div(class = "alert alert-success mt-2",
             icon("circle-check"), " Pairing looks valid.")
  })

  # ---------------------------------------------------------------------------
  # STEP 5 — Probability Threshold
  # ---------------------------------------------------------------------------
  .wiz_step5_ui <- function() {
    sw <- rv$dia_ptm_software %||% "spectronaut"
    if (sw == "spectronaut") {
      tagList(
        tags$h5("PTM Site Probability Threshold"),
        tags$p(class = "text-muted small",
          "Spectronaut reports a localization probability for each phosphosite. ",
          "Sites where the maximum probability across all samples falls below ",
          "this threshold are excluded from analysis. ",
          "Once a site passes, its intensity is collected from every matched sample, ",
          "including samples in which that site's own localization probability was ",
          "absent or below this threshold."),
        sliderInput("dia_ptm_wiz_prob_threshold",
                    "Minimum localization probability (%):",
                    min = 0, max = 100,
                    value = isolate(rv$dia_ptm_prob_threshold %||% 75), step = 1),
        tags$div(class = "alert alert-info",
          icon("circle-info"),
          " Recommended: 75%. This retains high-confidence phosphosites while ",
          "keeping reasonable coverage. Higher values increase specificity; ",
          "lower values increase sensitivity.")
      )
    } else if (sw == "spectronaut_pivot") {
      # Spectronaut applies localization probability filtering during export.
      # The pivot format typically lacks a per-sample probability column, so
      # sites are already pre-filtered — inform the user rather than warning.
      has_prob <- isTRUE(rv$dia_ptm_phospho_raw$has_probability)
      tagList(
        tags$h5("PTM Site Probability Threshold"),
        tags$div(class = "alert alert-info mb-3",
          icon("circle-info"),
          tags$strong(" Spectronaut PTM Pivot: "),
          "Localization probability filtering is typically applied inside Spectronaut before export. ",
          if (!has_prob)
            "No per-sample probability column was detected in your file \u2014 all sites will be treated as above threshold (pre-filtered by Spectronaut)."
          else
            "A probability column was detected and the threshold below will be applied."
        ),
        if (has_prob)
          sliderInput("dia_ptm_wiz_prob_threshold",
                      "Minimum localization probability (%):",
                      min = 0, max = 100,
                      value = isolate(rv$dia_ptm_prob_threshold %||% 75), step = 1)
        else
          sliderInput("dia_ptm_wiz_prob_threshold",
                      "Minimum localization probability (%) [inactive \u2014 no column present]:",
                      min = 0, max = 100, value = 0, step = 5)
      )
    } else {
      # MaxQuant — check whether probability column was detected
      has_prob <- isTRUE(rv$dia_ptm_phospho_raw$has_probability)
      if (has_prob) {
        tagList(
          tags$h5("PTM Site Probability Threshold"),
          tags$div(class = "alert alert-success mb-3",
            icon("circle-check"),
            tags$strong(" Probability column detected in your MaxQuant file."),
            " Sites below the threshold will be excluded."
          ),
          sliderInput("dia_ptm_wiz_prob_threshold",
                      "Minimum localization probability (%):",
                      min = 0, max = 100,
                      value = isolate(rv$dia_ptm_prob_threshold %||% 75), step = 1),
          tags$div(class = "alert alert-info",
            icon("circle-info"),
            " Recommended: 75%.")
        )
      } else {
        tagList(
          tags$h5("PTM Site Probability Threshold"),
          tags$div(class = "alert alert-warning",
            icon("triangle-exclamation"),
            tags$strong(" No localization probability column was found in your MaxQuant file."),
            tags$br(),
            "This typically means the data was pre-filtered before export. ",
            "As a result:",
            tags$ul(class = "mt-2 mb-1",
              tags$li("All sites in the file will be treated as above threshold."),
              tags$li("Multiplicity disambiguation (M1 vs M2) is not possible; ",
                      "each multiplicity entry is treated as an independent site."),
              tags$li("The probability threshold slider below has no effect.")
            ),
            tags$strong("Please confirm you understand before proceeding.")
          ),
          checkboxInput("dia_ptm_wiz_noprob_acknowledged",
                        "I understand — proceed without probability filtering.",
                        value = FALSE),
          sliderInput("dia_ptm_wiz_prob_threshold",
                      "Minimum localization probability (%) [inactive — no column present]:",
                      min = 0, max = 100, value = 0, step = 5)
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # STEP 6 — Summary & Summarize
  # ---------------------------------------------------------------------------
  .wiz_step6_ui <- function() {
    sw           <- rv$dia_ptm_software %||% "spectronaut"
    prob_thresh  <- rv$dia_ptm_prob_threshold %||%
                    input$dia_ptm_wiz_prob_threshold %||% 75
    n_phospho    <- rv$dia_ptm_phospho_raw$n_rows    %||% 0L
    n_global     <- rv$dia_ptm_global_raw$n_proteins %||% 0L
    n_samples    <- rv$dia_ptm_phospho_raw$n_samples %||% 0L
    n_grps       <- length(rv$dia_ptm_cond_groups)
    paired       <- isTRUE(input$dia_ptm_wiz_pairing_enabled)
    ptm_only     <- .wiz_ptm_only()
    tagList(
      tags$h5("Summary & Summarize"),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$tbody(
              tags$tr(tags$td(tags$strong("Mode")),
                      tags$td(if (ptm_only) "PTM only (no global proteome)" else "PTM with global proteomics")),
              tags$tr(tags$td(tags$strong("Software")),
                      tags$td(switch(sw,
                        "maxquant"          = "MaxQuant / Perseus (DDA)",
                        "spectronaut_pivot" = "Spectronaut PTM Pivot (DIA)",
                        "Spectronaut (DIA)"))),
              tags$tr(tags$td(tags$strong(if (sw %in% c("maxquant", "spectronaut_pivot")) "PTM file" else "Phospho file")),
                      tags$td(rv$dia_ptm_phospho_filename %||% "(none)")),
              if (!ptm_only)
                tags$tr(tags$td(tags$strong("Global file")),
                        tags$td(rv$dia_ptm_global_filename %||% "(none)")),
              tags$tr(tags$td(tags$strong(if (sw %in% c("maxquant", "spectronaut_pivot")) "PTM site rows" else "Phospho rows")),
                      tags$td(n_phospho)),
              if (sw %in% c("maxquant", "spectronaut_pivot") && isTRUE(rv$dia_ptm_phospho_raw$n_dropped_nonphospho > 0))
                tags$tr(tags$td(tags$strong("Non-Phospho rows dropped")),
                        tags$td(rv$dia_ptm_phospho_raw$n_dropped_nonphospho)),
              if (!ptm_only)
                tags$tr(tags$td(tags$strong("Global proteins")),
                        tags$td(n_global)),
              tags$tr(tags$td(tags$strong("Samples")),
                      tags$td(n_samples)),
              tags$tr(tags$td(tags$strong("Groups")),
                      tags$td(paste(names(rv$dia_ptm_cond_groups), collapse = ", "))),
              tags$tr(tags$td(tags$strong("Pairing")),
                      tags$td(if (paired) "Yes" else "No")),
              tags$tr(tags$td(tags$strong("Prob. threshold")),
                      tags$td(if (sw %in% c("maxquant", "spectronaut_pivot") && !isTRUE(rv$dia_ptm_phospho_raw$has_probability))
                                "N/A (no probability column \u2014 all sites pass)"
                              else paste0(prob_thresh, "%"))),
              tags$tr(tags$td(tags$strong("Summarization")),
                      tags$td(uiOutput("dia_ptm_wiz_summ_row", inline = TRUE))),
              tags$tr(tags$td(tags$strong("DE engine")),
                      tags$td(if (ptm_only) "limma (eBayes: robust=TRUE, trend=TRUE) on absolute PTM"
                              else "limma (eBayes: robust=TRUE; trend=TRUE for absolute PTM and protein, trend=FALSE for adjusted)"))
            )
          )
        )
      ),
      if (sw %in% c("spectronaut", "spectronaut_pivot")) .ptm_norm_question_ui(),
      actionButton("dia_ptm_wiz_summarize_btn", "Summarize",
                   class = "btn-success w-100", icon = icon("play")),
      uiOutput("dia_ptm_wiz_summarize_status")
    )
  }

  # ---------------------------------------------------------------------------
  # Step 6: were the Spectronaut exports already cross-run normalized?
  # One question per export (PTM and global), both default Yes (Spectronaut
  # normalizes across runs by default). dia_normalization_check() runs on each
  # uploaded matrix with the assigned groups and is shown next to its box.
  #   PTM yes    -> no median centering of the site matrix
  #   global yes -> no median centering of the global protein matrix
  #   both yes   -> no centering of the adjusted (site - protein) matrix either
  # Any "no" keeps that matrix's centering, and the adjusted-matrix centering.
  # PTM-only mode asks the PTM question alone.
  # MaxQuant (DDA) uploads are not affected.
  # ---------------------------------------------------------------------------
  .ptm_norm_question_ui <- function() {
    cache    <- rv$dia_ptm_master_cache
    ptm_only <- .wiz_ptm_only()
    tags$div(class = "card mb-3",
      tags$div(class = "card-body",
        tags$strong("Cross-run normalization in Spectronaut"),
        if (ptm_only)
          tags$p(class = "small mb-2 mt-1",
            "Spectronaut normalizes quantities across runs by default (Quantification > Cross Run ",
            "Normalization). Normalizing such data a second time can shift sites that did not change. ",
            "If the export was normalized, OSMSS uses the quantities as they are and does not apply its ",
            "own median centering.")
        else
          tags$p(class = "small mb-2 mt-1",
            "Spectronaut normalizes quantities across runs by default (Quantification > Cross Run ",
            "Normalization, set separately in the PTM-enriched and the global experiment). Normalizing ",
            "such data a second time can shift sites and proteins that did not change. For each export ",
            "that was normalized, OSMSS uses the quantities as they are and does not apply its own median ",
            "centering. When both were normalized, the site-minus-protein (stoichiometry-adjusted) values ",
            "are not re-centered either."),
        tags$strong(class = "small", "PTM-enriched export"),
        uiOutput("dia_ptm_wiz_norm_check_ptm"),
        checkboxInput("dia_ptm_wiz_ptm_normalized",
          "The PTM export was cross-run normalized (Spectronaut default)",
          value = if (is.null(cache)) TRUE else !isTRUE(cache$ptm_centering_applied %||% TRUE)),
        if (!ptm_only) tagList(
          tags$strong(class = "small", "Global protein export"),
          uiOutput("dia_ptm_wiz_norm_check_global"),
          checkboxInput("dia_ptm_wiz_global_normalized",
            "The global protein export was cross-run normalized (Spectronaut default)",
            value = if (is.null(cache)) TRUE else !isTRUE(cache$global_centering_applied %||% TRUE))),
        tags$p(class = "text-muted small mb-0",
          if (ptm_only)
            paste0("Untick the box only if you know cross-run normalization was switched off in ",
                   "Spectronaut. OSMSS then applies column-wise median centering (log2) to the site matrix.")
          else
            paste0("Untick a box only if you know cross-run normalization was switched off for that ",
                   "experiment in Spectronaut. OSMSS then applies column-wise median centering (log2) to ",
                   "that matrix and to the adjusted values."))
      )
    )
  }

  # Groups of the global file: PTM sample names, or positional matching when the
  # names differ (same rule as maxquant_ptm_summarize()).
  .ptm_global_groups <- function(groups, ptm_samples, global_samples) {
    if (sum(unlist(groups) %in% global_samples) >= 2L ||
        length(ptm_samples) != length(global_samples)) return(groups)
    map <- setNames(global_samples, ptm_samples)
    lapply(groups, function(s) unname(map[s][!is.na(map[s])]))
  }

  dia_ptm_norm_checks <- reactive({
    req(identical(rv$dia_ptm_wizard_step, 6L))
    sw <- rv$dia_ptm_software %||% "spectronaut"
    req(sw %in% c("spectronaut", "spectronaut_pivot"))
    ph <- rv$dia_ptm_phospho_raw; gl <- rv$dia_ptm_global_raw
    req(!is.null(ph), is.null(ph$error))
    if (!.wiz_ptm_only()) req(!is.null(gl), is.null(gl$error))
    groups <- rv$dia_ptm_cond_groups
    list(
      ptm    = tryCatch(dia_normalization_check(ptm_quant_matrix(ph), groups,
                          feature_label = if (sw == "spectronaut") "precursors" else "PTM sites"),
                        error = function(e) list(verdict = "inconclusive", message = conditionMessage(e))),
      global = if (!.wiz_ptm_only())
                 tryCatch(dia_normalization_check(gl$protein_matrix,
                            .ptm_global_groups(groups, ph$sample_names, gl$sample_names)),
                          error = function(e) list(verdict = "inconclusive", message = conditionMessage(e))))
  })

  .norm_check_alert <- function(chk) {
    cls <- switch(chk$verdict, normalized = "alert-success", borderline = "alert-warning",
                  not_normalized = "alert-danger", "alert-secondary")
    ic  <- switch(chk$verdict, normalized = "circle-check", inconclusive = "circle-info",
                  "triangle-exclamation")
    tags$div(class = paste("alert py-2 px-3 small mb-1", cls),
      icon(ic), tags$strong(" Data check: "), chk$message,
      if (chk$verdict %in% c("borderline", "not_normalized"))
        tags$span(" If cross-run normalization was switched off in Spectronaut, untick the box below."))
  }
  output$dia_ptm_wiz_norm_check_ptm    <- renderUI(.norm_check_alert(dia_ptm_norm_checks()$ptm))
  output$dia_ptm_wiz_norm_check_global <- renderUI(.norm_check_alert(dia_ptm_norm_checks()$global))

  # What the user answered (FALSE for MaxQuant uploads: earlier behaviour; global
  # FALSE in PTM-only mode, where there is no global export)
  .ptm_norm_flags <- function() {
    sw <- rv$dia_ptm_software %||% "spectronaut"
    sp <- sw %in% c("spectronaut", "spectronaut_pivot")
    list(ptm    = sp && isTRUE(input$dia_ptm_wiz_ptm_normalized    %||% TRUE),
         global = sp && !.wiz_ptm_only() && isTRUE(input$dia_ptm_wiz_global_normalized %||% TRUE))
  }

  output$dia_ptm_wiz_summ_row <- renderUI({
    f <- .ptm_norm_flags()
    paste0("Sum intensities/site → log₂ → ",
           if (f$ptm) "no PTM centering (normalized upstream)" else "PTM median-centering",
           if (.wiz_ptm_only()) " → no stoichiometric adj. (PTM only)"
           else paste0(" → stoichiometric adj.", if (f$ptm && f$global) " (not re-centered)" else ""),
           " (≥2 valid obs/condition retained)")
  })

  output$dia_ptm_wiz_summarize_status <- renderUI({
    if (is.null(rv$dia_ptm_master_cache)) {
      tags$p(class = "text-muted small mt-2",
             "Click 'Summarize' to process the data.")
    } else {
      m <- rv$dia_ptm_master_cache
      tags$div(class = "alert alert-success mt-2",
        icon("circle-check"), " Summarization complete: ",
        tags$strong(m$n_sites_output, " sites"),
        if (ptm_has_global(m))
          paste0(" (", m$n_sites_adjusted, " adjusted rows, including ",
                 m$n_shared_expanded, " from shared peptide expansion).")
        else " (PTM only: absolute sites, no stoichiometric adjustment).",
        if (isFALSE(m$ptm_centering_applied) || isFALSE(m$global_centering_applied))
          paste0(" Median centering skipped for: ",
                 paste(c(if (isFALSE(m$ptm_centering_applied)) "PTM matrix",
                         if (isFALSE(m$global_centering_applied)) "global matrix",
                         if (isFALSE(m$adj_centering_applied)) "adjusted matrix"), collapse = ", "),
                 " (normalized upstream).")
      )
    }
  })

  # The threshold slider is rendered inside a renderUI for step 5. When the
  # wizard advances, that UI leaves the DOM and Shiny sets the input to NULL.
  # Mirror it into rv the moment it changes so the value survives navigation.
  observeEvent(input$dia_ptm_wiz_prob_threshold, {
    v <- suppressWarnings(as.numeric(input$dia_ptm_wiz_prob_threshold))
    if (!is.na(v)) rv$dia_ptm_prob_threshold <- v
  }, ignoreInit = FALSE)

  observeEvent(input$dia_ptm_wiz_summarize_btn, {
    phospho_raw <- rv$dia_ptm_phospho_raw
    global_raw  <- rv$dia_ptm_global_raw
    if (is.null(phospho_raw) || !is.null(phospho_raw$error)) {
      showNotification("Please upload a valid phospho file first.", type = "warning")
      return()
    }
    ptm_only <- .wiz_ptm_only()
    if (!ptm_only && (is.null(global_raw) || !is.null(global_raw$error))) {
      showNotification("Please upload a valid global proteomics file first.", type = "warning")
      return()
    }
    if (length(rv$dia_ptm_cond_groups) < 2L) {
      showNotification("Define at least 2 condition groups first.", type = "warning")
      return()
    }

    # Collect subject pairing
    subj_pairing <- NULL
    if (isTRUE(input$dia_ptm_wiz_pairing_enabled)) {
      samples <- unlist(rv$dia_ptm_cond_groups)
      subj_pairing <- setNames(
        vapply(samples, function(sn)
          input[[paste0("dia_ptm_wiz_pair_", make.names(sn))]] %||% sn,
          character(1)),
        samples
      )
      rv$dia_ptm_subject_pairing <- subj_pairing
    } else {
      rv$dia_ptm_subject_pairing <- NULL
    }

    sw          <- rv$dia_ptm_software %||% "spectronaut"
    has_prob    <- isTRUE(phospho_raw$has_probability)
    prob_thresh <- if (sw %in% c("maxquant", "spectronaut_pivot") && !has_prob) 0 else {
      .pt <- rv$dia_ptm_prob_threshold
      if (is.null(.pt) || is.na(.pt)) {
        stop("Localization probability threshold is unset. This is a bug: the ",
             "slider lives in a renderUI on step 5 and its input is NULL once ",
             "that UI leaves the DOM. Refusing to silently substitute a default.")
      }
      as.numeric(.pt)
    }

    norm_flags <- .ptm_norm_flags()
    rv$dia_ptm_norm_check <- if (sw %in% c("spectronaut", "spectronaut_pivot"))
      tryCatch(isolate({
        list(ptm    = dia_normalization_check(ptm_quant_matrix(phospho_raw), rv$dia_ptm_cond_groups,
                       feature_label = if (sw == "spectronaut") "precursors" else "PTM sites"),
             global = if (!ptm_only)
                        dia_normalization_check(global_raw$protein_matrix,
                          .ptm_global_groups(rv$dia_ptm_cond_groups, phospho_raw$sample_names, global_raw$sample_names)))
      }), error = function(e) NULL) else NULL
    withProgress(message = "Summarizing PTM data...", value = 0.1, {
      incProgress(0.1, detail = "Parsing observations...")
      result <- if (sw %in% c("maxquant", "spectronaut_pivot")) {
        tryCatch(
          maxquant_ptm_summarize(
            ptm_parsed      = phospho_raw,
            global_parsed   = global_raw,
            cond_groups     = rv$dia_ptm_cond_groups,
            prob_threshold  = prob_thresh,
            subject_pairing = subj_pairing,
            acquisition_mode = if (sw == "maxquant") (rv$dia_ptm_acquisition %||% "dda") else "dia",
            ptm_input_normalized    = norm_flags$ptm,
            global_input_normalized = norm_flags$global,
            logger          = logger
          ),
          error = function(e) list(error = conditionMessage(e))
        )
      } else {
        tryCatch(
          dia_ptm_summarize(
            phospho_parsed  = phospho_raw,
            global_parsed   = if (ptm_only) NULL else global_raw,
            cond_groups     = rv$dia_ptm_cond_groups,
            prob_threshold  = prob_thresh,
            subject_pairing = subj_pairing,
            ptm_input_normalized    = norm_flags$ptm,
            global_input_normalized = norm_flags$global,
            logger          = logger
          ),
          error = function(e) list(error = conditionMessage(e))
        )
      }
      incProgress(0.8, detail = "Done.")
    })

    if (!is.null(result$error)) {
      showNotification(paste0("Summarization failed: ", result$error), type = "error",
                       duration = 10)
      return()
    }

    rv$dia_ptm_master_cache <- result
    logger(paste0("[DIA PTM] Summarization stored (",
                  toupper(sw), if (ptm_only) ", PTM only" else "", "): ",
                  result$n_sites_output, " sites",
                  if (ptm_only) "." else paste0(", ", result$n_sites_adjusted, " adjusted rows.")))
  })

  # ---------------------------------------------------------------------------
  # STEP 7 — Review
  # ---------------------------------------------------------------------------
  .wiz_step7_ui <- function() {
    cache <- rv$dia_ptm_master_cache
    sw    <- rv$dia_ptm_software %||% "spectronaut"
    if (is.null(cache)) {
      return(tags$div(class = "alert alert-warning",
        icon("triangle-exclamation"),
        " Summarization has not been run yet. Please go back to Step 6."))
    }
    tagList(
      tags$h5("Review & Finalize"),
      tags$div(class = "alert alert-success",
        icon("circle-check"),
        " Summarization complete. Click 'Done' to finalize and route data to PTM modules."
      ),
      if (sw == "maxquant" && !isTRUE(cache$has_probability))
        tags$div(class = "alert alert-warning",
          icon("triangle-exclamation"),
          tags$strong(" Note: "), "No probability column was present in the MaxQuant file. ",
          "All sites were treated as above threshold. Multiplicity entries (M1, M2, etc.) ",
          "are treated as independent sites."),
      if (sw == "spectronaut_pivot" && !isTRUE(cache$has_probability))
        tags$div(class = "alert alert-info",
          icon("circle-info"),
          tags$strong(" Note: "), "No per-sample probability column was present. ",
          "Sites are assumed pre-filtered by Spectronaut during export."),
      tags$div(class = "card",
        tags$div(class = "card-body",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$tbody(
              tags$tr(tags$td(tags$strong("Mode")),
                      tags$td(if (ptm_has_global(cache)) "PTM with global proteomics"
                              else "PTM only (no global proteome)")),
              tags$tr(tags$td(tags$strong("Input type")),
                      tags$td(switch(sw,
                        "maxquant"          = "DDA \u2014 MaxQuant / Perseus",
                        "spectronaut_pivot" = "DIA \u2014 Spectronaut PTM Pivot",
                        "DIA \u2014 Spectronaut (Precursor-Level)"))),
              tags$tr(tags$td(tags$strong("Analysis pipeline")),
                      tags$td(switch(sw,
                        "maxquant"          = "MinProb imputation \u2192 stoichiometric adj. \u2192 limma eBayes",
                        "spectronaut_pivot" = "kNN imputation \u2192 stoichiometric adj. \u2192 limma eBayes",
                        if (ptm_has_global(cache)) "TMP summarization \u2192 stoichiometric adj. \u2192 limma eBayes"
                        else "TMP summarization \u2192 limma eBayes (absolute PTM)"))),
              tags$tr(tags$td(tags$strong("Input rows")),
                      tags$td(cache$n_sites_input)),
              tags$tr(tags$td(tags$strong("Sites output")),
                      tags$td(cache$n_sites_output)),
              if (ptm_has_global(cache))
                tags$tr(tags$td(tags$strong("Adjusted rows (stoichiometry)")),
                        tags$td(cache$n_sites_adjusted)),
              if (ptm_has_global(cache))
                tags$tr(tags$td(tags$strong("Shared peptide rows")),
                        tags$td(cache$n_shared_expanded)),
              tags$tr(tags$td(tags$strong("Samples")),
                      tags$td(length(cache$matched_samples))),
              tags$tr(tags$td(tags$strong("Conditions")),
                      tags$td(paste(names(cache$cond_groups), collapse = ", "))),
              tags$tr(tags$td(tags$strong("Prob. threshold used")),
                      tags$td(if (!isTRUE(cache$has_probability) && sw == "maxquant")
                                "N/A (pre-filtered)" else paste0(cache$prob_threshold, "%")))
            )
          )
        )
      ),
      if (!ptm_has_global(cache))
        tags$div(class = "alert alert-info",
          icon("circle-info"), tags$strong(" PTM only: "),
          "every PTM module works on the absolute site abundances. A site's fold change then ",
          "includes any change in its protein's abundance."),
      tags$p(class = "text-muted small mt-2",
        "After clicking 'Done', navigate to PTM DE, PTM Trends, or PTM Single ",
        "to begin analysis. All three modules will recognize this dataset.")
    )
  }

  # ---------------------------------------------------------------------------
  # Footer
  # ---------------------------------------------------------------------------
  output$dia_ptm_wiz_footer <- renderUI({
    step <- rv$dia_ptm_wizard_step %||% 1L
    back_btn <- if (step > 1L)
      actionButton("dia_ptm_wiz_back", "\u2190 Back", class = "btn-secondary")
    else NULL
    fwd_btn <- if (step < N_STEPS)
      actionButton("dia_ptm_wiz_next", "Next \u2192", class = "btn-primary")
    else
      actionButton("dia_ptm_wiz_done", "\u2714 Done", class = "btn-success")
    cancel_btn <- modalButton("Cancel")
    tagList(
      tags$div(class = "d-flex justify-content-between",
        tags$div(cancel_btn),
        tags$div(back_btn, fwd_btn)
      )
    )
  })

  # ---------------------------------------------------------------------------
  # Navigation: Next
  # ---------------------------------------------------------------------------
  observeEvent(input$dia_ptm_wiz_next, {
    step <- rv$dia_ptm_wizard_step %||% 1L
    if (step == 1L) {
      # Step 1: software selection — just require a choice was made
      if (is.null(rv$dia_ptm_software) || !identical(rv$dia_ptm_software, "spectronaut")) {
        showNotification("Please select an input software first.", type = "warning"); return()
      }
      rv$dia_ptm_wizard_step <- 2L

    } else if (step == 2L) {
      phospho_ok <- !is.null(rv$dia_ptm_phospho_raw) && is.null(rv$dia_ptm_phospho_raw$error)
      global_ok  <- !is.null(rv$dia_ptm_global_raw)  && is.null(rv$dia_ptm_global_raw$error)
      if (!phospho_ok) { showNotification("Upload a valid PTM file first.", type = "warning"); return() }
      if (!global_ok && !.wiz_ptm_only()) {
        showNotification("Upload a valid global proteomics file first.", type = "warning"); return()
      }
      rv$dia_ptm_wizard_step <- 3L

    } else if (step == 3L) {
      grps <- rv$dia_ptm_cond_groups
      if (length(grps) < 2L) {
        showNotification("Define at least 2 condition groups.", type = "warning"); return()
      }
      if (!any(vapply(grps, length, integer(1)) >= 2L)) {
        showNotification("At least one group must have >=2 samples.", type = "warning"); return()
      }
      rv$dia_ptm_wizard_step <- 4L

    } else if (step == 4L) {
      rv$dia_ptm_wizard_step <- 5L

    } else if (step == 5L) {
      # Probability step: for MaxQuant without prob column, require acknowledgement
      sw       <- rv$dia_ptm_software %||% "spectronaut"
      has_prob <- isTRUE(rv$dia_ptm_phospho_raw$has_probability)
      # Only require explicit acknowledgement for MaxQuant DDA (where missing prob column is unexpected)
      if (sw == "maxquant" && !has_prob && !isTRUE(input$dia_ptm_wiz_noprob_acknowledged)) {
        showNotification("Please check the acknowledgement box before proceeding.", type = "warning"); return()
      }
      rv$dia_ptm_wizard_step <- 6L

    } else if (step == 6L) {
      if (is.null(rv$dia_ptm_master_cache)) {
        showNotification("Run Summarize before proceeding.", type = "warning"); return()
      }
      # A normalization answer changed after summarizing: the cache no longer matches it
      .c <- rv$dia_ptm_master_cache; .f <- .ptm_norm_flags()
      if (!identical(isTRUE(.c$ptm_centering_applied %||% TRUE), !.f$ptm) ||
          !identical(isTRUE(.c$global_centering_applied %||% TRUE), !.f$global)) {
        showNotification("A normalization setting changed after summarizing. Please click Summarize again.",
                         type = "warning"); return()
      }
      rv$dia_ptm_wizard_step <- 7L
    }
    showModal(.wiz_build_modal())
  })

  # ---------------------------------------------------------------------------
  # Navigation: Back
  # ---------------------------------------------------------------------------
  observeEvent(input$dia_ptm_wiz_back, {
    step <- rv$dia_ptm_wizard_step %||% 1L
    rv$dia_ptm_wizard_step <- max(1L, step - 1L)
    showModal(.wiz_build_modal())
  })

  # ---------------------------------------------------------------------------
  # Done — finalize; set rv fields; trigger PTM module sync
  # ---------------------------------------------------------------------------
  observeEvent(input$dia_ptm_wiz_done, {
    cache <- rv$dia_ptm_master_cache
    if (is.null(cache)) {
      showNotification("Run Summarize first.", type = "warning")
      return()
    }

    # ── Set active data source ──
    rv$data_source              <- "dia_ptm"
    rv$dia_ptm_wizard_completed <- TRUE
    rv$dia_ptm_wizard_modal_open <- FALSE
    removeModal()

    # ── Populate rv$ptm_data for cross-module sync ──
    # These fields are read by server_ptm_de.R, server_ptm_single.R,
    # server_ptm_trends.R in their sync observers.
    # Reset source to NULL first so re-uploads (same mode or mode switch)
    # always produce a NULL → "dia_ptm_upload" transition, guaranteeing
    # the downstream observeEvent(rv$ptm_data$source) fires.
    rv$ptm_data$source          <- NULL
    rv$ptm_data$design_option   <- "dia_ptm"
    rv$ptm_data$conditions      <- sort(unique(cache$design_df$condition))
    rv$ptm_data$data_format     <- "dia_ptm"
    rv$ptm_data$n_files         <- if (ptm_has_global(cache)) 2L else 1L   # phospho (+ global)
    rv$ptm_data$gene_lookup     <- NULL
    rv$ptm_data$protein_headers <- NULL
    rv$ptm_data$protein_df      <- NULL
    rv$ptm_data$annotation_map  <- NULL
    rv$ptm_data$ptm_long        <- NULL  # no MSstatsPTM intermediate tables
    rv$ptm_data$protein_long    <- NULL

    # Setting $source last triggers the observeEvent(rv$ptm_data$source, ...) in each module
    rv$ptm_data$source          <- "dia_ptm_upload"

    sw  <- rv$dia_ptm_software    %||% "spectronaut"
    acq <- rv$dia_ptm_acquisition %||% "dia"
    has_global <- ptm_has_global(cache)

    # ── Provenance ──
    rv$provenance$dia_ptm <- list(
      ptm_mode            = cache$ptm_mode %||% "with_global",
      phospho_file        = rv$dia_ptm_phospho_filename,
      global_file         = if (has_global) rv$dia_ptm_global_filename,
      software            = switch(sw,
                            "maxquant"          = "MaxQuant / Perseus",
                            "spectronaut_pivot" = "Spectronaut PTM Pivot (DIA)",
                            "Spectronaut"),
      acquisition_mode    = if (acq == "dda") "DDA" else "DIA",
      data_level          = switch(sw, "maxquant" = "site", "spectronaut_pivot" = "site (pivot)", "precursor"),
      has_probability     = if (sw %in% c("maxquant", "spectronaut_pivot")) isTRUE(cache$has_probability) else TRUE,
      n_sites_input       = cache$n_sites_input,
      n_sites_passing_prob = cache$n_sites_passing_prob,
      n_sites_output      = cache$n_sites_output,
      n_sites_adjusted    = cache$n_sites_adjusted,
      n_shared_expanded   = cache$n_shared_expanded,
      n_no_parent_protein = cache$n_no_parent_protein,
      n_sites_multi_peptidoform_all      = cache$n_sites_multi_peptidoform_all,
      max_peptidoforms_per_site_all      = cache$max_peptidoforms_per_site_all,
      n_sites_multi_peptidoform_retained = cache$n_sites_multi_peptidoform_retained,
      max_peptidoforms_per_site_retained = cache$max_peptidoforms_per_site_retained,
      site_summarization  = cache$site_summarization,
      n_medpolish_fallback = cache$n_medpolish_fallback,
      n_obs_collected           = cache$n_obs_collected,
      n_obs_no_own_prob         = cache$n_obs_no_own_prob,
      n_obs_below_own_threshold = cache$n_obs_below_own_threshold,
      min_valid_per_cond  = cache$min_valid_per_cond,
      min_valid_per_cond_effective = cache$min_valid_per_cond_effective,
      prob_threshold      = if (sw %in% c("maxquant", "spectronaut_pivot") && !isTRUE(cache$has_probability)) "N/A"
                            else suppressWarnings(as.numeric(cache$prob_threshold)),
      prob_threshold_label = if (sw %in% c("maxquant", "spectronaut_pivot") && !isTRUE(cache$has_probability)) "N/A"
                            else paste0(cache$prob_threshold, "%"),
      n_samples           = length(cache$matched_samples),
      n_conditions        = length(cache$cond_groups),
      condition_names     = names(cache$cond_groups),
      imputation          = paste0(
                              "No imputation is performed at summarization on the DIA acquisition paths; ",
                              "the DDA MaxQuant path applies MinProb at summarization before stoichiometric ",
                              "subtraction. All DE-stage imputation is recorded in the DE provenance block ",
                              "(rv$provenance$ptm_de: imputation_scope, imputation_eligibility_note, ",
                              "knn_max_missing_per_cond) - that block is the single owner of those parameters."),
      normalization       = switch(sw,
                              "maxquant"          = "Intensity summation (log2) + median centering (global) + stoichiometric adjustment",
                              "spectronaut_pivot" = "Intensity summation (log2) + median centering (global) + stoichiometric adjustment",
                              paste0(
                                if (identical(cache$site_summarization, "medpolish"))
                                  "Charge-state summation within peptidoform (linear) + log2 + Tukey median polish over the peptidoform x sample matrix per site"
                                else
                                  "Raw precursor intensity summation per site (log2)",
                                if (has_global)
                                  " + column median centering (PTM and global matrices, independently) + stoichiometric adjustment"
                                else
                                  " + column median centering of the site matrix (PTM only: no stoichiometric adjustment)")),
      stoichiometry       = if (has_global) "adjusted = PTM_log2 - global_log2"
                            else "none (PTM-only mode: no global proteome)",
      adj_centering_delta = cache$adj_centering_delta,
      # Spectronaut cross-run normalization answers and data checks (Step 6); NULL for MaxQuant,
      # and the global answer NULL in PTM-only mode
      ptm_input_normalized     = if (sw %in% c("spectronaut", "spectronaut_pivot")) !isTRUE(cache$ptm_centering_applied %||% TRUE) else NULL,
      global_input_normalized  = if (has_global && sw %in% c("spectronaut", "spectronaut_pivot")) !isTRUE(cache$global_centering_applied %||% TRUE) else NULL,
      ptm_centering_applied    = cache$ptm_centering_applied,
      global_centering_applied = cache$global_centering_applied,
      adj_centering_applied    = cache$adj_centering_applied,
      ptm_centering_max_shift_log2    = cache$ptm_centering_max_shift_log2,
      global_centering_max_shift_log2 = cache$global_centering_max_shift_log2,
      norm_check_ptm_verdict      = rv$dia_ptm_norm_check$ptm$verdict,
      norm_check_ptm_within_sd    = rv$dia_ptm_norm_check$ptm$within_sd,
      norm_check_ptm_n_complete   = rv$dia_ptm_norm_check$ptm$n_complete,
      norm_check_global_verdict   = rv$dia_ptm_norm_check$global$verdict,
      norm_check_global_within_sd = rv$dia_ptm_norm_check$global$within_sd,
      norm_check_global_n_complete = rv$dia_ptm_norm_check$global$n_complete,
      normalization_basis      = cache$normalization_basis,
      normalization_delta_log2 = cache$normalization_delta_log2,
      site_localization   = switch(sw,
                              "maxquant"          = "pre-resolved (PTM.CollapseKey)",
                              "spectronaut_pivot" = "pre-resolved (C: PTM.CollapseKey, Spectronaut pivot)",
                              "consensus (max prob across all samples per candidate position); once a site passes, its intensity is collected from every matched sample, including samples whose own EG.PTMProbabilities value for that row was absent or below threshold"),
      engine              = if (has_global)
                              "limma (eBayes: robust=TRUE throughout; trend=TRUE for the absolute PTM and global protein fits, trend=FALSE for the stoichiometry-adjusted fit)"
                            else
                              "limma (eBayes: robust=TRUE, trend=TRUE) on the absolute PTM matrix",
      pairing             = !is.null(rv$dia_ptm_subject_pairing),
      timestamp           = as.character(Sys.time())
    )

    logger("[Label Free PTM] Wizard completed. rv$ptm_data$source set to 'dia_ptm_upload'.")
  })

  # ---------------------------------------------------------------------------
  # Sidebar summary badge
  # ---------------------------------------------------------------------------
  output$dia_ptm_wiz_sidebar_badges <- renderUI({
    req(isTRUE(rv$dia_ptm_wizard_completed))
    cache <- rv$dia_ptm_master_cache
    tagList(
      tags$div(class = "alert alert-success py-2 px-3 small",
        icon("circle-check"), " Wizard complete.",
        tags$br(),
        if (!is.null(cache)) {
          tagList(
            tags$span(class = "badge bg-success me-1",
                      cache$n_sites_output, " sites"),
            tags$span(class = "badge bg-primary me-1",
                      length(cache$matched_samples), " samples"),
            tags$span(class = "badge bg-secondary me-1",
                      length(cache$cond_groups), " groups")
          )
        }
      )
    )
  })

  # ---------------------------------------------------------------------------
  # EXPORT SUMMARIZATION (.zip with full .rds + .preview.rds)
  # Bundle type = "ptm_summarization", route = "dia_ptm" to distinguish from
  # the MSstatsPTM (TMT) path.  Mirrors the DIA upload export pattern.
  # ---------------------------------------------------------------------------
  output$dia_ptm_upload_export_summ_btn <- downloadHandler(
    filename = function() {
      paste0("DIA_PTM_summarization_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      cache <- rv$dia_ptm_master_cache
      if (is.null(cache)) {
        showNotification("No summarization to export. Complete the wizard first.",
                         type = "warning")
        return()
      }
      # Strip the DE-stage k-NN memo before writing. These are session scratch:
      # they are large, and a bundle restored in a later session could otherwise
      # carry imputed matrices built under a cap that session does not use.
      cache$ptm_abs_matrix_imputed <- NULL
      cache$global_matrix_imputed  <- NULL
      cache$impute_stats           <- NULL
      cache$ptm_impute_cap         <- NULL

      ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
      base <- paste0("DIA_PTM_summarization_", ts)
      sw   <- rv$dia_ptm_software %||% "spectronaut"

      bundle <- list(
        type             = "ptm_summarization",
        route            = "dia_ptm",
        dia_ptm_cache    = cache,
        cond_groups      = rv$dia_ptm_cond_groups,
        conditions       = sort(unique(cache$design_df$condition)),
        software         = sw,
        subject_pairing  = rv$dia_ptm_subject_pairing,
        phospho_filename = rv$dia_ptm_phospho_filename %||% "(unknown)",
        global_filename  = if (ptm_has_global(cache)) rv$dia_ptm_global_filename %||% "(unknown)",
        design_option    = "dia_ptm",
        data_format      = "dia_ptm",
        meta             = list(
          software       = sw,
          pipeline_label = switch(sw,
            "maxquant"          = "LF PTM \u2014 MaxQuant / Perseus (DDA)",
            "spectronaut_pivot" = "LF PTM \u2014 Spectronaut PTM Pivot (DIA)",
            if (ptm_has_global(cache)) "LF PTM \u2014 Spectronaut (DIA)"
            else "LF PTM \u2014 Spectronaut (DIA, PTM only)"),
          n_sites        = cache$n_sites_output,
          n_samples      = length(cache$matched_samples)
        ),
        provenance       = reactiveValuesToList(rv$provenance),
        exported_at      = Sys.time()
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
      logger("[DIA PTM] Summarization exported to ZIP")
    }
  )

  # ---------------------------------------------------------------------------
  # PRERUN LOAD — restore a previously exported DIA PTM bundle.
  # Triggered by server_prerun.R setting rv$prerun_load_dia_ptm_bundle.
  # Replicates the same rv assignments as the wizard-done observer so all four
  # downstream PTM modules (PCA, DE, Trends, Single) sync correctly.
  # ---------------------------------------------------------------------------
  observeEvent(rv$prerun_load_dia_ptm_bundle, {
    bundle <- rv$prerun_load_dia_ptm_bundle
    cache  <- bundle$dia_ptm_cache
    if (is.null(cache)) {
      showNotification("Bundle is missing PTM cache data.", type = "error")
      rv$prerun_load_dia_ptm_bundle <- NULL
      return()
    }
    # The DDA MaxQuant / Perseus and the Spectronaut PTM pivot routes were removed (untested; to be rebuilt). A bundle
    # saved from either is refused rather than restored. Both build their cache with maxquant_ptm_summarize(), so
    # cache$software "maxquant" identifies them; the precursor route's cache never carries it.
    if (isTRUE(bundle$software %in% c("maxquant", "spectronaut_pivot")) || identical(cache$software, "maxquant")) {
      showNotification(paste0("This bundle was made with a PTM input route that has been removed pending ",
                              "re-validation (DDA MaxQuant / Perseus or Spectronaut PTM pivot). It cannot be loaded. ",
                              "Only the Spectronaut precursor-level report is supported."), type = "error", duration = 15)
      rv$prerun_load_dia_ptm_bundle <- NULL
      return()
    }

    withProgress(message = "Loading DIA PTM experiment\u2026", value = 0.2, {
      incProgress(0.5, detail = "Restoring session\u2026")

      rv$dia_ptm_master_cache     <- cache
      rv$dia_ptm_cond_groups      <- bundle$cond_groups
      rv$dia_ptm_software         <- bundle$software %||% "spectronaut"
      rv$dia_ptm_subject_pairing  <- bundle$subject_pairing
      rv$dia_ptm_mode             <- cache$ptm_mode %||% "with_global"
      rv$dia_ptm_phospho_filename <- bundle$phospho_filename %||% "(imported)"
      rv$dia_ptm_global_filename  <- if (ptm_has_global(cache)) bundle$global_filename %||% "(imported)"
      rv$dia_ptm_wizard_completed  <- TRUE
      rv$data_source               <- "dia_ptm"

      # Restore provenance if present
      if (!is.null(bundle$provenance$dia_ptm))
        rv$provenance$dia_ptm <- bundle$provenance$dia_ptm

      # Populate rv$ptm_data for cross-module sync.
      # NULL -> "dia_ptm_upload" transition guarantees observeEvent fires.
      rv$ptm_data$source         <- NULL
      rv$ptm_data$design_option  <- "dia_ptm"
      rv$ptm_data$conditions     <- sort(unique(cache$design_df$condition))
      rv$ptm_data$data_format    <- "dia_ptm"
      rv$ptm_data$n_files        <- if (ptm_has_global(cache)) 2L else 1L
      rv$ptm_data$gene_lookup    <- NULL
      rv$ptm_data$annotation_map <- NULL
      rv$ptm_data$ptm_long       <- NULL
      rv$ptm_data$protein_long   <- NULL

      incProgress(1, detail = "Done")
    })

    # Set source last — triggers all downstream observeEvent(rv$ptm_data$source)
    rv$ptm_data$source            <- "dia_ptm_upload"
    rv$prerun_load_dia_ptm_bundle <- NULL
    showNotification("DIA PTM experiment loaded.", type = "message", duration = 4)
    logger("[DIA PTM] Bundle loaded from Prerun Experiments")
  }, ignoreNULL = TRUE)

}
