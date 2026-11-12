# server_de.R ==================================================================
# Consolidated Differential Expression module (Label-Free DDA).
# Primary analysis modes:
#
#   (a) PD / MaxQuant - ProDA (+ QFeatures for PD peptides):  probabilistic dropout-aware DE
#
# Called from server.R as: server_de(input, output, session, rv, logger)
# ==============================================================================

server_de <- function(input, output, session, rv, logger) {

  # ===========================================================================
  # Reactive storage
  # ===========================================================================
  de_rv <- reactiveValues(
    analysis_mode     = NULL,   # "proda" | "saint" | "welch" | "dia" | "tmt"
    # --- Column selections ---
    ctrl_cols         = NULL,   # selected control abundance columns
    cond_cols         = NULL,   # selected condition abundance columns
    # --- Results ---
    result_df         = NULL,   # scored data.frame
    plot_data         = NULL,   # data.frame ready for plotting
    protein_df        = NULL,   # raw protein data.frame (for filter lookups)
    raw_abund_df      = NULL,   # raw abundance per protein (for CSV download)
    bait_name         = NULL,   # optional bait label
    # --- Plot state ---
    highlight_gene    = NULL,   # gene name to highlight
    render_count      = 0,      # bump to force re-render
    # --- Filters ---
    filter_conditions = list(),
    filter_counter    = 0L,
    # --- Highlight rules ---
    highlight_rules   = list(),
    highlight_counter = 0L,
    # --- DIA-specific ---
    dia_use_lod       = FALSE,   # Perseus LOD imputation toggle
    dia_nmar_result   = NULL,    # data.frame from dia_nmar_fisher()
    dia_nmar_test     = NULL,    # condition labels for NMAR plot
    dia_nmar_ref      = NULL,
    # --- TMT-specific ---
    tmt_one_arm_df    = NULL,    # proteins with adequate coverage in only one of the
                                 # two compared conditions: tested and flagged low
                                 # coverage. Descriptive table (see the tab text).
    tmt_filtered_df   = NULL     # proteins NOT tested: removed by a coverage filter,
                                 # with their per-condition / per-plex coverage
  )
  register_experiment_reset(rv, de_rv)

  # ProDA DE async task --------------------------------------------------------
  proda_de_task <- ExtendedTask$new(function(
      cache, ctrl_label, cond_label, cond_groups, factorial_design,
      batch_for_model, batch_groups, batch_filter_enabled,
      subject_pairing, expect_mnar, ip_sub, igg_conditions) {
    promises::future_promise({
      compute_proda_de(
        cache                = cache,
        ctrl_label           = ctrl_label,
        cond_label           = cond_label,
        cond_groups          = cond_groups,
        factorial_design     = factorial_design,
        batch_for_model      = batch_for_model,
        batch_groups         = batch_groups,
        batch_filter_enabled = batch_filter_enabled,
        subject_pairing      = subject_pairing,
        expect_mnar          = expect_mnar,
        ip_sub               = ip_sub,
        igg_conditions       = igg_conditions
      )
    }, seed = TRUE)
  })
  proda_cancelled <- reactiveVal(FALSE)

  # ProDA DE result observer ---------------------------------------------------
  observeEvent(proda_de_task$result(), {
    hide_async_overlay(session)
    if (isTRUE(proda_cancelled())) {
      proda_cancelled(FALSE)
      return()
    }
    result <- proda_de_task$result()
    for (msg in result$messages %||% character(0)) {
      showNotification(msg, type = "warning", duration = 8)
    }
    if (!isTRUE(result$success)) {
      showNotification(paste0("ProDA DE failed: ", result$error),
                       type = "error", duration = 12)
      return()
    }
    de_rv$result_df     <- result$result_df
    de_rv$plot_data     <- result$plot_data
    de_rv$acc_to_gene   <- result$acc_to_gene
    de_rv$analysis_mode <- "proda"
    de_rv$render_count  <- isolate(de_rv$render_count) + 1L
    rv$provenance$de    <- result$provenance
    showNotification(
      paste0("ProDA DE complete: ", nrow(result$result_df), " proteins (",
             sum(result$result_df$LowEvidence %in% TRUE),
             " low evidence)"),
      type = "message", duration = 5)
  })

  # ProDA cancel observer ------------------------------------------------------
  observeEvent(rv$async_cancel_requested, {
    if (isTRUE(rv$async_cancel_requested) &&
        identical(proda_de_task$status(), "running")) {
      proda_cancelled(TRUE)
      rv$async_cancel_requested <- FALSE
    }
  })

  # SAINTexpress async task ----------------------------------------------------
  saint_de_task <- ExtendedTask$new(function(
      cache, ctrl_label, cond_label, ctrl_cols, cond_cols,
      bait_gene, cond_groups, saint_binary_path, logs_dir) {
    promises::future_promise({
      compute_saint_de(
        cache             = cache,
        ctrl_label        = ctrl_label,
        cond_label        = cond_label,
        ctrl_cols         = ctrl_cols,
        cond_cols         = cond_cols,
        bait_gene         = bait_gene,
        cond_groups       = cond_groups,
        saint_binary_path = saint_binary_path,
        logs_dir          = logs_dir
      )
    }, seed = TRUE)
  })
  saint_cancelled <- reactiveVal(FALSE)

  # SAINTexpress result observer -----------------------------------------------
  observeEvent(saint_de_task$result(), {
    hide_async_overlay(session)
    if (isTRUE(saint_cancelled())) {
      saint_cancelled(FALSE)
      return()
    }
    result <- saint_de_task$result()
    for (msg in result$messages %||% character(0)) {
      showNotification(msg, type = "warning", duration = 8)
    }
    if (!isTRUE(result$success)) {
      showNotification(paste0("SAINTexpress failed: ", result$error),
                       type = "error", duration = 20)
      return()
    }
    de_rv$result_df     <- result$result_df
    de_rv$plot_data     <- result$plot_data
    de_rv$protein_df    <- result$result_df  # enables filtering by BFDR/SaintScore/Confidence
    de_rv$analysis_mode <- "saint"
    de_rv$bait_name     <- result$bait_gene
    de_rv$ctrl_cols     <- result$ctrl_cols
    de_rv$cond_cols     <- result$cond_cols
    de_rv$acc_to_gene   <- result$acc_to_gene
    de_rv$raw_abund_df  <- result$raw_abund_df
    de_rv$render_count  <- isolate(de_rv$render_count) + 1L
    rv$provenance$de    <- result$provenance
    showNotification(
      paste0("SAINTexpress complete: ", nrow(result$result_df), " interactions scored."),
      type = "message", duration = 8)
  })

  # SAINTexpress cancel observer -----------------------------------------------
  observeEvent(rv$async_cancel_requested, {
    if (isTRUE(rv$async_cancel_requested) &&
        identical(saint_de_task$status(), "running")) {
      saint_cancelled(TRUE)
      rv$async_cancel_requested <- FALSE
    }
  })

  # TMT DE async task ----------------------------------------------------------
  tmt_de_task <- ExtendedTask$new(function(
      tmt_master_cache, design_opt, cond_test, cond_ref,
      tmt_analysis_route, acc_to_gene) {
    promises::future_promise({
      compute_tmt_de(
        tmt_master_cache   = tmt_master_cache,
        design_opt         = design_opt,
        cond_test          = cond_test,
        cond_ref           = cond_ref,
        tmt_analysis_route = tmt_analysis_route,
        acc_to_gene        = acc_to_gene
      )
    }, seed = TRUE)
  })
  tmt_de_cancelled <- reactiveVal(FALSE)

  # TMT DE result observer -----------------------------------------------------
  observeEvent(tmt_de_task$result(), {
    shinyjs::enable("de_run_btn")
    hide_async_overlay(session)
    if (isTRUE(tmt_de_cancelled())) {
      tmt_de_cancelled(FALSE)
      return()
    }
    res <- tmt_de_task$result()
    for (msg in res$messages %||% character(0)) {
      showNotification(msg, type = "warning", duration = 8)
    }
    if (!isTRUE(res$success)) {
      showNotification(paste0("TMT DE failed: ", res$error), type = "error", duration = 15)
      return()
    }
    de_rv$result_df      <- res$result_df
    de_rv$plot_data      <- res$plot_data
    de_rv$analysis_mode  <- res$analysis_mode
    de_rv$cond_test      <- res$cond_test
    de_rv$cond_ref       <- res$cond_ref
    de_rv$tmt_one_arm_df <- res$one_arm_df
    de_rv$tmt_filtered_df <- res$filtered_df
    de_rv$render_count   <- isolate(de_rv$render_count) + 1L
    rv$provenance$tmt_volcano <- res$provenance_tmt_volcano
    rv$provenance$de          <- res$provenance_de
  })

  # TMT DE cancel observer -----------------------------------------------------
  observeEvent(rv$async_cancel_requested, {
    if (isTRUE(rv$async_cancel_requested) &&
        identical(tmt_de_task$status(), "running")) {
      tmt_de_cancelled(TRUE)
      rv$async_cancel_requested <- FALSE
    }
  })

  # ===========================================================================
  # Methods Text Module  -" renders provenance-based publication text
  # ===========================================================================
  methodsTextServer("de_methods", provenance = rv$provenance, section = "de")

  # ===========================================================================
  # DIA  -" expose reactive flags for conditionalPanel
  # ===========================================================================
  output$de_is_dia <- reactive({
    identical(rv$data_source, "dia")
  })
  outputOptions(output, "de_is_dia", suspendWhenHidden = FALSE)

  output$de_lod_active <- reactive({
    isTRUE(de_rv$dia_use_lod)
  })
  outputOptions(output, "de_lod_active", suspendWhenHidden = FALSE)

  # Show / hide NMAR tab depending on whether DIA data is loaded
  observe({
    is_dia <- identical(rv$data_source, "dia")
    if (is_dia) {
      shinyjs::runjs('var el = document.querySelector("[data-value=\'de_nmar_tab\']"); if (el) el.parentElement.style.display = "";')
    } else {
      shinyjs::runjs('var el = document.querySelector("[data-value=\'de_nmar_tab\']"); if (el) el.parentElement.style.display = "none";')
    }
  })

  # Show / hide the TMT one-condition-missing tab. Same mechanism as the NMAR
  # tab above, keyed on the TMT data source instead of DIA.
  observe({
    is_tmt <- identical(rv$data_source, "tmt")
    if (is_tmt) {
      shinyjs::runjs('var el = document.querySelector("[data-value=\'de_tmt_onearm_tab\']"); if (el) el.parentElement.style.display = "";')
    } else {
      shinyjs::runjs('var el = document.querySelector("[data-value=\'de_tmt_onearm_tab\']"); if (el) el.parentElement.style.display = "none";')
    }
  })

  # ---------------------------------------------------------------------------
  # NMAR tab content  -" single renderUI handles all states
  # ---------------------------------------------------------------------------
  output$de_nmar_content_ui <- renderUI({
    # State 1: LOD imputation is active  -" analysis not applicable
    if (isTRUE(de_rv$dia_use_lod)) {
      return(tags$div(
        class = "alert alert-warning d-flex gap-2 align-items-start my-3",
        style = "border-left: 4px solid #fd7e14;",
        icon("wand-magic-sparkles", class = "mt-1 flex-shrink-0"),
        tags$div(
          tags$b("Perseus LOD Imputation is active."),
          tags$br(),
          "When LOD imputation is enabled, missing values are filled for all proteins ",
          "before the DE test  -  so there are no \"excluded\" proteins to analyse here. ",
          "The presence/absence distinction that NMAR analysis relies on is absorbed ",
          "into the imputed fold changes in the main DE result.",
          tags$br(), tags$br(),
          tags$span(class = "text-muted small",
            "Deactivate Perseus LOD Mode and re-run DE to use this analysis."
          )
        )
      ))
    }

    df <- de_rv$dia_nmar_result

    # State 2: DIA loaded but DE not yet run
    if (is.null(df)) {
      return(tags$div(
        class = "alert alert-secondary my-3",
        icon("clock"), " Run DIA DE Analysis to populate this tab."
      ))
    }

    # State 3: Results available  -" render info box + plotly + summary
    cond_test      <- de_rv$dia_nmar_test %||% "Test"
    cond_ref       <- de_rv$dia_nmar_ref  %||% "Ref"
    n_total        <- nrow(df)
    n_sig          <- sum(df$fisher_adj_p < 0.05, na.rm = TRUE)
    n_drop_total   <- sum(df$imp_class == "dropout_one_cond", na.rm = TRUE)  # all triangles
    n_drop_sig     <- sum(df$imp_class == "dropout_one_cond" & df$fisher_adj_p < 0.05, na.rm = TRUE)

    # Fisher's exact test cannot reach BH-adjusted significance at small
    # replicate counts. When that is the case the p-values are suppressed and
    # the proteins are presented as a ranked list instead.
    has_power  <- !("fisher_has_power" %in% names(df)) || isTRUE(df$fisher_has_power[1])
    min_ach_p  <- if ("min_achievable_p" %in% names(df)) df$min_achievable_p[1] else NA_real_
    n_samp_t   <- if ("n_total_test" %in% names(df)) df$n_total_test[1] else NA_integer_
    n_samp_r   <- if ("n_total_ref"  %in% names(df)) df$n_total_ref[1]  else NA_integer_

    tagList(
      tags$div(
        class = "alert alert-info d-flex gap-2 align-items-start my-3",
        style = "border-left: 4px solid #0d6efd;",
        icon("circle-info", class = "mt-1 flex-shrink-0"),
        tags$div(
          tags$b("Presence/Absence Analysis \u2014 Proteins Excluded from Statistical Testing"),
          tags$br(),
          "These proteins were not testable by limma because they lack sufficient ",
          "detection across both groups. Rather than discarding them, this analysis ",
          "asks whether their ",
          tags$b("absence is itself non-random"),
          " \u2014 which can reflect biological regulation rather than stochastic dropout. ",
          "Significance (Fisher\u2019s exact, BH-corrected) reflects enriched detection ",
          "in one condition; fold change is estimated from LOD-imputed values to ",
          "give a directional signal even for proteins completely absent in one group.",
          tags$br(), tags$br(),
          tags$b("Benjamini-Hochberg correction for this table is computed within the ",
                 "excluded protein set only, independently of the main limma results. ",
                 "Do not combine the two tables into a single FDR-controlled list."),
          tags$br(), tags$br(),
          tags$span(class = "text-muted small",
            icon("triangle-exclamation"),
            " Triangle markers: detected exclusively in one condition. ",
            "Treat these fold change estimates as directional indicators, not precise measurements."
          )
        )
      ),
      if (!has_power)
        tags$div(
          class = "alert alert-warning d-flex gap-2 align-items-start my-3",
          style = "border-left: 4px solid #fd7e14;",
          icon("circle-exclamation", class = "mt-1 flex-shrink-0"),
          tags$div(
            paste0("This design (", n_samp_t, " vs ", n_samp_r,
                   ") cannot reach significance with Fisher's exact test ",
                   "(minimum achievable p = ", signif(min_ach_p, 3),
                   "). Proteins are ranked by detection difference.")
          )
        )
      else NULL,
      plotlyOutput("de_nmar_plotly", height = "600px"),
      br(),
      tags$div(class = "alert alert-secondary py-2 px-3 small",
        icon("table"), " ",
        tags$b(n_total), " proteins excluded from DE tested for non-random missingness. ",
        tags$b(n_drop_total), " detected exclusively in one condition (\u25b2 triangles in plot). ",
        if (has_power)
          tagList(tags$b(n_sig), " significant overall (BH-adj Fisher p\u00a0<\u00a00.05).")
        else
          "Fisher p-values are suppressed for this design; proteins are ranked by detection difference.",
        if (has_power && n_sig == 0L && n_total > 0L)
          tags$span(class = "text-muted d-block mt-1",
            icon("circle-info"), " ",
            "Fisher\u2019s exact test has limited power at typical DIA group sizes (n\u202f<\u202f5 per condition): ",
            "even perfect absence in one group (e.g. 3 vs. 0 detections) gives p\u202f=\u202f0.10 uncorrected. ",
            "The plot still ranks proteins directionally by their absence-pattern strength."
          )
      )
    )
  })

  output$de_nmar_plotly <- renderPlotly({
    df <- de_rv$dia_nmar_result
    req(!is.null(df), nrow(df) > 0, !isTRUE(de_rv$dia_use_lod))

    cond_test <- de_rv$dia_nmar_test %||% "Test"
    cond_ref  <- de_rv$dia_nmar_ref  %||% "Ref"

    # When the design cannot reach BH-adjusted significance with Fisher's exact
    # test the p-values are not interpretable, so they are suppressed and the
    # Y-axis becomes the absolute detection difference (a ranked list).
    .has_power <- !("fisher_has_power" %in% names(df)) || isTRUE(df$fisher_has_power[1])
    df$det_diff    <- abs(df$pct_test - df$pct_ref)
    if (.has_power) {
      df$neg_log_p <- -log10(pmax(df$fisher_adj_p, 1e-300))
      df$sig       <- df$fisher_adj_p < 0.05
    } else {
      df$neg_log_p <- df$det_diff
      df$sig       <- rep(FALSE, nrow(df))
    }
    fc_finite      <- df$log2fc[is.finite(df$log2fc)]
    fc_cap         <- if (length(fc_finite) > 0) max(abs(fc_finite), na.rm = TRUE) + 0.5 else 5
    df$log2fc_disp <- pmin(pmax(df$log2fc, -fc_cap), fc_cap)
    df$colour      <- if (.has_power)
                        ifelse(!df$sig, "#adb5bd",
                          ifelse(df$log2fc > 0, "#e74c3c", "#3498db"))
                      else
                        ifelse(df$log2fc > 0, "#e74c3c", "#3498db")
    df$symbol      <- ifelse(df$imp_class == "dropout_one_cond", "triangle-up", "circle")
    df$hover <- paste0(
      "<b>", df$Gene, "</b> (", df$Protein, ")<br>",
      "Log2FC (LOD-est): ", round(df$log2fc, 2), "<br>",
      cond_test, ": ", df$n_obs_test, "/", df$n_total_test, " (", df$pct_test, "%) detected<br>",
      cond_ref,  ": ", df$n_obs_ref,  "/", df$n_total_ref,  " (", df$pct_ref,  "%) detected<br>",
      if (.has_power)
        paste0("Fisher raw p: ", signif(df$fisher_p, 3), "<br>",
               "Fisher adj.p (BH): ", signif(df$fisher_adj_p, 3), "<br>")
      else
        paste0("Detection difference: ", round(df$det_diff, 1),
               " percentage points<br>",
               "Fisher p: suppressed (design has no resolution)<br>"),
      "Class: ", df$imp_class
    )
    # Ensure the y-axis always extends past the significance threshold so the
    # dashed line and annotation are visible even when all adj.p ≈ 1.
    y_thresh   <- -log10(0.05)          # ≈ 1.30
    y_max_plot <- if (.has_power)
      max(max(df$neg_log_p, na.rm = TRUE) * 1.15, y_thresh * 1.6)
    else
      max(max(df$neg_log_p, na.rm = TRUE) * 1.15, 1)
    .nmar_p <- plotly::plot_ly(df,
      x = ~log2fc_disp, y = ~neg_log_p,
      type = "scatter", mode = "markers",
      marker = list(size = 7, opacity = 0.75, color = ~colour, symbol = ~symbol,
                    line = list(width = 0.5, color = "white")),
      text = ~hover, hoverinfo = "text"
    ) %>%
      plotly::layout(
        title = list(
          text = paste0("NMAR: ", cond_test, " vs ", cond_ref,
                        "  <span style='font-size:12px;color:#6c757d;'>(proteins excluded from limma)</span>"),
          font = list(size = 14)
        ),
        xaxis = list(title = paste0("Log2 FC (LOD-imputed): ", cond_test, " / ", cond_ref),
                     zeroline = TRUE, zerolinecolor = "#ccc"),
        yaxis = list(title = if (.has_power) "\u2212log\u2081\u2080 (BH-adj Fisher p)"
                             else "|detection difference| (percentage points)",
                     range = c(0, y_max_plot)),
        shapes = if (.has_power) list(
          list(type = "line", x0 = -fc_cap, x1 = fc_cap,
               y0 = y_thresh, y1 = y_thresh,
               line = list(dash = "dot", color = "#e74c3c", width = 1))
        ) else list()
      )
    if (.has_power)
      .nmar_p <- plotly::add_annotations(
        .nmar_p,
        x = fc_cap, y = y_thresh, xref = "x", yref = "y",
        text = "adj.p = 0.05", showarrow = FALSE,
        font = list(size = 10, color = "#e74c3c"), xanchor = "right"
      )
    .nmar_p
  })

  # ---------------------------------------------------------------------------
  # TMT "Not Tested & Low Coverage" tab  -  descriptive tables, no inference
  #   1. proteins NOT tested: removed by a coverage filter (plex coverage at
  #      summarization on the multi-plex protein route; the missingness filter
  #      before DE on every route), with their coverage and the reason
  #   2. proteins tested with adequate coverage in only one of the two compared
  #      conditions (flagged LowCoverage in the results table)
  # ---------------------------------------------------------------------------
  output$de_tmt_onearm_content_ui <- renderUI({
    # TMT loaded but DE not yet run
    if (is.null(de_rv$result_df) || !identical(de_rv$analysis_mode, "tmt")) {
      return(tags$div(
        class = "alert alert-secondary my-3",
        icon("clock"), " Run TMT DE Analysis to populate this tab."
      ))
    }
    cond_test <- de_rv$cond_test %||% "Test"
    cond_ref  <- de_rv$cond_ref  %||% "Ref"
    fdf <- de_rv$tmt_filtered_df
    odf <- de_rv$tmt_one_arm_df
    n_f <- if (is.null(fdf)) 0L else nrow(fdf)
    n_o <- if (is.null(odf)) 0L else nrow(odf)
    n_f_pc <- if (n_f > 0L) sum(fdf$Stage == "Summarization (plex coverage)") else 0L
    n_f_nv <- if (n_f > 0L) sum(fdf$Stage == TMT_STAGE_NO_VALUES) else 0L       # never quantified
    n_f_mf <- n_f - n_f_pc - n_f_nv                                              # missingness filter
    .dopt  <- as.integer(rv$tmt_limma_cache$design_option %||% rv$tmt_design_option %||% 1L)
    odf_disp <- tmt_low_coverage_display(odf, de_rv$result_df)
    n_o_yes  <- if (n_o > 0L) sum(startsWith(odf_disp[["In DE table & volcano"]], "Yes")) else 0L

    not_tested_ui <- if (n_f == 0L) {
      tags$div(class = "alert alert-success my-3", icon("circle-check"),
        tags$b(" No proteins were removed by the coverage filters."),
        " Every quantified protein entered the model.")
    } else tagList(
      tags$div(
        class = "alert alert-info d-flex gap-2 align-items-start my-3",
        style = "border-left: 4px solid #0d6efd;",
        icon("circle-info", class = "mt-1 flex-shrink-0"),
        tags$div(
          tags$b(paste0(n_f, " protein(s) were not tested.")),
          tags$br(),
          "They are not in the DE results or the volcano plot: they have no p-value and ",
          "are not part of the FDR correction. They are listed here, with the channels ",
          "they were observed in, so you can judge them yourself.",
          tags$br(), tags$br(),
          tags$b("Why (Stage column). "),
          if (n_f_pc > 0L) tagList(tags$b("Summarization (plex coverage): "),
            "on a multi-plex design, a protein must be quantified in at least half of the ",
            "plexes (", n_f_pc, " removed). In isobaric labelling a protein is measured ",
            "in all channels of a plex or in none, so a protein missing from most plexes ",
            "rests on a minority of the experiment. ", tags$br()),
          if (n_f_mf > 0L) tagList(tags$b("Before DE (missingness filter): "),
            "at least 50% of a condition's channels must be observed in at least one ",
            "condition", if (.dopt > 1L) " (on a balanced multi-plex design, about half of the plexes)" else "",
            " (", n_f_mf, " removed). ", tags$br()),
          if (n_f_nv > 0L) tagList(tags$b(paste0(TMT_STAGE_NO_VALUES, ": ")),
            "in the uploaded file but without a reporter value in any channel - identified, ",
            "never quantified (", n_f_nv, "; listed last). Proteome Discoverer writes such rows ",
            "when it excludes a protein's PSMs from quantification (co-isolation or S/N ",
            "thresholds, missing reporter ions), or when its abundances are all zero. ",
            "Nothing was removed by OSMSS here: there was nothing to test."),
          tags$br(), tags$br(),
          tags$span(class = "text-muted small",
            tags$b("Columns: "), "Plexes_quantified = plexes with at least one observed channel; ",
            "<condition> / <plex> (channels observed) = observed of annotated channels. ",
            "No statistical test is applied to this table. Treat these proteins as ",
            "candidates to inspect, not as results.")
        )
      ),
      DT::DTOutput("de_tmt_filtered_table"),
      br(),
      downloadButton("de_tmt_filtered_download", "Download CSV (not tested)",
                     class = "btn-outline-secondary btn-sm", icon = icon("download"))
    )

    low_cov_ui <- if (n_o == 0L) {
      tags$div(class = "alert alert-success my-3", icon("circle-check"),
        tags$b(" No tested protein has one-sided coverage."),
        " Every tested protein had at least 50% of its channels observed in both ",
        tags$b(cond_test), " and ", tags$b(cond_ref), ".")
    } else tagList(
      tags$div(
        class = "alert alert-info d-flex gap-2 align-items-start my-3",
        style = "border-left: 4px solid #0d6efd;",
        icon("circle-info", class = "mt-1 flex-shrink-0"),
        tags$div(
          tags$b(paste0(n_o, " protein(s) have at least 50% of their channels observed in only one of ",
                        cond_test, " and ", cond_ref, ".")),
          tags$br(),
          tags$b(paste0(n_o_yes, " are in the DE results and the volcano plot")),
          " (column ", tags$i("In DE table & volcano"), " = Yes): they have at least one value ",
          "in both conditions, were tested normally with a real p-value, and are part of the ",
          "FDR correction; their thinner evidence shows as a larger standard error.",
          if (n_o - n_o_yes > 0L) tagList(" ", tags$b(paste0(n_o - n_o_yes, " are not")),
            " (= No): they have no value at all in one of the two conditions, so the ",
            "difference cannot be estimated. They appear in the DE table without a p-value ",
            "and are not plotted."),
          tags$br(), tags$br(),
          tags$span(class = "text-muted small",
            tags$b("Columns: "), "channels observed = channels with a value / channels of that ",
            "condition; Adequately covered in = the condition with at least 50% observed; ",
            "Zero values in input = gaps where the file held a zero (reporter ion present but at ",
            "the noise floor) rather than a blank; Peptides = peptides behind the protein ",
            "(1 = single-peptide protein). ",
            "In isobaric labelling every channel is quantified from the same fragment ",
            "spectrum, so a gap confined to one condition usually means the reporter ion ",
            "fell to the noise floor", if (.dopt > 1L) ", or, across plexes, that the protein was not identified in some plexes" else "",
            ". It can also be mechanical: zeros written by Proteome Discoverer, peptides ",
            "lost asymmetrically to the co-isolation filter, or a protein resting on a ",
            "single peptide.")
        )
      ),
      DT::DTOutput("de_tmt_onearm_table"),
      br(),
      downloadButton("de_tmt_onearm_download", "Download CSV (low coverage)",
                     class = "btn-outline-secondary btn-sm", icon = icon("download"))
    )

    tagList(
      tags$h5(class = "mt-3", icon("filter-circle-xmark"), " Not tested (removed by coverage filters)"),
      not_tested_ui,
      tags$hr(),
      tags$h5(icon("circle-half-stroke"), paste0(" Tested with low coverage (", cond_test, " vs ", cond_ref, ")")),
      low_cov_ui
    )
  })

  output$de_tmt_filtered_table <- DT::renderDT({
    df <- de_rv$tmt_filtered_df
    req(!is.null(df), nrow(df) > 0L)
    DT::datatable(df, rownames = FALSE, selection = "none", filter = "top",
                  options = list(pageLength = 25, scrollX = TRUE, order = list(list(0, "asc"))))
  })

  output$de_tmt_filtered_download <- downloadHandler(
    filename = function() paste0("TMT_not_tested_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      # The reactive value persists across a switch to a non-TMT data source.
      df <- if (identical(de_rv$analysis_mode, "tmt")) de_rv$tmt_filtered_df else NULL
      utils::write.csv(if (is.null(df)) data.frame() else df, file, row.names = FALSE)
    }
  )

  output$de_tmt_onearm_table <- DT::renderDT({
    df <- tmt_low_coverage_display(de_rv$tmt_one_arm_df, de_rv$result_df)
    req(!is.null(df), nrow(df) > 0L)
    DT::datatable(
      df,
      rownames  = FALSE,
      selection = "none",
      options   = list(pageLength = 25, scrollX = TRUE,
                       order = list(list(0, "asc")))
    )
  })

  output$de_tmt_onearm_download <- downloadHandler(
    filename = function() {
      paste0("TMT_low_coverage_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      # The reactive value persists across a switch to a non-TMT data source,
      # and this handler is reachable independently of the tab's visibility.
      df <- if (identical(de_rv$analysis_mode, "tmt"))
              tmt_low_coverage_display(de_rv$tmt_one_arm_df, de_rv$result_df) else NULL
      utils::write.csv(if (is.null(df)) data.frame() else df,
                       file, row.names = FALSE)
    }
  )

  # Perseus LOD toggle
  observeEvent(input$de_lod_toggle, {
    de_rv$dia_use_lod <- !isTRUE(de_rv$dia_use_lod)
    if (de_rv$dia_use_lod) {
      updateActionButton(session, "de_lod_toggle",
                         label = "Perseus LOD Mode Active",
                         icon  = icon("wand-magic-sparkles"))
      shinyjs::runjs('document.getElementById("de_lod_toggle").classList.remove("btn-secondary"); document.getElementById("de_lod_toggle").classList.add("btn-success");')
    } else {
      updateActionButton(session, "de_lod_toggle",
                         label = "Activate Perseus LOD Mode",
                         icon  = icon("wand-magic-sparkles"))
      shinyjs::runjs('document.getElementById("de_lod_toggle").classList.remove("btn-success"); document.getElementById("de_lod_toggle").classList.add("btn-secondary");')
      updateCheckboxInput(session, "de_lod_ack", value = FALSE)
    }
  })

  # ===========================================================================
  # File status
  # ===========================================================================
  # Enable/disable Run button based on file upload
  observe({
    src       <- rv$data_source
    dia_ready <- identical(src, "dia") && !is.null(rv$dia_master_cache)
    tmt_ready <- identical(src, "tmt") && !is.null(rv$tmt_limma_cache)
    mq_ready       <- identical(rv$lf_input_mode, "maxquant") && !is.null(rv$lf_master_cache)
    # lf_cache_ready covers all LF paths (peptide, SAINT, ProDA) loaded via
    # the normal summarize button or via Prerun Experiments import.
    lf_cache_ready <- !is.null(rv$lf_master_cache)
    lod_blocks <- identical(src, "dia") && isTRUE(de_rv$dia_use_lod) && !isTRUE(input$de_lod_ack)
    if ((dia_ready || tmt_ready || !is.null(rv$pd_filename) ||
         mq_ready || lf_cache_ready) && !lod_blocks) {
      shinyjs::enable("de_run_btn")
    } else {
      shinyjs::disable("de_run_btn")
    }
  })

  # ===========================================================================
  # SUMMARIZATION CACHE  -" Summarize / Export / Import
  # ===========================================================================
  output$upload_export_summarization <- downloadHandler(
    filename = function() {
      paste0("LF_Summarization_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
    },
    content = function(file) {
      cache <- rv$lf_master_cache; meta <- rv$lf_master_cache_meta
      if (is.null(cache)) {
        showNotification("No summarization to export.", type = "warning")
        return()
      }
      ts   <- format(Sys.time(), "%Y%m%d_%H%M")
      base <- paste0("LF_Summarization_", ts)
      bundle <- list(
        type            = "lf_summarization",
        version         = 3L,
        cache           = cache,
        meta            = meta,
        cond_groups     = rv$lf_cond_groups,
        batch_groups    = rv$lf_batch_groups,
        batch_corrected = isTRUE(rv$lf_batch_corrected),
        input_mode      = rv$lf_input_mode %||% "peptide",
        # Recorded so an imported bundle reports the setting it was built with,
        # rather than whatever the current default happens to be.
        shared_peptides = cache$shared_peptide_stats$shared_peptide_handling %||%
                          rv$lf_shared_peptides
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
      logger("DE: Summarization exported")
    }
  )

  observe({
    if (is.null(rv$lf_master_cache)) {
      shinyjs::disable("de_export_summarization")
    } else {
      shinyjs::enable("de_export_summarization")
    }
  })

  # ===========================================================================
  # Batch Assignment  -" visualization only (PCA & Heatmap removeBatchEffect)
  # ===========================================================================

  # Status indicator for batch assignment
  output$upload_batch_status <- renderUI({
    bg <- rv$lf_batch_groups
    if (is.null(bg) || length(bg) == 0) {
      tags$p(class = "small text-muted mb-0", icon("circle-info"),
             " No batches assigned.")
    } else {
      n <- length(bg)
      tags$div(class = "alert alert-info py-1 px-2 small mb-0",
        icon("circle-check"),
        tags$b(paste0(n, " batch", if (n > 1) "es" else "", " assigned.")),
        tags$br(),
        paste(vapply(names(bg), function(nm)
          paste0(nm, ": ", length(bg[[nm]]), " runs"), character(1)),
          collapse = " | ")
      )
    }
  })

  # Expose whether batch groups are active (for conditionalPanel).
  # Hidden during Background Z-Test (welch) mode  -" batch is not applied there.
  output$upload_batch_active <- reactive({
    !is.null(rv$lf_batch_groups) && length(rv$lf_batch_groups) > 0 &&
      !identical(de_rv$analysis_mode, "welch")
  })
  outputOptions(output, "upload_batch_active", suspendWhenHidden = FALSE)

  # Clear button
  observeEvent(input$upload_batch_clear, {
    rv$lf_batch_groups    <- NULL
    rv$lf_batch_corrected <- FALSE
    logger("DE: Batch assignment cleared")
    showNotification("Batch assignment cleared.", type = "message", duration = 3)
  })

  # ── Prerun Experiments: LF bundle load trigger ──────────────────────────────
  observeEvent(rv$prerun_load_lf_bundle, {
    bundle <- rv$prerun_load_lf_bundle
    req(!is.null(bundle))
    meta <- bundle$meta
    rv$lf_import_in_progress <- TRUE
    session$onFlushed(function() { rv$lf_import_in_progress <- FALSE }, once = TRUE)
    withProgress(message = "Loading experiment...", value = 0, {
      rv$lf_input_mode <- bundle$input_mode %||% "peptide"
      # Bundles predating the shared-peptide setting carry no key. Those runs
      # kept every peptide, so they restore as "keep_all" - defaulting them to
      # the new default would misreport what they actually did.
      rv$lf_shared_peptides <- bundle$shared_peptides %||% "keep_all"
      incProgress(0.4, detail = "Restoring pipeline settings...")
      if (!is.null(meta$pipeline))    rv$lf_pipeline     <- meta$pipeline
      if (!is.null(meta$ip_sub))      rv$lf_ip_sub       <- meta$ip_sub
      if (!is.null(meta$skip_impute)) rv$lf_skip_impute  <- meta$skip_impute
      if (!is.null(meta$bait))        rv$lf_bait_protein <- meta$bait
      if (!is.null(meta$de_engine))   rv$lf_de_engine    <- meta$de_engine
      if (!is.null(meta$groups) && length(meta$groups) > 0)
        rv$lf_cond_groups <- meta$groups
      if (!is.null(bundle$cond_groups) && length(bundle$cond_groups) > 0)
        rv$lf_cond_groups <- bundle$cond_groups
      incProgress(0.7, detail = "Loading cache...")
      rv$lf_master_cache      <- bundle$cache
      rv$lf_master_cache_meta <- meta
      rv$lf_cache_stale       <- FALSE
      if (!is.null(bundle$batch_groups))    rv$lf_batch_groups    <- bundle$batch_groups
      if (!is.null(bundle$batch_corrected)) rv$lf_batch_corrected <- bundle$batch_corrected
      rv$wizard_completed <- TRUE
      rv$data_source      <- "lf_dda"
      incProgress(1, detail = "Done!")
    })
    rv$prerun_load_lf_bundle <- NULL
    label <- meta$pipeline_label %||% meta$pipeline %||% "LF"
    showNotification(paste0("Experiment loaded: ", label), type = "message", duration = 4)
    logger("DE: LF bundle loaded from Prerun Experiments")

    # SAINTexpress bundles exported before bait was persisted to meta will
    # have rv$lf_bait_protein == NULL.  Prompt the user to supply it now so
    # the SAINT run can proceed without re-running the wizard.
    if (identical(rv$lf_de_engine, "saint") &&
        (is.null(rv$lf_bait_protein) || !nzchar(rv$lf_bait_protein %||% ""))) {
      gene_choices <- if (!is.null(rv$lf_master_cache$acc_to_gene)) {
        sort(unique(rv$lf_master_cache$acc_to_gene[
          nchar(rv$lf_master_cache$acc_to_gene) > 0]))
      } else character(0)

      showModal(modalDialog(
        title = tagList(icon("bullseye"), " Enter Bait Protein"),
        tags$p(class = "small",
          "This experiment uses SAINTexpress but the bait protein name was not ",
          "saved in the bundle. Please enter it now to proceed."),
        if (length(gene_choices) > 0)
          selectizeInput("prerun_bait_entry", "Bait Gene Symbol:",
                         choices  = gene_choices,
                         selected = NULL,
                         width    = "100%",
                         options  = list(placeholder = "Type to search -¦"))
        else
          textInput("prerun_bait_entry", "Bait Gene Symbol:", width = "100%"),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("prerun_bait_confirm", "Confirm",
                       class = "btn-success", icon = icon("check"))
        ),
        size = "s", easyClose = FALSE
      ))
    }
  }, ignoreNULL = TRUE)

observeEvent(input$prerun_bait_confirm, {
    bait <- trimws(input$prerun_bait_entry %||% "")
    if (!nzchar(bait)) {
      showNotification("Please enter a bait gene name.", type = "warning")
      return()
    }
    rv$lf_bait_protein <- bait
    if (!is.null(rv$lf_master_cache_meta))
      rv$lf_master_cache_meta$bait <- bait
      
    # 1. Close the Bait Selection Modal
    removeModal()
    
    showNotification(paste0("Bait set to: ", bait), type = "message", duration = 3)
    logger(paste0("DE: Bait protein set from import prompt: ", bait))
    
    # =========================================================================
    # 2. SEQUENTIAL LAUNCHPAD (Fires immediately after bait is saved)
    # =========================================================================
    showModal(modalDialog(
      title = tagList(icon("circle-check", class = "text-success"), " Bait Saved & Data Ready!"),
      tags$p("Your IP experiment is fully configured. Where would you like to begin?"),
      tags$div(class = "d-grid gap-2",
        
        tags$button(class = "btn btn-primary btn-lg text-start",
          icon("braille"), " Explore Data: PCA Analysis",
          onclick = "
            $('.modal').modal('hide'); 
            $('[data-tab=\"PCA Analysis\"]').click();
            $('#submenu_discoverer, #submenu_protein, #submenu_explore_data').addClass('open').prev('.nav-group-btn').addClass('open');
          "
        ),
        
        tags$button(class = "btn btn-info btn-lg text-start text-white",
          icon("chart-area"), " Analyze: Differential Expression",
          onclick = "
            $('.modal').modal('hide'); 
            $('[data-tab=\"Differential Expression\"]').click();
            $('#submenu_discoverer, #submenu_protein, #submenu_explore_data').addClass('open').prev('.nav-group-btn').addClass('open');
          "
        ),
        
        tags$button(class = "btn btn-outline-secondary text-start mt-2",
          icon("file-lines"), " View Data Summary",
          `data-bs-dismiss` = "modal"
        )
      ),
      easyClose = TRUE, footer = NULL
    ))
  })

  # ===========================================================================
  # Auto-select recommended columns when CSV arrives (mode b)
  # ===========================================================================
  observeEvent(rv$pd_protein_headers, {
    req(rv$pd_protein_headers)
    de_rv$plot_data   <- NULL
    de_rv$result_df   <- NULL
    de_rv$protein_df  <- NULL
    logger("DE: Shared CSV data updated, columns auto-selected")
  })

  # ===========================================================================
  # Condition Group Manager (DE tab  -" shares rv$lf_cond_groups with PCA/HM/TCL)
  # ===========================================================================

  # ── Helper: read ctrl/cond labels from standard pickers ──
  .get_de_ctrl_cond <- function() {
    list(ctrl = isolate(input$de_ctrl_group) %||% "",
         cond = isolate(input$de_cond_group) %||% "")
  }

  # Returns the effective condition groups for DE pickers.
  .de_effective_groups <- function() {
    src <- rv$data_source
    if (identical(src, "dia")) {
      grps <- rv$dia_cond_groups
      if (length(grps) == 0) return(list())
      return(grps)
    }
    if (identical(src, "tmt")) {
      conds <- rv$tmt_conditions
      if (is.null(conds) || length(conds) == 0) return(list())
      grps <- as.list(setNames(conds, conds))
      return(grps)
    }
    # LF DDA
    grps <- rv$lf_cond_groups
    if (!is.null(rv$lf_master_cache$processed$ProteinLevelData)) {
      cache_grps <- unique(as.character(
        rv$lf_master_cache$processed$ProteinLevelData$GROUP))
      grps <- grps[names(grps) %in% cache_grps]
    }
    grps
  }

  # A previous selection is kept only if it names a group of the current
  # dataset; otherwise (new dataset) the default is used. Control defaults to
  # the first group, Condition to the first group other than the control.
  .de_picker_selection <- function(prev, nms, exclude = NULL) {
    if (!is.null(prev) && prev %in% nms && !identical(prev, exclude)) return(prev)
    rest <- setdiff(nms, exclude)
    if (length(rest) > 0) rest[1] else nms[1]
  }

  output$de_ctrl_group_picker <- renderUI({
    grps <- .de_effective_groups()
    if (length(grps) == 0) return(NULL)
    nms <- names(grps)
    selectInput("de_ctrl_group", "Control Group:",
                choices = nms,
                selected = .de_picker_selection(isolate(input$de_ctrl_group), nms),
                width = "100%")
  })

  output$de_cond_group_picker <- renderUI({
    grps <- .de_effective_groups()
    if (length(grps) == 0) return(NULL)
    nms  <- names(grps)
    ctrl <- .de_picker_selection(isolate(input$de_ctrl_group), nms)
    selectInput("de_cond_group", "Condition Group:",
                choices = nms,
                selected = .de_picker_selection(isolate(input$de_cond_group), nms,
                                                exclude = ctrl),
                width = "100%")
  })

  # ===========================================================================
  # Helper: Parse protein-level data from shared CSV
  # ===========================================================================
  parse_protein_data <- function() {
    header_line   <- rv$pd_raw_lines[1]
    protein_lines <- rv$pd_raw_lines[rv$pd_protein_line_idx]
    text_block    <- paste(c(header_line, protein_lines), collapse = "\n")
    protein_df <- tryCatch(
      read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, na.strings = c("", "NA")),
      error = function(e) { logger(paste0("!!! Protein CSV parse error: ", e$message)); NULL }
    )
    if (is.null(protein_df)) return(NULL)
    if (!isTRUE(rv$pd_has_gn) && !is.null(rv$pd_gene_col) &&
        rv$pd_gene_col %in% colnames(protein_df)) {
      gene_names <- as.character(protein_df[[rv$pd_gene_col]])
      gene_names[is.na(gene_names)] <- ""
    } else {
      descriptions <- protein_df[["Description"]]
      gene_names <- stringr::str_extract(descriptions, "GN=([^ ]+)")
      gene_names <- stringr::str_replace(gene_names, "^GN=", "")
      gene_names[is.na(gene_names)] <- ""
    }
    accessions <- as.character(protein_df[["Accession"]])
    list(protein_df = protein_df, gene_names = gene_names, accessions = accessions)
  }

  # ===========================================================================
  # Bait normalization acknowledgment
  # ===========================================================================
  de_bait_acknowledged <- reactiveVal(FALSE)

  observeEvent(rv$lf_ip_sub, {
    # Only show the acknowledgment modal when the wizard is NOT open.
    # If the wizard is open, rv$lf_ip_sub is being set by the wizard itself;
    # the inline info box in the wizard already explains bait-group requirements.
    if (identical(rv$lf_ip_sub, "bait") &&
        identical(rv$lf_pipeline, "ip") &&
        !de_bait_acknowledged() &&
        !isTRUE(rv$wizard_modal_open)) {
      showModal(modalDialog(
        title = tags$span(icon("triangle-exclamation"), " Bait Normalization Notice"),
        tags$div(class = "alert alert-warning",
          tags$p(
            tags$strong("Bait subtraction is applied to ALL summarized groups."),
            " The bait protein\u2019s log\u2082 intensity per run is subtracted from every ",
            "protein in every sample that went through summarization."
          ),
          tags$p(
            "This is only valid when the bait protein was ", tags$strong("physically immunoprecipitated"),
            " in ", tags$strong("every group"), " included in the analysis. For example:",
            tags$ul(
              tags$li("WT bait (Day 0) vs Mutant bait (Day 0) \u2014 \u2705 valid"),
              tags$li("Bait IP Day 0 vs Bait IP Day 10 \u2014 \u2705 valid"),
              tags$li("Bait IP vs IgG control \u2014 \u274c ", tags$strong("invalid"))
            )
          ),
          tags$p(class = "mb-0",
            "IgG or empty-bead control groups should be ", tags$strong("excluded"),
            " via the group checkboxes in the upload wizard (Step 5).")
        ),
        size = "m", easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("de_bait_norm_ack", "I understand, continue",
                       class = "btn-warning")
        )
      ))
    }
  })

  observeEvent(input$de_bait_norm_ack, {
    de_bait_acknowledged(TRUE)
    removeModal()
  })

  # Reset acknowledgment when switching away from bait normalization
  observeEvent(rv$lf_ip_sub, {
    if (!identical(rv$lf_ip_sub, "bait")) {
      de_bait_acknowledged(FALSE)
    }
  }, priority = 10)

  # Show the bait ack modal the first time the wizard closes with bait mode active.
  # The ip_sub observer above is suppressed while the wizard is open, so this
  # catches the case where the user selected bait inside the wizard and then closed it.
  observeEvent(rv$wizard_modal_open, {
    if (!isTRUE(rv$wizard_modal_open) &&
        identical(rv$lf_ip_sub, "bait") &&
        identical(rv$lf_pipeline, "ip") &&
        !de_bait_acknowledged()) {
      showModal(modalDialog(
        title = tags$span(icon("triangle-exclamation"), " Bait Normalization Notice"),
        tags$div(class = "alert alert-warning",
          tags$p(
            tags$strong("Bait subtraction is applied to ALL summarized groups."),
            " The bait protein\u2019s log\u2082 intensity per run is subtracted from every ",
            "protein in every sample that went through summarization."
          ),
          tags$p(
            "This is only valid when the bait protein was ", tags$strong("physically immunoprecipitated"),
            " in ", tags$strong("every group"), " included in the analysis. For example:",
            tags$ul(
              tags$li("WT bait (Day 0) vs Mutant bait (Day 0) \u2014 \u2705 valid"),
              tags$li("Bait IP Day 0 vs Bait IP Day 10 \u2014 \u2705 valid"),
              tags$li("Bait IP vs IgG control \u2014 \u274c ", tags$strong("invalid"))
            )
          ),
          tags$p(class = "mb-0",
            "IgG or empty-bead control groups should be ", tags$strong("excluded"),
            " via the group checkboxes in the upload wizard (Step 5).")
        ),
        size = "m", easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("de_bait_norm_ack", "I understand, continue",
                       class = "btn-warning")
        )
      ))
    }
  }, ignoreInit = TRUE)

  # ===========================================================================
  # RUN DE ANALYSIS  -" dispatcher
  # ===========================================================================
  .de_dispatch <- function() {
    pipeline  <- isolate(rv$lf_pipeline)
    ip_sub    <- isolate(rv$lf_ip_sub)
    de_engine <- isolate(rv$lf_de_engine) %||% "proda"

    # Capture ctrl/cond ONCE here, before any async work or potential
    # renderUI re-renders that could reset input$de_ctrl/cond_group.
    .cc <- .get_de_ctrl_cond()

    # ── Quality-flag gate: override to Welch if replicate structure is broken ──
    # Applies globally (all pipelines) because the math failure is pipeline-agnostic.
    rep_flag <- rv$lf_rep_quality_flag %||%
                rv$lf_master_cache_meta$rep_quality_flag
    if (identical(rep_flag$severity, "block") &&
        de_engine %in% c("proda", "saint")) {
      cache <- rv$lf_master_cache
      if (is.null(cache)) {
        showNotification("Run summarization first.", type = "warning"); return()
      }
      logger(paste0("[DE Dispatch] Quality flag BLOCK (", rep_flag$worst_samp,
                    " NA=", round(rep_flag$worst_rate * 100), "%)  -  routing to Peptide-Level Background Z-Test"))
      run_welch_from_cache(cache, ctrl_override = .cc$ctrl, cond_override = .cc$cond)
      return()
    }

    # ── SAINT path ──
    if (identical(de_engine, "saint")) {
      if (!isTRUE(.SAINT_AVAILABLE)) {
        showNotification(
          "SAINTexpress binary not found. Cannot run SAINT analysis.",
          type = "error", duration = 10)
        return()
      }
      cache <- rv$lf_master_cache
      if (is.null(cache)) {
        showNotification("Run summarization first.", type = "warning")
        return()
      }

      run_saint_from_cache(cache)
      return()
    }

    if (identical(de_engine, "proda") && isTRUE(.PRODA_AVAILABLE)) {
      # ── ProDA path  -" all MNAR and non-MNAR ProDA scenarios ──
      cache <- rv$lf_master_cache
      if (!is.null(cache) && !is.null(rv$lf_master_cache_meta) &&
          identical(rv$lf_master_cache_meta$de_engine, "proda")) {
        logger("DE [ProDA]: Using cached ProDA summarization")
        run_proda_from_cache(cache)
      } else {
        grps <- rv$lf_cond_groups
        if (length(grps) < 2) {
          showNotification(
            "Define at least 2 condition groups before running DE.",
            type = "warning")
          return()
        }
        .de_proda_summarize_and_run(igg_conds = rv$lf_igg_conditions %||% NULL)
      }
    }
  }

  # ===========================================================================
  # TMT DE (unified route -- called when rv$data_source == "tmt")
  # ===========================================================================
  .run_de_tmt <- function() {
    if (is.null(rv$tmt_limma_cache)) {
      showNotification("Please run Summarize Data first (Upload Data tab).", type = "warning")
      return()
    }
    cond_test <- isolate(input$de_cond_group)
    cond_ref  <- isolate(input$de_ctrl_group)
    if (is.null(cond_test) || is.null(cond_ref) ||
        nchar(cond_test) == 0 || nchar(cond_ref) == 0) {
      showNotification("Select test and reference conditions in the sidebar.",
                       type = "warning")
      return()
    }
    if (cond_test == cond_ref) {
      showNotification("Test and reference conditions must be different.", type = "warning")
      return()
    }
    if (identical(tmt_de_task$status(), "running")) {
      showNotification("TMT DE already running.", type = "warning")
      return()
    }
    shinyjs::disable("de_run_btn")
    acc_to_gene_snapshot <- rv$tmt_limma_cache$gene_names %||% NULL
    active_cache  <- rv$tmt_limma_cache
    active_design <- as.integer(rv$tmt_limma_cache$design_option %||% rv$tmt_design_option %||% 1L)
    logger(paste0("[DE TMT] Invoking async TMT DE: ", cond_test, " vs ", cond_ref))
    show_async_overlay(session,
                       label  = paste0("TMT DE: ", cond_test, " vs ", cond_ref),
                       detail = "Please wait...")
    tmt_de_cancelled(FALSE)
    tmt_de_task$invoke(
      tmt_master_cache   = active_cache,
      design_opt         = active_design,
      cond_test          = cond_test,
      cond_ref           = cond_ref,
      tmt_analysis_route = rv$tmt_analysis_route,
      acc_to_gene        = acc_to_gene_snapshot
    )
  }

  # ===========================================================================
  # DIA DE (unified route  -" called when rv$data_source == "dia")
  # ===========================================================================
  .run_de_dia <- function() {
    if (is.null(rv$dia_master_cache)) {
      showNotification("Please complete the DIA Upload Wizard first.", type = "warning")
      return()
    }
    cond_test <- isolate(input$de_cond_group)
    cond_ref  <- isolate(input$de_ctrl_group)
    if (is.null(cond_test) || is.null(cond_ref) ||
        nchar(cond_test) == 0 || nchar(cond_ref) == 0) {
      showNotification("Select test and reference conditions in the sidebar.",
                       type = "warning")
      return()
    }
    if (cond_test == cond_ref) {
      showNotification("Test and reference conditions must be different.", type = "warning")
      return()
    }

    shinyjs::disable("de_run_btn")
    on.exit(shinyjs::enable("de_run_btn"), add = TRUE)

    prog_msg <- if (isTRUE(de_rv$dia_use_lod))
      "Running DIA DE (Perseus LOD + permutation FDR)..."
    else
      "Running DIA DE (limma)..."
    result <- withProgress(message = prog_msg, value = 0.2, {
      res <- tryCatch(
        dia_limma_de(
          cache            = rv$dia_master_cache,
          cond_test        = cond_test,
          cond_ref         = cond_ref,
          batch_groups     = rv$dia_batch_groups,
          batch_for_model  = isTRUE(rv$dia_batch_for_model),
          subject_pairing  = rv$dia_subject_pairing,
          factorial_design = rv$dia_factorial_design,
          perseus_lod      = isTRUE(de_rv$dia_use_lod),
          # all conditions in one model unless they differ clearly in variability
          # (dia_fit_scope_guard(), threshold 1.5)
          fit_all_groups   = "auto",
          logger           = logger
        ),
        error = function(e) list(error = conditionMessage(e))
      )
      incProgress(0.8)
      res
    })

    if (!is.null(result$error)) {
      showNotification(paste0("DE failed: ", result$error), type = "error"); return()
    }

    rdf <- result$result_df
    n_res <- nrow(rdf)

    # --- Per-group detection counts from n_valid_* columns ---
    col_test <- paste0("n_valid_", make.names(cond_test))
    col_ref  <- paste0("n_valid_", make.names(cond_ref))
    grps     <- rv$dia_cond_groups
    n_cond   <- length(grps[[cond_test]])
    n_ctrl   <- length(grps[[cond_ref]])
    det_cond <- if (col_test %in% colnames(rdf)) rdf[[col_test]] else rep(n_cond, n_res)
    det_ctrl <- if (col_ref  %in% colnames(rdf)) rdf[[col_ref]]  else rep(n_ctrl, n_res)

    prot_cond_int <- rdf$AveExpr   # log2 average intensity across compared samples

    # Result columns (display/export)
    rdf$AvgIntensity  <- round(prot_cond_int, 2)
    rdf$CondDetected  <- paste0(det_cond, "/", n_cond)
    rdf$CtrlDetected  <- paste0(det_ctrl, "/", n_ctrl)
    # Canonical column aliases expected by render_proda_plot and GSEA handlers
    rdf$Accession     <- rdf$Protein
    rdf$GeneName      <- rdf$Gene
    rdf$Log2FC        <- rdf$log2FC
    rdf$AdjPvalue     <- rdf$adj.pvalue

    # In LOD mode the volcano Y-axis (and the p-value shown in its hover, which
    # is labelled "p (Welch, s0=0.1)") must stay on the RAW p-value. dia_limma_de
    # now returns that as raw_pvalue and puts the permutation FDR in adj.pvalue.
    # Significance in the LOD volcano is driven by PermFDR below, so the plot is
    # unchanged. In standard mode raw_pvalue is absent and this is a no-op.
    .vol_p <- if (!is.null(rdf$raw_pvalue)) rdf$raw_pvalue else rdf$adj.pvalue
    plot_df <- data.frame(
      Accession  = rdf$Protein,
      GeneName   = rdf$Gene,
      Log2FC     = rdf$log2FC,
      AdjPValue  = .vol_p,
      NegLog10P  = -log10(pmax(.vol_p, 1e-300)),
      PermFDR    = if (!is.null(rdf$perm_fdr)) rdf$perm_fdr else NA_real_,
      IsCapped   = FALSE,
      stringsAsFactors = FALSE
    )

    # Attach per-sample raw + imputed abundance columns to rdf for CSV export
    .add_abund_cols <- function(df, mat, prefix) {
      if (is.null(mat) || nrow(mat) == 0) return(df)
      mat_t <- as.data.frame(t(mat[df$Protein, , drop = FALSE]))
      colnames(mat_t) <- paste0(prefix, colnames(mat_t))
      mat_t$Protein <- colnames(mat)  # placeholder; need row-aligned merge
      # row-aligned attach (mat rows == result_df rows via Protein key)
      m <- mat[df$Protein, , drop = FALSE]
      col_df <- as.data.frame(m, check.names = FALSE)
      colnames(col_df) <- paste0(prefix, colnames(col_df))
      cbind(df, col_df)
    }
    rdf <- .add_abund_cols(rdf, result$raw_mat,     "Raw_")
    rdf <- .add_abund_cols(rdf, result$imputed_mat, "Imp_")

    de_rv$result_df      <- rdf
    de_rv$plot_data      <- plot_df
    de_rv$analysis_mode  <- "dia"
    de_rv$cond_test      <- cond_test
    de_rv$cond_ref       <- cond_ref
    de_rv$render_count   <- isolate(de_rv$render_count) + 1L

    # Store LOD boundary parameters for volcano curve
    if (isTRUE(de_rv$dia_use_lod) && !is.null(rdf$perm_fdr) && !is.null(rdf$t_stat)) {
      sig_t <- rdf$t_stat[!is.na(rdf$perm_fdr) & rdf$perm_fdr <= 0.05]
      de_rv$lod_delta_star <- if (length(sig_t) > 0) min(abs(sig_t)) else NA_real_
      # Use median per-protein Welch df (stored in df_welch column) so the
      # boundary curve accounts for heterogeneous within-group variance.
      # Falls back to pooled n-2 if df_welch is unavailable (legacy RDS).
      de_rv$lod_df_rep <- if (!is.null(rdf$df_welch) && any(!is.na(rdf$df_welch)))
        max(median(rdf$df_welch, na.rm = TRUE), 1)
      else
        max(n_cond + n_ctrl - 2, 1)
    } else {
      de_rv$lod_delta_star <- NULL
      de_rv$lod_df_rep     <- NULL
    }

    # For LOD mode: perm_fdr column holds significance; adj.pvalue = raw p (for Y-axis)
    .lod_n_sig <- if (!is.null(rdf$perm_fdr)) {
      sum(rdf$perm_fdr <= 0.05, na.rm = TRUE)
    } else {
      sum(rdf$adj.pvalue < 0.05, na.rm = TRUE)
    }
# With NAs retained, lmFit returns NA for proteins with too few observations
    # for the design. These are inside n_tested but excluded from n_sig by
    # na.rm = TRUE, so the denominator would otherwise overstate what was tested.
    n_untestable <- sum(is.na(rdf$adj.pvalue))
    if (n_untestable > 0)
      logger(paste0("[DE DIA] ", n_untestable,
                    " proteins had insufficient observations to fit and returned NA."))
    rv$provenance$de <- list(
      source        = "dia",
      n_untestable  = n_untestable,
      analysis_mode = "dia",
      engine        = "limma",
      software_used = "limma",
      perseus_lod   = isTRUE(de_rv$dia_use_lod),
      comparison    = result$comparison,
      n_tested      = result$n_tested,
      min_valid_frac           = result$min_valid_frac           %||% 0.5,
      require_both             = result$require_both             %||% TRUE,
      min_n_test               = result$min_n_test               %||% NA_integer_,
      min_n_ref                = result$min_n_ref                %||% NA_integer_,
      n_test_samples           = result$n_test_samples           %||% NA_integer_,
      n_ref_samples            = result$n_ref_samples            %||% NA_integer_,
      max_missing_per_cond     = result$max_missing_per_cond     %||% NA_real_,
      knn_k                    = result$knn_k                    %||% NA_integer_,
      n_imputed                = result$n_imputed                %||% 0L,
      n_knn_rows               = result$n_knn_rows               %||% 0L,
      n_left_na                = result$n_left_na                %||% 0L,
      fit_all_groups           = isTRUE(result$fit_all_groups),
      fit_scope_reason         = result$fit_scope_guard$reason    %||% NA_character_,
      fit_scope_ratio          = result$fit_scope_guard$ratio     %||% NA_real_,
      fit_scope_threshold      = result$fit_scope_guard$threshold %||% NA_real_,
      n_conditions_fitted      = result$n_conditions_fitted      %||% NA_integer_,
      n_samples_fitted         = result$n_samples_fitted         %||% NA_integer_,
      n_sig         = .lod_n_sig,
      n_nmar        = if (!is.null(de_rv$dia_nmar_result)) nrow(de_rv$dia_nmar_result) else 0L,
      batch_corrected = isTRUE(result$batch_in_model),
      subject_paired  = !is.null(rv$dia_subject_pairing),
      cond_group    = cond_test,
      ctrl_group    = cond_ref,
      timestamp     = as.character(Sys.time())
    )
    logger(paste0("[DE DIA] Done: ", result$n_tested, " proteins, ",
                  .lod_n_sig, " significant.",
                  if (!is.null(rv$dia_subject_pairing))
                    paste0(" Paired design (duplicateCorrelation, ",
                           length(unique(rv$dia_subject_pairing)), " blocks).")
                  else " Unpaired design."))

    # NMAR analysis  -" only in standard (non-LOD) mode
    # Uses the 50% filter independently to determine which proteins were excluded
    de_rv$dia_nmar_test <- cond_test
    de_rv$dia_nmar_ref  <- cond_ref
    if (isTRUE(de_rv$dia_use_lod)) {
      # LOD mode imputes everything → NMAR concept does not apply
      de_rv$dia_nmar_result   <- NULL
      rv$provenance$de$n_nmar <- 0L
    } else {
      # tested_proteins = those that passed the 50% filter (i.e. rdf$Protein)
      nmar_df <- tryCatch(
        dia_nmar_fisher(
          cache           = rv$dia_master_cache,
          cond_test       = cond_test,
          cond_ref        = cond_ref,
          tested_proteins = rdf$Protein,
          cond_groups     = rv$dia_cond_groups,
          logger          = logger
        ),
        error = function(e) { logger(paste0("[NMAR] Error: ", e$message)); NULL }
      )
      de_rv$dia_nmar_result <- nmar_df
      # Provenance was written above, before NMAR ran, so n_nmar was stale
      # (0 on first run, previous comparison's count thereafter). Correct it now.
      rv$provenance$de$n_nmar <- if (!is.null(nmar_df)) nrow(nmar_df) else 0L
      # dia_nmar_fisher() exposes these run-level values as constant data-frame
      # columns rather than as scalars, so take the first element. Left absent
      # (not NA) when no NMAR analysis ran, so the methods writer falls back to
      # its unconditional wording rather than reporting a value it does not have.
      if (!is.null(nmar_df) && nrow(nmar_df) > 0L) {
        if (!is.null(nmar_df$fisher_has_power))
          rv$provenance$de$fisher_has_power <- nmar_df$fisher_has_power[1]
        if (!is.null(nmar_df$min_achievable_p))
          rv$provenance$de$fisher_min_achievable_p <- nmar_df$min_achievable_p[1]
      }
      if (!is.null(nmar_df))
        logger(paste0("[NMAR] ", nrow(nmar_df), " excluded proteins analysed."))
    }
  }

  # ===========================================================================
  observeEvent(input$de_run_btn, {
    # Route to TMT or DIA path if appropriate data source is active
    src <- rv$data_source
    if (identical(src, "tmt")) {
      .run_de_tmt()
      return()
    } else if (identical(src, "dia")) {
      .run_de_dia()
      return()
    }

    # LF DDA path (unchanged)
    if (identical(rv$lf_input_mode, "maxquant")) {
      if (is.null(rv$lf_master_cache)) {
        showNotification("Please run summarization before DE analysis.", type = "warning")
        return()
      }
    } else if (!is.null(rv$lf_master_cache)) {
      # Cache already present (from summarize button or Saved Experiments import).
      # Raw lines are not needed - proceed directly.
    } else {
      if (is.null(rv$pd_raw_lines) || is.null(rv$pd_protein_line_idx)) {
        showNotification("Please upload and parse a data file before running DE analysis.", type = "warning")
        return()
      }
    }

    shinyjs::disable("de_run_btn")
    on.exit(shinyjs::enable("de_run_btn"), add = TRUE)

    pipeline  <- isolate(rv$lf_pipeline)
    ip_sub    <- isolate(rv$lf_ip_sub)
    de_engine <- isolate(rv$lf_de_engine) %||% "proda"
    de_rv$analysis_mode <- de_engine

    # Check if we need re-summarization and confirm with user
    cache <- rv$lf_master_cache
    meta  <- rv$lf_master_cache_meta
    needs_resummarize <- is.null(cache) || is.null(meta) ||
      !identical(meta$de_engine, de_engine) ||
      FALSE

    if (needs_resummarize && isTRUE(rv$lf_cache_stale)) {
      showModal(modalDialog(
        title = tags$span(icon("rotate"), " Re-Summarization Required"),
        tags$div(class = "alert alert-info mb-0",
          tags$p(
            tags$b("A new analysis pipeline has been selected."),
            tags$br(),
            "The data will be re-summarized with the current settings before running ",
            "the DE analysis. This may take a few minutes."
          ),
          tags$p(class = "text-muted small mb-0",
            "Engine: ", tags$b(de_engine),
            " / Pipeline: ", tags$b(pipeline),
            if (identical(pipeline, "ip")) paste0(" (", ip_sub, ")") else "")
        ),
        size = "m", easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("de_confirm_resummarize", "Continue",
                       class = "btn-primary", icon = icon("play"))
        )
      ))
      return()
    }

    # Reset previous results
    de_rv$result_df      <- NULL
    de_rv$plot_data      <- NULL
    de_rv$highlight_gene <- NULL
    de_rv$raw_abund_df   <- NULL

    .de_dispatch()
  })

  observeEvent(input$de_confirm_resummarize, {
    removeModal()
    rv$lf_cache_stale <- FALSE
    .de_dispatch()
  })


  # ===========================================================================
  # MODE CACHE (ProDA): Run ProDA DE from pre-built ProDA cache
  # Async invoker: validation runs synchronously; heavy compute runs in
  # proda_de_task (ExtendedTask). Results handled by the result observer.
  # ===========================================================================
  run_proda_from_cache <- function(cache) {
    .cc        <- .get_de_ctrl_cond()
    ctrl_label <- .cc$ctrl
    cond_label <- .cc$cond

    if (is.null(ctrl_label) || !nzchar(ctrl_label) ||
        is.null(cond_label) || !nzchar(cond_label)) {
      showNotification("Select both Control and Condition groups.", type = "warning")
      return()
    }

    all_conditions <- unique(cache$design_df$condition)
    check_group <- function(label, role) {
      if (label %in% all_conditions) return(TRUE)
      showNotification(
        paste0(role, " group \u2018", label, "\u2019 not found in ProDA cache. ",
               "Available: ", paste(all_conditions, collapse = ", ")),
        type = "error", duration = 10)
      FALSE
    }
    if (!check_group(ctrl_label, "Control"))   return()
    if (!check_group(cond_label, "Condition")) return()

    if (identical(proda_de_task$status(), "running")) {
      showNotification("ProDA DE is already running.", type = "warning")
      return()
    }

    proda_cancelled(FALSE)
    show_async_overlay(session,
                       "Running ProDA DE\u2026",
                       "Fitting probabilistic dropout model")

    proda_de_task$invoke(
      cache                = cache,
      ctrl_label           = ctrl_label,
      cond_label           = cond_label,
      cond_groups          = isolate(rv$lf_cond_groups),
      factorial_design     = isolate(rv$lf_factorial_design),
      batch_for_model      = isTRUE(isolate(rv$lf_batch_for_model)),
      batch_groups         = isolate(rv$lf_batch_groups),
      batch_filter_enabled = isTRUE(isolate(rv$lf_batch_filter_enabled)),
      subject_pairing      = isolate(rv$lf_subject_pairing),
      expect_mnar          = isTRUE(isolate(rv$lf_expect_mnar)),
      ip_sub               = isolate(rv$lf_ip_sub),
      igg_conditions       = isolate(rv$lf_igg_conditions) %||% character(0)
    )
  }

  # ===========================================================================
  # SAINT PATHWAY  -" run_saint_from_cache
  # Batch: SAINTexpress has no covariate support; batch is viz-only (PCA/Heatmap).
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # run_saint_from_cache()
  # Async invoker: validation + column resolution runs synchronously;
  # heavy compute runs in saint_de_task (ExtendedTask).
  # ---------------------------------------------------------------------------
  run_saint_from_cache <- function(cache) {
    .cc        <- .get_de_ctrl_cond()
    ctrl_label <- .cc$ctrl
    cond_label <- .cc$cond

    if (!nzchar(ctrl_label) || !nzchar(cond_label)) {
      showNotification("Select Control and Condition groups first.", type = "warning")
      return()
    }
    if (identical(ctrl_label, cond_label)) {
      showNotification("Control and Condition groups must be different.", type = "warning")
      return()
    }

    grps      <- rv$lf_cond_groups
    ctrl_cols <- grps[[ctrl_label]]
    cond_cols <- grps[[cond_label]]

    # Remap column names from lf_cond_groups through design_df (wizard rename)
    if (!is.null(cache$design_df) && "sample" %in% colnames(cache$design_df)) {
      col_map   <- setNames(cache$design_df$sample_name, cache$design_df$sample)
      ctrl_cols <- unname(ifelse(ctrl_cols %in% names(col_map),
                                 col_map[ctrl_cols], ctrl_cols))
      cond_cols <- unname(ifelse(cond_cols %in% names(col_map),
                                 col_map[cond_cols], cond_cols))
    }

    if (is.null(ctrl_cols) || length(ctrl_cols) < 2) {
      if (isTRUE(.PRODA_AVAILABLE)) {
        showNotification(
          tagList(
            tags$strong("SAINTexpress requires \u22652 control replicates \u2014 falling back to ProDA."),
            tags$br(),
            "With only 1 control replicate the Bayesian background model cannot be estimated. ",
            "Analysis will proceed with ", tags$b("ProDA (SAINTexpress input matrix, no normalization)"), "."
          ),
          type = "warning", duration = 12L
        )
        logger("[SAINT gate] < 2 ctrl replicates \u2014 falling back to ProDA")
        run_proda_from_cache(cache)
      } else {
        showNotification(
          "SAINTexpress requires \u22652 control replicates and ProDA is not available.",
          type = "error", duration = 15L)
      }
      return()
    }
    if (is.null(cond_cols) || length(cond_cols) < 2) {
      if (isTRUE(.PRODA_AVAILABLE)) {
        showNotification(
          tagList(
            tags$strong("SAINTexpress requires \u22652 bait replicates \u2014 falling back to ProDA."),
            tags$br(),
            "Analysis will proceed with ", tags$b("ProDA (SAINTexpress input matrix, no normalization)"), "."
          ),
          type = "warning", duration = 12L
        )
        logger("[SAINT gate] < 2 bait replicates \u2014 falling back to ProDA")
        run_proda_from_cache(cache)
      } else {
        showNotification(
          "SAINTexpress requires \u22652 bait replicates and ProDA is not available.",
          type = "error", duration = 10L)
      }
      return()
    }

    bait_gene <- isolate(rv$lf_bait_protein) %||% ""
    if (!nzchar(bait_gene)) {
      showNotification("No bait gene recorded. Re-run the wizard and enter a bait gene name.",
                       type = "error")
      return()
    }

    # Subset + zero-variance filter (synchronous, fast)
    protein_matrix <- cache$protein_matrix
    all_cols  <- c(ctrl_cols, cond_cols)
    avail     <- intersect(all_cols, colnames(protein_matrix))
    if (length(avail) < length(all_cols)) {
      missing <- setdiff(all_cols, colnames(protein_matrix))
      showNotification(paste0("Columns not found in matrix: ",
                               paste(missing, collapse = ", ")),
                       type = "error")
      return()
    }
    pm_sub   <- protein_matrix[, avail, drop = FALSE]
    zero_var <- apply(pm_sub, 1, function(r) {
      vals <- r[!is.na(r)]
      length(vals) < 2
    })
    if (sum(!zero_var) < 5) {
      showNotification("Too few variable proteins for SAINT analysis.", type = "error")
      return()
    }
    # Build a trimmed cache with the filtered matrix for the worker
    cache_sub <- cache
    cache_sub$protein_matrix    <- pm_sub[!zero_var, , drop = FALSE]
    cache_sub$pre_saint_input   <- nrow(pm_sub)
    cache_sub$pre_saint_dropped <- sum(zero_var, na.rm = TRUE)

    if (identical(saint_de_task$status(), "running")) {
      showNotification("SAINTexpress is already running.", type = "warning")
      return()
    }

    logger(paste0("[SAINT] Running for bait=", bait_gene,
                  " ctrl=", ctrl_label, "(", length(ctrl_cols), " reps)",
                  " cond=", cond_label, "(", length(cond_cols), " reps)",
                  " proteins=", sum(!zero_var)))

    # Resolve binary path and logs dir in main thread (absolute paths for worker)
    saint_bin  <- tryCatch(normalizePath(.saint_int_binary(), mustWork = FALSE),
                           error = function(e) .saint_int_binary())
    logs_dir   <- tryCatch(normalizePath(file.path(getwd(), "logs"), mustWork = FALSE),
                           error = function(e) NULL)
    if (!dir.exists(logs_dir)) logs_dir <- NULL

    saint_cancelled(FALSE)
    show_async_overlay(session,
                       "Running SAINTexpress\u2026",
                       paste0("Bait: ", bait_gene))

    saint_de_task$invoke(
      cache             = cache_sub,
      ctrl_label        = ctrl_label,
      cond_label        = cond_label,
      ctrl_cols         = ctrl_cols,
      cond_cols         = cond_cols,
      bait_gene         = bait_gene,
      cond_groups       = isolate(rv$lf_cond_groups),
      saint_binary_path = saint_bin,
      logs_dir          = logs_dir
    )
  }

  # ===========================================================================
  # PEPTIDE-LEVEL BACKGROUND Z-TEST  -" run_welch_from_cache()
  # Replicates Proteome Discoverer's "t-test (Background Based)" architecture:
  #    -¢ Flatten valid peptide-level log2 FCs (both groups must have ≥1 obs  -" no imputation)
  #    -¢ Estimate global background SD via iterative 3σ clipping on real paired peptides only
  #    -¢ Roll up to protein: z = mean_peptide_log2FC * sqrt(N_peptides) / background_SD
  #    -¢ FDR: Benjamini -"Hochberg on all proteins with N_peptides ≥ 1
  # Requires nested PD CSV with peptide data (rv$pd_raw_lines).
  # If only protein-level data is available, returns an advisory error.
  # Called automatically when lf_rep_quality_flag$severity == "block".
  # Named in UI as: "Peptide-Level Background Z-Test (Iterative 3σ Clipping)"
  # ===========================================================================

  .calc_sigma_clip_sd <- function(x, sigma = 3, max_iter = 10) {
    x  <- x[!is.na(x)]
    mu <- median(x)
    s  <- sd(x)
    for (i in seq_len(max_iter)) {
      keep <- abs(x - mu) < sigma * s
      ns   <- sd(x[keep])
      nm   <- median(x[keep])
      if (abs(s - ns) < 1e-4) { s <- ns; break }
      s  <- ns
      mu <- nm
    }
    s
  }

  run_welch_from_cache <- function(cache, ctrl_override = NULL, cond_override = NULL) {
    # Use explicitly-passed labels (captured at dispatch time) when available;
    # fall back to reading the pickers directly only for standalone calls.
    if (!is.null(ctrl_override) && nzchar(ctrl_override)) {
      ctrl_label <- ctrl_override
      cond_label <- cond_override %||% ""
    } else {
      .cc        <- .get_de_ctrl_cond()
      ctrl_label <- .cc$ctrl
      cond_label <- .cc$cond
    }

    if (is.null(ctrl_label) || !nzchar(ctrl_label) ||
        is.null(cond_label) || !nzchar(cond_label)) {
      showNotification("Select both Control and Condition groups.", type = "warning")
      return()
    }

    design_df      <- cache$design_df
    protein_matrix <- cache$protein_matrix
    acc_to_gene    <- cache$acc_to_gene

    all_conditions <- unique(design_df$condition)
    for (lbl in c(ctrl_label, cond_label)) {
      if (!lbl %in% all_conditions) {
        showNotification(
          paste0("Group '", lbl, "' not found in cache. Available: ",
                 paste(all_conditions, collapse = ", ")),
          type = "error", duration = 10)
        return()
      }
    }

    # ── Check peptide data availability ──────────────────────────────────────
    # Peptide-level z-test requires nested PD CSV; protein-only exports cannot
    # provide per-peptide variance estimates.
    all_lines   <- rv$pd_raw_lines %||% character(0)
    pep_hdr_idx <- grep("^,.*Annotated Sequence", all_lines)
    has_peptides <- length(pep_hdr_idx) > 0

    if (!has_peptides) {
      showModal(modalDialog(
        title = tags$span(icon("triangle-exclamation"), " Peptide Data Required"),
        tagList(
          tags$div(class = "alert alert-danger",
            tags$strong("The Peptide-Level Background Z-Test cannot run."),
            tags$br(),
            "Your uploaded file contains only protein-level abundance data. ",
            "The background z-test requires nested peptide-level data ",
            "(Annotated Sequence rows) to estimate global variance and compute ",
            "statistically powered protein-level p-values."
          ),
          tags$hr(),
          tags$p(tags$strong("What you can do:")),
          tags$ul(
            tags$li(
              tags$strong("Recommended: "),
              "Re-export from Proteome Discoverer with peptide groups included ",
              "(File \u2192 Export \u2192 Result Files \u2192 check both Proteins and Peptide Groups)."
            ),
            tags$li(
              tags$strong("Alternative: "),
              "Remove the low-quality replicate (", flag$worst_samp %||% "the flagged sample",
              ") and re-run DE analysis with the remaining replicates."
            )
          ),
          tags$p(class = "text-muted small",
            "Without peptide data, a per-protein Welch t-test with n=3 replicates ",
            "cannot achieve statistical significance when one replicate has systematic ",
            "intensity suppression \u2014 the variance inflation dominates the denominator."
          )
        ),
        footer = modalButton("OK"), easyClose = TRUE, size = "m"
      ))
      return()
    }

    withProgress(message = "DE Analysis \u2014 Peptide-Level Background Z-Test...", value = 0, {
      incProgress(0.10, detail = "Selecting group columns...")

      ctrl_samps <- design_df$sample_name[design_df$condition == ctrl_label]
      cond_samps <- design_df$sample_name[design_df$condition == cond_label]
      ctrl_samps <- intersect(ctrl_samps, colnames(protein_matrix))
      cond_samps <- intersect(cond_samps, colnames(protein_matrix))

      ctrl_mat <- protein_matrix[, ctrl_samps, drop = FALSE]  # log2 normalized values
      cond_mat <- protein_matrix[, cond_samps, drop = FALSE]
      proteins   <- rownames(protein_matrix)
      n_proteins <- length(proteins)

      incProgress(0.20, detail = "Parsing peptide-level data...")
      logger(paste0("[BgZ] Parsing peptides from nested CSV (",
                    length(all_lines), " lines)"))

      # ── Parse peptide table from nested PD CSV ──────────────────────────────
      all_nested_idx <- grep("^,", all_lines)
      pep_data_idx   <- setdiff(all_nested_idx, pep_hdr_idx)
      pep_block      <- paste(c(all_lines[pep_hdr_idx[1]], all_lines[pep_data_idx]),
                              collapse = "\n")
      pep_df <- tryCatch(
        read.csv(text = pep_block, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, na.strings = c("", "NA")),
        error = function(e) { logger(paste0("[BgZ] Peptide parse error: ", e$message)); NULL }
      )
      if (is.null(pep_df) || nrow(pep_df) < 5) {
        showNotification("Peptide data parse error.", type = "error"); return()
      }
      logger(paste0("[BgZ] Parsed ", nrow(pep_df), " peptides"))

      # ── Map peptide abundance columns to ctrl/cond ──────────────────────────
      pep_cn       <- colnames(pep_df)
      master_col   <- grep("Master Protein Accessions", pep_cn, value = TRUE)[1]
      # Match by sample name (ctrl_samps/cond_samps are column names in protein_matrix
      # which use short names; peptide CSV uses full "Abundance: F1: ..." style)
      # Strategy: match any column whose name contains the sample name as a substring
      find_pep_cols <- function(samp_names) {
        cols <- character(0)
        for (sn in samp_names) {
          # exact match first, then substring
          m <- pep_cn[pep_cn == sn]
          if (!length(m)) m <- pep_cn[grep(sn, pep_cn, fixed = TRUE)]
          if (!length(m)) m <- pep_cn[grep(paste0("\\b", sn, "\\b"), pep_cn)]
          if (length(m)) cols <- c(cols, m[1])
        }
        cols
      }
      # Also try matching by "F<n>" run tag embedded in column names
      design_ctrl <- design_df[design_df$condition == ctrl_label, ]
      design_cond <- design_df[design_df$condition == cond_label, ]

      # Extract run IDs from design (e.g. "F1", "F3", "F5")
      run_tag_ctrl <- design_ctrl$run_id %||% design_ctrl$sample_name
      run_tag_cond <- design_cond$run_id %||% design_cond$sample_name

      # Build ctrl/cond abundance columns in peptide data
      build_pep_abund_cols <- function(run_tags) {
        abund_pep_all <- pep_cn[grep("^Abundance:", pep_cn, ignore.case = TRUE)]
        abund_pep_all <- abund_pep_all[!grepl("ratio|p-value|variability|normalized|grouped|count",
                                               abund_pep_all, ignore.case = TRUE)]
        matched <- character(0)
        for (tag in run_tags) {
          m <- abund_pep_all[grep(paste0("\\b", tag, "\\b"), abund_pep_all)]
          if (!length(m)) m <- abund_pep_all[grep(tag, abund_pep_all, fixed = TRUE)]
          if (length(m)) matched <- c(matched, m[1])
        }
        matched
      }

      ctrl_cols_pep <- build_pep_abund_cols(run_tag_ctrl)
      cond_cols_pep <- build_pep_abund_cols(run_tag_cond)

      # Fallback: use protein matrix column names directly if they exist in peptide CSV
      if (!length(ctrl_cols_pep)) ctrl_cols_pep <- find_pep_cols(ctrl_samps)
      if (!length(cond_cols_pep)) cond_cols_pep <- find_pep_cols(cond_samps)

      if (!length(ctrl_cols_pep) || !length(cond_cols_pep)) {
        # Last-resort: auto-assign first n_ctrl / n_cond abundance cols in peptide CSV
        abund_fallback <- pep_cn[grep("^Abundance:", pep_cn, ignore.case = TRUE)]
        abund_fallback <- abund_fallback[!grepl("ratio|p-value|variability|normalized|grouped|count",
                                                abund_fallback, ignore.case = TRUE)]
        if (length(abund_fallback) >= length(ctrl_samps) + length(cond_samps)) {
          ctrl_cols_pep <- abund_fallback[seq_len(length(ctrl_samps))]
          cond_cols_pep <- abund_fallback[seq_len(length(cond_samps)) + length(ctrl_samps)]
          logger("[BgZ] Peptide column fallback: using positional assignment")
        } else {
          showNotification(
            "Cannot map condition groups to peptide abundance columns. Upload a full PD nested export.",
            type = "error", duration = 12)
          return()
        }
      }

      logger(paste0("[BgZ] Peptide ctrl cols: ", paste(ctrl_cols_pep, collapse = ", ")))
      logger(paste0("[BgZ] Peptide cond cols: ", paste(cond_cols_pep, collapse = ", ")))

      incProgress(0.35, detail = "Computing per-peptide log2 fold-changes...")

      # ── Per-peptide log2 FC (valid pairs only  -" no imputation) ──────────────
      pep_ctrl_mat <- log2(matrix(as.numeric(unlist(pep_df[, ctrl_cols_pep, drop = FALSE])),
                                  ncol = length(ctrl_cols_pep)))
      pep_cond_mat <- log2(matrix(as.numeric(unlist(pep_df[, cond_cols_pep, drop = FALSE])),
                                  ncol = length(cond_cols_pep)))

      pep_ctrl_mean <- rowMeans(pep_ctrl_mat, na.rm = TRUE)
      pep_cond_mean <- rowMeans(pep_cond_mat, na.rm = TRUE)
      pep_ctrl_mean[rowSums(!is.na(pep_ctrl_mat)) == 0] <- NA
      pep_cond_mean[rowSums(!is.na(pep_cond_mat)) == 0] <- NA

      # Valid peptide: must have at least 1 obs in BOTH groups
      pep_log2fc <- pep_cond_mean - pep_ctrl_mean
      valid_pep  <- !is.na(pep_log2fc)
      logger(paste0("[BgZ] Valid paired peptides: ", sum(valid_pep), " / ", nrow(pep_df)))

      if (sum(valid_pep) < 20) {
        showNotification(
          paste0("Only ", sum(valid_pep), " peptides have valid observations in both groups. ",
                 "Background z-test requires \u226520 paired peptides."),
          type = "error", duration = 10)
        return()
      }

      # ── Background SD via iterative 3σ clipping ─────────────────────────────
      incProgress(0.50, detail = "Estimating background variance (3\u03c3 clipping)...")
      bg_sd <- .calc_sigma_clip_sd(pep_log2fc[valid_pep], sigma = 3)
      logger(paste0("[BgZ] Background SD (3\u03c3 clipped): ", round(bg_sd, 4)))

      # ── Roll up to protein level ─────────────────────────────────────────────
      incProgress(0.60, detail = "Rolling up to protein level...")
      pep_protein <- sub(";.*", "", trimws(pep_df[[master_col]]))

      n_pep_valid  <- integer(n_proteins)
      log2fc_vec   <- rep(NA_real_, n_proteins)
      for (i in seq_len(n_proteins)) {
        rows <- pep_protein == proteins[i] & valid_pep
        n_pep_valid[i] <- sum(rows)
        if (n_pep_valid[i] > 0)
          log2fc_vec[i] <- mean(pep_log2fc[rows], na.rm = TRUE)
      }

      # ── Z-test and BH correction ─────────────────────────────────────────────
      incProgress(0.70, detail = "Applying Benjamini-Hochberg correction...")
      se_vec   <- bg_sd / sqrt(pmax(n_pep_valid, 1L))
      z_vec    <- log2fc_vec / se_vec
      pval_vec <- 2 * pnorm(-abs(z_vec))
      pval_vec[n_pep_valid == 0] <- NA

      adj_pval_vec <- rep(NA_real_, n_proteins)
      valid_p      <- which(!is.na(pval_vec))
      if (length(valid_p) > 0)
        adj_pval_vec[valid_p] <- p.adjust(pval_vec[valid_p], method = "BH")

      # For protein-level FC display: also compute directly from protein matrix
      # (normalized, use as a consistent FC estimate for the volcano)
      ctrl_mean_vec <- rowMeans(ctrl_mat, na.rm = TRUE)
      cond_mean_vec <- rowMeans(cond_mat, na.rm = TRUE)
      ctrl_n_vec    <- rowSums(!is.na(ctrl_mat))
      cond_n_vec    <- rowSums(!is.na(cond_mat))
      # Prefer peptide-level FC (more observations), fall back to protein-level
      log2fc_vec    <- ifelse(!is.na(log2fc_vec), log2fc_vec,
                              cond_mean_vec - ctrl_mean_vec)

      incProgress(0.80, detail = "Assembling results...")

      avg_intensity <- rowMeans(cbind(ctrl_mat, cond_mat), na.rm = TRUE)
      gene_vec      <- acc_to_gene[proteins]
      gene_vec[is.na(gene_vec)] <- proteins[is.na(gene_vec)]

      ctrl_rep_df <- as.data.frame(ctrl_mat)
      colnames(ctrl_rep_df) <- paste0("Ctrl_Rep", seq_len(ncol(ctrl_rep_df)))
      ctrl_rep_df$Accession <- proteins
      cond_rep_df <- as.data.frame(cond_mat)
      colnames(cond_rep_df) <- paste0("Cond_Rep", seq_len(ncol(cond_rep_df)))
      cond_rep_df$Accession <- proteins

      result_df <- data.frame(
        GeneName     = gene_vec,
        Accession    = proteins,
        AvgIntensity = round(avg_intensity, 3),
        Log2FC       = round(log2fc_vec, 4),
        NPeptides    = n_pep_valid,
        CtrlDetected = paste0(ctrl_n_vec, "/", length(ctrl_samps)),
        CondDetected = paste0(cond_n_vec, "/", length(cond_samps)),
        AdjPvalue    = round(adj_pval_vec, 6),
        RawPvalue    = round(pval_vec, 6),
        stringsAsFactors = FALSE
      )

      result_df <- merge(result_df, ctrl_rep_df, by = "Accession", all.x = TRUE)
      result_df <- merge(result_df, cond_rep_df, by = "Accession", all.x = TRUE)

      col_order <- c("GeneName", "Accession", "AvgIntensity", "Log2FC", "NPeptides",
                     sort(grep("^Ctrl_Rep|^Cond_Rep", colnames(result_df), value = TRUE)),
                     "CtrlDetected", "CondDetected", "AdjPvalue", "RawPvalue")
      col_order  <- intersect(col_order, colnames(result_df))
      result_df  <- result_df[, col_order, drop = FALSE]
      # Most significant first (untested proteins last), larger effect breaking ties
      result_df  <- result_df[order(result_df$AdjPvalue, -abs(result_df$Log2FC),
                                    na.last = TRUE), ]
      rownames(result_df) <- NULL

      plot_df <- data.frame(
        Accession    = result_df$Accession,
        GeneName     = result_df$GeneName,
        FC_linear    = 2 ^ result_df$Log2FC,
        Log2FC       = result_df$Log2FC,
        AdjPValue    = result_df$AdjPvalue,
        NegLog10P    = -log10(pmax(result_df$AdjPvalue, .Machine$double.xmin)),
        IsCapped     = FALSE,
        IsImputed    = FALSE,
        AvgIntensity = result_df$AvgIntensity,
        stringsAsFactors = FALSE
      )
      plot_df$NegLog10P[is.na(result_df$AdjPvalue)] <- 0

      de_rv$result_df     <- result_df
      de_rv$plot_data     <- plot_df
      de_rv$analysis_mode <- "welch"
      de_rv$acc_to_gene   <- acc_to_gene

      raw_wide           <- as.data.frame(protein_matrix)
      raw_wide$Accession <- rownames(raw_wide)
      de_rv$raw_abund_df <- raw_wide
      de_rv$render_count <- isolate(de_rv$render_count) + 1

      incProgress(1.0, detail = "Complete!")

      n_sig <- sum(!is.na(result_df$AdjPvalue) & result_df$AdjPvalue < 0.05, na.rm = TRUE)
      n_na  <- sum(is.na(result_df$AdjPvalue))
      logger(paste0("[BgZ] ", nrow(result_df), " proteins. adjP<0.05: ", n_sig,
                    "  N_pep=0 (no p-value): ", n_na,
                    "  background_SD=", round(bg_sd, 4)))

      rv$provenance$de <- list(
        pipeline             = rv$lf_pipeline %||% "lf",
        ip_sub               = rv$lf_ip_sub,
        de_engine            = "welch",
        analysis_mode        = "welch",
        pipeline_label       = cache$pipeline_label %||% "Label-Free",
        software_used        = "Peptide-Level Background Z-Test (Iterative 3\u03c3 Clipping)",
        pvalue_adjustment    = "BH (Benjamini-Hochberg)",
        normalization_method = cache$normalization_method %||% "none",
        ctrl_group           = ctrl_label,
        cond_group           = cond_label,
        n_prots_input        = cache$n_prots_input     %||% NA_integer_,
        n_prots_processed    = cache$n_prots_processed %||% NA_integer_,
        n_prots_output       = nrow(result_df),
        n_prots_pval_na      = n_na,
        n_peptides_used      = sum(valid_pep),
        background_sd        = round(bg_sd, 4),
        welch_contingency    = TRUE,
        rep_quality_flag     = rv$lf_rep_quality_flag,
        ebn_ref_proteins     = cache$ebn_ref_vector %||% character(0),
        aggregation_method   = cache$aggregation_method %||% "unknown",
        # backward-compat shims
        normalization        = cache$normalization_method %||% "none",
        experiment_type      = "lf_welch",
        skip_impute          = FALSE,
        lod_n_imputed        = 0L,
        lod_imputation       = FALSE,
        mb_impute            = FALSE,
        batch_corrected      = FALSE,
        batch_groups         = NULL,
        factorial_enabled    = FALSE,
        subject_pairing_used = FALSE,
        starting_proteins    = as.numeric(cache$n_prots_input %||% nrow(result_df)),
        dropped_proteins     = 0
      )
    })

    n_sig <- sum(!is.na(de_rv$result_df$AdjPvalue) & de_rv$result_df$AdjPvalue < 0.05,
                 na.rm = TRUE)
    batch_note <- if (!is.null(rv$lf_batch_groups) && length(rv$lf_batch_groups) > 0)
      tags$span(class = "text-muted",
        " Note: Batch assignment is not applied in Background Z-Test mode (visualization only).")
    else NULL
    showNotification(
      tagList(
        tags$b("Peptide-Level Background Z-Test complete (Contingency Mode): "),
        paste0(nrow(de_rv$result_df), " proteins. "),
        paste0(n_sig, " significant (adj. p < 0.05)."),
        batch_note
      ),
      type = "warning", duration = 12)
  }

  # ===========================================================================
  # AUTO SUMMARIZE + RUN (ProDA path)
  # Builds protein matrix via QFeatures / PD nested, runs lf_core_summarization_proda,
  # caches result, then calls run_proda_from_cache.
  # ===========================================================================
  .de_proda_summarize_and_run <- function(igg_conds = NULL) {
    pipeline  <- isolate(rv$lf_pipeline)
    ip_sub    <- isolate(rv$lf_ip_sub)
    grps      <- rv$lf_cond_groups

    logger(paste0("DE [ProDA Auto-Summarize]: pipeline=", pipeline,
                  " ip_sub=", ip_sub %||% "NULL",
                  " igg=", if (length(igg_conds)) paste(igg_conds, collapse = ",") else "none"))

    showNotification(
      "Building protein matrix for ProDA analysis. This may take a moment...",
      type = "message", duration = NULL, id = "de_proda_summ_wait"
    )

    withProgress(message = "ProDA \u2014 Building Protein Matrix...", value = 0, {
      incProgress(0.1, detail = "Extracting protein-level data...")

      mat_result <- build_proda_protein_matrix(
        rv          = rv,
        cond_groups = grps,
        # See the wizard call site: %||% guards NULL only, and this value can
        # arrive as NA_character_ from a bundle exported on a non-peptide route.
        shared_peptides = { .sp <- rv$lf_shared_peptides
                            if (length(.sp) != 1L || is.na(.sp) ||
                                !nzchar(.sp)) "unique_only" else .sp },
        logger      = logger
      )
      if (!is.null(mat_result$error)) {
        removeNotification(id = "de_proda_summ_wait")
        showNotification(mat_result$error, type = "error", duration = 12)
        return()
      }

      incProgress(0.4, detail = "Normalising and building ProDA cache...")

      cache_result <- lf_core_summarization_proda(
        protein_matrix = mat_result$protein_matrix,
        design_df      = mat_result$design_df,
        acc_to_gene    = mat_result$acc_to_gene,
        ip_sub         = ip_sub,
        igg_conditions = igg_conds %||% character(0),
        ebn_ref_vector = NULL,
        # FALSE only for input normalized upstream (FragPipe/IonQuant, per fragpipe.workflow)
        apply_normalization = lf_apply_normalization_effective(rv),
        logger         = logger
      )
      incProgress(0.9, detail = "Done!")
    })

    removeNotification(id = "de_proda_summ_wait")

    if (!is.null(cache_result$error)) {
      showNotification(cache_result$error, type = "error", duration = 12)
      return()
    }

    # Carry aggregation route from build step into cache (used by provenance)
    cache_result$aggregation_method <-
      mat_result$input_route %||% "robustSummary"

    # MaxQuant route stats (quant type, contaminant counts, gene source)
    if (!is.null(mat_result$mq_stats))
      cache_result$mq_stats <- mat_result$mq_stats

    # Same for shared-peptide stats on the peptide-aggregating route. Left NULL
    # on routes that aggregate no peptides, so the renderer stays silent there.
    if (!is.null(mat_result$shared_peptide_handling) &&
        !is.na(mat_result$shared_peptide_handling))
      cache_result$shared_peptide_stats <- list(
        shared_peptide_handling     = mat_result$shared_peptide_handling,
        n_peptides_total            = mat_result$n_peptides_total,
        n_peptides_shared           = mat_result$n_peptides_shared,
        n_peptides_removed_shared   = mat_result$n_peptides_removed_shared,
        shared_detection_method     = mat_result$shared_detection_method,
        shared_detection_applicable = mat_result$shared_detection_applicable,
        n_peptides_razor_assigned             = mat_result$n_peptides_razor_assigned,
        n_peptides_dropped_no_unique_evidence = mat_result$n_peptides_dropped_no_unique_evidence,
        razor_basis                           = mat_result$razor_basis)

    # Store in shared cache (ProDA format)
    rv$lf_master_cache <- cache_result
    rv$lf_master_cache_meta <- list(
      pipeline    = pipeline,
      ip_sub      = ip_sub,
      skip_impute = FALSE,
      bait        = NULL,
      igg         = igg_conds,
      groups      = grps,
      pipeline_label = cache_result$pipeline_label,
      n_proteins  = cache_result$n_prots_processed,
      de_engine   = "proda",
      format      = "proda",
      input_mode  = rv$lf_input_mode,
      timestamp   = Sys.time()
    )
    rv$lf_cache_stale <- FALSE

    logger(paste0("DE [ProDA Auto-Summarize]: Cache built \u2014 ",
                  cache_result$n_prots_processed, " proteins. Running ProDA DE..."))
    showNotification(
      paste0("ProDA summarization complete: ", cache_result$n_prots_processed,
             " proteins. Running DE..."),
      type = "message", duration = 5)

    run_proda_from_cache(cache_result)
  }

  # ===========================================================================
  # DATA FILTERS  -" add / remove / clear / render
  # ===========================================================================
  observeEvent(input$de_open_filters, {
    # Resolve available column headers for any data source.
    # For SAINT: use analysis result columns (BFDR, SaintScore, Confidence, etc.).
    # For other modes: use raw PD headers supplemented with any result columns.
    hdrs <- if (identical(de_rv$analysis_mode, "saint") && !is.null(de_rv$result_df)) {
      colnames(de_rv$result_df)
    } else if (!is.null(rv$pd_protein_headers)) {
      result_extra <- if (!is.null(de_rv$result_df))
        setdiff(colnames(de_rv$result_df), rv$pd_protein_headers) else character(0)
      unique(c(rv$pd_protein_headers, result_extra))
    } else if (!is.null(de_rv$plot_data)) {
      colnames(de_rv$plot_data)
    } else if (!is.null(rv$dia_master_cache$protein_matrix) ||
               !is.null(rv$tmt_limma_cache$protein_matrix)) {
      c("GeneName", "Accession", "Log2FC", "AdjPValue", "NegLog10P")
    } else {
      showNotification("Run DE analysis first to enable filters.", type = "warning")
      return()
    }
    showModal(build_filter_modal("de_f", hdrs, de_rv$filter_conditions))
  })

  observeEvent(input$de_f_confirm, {
    col_name  <- input$de_f_column
    op        <- input$de_f_operator
    val       <- trimws(input$de_f_value)
    grp       <- input$de_f_group
    grp_logic <- input$de_f_group_logic
    if (is.null(grp_logic) || !grp_logic %in% c("OR", "AND")) grp_logic <- "OR"
    if (op %in% c("isNA", "notNA")) val <- ""
    if (is.null(col_name) || col_name == "") {
      showNotification("Please select a column.", type = "warning"); return()
    }
    if (!op %in% c("isNA", "notNA") && val == "") {
      showNotification("Please enter a value.", type = "warning"); return()
    }
    # Update logic for all existing conditions in the same group
    if (length(de_rv$filter_conditions) > 0) {
      de_rv$filter_conditions <- lapply(de_rv$filter_conditions, function(f) {
        if (f$group == grp) f$logic <- grp_logic; f
      })
    }
    de_rv$filter_counter <- de_rv$filter_counter + 1L
    de_rv$filter_conditions <- c(de_rv$filter_conditions, list(list(
      id = de_rv$filter_counter, group = grp, column = col_name,
      op = op, value = val, logic = grp_logic
    )))
    op_label <- switch(op, "isNA" = "is NA", "notNA" = "is not NA",
                        "contains" = "contains", "not_contains" = "not contains", op)
    logger(paste0("DE Filter: Added [", grp, " (", grp_logic, ")] ",
                  col_name, " ", op_label,
                  if (nzchar(val)) paste0(" ", val) else ""))
    removeModal()
  })

  # Remove single filter condition
  observeEvent(lapply(seq_len(100), function(i) input[[paste0("de_frm_", i)]]), {
    for (i in seq_along(de_rv$filter_conditions)) {
      btn_id <- paste0("de_frm_", de_rv$filter_conditions[[i]]$id)
      btn_val <- input[[btn_id]]
      if (!is.null(btn_val) && btn_val > 0) {
        cond <- de_rv$filter_conditions[[i]]
        attr_name <- paste0(".de_frm_seen_", cond$id)
        last_seen <- isolate(de_rv[[attr_name]])
        if (is.null(last_seen) || btn_val > last_seen) {
          de_rv[[attr_name]] <- btn_val
          de_rv$filter_conditions <- de_rv$filter_conditions[
            vapply(de_rv$filter_conditions, function(f) f$id != cond$id, logical(1))
          ]
          logger(paste0("DE Filter: Removed #", cond$id))
          break
        }
      }
    }
  })

  observeEvent(input$de_clear_filters, {
    de_rv$filter_conditions <- list()
    logger("DE Filter: All cleared")
  })

  output$de_active_filters <- renderUI({
    conds <- de_rv$filter_conditions
    if (length(conds) == 0)
      return(tags$p(class = "text-muted small", "No active filters."))
    render_filter_badges(conds, "de_frm_")
  })

  # ===========================================================================
  # HIGHLIGHT RULES  -" add / remove / clear / render
  # Same visual style as filters but applied as highlighting, not filtering
  # ===========================================================================
  observeEvent(input$de_add_highlight, {
    # Only offer columns that exist in the result/plot data (highlighting
    # is evaluated against that data, so raw-CSV-only columns won't work).
    result_cols <- NULL
    if (!is.null(de_rv$plot_data)) {
      result_cols <- colnames(de_rv$plot_data)
    } else if (!is.null(de_rv$result_df)) {
      result_cols <- colnames(de_rv$result_df)
    }
    all_choices <- unique(result_cols)
    req(length(all_choices) > 0)

    showModal(modalDialog(
      title = tags$span(icon("highlighter"), " Add Highlight Rule"),
      tagList(
        helpText(class = "small text-muted mb-2",
          "Highlight proteins matching this condition on the plot. ",
          "Multiple rules are combined with OR."),
        tags$div(class = "mb-3",
          tags$label(class = "form-label fw-bold", "Column"),
          tags$div(style = "max-height: 250px; overflow-y: auto; border: 1px solid #eee; padding: 6px; border-radius: 4px;",
            radioButtons("de_hl_column", label = NULL,
                         choices  = setNames(all_choices, all_choices),
                         selected = all_choices[1])
          )
        ),
        fluidRow(
          column(6, selectInput("de_hl_operator", "Operator",
                    choices = c(">" = ">", ">=" = ">=", "<" = "<", "<=" = "<=",
                                "==" = "==", "!=" = "!=",
                                "contains" = "contains",
                                "not contains" = "not_contains"))),
          column(6, textInput("de_hl_value", "Value", placeholder = "e.g. 0.01"))
        )
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("de_hl_confirm", "Add Highlight", class = "btn-warning",
                     icon = icon("highlighter"))
      )
    ))
  })

  observeEvent(input$de_hl_confirm, {
    col_name <- input$de_hl_column
    op       <- input$de_hl_operator
    val      <- trimws(input$de_hl_value)
    if (is.null(col_name) || col_name == "" || val == "") {
      showNotification("Please select a column and enter a value.", type = "warning")
      return()
    }
    de_rv$highlight_counter <- de_rv$highlight_counter + 1L
    de_rv$highlight_rules <- c(de_rv$highlight_rules, list(list(
      id = de_rv$highlight_counter, column = col_name, op = op, value = val
    )))
    de_rv$render_count <- isolate(de_rv$render_count) + 1
    logger(paste0("DE Highlight: Added ", col_name, " ", op, " ", val))
    removeModal()
  })

  # Remove single highlight rule
  observeEvent(lapply(seq_len(100), function(i) input[[paste0("de_hlrm_", i)]]), {
    for (i in seq_along(de_rv$highlight_rules)) {
      btn_id <- paste0("de_hlrm_", de_rv$highlight_rules[[i]]$id)
      btn_val <- input[[btn_id]]
      if (!is.null(btn_val) && btn_val > 0) {
        rule <- de_rv$highlight_rules[[i]]
        attr_name <- paste0(".de_hlrm_seen_", rule$id)
        last_seen <- isolate(de_rv[[attr_name]])
        if (is.null(last_seen) || btn_val > last_seen) {
          de_rv[[attr_name]] <- btn_val
          de_rv$highlight_rules <- de_rv$highlight_rules[
            vapply(de_rv$highlight_rules, function(r) r$id != rule$id, logical(1))
          ]
          de_rv$render_count <- isolate(de_rv$render_count) + 1
          logger(paste0("DE Highlight: Removed #", rule$id))
          break
        }
      }
    }
  })

  observeEvent(input$de_clear_highlights, {
    de_rv$highlight_rules <- list()
    de_rv$render_count <- isolate(de_rv$render_count) + 1
    logger("DE Highlight: All cleared")
  })

  output$de_active_highlights <- renderUI({
    rules <- de_rv$highlight_rules
    if (length(rules) == 0)
      return(tags$p(class = "text-muted small", "No active highlight rules."))

    rule_tags <- lapply(rules, function(r) {
      op_label <- if (r$op == "contains") "contains"
                  else if (r$op == "not_contains") "not contains"
                  else r$op
      tags$div(class = "d-flex align-items-center gap-1 mb-1",
        tags$span(class = "badge bg-warning text-dark",
                  style = "font-size: 0.8em; word-break: break-all;",
                  paste0(r$column, " ", op_label, " ", r$value)),
        actionButton(paste0("de_hlrm_", r$id), label = NULL,
                     icon = icon("xmark"),
                     class = "btn-outline-danger btn-sm",
                     style = "padding: 1px 5px; font-size: 0.7em;")
      )
    })

    # Show OR between rules
    or_separated <- list()
    for (j in seq_along(rule_tags)) {
      or_separated <- c(or_separated, list(rule_tags[[j]]))
      if (j < length(rule_tags))
        or_separated <- c(or_separated, list(
          tags$span(class = "badge bg-warning text-dark",
                    style = "font-size: 0.75em;", "OR")
        ))
    }

    tags$div(
      style = "background: #fef5e7; border-left: 3px solid #f39c12; padding: 6px 8px; border-radius: 4px;",
      tagList(or_separated)
    )
  })

  # ===========================================================================
  # Gene search highlight
  # ===========================================================================
  observeEvent(input$de_gene_search_btn, {
    gene <- trimws(input$de_gene_search)
    de_rv$highlight_gene <- if (nchar(gene) == 0) NULL else gene
    de_rv$render_count <- isolate(de_rv$render_count) + 1
  })

  # ===========================================================================
  # Render DE Plot (branches by mode)
  # ===========================================================================
  output$de_plot <- renderPlotly({
    mode <- de_rv$analysis_mode
    de_rv$render_count   # reactive dependency
    req(mode)

    if (identical(mode, "proda")) {
      render_proda_plot()
    } else {
      render_volcano_plot()
    }
  })

  # ---------------------------------------------------------------------------
  # ProDA plot: classic DE volcano
  # ---------------------------------------------------------------------------
  render_proda_plot <- function() {
    req(de_rv$result_df)
    df <- de_rv$result_df

    bait      <- de_rv$bait_name

    # Apply data filters  -" supplement raw PD data with analysis result columns so
    # users can filter by AdjPvalue, LowEvidence, etc. as well as raw upload columns.
    filter_conds <- de_rv$filter_conditions
    if (length(filter_conds) > 0) {
      pdata <- parse_protein_data()
      if (!is.null(pdata)) {
        filter_df <- pdata$protein_df
        if (!is.null(de_rv$result_df) && "Accession" %in% colnames(de_rv$result_df)) {
          res_extra <- setdiff(colnames(de_rv$result_df), colnames(filter_df))
          if (length(res_extra) > 0)
            filter_df <- merge(filter_df,
              de_rv$result_df[, c("Accession", res_extra), drop = FALSE],
              by = "Accession", all.x = TRUE)
        }
        keep_mask <- apply_filter_conditions(filter_df, filter_conds, "and_or")
        keep_acc  <- filter_df[["Accession"]][keep_mask]
        df <- df[df$Accession %in% keep_acc, ]
      }
    }

    # Show/hide capped
    show_capped <- isTRUE(input$de_show_capped)
    if (!show_capped) {
      df <- df[df$AdjPvalue != 0 & !is.na(df$AdjPvalue), ]
    }

    # Y-axis
    is_lod_mode <- identical(de_rv$analysis_mode, "dia") && isTRUE(de_rv$dia_use_lod)
    adj_p_vals <- df$AdjPvalue
    adj_p_vals[is.na(adj_p_vals)] <- 1
    df$y_val <- -log10(pmax(adj_p_vals, 1e-300))
    y_title  <- if (is_lod_mode) "-Log10(p) [Welch t-test, s0=0.1]" else "-Log10(Adj. P-value) [BH]"

    # Size by condition intensity
    # SAINT result_df uses Log2AvgIntensity; other modes use AvgIntensity.
    is_saint_mode <- identical(de_rv$analysis_mode, "saint")
    int_col      <- if (is_saint_mode && "Log2AvgIntensity" %in% colnames(df))
                      "Log2AvgIntensity" else "AvgIntensity"
    int_vals_raw <- df[[int_col]]
    valid_int    <- int_vals_raw[is.finite(int_vals_raw)]
    if (length(valid_int) > 0) {
      int_range <- range(valid_int, na.rm = TRUE)
      if (diff(int_range) > 0) {
        df$pt_size <- 0.75 * (5 + 15 *
          (pmin(pmax(int_vals_raw, int_range[1]), int_range[2]) -
           int_range[1]) / diff(int_range))
      } else df$pt_size <- 7.5
    } else df$pt_size <- 7.5
    df$pt_size[is.na(df$pt_size)] <- 6

    # Hover text
    pval_label  <- if (is_saint_mode) "BFDR" else if (is_lod_mode) "p-value (Welch)" else "Adj. P-value (BH)"
    pval_vals   <- if (is_saint_mode && "BFDR" %in% colnames(df))
                      signif(df$BFDR, 3) else signif(df$AdjPvalue, 3)
    saint_hover_extra <- if (is_saint_mode && "SaintScore" %in% colnames(df))
                           paste0("<br>SAINT Score: ", round(df$SaintScore, 3)) else ""
    lod_hover_extra <- if (is_lod_mode && "perm_fdr" %in% colnames(df))
                         paste0("<br>Perm. FDR: ", signif(df$perm_fdr, 3)) else ""
    pi_hover_extra <- if ("PiScore" %in% colnames(df))
                        paste0("<br>\u03c0-score: ", round(df$PiScore, 2)) else ""
    df$hover <- paste0(
      "<b>", df$GeneName, "</b> (", df$Accession, ")",
      "<br>Log2 FC: ", df$Log2FC,
      "<br>", pval_label, ": ", pval_vals,
      lod_hover_extra,
      saint_hover_extra,
      pi_hover_extra,
      "<br>Intensity: ", round(int_vals_raw, 1),
      "<br>Condition: ", df$CondDetected,
      " | Control: ", df$CtrlDetected
    )
    # LowEvidence (proDA runs only): fewer than two observations in at least one arm
    le_mask <- if ("LowEvidence" %in% colnames(df))
                 df$LowEvidence %in% c(TRUE, "TRUE") else rep(FALSE, nrow(df))
    df$hover[le_mask] <- paste0(df$hover[le_mask],
      "<br><i>Low evidence: &lt;2 values in a group</i>")

    # --- Determine highlight mode ---
    highlight_gene  <- de_rv$highlight_gene
    exact_match     <- isTRUE(input$de_exact_match)
    highlight_rules <- de_rv$highlight_rules

    # Compute highlight mask
    hl_mask <- rep(FALSE, nrow(df))
    if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      if (exact_match) hl_mask <- hl_mask | (toupper(df$GeneName) == toupper(highlight_gene))
      else hl_mask <- hl_mask | grepl(highlight_gene, df$GeneName, ignore.case = TRUE)
    }
    if (length(highlight_rules) > 0) {
      hl_mask <- hl_mask | apply_filter_conditions(df, highlight_rules, "or")
    }

    has_highlights <- any(hl_mask)
    p <- plot_ly()

    if (has_highlights) {
      # Highlight mode: faded background, prominent highlighted
      non_match <- df[!hl_mask, ]
      if (nrow(non_match) > 0)
        p <- add_trace(p, data = non_match, x = ~Log2FC, y = ~y_val,
          type = "scatter", mode = "markers", name = "Other",
          marker = list(size = 6, color = "#d4d4d4", opacity = 0.35,
                        line = list(width = 0.5, color = "#999")),
          hoverinfo = "text", text = ~hover)

      match_df <- df[hl_mask, ]
      if (nrow(match_df) > 0) {
        p <- add_trace(p, data = match_df, x = ~Log2FC, y = ~y_val,
          type = "scatter", mode = "markers", name = "Highlighted",
          marker = list(size = 12, color = "#e74c3c", opacity = 1,
                        line = list(width = 1, color = "#333")),
          hoverinfo = "text", text = ~hover)
        p <- add_trace(p, data = match_df, x = ~Log2FC, y = ~y_val,
          type = "scatter", mode = "text", text = ~GeneName,
          textposition = "top center",
          textfont = list(size = 11, color = "#c0392b", family = "Arial Black"),
          showlegend = FALSE, hoverinfo = "none")
      }
    } else {
      # Classic DE volcano: simple red/blue/gray by FC + adj.p thresholds.
      # LowEvidence proteins get their own grey-square group wherever they land,
      # so a red/blue point is never one whose estimate rests on <2 values in a group.
      adj_p_for_col <- df$AdjPvalue; adj_p_for_col[is.na(adj_p_for_col)] <- 1
      up_sig  <- adj_p_for_col < 0.05 & df$Log2FC >  1 & !le_mask
      dn_sig  <- adj_p_for_col < 0.05 & df$Log2FC < -1 & !le_mask
      ns_mask <- !up_sig & !dn_sig & !le_mask
      vol_groups <- list(
        list(mask = up_sig,  col = "#e74c3c", name = "Up (log\u2082FC\u00a0>\u00a01, adj.p\u00a0<\u00a00.05)",   opac = 0.85, symbol = "circle"),
        list(mask = dn_sig,  col = "#2980b9", name = "Down (log\u2082FC\u00a0<\u00a0\u22121, adj.p\u00a0<\u00a00.05)", opac = 0.85, symbol = "circle"),
        list(mask = ns_mask, col = "#bdc3c7", name = "Not significant",                        opac = 0.35, symbol = "circle"),
        list(mask = le_mask, col = "#bdc3c7", name = "Low evidence (<2 values in a group)",    opac = 0.6,  symbol = "square")
      )
      for (vg in vol_groups) {
        sub <- df[vg$mask, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- add_trace(p, data = sub, x = ~Log2FC, y = ~y_val,
          type = "scatter", mode = "markers", name = vg$name,
          marker = list(size = 7, color = vg$col, opacity = vg$opac, symbol = vg$symbol,
                        line = list(width = 0.5, color = "#555")),
          hoverinfo = "text", text = ~hover)
      }
    }

    # Threshold lines
    {
      x_max <- max(abs(df$Log2FC[is.finite(df$Log2FC)]), na.rm = TRUE) * 1.15
      x_max <- max(x_max, 2)
      y_max <- max(df$y_val[is.finite(df$y_val)], na.rm = TRUE) * 1.05

      # Perseus LOD mode: horizontal PermFDR = 0.05 threshold line.
      # The Y-axis is -log10(p_s0) where p_s0 comes from t_s0adj.
      # Significance is determined solely by |t_s0adj| >= delta_star, so the
      # boundary in this coordinate space is a horizontal line at the
      # -log10(p_s0) value corresponding to t_s0adj = delta_star.
      delta_star <- de_rv$lod_delta_star
      df_rep     <- de_rv$lod_df_rep
      if (is_lod_mode && !is.null(delta_star) && is.finite(delta_star) && delta_star > 0) {
        y_thresh <- -log10(2 * pt(-delta_star, df = max(df_rep, 1)))
        if (is.finite(y_thresh)) {
          ref_shapes <- list(
            list(type = "line", x0 = 0, x1 = 1, xref = "paper",
                 y0 = y_thresh, y1 = y_thresh,
                 line = list(dash = "dash", color = "black", width = 2))
          )
          ref_annot <- list(
            list(x = 1, y = y_thresh, xref = "paper", yref = "y",
                 text = "PermFDR = 0.05", showarrow = FALSE,
                 font = list(size = 10, color = "black"), xanchor = "right",
                 yanchor = "bottom")
          )
        } else {
          ref_shapes <- list()
          ref_annot  <- list()
        }
      } else {
        # Standard classic DE: horizontal adj.p=0.05 line + vertical FC=±1 lines
        y_sig   <- -log10(0.05)
        ref_shapes <- list(
          list(type = "line", x0 = 0, x1 = 1, xref = "paper",
               y0 = y_sig, y1 = y_sig,
               line = list(dash = "dash", color = "#555", width = 1.5)),
          list(type = "line", x0 = 1, x1 = 1, yref = "paper", y0 = 0, y1 = 1,
               line = list(dash = "dot", color = "#555", width = 1)),
          list(type = "line", x0 = -1, x1 = -1, yref = "paper", y0 = 0, y1 = 1,
               line = list(dash = "dot", color = "#555", width = 1))
        )
        ref_annot <- list(
          list(x = 1, y = y_sig, xref = "paper", yref = "y",
               text = "adj.p\u00a0=\u00a00.05", showarrow = FALSE,
               font = list(size = 10, color = "#555"), xanchor = "right",
               yanchor = "bottom")
        )
      }
    }

    type_label <- "Classic DE"
    engine_label <- switch(de_rv$analysis_mode,
                           "proda" = "ProDA", "saint" = "SAINTexpress",
                           "welch" = "Peptide BG Z-Test", "Unknown")
    p <- layout(p,
      title = paste0("Differential Expression (", engine_label, ", ", type_label, ")",
                     if (!is.null(bait) && bait != "Bait") paste0(" \u2014 ", bait) else ""),
      xaxis = list(title = "Log2 Fold Change"),
      yaxis = list(title = y_title),
      shapes = ref_shapes, annotations = ref_annot,
      legend = list(orientation = "h", y = -0.15),
      hovermode = "closest"
    )
    p
  }

  # ---------------------------------------------------------------------------
  # MODES B/C plot: Standard volcano
  # ---------------------------------------------------------------------------
  render_volcano_plot <- function() {
    req(de_rv$plot_data)
    df <- de_rv$plot_data

    # Detect Perseus LOD mode
    is_lod_vol <- identical(de_rv$analysis_mode, "dia") &&
                  !is.null(de_rv$lod_delta_star) && is.finite(de_rv$lod_delta_star)

    # Apply data filters
    filter_conds <- de_rv$filter_conditions
    if (length(filter_conds) > 0 && !is.null(de_rv$protein_df)) {
      keep_mask <- apply_filter_conditions(de_rv$protein_df, filter_conds, "and_or")
      keep_acc  <- de_rv$protein_df[["Accession"]][keep_mask]
      df <- df[df$Accession %in% keep_acc, ]
    } else if (length(filter_conds) > 0) {
      keep_mask <- apply_filter_conditions(df, filter_conds, "and_or")
      df <- df[keep_mask, ]
    }

    show_capped <- isTRUE(input$de_show_capped)
    if (!show_capped) {
      df <- df[!(df$IsCapped | (!is.na(df$AdjPValue) & df$AdjPValue == 0)), ]
    } else {
      capped_idx <- which(df$IsCapped | (!is.na(df$AdjPValue) & df$AdjPValue == 0))
      if (length(capped_idx) > 0) {
        set.seed(42)
        fc_range <- diff(range(df$Log2FC, na.rm = TRUE))
        p_range  <- diff(range(df$NegLog10P, na.rm = TRUE))
        df$Log2FC[capped_idx]    <- df$Log2FC[capped_idx]    + runif(length(capped_idx), -fc_range * 0.02, fc_range * 0.02)
        df$NegLog10P[capped_idx] <- df$NegLog10P[capped_idx] + runif(length(capped_idx), -p_range * 0.02,  p_range * 0.02)
      }
    }

    # Highlight logic
    highlight_gene  <- de_rv$highlight_gene
    exact_match     <- isTRUE(input$de_exact_match)
    highlight_rules <- de_rv$highlight_rules

    hl_mask <- rep(FALSE, nrow(df))
    if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      if (exact_match) hl_mask <- hl_mask | (toupper(df$GeneName) == toupper(highlight_gene))
      else hl_mask <- hl_mask | grepl(highlight_gene, df$GeneName, ignore.case = TRUE)
    }
    if (length(highlight_rules) > 0)
      hl_mask <- hl_mask | apply_filter_conditions(df, highlight_rules, "or")

    has_highlights <- any(hl_mask)

    # Mark imputed
    imputed_mask <- if ("IsImputed" %in% colnames(df))
      !is.na(df$IsImputed) & df$IsImputed else rep(FALSE, nrow(df))

    is_proda_vol <- identical(de_rv$analysis_mode, "proda")
    df$pt_symbol <- ifelse(imputed_mask, "triangle-up", "circle")

    # Hover text
    is_saint_vol <- identical(de_rv$analysis_mode, "saint")
    pval_hover_label <- if (is_saint_vol) "BFDR"
                        else if (is_lod_vol) "p (Welch, s0=0.1)"
                        else "Adj. P-Value (BH)"
    saint_score_line <- if (is_saint_vol && "Confidence" %in% colnames(df)) {
      ss_col <- if (!is.null(de_rv$result_df) && "SaintScore" %in% colnames(de_rv$result_df))
        de_rv$result_df$SaintScore[match(df$Accession, de_rv$result_df$Accession)]
      else rep(NA_real_, nrow(df))
      paste0("<br><b>SAINT Score:</b> ", round(ss_col, 3))
    } else rep("", nrow(df))
    perm_fdr_line <- if (is_lod_vol && "PermFDR" %in% colnames(df))
      paste0("<br><b>Perm. FDR:</b> ", signif(df$PermFDR, 3))
    else rep("", nrow(df))
    df$hover <- paste0(
      "<b>Gene:</b> ", df$GeneName,
      ifelse(df$IsCapped, " <i>(capped by PD)</i>",
             ifelse(!is.na(df$AdjPValue) & df$AdjPValue == 0, " <i>(p-value floored at 0)</i>", "")),
      ifelse(imputed_mask,
             if (is_proda_vol)
               " <i>(one group not detected \u2014 dropout-model estimate)</i>"
             else " <i>(imputed)</i>", ""),
      "<br><b>Accession:</b> ", df$Accession,
      "<br><b>", pval_hover_label, ":</b> ", signif(df$AdjPValue, 4),
      perm_fdr_line,
      saint_score_line,
      "<br><b>Log2 FC:</b> ", round(df$Log2FC, 3)
    )

    p <- plot_ly()

    if (has_highlights) {
      non_match <- df[!hl_mask, ]
      if (nrow(non_match) > 0)
        p <- add_trace(p, data = non_match, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "markers", name = "Other",
          marker = list(size = 6, color = "#d4d4d4", opacity = 0.35),
          hoverinfo = "text", text = ~hover)

      match_df <- df[hl_mask, ]
      if (nrow(match_df) > 0) {
        p <- add_trace(p, data = match_df, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "markers", name = "Highlighted",
          marker = list(size = 12, color = "#e74c3c", opacity = 1,
                        line = list(width = 1, color = "#333")),
          hoverinfo = "text", text = ~hover)
        p <- add_trace(p, data = match_df, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "text", text = ~GeneName,
          textposition = "top center",
          textfont = list(size = 11, color = "#c0392b", family = "Arial Black"),
          showlegend = FALSE, hoverinfo = "none")
      }
    } else if (identical(de_rv$analysis_mode, "saint") && "Confidence" %in% colnames(df)) {
      # SAINT: color dots by SaintScore-derived Confidence tier
      conf_colors  <- c("Exclusive" = "#9b59b6", "High" = "#e74c3c",
                        "Medium" = "#f39c12", "Low" = "#bdc3c7")
      conf_opacity <- c("Exclusive" = 0.9, "High" = 0.8, "Medium" = 0.7, "Low" = 0.4)
      # Marker class for T12 traces; control-absent wins over SingleReplicate
      vol_ctrl_absent <- if ("FC_Origin" %in% colnames(df))
        !is.na(df$FC_Origin) & df$FC_Origin == "control-absent"
      else rep(FALSE, nrow(df))
      vol_single_rep  <- if ("SingleReplicate" %in% colnames(df))
        !is.na(df$SingleReplicate) & df$SingleReplicate & !vol_ctrl_absent
      else rep(FALSE, nrow(df))
      df$MarkerClass <- ifelse(vol_ctrl_absent, "control-absent",
                        ifelse(vol_single_rep, "single-replicate", "measured"))
      for (conf in c("Exclusive", "High", "Medium", "Low")) {
        sub <- df[df$Confidence == conf & df$MarkerClass == "measured", ]
        if (nrow(sub) == 0) next
        p <- add_trace(p, data = sub, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "markers",
          name = switch(conf,
            "Exclusive" = "Exclusive (bait-only)",
            "High"      = "High (SAINT \u2265 0.80)",
            "Medium"    = "Medium (SAINT 0.50\u20130.80)",
            "Low"       = "Low (SAINT < 0.50)"),
          marker = list(size = 7, opacity = conf_opacity[conf],
                        color = conf_colors[conf]),
          hoverinfo = "text", text = ~hover)
      }
      # Separate traces so the legend can toggle them off independently
      for (mc in c("control-absent", "single-replicate")) {
        sub <- df[df$MarkerClass == mc, ]
        if (nrow(sub) == 0) next
        p <- add_trace(p, data = sub, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "markers",
          name = if (identical(mc, "control-absent"))
                   "Floor-substituted (control-absent)"
                 else "Single replicate (bait arm)",
          marker = list(
            size    = if (identical(mc, "control-absent")) 9 else 8,
            symbol  = if (identical(mc, "control-absent")) "triangle-up" else "star",
            opacity = 0.9,
            color   = unname(conf_colors[sub$Confidence]),
            line    = list(width = 1, color = "#333")),
          hoverinfo = "text", text = ~hover)
      }
    } else {
      # Normal significance-based coloring
      # LOD mode: use PermFDR <= 0.05; standard mode: use AdjPValue < 0.05
      if (is_lod_vol && "PermFDR" %in% colnames(df)) {
        sig_mask <- !is.na(df$PermFDR) & df$PermFDR <= 0.05
      } else {
        sig_mask <- !is.na(df$AdjPValue) & df$AdjPValue < 0.05
      }
      df$Category <- ifelse(
        sig_mask & df$Log2FC >  1, "Up-regulated (sig)",
        ifelse(sig_mask & df$Log2FC < -1, "Down-regulated (sig)",
               ifelse(abs(df$Log2FC) >= 1, "Log2FC > 1 only",
                      ifelse(sig_mask, "Sig. only", "Not significant"))))

      cat_colors <- c(
        "Up-regulated (sig)"   = "#e74c3c",
        "Down-regulated (sig)" = "#3498db",
        "Log2FC > 1 only"      = "#f39c12",
        "P < 0.05 only"        = "#95a5a6",
        "Sig. only"            = "#95a5a6",
        "Not significant"      = "#d4d4d4"
      )

      for (cat in unique(df$Category)) {
        # ProDA: marked rows are drawn in their own trace below (SAINT
        # convention), so the legend can toggle them. Every other engine gets
        # the same rows and marker list it got before.
        sub <- df[df$Category == cat & !(is_proda_vol & imputed_mask), ]
        if (is_proda_vol && nrow(sub) == 0) next
        .vol_marker <- list(size = 7, opacity = 0.7, color = cat_colors[cat])
        if (is_proda_vol) .vol_marker$symbol <- sub$pt_symbol
        p <- add_trace(p, data = sub, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "markers", name = cat,
          marker = .vol_marker,
          hoverinfo = "text", text = ~hover)
      }
      # Zero-detection rows: one toggleable triangle trace, coloured by category.
      if (is_proda_vol && any(imputed_mask)) {
        sub <- df[imputed_mask, ]
        p <- add_trace(p, data = sub, x = ~Log2FC, y = ~NegLog10P,
          type = "scatter", mode = "markers",
          name = "\u25b2 One group not detected",
          marker = list(size = 8, symbol = "triangle-up", opacity = 0.9,
                        color = unname(cat_colors[sub$Category]),
                        line = list(width = 1, color = "#333")),
          hoverinfo = "text", text = ~hover)
      }
    }

    mode_label <- if (identical(de_rv$analysis_mode, "imputation")) "Imputation + limma"
                  else if (identical(de_rv$analysis_mode, "saint"))  "SAINTexpress"
                  else if (identical(de_rv$analysis_mode, "welch"))  "Peptide-Level Background Z-Test (Iterative 3\u03c3 Clipping)"
                  else if (is_lod_vol) "Perseus LOD + Permutation FDR"
                  else "Proteome Discoverer"

    y_axis_title <- if (identical(de_rv$analysis_mode, "saint"))
      "\u2212log\u2081\u2080 (BFDR)"
      else if (is_lod_vol)
        "-Log10(p) [Welch t-test, s0=0.1]"
      else
        "\u2212log\u2081\u2080 (Adj. P\u2010Value [BH])"

    # Vertical FC threshold lines (always shown)
    vol_shapes <- list(
      list(type = "line", x0 = -1, x1 = -1, y0 = 0, y1 = 1, yref = "paper",
           line = list(dash = "dash", color = "#999", width = 1)),
      list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = 1, yref = "paper",
           line = list(dash = "dash", color = "#999", width = 1))
    )
    if (!is_lod_vol) {
      # Standard mode: horizontal Adj.P = 0.05 line
      vol_shapes[[3]] <- list(
        type = "line", x0 = 0, x1 = 1,
        y0 = -log10(0.05), y1 = -log10(0.05), xref = "paper",
        line = list(dash = "dash", color = "#999", width = 1)
      )
    }

    p <- layout(p,
      title = list(text = paste0("Volcano Plot (", mode_label, ")"), font = list(size = 18)),
      xaxis = list(title = "Log2 Fold Change",
                   zeroline = TRUE, zerolinecolor = "#999", zerolinewidth = 1),
      yaxis = list(title = y_axis_title, zeroline = FALSE),
      shapes = vol_shapes,
      legend = list(orientation = "h", y = -0.15),
      hovermode = "closest"
    )

    # LOD mode: horizontal PermFDR = 0.05 threshold line.
    # Y-axis = -log10(p_s0adj), significance boundary = -log10(p at |t_s0| = delta_star).
    if (is_lod_vol) {
      delta_star <- de_rv$lod_delta_star
      df_rep     <- de_rv$lod_df_rep
      if (!is.null(delta_star) && is.finite(delta_star) && delta_star > 0 && !is.null(df_rep)) {
        y_thresh <- -log10(2 * pt(-delta_star, df = max(df_rep, 1)))
        if (is.finite(y_thresh)) {
          p <- layout(p, shapes = c(
            p$x$layout$shapes,
            list(list(
              type = "line", x0 = 0, x1 = 1, xref = "paper",
              y0 = y_thresh, y1 = y_thresh,
              line = list(dash = "dash", color = "black", width = 2)
            ))
          ),
          annotations = list(list(
            x = 1, y = y_thresh, xref = "paper", yref = "y",
            text = "PermFDR = 0.05", showarrow = FALSE,
            font = list(size = 10, color = "black"),
            xanchor = "right", yanchor = "bottom"
          )))
        }
      }
    }

    # ── Welch contingency subtitle ────────────────────────────────────────────
    if (identical(de_rv$analysis_mode, "welch")) {
      p <- p %>% layout(
        annotations = list(list(
          text      = paste0("Statistical Engine: Peptide-Level Background Z-Test ",
                             "(Iterative 3\u03c3 Clipping \u2014 Contingency Mode)"),
          showarrow = FALSE,
          xref = "paper", yref = "paper",
          x = 0.5, y = -0.11,
          xanchor = "center", yanchor = "top",
          font = list(size = 10, color = "#e67e22")
        ))
      )
    }

    p
  }

  # ===========================================================================
  # Data table
  # ===========================================================================
  output$de_data_table <- DT::renderDT({
    mode <- de_rv$analysis_mode

    if (mode %in% c("proda", "saint", "welch")) {
      req(de_rv$result_df)
      df <- de_rv$result_df

      # Apply filters  -" for SAINT: de_rv$protein_df = result_df (BFDR, SaintScore, etc.)
      # For other modes: supplement raw PD data with analysis result columns.
      filter_conds <- de_rv$filter_conditions
      if (length(filter_conds) > 0) {
        if (!is.null(de_rv$protein_df)) {
          keep_mask <- apply_filter_conditions(de_rv$protein_df, filter_conds, "and_or")
          keep_acc  <- de_rv$protein_df[["Accession"]][keep_mask]
          df <- df[df$Accession %in% keep_acc, ]
        } else {
          pdata <- parse_protein_data()
          if (!is.null(pdata)) {
            filter_df <- pdata$protein_df
            if (!is.null(de_rv$result_df) && "Accession" %in% colnames(de_rv$result_df)) {
              res_extra <- setdiff(colnames(de_rv$result_df), colnames(filter_df))
              if (length(res_extra) > 0)
                filter_df <- merge(filter_df,
                  de_rv$result_df[, c("Accession", res_extra), drop = FALSE],
                  by = "Accession", all.x = TRUE)
            }
            keep_mask <- apply_filter_conditions(filter_df, filter_conds, "and_or")
            keep_acc  <- filter_df[["Accession"]][keep_mask]
            df <- df[df$Accession %in% keep_acc, ]
          }
        }
      }

      # Column names differ between SAINT (Log2AvgIntensity, BFDR, actual rep names)
      # and other modes (AvgIntensity, AdjPvalue, generic Ctrl_RepN/Cond_RepN).
      if (identical(mode, "saint")) {
        rep_cols <- intersect(
          c(de_rv$ctrl_cols %||% character(0), de_rv$cond_cols %||% character(0)),
          colnames(df))
        show_cols <- c("GeneName", "Accession", "Log2FC", "SaintScore", "BFDR",
                       rep_cols, "CondDetected", "CtrlDetected",
                       "Log2AvgIntensity", "PiScore", "Confidence")
      } else {
        rep_cols  <- sort(grep("^(Ctrl_Rep|Cond_Rep)", colnames(df), value = TRUE))
        show_cols <- c("GeneName", "Accession", "Log2FC", "AdjPvalue",
                       rep_cols, "CondDetected", "CtrlDetected",
                       "AvgIntensity", "PiScore", "LowEvidence",
                       "NObs", "NApprox")  # ProDA only; absent (dropped by intersect) on welch
      }
      df_show <- df[, intersect(show_cols, colnames(df)), drop = FALSE]
      names(df_show)[names(df_show) == "AdjPvalue"]        <- "Adj.Pvalue"
      names(df_show)[names(df_show) == "Log2AvgIntensity"] <- "Log2 Avg Intensity"

      # Default order: SAINT by SaintScore (high first); others by adj. p (low first)
      sort_col <- if (identical(mode, "saint")) "SaintScore" else "Adj.Pvalue"
      sort_dir <- if (identical(mode, "saint")) "desc" else "asc"
      sort_idx <- which(colnames(df_show) == sort_col) - 1
      if (length(sort_idx) == 0) { sort_idx <- 2; sort_dir <- "desc" }

      dt <- DT::datatable(df_show, rownames = FALSE, filter = "top",
        options = list(pageLength = 25, scrollX = TRUE,
                       order = list(list(sort_idx, sort_dir)),
                       dom = "Bfrtip", buttons = c("csv", "excel")),
        extensions = "Buttons", class = "compact stripe hover"
      )
      if ("Confidence" %in% colnames(df_show))   # SAINT only (SaintScore tiers)
        dt <- DT::formatStyle(dt, "Confidence",
          backgroundColor = DT::styleEqual(
            c("Exclusive", "High", "Medium", "Low"),
            c("#ede7f6", "#fde8e8", "#fef5e7", "#f5f5f5")),
          fontWeight = "bold"
        )
      dt

    } else {
      req(de_rv$plot_data)
      df <- de_rv$plot_data

      # Apply filters
      filter_conds <- de_rv$filter_conditions
      if (length(filter_conds) > 0 && !is.null(de_rv$protein_df)) {
        keep_mask <- apply_filter_conditions(de_rv$protein_df, filter_conds, "and_or")
        keep_acc  <- de_rv$protein_df[["Accession"]][keep_mask]
        df <- df[df$Accession %in% keep_acc, ]
      }

      # Add regulation column
      df$Regulation <- ifelse(
        !is.na(df$Log2FC) & df$AdjPValue < 0.05 & df$Log2FC >  1, "Up-regulated",
        ifelse(!is.na(df$Log2FC) & df$AdjPValue < 0.05 & df$Log2FC < -1, "Down-regulated",
               "Not significant"))

      show_cols <- c("GeneName", "Accession", "Log2FC", "AdjPValue", "NegLog10P",
                     "Regulation")
      if ("IsCapped" %in% colnames(df))
        show_cols <- c(show_cols, "IsCapped")
      if ("IsImputed" %in% colnames(df))
        show_cols <- c(show_cols, "IsImputed")

      df_show <- df[, intersect(show_cols, colnames(df)), drop = FALSE]
      df_show <- df_show[order(df_show$AdjPValue), ]

      DT::datatable(df_show, rownames = FALSE, filter = "top",
        options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
                       buttons = c("csv", "excel")),
        extensions = "Buttons", class = "compact stripe hover"
      )
    }
  })

  # ===========================================================================
  # Status summary
  # ===========================================================================
  output$de_status <- renderPrint({
    mode <- de_rv$analysis_mode

    if (mode == "proda" && !is.null(de_rv$result_df)) {
      df   <- de_rv$result_df
      bait <- de_rv$bait_name
      # Same categories as the Classic DE volcano
      le   <- df$LowEvidence %in% TRUE
      adjp <- df$AdjPvalue; adjp[is.na(adjp)] <- 1
      n_up   <- sum(!le & adjp < 0.05 & df$Log2FC >  1, na.rm = TRUE)
      n_down <- sum(!le & adjp < 0.05 & df$Log2FC < -1, na.rm = TRUE)

      cat(paste0("Differential Expression (ProDA)",
                 if (!is.null(bait) && bait != "Bait") paste0(" \u2014 Bait: ", bait) else "",
                 "\n"))
      cat(paste0("Total proteins tested: ", nrow(df), "\n\n"))
      cat(paste0("  Up   (log2FC > 1, adj.p < 0.05):  ", n_up, "\n"))
      cat(paste0("  Down (log2FC < -1, adj.p < 0.05): ", n_down, "\n"))
      cat(paste0("  Not significant:                  ", nrow(df) - n_up - n_down - sum(le), "\n"))
      cat(paste0("  Low evidence (<2 values in a group): ", sum(le), "\n"))

    } else if (mode == "saint" && !is.null(de_rv$result_df)) {
      df   <- de_rv$result_df
      bait <- de_rv$bait_name
      sp   <- as.numeric(df$SaintScore)
      sp[is.na(sp)] <- 0
      n_saint_high <- sum(sp >= 0.80)
      n_saint_mid  <- sum(sp >= 0.50 & sp < 0.80)
      n_saint_low  <- sum(sp <  0.50)
      n_excl <- sum(df$Confidence == "Exclusive", na.rm = TRUE)
      n_high <- sum(df$Confidence == "High",      na.rm = TRUE)
      n_med  <- sum(df$Confidence == "Medium",    na.rm = TRUE)
      n_low  <- sum(df$Confidence == "Low",       na.rm = TRUE)
      n_bfdr <- sum(as.numeric(df$BFDR) < 0.05,  na.rm = TRUE)

      cat(paste0("Differential Expression (SAINTexpress",
                 if (!is.null(bait) && bait != "Bait") paste0(" \u2014 Bait: ", bait) else "",
                 ")\n"))
      cat(paste0("Total proteins scored: ", nrow(df), "\n\n"))
      cat("SAINT Score tiers:\n")
      cat(paste0("  High   (\u2265 0.80):      ", n_saint_high, "\n"))
      cat(paste0("  Medium (0.50\u20130.80): ", n_saint_mid,  "\n"))
      cat(paste0("  Low    (< 0.50):      ", n_saint_low,  "\n\n"))
      cat("Confidence (includes bait-only bonus):\n")
      cat(paste0("  Exclusive: ", n_excl, "\n"))
      cat(paste0("  High:      ", n_high, "\n"))
      cat(paste0("  Medium:    ", n_med,  "\n"))
      cat(paste0("  Low:       ", n_low,  "\n\n"))
      cat(paste0("BFDR < 0.05: ", n_bfdr, " proteins\n"))

    } else if (!is.null(de_rv$plot_data)) {
      df <- de_rv$plot_data
      n_up   <- sum(df$Log2FC >  1 & df$AdjPValue < 0.05, na.rm = TRUE)
      n_down <- sum(df$Log2FC < -1 & df$AdjPValue < 0.05, na.rm = TRUE)
      n_capped <- sum(df$IsCapped, na.rm = TRUE)

      mode_label <- if (identical(mode, "imputation")) "Imputation + limma"
                    else "Proteome Discoverer"
      cat(paste0("Differential Expression (", mode_label, ")\n"))
      cat(paste0("Total proteins: ", nrow(df), "\n"))
      cat(paste0("  Capped/sentinel: ", n_capped, "\n"))
      cat(paste0("Significant (Adj.P < 0.05):\n"))
      cat(paste0("  Up-regulated   (Log2FC > 1):  ", n_up, "\n"))
      cat(paste0("  Down-regulated (Log2FC < -1): ", n_down, "\n"))

    } else {
      cat("Upload a CSV and run DE analysis to see results.\n")
    }
  })

  # ===========================================================================
  # Download CSV (Step 7)
  # ===========================================================================
  output$de_download_csv <- downloadHandler(
    filename = function() {
      paste0("DE_Table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      mode <- de_rv$analysis_mode

      if (identical(mode, "dia") && !is.null(de_rv$result_df)) {
        out_df <- de_rv$result_df
        is_lod <- !is.null(out_df$perm_fdr)

        # ── Remove duplicate alias columns created for the shared plot/export code ────────
        # (Accession=Protein, GeneName=Gene, Log2FC=log2FC, AdjPvalue=adj.pvalue)
        drop_dupe <- c("Accession", "GeneName", "Log2FC", "AdjPvalue", "df_welch")
        out_df <- out_df[, !names(out_df) %in% drop_dupe, drop = FALSE]

        if (is_lod) {
          # Perseus LOD mode: adj.pvalue holds the permutation FDR (the value
          # that determines significance) and raw_pvalue holds the s0-modified
          # Welch p. Both are exported alongside perm_fdr so a consumer
          # thresholding adj.pvalue gets an FDR, not a raw p. `pvalue`
          # duplicates raw_pvalue exactly and is dropped.
          out_df <- out_df[, names(out_df) != "pvalue", drop = FALSE]
          names(out_df)[names(out_df) == "log2FC"]   <- "Log2FC"
          names(out_df)[names(out_df) == "t_stat"]   <- "t_stat_s0adj"
        } else {
          # Standard limma mode: rename clearly
          names(out_df)[names(out_df) == "log2FC"]     <- "Log2FC"
          names(out_df)[names(out_df) == "pvalue"]     <- "Pvalue"
          names(out_df)[names(out_df) == "adj.pvalue"] <- "AdjP_BH"
        }

        # The limma table and the NMAR (Fisher) table receive INDEPENDENT BH
        # corrections. Concatenating them and filtering at FDR 5% gives an
        # invalid combined list, so the family each row belongs to travels with
        # the export.
        out_df$FDR_family <- "limma_tested"

        # Insert blank separator before Raw_ abundance columns
        raw_cols <- grep("^Raw_", names(out_df), value = TRUE)
        if (length(raw_cols) > 0) {
          first_raw <- which(names(out_df) == raw_cols[1])
          out_df <- cbind(
            out_df[, seq_len(first_raw - 1), drop = FALSE],
            `---` = NA_character_,
            out_df[, first_raw:ncol(out_df), drop = FALSE]
          )
        }
        write.csv(out_df, file, row.names = FALSE)

      } else if (mode == "proda" && !is.null(de_rv$result_df)) {
        out_df <- de_rv$result_df
        # Merge raw abundances if available
        if (!is.null(de_rv$raw_abund_df)) {
          out_df <- merge(out_df, de_rv$raw_abund_df, by = "Accession", all.x = TRUE)
        }
        names(out_df)[names(out_df) == "AdjPvalue"]  <- "Adj.Pvalue"
        write.csv(out_df, file, row.names = FALSE)

      } else if (mode == "saint" && !is.null(de_rv$result_df)) {
        out_df <- de_rv$result_df
        # Strip internal columns (AdjPvalue = floor-corrected BFDR for plot computation;
        # Confidence = colour-scheme tier; the clean BFDR column is already exported)
        out_df <- out_df[, !names(out_df) %in% c("AdjPvalue", "Confidence"), drop = FALSE]
        write.csv(out_df, file, row.names = FALSE)

      } else if (!is.null(de_rv$plot_data)) {
        out_df <- de_rv$plot_data
        # Add regulation
        out_df$Regulation <- ifelse(
          out_df$AdjPValue < 0.05 & out_df$Log2FC >  1, "Up-regulated",
          ifelse(out_df$AdjPValue < 0.05 & out_df$Log2FC < -1, "Down-regulated",
                 "Not significant"))
        # Merge raw abundances
        if (!is.null(de_rv$raw_abund_df))
          out_df <- merge(out_df, de_rv$raw_abund_df, by = "Accession", all.x = TRUE)
        # Merge normalized abundances (mode c)
        if (!is.null(de_rv$norm_abund_df))
          out_df <- merge(out_df, de_rv$norm_abund_df, by = "Accession", all.x = TRUE)
        write.csv(out_df, file, row.names = FALSE)
      }

      logger("DE: CSV downloaded")
    }
  )

  # ===========================================================================
  # Download NMAR Analysis (DIA only)
  # Proteins excluded from limma by the detection filter and tested for
  # presence/absence with Fisher's exact test (dia_nmar_fisher()).
  # ===========================================================================
  # dia_nmar_result persists across a later TMT/LF run, so gate on the mode too.
  .nmar_export_df <- reactive({
    df <- de_rv$dia_nmar_result
    if (!identical(de_rv$analysis_mode, "dia") || is.null(df) || nrow(df) == 0L)
      return(NULL)
    df
  })

  output$de_nmar_download_ready <- reactive({
    !is.null(.nmar_export_df())
  })
  outputOptions(output, "de_nmar_download_ready", suspendWhenHidden = FALSE)

  output$de_download_nmar <- downloadHandler(
    filename = function() {
      .safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x %||% "")
      paste0("NMAR_Analysis_", .safe(de_rv$dia_nmar_test), "_vs_",
             .safe(de_rv$dia_nmar_ref), "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      out_df <- .nmar_export_df()
      if (is.null(out_df)) {
        utils::write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      out_df$Comparison <- paste0(de_rv$dia_nmar_test, " vs ", de_rv$dia_nmar_ref)
      # Fisher BH is applied within this set only, independent of the limma
      # table (see de_download_csv), so the family travels with the export.
      out_df$FDR_family <- "nmar_fisher"

      # Per-sample abundances the Fisher test counted detections from
      # (NA = not detected), test samples first, after a blank separator.
      mat   <- rv$dia_master_cache$protein_matrix
      grps  <- rv$dia_cond_groups
      samps <- if (!is.null(mat) && !is.null(grps))
        intersect(c(grps[[de_rv$dia_nmar_test]], grps[[de_rv$dia_nmar_ref]]),
                  colnames(mat))
      else character(0)
      if (length(samps) > 0L) {
        abund <- as.data.frame(
          mat[match(out_df$Protein, rownames(mat)), samps, drop = FALSE],
          check.names = FALSE)
        colnames(abund) <- paste0("Raw_", samps)
        out_df <- cbind(out_df, `---` = NA_character_, abund)
      }

      utils::write.csv(out_df, file, row.names = FALSE)
      logger("DE: NMAR analysis CSV downloaded")
    }
  )

  # ===========================================================================
  # GSEA Export  -" visibility flags
  # ===========================================================================
  output$de_gsea_rnk_ready <- reactive({
    de_rv$analysis_mode %in% c("proda", "saint", "welch", "dia", "tmt") && !is.null(de_rv$result_df)
  })
  outputOptions(output, "de_gsea_rnk_ready", suspendWhenHidden = FALSE)

  output$de_gsea_gct_ready <- reactive({
    has_fc_table <- de_rv$analysis_mode %in% c("proda", "saint", "welch", "dia", "tmt") &&
      !is.null(de_rv$result_df)
    has_imputation <- identical(de_rv$analysis_mode, "imputation") &&
      !is.null(de_rv$norm_abund_df)
    has_fc_table || has_imputation
  })
  outputOptions(output, "de_gsea_gct_ready", suspendWhenHidden = FALSE)

  # ===========================================================================
  # GSEA Export  -" .rnk (Ranked List)
  # Format: tab-delimited, no header, GeneName \t SignedScore
  # Metric: signed pi-score = sign(Log2FC) * |Log2FC| * -log10(adj. p), all modes
  # (adj. p = BH for ProDA/limma, BFDR for SAINT)
  # ===========================================================================
  output$de_download_rnk <- downloadHandler(
    filename = function() {
      bait <- if (is.null(de_rv$bait_name)) "analysis" else de_rv$bait_name
      paste0("GSEA_Ranked_", bait, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rnk")
    },
    content = function(file) {
      req(de_rv$result_df)
      df <- de_rv$result_df

      # Signed ranking metric: fold change x significance (pi-score), every mode.
      # PiScore is |Log2FC| * -log10(adj. p) where the engine provides it
      # (ProDA, SAINT, welch); otherwise it is computed the same way here.
      pi <- if (!is.null(df$PiScore)) df$PiScore else
            abs(df$log2FC %||% df$Log2FC) *
              (-log10(pmax(df$adj.pvalue %||% df$AdjPvalue %||% df$AdjPValue, 1e-300)))
      df$RankMetric <- sign(df$log2FC %||% df$Log2FC) * pi

      # Normalise gene-name column
      if (is.null(df$GeneName) && !is.null(df$Gene)) df$GeneName <- df$Gene

      # Keep only rows with valid gene names
      df <- df[!is.na(df$GeneName) & nchar(trimws(df$GeneName)) > 0, ]

      # Deduplicate gene names: keep the row with the highest absolute score
      df <- df[order(-abs(df$RankMetric)), ]
      df <- df[!duplicated(df$GeneName), ]

      # Sort descending by rank metric (GSEA convention)
      df <- df[order(-df$RankMetric), ]

      # Write tab-delimited, no header
      write.table(
        data.frame(df$GeneName, round(df$RankMetric, 6)),
        file      = file,
        sep       = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote     = FALSE
      )

      logger(paste0("DE GSEA: Exported .rnk with ", nrow(df), " genes"))
    }
  )

  # ===========================================================================
  # GSEA Export  -" .gct (Expression Matrix, GCT v1.2)
  # Header lines: #1.2 / nGenes \t nSamples
  # Columns: NAME \t Description \t sample1 \t sample2 ...
  # ===========================================================================
  output$de_download_gct <- downloadHandler(
    filename = function() {
      bait <- if (is.null(de_rv$bait_name)) "analysis" else de_rv$bait_name
      paste0("GSEA_Expression_", bait, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".gct")
    },
    content = function(file) {
      mode <- de_rv$analysis_mode

      if (mode %in% c("proda", "saint", "welch")) {
        req(de_rv$result_df)
        df <- de_rv$result_df
        if (identical(mode, "saint")) {
          # SAINT names replicate columns after the actual samples (log2 values)
          ctrl_rep_cols <- intersect(de_rv$ctrl_cols %||% character(0), colnames(df))
          cond_rep_cols <- intersect(de_rv$cond_cols %||% character(0), colnames(df))
        } else {
          ctrl_rep_cols <- sort(grep("^Ctrl_Rep", colnames(df), value = TRUE))
          cond_rep_cols <- sort(grep("^Cond_Rep", colnames(df), value = TRUE))
        }
        sample_cols <- c(ctrl_rep_cols, cond_rep_cols)

        if (length(sample_cols) == 0) {
          showNotification("No per-replicate abundance columns found in results.",
                           type = "error")
          return()
        }

        gct_df <- data.frame(
          NAME        = df$GeneName,
          Description = "na",
          stringsAsFactors = FALSE
        )
        for (sc in sample_cols) gct_df[[sc]] <- df[[sc]]

      } else if (identical(mode, "imputation")) {
        req(de_rv$norm_abund_df, de_rv$result_df)
        norm_df <- de_rv$norm_abund_df
        res_df  <- de_rv$result_df

        # result_df may not exist for imputation; use plot_data for gene names
        if (is.null(res_df)) {
          req(de_rv$plot_data)
          res_df <- de_rv$plot_data
        }

        # Match gene names by accession
        gene_lookup <- setNames(res_df$GeneName, res_df$Accession)
        norm_df$NAME <- gene_lookup[norm_df$Accession]
        norm_df$NAME[is.na(norm_df$NAME) | norm_df$NAME == ""] <- norm_df$Accession[
          is.na(norm_df$NAME) | norm_df$NAME == ""]

        sample_cols <- grep("^Norm_", colnames(norm_df), value = TRUE)
        if (length(sample_cols) == 0) {
          showNotification("No normalized abundance columns found.", type = "error")
          return()
        }

        # Clean up column names (remove Norm_ prefix)
        gct_df <- data.frame(
          NAME        = norm_df$NAME,
          Description = "na",
          stringsAsFactors = FALSE
        )
        for (sc in sample_cols) {
          clean_name <- sub("^Norm_", "", sc)
          gct_df[[clean_name]] <- norm_df[[sc]]
        }
        sample_cols <- sub("^Norm_", "", sample_cols)

      } else if (identical(mode, "dia")) {
        req(de_rv$result_df, de_rv$cond_test, de_rv$cond_ref)
        rdf   <- de_rv$result_df
        cache <- rv$dia_master_cache
        req(cache, cache$protein_matrix)
        mat  <- cache$protein_matrix
        grps <- rv$dia_cond_groups
        ref_samples  <- intersect(grps[[de_rv$cond_ref]],  colnames(mat))
        cond_samples <- intersect(grps[[de_rv$cond_test]], colnames(mat))
        sample_cols  <- c(ref_samples, cond_samples)
        if (length(sample_cols) == 0) {
          showNotification("No sample columns found in DIA protein matrix.", type = "error")
          return()
        }
        sub_mat    <- mat[rdf$Protein, sample_cols, drop = FALSE]
        gene_names <- rdf$GeneName %||% rdf$Gene
        gene_names[is.na(gene_names) | gene_names == ""] <-
          rdf$Protein[is.na(gene_names) | gene_names == ""]
        gct_df <- data.frame(NAME = gene_names, Description = "na",
                             stringsAsFactors = FALSE)
        for (sc in sample_cols) gct_df[[sc]] <- sub_mat[, sc]

      } else if (identical(mode, "tmt")) {
        req(de_rv$result_df, de_rv$cond_test, de_rv$cond_ref)
        lc  <- rv$tmt_limma_cache
        req(!is.null(lc), !is.null(lc$protein_matrix))
        rdf <- de_rv$result_df
        ann <- lc$annotation_map
        pm  <- lc$protein_matrix
        keep_conds  <- c(de_rv$cond_ref, de_rv$cond_test)
        use_chans   <- ann$ColumnName[ann$Condition %in% keep_conds]
        use_chans   <- intersect(use_chans, colnames(pm))
        if (length(use_chans) == 0) {
          showNotification("No TMT channels found for the selected conditions.", type = "error")
          return()
        }
        sub_mat     <- pm[, use_chans, drop = FALSE]
        sample_cols <- colnames(sub_mat)
        keep_prots  <- intersect(as.character(rdf$Protein), rownames(sub_mat))
        sub_mat     <- sub_mat[keep_prots, , drop = FALSE]
        gene_lookup <- setNames(as.character(rdf$GeneName %||% rdf$Gene),
                                as.character(rdf$Protein))
        gene_names  <- gene_lookup[rownames(sub_mat)]
        gene_names[is.na(gene_names) | !nzchar(gene_names)] <-
          rownames(sub_mat)[is.na(gene_names) | !nzchar(gene_names)]
        gct_df <- data.frame(NAME = gene_names, Description = "na",
                             stringsAsFactors = FALSE)
        for (sc in sample_cols) gct_df[[sc]] <- sub_mat[, sc]

      } else {
        showNotification("GCT export requires ProDA, DIA, or Imputation mode.",
                         type = "error")
        return()
      }

      # Pull description from protein_df if available
      if (!is.null(de_rv$protein_df) && "Description" %in% colnames(de_rv$protein_df)) {
        acc_to_desc <- setNames(
          as.character(de_rv$protein_df[["Description"]]),
          as.character(de_rv$protein_df[["Accession"]])
        )
        if (mode %in% c("proda", "saint", "welch") && !is.null(de_rv$result_df)) {
          desc_vals <- acc_to_desc[de_rv$result_df$Accession]
        } else if (!is.null(de_rv$norm_abund_df)) {
          desc_vals <- acc_to_desc[de_rv$norm_abund_df$Accession]
        } else {
          desc_vals <- rep(NA, nrow(gct_df))
        }
        desc_vals[is.na(desc_vals) | desc_vals == ""] <- "na"
        gct_df$Description <- desc_vals
      }

      # Remove rows with empty/NA gene names
      gct_df <- gct_df[!is.na(gct_df$NAME) & nchar(trimws(gct_df$NAME)) > 0, ]

      # Deduplicate gene names: keep the row with highest mean abundance
      abund_cols <- setdiff(colnames(gct_df), c("NAME", "Description"))
      row_means <- rowMeans(gct_df[, abund_cols, drop = FALSE], na.rm = TRUE)
      gct_df$.row_mean <- row_means
      gct_df <- gct_df[order(-gct_df$.row_mean), ]
      gct_df <- gct_df[!duplicated(gct_df$NAME), ]
      gct_df$.row_mean <- NULL

      n_genes   <- nrow(gct_df)
      n_samples <- length(abund_cols)

      # Write GCT v1.2 format
      con <- file(file, "w")
      writeLines("#1.2", con)
      writeLines(paste0(n_genes, "\t", n_samples), con)
      close(con)

      write.table(
        gct_df,
        file      = file,
        sep       = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote     = FALSE,
        append    = TRUE
      )

      logger(paste0("DE GSEA: Exported .gct with ", n_genes, " genes, ",
                    n_samples, " samples"))
    }
  )

  # ===========================================================================
  # GSEA Export  -" .cls (Phenotype Labels)
  # Format: 3 lines
  #   <nSamples> 2 1
  #   # <label1> <label2>
  #   0 0 0 1 1 1
  # ===========================================================================
  output$de_download_cls <- downloadHandler(
    filename = function() {
      bait <- if (is.null(de_rv$bait_name)) "analysis" else de_rv$bait_name
      paste0("GSEA_Phenotype_", bait, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".cls")
    },
    content = function(file) {
      mode <- de_rv$analysis_mode
      bait <- if (is.null(de_rv$bait_name) || de_rv$bait_name == "Bait")
                "Condition" else de_rv$bait_name
      # Sanitize label (no spaces)
      bait <- gsub("\\s+", "_", bait)

      if (mode %in% c("proda", "saint", "welch")) {
        req(de_rv$result_df)
        if (identical(mode, "saint")) {
          # same sample columns, in the same order, as the SAINT .gct above
          ctrl_rep_cols <- intersect(de_rv$ctrl_cols %||% character(0), colnames(de_rv$result_df))
          cond_rep_cols <- intersect(de_rv$cond_cols %||% character(0), colnames(de_rv$result_df))
        } else {
          ctrl_rep_cols <- sort(grep("^Ctrl_Rep", colnames(de_rv$result_df), value = TRUE))
          cond_rep_cols <- sort(grep("^Cond_Rep", colnames(de_rv$result_df), value = TRUE))
        }
        n_ctrl <- length(ctrl_rep_cols)
        n_cond <- length(cond_rep_cols)
        if (n_ctrl == 0 || n_cond == 0) {
          showNotification("Could not determine sample counts for CLS export.", type = "error")
          return()
        }
      } else if (identical(mode, "imputation")) {
        n_ctrl <- length(de_rv$ctrl_cols)
        n_cond <- length(de_rv$cond_cols)
      } else if (mode %in% c("dia", "tmt")) {
        # Use stored comparison groups to get replicate counts
        if (identical(mode, "dia")) {
          grps   <- rv$dia_cond_groups
          n_ctrl <- length(grps[[de_rv$cond_ref]])
          n_cond <- length(grps[[de_rv$cond_test]])
        } else {
          # TMT: count channels per condition from annotation map
          ann_cls <- rv$tmt_annotation_map %||% rv$tmt_limma_cache$annotation_map
          if (!is.null(ann_cls)) {
            bio_ann <- ann_cls[!ann_cls$Condition %in% "Norm" &
                               !is.na(ann_cls$Condition), , drop = FALSE]
            n_ctrl <- sum(bio_ann$Condition == de_rv$cond_ref, na.rm = TRUE)
            n_cond <- sum(bio_ann$Condition == de_rv$cond_test, na.rm = TRUE)
          } else {
            n_ctrl <- 0L; n_cond <- 0L
          }
        }
        if (n_ctrl == 0 || n_cond == 0) {
          showNotification("Could not determine sample counts for CLS export.", type = "error")
          return()
        }
      } else {
        showNotification("CLS export requires ProDA or Imputation mode.",
                         type = "error")
        return()
      }

      n_samples <- n_ctrl + n_cond
      labels <- c(rep(0, n_ctrl), rep(1, n_cond))

      con <- file(file, "w")
      writeLines(paste(n_samples, "2", "1"), con)
      writeLines(paste("#", "Control", bait), con)
      writeLines(paste(labels, collapse = " "), con)
      close(con)

      logger(paste0("DE GSEA: Exported .cls (", n_ctrl, " ctrl, ", n_cond, " cond)"))
    }
  )

}
