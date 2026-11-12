# server_ptm_de.R
# ==============================================================================
# PTM Differential Expression (TMT) module
# Primary: limma-based stoichiometric PTM DE (DIA and TMT PTM)
# No fallback to MSstatsTMT in this module.
# ==============================================================================

server_ptm_de <- function(input, output, session, rv, logger) {

  ptm_rv <- reactiveValues(
    parsed_files    = list(),
    files           = list(),
    n_files         = 0,
    # Design grid state (required by tmt_register_design_grid)
    all_channels    = character(0),
    n_plexes        = 1L,
    has_bridge      = FALSE,
    design_conds    = character(0),
    design_option   = NULL,         # 1, 2, or 3
    annotation_map  = NULL,         # data.frame from design grid
    # PTM type selection
    detected_ptms   = character(0),
    selected_ptms   = "Phospho",
    # Data tables — stored in rv$ptm_data (shared store); gene_lookup kept local
    gene_lookup_df  = data.frame(Accession = character(0), Gene = character(0),
                                stringsAsFactors = FALSE),
    protein_headers = character(0),
    protein_df      = NULL,
    conditions      = NULL,
    de_result       = NULL,
    de_result_adjusted   = NULL,    # ADJUSTED.Model (stoichiometry)
    de_result_unadjusted = NULL,    # PTM.Model (absolute PTM)
    de_result_protein    = NULL,    # PROTEIN.Model (protein-level)
    de_master       = NULL,         # merged master table (all 3 models)
    de_view_mode    = "adjusted",   # "adjusted" | "unadjusted" | "protein"
    method_used     = NULL,
    comparison      = NULL,
    highlight_gene  = NULL,
    highlight_ptm   = NULL,
    highlight_rules    = list(),
    highlight_counter  = 0L,
    filter_conditions  = list(),
    filter_counter     = 0L,
    render_count     = 0L,
    de_result_prepped = NULL,  # de_result with neg_log10_p and hover pre-computed
    # Hybrid limma fork
    pending_hybrid_comparison = NULL,  # stored comparison while modal is open
    bypass_hybrid_check       = FALSE, # TRUE = skip fork modal
    force_standard_pipeline   = FALSE, # TRUE = user chose standard over hybrid
    status          = "Use the PTM upload wizard to load data.",
    data_format     = NULL,
    fc_scatter_data = NULL       # stores computed FC x FC data.frame + meta
  )
  methodsTextServer("ptm_de_methods", provenance = rv$provenance,
                    section = reactive({
                      if (isTRUE(ptm_rv$design_option %in% c("dia_ptm", "tmt_ptm")))
                        "dia_ptm_de"
                      else
                        "ptm_de"
                    }))

  # Offer only the views this cache has (ptm_layers(): a PTM-only cache has the
  # absolute PTM view alone) and select its primary one.
  .sync_view_choices <- function(cache) {
    lay <- ptm_layers(cache)
    ptm_rv$de_view_mode <- lay[[1L]]
    updateRadioButtons(session, "ptm_view_mode",
                       choices = setNames(lay, PTM_LAYER_LABELS[lay]), selected = lay[[1L]])
  }

  output$ptm_view_mode_help <- renderUI({
    lay  <- ptm_layers(rv$dia_ptm_master_cache)
    name <- c(adjusted = "Adjusted:", unadjusted = "Absolute PTM:", protein = "Protein:")
    help <- c(adjusted   = " PTM change minus protein change (stoichiometry).",
              unadjusted = " Raw PTM fold change (not protein-corrected).",
              protein    = " Total protein fold change only.")
    helpText(class = "small text-muted mb-2",
      if (!ptm_has_global(rv$dia_ptm_master_cache))
        tagList(tags$b("PTM only:"), " no global proteome was uploaded, so a site's fold change ",
                "includes any change in its protein's abundance.", tags$br()),
      lapply(seq_along(lay), function(i)
        tagList(tags$b(name[[lay[i]]]), help[[lay[i]]], if (i < length(lay)) tags$br())))
  })

  # Sync from other PTM modules when source changes
  observeEvent(rv$ptm_data$source, {
    if (!isTRUE(rv$ptm_data$source %in% c("ptm_trends", "ptm_single", "dia_ptm_upload", "tmt_ptm_upload"))) return()
    has_tables <- !is.null(rv$ptm_data$ptm_long) && nrow(rv$ptm_data$ptm_long) > 0
    is_dia_ptm_src <- isTRUE(rv$ptm_data$source %in% c("dia_ptm_upload", "tmt_ptm_upload"))
    if (!has_tables && !is_dia_ptm_src) return()

    # Guard against annotation_map observer wiping cache during sync
    ptm_rv$importing <- TRUE
    session$onFlushed(function() { ptm_rv$importing <- FALSE }, once = TRUE)

    # Tables are read directly from rv$ptm_data — no local copies needed
    # Always sync gene_lookup_df (even without ptm_long, e.g. cache-only import)
    if (!is.null(rv$ptm_data$gene_lookup))
      ptm_rv$gene_lookup_df  <- rv$ptm_data$gene_lookup
    if (!is.null(rv$ptm_data$data_format))
      ptm_rv$data_format    <- rv$ptm_data$data_format
    if (!is.null(rv$ptm_data$protein_headers))
      ptm_rv$protein_headers <- rv$ptm_data$protein_headers
    if (!is.null(rv$ptm_data$protein_df))
      ptm_rv$protein_df      <- rv$ptm_data$protein_df
    if (!is.null(rv$ptm_data$n_files))
      ptm_rv$n_files         <- rv$ptm_data$n_files
    if (!is.null(rv$ptm_data$annotation_map))
      ptm_rv$annotation_map  <- rv$ptm_data$annotation_map
    if (!is.null(rv$ptm_data$design_option))
      ptm_rv$design_option   <- rv$ptm_data$design_option
    ptm_rv$conditions      <- sort(unique(as.character(rv$ptm_data$conditions)))
    ptm_rv$de_result       <- NULL

    # Populate Group A / Group B dropdowns
    conds <- ptm_rv$conditions
    if (length(conds) >= 2) {
      updateSelectInput(session, "ptm_group_a", choices = conds, selected = conds[1])
      updateSelectInput(session, "ptm_group_b", choices = conds, selected = conds[2])
    } else if (length(conds) == 1) {
      updateSelectInput(session, "ptm_group_a", choices = conds)
      updateSelectInput(session, "ptm_group_b", choices = conds)
    }

    if (has_tables) {
      n_ptm  <- length(unique(rv$ptm_data$ptm_long$ProteinName))
      n_prot <- length(unique(rv$ptm_data$protein_long$ProteinName))
      src    <- rv$ptm_data$source
      ptm_rv$status <- paste0("Synced from ", src, ": ", n_ptm, " PTM sites, ", n_prot, " proteins.")
      logger(paste0("PTM DE: Synced from ", src, " — ", n_ptm, " PTM sites"))
    } else if (is_dia_ptm_src) {
      cache   <- rv$dia_ptm_master_cache
      n_sites <- if (!is.null(cache)) cache$n_sites_output else "?"
      src_tag <- rv$ptm_data$source
      .sw_sync_lbl <- if (identical(src_tag, "tmt_ptm_upload")) {
        "TMT PTM (limma)"
      } else {
        switch(rv$dia_ptm_software %||% "spectronaut",
          "maxquant"          = "Label-Free DDA (MaxQuant)",
          "spectronaut_pivot" = "DIA (Spectronaut Pivot)",
          "DIA (Spectronaut)")
      }
      ptm_rv$status <- paste0(.sw_sync_lbl, " PTM data synced: ", n_sites,
                              " phosphosites across ",
                              length(ptm_rv$conditions), " conditions",
                              if (!ptm_has_global(cache)) " (PTM only: absolute PTM view)" else "",
                              ". Select groups and click Run DE.")
      logger(paste0("PTM DE: ", .sw_sync_lbl, " PTM synced — ", n_sites, " sites"))
      .sync_view_choices(cache)
      .enable_run_btn()
    } else {
      src <- rv$ptm_data$source
      ptm_rv$status <- paste0("Design synced from ", src, " (summarization cache imported).")
      logger(paste0("PTM DE: Design synced from ", src, " (import mode, no intermediate tables)"))
    }
  }, ignoreInit = TRUE)


  # Helper: enable/disable Run button with colour swap
  # ---------------------------------------------------------------------------
  .enable_run_btn <- function() {
    shinyjs::enable("ptm_run_btn")
    shinyjs::removeClass("ptm_run_btn", "btn-secondary")
    shinyjs::addClass("ptm_run_btn", "btn-success")
  }

  .disable_run_btn <- function() {
    shinyjs::disable("ptm_run_btn")
    shinyjs::removeClass("ptm_run_btn", "btn-success")
    shinyjs::addClass("ptm_run_btn", "btn-secondary")
  }

  # Also enable when DIA PTM cache is populated
  observeEvent(rv$dia_ptm_master_cache, {
    if (!is.null(rv$dia_ptm_master_cache)) {
      .enable_run_btn()
    }
  }, ignoreInit = TRUE)

  # ---------------------------------------------------------------------------
  # Run analysis (group comparison only — requires summarization first)
  # ---------------------------------------------------------------------------
  observeEvent(input$ptm_run_btn, {

    is_dia_ptm <- isTRUE(ptm_rv$design_option %in% c("dia_ptm", "tmt_ptm"))
    if (!is_dia_ptm) return()

    group_a <- input$ptm_group_a
    group_b <- input$ptm_group_b
    if (is.null(group_a) || is.null(group_b) || group_a == "" || group_b == "") {
      showNotification("Please select both Group A and Group B.", type = "warning")
      return()
    }
    if (group_a == group_b) {
      showNotification("Group A and Group B must be different.", type = "warning")
      return()
    }
    cache <- rv$dia_ptm_master_cache
    if (is.null(cache)) {
      showNotification("DIA PTM data not found. Please complete the DIA PTM Upload Wizard.", type = "warning")
      return()
    }

    comparison <- paste0(group_a, " vs ", group_b)
    ptm_rv$comparison           <- comparison
    ptm_rv$highlight_gene       <- NULL
    ptm_rv$highlight_ptm        <- NULL
    ptm_rv$de_result_adjusted   <- NULL
    ptm_rv$de_result_unadjusted <- NULL
    ptm_rv$de_result_protein    <- NULL
    ptm_rv$de_result            <- NULL
    ptm_rv$fc_scatter_data      <- NULL
    .sw_run_lbl <- if (isTRUE(ptm_rv$design_option == "tmt_ptm")) {
      "TMT PTM"
    } else {
      switch(rv$dia_ptm_software %||% "spectronaut",
        "maxquant"          = "Label-Free DDA",
        "spectronaut_pivot" = "DIA (Spectronaut Pivot)",
        "DIA")
    }
    logger(paste0("PTM DE [", .sw_run_lbl, " PTM]: Running comparison ", comparison))

    withProgress(message = paste0("Running ", .sw_run_lbl, " PTM comparison..."), value = 0.1, {
      # For TMT PTM with >=2 plexes and non-IRS normalization (design_opt != 2),
      # include plex as a batch covariate so inter-plex loading differences do
      # not inflate within-condition variance and deflate t-statistics.
      # IRS (design_opt=2) already equalises plexes; batch term is not needed.
      .batch_grps <- ptm_de_batch_groups(cache, ptm_rv$design_option)
      .batch_for_model <- !is.null(.batch_grps)
      if (.batch_for_model)
        logger(paste0("PTM DE [TMT PTM]: plex batch covariate included (",
                      length(.batch_grps), " plexes)."))
      res <- tryCatch(
        dia_ptm_limma_de(
          cache           = cache,
          cond_test       = group_a,
          cond_ref        = group_b,
          de_min_valid_per_cond = as.integer(
            cache$min_valid_per_cond %||% formals(dia_ptm_limma_de)$de_min_valid_per_cond),
          batch_groups    = .batch_grps,
          batch_for_model = .batch_for_model,
          subject_pairing = rv$dia_ptm_subject_pairing,
          logger          = logger
        ),
        error = function(e) list(error = conditionMessage(e))
      )
      incProgress(0.8)
    })

    if (!is.null(res$error)) {
      showNotification(paste0("DIA PTM DE failed: ", res$error), type = "error", duration = 10)
      ptm_rv$status <- paste0("Error: ", res$error)
      return()
    }

    # Memoise the cache-scoped k-NN matrices so every later contrast reuses the
    # same donor pool and the same filled values (computed lazily on first run).
    # These three fields are session-scratch: the export handler in
    # server_dia_ptm_upload.R strips them before writing the prerun bundle, so
    # they never reach disk. ptm_impute_cap records the cap they were built
    # under, and dia_ptm_limma_de() invalidates the memo when it differs.
    if (!is.null(res$imputed_cache) && is.null(cache$ptm_abs_matrix_imputed)) {
      rv$dia_ptm_master_cache$ptm_abs_matrix_imputed <-
        res$imputed_cache$ptm_abs_matrix_imputed
      rv$dia_ptm_master_cache$global_matrix_imputed  <-
        res$imputed_cache$global_matrix_imputed
      rv$dia_ptm_master_cache$impute_stats           <-
        res$imputed_cache$impute_stats
      rv$dia_ptm_master_cache$ptm_impute_cap         <-
        res$imputed_cache$ptm_impute_cap
    }

    # ── Annotate all three DIA PTM results with display columns ──
    gene_map <- cache$site_to_gene    # abs_site_id -> gene label(s)
    acc_gene <- cache$acc_to_gene     # accession   -> gene name
    adj_map  <- cache$adj_to_abs_map  # adj_row_name -> abs_site_id

    # The three annotation closures that used to live here are now in global.R
    # (ptm_de_annotate_views), because PTM Single needs them too: a contrast the
    # DE tab has not run has no numbers in Single, and Single could not fit it
    # itself while the helpers were locals of this observer. The bodies moved
    # unchanged; gene_map / acc_gene / adj_map are read from the cache there.
    .views <- ptm_de_annotate_views(res, cache)
    ptm_rv$de_result_adjusted   <- .views$adjusted
    ptm_rv$de_result_unadjusted <- .views$unadjusted
    ptm_rv$de_result_protein    <- .views$protein
    # Default view = the cache's primary layer (adjusted; absolute PTM in PTM-only mode)
    .primary                    <- ptm_primary_layer(cache)
    ptm_rv$de_result            <- .views[[.primary]]
    .sync_view_choices(cache)
    ptm_rv$nmar_table           <- res$nmar_table
    # Recorded beside the table, from the same result object, so the saved bundle
    # cannot pair one run's table with another run's power flags.
    ptm_rv$nmar_has_power       <- isTRUE(res$nmar_has_power)
    ptm_rv$nmar_min_achievable_p <- as.numeric(res$nmar_min_achievable_p %||% NA_real_)

    # ── Publish to the shared store so single-site views report these numbers ──
    # PTM Single must never refit a model; it reads the result for the contrast it
    # is displaying, so a site's statistics are identical in both tabs.  Keyed by
    # comparison label ("<test> vs <ref>", the same string dia_ptm_limma_de() puts
    # in the Label column) so results accumulate across contrasts instead of the
    # most recent run erasing the previous one.
    #
    # The ANNOTATED tables are published, not the raw res$ tables: .annotate_*_view()
    # moves the matrix row key to PTM_ID and overwrites Protein with the parsed
    # accession, so PTM_ID is the join key for all three - row key for the two site
    # tables, accession for the protein table.
    if (is.null(rv$dia_ptm_de_results)) rv$dia_ptm_de_results <- list()
    .store <- rv$dia_ptm_de_results
    # Built by the shared constructor in global.R so a contrast fitted here and the
    # same contrast fitted on demand from PTM Single cannot differ in shape.
    .store[[as.character(comparison)]] <- ptm_de_store_entry(
      views = list(adjusted   = ptm_rv$de_result_adjusted,
                   unadjusted = ptm_rv$de_result_unadjusted,
                   protein    = ptm_rv$de_result_protein),
      res       = res,
      cache     = cache,
      cond_test = group_a,
      cond_ref  = group_b)
    rv$dia_ptm_de_results <- .store
    logger(paste0("PTM DE: published '", comparison, "' to the shared store (",
                  length(.store), " contrast(s) now available to single-site views)."))
    .sw_done_lbl <- switch(rv$dia_ptm_software %||% "spectronaut",
                      "maxquant"          = "Label-Free DDA",
                      "spectronaut_pivot" = "DIA (Spectronaut Pivot)",
                      "DIA")
    ptm_rv$status <- paste0(.sw_done_lbl, " PTM DE complete",
                             if (.primary == "unadjusted") " (PTM only, absolute PTM)" else "", ": ",
                             res[[paste0("n_tested_", .primary)]], " sites tested, ",
                             res[[paste0("n_estimable_", .primary)]] %||% NA, " estimable (",
                             comparison, ")")
    logger(paste0("PTM DE [", .sw_done_lbl, " PTM]: Done — ",
                  paste(vapply(ptm_layers(cache), function(l)
                    paste0(res[[paste0("n_tested_", l)]], " ",
                           c(adjusted = "adjusted", unadjusted = "absolute", protein = "protein")[[l]],
                           " (", res[[paste0("n_estimable_", l)]] %||% NA, " estimable)"),
                    character(1)), collapse = ", "), "."))

    # Build master table (for DE table tab + CSV download)
    ptm_rv$de_master <- .build_master_table(
      ptm_rv$de_result_adjusted,
      ptm_rv$de_result_unadjusted,
      ptm_rv$de_result_protein,
      adj_map,
      primary = .primary
    )
    if (!is.null(ptm_rv$de_master))
      logger(paste0("PTM DE [", .sw_done_lbl, " PTM]: Master table ",
                    nrow(ptm_rv$de_master), " rows."))

    # ── Provenance write (DE-specific) ──
    is_tmt_ptm_de <- isTRUE(ptm_rv$design_option == "tmt_ptm")
    .ptm_only_de <- !ptm_has_global(cache)
    rv$provenance$ptm_de <- list(
      ptm_mode          = cache$ptm_mode %||% "with_global",
      analysis_mode     = if (is_tmt_ptm_de) "tmt_ptm_limma"
                          else switch(rv$dia_ptm_software %||% "spectronaut",
                            "maxquant"          = "lf_dda_limma",
                            "spectronaut_pivot" = "dia_pivot_limma",
                            "dia_limma"),
      data_route        = if (is_tmt_ptm_de) "TMT PTM (Proteome Discoverer, limma eBayes)"
                          else switch(rv$dia_ptm_software %||% "spectronaut",
                            "maxquant"          = "Label-Free DDA (MaxQuant)",
                            "spectronaut_pivot" = "DIA (Spectronaut Pivot)",
                            "DIA (Spectronaut)"),
      comparison        = as.character(comparison),
      # Counts of the primary table: adjusted, or absolute PTM in PTM-only mode
      n_ptm_sites             = as.numeric(res[[paste0("n_tested_", .primary)]] %||% 0),
      n_ptm_sites_estimable   = as.numeric(res[[paste0("n_estimable_", .primary)]] %||% NA_real_),
      n_ptm_sites_unadjusted_estimable = as.numeric(res$n_estimable_unadjusted %||% NA_real_),
      n_proteins_estimable    = as.numeric(res$n_estimable_protein %||% NA_real_),
      estimability_note       = "Rows are 'tested' if they pass the per-contrast missingness filter; a row is 'estimable' only if the fit returned a non-NA coefficient. Non-estimable rows appear in the CSV with blank statistics and are absent from the volcano.",
      de_min_valid_per_cond         = as.numeric(res$de_min_valid_per_cond %||% NA_real_),
      de_min_valid_per_cond_effective = res$de_min_valid_per_cond_effective,
      min_valid_rule                = if (isTRUE(res$de_require_both_arms))
        "A site is retained for testing when BOTH compared conditions carry at least de_min_valid_per_cond valid values. The threshold is a flat count and is not scaled by replicate number. Sites meeting it in only one condition are excluded from the fitted model and analysed instead by the presence/absence (Fisher's exact) route."
      else
        "A site is retained for testing when at least one condition carries a valid value in at least 50% of that condition's samples, subject to a floor of de_min_valid_per_cond. The threshold is therefore computed per condition from that condition's own replicate count: effective = max(floor, ceiling(0.5 * n_samples_in_condition)). This is the same >=50%-valid-in-at-least-one-group rule the global DIA route applies.",
      knn_max_missing_per_cond      = as.numeric(res$knn_max_missing_per_cond %||% NA_real_),
      n_knn_rows_skipped_unadjusted = as.numeric(res$n_knn_rows_skipped_unadjusted %||% NA_real_),
      n_knn_rows_skipped_protein    = as.numeric(res$n_knn_rows_skipped_protein %||% NA_real_),
      n_knn_values_unadjusted       = as.numeric(res$n_knn_values_unadjusted %||% NA_real_),
      n_knn_values_protein          = as.numeric(res$n_knn_values_protein %||% NA_real_),
      imputation_eligibility_note   = if (is_tmt_ptm_de)
        "No imputation is performed on any matrix on this route: dia_ptm_limma_de() disables k-NN for pd_tmt caches. All three models are fitted by lmFit on the observed values only."
      else paste0(
        "kNN (k = 10, rowmax = 1, colmax = 0.8, maxp = Inf, seed 42) is applied only to rows with at most knn_max_missing_per_cond missing values in every condition in scope. The cap is a count of missing values, not a fraction. Rows above that cap keep their NAs and are fitted by lmFit on the observations that exist.",
        if (.ptm_only_de) "" else " The adjusted matrix is never imputed."),
      knn_cap_degenerate            = isTRUE(res$knn_cap_degenerate),
      knn_cap_exhaustive            = isTRUE(res$knn_cap_exhaustive),
      n_nmar                        = as.numeric(res$n_nmar %||% 0),
      nmar_has_power                = res$nmar_has_power,
      nmar_min_achievable_p         = as.numeric(res$nmar_min_achievable_p %||% NA_real_),
      nmar_note                     = "Sites failing the valid-value filter are not testable by limma and were instead assessed for a non-random detection pattern by Fisher's exact test on detected/not-detected counts, BH-corrected within this excluded set only. This correction is independent of the correction applied to the limma-tested set, so no joint false discovery rate is claimed across the two. Where the replicate count makes the smallest attainable two-sided Fisher p-value exceed 0.05, p-values are reported as NA and the set is presented as a ranked descriptive list of detection frequencies.",
      imputation_scope              = as.character(res$knn_scope %||% "unknown"),
      imputation_scope_note         = "cache = imputation computed once over every sample in the cache, so a site's filled values and p-value do not depend on which contrast is running; eligibility is evaluated in every condition of the cache. contrast = legacy fallback, imputed within the two contrasted groups only. The DE missingness filter stays contrast-scoped either way.",
      n_knn_rows_skipped_cache_ptm     = as.numeric(res$imputed_cache$impute_stats$n_rows_skipped_ptm %||% NA_real_),
      n_knn_rows_skipped_cache_protein = as.numeric(res$imputed_cache$impute_stats$n_rows_skipped_protein %||% NA_real_),
      n_knn_values_cache_ptm           = as.numeric(res$imputed_cache$impute_stats$n_values_filled_ptm %||% NA_real_),
      n_knn_values_cache_protein       = as.numeric(res$imputed_cache$impute_stats$n_values_filled_protein %||% NA_real_),
      n_significant           = as.numeric(sum(
                                  .views[[.primary]]$adj.pvalue < 0.05, na.rm = TRUE)),
      n_ptm_sites_unadjusted  = as.numeric(res$n_tested_unadjusted %||% 0),
      n_significant_unadjusted = if (!is.null(ptm_rv$de_result_unadjusted))
                                   as.numeric(sum(
                                     ptm_rv$de_result_unadjusted$adj.pvalue < 0.05,
                                     na.rm = TRUE))
                                 else numeric(),
      n_proteins              = as.numeric(res$n_tested_protein %||% 0),
      prob_threshold          = as.numeric(cache$prob_threshold %||% 75),
      design_option           = as.numeric(cache$design_opt %||% 1L),
      software_used           = if (.ptm_only_de)
        "limma (eBayes with robust = TRUE and an intensity-variance trend (trend = TRUE)) on the absolute PTM matrix"
      else
        "limma (eBayes with robust = TRUE throughout; an intensity-variance trend (trend = TRUE) for the absolute PTM and global protein fits, and trend = FALSE for the stoichiometry-adjusted fit, whose row means are log-ratios rather than abundances)",
      pvalue_adjustment       = "BH",
      result_tables           = if (.ptm_only_de)
        "One limma fit: Absolute PTM. No global proteome was supplied (PTM-only mode), so there is no stoichiometry-adjusted or protein table, and a site's fold change includes any change in its parent protein's abundance."
      else if (is_tmt_ptm_de)
        "Three limma fits on the same samples: Adjusted (per-sample PTM log2 minus protein log2), Absolute PTM, and Protein. All three are built from the same measurements and are not independent evidence. No matrix is imputed on this route. Adjusted still does not equal Absolute minus Protein numerically, because the adjusted matrix is column-median-centred again after the subtraction and because each fit applies its own per-comparison missingness filter."
      else
        "Three limma fits on the same samples: Adjusted (per-sample PTM log2 minus protein log2), Absolute PTM, and Protein. All three are built from the same measurements and are not independent evidence. They are not arithmetically related either: the Absolute PTM and Protein matrices are kNN-imputed before fitting while the adjusted matrix is fitted on observed values only, so Adjusted does not equal Absolute minus Protein numerically.",
      fold_change_sources     = if (.ptm_only_de)
        "DE table log2FC comes from the absolute PTM limma fit."
      else
        "DE table log2FC comes from the limma fit for the selected view. The FC x FC scatter log2FC is computed separately from Perseus-style LOD-imputed values and will not match.",
      bh_family               = if (.ptm_only_de)
        "BH applied within the absolute PTM table. No correction across comparisons."
      else
        "BH applied independently within each of the three result tables. No correction across tables or across comparisons.",
      n_distinct_measured_sites = as.numeric({
        .ids <- .views[[.primary]]$PTM_ID
        .abs <- if (!is.null(adj_map)) unname(adj_map[.ids]) else rep(NA_character_, length(.ids))
        .abs[is.na(.abs)] <- .ids[is.na(.abs)]
        length(unique(.abs))
      }),
      shared_expansion_note   = if (.ptm_only_de)
        "Shared-peptide sites are NOT expanded in PTM-only mode (expansion is part of the stoichiometric adjustment). A site whose peptide maps to several proteins is tested once, as one row keyed by all of them, and flagged Is_Shared = TRUE."
      else if (is_tmt_ptm_de)
        "Shared-peptide sites are NOT expanded on this route. A peptide matching several master proteins is assigned to one of them by the razor rule (greatest total PSM support, summed across every phospho file), so each measured site contributes exactly one tested row and the tested row count equals the number of distinct measured sites. The remaining master proteins are retained in Also_Maps_To and the site is flagged Is_Shared = TRUE. There is no single agreed field standard for shared-peptide handling; the alternatives (assign to every matching protein, or merge them into one combined site) would each produce a different tested row count from the same measurements."
      else
        "Shared-peptide sites are expanded into one adjusted row per constituent protein (Is_Shared = TRUE). Expanded rows share a single measured intensity vector, so BH was applied over a row set larger than the number of distinct measured sites. BH remains valid under positive dependence, but the test count is not a count of independent hypotheses.",
      plex_covariate_used     = isTRUE(.batch_for_model),
      n_plex_batches          = if (isTRUE(.batch_for_model))
                                  as.integer(length(.batch_grps)) else NA_integer_,
      n_absolute_key_unmapped = as.numeric(
                                  attr(ptm_rv$de_master, "n_absolute_key_unmapped") %||% NA_real_),
      # No FC x FC scatter in PTM-only mode: it plots site against protein fold change
      fc_scatter_thresholds   = if (!.ptm_only_de) c(ptm = 1.0, protein = 0.5),
      fc_scatter_basis        = if (!.ptm_only_de) "FC x FC scatter categories are fold-change only, computed on Perseus-style LOD-imputed values (PTM |log2FC| >= 1.0, protein |log2FC| >= 0.5). No p-value or significance test enters the categorisation."
    )

    ptm_rv$method_used <- if (is_tmt_ptm_de) {
      "TMT PTM (Proteome Discoverer) \u2014 limma eBayes"
    } else {
      switch(rv$dia_ptm_software %||% "spectronaut",
        "maxquant"          = "Label-Free DDA (MaxQuant) \u2014 limma eBayes",
        "spectronaut_pivot" = "DIA (Spectronaut Pivot) \u2014 limma eBayes",
        "DIA (Spectronaut) \u2014 limma eBayes")
    }

    # PTM-only mode has no protein fold changes to plot sites against: no FC x FC scatter
    if (.ptm_only_de) {
      .enable_run_btn()
      return()
    }

    # -- STEP 2: FC x FC scatter (Perseus LOD imputation, same groups/cache) --
    LOD_DOWNSHIFT <- 1.8
    LOD_WIDTH     <- 0.3
    withProgress(message = "Computing FC x FC (Perseus LOD imputation)...", value = 0.1, {
      fc_df <- tryCatch({
        design_df <- cache$design_df
        samp_test <- design_df$sample_name[design_df$condition == group_a]
        samp_ref  <- design_df$sample_name[design_df$condition == group_b]
        ptm_fc_df <- lod_impute_fc_scatter(cache$ptm_abs_matrix, samp_test, samp_ref,
                                           LOD_DOWNSHIFT, LOD_WIDTH)
        glb_fc_df <- lod_impute_fc_scatter(cache$global_matrix, samp_test, samp_ref,
                                           LOD_DOWNSHIFT, LOD_WIDTH)
        if (is.null(ptm_fc_df) || is.null(glb_fc_df))
          stop("Insufficient data for selected comparison.")
        ptm_fc_df$acc <- sub("_.*$", "", ptm_fc_df$feature)
        # For a shared absolute site ID (ACC1_S44_Phospho_M1;ACC2_...) the substring
        # before the first underscore is the FIRST accession only, so the plotted
        # global fold change is that protein's regardless of which protein the site
        # belongs to.  Flagged here; which accession is used is left unchanged.
        ptm_fc_df$acc_ambiguous <- grepl(";", ptm_fc_df$feature, fixed = TRUE)
        colnames(glb_fc_df)[colnames(glb_fc_df) == "log2fc"] <- "global_log2fc"
        df_out <- merge(ptm_fc_df, glb_fc_df[, c("feature", "global_log2fc"), drop = FALSE],
                        by.x = "acc", by.y = "feature", all.x = TRUE)
        .n_before_glb  <- nrow(df_out)
        df_out <- df_out[!is.na(df_out$global_log2fc), , drop = FALSE]
        .n_dropped_glb <- .n_before_glb - nrow(df_out)
        logger(paste0("PTM DE [", .sw_run_lbl, " PTM]: FC x FC scatter dropped ",
                      .n_dropped_glb, " of ", .n_before_glb,
                      " site(s) with no matched protein fold change; ",
                      sum(df_out$acc_ambiguous),
                      " retained site(s) have an ambiguous parent accession."))
        attr(df_out, "n_dropped_no_global") <- .n_dropped_glb
        # Label: Gene_SitePositions (e.g. KRT16_S44+Y46)
        site_gene <- cache$site_to_gene
        df_out$PTM_Label <- vapply(df_out$feature, function(sid) {
          g <- if (!is.null(site_gene)) site_gene[sid] else NA_character_
          if (!is.na(g) && nzchar(g)) {
            gene       <- strsplit(g, ";", fixed = TRUE)[[1]][1]
            first_part <- strsplit(sid, ";", fixed = TRUE)[[1]][1]
            # Spectronaut format: ACC_S40_Phospho[...] captures "S40"
            sites_part <- sub("^[^_]+_(.+)_Phospho.*$", "\\1", first_part)
            if (identical(sites_part, first_part)) {
              # MaxQuant format: ACC_S40_M1 captures "S40"
              sites_part <- sub("^.*_([A-Za-z][0-9]+(?:[+][A-Za-z][0-9]+)*)_M[0-9]+$",
                                "\\1", first_part)
            }
            if (nzchar(sites_part) && !identical(sites_part, first_part))
              paste0(gene, "_", sites_part) else gene
          } else sid
        }, character(1))
        df_out
      }, error = function(e) {
        showNotification(paste0("FC x FC scatter failed: ", e$message), type = "warning", duration = 8)
        logger(paste0("PTM DE [", .sw_run_lbl, " PTM]: FC x FC scatter failed \u2014 ", e$message))
        NULL
      })
    })

    if (!is.null(fc_df)) {
      ptm_rv$fc_scatter_data <- list(
        df         = fc_df,
        comp_label = comparison,
        method_lbl = paste0(.sw_run_lbl, " PTM (Perseus LOD)"),
        n_dropped_no_global = attr(fc_df, "n_dropped_no_global") %||% 0L,
        n_acc_ambiguous     = sum(fc_df$acc_ambiguous, na.rm = TRUE)
      )
    }

    .enable_run_btn()
  })
  # ---------------------------------------------------------------------------
  # View mode toggle: switch between Adjusted / Absolute PTM / Protein
  # ---------------------------------------------------------------------------
  observeEvent(input$ptm_view_mode, {
    mode <- input$ptm_view_mode
    ptm_rv$de_view_mode <- mode

    chosen <- switch(mode,
      "adjusted"   = ptm_rv$de_result_adjusted,
      "unadjusted" = ptm_rv$de_result_unadjusted,
      "protein"    = ptm_rv$de_result_protein,
      ptm_rv$de_result_adjusted
    )
    mode_label <- PTM_LAYER_LABELS[[mode]]

    if (is.null(chosen)) {
      showNotification(
        paste0("No data available for view mode: '", mode_label, "'. Run PTM DE first."),
        type = "warning", duration = 6
      )
      return()
    }

    ptm_rv$de_result <- chosen

    n_sig <- sum(chosen$adj.pvalue < 0.05, na.rm = TRUE)
    ptm_rv$status <- paste0(
      "[", mode_label, "] ", ptm_rv$method_used, ": ",
      nrow(chosen), " PTM sites, ",
      n_sig, " significant (adj.p < 0.05)"
    )
  }, ignoreInit = TRUE)

  # ---------------------------------------------------------------------------
  # Pre-compute hover text whenever de_result changes (avoids O(n×k) in render)
  # ---------------------------------------------------------------------------
  observeEvent(ptm_rv$de_result, {
    df <- ptm_rv$de_result
    if (is.null(df) || nrow(df) == 0) {
      ptm_rv$de_result_prepped <- NULL
      return()
    }
    df$adj.pvalue  <- pmax(df$adj.pvalue, .Machine$double.xmin)
    df$neg_log10_p <- -log10(df$adj.pvalue)
    site_h <- ifelse(!is.na(df$Site) & nzchar(df$Site),
                     paste0("<br>Site: ", df$Site), "")
    mod_h  <- ifelse(!is.na(df$Modification) & nzchar(df$Modification),
                     paste0("<br>Modification: ", df$Modification), "")
    df$hover <- paste0(
      "<b>", df$PTM_Label, "</b>",
      "<br>Protein: ", df$Protein,
      site_h, mod_h,
      "<br>log2FC: ", round(df$log2FC, 3),
      "<br>adj.p: ", signif(df$adj.pvalue, 3),
      "<br>-log10(adj.p): ", round(df$neg_log10_p, 2)
    )
    # Razor assignment is visible in the results table but was invisible in the
    # plot: a shared site is tested once, under its razor protein, and the other
    # master proteins live in Also_Maps_To. Without this a reader hovering the
    # point cannot tell the site is shared, or with what. Guarded on the column,
    # not on the route: .annotate_prot_view() does not create it.
    if ("Also_Maps_To" %in% colnames(df)) {
      .alt_h <- ifelse(!is.na(df$Also_Maps_To) & nzchar(df$Also_Maps_To),
                       paste0("<br>Also maps to: ", df$Also_Maps_To,
                              "<br><i>shared peptide - one measurement, protein assignment is an inference</i>"),
                       "")
      df$hover <- paste0(df$hover, .alt_h)
    }
    det_cols <- grep("^Detected_", colnames(df), value = TRUE)
    if (length(det_cols) > 0) {
      det_text <- character(nrow(df))
      for (dc in det_cols) {
        cond_name <- sub("^Detected_", "", dc)
        det_text <- paste0(det_text, "<br>Detected (", cond_name, "): ", df[[dc]])
      }
      df$hover <- paste0(df$hover, det_text)
    }
    # DIA PTM path: Found_<group> columns (obs/total per group, before imputation)
    found_cols <- grep("^Found_", colnames(df), value = TRUE)
    if (length(found_cols) > 0) {
      found_text <- character(nrow(df))
      for (fc_col in found_cols) {
        cond_name  <- sub("^Found_", "", fc_col)
        found_text <- paste0(found_text, "<br>Found (", cond_name, "): ", df[[fc_col]])
      }
      df$hover <- paste0(df$hover, found_text)
    }
    # High missingness warning (>75% missing in any group)
    if ("high_missingness" %in% colnames(df)) {
      hi_flag <- !is.na(df$high_missingness) & df$high_missingness
      df$hover[hi_flag] <- paste0(
        df$hover[hi_flag],
        "<br><b>\u26a0 High Missingness (&gt;75% missing in a group).",
        "<br>Validate using FC\u00d7FC graph and single PTM analysis.</b>"
      )
    }
    ptm_rv$de_result_prepped <- df
  }, ignoreNULL = FALSE)

  # ---------------------------------------------------------------------------
  # Volcano plot
  # ---------------------------------------------------------------------------
  output$ptm_volcano_plotly <- renderPlotly({
    req(ptm_rv$de_result_prepped)
    req(isTRUE(input$ptm_de_tabs == "volcano"))

    df <- ptm_rv$de_result_prepped

    # --- Apply data filters ---
    filter_conds <- ptm_rv$filter_conditions
    if (length(filter_conds) > 0) {
      keep_mask <- apply_filter_conditions(df, filter_conds, "and_or")
      df <- df[keep_mask, , drop = FALSE]
    }

    fc_thresh <- input$ptm_fc_thresh
    p_thresh  <- input$ptm_p_thresh

    df$Significance <- "Not Significant"
    is_sig <- abs(df$log2FC) >= fc_thresh & df$adj.pvalue < p_thresh
    df$Significance[is_sig & df$log2FC > 0] <- "Up"
    df$Significance[is_sig & df$log2FC < 0] <- "Down"

    n_label <- input$ptm_n_labels
    sig_df <- df[df$Significance %in% c("Up", "Down"), , drop = FALSE]
    sig_df <- sig_df[order(sig_df$adj.pvalue), , drop = FALSE]
    label_targets <- head(sig_df$PTM_Label, n_label)

    # --- Build combined highlight mask ---
    highlight_gene  <- ptm_rv$highlight_gene
    highlight_ptm   <- ptm_rv$highlight_ptm
    exact_match     <- isTRUE(input$ptm_exact_match)
    highlight_rules <- ptm_rv$highlight_rules

    hl_mask <- rep(FALSE, nrow(df))
    if (!is.null(highlight_ptm) && nchar(highlight_ptm) > 0) {
      if (exact_match) hl_mask <- hl_mask | (toupper(df$PTM_Label) == toupper(highlight_ptm))
      else hl_mask <- hl_mask | grepl(highlight_ptm, df$PTM_Label, ignore.case = TRUE)
    } else if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      # Match the alternates as well as the razor winner. A shared site is stored
      # under one protein with the others in Also_Maps_To, so searching for an
      # alternate previously highlighted nothing and looked like an absence.
      # Exact mode compares against each comma-separated alternate individually,
      # so "KRT16" matches "KRT16, KRT17" without matching a longer name.
      .alt_v <- if ("Also_Maps_To" %in% colnames(df)) {
        .a <- as.character(df$Also_Maps_To); .a[is.na(.a)] <- ""; .a
      } else rep("", nrow(df))
      if (exact_match) {
        .alt_hit <- vapply(strsplit(.alt_v, ",", fixed = TRUE), function(v)
          any(toupper(trimws(v)) == toupper(highlight_gene)), logical(1))
        hl_mask <- hl_mask | (toupper(df$Gene) == toupper(highlight_gene)) | .alt_hit
      } else {
        hl_mask <- hl_mask | grepl(highlight_gene, df$Gene, ignore.case = TRUE) |
                   grepl(highlight_gene, .alt_v, ignore.case = TRUE)
      }
    }
    if (length(highlight_rules) > 0) {
      hl_mask <- hl_mask | apply_filter_conditions(df, highlight_rules, "or")
    }

    # High-missingness flag (>75% of samples missing in any group, computed pre-imputation;
    # these rows are above the k-NN eligibility cap and are never imputed)
    df$heavy_flag <- if ("high_missingness" %in% colnames(df))
      !is.na(df$high_missingness) & df$high_missingness
    else
      rep(FALSE, nrow(df))

    p <- plot_ly()

    if (any(hl_mask)) {

      non_match <- df[!hl_mask & !df$heavy_flag, , drop = FALSE]
      if (nrow(non_match) > 0) {
        p <- add_trace(
          p, data = non_match,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "markers",
          name = "Other",
          marker = list(size = 5, color = "#d4d4d4", opacity = 0.35,
                        line = list(width = 0.5, color = "#999")),
          hoverinfo = "text", text = ~hover
        )
      }

      hit_df <- df[hl_mask & !df$heavy_flag, , drop = FALSE]
      if (nrow(hit_df) > 0) {
        hl_label <- if (!is.null(highlight_ptm) && nchar(highlight_ptm) > 0)
                      paste0("Highlighted (incl. ", highlight_ptm, ")")
                    else if (!is.null(highlight_gene) && nchar(highlight_gene) > 0)
                      paste0("Highlighted (incl. ", highlight_gene, ")")
                    else "Highlighted"
        p <- add_trace(
          p, data = hit_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "markers",
          name = hl_label,
          marker = list(size = 11, color = "#e74c3c", opacity = 1,
                        line = list(width = 1, color = "#333")),
          hoverinfo = "text", text = ~hover
        )

        p <- add_trace(
          p, data = hit_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "text",
          text = ~PTM_Label,
          textposition = "top center",
          textfont = list(size = 11, color = "#c0392b", family = "Arial Black"),
          showlegend = FALSE, hoverinfo = "none"
        )
      }

      # High Misngness points (gray triangles, regardless of highlight state)
      heavy_df <- df[df$heavy_flag, , drop = FALSE]
      if (nrow(heavy_df) > 0) {
        p <- add_trace(
          p, data = heavy_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "markers",
          name = "High Missingness",
          marker = list(size = 7, color = "#999999", symbol = "triangle-up",
                        opacity = 0.75, line = list(width = 0.5, color = "#555")),
          hoverinfo = "text", text = ~hover
        )
      }
      # Text labels for heavy-imputed highlighted points
      heavy_hit_df <- df[df$heavy_flag & hl_mask, , drop = FALSE]
      if (nrow(heavy_hit_df) > 0) {
        p <- add_trace(
          p, data = heavy_hit_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "text",
          text = ~PTM_Label,
          textposition = "top center",
          textfont = list(size = 11, color = "#888", family = "Arial Black"),
          showlegend = FALSE, hoverinfo = "none"
        )
      }

    } else {
      sig_colors <- c("Not Significant" = "#BBBBBB", "Down" = "#3498DB", "Up" = "#E74C3C")
      for (grp in c("Not Significant", "Down", "Up")) {
        grp_df <- df[df$Significance == grp & !df$heavy_flag, , drop = FALSE]
        if (nrow(grp_df) == 0) next

        p <- add_trace(
          p, data = grp_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "markers",
          name = grp,
          marker = list(size = 5, color = sig_colors[[grp]], opacity = 0.6,
                        line = list(width = 0.3, color = "#555")),
          hoverinfo = "text", text = ~hover
        )
      }

      # Heavy imputation points: gray triangles regardless of significance
      heavy_df <- df[df$heavy_flag, , drop = FALSE]
      if (nrow(heavy_df) > 0) {
        p <- add_trace(
          p, data = heavy_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "markers",
          name = "High Missingness",
          marker = list(size = 7, color = "#999999", symbol = "triangle-up",
                        opacity = 0.75, line = list(width = 0.5, color = "#555")),
          hoverinfo = "text", text = ~hover
        )
      }

      lab_df <- df[df$PTM_Label %in% label_targets & !df$heavy_flag, , drop = FALSE]
      if (nrow(lab_df) > 0) {
        p <- add_trace(
          p, data = lab_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "text",
          text = ~PTM_Label,
          textposition = "top center",
          textfont = list(size = 10, color = "#333", family = "Arial"),
          showlegend = FALSE, hoverinfo = "none"
        )
      }
      # Labels for heavy-imputed points that rank in top-N by significance
      lab_heavy_df <- df[df$PTM_Label %in% label_targets & df$heavy_flag, , drop = FALSE]
      if (nrow(lab_heavy_df) > 0) {
        p <- add_trace(
          p, data = lab_heavy_df,
          x = ~log2FC, y = ~neg_log10_p,
          type = "scatter", mode = "text",
          text = ~PTM_Label,
          textposition = "top center",
          textfont = list(size = 10, color = "#999", family = "Arial"),
          showlegend = FALSE, hoverinfo = "none"
        )
      }
    }

    # Build directional arrow annotations from comparison string
    comp_parts <- strsplit(ptm_rv$comparison, " vs ")[[1]]
    cond_test <- if (length(comp_parts) >= 1) trimws(comp_parts[1]) else "Test"
    cond_ref  <- if (length(comp_parts) >= 2) trimws(comp_parts[2]) else "Reference"

    direction_annots <- list(
      list(
        x = 0.01, y = -0.08, xref = "paper", yref = "paper",
        text = paste0("\u2190 Higher in <b>", cond_ref, "</b>"),
        showarrow = FALSE,
        font = list(size = 12, color = "#3498DB"),
        xanchor = "left"
      ),
      list(
        x = 0.99, y = -0.08, xref = "paper", yref = "paper",
        text = paste0("Higher in <b>", cond_test, "</b> \u2192"),
        showarrow = FALSE,
        font = list(size = 12, color = "#E74C3C"),
        xanchor = "right"
      )
    )

    # Dynamic title and axis label based on view mode
    view_mode <- ptm_rv$de_view_mode
    view_suffix <- switch(view_mode,
      "adjusted"   = "Adjusted to Total Protein",
      "unadjusted" = "Absolute PTM",
      "protein"    = "Protein View",
      "Adjusted to Total Protein"
    )
    x_label <- switch(view_mode,
      "adjusted"   = "PTM log2 Fold Change (adjusted for protein)",
      "unadjusted" = "PTM log2 Fold Change (absolute, unadjusted)",
      "protein"    = "Protein log2 Fold Change",
      "PTM log2 Fold Change (adjusted for protein)"
    )

    p <- layout(
      p,
      title = paste0("PTM Volcano: ", ptm_rv$comparison,
                      " \u2014 ", view_suffix,
                      " (", ptm_rv$method_used, ")",
                      ""),
      xaxis = list(title = x_label),
      yaxis = list(title = "-log10(Adj. P-Value)"),
      shapes = list(
        list(type = "line", x0 = -fc_thresh, x1 = -fc_thresh,
             y0 = 0, y1 = 1, yref = "paper",
             line = list(dash = "dash", color = "grey50", width = 1)),
        list(type = "line", x0 = fc_thresh, x1 = fc_thresh,
             y0 = 0, y1 = 1, yref = "paper",
             line = list(dash = "dash", color = "grey50", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = -log10(p_thresh), y1 = -log10(p_thresh),
             line = list(dash = "dash", color = "grey50", width = 1))
      ),
      annotations = direction_annots,
      legend = list(orientation = "h", x = 0.3, y = -0.2),
      margin = list(b = 100)
    )

    p
  })

  # ---------------------------------------------------------------------------
  # FC x FC scatter — triggered by "Run FC x FC Analysis" button.
  # Data computation (LOD imputation, labelling, classification) happens in
  # observeEvent(input$ptm_run_btn) FC branch; result stored in ptm_rv$fc_scatter_data.
  # ---------------------------------------------------------------------------
  output$ptm_fc_scatter_plotly <- renderPlotly({
    req(!is.null(ptm_rv$fc_scatter_data), !is.null(ptm_rv$fc_scatter_data$df))

    fc_data    <- ptm_rv$fc_scatter_data
    df         <- fc_data$df
    comp_label <- fc_data$comp_label %||% ""
    method_lbl <- fc_data$method_lbl %||% ""

    PTM_THRESH <- 1.0
    GLB_THRESH <- 0.5

    # ── 9-category color scheme ─────────────────────────────────────────────
    # PTM threshold: ±1.0 log2FC  |  Global threshold: ±0.5 log2FC
    ptm_up   <- df$log2fc >  PTM_THRESH
    ptm_down <- df$log2fc < -PTM_THRESH
    glb_up   <- df$global_log2fc >  GLB_THRESH
    glb_down <- df$global_log2fc < -GLB_THRESH

    df$category <- "NS in Both"
    df$category[ptm_down & glb_down]                  <- "Both Down"
    df$category[ptm_down & !glb_up & !glb_down]       <- "PTM Down / Global NS"
    df$category[!ptm_up & !ptm_down & glb_up]         <- "PTM NS / Global Up"
    df$category[ptm_up  & !glb_up & !glb_down]        <- "PTM Up / Global NS"
    df$category[ptm_up  & glb_up]                     <- "Both Up"
    df$category[ptm_down & glb_up]                    <- "PTM Down / Global Up"
    df$category[!ptm_up & !ptm_down & glb_down]       <- "PTM NS / Global Down"
    df$category[ptm_up  & glb_down]                   <- "PTM Up / Global Down"

    cat_colors <- c(
      "Both Down"             = "#1565C0",
      "PTM Down / Global NS"  = "#64B5F6",
      "NS in Both"            = "#BDBDBD",
      "PTM NS / Global Up"    = "#FFCC80",
      "PTM Up / Global NS"    = "#EF6C00",
      "Both Up"               = "#B71C1C",
      "PTM Down / Global Up"  = "#00897B",
      "PTM NS / Global Down"  = "#B0BEC5",
      "PTM Up / Global Down"  = "#8E24AA"
    )

    # ── Shape encoding: imputation status ───────────────────────────────────
    shape_map <- c(
      "none"            = "circle",
      "partial"         = "triangle-up",
      "dropout_one_cond" = "star"
    )
    df$marker_symbol <- shape_map[df$imp_class]
    df$marker_symbol[is.na(df$marker_symbol)] <- "circle"

    # ── Hover text ───────────────────────────────────────────────────────────
    imp_label <- c(
      "none"             = "None",
      "partial"          = "Partial LOD imputation",
      "dropout_one_cond" = "Condition-wide dropout (LOD filled)"
    )
    df$hover <- paste0(
      "<b>", df$PTM_Label, "</b>",
      "<br>PTM log2FC: ",    round(df$log2fc, 3),
      "<br>Global log2FC: ", round(df$global_log2fc, 3),
      "<br>Category: ",      df$category,
      "<br>Imputation: ",    imp_label[df$imp_class],
      "<br>Obs (test): ",    df$n_obs_test, "/", df$n_total_test,
      "<br>Obs (ref): ",     df$n_obs_ref,  "/", df$n_total_ref,
      ifelse(if (is.null(df$acc_ambiguous)) FALSE
             else !is.na(df$acc_ambiguous) & df$acc_ambiguous,
             "<br><i>Shared site: global log2FC is that of the first listed accession</i>",
             "")
    )

    # ── Highlight logic ──────────────────────────────────────────────────────
    highlight_gene  <- ptm_rv$highlight_gene
    highlight_ptm   <- ptm_rv$highlight_ptm
    exact_match     <- isTRUE(input$ptm_exact_match)
    highlight_rules <- ptm_rv$highlight_rules

    hl_mask <- rep(FALSE, nrow(df))
    if (!is.null(highlight_ptm) && nzchar(highlight_ptm %||% "")) {
      hl_mask <- if (exact_match) toupper(df$PTM_Label) == toupper(highlight_ptm)
                 else grepl(highlight_ptm, df$PTM_Label, ignore.case = TRUE)
    } else if (!is.null(highlight_gene) && nzchar(highlight_gene %||% "")) {
      hl_mask <- if (exact_match) toupper(df$PTM_Label) == toupper(highlight_gene)
                 else grepl(highlight_gene, df$PTM_Label, ignore.case = TRUE)
    }
    if (length(highlight_rules) > 0)
      hl_mask <- hl_mask | apply_filter_conditions(df, highlight_rules, "or")

    # ── Label top-N non-NS points ────────────────────────────────────────────
    n_label      <- input$ptm_n_labels %||% 15L
    label_df     <- df[df$category != "NS in Both", , drop = FALSE]
    label_df     <- label_df[order(abs(label_df$log2fc), decreasing = TRUE), ]
    label_targets <- head(label_df$PTM_Label, n_label)

    # ── Plot ─────────────────────────────────────────────────────────────────
    p <- plot_ly()

    if (any(hl_mask)) {
      non_hl <- df[!hl_mask, , drop = FALSE]
      if (nrow(non_hl) > 0)
        p <- add_trace(p, data = non_hl,
                       x = ~global_log2fc, y = ~log2fc,
                       type = "scatter", mode = "markers",
                       name = "Other",
                       marker = list(size = 6, color = "#d4d4d4", opacity = 0.35,
                                     symbol = ~marker_symbol,
                                     line = list(width = 0.5, color = "#999")),
                       hoverinfo = "text", text = ~hover)
      hit_df <- df[hl_mask, , drop = FALSE]
      if (nrow(hit_df) > 0) {
        hl_lbl <- paste0("Highlighted (",
                         if (!is.null(highlight_ptm) && nzchar(highlight_ptm %||% ""))
                           highlight_ptm else highlight_gene, ")")
        p <- add_trace(p, data = hit_df,
                       x = ~global_log2fc, y = ~log2fc,
                       type = "scatter", mode = "markers",
                       name = hl_lbl,
                       marker = list(size = 12, color = "#e74c3c", opacity = 1,
                                     symbol = ~marker_symbol,
                                     line = list(width = 1, color = "#333")),
                       hoverinfo = "text", text = ~hover)
        p <- add_trace(p, data = hit_df,
                       x = ~global_log2fc, y = ~log2fc,
                       type = "scatter", mode = "text",
                       text = ~PTM_Label, textposition = "top center",
                       textfont = list(size = 11, color = "#c0392b", family = "Arial Black"),
                       showlegend = FALSE, hoverinfo = "none")
      }
    } else {
      # Plot by category (color) with shape per imputation class within each
      cat_order <- c("Both Down", "PTM Down / Global NS", "NS in Both",
                     "PTM NS / Global Up", "PTM Up / Global NS", "Both Up",
                     "PTM Down / Global Up", "PTM NS / Global Down", "PTM Up / Global Down")
      for (cat in cat_order) {
        cat_df <- df[df$category == cat, , drop = FALSE]
        if (nrow(cat_df) == 0) next
        p <- add_trace(p, data = cat_df,
                       x = ~global_log2fc, y = ~log2fc,
                       type = "scatter", mode = "markers",
                       name = cat,
                       marker = list(
                         size   = 7,
                         color  = cat_colors[[cat]],
                         symbol = ~marker_symbol,
                         opacity = if (cat == "NS in Both") 0.4 else 0.8,
                         line   = list(width = 0.5, color = "#555")
                       ),
                       hoverinfo = "text", text = ~hover)
      }
      lab_df <- df[df$PTM_Label %in% label_targets, , drop = FALSE]
      if (nrow(lab_df) > 0)
        p <- add_trace(p, data = lab_df,
                       x = ~global_log2fc, y = ~log2fc,
                       type = "scatter", mode = "text",
                       text = ~PTM_Label, textposition = "top center",
                       textfont = list(size = 10, color = "#222", family = "Arial"),
                       showlegend = FALSE, hoverinfo = "none")
    }

    # ── Threshold reference lines and axis ────────────────────────────────────
    all_fc    <- c(df$global_log2fc, df$log2fc)
    fc_range  <- range(all_fc, na.rm = TRUE)
    pad       <- diff(fc_range) * 0.1
    ax_lims   <- c(fc_range[1] - pad, fc_range[2] + pad)

    p <- layout(p,
      title = paste0("FC Graph (with Perseus LOD Imputation) \u2014 ",
                     comp_label, " (", method_lbl, ")"),
      xaxis = list(title = paste0("Global Protein log2FC  [threshold \u00b1", GLB_THRESH, "]"),
                   zeroline = TRUE, range = ax_lims),
      yaxis = list(title = paste0("Absolute PTM log2FC  [threshold \u00b1", PTM_THRESH, "]"),
                   zeroline = TRUE, range = ax_lims),
      shapes = list(
        list(type = "line", x0 = ax_lims[1], x1 = ax_lims[2],
             y0 = -PTM_THRESH, y1 = -PTM_THRESH,
             line = list(dash = "dot", color = "#888", width = 1)),
        list(type = "line", x0 = ax_lims[1], x1 = ax_lims[2],
             y0 =  PTM_THRESH, y1 =  PTM_THRESH,
             line = list(dash = "dot", color = "#888", width = 1)),
        list(type = "line", x0 = -GLB_THRESH, x1 = -GLB_THRESH,
             y0 = ax_lims[1], y1 = ax_lims[2],
             line = list(dash = "dot", color = "#888", width = 1)),
        list(type = "line", x0 =  GLB_THRESH, x1 =  GLB_THRESH,
             y0 = ax_lims[1], y1 = ax_lims[2],
             line = list(dash = "dot", color = "#888", width = 1))
      ),
      annotations = list(
        list(x = ax_lims[2] * 0.98, y = ax_lims[1] * 0.98,
             text = paste0("\u25cf No imputation  \u25b2 Partial LOD  \u2605 Condition dropout",
                           "<br>Categories: |PTM log2FC| > ", PTM_THRESH,
                           ", |protein log2FC| > ", GLB_THRESH,
                           " \u2014 fold change only \u2014 no significance test",
                           "<br>", fc_data$n_dropped_no_global %||% 0,
                           " site(s) omitted (no matched protein log2FC); ",
                           fc_data$n_acc_ambiguous %||% 0,
                           " shared site(s) use the first listed accession"),
             showarrow = FALSE, xanchor = "right", yanchor = "bottom",
             font = list(size = 10, color = "#555"))
      ),
      legend = list(orientation = "h", x = 0, y = -0.18, font = list(size = 10))
    )
    p
  })

  # ---------------------------------------------------------------------------
  # Build master table: merge adjusted, absolute, protein into one wide table
  # Shared cols appear once; model-specific cols get prefixed.
  # primary = ptm_primary_layer() of the cache: "adjusted" (rows keyed on the
  # adjusted table), or "unadjusted" for a PTM-only cache, whose master table is
  # the absolute table alone (Absolute_* columns; no Adjusted_ / ProteinOnly_).
  # ---------------------------------------------------------------------------
  .build_master_table <- function(adj_df, unadj_df, prot_df, adj_map = NULL,
                                  primary = "adjusted") {
    key_df <- if (identical(primary, "unadjusted")) unadj_df else adj_df
    if (is.null(key_df) || nrow(key_df) == 0) return(NULL)

    shared_cols <- c("PTM_ID", "Gene", "Also_Maps_To", "Is_Shared", "Protein", "Site", "Modification",
                     "PTM_Label", "Comparison")
    loc_cols <- intersect(c("Localization_Prob", "Site_ID"), colnames(key_df))
    det_cols <- grep("^Detected_", colnames(key_df), value = TRUE)
    shared_cols <- c(shared_cols, loc_cols, det_cols)
    shared_cols <- intersect(shared_cols, colnames(key_df))

    model_cols <- c("log2FC", "SE", "DF", "pvalue", "adj.pvalue")

    if (identical(primary, "unadjusted")) {
      master <- key_df[, intersect(c(shared_cols, model_cols), colnames(key_df)), drop = FALSE]
      abs_rename <- intersect(model_cols, colnames(master))
      names(master)[names(master) %in% abs_rename] <-
        paste0("Absolute_", names(master)[names(master) %in% abs_rename])
      master <- master[order(master$Absolute_adj.pvalue), , drop = FALSE]
      rownames(master) <- NULL
      return(master)
    }

    # Start with adjusted (primary key = PTM_ID)
    master <- adj_df[, intersect(c(shared_cols, model_cols), colnames(adj_df)), drop = FALSE]
    adj_rename <- intersect(model_cols, colnames(master))
    names(master)[names(master) %in% adj_rename] <-
      paste0("Adjusted_", names(master)[names(master) %in% adj_rename])

    # Merge unadjusted (absolute PTM) on the ABSOLUTE site ID.
    # The adjusted PTM_ID is the adjusted-matrix row name (for shared sites
    # "<per-protein site>_SharedInference") while the absolute PTM_ID is the
    # semicolon-joined absolute site ID, so the two never match directly.
    # Resolve through adj_to_abs_map (falling back to PTM_ID when unmapped).
    n_key_unmapped <- NA_integer_
    if (!is.null(unadj_df) && nrow(unadj_df) > 0) {
      master$Absolute_key <- if (!is.null(adj_map)) {
        .res <- unname(adj_map[master$PTM_ID])
        # Count the silent fallbacks. If the map is ever incomplete these rows
        # revert to NA Absolute_* statistics - the defect this merge key fixed,
        # reappearing without a signal. Behaviour is unchanged; only counted.
        n_key_unmapped <- sum(is.na(.res))
        if (n_key_unmapped > 0L)
          logger(paste0("PTM DE: ", n_key_unmapped, " of ", nrow(master),
                        " adjusted row(s) did not resolve through adj_to_abs_map ",
                        "- fell back to PTM_ID for the Absolute merge key."))
        ifelse(is.na(.res), master$PTM_ID, .res)
      } else {
        master$PTM_ID
      }
      unadj_sub <- unadj_df[, intersect(c("PTM_ID", model_cols), colnames(unadj_df)), drop = FALSE]
      unadj_mcols <- intersect(model_cols, colnames(unadj_sub))
      names(unadj_sub)[names(unadj_sub) %in% unadj_mcols] <-
        paste0("Absolute_", names(unadj_sub)[names(unadj_sub) %in% unadj_mcols])
      names(unadj_sub)[names(unadj_sub) == "PTM_ID"] <- "Absolute_key"
      unadj_sub <- unadj_sub[!duplicated(unadj_sub$Absolute_key), , drop = FALSE]
      master <- merge(master, unadj_sub, by = "Absolute_key", all.x = TRUE, sort = FALSE)
      master$Absolute_key <- NULL
    } else {
      for (mc in model_cols) master[[paste0("Absolute_", mc)]] <- NA_real_
    }

    # Merge protein on Protein (accession) — duplicates across PTM sites
    if (!is.null(prot_df) && nrow(prot_df) > 0) {
      prot_sub <- prot_df[, intersect(c("Protein", model_cols), colnames(prot_df)), drop = FALSE]
      # Protein model PTM_ID = accession, map via Protein column
      prot_mcols <- intersect(model_cols, colnames(prot_sub))
      names(prot_sub)[names(prot_sub) %in% prot_mcols] <-
        paste0("ProteinOnly_", names(prot_sub)[names(prot_sub) %in% prot_mcols])
      # De-dup protein rows (one per accession)
      prot_sub <- prot_sub[!duplicated(prot_sub$Protein), , drop = FALSE]
      master <- merge(master, prot_sub, by = "Protein", all.x = TRUE, sort = FALSE)
    } else {
      for (mc in model_cols) master[[paste0("ProteinOnly_", mc)]] <- NA_real_
    }

    # Order: shared cols first, then Adjusted_, Absolute_, ProteinOnly_
    adj_cols  <- grep("^Adjusted_", colnames(master), value = TRUE)
    abs_cols  <- grep("^Absolute_", colnames(master), value = TRUE)
    prot_cols <- grep("^ProteinOnly_", colnames(master), value = TRUE)
    ordered <- c(intersect(shared_cols, colnames(master)), adj_cols, abs_cols, prot_cols)
    extra <- setdiff(colnames(master), ordered)
    master <- master[, c(ordered, extra), drop = FALSE]
    master <- master[order(master$Adjusted_adj.pvalue), , drop = FALSE]
    rownames(master) <- NULL
    attr(master, "n_absolute_key_unmapped") <- n_key_unmapped
    master
  }

  # ---------------------------------------------------------------------------
  # Table and status
  # ---------------------------------------------------------------------------
  output$ptm_de_table <- DT::renderDT({
    req(ptm_rv$de_master)
    req(isTRUE(input$ptm_de_tabs == "de_table"))

    display_df <- ptm_rv$de_master

    # Round numeric columns for display
    num_cols <- sapply(display_df, is.numeric)
    for (col in names(display_df)[num_cols]) {
      if (grepl("pvalue", col, ignore.case = TRUE)) {
        display_df[[col]] <- signif(display_df[[col]], 3)
      } else {
        display_df[[col]] <- round(display_df[[col]], 4)
      }
    }

    DT::datatable(
      display_df,
      filter = "top",
      rownames = FALSE,
      caption = NULL,
      options = list(pageLength = 25, scrollX = TRUE),
      class = "display compact"
    )
  })


  output$ptm_status <- renderText({
    ptm_rv$status
  })

  # ---------------------------------------------------------------------------
  # Highlight
  # ---------------------------------------------------------------------------
  observeEvent(input$ptm_gene_search_btn, {
    req(!is.null(ptm_rv$de_result) || !is.null(ptm_rv$fc_scatter_data))

    gene <- trimws(if (is.null(input$ptm_gene_search)) "" else input$ptm_gene_search)
    if (nchar(gene) == 0) {
      ptm_rv$highlight_gene <- NULL
      ptm_rv$highlight_ptm  <- NULL
      ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
      showNotification("Highlight cleared.", type = "message")
      return()
    }

    ptm_rv$highlight_ptm <- NULL          # clear PTM-level search when gene search is used

    # Search in whichever result is available; FC mode: search PTM_Label
    search_df  <- if (!is.null(ptm_rv$de_result)) ptm_rv$de_result else ptm_rv$fc_scatter_data$df
    search_col <- if (!is.null(ptm_rv$de_result) && "Gene" %in% colnames(ptm_rv$de_result)) "Gene" else "PTM_Label"
    if (search_col %in% colnames(search_df)) {
      hits <- if (isTRUE(input$ptm_exact_match))
        which(tolower(search_df[[search_col]]) == tolower(gene))
      else
        grep(gene, search_df[[search_col]], ignore.case = TRUE)
    } else {
      hits <- integer(0)
    }

    if (length(hits) == 0) {
      showNotification(paste0("Gene '", gene, "' not found."), type = "warning")
      ptm_rv$highlight_gene <- NULL
    } else {
      ptm_rv$highlight_gene <- gene
      showNotification(paste0("Found ", length(hits), " match(es) for '", gene, "'."), type = "message")
      ptm_rv$status <- paste0("Highlighted: ", gene, " (", length(hits), " PTM rows)")
    }

    ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
  })

  # PTM-level highlight
  observeEvent(input$ptm_ptm_search_btn, {
    req(!is.null(ptm_rv$de_result) || !is.null(ptm_rv$fc_scatter_data))

    ptm_q <- trimws(if (is.null(input$ptm_ptm_search)) "" else input$ptm_ptm_search)
    if (nchar(ptm_q) == 0) {
      ptm_rv$highlight_ptm  <- NULL
      ptm_rv$highlight_gene <- NULL
      ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
      showNotification("PTM highlight cleared.", type = "message")
      return()
    }

    ptm_rv$highlight_gene <- NULL          # clear gene-level search when PTM search is used

    # Search PTM_Label in whichever result is available
    search_df <- if (!is.null(ptm_rv$de_result)) ptm_rv$de_result else ptm_rv$fc_scatter_data$df
    if ("PTM_Label" %in% colnames(search_df)) {
      hits <- if (isTRUE(input$ptm_exact_match))
        which(tolower(search_df$PTM_Label) == tolower(ptm_q))
      else
        grep(ptm_q, search_df$PTM_Label, ignore.case = TRUE)
    } else {
      hits <- integer(0)
    }

    if (length(hits) == 0) {
      showNotification(paste0("PTM '", ptm_q, "' not found."), type = "warning")
      ptm_rv$highlight_ptm <- NULL
    } else {
      ptm_rv$highlight_ptm <- ptm_q
      showNotification(paste0("Found ", length(hits), " match(es) for PTM '", ptm_q, "'."), type = "message")
      ptm_rv$status <- paste0("Highlighted PTM: ", ptm_q, " (", length(hits), " rows)")
    }

    ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
  })

  # =========================================================================
  # HIGHLIGHT RULES — add / remove / clear / render  (PTM-DE)
  # =========================================================================
  observeEvent(input$ptm_de_add_highlight, {
    result_cols <-
      if (!is.null(ptm_rv$de_result))          colnames(ptm_rv$de_result)
      else if (!is.null(ptm_rv$fc_scatter_data)) colnames(ptm_rv$fc_scatter_data$df)
      else                                       NULL
    all_choices <- unique(result_cols)
    req(length(all_choices) > 0)

    showModal(modalDialog(
      title = tags$span(icon("highlighter"), " Add Highlight Rule"),
      tagList(
        helpText(class = "small text-muted mb-2",
          "Highlight PTMs matching this condition on the plot. ",
          "Multiple rules are combined with OR."),
        tags$div(class = "mb-3",
          tags$label(class = "form-label fw-bold", "Column"),
          tags$div(style = "max-height: 250px; overflow-y: auto; border: 1px solid #eee; padding: 6px; border-radius: 4px;",
            radioButtons("ptm_de_hl_column", label = NULL,
                         choices  = setNames(all_choices, all_choices),
                         selected = all_choices[1])
          )
        ),
        fluidRow(
          column(6, selectInput("ptm_de_hl_operator", "Operator",
                    choices = c(">" = ">", ">=" = ">=", "<" = "<", "<=" = "<=",
                                "==" = "==", "!=" = "!=",
                                "contains" = "contains",
                                "not contains" = "not_contains"))),
          column(6, textInput("ptm_de_hl_value", "Value", placeholder = "e.g. 0.01"))
        )
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ptm_de_hl_confirm", "Add Highlight", class = "btn-warning",
                     icon = icon("highlighter"))
      )
    ))
  })

  observeEvent(input$ptm_de_hl_confirm, {
    col_name <- input$ptm_de_hl_column
    op       <- input$ptm_de_hl_operator
    val      <- trimws(input$ptm_de_hl_value)
    if (is.null(col_name) || col_name == "" || val == "") {
      showNotification("Please select a column and enter a value.", type = "warning")
      return()
    }
    ptm_rv$highlight_counter <- ptm_rv$highlight_counter + 1L
    ptm_rv$highlight_rules <- c(ptm_rv$highlight_rules, list(list(
      id = ptm_rv$highlight_counter, column = col_name, op = op, value = val
    )))
    ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
    logger(paste0("PTM DE Highlight: Added ", col_name, " ", op, " ", val))
    removeModal()
  })

  observeEvent(lapply(seq_len(100), function(i) input[[paste0("ptm_de_hlrm_", i)]]), {
    for (i in seq_along(ptm_rv$highlight_rules)) {
      btn_id <- paste0("ptm_de_hlrm_", ptm_rv$highlight_rules[[i]]$id)
      btn_val <- input[[btn_id]]
      if (!is.null(btn_val) && btn_val > 0) {
        rule <- ptm_rv$highlight_rules[[i]]
        attr_name <- paste0(".ptm_de_hlrm_seen_", rule$id)
        last_seen <- isolate(ptm_rv[[attr_name]])
        if (is.null(last_seen) || btn_val > last_seen) {
          ptm_rv[[attr_name]] <- btn_val
          ptm_rv$highlight_rules <- ptm_rv$highlight_rules[
            vapply(ptm_rv$highlight_rules, function(r) r$id != rule$id, logical(1))
          ]
          ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
          logger(paste0("PTM DE Highlight: Removed #", rule$id))
          break
        }
      }
    }
  })

  observeEvent(input$ptm_de_clear_highlights, {
    ptm_rv$highlight_rules <- list()
    ptm_rv$render_count <- isolate(ptm_rv$render_count) + 1L
    logger("PTM DE Highlight: All cleared")
  })

  output$ptm_de_active_highlights <- renderUI({
    rules <- ptm_rv$highlight_rules
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
        actionButton(paste0("ptm_de_hlrm_", r$id), label = NULL,
                     icon = icon("xmark"),
                     class = "btn-outline-danger btn-sm",
                     style = "padding: 1px 5px; font-size: 0.7em;")
      )
    })
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

  # =========================================================================
  # DATA FILTERS — add / remove / clear / render  (PTM-DE)
  # =========================================================================
  observeEvent(input$ptm_de_open_filters, {
    filter_headers <-
      if (!is.null(rv$pd_protein_headers) && length(rv$pd_protein_headers) > 0)
        rv$pd_protein_headers
      else if (!is.null(ptm_rv$de_master))     colnames(ptm_rv$de_master)
      else if (!is.null(ptm_rv$de_result))     colnames(ptm_rv$de_result)
      else if (!is.null(ptm_rv$fc_scatter_data)) colnames(ptm_rv$fc_scatter_data$df)
      else NULL
    req(!is.null(filter_headers))
    showModal(build_filter_modal("ptm_de_f", filter_headers, ptm_rv$filter_conditions))
  })

  observeEvent(input$ptm_de_f_confirm, {
    col_name  <- input$ptm_de_f_column
    op        <- input$ptm_de_f_operator
    val       <- trimws(input$ptm_de_f_value)
    grp       <- input$ptm_de_f_group
    grp_logic <- input$ptm_de_f_group_logic
    if (is.null(grp_logic) || !grp_logic %in% c("OR", "AND")) grp_logic <- "OR"
    if (op %in% c("isNA", "notNA")) val <- ""
    if (is.null(col_name) || col_name == "") {
      showNotification("Please select a column.", type = "warning"); return()
    }
    if (!op %in% c("isNA", "notNA") && val == "") {
      showNotification("Please enter a value.", type = "warning"); return()
    }
    if (length(ptm_rv$filter_conditions) > 0) {
      ptm_rv$filter_conditions <- lapply(ptm_rv$filter_conditions, function(f) {
        if (f$group == grp) f$logic <- grp_logic; f
      })
    }
    ptm_rv$filter_counter <- ptm_rv$filter_counter + 1L
    ptm_rv$filter_conditions <- c(ptm_rv$filter_conditions, list(list(
      id = ptm_rv$filter_counter, group = grp, column = col_name,
      op = op, value = val, logic = grp_logic
    )))
    op_label <- switch(op, "isNA" = "is NA", "notNA" = "is not NA",
                        "contains" = "contains", "not_contains" = "not contains", op)
    logger(paste0("PTM DE Filter: Added [", grp, " (", grp_logic, ")] ",
                  col_name, " ", op_label,
                  if (nzchar(val)) paste0(" ", val) else ""))
    removeModal()
  })

  observeEvent(lapply(seq_len(100), function(i) input[[paste0("ptm_de_frm_", i)]]), {
    for (i in seq_along(ptm_rv$filter_conditions)) {
      btn_id <- paste0("ptm_de_frm_", ptm_rv$filter_conditions[[i]]$id)
      btn_val <- input[[btn_id]]
      if (!is.null(btn_val) && btn_val > 0) {
        cond <- ptm_rv$filter_conditions[[i]]
        attr_name <- paste0(".ptm_de_frm_seen_", cond$id)
        last_seen <- isolate(ptm_rv[[attr_name]])
        if (is.null(last_seen) || btn_val > last_seen) {
          ptm_rv[[attr_name]] <- btn_val
          ptm_rv$filter_conditions <- ptm_rv$filter_conditions[
            vapply(ptm_rv$filter_conditions, function(f) f$id != cond$id, logical(1))
          ]
          logger(paste0("PTM DE Filter: Removed #", cond$id))
          break
        }
      }
    }
  })

  observeEvent(input$ptm_de_clear_filters, {
    ptm_rv$filter_conditions <- list()
    logger("PTM DE Filter: All cleared")
  })

  output$ptm_de_active_filters <- renderUI({
    conds <- ptm_rv$filter_conditions
    if (length(conds) == 0)
      return(tags$p(class = "text-muted small", "No active filters."))
    render_filter_badges(conds, "ptm_de_frm_")
  })

  # ---------------------------------------------------------------------------
  # Downloads and RDS save/load
  # ---------------------------------------------------------------------------
  output$ptm_fc_data_ready <- reactive({ !is.null(ptm_rv$fc_scatter_data) })
  outputOptions(output, "ptm_fc_data_ready", suspendWhenHidden = FALSE)

  output$ptm_download_csv <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_DE_Master_", comp_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(ptm_rv$de_master)
      write.csv(ptm_rv$de_master, file, row.names = FALSE)
      logger("PTM DE: Master CSV downloaded")
    }
  )

  output$ptm_download_fc_csv <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_FCxFC_", comp_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(ptm_rv$fc_scatter_data)
      fc_df <- ptm_rv$fc_scatter_data$df
      out_cols <- intersect(c("PTM_Label", "feature", "log2fc", "global_log2fc",
                              "category", "imp_class",
                              "n_obs_test", "n_total_test", "n_obs_ref", "n_total_ref"),
                            colnames(fc_df))
      write.csv(fc_df[, out_cols, drop = FALSE], file, row.names = FALSE)
      logger("PTM FC: FC x FC CSV downloaded")
    }
  )

  output$ptm_save_rds <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_Analysis_", comp_safe, "_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(ptm_rv$de_result)
      save_data <- list(
        de_result = ptm_rv$de_result,
        de_result_adjusted   = ptm_rv$de_result_adjusted,
        de_result_unadjusted = ptm_rv$de_result_unadjusted,
        de_result_protein    = ptm_rv$de_result_protein,
        de_master  = ptm_rv$de_master,
        de_view_mode = ptm_rv$de_view_mode,
        fc_scatter_data = ptm_rv$fc_scatter_data,
        nmar_table = ptm_rv$nmar_table,
        nmar_has_power = ptm_rv$nmar_has_power,
        nmar_min_achievable_p = ptm_rv$nmar_min_achievable_p,
        comparison = ptm_rv$comparison,
        method_used = ptm_rv$method_used,
        n_files = ptm_rv$n_files,
        conditions = ptm_rv$conditions,
        cache_fingerprint = ptm_cache_fingerprint(rv$dia_ptm_master_cache),
        timestamp = Sys.time(),
        provenance = reactiveValuesToList(rv$provenance)
      )
      saveRDS(save_data, file)
      logger("PTM DE: RDS saved (DE + FC x FC)")
    }
  )

  # ===========================================================================
  # DETECTION TEST TABLE (Tier C: sites not testable by limma)
  # ===========================================================================
  # One table, two populations. These were previously two tabs that looked like
  # duplicates but are disjoint:
  #   (a) EXCLUDED - failed the per-comparison valid-value filter, so limma could
  #       not fit them; assessed instead by Fisher's exact test on detection
  #       counts, BH-corrected within this set only.
  #   (b) MODELLED - passed the filter and were fitted, but have a condition with
  #       zero detections or with more than 75% of samples missing.
  # On the DIA route the any-arm rule puts zero-detection sites in (b). On the TMT
  # route the both-arms rule puts them all in (a), leaving (b) reachable only via
  # the high-missingness flag - and unreachable at four replicates, where two valid
  # of four is 50% missing. The old Condition Dropout tab was therefore empty by
  # construction on TMT and read as "nothing to see" when everything was next door.
  # Fisher and limma statistics are kept in SEPARATE columns: they are different
  # tests on different populations and must never share a column.
  output$ptm_nmar_table <- DT::renderDT({
    .empty_msg <- function(m) DT::datatable(
      data.frame(Message = m), options = list(dom = "t"), rownames = FALSE)

    .cols <- c("Site", "Gene", "Assessment", "Detection_diff",
               "Fisher_p", "Fisher_adj_p", "limma_log2FC", "limma_adj_p")
    .blank <- function(n) {
      d <- data.frame(Site = character(n), Gene = character(n),
                      Assessment = character(n), Detection_diff = rep(NA_real_, n),
                      Fisher_p = rep(NA_real_, n), Fisher_adj_p = rep(NA_real_, n),
                      limma_log2FC = rep(NA_real_, n), limma_adj_p = rep(NA_real_, n),
                      stringsAsFactors = FALSE)
      d
    }
    parts <- list()

    # (a) excluded from the model - Fisher detection test
    fx <- ptm_rv$nmar_table
    if (!is.null(fx) && nrow(fx) > 0) {
      d <- .blank(nrow(fx))
      d$Site       <- if ("Protein" %in% colnames(fx)) as.character(fx$Protein)
                      else as.character(fx[[1]])
      d$Gene       <- if ("Gene" %in% colnames(fx)) as.character(fx$Gene) else NA_character_
      d$Assessment <- "Excluded from model - Fisher detection test"
      if ("Detection_diff" %in% colnames(fx)) d$Detection_diff <- round(fx$Detection_diff, 3)
      if ("pvalue"         %in% colnames(fx)) d$Fisher_p        <- signif(fx$pvalue, 3)
      if ("adj.pvalue"     %in% colnames(fx)) d$Fisher_adj_p    <- signif(fx$adj.pvalue, 3)
      for (cn in grep("^Found_", colnames(fx), value = TRUE)) d[[cn]] <- as.character(fx[[cn]])
      parts[["excluded"]] <- d
    }

    # (b) modelled, but with a zero-detection or heavily-missing condition
    # (primary table: adjusted, or absolute PTM for a PTM-only result)
    md <- ptm_rv$de_result_adjusted %||% ptm_rv$de_result_unadjusted
    if (!is.null(md) && nrow(md) > 0) {
      fc_cols <- grep("^Found_", colnames(md), value = TRUE)
      zero_hit <- if (length(fc_cols) > 0)
        apply(md[, fc_cols, drop = FALSE], 1, function(r) any(grepl("^0/", r)))
      else rep(FALSE, nrow(md))
      heavy_hit <- if ("high_missingness" %in% colnames(md))
        !is.na(md$high_missingness) & md$high_missingness else rep(FALSE, nrow(md))
      keep <- zero_hit | heavy_hit
      if (any(keep)) {
        m <- md[keep, , drop = FALSE]
        d <- .blank(nrow(m))
        d$Site       <- if ("PTM_Label" %in% colnames(m)) as.character(m$PTM_Label)
                        else if ("PTM_ID" %in% colnames(m)) as.character(m$PTM_ID)
                        else as.character(m[[1]])
        d$Gene       <- if ("Gene" %in% colnames(m)) as.character(m$Gene) else NA_character_
        d$Assessment <- ifelse(zero_hit[keep],
                               "Modelled by limma - condition with zero detections",
                               "Modelled by limma - condition >75% missing")
        if ("log2FC"     %in% colnames(m)) d$limma_log2FC <- round(m$log2FC, 4)
        if ("adj.pvalue" %in% colnames(m)) d$limma_adj_p  <- signif(m$adj.pvalue, 3)
        for (cn in fc_cols) d[[cn]] <- as.character(m[[cn]])
        parts[["modelled"]] <- d
      }
    }

    if (length(parts) == 0)
      return(.empty_msg(paste0(
        "No site in this comparison has a detection concern: none failed the ",
        "valid-value filter, and no modelled site has a condition with zero ",
        "detections. If PTM DE has not been run yet, run it first.")))

    # Union the Found_* columns so rbind succeeds when the two parts differ.
    all_found <- unique(unlist(lapply(parts, function(d)
      grep("^Found_", colnames(d), value = TRUE))))
    parts <- lapply(parts, function(d) {
      for (cn in setdiff(all_found, colnames(d))) d[[cn]] <- NA_character_
      d[, c(.cols, all_found), drop = FALSE]
    })
    out <- do.call(rbind, parts)
    out <- out[order(out$Assessment,
                     ifelse(is.na(out$Fisher_adj_p), out$limma_adj_p, out$Fisher_adj_p)),
               , drop = FALSE]

    .n_ex <- if (!is.null(parts[["excluded"]])) nrow(parts[["excluded"]]) else 0L
    .n_md <- if (!is.null(parts[["modelled"]])) nrow(parts[["modelled"]]) else 0L
    DT::datatable(out,
      options = list(pageLength = 25, scrollX = TRUE, scrollY = "400px", dom = "frtip"),
      filter = "top", rownames = FALSE, class = "compact stripe hover",
      caption = paste0(
        "Detection concerns - ", .n_ex, " excluded from the model (Fisher test), ",
        .n_md, " modelled but with a zero-detection or heavily-missing condition. ",
        "Fisher and limma statistics are separate columns: they are different tests ",
        "on different sites, and no joint false discovery rate is claimed across them."))
  })

  observeEvent(input$ptm_load_rds, {
    req(input$ptm_load_rds)
    tryCatch({
      loaded <- readRDS(input$ptm_load_rds$datapath)

      # Restore all 3 models if available (backward-compatible with old RDS files,
      # whose single de_result was the adjusted table; a PTM-only save has an
      # absolute table and no adjusted one)
      ptm_rv$de_result_adjusted   <- if (!is.null(loaded$de_result_adjusted)) loaded$de_result_adjusted
                                     else if (is.null(loaded$de_result_unadjusted)) loaded$de_result
      ptm_rv$de_result_unadjusted <- loaded$de_result_unadjusted
      ptm_rv$de_result_protein    <- loaded$de_result_protein
      .loaded_primary             <- if (is.null(ptm_rv$de_result_adjusted)) "unadjusted" else "adjusted"
      ptm_rv$de_view_mode         <- loaded$de_view_mode %||% .loaded_primary
      ptm_rv$de_result            <- loaded$de_result
      # Restore or rebuild master table
      if (!is.null(loaded$de_master)) {
        ptm_rv$de_master <- loaded$de_master
      } else {
        ptm_rv$de_master <- .build_master_table(
          ptm_rv$de_result_adjusted,
          ptm_rv$de_result_unadjusted,
          ptm_rv$de_result_protein,
          rv$dia_ptm_master_cache$adj_to_abs_map,
          primary = .loaded_primary
        )
      }
      ptm_rv$fc_scatter_data <- loaded$fc_scatter_data
      # Read from the bundle, never from ptm_rv: this loader does not otherwise
      # reset nmar_table, so falling back to the reactive would carry the PREVIOUS
      # run's presence/absence table into a restored contrast. NULL for a bundle
      # saved before this field existed, which is the correct answer.
      ptm_rv$nmar_table <- loaded$nmar_table
      # Restored beside the table, from the same bundle. Without these two the
      # reactives keep the PREVIOUS run's power flags, and a save made after a
      # load would pair a restored table with another run's flags.
      ptm_rv$nmar_has_power        <- isTRUE(loaded$nmar_has_power)
      ptm_rv$nmar_min_achievable_p <- as.numeric(loaded$nmar_min_achievable_p %||% NA_real_)
      ptm_rv$comparison <- loaded$comparison
      ptm_rv$method_used <- loaded$method_used

      # ── Republish to the shared store so single-site views see a restored run ──
      # The run observer publishes; this loader previously bypassed it entirely,
      # leaving a populated DE tab and an empty PTM Single.  The cache-watching
      # observer in server.R cannot help here - the cache has not changed - so the
      # bundle carries its own fingerprint and is published only when it matches the
      # cache currently loaded.  A bundle saved against a different dataset, or one
      # predating the fingerprint, is not published.
      .fp_now    <- ptm_cache_fingerprint(rv$dia_ptm_master_cache)
      .fp_loaded <- loaded$cache_fingerprint
      .can_publish <- !is.null(loaded$comparison) && nzchar(as.character(loaded$comparison)) &&
        !is.null(.fp_loaded) && !is.na(.fp_loaded) && !is.na(.fp_now) &&
        identical(as.character(.fp_loaded), as.character(.fp_now))
      if (isTRUE(.can_publish)) {
        if (is.null(rv$dia_ptm_de_results)) rv$dia_ptm_de_results <- list()
        .store <- rv$dia_ptm_de_results
        .parts <- strsplit(as.character(loaded$comparison), " vs ", fixed = TRUE)[[1]]
        .store[[as.character(loaded$comparison)]] <- list(
          adjusted   = ptm_rv$de_result_adjusted,
          unadjusted = ptm_rv$de_result_unadjusted,
          protein    = ptm_rv$de_result_protein,
          nmar       = list(
            table            = loaded$nmar_table,
            has_power        = isTRUE(loaded$nmar_has_power),
            min_achievable_p = as.numeric(loaded$nmar_min_achievable_p %||% NA_real_),
            key              = "abs_id"
          ),
          comparison = as.character(loaded$comparison),
          cond_test  = if (length(.parts) >= 1) .parts[1] else NA_character_,
          cond_ref   = if (length(.parts) >= 2) .parts[2] else NA_character_,
          join_key   = "PTM_ID",
          engine     = "dia_ptm_limma_de (restored from .rds)",
          cache_fingerprint = .fp_now,
          run_at     = loaded$timestamp %||% Sys.time()
        )
        rv$dia_ptm_de_results <- .store
        logger(paste0("PTM DE: restored '", loaded$comparison,
                      "' from .rds and published to the shared store."))
      } else {
        logger(paste0("PTM DE: restored bundle was not published to the shared store ",
                      "(no fingerprint, or it does not match the loaded summarization ",
                      "cache). Single-site statistics will be unavailable for this ",
                      "contrast until PTM DE is re-run."))
      }
      ptm_rv$n_files <- if (!is.null(loaded$n_files)) loaded$n_files else ptm_rv$n_files
      ptm_rv$conditions <- loaded$conditions

      # Update view mode radio to match: the views this result has, its saved view selected
      .lay <- if (.loaded_primary == "unadjusted") "unadjusted" else c("adjusted", "unadjusted", "protein")
      updateRadioButtons(session, "ptm_view_mode",
                         choices = setNames(.lay, PTM_LAYER_LABELS[.lay]), selected = ptm_rv$de_view_mode)

      # Restore provenance if present
      if (!is.null(loaded$provenance) && !is.null(loaded$provenance$ptm_de)) {
        rv$provenance$ptm_de <- loaded$provenance$ptm_de
      }

      ptm_rv$status <- paste0("Loaded saved PTM analysis: ", loaded$comparison,
                              " (", loaded$method_used, ", ",
                              nrow(loaded$de_result), " rows)")

      if (!is.null(loaded$conditions) && length(loaded$conditions) >= 2) {
        comp_parts <- strsplit(loaded$comparison, " vs ")[[1]]
        sel_a <- if (length(comp_parts) >= 1 && comp_parts[1] %in% loaded$conditions) comp_parts[1] else loaded$conditions[2]
        sel_b <- if (length(comp_parts) >= 2 && comp_parts[2] %in% loaded$conditions) comp_parts[2] else loaded$conditions[1]
        updateSelectInput(session, "ptm_group_a", choices = loaded$conditions, selected = sel_a)
        updateSelectInput(session, "ptm_group_b", choices = loaded$conditions, selected = sel_b)
      }

      showNotification("PTM analysis loaded.", type = "message")
      logger(paste0("PTM DE: Loaded RDS - ", loaded$comparison))
    }, error = function(e) {
      logger(paste0("!!! [PTM_DE] Failed to load RDS: ", e$message))
      showNotification(paste0("Error loading RDS: ", e$message), type = "error")
    })
  })

  # ===========================================================================
  # Post-analysis export: PhosR & KSEA
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  # Lazy FASTA loader — cached per session, keyed by species.
  # species: "human" (default) or "mouse"
  .fasta_env <- new.env(parent = emptyenv())

  # The UniProt reference proteome (one sequence per gene) is not shipped with
  # OSMSS. On first use it is downloaded from UniProt (CC BY 4.0), about 8 MB,
  # and kept in data/KSEA/; later sessions read the local copy. A plain .fasta
  # placed there by the user is used in preference to the download.
  .ptm_proteome_ids <- list(
    human = c(proteome = "UP000005640", taxon = "9606"),
    mouse = c(proteome = "UP000000589", taxon = "10090"))

  .ptm_fasta_path <- function(species) {
    ids  <- .ptm_proteome_ids[[if (identical(species, "mouse")) "mouse" else "human"]]
    stem <- paste0(ids[["proteome"]], "_", ids[["taxon"]], ".fasta")
    dir  <- file.path(DATA_DIR, "KSEA")
    for (p in file.path(dir, c(stem, paste0(stem, ".gz"))))
      if (file.exists(p) && file.info(p)$size > 0) return(p)
    # Not present yet: download the gzipped FASTA once.
    url  <- paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/",
                   "knowledgebase/reference_proteomes/Eukaryota/", ids[["proteome"]], "/",
                   stem, ".gz")
    dest <- file.path(dir, paste0(stem, ".gz"))
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    nid <- showNotification(paste0("Downloading the UniProt reference proteome (",
                                   ids[["proteome"]], ", about 8 MB). This happens once."),
                            type = "message", duration = NULL)
    on.exit(removeNotification(nid), add = TRUE)
    logger(paste0("[Export] Downloading UniProt reference proteome: ", url))
    tmp <- paste0(dest, ".part")
    old_timeout <- getOption("timeout"); options(timeout = max(300, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    ok <- tryCatch({
      utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
      file.exists(tmp) && file.info(tmp)$size > 1e6
    }, error = function(e) { logger(paste0("[Export] Proteome download failed: ", conditionMessage(e))); FALSE },
       warning = function(w) { logger(paste0("[Export] Proteome download failed: ", conditionMessage(w))); FALSE })
    if (!isTRUE(ok)) {
      if (file.exists(tmp)) unlink(tmp)
      showNotification(paste0("Could not download the UniProt reference proteome. The export needs it ",
                              "for sequence windows. Check the internet connection, or place ", stem,
                              " (from uniprot.org, proteome ", ids[["proteome"]], ") in data/KSEA/."),
                       type = "error", duration = 15)
      return(NULL)
    }
    file.rename(tmp, dest)
    logger(paste0("[Export] Proteome saved: ", dest))
    dest
  }

  .ptm_load_fasta <- function(species = "human") {
    cache_key  <- paste0("proteome_", species)
    if (!is.null(.fasta_env[[cache_key]])) return(.fasta_env[[cache_key]])
    fasta_path <- .ptm_fasta_path(species)
    if (is.null(fasta_path)) {
      logger("[Export] FASTA not available")
      return(NULL)
    }
    tryCatch({
      prot <- Biostrings::readAAStringSet(fasta_path)   # reads .fasta and .fasta.gz
      # Normalise names: sp|ACCESSION|NAME or tr|ACCESSION|NAME -> ACCESSION
      nms <- sub("^(sp|tr)\\|([A-Za-z0-9_-]+)\\|.*", "\\2", names(prot))
      names(prot) <- nms
      .fasta_env[[cache_key]] <- prot
      logger(paste0("[Export] FASTA loaded (", species, "): ", length(prot), " proteins"))
      prot
    }, error = function(e) {
      logger(paste0("[Export] FASTA load error: ", e$message))
      NULL
    })
  }

  # ---------------------------------------------------------------------------
  # .ptm_plex_metadata_row()
  # For TMT exports only: builds a single-row data.frame with PhosR_ID =
  # "__plex__" whose sample columns contain the plex (Mixture) label for each
  # sample.  This row is prepended to the PhosR CSV so server_ks_scoring.R can
  # detect it, extract batch labels, and apply RUVphospho() before scoring.
  # Returns NULL when annotation_map is unavailable or has no Mixture column.
  # ---------------------------------------------------------------------------
  .ptm_plex_metadata_row <- function(out_df) {
    ann <- ptm_rv$annotation_map
    if (is.null(ann) || !"Mixture" %in% colnames(ann))
      return(NULL)

    # Sample columns in out_df (exclude fixed metadata and DE stat columns)
    sample_cols <- setdiff(colnames(out_df),
                           c("PhosR_ID", "Gene", "Protein", "Residue",
                             "Flanking_Sequence_15mer", "log2FC", "adj.pvalue"))
    if (length(sample_cols) == 0L) return(NULL)

    # Build sample_label -> Mixture map.
    # Two possible naming conventions depending on how the matrix was built:
    #   (a) DIA / TMT PTM prerun: columns are raw PD ColumnName values
    #       (e.g. "Abundance: F4: 126, Sample, TPA 0 hr, 1")
    #   (b) MSstatsTMT path (.ptm_tmt_to_matrices): columns are Condition_BioReplicate
    # Try (a) first via ann$ColumnName; fall back to (b).
    plex_map <- NULL
    if ("ColumnName" %in% colnames(ann)) {
      plex_map <- setNames(ann$Mixture, ann$ColumnName)
    }
    # If none of the sample_cols matched (a), try (b)
    if (is.null(plex_map) || !any(sample_cols %in% names(plex_map))) {
      if (all(c("Condition", "BioReplicate") %in% colnames(ann))) {
        plex_map <- setNames(ann$Mixture, paste0(ann$Condition, "_", ann$BioReplicate))
      } else {
        return(NULL)
      }
    }

    meta <- out_df[1L, , drop = FALSE]
    meta[1L, c("PhosR_ID", "Gene", "Protein", "Residue",
               "Flanking_Sequence_15mer")] <- "__plex__"
    if ("log2FC"     %in% colnames(meta)) meta[1L, "log2FC"]     <- NA_real_
    if ("adj.pvalue" %in% colnames(meta)) meta[1L, "adj.pvalue"] <- NA_real_
    for (sc in sample_cols)
      meta[1L, sc] <- if (!is.na(plex_map[sc])) plex_map[sc] else NA_character_
    meta
  }

  # ---------------------------------------------------------------------------
  # .ptm_condition_metadata_row()
  # Same structure as .ptm_plex_metadata_row() but encodes the Condition for
  # each sample column.  Written as PhosR_ID == "__condition__" so the KS
  # scoring module can build the RUVphospho M matrix (biology to preserve).
  # ---------------------------------------------------------------------------
  .ptm_condition_metadata_row <- function(out_df) {
    ann <- ptm_rv$annotation_map
    if (is.null(ann) || !"Condition" %in% colnames(ann))
      return(NULL)

    sample_cols <- setdiff(colnames(out_df),
                           c("PhosR_ID", "Gene", "Protein", "Residue",
                             "Flanking_Sequence_15mer", "log2FC", "adj.pvalue"))
    if (length(sample_cols) == 0L) return(NULL)

    # Build sample_label -> Condition map (same two-path logic as plex helper)
    cond_map <- NULL
    if ("ColumnName" %in% colnames(ann)) {
      cond_map <- setNames(ann$Condition, ann$ColumnName)
    }
    if (is.null(cond_map) || !any(sample_cols %in% names(cond_map))) {
      if (all(c("Condition", "BioReplicate") %in% colnames(ann))) {
        cond_map <- setNames(ann$Condition, paste0(ann$Condition, "_", ann$BioReplicate))
      } else {
        return(NULL)
      }
    }

    meta <- out_df[1L, , drop = FALSE]
    meta[1L, c("PhosR_ID", "Gene", "Protein", "Residue",
               "Flanking_Sequence_15mer")] <- "__condition__"
    if ("log2FC"     %in% colnames(meta)) meta[1L, "log2FC"]     <- NA_real_
    if ("adj.pvalue" %in% colnames(meta)) meta[1L, "adj.pvalue"] <- NA_real_
    for (sc in sample_cols)
      meta[1L, sc] <- if (!is.na(cond_map[sc])) cond_map[sc] else NA_character_
    meta
  }


  # Extract a ±7 aa flanking window (15-mer) for a given protein + position.
  # Returns a 15-character string padded with "_" at sequence edges.
  .ptm_get_flanking_seq <- function(acc, pos_int, proteome, window = 7L) {
    if (is.null(proteome) || !(acc %in% names(proteome))) return("")
    if (is.na(pos_int) || pos_int < 1L) return("")
    prot_seq  <- as.character(proteome[[acc]])
    prot_len  <- nchar(prot_seq)
    if (pos_int > prot_len) return("")
    start     <- max(1L, pos_int - window)
    end       <- min(prot_len, pos_int + window)
    flank     <- substr(prot_seq, start, end)
    left_pad  <- paste(rep("_", window - (pos_int - start)), collapse = "")
    right_pad <- paste(rep("_", window - (end - pos_int)),   collapse = "")
    paste0(left_pad, flank, right_pad)
  }

  # Build a KSEA-format data.frame from a formatted PTM DE result.
  # de_df columns used: Protein (accession), Gene, Site (e.g. "S14" or "S14+S18"),
  #                     log2FC, adj.pvalue
  # FC is converted to LINEAR scale (2^log2FC) as required by KSEA.
  # Multi-site residues ("S14+S18") are kept as one row with ";" separator.
  .ptm_build_ksea_export <- function(de_df, proteome) {
    if (is.null(de_df) || nrow(de_df) == 0L) return(NULL)
    rows <- vector("list", nrow(de_df))
    for (i in seq_len(nrow(de_df))) {
      acc      <- as.character(de_df$Protein[i])
      gene     <- as.character(de_df$Gene[i])
      site_str <- as.character(de_df$Site[i])
      log2fc   <- suppressWarnings(as.numeric(de_df$log2FC[i]))
      adjp     <- suppressWarnings(as.numeric(de_df$adj.pvalue[i]))

      # Convert "+" -> ";" for KSEA Residue.Both (one row per peptide entry)
      residue <- gsub("+", ";", site_str, fixed = TRUE)

      # Flanking sequence from the first site in the peptide
      pep        <- ""
      first_site <- trimws(strsplit(site_str, "+", fixed = TRUE)[[1]])[1]
      pos        <- suppressWarnings(as.integer(gsub("[^0-9]", "", first_site)))
      if (!is.na(pos) && !is.null(proteome)) {
        pep <- .ptm_get_flanking_seq(acc, pos, proteome)
      }

      # Linear FC — KSEA requires linear fold change, NOT log2
      fc_lin <- if (!is.na(log2fc)) 2^log2fc else NA_real_

      rows[[i]] <- data.frame(
        Protein      = acc,
        Gene         = gene,
        Peptide      = pep,
        Residue.Both = residue,
        p            = adjp,
        FC           = fc_lin,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  # Build a PhosR-format data.frame from a per-sample abundance matrix.
  # mat       : rows = site IDs (OSMSS format), cols = sample names, values = log2 intensities
  # site_to_gene : named char vector site_id -> gene_name (from dia_ptm_master_cache)
  # Multi-phospho peptides (e.g. "P60709_S14+S18_Phospho_M2") are EXPANDED into
  # separate single-site rows — both rows carry the same per-sample values since
  # PhosR cannot distinguish individual site contributions within one peptide.
  # Row ID format follows PhosR convention: "GENE;RESIDUE;" (e.g. "KRT16;S44;")
  .ptm_build_phosr_export <- function(mat, site_to_gene, proteome) {
    if (is.null(mat) || nrow(mat) == 0L) return(NULL)
    site_ids <- rownames(mat)
    samples  <- colnames(mat)
    rows     <- list()

    for (i in seq_along(site_ids)) {
      raw_id   <- site_ids[i]
      row_vals <- as.numeric(mat[i, ])

      # Strip SharedInference suffix added by the DIA PTM pipeline for shared peptides
      clean_id <- sub("_SharedInference$", "", raw_id)

      # Resolve gene name: prefer cache map, fall back to first token of ID
      gene <- if (!is.null(site_to_gene) && raw_id %in% names(site_to_gene)) {
        site_to_gene[[raw_id]]
      } else {
        sub("_.*", "", clean_id)
      }
      # site_to_gene may carry every mapped protein, semicolon-joined. PhosR's
      # row-name format is itself semicolon-delimited ("GENE;RESIDUE;"), so a
      # joined label would be parsed as gene + residue + stray token. Take the
      # primary gene, which is the razor winner on TMT and the first listed
      # protein on DIA.
      gene <- strsplit(gene, ";", fixed = TRUE)[[1]][1L]
      # strsplit("", ";") returns character(0), so an empty gene label would make
      # the row name "NA;S35;". Fall back to the accession, as the else-branch
      # above already does.
      if (is.na(gene) || !nzchar(gene)) gene <- sub("_.*", "", clean_id)

      # Parse accession and residue(s) from site ID
      # Expected pattern: ACCESSION_RESIDUES_Phospho... where RESIDUES may contain "+"
      m <- regmatches(clean_id,
                      regexec("^([A-Za-z0-9]+)_([A-Za-z][0-9]+(\\+[A-Za-z][0-9]+)*)_Phospho",
                              clean_id, perl = TRUE))[[1]]
      if (length(m) >= 3L) {
        acc       <- m[2L]
        sites_str <- m[3L]
      } else {
        # Fallback: split on "_"
        parts     <- strsplit(clean_id, "_", fixed = TRUE)[[1]]
        acc       <- if (length(parts) >= 1L) parts[1L] else clean_id
        sites_str <- if (length(parts) >= 2L) parts[2L] else ""
      }

      # Expand multi-site (e.g. "S14+S18" -> c("S14", "S18"))
      individual_sites <- trimws(strsplit(sites_str, "+", fixed = TRUE)[[1]])
      individual_sites <- individual_sites[nzchar(individual_sites)]
      if (length(individual_sites) == 0L) individual_sites <- sites_str

      for (site in individual_sites) {
        if (!nzchar(site)) next

        # PhosR standard row-name format: GENE;RESIDUE;
        phosr_id <- paste0(gene, ";", site, ";")

        # Flanking sequence (15-mer) for the kinaseSubstrateScore() seqs argument
        pos   <- suppressWarnings(as.integer(gsub("[^0-9]", "", site)))
        flank <- if (!is.na(pos) && !is.null(proteome)) {
          .ptm_get_flanking_seq(acc, pos, proteome)
        } else ""

        row_df          <- data.frame(PhosR_ID                = phosr_id,
                                      Gene                    = gene,
                                      Protein                 = acc,
                                      Residue                 = site,
                                      Flanking_Sequence_15mer = flank,
                                      stringsAsFactors        = FALSE)
        row_df[samples] <- as.list(row_vals)
        rows[[length(rows) + 1L]] <- row_df
      }
    }

    if (length(rows) == 0L) return(NULL)
    do.call(rbind, rows)
  }

  # ---------------------------------------------------------------------------
  # UI renderers — show buttons when data is available, placeholder otherwise
  # ---------------------------------------------------------------------------


  output$ptm_phosr_export_ui <- renderUI({
    cache  <- rv$dia_ptm_master_cache
    is_dia <- !is.null(cache) && !is.null(cache$ptm_abs_matrix)

    if (!is_dia) {
      return(tags$p(class = "text-muted small fst-italic",
                    icon("triangle-exclamation"),
                    " Complete the DIA PTM Upload Wizard first."))
    }

    tagList(
      radioButtons(
        "ptm_phosr_species", label = NULL,
        choices  = c("Human" = "human", "Mouse" = "mouse"),
        selected = "human", inline = TRUE),
      # PTM-only caches have no adjusted matrix
      if ("adjusted" %in% ptm_layers(cache))
        downloadButton(
          "ptm_phosr_adj_csv", "Adjusted (Stoichiometry)",
          class = "btn-outline-primary btn-sm w-100 mb-2", icon = icon("flask")),
      downloadButton(
        "ptm_phosr_abs_csv", "Absolute PTM",
        class = "btn-outline-primary btn-sm w-100 mb-2", icon = icon("flask")))
  })

  output$ptm_ksea_export_ui <- renderUI({
    has_adj <- !is.null(ptm_rv$de_result_adjusted)  && nrow(ptm_rv$de_result_adjusted)  > 0L
    has_abs <- !is.null(ptm_rv$de_result_unadjusted) && nrow(ptm_rv$de_result_unadjusted) > 0L
    if (!has_adj && !has_abs) {
      return(tags$p(class = "text-muted small fst-italic",
                    icon("triangle-exclamation"),
                    " Run PTM DE analysis first."))
    }
    btns <- list()
    if (has_adj)
      btns[[length(btns) + 1L]] <- downloadButton(
        "ptm_ksea_adj_csv", "Adjusted (Stoichiometry)",
        class = "btn-outline-success btn-sm w-100 mb-2", icon = icon("chart-bar"))
    if (has_abs)
      btns[[length(btns) + 1L]] <- downloadButton(
        "ptm_ksea_abs_csv", "Absolute PTM",
        class = "btn-outline-success btn-sm w-100 mb-2", icon = icon("chart-bar"))
    do.call(tagList, btns)
  })

  # ---------------------------------------------------------------------------
  # Download handlers
  # ---------------------------------------------------------------------------

  output$ptm_phosr_adj_csv <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_PhosR_Adjusted_", comp_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      proteome <- .ptm_load_fasta(species = isolate(input$ptm_phosr_species) %||% "human")

      # Detect source: DIA or TMT
      if (!is.null(rv$dia_ptm_master_cache) &&
          !is.null(rv$dia_ptm_master_cache$ptm_adj_matrix)) {
        # DIA path (also used by TMT PTM prerun bundles)
        cache  <- rv$dia_ptm_master_cache
        out_df <- .ptm_build_phosr_export(cache$ptm_adj_matrix,
                                          cache$site_to_gene, proteome)
        src_label <- if (isTRUE(rv$data_source == "tmt_ptm")) "TMT adjusted" else "DIA adjusted"
      } else {
        showNotification("No PTM data available for PhosR adjusted export.", type = "warning")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      if (is.null(out_df)) {
        showNotification("No data available for PhosR adjusted export.", type = "warning")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      # Attach DE stats (log2FC, adj.pvalue) so run_ks_scoring.R can use DE-guided
      # site selection. Multi-site entries ("S14+S18") are expanded to individual rows.
      de_for_join <- ptm_rv$de_result_adjusted
      if (!is.null(de_for_join) && nrow(de_for_join) > 0L &&
          all(c("Gene", "Site", "log2FC", "adj.pvalue") %in% colnames(de_for_join))) {
        de_exp <- do.call(rbind, lapply(seq_len(nrow(de_for_join)), function(j) {
          sites <- trimws(strsplit(as.character(de_for_join$Site[j]), "+", fixed = TRUE)[[1]])
          data.frame(
            PhosR_ID   = paste0(de_for_join$Gene[j], ";", sites, ";"),
            log2FC     = suppressWarnings(as.numeric(de_for_join$log2FC[j])),
            adj.pvalue = suppressWarnings(as.numeric(de_for_join$adj.pvalue[j])),
            stringsAsFactors = FALSE
          )
        }))
        m <- match(out_df$PhosR_ID, de_exp$PhosR_ID)
        out_df$log2FC     <- de_exp$log2FC[m]
        out_df$adj.pvalue <- de_exp$adj.pvalue[m]
        logger(paste0("[Export] DE stats joined (adj): ",
                      sum(!is.na(out_df$log2FC)), " / ", nrow(out_df), " sites matched"))
      }

      # For TMT: prepend plex + condition metadata rows so KS scoring module
      # can apply RUVphospho() with a proper M matrix (biology to preserve).
      if (src_label == "TMT adjusted") {
        plex_row <- .ptm_plex_metadata_row(out_df)
        cond_row <- .ptm_condition_metadata_row(out_df)
        if (!is.null(plex_row)) {
          rows_to_prepend <- if (!is.null(cond_row)) rbind(plex_row, cond_row) else plex_row
          out_df <- rbind(rows_to_prepend, out_df)
          logger(paste0("[Export] Plex", if (!is.null(cond_row)) " + condition" else "",
                        " metadata row(s) prepended to PhosR adjusted CSV"))
        }
      }

      write.csv(out_df, file, row.names = FALSE)
      logger(paste0("[Export] PhosR adjusted CSV (", src_label, "): ",
                    nrow(out_df), " rows, ", ncol(out_df) - 7L, " samples"))
    }
  )

  output$ptm_phosr_abs_csv <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_PhosR_Absolute_", comp_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      proteome <- .ptm_load_fasta(species = isolate(input$ptm_phosr_species) %||% "human")

      # Detect source: DIA or TMT
      if (!is.null(rv$dia_ptm_master_cache) &&
          !is.null(rv$dia_ptm_master_cache$ptm_abs_matrix)) {
        # DIA path (also used by TMT PTM prerun bundles)
        cache  <- rv$dia_ptm_master_cache
        out_df <- .ptm_build_phosr_export(cache$ptm_abs_matrix,
                                          cache$site_to_gene, proteome)
        src_label <- if (isTRUE(rv$data_source == "tmt_ptm")) "TMT absolute" else "DIA absolute"
      } else {
        showNotification("No PTM data available for PhosR absolute export.", type = "warning")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      if (is.null(out_df)) {
        showNotification("No data available for PhosR absolute export.", type = "warning")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }

      # Attach DE stats from unadjusted (absolute) results
      de_for_join <- ptm_rv$de_result_unadjusted
      if (!is.null(de_for_join) && nrow(de_for_join) > 0L &&
          all(c("Gene", "Site", "log2FC", "adj.pvalue") %in% colnames(de_for_join))) {
        de_exp <- do.call(rbind, lapply(seq_len(nrow(de_for_join)), function(j) {
          sites <- trimws(strsplit(as.character(de_for_join$Site[j]), "+", fixed = TRUE)[[1]])
          data.frame(
            PhosR_ID   = paste0(de_for_join$Gene[j], ";", sites, ";"),
            log2FC     = suppressWarnings(as.numeric(de_for_join$log2FC[j])),
            adj.pvalue = suppressWarnings(as.numeric(de_for_join$adj.pvalue[j])),
            stringsAsFactors = FALSE
          )
        }))
        m <- match(out_df$PhosR_ID, de_exp$PhosR_ID)
        out_df$log2FC     <- de_exp$log2FC[m]
        out_df$adj.pvalue <- de_exp$adj.pvalue[m]
        logger(paste0("[Export] DE stats joined (abs): ",
                      sum(!is.na(out_df$log2FC)), " / ", nrow(out_df), " sites matched"))
      }

      # For TMT: prepend plex + condition metadata rows so KS scoring module
      # can apply RUVphospho() with a proper M matrix (biology to preserve).
      if (src_label == "TMT absolute") {
        plex_row <- .ptm_plex_metadata_row(out_df)
        cond_row <- .ptm_condition_metadata_row(out_df)
        if (!is.null(plex_row)) {
          rows_to_prepend <- if (!is.null(cond_row)) rbind(plex_row, cond_row) else plex_row
          out_df <- rbind(rows_to_prepend, out_df)
          logger(paste0("[Export] Plex", if (!is.null(cond_row)) " + condition" else "",
                        " metadata row(s) prepended to PhosR absolute CSV"))
        }
      }

      write.csv(out_df, file, row.names = FALSE)
      logger(paste0("[Export] PhosR absolute CSV (", src_label, "): ",
                    nrow(out_df), " rows, ", ncol(out_df) - 7L, " samples"))
    }
  )

  output$ptm_ksea_adj_csv <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_KSEA_Adjusted_", comp_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(ptm_rv$de_result_adjusted)
      proteome <- .ptm_load_fasta()
      out_df   <- .ptm_build_ksea_export(ptm_rv$de_result_adjusted, proteome)
      if (is.null(out_df)) {
        showNotification("No data available for KSEA adjusted export.", type = "warning")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }
      write.csv(out_df, file, row.names = FALSE)
      logger(paste0("[Export] KSEA adjusted CSV: ", nrow(out_df), " rows"))
    }
  )

  output$ptm_ksea_abs_csv <- downloadHandler(
    filename = function() {
      comp_safe <- gsub("[^a-zA-Z0-9_]", "_",
                        if (is.null(ptm_rv$comparison)) "results" else ptm_rv$comparison)
      paste0("PTM_KSEA_Absolute_", comp_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(ptm_rv$de_result_unadjusted)
      proteome <- .ptm_load_fasta()
      out_df   <- .ptm_build_ksea_export(ptm_rv$de_result_unadjusted, proteome)
      if (is.null(out_df)) {
        showNotification("No data available for KSEA absolute export.", type = "warning")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }
      write.csv(out_df, file, row.names = FALSE)
      logger(paste0("[Export] KSEA absolute CSV: ", nrow(out_df), " rows"))
    }
  )
}
