# server_ptm_trends.R
# ==============================================================================
# PTM Trends — fuzzy c-means clustering of stoichiometry-adjusted profiles
#              (DIA and TMT PTM).  No model is fitted and no value is imputed.
#
# Pipeline:
#   1. Sync from PTM DE / PTM Single OR DIA PTM Upload Wizard
#   2. Fuzzy c-means clustering on stoichiometry-adjusted log2 values (absolute
#      PTM values for a PTM-only DIA cache, which has no global proteome)
#   3. Trend visualization with optional gene highlight
#
# Shared globals used:
#   tmt_register_design_grid(), detect_gene_column()
# ==============================================================================
#
# Called from server.R as:
#   server_ptm_trends(input, output, session, rv, logger)
# ==============================================================================

server_ptm_trends <- function(input, output, session, rv, logger) {

  # ── module-local reactive state ──
  ptrnd <- reactiveValues(
    # Raw parsed files
    parsed_files    = list(),
    files           = list(),
    n_files         = 0L,
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
    # Built data tables — stored in rv$ptm_data (shared store); gene_lookup kept local
    gene_lookup_df  = NULL,
    protein_headers = character(0),
    protein_df      = NULL,
    n_files_built   = 0L,
    conditions      = NULL,
    cond_order      = NULL,
    # Analysis outputs
    zscore_df       = NULL,       # long: PTM_ID, Gene, Site, Condition, ZScore, Cluster, Membership
    cluster_means   = NULL,       # long: Cluster, Condition, CondIdx, MeanZ
    cluster_info    = NULL,       # data.frame: Cluster, n_sites
    n_clusters      = 6L,
    membership_thresh = 0.7,
    stringency      = "standard", # "standard" (>=50% valid per condition) or "strict" (100%)
    site_gate       = "none",     # "none" or "de_sig": restrict clustering to DE hits
    site_gate_applied = NULL,     # what the gate actually did on the last run
    stringency_effective = NULL,  # named integer vector: realised threshold per condition
    excluded_df     = NULL,       # sites dropped by the eligibility rule, with reasons
    memberships     = NULL,       # full membership matrix
    render_count    = 0L,
    highlight_gene  = NULL,
    highlight_ptm   = NULL,
    status          = "Upload nested PD CSV file(s) to begin.",
    # For export
    zscore_wide     = NULL,
    avg_wide        = NULL,
    n_omitted       = 0L,
    n_imputed       = 0L,
    lod_value       = NA_real_,
    # Total PTM trends (separate clustering)
    total_zscore_df     = NULL,
    total_cluster_means = NULL,
    total_cluster_info  = NULL,
    total_zscore_wide   = NULL,
    total_avg_wide      = NULL,
    total_memberships   = NULL
  )

  # Stringency control.  Default "standard" when the input is absent, so the module
  # behaves identically whether or not the UI control has been wired.
  observeEvent(input$ptrnd_stringency, {
    val <- input$ptrnd_stringency
    if (!is.null(val) && val %in% c("standard", "strict")) ptrnd$stringency <- val
  }, ignoreInit = TRUE)

  observeEvent(input$ptrnd_site_gate, {
    val <- input$ptrnd_site_gate
    if (!is.null(val) && val %in% c("none", "de_sig")) ptrnd$site_gate <- val
  }, ignoreInit = TRUE)

  # The compute path below is gated on %in% c("dia_ptm", "tmt_ptm"); this selector
  # was not widened with it, so a TMT run fell through to the legacy "ptm_trends"
  # renderer, which reads prov$ptm_trends - a slot the TMT run never writes. The
  # result was no methods text at all, or, in a session where a legacy label-free
  # PTM analysis had run earlier, that run's text rendered against TMT figures.
  # render_dia_ptm_trends_methods() already distinguishes TMT from DIA internally.
  methodsTextServer("ptm_trnd_methods", provenance = rv$provenance,
                    section = function() {
                      if (isTRUE(ptrnd$design_option %in% c("dia_ptm", "tmt_ptm")))
                        "dia_ptm_trends"
                      else "ptm_trends"
                    })

  # ── Local helpers ──────────────────────────────────────────────────────────

  # Schwaemmle-Jensen fuzzifier estimate: the published empirical formula of
  # Schwaemmle and Jensen, Bioinformatics 26(22), 2841-2848, 2010, for a matrix
  # of N rows and D columns:
  #     m = 1 + (1418/N + 22.05) * D^(-2) + (12.33/N + 0.243) * D^(-0.0406*ln(N) - 0.1134)
  # ln is the natural logarithm (R's log()). The estimate depends only on the
  # SHAPE of the matrix, never on its values. It gives the same value as
  # Mfuzz::mestimate(); the formula is written out here so that no package
  # dependency is introduced.
  .sj_fuzzifier <- function(mat) {
    if (is.null(mat)) return(NA_real_)
    N <- nrow(mat); D <- ncol(mat)
    if (is.null(N) || is.null(D) || N < 2L || D < 2L) return(NA_real_)
    1 + (1418/N + 22.05) * D^(-2) +
        (12.33/N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
  }
  .parse_ptm_id <- function(ptm_id) {
    s <- as.character(ptm_id)
    # Match _<Site>_<Modification> from the right (site = letter + digits)
    m <- regmatches(s, regexpr("_[A-Za-z][0-9]+_[^_]+$", s))
    if (length(m) == 1L && nchar(m) > 0L) {
      acc <- sub(m, "", s, fixed = TRUE)
      parts <- strsplit(substring(m, 2L), "_", fixed = TRUE)[[1]]
      list(accession = acc, site = parts[1], modification = paste(parts[-1], collapse = "_"))
    } else {
      # Fallback for IDs without the expected pattern
      parts <- strsplit(s, "_", fixed = TRUE)[[1]]
      if (length(parts) >= 3) {
        list(accession = parts[1],
             site = parts[2],
             modification = paste(parts[3:length(parts)], collapse = "_"))
      } else if (length(parts) == 2) {
        list(accession = parts[1], site = parts[2], modification = "")
      } else {
        list(accession = ptm_id, site = "", modification = "")
      }
    }
  }

  # ===========================================================================
  # SYNC: Inherit from PTM DE (parsed files + design + built tables)
  # ===========================================================================
  observeEvent(rv$ptm_data$source, {
    if (!isTRUE(rv$ptm_data$source %in% c("ptm_de", "ptm_single", "dia_ptm_upload", "tmt_ptm_upload"))) return()
    has_tables <- !is.null(rv$ptm_data$ptm_long) && nrow(rv$ptm_data$ptm_long) > 0
    is_dia_ptm_src <- isTRUE(rv$ptm_data$source %in% c("dia_ptm_upload", "tmt_ptm_upload"))
    if (!has_tables && !is_dia_ptm_src) return()

    # Guard against annotation_map observer wiping cache during sync
    ptrnd$importing <- TRUE
    session$onFlushed(function() { ptrnd$importing <- FALSE }, once = TRUE)

    # Tables are read directly from rv$ptm_data — no local copies needed
    # Always sync gene_lookup_df (even without ptm_long, e.g. cache-only import)
    if (!is.null(rv$ptm_data$gene_lookup))
      ptrnd$gene_lookup_df  <- rv$ptm_data$gene_lookup
    if (!is.null(rv$ptm_data$data_format))
      ptrnd$data_format    <- rv$ptm_data$data_format
    if (!is.null(rv$ptm_data$annotation_map))
      ptrnd$annotation_map <- rv$ptm_data$annotation_map
    if (!is.null(rv$ptm_data$design_option))
      ptrnd$design_option <- rv$ptm_data$design_option
    if (!is.null(rv$ptm_data$selected_ptms))
      ptrnd$selected_ptms <- rv$ptm_data$selected_ptms
    if (!is.null(rv$ptm_data$protein_headers))
      ptrnd$protein_headers <- rv$ptm_data$protein_headers
    if (!is.null(rv$ptm_data$protein_df))
      ptrnd$protein_df      <- rv$ptm_data$protein_df
    if (!is.null(rv$ptm_data$n_files)) {
      ptrnd$n_files         <- rv$ptm_data$n_files
      ptrnd$n_files_built   <- rv$ptm_data$n_files
    }

    all_conditions <- sort(unique(as.character(rv$ptm_data$conditions)))
    ptrnd$conditions <- all_conditions
    ptrnd$cond_order <- if (!is.null(rv$ptm_data$cond_order) &&
                            setequal(rv$ptm_data$cond_order, all_conditions))
                          rv$ptm_data$cond_order else all_conditions

    # Extract design grid fields for display
    ann <- rv$ptm_data$annotation_map
    if (!is.null(ann)) {
      ptrnd$all_channels <- unique(ann$ColumnName)
      ptrnd$n_plexes     <- length(unique(ann$Mixture))
      ptrnd$has_bridge   <- any(ann$Condition == "Norm")
      ptrnd$design_conds <- sort(unique(ann$Condition[ann$Condition != "Norm"]))
    }

    src <- rv$ptm_data$source
    if (has_tables) {
      n_ptm <- length(unique(rv$ptm_data$ptm_long$ProteinName))
      ptrnd$status <- paste0("Synced from ", src, ": ", n_ptm,
                             " PTM sites, ", length(all_conditions), " conditions.")
      logger(paste0("PTM Trends: Synced from ", src, " — ", n_ptm, " PTM sites, ",
                    "conditions: ", paste(all_conditions, collapse = ", ")))
    } else if (is_dia_ptm_src) {
      cache <- rv$dia_ptm_master_cache
      n_sites <- if (!is.null(cache)) cache$n_sites_output else "?"
      ptrnd$status <- paste0("PTM data synced: ", n_sites,
                             " PTM sites across ",
                             length(all_conditions), " conditions. ",
                             "Click Run Trends.")
      logger(paste0("PTM Trends: DIA PTM synced — ", n_sites, " sites, ",
                    "conditions: ", paste(all_conditions, collapse = ", ")))
      .enable_run_btn()
    } else {
      ptrnd$status <- paste0("Design synced from ", src, " (summarization cache imported). ",
                             length(all_conditions), " conditions.")
      logger(paste0("PTM Trends: Design synced from ", src, " (import mode), ",
                    "conditions: ", paste(all_conditions, collapse = ", ")))
    }
  }, ignoreInit = TRUE)

  observeEvent(ptrnd$annotation_map, {
    if (isTRUE(ptrnd$importing)) return()

  }, ignoreInit = TRUE)

  # ===========================================================================
  # CONDITION ORDER MANAGEMENT
  # ===========================================================================
  output$ptm_trends_cond_label <- renderUI({
    conds <- ptrnd$cond_order
    if (is.null(conds) || length(conds) == 0) {
      tags$p(class = "text-muted small", "Upload files to detect conditions.")
    } else {
      tags$div(class = "small",
        tags$span(class = "badge bg-info", paste0(length(conds), " conditions")),
        tags$br(),
        tags$span(class = "text-muted", paste(conds, collapse = " \u2192 "))
      )
    }
  })

  observeEvent(input$ptm_trends_reorder_btn, {
    req(ptrnd$conditions)
    conds <- ptrnd$cond_order
    if (is.null(conds)) conds <- ptrnd$conditions

    showModal(modalDialog(
      title = tags$span(icon("arrows-up-down"), " Reorder Conditions"),
      tagList(
        helpText(class = "small text-muted",
          "Drag conditions to set the X-axis order (left to right). ",
          "Typically chronological or dose-increasing."),
        sortable::rank_list(
          text       = "Condition Order (drag to reorder):",
          labels     = conds,
          input_id   = "ptm_trends_cond_order_sorted",
          options    = sortable::sortable_options(multiDrag = FALSE),
          class      = "default-sortable"
        )
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ptm_trends_cond_order_save", "Apply Order",
                     class = "btn-success", icon = icon("check"))
      )
    ))
  })

  observeEvent(input$ptm_trends_cond_order_save, {
    new_order <- input$ptm_trends_cond_order_sorted
    if (is.null(new_order) || length(new_order) == 0) {
      showNotification("No conditions to order.", type = "warning")
      return()
    }
    ptrnd$cond_order       <- new_order
    rv$ptm_data$cond_order <- new_order
    removeModal()
    logger(paste0("PTM Trends: Condition order set to: ",
                  paste(new_order, collapse = " -> ")))
  })

  # ===========================================================================
  # Helper: enable/disable Run button with colour swap
  # ===========================================================================
  .enable_run_btn <- function() {
    shinyjs::enable("ptm_trends_run_btn")
    shinyjs::removeClass("ptm_trends_run_btn", "btn-secondary")
    shinyjs::addClass("ptm_trends_run_btn", "btn-success")
  }

  .disable_run_btn <- function() {
    shinyjs::disable("ptm_trends_run_btn")
    shinyjs::removeClass("ptm_trends_run_btn", "btn-success")
    shinyjs::addClass("ptm_trends_run_btn", "btn-secondary")
  }

  # Sync Run button state when shared cache changes (e.g. from PTM DE)
  # Also enable when DIA PTM cache is populated
  observeEvent(rv$dia_ptm_master_cache, {
    if (!is.null(rv$dia_ptm_master_cache)) {
      .enable_run_btn()
    }
  }, ignoreInit = TRUE)



  # ===========================================================================
  # RUN ANALYSIS (trend clustering — requires summarization first)
  # ===========================================================================
  observeEvent(input$ptm_trends_run_btn, {
    # ── DIA PTM BRANCH (also handles tmt_ptm which uses the same cache path) ──
    is_dia_ptm <- isTRUE(ptrnd$design_option %in% c("dia_ptm", "tmt_ptm"))
    if (is_dia_ptm) {
      cache <- rv$dia_ptm_master_cache
      if (is.null(cache)) {
        showNotification("PTM summarization cache not found. Complete the PTM Upload Wizard.", type = "warning")
        return()
      }

      cond_order <- ptrnd$cond_order
      if (is.null(cond_order)) cond_order <- cache$design_df$condition |> unique() |> sort()

      if (length(cond_order) < 3) {
        showNotification("Trend analysis requires at least 3 conditions.", type = "warning")
        return()
      }
      if (!requireNamespace("e1071", quietly = TRUE)) {
        showNotification("e1071 package is required.", type = "error", duration = 10)
        return()
      }

      k <- as.integer(if (!is.null(input$ptm_trends_n_clusters)) input$ptm_trends_n_clusters else 6L)
      k <- max(2L, min(20L, k))
      ptrnd$n_clusters        <- k
      ptrnd$membership_thresh <- if (!is.null(input$ptm_trends_membership)) input$ptm_trends_membership else 0.7
      # The fuzzifier is decided below, once the matrix that will actually be
      # clustered is known, because the Schwaemmle-Jensen estimate is a function of
      # that matrix's shape. Estimation is the default; the manual slider applies
      # only when the user unticks the estimate box. A NULL checkbox (control not
      # yet rendered) selects estimation, matching the checkbox's own default, so
      # the two can never disagree.
      fuzziness_input         <- if (!is.null(input$ptrnd_fuzziness))
                                   suppressWarnings(as.numeric(input$ptrnd_fuzziness))
                                 else NA_real_
      .m_auto                 <- is.null(input$ptrnd_fuzz_auto) ||
                                 isTRUE(input$ptrnd_fuzz_auto)

      logger(paste0("PTM Trends [DIA PTM]: Running k=", k, ", conds=",
                    paste(cond_order, collapse = " -> ")))
      ptrnd$highlight_gene <- NULL
      ptrnd$highlight_ptm  <- NULL

      withProgress(message = "Running DIA PTM Trend Analysis...", value = 0.1, {
        # Cluster the cache's primary layer: the stoichiometry-adjusted matrix, or the
        # absolute PTM matrix for a PTM-only cache (ptm_primary_layer()).  In PTM-only
        # mode the row keys are absolute site IDs, so adj_to_abs is the identity, and
        # the separate Total PTM clustering below is skipped as it would repeat this one.
        .primary   <- ptm_primary_layer(cache)
        adj_mat    <- if (.primary == "adjusted") cache$ptm_adj_matrix else cache$ptm_abs_matrix
        design_df  <- cache$design_df
        site_gene  <- cache$site_to_gene
        adj_to_abs <- cache$adj_to_abs_map %||% setNames(rownames(adj_mat), rownames(adj_mat))

        # Keep only samples in cond_order
        use_samples <- design_df$sample_name[design_df$condition %in% cond_order]
        use_samples <- intersect(use_samples, colnames(adj_mat))
        adj_mat     <- adj_mat[, use_samples, drop = FALSE]

        # ── Site eligibility: per-replicate valid-value rule, NO imputation ──
        # Clustering needs a complete condition profile, so the rule is applied to
        # EVERY condition.  This is deliberately stricter than the DE filter, which
        # requires only one of the two contrasted groups to pass, because limma can
        # fit on available observations and cmeans cannot.
        #   standard: min(max(2, ceiling(0.5 * n)), n) valid values per condition
        #   strict  : n valid values per condition (complete cases)
        # The inner min() is a degeneracy guard: at n = 1 the floor of 2 is
        # unsatisfiable and would exclude every site.
        .strict_mode <- identical(ptrnd$stringency, "strict")
        .cond_samps_list <- lapply(cond_order, function(cond) {
          intersect(design_df$sample_name[design_df$condition == cond], use_samples)
        })
        names(.cond_samps_list) <- cond_order
        .min_valid_eff <- vapply(cond_order, function(cond) {
          n_s <- length(.cond_samps_list[[cond]])
          if (n_s == 0L) return(NA_integer_)
          if (.strict_mode) as.integer(n_s)
          else min(max(2L, as.integer(ceiling(0.5 * n_s))), as.integer(n_s))
        }, integer(1))
        names(.min_valid_eff) <- cond_order

        .degenerate <- cond_order[!is.na(.min_valid_eff) &
                                    vapply(cond_order, function(cond)
                                      length(.cond_samps_list[[cond]]) < 2L, logical(1))]
        if (length(.degenerate) > 0)
          logger(paste0("PTM Trends [DIA PTM]: WARNING: condition(s) ",
                        paste(.degenerate, collapse = ", "),
                        " have fewer than 2 samples, so the valid-value threshold ",
                        "degenerates to 1 for them."))

        .empty_conds <- cond_order[is.na(.min_valid_eff)]
        if (length(.empty_conds) > 0) {
          showNotification(paste0("Condition(s) with no samples in the cache: ",
                                  paste(.empty_conds, collapse = ", "),
                                  ". Adjust the condition order or re-run the wizard."),
                           type = "error", duration = 10)
          return()
        }

        .cond_valid <- vapply(cond_order, function(cond) {
          samps <- .cond_samps_list[[cond]]
          rowSums(!is.na(adj_mat[, samps, drop = FALSE])) >= .min_valid_eff[[cond]]
        }, logical(nrow(adj_mat)))
        keep_rows <- if (is.matrix(.cond_valid)) apply(.cond_valid, 1, all)
                     else all(.cond_valid)
        names(keep_rows) <- rownames(adj_mat)

        # ── Excluded-site ledger: sparse sites leave the figure, not the app ──
        .excl_ids <- rownames(adj_mat)[!keep_rows]
        if (length(.excl_ids) > 0) {
          .excl_detected <- vapply(cond_order, function(cond) {
            rowSums(!is.na(adj_mat[.excl_ids, .cond_samps_list[[cond]], drop = FALSE]))
          }, numeric(length(.excl_ids)))
          if (!is.matrix(.excl_detected))
            .excl_detected <- matrix(.excl_detected, nrow = length(.excl_ids),
                                     dimnames = list(.excl_ids, cond_order))
          .excl_df <- data.frame(PTM_ID = .excl_ids, stringsAsFactors = FALSE)
          for (cond in cond_order)
            .excl_df[[paste0("Detected_", cond)]] <-
              paste0(.excl_detected[, cond], "/", length(.cond_samps_list[[cond]]))
          .excl_df$Reason <- paste0(
            "below ", if (.strict_mode) "strict" else "standard",
            " valid-value threshold in >= 1 condition")
          ptrnd$excluded_df <- .excl_df
        } else {
          ptrnd$excluded_df <- data.frame()
        }
        ptrnd$n_omitted   <- length(.excl_ids)
        ptrnd$n_imputed   <- 0L
        ptrnd$lod_value   <- NA_real_
        ptrnd$stringency_effective <- .min_valid_eff
        logger(paste0("PTM Trends [DIA PTM]: eligibility (",
                      if (.strict_mode) "strict" else "standard", ") = ",
                      paste0(names(.min_valid_eff), ":", as.integer(.min_valid_eff),
                             collapse = ", "),
                      "; ", sum(keep_rows), " of ", nrow(adj_mat),
                      " sites retained, ", length(.excl_ids), " excluded."))

        adj_mat <- adj_mat[keep_rows, , drop = FALSE]

        # ── Optional gate: cluster only sites with a DE hit ────────────────
        # Row z-scoring forces every profile to unit amplitude, so without a gate a
        # site moving 0.05 log2 across the experiment is drawn with the same shape
        # and weight as one moving 3.  Off by default; when on, a site qualifies if
        # it reaches adj.p < 0.05 in the table being clustered (ADJUSTED, or ABSOLUTE
        # PTM for a PTM-only cache) of any stored contrast.
        # Applied AFTER the eligibility filter so the excluded-site ledger keeps
        # meaning "failed the valid-value rule" and is not polluted with
        # "not significant".
        ptrnd$site_gate_applied <- NULL
        if (identical(ptrnd$site_gate, "de_sig")) {
          de_store <- rv$dia_ptm_de_results
          if (is.null(de_store) || length(de_store) == 0L) {
            showNotification(paste0("Site gate is set to 'DE-significant only' but no PTM DE ",
                                    "result is stored. Clustering all eligible sites instead."),
                             type = "warning", duration = 8)
            ptrnd$site_gate_applied <- "requested de_sig; no DE result stored; gate not applied"
            logger("PTM Trends [DIA PTM]: de_sig gate requested but no DE result stored - not applied.")
          } else {
            sig_ids <- unique(unlist(lapply(de_store, function(res) {
              tb <- res[[.primary]]
              if (is.null(tb) || !all(c("PTM_ID", "adj.pvalue") %in% colnames(tb)))
                return(character(0))
              as.character(tb$PTM_ID[!is.na(tb$adj.pvalue) & tb$adj.pvalue < 0.05])
            })))
            n_before <- nrow(adj_mat)
            hit      <- rownames(adj_mat) %in% sig_ids
            if (sum(hit) < k) {
              showNotification(paste0("Only ", sum(hit), " DE-significant site(s) - fewer than ",
                                      "k = ", k, ". Clustering all eligible sites instead."),
                               type = "warning", duration = 8)
              ptrnd$site_gate_applied <- paste0("requested de_sig; only ", sum(hit),
                                                " hits (< k); gate not applied")
              logger(paste0("PTM Trends [DIA PTM]: de_sig gate found ", sum(hit),
                            " hits, fewer than k = ", k, " - not applied."))
            } else {
              adj_mat <- adj_mat[hit, , drop = FALSE]
              ptrnd$site_gate_applied <- paste0("de_sig across ", length(de_store),
                                                " contrast(s): ", n_before, " -> ", nrow(adj_mat),
                                                " sites")
              logger(paste0("PTM Trends [DIA PTM]: de_sig gate ", n_before, " -> ",
                            nrow(adj_mat), " sites (adj.p < 0.05 in >= 1 stored contrast)."))
            }
          }
        }

        # Per-condition mean (for trend profile)
        cond_means_list <- lapply(cond_order, function(cond) {
          cond_samps <- .cond_samps_list[[cond]]
          apply(adj_mat[, cond_samps, drop = FALSE], 1, mean, na.rm = TRUE)
        })
        mean_mat  <- do.call(cbind, cond_means_list)
        colnames(mean_mat) <- cond_order
        rownames(mean_mat) <- rownames(adj_mat)
        # Drop rows with no between-condition variation, matching the total-PTM
        # path below. Such a row z-scores to all zeros (the sd == 0 branch of the
        # z-score step) and would otherwise sit at the origin and anchor a cluster
        # centroid while carrying no trend information.
        .adj_rv <- apply(mean_mat, 1, var)
        .adj_rv[is.na(.adj_rv)] <- 0
        n_zero_var_adj <- sum(.adj_rv <= 0)
        if (n_zero_var_adj > 0L) {
          mean_mat <- mean_mat[.adj_rv > 0, , drop = FALSE]
          logger(paste0("PTM Trends [DIA PTM]: ", n_zero_var_adj,
                        " site(s) had zero between-condition variance and were ",
                        "excluded from clustering."))
        }
        n_sites_cl <- nrow(mean_mat)

        if (n_sites_cl < k) {
          showNotification(paste0("Too few sites (", n_sites_cl, ") for ", k, " clusters."),
                           type = "error", duration = 10)
          return()
        }

        # Z-score rows (always done for trend clustering).  Every retained row is
        # complete by construction, so no NA can survive to the clustering step and
        # nothing is filled in.
        mean_mat_z <- t(apply(mean_mat, 1, function(r) {
          s <- sd(r)
          if (is.na(s) || s == 0) return(rep(0, length(r)))
          (r - mean(r)) / s
        }))
        colnames(mean_mat_z) <- cond_order
        rownames(mean_mat_z) <- rownames(mean_mat)

        .m_est_adj  <- .sj_fuzzifier(mean_mat_z)
        .m_manual   <- !.m_auto && !is.na(fuzziness_input)
        fuzziness_m <- if (.m_manual) fuzziness_input
                       else if (is.finite(.m_est_adj) && .m_est_adj > 1) .m_est_adj
                       else 2.0
        fuzzifier_source <- if (.m_manual) "user-specified"
                            else if (is.finite(.m_est_adj) && .m_est_adj > 1)
                              "estimated (Schwaemmle-Jensen)"
                            else "fallback 2.0 (estimate unavailable)"
        logger(paste0("PTM Trends [DIA PTM]: fuzzifier m = ", signif(fuzziness_m, 4),
                      " (", fuzzifier_source, "; estimate = ",
                      if (is.finite(.m_est_adj)) signif(.m_est_adj, 4) else "NA", ")."))

        incProgress(0.4, detail = "Fuzzy c-means clustering...")
        # cmeans initialises by sampling k rows at random.  Seeded so that a re-run
        # on the same cache reproduces the same partition; matches the seed used by
        # the DE route's impute.knn call.
        set.seed(42L)
        cm_res <- tryCatch(
          e1071::cmeans(mean_mat_z, centers = k, m = fuzziness_m,
                        iter.max = 200, method = "cmeans"),
          error = function(e) NULL
        )
        if (is.null(cm_res)) {
          showNotification("Fuzzy c-means failed.", type = "error"); return()
        }

        # Build gene annotation
        gene_vec <- vapply(rownames(mean_mat_z), function(adj_id) {
          abs_id <- adj_to_abs[adj_id]
          if (!is.na(abs_id) && abs_id %in% names(site_gene)) site_gene[[abs_id]]
          # Modification-agnostic fallback: the previous expression was anchored on
          # "_Phospho" and returned the key unchanged for any other modification.
          else .parse_ptm_id(adj_id)$accession
        }, character(1))
        # site_to_gene may carry every mapped protein, semicolon-joined. Display the
        # primary gene, as the DE and Single modules do.
        gene_vec <- vapply(strsplit(gene_vec, ";", fixed = TRUE),
                           `[`, character(1), 1L)

        cluster_assign <- cm_res$cluster
        membership_mat <- cm_res$membership
        membership_max <- apply(membership_mat, 1, max)

        trend_result <- data.frame(
          PTM_ID          = rownames(mean_mat_z),
          Gene            = gene_vec,
          Cluster         = cluster_assign,
          Max_Membership  = membership_max,
          stringsAsFactors = FALSE
        )
        trend_result <- cbind(trend_result, as.data.frame(mean_mat_z))

        ptrnd$trend_result    <- trend_result
        ptrnd$trend_mat_z     <- mean_mat_z
        ptrnd$cluster_assign  <- cluster_assign
        ptrnd$membership_mat  <- membership_mat
        ptrnd$cond_order      <- cond_order
        ptrnd$n_sites_clustered <- n_sites_cl

        # ── Build long_df / cluster_means in the same format the renders expect ──
        {
          adj_ids    <- rownames(mean_mat_z)
          n_sites_d  <- length(adj_ids)
          n_conds_d  <- length(cond_order)

          parsed_d    <- lapply(adj_ids, function(aid) {
            abs_id <- adj_to_abs[aid]
            if (!is.na(abs_id) && abs_id != "") .parse_ptm_id(abs_id)
            else .parse_ptm_id(aid)
          })
          accessions_d <- vapply(parsed_d, function(x) x$accession,   character(1))
          sites_d      <- vapply(parsed_d, function(x) x$site,        character(1))
          mods_d       <- vapply(parsed_d, function(x) x$modification, character(1))
          ptm_labels_d <- paste0(gene_vec, "_", sites_d)

          long_df <- data.frame(
            PTM_ID       = rep(adj_ids,      each  = n_conds_d),
            PTM_Label    = rep(ptm_labels_d, each  = n_conds_d),
            Gene         = rep(gene_vec,     each  = n_conds_d),
            Accession    = rep(accessions_d, each  = n_conds_d),
            Site         = rep(sites_d,      each  = n_conds_d),
            Modification = rep(mods_d,       each  = n_conds_d),
            Condition    = rep(cond_order,   times = n_sites_d),
            CondIdx      = rep(seq_len(n_conds_d), times = n_sites_d),
            ZScore       = as.vector(t(mean_mat_z)),
            Cluster      = rep(cluster_assign, each = n_conds_d),
            Membership   = rep(membership_max, each = n_conds_d),
            stringsAsFactors = FALSE
          )

          cluster_means_d <- aggregate(ZScore ~ Cluster + Condition + CondIdx,
                                       data = long_df, FUN = mean)
          colnames(cluster_means_d)[colnames(cluster_means_d) == "ZScore"] <- "MeanZ"

          cluster_info_d <- as.data.frame(table(Cluster = cluster_assign))
          colnames(cluster_info_d) <- c("Cluster", "n_sites")
          cluster_info_d$Cluster <- as.integer(as.character(cluster_info_d$Cluster))

          # Sort clusters by mean at last condition
          last_row <- cluster_means_d[cluster_means_d$CondIdx == n_conds_d, ]
          cl_order <- last_row$Cluster[order(-last_row$MeanZ)]
          cl_map   <- setNames(seq_along(cl_order), cl_order)

          long_df$Cluster         <- cl_map[as.character(long_df$Cluster)]
          cluster_means_d$Cluster <- cl_map[as.character(cluster_means_d$Cluster)]
          cluster_info_d$Cluster  <- cl_map[as.character(cluster_info_d$Cluster)]

          cl_labels <- paste0(
            "Cluster ", cluster_info_d$Cluster[order(cluster_info_d$Cluster)],
            " (n=", cluster_info_d$n_sites[order(cluster_info_d$Cluster)], ")")
          long_df$ClusterLabel <- factor(
            paste0("Cluster ", long_df$Cluster,
                   " (n=", cluster_info_d$n_sites[
                     match(long_df$Cluster, cluster_info_d$Cluster)], ")"),
            levels = cl_labels)
          cluster_means_d$ClusterLabel <- factor(
            paste0("Cluster ", cluster_means_d$Cluster,
                   " (n=", cluster_info_d$n_sites[
                     match(cluster_means_d$Cluster, cluster_info_d$Cluster)], ")"),
            levels = cl_labels)

          long_df$Condition         <- factor(long_df$Condition,         levels = cond_order)
          cluster_means_d$Condition <- factor(cluster_means_d$Condition, levels = cond_order)

          zscore_wide_d <- as.data.frame(mean_mat_z)
          colnames(zscore_wide_d) <- paste0("ZScore_", cond_order)
          zscore_wide_d$PTM_ID       <- adj_ids
          zscore_wide_d$PTM_Label    <- ptm_labels_d
          zscore_wide_d$Gene         <- gene_vec
          zscore_wide_d$Accession    <- accessions_d
          zscore_wide_d$Site         <- sites_d
          zscore_wide_d$Modification <- mods_d
          zscore_wide_d$Cluster      <- cl_map[as.character(cluster_assign)]
          zscore_wide_d$MaxMembership <- membership_max

          ptrnd$zscore_df           <- long_df
          ptrnd$cluster_means       <- cluster_means_d
          ptrnd$cluster_info        <- cluster_info_d
          avg_wide_d <- as.data.frame(mean_mat)
          colnames(avg_wide_d) <- paste0("MeanLog2_", cond_order)

          ptrnd$zscore_wide         <- zscore_wide_d
          ptrnd$avg_wide            <- avg_wide_d
          ptrnd$memberships         <- membership_mat

          # ── Total PTM trends: cluster on ptm_abs_matrix (unadjusted absolute abundance) ──
          # Not in PTM-only mode, where the clustering above already used this matrix.
          ptrnd$total_zscore_df     <- NULL
          ptrnd$total_cluster_means <- NULL
          tryCatch({
            abs_mat <- cache$ptm_abs_matrix
            if (.primary == "adjusted" && !is.null(abs_mat) && nrow(abs_mat) > 0) {
              t_use_samps <- intersect(use_samples, colnames(abs_mat))
              abs_mat_f   <- abs_mat[, t_use_samps, drop = FALSE]

              # Per-condition means
              t_means_list <- lapply(cond_order, function(cond) {
                cs <- design_df$sample_name[design_df$condition == cond]
                cs <- intersect(cs, t_use_samps)
                if (length(cs) == 0) return(rep(NA_real_, nrow(abs_mat_f)))
                apply(abs_mat_f[, cs, drop = FALSE], 1, mean, na.rm = TRUE)
              })
              t_mean_mat <- do.call(cbind, t_means_list)
              colnames(t_mean_mat) <- cond_order
              rownames(t_mean_mat) <- rownames(abs_mat_f)

              # Same per-replicate eligibility rule as the adjusted path, evaluated
              # on the absolute matrix.  No imputation.
              t_cond_samps <- lapply(cond_order, function(cond)
                intersect(design_df$sample_name[design_df$condition == cond], t_use_samps))
              names(t_cond_samps) <- cond_order
              t_min_valid <- vapply(cond_order, function(cond) {
                n_s <- length(t_cond_samps[[cond]])
                if (n_s == 0L) return(NA_integer_)
                if (.strict_mode) as.integer(n_s)
                else min(max(2L, as.integer(ceiling(0.5 * n_s))), as.integer(n_s))
              }, integer(1))
              names(t_min_valid) <- cond_order
              t_cv <- vapply(cond_order, function(cond) {
                samps <- t_cond_samps[[cond]]
                if (length(samps) == 0L) return(rep(FALSE, nrow(abs_mat_f)))
                rowSums(!is.na(abs_mat_f[, samps, drop = FALSE])) >= t_min_valid[[cond]]
              }, logical(nrow(abs_mat_f)))
              t_keep <- if (is.matrix(t_cv)) apply(t_cv, 1, all) else all(t_cv)
              t_mean_mat <- t_mean_mat[t_keep, , drop = FALSE]
              logger(paste0("PTM Trends [DIA PTM Total]: ", sum(t_keep), " of ",
                            length(t_keep), " absolute sites passed eligibility."))

              if (nrow(t_mean_mat) >= k) {
                # Variance filter: drop rows with no between-condition variation at
                # all.  Retained because a zero-variance row z-scores to all zeros
                # and would otherwise anchor a cluster centroid.
                t_rv <- apply(t_mean_mat, 1, var)
                t_rv[is.na(t_rv)] <- 0
                t_mean_mat <- t_mean_mat[t_rv > 0, , drop = FALSE]

                # Z-score.  Rows are complete by construction after the eligibility
                # filter above; nothing is filled in.
                t_zscore_mat <- t(apply(t_mean_mat, 1, function(r) {
                  s <- sd(r)
                  if (is.na(s) || s == 0) return(rep(0, length(r)))
                  (r - mean(r)) / s
                }))
                colnames(t_zscore_mat) <- cond_order
                rownames(t_zscore_mat) <- rownames(t_mean_mat)
                t_abs_ids    <- rownames(t_zscore_mat)
                t_n_sites    <- nrow(t_zscore_mat)
                t_k          <- max(2L, min(k, t_n_sites - 1L))
                if (!identical(as.integer(t_k), as.integer(k)))
                  logger(paste0("PTM Trends [DIA PTM Total]: k reduced from ", k,
                                " to ", t_k, " (", t_n_sites,
                                " sites available after eligibility)."))

                if (t_n_sites >= 2) {
                  # Own estimate: the total matrix has a different row count, and
                  # the estimate is a function of matrix shape.
                  .m_est_tot <- .sj_fuzzifier(t_zscore_mat)
                  t_fuzz     <- if (!.m_auto && !is.na(fuzziness_input)) fuzziness_input
                                else if (is.finite(.m_est_tot) && .m_est_tot > 1) .m_est_tot
                                else 2.0
                  logger(paste0("PTM Trends [DIA PTM Total]: fuzzifier m = ",
                                signif(t_fuzz, 4), "."))
                  set.seed(42L)
                  t_cl_res   <- e1071::cmeans(t_zscore_mat, centers = t_k, m = t_fuzz,
                                              method = "cmeans", iter.max = 200)
                  t_clusters <- t_cl_res$cluster
                  t_mems     <- t_cl_res$membership
                  t_mem_max  <- apply(t_mems, 1, max)

                  # Gene annotation via site_gene (abs_ids are like P12345_S55_Phospho)
                  t_parsed <- lapply(t_abs_ids, .parse_ptm_id)
                  t_accs   <- vapply(t_parsed, function(x) x$accession,    character(1))
                  t_sites  <- vapply(t_parsed, function(x) x$site,         character(1))
                  t_mods   <- vapply(t_parsed, function(x) x$modification, character(1))
                  t_genes  <- vapply(t_abs_ids, function(aid) {
                    g <- site_gene[aid]
                    if (!is.na(g) && nzchar(g)) g
                    # Modification-agnostic fallback, as in the adjusted path above.
                    else .parse_ptm_id(aid)$accession
                  }, character(1))
                  # site_to_gene may carry every mapped protein, semicolon-joined.
                  # Display the primary gene, as the DE and Single modules do.
                  t_genes  <- vapply(strsplit(t_genes, ";", fixed = TRUE),
                                     `[`, character(1), 1L)
                  t_lbls   <- paste0(t_genes, "_", t_sites)
                  t_n_conds <- length(cond_order)

                  t_long_df <- data.frame(
                    PTM_ID       = rep(t_abs_ids,  each  = t_n_conds),
                    PTM_Label    = rep(t_lbls,     each  = t_n_conds),
                    Gene         = rep(t_genes,    each  = t_n_conds),
                    Accession    = rep(t_accs,     each  = t_n_conds),
                    Site         = rep(t_sites,    each  = t_n_conds),
                    Modification = rep(t_mods,     each  = t_n_conds),
                    Condition    = rep(cond_order, times = t_n_sites),
                    CondIdx      = rep(seq_len(t_n_conds), times = t_n_sites),
                    ZScore       = as.vector(t(t_zscore_mat)),
                    Cluster      = rep(t_clusters, each  = t_n_conds),
                    Membership   = rep(t_mem_max,  each  = t_n_conds),
                    stringsAsFactors = FALSE
                  )

                  t_cmeans <- aggregate(ZScore ~ Cluster + Condition + CondIdx,
                                        data = t_long_df, FUN = mean)
                  colnames(t_cmeans)[colnames(t_cmeans) == "ZScore"] <- "MeanZ"

                  t_cinfo <- as.data.frame(table(Cluster = t_clusters))
                  colnames(t_cinfo) <- c("Cluster", "n_sites")
                  t_cinfo$Cluster   <- as.integer(as.character(t_cinfo$Cluster))

                  # Sort clusters by mean at last condition
                  t_last  <- t_cmeans[t_cmeans$CondIdx == t_n_conds, ]
                  t_ord   <- t_last$Cluster[order(-t_last$MeanZ)]
                  t_cmap  <- setNames(seq_along(t_ord), t_ord)

                  t_long_df$Cluster <- t_cmap[as.character(t_long_df$Cluster)]
                  t_cmeans$Cluster  <- t_cmap[as.character(t_cmeans$Cluster)]
                  t_cinfo$Cluster   <- t_cmap[as.character(t_cinfo$Cluster)]

                  t_cl_lbls <- paste0(
                    "Cluster ", t_cinfo$Cluster[order(t_cinfo$Cluster)],
                    " (n=",    t_cinfo$n_sites[order(t_cinfo$Cluster)], ")")
                  t_long_df$ClusterLabel <- factor(
                    paste0("Cluster ", t_long_df$Cluster,
                           " (n=", t_cinfo$n_sites[match(t_long_df$Cluster, t_cinfo$Cluster)], ")"),
                    levels = t_cl_lbls)
                  t_cmeans$ClusterLabel <- factor(
                    paste0("Cluster ", t_cmeans$Cluster,
                           " (n=", t_cinfo$n_sites[match(t_cmeans$Cluster, t_cinfo$Cluster)], ")"),
                    levels = t_cl_lbls)

                  t_long_df$Condition  <- factor(t_long_df$Condition,  levels = cond_order)
                  t_cmeans$Condition   <- factor(t_cmeans$Condition,   levels = cond_order)

                  t_zscore_wide_d <- as.data.frame(t_zscore_mat)
                  colnames(t_zscore_wide_d) <- paste0("ZScore_", cond_order)
                  t_zscore_wide_d$PTM_ID        <- t_abs_ids
                  t_zscore_wide_d$PTM_Label     <- t_lbls
                  t_zscore_wide_d$Gene          <- t_genes
                  t_zscore_wide_d$Accession     <- t_accs
                  t_zscore_wide_d$Site          <- t_sites
                  t_zscore_wide_d$Modification  <- t_mods
                  t_zscore_wide_d$Cluster       <- t_cmap[as.character(t_clusters)]
                  t_zscore_wide_d$MaxMembership <- t_mem_max

                  t_avg_wide_d <- as.data.frame(t_mean_mat[t_abs_ids, , drop = FALSE])
                  colnames(t_avg_wide_d) <- paste0("MeanLog2_", cond_order)

                  ptrnd$total_zscore_df     <- t_long_df
                  ptrnd$total_cluster_means <- t_cmeans
                  ptrnd$total_cluster_info  <- t_cinfo
                  ptrnd$total_zscore_wide   <- t_zscore_wide_d
                  ptrnd$total_avg_wide      <- t_avg_wide_d
                  ptrnd$total_memberships   <- t_mems
                  logger(paste0("PTM Trends [DIA PTM Total]: Done — ",
                                t_n_sites, " abs sites, ", t_k, " clusters"))
                }
              } else {
                logger(paste0("PTM Trends [DIA PTM Total]: Only ", nrow(t_mean_mat),
                              " sites — fewer than k=", k, ", skipping."))
              }
            }
          }, error = function(e) {
            logger(paste0("PTM Trends [DIA PTM Total]: failed — ", conditionMessage(e)))
            ptrnd$total_zscore_df <- NULL
          })
        }

        rv$provenance$dia_ptm_trends_run <- list(
          ptm_mode          = cache$ptm_mode %||% "with_global",
          matrix_clustered  = if (.primary == "adjusted") "stoichiometry-adjusted" else "absolute PTM",
          stringency        = if (identical(ptrnd$stringency, "strict")) "strict" else "standard",
          min_valid_per_cond = .min_valid_eff,
          n_sites_input     = nrow(if (.primary == "adjusted") cache$ptm_adj_matrix else cache$ptm_abs_matrix),
          n_sites_clustered = n_sites_cl,
          n_sites_excluded  = ptrnd$n_omitted,
          n_imputed         = 0L,
          n_clusters        = k,
          k_selection       = "user-specified; not determined by a cluster-validity criterion",
          fuzzifier         = fuzziness_m,
          fuzzifier_source  = fuzzifier_source,
          fuzzifier_estimate = .m_est_adj,
          n_sites_zero_variance = n_zero_var_adj,
          seed              = 42L,
          membership_thresh = ptrnd$membership_thresh,
          condition_names   = cond_order,
          cluster_order_rule = "descending mean z-score at the last condition",
          site_gate          = ptrnd$site_gate,
          site_gate_applied  = ptrnd$site_gate_applied
        )
        ptrnd$status <- paste0("DIA PTM Trends: ", n_sites_cl, " sites, ", k, " clusters.")
        logger(paste0("PTM Trends [DIA PTM]: Done — ", n_sites_cl, " sites, ", k, " clusters."))
        incProgress(0.5, detail = "Done.")
      })
      return()
    }
    logger(paste0("PTM Trends: Run clicked but design_option is '",
                  as.character(ptrnd$design_option %||% "NULL"),
                  "' - this handler serves the DIA/TMT PTM cache route only; no action taken."))
  })

  # ===========================================================================
  # HIGHLIGHT CONTROL
  # ===========================================================================
  observeEvent(input$ptm_trends_gene_search_btn, {
    gene <- trimws(input$ptm_trends_gene_search)
    ptrnd$highlight_ptm  <- NULL
    ptrnd$highlight_gene <- if (nchar(gene) == 0) NULL else gene
    ptrnd$render_count <- isolate(ptrnd$render_count) + 1L
  })

  observeEvent(input$ptm_trends_ptm_search_btn, {
    ptm_q <- trimws(input$ptm_trends_ptm_search)
    ptrnd$highlight_gene <- NULL
    ptrnd$highlight_ptm  <- if (nchar(ptm_q) == 0) NULL else ptm_q
    ptrnd$render_count <- isolate(ptrnd$render_count) + 1L
  })

  observeEvent(input$ptm_trends_membership, {
    ptrnd$membership_thresh <- input$ptm_trends_membership
    if (!is.null(ptrnd$zscore_df)) {
      ptrnd$render_count <- isolate(ptrnd$render_count) + 1L
    }
  }, ignoreInit = TRUE)

  # ===========================================================================
  # PROFILE PLOT — renderPlotly
  # ===========================================================================
  output$ptm_trends_profile_plotly <- renderPlotly({
    req(ptrnd$zscore_df, ptrnd$cluster_means)
    ptrnd$render_count

    long_df       <- ptrnd$zscore_df
    cluster_means <- ptrnd$cluster_means
    cluster_info  <- ptrnd$cluster_info
    highlight     <- ptrnd$highlight_gene
    highlight_ptm <- ptrnd$highlight_ptm
    exact_match   <- isTRUE(input$ptm_trends_exact_match)
    col_labels    <- levels(long_df$Condition)
    k             <- ptrnd$n_clusters
    mem_thresh    <- ptrnd$membership_thresh

    clusters_sorted <- sort(unique(long_df$Cluster))

    plot_list <- lapply(clusters_sorted, function(cl) {
      cl_data  <- long_df[long_df$Cluster == cl, ]
      cl_means <- cluster_means[cluster_means$Cluster == cl, ]
      cl_means <- cl_means[order(cl_means$CondIdx), ]
      n_sites  <- cluster_info$n_sites[cluster_info$Cluster == cl]

      # Membership threshold filter
      if (mem_thresh > 0) {
        cl_data <- cl_data[cl_data$Membership >= mem_thresh, ]
      }
      n_displayed <- length(unique(cl_data$PTM_ID))

      p <- plot_ly()

      # Consolidated PTM site traces (single trace with NA breaks)
      unique_ptm <- unique(cl_data$PTM_ID)
      n_conds    <- length(col_labels)
      all_x    <- vector("numeric",   length(unique_ptm) * (n_conds + 1))
      all_y    <- vector("numeric",   length(unique_ptm) * (n_conds + 1))
      all_text <- vector("character", length(unique_ptm) * (n_conds + 1))
      idx <- 1L

      for (ptm_id in unique_ptm) {
        ptm_df <- cl_data[cl_data$PTM_ID == ptm_id, ]
        ptm_df <- ptm_df[order(ptm_df$CondIdx), ]
        lbl    <- ptm_df$PTM_Label[1]
        n_r    <- nrow(ptm_df)
        rng    <- idx:(idx + n_r - 1L)
        all_x[rng] <- ptm_df$CondIdx
        all_y[rng] <- ptm_df$ZScore
        all_text[rng] <- paste0(
          "<b>PTM:</b> ", lbl,
          "<br><b>Gene:</b> ", ptm_df$Gene,
          "<br><b>Accession:</b> ", ptm_df$Accession,
          "<br><b>Site:</b> ", ptm_df$Site,
          "<br><b>Condition:</b> ", ptm_df$Condition,
          "<br><b>Z-Score:</b> ", round(ptm_df$ZScore, 3),
          "<br><b>Membership:</b> ", round(ptm_df$Membership, 3)
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
        "<b>Cluster ", cl, " Average</b>",
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

      # Highlight: PTM-level takes priority over gene-level
      hl_active <- NULL
      if (!is.null(highlight_ptm) && nchar(highlight_ptm) > 0) {
        if (exact_match) {
          hl_df <- cl_data[toupper(cl_data$PTM_Label) == toupper(highlight_ptm), ]
        } else {
          hl_df <- cl_data[grepl(highlight_ptm, cl_data$PTM_Label, ignore.case = TRUE), ]
        }
        hl_active <- highlight_ptm
      } else if (!is.null(highlight) && nchar(highlight) > 0) {
        if (exact_match) {
          hl_df <- cl_data[toupper(cl_data$Gene) == toupper(highlight), ]
        } else {
          hl_df <- cl_data[grepl(highlight, cl_data$Gene, ignore.case = TRUE), ]
        }
        hl_active <- highlight
      } else {
        hl_df <- cl_data[0, ]
      }

      if (!is.null(hl_active) && nrow(hl_df) > 0) {
        hl_ptms <- unique(hl_df$PTM_ID)
        for (hp in hl_ptms) {
          h_df <- hl_df[hl_df$PTM_ID == hp, ]
          h_df <- h_df[order(h_df$CondIdx), ]
          hl_hover <- paste0(
            "<b>\u2b50 ", h_df$PTM_Label[1], "</b>",
            "<br><b>Gene:</b> ", h_df$Gene[1],
            "<br><b>Accession:</b> ", h_df$Accession[1],
            "<br><b>Site:</b> ", h_df$Site[1],
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

      # Per-subplot layout
      anno_text <- if (n_displayed < n_sites)
        paste0("<b>Cluster ", cl, "</b> (", n_displayed, "/", n_sites, ")")
      else
        paste0("<b>Cluster ", cl, "</b> (n=", n_sites, ")")

      p <- layout(p,
        xaxis = list(
          tickvals = seq_along(col_labels),
          ticktext = col_labels,
          tickangle = 45,
          title = ""
        ),
        yaxis = list(title = if (cl == clusters_sorted[1]) "Z-Score" else ""),
        annotations = list(
          list(
            text = anno_text,
            xref = "paper", yref = "paper",
            x = 0.5, y = 1.08,
            showarrow = FALSE,
            font = list(size = 12, color = "#2c3e50")
          )
        )
      )
      p
    })

    n_cols <- min(3, length(plot_list))
    n_rows <- ceiling(length(plot_list) / n_cols)

    fig <- subplot(plot_list,
                   nrows = n_rows,
                   shareX = TRUE, shareY = TRUE,
                   titleY = TRUE,
                   margin = 0.06)

    fig <- layout(fig,
      title = list(
        text = paste0(if (.main_is_adjusted()) "PTM Stoichiometric" else "Absolute PTM",
                      " Trend Profiles \u2014 ",
                      "Fuzzy c-Means (k = ", k, ", m = ",
                      round(if (!is.null(input$ptrnd_fuzziness)) input$ptrnd_fuzziness else 2, 1),
                      ", seed 42)<br>",
                      "<span style='font-size:11px'>Clusters numbered by descending mean ",
                      "z-score at the last condition. Traces shown for sites with max ",
                      "membership \u2265 ", round(mem_thresh, 2),
                      "; counts are displayed/total.</span>"),
        font = list(size = 16, color = "#2c3e50")
      ),
      showlegend = FALSE,
      hovermode = "closest",
      margin = list(t = 80, b = 80)
    )

    fig
  })

  # ===========================================================================
  # CLUSTER MEMBERSHIP TABLE
  # ===========================================================================
  output$ptm_trends_cluster_table <- renderDT({
    req(ptrnd$zscore_df)
    long_df <- ptrnd$zscore_df

    # One row per PTM site
    ptm_df <- long_df[!duplicated(long_df$PTM_ID),
                       c("PTM_ID", "PTM_Label", "Gene", "Accession", "Site",
                         "Modification", "Cluster", "Membership")]
    ptm_df <- ptm_df[order(ptm_df$Cluster, ptm_df$Gene, ptm_df$Site), ]
    rownames(ptm_df) <- NULL
    ptm_df$Membership <- round(ptm_df$Membership, 3)

    # Append Z-scores per condition for the table
    zw <- ptrnd$zscore_wide
    if (!is.null(zw)) {
      zcols <- grep("^ZScore_", colnames(zw), value = TRUE)
      z_sub <- zw[, c("PTM_ID", zcols), drop = FALSE]
      for (zc in zcols) z_sub[[zc]] <- round(z_sub[[zc]], 3)
      ptm_df <- merge(ptm_df, z_sub, by = "PTM_ID", all.x = TRUE, sort = FALSE)
      ptm_df <- ptm_df[order(ptm_df$Cluster, ptm_df$Gene, ptm_df$Site), ]
    }

    datatable(ptm_df,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        scrollY = "400px",
        dom = "frtip",
        order = list(list(6, "asc"))
      ),
      filter = "top",
      rownames = FALSE,
      class = "compact stripe hover"
    )
  })

  # ===========================================================================
  # TOTAL PTM PROFILE PLOT — renderPlotly
  # ===========================================================================
  output$ptm_total_profile_plotly <- renderPlotly({
    req(ptrnd$total_zscore_df, ptrnd$total_cluster_means)
    ptrnd$render_count

    long_df       <- ptrnd$total_zscore_df
    cluster_means <- ptrnd$total_cluster_means
    cluster_info  <- ptrnd$total_cluster_info
    highlight     <- ptrnd$highlight_gene
    highlight_ptm <- ptrnd$highlight_ptm
    exact_match   <- isTRUE(input$ptm_trends_exact_match)
    col_labels    <- levels(long_df$Condition)
    k             <- length(unique(long_df$Cluster))
    mem_thresh    <- ptrnd$membership_thresh

    clusters_sorted <- sort(unique(long_df$Cluster))

    plot_list <- lapply(clusters_sorted, function(cl) {
      cl_data  <- long_df[long_df$Cluster == cl, ]
      cl_means <- cluster_means[cluster_means$Cluster == cl, ]
      cl_means <- cl_means[order(cl_means$CondIdx), ]
      n_sites  <- cluster_info$n_sites[cluster_info$Cluster == cl]

      if (mem_thresh > 0) cl_data <- cl_data[cl_data$Membership >= mem_thresh, ]
      n_displayed <- length(unique(cl_data$PTM_ID))

      p <- plot_ly()

      unique_ptm <- unique(cl_data$PTM_ID)
      n_conds    <- length(col_labels)
      all_x <- vector("numeric",   length(unique_ptm) * (n_conds + 1))
      all_y <- vector("numeric",   length(unique_ptm) * (n_conds + 1))
      all_text <- vector("character", length(unique_ptm) * (n_conds + 1))
      idx <- 1L
      for (ptm_id in unique_ptm) {
        ptm_df <- cl_data[cl_data$PTM_ID == ptm_id, ]
        ptm_df <- ptm_df[order(ptm_df$CondIdx), ]
        lbl    <- ptm_df$PTM_Label[1]
        n_r    <- nrow(ptm_df)
        rng    <- idx:(idx + n_r - 1L)
        all_x[rng] <- ptm_df$CondIdx
        all_y[rng] <- ptm_df$ZScore
        all_text[rng] <- paste0(
          "<b>PTM:</b> ", lbl,
          "<br><b>Gene:</b> ", ptm_df$Gene,
          "<br><b>Condition:</b> ", ptm_df$Condition,
          "<br><b>Z-Score:</b> ", round(ptm_df$ZScore, 3),
          "<br><b>Membership:</b> ", round(ptm_df$Membership, 3))
        idx <- idx + n_r
        all_x[idx] <- NA_real_; all_y[idx] <- NA_real_; all_text[idx] <- ""
        idx <- idx + 1L
      }
      all_x <- all_x[seq_len(idx - 1L)]
      all_y <- all_y[seq_len(idx - 1L)]
      all_text <- all_text[seq_len(idx - 1L)]

      p <- add_trace(p, x = all_x, y = all_y,
        type = "scatter", mode = "lines",
        line = list(color = "rgba(170,170,170,0.15)", width = 0.8),
        hoverinfo = "text", text = all_text,
        showlegend = FALSE, connectgaps = FALSE)

      mean_hover <- paste0(
        "<b>Cluster ", cl, " Average</b>",
        "<br><b>Condition:</b> ", cl_means$Condition,
        "<br><b>Mean Z-Score:</b> ", round(cl_means$MeanZ, 3))
      p <- add_trace(p, x = cl_means$CondIdx, y = cl_means$MeanZ,
        type = "scatter", mode = "lines+markers",
        line = list(color = "#e74c3c", width = 3),
        marker = list(color = "#e74c3c", size = 7),
        hoverinfo = "text", text = mean_hover,
        name = paste0("Cluster ", cl, " mean"), showlegend = FALSE)

      hl_active <- NULL
      if (!is.null(highlight_ptm) && nchar(highlight_ptm) > 0) {
        if (exact_match) hl_df <- cl_data[toupper(cl_data$PTM_Label) == toupper(highlight_ptm), ]
        else hl_df <- cl_data[grepl(highlight_ptm, cl_data$PTM_Label, ignore.case = TRUE), ]
        hl_active <- highlight_ptm
      } else if (!is.null(highlight) && nchar(highlight) > 0) {
        if (exact_match) hl_df <- cl_data[toupper(cl_data$Gene) == toupper(highlight), ]
        else hl_df <- cl_data[grepl(highlight, cl_data$Gene, ignore.case = TRUE), ]
        hl_active <- highlight
      } else { hl_df <- cl_data[0, ] }

      if (!is.null(hl_active) && nrow(hl_df) > 0) {
        for (hp in unique(hl_df$PTM_ID)) {
          h_df <- hl_df[hl_df$PTM_ID == hp, ]
          h_df <- h_df[order(h_df$CondIdx), ]
          hl_hover <- paste0(
            "<b>\u2b50 ", h_df$PTM_Label[1], "</b>",
            "<br><b>Condition:</b> ", h_df$Condition,
            "<br><b>Z-Score:</b> ", round(h_df$ZScore, 3))
          p <- add_trace(p, x = h_df$CondIdx, y = h_df$ZScore,
            type = "scatter", mode = "lines+markers",
            line = list(color = "#f39c12", width = 2.5),
            marker = list(color = "#f39c12", size = 6),
            hoverinfo = "text", text = hl_hover, showlegend = FALSE)
        }
      }

      anno_text <- if (n_displayed < n_sites)
        paste0("<b>Cluster ", cl, "</b> (", n_displayed, "/", n_sites, ")")
      else paste0("<b>Cluster ", cl, "</b> (n=", n_sites, ")")

      p <- layout(p,
        xaxis = list(tickvals = seq_along(col_labels), ticktext = col_labels,
                     tickangle = 45, title = ""),
        yaxis = list(title = if (cl == clusters_sorted[1]) "Z-Score" else ""),
        annotations = list(list(text = anno_text, xref = "paper", yref = "paper",
                                x = 0.5, y = 1.08, showarrow = FALSE,
                                font = list(size = 12, color = "#2c3e50"))))
      p
    })

    n_cols <- min(3, length(plot_list))
    n_rows <- ceiling(length(plot_list) / n_cols)
    fig <- subplot(plot_list, nrows = n_rows, shareX = TRUE, shareY = TRUE,
                   titleY = TRUE, margin = 0.06)
    fig <- layout(fig,
      title = list(
        text = paste0("Total PTM Trend Profiles \u2014 ",
                      "Fuzzy c-Means (k = ", k, ", m = ",
                      round(if (!is.null(input$ptrnd_fuzziness)) input$ptrnd_fuzziness else 2, 1),
                      ", seed 42)<br>",
                      "<span style='font-size:11px'>Clusters numbered by descending mean ",
                      "z-score at the last condition. Traces shown for sites with max ",
                      "membership \u2265 ", round(mem_thresh, 2),
                      "; counts are displayed/total.</span>"),
        font = list(size = 16, color = "#2c3e50")),
      showlegend = FALSE, hovermode = "closest",
      margin = list(t = 80, b = 80))
    fig
  })

  # ===========================================================================
  # TOTAL PTM CLUSTER TABLE
  # ===========================================================================
  output$ptm_total_cluster_table <- renderDT({
    req(ptrnd$total_zscore_df)
    long_df <- ptrnd$total_zscore_df

    ptm_df <- long_df[!duplicated(long_df$PTM_ID),
                       c("PTM_ID", "PTM_Label", "Gene", "Accession", "Site",
                         "Modification", "Cluster", "Membership")]
    ptm_df <- ptm_df[order(ptm_df$Cluster, ptm_df$Gene, ptm_df$Site), ]
    rownames(ptm_df) <- NULL
    ptm_df$Membership <- round(ptm_df$Membership, 3)

    zw <- ptrnd$total_zscore_wide
    if (!is.null(zw)) {
      zcols <- grep("^ZScore_", colnames(zw), value = TRUE)
      z_sub <- zw[, c("PTM_ID", zcols), drop = FALSE]
      for (zc in zcols) z_sub[[zc]] <- round(z_sub[[zc]], 3)
      ptm_df <- merge(ptm_df, z_sub, by = "PTM_ID", all.x = TRUE, sort = FALSE)
      ptm_df <- ptm_df[order(ptm_df$Cluster, ptm_df$Gene, ptm_df$Site), ]
    }

    datatable(ptm_df,
      options = list(pageLength = 25, scrollX = TRUE, scrollY = "400px",
                     dom = "frtip", order = list(list(6, "asc"))),
      filter = "top", rownames = FALSE, class = "compact stripe hover")
  })

  # ===========================================================================
  # EXCLUDED SITES (failed the eligibility rule - shown, not silently dropped)
  # ===========================================================================
  output$ptm_trends_excluded_table <- renderDT({
    ex <- ptrnd$excluded_df
    if (is.null(ex) || nrow(ex) == 0)
      return(datatable(data.frame(Message = "No sites were excluded."),
                       options = list(dom = "t"), rownames = FALSE))
    datatable(ex,
      options = list(pageLength = 25, scrollX = TRUE, scrollY = "300px", dom = "frtip"),
      filter = "top", rownames = FALSE, class = "compact stripe hover")
  })

  output$ptm_trends_download_excluded <- downloadHandler(
    filename = function() {
      paste0("PTM_Trends_ExcludedSites_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      ex <- ptrnd$excluded_df
      if (is.null(ex)) ex <- data.frame()
      write.csv(ex, file, row.names = FALSE)
      logger(paste0("PTM Trends: excluded-site list exported — ", nrow(ex), " sites"))
    }
  )

  # ===========================================================================
  # CSV DOWNLOAD
  # ===========================================================================
  # The main clustering's layer follows the cache (ptm_primary_layer()): the
  # tab title, download label and file name name it.
  .main_is_adjusted <- reactive(ptm_has_global(rv$dia_ptm_master_cache))

  output$ptm_trends_main_tab_title <- renderText(
    if (.main_is_adjusted()) "Adjusted PTM Trends" else "Absolute PTM Trends")

  output$ptm_trends_download_ui <- renderUI(tagList(
    downloadButton("ptm_trends_download_csv",
                   if (.main_is_adjusted()) "Download Adjusted CSV" else "Download Absolute PTM CSV",
                   class = "btn-primary w-100 mb-2"),
    # PTM-only: the main clustering is the absolute one; there is no separate Total PTM run
    if (.main_is_adjusted())
      downloadButton("ptm_trends_download_total_csv", "Download Total PTM CSV",
                     class = "btn-outline-primary w-100 mb-2")))

  observe({
    if (.main_is_adjusted()) bslib::nav_show("ptm_trends_tabs", "total", session = session)
    else                     bslib::nav_hide("ptm_trends_tabs", "total", session = session)
  })

  output$ptm_trends_download_csv <- downloadHandler(
    filename = function() {
      paste0(if (.main_is_adjusted()) "PTM_Stoichiometric_TrendClusters_k" else "PTM_Absolute_TrendClusters_k",
             ptrnd$n_clusters, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(ptrnd$zscore_wide, ptrnd$avg_wide)

      zw <- ptrnd$zscore_wide
      aw <- ptrnd$avg_wide

      export_df <- data.frame(
        Cluster      = zw$Cluster,
        PTM_ID       = zw$PTM_ID,
        PTM_Label    = zw$PTM_Label,
        Gene         = zw$Gene,
        Accession    = zw$Accession,
        Site         = zw$Site,
        Modification = zw$Modification,
        MaxMembership = round(zw$MaxMembership, 4),
        stringsAsFactors = FALSE
      )

      for (col in colnames(aw)) export_df[[col]] <- aw[[col]]
      zcols <- grep("^ZScore_", colnames(zw), value = TRUE)
      for (col in zcols) export_df[[col]] <- zw[[col]]

      export_df <- export_df[order(export_df$Cluster, export_df$Gene, export_df$Site), ]
      write.csv(export_df, file, row.names = FALSE)
      logger(paste0("PTM Trends: CSV exported — ", nrow(export_df), " PTM sites"))
    }
  )

  # ===========================================================================
  # TOTAL PTM CSV DOWNLOAD
  # ===========================================================================
  output$ptm_trends_download_total_csv <- downloadHandler(
    filename = function() {
      paste0("PTM_Total_TrendClusters_k", ptrnd$n_clusters, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(ptrnd$total_zscore_wide, ptrnd$total_avg_wide)

      zw <- ptrnd$total_zscore_wide
      aw <- ptrnd$total_avg_wide

      export_df <- data.frame(
        Cluster      = zw$Cluster,
        PTM_ID       = zw$PTM_ID,
        PTM_Label    = zw$PTM_Label,
        Gene         = zw$Gene,
        Accession    = zw$Accession,
        Site         = zw$Site,
        Modification = zw$Modification,
        MaxMembership = round(zw$MaxMembership, 4),
        stringsAsFactors = FALSE
      )

      for (col in colnames(aw)) export_df[[col]] <- aw[[col]]
      zcols <- grep("^ZScore_", colnames(zw), value = TRUE)
      for (col in zcols) export_df[[col]] <- zw[[col]]

      export_df <- export_df[order(export_df$Cluster, export_df$Gene, export_df$Site), ]
      write.csv(export_df, file, row.names = FALSE)
      logger(paste0("PTM Trends [Total]: CSV exported \u2014 ", nrow(export_df), " PTM sites"))
    }
  )

  # ===========================================================================
  # RDS SAVE / LOAD
  # ===========================================================================
  output$ptm_trends_save_rds <- downloadHandler(
    filename = function() {
      paste0(if (.main_is_adjusted()) "PTM_Stoichiometric_TrendClusters_k" else "PTM_Absolute_TrendClusters_k",
             ptrnd$n_clusters, "_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(ptrnd$zscore_df)
      save_data <- list(
        zscore_df         = ptrnd$zscore_df,
        cluster_means     = ptrnd$cluster_means,
        cluster_info      = ptrnd$cluster_info,
        zscore_wide       = ptrnd$zscore_wide,
        avg_wide          = ptrnd$avg_wide,
        memberships       = ptrnd$memberships,
        n_clusters        = ptrnd$n_clusters,
        cond_order        = ptrnd$cond_order,
        conditions        = ptrnd$conditions,
        membership_thresh = ptrnd$membership_thresh,
        n_files           = ptrnd$n_files,
        annotation_map    = ptrnd$annotation_map,
        design_option     = ptrnd$design_option,
        n_omitted         = ptrnd$n_omitted,
        n_imputed         = ptrnd$n_imputed,
        lod_value         = ptrnd$lod_value,
        rule_version      = "2026-08-eligibility-v1",
        stringency        = ptrnd$stringency,
        stringency_effective = ptrnd$stringency_effective,
        site_gate         = ptrnd$site_gate,
        seed              = 42L,
        timestamp         = Sys.time(),
        provenance        = reactiveValuesToList(rv$provenance)
      )
      saveRDS(save_data, file)
      logger("PTM Trends: RDS saved")
    }
  )

  observeEvent(input$ptm_trends_load_rds, {
    req(input$ptm_trends_load_rds)

    tryCatch({
      loaded <- readRDS(input$ptm_trends_load_rds$datapath)

      ptrnd$zscore_df       <- loaded$zscore_df
      ptrnd$cluster_means   <- loaded$cluster_means
      ptrnd$cluster_info    <- loaded$cluster_info
      ptrnd$zscore_wide     <- loaded$zscore_wide
      ptrnd$avg_wide        <- loaded$avg_wide
      ptrnd$memberships     <- loaded$memberships
      ptrnd$n_clusters      <- loaded$n_clusters
      ptrnd$cond_order      <- loaded$cond_order
      ptrnd$conditions      <- loaded$conditions
      ptrnd$membership_thresh <- if (!is.null(loaded$membership_thresh))
                                   loaded$membership_thresh else 0.7
      if (!is.null(loaded$annotation_map)) ptrnd$annotation_map <- loaded$annotation_map
      if (!is.null(loaded$design_option))  ptrnd$design_option  <- loaded$design_option
      ptrnd$n_omitted <- if (!is.null(loaded$n_omitted)) loaded$n_omitted else 0L
      ptrnd$n_imputed <- if (!is.null(loaded$n_imputed)) loaded$n_imputed else 0L
      ptrnd$lod_value <- if (!is.null(loaded$lod_value)) loaded$lod_value else NA_real_
      # Bundles written before the eligibility rewrite carry clusters computed with
      # column-mean-filled z-scores.  Say so rather than presenting them as current.
      if (is.null(loaded$rule_version)) {
        ptrnd$stringency_effective <- NULL
        showNotification(paste0("This saved analysis predates the current missing-value rule. ",
                                "Its clusters were computed with column-mean imputation of ",
                                "missing conditions. Re-run to apply the current rule."),
                         type = "warning", duration = 12)
        logger("PTM Trends: loaded a pre-eligibility-rule bundle (rule_version absent).")
      } else {
        ptrnd$stringency           <- loaded$stringency %||% "standard"
        ptrnd$stringency_effective <- loaded$stringency_effective
      }

      n_total <- sum(ptrnd$cluster_info$n_sites)
      ptrnd$status <- paste0("Loaded saved analysis: Fuzzy c-means, k=",
                             ptrnd$n_clusters, ", ", n_total, " PTM sites")

      updateNumericInput(session, "ptm_trends_n_clusters", value = ptrnd$n_clusters)
      if (!is.null(loaded$membership_thresh))
        updateSliderInput(session, "ptm_trends_membership", value = loaded$membership_thresh)

      # Restore provenance if present. Both slots are restored because the legacy
      # label-free route writes $ptm_trends while the DIA/TMT route writes
      # $dia_ptm_trends_run; a bundle carries whichever its run produced.
      if (!is.null(loaded$provenance) && !is.null(loaded$provenance$ptm_trends)) {
        rv$provenance$ptm_trends <- loaded$provenance$ptm_trends
      }
      if (!is.null(loaded$provenance) && !is.null(loaded$provenance$dia_ptm_trends_run)) {
        rv$provenance$dia_ptm_trends_run <- loaded$provenance$dia_ptm_trends_run
      }

      ptrnd$render_count <- isolate(ptrnd$render_count) + 1L
      logger(paste0("PTM Trends: Loaded RDS (Fuzzy c-means, k=",
                    ptrnd$n_clusters, ", ", n_total, " PTM sites)"))
      showNotification("PTM trend analysis loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste0("Error loading RDS: ", e$message), type = "error")
      logger(paste0("PTM Trends: RDS load error — ", e$message))
    })
  })

  # ===========================================================================
  # STATUS
  # ===========================================================================
  output$ptm_trends_status <- renderPrint({
    if (is.null(ptrnd$zscore_df)) {
      cat("Upload nested PD CSV files, set design & condition order, and run.\n")
      cat(paste0("\nStatus: ", ptrnd$status, "\n"))
    } else {
      info <- ptrnd$cluster_info
      info <- info[order(info$Cluster), ]
      n_total <- sum(info$n_sites)

      cat(paste0("Method: Fuzzy c-means on ",
                 if (ptm_has_global(rv$dia_ptm_master_cache)) "stoichiometry-adjusted"
                 else "absolute PTM (PTM only, no global proteome)",
                 " z-scores\n"))
      cat(paste0("Missing-value handling: none - no imputation. Sites must meet the ",
                 if (identical(ptrnd$stringency, "strict")) "strict (100%)"
                 else "standard (>=50%)",
                 " valid-value rule in every condition.\n"))
      if (!is.null(ptrnd$stringency_effective))
        cat(paste0("Realised threshold per condition: ",
                   paste0(names(ptrnd$stringency_effective), ":",
                          as.integer(ptrnd$stringency_effective), collapse = ", "), "\n"))
      cat(paste0("Total PTM sites clustered: ", n_total, "\n"))
      cat(paste0("Clusters (k): ", ptrnd$n_clusters, "\n"))
      cat(paste0("Random seed: 42\n"))
      cat(paste0("Conditions (", length(levels(ptrnd$zscore_df$Condition)), "): ",
                 paste(levels(ptrnd$zscore_df$Condition), collapse = " \u2192 "), "\n"))

      if (ptrnd$n_omitted > 0)
        cat(paste0("PTM sites excluded by the eligibility rule: ", ptrnd$n_omitted,
                    " (listed in the excluded-sites table)\n"))
      if (!is.null(ptrnd$site_gate_applied))
        cat(paste0("Site gate: ", ptrnd$site_gate_applied, "\n"))

      cat("\n")
      for (i in seq_len(nrow(info))) {
        cat(paste0("  Cluster ", info$Cluster[i], ": ",
                    info$n_sites[i], " PTM sites\n"))
      }

      # Report highlight
      hl     <- ptrnd$highlight_gene
      hl_ptm <- ptrnd$highlight_ptm
      exact_match <- isTRUE(input$ptm_trends_exact_match)
      long_df <- ptrnd$zscore_df

      if (!is.null(hl_ptm) && nchar(hl_ptm) > 0) {
        if (exact_match) {
          matches <- long_df[toupper(long_df$PTM_Label) == toupper(hl_ptm), ]
        } else {
          matches <- long_df[grepl(hl_ptm, long_df$PTM_Label, ignore.case = TRUE), ]
        }
        unique_matches <- unique(matches[, c("PTM_ID", "PTM_Label", "Gene", "Cluster")])
        if (nrow(unique_matches) > 0) {
          cat(paste0("\nHighlighted PTM '", hl_ptm, "':\n"))
          for (j in seq_len(nrow(unique_matches))) {
            cat(paste0("  ", unique_matches$PTM_Label[j],
                        " (", unique_matches$PTM_ID[j],
                        ") -> Cluster ", unique_matches$Cluster[j], "\n"))
          }
        } else {
          cat(paste0("\nHighlighted PTM '", hl_ptm, "': no matches\n"))
        }
      } else if (!is.null(hl) && nchar(hl) > 0) {
        if (exact_match) {
          matches <- long_df[toupper(long_df$Gene) == toupper(hl), ]
        } else {
          matches <- long_df[grepl(hl, long_df$Gene, ignore.case = TRUE), ]
        }
        unique_matches <- unique(matches[, c("PTM_ID", "PTM_Label", "Gene", "Cluster")])
        if (nrow(unique_matches) > 0) {
          cat(paste0("\nHighlighted '", hl, "':\n"))
          for (j in seq_len(nrow(unique_matches))) {
            cat(paste0("  ", unique_matches$PTM_Label[j],
                        " (", unique_matches$PTM_ID[j],
                        ") -> Cluster ", unique_matches$Cluster[j], "\n"))
          }
        } else {
          cat(paste0("\nHighlighted '", hl, "': no matches\n"))
        }
      }
    }
  })

}
