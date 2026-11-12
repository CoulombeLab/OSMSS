# server_ptm_pca.R =============================================================
# PTM Principal Component Analysis module (v1.0)
#
# Primary tab  : prcomp on log2 absolute (unadjusted) PTM site intensities
# Secondary tab: prcomp on log2 global protein intensities (same samples)
# Third tab    : CV density curves by condition
#
#   TMT PTM  — rv$dia_ptm_master_cache (tmt_ptm_upload wizard via limma)
#   TMT  — rv$ptm_master_cache   (TMT PTM summarization output)
#             $PTM$ProteinLevelData   : Protein, GROUP, originalRun, Abundance
#             $PROTEIN$ProteinLevelData : same structure
#   DIA  — rv$dia_ptm_master_cache (dia_ptm_summarize / maxquant_ptm_summarize)
#             $ptm_abs_matrix   : sites × samples matrix (log2)
#             $global_matrix    : proteins × samples matrix (log2, median-centred)
#             $design_df        : data.frame(sample_name, condition)
#
# Missing-value strategy (user-controlled):
#   1. Drop sites/proteins where < pct_valid % of samples have a finite value.
#   2. Row-median impute any remaining NAs in surviving rows.
#
# Called from server.R as: server_ptm_pca(input, output, session, rv, logger)
# ==============================================================================

server_ptm_pca <- function(input, output, session, rv, logger) {

  # ============================================================================
  # Reactive storage
  # ============================================================================
  ptm_pca_rv <- reactiveValues(
    # PTM (unadjusted) PCA
    ptm_result   = NULL,   # prcomp() output
    ptm_scores   = NULL,   # data.frame: PC1, PC2, Sample, Group
    ptm_mat      = NULL,   # matrix: rows=samples, cols=sites  (post-filter)
    n_sites_used = NULL,

    # Global protein PCA
    prot_result  = NULL,
    prot_scores  = NULL,
    prot_mat     = NULL,
    n_prots_used = NULL,

    # CV groups (condition -> vector of sample names)
    cv_groups    = NULL,

    # Meta
    data_source    = NULL,  # "tmt" | "dia"
    batch_corrected = FALSE, # TRUE when limma::removeBatchEffect was applied
    status         = "Complete the PTM upload wizard and run Summarize Data to begin."
  )

  # ============================================================================
  # Internal helpers
  # ============================================================================

  # -- Detect which summarization cache is active -----------------------------
  .data_src <- function() {
    src <- rv$ptm_data$source
    if (isTRUE(src %in% c("dia_ptm_upload", "tmt_ptm_upload")) && !is.null(rv$dia_ptm_master_cache))
      return("dia")
    NULL
  }

  # -- Filter rows by % valid, then row-median impute remaining NAs -----------
  # mat : rows = features (sites / proteins), cols = samples
  .filter_impute_mat <- function(mat, valid_pct) {
    pct_valid <- rowMeans(!is.na(mat) & is.finite(mat))
    keep      <- pct_valid >= (valid_pct / 100)
    mat       <- mat[keep, , drop = FALSE]
    if (nrow(mat) == 0) return(NULL)

    # Row-median imputation of remaining NAs
    mat <- t(apply(mat, 1, function(r) {
      if (!any(is.na(r))) return(r)
      med <- median(r, na.rm = TRUE)
      r[is.na(r) | !is.finite(r)] <- if (is.finite(med)) med else 0
      r
    }))
    mat
  }

  # -- Run prcomp, return list(pca, mat_t [rows=samples], var_pct) ------------
  .run_prcomp <- function(mat) {
    mat_t <- t(mat)   # rows = samples, cols = features
    pca   <- prcomp(mat_t, center = TRUE, scale. = FALSE)
    sdev  <- pca$sdev
    var_p <- (sdev^2) / sum(sdev^2) * 100
    list(pca = pca, mat_t = mat_t, var_pct = var_p)
  }

  # -- Render a plotly PCA scatter (PC1 vs PC2) -------------------------------
  .render_scatter <- function(pca_out, scores_df, title_str, n_features) {
    var_ex  <- summary(pca_out)$importance[2, ] * 100
    pc1_var <- round(var_ex[1], 1)
    pc2_var <- round(var_ex[min(2, length(var_ex))], 1)

    scores_df$hover <- paste0(
      "<b>Sample:</b> ", scores_df$Sample,
      "<br><b>Group:</b> ", scores_df$Group,
      "<br><b>PC1:</b> ", round(scores_df$PC1, 2),
      "<br><b>PC2:</b> ", round(scores_df$PC2, 2)
    )

    group_colors <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12",
                      "#9b59b6", "#1abc9c", "#e67e22", "#34495e",
                      "#16a085", "#c0392b", "#2980b9", "#27ae60")
    n_groups <- length(unique(scores_df$Group))

    if (n_groups <= 1) {
      p <- plot_ly(
        scores_df, x = ~PC1, y = ~PC2, type = "scatter",
        mode = "markers+text",
        text = ~Sample, textposition = "top center",
        marker = list(size = 12, opacity = 0.8, color = "#3498db",
                      line = list(width = 1, color = "#2c3e50")),
        hoverinfo = "text", hovertext = ~hover
      )
    } else {
      p <- plot_ly(
        scores_df, x = ~PC1, y = ~PC2, type = "scatter",
        mode = "markers+text",
        text = ~Sample, textposition = "top center",
        color = ~Group,
        colors = group_colors[seq_len(min(n_groups, length(group_colors)))],
        marker = list(size = 12, opacity = 0.8,
                      line = list(width = 1, color = "#2c3e50")),
        hoverinfo = "text", hovertext = ~hover
      )
    }

    feat_label <- format(n_features %||% 0, big.mark = ",")
    p %>% layout(
      title = list(
        text = paste0(title_str, "<br><sub>", feat_label, " features used</sub>"),
        font = list(size = 15)
      ),
      xaxis = list(title = paste0("PC1 (", pc1_var, "% variance)"),
                   zeroline = TRUE, zerolinecolor = "#ccc"),
      yaxis = list(title = paste0("PC2 (", pc2_var, "% variance)"),
                   zeroline = TRUE, zerolinecolor = "#ccc"),
      hovermode = "closest"
    )
  }

  # -- Render a plotly scree plot (up to 15 PCs) ------------------------------
  .render_scree <- function(pca_out) {
    var_explained <- summary(pca_out)$importance[2, ] * 100
    cum_var       <- summary(pca_out)$importance[3, ] * 100
    n_pcs         <- min(length(var_explained), 15)

    scree_df <- data.frame(
      PC         = factor(paste0("PC", seq_len(n_pcs)),
                          levels = paste0("PC", seq_len(n_pcs))),
      Variance   = var_explained[seq_len(n_pcs)],
      Cumulative = cum_var[seq_len(n_pcs)]
    )

    plot_ly(scree_df) %>%
      add_bars(x = ~PC, y = ~Variance, name = "Individual",
               marker = list(color = "#3498db", opacity = 0.8)) %>%
      add_trace(x = ~PC, y = ~Cumulative, name = "Cumulative",
                type = "scatter", mode = "lines+markers",
                line   = list(color = "#e74c3c", width = 2),
                marker = list(color = "#e74c3c", size = 8)) %>%
      layout(
        title  = list(text = "Scree Plot \u2014 Variance Explained",
                      font = list(size = 14)),
        xaxis  = list(title = "Principal Component"),
        yaxis  = list(title = "% Variance Explained", range = c(0, 105)),
        legend = list(orientation = "h", y = -0.25),
        barmode = "overlay"
      )
  }

  # ============================================================================
  # Data status banner (sidebar)
  # ============================================================================
  output$ptm_pca_data_status <- renderUI({
    src <- .data_src()

    if (is.null(src)) {
      # Check whether upload happened but summarization is missing
      has_upload <- (!is.null(rv$ptm_data$ptm_long) &&
                       nrow(rv$ptm_data$ptm_long) > 0) ||
                    isTRUE(rv$dia_ptm_wizard_completed)

      if (has_upload) {
        tags$div(class = "alert alert-warning py-2 px-3 small mb-2",
          icon("triangle-exclamation"),
          tags$b(" Summarization required."),
          " Return to ", tags$b("PTM Upload Data"),
          " and run Summarize Data first."
        )
      } else {
        tags$div(class = "alert alert-warning py-2 px-3 small mb-2",
          icon("triangle-exclamation"),
          tags$b(" No PTM data loaded."),
          " Complete the PTM Upload Wizard to begin."
        )
      }


    } else {
      cache  <- rv$dia_ptm_master_cache
      n_sites <- if (!is.null(cache)) cache$n_sites_output else "?"
      n_conds <- if (!is.null(cache)) length(cache$cond_groups) else "?"
      sw_lbl  <- switch(rv$dia_ptm_software %||% "spectronaut",
                   "maxquant"          = "MaxQuant",
                   "spectronaut_pivot" = "Spectronaut Pivot",
                   "Spectronaut")
      tags$div(class = "alert alert-success py-2 px-3 small mb-2",
        icon("circle-check"),
        tags$b(paste0(" ", sw_lbl, " PTM data ready. ")),
        paste0(n_sites, " sites \u2022 ", n_conds, " conditions.")
      )
    }
  })

  # ============================================================================
  # Enable / disable run button
  # ============================================================================
  observe({
    if (!is.null(.data_src())) {
      shinyjs::enable("ptm_pca_run_btn")
    } else {
      shinyjs::disable("ptm_pca_run_btn")
    }
  })

  # ============================================================================
  # Reset on new data
  # ============================================================================
  observeEvent(rv$ptm_data$source, {
    ptm_pca_rv$ptm_result   <- NULL
    ptm_pca_rv$ptm_scores   <- NULL
    ptm_pca_rv$ptm_mat      <- NULL
    ptm_pca_rv$n_sites_used <- NULL
    ptm_pca_rv$prot_result  <- NULL
    ptm_pca_rv$prot_scores  <- NULL
    ptm_pca_rv$prot_mat     <- NULL
    ptm_pca_rv$n_prots_used <- NULL
    ptm_pca_rv$cv_groups    <- NULL
    ptm_pca_rv$data_source  <- NULL
    ptm_pca_rv$status <- "New data detected. Configure settings and run PCA."
    logger("PTM PCA: Data source changed, state reset.")
  }, ignoreInit = TRUE)

  # ============================================================================
  # RUN PCA button — route to correct pipeline
  # ============================================================================
  observeEvent(input$ptm_pca_run_btn, {
    src       <- .data_src()
    valid_pct <- isolate(input$ptm_pca_valid_pct) %||% 50

    if (is.null(src)) {
      showNotification(
        "No summarized PTM data found. Complete the Upload Wizard first.",
        type = "warning")
      return()
    }

    shinyjs::disable("ptm_pca_run_btn")
    on.exit(shinyjs::enable("ptm_pca_run_btn"), add = TRUE)

    # Clear previous results
    ptm_pca_rv$ptm_result  <- NULL
    ptm_pca_rv$prot_result <- NULL
    ptm_pca_rv$cv_groups   <- NULL

    .run_ptm_pca_dia(valid_pct)
    ptm_pca_rv$data_source <- src
  })

  # ============================================================================
  # DIA pipeline
  # ============================================================================
  .run_ptm_pca_dia <- function(valid_pct) {
    cache <- rv$dia_ptm_master_cache
    if (is.null(cache)) {
      showNotification("DIA PTM summarization cache not found.", type = "error")
      return()
    }

    ptm_mat    <- cache$ptm_abs_matrix   # sites × samples, already log2
    global_mat <- cache$global_matrix    # proteins × samples, already log2
    design_df  <- cache$design_df        # sample_name, condition
    cond_groups <- cache$cond_groups     # named list: cond -> c(sample_names)

    if (is.null(ptm_mat) || ncol(ptm_mat) == 0) {
      showNotification("PTM intensity matrix is empty.", type = "error")
      return()
    }

    # Build run-to-condition map
    if (!is.null(design_df) &&
        all(c("sample_name", "condition") %in% colnames(design_df))) {
      run_cond_map <- setNames(as.character(design_df$condition),
                               as.character(design_df$sample_name))
    } else if (!is.null(cond_groups)) {
      run_cond_map <- setNames(
        rep(names(cond_groups), vapply(cond_groups, length, integer(1))),
        unlist(cond_groups)
      )
    } else {
      run_cond_map <- setNames(colnames(ptm_mat), colnames(ptm_mat))
    }

    # ── Restrict matrices to assigned samples only ────────────────────────────
    # ptm_abs_matrix is built from ALL matched samples in the input file, but
    # design_df / cond_groups only contains samples the user actually assigned
    # to a condition. Unassigned samples must be excluded BEFORE PCA so they
    # neither contribute to PC calculation nor appear as spurious "Unknown" points.
    assigned_samples <- names(run_cond_map)
    keep_ptm <- intersect(assigned_samples, colnames(ptm_mat))
    if (length(keep_ptm) < ncol(ptm_mat)) {
      n_dropped <- ncol(ptm_mat) - length(keep_ptm)
      logger(paste0("PTM PCA: Excluding ", n_dropped,
                    " unassigned sample(s) from PTM matrix before PCA."))
      ptm_mat <- ptm_mat[, keep_ptm, drop = FALSE]
    }
    if (!is.null(global_mat) && ncol(global_mat) > 0) {
      keep_glob <- intersect(assigned_samples, colnames(global_mat))
      if (length(keep_glob) < ncol(global_mat)) {
        global_mat <- global_mat[, keep_glob, drop = FALSE]
      }
    }

    withProgress(message = "PTM PCA (DIA) \u2014 Computing...", value = 0, {

      incProgress(0.20, detail = "Filtering PTM sites...")
      ptm_mat_f <- .filter_impute_mat(ptm_mat, valid_pct)

      if (is.null(ptm_mat_f) || nrow(ptm_mat_f) == 0) {
        showNotification(
          paste0("No PTM sites survived the ", valid_pct,
                 "% valid-values filter. Try lowering the threshold."),
          type = "error")
        return()
      }

      n_sites <- nrow(ptm_mat_f)
      n_samps <- ncol(ptm_mat_f)

      if (n_samps < 3) {
        showNotification("Need at least 3 samples to run PCA.", type = "error")
        return()
      }

      # ── Multi-plex batch correction (TMT PTM with ≥2 plexes) ───────────────
      # removeBatchEffect removes the plex mean shift while preserving condition
      # effects, matching the behaviour of the TMT protein PCA (server_pca.R).
      did_bc_ptm  <- FALSE
      is_tmt_ptm  <- isTRUE(isolate(rv$ptm_data$source) == "tmt_ptm_upload")
      if (is_tmt_ptm && "mixture" %in% colnames(design_df)) {
        mix_map     <- setNames(as.character(design_df$mixture),
                                as.character(design_df$sample_name))
        mix_vec_ptm <- mix_map[colnames(ptm_mat_f)]
        n_plexes_bc <- length(unique(mix_vec_ptm[!is.na(mix_vec_ptm)]))
        if (n_plexes_bc >= 2 && requireNamespace("limma", quietly = TRUE)) {
          cond_bc <- run_cond_map[colnames(ptm_mat_f)]
          cond_bc[is.na(cond_bc)] <- "Unknown"
          ptm_mat_f <- tryCatch(
            limma::removeBatchEffect(
              ptm_mat_f,
              batch  = factor(mix_vec_ptm),
              design = model.matrix(~ factor(cond_bc))
            ),
            error = function(e) {
              logger(paste0("PTM PCA: removeBatchEffect failed: ", e$message))
              ptm_mat_f
            }
          )
          did_bc_ptm <- TRUE
          logger(paste0("PTM PCA [TMT]: plex batch correction applied (",
                        n_plexes_bc, " plexes)."))
        }
      }

      incProgress(0.45, detail = paste0("Running PTM PCA (", n_sites, " sites)..."))
      ptm_pca_res <- .run_prcomp(ptm_mat_f)

      cond_vec <- run_cond_map[colnames(ptm_mat_f)]
      cond_vec[is.na(cond_vec)] <- "Unknown"

      scores_df <- as.data.frame(ptm_pca_res$pca$x)
      scores_df$Sample <- colnames(ptm_mat_f)
      scores_df$Group  <- cond_vec

      ptm_pca_rv$ptm_result    <- ptm_pca_res$pca
      ptm_pca_rv$ptm_scores    <- scores_df
      ptm_pca_rv$ptm_mat       <- ptm_pca_res$mat_t
      ptm_pca_rv$n_sites_used  <- n_sites
      ptm_pca_rv$cv_groups     <- split(colnames(ptm_mat_f), cond_vec)
      ptm_pca_rv$batch_corrected <- did_bc_ptm

      # ── Global protein PCA ─────────────────────────────────────────────────
      incProgress(0.65, detail = "Running global protein PCA...")

      if (!is.null(global_mat) && ncol(global_mat) > 0) {
        global_mat_f <- .filter_impute_mat(global_mat, valid_pct)

        if (!is.null(global_mat_f) &&
            nrow(global_mat_f) >= 2 && ncol(global_mat_f) >= 3) {

          # Apply the same plex batch correction to the global protein matrix
          if (is_tmt_ptm && "mixture" %in% colnames(design_df)) {
            mix_map_g  <- setNames(as.character(design_df$mixture),
                                   as.character(design_df$sample_name))
            mix_vec_g  <- mix_map_g[colnames(global_mat_f)]
            n_plexes_g <- length(unique(mix_vec_g[!is.na(mix_vec_g)]))
            if (n_plexes_g >= 2 && requireNamespace("limma", quietly = TRUE)) {
              cond_g <- run_cond_map[colnames(global_mat_f)]
              cond_g[is.na(cond_g)] <- "Unknown"
              global_mat_f <- tryCatch(
                limma::removeBatchEffect(
                  global_mat_f,
                  batch  = factor(mix_vec_g),
                  design = model.matrix(~ factor(cond_g))
                ),
                error = function(e) {
                  logger(paste0("PTM PCA: global removeBatchEffect failed: ", e$message))
                  global_mat_f
                }
              )
            }
          }

          prt_pca_res  <- .run_prcomp(global_mat_f)
          prt_cond_vec <- run_cond_map[colnames(global_mat_f)]
          prt_cond_vec[is.na(prt_cond_vec)] <- "Unknown"

          prot_scores <- as.data.frame(prt_pca_res$pca$x)
          prot_scores$Sample <- colnames(global_mat_f)
          prot_scores$Group  <- prt_cond_vec

          ptm_pca_rv$prot_result  <- prt_pca_res$pca
          ptm_pca_rv$prot_scores  <- prot_scores
          ptm_pca_rv$prot_mat     <- prt_pca_res$mat_t
          ptm_pca_rv$n_prots_used <- nrow(global_mat_f)

        } else {
          logger("PTM PCA: Global protein data insufficient after filter — skipping.")
          ptm_pca_rv$prot_result <- NULL
        }
      }

      incProgress(1, detail = "Done!")
      pipeline_lbl <- if (is_tmt_ptm)
        paste0("TMT PTM", if (did_bc_ptm) " (plex batch-corrected)" else "")
      else {
        sw_lbl <- switch(rv$dia_ptm_software %||% "spectronaut",
                    "maxquant"          = "MaxQuant",
                    "spectronaut_pivot" = "Spectronaut Pivot",
                    "Spectronaut")
        paste0("DIA \u2014 ", sw_lbl)
      }
      ptm_pca_rv$status <- paste0(
        "PTM PCA complete (", pipeline_lbl, ").\n",
        "Samples: ", n_samps, "\n",
        "PTM sites used: ", format(n_sites, big.mark = ","),
        " (after \u2265", valid_pct, "% valid-values filter)\n",
        "PC1: ", round(ptm_pca_res$var_pct[1], 1),
        "% | PC2: ",
        round(ptm_pca_res$var_pct[min(2, length(ptm_pca_res$var_pct))], 1), "%\n",
        if (!is.null(ptm_pca_rv$prot_result))
          paste0("Global proteins used: ",
                 format(ptm_pca_rv$n_prots_used, big.mark = ","), "\n")
        else "Global protein PCA: insufficient data after filter.\n"
      )
      logger(paste0("PTM PCA [", pipeline_lbl, "]: Complete. ",
                    n_samps, " samples, ", n_sites, " PTM sites."))
    })
  }

  # ============================================================================
  # Outputs — PTM PCA scatter
  # ============================================================================
  output$ptm_pca_ptm_scatter <- renderPlotly({
    req(ptm_pca_rv$ptm_result, ptm_pca_rv$ptm_scores)
    title_str <- if (isTRUE(ptm_pca_rv$batch_corrected))
      "PTM PCA \u2014 Absolute Intensities (Plex Batch-Corrected)"
    else
      "PTM PCA \u2014 Absolute Intensities (Unadjusted)"
    .render_scatter(
      ptm_pca_rv$ptm_result,
      ptm_pca_rv$ptm_scores,
      title_str,
      ptm_pca_rv$n_sites_used
    )
  })

  # PTM scree plot
  output$ptm_pca_ptm_scree <- renderPlotly({
    req(ptm_pca_rv$ptm_result)
    .render_scree(ptm_pca_rv$ptm_result)
  })

  # ============================================================================
  # Outputs — Global protein PCA scatter
  # ============================================================================
  output$ptm_pca_prot_scatter <- renderPlotly({
    if (is.null(ptm_pca_rv$prot_result)) {
      plotly_empty() %>% layout(
        title = list(
          text = "Global protein PCA not available<br><sub>Run PTM PCA first, or check if global protein data is present and passes the valid-values filter.</sub>",
          font = list(size = 14)
        )
      )
    } else {
      .render_scatter(
        ptm_pca_rv$prot_result,
        ptm_pca_rv$prot_scores,
        "Global Protein PCA",
        ptm_pca_rv$n_prots_used
      )
    }
  })

  # Protein scree plot
  output$ptm_pca_prot_scree <- renderPlotly({
    if (is.null(ptm_pca_rv$prot_result)) return(plotly_empty())
    .render_scree(ptm_pca_rv$prot_result)
  })

  # ============================================================================
  # Outputs — CV analysis (density curves per condition group)
  # Uses the PTM mat (rows = samples, cols = sites).
  # CV is computed in linear space (2^log2 values) per group.
  # ============================================================================
  output$ptm_pca_cv_plot <- renderPlotly({
    req(ptm_pca_rv$cv_groups, ptm_pca_rv$ptm_mat)

    groups <- ptm_pca_rv$cv_groups
    mat    <- ptm_pca_rv$ptm_mat   # rows = samples, cols = features

    if (length(groups) == 0)
      return(plotly_empty() %>% layout(title = "No group data for CV analysis."))

    group_colors <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12",
                      "#9b59b6", "#1abc9c", "#e67e22", "#34495e",
                      "#16a085", "#c0392b", "#2980b9", "#27ae60")

    cv_data <- list()
    for (grp_name in names(groups)) {
      grp_samps <- groups[[grp_name]]
      grp_idx   <- which(rownames(mat) %in% grp_samps)
      if (length(grp_idx) < 2) next

      grp_mat    <- mat[grp_idx, , drop = FALSE]
      grp_linear <- 2^grp_mat  # back to linear for CV

      col_means <- colMeans(grp_linear, na.rm = TRUE)
      col_sds   <- apply(grp_linear, 2, sd, na.rm = TRUE)
      cv_vals   <- (col_sds / col_means) * 100
      cv_vals   <- cv_vals[is.finite(cv_vals) & !is.na(cv_vals)]

      if (length(cv_vals) > 0)
        cv_data[[grp_name]] <- data.frame(
          CV    = cv_vals,
          Group = grp_name,
          stringsAsFactors = FALSE
        )
    }

    if (length(cv_data) == 0)
      return(plotly_empty() %>% layout(title = "No valid CV data."))

    all_cv <- do.call(rbind, cv_data)
    cv_max <- max(all_cv$CV, na.rm = TRUE) * 1.1

    p <- plot_ly()
    for (i in seq_along(names(groups))) {
      grp    <- names(groups)[i]
      grp_cv <- all_cv[all_cv$Group == grp, ]
      if (nrow(grp_cv) < 2) next

      dens <- tryCatch(
        density(grp_cv$CV, from = 0, to = cv_max, n = 256, na.rm = TRUE),
        error = function(e) NULL
      )
      if (is.null(dens)) next

      col    <- group_colors[((i - 1) %% length(group_colors)) + 1]
      med_cv <- median(grp_cv$CV)
      med_y  <- tryCatch(approx(dens$x, dens$y, xout = med_cv)$y, error = function(e) 0)
      if (!is.finite(med_y)) med_y <- 0

      p <- p %>%
        add_trace(
          x = dens$x, y = dens$y,
          type = "scatter", mode = "lines",
          name      = paste0(grp, " (median=", round(med_cv, 1), "%)"),
          fill      = "tozeroy",
          fillcolor = paste0(col, "33"),
          line      = list(color = col, width = 2)
        ) %>%
        add_segments(
          x = med_cv, xend = med_cv,
          y = 0, yend = med_y,
          line      = list(color = col, width = 2, dash = "dash"),
          showlegend = FALSE,
          hoverinfo  = "text",
          text       = paste0(grp, " median: ", round(med_cv, 1), "%")
        )
    }

    p %>% layout(
      title  = list(text = "Coefficient of Variation \u2014 PTM Sites by Condition",
                    font = list(size = 15)),
      xaxis  = list(title = "CV (%)", rangemode = "tozero"),
      yaxis  = list(title = "Density"),
      legend = list(orientation = "h", y = -0.2),
      hovermode = "x unified"
    )
  })

  # ============================================================================
  # Status text
  # ============================================================================
  output$ptm_pca_status <- renderPrint({
    cat(ptm_pca_rv$status, "\n")
  })
}
