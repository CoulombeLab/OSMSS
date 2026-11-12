# server_trends.R ================================================================
# Server logic for Proteome Discoverer > Single-Shot > Trends Analysis module.
# Parses Proteome Discoverer CSV, lets user pick multiple Abundance columns
# (conditions/timepoints), log2-transforms, z-scores per protein,
# clusters with fuzzy c-means, and plots a faceted profile plot.
# Called from server.R as: server_trends(input, output, session, rv, logger)
# ==============================================================================

server_trends <- function(input, output, session, rv, logger) {

  # ---- Module-local reactive values ----
  trd <- reactiveValues(
    selected_cols   = NULL,    # character vector of user-selected Abundance columns
    protein_df      = NULL,    # raw protein data.frame (all CSV cols)
    zscore_df       = NULL,    # long-form data.frame: Accession, GeneName, Condition, ZScore, Cluster
    cluster_means   = NULL,    # long-form data.frame: Cluster, Condition, MeanZ
    cluster_info    = NULL,    # data.frame: Cluster, n_proteins
    n_clusters      = 6L,      # number of clusters for fuzzy c-means
    render_count    = 0L,      # counter to force re-render
    highlight_gene  = NULL,    # gene to highlight on the plot
    col_labels      = NULL,    # short labels for conditions (derived from column names)
    zscore_wide     = NULL     # wide-form Z-score matrix (rows=proteins, cols=conditions)
  )
  methodsTextServer("trd_methods", provenance = rv$provenance, section = "trends")

  # ===========================================================================
  # React to shared CSV data — auto-select recommended Abundance columns
  # ===========================================================================
  observeEvent(rv$ss_pd_protein_headers, {
    req(rv$ss_pd_protein_headers)
    hdr <- rv$ss_pd_protein_headers
    abundance_cols <- grep("^Abundance:", hdr, value = TRUE)
    abundance_cols <- abundance_cols[!grepl("Normalized|Count|Grouped", abundance_cols)]
    trd$selected_cols <- if (length(abundance_cols) >= 2) abundance_cols else NULL
    trd$zscore_df     <- NULL
    trd$cluster_means <- NULL
    trd$cluster_info  <- NULL
    trd$protein_df    <- NULL
    trd$highlight_gene <- NULL
    logger("Trends: Shared CSV data updated, columns auto-selected")
  })

  # ===========================================================================
  # File status banner — shown after CSV upload
  # ===========================================================================
  output$trd_file_status <- renderUI({
    hdr <- rv$ss_pd_protein_headers
    fname <- rv$ss_pd_filename
    n_prot <- length(rv$ss_pd_protein_line_idx)
    if (is.null(hdr)) return(NULL)
    tags$div(class = "alert alert-success py-2 px-3 small mt-1",
      icon("circle-check"), " ",
      tags$b(fname), " — ",
      paste0(n_prot, " proteins, ", length(hdr), " columns loaded.")
    )
  })

  # ===========================================================================
  # Column selection label
  # ===========================================================================
  output$trd_col_label <- renderUI({
    sel <- trd$selected_cols
    if (is.null(sel) || length(sel) == 0) {
      tags$p(class = "text-muted small", "None selected")
    } else {
      tags$p(class = "small",
             tags$span(class = "badge bg-info", paste0(length(sel), " columns")),
             " selected")
    }
  })

  # ===========================================================================
  # Multi-column picker modal (checkboxGroupInput, same pattern as imputation)
  # ===========================================================================
  observeEvent(input$trd_pick_cols, {
    req(rv$ss_pd_protein_headers)
    all_cols <- rv$ss_pd_protein_headers
    raw_abund_cols <- grep("^Abundance", all_cols, ignore.case = TRUE, value = TRUE)
    raw_abund_cols <- raw_abund_cols[!grepl("ratio|p-value|variability|normalized|normalised|count|grouped",
                                            raw_abund_cols, ignore.case = TRUE)]
    other_cols <- setdiff(all_cols, raw_abund_cols)
    current_sel <- trd$selected_cols

    showModal(modalDialog(
      title = tags$span(icon("list-check"), " Select Abundance Columns (Conditions)"),
      tagList(
        helpText(class = "small",
          "Select the ", tags$b("raw Abundance"), " columns representing your ",
          "conditions or timepoints. Each column will become one point on the X-axis."),
        tags$style(HTML(
          "#trd_col_checks .checkbox label { display: block; padding: 3px 8px;",
          " border-radius: 4px; margin-bottom: 1px; cursor: pointer; font-size: 0.88em; }",
          "#trd_col_checks .checkbox label:hover { background-color: #e8f6f3; }",
          ".trd-scroll { max-height: 350px; overflow-y: auto;",
          " border: 1px solid #eee; padding: 6px; border-radius: 4px; }"
        )),
        if (length(raw_abund_cols) > 0) tagList(
          tags$div(
            style = "background: #f0faf7; border-left: 3px solid #18bc9c; padding: 8px 10px; margin-bottom: 10px; border-radius: 4px;",
            tags$strong(tags$span(style = "color: #18bc9c;",
                                  "\u2b50 Recommended (Raw Abundance columns)"))
          )
        ),
        tags$div(class = "trd-scroll",
          checkboxGroupInput("trd_col_checks", label = NULL,
                             choices  = setNames(c(raw_abund_cols, other_cols),
                                                  c(raw_abund_cols, other_cols)),
                             selected = if (!is.null(current_sel))
                               intersect(current_sel, c(raw_abund_cols, other_cols))
                             else raw_abund_cols)
        ),
        tags$div(class = "d-flex gap-2 mt-2",
          actionButton("trd_col_select_all", "Select All Recommended",
                       class = "btn-outline-success btn-sm", icon = icon("check-double")),
          actionButton("trd_col_clear_all", "Clear All",
                       class = "btn-outline-secondary btn-sm", icon = icon("eraser"))
        )
      ),
      size = "l", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("trd_col_confirm", "Confirm Selection", class = "btn-success")
      )
    ))
  })

  # Select All / Clear All buttons inside the modal
  observeEvent(input$trd_col_select_all, {
    all_cols <- rv$ss_pd_protein_headers
    raw_abund_cols <- grep("^Abundance", all_cols, ignore.case = TRUE, value = TRUE)
    raw_abund_cols <- raw_abund_cols[!grepl("ratio|p-value|variability|normalized|normalised|count|grouped",
                                            raw_abund_cols, ignore.case = TRUE)]
    updateCheckboxGroupInput(session, "trd_col_checks", selected = raw_abund_cols)
  })

  observeEvent(input$trd_col_clear_all, {
    updateCheckboxGroupInput(session, "trd_col_checks", selected = character(0))
  })

  observeEvent(input$trd_col_confirm, {
    sel <- input$trd_col_checks
    if (length(sel) < 2) {
      showNotification("Please select at least 2 columns.", type = "warning")
      return()
    }
    trd$selected_cols <- sel
    logger(paste0("Trends: Selected ", length(sel), " abundance columns"))
    removeModal()
  })

  # ===========================================================================
  # GENERATE — parse, log2, z-score, fuzzy c-means
  # ===========================================================================
  observeEvent(input$trd_plot_btn, {
    req(rv$ss_pd_protein_headers, rv$ss_pd_raw_lines, rv$ss_pd_protein_line_idx)

    sel_cols <- trd$selected_cols
    if (is.null(sel_cols) || length(sel_cols) < 2) {
      showNotification("Select at least 2 Abundance columns first.", type = "warning")
      return()
    }

    k <- as.integer(if (!is.null(input$trd_n_clusters)) input$trd_n_clusters else 6L)
    if (k < 2) k <- 2L
    if (k > 20) k <- 20L
    trd$n_clusters <- k

    logger(paste0("Trends: Generating with k=", k, " clusters, ",
                  length(sel_cols), " conditions"))

    # Reset viewing state
    trd$highlight_gene <- NULL
    updateTextInput(session, "trd_gene_search", value = "")

    withProgress(message = "Building Trends Analysis...", value = 0, {

      # ── Parse protein data ──
      header_line   <- rv$ss_pd_raw_lines[1]
      protein_lines <- rv$ss_pd_raw_lines[rv$ss_pd_protein_line_idx]
      text_block    <- paste(c(header_line, protein_lines), collapse = "\n")

      incProgress(0.2, detail = "Parsing protein data...")
      protein_df <- tryCatch({
        read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, na.strings = c("", "NA"))
      }, error = function(e) {
        logger(paste0("!!! [Trends] CSV parse error: ", e$message))
        showNotification(paste("CSV parse error:", e$message), type = "error")
        return(NULL)
      })
      if (is.null(protein_df)) return()

      # ── Extract gene names ──
      accessions <- protein_df[["Accession"]]
      if (!isTRUE(rv$ss_pd_has_gn) && !is.null(rv$ss_pd_gene_col) &&
          rv$ss_pd_gene_col %in% colnames(protein_df)) {
        gene_names <- as.character(protein_df[[rv$ss_pd_gene_col]])
        gene_names[is.na(gene_names)] <- ""
      } else {
        descriptions <- protein_df[["Description"]]
        gene_names   <- stringr::str_extract(descriptions, "GN=([^ ]+)")
        gene_names   <- stringr::str_replace(gene_names, "^GN=", "")
        gene_names[is.na(gene_names)] <- ""
      }

      incProgress(0.4, detail = "Log2 transform & Z-score...")

      # ── Extract abundance matrix ──
      abund_mat <- as.data.frame(lapply(sel_cols, function(col) {
        as.numeric(protein_df[[col]])
      }))
      colnames(abund_mat) <- sel_cols

      # ── Log2 transform ──
      log2_mat <- log2(abund_mat)

      # ── Remove rows with any NA/Inf ──
      valid_rows <- apply(log2_mat, 1, function(r) all(is.finite(r)))
      log2_mat   <- log2_mat[valid_rows, ]
      accessions <- accessions[valid_rows]
      gene_names <- gene_names[valid_rows]
      protein_df_valid <- protein_df[valid_rows, ]

      if (nrow(log2_mat) < 2) {
        showNotification("Fewer than 2 proteins — cannot cluster.", type = "error")
        return()
      }
      if (nrow(log2_mat) <= k) {
        k <- max(2L, nrow(log2_mat) - 1L)
        showNotification(
          paste0("Too few proteins for requested k. Auto-reduced to k=", k, "."),
          type = "warning", duration = 6)
      }

      # ── Row-wise Z-score (with zero-variance purge) ──
      zs_result  <- rowwise_zscore(log2_mat, sel_cols)
      zscore_mat <- zs_result$z
      kept_idx   <- zs_result$kept_idx
      accessions       <- accessions[kept_idx]
      gene_names       <- gene_names[kept_idx]
      protein_df_valid <- protein_df_valid[kept_idx, ]

      if (nrow(zscore_mat) < 2) {
        showNotification("Fewer than 2 proteins with non-zero variance — cannot cluster.", type = "error")
        return()
      }
      if (nrow(zscore_mat) <= k) {
        k <- max(2L, nrow(zscore_mat) - 1L)
        showNotification(
          paste0("Too few proteins after variance filter. Auto-reduced to k=", k, "."),
          type = "warning", duration = 6)
      }

      incProgress(0.6, detail = paste0("Fuzzy c-means clustering (k=", k, ")..."))

      # ── Fuzzy c-means ──
      if (!requireNamespace("e1071", quietly = TRUE)) {
        showNotification("Package 'e1071' is required for fuzzy c-means clustering but is not installed.",
                         type = "error", duration = 10)
        return()
      }
      N <- nrow(zscore_mat); D <- ncol(zscore_mat)
      m_est <- 1 + (1418/N + 22.05) * D^(-2) +
               (12.33/N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
      memberships <- NULL
      clusters <- NULL
      tryCatch({
        cl <- e1071::cmeans(zscore_mat, centers = k, m = m_est,
                            method = "cmeans", dist = "euclidean")
        clusters <- cl$cluster
        memberships <- cl$membership
        logger(paste0("Trends: Fuzzy c-means complete (m=", round(m_est, 2), ", k=", k, ")"))
      }, error = function(e) {
        logger(paste0("Trends: cmeans failed: ", conditionMessage(e)))
        showNotification(paste0("Fuzzy c-means clustering failed: ", conditionMessage(e)),
                         type = "error", duration = 10)
      })
      if (is.null(clusters)) return()

      # ── Short condition labels (strip "Abundance: F..: Sample, n/a, " prefix) ──
      col_labels <- sub("^Abundance:\\s*F\\d+:\\s*Sample,\\s*n/a,\\s*", "", sel_cols)
      col_labels <- sub("^Abundance:\\s*", "", col_labels)
      col_labels <- trimws(col_labels)
      # If labels are still very long, truncate
      if (any(nchar(col_labels) > 30)) {
        col_labels <- substr(col_labels, 1, 30)
      }

      # ── Build long-form data + cluster sorting (shared helper) ──
      result <- build_cluster_long_df(zscore_mat, accessions, gene_names,
                                       col_labels, clusters, memberships)
      long_df       <- result$long_df
      cluster_means <- result$cluster_means
      cluster_info  <- result$cluster_info
      cluster_map   <- result$cluster_map

      # ── Wide-form Z-score for CSV export ──
      zscore_wide <- as.data.frame(zscore_mat)
      colnames(zscore_wide) <- paste0("ZScore_", col_labels)
      zscore_wide$Accession <- accessions
      zscore_wide$GeneName  <- gene_names
      zscore_wide$Cluster   <- clusters
      # Re-map clusters consistently
      zscore_wide$Cluster   <- cluster_map[as.character(zscore_wide$Cluster)]

      trd$zscore_df     <- long_df
      trd$cluster_means <- cluster_means
      trd$cluster_info  <- cluster_info
      trd$protein_df    <- protein_df_valid
      trd$col_labels    <- col_labels
      trd$zscore_wide   <- zscore_wide

      # --- Provenance write ---
      rv$provenance$trends <- list(
        n_proteins_input     = as.numeric(nrow(protein_df)),
        n_proteins_clustered = as.numeric(nrow(zscore_mat)),
        n_clusters           = as.numeric(k),
        clustering_method    = "fuzzy c-means (e1071::cmeans)",
        distance_metric      = "euclidean",
        fuzzifier            = round(m_est, 3),
        zscore_transform     = TRUE,
        condition_names      = as.character(col_labels),
        software_used        = "e1071"
      )

      n_proteins <- nrow(zscore_mat)
      incProgress(1, detail = "Done!")
      logger(paste0("Trends: Analysis complete - ", n_proteins,
                    " proteins clustered into ", k, " groups"))
    })
  })

  # ===========================================================================
  # Highlight control
  # ===========================================================================
  observeEvent(input$trd_gene_search_btn, {
    gene <- trimws(input$trd_gene_search)
    trd$highlight_gene <- if (nchar(gene) == 0) NULL else gene
    trd$render_count <- isolate(trd$render_count) + 1L
  })

  # ===========================================================================
  # Profile Plot — direct plotly rendering with subplots (faceted by cluster)
  # ===========================================================================
  output$trd_profile_plotly <- renderPlotly({
    req(trd$zscore_df, trd$cluster_means)
    trd$render_count   # force-refresh dependency

    render_trend_subplots(
      long_df       = trd$zscore_df,
      cluster_means = trd$cluster_means,
      cluster_info  = trd$cluster_info,
      col_labels    = trd$col_labels,
      k             = trd$n_clusters,
      highlight     = trd$highlight_gene,
      exact_match   = isTRUE(input$trd_exact_match),
      title         = paste0("Protein Trend Profiles (k = ", trd$n_clusters, ")"),
      method        = "fcm"
    )
  })

  # ===========================================================================
  # Cluster membership table — renderDT
  # ===========================================================================
  output$trd_cluster_table <- renderDT({
    req(trd$zscore_df)
    # One row per protein
    long_df <- trd$zscore_df
    prot_df <- long_df[!duplicated(long_df$Accession), c("Accession", "GeneName", "Cluster")]
    prot_df <- prot_df[order(prot_df$Cluster, prot_df$GeneName), ]
    rownames(prot_df) <- NULL
    colnames(prot_df) <- c("Accession", "Gene", "Cluster")

    datatable(prot_df,
      options = list(
        pageLength = 25,
        scrollY = "400px",
        dom = "frtip",
        order = list(list(2, "asc"))
      ),
      filter = "top",
      rownames = FALSE,
      class = "compact stripe hover"
    )
  })

  # ===========================================================================
  # CSV Download
  # ===========================================================================
  output$trd_download_csv <- downloadHandler(
    filename = function() {
      paste0("TrendsAnalysis_k", trd$n_clusters, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(trd$zscore_wide, trd$protein_df, trd$selected_cols)

      zw <- trd$zscore_wide
      pdf <- trd$protein_df
      sel_cols <- trd$selected_cols

      # Build export: Cluster | Accession | Gene | raw abundance cols | Z-score cols
      export_df <- data.frame(
        Cluster   = zw$Cluster,
        Accession = zw$Accession,
        Gene      = zw$GeneName,
        stringsAsFactors = FALSE
      )

      # Raw abundance columns
      for (col in sel_cols) {
        export_df[[col]] <- as.numeric(pdf[[col]])
      }

      # Z-score columns
      zcols <- grep("^ZScore_", colnames(zw), value = TRUE)
      for (col in zcols) {
        export_df[[col]] <- zw[[col]]
      }

      # Sort by cluster then gene
      export_df <- export_df[order(export_df$Cluster, export_df$Gene), ]

      write.csv(export_df, file, row.names = FALSE)
      logger(paste0("Trends: CSV exported — ", nrow(export_df), " proteins"))
    }
  )

  # ===========================================================================
  # Status summary
  # ===========================================================================
  output$trd_status <- renderPrint({
    if (is.null(trd$zscore_df)) {
      cat("Load data and generate the analysis to see results.\n")
    } else {
      info <- trd$cluster_info
      info <- info[order(info$Cluster), ]
      n_total <- sum(info$n_proteins)
      cat(paste0("Total proteins clustered: ", n_total, "\n"))
      cat(paste0("Number of clusters (k): ", trd$n_clusters, "\n"))
      cat(paste0("Conditions: ", length(trd$col_labels), "\n\n"))
      for (i in seq_len(nrow(info))) {
        cat(paste0("  Cluster ", info$Cluster[i], ": ",
                    info$n_proteins[i], " proteins\n"))
      }
      hl <- trd$highlight_gene
      if (!is.null(hl) && nchar(hl) > 0) {
        long_df <- trd$zscore_df
        matches <- long_df[grepl(hl, long_df$GeneName, ignore.case = TRUE), ]
        unique_matches <- unique(matches[, c("Accession", "GeneName", "Cluster")])
        if (nrow(unique_matches) > 0) {
          cat(paste0("\nHighlighted '", hl, "':\n"))
          for (j in seq_len(nrow(unique_matches))) {
            cat(paste0("  ", unique_matches$GeneName[j],
                        " (", unique_matches$Accession[j],
                        ") -> Cluster ", unique_matches$Cluster[j], "\n"))
          }
        } else {
          cat(paste0("\nHighlighted '", hl, "': no matches\n"))
        }
      }
    }
  })
}
