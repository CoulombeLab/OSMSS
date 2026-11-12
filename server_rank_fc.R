# server_rank_fc.R ================================================================
# Server logic for Proteome Discoverer > Single-Shot > Rank Fold Change module.
# Parses Proteome Discoverer CSV (nested table-within-table format),
# extracts protein-level data with two user-selected Abundance columns,
# log2-transforms, computes fold change & z-scores, and generates:
#   Tab 1: MA plot  (Average Abundance vs Log2 FC)
#   Tab 2: Diagonal scatter  (Condition A vs Condition B)
# Called from server.R as: server_rank_fc(input, output, session, rv, logger)
# ==============================================================================

server_rank_fc <- function(input, output, session, rv, logger) {

  # ---- Module-local reactive values ----
  rfc <- reactiveValues(
    selected_col_a    = NULL,   # user-selected Abundance column A
    selected_col_b    = NULL,   # user-selected Abundance column B
    plot_data         = NULL,   # data.frame ready for plotting
    protein_df        = NULL,   # raw protein data.frame (all CSV cols, for filters)
    highlight_gene    = NULL,   # gene name to highlight (or NULL)
    label_top_n       = 20L,    # how many top-FC proteins to label
    show_top_labels   = FALSE,  # whether to label top N
    render_count      = 0L,     # counter to force plot re-render
    filter_conditions = list(), # list of filter condition lists
    filter_counter    = 0L,     # unique ID counter for conditions
    fc_gutter         = 1,      # log2FC gutter width for diagonal plot colouring
    imputed_flags     = NULL    # logical vector: TRUE for rows with any imputed value
  )
  methodsTextServer("rfc_methods", provenance = rv$provenance, section = "rank_fc")

  # ===========================================================================
  # React to shared CSV data loaded — auto-select recommended Abundance columns
  # ===========================================================================
  observeEvent(rv$ss_pd_protein_headers, {
    req(rv$ss_pd_protein_headers)
    hdr <- rv$ss_pd_protein_headers

    # Recommend individual Abundance: columns (not Grouped, not Normalized, not Count)
    abundance_cols <- grep("^Abundance:", hdr, value = TRUE)
    abundance_cols <- abundance_cols[!grepl("Normalized|Count|Grouped", abundance_cols)]

    if (length(abundance_cols) >= 2) {
      rfc$selected_col_a <- abundance_cols[1]
      rfc$selected_col_b <- abundance_cols[2]
    } else {
      rfc$selected_col_a <- NULL
      rfc$selected_col_b <- NULL
    }
    rfc$plot_data      <- NULL
    rfc$protein_df     <- NULL
    rfc$imputed_flags  <- NULL
    rfc$highlight_gene <- NULL
    rfc$show_top_labels <- FALSE
    logger("RankFC: Shared CSV data updated, columns auto-selected")
  })

  # ===========================================================================
  # File status banner — shown after CSV upload
  # ===========================================================================
  output$rfc_file_status <- renderUI({
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
  # Column selection labels
  # ===========================================================================
  output$rfc_col_a_label <- renderUI({
    sel <- rfc$selected_col_a
    if (is.null(sel)) {
      tags$p(class = "text-muted small", "None selected")
    } else {
      tags$p(class = "small",
             tags$span(class = "badge bg-info", "A"), " ",
             tags$span(style = "word-break: break-all;", sel))
    }
  })

  output$rfc_col_b_label <- renderUI({
    sel <- rfc$selected_col_b
    if (is.null(sel)) {
      tags$p(class = "text-muted small", "None selected")
    } else {
      tags$p(class = "small",
             tags$span(class = "badge bg-warning", "B"), " ",
             tags$span(style = "word-break: break-all;", sel))
    }
  })

  # ===========================================================================
  # Reusable column-picker modal (same pattern as volcano)
  # ===========================================================================
  build_rfc_col_modal <- function(modal_title, radio_id, recommended, others,
                                  rec_label, current_selection) {
    rec_choices <- setNames(recommended, recommended)
    oth_choices <- setNames(others, others)
    other_id <- paste0(radio_id, "_other")

    # JS to deselect the other radio group when one is clicked
    deselect_js <- HTML(sprintf("
      $(document).off('change.rfc_col_%1$s').on('change.rfc_col_%1$s',
        '#%1$s input[type=radio]', function(){
          $('#%2$s input[type=radio]').prop('checked', false);
          Shiny.setInputValue('%2$s', '', {priority: 'event'});
      });
      $(document).off('change.rfc_col_%2$s').on('change.rfc_col_%2$s',
        '#%2$s input[type=radio]', function(){
          $('#%1$s input[type=radio]').prop('checked', false);
          Shiny.setInputValue('%1$s', '', {priority: 'event'});
      });
    ", radio_id, other_id))

    in_rec   <- !is.null(current_selection) && current_selection %in% recommended
    in_other <- !is.null(current_selection) && current_selection %in% others

    modal_body <- tagList(
      tags$style(HTML(paste0(
        "#", radio_id, " .radio label, #", other_id, " .radio label",
        " { display: block; padding: 5px 10px;",
        " border-radius: 4px; margin-bottom: 2px; cursor: pointer; font-size: 0.9em; }",
        "#", radio_id, " .radio label:hover, #", other_id, " .radio label:hover",
        " { background-color: #e8f6f3; }",
        ".rec-section { background: #f0faf7; border-left: 3px solid #18bc9c;",
        " padding: 8px 10px; margin-bottom: 10px; border-radius: 4px; }",
        ".other-section { max-height: 300px; overflow-y: auto;",
        " border: 1px solid #eee; padding: 6px; border-radius: 4px; }"
      ))),
      tags$script(deselect_js),
      if (length(recommended) > 0) tagList(
        tags$div(class = "rec-section",
          tags$strong(tags$span(style = "color: #18bc9c;", rec_label)),
          radioButtons(radio_id, label = NULL,
                       choices  = rec_choices,
                       selected = if (in_rec) current_selection
                                  else if (in_other) character(0)
                                  else recommended[1])
        )
      ),
      if (length(others) > 0) tagList(
        tags$p(class = "text-muted small mb-1", "Other columns:"),
        tags$div(class = "other-section",
          radioButtons(other_id, label = NULL,
                       choices  = oth_choices,
                       selected = if (in_other) current_selection
                                  else character(0))
        )
      )
    )

    modalDialog(
      title = modal_title,
      modal_body,
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton(paste0(radio_id, "_confirm"), "Confirm Selection",
                     class = "btn-success")
      )
    )
  }

  # Helper: classify columns into recommended abundance and others
  get_abundance_cols <- function() {
    hdr <- rv$ss_pd_protein_headers
    abundance_cols <- grep("^Abundance:", hdr, value = TRUE)
    abundance_cols <- abundance_cols[!grepl("Normalized|Count|Grouped", abundance_cols)]
    other_cols <- setdiff(hdr, abundance_cols)
    list(recommended = abundance_cols, others = other_cols)
  }

  # ---- Column A picker ----
  observeEvent(input$rfc_pick_col_a, {
    req(rv$ss_pd_protein_headers)
    cols <- get_abundance_cols()
    showModal(build_rfc_col_modal(
      modal_title       = "Select Abundance Column — Condition A",
      radio_id          = "rfc_cola_radio",
      recommended       = cols$recommended,
      others            = cols$others,
      rec_label         = "\u2b50 Recommended (Abundance columns)",
      current_selection = rfc$selected_col_a
    ))
  })

  observeEvent(input$rfc_cola_radio_confirm, {
    sel_rec <- input$rfc_cola_radio
    sel_oth <- input$rfc_cola_radio_other
    sel <- if (!is.null(sel_oth) && nzchar(sel_oth)) sel_oth
           else if (!is.null(sel_rec) && nzchar(sel_rec)) sel_rec
           else NULL
    if (!is.null(sel) && nzchar(sel)) {
      rfc$selected_col_a <- sel
      logger(paste0("RankFC: Column A set to: ", sel))
    }
    removeModal()
  })

  # ---- Column B picker ----
  observeEvent(input$rfc_pick_col_b, {
    req(rv$ss_pd_protein_headers)
    cols <- get_abundance_cols()
    showModal(build_rfc_col_modal(
      modal_title       = "Select Abundance Column — Condition B",
      radio_id          = "rfc_colb_radio",
      recommended       = cols$recommended,
      others            = cols$others,
      rec_label         = "\u2b50 Recommended (Abundance columns)",
      current_selection = rfc$selected_col_b
    ))
  })

  observeEvent(input$rfc_colb_radio_confirm, {
    sel_rec <- input$rfc_colb_radio
    sel_oth <- input$rfc_colb_radio_other
    sel <- if (!is.null(sel_oth) && nzchar(sel_oth)) sel_oth
           else if (!is.null(sel_rec) && nzchar(sel_rec)) sel_rec
           else NULL
    if (!is.null(sel) && nzchar(sel)) {
      rfc$selected_col_b <- sel
      logger(paste0("RankFC: Column B set to: ", sel))
    }
    removeModal()
  })

  # ===========================================================================
  # DATA FILTERS — mirroring volcano pattern (OR within group, AND across)
  # ===========================================================================

  observeEvent(input$rfc_open_filters, {
    req(rv$ss_pd_protein_headers)
    showModal(build_filter_modal("rfcf", rv$ss_pd_protein_headers, rfc$filter_conditions))
  })

  observeEvent(input$rfcf_confirm, {
    col_name  <- input$rfcf_column
    op        <- input$rfcf_operator
    val       <- trimws(input$rfcf_value)
    grp       <- input$rfcf_group
    grp_logic <- input$rfcf_group_logic
    if (is.null(grp_logic) || !grp_logic %in% c("OR", "AND")) grp_logic <- "OR"
    if (op %in% c("isNA", "notNA")) val <- ""
    if (is.null(col_name) || col_name == "") {
      showNotification("Please select a column.", type = "warning"); return()
    }
    if (!op %in% c("isNA", "notNA") && val == "") {
      showNotification("Please enter a value.", type = "warning"); return()
    }
    # Update logic for all existing conditions in the same group
    if (length(rfc$filter_conditions) > 0) {
      rfc$filter_conditions <- lapply(rfc$filter_conditions, function(f) {
        if (f$group == grp) f$logic <- grp_logic; f
      })
    }
    rfc$filter_counter <- rfc$filter_counter + 1L
    rfc$filter_conditions <- c(rfc$filter_conditions, list(list(
      id = rfc$filter_counter, group = grp, column = col_name,
      op = op, value = val, logic = grp_logic
    )))
    op_label <- switch(op, "isNA" = "is NA", "notNA" = "is not NA",
                        "contains" = "contains", "not_contains" = "not contains", op)
    logger(paste0("RankFC Filter: Added [", grp, " (", grp_logic, ")] ",
                  col_name, " ", op_label,
                  if (nzchar(val)) paste0(" ", val) else ""))
    removeModal()
  })

  observeEvent(input$rfc_clear_filters, {
    rfc$filter_conditions <- list()
    logger("RankFC Filter: All conditions cleared")
  })

  # Remove single condition
  observeEvent(lapply(seq_len(100), function(i) input[[paste0("rfcf_rm_", i)]]), {
    for (i in seq_along(rfc$filter_conditions)) {
      btn_id <- paste0("rfcf_rm_", rfc$filter_conditions[[i]]$id)
      btn_val <- input[[btn_id]]
      if (!is.null(btn_val) && btn_val > 0) {
        cond <- rfc$filter_conditions[[i]]
        attr_name <- paste0(".rfcf_rm_seen_", cond$id)
        last_seen <- isolate(rfc[[attr_name]])
        if (is.null(last_seen) || btn_val > last_seen) {
          rfc[[attr_name]] <- btn_val
          rfc$filter_conditions <- rfc$filter_conditions[
            vapply(rfc$filter_conditions, function(f) f$id != cond$id, logical(1))
          ]
          logger(paste0("RankFC Filter: Removed condition #", cond$id))
          break
        }
      }
    }
  })

  # Render active filters in sidebar
  output$rfc_active_filters <- renderUI({
    conds <- rfc$filter_conditions
    if (length(conds) == 0)
      return(tags$p(class = "text-muted small", "No active filters."))
    render_filter_badges(conds, "rfcf_rm_")
  })

  # ===========================================================================
  # GENERATE PLOT DATA
  # ===========================================================================
  observeEvent(input$rfc_plot_btn, {
    req(rv$ss_pd_protein_headers, rv$ss_pd_raw_lines, rv$ss_pd_protein_line_idx)

    col_a <- rfc$selected_col_a
    col_b <- rfc$selected_col_b
    if (is.null(col_a) || col_a == "" || is.null(col_b) || col_b == "") {
      showNotification("Please select both Abundance columns (A and B).",
                       type = "warning")
      return()
    }
    if (col_a == col_b) {
      showNotification("Columns A and B must be different.", type = "warning")
      return()
    }

    logger(paste0("RankFC: Generating plot  A='", col_a, "'  B='", col_b, "'"))

    # Reset viewing state
    rfc$highlight_gene  <- NULL
    rfc$show_top_labels <- FALSE
    updateTextInput(session, "rfc_gene_search", value = "")

    withProgress(message = "Building Rank FC analysis...", value = 0, {

      # ── Parse protein data from PD nested CSV ──
      header_line   <- rv$ss_pd_raw_lines[1]
      protein_lines <- rv$ss_pd_raw_lines[rv$ss_pd_protein_line_idx]
      text_block    <- paste(c(header_line, protein_lines), collapse = "\n")

      incProgress(0.3, detail = "Parsing protein data...")
      protein_df <- tryCatch({
        read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, na.strings = c("", "NA"))
      }, error = function(e) {
        showNotification(paste("CSV parse error:", e$message), type = "error")
        logger(paste0("RankFC: CSV parse error - ", e$message))
        return(NULL)
      })
      if (is.null(protein_df)) return()

      incProgress(0.5, detail = "Computing fold changes...")

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

      # ── Extract abundances ──
      abund_a <- suppressWarnings(as.numeric(as.character(protein_df[[col_a]])))
      abund_b <- suppressWarnings(as.numeric(as.character(protein_df[[col_b]])))

      # Treat zero and negative as missing
      abund_a[!is.na(abund_a) & abund_a <= 0] <- NA
      abund_b[!is.na(abund_b) & abund_b <= 0] <- NA

      # ── Log2 transform ──
      log2_a <- log2(abund_a)
      log2_b <- log2(abund_b)

      # ── Perseus-style imputation (optional) ──
      do_impute <- isTRUE(input$rfc_impute_toggle)
      was_imputed <- rep(FALSE, length(log2_a))

      if (do_impute) {
        set.seed(42)
        imp_downshift <- if (!is.null(input$rfc_impute_downshift)) {
          input$rfc_impute_downshift
        } else {
          1.8
        }
        imp_width <- if (!is.null(input$rfc_impute_width)) {
          input$rfc_impute_width
        } else {
          0.3
        }

        na_a <- is.na(log2_a)
        obs_a <- log2_a[!na_a]
        if (sum(na_a) > 0) {
          if (length(obs_a) >= 3) {
            mu_a <- mean(obs_a)
            sd_a <- sd(obs_a)
            if (is.na(sd_a) || sd_a < 0.1) sd_a <- 0.5
            log2_a[na_a] <- rnorm(
              sum(na_a),
              mean = mu_a - imp_downshift * sd_a,
              sd   = imp_width * sd_a
            )
            was_imputed[na_a] <- TRUE
          } else {
            logger("RankFC: Skipping imputation for column A (fewer than 3 observed values)")
          }
        }

        na_b <- is.na(log2_b)
        obs_b <- log2_b[!na_b]
        if (sum(na_b) > 0) {
          if (length(obs_b) >= 3) {
            mu_b <- mean(obs_b)
            sd_b <- sd(obs_b)
            if (is.na(sd_b) || sd_b < 0.1) sd_b <- 0.5
            log2_b[na_b] <- rnorm(
              sum(na_b),
              mean = mu_b - imp_downshift * sd_b,
              sd   = imp_width * sd_b
            )
            was_imputed[na_b] <- TRUE
          } else {
            logger("RankFC: Skipping imputation for column B (fewer than 3 observed values)")
          }
        }

        logger(paste0(
          "RankFC: Imputed ", sum(was_imputed), " proteins (",
          sum(na_a), " in A, ", sum(na_b), " in B) — downshift=",
          imp_downshift, " width=", imp_width
        ))
      }

      # ── Fold change (B / A) ──
      log2_fc <- log2_b - log2_a          # log2(B/A)
      avg_intensity <- (log2_a + log2_b) / 2   # M-A average

      # ── Z-score (row-wise across conditions) ──
      # For two conditions: Z = (value - row_mean) / row_sd
      row_mean <- (log2_a + log2_b) / 2
      row_sd   <- abs(log2_a - log2_b) / sqrt(2)  # sd of two values
      row_sd[row_sd == 0] <- NA  # avoid div by zero
      zscore_a <- (log2_a - row_mean) / row_sd
      zscore_b <- (log2_b - row_mean) / row_sd

      # ── Build plot data.frame ──
      plot_df <- data.frame(
        Accession     = accessions,
        GeneName      = gene_names,
        Abundance_A   = abund_a,
        Abundance_B   = abund_b,
        Log2_A        = log2_a,
        Log2_B        = log2_b,
        Log2FC        = log2_fc,
        AvgIntensity  = avg_intensity,
        ZScore_A      = zscore_a,
        ZScore_B      = zscore_b,
        AbsLog2FC     = abs(log2_fc),
        IsImputed     = was_imputed,
        stringsAsFactors = FALSE
      )

      # Remove rows with NA/Inf in essential columns
      valid <- !is.na(plot_df$Log2_A) & is.finite(plot_df$Log2_A) &
               !is.na(plot_df$Log2_B) & is.finite(plot_df$Log2_B) &
               !is.na(plot_df$Log2FC)  & is.finite(plot_df$Log2FC)
      plot_df <- plot_df[valid, ]

      # ── Rank by absolute FC ──
      plot_df$FC_Rank <- rank(-plot_df$AbsLog2FC, ties.method = "min")

      rfc$plot_data  <- plot_df
      rfc$protein_df <- protein_df
      rfc$imputed_flags <- plot_df$IsImputed

      # --- Provenance write ---
      rv$provenance$rank_fc <- list(
        n_proteins_input   = as.numeric(nrow(protein_df)),
        n_proteins_plotted = as.numeric(nrow(plot_df)),
        n_imputed          = as.numeric(sum(was_imputed, na.rm = TRUE)),
        condition_a        = col_a,
        condition_b        = col_b,
        fc_method          = "log2(B/A)",
        imputation_enabled = isTRUE(do_impute),
        imputation_method  = if (isTRUE(do_impute)) "Perseus-style random normal" else character(),
        imp_downshift      = if (isTRUE(do_impute)) as.numeric(imp_downshift) else numeric(),
        imp_width          = if (isTRUE(do_impute)) as.numeric(imp_width) else numeric(),
        software_used      = "plotly"
      )

      incProgress(1, detail = "Done!")
      if (do_impute) {
        logger(paste0("RankFC: Plot data ready - ", nrow(plot_df),
                      " proteins (", sum(plot_df$IsImputed, na.rm = TRUE), " imputed)"))
      } else {
        logger(paste0("RankFC: Plot data ready - ", nrow(plot_df),
                      " proteins with valid abundance values"))
      }
    })
  })

  # ===========================================================================
  # Highlight & Label controls
  # ===========================================================================
  observeEvent(input$rfc_gene_search_btn, {
    gene <- trimws(input$rfc_gene_search)
    rfc$highlight_gene <- if (nchar(gene) == 0) NULL else gene
    rfc$render_count <- isolate(rfc$render_count) + 1L
  })

  observeEvent(input$rfc_label_top_btn, {
    rfc$show_top_labels <- !isTRUE(rfc$show_top_labels)
    rfc$render_count <- isolate(rfc$render_count) + 1L
  })

  observeEvent(input$rfc_top_n, {
    n <- input$rfc_top_n
    if (!is.null(n) && is.numeric(n) && n >= 1) {
      rfc$label_top_n <- as.integer(n)
      if (isTRUE(rfc$show_top_labels)) {
        rfc$render_count <- isolate(rfc$render_count) + 1L
      }
    }
  })

  observeEvent(input$rfc_gutter_width, {
    val <- input$rfc_gutter_width
    if (!is.null(val) && is.numeric(val) && val >= 0) {
      rfc$fc_gutter <- val
      rfc$render_count <- isolate(rfc$render_count) + 1L
    }
  })

  # ===========================================================================
  # Helper: prepare data for rendering (filter + highlight + labels)
  # ===========================================================================
  prepare_rfc_data <- function() {
    req(rfc$plot_data)
    df       <- rfc$plot_data
    raw_prot <- rfc$protein_df
    rfc$render_count  # reactive dependency for force-refresh

    # Apply filters
    filter_conds <- rfc$filter_conditions
    if (length(filter_conds) > 0 && !is.null(raw_prot)) {
      keep_mask <- apply_filter_conditions(raw_prot, filter_conds, "and_or")
      keep_acc  <- raw_prot[["Accession"]][keep_mask]
      df <- df[df$Accession %in% keep_acc, ]
    }

    highlight_gene  <- rfc$highlight_gene
    show_top        <- rfc$show_top_labels
    top_n           <- rfc$label_top_n
    exact_match     <- isTRUE(input$rfc_exact_match)
    gutter          <- rfc$fc_gutter

    # Re-rank after filtering
    df$FC_Rank <- rank(-df$AbsLog2FC, ties.method = "min")

    list(df = df, highlight_gene = highlight_gene, show_top = show_top,
         top_n = top_n, exact_match = exact_match, gutter = gutter)
  }

  # ===========================================================================
  # Tab 1: MA Plot — renderPlotly
  # ===========================================================================
  output$rfc_ma_plot <- renderPlotly({
    d <- prepare_rfc_data()
    df <- d$df
    if (is.null(df) || nrow(df) == 0) return(plotly_empty())

    highlight_gene <- d$highlight_gene
    show_top       <- d$show_top
    top_n          <- d$top_n
    exact_match    <- d$exact_match
    gutter         <- d$gutter

    # ── Determine categories ──
    if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      # Highlight mode
      if (exact_match) {
        is_match <- toupper(df$GeneName) == toupper(highlight_gene)
      } else {
        is_match <- grepl(highlight_gene, df$GeneName, ignore.case = TRUE)
      }
      match_name <- paste0("Searched: ", highlight_gene)
      df$Category <- ifelse(is_match, match_name, "Other")
      cat_colors  <- c("Other" = "#d4d4d4")
      cat_colors[match_name] <- "#e74c3c"
      df$pt_size    <- ifelse(is_match, 12, 6)
      df$pt_opacity <- ifelse(is_match, 1, 0.35)
      df$pt_symbol  <- ifelse(df$IsImputed, "triangle-up", "circle")
    } else {
      # FC-based coloring
      df$Category <- ifelse(
        df$Log2FC >  gutter, "Enriched in B",
        ifelse(df$Log2FC < -gutter, "Enriched in A", "Unchanged"))
      cat_colors <- c(
        "Enriched in B" = "#e74c3c",
        "Enriched in A" = "#3498db",
        "Unchanged"     = "#d4d4d4"
      )
      df$pt_size    <- ifelse(df$Category == "Unchanged", 6, 8)
      df$pt_opacity <- ifelse(df$Category == "Unchanged", 0.4, 0.8)
      df$pt_symbol  <- ifelse(df$IsImputed, "triangle-up", "circle")
    }

    # ── Hover text ──
    col_a_short <- sub("^Abundance:\\s*", "", rfc$selected_col_a)
    col_b_short <- sub("^Abundance:\\s*", "", rfc$selected_col_b)
    df$hover <- paste0(
      "<b>Gene:</b> ", df$GeneName,
      "<br><b>Accession:</b> ", df$Accession,
      "<br><b>Log2 FC (B/A):</b> ", round(df$Log2FC, 3),
      "<br><b>Avg Log2 Intensity:</b> ", round(df$AvgIntensity, 2),
      "<br><b>", col_a_short, ":</b> ", signif(df$Abundance_A, 4),
      "<br><b>", col_b_short, ":</b> ", signif(df$Abundance_B, 4),
      "<br><b>FC Rank:</b> #", df$FC_Rank,
      ifelse(df$IsImputed, "<br><i style='color:#e67e22;'>⚠ Imputed</i>", "")
    )

    # ── Build traces per category ──
    categories <- unique(df$Category)
    p <- plot_ly()
    for (cat in categories) {
      sub <- df[df$Category == cat, ]
      p <- add_trace(p, data = sub,
        x = ~AvgIntensity, y = ~Log2FC,
        type = "scatter", mode = "markers", name = cat,
        marker = list(size = sub$pt_size, opacity = sub$pt_opacity,
                color = cat_colors[cat], symbol = sub$pt_symbol),
        hoverinfo = "text", text = ~hover)
    }

    # ── Labels: top N by FC ──
    if (show_top && !is.null(highlight_gene) == FALSE) {
      top_df <- df[df$FC_Rank <= top_n, ]
      if (nrow(top_df) > 0) {
        arrow_col <- ifelse(top_df$Log2FC > 0, "#c0392b", "#2471a3")
        p <- add_annotations(p,
          x = top_df$AvgIntensity, y = top_df$Log2FC,
          text = top_df$GeneName,
          showarrow = TRUE, arrowhead = 0, arrowsize = 0.5,
          ax = 25, ay = ifelse(top_df$Log2FC > 0, -15, 15),
          font = list(size = 9, color = arrow_col),
          arrowcolor = arrow_col)
      }
    }

    # ── Highlight labels ──
    if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      if (exact_match) {
        hits <- df[toupper(df$GeneName) == toupper(highlight_gene), ]
      } else {
        hits <- df[grepl(highlight_gene, df$GeneName, ignore.case = TRUE), ]
      }
      if (nrow(hits) > 0) {
        p <- add_annotations(p,
          x = hits$AvgIntensity, y = hits$Log2FC,
          text = hits$GeneName,
          showarrow = TRUE, arrowhead = 2, arrowsize = 1,
          ax = 30, ay = -25,
          font = list(size = 11, color = "#c0392b", family = "Arial Black"),
          arrowcolor = "#c0392b")
      }
    }

    # ── Layout ──
    p <- layout(p,
      title = list(text = "MA Plot — Single-Shot Fold Change",
                   font = list(size = 18)),
      xaxis = list(title = "Average Log2 Intensity (A)",
                   zeroline = FALSE),
      yaxis = list(title = "Log2 Fold Change (B / A)",
                   zeroline = TRUE, zerolinecolor = "#999", zerolinewidth = 1),
      shapes = list(
        list(type = "line", x0 = 0, x1 = 1, y0 =  gutter, y1 =  gutter,
             xref = "paper",
             line = list(dash = "dash", color = "#e74c3c", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, y0 = -gutter, y1 = -gutter,
             xref = "paper",
             line = list(dash = "dash", color = "#3498db", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, y0 = 0, y1 = 0,
             xref = "paper",
             line = list(dash = "dot", color = "#999", width = 1))
      ),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "closest"
    )
    p
  })

  # ===========================================================================
  # Tab 2: Diagonal Scatter Plot — renderPlotly
  # ===========================================================================
  output$rfc_diag_plot <- renderPlotly({
    d <- prepare_rfc_data()
    df <- d$df
    if (is.null(df) || nrow(df) == 0) return(plotly_empty())

    highlight_gene <- d$highlight_gene
    show_top       <- d$show_top
    top_n          <- d$top_n
    exact_match    <- d$exact_match
    gutter         <- d$gutter

    # ── Determine categories ──
    # Points outside the diagonal gutter are coloured; inside are grey
    if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      if (exact_match) {
        is_match <- toupper(df$GeneName) == toupper(highlight_gene)
      } else {
        is_match <- grepl(highlight_gene, df$GeneName, ignore.case = TRUE)
      }
      match_name <- paste0("Searched: ", highlight_gene)
      df$Category <- ifelse(is_match, match_name, "Other")
      cat_colors  <- c("Other" = "#d4d4d4")
      cat_colors[match_name] <- "#e74c3c"
      df$pt_size    <- ifelse(is_match, 12, 6)
      df$pt_opacity <- ifelse(is_match, 1, 0.35)
      df$pt_symbol  <- ifelse(df$IsImputed, "triangle-up", "circle")
    } else {
      # Distance from diagonal: Log2FC = Log2_B - Log2_A
      df$Category <- ifelse(
        df$Log2FC >  gutter, "Enriched in B",
        ifelse(df$Log2FC < -gutter, "Enriched in A", "Unchanged"))
      cat_colors <- c(
        "Enriched in B" = "#e74c3c",
        "Enriched in A" = "#3498db",
        "Unchanged"     = "#d4d4d4"
      )
      df$pt_size    <- ifelse(df$Category == "Unchanged", 6, 8)
      df$pt_opacity <- ifelse(df$Category == "Unchanged", 0.4, 0.8)
      df$pt_symbol  <- ifelse(df$IsImputed, "triangle-up", "circle")
    }

    # ── Hover text ──
    col_a_short <- sub("^Abundance:\\s*", "", rfc$selected_col_a)
    col_b_short <- sub("^Abundance:\\s*", "", rfc$selected_col_b)
    df$hover <- paste0(
      "<b>Gene:</b> ", df$GeneName,
      "<br><b>Accession:</b> ", df$Accession,
      "<br><b>Log2 A (", col_a_short, "):</b> ", round(df$Log2_A, 2),
      "<br><b>Log2 B (", col_b_short, "):</b> ", round(df$Log2_B, 2),
      "<br><b>Log2 FC:</b> ", round(df$Log2FC, 3),
      "<br><b>FC Rank:</b> #", df$FC_Rank,
      ifelse(df$IsImputed, "<br><i style='color:#e67e22;'>⚠ Imputed</i>", "")
    )

    # ── Build traces per category ──
    categories <- unique(df$Category)
    p <- plot_ly()
    for (cat in categories) {
      sub <- df[df$Category == cat, ]
      p <- add_trace(p, data = sub,
        x = ~Log2_A, y = ~Log2_B,
        type = "scatter", mode = "markers", name = cat,
        marker = list(size = sub$pt_size, opacity = sub$pt_opacity,
                color = cat_colors[cat], symbol = sub$pt_symbol),
        hoverinfo = "text", text = ~hover)
    }

    # ── Top N labels ──
    if (show_top && is.null(highlight_gene)) {
      top_df <- df[df$FC_Rank <= top_n, ]
      if (nrow(top_df) > 0) {
        arrow_col <- ifelse(top_df$Log2FC > 0, "#c0392b", "#2471a3")
        p <- add_annotations(p,
          x = top_df$Log2_A, y = top_df$Log2_B,
          text = top_df$GeneName,
          showarrow = TRUE, arrowhead = 0, arrowsize = 0.5,
          ax = 25, ay = -15,
          font = list(size = 9, color = arrow_col),
          arrowcolor = arrow_col)
      }
    }

    # ── Highlight labels ──
    if (!is.null(highlight_gene) && nchar(highlight_gene) > 0) {
      if (exact_match) {
        hits <- df[toupper(df$GeneName) == toupper(highlight_gene), ]
      } else {
        hits <- df[grepl(highlight_gene, df$GeneName, ignore.case = TRUE), ]
      }
      if (nrow(hits) > 0) {
        p <- add_annotations(p,
          x = hits$Log2_A, y = hits$Log2_B,
          text = hits$GeneName,
          showarrow = TRUE, arrowhead = 2, arrowsize = 1,
          ax = 30, ay = -25,
          font = list(size = 11, color = "#c0392b", family = "Arial Black"),
          arrowcolor = "#c0392b")
      }
    }

    # ── Diagonal guidelines ──
    data_range <- range(c(df$Log2_A, df$Log2_B), na.rm = TRUE)
    rng_min <- floor(data_range[1]) - 1
    rng_max <- ceiling(data_range[2]) + 1

    p <- layout(p,
      title = list(text = "Diagonal Scatter — Condition A vs B",
                   font = list(size = 18)),
      xaxis = list(title = paste0("Log2 Abundance — Condition A"),
                   range = c(rng_min, rng_max)),
      yaxis = list(title = paste0("Log2 Abundance — Condition B"),
                   scaleanchor = "x", scaleratio = 1,
                   range = c(rng_min, rng_max)),
      shapes = list(
        # y = x diagonal (identity)
        list(type = "line", x0 = rng_min, x1 = rng_max,
             y0 = rng_min, y1 = rng_max,
             line = list(dash = "dot", color = "#999", width = 1)),
        # y = x + gutter (upper threshold)
        list(type = "line", x0 = rng_min, x1 = rng_max,
             y0 = rng_min + gutter, y1 = rng_max + gutter,
             line = list(dash = "dash", color = "#e74c3c", width = 1)),
        # y = x - gutter (lower threshold)
        list(type = "line", x0 = rng_min, x1 = rng_max,
             y0 = rng_min - gutter, y1 = rng_max - gutter,
             line = list(dash = "dash", color = "#3498db", width = 1))
      ),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "closest"
    )
    p
  })

  # ===========================================================================
  # Status summary
  # ===========================================================================
  output$rfc_status <- renderPrint({
    if (is.null(rfc$plot_data)) {
      cat("Load data and generate the plot to see results.\n")
    } else {
      df <- rfc$plot_data
      gutter <- rfc$fc_gutter
      n_enriched_b <- sum(df$Log2FC >  gutter, na.rm = TRUE)
      n_enriched_a <- sum(df$Log2FC < -gutter, na.rm = TRUE)
      n_unchanged  <- sum(abs(df$Log2FC) <= gutter, na.rm = TRUE)
      n_imputed    <- sum(df$IsImputed, na.rm = TRUE)
      cat(paste0("Total proteins: ", nrow(df), "\n"))
      cat(paste0("  Enriched in B (Log2FC > ", gutter, "): ", n_enriched_b, "\n"))
      cat(paste0("  Enriched in A (Log2FC < -", gutter, "): ", n_enriched_a, "\n"))
      cat(paste0("  Unchanged:              ", n_unchanged, "\n"))
      if (n_imputed > 0) {
        cat(paste0("  Imputed values:         ", n_imputed, "\n"))
      }
      cat(paste0("\nCondition A: ", rfc$selected_col_a, "\n"))
      cat(paste0("Condition B: ", rfc$selected_col_b, "\n"))
      hl <- rfc$highlight_gene
      if (!is.null(hl) && nchar(hl) > 0) {
        n_match <- sum(grepl(hl, df$GeneName, ignore.case = TRUE))
        cat(paste0("\nHighlighted gene '", hl, "': ", n_match, " match(es)\n"))
      }
    }
  })

  # ===========================================================================
  # GSEA Export — visibility flag
  # ===========================================================================
  output$rfc_gsea_rnk_ready <- reactive({
    !is.null(rfc$plot_data) && nrow(rfc$plot_data) > 0
  })
  outputOptions(output, "rfc_gsea_rnk_ready", suspendWhenHidden = FALSE)

  # ===========================================================================
  # CSV Download — visibility flag + handler
  # ===========================================================================
  output$rfc_csv_ready <- reactive({
    !is.null(rfc$plot_data) && nrow(rfc$plot_data) > 0
  })
  outputOptions(output, "rfc_csv_ready", suspendWhenHidden = FALSE)

  output$rfc_download_csv <- downloadHandler(
    filename = function() {
      paste0("RankFC_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      prep <- prepare_rfc_data()
      df   <- prep$df

      out <- data.frame(
        Accession    = df$Accession,
        Gene         = df$GeneName,
        Abundance_A  = df$Abundance_A,
        Abundance_B  = df$Abundance_B,
        Log2_A       = round(df$Log2_A, 6),
        Log2_B       = round(df$Log2_B, 6),
        Log2FC       = round(df$Log2FC, 6),
        AbsLog2FC    = round(df$AbsLog2FC, 6),
        FC_Rank      = df$FC_Rank,
        AvgIntensity = round(df$AvgIntensity, 6),
        IsImputed    = df$IsImputed,
        stringsAsFactors = FALSE
      )
      out <- out[order(out$FC_Rank), ]

      write.csv(out, file, row.names = FALSE)
      logger(paste0("RankFC: Exported CSV with ", nrow(out), " proteins"))
    }
  )

  # ===========================================================================
  # GSEA Export — .rnk (Ranked List)
  # Format: tab-delimited, no header, GeneName \t Log2FC
  # ===========================================================================
  output$rfc_download_rnk <- downloadHandler(
    filename = function() {
      paste0("GSEA_Ranked_RankFC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rnk")
    },
    content = function(file) {
      req(rfc$plot_data)
      df <- rfc$plot_data

      # Keep only rows with valid gene names
      df <- df[!is.na(df$GeneName) & nchar(trimws(df$GeneName)) > 0, ]

      # Deduplicate gene names: keep the row with the highest absolute Log2FC
      df <- df[order(-abs(df$Log2FC)), ]
      df <- df[!duplicated(df$GeneName), ]

      # Sort descending by Log2FC (GSEA convention)
      df <- df[order(-df$Log2FC), ]

      # Write tab-delimited, no header
      write.table(
        data.frame(df$GeneName, round(df$Log2FC, 6)),
        file      = file,
        sep       = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote     = FALSE
      )

      logger(paste0("RankFC GSEA: Exported .rnk with ", nrow(df), " genes"))
    }
  )
}
