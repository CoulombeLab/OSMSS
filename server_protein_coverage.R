# server_protein_coverage.R =====================================================
# Server logic for Proteome Discoverer > Protein Coverage module.
# Visualises peptide coverage across a protein sequence using the nested
# PD CSV format. Peptides are extracted from the interleaved peptide rows
# and mapped onto the protein using "Positions in Master Proteins".
#
# - Unique peptides (mapped to one protein only) are shown in GREEN.
# - Shared peptides (mapped to multiple proteins) are shown in MAGENTA.
# - Hovering shows position, sequence, and all mapped accessions.
# - Protein length is obtained from UniProt (drawProteins).
# - Coverage % is either taken from the CSV or computed from peptide spans.
#
# Called from server.R as:
#   server_protein_coverage(input, output, session, rv, logger)
# ==============================================================================

server_protein_coverage <- function(input, output, session, rv, logger) {

  cov_rv <- reactiveValues(
    plot_data     = NULL,   # data.frame of peptide coverage records
    protein_len   = NULL,   # integer: protein length
    coverage_pct  = NULL,   # numeric: coverage %
    accession     = NULL,   # target accession
    domains_df    = NULL    # domain annotation from UniProt
  )
  register_experiment_reset(rv, cov_rv)

  # ===========================================================================
  # Domain overlay selector (populated after domains are fetched)
  # ===========================================================================
  output$cov_domain_selector <- renderUI({
    domains <- cov_rv$domains_df
    if (is.null(domains) || nrow(domains) == 0) return(NULL)
    domain_names <- domains$description
    tagList(
      hr(),
      h5(icon("layer-group"), " Protein Domains"),
      helpText("Select domains to overlay on the coverage map:"),
      div(class = "d-flex gap-2 mb-2",
        actionButton("cov_dom_all", "Select All",
                     class = "btn-sm btn-outline-primary"),
        actionButton("cov_dom_none", "Deselect All",
                     class = "btn-sm btn-outline-secondary")
      ),
      checkboxGroupInput("cov_selected_domains", NULL,
                         choices = domain_names, selected = domain_names)
    )
  })

  observeEvent(input$cov_dom_all, {
    domains <- cov_rv$domains_df
    if (!is.null(domains) && nrow(domains) > 0) {
      updateCheckboxGroupInput(session, "cov_selected_domains",
                               selected = domains$description)
    }
  })

  observeEvent(input$cov_dom_none, {
    updateCheckboxGroupInput(session, "cov_selected_domains",
                             selected = character(0))
  })

  # ===========================================================================
  # Column picker for "Positions in Master Proteins" (peptide-level column)
  # Primary source: peptide_table stored in cache (works after import without CSV)
  # Fallback: parse header row from raw CSV lines
  # ===========================================================================
  cov_pep_headers <- reactive({
    pt <- rv$lf_master_cache$peptide_table
    if (!is.null(pt) && ncol(pt) > 0) return(colnames(pt))
    req(rv$pd_raw_lines)
    pep_hdr_idx <- grep("^,.*Annotated Sequence", rv$pd_raw_lines)
    if (length(pep_hdr_idx) == 0) return(character(0))
    hdr_line <- rv$pd_raw_lines[pep_hdr_idx[1]]
    cols <- as.character(
      read.csv(text = hdr_line, header = FALSE, stringsAsFactors = FALSE,
               check.names = FALSE)
    )
    cols <- cols[!is.na(cols) & nchar(trimws(cols)) > 0]
    cols
  })

  # Auto-select "Positions in Master Proteins" column on data load
  observeEvent(cov_pep_headers(), {
    ph <- cov_pep_headers()
    if (length(ph) > 0) {
      pos_col <- grep("Positions in Master Proteins", ph, value = TRUE,
                       ignore.case = TRUE)[1]
      if (!is.na(pos_col)) {
        updateSelectInput(session, "cov_pos_col",
                          choices = ph, selected = pos_col)
      } else {
        updateSelectInput(session, "cov_pos_col", choices = ph)
      }
    }
  })

  # ===========================================================================
  # Generate coverage plot
  # ===========================================================================
  observeEvent(input$cov_generate_btn, {
    has_cache_pt <- !is.null(rv$lf_master_cache$peptide_table)
    has_raw      <- !is.null(rv$pd_raw_lines)
    if (!has_cache_pt && !has_raw) {
      showNotification(
        "No peptide data available. Upload CSV or run summarization with peptide data first.",
        type = "error", duration = 8)
      return()
    }

    accession <- trimws(input$cov_accession)
    pos_col   <- input$cov_pos_col

    if (nchar(accession) == 0) {
      showNotification("Please enter a protein accession.", type = "warning")
      return()
    }
    if (is.null(pos_col) || pos_col == "") {
      showNotification("Please select the Positions column.", type = "warning")
      return()
    }

    logger(paste0("Protein Coverage: Generating for ", accession,
                  " using column '", pos_col, "'"))

    withProgress(message = "Building protein coverage...", value = 0, {

      # ------- 1. Resolve peptide data frame -------
      incProgress(0.1, detail = "Parsing peptide data...")

      if (has_cache_pt) {
        pep_df <- rv$lf_master_cache$peptide_table
        logger(paste0("Protein Coverage: Using cached peptide_table (",
                      nrow(pep_df), " rows)"))
      } else {
        all_lines     <- rv$pd_raw_lines
        pep_hdr_idx   <- grep("^,.*Annotated Sequence", all_lines)
        all_nested_idx <- grep("^,", all_lines)
        pep_data_idx   <- setdiff(all_nested_idx, pep_hdr_idx)

        if (length(pep_data_idx) < 1) {
          showNotification("No peptide-level data found in this CSV.", type = "error")
          return()
        }

        pep_block <- paste(c(all_lines[pep_hdr_idx[1]], all_lines[pep_data_idx]),
                           collapse = "\n")
        pep_df <- tryCatch(
          read.csv(text = pep_block, header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = c("", "NA")),
          error = function(e) {
            logger(paste0("!!! [Coverage] Peptide CSV parse error: ", e$message))
            NULL
          }
        )

        if (is.null(pep_df) || nrow(pep_df) == 0) {
          showNotification("Failed to parse peptide data.", type = "error")
          return()
        }
      }

      # ------- 2. Find the positions column -------
      if (!(pos_col %in% colnames(pep_df))) {
        showNotification(paste0("Column '", pos_col, "' not found in peptide data."),
                         type = "error")
        return()
      }

      incProgress(0.3, detail = "Extracting peptide positions...")

      # Also find the Annotated Sequence column
      seq_col <- grep("Annotated Sequence", colnames(pep_df), value = TRUE)[1]
      master_col <- grep("Master Protein Accessions", colnames(pep_df),
                          value = TRUE)[1]
      psm_col <- grep("^#\\s*PSMs$|^X\\.\\.PSMs$", colnames(pep_df), value = TRUE)[1]

      # ------- 3. Filter for target accession and extract position ranges -------
      pos_data <- pep_df[[pos_col]]
      pos_data[is.na(pos_data)] <- ""

      # Each entry can be like "P60709 [40-50]; Q562R1 [41-51]; P63261 [40-50]"
      # Unique peptides may omit the accession prefix: "[40-50]"
      # Search the positions column first; fall back to Master Protein Accessions.
      has_target_pos <- grepl(accession, pos_data, fixed = TRUE)
      has_target_master <- logical(nrow(pep_df))
      if (!is.na(master_col) && master_col %in% colnames(pep_df)) {
        master_data <- pep_df[[master_col]]
        master_data[is.na(master_data)] <- ""
        has_target_master <- grepl(accession, master_data, fixed = TRUE)
      }
      has_target <- has_target_pos | has_target_master
      target_pep <- pep_df[has_target, , drop = FALSE]

      if (nrow(target_pep) == 0) {
        showNotification(paste0("No peptides found for accession '", accession, "'."),
                         type = "error")
        return()
      }

      logger(paste0("Protein Coverage: Found ", nrow(target_pep),
                    " peptide rows for ", accession))

      # ------- 3b. Identify peptide abundance columns per condition -------
      pep_abund_cols <- grep("^Abundance", colnames(pep_df), value = TRUE)
      pep_abund_cols <- pep_abund_cols[!grepl(
        "Ratio|P-Value|Variability|Normalized|Normalised|Count|Grouped",
        pep_abund_cols, ignore.case = TRUE)]

      pep_cond_map <- list()
      cond_groups  <- rv$lf_cond_groups
      if (length(cond_groups) > 0 && length(pep_abund_cols) > 0) {
        for (gname in names(cond_groups)) {
          matched <- intersect(cond_groups[[gname]], pep_abund_cols)
          if (length(matched) > 0) pep_cond_map[[gname]] <- matched
        }
      }
      # If no condition groups matched, fall back to raw abundance columns
      if (length(pep_cond_map) == 0 && length(pep_abund_cols) > 0) {
        short_names <- gsub("^Abundance:\\s*", "", pep_abund_cols)
        pep_cond_map <- setNames(as.list(pep_abund_cols), short_names)
      }

      # Parse each peptide's position for the target accession
      records <- list()
      for (i in seq_len(nrow(target_pep))) {
        pos_str <- as.character(target_pep[[pos_col]][i])
        seq_str <- if (!is.na(seq_col)) as.character(target_pep[[seq_col]][i]) else ""
        master_str <- if (!is.na(master_col)) as.character(target_pep[[master_col]][i]) else ""
        psm_val <- if (!is.na(psm_col)) suppressWarnings(as.integer(target_pep[[psm_col]][i])) else NA_integer_
        if (is.na(psm_val)) psm_val <- 0L

        # Clean annotated sequence (remove flanking residues)
        clean_seq <- gsub("^\\[.*?\\]\\.", "", seq_str)
        clean_seq <- gsub("\\.\\[.*?\\]$", "", clean_seq)

        # Skip if pos_str is NA or empty
        if (is.na(pos_str) || nchar(pos_str) == 0) next

        # Normalize en-dash / em-dash → ASCII hyphen (PD Windows export artefact)
        pos_str <- chartr("\u2013\u2014", "--", pos_str)

        # Split position entries by semicolon
        entries <- strsplit(pos_str, ";")[[1]]
        entries <- trimws(entries)
        entries <- entries[!is.na(entries) & nchar(entries) > 0]

        # Find the entry for our target accession.
        # Unique PD peptides may store only "[N-M]" without an accession prefix;
        # in that case use all entries (row was already confirmed via master_col).
        target_entry <- grep(paste0("\\b", accession, "\\b"), entries, value = TRUE)
        if (length(target_entry) == 0) {
          # Fallback: unique-peptide format "[N-M]" with no accession prefix
          target_entry <- entries[grepl("^\\s*\\[\\s*\\d+", entries)]
        }
        all_accessions <- gsub("\\s*\\[.*\\]", "", entries)
        all_accessions <- trimws(all_accessions)
        all_accessions <- all_accessions[nchar(all_accessions) > 0]

        # Determine uniqueness: is this peptide mapped to only one protein?
        is_unique <- length(all_accessions) == 1

        # Compute abundance text for this peptide row
        abund_text <- ""
        if (length(pep_cond_map) > 0) {
          abund_parts <- vapply(names(pep_cond_map), function(gname) {
            cols <- pep_cond_map[[gname]]
            vals <- suppressWarnings(as.numeric(target_pep[i, cols, drop = TRUE]))
            if (all(is.na(vals))) {
              paste0(gname, ": N/A")
            } else {
              mean_val <- mean(vals, na.rm = TRUE)
              paste0(gname, ": ", format(round(mean_val, 1), big.mark = ","))
            }
          }, character(1))
          abund_text <- paste(abund_parts, collapse = "<br>")
        }

        for (te in target_entry) {
          # Extract [start-end]: accept optional spaces around hyphen, e.g. [40-50] or [40 - 50]
          range_match <- regmatches(te, regexpr("\\[(\\d+)\\s*-\\s*(\\d+)\\]", te))
          if (length(range_match) == 1) {
            nums <- as.integer(regmatches(range_match, gregexpr("\\d+", range_match))[[1]])
            start_pos <- nums[1]
            end_pos   <- nums[2]

            records[[length(records) + 1]] <- data.frame(
              Start         = start_pos,
              End           = end_pos,
              Sequence      = clean_seq,
              IsUnique      = is_unique,
              AllProteins   = paste(all_accessions, collapse = ", "),
              PSMs          = psm_val,
              AbundanceText = abund_text,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(records) == 0) {
        showNotification(paste0("Could not parse positions for '", accession, "'."),
                         type = "error")
        return()
      }

      plot_df <- do.call(rbind, records)
      # De-duplicate identical ranges
      plot_df <- unique(plot_df)

      incProgress(0.5, detail = "Fetching protein length from UniProt...")

      # ------- 4. Get protein length -------
      protein_len <- tryCatch({
        get_uniprot_length(accession)
      }, error = function(e) { logger(paste0("!!! [Coverage] UniProt length fetch error: ", e$message)); NA_integer_ })

      if (is.na(protein_len) || protein_len < 1) {
        # Fall back to max end position from data
        fallback <- suppressWarnings(max(plot_df$End, na.rm = TRUE))
        if (is.finite(fallback) && fallback >= 1) {
          protein_len <- as.integer(fallback)
          logger(paste0("Protein Coverage: UniProt length unavailable, using max peptide end = ",
                        protein_len))
        } else {
          showNotification(
            "Could not determine protein length from UniProt or peptide data.",
            type = "error", duration = 8)
          logger("Protein Coverage: No valid protein length from API or peptide data")
          return()
        }
      } else {
        logger(paste0("Protein Coverage: Protein length = ", protein_len, " aa"))
      }

      # ------- 5. Calculate coverage % -------
      incProgress(0.6, detail = "Computing coverage...")

      # Try to get coverage from the CSV first
      csv_coverage <- NA_real_
      if (!is.null(rv$pd_protein_headers) && !is.null(rv$pd_raw_lines)) {
        header_line   <- rv$pd_raw_lines[1]
        protein_lines <- rv$pd_raw_lines[rv$pd_protein_line_idx]
        text_block    <- paste(c(header_line, protein_lines), collapse = "\n")
        protein_df <- tryCatch(
          read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = c("", "NA")),
          error = function(e) { logger(paste0("!!! [Coverage] Protein CSV re-parse error: ", e$message)); NULL }
        )
        if (!is.null(protein_df)) {
          cov_col <- grep("Coverage \\[%\\]", colnames(protein_df), value = TRUE)[1]
          acc_col <- grep("^Accession$", colnames(protein_df), value = TRUE)[1]
          if (!is.na(cov_col) && !is.na(acc_col)) {
            match_row <- which(protein_df[[acc_col]] == accession)
            if (length(match_row) > 0) {
              csv_coverage <- suppressWarnings(
                as.numeric(protein_df[[cov_col]][match_row[1]])
              )
            }
          }
        }
      }

      if (!is.na(csv_coverage)) {
        coverage_pct <- csv_coverage
      } else {
        # Calculate from peptide positions
        covered <- rep(FALSE, protein_len)
        for (r in seq_len(nrow(plot_df))) {
          s <- max(1, plot_df$Start[r])
          e <- min(protein_len, plot_df$End[r])
          if (s <= e) covered[s:e] <- TRUE
        }
        coverage_pct <- round(sum(covered) / protein_len * 100, 1)
      }

      incProgress(0.7, detail = "Fetching domain annotations...")

      # ------- 6. Get domain annotations -------
      domains_df <- tryCatch({
        fetch_uniprot_domains(accession)
      }, error = function(e) {
        logger(paste0("!!! [Coverage] Domain fetch error: ", e$message))
        data.frame(type = character(0), begin = integer(0), end = integer(0),
                   description = character(0), label = character(0))
      })

      # Store results
      cov_rv$plot_data   <- plot_df
      cov_rv$protein_len <- protein_len
      cov_rv$coverage_pct <- coverage_pct
      cov_rv$accession   <- accession
      cov_rv$domains_df  <- domains_df

      incProgress(1, detail = "Done!")

      n_unique <- sum(plot_df$IsUnique == TRUE, na.rm = TRUE)
      n_shared <- sum(plot_df$IsUnique != TRUE, na.rm = TRUE)
      logger(paste0("Protein Coverage: ", nrow(plot_df), " peptide ranges (",
                    n_unique, " unique, ", n_shared, " shared), ",
                    coverage_pct, "% coverage"))

      showNotification(
        paste0("Coverage: ", coverage_pct, "% | ",
               nrow(plot_df), " peptides (", n_unique, " unique, ",
               n_shared, " shared)"),
        type = "message", duration = 6)
    })
  })

  # ===========================================================================
  # Render coverage plot (plotly)
  # ===========================================================================
  output$cov_plot <- renderPlotly({
    req(cov_rv$plot_data, cov_rv$protein_len)

    df          <- cov_rv$plot_data
    protein_len <- cov_rv$protein_len
    accession   <- cov_rv$accession
    coverage    <- cov_rv$coverage_pct
    domains     <- cov_rv$domains_df

    # Filter domains by sidebar checkbox selection
    if (!is.null(domains) && nrow(domains) > 0) {
      sel <- input$cov_selected_domains
      if (is.null(sel)) sel <- character(0)
      domains <- domains[domains$description %in% sel, , drop = FALSE]
    }

    # ---- Layout constants ----
    protein_y     <- 0       # centre of protein bar
    protein_h     <- 0.4     # half-height of protein bar
    domain_h      <- 0.35    # half-height of domain bars
    peptide_h     <- 0.15    # half-height of peptide tracks
    peptide_gap   <- 0.38    # vertical gap between peptide rows
    peptide_start <- -1.0    # y-offset for first peptide track row

    # ---- Assign peptide tracks (pack non-overlapping into rows) ----
    df <- df[order(df$Start, df$End), ]
    # Remove rows with NA positions
    df <- df[!is.na(df$Start) & !is.na(df$End), , drop = FALSE]

    if (nrow(df) == 0) {
      return(
        plot_ly(x = 0, y = 0, type = "scatter", mode = "markers",
                marker = list(size = 0, opacity = 0),
                hoverinfo = "none", showlegend = FALSE) %>%
          layout(title = list(
            text = paste0("No valid peptide positions for ", accession),
            font = list(size = 14)
          ))
      )
    }

    df$Track <- NA_integer_
    track_ends <- numeric(0)

    for (i in seq_len(nrow(df))) {
      placed <- FALSE
      for (t in seq_along(track_ends)) {
        if (df$Start[i] > track_ends[t] + 1) {
          df$Track[i] <- t
          track_ends[t] <- df$End[i]
          placed <- TRUE
          break
        }
      }
      if (!placed) {
        track_ends <- c(track_ends, df$End[i])
        df$Track[i] <- length(track_ends)
      }
    }

    # Compute y positions for each track
    df$y_center <- peptide_start - (df$Track - 1) * peptide_gap

    # ---- Build shapes ----
    shapes <- list()

    # Protein bar (full length)
    shapes[[length(shapes) + 1]] <- list(
      type = "rect",
      x0 = 1, x1 = protein_len,
      y0 = protein_y - protein_h, y1 = protein_y + protein_h,
      fillcolor = "rgba(189, 195, 199, 0.7)",
      line = list(color = "#7f8c8d", width = 1.5)
    )

    # Domain annotations
    if (!is.null(domains) && nrow(domains) > 0) {
      dom_colors <- DOMAIN_COLORS
      for (d in seq_len(nrow(domains))) {
        col_idx <- ((d - 1) %% length(dom_colors)) + 1
        shapes[[length(shapes) + 1]] <- list(
          type = "rect",
          x0 = domains$begin[d], x1 = domains$end[d],
          y0 = protein_y - domain_h, y1 = protein_y + domain_h,
          fillcolor = paste0(dom_colors[col_idx], "CC"),
          line = list(color = "#333", width = 1)
        )
      }
    }

    # Peptide rectangles
    for (i in seq_len(nrow(df))) {
      fill_col <- if (isTRUE(df$IsUnique[i])) "rgba(39, 174, 96, 0.75)"
                  else "rgba(192, 57, 183, 0.75)"
      border_col <- if (isTRUE(df$IsUnique[i])) "#1a9641" else "#8e24aa"

      shapes[[length(shapes) + 1]] <- list(
        type = "rect",
        x0 = df$Start[i], x1 = df$End[i],
        y0 = df$y_center[i] - peptide_h,
        y1 = df$y_center[i] + peptide_h,
        fillcolor = fill_col,
        line = list(color = border_col, width = 1)
      )
    }

    # ---- Build hover traces (invisible scatter for hover) ----
    psm_label <- if ("PSMs" %in% colnames(df)) df$PSMs else rep(NA_integer_, nrow(df))
    abund_label <- if ("AbundanceText" %in% colnames(df)) df$AbundanceText else rep("", nrow(df))
    df$hover <- paste0(
      "<b>Position:</b> [", df$Start, "-", df$End, "]",
      "<br><b>Sequence:</b> ", df$Sequence,
      "<br><b># PSMs:</b> ", ifelse(is.na(psm_label), "?", psm_label),
      "<br><b>Proteins:</b> ", df$AllProteins,
      "<br><b>Type:</b> ", ifelse(df$IsUnique == TRUE & !is.na(df$IsUnique), "Unique", "Shared"),
      ifelse(nchar(abund_label) > 0,
             paste0("<br><br><b>Abundance:</b><br>", abund_label),
             "")
    )

    # Midpoint x for hover targets
    df$x_mid <- (df$Start + df$End) / 2

    # Anchor trace: invisible markers at protein endpoints to force correct axis range
    p <- plot_ly(x = c(1, protein_len), y = c(protein_y, protein_y),
                 type = "scatter", mode = "markers",
                 marker = list(size = 0.1, opacity = 0),
                 hoverinfo = "none", showlegend = FALSE)

    # Unique peptide traces (green)
    unique_df <- df[which(df$IsUnique == TRUE), , drop = FALSE]
    if (nrow(unique_df) > 0) {
      p <- add_trace(p,
        data = unique_df,
        x = ~x_mid, y = ~y_center,
        type = "scatter", mode = "markers",
        name = "Unique peptide",
        marker = list(size = 1, color = "rgba(39, 174, 96, 0)", opacity = 0),
        hoverinfo = "text", text = ~hover
      )
    }

    # Shared peptide traces (magenta)
    shared_df <- df[which(df$IsUnique != TRUE), , drop = FALSE]
    if (nrow(shared_df) > 0) {
      p <- add_trace(p,
        data = shared_df,
        x = ~x_mid, y = ~y_center,
        type = "scatter", mode = "markers",
        name = "Shared peptide",
        marker = list(size = 1, color = "rgba(192, 57, 183, 0)", opacity = 0),
        hoverinfo = "text", text = ~hover
      )
    }

    # Domain hover annotations
    if (!is.null(domains) && nrow(domains) > 0) {
      dom_hover <- paste0(
        "<b>", domains$description, "</b>",
        "<br>[", domains$begin, "-", domains$end, "]"
      )
      dom_x_mid <- (domains$begin + domains$end) / 2
      p <- add_trace(p,
        x = dom_x_mid, y = rep(protein_y, nrow(domains)),
        type = "scatter", mode = "markers",
        name = "Domains",
        marker = list(size = 1, opacity = 0),
        hoverinfo = "text", text = dom_hover,
        showlegend = FALSE
      )
    }

    # ---- Annotations ----
    annotations <- list()

    # Coverage % above the protein bar
    annotations[[1]] <- list(
      x = protein_len / 2, y = protein_y + protein_h + 0.35,
      text = paste0("<b>", accession, "</b>  |  Coverage: <b>",
                    coverage, "%</b>  |  Length: ", protein_len, " aa"),
      showarrow = FALSE,
      font = list(size = 14, color = "#2c3e50"),
      xanchor = "center"
    )

    # Domain labels
    if (!is.null(domains) && nrow(domains) > 0) {
      for (d in seq_len(nrow(domains))) {
        annotations[[length(annotations) + 1]] <- list(
          x = (domains$begin[d] + domains$end[d]) / 2,
          y = protein_y + protein_h + 0.12,
          text = domains$description[d],
          showarrow = FALSE,
          font = list(size = 8, color = "#555"),
          xanchor = "center"
        )
      }
    }

    # Legend markers (manual since shapes don't have legend)
    # Add invisible legend-only traces with valid coordinates
    p <- add_trace(p,
      x = c(-999), y = c(-999),
      type = "scatter", mode = "markers",
      name = "Unique peptide",
      marker = list(size = 12, color = "rgba(39, 174, 96, 0.75)",
                    symbol = "square"),
      showlegend = TRUE, hoverinfo = "none"
    )
    p <- add_trace(p,
      x = c(-999), y = c(-999),
      type = "scatter", mode = "markers",
      name = "Shared peptide",
      marker = list(size = 12, color = "rgba(192, 57, 183, 0.75)",
                    symbol = "square"),
      showlegend = TRUE, hoverinfo = "none"
    )

    # Y-axis range
    min_y <- min(df$y_center) - 0.5
    max_y <- protein_y + protein_h + 0.7

    p <- layout(p,
      title = list(
        text = paste0("Protein Coverage Map \u2014 ", accession),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Residue Position",
        range = c(0, protein_len + 5),
        autorange = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "",
        range = c(min_y, max_y),
        autorange = FALSE,
        showticklabels = FALSE,
        zeroline = FALSE,
        showgrid = FALSE
      ),
      shapes = shapes,
      annotations = annotations,
      legend = list(orientation = "h", y = -0.1),
      hovermode = "closest",
      margin = list(t = 60)
    )

    p
  })

  # ===========================================================================
  # Status summary
  # ===========================================================================
  output$cov_status <- renderPrint({
    if (is.null(cov_rv$plot_data)) {
      cat("Upload a Proteome Discoverer CSV, enter a protein accession,\n")
      cat("and click Generate to see the coverage map.\n")
    } else {
      df <- cov_rv$plot_data
      n_unique <- sum(df$IsUnique == TRUE, na.rm = TRUE)
      n_shared <- sum(df$IsUnique != TRUE, na.rm = TRUE)
      cat(paste0("Protein: ", cov_rv$accession, "\n"))
      cat(paste0("Length:  ", cov_rv$protein_len, " amino acids\n"))
      cat(paste0("Coverage: ", cov_rv$coverage_pct, "%\n\n"))
      cat(paste0("Peptide ranges: ", nrow(df), "\n"))
      cat(paste0("  Unique (green):  ", n_unique, "\n"))
      cat(paste0("  Shared (magenta): ", n_shared, "\n"))
    }
  })
}
