# server_ptm_cross_exp.R =======================================================
# PTM Cross-Experiment Comparison module.
#
# Accepts 2+ PTM_DE_Master_*.csv files (from PTM > Pairwise DE > Download CSV)
# and computes regulated site/gene sets across four analytical layers:
#
#   Tab 1  Adjusted PTM            — site-level, PTM_Label identifier
#   Tab 2  Absolute PTM            — site-level, PTM_Label identifier
#   Tab 3  Protein Only            — protein-level, Gene identifier (de-duped)
#   Tab 4  Adjusted (Collapsed)    — Adjusted FC/p, collapsed to Gene level;
#                                    a gene is regulated if ANY site qualifies
#
# Thresholds: |Log2FC| > 1 AND adj.p / BFDR < 0.05 (matching DE volcano).
#
# Called from server.R as: server_ptm_cross_exp(input, output, session, rv, logger)
# ==============================================================================

server_ptm_cross_exp <- function(input, output, session, rv, logger) {

  CPTM_FC_UP    <-  1
  CPTM_FC_DOWN  <- -1
  CPTM_P_THRESH <- 0.05

  # Ordered layer definitions — drives all downstream helpers
  .LAYERS <- list(
    adj = list(
      label    = "Adjusted PTM",
      id_col   = "PTM_Label",
      fc_col   = "Adjusted_log2FC",
      p_col    = "Adjusted_adj.pvalue",
      dedup    = FALSE,
      id_type  = "PTM site"
    ),
    abs = list(
      label    = "Absolute PTM",
      id_col   = "PTM_Label",
      fc_col   = "Absolute_log2FC",
      p_col    = "Absolute_adj.pvalue",
      dedup    = FALSE,
      id_type  = "PTM site"
    ),
    prot = list(
      label    = "Protein Only",
      id_col   = "Gene",
      fc_col   = "ProteinOnly_log2FC",
      p_col    = "ProteinOnly_adj.pvalue",
      dedup    = TRUE,
      id_type  = "protein"
    ),
    collprot = list(
      label    = "Adjusted \u2014 Collapsed by Protein",
      id_col   = "Gene",
      fc_col   = "Adjusted_log2FC",
      p_col    = "Adjusted_adj.pvalue",
      dedup    = TRUE,
      id_type  = "protein"
    )
  )

  # ---------------------------------------------------------------------------
  # Reactive storage
  # ---------------------------------------------------------------------------
  cptm_rv <- reactiveValues(
    files_parsed  = list(),
    parse_errors  = character(0),
    n_files       = 0L
  )

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  #' Check which layers are available in a data.frame
  available_layers <- function(df) {
    nms <- names(df)
    Filter(function(lk) {
      l <- .LAYERS[[lk]]
      all(c(l$id_col, l$fc_col, l$p_col) %in% nms)
    }, names(.LAYERS))
  }

  #' Classify items in one layer into up / down / nonsig sets.
  #' For dedup layers, collapses to unique Gene IDs;
  #' a gene is up/down if ANY row for that gene meets the threshold.
  classify_layer <- function(df, layer) {
    id_col <- layer$id_col
    fc_col <- layer$fc_col
    p_col  <- layer$p_col
    dedup  <- layer$dedup

    id  <- trimws(as.character(df[[id_col]]))
    fc  <- suppressWarnings(as.numeric(df[[fc_col]]))
    sig <- suppressWarnings(as.numeric(df[[p_col]]))

    sig[is.na(sig)]                  <- 1
    id[id == "" | id == "NA"]        <- NA_character_

    valid     <- !is.na(id) & !is.na(fc) & is.finite(fc)
    up_mask   <- valid & fc >  CPTM_FC_UP   & sig < CPTM_P_THRESH
    down_mask <- valid & fc <  CPTM_FC_DOWN & sig < CPTM_P_THRESH

    up_ids <- unique(id[up_mask])
    dn_ids <- unique(id[down_mask])

    if (dedup) {
      # For collapsed layers: total = unique valid genes,
      # nonsig = those with no site meeting either threshold
      all_valid_ids <- unique(id[valid])
      ns_ids  <- setdiff(all_valid_ids, union(up_ids, dn_ids))
      n_total <- length(all_valid_ids)
    } else {
      ns_ids  <- unique(id[valid & !up_mask & !down_mask])
      n_total <- length(unique(id[valid]))
    }

    list(up = up_ids, down = dn_ids, nonsig = ns_ids, n_total = n_total)
  }

  # ---------------------------------------------------------------------------
  # File parsing
  # ---------------------------------------------------------------------------
  observeEvent(input$cptm_files, {
    req(input$cptm_files)
    finfo <- input$cptm_files
    n     <- nrow(finfo)

    parsed <- vector("list", n)
    errors <- character(0)

    for (i in seq_len(n)) {
      fname <- finfo$name[i]
      fpath <- finfo$datapath[i]

      df <- tryCatch(
        read.csv(fpath, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )

      if (is.null(df) || nrow(df) == 0) {
        errors <- c(errors, paste0(fname, ": Failed to read or empty CSV."))
        next
      }

      avail <- available_layers(df)

      if (length(avail) == 0) {
        errors <- c(errors, paste0(
          fname, ": Not recognised as a PTM DE Master CSV. ",
          "Expected columns like Adjusted_log2FC, Absolute_log2FC, ProteinOnly_log2FC."
        ))
        next
      }

      # Pre-compute all available layers
      layer_sets <- lapply(avail, function(lk) {
        classify_layer(df, .LAYERS[[lk]])
      })
      names(layer_sets) <- avail

      # Auto-name: strip OSMSS default prefix and date
      exp_name <- tools::file_path_sans_ext(basename(fname))
      exp_name <- sub("^PTM_DE_Master_", "", exp_name)
      exp_name <- sub("_\\d{4}-\\d{2}-\\d{2}$", "", exp_name)  # strip trailing date
      exp_name <- gsub("_", " ", exp_name)

      parsed[[i]] <- list(
        fname      = fname,
        exp_name   = exp_name,
        avail_layers = avail,
        layer_sets = layer_sets
      )
    }

    valid_parsed <- Filter(Negate(is.null), parsed)

    cptm_rv$files_parsed <- valid_parsed
    cptm_rv$n_files      <- length(valid_parsed)
    cptm_rv$parse_errors <- errors

    for (e in errors) showNotification(e, type = "warning", duration = 12)

    if (length(valid_parsed) >= 2) {
      showNotification(
        paste0(length(valid_parsed), " PTM DE file(s) loaded successfully."),
        type = "message", duration = 5
      )
    } else if (length(valid_parsed) == 1) {
      showNotification(
        "Only 1 valid file — upload at least 2 to compare.",
        type = "warning", duration = 8
      )
    }
  })

  # ---------------------------------------------------------------------------
  # Experiment name inputs
  # ---------------------------------------------------------------------------
  output$cptm_exp_name_inputs <- renderUI({
    req(cptm_rv$n_files >= 2)
    tagList(lapply(seq_len(cptm_rv$n_files), function(i) {
      textInput(
        inputId     = paste0("cptm_exp_name_", i),
        label       = NULL,
        value       = cptm_rv$files_parsed[[i]]$exp_name,
        width       = "100%",
        placeholder = paste0("Experiment ", i)
      )
    }))
  })

  get_exp_names <- reactive({
    req(cptm_rv$n_files >= 2)
    vapply(seq_len(cptm_rv$n_files), function(i) {
      nm <- input[[paste0("cptm_exp_name_", i)]]
      if (is.null(nm) || nchar(trimws(nm)) == 0) paste0("Exp ", i)
      else trimws(nm)
    }, character(1))
  })

  # ---------------------------------------------------------------------------
  # Ready flag
  # ---------------------------------------------------------------------------
  output$cptm_ready <- reactive({ cptm_rv$n_files >= 2 })
  outputOptions(output, "cptm_ready", suspendWhenHidden = FALSE)

  # ---------------------------------------------------------------------------
  # Highlight search state (press-to-render)
  # ---------------------------------------------------------------------------
  cptm_hl_gene <- reactiveVal(NULL)

  observeEvent(input$cptm_highlight_btn, {
    hl <- if (is.null(input$cptm_search_input)) "" else trimws(input$cptm_search_input)
    cptm_hl_gene(hl)
  })

  observeEvent(input$cptm_clear_btn, {
    cptm_hl_gene(NULL)
    updateTextInput(session, "cptm_search_input", value = "")
  })

  output$cptm_search_result_ui <- renderUI({
    hl <- cptm_hl_gene()
    req(!is.null(hl) && nchar(hl) > 0 && cptm_rv$n_files >= 2)
    fps       <- cptm_rv$files_parsed
    exp_names <- get_exp_names()
    badges <- lapply(seq_along(fps), function(i) {
      fp      <- fps[[i]]
      in_up   <- FALSE; in_down <- FALSE
      for (lk in fp$avail_layers) {
        s <- fp$layer_sets[[lk]]
        if (any(grepl(hl, s$up,   ignore.case = TRUE))) in_up   <- TRUE
        if (any(grepl(hl, s$down, ignore.case = TRUE))) in_down <- TRUE
      }
      info <- if (in_up && in_down)  list(label = "Discordant", cls = "warning")
              else if (in_up)        list(label = "\u2191 Up",   cls = "danger")
              else if (in_down)      list(label = "\u2193 Down", cls = "primary")
              else                   list(label = "NS",          cls = "secondary")
      tags$div(class = "mb-1",
        tags$strong(exp_names[i], ": "),
        tags$span(class = paste0("badge text-bg-", info$cls), info$label)
      )
    })
    tags$div(class = "small mt-2 p-2 bg-light rounded",
             tags$strong(style = "color:#b36a00;", "\u25BA ", hl),
             tags$div(class = "mt-1", tagList(badges)))
  })

  # ---------------------------------------------------------------------------
  # Sidebar status
  # ---------------------------------------------------------------------------
  output$cptm_status_msg <- renderUI({
    if (cptm_rv$n_files == 0) {
      return(tags$p(class = "text-muted small mt-2",
                    icon("circle-info"),
                    " Upload \u2265 2 PTM DE Master CSV files to begin."))
    }
    if (cptm_rv$n_files < 2) {
      return(tags$p(class = "text-danger small mt-2",
                    icon("triangle-exclamation"),
                    " Need at least 2 valid CSV files."))
    }

    has_venn <- .pkg_installed("ggvenn")
    venn_warn <- if (!has_venn) {
      tags$div(
        class = "alert alert-warning py-2 px-3 small mt-2",
        icon("triangle-exclamation"),
        " Venn diagrams require the ",
        tags$code("ggvenn"),
        " package. Install with:",
        tags$br(),
        tags$code('install.packages("ggvenn")')
      )
    }

    # Warn about layers that are missing from some files
    all_avail <- lapply(cptm_rv$files_parsed, `[[`, "avail_layers")
    shared <- Reduce(intersect, all_avail)
    missing_note <- if (length(shared) < length(.LAYERS)) {
      missing_lks <- setdiff(names(.LAYERS), shared)
      missing_labels <- vapply(missing_lks, function(k) .LAYERS[[k]]$label, character(1))
      tags$p(class = "text-warning small mt-1",
             icon("triangle-exclamation"),
             paste0(" Some tabs may be unavailable: ",
                    paste(missing_labels, collapse = ", "), "."))
    }

    tagList(
      tags$p(class = "text-success small mt-2",
             icon("circle-check"),
             paste0(" ", cptm_rv$n_files, " files ready.")),
      missing_note,
      venn_warn
    )
  })

  # ---------------------------------------------------------------------------
  # Panel header (shown before files are loaded)
  # ---------------------------------------------------------------------------
  output$cptm_panel_header <- renderUI({
    if (cptm_rv$n_files < 2) {
      tags$div(
        class = "alert alert-info py-3 px-4",
        tags$h5(icon("atom"), " PTM Cross-Experiment Comparison"),
        tags$p(class = "mb-1",
          "Upload two or more ", tags$b("PTM DE Master CSV"), " files from the ",
          tags$b("PTM \u2192 Pairwise DE \u2192 Download CSV"), " module using the sidebar."
        ),
        tags$ul(class = "mb-1 small",
          tags$li(tags$b("Tab 1 \u2014 Adjusted PTM:"), " site-level, occupancy-corrected FC"),
          tags$li(tags$b("Tab 2 \u2014 Absolute PTM:"), " site-level, unadjusted FC"),
          tags$li(tags$b("Tab 3 \u2014 Protein Only:"), " protein-level model (unique genes)"),
          tags$li(tags$b("Tab 4 \u2014 Adjusted (Collapsed):"),
                  " adjusted FC collapsed to gene level \u2014 regulated if any site qualifies")
        ),
        tags$p(class = "small mb-0 text-muted",
               "Thresholds: |Log\u2082FC| > 1 \u2227 adj. p-value < 0.05")
      )
    }
  })

  # ---------------------------------------------------------------------------
  # Venn diagram renderer  (ggvenn · static ggplot · press-to-render highlight)
  # ---------------------------------------------------------------------------
  # Fill palettes: up = warm colours, down = cool colours (4 experiments max)
  .VF_UP   <- c("#e74c3c", "#e67e22", "#8e44ad", "#d35400")
  .VF_DOWN <- c("#2980b9", "#1abc9c", "#2c3e50", "#7f8c8d")

  render_venn <- function(gene_sets, title, fill_colors, highlight_gene = NULL) {
    n <- length(gene_sets)
    txt_plot <- function(msg, col = "#888888") {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                          label = msg, size = 4.5, colour = col, lineheight = 1.4) +
        ggplot2::theme_void() +
        ggplot2::labs(title = title) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          size = 12, face = "bold", hjust = 0.5))
    }
    if (n < 2) return(txt_plot("Need \u2265 2 experiments"))
    if (n > 4) {
      return(txt_plot(paste0(
        n, " experiments: Venn diagrams support \u2264 4 sets.\n",
        "See the membership table below."
      )))
    }
    if (!.pkg_installed("ggvenn")) {
      return(txt_plot(
        "Venn diagrams require the 'ggvenn' package.\nRun:  install.packages(\"ggvenn\")",
        col = "#cc6600"
      ))
    }
    require_pkg("ggvenn")
    hl_sub <- NULL
    hl <- if (!is.null(highlight_gene)) trimws(highlight_gene) else ""
    if (nchar(hl) > 0) {
      membership <- which(vapply(gene_sets,
                                 function(s) any(grepl(hl, s, ignore.case = TRUE)),
                                 logical(1)))
      hl_sub <- if (length(membership) > 0) {
        paste0("\u25BA ", hl, "  \u2208  ",
               paste(names(gene_sets)[membership], collapse = " \u2229 "))
      } else {
        paste0("\u25BA ", hl, "  \u2014  not found in this direction")
      }
    }
    p <- ggvenn::ggvenn(
      gene_sets,
      fill_color      = fill_colors[seq_len(n)],
      fill_alpha      = 0.55,
      stroke_size     = 0.5,
      set_name_size   = 3.8,
      text_size       = 4.0,
      show_percentage = FALSE
    ) +
      ggplot2::labs(title = title, subtitle = hl_sub) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 12, face = "bold",
                                              hjust = 0.5,
                                              margin = ggplot2::margin(b = 3)),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "#b36a00",
                                              hjust = 0.5, face = "bold",
                                              margin = ggplot2::margin(b = 3))
      )
    p
  }

  # ---------------------------------------------------------------------------
  # Per-layer output factory
  # (uses local() to capture lk in each iteration)
  # ---------------------------------------------------------------------------
  for (.lk in names(.LAYERS)) {
    local({
      lk    <- .lk
      layer <- .LAYERS[[lk]]

      get_sets <- function() {
        req(cptm_rv$n_files >= 2)
        fps <- cptm_rv$files_parsed
        avail_in_all <- all(vapply(fps, function(fp) lk %in% fp$avail_layers, logical(1)))
        if (!avail_in_all) return(NULL)
        lapply(fps, function(fp) fp$layer_sets[[lk]])
      }

      # -- Venn: up -----------------------------------------------------------
      output[[paste0("cptm_venn_", lk, "_up")]] <- renderPlot({
        sets <- get_sets()
        shiny::validate(need(!is.null(sets),
                      paste0(layer$label, " data not available in all uploaded files.")))
        gene_sets <- lapply(sets, `[[`, "up")
        names(gene_sets) <- get_exp_names()
        render_venn(gene_sets,
                    title          = paste0("Up-regulated ", layer$id_type,
                                            "s  (Log\u2082FC > 1, adj.p < 0.05)"),
                    fill_colors    = .VF_UP,
                    highlight_gene = cptm_hl_gene())
      }, res = 96)

      # -- Venn: down ---------------------------------------------------------
      output[[paste0("cptm_venn_", lk, "_down")]] <- renderPlot({
        sets <- get_sets()
        shiny::validate(need(!is.null(sets),
                      paste0(layer$label, " data not available in all uploaded files.")))
        gene_sets <- lapply(sets, `[[`, "down")
        names(gene_sets) <- get_exp_names()
        render_venn(gene_sets,
                    title          = paste0("Down-regulated ", layer$id_type,
                                            "s  (Log\u2082FC < \u22121, adj.p < 0.05)"),
                    fill_colors    = .VF_DOWN,
                    highlight_gene = cptm_hl_gene())
      }, res = 96)

      # -- Summary table ------------------------------------------------------
      output[[paste0("cptm_summary_", lk)]] <- DT::renderDT({
        sets <- get_sets()
        shiny::validate(need(!is.null(sets),
                      paste0(layer$label, " data not available in all uploaded files.")))
        exp_names <- get_exp_names()
        fps       <- cptm_rv$files_parsed
        rows <- lapply(seq_along(fps), function(i) {
          s <- sets[[i]]
          other_up <- unique(unlist(lapply(sets[-i], `[[`, "up")))
          other_dn <- unique(unlist(lapply(sets[-i], `[[`, "down")))
          data.frame(
            Experiment      = exp_names[i],
            Source_File     = fps[[i]]$fname,
            Total_Detected  = s$n_total,
            Upregulated     = length(s$up),
            Up_Unique       = sum(!s$up   %in% other_up),
            Downregulated   = length(s$down),
            Down_Unique     = sum(!s$down %in% other_dn),
            Non_Significant = length(s$nonsig),
            stringsAsFactors = FALSE
          )
        })
        df <- do.call(rbind, rows)
        DT::datatable(
          df,
          rownames = FALSE,
          colnames = c("Experiment", "Source File",
                       paste0("Total ", tools::toTitleCase(layer$id_type), "s"),
                       "Upregulated", "Unique Up\u2020",
                       "Downregulated", "Unique Down\u2020",
                       "Non-Significant"),
          caption = htmltools::tags$caption(
            style = "caption-side:bottom; text-align:left; font-size:0.85em; color:#666;",
            "\u2020 Unique = present in this experiment but absent (same direction) from every other experiment."
          ),
          options = list(pageLength = 15, dom = "t", scrollX = TRUE),
          class   = "compact stripe hover"
        ) %>%
          DT::formatStyle("Upregulated",
            color = "white", backgroundColor = "#dc3545", fontWeight = "bold") %>%
          DT::formatStyle("Downregulated",
            color = "white", backgroundColor = "#0d6efd", fontWeight = "bold")
      })

      # -- Overlap data reactive (shared: DT render + download) ---------------
      overlap_data <- reactive({
        sets <- get_sets()
        req(!is.null(sets))
        exp_names   <- get_exp_names()
        exp_letters <- LETTERS[seq_along(exp_names)]
        all_up  <- unique(unlist(lapply(sets, `[[`, "up")))
        all_dn  <- unique(unlist(lapply(sets, `[[`, "down")))
        all_reg <- union(all_up, all_dn)
        if (length(all_reg) == 0) return(data.frame())
        mat <- data.frame(
          ID        = all_reg,
          Consensus = ifelse(
            all_reg %in% all_up & all_reg %in% all_dn, "Discordant",
            ifelse(all_reg %in% all_up, "Up", "Down")
          ),
          stringsAsFactors = FALSE
        )
        for (i in seq_along(sets)) {
          s <- sets[[i]]
          mat[[exp_names[i]]] <- ifelse(
            all_reg %in% s$up,   "Up",
            ifelse(all_reg %in% s$down, "Down", "")
          )
        }
        mat$Up_Group <- vapply(all_reg, function(g) {
          present <- which(vapply(sets, function(s) g %in% s$up, logical(1)))
          if (length(present) == 0) "\u2014"
          else paste(exp_letters[present], collapse = "\u2229")
        }, character(1))
        mat$Down_Group <- vapply(all_reg, function(g) {
          present <- which(vapply(sets, function(s) g %in% s$down, logical(1)))
          if (length(present) == 0) "\u2014"
          else paste(exp_letters[present], collapse = "\u2229")
        }, character(1))
        mat$N_Exp_Regulated <- rowSums(
          mat[, exp_names, drop = FALSE] != "", na.rm = TRUE
        )
        mat <- mat[order(-mat$N_Exp_Regulated, mat$Consensus, mat$ID), ]
        mat
      })

      # -- Overlap table ------------------------------------------------------
      output[[paste0("cptm_overlap_", lk)]] <- DT::renderDT({
        mat       <- overlap_data()
        shiny::validate(need(nrow(mat) > 0,
                      "No regulated entries found across the uploaded experiments."))
        exp_names   <- get_exp_names()
        exp_letters <- LETTERS[seq_along(exp_names)]
        key_str     <- paste(paste0(exp_letters[seq_along(exp_names)], " = ", exp_names),
                             collapse = ";  ")
        id_label <- if (layer$id_type == "PTM site") "PTM Site" else "Gene"
        DT::datatable(
          mat, rownames = FALSE,
          colnames = c(id_label, "Consensus", exp_names,
                       "Up Venn Group", "Down Venn Group", "# Experiments"),
          filter   = "top",
          caption  = htmltools::tags$caption(
            style = "caption-side:bottom;text-align:left;font-size:0.85em;color:#666;",
            paste0("Set key:  ", key_str)
          ),
          options = list(pageLength = 25, scrollX = TRUE, dom = "ftp"),
          class   = "compact stripe hover"
        ) %>%
          DT::formatStyle(
            exp_names,
            backgroundColor = DT::styleEqual(c("Up","Down",""), c("#ffe0e0","#dde8ff","#ffffff")),
            color           = DT::styleEqual(c("Up","Down",""), c("#9b0000","#003580","#888888")),
            fontWeight      = DT::styleEqual(c("Up","Down",""), c("bold","bold","normal"))
          ) %>%
          DT::formatStyle(
            "Consensus",
            backgroundColor = DT::styleEqual(
              c("Up","Down","Discordant"),
              c("#ffe0e0","#dde8ff","#fff3cd")
            )
          )
      })

      # -- Overlap download ---------------------------------------------------
      output[[paste0("cptm_dl_overlap_", lk)]] <- downloadHandler(
        filename = function() {
          paste0("PTMCrossExp_", gsub("[^A-Za-z0-9]", "_", layer$label),
                 "_Overlap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
          mat <- tryCatch(overlap_data(), error = function(e) data.frame())
          utils::write.csv(mat, file, row.names = FALSE)
        }
      )

    }) # end local()
  } # end for loop

} # end server_ptm_cross_exp
