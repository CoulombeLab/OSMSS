# server_cross_exp.R ===========================================================
# Cross-Experiment Comparison module.
# Accepts 2+ Protein > DE result CSV files, auto-detects gene name / Log2FC /
# significance columns, and produces:
#   - Two Venn diagrams (up-regulated and down-regulated gene sets)
#   - Per-experiment summary statistics table
#   - Gene-level overlap table
#
# Thresholds (hard-coded, matching volcano plot convention):
#   Up-regulated  : Log2FC > 1  AND  adj. p-value (or BFDR) < 0.05
#   Down-regulated: Log2FC < -1 AND  adj. p-value (or BFDR) < 0.05
#
# Called from server.R as: server_cross_exp(input, output, session, rv, logger)
# ==============================================================================

server_cross_exp <- function(input, output, session, rv, logger) {

  CE_FC_UP    <-  1
  CE_FC_DOWN  <- -1
  CE_P_THRESH <- 0.05

  # ---------------------------------------------------------------------------
  # Reactive storage
  # ---------------------------------------------------------------------------
  ce_rv <- reactiveValues(
    files_parsed  = list(),       # list of structs (one per valid file)
    parse_errors  = character(0),
    n_files       = 0L
  )

  # ---------------------------------------------------------------------------
  # Column detection helpers
  # ---------------------------------------------------------------------------

  #' Detect the gene-name, Log2FC, and significance columns in an OSMSS DE CSV.
  #' Priority lists are ordered from most-specific to most-generic to avoid
  #' false positives across MSstats / SAINT / DIA-limma / ProDA output formats.
  detect_de_cols <- function(df) {
    nms <- names(df)

    # Gene name column
    gene_col <- NA_character_
    for (cn in c("GeneName", "Gene", "Prey", "Protein")) {
      if (cn %in% nms) { gene_col <- cn; break }
    }

    # Log2 fold-change column
    fc_col <- NA_character_
    for (cn in c("Log2FC", "log2FC", "log2fc")) {
      if (cn %in% nms) { fc_col <- cn; break }
    }

    # Significance column
    # Covers: MSstats (Adj.Pvalue), DIA limma (AdjP_BH), ProDA/plot_data
    # (AdjPValue), generic (AdjPvalue), SAINT (BFDR), DIA LOD (PermFDR)
    sig_col   <- NA_character_
    sig_label <- NA_character_
    for (cn in c("Adj.Pvalue", "AdjP_BH", "AdjPValue", "AdjPvalue",
                 "BFDR", "PermFDR")) {
      if (cn %in% nms) {
        sig_col   <- cn
        sig_label <- switch(cn,
          "Adj.Pvalue" = "Adj. P-value (BH)",
          "AdjP_BH"    = "Adj. P-value (BH)",
          "AdjPValue"  = "Adj. P-value",
          "AdjPvalue"  = "Adj. P-value",
          "BFDR"       = "BFDR (SAINTexpress)",
          "PermFDR"    = "Permutation FDR",
          cn
        )
        break
      }
    }

    list(gene_col  = gene_col,
         fc_col    = fc_col,
         sig_col   = sig_col,
         sig_label = sig_label)
  }

  #' Classify genes in a parsed data.frame into up/down/non-significant sets.
  classify_genes <- function(df, gene_col, fc_col, sig_col) {
    fc   <- suppressWarnings(as.numeric(df[[fc_col]]))
    sig  <- suppressWarnings(as.numeric(df[[sig_col]]))
    gene <- trimws(as.character(df[[gene_col]]))

    sig[is.na(sig)]                   <- 1          # treat missing p as ns
    gene[gene == "" | gene == "NA"]   <- NA_character_

    valid     <- !is.na(gene) & !is.na(fc) & is.finite(fc)
    up_mask   <- valid & fc >  CE_FC_UP   & sig < CE_P_THRESH
    down_mask <- valid & fc <  CE_FC_DOWN & sig < CE_P_THRESH

    list(
      up      = unique(gene[up_mask]),
      down    = unique(gene[down_mask]),
      nonsig  = unique(gene[valid & !up_mask & !down_mask]),
      n_total = sum(valid)
    )
  }

  # ---------------------------------------------------------------------------
  # File parsing
  # ---------------------------------------------------------------------------
  observeEvent(input$ce_files, {
    req(input$ce_files)
    finfo <- input$ce_files
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

      cols <- detect_de_cols(df)

      if (is.na(cols$gene_col)) {
        errors <- c(errors, paste0(fname, ": Could not find gene name column ",
                                   "(expected GeneName, Gene, Prey, or Protein)."))
        next
      }
      if (is.na(cols$fc_col)) {
        errors <- c(errors, paste0(fname, ": Could not find Log2FC column."))
        next
      }
      if (is.na(cols$sig_col)) {
        errors <- c(errors,
          paste0(fname, ": Could not find a significance column. ",
                 "Expected one of: Adj.Pvalue, AdjP_BH, AdjPValue, AdjPvalue, BFDR, PermFDR."))
        next
      }

      classified <- classify_genes(df, cols$gene_col, cols$fc_col, cols$sig_col)

      # Auto-name: strip OSMSS default prefix and extension
      exp_name <- tools::file_path_sans_ext(basename(fname))
      exp_name <- sub("^DE_Table_\\d{8}_\\d{6}$", "DE Table", exp_name)
      exp_name <- sub("^DE_Table_", "", exp_name)

      parsed[[i]] <- list(
        fname     = fname,
        exp_name  = exp_name,
        gene_col  = cols$gene_col,
        fc_col    = cols$fc_col,
        sig_col   = cols$sig_col,
        sig_label = cols$sig_label,
        up        = classified$up,
        down      = classified$down,
        nonsig    = classified$nonsig,
        n_total   = classified$n_total
      )
    }

    valid_parsed <- Filter(Negate(is.null), parsed)

    # Deduplicate experiment names — duplicate names (e.g. multiple OSMSS
    # default-named files all auto-labelled "DE Table") collapse to a single
    # column inside ggvenn's list_to_data_frame(), making n_sets < 2.
    if (length(valid_parsed) > 0) {
      raw_names <- vapply(valid_parsed, `[[`, character(1), "exp_name")
      dedup     <- make.unique(raw_names, sep = " ")
      for (k in seq_along(valid_parsed))
        valid_parsed[[k]]$exp_name <- dedup[k]
    }

    ce_rv$files_parsed <- valid_parsed
    ce_rv$n_files      <- length(valid_parsed)
    ce_rv$parse_errors <- errors

    for (e in errors) showNotification(e, type = "warning", duration = 12)

    if (length(valid_parsed) >= 2) {
      showNotification(
        paste0(length(valid_parsed), " DE file(s) loaded successfully."),
        type = "message", duration = 5
      )
    } else if (length(valid_parsed) == 1) {
      showNotification("Only 1 valid file found — upload at least 2 to compare.",
                       type = "warning", duration = 8)
    }
  })

  # ---------------------------------------------------------------------------
  # Experiment name inputs (dynamic text fields, one per valid file)
  # ---------------------------------------------------------------------------
  output$ce_exp_name_inputs <- renderUI({
    req(ce_rv$n_files >= 2)
    tagList(
      lapply(seq_len(ce_rv$n_files), function(i) {
        textInput(
          inputId   = paste0("ce_exp_name_", i),
          label     = NULL,
          value     = ce_rv$files_parsed[[i]]$exp_name,
          width     = "100%",
          placeholder = paste0("Experiment ", i)
        )
      })
    )
  })

  #' Collect current experiment labels from dynamic text inputs
  get_exp_names_now <- reactive({
    req(ce_rv$n_files >= 2)
    vapply(seq_len(ce_rv$n_files), function(i) {
      nm <- input[[paste0("ce_exp_name_", i)]]
      if (is.null(nm) || nchar(trimws(nm)) == 0) paste0("Exp ", i)
      else trimws(nm)
    }, character(1))
  })
  # Debounced: plots and tables update once typing pauses, not on every keystroke
  # (per-keystroke re-renders raced the DataTables requests and raised
  # "column name ... is not found in data" alerts).
  get_exp_names <- debounce(get_exp_names_now, 700)

  # ---------------------------------------------------------------------------
  # Ready flag (used by conditionalPanel in ui.R)
  # ---------------------------------------------------------------------------
  output$ce_ready <- reactive({ ce_rv$n_files >= 2 })
  outputOptions(output, "ce_ready", suspendWhenHidden = FALSE)

  # ---------------------------------------------------------------------------
  # Highlight search state (press-to-render)
  # ---------------------------------------------------------------------------
  ce_hl_gene <- reactiveVal(NULL)

  observeEvent(input$ce_highlight_btn, {
    hl <- if (is.null(input$ce_search_input)) "" else trimws(input$ce_search_input)
    ce_hl_gene(hl)
  })

  observeEvent(input$ce_clear_btn, {
    ce_hl_gene(NULL)
    updateTextInput(session, "ce_search_input", value = "")
  })

  output$ce_search_result_ui <- renderUI({
    hl <- ce_hl_gene()
    req(!is.null(hl) && nchar(hl) > 0 && ce_rv$n_files >= 2)
    fps       <- ce_rv$files_parsed
    exp_names <- get_exp_names()
    badges <- lapply(seq_along(fps), function(i) {
      fp <- fps[[i]]
      in_up   <- any(grepl(hl, fp$up,   ignore.case = TRUE))
      in_down <- any(grepl(hl, fp$down, ignore.case = TRUE))
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
  # Status message shown in the sidebar
  # ---------------------------------------------------------------------------
  output$ce_status_msg <- renderUI({
    if (ce_rv$n_files == 0) {
      return(tags$p(class = "text-muted small mt-2",
                    icon("circle-info"),
                    " Upload \u2265 2 OSMSS DE result CSV files to begin."))
    }
    if (ce_rv$n_files < 2) {
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

    tagList(
      tags$p(class = "text-success small mt-2",
             icon("circle-check"),
             paste0(" ", ce_rv$n_files, " files ready.")),
      venn_warn
    )
  })

  # ---------------------------------------------------------------------------
  # Venn diagram renderer  (ggvenn · static ggplot · press-to-render highlight)
  # ---------------------------------------------------------------------------
  # Fill palettes: up = warm colours, down = cool colours (4 experiments max)
  .VF_UP   <- c("#e74c3c", "#e67e22", "#8e44ad", "#d35400")
  .VF_DOWN <- c("#2980b9", "#1abc9c", "#2c3e50", "#7f8c8d")

  #' Build a static ggvenn Venn diagram.
  #' Supports 2 sets (circles), 3 sets (circles), 4 sets (ellipses).
  #' 5+ sets: returns a text-only placeholder; use the membership table.
  #' highlight_gene: if non-empty, a subtitle names the sets containing it.
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

    # Build highlight subtitle
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
  # Output: Venn — upregulated
  # ---------------------------------------------------------------------------
  output$ce_venn_up <- renderPlot({
    req(ce_rv$n_files >= 2)
    gene_sets <- lapply(ce_rv$files_parsed, `[[`, "up")
    names(gene_sets) <- get_exp_names()
    render_venn(gene_sets,
                title          = "Upregulated  (Log\u2082FC > 1, FDR < 0.05)",
                fill_colors    = .VF_UP,
                highlight_gene = ce_hl_gene())
  }, res = 96)

  # ---------------------------------------------------------------------------
  # Output: Venn — downregulated
  # ---------------------------------------------------------------------------
  output$ce_venn_down <- renderPlot({
    req(ce_rv$n_files >= 2)
    gene_sets <- lapply(ce_rv$files_parsed, `[[`, "down")
    names(gene_sets) <- get_exp_names()
    render_venn(gene_sets,
                title          = "Downregulated  (Log\u2082FC < \u22121, FDR < 0.05)",
                fill_colors    = .VF_DOWN,
                highlight_gene = ce_hl_gene())
  }, res = 96)

  # ---------------------------------------------------------------------------
  # Output: per-experiment summary table
  # ---------------------------------------------------------------------------
  output$ce_summary_table <- DT::renderDT({
    req(ce_rv$n_files >= 2)
    exp_names <- get_exp_names()
    fps       <- ce_rv$files_parsed

    rows <- lapply(seq_along(fps), function(i) {
      fp <- fps[[i]]

      # Genes exclusive to this experiment (not found in any other experiment)
      other_up   <- unique(unlist(lapply(fps[-i], `[[`, "up")))
      other_down <- unique(unlist(lapply(fps[-i], `[[`, "down")))

      n_up_excl   <- sum(!fp$up   %in% other_up)
      n_down_excl <- sum(!fp$down %in% other_down)

      data.frame(
        Experiment       = exp_names[i],
        Source_File      = fp$fname,
        Total_Detected   = fp$n_total,
        Upregulated      = length(fp$up),
        Up_Unique        = n_up_excl,
        Downregulated    = length(fp$down),
        Down_Unique      = n_down_excl,
        Non_Significant  = length(fp$nonsig),
        Sig_Column       = fp$sig_label,
        stringsAsFactors = FALSE
      )
    })

    df <- do.call(rbind, rows)

    DT::datatable(
      df,
      rownames = FALSE,
      colnames = c(
        "Experiment", "Source File", "Total Detected",
        "Upregulated", "Unique Up\u2020",
        "Downregulated", "Unique Down\u2020",
        "Non-Significant", "Sig. Column"
      ),
      caption  = htmltools::tags$caption(
        style = "caption-side: bottom; text-align: left; font-size: 0.85em; color: #666;",
        "\u2020 Unique = gene found in this experiment but not in any other uploaded experiment (for the same direction)."
      ),
      options  = list(pageLength = 15, dom = "t", scrollX = TRUE),
      class    = "compact stripe hover"
    ) %>%
      DT::formatStyle(
        "Upregulated",
        color           = "white",
        backgroundColor = "#dc3545",
        fontWeight      = "bold"
      ) %>%
      DT::formatStyle(
        "Downregulated",
        color           = "white",
        backgroundColor = "#0d6efd",
        fontWeight      = "bold"
      )
  }, server = FALSE)  # small table: render in the browser, no per-request server lookups

  # ---------------------------------------------------------------------------
  # Gene-level overlap data reactive (shared by DT render and download handler)
  # ---------------------------------------------------------------------------
  ce_overlap_data <- reactive({
    req(ce_rv$n_files >= 2)
    exp_names   <- get_exp_names()
    fps         <- ce_rv$files_parsed
    exp_letters <- LETTERS[seq_along(exp_names)]

    all_up  <- unique(unlist(lapply(fps, `[[`, "up")))
    all_dn  <- unique(unlist(lapply(fps, `[[`, "down")))
    all_reg <- union(all_up, all_dn)
    if (length(all_reg) == 0) return(data.frame())

    mat <- data.frame(
      Gene      = all_reg,
      Consensus = ifelse(
        all_reg %in% all_up & all_reg %in% all_dn, "Discordant",
        ifelse(all_reg %in% all_up, "Up", "Down")
      ),
      stringsAsFactors = FALSE
    )

    # Per-experiment direction columns
    for (i in seq_along(fps)) {
      mat[[exp_names[i]]] <- ifelse(
        all_reg %in% fps[[i]]$up,   "Up",
        ifelse(all_reg %in% fps[[i]]$down, "Down", "")
      )
    }

    # Venn group: which experiments contain this gene per direction
    mat$Up_Group <- vapply(all_reg, function(g) {
      present <- which(vapply(fps, function(fp) g %in% fp$up, logical(1)))
      if (length(present) == 0) "\u2014"
      else paste(exp_letters[present], collapse = "\u2229")
    }, character(1))

    mat$Down_Group <- vapply(all_reg, function(g) {
      present <- which(vapply(fps, function(fp) g %in% fp$down, logical(1)))
      if (length(present) == 0) "\u2014"
      else paste(exp_letters[present], collapse = "\u2229")
    }, character(1))

    mat$N_Exp_Regulated <- rowSums(
      mat[, exp_names, drop = FALSE] != "", na.rm = TRUE
    )
    mat <- mat[order(-mat$N_Exp_Regulated, mat$Consensus, mat$Gene), ]
    mat
  })

  # ---------------------------------------------------------------------------
  # Output: gene-level overlap table
  # ---------------------------------------------------------------------------
  output$ce_overlap_table <- DT::renderDT({
    mat       <- ce_overlap_data()
    shiny::validate(need(nrow(mat) > 0,
                  "No regulated genes found across the uploaded experiments."))
    exp_names   <- get_exp_names()
    exp_letters <- LETTERS[seq_along(exp_names)]
    key_str     <- paste(paste0(exp_letters[seq_along(exp_names)], " = ", exp_names),
                         collapse = ";  ")

    DT::datatable(
      mat,
      rownames = FALSE,
      colnames = c("Gene", "Consensus", exp_names,
                   "Up Venn Group", "Down Venn Group", "# Experiments"),
      filter   = "top",
      caption  = htmltools::tags$caption(
        style = "caption-side:bottom;text-align:left;font-size:0.85em;color:#666;",
        paste0("Set key:  ", key_str)
      ),
      options  = list(
        pageLength = 25, scrollX = TRUE, dom = "ftp",
        columnDefs = list(
          list(targets = 0, className = "dt-left"),
          list(targets = 1, className = "dt-center")
        )
      ),
      class = "compact stripe hover"
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
  }, server = FALSE)  # render in the browser: renaming an experiment can't race a server-side data request

  # Download handler for the overlap table
  output$ce_dl_overlap <- downloadHandler(
    filename = function() {
      paste0("CrossExp_Overlap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      utils::write.csv(ce_overlap_data(), file, row.names = FALSE)
    }
  )

  # ---------------------------------------------------------------------------
  # Output: content area header / description (rendered when no files yet)
  # ---------------------------------------------------------------------------
  output$ce_panel_header <- renderUI({
    if (ce_rv$n_files < 2) {
      tags$div(
        class = "alert alert-info py-3 px-4",
        tags$h5(icon("circle-nodes"), " Cross-Experiment Comparison"),
        tags$p(class = "mb-1",
          "Upload two or more ", tags$b("DE Result CSV"), " files from the ",
          tags$b("Protein \u2192 Differential Expression"), " module using the sidebar."
        ),
        tags$ul(class = "mb-0 small",
          tags$li(tags$b("Upregulated: "), "Log\u2082FC > 1 and FDR / BFDR < 0.05"),
          tags$li(tags$b("Downregulated: "), "Log\u2082FC < \u22121 and FDR / BFDR < 0.05"),
          tags$li("Compatible with MSstats, SAINTexpress, ProDA, and DIA-limma outputs")
        )
      )
    }
  })

}
