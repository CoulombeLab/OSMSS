# server_local_ms.R =============================================================
# Server logic for the "Analyze Local MS" module:
#   - FragPipe PTM Analysis (Modification Viewer + Data Table)
# Database management (import / metadata) is handled by server_fragpipe_db.R.
# Called from server.R as: server_local_ms(input, output, session, rv, logger)
# All helper functions (build_modification_plot, etc.)
# are defined in global.R to avoid duplication.
# ==============================================================================

server_local_ms <- function(input, output, session, rv, logger) {

  # ---- FragpipeDB reader from session cache ----
  get_fragpipe_db <- function() {
    rv$fpdb_fragpipe_db
  }

  # ---- Helper: get allowed ExperimentIDs for the selected organism ----
  get_organism_experiment_ids <- function(organism) {
    exp_db <- rv$fpdb_experiment_db
    if (is.null(exp_db) || nrow(exp_db) == 0) return(NULL)
    if (!"Organism" %in% colnames(exp_db)) return(NULL)
    # Filter experiments matching the organism
    matched <- exp_db[exp_db$Organism == organism, , drop = FALSE]
    if (nrow(matched) == 0) return(character(0))
    # Extract all filenames (ExperimentIDs) from matching experiments
    all_fnames <- unlist(lapply(matched[["Order of Filenames"]], function(x) {
      f <- trimws(unlist(strsplit(as.character(x), "\n")))
      f[nchar(f) > 0]
    }))
    unique(all_fnames)
  }

  filter_db_by_organism <- function(db, organism) {
    allowed <- get_organism_experiment_ids(organism)
    if (is.null(allowed)) return(db)   # no exp DB or no Organism col — pass through
    if (length(allowed) == 0) return(db[0, , drop = FALSE])  # no matching experiments
    db %>% filter(ExperimentID %in% allowed)
  }

  # ---- Plot logic ----
  ms_plot_result <- reactiveVal(NULL)
  ms_plot_key <- reactiveVal(0L)  # unique key to force full widget rebuild
  ms_domains_all <- reactiveVal(NULL)  # full domain list for checkbox
  
  # Domain select all / deselect all
  observeEvent(input$ms_dom_all, {
    all_dom <- ms_domains_all()
    if (!is.null(all_dom) && nrow(all_dom) > 0) {
      updateCheckboxGroupInput(session, "ms_selected_domains",
                               selected = unique(all_dom$label))
    }
  })
  observeEvent(input$ms_dom_none, {
    updateCheckboxGroupInput(session, "ms_selected_domains", selected = character(0))
  })
  
  # When UniProt ID changes or Generate Plot is clicked, fetch domains
  observeEvent(input$ms_plot_btn, {
    ms_plot_result(NULL)
    ms_plot_key(isolate(ms_plot_key()) + 1L)  # bump key to force fresh widget
    
    if (is.null(rv$fpdb_fragpipe_db) || nrow(rv$fpdb_fragpipe_db) == 0) {
      ms_plot_result(list(plot = NULL, error = "No database loaded. Go to Load Database to import or create one."))
      return()
    }
    
    gene <- trimws(input$ms_gene_name)
    uniprot <- trimws(input$ms_uniprot_id)
    if (nchar(gene) == 0 || nchar(uniprot) == 0) {
      ms_plot_result(list(plot = NULL, error = "Please enter both Gene Name and UniProt ID."))
      return()
    }
    
    withProgress(message = "Building plot...", {
      incProgress(0.1, detail = "Fetching domains...")
      all_dom <- fetch_uniprot_domains(uniprot)
      ms_domains_all(all_dom)
      
      # Update domain checkbox choices — default to none selected
      if (nrow(all_dom) > 0) {
        choices <- unique(all_dom$label)
        updateCheckboxGroupInput(session, "ms_selected_domains",
                                 choices = choices, selected = character(0))
      } else {
        updateCheckboxGroupInput(session, "ms_selected_domains",
                                 choices = character(0), selected = character(0))
      }
      
      incProgress(0.3, detail = "Loading data...")
      db <- get_fragpipe_db()
      if (is.null(db)) { ms_plot_result(list(plot = NULL, error = "Failed to read FragpipeDB.")); return() }
      db <- filter_db_by_organism(db, input$ms_organism)
      
      # No domains selected by default
      sel_domains <- NULL
      
      incProgress(0.3, detail = "Rendering plot...")
      res <- build_modification_plot(db, gene, uniprot,
                                      y_mode = input$ms_view_mode,
                                      domains_df = sel_domains)
      ms_plot_result(res)
    })
  })
  
  # Refresh plot when domain selection changes
  observeEvent(input$ms_selected_domains, {
    req(ms_domains_all())
    req(ms_plot_result())  # only if we already have a plot
    
    gene <- trimws(input$ms_gene_name)
    uniprot <- trimws(input$ms_uniprot_id)
    if (nchar(gene) == 0 || nchar(uniprot) == 0) return()
    if (is.null(rv$fpdb_fragpipe_db) || nrow(rv$fpdb_fragpipe_db) == 0) return()
    
    db <- get_fragpipe_db()
    if (is.null(db)) return()
    db <- filter_db_by_organism(db, input$ms_organism)
    all_dom <- ms_domains_all()
    
    # Filter to selected labels
    if (nrow(all_dom) > 0 && length(input$ms_selected_domains) > 0) {
      sel_domains <- all_dom %>% filter(label %in% input$ms_selected_domains)
    } else {
      sel_domains <- NULL
    }
    
    res <- build_modification_plot(db, gene, uniprot,
                                    y_mode = input$ms_view_mode,
                                    domains_df = sel_domains)
    ms_plot_result(res)
  }, ignoreInit = TRUE)
  
  # Refresh plot when view mode toggle changes
  observeEvent(input$ms_view_mode, {
    req(ms_plot_result())  # only if we already have a plot
    
    gene <- trimws(input$ms_gene_name)
    uniprot <- trimws(input$ms_uniprot_id)
    if (nchar(gene) == 0 || nchar(uniprot) == 0) return()
    if (is.null(rv$fpdb_fragpipe_db) || nrow(rv$fpdb_fragpipe_db) == 0) return()
    
    db <- get_fragpipe_db()
    if (is.null(db)) return()
    db <- filter_db_by_organism(db, input$ms_organism)
    all_dom <- ms_domains_all()
    sel_domains <- NULL
    if (!is.null(all_dom) && nrow(all_dom) > 0 && length(input$ms_selected_domains) > 0) {
      sel_domains <- all_dom %>% filter(label %in% input$ms_selected_domains)
    }
    
    res <- build_modification_plot(db, gene, uniprot,
                                    y_mode = input$ms_view_mode,
                                    domains_df = sel_domains)
    ms_plot_result(res)
  }, ignoreInit = TRUE)
  
  output$ms_plot_container <- renderUI({
    res <- ms_plot_result()
    key <- ms_plot_key()
    # Dynamic plotlyOutput with unique ID forces full teardown/rebuild
    plot_id <- paste0("ms_mod_plot_", key)
    tagList(
      plotlyOutput(plot_id, height = "800px")
    )
  })

  observe({
    res <- ms_plot_result()
    key <- ms_plot_key()
    plot_id <- paste0("ms_mod_plot_", key)
    output[[plot_id]] <- renderPlotly({
      if (is.null(res) || is.null(res$plot)) return(plotly_empty())
      res$plot
    })
  })
  
  output$ms_plot_status <- renderPrint({
    res <- ms_plot_result()
    if (is.null(res)) {
      cat("Enter a gene name and UniProt ID, then click Generate Plot.\n")
    } else if (!is.null(res$error)) {
      cat(res$error, "\n")
    } else {
      cat("Plot generated successfully.\n")
    }
  })
  
  # ---- Data Table logic ----
  ms_table_data <- reactiveVal(NULL)
  
  observeEvent(input$ms_table_btn, {
    ms_table_data(NULL)
    if (!isTRUE(rv$fpdb_db_loaded) || is.null(rv$fpdb_fragpipe_db)) {
      ms_table_data("NO_DB")
      return()
    }
    gene <- trimws(input$ms_table_gene)
    if (nchar(gene) == 0) {
      ms_table_data("NO_GENE")
      return()
    }
    db <- get_fragpipe_db()
    if (is.null(db)) { ms_table_data("NO_DB"); return() }
    db <- filter_db_by_organism(db, input$ms_table_organism)
    gene_rows <- db %>% filter(Gene == gene)
    if (nrow(gene_rows) == 0) {
      ms_table_data("NO_ROWS")
    } else {
      ms_table_data(gene_rows)
    }
  })
  
  output$ms_gene_table <- renderDT({
    td <- ms_table_data()
    if (is.null(td) || is.character(td)) return(NULL)
    datatable(td, filter = "top", rownames = FALSE,
              options = list(pageLength = 25, scrollX = TRUE,
                             dom = "lftipr",
                             order = list(list(0, "asc"))))
  })
  
  output$ms_table_status <- renderPrint({
    td <- ms_table_data()
    if (is.null(td)) {
      cat("Enter a gene name and click Load Table.\n")
    } else if (identical(td, "NO_DB")) {
      cat("No FragpipeDB.csv found. Import data first.\n")
    } else if (identical(td, "NO_GENE")) {
      cat("Please enter a gene name.\n")
    } else if (identical(td, "NO_ROWS")) {
      cat("No entries found for this gene.\n")
    } else {
      cat(paste0("Showing ", nrow(td), " entries for gene '", input$ms_table_gene, "'.\n"))
    }
  })
}
