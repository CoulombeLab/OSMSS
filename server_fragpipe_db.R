# server_fragpipe_db.R ==========================================================
# Server logic for FragPipe > Manage Database module.
# Provides:
#   - Compact list view of experiments from MS_Experiment_DB.csv
#   - Table view with CSV download
#   - Add Experiment (with or without data)
#   - Edit / Delete experiment metadata and data
#   - Upload data for experiments missing FragpipeDB entries
#
# Called from server.R as: server_fragpipe_db(input, output, session, rv, logger)
# ==============================================================================

server_fragpipe_db <- function(input, output, session, rv, logger) {

  # ---- Reactive storage ----
  fpdb <- reactiveValues(
    exp_db         = NULL,   # data.frame from MS_Experiment_DB (rv-backed)
    db_version     = 0L,     # bump to trigger refresh
    selected_exp   = NULL,   # Experiment ID currently selected
    detail_mode    = "view"  # "view" or "edit"
  )

  # ---- Load / reload experiment DB from session cache ----
  observe({
    fpdb$db_version  # dependency
    fpdb$exp_db <- rv$fpdb_experiment_db
  })

  # ---- Cached FragpipeDB (column-subset for stats/summary) ----
  frag_db_cache <- reactive({
    fpdb$db_version
    db <- rv$fpdb_fragpipe_db
    if (is.null(db)) return(NULL)
    db[, intersect(c("ExperimentID", "Gene", "Modification", "UserID"), colnames(db)), drop = FALSE]
  })

  # ---- Show/hide loading overlay based on db_loaded state ----
  observe({
    if (isTRUE(rv$fpdb_db_loaded)) {
      shinyjs::runjs("$('#fpdb_init_loading').fadeOut(300);")
    }
  })

  # =========================================================================
  # LOAD DATABASE — Import observer
  # =========================================================================
  observeEvent(input$fpdb_import_btn, {
    req(input$fpdb_import_fragpipe, input$fpdb_import_experiments)

    # Read FragpipeDB
    fdb <- tryCatch(
      readr::read_csv(input$fpdb_import_fragpipe$datapath, show_col_types = FALSE),
      error = function(e) { showNotification(paste0("Error reading PTM dataset: ", e$message), type = "error", duration = 12); NULL }
    )
    if (is.null(fdb)) return()
    required_fdb_cols <- c("ExperimentID", "Gene", "Modification")
    missing_fdb <- setdiff(required_fdb_cols, colnames(fdb))
    if (length(missing_fdb) > 0) {
      showNotification(paste0("FragpipeDB.csv is missing required columns: ",
                              paste(missing_fdb, collapse = ", ")),
                       type = "error", duration = 12)
      return()
    }

    # Read MS_Experiment_DB
    edb <- tryCatch(
      read.csv(input$fpdb_import_experiments$datapath,
               stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) { showNotification(paste0("Error reading experiment dataset: ", e$message), type = "error", duration = 12); NULL }
    )
    if (is.null(edb)) return()
    if (!"Experiment ID" %in% colnames(edb)) {
      showNotification("MS_Experiment_DB.csv is missing the 'Experiment ID' column.",
                       type = "error", duration = 12)
      return()
    }
    colnames(edb) <- trimws(colnames(edb))
    if (!"Organism" %in% colnames(edb)) edb$Organism <- "Human"

    rv$fpdb_fragpipe_db   <- fdb
    rv$fpdb_experiment_db <- edb
    rv$fpdb_db_loaded     <- TRUE
    fpdb$db_version       <- fpdb$db_version + 1L
    logger(paste0("FragPipe DB: Imported ", nrow(edb), " experiments, ",
                  nrow(fdb), " PTM entries."))
    showNotification(
      paste0("Database loaded: ", nrow(edb), " experiments, ", nrow(fdb), " PTM entries."),
      type = "message", duration = 8
    )
  })

  # =========================================================================
  # LOAD DATABASE — New empty database observer
  # =========================================================================
  observeEvent(input$fpdb_new_db_btn, {
    edb <- data.frame(
      Lab                   = character(0),
      Investigator          = character(0),
      `Experiment Year`     = character(0),
      `Experiment ID`       = character(0),
      `Cell/Tissue type`    = character(0),
      `Experiment Type`     = character(0),
      `Experimental Design` = character(0),
      Samples               = character(0),
      Notes                 = character(0),
      `Order of conditions` = character(0),
      `Order of Filenames`  = character(0),
      Enzyme                = character(0),
      `DIA or DDA`          = character(0),
      Organism              = character(0),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    fdb <- data.frame(
      Gene                             = character(0),
      Peptide                          = character(0),
      Modification                     = character(0),
      Residue_Location                 = character(0),
      Protein_Start                    = integer(0),
      Protein_End                      = integer(0),
      Modification_Location_in_Protein = numeric(0),
      Delta_Mass                       = numeric(0),
      Additional_Gene_Hits             = character(0),
      Best_PeptideProphet              = numeric(0),
      Best_Expectation                 = numeric(0),
      Best_Hyperscore                  = numeric(0),
      Best_Score_With_DM               = numeric(0),
      Best_Score_Without_DM            = numeric(0),
      Is_Unique                        = logical(0),
      Min_Missed_Cleavages             = integer(0),
      Has_Localization                 = logical(0),
      Mass_Error                       = numeric(0),
      In_Peptide_File                  = logical(0),
      Total_Peptide_Identified         = integer(0),
      Times_Mod_Identified             = integer(0),
      Mod_Type_Class                   = character(0),
      Confidence_Score                 = numeric(0),
      Biological_Confidence            = character(0),
      UserID                           = character(0),
      ExperimentID                     = character(0),
      stringsAsFactors = FALSE
    )

    rv$fpdb_fragpipe_db   <- fdb
    rv$fpdb_experiment_db <- edb
    rv$fpdb_db_loaded     <- TRUE
    fpdb$db_version       <- fpdb$db_version + 1L
    logger("FragPipe DB: New empty database created.")
    showNotification("Empty database created. Add experiments in Manage Database.",
                     type = "message", duration = 8)
  })

  # ---- Load status output ----
  output$fpdb_load_status <- renderUI({
    if (!isTRUE(rv$fpdb_db_loaded)) {
      tags$div(class = "alert alert-info text-center py-3",
        icon("circle-info"), " No database loaded yet.",
        tags$br(),
        tags$small(class = "text-muted", "Import or create a database above to get started.")
      )
    } else {
      tags$div(class = "alert alert-success py-3",
        icon("circle-check"), tags$strong(" Database loaded."),
        tags$br(),
        tags$small(
          nrow(rv$fpdb_experiment_db), " experiment(s)  \u2022  ",
          nrow(rv$fpdb_fragpipe_db), " PTM entries"
        )
      )
    }
  })

  # ---- No-DB banner for Manage Database ----
  output$fpdb_no_db_banner <- renderUI({
    if (!isTRUE(rv$fpdb_db_loaded)) {
      tags$div(
        class = "alert alert-warning text-center py-4 my-3",
        icon("triangle-exclamation"), tags$strong(" No database loaded."),
        tags$br(), tags$br(),
        "Go to ", tags$b("Load Database"), " to import an existing database or create a new one."
      )
    } else {
      NULL
    }
  })

  # ---- Pre-compute stats for ALL experiments in one pass ----
  all_exp_stats <- reactive({
    exp_db  <- fpdb$exp_db
    frag_db <- frag_db_cache()
    if (is.null(exp_db) || nrow(exp_db) == 0) return(list())

    empty <- list(has_data = FALSE, n_entries = 0, n_genes = 0, n_mods = 0,
                  user_ids = character(0), file_ids = character(0))

    stats <- lapply(seq_len(nrow(exp_db)), function(i) {
      fnames_str <- as.character(exp_db[i, "Order of Filenames"])
      fnames <- trimws(unlist(strsplit(fnames_str, "\n")))
      fnames <- fnames[nchar(fnames) > 0]
      if (is.null(frag_db) || nrow(frag_db) == 0 || length(fnames) == 0)
        return(empty)

      matched <- frag_db %>% filter(ExperimentID %in% fnames)
      if (nrow(matched) == 0) return(empty)

      list(has_data  = TRUE,
           n_entries = nrow(matched),
           n_genes   = n_distinct(matched$Gene),
           n_mods    = n_distinct(matched$Modification),
           user_ids  = unique(matched$UserID),
           file_ids  = unique(matched$ExperimentID))
    })
    names(stats) <- exp_db[["Experiment ID"]]
    stats
  })

  get_exp_stats <- function(exp_id) {
    s <- all_exp_stats()
    if (exp_id %in% names(s)) s[[exp_id]]
    else list(has_data = FALSE, n_entries = 0, n_genes = 0, n_mods = 0,
              user_ids = character(0), file_ids = character(0))
  }

  # ===========================================================================
  # CARD VIEW — compact experiment list
  # ===========================================================================
  output$fpdb_card_grid <- renderUI({
    # Show nothing until a database is loaded.
    if (!isTRUE(rv$fpdb_db_loaded)) return(NULL)
    fpdb$db_version
    exp_db   <- fpdb$exp_db
    selected <- fpdb$selected_exp

    if (is.null(exp_db) || nrow(exp_db) == 0) {
      return(tags$div(class = "text-center text-muted py-5",
        icon("database", class = "fa-3x mb-3"), tags$br(),
        tags$h5("No experiments found."),
        tags$p("Click ", tags$b("Add Experiment"), " to get started.")
      ))
    }

    stats_all <- all_exp_stats()

    rows <- lapply(seq_len(nrow(exp_db)), function(i) {
      row    <- exp_db[i, ]
      exp_id <- as.character(row[["Experiment ID"]])
      stats  <- if (exp_id %in% names(stats_all)) stats_all[[exp_id]]
                else list(has_data = FALSE, n_entries = 0)
      is_sel <- identical(exp_id, selected)

      status_icon <- if (stats$has_data) {
        tags$span(class = "text-success me-2", icon("circle-check"))
      } else {
        tags$span(class = "text-warning me-2", icon("circle-exclamation"))
      }

      data_label <- if (stats$has_data) {
        tags$span(class = "badge bg-success",
          paste0(stats$n_entries, " entries"))
      } else {
        tags$span(class = "badge bg-secondary", "No data")
      }

      organism <- if ("Organism" %in% colnames(row))
        as.character(row[["Organism"]]) else "Human"

      tags$div(
        class = paste0(
          "list-group-item list-group-item-action py-2 px-3 fpdb-list-item",
          if (is_sel) " fpdb-list-active" else ""
        ),
        style = "cursor: pointer; border-left: 3px solid transparent;",
        onclick = paste0("Shiny.setInputValue('fpdb_select_exp','", exp_id,
                         "',{priority:'event'})"),
        # Line 1: status icon + experiment ID
        tags$div(class = "d-flex align-items-center",
          status_icon,
          tags$span(class = "fw-bold", style = "font-size: 0.92em;", exp_id)
        ),
        # Line 2: metadata summary + data badge
        tags$div(class = "d-flex justify-content-between align-items-center mt-1",
          tags$span(class = "text-muted", style = "font-size: 0.78em;",
            paste0(row[["Lab"]], " / ", row[["Investigator"]],
                   " \u00b7 ", row[["Experiment Type"]],
                   " \u00b7 ", organism,
                   " \u00b7 ", row[["Experiment Year"]])),
          data_label
        )
      )
    })

    tags$div(class = "list-group list-group-flush border rounded", rows)
  })

  # ===========================================================================
  # DETAIL PANEL
  # ===========================================================================
  observeEvent(input$fpdb_select_exp, {
    fpdb$selected_exp <- input$fpdb_select_exp
    fpdb$detail_mode  <- "view"
  })

  observeEvent(input$fpdb_close_detail, {
    fpdb$selected_exp <- NULL
    fpdb$detail_mode  <- "view"
  })

  observeEvent(input$fpdb_edit_exp,   { fpdb$detail_mode <- "edit" })
  observeEvent(input$fpdb_cancel_edit, { fpdb$detail_mode <- "view" })

  # ---- Save edits ----
  observeEvent(input$fpdb_save_edit, {
    req(fpdb$exp_db, fpdb$selected_exp)
    exp_id <- fpdb$selected_exp
    idx <- which(fpdb$exp_db[["Experiment ID"]] == exp_id)
    if (length(idx) == 0) return()
    idx <- idx[1]

    fpdb$exp_db[idx, "Lab"]                 <- input$fpdb_edit_lab
    fpdb$exp_db[idx, "Investigator"]        <- input$fpdb_edit_investigator
    fpdb$exp_db[idx, "Experiment Year"]     <- input$fpdb_edit_year
    fpdb$exp_db[idx, "Cell/Tissue type"]    <- input$fpdb_edit_celltype
    fpdb$exp_db[idx, "Experiment Type"]     <- input$fpdb_edit_exptype
    fpdb$exp_db[idx, "Experimental Design"] <- input$fpdb_edit_design
    fpdb$exp_db[idx, "Samples"]             <- input$fpdb_edit_samples
    fpdb$exp_db[idx, "Notes"]               <- input$fpdb_edit_notes
    fpdb$exp_db[idx, "Order of conditions"] <- input$fpdb_edit_conditions
    fpdb$exp_db[idx, "Enzyme"]              <- input$fpdb_edit_enzyme
    fpdb$exp_db[idx, "DIA or DDA"]          <- input$fpdb_edit_diadda
    fpdb$exp_db[idx, "Organism"]            <- input$fpdb_edit_organism

    rv$fpdb_experiment_db <- fpdb$exp_db
    fpdb$db_version  <- fpdb$db_version + 1L
    fpdb$detail_mode <- "view"
    logger(paste0("FragPipe DB: Updated metadata for ", exp_id))
    showNotification(paste0("Updated ", exp_id), type = "message")
  })

  # ---- Delete experiment ----
  observeEvent(input$fpdb_delete_exp, {
    req(fpdb$selected_exp)
    exp_id <- fpdb$selected_exp
    showModal(modalDialog(
      title = tags$span(class = "text-danger",
        icon("triangle-exclamation"), " Delete Experiment"),
      tags$div(
        tags$p("Are you sure you want to delete ",
               tags$strong(tags$code(exp_id)), "?"),
        tags$div(class = "alert alert-danger py-2",
          icon("circle-exclamation"),
          tags$strong(" Warning: "),
          "This will delete ", tags$b("BOTH"), " information and data ",
          "associated with this experiment. This action cannot be undone."
        )
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("fpdb_confirm_delete", "Delete Permanently",
                     class = "btn-danger", icon = icon("trash"))
      )
    ))
  })

  observeEvent(input$fpdb_confirm_delete, {
    req(fpdb$selected_exp)
    exp_id <- fpdb$selected_exp
    exp_db <- fpdb$exp_db
    idx <- which(exp_db[["Experiment ID"]] == exp_id)
    if (length(idx) == 0) { removeModal(); return() }

    row <- exp_db[idx[1], , drop = FALSE]
    fnames_str <- as.character(row[["Order of Filenames"]])
    fnames <- trimws(unlist(strsplit(fnames_str, "\n")))
    fnames <- fnames[nchar(fnames) > 0]

    # 1. Update experiment DB in session cache
    exp_db <- exp_db[-idx, , drop = FALSE]
    rv$fpdb_experiment_db <- exp_db

    # 2. Remove matching rows from FragpipeDB session cache
    if (!is.null(rv$fpdb_fragpipe_db) && nrow(rv$fpdb_fragpipe_db) > 0 && length(fnames) > 0) {
      n_before <- nrow(rv$fpdb_fragpipe_db)
      rv$fpdb_fragpipe_db <- rv$fpdb_fragpipe_db %>% filter(!ExperimentID %in% fnames)
      logger(paste0("FragPipe DB: Removed ", n_before - nrow(rv$fpdb_fragpipe_db),
                    " data entries for ", exp_id))
    }

    fpdb$selected_exp <- NULL
    fpdb$detail_mode  <- "view"
    fpdb$db_version   <- fpdb$db_version + 1L
    logger(paste0("FragPipe DB: Deleted experiment '", exp_id, "'"))
    showNotification(paste0("Experiment '", exp_id, "' deleted."),
                     type = "warning", duration = 8)
    removeModal()
  })

  # ---- Detail panel render ----
  output$fpdb_detail_panel <- renderUI({
    if (!isTRUE(rv$fpdb_db_loaded)) return(NULL)
    exp_id <- fpdb$selected_exp
    if (is.null(exp_id)) {
      return(tags$div(class = "text-center text-muted py-5",
        icon("hand-pointer", class = "fa-2x mb-3"), tags$br(),
        tags$p("Click on an experiment to view its details.")
      ))
    }

    exp_db <- fpdb$exp_db
    if (is.null(exp_db)) return(NULL)
    idx <- which(exp_db[["Experiment ID"]] == exp_id)
    if (length(idx) == 0) return(tags$p(class = "text-danger", "Experiment not found."))
    row <- exp_db[idx[1], , drop = FALSE]

    fnames_str <- as.character(row[["Order of Filenames"]])
    stats <- get_exp_stats(exp_id)

    conditions <- trimws(unlist(strsplit(as.character(row[["Order of conditions"]]), "\n")))
    filenames  <- trimws(unlist(strsplit(fnames_str, "\n")))
    conditions <- conditions[nchar(conditions) > 0]
    filenames  <- filenames[nchar(filenames) > 0]

    max_len <- max(length(conditions), length(filenames))
    if (length(conditions) < max_len)
      conditions <- c(conditions, rep("", max_len - length(conditions)))
    if (length(filenames) < max_len)
      filenames <- c(filenames, rep("", max_len - length(filenames)))

    organism <- if ("Organism" %in% colnames(row))
      as.character(row[["Organism"]]) else "Human"

    if (fpdb$detail_mode == "edit") {
      # ===== EDIT MODE =====
      return(tags$div(
        tags$div(class = "d-flex justify-content-between align-items-center mb-3",
          tags$h5(class = "mb-0", icon("pen-to-square"), " Editing: ", tags$code(exp_id)),
          tags$div(
            actionButton("fpdb_save_edit", "Save", class = "btn-success btn-sm me-2",
                         icon = icon("floppy-disk")),
            actionButton("fpdb_cancel_edit", "Cancel",
                         class = "btn-outline-secondary btn-sm", icon = icon("xmark"))
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          textInput("fpdb_edit_lab", "Lab", value = as.character(row[["Lab"]])),
          textInput("fpdb_edit_investigator", "Investigator",
                    value = as.character(row[["Investigator"]]))
        ),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          textInput("fpdb_edit_year", "Year",
                    value = as.character(row[["Experiment Year"]])),
          textInput("fpdb_edit_celltype", "Cell/Tissue Type",
                    value = as.character(row[["Cell/Tissue type"]])),
          selectInput("fpdb_edit_exptype", "Experiment Type",
                      choices = c("IP/MS", "LC/MS", "TMT", "DIA", "Other"),
                      selected = as.character(row[["Experiment Type"]])),
          selectizeInput("fpdb_edit_organism", "Organism",
                         choices = c("Human", "Mouse", "Rat"),
                         selected = organism,
                         options = list(create = TRUE))
        ),
        textAreaInput("fpdb_edit_design", "Experimental Design",
                      value = as.character(row[["Experimental Design"]]),
                      rows = 3, width = "100%"),
        layout_columns(
          col_widths = c(4, 4, 4),
          numericInput("fpdb_edit_samples", "Samples",
                       value = as.integer(row[["Samples"]]), min = 1),
          textInput("fpdb_edit_enzyme", "Enzyme",
                    value = as.character(row[["Enzyme"]])),
          selectInput("fpdb_edit_diadda", "DIA or DDA",
                      choices = c("DDA", "DIA"),
                      selected = as.character(row[["DIA or DDA"]]))
        ),
        textAreaInput("fpdb_edit_notes", "Notes",
                      value = as.character(row[["Notes"]]),
                      rows = 2, width = "100%"),
        textAreaInput("fpdb_edit_conditions", "Order of Conditions (one per line)",
                      value = as.character(row[["Order of conditions"]]),
                      rows = max(2, length(conditions)), width = "100%"),
        tags$div(class = "alert alert-secondary mt-2",
          tags$strong(icon("lock"), " Filenames (read-only):"),
          tags$pre(class = "mb-0 mt-1", style = "font-size: 0.85em;", fnames_str)
        )
      ))
    }

    # ===== VIEW MODE =====
    mapping_table <- tags$table(
      class = "table table-sm table-bordered mb-0",
      style = "font-size: 0.85em;",
      tags$thead(
        tags$tr(
          tags$th(style = "width: 5%;", "#"),
          tags$th(style = "width: 47.5%;", "Condition"),
          tags$th(style = "width: 47.5%;", "Filename")
        )
      ),
      tags$tbody(
        lapply(seq_len(max_len), function(j) {
          fname <- filenames[j]
          has_file_data <- fname != "" && stats$has_data && fname %in% stats$file_ids
          file_icon <- if (has_file_data) {
            tags$span(class = "text-success", icon("circle-check"))
          } else if (fname != "") {
            tags$span(class = "text-warning", icon("circle-xmark"))
          } else NULL

          tags$tr(
            tags$td(class = "text-center text-muted", j),
            tags$td(conditions[j]),
            tags$td(tags$span(file_icon, " ", tags$code(fname)))
          )
        })
      )
    )

    data_status_card <- if (stats$has_data) {
      tags$div(class = "alert alert-success py-2 px-3 mb-3",
        tags$div(class = "d-flex justify-content-between",
          tags$div(
            icon("database"), tags$strong(" FragpipeDB Data Available"),
            tags$ul(class = "mb-0 mt-1 small",
              tags$li(paste0("Entries: ", stats$n_entries)),
              tags$li(paste0("Unique genes: ", stats$n_genes)),
              tags$li(paste0("Unique modifications: ", stats$n_mods)),
              tags$li(paste0("User IDs: ", paste(stats$user_ids, collapse = ", ")))
            )
          )
        )
      )
    } else {
      tags$div(class = "alert alert-warning py-2 px-3 mb-3",
        icon("triangle-exclamation"),
        tags$strong(" No FragpipeDB data found for this experiment."),
        tags$br(),
        tags$span(class = "small",
          "You can upload PSM and peptide files to populate the database."),
        tags$div(class = "mt-2",
          actionButton("fpdb_upload_for_exp", "Upload Experiment Data",
                       class = "btn-primary btn-sm", icon = icon("upload"))
        )
      )
    }

    tags$div(
      # Header with Edit / Delete / Close buttons
      tags$div(class = "d-flex justify-content-between align-items-center mb-3",
        tags$h5(class = "mb-0", icon("flask"), " ", tags$code(exp_id)),
        tags$div(
          actionButton("fpdb_edit_exp", "Edit",
                       class = "btn-outline-primary btn-sm me-2",
                       icon = icon("pen-to-square")),
          actionButton("fpdb_delete_exp", "Delete",
                       class = "btn-outline-danger btn-sm me-2",
                       icon = icon("trash")),
          actionButton("fpdb_close_detail", "Close",
                       class = "btn-outline-secondary btn-sm",
                       icon = icon("xmark"))
        )
      ),

      # Meta info (with Organism)
      tags$div(class = "row mb-3",
        tags$div(class = "col-md-6",
          tags$table(class = "table table-sm table-borderless mb-0",
                     style = "font-size: 0.9em;",
            tags$tr(tags$td(class = "fw-bold text-muted", "Lab"),
                    tags$td(as.character(row[["Lab"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "Investigator"),
                    tags$td(as.character(row[["Investigator"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "Year"),
                    tags$td(as.character(row[["Experiment Year"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "Cell/Tissue"),
                    tags$td(as.character(row[["Cell/Tissue type"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "Organism"),
                    tags$td(organism))
          )
        ),
        tags$div(class = "col-md-6",
          tags$table(class = "table table-sm table-borderless mb-0",
                     style = "font-size: 0.9em;",
            tags$tr(tags$td(class = "fw-bold text-muted", "Type"),
                    tags$td(as.character(row[["Experiment Type"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "Samples"),
                    tags$td(as.character(row[["Samples"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "Enzyme"),
                    tags$td(as.character(row[["Enzyme"]]))),
            tags$tr(tags$td(class = "fw-bold text-muted", "DIA/DDA"),
                    tags$td(as.character(row[["DIA or DDA"]])))
          )
        )
      ),

      # Experimental Design
      tags$div(class = "mb-3",
        tags$label(class = "fw-bold text-muted small", "Experimental Design"),
        tags$div(class = "border rounded p-2",
                 style = "background: #f8f9fa; font-size: 0.9em; white-space: pre-wrap;",
          as.character(row[["Experimental Design"]]))
      ),

      # Notes
      if (nchar(trimws(as.character(row[["Notes"]]))) > 0) {
        tags$div(class = "mb-3",
          tags$label(class = "fw-bold text-muted small", "Notes"),
          tags$div(class = "border rounded p-2",
                   style = "background: #fff3cd; font-size: 0.85em; white-space: pre-wrap;",
            as.character(row[["Notes"]])))
      },

      # Condition / Filename mapping
      tags$div(class = "mb-3",
        tags$label(class = "fw-bold text-muted small", "Sample Mapping"),
        mapping_table
      ),

      # FragpipeDB status
      data_status_card
    )
  })

  # ===========================================================================
  # Upload data for an existing experiment (no FragpipeDB data)
  # ===========================================================================
  observeEvent(input$fpdb_upload_for_exp, {
    req(fpdb$selected_exp)
    exp_id <- fpdb$selected_exp
    exp_db <- fpdb$exp_db
    idx <- which(exp_db[["Experiment ID"]] == exp_id)
    if (length(idx) == 0) return()
    row <- exp_db[idx[1], , drop = FALSE]

    fnames <- trimws(unlist(strsplit(as.character(row[["Order of Filenames"]]), "\n")))
    fnames <- fnames[nchar(fnames) > 0]
    conditions <- trimws(unlist(strsplit(as.character(row[["Order of conditions"]]), "\n")))
    conditions <- conditions[nchar(conditions) > 0]

    showModal(modalDialog(
      title = tags$span(icon("upload"), " Upload Data for ", tags$code(exp_id)),
      tagList(
        tags$div(class = "alert alert-info py-2 small",
          icon("circle-info"),
          " Upload the FragPipe output files for this experiment.",
          " Each pair of PSM.tsv + peptide.tsv corresponds to one filename/condition."
        ),
        tags$p(class = "fw-bold small", "Expected files:"),
        tags$ul(class = "small",
          lapply(seq_along(fnames), function(j) {
            cond_label <- if (j <= length(conditions)) conditions[j] else ""
            tags$li(tags$code(fnames[j]),
                    if (nchar(cond_label) > 0) paste0(" \u2014 ", cond_label))
          })
        ),
        hr(),
        textInput("fpdb_upload_user_id", "User ID / Label",
                  placeholder = "e.g. CellLine_Condition"),
        fileInput("fpdb_upload_psm", "PSM.tsv file", accept = c(".tsv", ".txt")),
        fileInput("fpdb_upload_peptide", "peptide.tsv file", accept = c(".tsv", ".txt")),
        tags$p(class = "text-muted small",
          icon("circle-info"),
          " Import one file pair at a time. Repeat for each condition/filename.")
      ),
      size = "l", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("fpdb_upload_run", "Import", class = "btn-primary",
                     icon = icon("upload"))
      )
    ))
  })

  observeEvent(input$fpdb_upload_run, {
    if (is.null(input$fpdb_upload_psm)) {
      showNotification("Please select a PSM.tsv file.", type = "error"); return()
    }
    if (is.null(input$fpdb_upload_peptide)) {
      showNotification("Please select a peptide.tsv file.", type = "error"); return()
    }
    user_id <- trimws(input$fpdb_upload_user_id)
    if (nchar(user_id) == 0) {
      showNotification("Please enter a User ID.", type = "error"); return()
    }

    shinyjs::disable("fpdb_upload_run")

    withProgress(message = "Importing FragPipe data...", value = 0, {
      incProgress(0.1, detail = "Reading files...")
      result <- tryCatch({
        update_fragpipe_db(
          psm_path     = input$fpdb_upload_psm$datapath,
          peptide_path = input$fpdb_upload_peptide$datapath,
          user_id      = user_id,
          existing_db  = rv$fpdb_fragpipe_db,
          progress_fn  = function(msg) incProgress(0, detail = msg)
        )
      }, error = function(e) {
        logger(paste0("!!! [FragPipeDB] Import error (existing exp): ", e$message))
        list(success = FALSE, message = paste0("Error: ", e$message), n_entries = 0)
      })
      incProgress(0.9, detail = "Done")
    })

    shinyjs::enable("fpdb_upload_run")

    if (result$success) {
      rv$fpdb_fragpipe_db <- if (is.null(rv$fpdb_fragpipe_db)) result$data
                             else dplyr::bind_rows(rv$fpdb_fragpipe_db, result$data)
      fpdb$db_version <- fpdb$db_version + 1L
      logger(paste0("FragPipe DB: Imported ", result$n_entries,
                    " entries for experiment ", fpdb$selected_exp))
      showNotification(result$message, type = "message", duration = 8)
      removeModal()
    } else {
      showNotification(result$message, type = "error", duration = 10)
    }
  })

  # ===========================================================================
  # TABLE VIEW
  # ===========================================================================
  output$fpdb_exp_table <- renderDT({
    if (!isTRUE(rv$fpdb_db_loaded)) return(NULL)
    fpdb$db_version
    exp_db <- fpdb$exp_db
    if (is.null(exp_db) || nrow(exp_db) == 0) return(NULL)

    stats_all <- all_exp_stats()

    exp_db$`Data Status` <- sapply(seq_len(nrow(exp_db)), function(i) {
      eid <- as.character(exp_db[i, "Experiment ID"])
      s <- if (eid %in% names(stats_all)) stats_all[[eid]]
           else list(has_data = FALSE)
      if (s$has_data) paste0("Loaded (", s$n_entries, " entries)") else "No data"
    })

    display_cols <- c("Experiment ID", "Lab", "Investigator", "Experiment Year",
                      "Cell/Tissue type", "Experiment Type", "Organism",
                      "Samples", "Experimental Design", "Notes",
                      "Order of conditions", "Order of Filenames",
                      "Enzyme", "DIA or DDA", "Data Status")
    display_cols <- intersect(display_cols, colnames(exp_db))

    datatable(
      exp_db[, display_cols, drop = FALSE],
      filter = "top", rownames = FALSE,
      options = list(
        pageLength = 25, scrollX = TRUE,
        dom = "Blftipr",
        buttons = list(
          list(extend = "csv", text = "Download CSV",
               filename = "MS_Experiment_DB")
        ),
        columnDefs = list(
          list(targets = which(display_cols == "Experimental Design") - 1,
               render = JS("function(data) {
                 if (data && data.length > 80) return data.substr(0, 80) + '...';
                 return data;
               }"))
        )
      ),
      extensions = "Buttons",
      class = "compact stripe"
    )
  })

  # ===========================================================================
  # ADD EXPERIMENT — modal
  # ===========================================================================
  observeEvent(input$fpdb_add_exp_btn, {
    showModal(modalDialog(
      title = tags$span(icon("plus-circle"), " Add New Experiment"),
      tagList(
        tags$p(class = "text-muted small",
          "Choose whether to add experiment documentation with data ",
          "(PSM + peptide files), or documentation only."
        ),
        tags$div(class = "d-grid gap-3",
          actionButton("fpdb_add_with_data",
                       tags$span(icon("upload"), " Add Experiment with Data"),
                       class = "btn-primary btn-lg w-100"),
          tags$p(class = "small text-muted text-center mb-0",
            "Provide experiment metadata AND upload PSM.tsv + peptide.tsv files."),
          hr(),
          actionButton("fpdb_add_doc_only",
                       tags$span(icon("file-lines"), " Add Documentation Only"),
                       class = "btn-outline-secondary btn-lg w-100"),
          tags$p(class = "small text-muted text-center mb-0",
            "Record experiment metadata without uploading any data files.")
        )
      ),
      size = "m", easyClose = TRUE,
      footer = modalButton("Cancel")
    ))
  })

  # ---- Add with data form ----
  observeEvent(input$fpdb_add_with_data, {
    removeModal()
    showModal(modalDialog(
      title = tags$span(icon("flask"), " New Experiment (with Data)"),
      tagList(
        tags$h6(icon("clipboard"), " Experiment Information"),
        layout_columns(
          col_widths = c(6, 6),
          textInput("fpdb_new_lab", "Lab", placeholder = "e.g. Coulombe"),
          textInput("fpdb_new_investigator", "Investigator", placeholder = "e.g. EC")
        ),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          textInput("fpdb_new_year", "Year", value = format(Sys.Date(), "%Y")),
          textInput("fpdb_new_expid", "Experiment ID",
                    placeholder = "e.g. PRF_Q_2024_P_LAB_123"),
          textInput("fpdb_new_celltype", "Cell/Tissue Type",
                    placeholder = "e.g. A431"),
          selectizeInput("fpdb_new_organism", "Organism",
                         choices = c("Human", "Mouse", "Rat"),
                         selected = "Human",
                         options = list(create = TRUE))
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          selectInput("fpdb_new_exptype", "Experiment Type",
                      choices = c("IP/MS", "LC/MS", "TMT", "DIA", "Other")),
          numericInput("fpdb_new_samples", "Samples", value = 2, min = 1),
          selectInput("fpdb_new_diadda", "DIA or DDA",
                      choices = c("DDA", "DIA"))
        ),
        textAreaInput("fpdb_new_design", "Experimental Design",
                      rows = 3, width = "100%",
                      placeholder = "Describe the experimental setup..."),
        textAreaInput("fpdb_new_notes", "Notes", rows = 2, width = "100%",
                      placeholder = "Any relevant notes..."),
        textInput("fpdb_new_enzyme", "Enzyme", value = "Trypsin"),
        textAreaInput("fpdb_new_conditions", "Order of Conditions (one per line)",
                      rows = 3, width = "100%",
                      placeholder = "Condition 1\nCondition 2\n..."),
        textAreaInput("fpdb_new_filenames", "Order of Filenames (one per line)",
                      rows = 3, width = "100%",
                      placeholder = "PRF_Q_2024_P_LAB_123_00001\n..."),
        hr(),
        tags$h6(icon("upload"), " FragPipe Data Files"),
        textInput("fpdb_new_userid", "User ID / Label",
                  placeholder = "e.g. CellLine_Condition"),
        fileInput("fpdb_new_psm", "PSM.tsv file", accept = c(".tsv", ".txt")),
        fileInput("fpdb_new_peptide", "peptide.tsv file", accept = c(".tsv", ".txt"))
      ),
      size = "l", easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("fpdb_new_submit_data", "Save & Import",
                     class = "btn-success", icon = icon("floppy-disk"))
      )
    ))
  })

  observeEvent(input$fpdb_new_submit_data, {
    exp_id <- trimws(input$fpdb_new_expid)
    if (nchar(exp_id) == 0) {
      showNotification("Experiment ID is required.", type = "error"); return()
    }
    if (!is.null(fpdb$exp_db) && exp_id %in% fpdb$exp_db[["Experiment ID"]]) {
      showNotification("This Experiment ID already exists.", type = "error"); return()
    }
    if (is.null(input$fpdb_new_psm)) {
      showNotification("Please select a PSM.tsv file.", type = "error"); return()
    }
    if (is.null(input$fpdb_new_peptide)) {
      showNotification("Please select a peptide.tsv file.", type = "error"); return()
    }
    user_id <- trimws(input$fpdb_new_userid)
    if (nchar(user_id) == 0) {
      showNotification("Please enter a User ID.", type = "error"); return()
    }

    shinyjs::disable("fpdb_new_submit_data")

    new_row <- data.frame(
      Lab                  = trimws(input$fpdb_new_lab),
      Investigator         = trimws(input$fpdb_new_investigator),
      `Experiment Year`    = trimws(input$fpdb_new_year),
      `Experiment ID`      = exp_id,
      `Cell/Tissue type`   = trimws(input$fpdb_new_celltype),
      `Experiment Type`    = input$fpdb_new_exptype,
      `Experimental Design`= trimws(input$fpdb_new_design),
      Samples              = input$fpdb_new_samples,
      Notes                = trimws(input$fpdb_new_notes),
      `Order of conditions`= trimws(input$fpdb_new_conditions),
      `Order of Filenames` = trimws(input$fpdb_new_filenames),
      Enzyme               = trimws(input$fpdb_new_enzyme),
      `DIA or DDA`         = input$fpdb_new_diadda,
      Organism             = input$fpdb_new_organism,
      stringsAsFactors = FALSE, check.names = FALSE
    )

    combined <- if (is.null(fpdb$exp_db)) new_row
                else rbind(fpdb$exp_db, new_row)
    rv$fpdb_experiment_db <- combined

    withProgress(message = "Importing FragPipe data...", value = 0, {
      incProgress(0.1, detail = "Reading files...")
      result <- tryCatch({
        update_fragpipe_db(
          psm_path     = input$fpdb_new_psm$datapath,
          peptide_path = input$fpdb_new_peptide$datapath,
          user_id      = user_id,
          existing_db  = rv$fpdb_fragpipe_db,
          progress_fn  = function(msg) incProgress(0, detail = msg)
        )
      }, error = function(e) {
        logger(paste0("!!! [FragPipeDB] Import error (new exp): ", e$message))
        list(success = FALSE, message = paste0("Error: ", e$message), n_entries = 0)
      })
      incProgress(0.9, detail = "Done")
    })

    shinyjs::enable("fpdb_new_submit_data")

    if (result$success) {
      rv$fpdb_fragpipe_db <- if (is.null(rv$fpdb_fragpipe_db)) result$data
                             else dplyr::bind_rows(rv$fpdb_fragpipe_db, result$data)
      fpdb$db_version <- fpdb$db_version + 1L
      logger(paste0("FragPipe DB: Added experiment '", exp_id,
                    "' with ", result$n_entries, " data entries."))
      showNotification(
        paste0("Experiment added! ", result$n_entries, " entries imported."),
        type = "message", duration = 8)
      removeModal()
    } else {
      fpdb$db_version <- fpdb$db_version + 1L
      showNotification(
        paste0("Metadata saved, but data import failed: ", result$message),
        type = "warning", duration = 10)
    }
  })

  # ---- Add documentation only form ----
  observeEvent(input$fpdb_add_doc_only, {
    removeModal()
    showModal(modalDialog(
      title = tags$span(icon("file-lines"), " New Experiment (Documentation Only)"),
      tagList(
        tags$div(class = "alert alert-info py-2 small",
          icon("circle-info"),
          " You are adding experiment documentation ", tags$b("without"),
          " data files. You can upload data later from the experiment detail view."
        ),
        tags$h6(icon("clipboard"), " Experiment Information"),
        layout_columns(
          col_widths = c(6, 6),
          textInput("fpdb_doc_lab", "Lab", placeholder = "e.g. Coulombe"),
          textInput("fpdb_doc_investigator", "Investigator", placeholder = "e.g. EC")
        ),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          textInput("fpdb_doc_year", "Year", value = format(Sys.Date(), "%Y")),
          textInput("fpdb_doc_expid", "Experiment ID",
                    placeholder = "e.g. PRF_Q_2024_P_LAB_123"),
          textInput("fpdb_doc_celltype", "Cell/Tissue Type",
                    placeholder = "e.g. A431"),
          selectizeInput("fpdb_doc_organism", "Organism",
                         choices = c("Human", "Mouse", "Rat"),
                         selected = "Human",
                         options = list(create = TRUE))
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          selectInput("fpdb_doc_exptype", "Experiment Type",
                      choices = c("IP/MS", "LC/MS", "TMT", "DIA", "Other")),
          numericInput("fpdb_doc_samples", "Samples", value = 2, min = 1),
          selectInput("fpdb_doc_diadda", "DIA or DDA",
                      choices = c("DDA", "DIA"))
        ),
        textAreaInput("fpdb_doc_design", "Experimental Design",
                      rows = 3, width = "100%",
                      placeholder = "Describe the experimental setup..."),
        textAreaInput("fpdb_doc_notes", "Notes", rows = 2, width = "100%",
                      placeholder = "Any relevant notes..."),
        textInput("fpdb_doc_enzyme", "Enzyme", value = "Trypsin"),
        textAreaInput("fpdb_doc_conditions", "Order of Conditions (one per line)",
                      rows = 3, width = "100%",
                      placeholder = "Condition 1\nCondition 2\n..."),
        textAreaInput("fpdb_doc_filenames", "Order of Filenames (one per line)",
                      rows = 3, width = "100%",
                      placeholder = "PRF_Q_2024_P_LAB_123_00001\n..."),
        hr(),
        textAreaInput("fpdb_doc_reason", "Reason for no data upload",
                      rows = 2, width = "100%",
                      placeholder = "e.g. Data files not yet available, awaiting core facility...")
      ),
      size = "l", easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("fpdb_doc_submit", "Save Documentation",
                     class = "btn-success", icon = icon("floppy-disk"))
      )
    ))
  })

  observeEvent(input$fpdb_doc_submit, {
    exp_id <- trimws(input$fpdb_doc_expid)
    if (nchar(exp_id) == 0) {
      showNotification("Experiment ID is required.", type = "error"); return()
    }
    if (!is.null(fpdb$exp_db) && exp_id %in% fpdb$exp_db[["Experiment ID"]]) {
      showNotification("This Experiment ID already exists.", type = "error"); return()
    }
    reason <- trimws(input$fpdb_doc_reason)
    if (nchar(reason) == 0) {
      showNotification("Please provide a reason for not uploading data.",
                       type = "error"); return()
    }

    notes_text <- trimws(input$fpdb_doc_notes)
    notes_text <- if (nchar(notes_text) > 0) paste0(notes_text, " | NO DATA: ", reason)
                  else paste0("NO DATA: ", reason)

    new_row <- data.frame(
      Lab                  = trimws(input$fpdb_doc_lab),
      Investigator         = trimws(input$fpdb_doc_investigator),
      `Experiment Year`    = trimws(input$fpdb_doc_year),
      `Experiment ID`      = exp_id,
      `Cell/Tissue type`   = trimws(input$fpdb_doc_celltype),
      `Experiment Type`    = input$fpdb_doc_exptype,
      `Experimental Design`= trimws(input$fpdb_doc_design),
      Samples              = input$fpdb_doc_samples,
      Notes                = notes_text,
      `Order of conditions`= trimws(input$fpdb_doc_conditions),
      `Order of Filenames` = trimws(input$fpdb_doc_filenames),
      Enzyme               = trimws(input$fpdb_doc_enzyme),
      `DIA or DDA`         = input$fpdb_doc_diadda,
      Organism             = input$fpdb_doc_organism,
      stringsAsFactors = FALSE, check.names = FALSE
    )

    combined <- if (is.null(fpdb$exp_db)) new_row
                else rbind(fpdb$exp_db, new_row)
    rv$fpdb_experiment_db <- combined

    fpdb$db_version <- fpdb$db_version + 1L
    logger(paste0("FragPipe DB: Added documentation for '", exp_id, "' (no data)"))
    showNotification(paste0("Experiment '", exp_id, "' documented."),
                     type = "message", duration = 6)
    removeModal()
  })

  # ===========================================================================
  # Database summary (uses cached FragpipeDB)
  # ===========================================================================
  output$fpdb_db_summary <- renderPrint({
    fpdb$db_version
    if (!isTRUE(rv$fpdb_db_loaded)) { cat("No database loaded.\n"); return(invisible(NULL)) }
    n_exp <- if (!is.null(fpdb$exp_db)) nrow(fpdb$exp_db) else 0
    cat("Experiment Metadata: ", n_exp, " experiments documented\n")

    frag_db <- frag_db_cache()
    if (!is.null(frag_db) && nrow(frag_db) > 0) {
      cat("FragpipeDB: ", nrow(frag_db), " entries, ",
          n_distinct(frag_db$Gene), " genes, ",
          n_distinct(frag_db$Modification), " modifications\n")
      if ("ExperimentID" %in% colnames(frag_db))
        cat("Data files loaded: ", n_distinct(frag_db$ExperimentID), "\n")
    } else {
      cat("FragpipeDB: No data file found.\n")
    }
  })

  # ===========================================================================
  # Export Database (downloads both files as a ZIP)
  # ===========================================================================
  output$fpdb_export_btn <- downloadHandler(
    filename = function() {
      paste0(format(Sys.Date(), "%Y%m%d"), "_FragpipeDB_Export.zip")
    },
    content = function(file) {
      req(rv$fpdb_db_loaded)
      tmpdir <- tempdir()
      date_prefix <- format(Sys.Date(), "%Y%m%d")
      fdb_path <- file.path(tmpdir, paste0(date_prefix, "_FragpipeDB.csv"))
      edb_path <- file.path(tmpdir, paste0(date_prefix, "_MS_Experiment_DB.csv"))
      readr::write_csv(if (!is.null(rv$fpdb_fragpipe_db)) rv$fpdb_fragpipe_db
                       else data.frame(), fdb_path, na = "")
      write.csv(if (!is.null(rv$fpdb_experiment_db)) rv$fpdb_experiment_db
                else data.frame(), edb_path, row.names = FALSE)
      zip(file, c(fdb_path, edb_path), flags = "-j")
    },
    contentType = "application/zip"
  )
}
