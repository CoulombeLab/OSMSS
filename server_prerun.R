# server_prerun.R ─ Prerun Experiments module
# ===========================================================================
# FOLDER-SCAN APPROACH (client-side, preview files only)
#
# The user clicks "Select Folder" which triggers a webkitdirectory input.
# JavaScript:
#   • filters the selection to *.preview.rds files and relays them to the
#     hidden Shiny fileInput (tiny uploads, grid renders instantly).
#   • stores ALL non-preview *.rds files in window.prerunFullRdsMap.
#   • sends the relative paths of each preview as prerun_preview_paths.
#
# When the user clicks "Load" on a tile (plain button, class prerun-load-btn):
#   • JS looks for the matching full .rds in prerunFullRdsMap.
#   • FOUND   → JS injects it into the hidden #prerun_auto_rds_input fileInput
#               and fires prerun_auto_load_request → R auto-loads without modal.
#   • NOT FOUND → JS fires prerun_manual_load_request → R shows a browse modal.
#
# Reactive values:
#   prerun_rv$previews          — lightweight preview objects
#   prerun_rv$preview_rel_paths — webkitRelativePath for each preview (from JS)
#   prerun_rv$pending_idx       — tile index awaiting the manual fallback upload
#   prerun_rv$pending_auto_idx  — tile index awaiting the auto-injected upload
# ===========================================================================

server_prerun <- function(input, output, session, rv, logger) {

  prerun_rv <- reactiveValues(
    previews          = list(),
    scan_error        = NULL,
    pending_idx       = NULL,
    pending_auto_idx  = NULL,
    preview_rel_paths = list()
  )

  valid_types <- c("lf_summarization", "dia_summarization",
                   "tmt_summarization", "ptm_summarization")

  # ── Folder status display ──────────────────────────────────────────────────
  output$prerun_folder_status <- renderUI({
    info <- input$prerun_folder_info
    if (is.null(info)) return(NULL)
    n        <- info$n        %||% 0L
    nFolders <- info$nFolders %||% 0L
    folder   <- info$folder   %||% ""
    tags$div(
      class = "small",
      icon("folder"), " ",
      tags$b(folder),
      tags$span(class = "text-muted ms-2",
        paste0(n, " preview", if (n != 1) "s" else "",
               if (nFolders > 1) paste0(" across ", nFolders, " sub-folders") else ""))
    )
  })

  # ── No previews found notification ────────────────────────────────────────
  observeEvent(input$prerun_no_previews, {
    showNotification(
      tagList(icon("triangle-exclamation"),
              " No .preview.rds files found in the selected folder. ",
              "Make sure you have extracted the export ZIPs first."),
      type = "warning", duration = 6)
  }, ignoreNULL = TRUE)

  # ── Capture relative paths sent by JS before the file upload completes ────
  observeEvent(input$prerun_preview_paths, {
    prerun_rv$preview_rel_paths <- input$prerun_preview_paths
  }, ignoreNULL = TRUE)

  # ── Process uploaded preview files ────────────────────────────────────────
  process_uploaded_files <- function(file_df) {
    if (is.null(file_df) || nrow(file_df) == 0) return()

    prerun_rv$scan_error <- NULL
    prerun_rv$previews   <- list()

    # Only process .preview.rds rows (JS pre-filters, but be defensive)
    preview_df <- file_df[grepl("\\.preview\\.rds$", file_df$name,
                                ignore.case = TRUE), , drop = FALSE]

    if (nrow(preview_df) == 0) {
      prerun_rv$scan_error <- "No .preview.rds files received."
      return()
    }

    previews <- list()
    rel_paths <- prerun_rv$preview_rel_paths %||% list()
    withProgress(message = "Reading experiment metadata...", value = 0, {
      for (i in seq_len(nrow(preview_df))) {
        incProgress(1 / nrow(preview_df), detail = preview_df$name[i])

        bundle <- tryCatch(readRDS(preview_df$datapath[i]), error = function(e) NULL)
        if (is.null(bundle) || !is.list(bundle))           next
        if (!isTRUE(bundle$type %in% valid_types))         next

        rel_path <- if (i <= length(rel_paths)) rel_paths[[i]] %||% "" else ""

        # Title: use relative path (e.g. "SAMIP/LF_Summarization_20260520_1412")
        # so the user can see which sub-folder the experiment lives in.
        display_name <- if (nzchar(rel_path)) {
          sub("\\.preview\\.rds$", "", rel_path, ignore.case = TRUE)
        } else {
          tools::file_path_sans_ext(
            sub("\\.preview$", "",
                tools::file_path_sans_ext(preview_df$name[i])))
        }

        # Expected full .rds filename (basename only — used for JS auto-lookup)
        expected_full <- sub("\\.preview\\.rds$", ".rds",
                             basename(if (nzchar(rel_path)) rel_path
                                      else preview_df$name[i]),
                             ignore.case = TRUE)

        preview <- list(
          type              = bundle$type,
          meta              = bundle$meta,
          route             = bundle$route,
          pca_thumbnail     = bundle$pca_thumbnail,
          n_features        = bundle$n_features,
          cond_groups       = bundle$cond_groups,
          conditions        = bundle$conditions,
          annotation_map    = bundle$annotation_map,
          data_format       = bundle$data_format,
          design_option     = bundle$design_option,
          .file_title       = display_name,
          .expected_full    = expected_full,
          .relative_path    = rel_path,
          .source_filename  = bundle$.source_filename %||% NULL
        )
        previews[[ length(previews) + 1 ]] <- preview
        rm(bundle)
      }
    })
    gc(verbose = FALSE)

    prerun_rv$previews <- previews
    if (length(previews) == 0)
      prerun_rv$scan_error <- "No valid summarization previews found."
  }

  # ── Observe file upload (triggered by JS) ─────────────────────────────────
  observeEvent(input$prerun_files, {
    process_uploaded_files(input$prerun_files)
  })

  # ── Refresh ────────────────────────────────────────────────────────────────
  observeEvent(input$prerun_refresh_btn, {
    process_uploaded_files(input$prerun_files)
  })

  # ── Render grid ────────────────────────────────────────────────────────────
  output$prerun_grid <- renderUI({
    previews <- prerun_rv$previews

    if (!is.null(prerun_rv$scan_error) && length(previews) == 0) {
      return(tags$div(class = "alert alert-warning mt-3",
                      icon("triangle-exclamation"), " ", prerun_rv$scan_error))
    }

    if (length(previews) == 0) {
      return(tags$div(
        class = "text-center mt-5 text-muted",
        style = "padding:60px 0;",
        tags$div(style = "font-size:3.5em; margin-bottom:16px; color:#dee2e6;",
                 icon("folder-open")),
        tags$h5("No folder selected"),
        tags$p(class = "small",
          "Click ", tags$b("Select Folder"), " and choose the folder where you ",
          "extracted your experiment ZIPs. The app will scan all sub-folders ",
          "for ", tags$code(".preview.rds"), " files automatically.")
      ))
    }

    tiles <- lapply(seq_along(previews), function(i) {
      .make_prerun_tile(previews[[i]], i)
    })

    tagList(
      tags$p(class = "text-muted small mb-3",
        tags$b(length(previews)), " experiment(s) \u2014 click ",
        tags$b("Load"), " on any tile to load its full .rds file."),
      layout_columns(col_widths = 4, !!!tiles)
    )
  })

  # ── Auto-load: JS found the full .rds and is injecting it ─────────────────
  # prerun_auto_load_request fires (fast, setInputValue) before the file upload
  # completes, so we store the pending index here and consume it when the file
  # arrives via prerun_auto_rds_input.
  observeEvent(input$prerun_auto_load_request, {
    prerun_rv$pending_auto_idx <- input$prerun_auto_load_request$idx
  }, ignoreNULL = TRUE)

  observeEvent(input$prerun_auto_rds_input, {
    idx <- prerun_rv$pending_auto_idx
    req(!is.null(idx))
    previews <- prerun_rv$previews
    if (idx > length(previews)) { prerun_rv$pending_auto_idx <- NULL; return() }
    prev <- previews[[idx]]

    fp <- input$prerun_auto_rds_input$datapath
    withProgress(message = paste0("Loading \u2018", basename(prev$.file_title), "\u2019\u2026"), value = 0.2, {
      bundle <- tryCatch(readRDS(fp), error = function(e) NULL)
      if (is.null(bundle)) {
        showNotification("Failed to read the .rds file.", type = "error")
        prerun_rv$pending_auto_idx <- NULL; return()
      }
      if (!identical(bundle$type, prev$type)) {
        showNotification(
          paste0("Wrong file type: expected '", prev$type,
                 "', got '", bundle$type %||% "unknown", "'."),
          type = "error")
        prerun_rv$pending_auto_idx <- NULL; return()
      }
      incProgress(0.6, detail = "Restoring session\u2026")
      .load_bundle_into_session(bundle, rv, session, input, logger)
      incProgress(1, detail = "Done")
    })
    prerun_rv$pending_auto_idx <- NULL
    logger(paste0("Prerun: auto-loaded \u2018", prev$.file_title, "\u2019"))
  })

  # ── Manual load: JS could not find the full .rds → show browse modal ──────
  observeEvent(input$prerun_manual_load_request, {
    idx <- input$prerun_manual_load_request$idx
    previews <- prerun_rv$previews
    if (is.null(idx) || idx > length(previews)) return()
    prev <- previews[[idx]]

    prerun_rv$pending_idx <- idx
    showModal(modalDialog(
      title = tagList(icon("file-arrow-up"), " Load Experiment"),
      tags$p(class = "small",
        "The full ", tags$code(".rds"), " file was not found automatically. ",
        "Please locate it manually:"),
      tags$p(class = "text-muted small mb-1",
        icon("folder"), " Preview path: ",
        tags$code(prev$.relative_path %||% prev$.file_title)),
      tags$p(class = "text-muted small mb-2",
        icon("file"), " Expected file: ",
        tags$code(prev$.expected_full)),
      fileInput("prerun_load_fallback_file",
                label  = NULL,
                accept = ".rds",
                width  = "100%",
                buttonLabel = tagList(icon("folder-open"), " Browse\u2026")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("prerun_load_fallback_go", "Load",
                     class = "btn-success", icon = icon("upload"))
      ),
      easyClose = TRUE, size = "m"
    ))
  }, ignoreNULL = TRUE)

  # ── Fallback load handler (manual modal confirm) ───────────────────────────
  observeEvent(input$prerun_load_fallback_go, {
    idx <- prerun_rv$pending_idx
    req(!is.null(idx), !is.null(input$prerun_load_fallback_file))

    previews <- prerun_rv$previews
    if (idx > length(previews)) { removeModal(); return() }
    prev <- previews[[idx]]

    fp <- input$prerun_load_fallback_file$datapath
    withProgress(message = paste0("Loading '", prev$.file_title, "'..."), value = 0.2, {
      bundle <- tryCatch(readRDS(fp), error = function(e) NULL)
      if (is.null(bundle)) {
        showNotification("Failed to read the .rds file.", type = "error")
        removeModal(); return()
      }
      if (!identical(bundle$type, prev$type)) {
        showNotification(
          paste0("Wrong file: expected type '", prev$type,
                 "', got '", bundle$type %||% "unknown", "'."),
          type = "error"); removeModal(); return()
      }
      incProgress(0.6, detail = "Restoring session...")
      .load_bundle_into_session(bundle, rv, session, input, logger)
      incProgress(1, detail = "Done")
    })
    prerun_rv$pending_idx <- NULL
    removeModal()
    logger(paste0("Prerun: loaded '", prev$.file_title, "'"))
  })

# ── helpers ────────────────────────────────────────────────────────────────
  .load_bundle_into_session <- function(bundle, rv, session, input, logger) {
    
    # =========================================================================
    # 0. RESET PREVIOUS SESSION STATE
    # =========================================================================
    # Wipe all analytical state back to factory defaults using our snapshot manager.
    # This prevents old dataset settings (like Bait proteins) from bleeding over!
    rv$reset_experiment_state()

    type <- bundle$type %||% "unknown"
    is_ptm <- FALSE
    
    # Check if this is an IP / SAINT pipeline experiment
    is_ip <- identical(bundle$meta$pipeline, "ip")

    # 1. Route the data to the correct module
    if (identical(type, "lf_summarization")) {
      rv$prerun_load_lf_bundle <- bundle
      updateTabsetPanel(session, "main_content_area", selected = "Upload Data")
      
    } else if (identical(type, "dia_summarization")) {
      rv$prerun_load_dia_bundle <- bundle
      updateTabsetPanel(session, "main_content_area", selected = "DIA Upload Data")
      
    } else if (identical(type, "tmt_summarization")) {
      rv$prerun_load_tmt_bundle <- bundle
      updateTabsetPanel(session, "main_content_area", selected = "TMT Upload Data")
      
    } else if (identical(type, "ptm_summarization")) {
      is_ptm <- TRUE
      if (identical(bundle$route, "dia_ptm")) {
        rv$prerun_load_dia_ptm_bundle <- bundle
        updateTabsetPanel(session, "main_content_area", selected = "Label Free PTM Upload")
      } else if (identical(bundle$route, "tmt_ptm")) {
        rv$prerun_load_tmt_ptm_bundle <- bundle
        updateTabsetPanel(session, "main_content_area", selected = "PTM Upload Data")
      } else {
        showNotification(
          tagList(icon("clock"), " Session loading for this PTM route is coming in a future update."),
          type = "warning", duration = 5)
        return() 
      }
    } else {
      showNotification("Unknown summarization type.", type = "error")
      return()
    }

    # =========================================================================
    # 2. THE LAUNCHPAD UX 
    # =========================================================================
    
    # COLLISION CHECK: If this is a SAINT experiment missing a bait, 
    # the system will ask for it. We must exit early here so we don't 
    # overwrite the Bait Selection modal!
    is_saint_missing_bait <- identical(bundle$meta$de_engine, "saint") &&
                             (is.null(bundle$meta$bait) || !nzchar(bundle$meta$bait))       
    if (is_saint_missing_bait) {
      return() # Exit silently. The Launchpad will fire after they type the bait.
    }
    
    if (is_ptm) {
      # --- PTM EXPERIMENT POPUP ---
      showModal(modalDialog(
        title = tagList(icon("circle-check", class = "text-success"), " PTM Experiment Loaded!"),
        tags$p("Your PTM data is successfully loaded. Where would you like to begin?"),
        tags$div(class = "d-grid gap-2",
          tags$button(class = "btn btn-primary btn-lg text-start",
            icon("braille"), " Explore Data: PTM PCA",
            onclick = "
              $('.modal').modal('hide'); 
              $('[data-tab=\"PTM PCA\"]').click();
              $('#submenu_discoverer, #submenu_pd_tmt_ptm, #submenu_ptm_explore').addClass('open').prev('.nav-group-btn').addClass('open');
            "
          ),
          tags$button(class = "btn btn-info btn-lg text-start text-white",
            icon("volcano"), " Analyze: PTM Differential Expression",
            onclick = "
              $('.modal').modal('hide'); 
              $('[data-tab=\"PTM DE\"]').click();
              $('#submenu_discoverer, #submenu_pd_tmt_ptm, #submenu_ptm_explore').addClass('open').prev('.nav-group-btn').addClass('open');
            "
          ),
          tags$button(class = "btn btn-outline-secondary text-start mt-2",
            icon("file-lines"), " View Data Summary",
            `data-bs-dismiss` = "modal"
          )
        ),
        easyClose = TRUE, footer = NULL
      ))
      
    } else {
      # --- PROTEIN EXPERIMENT POPUP ---
      showModal(modalDialog(
        title = tagList(icon("circle-check", class = "text-success"), " Protein Experiment Loaded!"),
        tags$p("Your Protein data is successfully loaded. Where would you like to begin?"),
        tags$div(class = "d-grid gap-2",
          tags$button(class = "btn btn-primary btn-lg text-start",
            icon("braille"), " Explore Data: PCA Analysis",
            onclick = "
              $('.modal').modal('hide'); 
              $('[data-tab=\"PCA Analysis\"]').click();
              $('#submenu_discoverer, #submenu_protein, #submenu_explore_data').addClass('open').prev('.nav-group-btn').addClass('open');
            "
          ),
          tags$button(class = "btn btn-info btn-lg text-start text-white",
            icon("chart-area"), " Analyze: Differential Expression",
            onclick = "
              $('.modal').modal('hide'); 
              $('[data-tab=\"Differential Expression\"]').click();
              $('#submenu_discoverer, #submenu_protein, #submenu_explore_data').addClass('open').prev('.nav-group-btn').addClass('open');
            "
          ),
          tags$button(class = "btn btn-outline-secondary text-start mt-2",
            icon("file-lines"), " View Data Summary",
            `data-bs-dismiss` = "modal"
          )
        ),
        easyClose = TRUE, footer = NULL
      ))
    }
  }

}  # end server_prerun

# ---------------------------------------------------------------------------
# .make_prerun_tile()
# ---------------------------------------------------------------------------
.make_prerun_tile <- function(prev, idx, has_full = FALSE) {
  type  <- prev$type %||% "unknown"
  title <- prev$.file_title %||% paste0("Experiment ", idx)

  type_cfg <- switch(type,
    lf_summarization  = list(bg = "primary",   label = "LF DDA"),
    dia_summarization = list(bg = "info",       label = "DIA"),
    tmt_summarization = list(bg = "success",    label = "TMT"),
    ptm_summarization = list(bg = "warning",    label = "PTM"),
    list(bg = "secondary", label = "Unknown")
  )

  n_feat     <- prev$n_features
  feat_label <- if (identical(type, "ptm_summarization")) "PTM sites" else "proteins"
  n_str      <- if (!is.null(n_feat) && !is.na(n_feat))
                  paste0(format(n_feat, big.mark = ","), " ", feat_label)
                else "\u2014"

  thumb    <- prev$pca_thumbnail
  thumb_ui <- if (!is.null(thumb) && nzchar(thumb)) {
    tags$img(src   = thumb,
             style = "width:100%; border-radius:4px; margin:6px 0; display:block;")
  } else {
    tags$div(
      style = paste0("width:100%; height:160px; display:flex; align-items:center;",
                     " justify-content:center; background:#f8f9fa; border-radius:4px;",
                     " margin:6px 0; color:#adb5bd; font-size:0.82em;",
                     " flex-direction:column; gap:6px;"),
      icon("chart-line"), tags$span("No preview available"))
  }

  can_load <- !identical(type, "ptm_summarization") ||
              identical(prev$route, "dia_ptm") ||
              identical(prev$route, "tmt_ptm")
  load_ui  <- if (can_load) {
    tags$button(
      class        = "btn btn-outline-success btn-sm w-100 prerun-load-btn",
      `data-idx`   = idx,
      tagList(icon("file-arrow-up"), " Load into Session\u2026")
    )
  } else {
    tags$div(class = "text-center small text-muted",
             icon("circle-info"), " This PTM bundle cannot be reloaded in this version")
  }

  card(
    class = "h-100",
    card_header(
      class = "py-2 px-3",
      tags$div(
        class = "d-flex align-items-start justify-content-between gap-1",
        tags$span(class = "fw-bold lh-sm",
                  style = "font-size:0.88em; word-break:break-word; flex:1;",
                  title),
        tags$span(class = paste0("badge text-white bg-", type_cfg$bg),
                  style = "font-size:0.7em; flex-shrink:0; margin-top:1px;",
                  type_cfg$label)
      )
    ),
    card_body(
      class = "p-2",
      thumb_ui,
      tags$div(
        class = "small mt-1",
        tags$div(class = "mb-1",
          tags$span(class = "fw-semibold text-muted", "Route: "),
          tags$span(.get_route_label(prev))),
        tags$div(class = "mb-1",
          tags$span(class = "fw-semibold text-muted", "Features: "),
          tags$span(n_str)),
        .make_groups_display(prev)
      )
    ),
    card_footer(class = "p-2", load_ui)
  )
}

# ---------------------------------------------------------------------------
# .get_route_label()
# ---------------------------------------------------------------------------
.get_route_label <- function(prev) {
  type <- prev$type %||% ""
  if (identical(type, "lf_summarization")) {
    prev$meta$pipeline_label %||% prev$meta$pipeline %||% "Label-Free DDA"
  } else if (identical(type, "dia_summarization")) {
    prev$meta$pipeline_label %||% "Label-Free DIA"
  } else if (identical(type, "tmt_summarization")) {
    switch(prev$route %||% "TMT",
      limma_protein = "Limma (Protein)",
      limma_qf      = "Limma + QFeatures",
      prev$route %||% "TMT")
  } else if (identical(type, "ptm_summarization")) {
    if (identical(prev$route %||% "", "dia_ptm")) {
      sw <- prev$meta$software %||% ""
      switch(sw,
        "maxquant"          = "LF PTM \u2014 MaxQuant / Perseus (DDA)",
        "spectronaut_pivot" = "LF PTM \u2014 Spectronaut PTM Pivot (DIA)",
        "LF PTM \u2014 Spectronaut (DIA)")
    } else if (identical(prev$route %||% "", "tmt_ptm")) {
      "TMT PTM \u2014 Proteome Discoverer (limma)"
    } else {
      paste0("PTM \u2014 limma (",
             if (identical(prev$data_format %||% "pd_nested", "flat_peptide"))
               "Flat Peptide" else "PD Nested", ")")
    }
  } else "Unknown"
}

# ---------------------------------------------------------------------------
# .make_groups_display()
# ---------------------------------------------------------------------------
.make_groups_display <- function(prev) {
  type      <- prev$type %||% ""
  group_map <- NULL

  if (type %in% c("lf_summarization", "dia_summarization")) {
    group_map <- prev$cond_groups %||% prev$meta$groups
  } else if (type %in% c("tmt_summarization", "ptm_summarization")) {
    if (!is.null(prev$cond_groups) && length(prev$cond_groups) > 0L) {
      # DIA PTM path stores groups as cond_groups (same structure as lf/dia)
      group_map <- prev$cond_groups
    } else {
      ann   <- prev$annotation_map
      conds <- prev$conditions
      if (!is.null(ann) && !is.null(conds))
        group_map <- lapply(setNames(conds, conds),
                            function(cond) ann$Run[ann$Condition == cond])
    }
  }

  if (is.null(group_map) || length(group_map) == 0)
    return(tags$div(class = "text-muted small", "Groups: \u2014"))

  items <- lapply(names(group_map), function(grp) {
    n <- length(group_map[[grp]])
    tags$div(
      class = "d-flex justify-content-between align-items-center py-0",
      style = "line-height:1.5;",
      tags$span(class = "fw-semibold", style = "color:#2980b9; font-size:0.82em;", grp),
      tags$span(class = "badge bg-light text-dark border",
                style = "font-size:0.7em;", paste0("n=", n))
    )
  })

  tagList(
    tags$div(class = "fw-semibold text-muted mb-1", style = "font-size:0.82em;",
             "Groups:"),
    tags$div(class = "border rounded px-2 py-1", style = "background:#f8f9fa;",
             !!!items)
  )
}
