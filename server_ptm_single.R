# server_ptm_single.R
# ==============================================================================
# Single PTM Inspector — measured-value viewer
#
# This module fits no model.  It displays measured abundances for one site and
# reports the moderated statistics computed by dia_ptm_limma_de() for the
# selected contrast.  Nothing is imputed and no site is hidden.
#
# Pipeline:
#   1. Sync from PTM DE / PTM Trends OR DIA PTM Upload Wizard
#   2. User selects PTM_ID + Reference Condition
#   3. Targeted extraction — NO imputation
#   4. Three-panel visualization:
#      (A) Summarized Abundance (PTM + Protein boxplot/jitter)
#      (B) Stoichiometric Shape (unscaled log2 condition means)
#      (C) Mechanistic Fold Changes (grouped bar chart)
#
# Shared globals used:
#   tmt_register_design_grid(), detect_gene_column(), parse_ptm_mmp_sites()
# ==============================================================================

server_ptm_single <- function(input, output, session, rv, logger) {

  psg <- reactiveValues(
    # Raw parsed files
    parsed_files    = list(),
    files           = list(),
    n_files         = 0L,
    # Design grid state (required by tmt_register_design_grid)
    all_channels    = character(0),
    n_plexes        = 1L,
    has_bridge      = FALSE,
    design_conds    = character(0),
    design_option   = NULL,
    annotation_map  = NULL,
    # PTM type selection
    detected_ptms   = character(0),
    selected_ptms   = "Phospho",
    # Built data tables \u2014 stored in rv$ptm_data (shared store); gene_lookup kept local
    gene_lookup_df  = data.frame(Accession = character(0), Gene = character(0),
                                 stringsAsFactors = FALSE),
    protein_headers = character(0),
    protein_df      = NULL,
    conditions      = NULL,
    cond_order      = NULL,
    # PTM choices and selection
    ptm_choices     = NULL,
    selected_ptm    = NULL,
    # Plot data
    plot_data       = NULL,
    render_count    = 0L,
    status          = "Upload nested PD CSV file(s) or sync from PTM DE / Trends.",
    data_format     = NULL
  )
  # The compute path in this module is gated on %in% c("dia_ptm", "tmt_ptm"); this
  # selector was not widened with it, so a TMT run fell through to the legacy
  # "ptm_single" renderer, which reads prov$ptm_single - a slot the TMT run never
  # writes. The result was no methods text, or, in a session where a legacy
  # label-free PTM analysis had run earlier, that run's text against TMT figures.
  # render_dia_ptm_single_methods() distinguishes TMT from DIA internally.
  methodsTextServer("ptm_sg_methods", provenance = rv$provenance,
                    section = function() {
                      if (isTRUE(psg$design_option %in% c("dia_ptm", "tmt_ptm")))
                        "dia_ptm_single"
                      else "ptm_single"
                    })

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
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
             site      = parts[2],
             modification = paste(parts[3:length(parts)], collapse = "_"))
      } else if (length(parts) == 2) {
        list(accession = parts[1], site = parts[2], modification = "PTM")
      } else {
        list(accession = s, site = "", modification = "")
      }
    }
  }

  .build_ptm_choices <- function(ptm_long, gene_lookup_vec = NULL) {
    ptm_ids <- sort(unique(as.character(ptm_long$ProteinName)))
    labels <- vapply(ptm_ids, function(pid) {
      parsed <- .parse_ptm_id(pid)
      gene <- ""
      gene_vals <- unique(as.character(ptm_long$Gene[ptm_long$ProteinName == pid]))
      gene_vals <- gene_vals[!is.na(gene_vals) & gene_vals != "" & gene_vals != parsed$accession]
      if (length(gene_vals) > 0) gene <- gene_vals[1]
      if ((is.na(gene) || gene == "" || gene == parsed$accession) && !is.null(gene_lookup_vec)) {
        g2 <- gene_lookup_vec[parsed$accession]
        if (!is.na(g2) && g2 != "" && g2 != parsed$accession) gene <- g2
      }
      if (is.na(gene) || gene == "") gene <- parsed$accession
      paste0(gene, "_", parsed$site, " (", parsed$modification, ")")
    }, character(1))
    setNames(ptm_ids, labels)
  }

  .build_gene_lookup_vec <- function() {
    lk <- psg$gene_lookup_df
    if (is.null(lk) || nrow(lk) == 0) return(setNames(character(0), character(0)))
    lk <- lk[!is.na(lk$Accession) & lk$Accession != "", , drop = FALSE]
    lk <- lk[!duplicated(lk$Accession), , drop = FALSE]
    setNames(as.character(lk$Gene), as.character(lk$Accession))
  }

  .extract_col <- function(df, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  # Parent accession from a site key.  The summarizers build site IDs differently
  # per route, so the pattern must match the one that built the key.  Mirrors the
  # route-aware logic in dia_ptm_limma_de(); duplicated here rather than shared to
  # avoid touching global.R.  The Spectronaut branch is '+'-aware, which the
  # previous single-site pattern was not: multi-site keys such as
  # P08779_S44+S48_Phospho_M2 were parsed to "P08779_S44+" and never matched a
  # protein row.
  .psg_parent_acc_of <- function(keys, software = NULL) ptm_parent_acc_of(keys, software)

  # Residue set from a site key, e.g. "S44" or "S44+S48".  Non-greedy anchor so a
  # multi-site key keeps every residue instead of only the last one.
  .psg_site_part_of <- function(key) ptm_site_part_of(key)

  # Plex provenance for the measured-value panels.
  # TMT carries a between-plex variance component the DIA route does not have.
  # IRS or the plex covariate removes it at the model level, but the values shown
  # in these panels are the measured ones, so the offset is still in them and the
  # reader needs to know which plex each point came from.
  # build_pld_from_matrix() returns originalRUN (the channel column name), which
  # is the key into design_df, so the plex is recovered here without touching that
  # shared helper - the DIA route is unaffected.
  # Relabelling happens ONLY when two or more plexes are actually present, so a
  # single-plex TMT run and every DIA run keep their existing labels exactly.
  .psg_add_plex <- function(pld, design_df) {
    if (is.null(pld) || nrow(pld) == 0L) return(pld)
    if (is.null(design_df) || is.null(design_df$mixture) ||
        is.null(pld$originalRUN)) return(pld)
    mx <- as.character(design_df$mixture)[match(pld$originalRUN,
                                                design_df$sample_name)]
    if (all(is.na(mx))) return(pld)
    pld$Plex <- mx
    if (length(unique(mx[!is.na(mx)])) < 2L) return(pld)
    # A Mixture label that already reads as a plex name is used verbatim; a bare
    # value such as "1" or "F1" gets a "plex" prefix so the label is legible.
    pfx <- ifelse(grepl("^(plex|mix)", mx, ignore.case = TRUE), mx,
                  paste0("plex", mx))
    pld$BioReplicate <- ifelse(is.na(mx), as.character(pld$BioReplicate),
                               paste0(pfx, "_", as.character(pld$BioReplicate)))
    pld
  }

  # ── Statistics accessor: read the stored limma result, never recompute ─────
  # Single PTM does not fit its own model.  It reports the moderated statistics
  # produced by dia_ptm_limma_de() for the requested contrast, so that a site's
  # numbers are identical here and in the PTM DE tab.  Returns NULL when no DE
  # result exists for that contrast; the caller renders an empty state and the
  # data panels are unaffected.
  #
  # INTEGRATION POINT (Task 0 / D1): .psg_de_store() must return either a single
  # dia_ptm_limma_de() result list, or a list of them keyed by "<test> vs <ref>".
  .psg_de_store <- function() {
    rv$dia_ptm_de_results
  }

  .psg_de_lookup <- function(site_key, ref_cond, test_cond) {
    store <- tryCatch(.psg_de_store(), error = function(e) NULL)
    if (is.null(store) || length(store) == 0L) return(NULL)
    want <- paste0(test_cond, " vs ", ref_cond)
    res  <- if (!is.null(store$adjusted)) {
      # Single unkeyed result (legacy shape) - accept only if it is this contrast.
      if (identical(store$comparison, want)) store else NULL
    } else if (!is.null(store[[want]])) store[[want]] else NULL
    if (is.null(res)) return(NULL)

    # Join on PTM_ID, not Protein.  server_ptm_de.R's .annotate_ptm_view() moves the
    # matrix row key to PTM_ID and overwrites Protein with the parsed accession;
    # .annotate_prot_view() sets PTM_ID to the accession.  So PTM_ID is the correct
    # key for all three tables and Protein is the correct key for none of them.
    .pull <- function(tbl, key) {
      if (is.null(tbl) || is.null(key) || length(key) != 1L || is.na(key)) return(NULL)
      kcol <- if ("PTM_ID" %in% colnames(tbl)) "PTM_ID"
              else if ("Protein" %in% colnames(tbl)) "Protein"
              else return(NULL)
      i <- which(tbl[[kcol]] == key)
      if (length(i) == 0L) return(NULL)
      tbl[i[1], , drop = FALSE]
    }

    cache_ref <- rv$dia_ptm_master_cache
    abs_key <- if (!is.null(cache_ref$adj_to_abs_map))
                 unname(cache_ref$adj_to_abs_map[site_key]) else NA_character_
    if (length(abs_key) != 1L || is.na(abs_key)) abs_key <- site_key
    parent  <- .psg_parent_acc_of(abs_key, cache_ref$software)
    list(
      adjusted   = .pull(res$adjusted,   site_key),
      unadjusted = .pull(res$unadjusted, abs_key),
      protein    = .pull(res$protein,    parent),
      comparison = want
    )
  }

  # Presence/absence (Fisher) lookup.
  # A site failing the per-comparison valid-value filter cannot be fitted by limma,
  # so dia_ptm_limma_de() assesses it for a non-random detection pattern instead.
  # That result exists but was unreachable from this module, so such a site was
  # reported here as "no DE result stored" - false. It matters most on the TMT
  # route, whose both-arms rule sends every single-arm site down this path.
  #
  # KEYED ON abs_key, NOT on the adjusted row name: .presence_absence() runs on
  # cache$ptm_abs_matrix, so its Protein column holds absolute-matrix row names.
  # The detection columns are named Found_<condition>, so they are constructed from
  # the contrast rather than taken by position.
  .psg_nmar_lookup <- function(abs_key, ref_cond, test_cond) {
    store <- tryCatch(.psg_de_store(), error = function(e) NULL)
    if (is.null(store) || length(store) == 0L) return(NULL)
    want <- paste0(test_cond, " vs ", ref_cond)
    # Accept the same two store shapes .psg_de_lookup accepts: a keyed store, or a
    # single unkeyed result, which is taken only when it is this contrast. Current
    # code always writes a keyed store, so this is defensive - but the two helpers
    # sit side by side and should not disagree about what a store can look like.
    res  <- if (!is.null(store$adjusted)) {
      if (identical(store$comparison, want)) store else NULL
    } else if (!is.null(store[[want]])) store[[want]] else NULL
    if (is.null(res) || is.null(res$nmar) || is.null(res$nmar$table)) return(NULL)
    tb <- res$nmar$table
    if (!("Protein" %in% colnames(tb)) || nrow(tb) == 0L) return(NULL)
    i <- which(as.character(tb$Protein) == as.character(abs_key))
    if (length(i) == 0L) return(NULL)
    row <- tb[i[1], , drop = FALSE]
    .found <- function(cond) {
      cn <- paste0("Found_", cond)
      if (cn %in% colnames(row)) as.character(row[[cn]][1]) else NA_character_
    }
    list(
      comparison       = want,
      detection_diff   = row$Detection_diff[1],
      pvalue           = if ("pvalue"     %in% colnames(row)) row$pvalue[1]     else NA_real_,
      adj_pvalue       = if ("adj.pvalue" %in% colnames(row)) row$adj.pvalue[1] else NA_real_,
      found_test       = .found(test_cond),
      found_ref        = .found(ref_cond),
      has_power        = isTRUE(res$nmar$has_power),
      min_achievable_p = as.numeric(res$nmar$min_achievable_p %||% NA_real_)
    )
  }

  # Assemble one FC table per layer, in the shape the render functions expect.
  .psg_fc_from_de <- function(site_key, ref_cond, test_conds, layer) {
    if (length(test_conds) == 0) return(NULL)
    rows <- lapply(test_conds, function(tc) {
      hit <- .psg_de_lookup(site_key, ref_cond, tc)
      r   <- if (is.null(hit)) NULL else hit[[layer]]
      # Two different reasons for an NA, previously indistinguishable on screen.
      # .psg_de_lookup() returns NULL when the CONTRAST is not in the store - it
      # has not been run in the PTM DE tab. It returns a list whose [[layer]] is
      # NULL when the contrast ran but this SITE is not in that table. A column of
      # NAs on a site detected in every replicate means the former, and reads as
      # the latter unless it is labelled.
      data.frame(
        Label      = paste0(tc, " vs ", ref_cond),
        log2FC     = if (is.null(r)) NA_real_ else r$log2FC[1],
        SE         = if (is.null(r)) NA_real_ else r$SE[1],
        pvalue     = if (is.null(r)) NA_real_ else r$pvalue[1],
        adj.pvalue = if (is.null(r)) NA_real_ else r$adj.pvalue[1],
        Status     = if (is.null(hit)) "not run"
                     else if (is.null(r)) "no estimate"
                     else "ok",
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    if (all(is.na(out$log2FC))) NULL else out
  }

  # Populate PTM selector and condition inputs after data is available
  .update_ui_after_data <- function() {
    gene_lookup_vec <- .build_gene_lookup_vec()

    choices <- .build_ptm_choices(rv$ptm_data$ptm_long, gene_lookup_vec)

    # DIA PTM: always populate choices from the cache, and never restrict them to
    # sites that passed a downstream filter.  This module is an inspector: a site
    # too sparse to be tested or clustered must still be selectable and viewable.
    # Sites are the rows of the cache's primary matrix: the adjusted matrix, or the
    # absolute PTM matrix for a PTM-only cache (whose keys map to themselves).
    .cache    <- rv$dia_ptm_master_cache
    .site_mat <- if (ptm_has_global(.cache)) .cache$ptm_adj_matrix else .cache$ptm_abs_matrix
    if (isTRUE(psg$design_option %in% c("dia_ptm", "tmt_ptm")) && !is.null(.site_mat)) {
      adj_ids     <- rownames(.site_mat)
      site_map    <- .cache$adj_to_abs_map %||% setNames(adj_ids, adj_ids)
      site_gene   <- .cache$site_to_gene
      site_alt <- rv$dia_ptm_master_cache$site_alt_accessions
      labels <- vapply(adj_ids, function(adj_id) {
        abs_id   <- site_map[adj_id]
        gene_lbl <- if (!is.na(abs_id) && abs_id %in% names(site_gene))
          site_gene[[abs_id]]
          # Modification-agnostic fallback. The previous expression was anchored on
          # "_Phospho" and returned the key unchanged both for other modifications
          # and for multi-residue keys such as P08779_S44+S48_Phospho, where the "+"
          # defeats the pattern - displaying a full site ID in place of a gene name.
          else .parse_ptm_id(adj_id)$accession
        # Strip _SharedInference suffix for display
        clean_id <- sub("_SharedInference$", "", adj_id)
        # Alternates appear in the label so the selector's search finds a shared
        # site under any of its proteins, not only the razor winner.
        alt_txt <- ""
        if (!is.null(site_alt) && !is.na(abs_id) && abs_id %in% names(site_alt) &&
            nzchar(site_alt[[abs_id]]))
          alt_txt <- paste0("  [also ", gsub(";", ", ", site_alt[[abs_id]],
                                             fixed = TRUE), "]")
        paste0(gsub(";", " / ", gene_lbl, fixed = TRUE), " | ", clean_id, alt_txt)
      }, character(1))
      choices <- setNames(adj_ids, labels)
      logger(paste0("PTM Single: choices derived from ",
                    if (ptm_has_global(.cache)) "adj" else "absolute (PTM-only)", " matrix"))
    }

    psg$ptm_choices <- choices

    updateSelectizeInput(session, "psg_ptm_select",
                         choices = choices, server = TRUE)

    all_conditions <- psg$conditions
    if (!is.null(all_conditions) && length(all_conditions) > 0) {
      updateSelectInput(session, "psg_reference",
                        choices = all_conditions, selected = all_conditions[1])
    }
  }

  # ===========================================================================
  # SYNC: Inherit from PTM DE or PTM Trends
  # ===========================================================================
  observeEvent(rv$ptm_data$source, {
    src <- rv$ptm_data$source
    if (!isTRUE(src %in% c("ptm_de", "ptm_trends", "dia_ptm_upload", "tmt_ptm_upload"))) return()
    has_tables <- !is.null(rv$ptm_data$ptm_long) && nrow(rv$ptm_data$ptm_long) > 0
    is_dia_ptm_src <- isTRUE(src %in% c("dia_ptm_upload", "tmt_ptm_upload"))
    if (!has_tables && !is_dia_ptm_src) return()

    # Guard against design observer wiping cache during sync
    psg$importing <- TRUE
    session$onFlushed(function() { psg$importing <- FALSE }, once = TRUE)

    # Tables are read directly from rv$ptm_data — no local copies needed
    # Always sync gene_lookup_df (even without ptm_long, e.g. cache-only import)
    if (!is.null(rv$ptm_data$gene_lookup))
      psg$gene_lookup_df  <- rv$ptm_data$gene_lookup
    if (!is.null(rv$ptm_data$data_format))
      psg$data_format    <- rv$ptm_data$data_format
    if (!is.null(rv$ptm_data$annotation_map))
      psg$annotation_map <- rv$ptm_data$annotation_map
    if (!is.null(rv$ptm_data$design_option))
      psg$design_option <- rv$ptm_data$design_option
    if (!is.null(rv$ptm_data$selected_ptms))
      psg$selected_ptms <- rv$ptm_data$selected_ptms
    if (!is.null(rv$ptm_data$protein_headers))
      psg$protein_headers <- rv$ptm_data$protein_headers
    if (!is.null(rv$ptm_data$protein_df))
      psg$protein_df      <- rv$ptm_data$protein_df
    if (!is.null(rv$ptm_data$n_files))
      psg$n_files         <- rv$ptm_data$n_files

    all_conditions <- sort(unique(as.character(rv$ptm_data$conditions)))
    psg$conditions <- all_conditions
    psg$cond_order <- if (!is.null(rv$ptm_data$cond_order) &&
                          setequal(rv$ptm_data$cond_order, all_conditions))
                        rv$ptm_data$cond_order else all_conditions

    # Extract design grid fields
    ann <- rv$ptm_data$annotation_map
    if (!is.null(ann)) {
      psg$all_channels <- unique(ann$ColumnName)
      psg$n_plexes     <- length(unique(ann$Mixture))
      psg$has_bridge   <- any(ann$Condition == "Norm")
      psg$design_conds <- sort(unique(ann$Condition[ann$Condition != "Norm"]))
    }

    if (has_tables) {
      psg$gc_result       <- NULL
      psg$last_reference  <- NULL
      psg$plot_data       <- NULL

      n_ptm <- length(unique(rv$ptm_data$ptm_long$ProteinName))
      psg$status <- paste0("Synced from ", src, ": ", n_ptm,
                           " PTM sites, ", length(all_conditions), " conditions. ",
                           "Run Summarization to proceed.")
      logger(paste0("PTM Single: Synced from ", src, " \u2014 ", n_ptm, " PTM sites"))
    } else if (is_dia_ptm_src) {
      psg$gc_result       <- NULL
      psg$last_reference  <- NULL
      psg$plot_data       <- NULL
      cache <- rv$dia_ptm_master_cache
      n_sites <- if (is.null(cache)) "?"
                 else if (ptm_has_global(cache)) cache$n_sites_adjusted
                 else cache$n_sites_output
      .sw_single_lbl <- if (identical(src, "tmt_ptm_upload")) {
        "TMT PTM (limma)"
      } else {
        switch(rv$dia_ptm_software %||% "spectronaut",
          "maxquant"          = "Label-Free DDA (MaxQuant)",
          "spectronaut_pivot" = "DIA (Spectronaut Pivot)",
          "DIA (Spectronaut)")
      }
      psg$status <- paste0(.sw_single_lbl, " data synced: ", n_sites,
                           " phosphosite rows across ", length(all_conditions),
                           " conditions. Select a site and click Plot.")
      logger(paste0("PTM Single: DIA PTM synced \u2014 ", n_sites, " adj rows"))
    } else {
      psg$gc_result       <- NULL
      psg$last_reference  <- NULL
      psg$plot_data       <- NULL
      psg$status <- paste0("Design synced from ", src,
                           " (summarization cache imported). ",
                           length(all_conditions), " conditions.")
      logger(paste0("PTM Single: Design synced from ", src, " (import mode)"))
    }

    .update_ui_after_data()
  }, ignoreInit = TRUE)

  # Condition order is a property of the experiment, not of a tab.  Mirror changes
  # made in PTM Trends without clobbering a local reorder of the same set.
  observeEvent(rv$ptm_data$cond_order, {
    new_order <- rv$ptm_data$cond_order
    if (is.null(new_order) || length(new_order) == 0) return()
    if (!is.null(psg$conditions) && setequal(new_order, psg$conditions) &&
        !identical(new_order, psg$cond_order))
      psg$cond_order <- new_order
  }, ignoreInit = TRUE)

  # ── Stale cache warning on design changes ──
  observeEvent(psg$design_option, {
    if (isTRUE(psg$importing)) return()

  }, ignoreInit = TRUE)

  # ===========================================================================
  # CONDITION ORDER MANAGEMENT
  # ===========================================================================
  output$psg_cond_label <- renderUI({
    conds <- psg$cond_order
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

  observeEvent(input$psg_reorder_btn, {
    req(psg$conditions)
    conds <- psg$cond_order
    if (is.null(conds)) conds <- psg$conditions

    showModal(modalDialog(
      title = tags$span(icon("arrows-up-down"), " Reorder Conditions"),
      tagList(
        helpText(class = "small text-muted",
          "Drag conditions to set the X-axis order (left to right). ",
          "Typically chronological or dose-increasing."),
        sortable::rank_list(
          text       = "Condition Order (drag to reorder):",
          labels     = conds,
          input_id   = "psg_cond_order_sorted",
          options    = sortable::sortable_options(multiDrag = FALSE),
          class      = "default-sortable"
        )
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("psg_cond_order_save", "Apply Order",
                     class = "btn-success", icon = icon("check"))
      )
    ))
  })

  observeEvent(input$psg_cond_order_save, {
    new_order <- input$psg_cond_order_sorted
    if (is.null(new_order) || length(new_order) == 0) {
      showNotification("No conditions to order.", type = "warning")
      return()
    }
    psg$cond_order         <- new_order
    rv$ptm_data$cond_order <- new_order
    removeModal()
    logger(paste0("PTM Single: Condition order set to: ",
                  paste(new_order, collapse = " -> ")))
  })

  # ===========================================================================
  # EXPORT SUMMARIZATION (.rds)


  # ===========================================================================
  # Check if summarization is cached (for disabling plot button)
  # ===========================================================================
  output$psg_cache_ready <- reactive({
    !is.null(rv$dia_ptm_master_cache)
  })
  outputOptions(output, "psg_cache_ready", suspendWhenHidden = FALSE)

  # ===========================================================================
  # PLOT PTM (targeted extraction + DE lookup + visualization)
  # ===========================================================================
  observeEvent(input$psg_plot_btn, {
    # ── DIA PTM BRANCH (also handles tmt_ptm which uses the same cache path) ──
    is_dia_ptm <- isTRUE(psg$design_option %in% c("dia_ptm", "tmt_ptm"))
    if (is_dia_ptm) {
      cache <- rv$dia_ptm_master_cache
      if (is.null(cache)) {
        showNotification("PTM summarization cache not found. Complete the PTM Upload Wizard.", type = "warning")
        return()
      }
      ptm_id   <- input$psg_ptm_select
      ref_cond <- input$psg_reference
      if (is.null(ptm_id) || ptm_id == "") {
        showNotification("Select a PTM site first.", type = "warning")
        return()
      }

      # Build ProteinLevelData from the absolute (unadjusted) matrix for the selected site
      # Use adj_to_abs_map to handle _SharedInference suffixed IDs
      abs_id   <- if (!is.null(cache$adj_to_abs_map))
                    unname(cache$adj_to_abs_map[ptm_id]) else NA_character_
      if (length(abs_id) != 1L || is.na(abs_id)) abs_id <- ptm_id   # fallback: treat as direct

      # Build PLDs using build_pld_from_matrix
      ptm_pld  <- build_pld_from_matrix(
        cache$ptm_abs_matrix[abs_id, , drop = FALSE],
        cache$design_df
      )
      # Parse metadata from the primary (first) component of abs_id
      primary_abs <- strsplit(abs_id, ";", fixed = TRUE)[[1]][1]
      parent_acc  <- .psg_parent_acc_of(primary_abs, cache$software)
      # Residue set, e.g. "S100" or "S100+S104".  ptm_site_part_of() is anchored on
      # "_Phospho" and returns "" for any other modification, so the fallback must
      # not be anchored the same way.  .parse_ptm_id() is modification-agnostic and
      # also handles multi-residue keys; on a Phospho key the two agree exactly.
      .parsed_abs  <- .parse_ptm_id(primary_abs)
      site_id_part <- .psg_site_part_of(primary_abs)
      if (!nzchar(site_id_part)) site_id_part <- .parsed_abs$site
      .mod_lbl <- .parsed_abs$modification
      # "_SharedInference" is a row-expansion suffix on the DIA route, not part of
      # the modification. .parse_ptm_id() cannot separate it (the tail contains an
      # underscore, so its regex falls through to the split branch and glues the
      # suffix onto the modification), and it reaches here whenever the
      # adj_to_abs_map lookup misses and abs_id falls back to the adjusted id.
      .mod_lbl <- sub("_SharedInference$", "", .mod_lbl)
      if (is.null(.mod_lbl) || is.na(.mod_lbl) || !nzchar(.mod_lbl))
        .mod_lbl <- "Phospho"

      prot_pld <- if (!is.null(cache$global_matrix) &&
                       parent_acc %in% rownames(cache$global_matrix)) {
        build_pld_from_matrix(
          cache$global_matrix[parent_acc, , drop = FALSE],
          cache$design_df
        )
      } else NULL
      # Tag both layers with their plex before anything reads them.
      ptm_pld  <- .psg_add_plex(ptm_pld,  cache$design_df)
      prot_pld <- .psg_add_plex(prot_pld, cache$design_df)
      # A PTM-only cache has no protein layer by design: nothing to warn about
      if (is.null(prot_pld) && ptm_has_global(cache)) {
        showNotification(paste0("Parent protein '", parent_acc,
                                "' was not quantified in the global matrix - protein and ",
                                "adjusted layers are unavailable for this site."),
                         type = "warning", duration = 8)
        logger(paste0("PTM Single [DIA PTM]: parent protein not found for ", primary_abs,
                      " (parsed as '", parent_acc, "')"))
      }

      # Gene name
      gene_lbl  <- ""
      abs_gene  <- cache$site_to_gene[abs_id]
      if (!is.na(abs_gene) && nzchar(abs_gene)) {
        gene_lbl <- strsplit(abs_gene, ";", fixed = TRUE)[[1]][1]
      } else if (!is.null(cache$acc_to_gene) && parent_acc %in% names(cache$acc_to_gene)) {
        gene_lbl <- cache$acc_to_gene[[parent_acc]]
      }
      if (is.na(gene_lbl) || !nzchar(gene_lbl)) gene_lbl <- parent_acc
      ptm_label <- paste0(gene_lbl, "_", site_id_part, "_", .mod_lbl)

      # Condition order
      cond_order_dia <- if (!is.null(psg$cond_order) && length(psg$cond_order) > 0)
        psg$cond_order else psg$conditions
      ref_cond_final  <- if (!is.null(ref_cond) && nzchar(ref_cond))
        ref_cond else cond_order_dia[1]
      non_ref_conds_dia <- setdiff(cond_order_dia, ref_cond_final)

      # Missing conditions
      ptm_conds_present <- unique(as.character(
        ptm_pld$Condition[!is.na(ptm_pld$Abundance)]))
      missing_conds_dia <- setdiff(cond_order_dia, ptm_conds_present)

      # ── Line profile (Panel B) ────────────────────────────────────────────
      # Each line is aggregated independently from its own source, so a line is
      # never silently restricted to the samples where a DIFFERENT measurement
      # also exists.  The adjusted line is read from cache$ptm_adj_matrix - the
      # canonical, twice-centred matrix that PTM Trends and PTM DE both use - and
      # is NOT re-derived as (PTM - protein) here; re-deriving it drops the
      # post-subtraction median centering and yields a different fold change.
      adj_pld <- if (!is.null(cache$ptm_adj_matrix) &&
                     ptm_id %in% rownames(cache$ptm_adj_matrix)) {
        build_pld_from_matrix(
          cache$ptm_adj_matrix[ptm_id, , drop = FALSE],
          cache$design_df
        )
      } else NULL

      .cond_mean_of <- function(pld) {
        if (is.null(pld) || nrow(pld) == 0) return(NULL)
        v <- pld[!is.na(pld$Abundance), , drop = FALSE]
        if (nrow(v) == 0) return(NULL)
        a <- aggregate(Abundance ~ Condition, data = v, FUN = mean)
        setNames(a$Abundance, as.character(a$Condition))
      }
      .cond_n_of <- function(pld) {
        if (is.null(pld) || nrow(pld) == 0) return(setNames(integer(0), character(0)))
        v <- pld[!is.na(pld$Abundance), , drop = FALSE]
        if (nrow(v) == 0) return(setNames(integer(0), character(0)))
        tb <- table(as.character(v$Condition))
        setNames(as.integer(tb), names(tb))
      }

      .m_ptm  <- .cond_mean_of(ptm_pld)
      .m_prot <- .cond_mean_of(prot_pld)
      .m_adj  <- .cond_mean_of(adj_pld)
      .n_ptm  <- .cond_n_of(ptm_pld)
      .n_prot <- .cond_n_of(prot_pld)
      .n_adj  <- .cond_n_of(adj_pld)

      stoich_df_dia <- NULL
      if (!is.null(.m_ptm)) {
        stoich_df_dia <- data.frame(
          Condition      = cond_order_dia,
          Mean_Adjusted  = if (is.null(.m_adj))  NA_real_ else unname(.m_adj[cond_order_dia]),
          Mean_Protein   = if (is.null(.m_prot)) NA_real_ else unname(.m_prot[cond_order_dia]),
          Mean_Total_PTM = unname(.m_ptm[cond_order_dia]),
          N_PTM          = unname(ifelse(is.na(.n_ptm[cond_order_dia]),  0L, .n_ptm[cond_order_dia])),
          N_Protein      = unname(ifelse(is.na(.n_prot[cond_order_dia]), 0L, .n_prot[cond_order_dia])),
          N_Adjusted     = unname(ifelse(is.na(.n_adj[cond_order_dia]),  0L, .n_adj[cond_order_dia])),
          stringsAsFactors = FALSE
        )
        stoich_df_dia$Condition <- factor(stoich_df_dia$Condition, levels = cond_order_dia)
        stoich_df_dia <- stoich_df_dia[order(stoich_df_dia$Condition), , drop = FALSE]
        stoich_df_dia <- stoich_df_dia[!is.na(stoich_df_dia$Mean_Total_PTM), , drop = FALSE]
      }

      # ── Fit any contrast against this reference that has not been fitted yet ──
      # PTM Single reports the DE tab's moderated statistics and never fits its own
      # model, so a contrast the DE tab has not run has no numbers here. That is
      # what produced a site detected in 4/4 replicates of every condition showing
      # NA fold changes for every contrast except the one that happened to be
      # loaded. Every missing contrast against this reference is fitted now,
      # through the same helper the DE tab uses and published under the same key,
      # so a site's statistics are identical whichever tab produced them.
      # Contrasts already in the store are never refitted, so a second Plot on the
      # same reference costs nothing.
      .bf_store <- rv$dia_ptm_de_results %||% list()
      .bf_todo  <- non_ref_conds_dia[
        !(paste0(non_ref_conds_dia, " vs ", ref_cond_final) %in% names(.bf_store))]
      if (length(.bf_todo) > 0) {
        .bf_n <- length(.bf_todo)
        .bf <- NULL
        withProgress(
          message = paste0("Fitting ", .bf_n, " contrast(s) vs ", ref_cond_final, "..."),
          value = 0, {
            .bf <- ptm_de_backfill_contrasts(
              cache           = cache,
              ref_cond        = ref_cond_final,
              test_conds      = non_ref_conds_dia,
              store           = .bf_store,
              design_option   = psg$design_option,
              subject_pairing = rv$dia_ptm_subject_pairing,
              logger          = logger,
              progress        = function(lbl)
                incProgress(1 / .bf_n, detail = lbl)
            )
          })
        if (!is.null(.bf)) {
          # The k-NN memo the DE tab would also have written. Cache-scoped and
          # seeded, so the filled values do not depend on which tab computed them.
          if (!is.null(.bf$cache_update)) {
            rv$dia_ptm_master_cache$ptm_abs_matrix_imputed <- .bf$cache_update$ptm_abs_matrix_imputed
            rv$dia_ptm_master_cache$global_matrix_imputed  <- .bf$cache_update$global_matrix_imputed
            rv$dia_ptm_master_cache$impute_stats           <- .bf$cache_update$impute_stats
            rv$dia_ptm_master_cache$ptm_impute_cap         <- .bf$cache_update$ptm_impute_cap
          }
          rv$dia_ptm_de_results <- .bf$store
          if (length(.bf$ran) > 0)
            logger(paste0("PTM Single [DIA PTM]: fitted and published ",
                          length(.bf$ran), " contrast(s) on demand - ",
                          paste(.bf$ran, collapse = "; ")))
          if (length(.bf$failed) > 0) {
            showNotification(
              paste0("Some contrasts could not be fitted: ",
                     paste(.bf$failed, collapse = "; "),
                     ". Their fold changes will be blank."),
              type = "warning", duration = 10)
            logger(paste0("PTM Single [DIA PTM]: ", length(.bf$failed),
                          " contrast(s) failed - ", paste(.bf$failed, collapse = "; ")))
          }
        }
      }

      # ── Fold changes: read from the stored limma fit, never recomputed here ──
      # Single PTM reports the same moderated statistics the PTM DE tab reports.
      # When no DE result exists for a contrast every field is NA and the panels
      # render an explicit empty state; the measured-data panels are unaffected.
      fc_unadjusted_dia <- .psg_fc_from_de(ptm_id, ref_cond_final,
                                           non_ref_conds_dia, "unadjusted")
      fc_protein_dia    <- .psg_fc_from_de(ptm_id, ref_cond_final,
                                           non_ref_conds_dia, "protein")
      fc_adjusted_dia   <- .psg_fc_from_de(ptm_id, ref_cond_final,
                                           non_ref_conds_dia, "adjusted")
      de_available_dia  <- !is.null(fc_adjusted_dia) || !is.null(fc_unadjusted_dia)
      # No limma result does not mean the site was not analysed. Check the
      # presence/absence set before reporting an absence of evidence.
      nmar_dia <- list()
      if (!de_available_dia && length(non_ref_conds_dia) > 0L) {
        nmar_dia <- lapply(non_ref_conds_dia, function(tc)
          .psg_nmar_lookup(abs_id, ref_cond_final, tc))
        names(nmar_dia) <- non_ref_conds_dia
        nmar_dia <- nmar_dia[!vapply(nmar_dia, is.null, logical(1))]
      }
      nmar_available_dia <- length(nmar_dia) > 0L
      if (!de_available_dia)
        logger(paste0("PTM Single [DIA PTM]: no stored limma result for reference '",
                      ref_cond_final, "' - ",
                      if (nmar_available_dia)
                        paste0("site is in the presence/absence set for ",
                               length(nmar_dia), " contrast(s); reporting detection counts.")
                      else "statistics panel will show an empty state."))

      psg$selected_ptm <- ptm_id
      psg$gc_result    <- NULL
      psg$plot_data    <- list(
        ptm_label     = ptm_label,
        ptm_id        = ptm_id,
        abs_id        = abs_id,
        accession     = parent_acc,
        site          = site_id_part,
        modification  = .mod_lbl,
        gene          = gene_lbl,
        ref_cond      = ref_cond_final,
        cond_order    = cond_order_dia,
        ptm_pld       = ptm_pld,
        prot_pld      = prot_pld,
        stoich_df     = stoich_df_dia,
        fc_adjusted   = fc_adjusted_dia,
        fc_unadjusted = fc_unadjusted_dia,
        fc_protein    = fc_protein_dia,
        non_ref_conds = non_ref_conds_dia,
        has_protein   = !is.null(prot_pld) && nrow(prot_pld) > 0,
        # FALSE for a PTM-only cache: no adjusted or protein layer exists at all
        has_global    = ptm_has_global(cache),
        missing_conds = missing_conds_dia,
        de_available  = de_available_dia,
        nmar          = nmar_dia,
        nmar_available = nmar_available_dia,
        source        = "dia_ptm"
      )
      logger(paste0("PTM Single [DIA PTM]: Plot data built for ", ptm_label,
                    " (", ptm_id, "), ref=", ref_cond_final))
      rv$provenance$dia_ptm_single_run <- list(
        ptm_id        = ptm_id,
        abs_id        = abs_id,
        gene          = gene_lbl,
        site          = site_id_part,
        modification  = .mod_lbl,
        reference_condition = ref_cond_final,
        n_contrasts   = length(non_ref_conds_dia),
        has_protein   = !is.null(prot_pld) && nrow(prot_pld) > 0,
        de_available  = de_available_dia,
        nmar_available = nmar_available_dia,
        stats_source  = "dia_ptm_limma_de",
        imputation    = "none",
        missing_conds = missing_conds_dia
      )
      psg$status <- paste0("Plotted: ", ptm_label, " (", ptm_id, ")")
      return()
    }
    logger(paste0("PTM Single: Plot clicked but design_option is '",
                  as.character(psg$design_option %||% "NULL"),
                  "' - this handler serves the DIA/TMT PTM cache route only; no action taken."))
  })

  # ===========================================================================
  # PANEL A: Summarized Abundance (PTM + Protein boxplot/jitter)
  # ===========================================================================
  output$psg_abundance_plot <- renderPlotly({
    req(psg$plot_data)
    pd <- psg$plot_data
    cond_order <- pd$cond_order

    # Build combined data frame for PTM and Protein
    plot_rows <- list()

    if (!is.null(pd$ptm_pld) && nrow(pd$ptm_pld) > 0) {
      ptm_d <- pd$ptm_pld[!is.na(pd$ptm_pld$Abundance), , drop = FALSE]
      if (nrow(ptm_d) > 0) {
        plot_rows[["PTM"]] <- data.frame(
          Condition = as.character(ptm_d$Condition),
          BioReplicate = as.character(ptm_d$BioReplicate),
          Plex = if (!is.null(ptm_d$Plex)) as.character(ptm_d$Plex) else NA_character_,
          Abundance = ptm_d$Abundance,
          Layer = "PTM Abundance",
          stringsAsFactors = FALSE
        )
      }
    }

    if (!is.null(pd$prot_pld) && nrow(pd$prot_pld) > 0) {
      prot_d <- pd$prot_pld[!is.na(pd$prot_pld$Abundance), , drop = FALSE]
      if (nrow(prot_d) > 0) {
        plot_rows[["Protein"]] <- data.frame(
          Condition = as.character(prot_d$Condition),
          BioReplicate = as.character(prot_d$BioReplicate),
          Plex = if (!is.null(prot_d$Plex)) as.character(prot_d$Plex) else NA_character_,
          Abundance = prot_d$Abundance,
          Layer = "Global Protein",
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(plot_rows) == 0) {
      return(plotly_empty() %>% layout(title = "No abundance data available"))
    }

    combined <- do.call(rbind, plot_rows)
    combined$Condition <- factor(combined$Condition, levels = cond_order)
    combined$Layer <- factor(combined$Layer, levels = c("PTM Abundance", "Global Protein"))

    ptm_color  <- "#E74C3C"
    prot_color <- "#3498DB"

    p <- plot_ly()

    for (layer in levels(combined$Layer)) {
      layer_data <- combined[combined$Layer == layer, , drop = FALSE]
      if (nrow(layer_data) == 0) next

      col <- if (layer == "PTM Abundance") ptm_color else prot_color

      # Boxplot
      p <- add_trace(p, data = layer_data,
                     x = ~Condition, y = ~Abundance,
                     type = "box",
                     name = layer,
                     marker = list(color = col),
                     line = list(color = col),
                     fillcolor = paste0(col, "33"),
                     boxpoints = FALSE,
                     legendgroup = layer)

      # Jittered points
      # Build detection fraction text from missingness map
      det_label <- character(nrow(layer_data))
      miss_map <- rv$dia_ptm_master_cache$missingness_map
      site_key_psg <- pd$abs_id
      if (!is.null(miss_map) && nrow(miss_map) > 0 &&
          all(c("Site_ID", "Condition", "Detected", "Total") %in% colnames(miss_map)) &&
          layer == "PTM Abundance") {
        for (ri in seq_len(nrow(layer_data))) {
          cond <- as.character(layer_data$Condition[ri])
          mm_row <- miss_map[miss_map$Site_ID == site_key_psg &
                             miss_map$Condition == cond, , drop = FALSE]
          if (nrow(mm_row) > 0) {
            det_label[ri] <- paste0("<br>Detected: ", mm_row$Detected[1],
                                    "/", mm_row$Total[1])
          }
        }
      }

      .plex_txt <- if (!is.null(layer_data$Plex) && !all(is.na(layer_data$Plex)))
        paste0("<br>Plex: ", layer_data$Plex) else ""
      hover <- paste0(
        "<b>", layer, "</b>",
        "<br>Condition: ", layer_data$Condition,
        "<br>BioReplicate: ", layer_data$BioReplicate,
        .plex_txt,
        "<br>log2(Abundance): ", round(layer_data$Abundance, 3),
        det_label
      )

      # Symbol by plex where more than one is present, so between-plex spread is
      # visible without reading the hover. Falls back to a single symbol otherwise,
      # which is what every DIA run and every single-plex TMT run gets.
      .plex_lv  <- if (!is.null(layer_data$Plex)) unique(layer_data$Plex[!is.na(layer_data$Plex)]) else character(0)
      .plex_sym <- if (length(.plex_lv) >= 2L) {
        .syms <- c("circle", "square", "diamond", "cross", "x",
                   "triangle-up", "triangle-down", "pentagon")
        unname(setNames(rep_len(.syms, length(.plex_lv)),
                        .plex_lv)[as.character(layer_data$Plex)])
      } else "circle"
      p <- add_trace(p, data = layer_data,
                     x = ~Condition, y = ~Abundance,
                     type = "scatter", mode = "markers",
                     marker = list(size = 8, color = col, symbol = .plex_sym,
                                   line = list(width = 1, color = "#333"),
                                   opacity = 0.7),
                     name = paste0(layer, " (points)"),
                     legendgroup = layer,
                     showlegend = FALSE,
                     hoverinfo = "text", text = hover,
                     jitter = 0.3)
    }

    p <- layout(p,
      title = paste0("Summarized Abundance: ", pd$ptm_label,
                      " (", pd$modification, ")"),
      xaxis = list(title = "Condition", categoryorder = "array",
                   categoryarray = cond_order),
      yaxis = list(title = "log2(Summarized Intensity)"),
      boxmode = "group",
      legend = list(orientation = "h", x = 0, y = -0.15)
    )
    p
  })

  # ===========================================================================
  # PANEL C: Mechanistic Fold Changes — line profiles across contrasts
  # (three traces: Adjusted (PTM - Protein), Unadjusted PTM, Global Protein)
  # ===========================================================================
  output$psg_fc_plot <- renderPlotly({
    req(psg$plot_data)
    pd <- psg$plot_data

    non_ref_conds <- pd$non_ref_conds
    if (length(non_ref_conds) == 0) {
      return(plotly_empty() %>% layout(title = "No contrasts to display"))
    }
    if (!isTRUE(pd$de_available)) {
      .avail <- names(rv$dia_ptm_de_results)
      # Three distinct states, previously collapsed into one misleading message:
      #   (a) the site was excluded from the model and assessed by the
      #       presence/absence route - a result exists, it is just not a fold change
      #   (b) no PTM DE run exists for this reference condition at all
      if (isTRUE(pd$nmar_available)) {
        .nm  <- pd$nmar
        # Power is a property of the CONTRAST, not of the site: min_achievable_p is
        # fisher.test(matrix(c(n_test, 0, 0, n_ref))), so it depends only on the two
        # group sizes. In an unbalanced design (say 4/4/2) one contrast against the
        # same reference can have power while another does not, so the verdict is
        # stated on each contrast's own line rather than once for all of them.
        .txt <- vapply(names(.nm), function(tc) {
          x <- .nm[[tc]]
          paste0(tc, " vs ", pd$ref_cond, ": detected ", x$found_test, " vs ",
                 x$found_ref,
                 if (isTRUE(x$has_power) && !is.na(x$adj_pvalue))
                   paste0(", Fisher adj.p = ", signif(x$adj_pvalue, 3))
                 else if (!isTRUE(x$has_power))
                   paste0(" (no attainable significance at these replicate counts; ",
                          "smallest possible two-sided p = ",
                          signif(x$min_achievable_p, 3), ", so p is reported as NA)")
                 else " (p not available for this site)")
        }, character(1))
        return(plotly_empty() %>% layout(
          title = paste0("This site was not fitted by the linear model.<br>",
                         "<span style='font-size:12px'>It did not meet the ",
                         "valid-value requirement in both compared conditions, so it ",
                         "was assessed for a non-random detection pattern instead ",
                         "(Fisher's exact test on detected / not-detected counts, ",
                         "corrected within that set only).<br>",
                         paste(.txt, collapse = "<br>"),
                         "</span>")))
      }
      # Two different situations, and after the on-demand fitting in PTM Single
      # only one of them is "nothing has been run". If contrasts against this
      # reference are in the store, the site is simply absent from all of them -
      # saying "no result stored" would be false. Decide from the store, not from
      # the fact that no fold change came back.
      .avail_ref <- .avail[grepl(paste0(" vs ", pd$ref_cond, "$"), .avail)]
      return(plotly_empty() %>% layout(
        title = paste0("No fold change is available for this site against '",
                       pd$ref_cond, "'.<br>",
                       "<span style='font-size:12px'>Measured abundances are shown in the ",
                       "panels above. ",
                       if (length(.avail_ref) > 0)
                         paste0("The comparison(s) against this reference have been ",
                                "run, so this site was not returned by any of them - ",
                                "it did not survive summarization, or it carries no ",
                                "estimable coefficient in this contrast set.")
                       else
                         paste0("No comparison against this reference has been run ",
                                "yet. Run PTM DE for it, or re-open this site, to ",
                                "populate fold changes and moderated statistics."),
                       if (length(.avail) > 0)
                         paste0("<br>Available: ", paste(.avail, collapse = "; "))
                       else "",
                       "</span>")))
    }

    contrast_labels <- c(paste0(pd$ref_cond, " (ref)"), paste0(non_ref_conds, " vs ", pd$ref_cond))

    # Match FC data by label ordering
    .match_fc <- function(fc_df, contrast_labels) {
      out <- rep(NA_real_, length(contrast_labels))
      out_se <- rep(NA_real_, length(contrast_labels))
      out_p  <- rep(NA_real_, length(contrast_labels))
      if (!is.null(fc_df) && nrow(fc_df) > 0) {
        for (i in seq_along(contrast_labels)) {
          idx <- which(fc_df$Label == contrast_labels[i])
          if (length(idx) > 0) {
            out[i]    <- fc_df$log2FC[idx[1]]
            out_se[i] <- fc_df$SE[idx[1]]
            out_p[i]  <- fc_df$adj.pvalue[idx[1]]
          }
        }
      }
      list(fc = out, se = out_se, p = out_p)
    }

    fc_labels <- paste0(non_ref_conds, " vs ", pd$ref_cond)
    adj_vals   <- .match_fc(pd$fc_adjusted, fc_labels)
    unadj_vals <- .match_fc(pd$fc_unadjusted, fc_labels)
    prot_vals  <- .match_fc(pd$fc_protein, fc_labels)

    # Prepend reference condition at (0, 0)
    adj_vals$fc   <- c(0, adj_vals$fc);   adj_vals$se   <- c(0, adj_vals$se);   adj_vals$p   <- c(NA, adj_vals$p)
    unadj_vals$fc <- c(0, unadj_vals$fc); unadj_vals$se <- c(0, unadj_vals$se); unadj_vals$p <- c(NA, unadj_vals$p)
    prot_vals$fc  <- c(0, prot_vals$fc);  prot_vals$se  <- c(0, prot_vals$se);  prot_vals$p  <- c(NA, prot_vals$p)

    p <- plot_ly()

    # Build per-condition detection labels from missingness map
    miss_map <- rv$dia_ptm_master_cache$missingness_map
    site_key_psg <- pd$abs_id
    det_by_cond <- character(length(contrast_labels))
    if (!is.null(miss_map) && nrow(miss_map) > 0 &&
        all(c("Site_ID", "Condition", "Detected", "Total") %in% colnames(miss_map))) {
      all_test_conds <- c(pd$ref_cond, non_ref_conds)
      for (ci in seq_along(all_test_conds)) {
        cond <- all_test_conds[ci]
        mm_row <- miss_map[miss_map$Site_ID == site_key_psg &
                           miss_map$Condition == cond, , drop = FALSE]
        if (nrow(mm_row) > 0) {
          det_by_cond[ci] <- paste0("<br>Detected: ", mm_row$Detected[1],
                                    "/", mm_row$Total[1])
        }
      }
    }

    # PTM-only cache: the absolute PTM line alone (no adjusted or protein layer).
    # plot_data saved before has_global existed always had a global proteome.
    .with_global <- !isFALSE(pd$has_global)

    # Adjusted (red)
    if (.with_global) {
      hover_adj <- paste0(
        "<b>Adjusted (Stoichiometric)</b>",
        "<br>Contrast: ", contrast_labels,
        "<br>log2FC: ", round(adj_vals$fc, 3),
        "<br>SE: \u00b1", round(adj_vals$se, 3),
        "<br>adj.p: ", signif(adj_vals$p, 3),
        det_by_cond
      )
      p <- add_trace(p, x = contrast_labels, y = adj_vals$fc,
                     type = "scatter", mode = "lines+markers",
                     line = list(color = "#E74C3C", width = 3),
                     marker = list(size = 9, color = "#E74C3C",
                                   line = list(width = 1, color = "#333")),
                     hoverinfo = "text", text = hover_adj,
                     name = "Adjusted (PTM \u2212 Protein)")
    }

    # Unadjusted PTM (blue)
    hover_unadj <- paste0(
      "<b>Unadjusted PTM</b>",
      "<br>Contrast: ", contrast_labels,
      "<br>log2FC: ", round(unadj_vals$fc, 3),
      "<br>SE: \u00b1", round(unadj_vals$se, 3),
      "<br>adj.p: ", signif(unadj_vals$p, 3),
      det_by_cond
    )
    p <- add_trace(p, x = contrast_labels, y = unadj_vals$fc,
                   type = "scatter", mode = "lines+markers",
                   line = list(color = "#3498DB", width = 3),
                   marker = list(size = 9, color = "#3498DB",
                                 line = list(width = 1, color = "#333")),
                   hoverinfo = "text", text = hover_unadj,
                   name = "Unadjusted PTM")

    # Global Protein (grey)
    if (.with_global) {
      hover_prot <- paste0(
        "<b>Global Protein</b>",
        "<br>Contrast: ", contrast_labels,
        "<br>log2FC: ", round(prot_vals$fc, 3),
        "<br>SE: \u00b1", round(prot_vals$se, 3),
        "<br>adj.p: ", signif(prot_vals$p, 3)
      )
      p <- add_trace(p, x = contrast_labels, y = prot_vals$fc,
                     type = "scatter", mode = "lines+markers",
                     line = list(color = "#95A5A6", width = 3),
                     marker = list(size = 9, color = "#95A5A6",
                                   line = list(width = 1, color = "#333")),
                     hoverinfo = "text", text = hover_prot,
                     name = "Global Protein")
    }

    p <- layout(p,
      title = paste0("Mechanistic Fold Changes: ", pd$ptm_label,
                      " (vs ", pd$ref_cond, ")"),
      xaxis = list(title = "Contrast", categoryorder = "array",
                   categoryarray = contrast_labels),
      yaxis = list(title = "log2 Fold Change"),
      legend = list(orientation = "h", x = 0, y = -0.15),
      shapes = list(
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = 0, y1 = 0,
             line = list(dash = "dash", color = "grey50", width = 1))
      )
    )
    p
  })

  # ===========================================================================
  # STATUS OUTPUT
  # ===========================================================================
  output$psg_status <- renderText({
    psg$status
  })

  # -- NA condition warning banner --
  output$psg_na_warning <- renderUI({
    req(psg$plot_data)
    mc <- psg$plot_data$missing_conds
    if (length(mc) == 0) return(NULL)
    tags$div(
      class = "alert alert-warning py-2 px-3 small mb-2",
      icon("triangle-exclamation"),
      tags$b(paste0(length(mc), " condition(s) have no PTM data: ")),
      paste(mc, collapse = ", "),
      " \u2014 these conditions are absent from the plots below."
    )
  })

  output$psg_plot_status <- renderText({
    req(psg$plot_data)
    pd <- psg$plot_data

    .safe_fc <- function(tbl, i) {
      if (is.null(tbl) || nrow(tbl) < i) return(NA_real_)
      tbl$log2FC[i]
    }
    # Lines follow the primary table: adjusted, or absolute PTM for a PTM-only cache
    .with_global <- !isFALSE(pd$has_global)
    .fc_main     <- if (.with_global) pd$fc_adjusted else pd$fc_unadjusted
    fc_lines <- character(0)
    if (!is.null(.fc_main) && nrow(.fc_main) > 0) {
      for (i in seq_len(nrow(.fc_main))) {
        .st <- if ("Status" %in% colnames(.fc_main))
                 as.character(.fc_main$Status[i]) else "ok"
        fc_lines <- c(fc_lines, if (identical(.st, "not run")) paste0(
          "  ", .fc_main$Label[i],
          ": comparison not yet run in PTM DE - run it there to populate ",
          "fold changes and p-values for this contrast."
        ) else paste0(
          "  ", .fc_main$Label[i],
          if (.with_global)
            paste0(": adj.FC=", round(.fc_main$log2FC[i], 3),
                   ", unadj.FC=", round(.safe_fc(pd$fc_unadjusted, i), 3),
                   ", prot.FC=", round(.safe_fc(pd$fc_protein, i), 3))
          else
            paste0(": FC=", round(.fc_main$log2FC[i], 3), " (absolute PTM)"),
          ", adj.p=", signif(.fc_main$adj.pvalue[i], 3),
          if (identical(.st, "no estimate"))
            " (site not present in this contrast's result table)" else ""
        ))
      }
    }

    det_lines <- character(0)
    mm <- rv$dia_ptm_master_cache$missingness_map
    if (!is.null(mm) && nrow(mm) > 0 &&
        all(c("Site_ID", "Condition", "Detected", "Total") %in% colnames(mm))) {
      mrows <- mm[mm$Site_ID == pd$abs_id, , drop = FALSE]
      if (nrow(mrows) > 0)
        det_lines <- paste0("  ", mrows$Condition, ": ",
                            mrows$Detected, "/", mrows$Total, " replicates detected")
    }

    paste0(
      "PTM: ", pd$ptm_label, " (", pd$ptm_id, ")\n",
      "Gene: ", if (pd$gene != "") pd$gene else "N/A",
      "  |  Accession: ", pd$accession,
      {
        .alt <- rv$dia_ptm_master_cache$site_alt_accessions
        if (!is.null(.alt) && !is.null(pd$abs_id) && pd$abs_id %in% names(.alt) &&
            nzchar(.alt[[pd$abs_id]]))
          paste0(" (shared peptide; also maps to ",
                 gsub(";", ", ", .alt[[pd$abs_id]], fixed = TRUE),
                 " - quantities are identical, the protein assignment is an inference)")
        else ""
      },
      "  |  Site: ", pd$site,
      "  |  Modification: ", pd$modification, "\n",
      "Reference: ", pd$ref_cond,
      "  |  Protein data: ", if (!.with_global) "none (PTM only)" else if (pd$has_protein) "yes" else "no", "\n",
      "Abundances shown are measured values only - nothing is imputed in this view.\n",
      "Statistics: ",
      if (isTRUE(pd$de_available))
        "from the stored limma fit (empirical Bayes moderated, BH across all tested sites)."
      else if (isTRUE(pd$nmar_available))
        paste0("this site was not fitted by the linear model - it did not meet the ",
               "valid-value requirement in both compared conditions - and was assessed ",
               "for a non-random detection pattern instead (Fisher's exact test, ",
               "corrected within that set only).\n",
               paste(vapply(names(pd$nmar), function(tc) {
                 x <- pd$nmar[[tc]]
                 paste0("  ", tc, " vs ", pd$ref_cond, ": detected ", x$found_test,
                        " vs ", x$found_ref,
                        if (isTRUE(x$has_power) && !is.na(x$adj_pvalue))
                          paste0(", Fisher adj.p = ", signif(x$adj_pvalue, 3))
                        else if (!isTRUE(x$has_power))
                          paste0(" (no attainable significance at these replicate counts; ",
                                 "smallest possible two-sided p = ",
                                 signif(x$min_achievable_p, 3), ", so p is reported as NA)")
                        else " (p not available for this site)")
               }, character(1)), collapse = "\n"))
      else
        paste0("not available - no PTM DE result is stored for this reference condition. ",
               if (length(names(rv$dia_ptm_de_results)) > 0)
                 paste0("Contrasts currently available: ",
                        paste(names(rv$dia_ptm_de_results), collapse = "; "), ".")
               else "No PTM DE run has been completed for this dataset yet."),
      "\n",
      if (length(det_lines) > 0)
        paste0("\nDetection:\n", paste(det_lines, collapse = "\n"), "\n") else "",
      if (length(fc_lines) > 0) paste0("\nFold Changes:\n", paste(fc_lines, collapse = "\n")) else ""
    )
  })

  # ===========================================================================
  # EXPORT
  # ===========================================================================
  output$psg_download_csv <- downloadHandler(
    filename = function() {
      ptm_safe <- gsub("[^a-zA-Z0-9_]", "_", psg$plot_data$ptm_label)
      paste0("SinglePTM_", ptm_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(psg$plot_data)
      pd <- psg$plot_data

      lines_out <- character(0)

      # Header
      lines_out <- c(lines_out,
        paste0("# Single PTM Report: ", pd$ptm_label, " (", pd$ptm_id, ")"),
        paste0("# Gene: ", if (pd$gene != "") pd$gene else "N/A",
               "  |  Accession: ", pd$accession,
               "  |  Site: ", pd$site,
               "  |  Modification: ", pd$modification),
        paste0("# Reference Condition: ", pd$ref_cond),
        paste0("# Generated: ", Sys.time()),
        ""
      )

      # Panel A: Summarized Abundance
      lines_out <- c(lines_out, "## Panel A: Summarized Abundance", "")
      if (!is.null(pd$ptm_pld) && nrow(pd$ptm_pld) > 0) {
        lines_out <- c(lines_out, "# PTM Abundance")
        ptm_exp <- data.frame(
          Condition = as.character(pd$ptm_pld$Condition),
          BioReplicate = as.character(pd$ptm_pld$BioReplicate),
          Plex = if (!is.null(pd$ptm_pld$Plex)) as.character(pd$ptm_pld$Plex) else NA_character_,
          Abundance = round(pd$ptm_pld$Abundance, 6),
          Layer = "PTM",
          stringsAsFactors = FALSE
        )
        tc <- textConnection("tmp_lines", "w")
        write.csv(ptm_exp, tc, row.names = FALSE)
        close(tc)
        lines_out <- c(lines_out, tmp_lines, "")
      }
      if (!is.null(pd$prot_pld) && nrow(pd$prot_pld) > 0) {
        lines_out <- c(lines_out, "# Protein Abundance")
        prot_exp <- data.frame(
          Condition = as.character(pd$prot_pld$Condition),
          BioReplicate = as.character(pd$prot_pld$BioReplicate),
          Plex = if (!is.null(pd$prot_pld$Plex)) as.character(pd$prot_pld$Plex) else NA_character_,
          Abundance = round(pd$prot_pld$Abundance, 6),
          Layer = "Protein",
          stringsAsFactors = FALSE
        )
        tc <- textConnection("tmp_lines", "w")
        write.csv(prot_exp, tc, row.names = FALSE)
        close(tc)
        lines_out <- c(lines_out, tmp_lines, "")
      }

      # PTM-only cache: no protein layer, so no stoichiometric shape and only the
      # absolute PTM fold changes
      .with_global <- !isFALSE(pd$has_global)

      # Stoichiometric shape — data only; this view has no plot panel in the UI
      if (.with_global) {
        lines_out <- c(lines_out, "## Stoichiometric Shape (PTM - Protein)", "")
        if (!is.null(pd$stoich_df) && nrow(pd$stoich_df) > 0) {
          tc <- textConnection("tmp_lines", "w")
          write.csv(pd$stoich_df, tc, row.names = FALSE)
          close(tc)
          lines_out <- c(lines_out, tmp_lines, "")
        }
      }

      # Panel C: Fold Changes
      lines_out <- c(lines_out, "## Panel C: Mechanistic Fold Changes", "")
      .fc_main <- if (.with_global) pd$fc_adjusted else pd$fc_unadjusted
      if (is.null(.fc_main) || nrow(.fc_main) == 0) {
        lines_out <- c(lines_out,
          "# No differential expression result is stored for this contrast.",
          "# Run PTM DE for this comparison, then re-export.", "")
        writeLines(lines_out, file)
        logger(paste0("PTM Single: CSV exported (no DE result) for ", pd$ptm_label))
        return(invisible(NULL))
      }
      .fc_col <- function(tbl, col, n) {
        if (is.null(tbl) || !(col %in% colnames(tbl))) return(rep(NA_real_, n))
        tbl[[col]]
      }
      .n_fc <- nrow(.fc_main)
      fc_combined <- if (.with_global) data.frame(
        Contrast = pd$fc_adjusted$Label,
        Adjusted_log2FC   = round(pd$fc_adjusted$log2FC, 6),
        Adjusted_SE       = round(pd$fc_adjusted$SE, 6),
        Adjusted_adj.p    = signif(pd$fc_adjusted$adj.pvalue, 6),
        Unadjusted_log2FC = round(.fc_col(pd$fc_unadjusted, "log2FC", .n_fc), 6),
        Unadjusted_SE     = round(.fc_col(pd$fc_unadjusted, "SE", .n_fc), 6),
        Unadjusted_adj.p  = signif(.fc_col(pd$fc_unadjusted, "adj.pvalue", .n_fc), 6),
        Protein_log2FC    = round(.fc_col(pd$fc_protein, "log2FC", .n_fc), 6),
        Protein_SE        = round(.fc_col(pd$fc_protein, "SE", .n_fc), 6),
        Protein_adj.p     = signif(.fc_col(pd$fc_protein, "adj.pvalue", .n_fc), 6),
        stringsAsFactors = FALSE
      ) else data.frame(
        Contrast = .fc_main$Label,
        Unadjusted_log2FC = round(.fc_main$log2FC, 6),
        Unadjusted_SE     = round(.fc_main$SE, 6),
        Unadjusted_adj.p  = signif(.fc_main$adj.pvalue, 6),
        stringsAsFactors = FALSE
      )
      tc <- textConnection("tmp_lines", "w")
      write.csv(fc_combined, tc, row.names = FALSE)
      close(tc)
      lines_out <- c(lines_out, tmp_lines)

      writeLines(lines_out, file)
      logger(paste0("PTM Single: CSV exported for ", pd$ptm_label))
    }
  )

  output$psg_save_rds <- downloadHandler(
    filename = function() {
      ptm_safe <- gsub("[^a-zA-Z0-9_]", "_", psg$plot_data$ptm_label)
      paste0("SinglePTM_", ptm_safe, "_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(psg$plot_data)
      saveRDS(list(
        plot_data  = psg$plot_data,
        gc_result  = psg$gc_result,
        reference  = psg$last_reference,
        conditions = psg$conditions,
        cond_order = psg$cond_order,
        timestamp  = Sys.time(),
        provenance = reactiveValuesToList(rv$provenance)
      ), file)
      logger("PTM Single: RDS exported")
    }
  )

}
