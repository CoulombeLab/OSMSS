# server_heatmap.R ==============================================================
# Server logic for Heatmap module.
# Two LF DDA pipelines (reads from wizard cache — rv$lf_master_cache):
#   1. IP (Immunoprecipitation) - invariant-protein median normalization via ProDA
#   2. Global (Whole Tissue)    - ProDA median normalization
# TMT: reads rv$tmt_limma_cache (limma route) or rv$tmt_master_cache (legacy).
# DIA: reads rv$dia_master_cache.
# Output: protein x sample matrix — row z-score + hierarchical clustering.
# Called from server.R as: server_heatmap(input, output, session, rv, logger)
# ==============================================================================

server_heatmap <- function(input, output, session, rv, logger) {

  # =========================================================================
  # Reactive storage
  # =========================================================================
  hm <- reactiveValues(
    # Results
    abundance_mat    = NULL,   # numeric matrix (rows = proteins, cols = samples)
    stats_df         = NULL,   # data.frame with GeneName, Accession, Variance
    gene_names       = NULL,   # character vector aligned with abundance_mat rows
    accessions       = NULL,
    all_genes        = NULL,   # unique gene names for selectize
    heatmap_obj      = NULL,   # list(plot, height, n_genes, n_samples, mode)
    render_count     = 0L,
    # Invariant reference proteins used for normalization
    ebn_refs         = NULL,   # character vector of accessions, or NULL
    ebn_acc_to_gene  = NULL,   # named character: accession -> gene name (for display)
    # Custom column order (NULL = default clustering order)
    col_order        = NULL,   # character vector of column names in user-chosen order
    # Column-to-group mapping for averaging (named char vec, or NULL)
    col_to_group     = NULL,   # e.g. c("Rep1" = "Cond1", "Rep2" = "Cond1", ...)
    # Pre-imputation matrix (after channel selection and batch correction), or
    # NULL. Used for variance ranking only -- never displayed or exported.
    observed_mat     = NULL
  )
  register_experiment_reset(rv, hm)
  methodsTextServer("hm_methods", provenance = rv$provenance, section = "heatmap")

  # =========================================================================
  # Shared helpers
  # =========================================================================

  # Parse protein-level data from PD CSV
  .parse_protein_df <- function() {
    req(rv$pd_raw_lines, rv$pd_protein_line_idx)
    header_line   <- rv$pd_raw_lines[1]
    protein_lines <- rv$pd_raw_lines[rv$pd_protein_line_idx]
    text_block    <- paste(c(header_line, protein_lines), collapse = "\n")
    tryCatch(
      read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, na.strings = c("", "NA")),
      error = function(e) {
        logger(paste0("!!! [Heatmap] Protein CSV parse error: ", e$message))
        NULL
      }
    )
  }

  # Build row labels: GeneName fallback to Accession
  .row_labels <- function(gene_names, accessions) {
    labels <- ifelse(gene_names != "", gene_names, accessions)
    make.unique(labels, sep = "_")
  }

  # Require a package, showing notification if missing
  require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      showNotification(
        paste0("Package '", pkg, "' is not installed. ",
               "Install via: BiocManager::install('", pkg, "')"),
        type = "error", duration = 10)
      return(FALSE)
    }
    TRUE
  }

  # Build plotly heatmap -- delegates to shared build_heatmap_plotly() in global.R
  .build_heatmap <- function(mat, title_txt = "Heatmap", mode_label = "",
                             col_order_mode = "cluster", scale_mode = "zscore",
                             na_display = "zero") {
    build_heatmap_plotly(mat, title_txt, mode_label,
                        col_order_mode = col_order_mode, scale_mode = scale_mode,
                        na_display = na_display,
                        log_fn = logger)
  }

  # =========================================================================
  #  GENERATE HEATMAP -- dispatcher
  # =========================================================================
  .hm_dispatch <- function() {
    pipeline <- isolate(rv$lf_pipeline)
    if (is.null(pipeline)) pipeline <- "ip"
    logger(paste0("Heatmap: Generating (pipeline = ", pipeline, ")..."))

    grps <- rv$lf_cond_groups
    if (length(grps) < 2) {
      showNotification("Define at least 2 condition groups first (Step 2).",
                       type = "warning"); return()
    }
    cache <- rv$lf_master_cache
    meta  <- rv$lf_master_cache_meta
    if (is.null(cache) || is.null(meta)) {
      showNotification(
        "Please complete the Upload Wizard (Summarize Data step) before generating the heatmap.",
        type = "warning", duration = 8)
      return()
    }
    logger(paste0("HM [", pipeline, "]: Using cached summarization"))
    .run_heatmap_from_cache(cache)
  }

  # =========================================================================
  # TMT Heatmap (unified route  -" called when rv$data_source == "tmt")
  # =========================================================================
  .run_hm_tmt <- function() {
    if (is.null(rv$tmt_limma_cache)) {
      showNotification("Please run TMT Summarize Data first (Upload Data tab).", type = "warning")
      return()
    }

    missing_mode <- isolate(input$hm_missing_mode) %||% "omit"
    logger(paste0("Heatmap [TMT]: Starting (missing=", missing_mode,
                  ", route=", rv$tmt_analysis_route %||% "legacy", ")"))

    shinyjs::disable("hm_generate")
    on.exit(shinyjs::enable("hm_generate"), add = TRUE)

    lc  <- rv$tmt_limma_cache
    ann <- lc$annotation_map
    pm  <- lc$protein_matrix   # proteins × channels

    withProgress(message = "Running TMT Heatmap (Limma)...", value = 0, {
      incProgress(0.30, detail = "Building matrix...")

      use_chans <- ann$ColumnName[!ann$Condition %in% "Norm"]
      use_chans <- intersect(use_chans, colnames(pm))
      hm_mat    <- pm[, use_chans, drop = FALSE]   # proteins × samples

      proteins  <- rownames(hm_mat)
      n_omitted <- 0L; n_imputed <- 0L
      # render_tmt_heatmap_methods guards its LOD sentence on lod_value, and the
      # computed floor used to be discarded. Recorded on both paths.
      hm_lod_value <- NA_real_

      incProgress(0.50, detail = "Batch correction...")
      sample_groups <- ann$Condition[match(colnames(hm_mat), ann$ColumnName)]
      mixture_vec   <- ann$Mixture[match(colnames(hm_mat), ann$ColumnName)]
      n_plexes      <- length(unique(mixture_vec[!is.na(mixture_vec)]))
      did_batch_hm  <- FALSE
      n_bc_not_est_hm <- 0L
      n_single_plex_hm <- 0L
      # See tmt_plex_correctable(): with the Norm channels dropped, a plex holding
      # one sample cannot have its offset separated from that sample's own value.
      plex_ok_hm <- tmt_plex_correctable(mixture_vec)
      if (n_plexes >= 2 && !plex_ok_hm)
        logger(paste0("Heatmap [TMT]: plex batch correction skipped  -  at least one ",
                      "plex holds a single non-reference channel, so its offset cannot ",
                      "be separated from that sample (correction would set each such ",
                      "sample to its condition mean). Values are shown as summarized."))

      if (n_plexes >= 2 && plex_ok_hm && requireNamespace("limma", quietly = TRUE)) {
        # The flag is set on the SUCCESS path only. It used to be set
        # unconditionally after the tryCatch, so a caught failure still recorded
        # batch_corrected = TRUE in provenance and the methods text asserted a
        # correction that never happened.
        .hm_uncorrected <- hm_mat
        hm_mat <- tryCatch({
          .hm_bc <- limma::removeBatchEffect(
            hm_mat,
            batch  = factor(mixture_vec),
            design = model.matrix(~ factor(sample_groups))
          )
          did_batch_hm <- TRUE
          .hm_bc
        },
          error = function(e) {
            logger(paste0("Heatmap [TMT limma]: removeBatchEffect failed: ", e$message,
                          "  -  heatmap shown WITHOUT plex batch correction."))
            did_batch_hm <<- FALSE
            hm_mat
          }
        )
        if (isTRUE(did_batch_hm)) {
          # A protein seen in only one plex has no estimable plex coefficient
          # and removeBatchEffect can return NA for its whole row. Correcting
          # before imputation (the order this block now runs in) exposes that,
          # so restore those rows uncorrected rather than losing them.
          .bc_lost <- rowSums(!is.na(hm_mat)) == 0L &
                      rowSums(!is.na(.hm_uncorrected)) > 0L
          n_bc_not_est_hm <- sum(.bc_lost)
          if (any(.bc_lost)) {
            hm_mat[.bc_lost, ] <- .hm_uncorrected[.bc_lost, ]
            logger(paste0("Heatmap [TMT]: ", n_bc_not_est_hm,
                          " proteins were not estimable for the plex term ",
                          "(present in a single plex) and are shown uncorrected."))
          }
          # A protein observed in fewer than two plexes carries no between-plex
          # contrast, so its batch coefficient is zero and the row passes
          # through unadjusted rather than being voided. The all-NA guard above
          # does not catch that, so count it separately: after correction the
          # matrix holds two populations on different footings, and the methods
          # text should not imply the adjustment reached every row.
          .hm_plex_n <- Reduce(`+`, lapply(
            unique(mixture_vec[!is.na(mixture_vec)]),
            function(px) as.integer(rowSums(
              !is.na(.hm_uncorrected[, which(mixture_vec == px), drop = FALSE])) > 0L)))
          n_single_plex_hm <- sum(.hm_plex_n < 2L)
          if (n_single_plex_hm > 0)
            logger(paste0("Heatmap [TMT]: ", n_single_plex_hm,
                          " proteins were observed in a single plex and therefore ",
                          "carry no plex adjustment."))
          logger(paste0("Heatmap [TMT limma]: plex batch correction applied (",
                        n_plexes, " plexes)."))
        }
      }

      # The matrix as it stands after channel selection and batch correction but
      # before any missingness handling. .finalize_heatmap ranks on these values
      # so that a protein absent from a whole condition is not promoted into the
      # display by the spread its substituted values create. Handing over the
      # values rather than a missingness mask means the condition-averaged case
      # needs no threshold on a collapsed mask.
      hm$observed_mat <- hm_mat

      if (missing_mode == "lod") {
        sample_conds <- ann$Condition[match(colnames(hm_mat), ann$ColumnName)]
        keep_rows    <- rep(FALSE, nrow(hm_mat))
        for (cond in unique(sample_conds)) {
          cond_cols <- which(sample_conds == cond)
          if (length(cond_cols) >= 2) {
            n_valid   <- rowSums(!is.na(hm_mat[, cond_cols, drop = FALSE]))
            keep_rows <- keep_rows | (n_valid >= 2)
          }
        }
        n_omitted      <- sum(!keep_rows)
        hm_mat         <- hm_mat[keep_rows, , drop = FALSE]
        proteins       <- proteins[keep_rows]
        # Per-COLUMN (per-sample) 1st-percentile substitution, mirroring the DIA
        # route. The 1st percentile of a channel's own observed distribution is
        # that channel's detection floor; a single floor over the whole matrix
        # has no standard equivalent and cannot be named in a methods section.
        # This is imputeLCMD::impute.MinDet with q = 0.01.
        n_imputed      <- sum(is.na(hm_mat))
        lod_vec        <- apply(hm_mat, 2, quantile, probs = 0.01, na.rm = TRUE)
        names(lod_vec) <- colnames(hm_mat)
        # A channel observed nowhere yields a non-finite floor; substituting it
        # would leave NAs in a matrix the plotting path expects to be complete.
        if (any(!is.finite(lod_vec)))
          lod_vec[!is.finite(lod_vec)] <- median(lod_vec[is.finite(lod_vec)])
        for (ci in seq_len(ncol(hm_mat))) {
          na_m <- is.na(hm_mat[, ci])
          if (any(na_m)) hm_mat[na_m, ci] <- lod_vec[ci]
        }
        hm_lod_value   <- as.numeric(median(lod_vec))
        # Diagnostic. Within a plex all channels are quantified from the same
        # MS2 scans, so after median centering the per-channel floors should sit
        # close together; a wide spread is information the user wants.
        logger(paste0("Heatmap [TMT]: per-sample LOD floors  -  min=",
                      round(min(lod_vec), 3), ", median=", round(median(lod_vec), 3),
                      ", max=", round(max(lod_vec), 3), ", max-min=",
                      round(max(lod_vec) - min(lod_vec), 3), " log2 units (ratio=",
                      round(max(lod_vec) / min(lod_vec), 4), ")"))
      } else {
        complete_rows <- complete.cases(hm_mat)
        n_omitted     <- sum(!complete_rows)
        hm_mat        <- hm_mat[complete_rows, , drop = FALSE]
        proteins      <- proteins[complete_rows]
        # Nothing is substituted in omit mode, so ranking on the post-filter
        # matrix and on the observed values are the same computation.
        hm$observed_mat <- NULL
      }

      if (nrow(hm_mat) < 2) {
        showNotification("Too few proteins after filtering.", type = "error"); return()
      }

      incProgress(0.70, detail = "Gene names...")
      gene_names <- tmt_cache_gene_labels(lc, proteins)

      cond_names <- sort(unique(sample_groups[!is.na(sample_groups)]))
      hm$col_to_group <- setNames(
        ann$Condition[match(colnames(hm_mat), ann$ColumnName)],
        colnames(hm_mat)
      )

      rv$provenance$heatmap <- list(
        source           = "tmt",
        pipeline         = "tmt",
        pipeline_label   = "TMT quantitative proteomics pipeline",
        route            = rv$tmt_analysis_route,
        n_proteins_input = nrow(hm_mat) + n_omitted,
        missing_mode     = missing_mode,
        lod_value        = hm_lod_value,
        lod_scope        = if (missing_mode == "lod") "per_sample_1st_percentile" else NA_character_,
        n_omitted        = n_omitted,
        n_imputed        = n_imputed,
        batch_corrected  = did_batch_hm,
        n_proteins_batch_not_estimable = n_bc_not_est_hm,
        n_proteins_single_plex = n_single_plex_hm,
        condition_names  = cond_names,
        software_used    = "limma + plotly"
      )

      incProgress(0.90, detail = "Building heatmap...")
      hm$na_display <- "zero"   # unchanged TMT behaviour
      n_conds    <- length(cond_names)
      mode_label <- paste0("TMT Limma (", n_conds, " conditions)")
      .finalize_heatmap(hm_mat, gene_names, proteins, mode_label = mode_label)
      incProgress(1, detail = "Done!")
      logger(paste0("Heatmap [TMT limma]: Done  -  ", nrow(hm_mat), " proteins x ",
                    ncol(hm_mat), " samples"))
    })
  }


  # =========================================================================
  # DIA Heatmap (unified route  -" called when rv$data_source == "dia")
  # =========================================================================
  .run_hm_dia <- function() {
    if (is.null(rv$dia_master_cache)) {
      showNotification("Please complete the DIA Upload Wizard first.", type = "warning")
      return()
    }

    shinyjs::disable("hm_generate")
    on.exit(shinyjs::enable("hm_generate"), add = TRUE)

    mat  <- rv$dia_master_cache$protein_matrix
    grps <- rv$dia_cond_groups
    cond_map <- setNames(
      rep(names(grps), vapply(grps, length, integer(1))),
      unlist(grps)
    )
    use_samps <- intersect(colnames(mat), unlist(grps))
    mat <- mat[, use_samps, drop = FALSE]
    mat_full <- mat  # keep full matrix for gene picker

    withProgress(message = "Building DIA Heatmap...", value = 0.1, {
      missing_mode <- isolate(input$hm_missing_mode) %||% "omit"
      n_omitted <- 0L; n_imputed <- 0L
      if (missing_mode == "omit") {
        keep      <- complete.cases(mat)
        n_omitted <- sum(!keep)
        mat       <- mat[keep, , drop = FALSE]
      } else {
        # Per-COLUMN (per-sample) 1st-percentile substitution. The 1st percentile
        # of a sample's own observed distribution is a sample-level detection
        # floor; the previous per-row version substituted each protein's own
        # minimum, which is not a detection floor.
        n_imputed <- sum(is.na(mat))
        for (ci in seq_len(ncol(mat))) {
          na_m <- is.na(mat[, ci])
          if (any(na_m)) {
            lod_v <- quantile(mat[, ci], 0.01, na.rm = TRUE)
            mat[na_m, ci] <- lod_v
          }
        }
      }
      incProgress(0.2)

      if (isTRUE(rv$dia_batch_corrected) && !is.null(rv$dia_batch_groups) &&
          requireNamespace("limma", quietly = TRUE)) {
        bat_map   <- setNames(
          rep(names(rv$dia_batch_groups), vapply(rv$dia_batch_groups, length, integer(1))),
          unlist(rv$dia_batch_groups)
        )
        batch_vec <- bat_map[colnames(mat)]
        cond_vec  <- cond_map[colnames(mat)]
        cond_vec[is.na(cond_vec)] <- "Unknown"
        mat <- tryCatch(
          limma::removeBatchEffect(mat, batch = batch_vec,
                                   design = model.matrix(~ factor(cond_vec))),
          error = function(e) mat)
      }
      incProgress(0.2)

      top_n    <- as.integer(isolate(input$hm_n_top) %||% 50L)
      # Rank on observed cells only. Imputed cells are masked back to NA first,
      # so a protein's rank reflects measured variation rather than the spread
      # introduced by substituting a constant for its missing values.
      .obs_na  <- is.na(mat_full[rownames(mat), colnames(mat), drop = FALSE])
      .obs_mat <- mat
      .obs_mat[.obs_na] <- NA_real_
      row_vars <- apply(.obs_mat, 1, var, na.rm = TRUE)
      row_vars[is.na(row_vars)] <- 0
      top_rows <- names(sort(row_vars, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))])
      mat_sel  <- mat[top_rows, , drop = FALSE]

      # If user has selected specific genes AND LOD imputation is active, include
      # them even if they were dropped by the complete.cases() filter.
      # In Omit-NA mode, proteins with missing values are intentionally excluded.
      .sel_genes_now <- isolate(input$hm_gene_filter)
      if (length(.sel_genes_now) > 0 && missing_mode != "omit") {
        acc_to_gene_full <- rv$dia_master_cache$acc_to_gene %||% character(0)
        full_gn <- acc_to_gene_full[rownames(mat_full)]
        full_gn[is.na(full_gn)] <- rownames(mat_full)[is.na(full_gn)]
        sel_accs <- rownames(mat_full)[full_gn %in% .sel_genes_now]
        missing_accs <- setdiff(sel_accs, rownames(mat_sel))
        if (length(missing_accs) > 0) {
          add_mat <- mat_full[missing_accs, use_samps, drop = FALSE]
          # Column LODs must come from the FULL matrix: add_mat holds only the
          # handful of rescued proteins, and a 1st percentile over 2-3 values is
          # meaningless.
          .rescue_lod <- apply(mat_full[, use_samps, drop = FALSE], 2,
                               quantile, probs = 0.01, na.rm = TRUE)
          for (ci in seq_len(ncol(add_mat))) {
            na_m <- is.na(add_mat[, ci])
            if (any(na_m)) add_mat[na_m, ci] <- .rescue_lod[ci]
          }
          mat_sel <- rbind(mat_sel, add_mat)
        }
      }
      incProgress(0.2)

      acc_to_gene <- rv$dia_master_cache$acc_to_gene %||% character(0)
      gene_names  <- acc_to_gene[rownames(mat_sel)]
      gene_names[is.na(gene_names)] <- rownames(mat_sel)[is.na(gene_names)]
      proteins    <- rownames(mat_sel)

      hm$col_to_group <- setNames(cond_map[colnames(mat_sel)], colnames(mat_sel))
      hm$full_mat     <- mat_full   # full matrix for gene-picker rescue in Refresh View
      hm$acc_to_gene  <- rv$dia_master_cache$acc_to_gene %||% character(0)

      rv$provenance$heatmap <- list(
        source           = "dia",
        n_proteins_input = nrow(mat),
        top_n            = nrow(mat_sel),
        missing_mode     = missing_mode,
        n_omitted        = n_omitted,
        n_imputed        = n_imputed,
        batch_corrected  = isTRUE(rv$dia_batch_corrected),
        condition_names  = names(grps),
        software_used    = "plotly"
      )

      incProgress(0.2)
      hm$na_display <- "zero"   # unchanged DIA behaviour
      mode_label <- "DIA (log2 Abundance)"
      .finalize_heatmap(mat_sel, gene_names, proteins,
                        mode_label = mode_label)
      # Override all_genes with full protein set (pre-missing-filter) for gene picker
      all_acc_genes <- acc_to_gene[rownames(mat_full)]
      all_acc_genes[is.na(all_acc_genes)] <- rownames(mat_full)[is.na(all_acc_genes)]
      hm$all_genes <- sort(unique(all_acc_genes[nzchar(all_acc_genes)]))
      updateSelectizeInput(session, "hm_gene_filter",
                           choices = hm$all_genes, selected = character(0), server = TRUE)
      incProgress(0.2)
      logger(paste0("[Heatmap DIA] Done: ", nrow(mat_sel), " proteins x ",
                    ncol(mat_sel), " samples."))
    })
  }

  # =========================================================================
  # GENERATE HEATMAP button handler  -" source router + LF DDA path
  # =========================================================================
  observeEvent(input$hm_generate, {
    # Route to TMT or DIA path if appropriate data source is active
    src <- rv$data_source
    if (identical(src, "tmt")) {
      .run_hm_tmt()
      return()
    } else if (identical(src, "dia")) {
      .run_hm_dia()
      return()
    }

    # LF DDA path (unchanged)
    if (is.null(rv$lf_master_cache) && is.null(rv$pd_raw_lines)) {
      showNotification("Please run the Upload Wizard first to summarize data.", type = "warning")
      return()
    }

    shinyjs::disable("hm_generate")
    on.exit(shinyjs::enable("hm_generate"), add = TRUE)

    pipeline <- isolate(rv$lf_pipeline)
    if (is.null(pipeline)) pipeline <- "ip"

    # Check if re-summarization needed and confirm with user
    cache <- rv$lf_master_cache; meta <- rv$lf_master_cache_meta
    needs_resummarize <- is.null(cache) || is.null(meta) ||
      !identical(meta$pipeline, pipeline) ||
      (identical(pipeline, "ip") && !identical(meta$ip_sub, isolate(rv$lf_ip_sub)))

    if (needs_resummarize && isTRUE(rv$lf_cache_stale)) {
        showModal(modalDialog(
          title = tags$span(icon("rotate"), " Re-Summarization Required"),
          tags$div(class = "alert alert-info mb-0",
            tags$p(
              tags$b("A new analysis pipeline has been selected."),
              tags$br(),
              "The data will be re-summarized with the current settings before ",
              "generating the heatmap. This may take a few minutes."
            ),
            tags$p(class = "text-muted small mb-0",
              "Pipeline: ", tags$b(pipeline),
              if (identical(pipeline, "ip")) paste0(" (", isolate(rv$lf_ip_sub), ")") else "")
          ),
          size = "m", easyClose = FALSE,
          footer = tagList(
            modalButton("Cancel"),
            actionButton("hm_confirm_resummarize", "Continue",
                         class = "btn-primary", icon = icon("play"))
          )
        ))
        return()
    }

    .hm_dispatch()
  })

  observeEvent(input$hm_confirm_resummarize, {
    removeModal()
    rv$lf_cache_stale <- FALSE
    .hm_dispatch()
  })


  # =========================================================================
  #  Select top genes & build heatmap from results
  # =========================================================================
  .finalize_heatmap <- function(abundance_mat, gene_names, accessions,
                                mode_label = "",
                                update_gene_list = TRUE) {
    if (is.null(abundance_mat) || nrow(abundance_mat) < 2) {
      showNotification("Too few proteins with valid data.", type = "error")
      return()
    }

    # Store full (per-replicate) results
    row_labels <- .row_labels(gene_names, accessions)
    rownames(abundance_mat) <- row_labels
    hm$abundance_mat <- abundance_mat
    hm$gene_names    <- gene_names
    hm$accessions    <- accessions
    # Only overwrite the full gene list on initial generation (not on toggle/refresh),
    # so that pipelines that set a broader hm$all_genes after the call (e.g. DIA
    # overrides with mat_full) are not clobbered when the user re-renders.
    if (update_gene_list) {
      hm$all_genes <- sort(unique(row_labels[row_labels != ""]))
    }

    # --- Averaging by condition (if toggled ON) ---
    avg_mode <- isTRUE(isolate(input$hm_avg_by_cond))
    ctg      <- hm$col_to_group

    if (avg_mode && !is.null(ctg) && length(ctg) > 0) {
      work_mat <- heatmap_group_average(abundance_mat, ctg)
      if (is.null(work_mat) || ncol(work_mat) < 2) {
        showNotification("Cannot average: fewer than 2 groups with data.",
                         type = "warning")
        work_mat <- abundance_mat
      }
    } else {
      work_mat <- abundance_mat
    }

    # Build enriched stats_df from the *working* matrix
    row_vars <- apply(work_mat, 1, var, na.rm = TRUE)
    row_vars[is.na(row_vars)] <- 0

    # Selection ranks on the OBSERVED, pre-imputation values where the caller
    # supplied them (TMT LOD mode). Ranking on the post-imputation matrix
    # promotes a protein absent from a whole condition on the spread its
    # substituted values create. Taking the values rather than a missingness
    # mask means the condition-averaged case needs no threshold on a collapsed
    # mask: each condition mean is simply formed from its observed channels.
    # stats_df, the displayed cell values, the colour scale and the export path
    # all keep the post-imputation matrix; this changes which proteins are
    # selected and nothing else. A protein with one imputed channel in a
    # condition still displays that condition's mean; it just does not gain
    # rank from the imputer's spread.
    # Absent observed_mat (LF, IP/SAINT, DIA, or omit mode) -> unchanged.
    #
    # A condition mean used for variance ranking is only formed when at least
    # 50% of that condition's channels were observed. Matches .tcl_min_rep_frac
    # on the DIA and TMT trends routes, so one coverage rule governs the app.
    .hm_min_obs_frac <- 0.5
    .hm_rank_vars    <- row_vars
    .hm_rank_basis   <- NA_character_
    .hm_n_unrankable <- NA_integer_
    .hm_obs <- hm$observed_mat
    if (!is.null(.hm_obs) && is.matrix(.hm_obs) &&
        identical(rownames(.hm_obs), accessions) &&
        identical(colnames(.hm_obs), colnames(abundance_mat))) {
      .rank_mat <- .hm_obs
      rownames(.rank_mat) <- row_labels
      if (avg_mode && !is.null(ctg) && length(ctg) > 0 &&
          !identical(colnames(.rank_mat), colnames(work_mat))) {
        # Condition means from observed channels only. Below the coverage rule
        # the mean is NA rather than being formed from one surviving channel.
        .rk_grp <- ctg[colnames(.rank_mat)]
        .rk_grp[is.na(.rk_grp)] <- "Ungrouped"
        .rk_avg <- matrix(NA_real_, nrow = nrow(.rank_mat), ncol = ncol(work_mat),
                          dimnames = list(rownames(.rank_mat), colnames(work_mat)))
        for (.g in colnames(work_mat)) {
          .g_cols <- colnames(.rank_mat)[.rk_grp == .g]
          if (length(.g_cols) == 0) next
          .g_sub  <- .rank_mat[, .g_cols, drop = FALSE]
          .g_nobs <- rowSums(!is.na(.g_sub))
          .g_mu   <- rowMeans(.g_sub, na.rm = TRUE)
          .g_mu[(.g_nobs / length(.g_cols)) < .hm_min_obs_frac] <- NA_real_
          .rk_avg[, .g] <- .g_mu
        }
        .rank_mat <- .rk_avg
      }
      if (identical(dim(.rank_mat), dim(work_mat))) {
        # Fewer than two observed values yields no variance: rank last rather
        # than error, and count how many proteins that applies to.
        .rk_nobs         <- rowSums(!is.na(.rank_mat))
        .hm_n_unrankable <- sum(.rk_nobs < 2L)
        .hm_rank_vars    <- apply(.rank_mat, 1, var, na.rm = TRUE)
        .hm_rank_vars[is.na(.hm_rank_vars) | .rk_nobs < 2L] <- 0
        .hm_rank_basis   <- "observed_pre_imputation"
      }
    }

    stats_df <- data.frame(
      GeneName  = gene_names,
      Accession = accessions,
      Variance  = round(row_vars, 4),
      stringsAsFactors = FALSE
    )

    # Append raw (pre-z-score) abundance values as additional columns
    abundance_df <- as.data.frame(work_mat, stringsAsFactors = FALSE)
    colnames(abundance_df) <- colnames(work_mat)
    stats_df <- cbind(stats_df, abundance_df)

    # hm$stats_df is set below after sub_mat and scale mode are finalized

    # Update gene filter choices
    updateSelectizeInput(session, "hm_gene_filter",
                         choices = hm$all_genes, selected = character(0), server = TRUE)

    # Select genes to display  -" ranked by working-matrix variance
    n_top <- isolate(input$hm_n_top)
    if (is.null(n_top) || is.na(n_top) || n_top < 1) n_top <- 50
    selected_genes <- isolate(input$hm_gene_filter)

    # Rescue: if user selected specific genes not present in abundance_mat,
    # pull them from hm$full_mat (if available) and impute NAs at per-row LOD.
    .hm_mm_now <- isolate(input$hm_missing_mode) %||% "omit"
    if (length(selected_genes) > 0 && !is.null(hm$full_mat) &&
        !identical(.hm_mm_now, "omit")) {
      .lookup_gn <- if (!is.null(hm$acc_to_gene) && length(hm$acc_to_gene) > 0) {
        hm$acc_to_gene
      } else if (!is.null(hm$ebn_acc_to_gene) && length(hm$ebn_acc_to_gene) > 0) {
        hm$ebn_acc_to_gene
      } else {
        setNames(hm$gene_names, hm$accessions)
      }
      full_gn_mapped <- .lookup_gn[rownames(hm$full_mat)]
      full_gn_mapped[is.na(full_gn_mapped)] <- rownames(hm$full_mat)[is.na(full_gn_mapped)]
      missing_genes <- setdiff(selected_genes, row_labels)
      if (length(missing_genes) > 0) {
        rescue_accs <- rownames(hm$full_mat)[full_gn_mapped %in% missing_genes]
        if (length(rescue_accs) > 0) {
          add_cols <- intersect(colnames(hm$full_mat), colnames(abundance_mat))
          add_mat  <- hm$full_mat[rescue_accs, add_cols, drop = FALSE]
          .rescue_lod <- apply(hm$full_mat[, add_cols, drop = FALSE], 2,
                               quantile, probs = 0.01, na.rm = TRUE)
          for (ci in seq_len(ncol(add_mat))) {
            na_m <- is.na(add_mat[, ci])
            if (any(na_m)) add_mat[na_m, ci] <- .rescue_lod[ci]
          }
          add_gn   <- full_gn_mapped[rescue_accs]
          abundance_mat <- rbind(abundance_mat, add_mat)
          gene_names    <- c(gene_names, add_gn)
          accessions    <- c(accessions, rescue_accs)
          row_labels    <- .row_labels(gene_names, accessions)
          rownames(abundance_mat) <- row_labels
          # Recompute work_mat to include rescued rows
          if (avg_mode && !is.null(ctg) && length(ctg) > 0) {
            work_mat <- heatmap_group_average(abundance_mat, ctg)
            if (is.null(work_mat) || ncol(work_mat) < 2) work_mat <- abundance_mat
          } else {
            work_mat <- abundance_mat
          }
          row_vars <- apply(work_mat, 1, var, na.rm = TRUE)
          row_vars[is.na(row_vars)] <- 0
          # Rescued rows change the row set, so observed_mat no longer aligns:
          # fall back to post-imputation ranking for this call.
          .hm_rank_vars    <- row_vars
          .hm_rank_basis   <- NA_character_
          .hm_n_unrankable <- NA_integer_
        }
      }
    }

    if (length(selected_genes) > 0) {
      keep_idx <- which(row_labels %in% selected_genes)
    } else {
      ranked   <- order(.hm_rank_vars, decreasing = TRUE)
      keep_idx <- head(ranked, n_top)
    }

    if (length(keep_idx) == 0) {
      showNotification("No genes match the filter.", type = "warning")
      return()
    }
    if (length(keep_idx) < 2) keep_idx <- rep(keep_idx, 2)  # plotly needs >=2 rows

    sub_mat <- work_mat[keep_idx, , drop = FALSE]

    # Remove rows that are entirely NA
    all_na <- apply(sub_mat, 1, function(r) all(is.na(r)))
    sub_mat <- sub_mat[!all_na, , drop = FALSE]

    if (nrow(sub_mat) < 2) {
      showNotification("Too few proteins with data after filtering.", type = "error")
      return()
    }

    # --- Heatmap scale mode ---
    hm_scale_mode <- isolate(input$hm_scale_mode) %||% "zscore"

    # --- When averaging by condition AND z-score, z-score per-replicate FIRST, then average ---
    if (identical(hm_scale_mode, "zscore") && avg_mode && !is.null(ctg) && length(ctg) > 0) {
      full_sub <- abundance_mat[rownames(sub_mat), , drop = FALSE]
      rm_vals  <- rowMeans(full_sub, na.rm = TRUE)
      rs_vals  <- apply(full_sub, 1, sd, na.rm = TRUE)
      ok       <- !is.na(rs_vals) & rs_vals > 0
      if (sum(ok) >= 2) {
        z_full <- (full_sub[ok, , drop = FALSE] - rm_vals[ok]) / rs_vals[ok]
        avg_z  <- heatmap_group_average(z_full, ctg)
        if (!is.null(avg_z) && nrow(avg_z) >= 2) {
          sub_mat       <- avg_z
          hm_scale_mode <- "zscore_pre"   # z-score already applied; keep z-score display styling
        }
      }
    } else if (identical(hm_scale_mode, "abundance")) {
      hm_scale_mode <- "abundance"   # pass through as-is (no z-score)
    }

    # Which z-score ordering actually ran. These are different quantities:
    # averaging first shrinks the SD and enlarges the resulting z-scores.
    .hm_z_order <- if (identical(hm_scale_mode, "zscore_pre")) {
      "z-score per replicate, then averaged by condition"
    } else if (identical(hm_scale_mode, "abundance")) {
      NA_character_
    } else if (isTRUE(avg_mode) && !is.null(ctg) && length(ctg) > 0) {
      "averaged by condition, then z-scored"
    } else {
      "z-score per sample column"
    }

    n_showing <- nrow(sub_mat)

    # Determine column ordering mode
    user_col_order <- hm$col_order

    # --- Update provenance display fields ---
    sel_genes_ui <- isolate(input$hm_gene_filter)
    .hm_sel_method  <- if (length(sel_genes_ui) > 0) "gene_filter" else "top_variance"
    # Mirror what is actually passed to build_heatmap_plotly() below: only
    # "manual" or "cluster" are ever used, so a "grouped" report was false.
    .hm_col_mode    <- if (!is.null(user_col_order) && length(user_col_order) > 0 &&
                           length(intersect(user_col_order, colnames(sub_mat))) > 0) {
      "manual"
    } else {
      "hierarchical (ward.D2, euclidean)"
    }
    existing_prov <- tryCatch(
      if (is.reactivevalues(rv$provenance$heatmap))
        isolate(reactiveValuesToList(rv$provenance$heatmap))
      else
        isolate(rv$provenance$heatmap),
      error = function(e) list()
    )
    if (length(existing_prov) > 0) {
      rv$provenance$heatmap$n_proteins_displayed <- as.numeric(n_showing)
      rv$provenance$heatmap$selection_method     <- .hm_sel_method
      rv$provenance$heatmap$top_n                <- if (.hm_sel_method == "top_variance") as.numeric(n_top) else numeric()
      rv$provenance$heatmap$gene_filter          <- if (.hm_sel_method == "gene_filter") sel_genes_ui else character()
      rv$provenance$heatmap$avg_by_condition     <- isTRUE(avg_mode)
      # Only written when ranking actually used the observed values, so the LF,
      # IP/SAINT and DIA ledger entries are left exactly as they were.
      if (!is.na(.hm_rank_basis)) {
        rv$provenance$heatmap$rank_basis            <- .hm_rank_basis
        rv$provenance$heatmap$rank_min_obs_frac     <- .hm_min_obs_frac
        rv$provenance$heatmap$n_proteins_unrankable <- .hm_n_unrankable
      }
      rv$provenance$heatmap$row_scaling           <- if (identical(hm_scale_mode, "abundance")) "normalized-abundance" else "z-score"
      rv$provenance$heatmap$z_score_order         <- .hm_z_order
      rv$provenance$heatmap$col_clustering       <- .hm_col_mode
      # Written on every route, not just DIA. On TMT the key used to be absent
      # and render_tmt_heatmap_methods fell back to asserting hierarchical
      # clustering unconditionally, including when it had been skipped. Refined
      # below from what build_heatmap_plotly actually did.
      rv$provenance$heatmap$row_clustering     <- if (nrow(sub_mat) > 2)
        "hierarchical (ward.D2, euclidean)" else "none (fewer than three rows)"
      # rv$lf_cond_groups is NULL on the DIA and TMT routes, so assigning
      # unconditionally overwrote the condition names those routes had just
      # written with an empty vector. Only assign when there is something to
      # assign; the LF route always populates this reactive, so it is unaffected.
      .hm_lf_conds <- names(rv$lf_cond_groups)
      if (length(.hm_lf_conds) > 0)
        rv$provenance$heatmap$condition_names    <- .hm_lf_conds
    }

    # --- Build title ---
    scale_label <- if (identical(hm_scale_mode, "abundance")) "log2 Abundance" else "Row Z-score"
    avg_tag <- if (avg_mode && !is.null(ctg)) " [Averaged by Condition]" else ""
    title_txt <- paste0("Top ", n_showing, " Genes \u2014 ", mode_label,
                        " (", scale_label, ", Hierarchical Clustering)", avg_tag)

    col_mode <- "cluster"
    if (!is.null(user_col_order) && length(user_col_order) > 0) {
      matched <- intersect(user_col_order, colnames(sub_mat))
      extra   <- setdiff(colnames(sub_mat), user_col_order)
      if (length(matched) > 0) {
        sub_mat <- sub_mat[, c(matched, extra), drop = FALSE]
        col_mode <- "manual"
      }
    }

    # Missing-cell display is set by the route that built the matrix. Only the
    # LF-cache ProDA subroute asks for gaps; every other route takes "zero",
    # which is the pre-existing behaviour.
    .hm_built <- .build_heatmap(sub_mat, title_txt, mode_label,
                                col_order_mode = col_mode,
                                scale_mode = hm_scale_mode,
                                na_display = hm$na_display %||% "zero")
    hm$heatmap_obj <- .hm_built
    hm$render_count <- isolate(hm$render_count) + 1

    # A row count is not proof the clustering succeeded: build_heatmap_plotly
    # catches hclust/dist errors and returns the matrix unordered. Prefer the
    # flag it now reports over the row-count heuristic written above.
    .hm_rc_applied <- .hm_built$row_clustering_applied
    if (length(existing_prov) > 0 && !is.null(.hm_rc_applied)) {
      rv$provenance$heatmap$row_clustering <- if (isTRUE(.hm_rc_applied)) {
        "hierarchical (ward.D2, euclidean)"
      } else if (nrow(sub_mat) > 2) {
        "none (clustering failed)"
      } else {
        "none (fewer than three rows)"
      }
    }

    # Same refinement for columns. A manual order is left reported as "manual":
    # no clustering is attempted on that path, so there is no flag to consult.
    .hm_cc_applied <- .hm_built$col_clustering_applied
    if (length(existing_prov) > 0 && !is.null(.hm_cc_applied) &&
        !identical(.hm_col_mode, "manual")) {
      rv$provenance$heatmap$col_clustering <- if (isTRUE(.hm_cc_applied)) {
        "hierarchical (ward.D2, euclidean)"
      } else if (ncol(sub_mat) > 2) {
        "none (clustering failed)"
      } else {
        "none (fewer than three columns)"
      }
    }

    # Build stats_df from the displayed sub_mat with the correct scale.
    # For z-score modes: compute row z-scores; for abundance: use as-is.
    is_zscore_mode <- !identical(hm_scale_mode, "abundance")
    if (is_zscore_mode && !identical(hm_scale_mode, "zscore_pre")) {
      rm_export <- rowMeans(sub_mat, na.rm = TRUE)
      rs_export <- apply(sub_mat, 1, sd, na.rm = TRUE)
      rs_export[is.na(rs_export) | rs_export == 0] <- NA
      export_mat <- (sub_mat - rm_export) / rs_export
    } else {
      export_mat <- sub_mat
    }
    export_vars <- apply(export_mat, 1, var, na.rm = TRUE)
    export_vars[is.na(export_vars)] <- 0
    displayed_rows <- rownames(sub_mat)
    stats_df_export <- data.frame(
      GeneName  = gene_names[match(displayed_rows, row_labels)],
      Accession = accessions[match(displayed_rows, row_labels)],
      Variance  = round(export_vars, 4),
      stringsAsFactors = FALSE
    )
    export_df <- as.data.frame(export_mat, stringsAsFactors = FALSE)
    colnames(export_df) <- colnames(export_mat)
    hm$stats_df <- cbind(stats_df_export, export_df)

    n_display_cols <- ncol(sub_mat)
    logger(paste0("Heatmap: Built with ", n_showing, " genes x ",
                  n_display_cols, " samples (", mode_label,
                  if (avg_mode && !is.null(ctg)) ", averaged" else "", ")"))
    showNotification(
      paste0("Heatmap generated: ", n_showing, " genes x ",
             n_display_cols, " samples."),
      type = "message", duration = 5)
  }

  # =========================================================================
  #  PIPELINE CACHE: Run Heatmap from pre-built lf_master_cache
  # =========================================================================
  .run_heatmap_from_cache <- function(cache) {
    # Single-trunk: always use cache$processed; protein-level removeBatchEffect
    # applied below when batch correction is active (visualization only).
    processed   <- cache$processed
    acc_to_gene <- cache$acc_to_gene
    pld         <- processed$ProteinLevelData
    grps        <- rv$lf_cond_groups

    withProgress(message = "Heatmap \u2014 From Cache...", value = 0, {
      incProgress(0.2, detail = "Pivoting protein-level data...")

      pld_wide <- tryCatch(
        tidyr::pivot_wider(pld[, c("Protein", "RUN", "LogIntensities")],
                           names_from = "RUN", values_from = "LogIntensities"),
        error = function(e) NULL
      )
      if (is.null(pld_wide)) {
        showNotification("Cache pivot failed  -  try re-summarizing.", type = "error"); return()
      }

      proteins_v <- pld_wide$Protein
      mat <- as.matrix(pld_wide[, -1]); rownames(mat) <- proteins_v
      genes_v <- acc_to_gene[proteins_v]
      n_mapped <- sum(!is.na(genes_v) & genes_v != "")
      genes_v[is.na(genes_v) | genes_v == ""] <- proteins_v[is.na(genes_v) | genes_v == ""]

      logger(paste0("HM [Cache]: ", length(proteins_v), " proteins in PLD, ",
                     n_mapped, " mapped to gene names, ",
                     length(acc_to_gene), " entries in acc_to_gene, ",
                     sum(is.na(mat)), " NAs in matrix (",
                     ncol(mat), " runs)"))

      # Build run-to-label & condition maps from PLD metadata
      run_label_map <- unique(pld[, c("RUN", "GROUP", "SUBJECT")])
      run_to_cond  <- setNames(as.character(run_label_map$GROUP),
                               as.character(run_label_map$RUN))
      run_to_label <- setNames(as.character(run_label_map$SUBJECT),
                               as.character(run_label_map$RUN))
      cn             <- colnames(mat)
      matched_conds  <- run_to_cond[cn]
      matched_conds[is.na(matched_conds)] <- "Unknown"
      col_labels     <- run_to_label[cn]
      col_labels[is.na(col_labels)] <- cn[is.na(col_labels)]
      colnames(mat)  <- col_labels

      incProgress(0.5, detail = "Computing Z-scores...")

      # ── Missing value strategy (user-selected) ────────────────────────────
      missing_mode  <- isolate(input$hm_missing_mode) %||% "omit"
      n_omitted_mv  <- 0L
      n_lod_imputed <- 0L
      if (missing_mode == "omit") {
        complete_rows <- complete.cases(mat)
        n_omitted_mv  <- sum(!complete_rows)
        mat           <- mat[complete_rows, , drop = FALSE]
        proteins_v    <- proteins_v[complete_rows]
        genes_v       <- genes_v[complete_rows]
        if (nrow(mat) < 2) {
          showNotification(
            paste0("Too few complete proteins after omitting rows with missing values (",
                   n_omitted_mv, " omitted). Try LOD imputation instead."),
            type = "error"); return()
        }
      } else if (missing_mode == "lod") {
        for (ci in seq_len(ncol(mat))) {
          na_m <- is.na(mat[, ci])
          if (any(na_m)) {
            lod_v          <- quantile(mat[, ci], 0.01, na.rm = TRUE)
            mat[na_m, ci]  <- lod_v
            n_lod_imputed  <- n_lod_imputed + sum(na_m)
          }
        }
      }
      # ── End missing value handling ────────────────────────────────────────

      # Batch correction for visualization (limma::removeBatchEffect)
      batch_applied <- FALSE
      # Cells fabricated by the pre-batch fill below. Declared here so
      # provenance reports 0 when batch correction does not run at all.
      .n_prebatch_imputed     <- 0L
      .batch_remap_ambiguous  <- FALSE
      # Rows limma could not estimate a batch coefficient for. These pass
      # through removeBatchEffect uncorrected and are displayed alongside
      # corrected rows, so the count is captured and reported rather than
      # left in limma's console warning where the app log never sees it.
      .n_partial_na_rows      <- 0L
      bg <- rv$lf_batch_groups
      if (isTRUE(rv$lf_batch_corrected) && !is.null(bg) && length(bg) >= 2) {
        require_pkg("limma")
        cn_mat <- colnames(mat)   # post-relabel: SUBJECT values

        # Remap batch groups: user assigned original PD column names; ProDA PLD
        # uses bioreplicate (SUBJECT after relabeling) or sample_name (RUN).
        # Use cache$design_df to build original_col → bioreplicate map.
        bg_remapped <- bg
        design_map  <- cache$design_df
        if (!is.null(design_map) &&
            all(c("sample", "bioreplicate") %in% colnames(design_map))) {
          col_to_subj <- setNames(design_map$bioreplicate, design_map$sample)
          # Also try via sample_name in case mat was relabeled to sample_name
          col_to_sname <- setNames(design_map$sample_name, design_map$sample)
          bg_remapped <- lapply(bg, function(runs) {
            via_subj  <- col_to_subj[runs];  via_sname <- col_to_sname[runs]
            # A run that resolves through BOTH maps to DIFFERENT columns is
            # ambiguous: the ifelse below silently prefers the bioreplicate
            # mapping and the loser is never mentioned. Which value wins is
            # deliberately left alone -- changing it would re-batch existing
            # analyses -- but the ambiguity is now reported.
            .amb <- !is.na(via_subj) & !is.na(via_sname) &
                    via_subj %in% cn_mat & via_sname %in% cn_mat &
                    via_subj != via_sname
            if (any(.amb)) {
              .batch_remap_ambiguous <<- TRUE
              for (.k in which(.amb))
                logger(paste0("HM [Cache]: WARNING ambiguous batch remap for run '",
                              runs[.k], "' \u2014 resolves to both '",
                              via_subj[.k], "' (bioreplicate) and '",
                              via_sname[.k], "' (sample_name); using '",
                              via_subj[.k], "'"))
            }
            # Prefer whichever mapping produces hits in cn_mat
            out <- ifelse(!is.na(via_subj)  & via_subj  %in% cn_mat, via_subj,
                   ifelse(!is.na(via_sname) & via_sname %in% cn_mat, via_sname,
                          runs))
            unique(out)
          })
        }

        # Map run IDs (post-relabel) → batch labels
        batch_vec <- rep(NA_character_, length(cn_mat))
        for (bl in names(bg_remapped)) {
          batch_vec[cn_mat %in% bg_remapped[[bl]]] <- bl
        }
        if (!any(is.na(batch_vec)) && length(unique(batch_vec)) >= 2) {
          # Batch correction runs on REAL VALUES ONLY. limma fits each row from
          # its available observations, so no fill is needed beforehand, and
          # filling first was circular: a per-sample median is close to the very
          # quantity a batch correction estimates, so the fabricated cells pulled
          # the estimated coefficients toward the fill values and those
          # coefficients were then applied to every real value. A row with too
          # few observations to estimate a coefficient passes through
          # uncorrected rather than being dropped.
          mat <- withCallingHandlers(
            limma::removeBatchEffect(mat, batch = batch_vec),
            warning = function(w) {
              .wm <- conditionMessage(w)
              logger(paste0("HM [Cache]: limma removeBatchEffect: ", .wm))
              if (grepl("Partial NA coefficient", .wm, fixed = TRUE)) {
                .wn <- suppressWarnings(as.integer(
                  sub("^[^0-9]*([0-9]+).*$", "\\1", .wm)))
                if (!is.na(.wn)) .n_partial_na_rows <<- .wn
              }
              invokeRestart("muffleWarning")
            }
          )
          batch_applied <- TRUE
          if (.n_partial_na_rows > 0L)
            logger(paste0("HM [Cache]: ", .n_partial_na_rows, " protein(s) had ",
                          "too few observations to estimate a batch coefficient ",
                          "and are displayed UNCORRECTED"))
          logger(paste0("HM [Cache]: removeBatchEffect applied (",
                        length(unique(batch_vec)), " batches)"))
          # Fill with the PER-SAMPLE (column) median, now AFTER correction, so
          # these values affect display only. This is a column median, not a row
          # median: median(mat[!na_m, ci]) is taken over the observed values of
          # column ci. Counted at the point the loop actually runs.
          .n_prebatch_imputed <- as.integer(sum(is.na(mat)))
          for (ci in seq_len(ncol(mat))) {
            na_m <- is.na(mat[, ci])
            if (any(na_m)) mat[na_m, ci] <- median(mat[!na_m, ci], na.rm = TRUE)
          }
          if (.n_prebatch_imputed > 0L)
            logger(paste0("HM [Cache]: column-median fill replaced ",
                          .n_prebatch_imputed, " missing cell(s) AFTER ",
                          "removeBatchEffect (display only)"))
        } else {
          logger(paste0("HM [Cache]: batch correction skipped (",
                        sum(is.na(batch_vec)), " runs unmatched to any batch)"))
        }
      }

      n_before_zs <- nrow(mat)
      kept_idx   <- rowwise_zscore(mat, col_labels)$kept_idx
      proteins_v <- proteins_v[kept_idx]
      genes_v    <- genes_v[kept_idx]
      mat_raw    <- mat[kept_idx, , drop = FALSE]
      rownames(mat_raw) <- genes_v

      logger(paste0("HM [Cache]: Variance filter kept ", length(kept_idx), "/",
                     n_before_zs, " proteins (",
                     n_before_zs - length(kept_idx), " dropped)"))

      hm$acc_to_gene  <- acc_to_gene
      hm$col_to_group <- setNames(matched_conds, col_labels)
      # On the ProDA subroute nothing is imputed, so a zero-filled cell in a
      # z-score display is indistinguishable from a protein sitting exactly at
      # its own mean. Render those cells as gaps instead. SAINT and every
      # other route keep the zero fill.
      hm$na_display <- if (identical(isolate(rv$lf_master_cache_meta$de_engine),
                                     "proda")) "gap" else "zero"

      # Provenance
      pipeline <- isolate(rv$lf_pipeline); ip_sub <- isolate(rv$lf_ip_sub)
      pipeline_label <- cache$pipeline_label %||% pipeline
      rv$provenance$heatmap <- list(
        pipeline         = pipeline, ip_sub = ip_sub,
        pipeline_label   = pipeline_label,
        n_proteins       = nrow(mat_raw), n_samples = ncol(mat_raw),
        condition_names  = names(rv$lf_cond_groups),
        batch_corrected  = batch_applied,
        missing_mode     = missing_mode,
        n_omitted        = n_omitted_mv,
        lod_n_imputed    = n_lod_imputed,
        lod_scope        = if (missing_mode == "lod") "per_sample_1st_percentile" else NA_character_,
        prebatch_imputation_method   = "column_median",
        prebatch_n_imputed           = .n_prebatch_imputed,
        prebatch_imputation_position = "after_batch_correction",
        batch_remap_ambiguous      = .batch_remap_ambiguous,
        batch_partial_na_rows      = .n_partial_na_rows,
        ebn_ref_proteins = cache$ebn_ref_vector %||% character(0),
        igg_conditions   = rv$lf_igg_conditions %||% character(0),
        from_cache       = TRUE
      )

      incProgress(0.8, detail = "Rendering heatmap...")
      .finalize_heatmap(mat_raw, genes_v, proteins_v,
                        mode_label = pipeline_label)

      incProgress(1, detail = "Done!")
      logger(paste0("HM [Cache]: Complete. ", nrow(mat_raw), " proteins x ",
                    ncol(mat_raw), " samples."))
    })
  }


  # =========================================================================
  #  RDS SAVE / LOAD
  # =========================================================================
  output$hm_download_csv <- downloadHandler(
    filename = function() {
      mode_safe <- gsub("[^a-zA-Z0-9_]", "_",
                         if (!is.null(hm$heatmap_obj)) hm$heatmap_obj$mode else "heatmap")
      paste0("Heatmap_", mode_safe, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(hm$stats_df)
      write.csv(hm$stats_df, file, row.names = FALSE)
      logger("Heatmap: CSV downloaded")
    }
  )

  # =========================================================================
  #  Column Reorder -- drag-and-drop modal (sortable::rank_list)
  # =========================================================================

  # Label showing current column order status
  output$hm_col_order_label <- renderUI({
    ord <- hm$col_order
    if (is.null(ord) || length(ord) == 0) {
      tags$p(class = "text-muted small", "Default (hierarchical clustering)")
    } else {
      tags$p(class = "small",
             tags$span(class = "badge bg-info", "Custom"),
             " ", paste(ord, collapse = " \u2192 "))
    }
  })

  # Open reorder modal
  observeEvent(input$hm_reorder_cols, {
    req(hm$abundance_mat)
    current_cols <- colnames(hm$abundance_mat)

    # If user already has a custom order, use that; otherwise matrix order
    if (!is.null(hm$col_order) && length(hm$col_order) > 0) {
      matched <- intersect(hm$col_order, current_cols)
      extras  <- setdiff(current_cols, hm$col_order)
      current_cols <- c(matched, extras)
    }

    showModal(modalDialog(
      title = tags$span(icon("arrows-up-down"), " Reorder Heatmap Columns"),
      tagList(
        helpText(class = "small text-muted",
          "Drag columns to set the display order (left to right). ",
          "Click Apply to use this order; the heatmap will re-render."),
        sortable::rank_list(
          text       = "Column Order (drag to reorder):",
          labels     = current_cols,
          input_id   = "hm_col_order_sorted",
          options    = sortable::sortable_options(multiDrag = FALSE),
          class      = "default-sortable"
        ),
        tags$div(class = "mt-2",
          actionButton("hm_col_order_reset", "Reset to Default",
                       class = "btn-outline-secondary btn-sm",
                       icon = icon("rotate-left"))
        )
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("hm_col_order_save", "Apply Order",
                     class = "btn-success", icon = icon("check"))
      )
    ))
  })

  # Save column order
  observeEvent(input$hm_col_order_save, {
    new_order <- input$hm_col_order_sorted
    if (is.null(new_order) || length(new_order) == 0) {
      showNotification("No columns to order.", type = "warning")
      return()
    }
    hm$col_order <- new_order
    removeModal()
    logger(paste0("Heatmap: Column order set to: ",
                  paste(new_order, collapse = " -> ")))

    # Auto-refresh the heatmap with new order
    if (!is.null(hm$abundance_mat)) {
      mode_label <- if (!is.null(hm$heatmap_obj)) hm$heatmap_obj$mode else ""
      .finalize_heatmap(hm$abundance_mat, hm$gene_names, hm$accessions,
                        mode_label = mode_label, update_gene_list = FALSE)
    }
  })

  # Reset column order to default (clustering)
  observeEvent(input$hm_col_order_reset, {
    hm$col_order <- NULL
    logger("Heatmap: Column order reset to default (hierarchical clustering)")
    showNotification("Column order reset to default.", type = "message", duration = 3)
  })

  # =========================================================================
  #  Re-render when gene filter or N changes (without re-running pipeline)
  # =========================================================================
  observeEvent(input$hm_refresh, {
    req(hm$abundance_mat)

    mode_label <- if (!is.null(hm$heatmap_obj)) hm$heatmap_obj$mode else ""
    .finalize_heatmap(hm$abundance_mat, hm$gene_names, hm$accessions,
                      mode_label = mode_label, update_gene_list = FALSE)
  })

  # =========================================================================
  #  Heatmap Scale toggle
  # =========================================================================
  observeEvent(input$hm_scale_mode, {
    req(hm$abundance_mat)
    mode_label <- if (!is.null(hm$heatmap_obj)) hm$heatmap_obj$mode else ""
    .finalize_heatmap(hm$abundance_mat, hm$gene_names, hm$accessions,
                      mode_label = mode_label, update_gene_list = FALSE)
  })

  # =========================================================================
  #  Average by condition toggle
  # =========================================================================
  observeEvent(input$hm_avg_by_cond, {
    req(hm$abundance_mat)

    if (isTRUE(input$hm_avg_by_cond) && is.null(hm$col_to_group)) {
      showNotification(
        "Average by condition is only available for IP/Global pipelines (condition groups required).",
        type = "warning", duration = 6)
      updateCheckboxInput(session, "hm_avg_by_cond", value = FALSE)
      return()
    }

    # Reset manual column order when toggling (columns change)
    hm$col_order <- NULL

    mode_label <- if (!is.null(hm$heatmap_obj)) hm$heatmap_obj$mode else ""
    .finalize_heatmap(hm$abundance_mat, hm$gene_names, hm$accessions,
                      mode_label = mode_label, update_gene_list = FALSE)
  })

  # =========================================================================
  #  Render heatmap
  # =========================================================================
  output$hm_plotly <- renderPlotly({
    req(hm$heatmap_obj)
    hm$render_count  # dependency
    hm$heatmap_obj$plot
  })

  output$hm_plot_container <- renderUI({
    h <- if (!is.null(hm$heatmap_obj)) hm$heatmap_obj$height else 600
    plotlyOutput("hm_plotly", height = paste0(h, "px"))
  })

  # =========================================================================
  #  Invariant Reference Protein Info (displayed below heatmap)
  # =========================================================================
  output$hm_ebn_info <- renderUI({
    refs <- hm$ebn_refs
    if (is.null(refs) || length(refs) == 0) return(NULL)

    # Build display labels: GeneName (Accession) when gene name available
    acc_map <- hm$ebn_acc_to_gene
    labels <- if (!is.null(acc_map)) {
      sapply(refs, function(acc) {
        gn <- acc_map[acc]
        if (!is.na(gn) && nzchar(gn) && gn != acc) {
          paste0(gn, " (", acc, ")")
        } else {
          acc
        }
      }, USE.NAMES = FALSE)
    } else {
      refs
    }

    tags$div(
      class = "alert alert-info py-2 px-3 mt-3 mb-0",
      style = "font-size: 0.88em;",
      tags$strong(icon("dna"), " Invariant Reference Proteins"),
      tags$span(class = "badge bg-primary ms-2", length(refs)),
      tags$p(class = "mb-0 mt-1",
        "The following proteins were used as endogenous normalization standards ",
        "(invariant-protein median normalization via ProDA):"),
      tags$div(class = "mt-1",
        style = "max-height: 120px; overflow-y: auto;",
        tags$small(paste(labels, collapse = " \u2022 "))
      )
    )
  })

  # =========================================================================
  #  Statistics table
  # =========================================================================
  output$hm_stats_table <- DT::renderDT({
    req(hm$stats_df)
    df <- hm$stats_df

    # Identify numeric columns beyond GeneName/Accession for rounding
    numeric_cols <- names(df)[sapply(df, is.numeric)]

    DT::datatable(
      df,
      options = list(pageLength = 25, scrollX = TRUE,
                     order = list(list(2, "desc"))),
      rownames = FALSE,
      filter = "top"
    ) %>% DT::formatRound(columns = numeric_cols, digits = 3)
  })

  # =========================================================================
  #  Status
  # =========================================================================
  output$hm_status <- renderPrint({
    if (is.null(hm$heatmap_obj)) {
      cat("Upload a CSV, define condition groups, and generate the heatmap.\n")
      cat("Pipelines:\n")
      cat("  IP (Immunoprecipitation) - invariant-protein median normalization (ProDA)\n")
      cat("  Global (Whole Tissue)    - equalizeMedians \u2192 TMP \u2192 AFT imputation\n")
      cat("  PD Normalized            - Direct use of PD-normalized abundances\n")
      cat("\nIP and Global read from wizard cache (ProDA + QFeatures summarization).\n")
      cat("Proteins ranked by variance. Row Z-score applied for display.\n")
      cat("Row/column ordering: Hierarchical clustering (Euclidean distance, Ward.D2 linkage)\n")
    } else {
      obj <- hm$heatmap_obj
      # If DE was run via SAINT, annotate the mode label accordingly
      de_prov <- tryCatch(isolate(rv$provenance$de), error = function(e) NULL)
      is_saint_de <- !is.null(de_prov) &&
        identical(de_prov$analysis_mode %||% de_prov$de_engine, "saint")
      mode_display <- if (is_saint_de) {
        bait <- de_prov$bait_protein %||% "unknown"
        paste0("IP (SAINTexpress)  -  Bait: ", bait)
      } else {
        obj$mode
      }
      cat(paste0("Heatmap: ", obj$n_genes, " genes x ", obj$n_samples,
                 " samples\n"))
      cat(paste0("Mode: ", mode_display, "\n"))
      cat(paste0("Row scaling: ",
                 if (identical(isolate(input$hm_scale_mode) %||% "zscore", "abundance"))
                   "none (log2 abundance)" else "Row Z-score", "\n"))
      cat("Row/column ordering: Hierarchical clustering (Euclidean distance, Ward.D2 linkage)\n")
      if (!is.null(hm$ebn_refs) && length(hm$ebn_refs) > 0) {
        cat(paste0("Invariant reference proteins: ", length(hm$ebn_refs), "\n"))
      }
      if (!is.null(hm$stats_df)) {
        cat(paste0("Total proteins: ", nrow(hm$stats_df), "\n"))
      }
    }
  })
}
