# server_trend_clusters.R ======================================================
# Server logic for Proteome Discoverer > Replicated Quantitation >
# Trend Clusters module.
# MSstats-based pipeline for replicated multi-condition data:
#   - User creates named condition groups (e.g., "0h", "6h", "24h")
#   - Assigns protein-level abundance columns to each group (replicates)
#   - Peptide-level data parsed → MSstats long format → dataProcess → imputation
#   - Protein-level summaries averaged per condition → Z-score → Fuzzy c-means
#   - Two pipeline modes:
#       IP:     invariant-protein median normalization → TMP → LOD imputation (0.5× run min)
#       Global: equalizeMedians  → TMP → AFT imputation (MBimpute=TRUE)
# Called from server.R as: server_trend_clusters(input, output, session, rv, logger)
# ==============================================================================

server_trend_clusters <- function(input, output, session, rv, logger) {

  # ---- Module-local reactive values ----
  tcl <- reactiveValues(
    cond_groups       = list(),   # named list: list(name = c(col1, col2, ...))
    protein_df        = NULL,     # raw protein data.frame
    zscore_df         = NULL,     # long-form: Accession, GeneName, Condition, ZScore, Cluster
    cluster_means     = NULL,     # long-form: Cluster, Condition, MeanZ
    cluster_info      = NULL,     # data.frame: Cluster, n_proteins
    n_clusters        = 6L,       # k for clustering
    render_count      = 0L,       # force re-render counter
    highlight_gene    = NULL,     # gene to highlight
    col_labels        = NULL,     # condition labels (group names)
    pipeline          = "ip",     # pipeline mode ("ip" or "global")
    memberships       = NULL,     # full membership matrix (n_proteins x k) for fuzzy c-means
    membership_thresh = 0.5,      # fuzzy c-means membership threshold for display
    zscore_wide       = NULL,     # wide-form Z-score matrix (one col per condition)
    avg_wide          = NULL,     # wide-form averaged abundances (one col per condition)
    cond_order        = NULL,     # custom condition order set by user (NULL = default)
    # Pending state while waiting for IgG modal response
    pending           = NULL      # list(msstats_input, acc_to_gene, cond_names, protein_df, k, pipeline_label)
  )
  register_experiment_reset(rv, tcl)
  methodsTextServer("tcl_methods", provenance = rv$provenance, section = "trend_clusters")

  # =========================================================================
  # SUMMARIZATION CACHE  -" Summarize / Export / Import
  # =========================================================================
  output$tcl_export_summarization <- downloadHandler(
    filename = function() paste0("LF_Summarization_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"),
    content = function(file) {
      cache <- rv$lf_master_cache
      if (is.null(cache)) {
        showNotification("No summarization to export.", type = "warning")
        saveRDS(list(type = "lf_summarization_empty", error = "No cache"), file)
        return()
      }
      saveRDS(list(type = "lf_summarization", version = 1L,
                   cache = cache, meta = rv$lf_master_cache_meta), file)
      logger("TCL: Summarization exported")
    }
  )

  observe({
    if (is.null(rv$lf_master_cache)) {
      shinyjs::disable("tcl_export_summarization")
    } else {
      shinyjs::enable("tcl_export_summarization")
    }
  })

  # =========================================================================
  # Helper: Parse protein-level data from shared CSV
  # =========================================================================
  parse_protein_data <- function() {
    header_line   <- rv$pd_raw_lines[1]
    protein_lines <- rv$pd_raw_lines[rv$pd_protein_line_idx]
    text_block    <- paste(c(header_line, protein_lines), collapse = "\n")
    protein_df <- tryCatch(
      read.csv(text = text_block, header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, na.strings = c("", "NA")),
      error = function(e) { logger(paste0("!!! Protein CSV parse error: ", e$message)); NULL }
    )
    if (is.null(protein_df)) return(NULL)
    if (!isTRUE(rv$pd_has_gn) && !is.null(rv$pd_gene_col) &&
        rv$pd_gene_col %in% colnames(protein_df)) {
      gene_names <- as.character(protein_df[[rv$pd_gene_col]])
      gene_names[is.na(gene_names)] <- ""
    } else {
      descriptions <- protein_df[["Description"]]
      gene_names <- stringr::str_extract(descriptions, "GN=([^ ]+)")
      gene_names <- stringr::str_replace(gene_names, "^GN=", "")
      gene_names[is.na(gene_names)] <- ""
    }
    accessions <- as.character(protein_df[["Accession"]])
    list(protein_df = protein_df, gene_names = gene_names, accessions = accessions)
  }

  # ===========================================================================
  # React to shared CSV data  -" reset state

  # ===========================================================================
  observeEvent(rv$pd_protein_headers, {
    req(rv$pd_protein_headers)
    tcl$zscore_df      <- NULL
    tcl$cluster_means  <- NULL
    tcl$cluster_info   <- NULL
    tcl$protein_df     <- NULL
    tcl$highlight_gene <- NULL
    tcl$col_labels     <- NULL
    logger("TrendClusters: Shared CSV data updated, state reset")
  })

  # ===========================================================================
  # CONDITION ORDER  -" label + reorder modal
  # ===========================================================================

  # Helper: get the default condition list for current data source
  .tcl_default_conds <- function() {
    src <- rv$data_source
    if (identical(src, "tmt")) {
      rv$tmt_conditions
    } else if (identical(src, "dia")) {
      names(rv$dia_cond_groups)
    } else {
      names(rv$lf_cond_groups)
    }
  }

  output$tcl_cond_order_label <- renderUI({
    ord <- tcl$cond_order
    defaults <- .tcl_default_conds()
    if (is.null(ord) || length(ord) == 0) {
      if (is.null(defaults) || length(defaults) == 0) {
        tags$p(class = "text-muted small", "No conditions loaded yet.")
      } else {
        tags$p(class = "small",
               tags$span(class = "badge bg-secondary", "Default"),
               " ", paste(defaults, collapse = " \u2192 "))
      }
    } else {
      tags$p(class = "small",
             tags$span(class = "badge bg-info", "Custom"),
             " ", paste(ord, collapse = " \u2192 "))
    }
  })

  observeEvent(input$tcl_reorder_btn, {
    defaults <- .tcl_default_conds()
    if (is.null(defaults) || length(defaults) == 0) {
      showNotification("No conditions available to reorder. Load data first.", type = "warning")
      return()
    }
    current <- if (!is.null(tcl$cond_order) && length(tcl$cond_order) > 0) {
      matched <- intersect(tcl$cond_order, defaults)
      extras  <- setdiff(defaults, tcl$cond_order)
      c(matched, extras)
    } else {
      defaults
    }
    showModal(modalDialog(
      title = tags$span(icon("arrows-up-down"), " Reorder Conditions"),
      tagList(
        helpText(class = "small text-muted",
          "Drag conditions to set the X-axis order (left to right). ",
          "Click Apply to use this order."),
        sortable::rank_list(
          text       = "Condition Order (drag to reorder):",
          labels     = current,
          input_id   = "tcl_cond_order_sorted",
          options    = sortable::sortable_options(multiDrag = FALSE),
          class      = "default-sortable"
        ),
        tags$div(class = "mt-2",
          actionButton("tcl_cond_order_reset", "Reset to Default",
                       class = "btn-outline-secondary btn-sm",
                       icon = icon("rotate-left")))
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("tcl_cond_order_save", "Apply Order",
                     class = "btn-success", icon = icon("check"))
      )
    ))
  })

  observeEvent(input$tcl_cond_order_save, {
    new_order <- input$tcl_cond_order_sorted
    if (is.null(new_order) || length(new_order) == 0) {
      showNotification("No conditions to order.", type = "warning")
      return()
    }
    tcl$cond_order <- new_order
    removeModal()
    logger(paste0("TrendClusters: Condition order set to: ",
                  paste(new_order, collapse = " -> ")))
  })

  observeEvent(input$tcl_cond_order_reset, {
    tcl$cond_order <- NULL
    removeModal()
    logger("TrendClusters: Condition order reset to default.")
  })

  # ===========================================================================

  # TMT Trend Clusters (unified route  -" called when rv$data_source == "tmt")
  # ===========================================================================
  .run_tcl_tmt <- function() {
    if (is.null(rv$tmt_master_cache)) {
      showNotification("Please run Summarize Data first (Upload Data tab).", type = "warning")
      return()
    }
    default_conds <- rv$tmt_conditions
    cond_order <- if (!is.null(tcl$cond_order) && length(tcl$cond_order) >= 3) {
      intersect(tcl$cond_order, default_conds)
    } else {
      default_conds
    }
    if (is.null(cond_order) || length(cond_order) < 3) {
      showNotification("Trend Clusters requires at least 3 conditions.", type = "warning")
      return()
    }

    k            <- max(2L, min(20L, as.integer(isolate(input$tcl_n_clusters) %||% 6L)))
    tcl$n_clusters      <- k
    tcl$membership_thresh <- isolate(input$tcl_membership) %||% 0.5
    missing_mode <- isolate(input$tcl_missing_mode) %||% "omit"
    tcl$method   <- "fcm"

    logger(paste0("TrendClusters [TMT]: k=", k,
                  ", conditions=", paste(cond_order, collapse = " -> ")))

    withProgress(message = "Running TMT Trend Clusters...", value = 0.1, {
      pld <- rv$tmt_master_cache$ProteinLevelData
      if (is.null(pld) || nrow(pld) == 0) {
        showNotification("No protein data in TMT cache.", type = "error"); return()
      }
      pld$Protein      <- as.character(pld$Protein)
      pld$Condition    <- as.character(pld$Condition)
      # Channels per condition in the FULL annotation, counted before the
      # !is.na(Abundance) filter below. Counting after it would shrink the
      # denominator of the coverage rule and let a thin mean pass.
      .tcl_chan_key <- if ("Run" %in% colnames(pld)) {
        paste(pld$Condition, pld$BioReplicate, pld$Run, sep = "||")
      } else {
        paste(pld$Condition, pld$BioReplicate, sep = "||")
      }
      .tcl_cond_nchan <- tapply(.tcl_chan_key, pld$Condition,
                                function(x) length(unique(x)))
      pld <- pld[pld$Condition %in% cond_order & !is.na(pld$Abundance), ]
      mat_cols <- intersect(cond_order, unique(pld$Condition))
      if (length(mat_cols) < 3) {
        showNotification("Need at least 3 conditions with data.", type = "error"); return()
      }

      # -- Plex batch correction (display only) --------------------------------
      # Trend clusters was the only TMT display consumer that did not remove plex
      # effects; the heatmap, the PCA and the PCA thumbnail all correct
      # themselves. Condition-averaging cancels the offsets under a balanced
      # design but not under an unbalanced one, and this is a general-purpose
      # pipeline. Applied locally rather than in the shared ProteinLevelData
      # cache, because the heatmap and the PCA read that same cache and would
      # then be corrected twice. Mirrors the heatmap's TMT limma branch: same
      # batch factor, same ~ factor(condition) design, same >=2 plex guard, same
      # log-and-proceed on failure. Never used for testing.
      tcl_batch_corrected <- FALSE
      tcl_n_plexes        <- NA_integer_
      tcl_n_bc_not_est    <- 0L
      tcl_n_single_plex   <- 0L
      if ("Run" %in% colnames(pld)) {
        .run_v       <- as.character(pld$Run)
        tcl_n_plexes <- length(unique(.run_v[!is.na(.run_v)]))
        if (tcl_n_plexes >= 2 && requireNamespace("limma", quietly = TRUE)) {
          .ck    <- paste(pld$Condition, pld$BioReplicate, .run_v, sep = "||")
          .do_bc <- TRUE
          if (anyNA(.run_v)) {
            logger(paste0("TrendClusters [TMT]: plex membership could not be resolved ",
                          "for every channel  -  batch correction skipped, ",
                          "proceeding uncorrected."))
            .do_bc <- FALSE
          } else if (anyDuplicated(paste(pld$Protein, .ck)) > 0L) {
            # Condition + BioReplicate + Run must identify a channel uniquely,
            # otherwise the long-to-wide round-trip below would silently collapse
            # two channels into one cell and write the same value back to both.
            # Skip rather than corrupt the matrix.
            logger(paste0("TrendClusters [TMT]: channel identity is not unique on ",
                          "Condition/BioReplicate/Run  -  batch correction skipped, ",
                          "proceeding uncorrected."))
            .do_bc <- FALSE
          } else if (!tmt_plex_correctable(.run_v[match(unique(.ck), .ck)])) {
            # See tmt_plex_correctable(): a plex holding one channel cannot have its
            # offset separated from that channel's own value.
            logger(paste0("TrendClusters [TMT]: at least one plex holds a single ",
                          "channel  -  batch correction skipped, proceeding uncorrected."))
            .do_bc <- FALSE
          }
          if (isTRUE(.do_bc)) {
            .rn   <- unique(pld$Protein)
            .cn   <- unique(.ck)
            .bm   <- matrix(NA_real_, nrow = length(.rn), ncol = length(.cn),
                            dimnames = list(.rn, .cn))
            .bm[cbind(match(pld$Protein, .rn), match(.ck, .cn))] <- pld$Abundance
            .col_cond <- pld$Condition[match(.cn, .ck)]
            .col_run  <- .run_v[match(.cn, .ck)]
            .bm_c <- tryCatch(
              limma::removeBatchEffect(
                .bm,
                batch  = factor(.col_run),
                design = model.matrix(~ factor(.col_cond))
              ),
              error = function(e) {
                logger(paste0("TrendClusters [TMT]: removeBatchEffect failed: ",
                              e$message, "  -  proceeding uncorrected."))
                NULL
              }
            )
            if (!is.null(.bm_c)) {
              # A protein seen in only one plex has no estimable plex
              # coefficient and removeBatchEffect can return NA for its whole
              # row; this route corrects before imputation, so such proteins
              # were being silently lost. Restore them uncorrected.
              .bc_lost <- rowSums(!is.na(.bm_c)) == 0L &
                          rowSums(!is.na(.bm)) > 0L
              tcl_n_bc_not_est <- sum(.bc_lost)
              if (any(.bc_lost)) {
                .bm_c[.bc_lost, ] <- .bm[.bc_lost, ]
                logger(paste0("TrendClusters [TMT]: ", tcl_n_bc_not_est,
                              " proteins were not estimable for the plex term ",
                              "(present in a single plex) and are shown uncorrected."))
              }
              pld$Abundance       <- .bm_c[cbind(match(pld$Protein, .rn),
                                                 match(.ck, .cn))]
              tcl_batch_corrected <- TRUE
              # See .run_hm_tmt: a protein observed in fewer than two plexes has
              # no between-plex contrast, so limma leaves it unadjusted rather
              # than voiding it. Counted so the methods text does not imply the
              # correction reached every protein.
              .tcl_plex_n <- Reduce(`+`, lapply(
                unique(.col_run[!is.na(.col_run)]),
                function(px) as.integer(rowSums(
                  !is.na(.bm[, which(.col_run == px), drop = FALSE])) > 0L)))
              tcl_n_single_plex <- sum(.tcl_plex_n < 2L)
              if (tcl_n_single_plex > 0)
                logger(paste0("TrendClusters [TMT]: ", tcl_n_single_plex,
                              " proteins were observed in a single plex and ",
                              "therefore carry no plex adjustment."))
              logger(paste0("TrendClusters [TMT]: plex batch correction applied (",
                            tcl_n_plexes, " plexes)."))
            }
          }
        }
      }
      incProgress(0.2)

      # Missing-data provenance counters. Initialised on both paths so the TMT
      # trends methods text always gets its missing-data and protein-count
      # sentences; render_tmt_trends_methods drops both when a key is absent.
      tcl_n_input   <- NA_integer_
      tcl_n_omitted <- 0L
      tcl_n_imputed <- 0L
      tcl_lod_value <- NA_real_
      tcl_n_means_suppressed <- 0L
      tcl_biorep_disambig    <- FALSE
      .tcl_excluded <- character(0)
      tcl$excluded_df <- NULL

      # A condition mean is only formed when at least 50% of that condition's
      # channels are observed. Below that the mean is set NA rather than being
      # built from one surviving replicate and treated as equal evidence.
      # Matches the DIA route's .tcl_min_rep_frac. Shared by both missing modes.
      .tcl_min_rep_frac <- 0.5

      if (missing_mode == "lod") {
        # Condition + BioReplicate must identify a channel uniquely, otherwise
        # the long-to-wide reshape below silently collapses two channels into
        # one column and every downstream count, floor and cluster is built on
        # merged data. BioReplicate is commonly numbered within plex in TMT, so
        # plex 1 replicate 1 and plex 2 replicate 1 both yield "C1_1". Same test
        # as the batch-correction block above, but this branch disambiguates
        # instead of skipping, since there is no uncorrected path to fall back
        # to here.
        if (anyDuplicated(paste(pld$Protein, pld$Condition, pld$BioReplicate)) > 0L) {
          if ("Run" %in% colnames(pld) && !anyNA(pld$Run)) {
            pld$BioRep_Label    <- paste0(pld$Condition, "_", pld$BioReplicate,
                                          "_", as.character(pld$Run))
            tcl_biorep_disambig <- TRUE
            logger(paste0("TrendClusters [TMT]: channel identity is not unique on ",
                          "Condition/BioReplicate  -  plex (Run) added to the ",
                          "channel label so every channel keeps its own column."))
          } else {
            logger(paste0("TrendClusters [TMT]: channel identity is not unique on ",
                          "Condition/BioReplicate and Run is unavailable  -  LOD ",
                          "branch aborted rather than building a collapsed matrix."))
            showNotification(
              paste0("Channel identity is ambiguous: BioReplicate repeats within a ",
                     "condition and plex (Run) is not available to separate them. ",
                     "LOD mode cannot build the channel matrix without merging ",
                     "channels. Use Omit-NA mode, or supply a Run column."),
              type = "error", duration = 12)
            return()
          }
        } else {
          pld$BioRep_Label <- paste0(pld$Condition, "_", pld$BioReplicate)
        }
        wide_bio <- reshape(pld[, c("Protein", "BioRep_Label", "Abundance")],
                            idvar = "Protein", timevar = "BioRep_Label",
                            direction = "wide")
        colnames(wide_bio) <- sub("^Abundance\\.", "", colnames(wide_bio))
        all_proteins <- wide_bio$Protein
        bio_cols <- setdiff(colnames(wide_bio), "Protein")
        bio_mat  <- as.matrix(wide_bio[, bio_cols, drop = FALSE])
        rownames(bio_mat) <- all_proteins

        tcl_n_input <- nrow(bio_mat)

        keep_rows <- rep(FALSE, nrow(bio_mat))
        for (cond in mat_cols) {
          cond_cols_bc <- bio_cols[startsWith(bio_cols, paste0(cond, "_"))]
          if (length(cond_cols_bc) >= 2) {
            keep_rows <- keep_rows | (rowSums(!is.na(bio_mat[, cond_cols_bc, drop = FALSE])) >= 2)
          }
        }
        tcl_n_omitted <- sum(!keep_rows)
        bio_mat      <- bio_mat[keep_rows, , drop = FALSE]
        all_proteins <- all_proteins[keep_rows]

        # Proteins with no supported mean in at least one condition are routed
        # out of clustering rather than imputed into it. Every remaining cell
        # would otherwise be filled at the detection floor, and the protein
        # would z-score into a stereotyped step profile that clusters with
        # every other on/off protein - a cluster produced by the substitution,
        # not by biology. Same coverage rule as the omit branch, and the same
        # treatment as the DIA route: reported descriptively, never tested.
        .tcl_cond_nobs <- matrix(0L, nrow = nrow(bio_mat), ncol = length(mat_cols),
                                 dimnames = list(rownames(bio_mat), mat_cols))
        .tcl_cond_reps <- setNames(integer(length(mat_cols)), mat_cols)
        .tcl_supported <- matrix(FALSE, nrow = nrow(bio_mat), ncol = length(mat_cols),
                                 dimnames = list(rownames(bio_mat), mat_cols))
        for (ci in seq_along(mat_cols)) {
          cond_cols_bc <- bio_cols[startsWith(bio_cols, paste0(mat_cols[ci], "_"))]
          .tcl_cond_reps[ci] <- length(cond_cols_bc)
          if (length(cond_cols_bc) > 0) {
            n_ok <- rowSums(!is.na(bio_mat[, cond_cols_bc, drop = FALSE]))
            .tcl_cond_nobs[, ci] <- n_ok
            .tcl_supported[, ci] <- (n_ok / length(cond_cols_bc)) >= .tcl_min_rep_frac
          }
        }
        .tcl_excluded <- rownames(bio_mat)[rowSums(!.tcl_supported) > 0]

        if (length(.tcl_excluded) > 0) {
          .exc_gene <- tmt_cache_gene_labels(rv$tmt_limma_cache, .tcl_excluded)
          .exc_df <- data.frame(Accession = .tcl_excluded, Gene = as.character(.exc_gene),
                                stringsAsFactors = FALSE)
          for (ci in seq_along(mat_cols))
            .exc_df[[paste0("Detected_", mat_cols[ci])]] <-
              paste0(.tcl_cond_nobs[.tcl_excluded, ci], "/", .tcl_cond_reps[ci])
          .exc_df$n_conditions_supported <-
            rowSums(.tcl_supported[.tcl_excluded, , drop = FALSE])
          rownames(.exc_df) <- NULL
          tcl$excluded_df <- .exc_df

          .keep_cl     <- !(rownames(bio_mat) %in% .tcl_excluded)
          bio_mat      <- bio_mat[.keep_cl, , drop = FALSE]
          all_proteins <- all_proteins[.keep_cl]
        }
        logger(paste0("TrendClusters [TMT]: ", length(.tcl_excluded),
                      " proteins routed out of clustering (no supported mean in ",
                      "every condition); ", nrow(bio_mat), " retained."))
        if (nrow(bio_mat) < 2) {
          showNotification(
            paste0("All proteins lack a supported mean in at least one condition ",
                   "(>=", round(.tcl_min_rep_frac * 100), "% of channels observed). ",
                   "Nothing left to cluster."), type = "error", duration = 10)
          return()
        }
        # Per-COLUMN (per-channel) 1st-percentile substitution, mirroring the
        # DIA route and the TMT heatmap. See .run_hm_tmt for the rationale.
        tcl_n_imputed <- sum(is.na(bio_mat))
        lod_vec       <- apply(bio_mat, 2, quantile, probs = 0.01, na.rm = TRUE)
        names(lod_vec) <- colnames(bio_mat)
        if (any(!is.finite(lod_vec)))
          lod_vec[!is.finite(lod_vec)] <- median(lod_vec[is.finite(lod_vec)])
        for (ci in seq_len(ncol(bio_mat))) {
          na_m <- is.na(bio_mat[, ci])
          if (any(na_m)) bio_mat[na_m, ci] <- lod_vec[ci]
        }
        lod_value     <- median(lod_vec)
        tcl_lod_value <- as.numeric(lod_value)
        logger(paste0("TrendClusters [TMT]: per-sample LOD floors  -  min=",
                      round(min(lod_vec), 3), ", median=", round(median(lod_vec), 3),
                      ", max=", round(max(lod_vec), 3), ", max-min=",
                      round(max(lod_vec) - min(lod_vec), 3), " log2 units (ratio=",
                      round(max(lod_vec) / min(lod_vec), 4), ")"))

        .tcl_empty_conds <- character(0)
        cond_avg_list <- lapply(mat_cols, function(cond) {
          cond_cols_bc <- bio_cols[startsWith(bio_cols, paste0(cond, "_"))]
          if (length(cond_cols_bc) == 0) {
            # Filling the whole condition at the detection floor fabricates a
            # timepoint. Return NA and let complete.cases() on the z-score
            # matrix decide instead.
            .tcl_empty_conds <<- c(.tcl_empty_conds, cond)
            logger(paste0("TrendClusters [TMT]: condition '", cond,
                          "' has no replicate columns  -  left missing rather ",
                          "than filled at the detection floor."))
            rep(NA_real_, nrow(bio_mat))
          }
          else rowMeans(bio_mat[, cond_cols_bc, drop = FALSE])
        })
        if (length(.tcl_empty_conds) > 0) {
          # Every protein would now be dropped by the z-score complete.cases,
          # so say why instead of reporting "too few proteins" further down.
          showNotification(
            paste0("No replicate columns for condition(s): ",
                   paste(.tcl_empty_conds, collapse = ", "),
                   ". No protein can be clustered across the selected ",
                   "conditions - deselect these and re-run."),
            type = "error", duration = 12)
          return()
        }
        abund_mat <- do.call(cbind, cond_avg_list)
        colnames(abund_mat) <- mat_cols; rownames(abund_mat) <- all_proteins
        proteins <- all_proteins
      } else {
        avg_df  <- aggregate(Abundance ~ Protein + Condition, data = pld, FUN = mean, na.rm = TRUE)
        # Suppress means built from too few observed channels (see
        # .tcl_min_rep_frac above); the complete.cases() at the end of this
        # branch then drops those proteins, as on the DIA route.
        .nobs_df <- aggregate(Abundance ~ Protein + Condition, data = pld, FUN = length)
        names(.nobs_df)[names(.nobs_df) == "Abundance"] <- "NObs"
        avg_df   <- merge(avg_df, .nobs_df, by = c("Protein", "Condition"),
                          sort = FALSE)
        .cond_tot <- as.numeric(.tcl_cond_nchan[avg_df$Condition])
        .suppress <- !is.na(.cond_tot) & .cond_tot > 0 &
                     (avg_df$NObs / .cond_tot) < .tcl_min_rep_frac
        tcl_n_means_suppressed <- sum(.suppress)
        avg_df$Abundance[.suppress] <- NA_real_
        avg_df   <- avg_df[, c("Protein", "Condition", "Abundance")]
        if (tcl_n_means_suppressed > 0)
          logger(paste0("TrendClusters [TMT]: ", tcl_n_means_suppressed,
                        " condition means suppressed for coverage below ",
                        round(.tcl_min_rep_frac * 100), "% of that condition's channels."))
        wide_df <- reshape(avg_df, idvar = "Protein", timevar = "Condition", direction = "wide")
        colnames(wide_df) <- sub("^Abundance\\.", "", colnames(wide_df))
        proteins  <- wide_df$Protein
        abund_mat <- as.matrix(wide_df[, mat_cols, drop = FALSE])
        rownames(abund_mat) <- proteins
        tcl_n_input <- nrow(abund_mat)
        ok        <- complete.cases(abund_mat)
        tcl_n_omitted <- sum(!ok)
        abund_mat <- abund_mat[ok, , drop = FALSE]; proteins <- proteins[ok]
      }

      row_var <- apply(abund_mat, 1, var, na.rm = TRUE); row_var[is.na(row_var)] <- 0
      abund_mat <- abund_mat[row_var > 0, , drop = FALSE]
      proteins  <- proteins[row_var > 0]
      zscore_mat <- t(apply(abund_mat, 1, function(r) { m <- mean(r); s <- sd(r); (r - m)/s }))
      colnames(zscore_mat) <- mat_cols
      zscore_mat <- zscore_mat[complete.cases(zscore_mat), , drop = FALSE]
      abund_mat  <- abund_mat[rownames(zscore_mat), , drop = FALSE]
      # abund_mat has already been subset to rownames(zscore_mat) on the line
      # above, so match() returned 1:n while `proteins` was still the pre-filter
      # vector: any row dropped by complete.cases() reassigned accessions to the
      # wrong z-score profiles. Row names are the accessions by construction.
      proteins   <- rownames(zscore_mat)

      if (nrow(zscore_mat) < 2) {
        showNotification("Too few proteins after z-scoring.", type = "error"); return()
      }
      if (nrow(zscore_mat) <= k) k <- max(2L, nrow(zscore_mat) - 1L)
      incProgress(0.4)

      N <- nrow(zscore_mat); D <- ncol(zscore_mat)
      m_est <- 1 + (1418/N + 22.05) * D^(-2) +
                (12.33/N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
      # cmeans picks initial centres at random from the data rows. Seed it so a
      # given dataset always yields the same partition and the same cluster IDs.
      .tcl_seed <- 42L
      set.seed(.tcl_seed)
      cl <- tryCatch(
        e1071::cmeans(zscore_mat, centers = k, iter.max = 200, m = m_est,
                      method = "cmeans", dist = "euclidean"),
        error = function(e) {
          showNotification(paste0("Clustering failed: ", e$message), type = "error"); NULL
        }
      )
      if (is.null(cl)) return()
      incProgress(0.2)

      # Same gene map DE and Heatmap use: built at summarization time and
      # restored with saved experiments, so it never depends on the raw files.
      gene_names_v <- tmt_cache_gene_labels(rv$tmt_limma_cache, proteins)
      result <- build_cluster_long_df(zscore_mat, proteins, gene_names_v,
                                       mat_cols, cl$cluster, cl$membership)
      zscore_wide <- as.data.frame(zscore_mat)
      colnames(zscore_wide) <- paste0("ZScore_", mat_cols)
      zscore_wide$Accession <- proteins; zscore_wide$GeneName <- gene_names_v
      zscore_wide$Cluster   <- result$cluster_map[as.character(cl$cluster)]
      zscore_wide$MaxMembership <- apply(cl$membership, 1, max)

      avg_wide <- as.data.frame(abund_mat)
      colnames(avg_wide) <- paste0("Log2_", mat_cols)
      prot_df_v <- data.frame(Accession = proteins, Gene = gene_names_v, stringsAsFactors = FALSE)
      for (ci in seq_along(mat_cols))
        prot_df_v[[paste0("Mean_", mat_cols[ci])]] <- round(abund_mat[, ci], 3)

      tcl$zscore_df     <- result$long_df
      tcl$cluster_means <- result$cluster_means
      tcl$cluster_info  <- result$cluster_info
      tcl$protein_df    <- prot_df_v
      tcl$col_labels    <- mat_cols
      tcl$zscore_wide   <- zscore_wide
      tcl$avg_wide      <- avg_wide
      tcl$memberships   <- result$memberships
      tcl$acc_to_gene   <- setNames(gene_names_v, proteins)
      tcl$pipeline      <- "tmt"
      tcl$render_count  <- isolate(tcl$render_count) + 1L

      rv$provenance$trend_clusters <- list(
        source               = "tmt",
        pipeline             = "tmt",
        pipeline_label       = "TMT quantitative proteomics pipeline",
        missing_mode         = missing_mode,
        n_proteins_input     = tcl_n_input,
        n_omitted            = tcl_n_omitted,
        n_imputed            = tcl_n_imputed,
        lod_value            = tcl_lod_value,
        lod_scope            = if (missing_mode == "lod") "per_sample_1st_percentile" else NA_character_,
        min_replicate_frac   = .tcl_min_rep_frac,
        n_condition_means_suppressed = tcl_n_means_suppressed,
        bioreplicate_disambiguated = tcl_biorep_disambig,
        n_proteins_routed_out = length(.tcl_excluded),
        inclusion_rule       = paste0(
          "a condition mean requires >=", round(.tcl_min_rep_frac * 100),
          "% valid replicates; clustering requires a supported mean in every condition"),
        membership_threshold = tcl$membership_thresh,
        zscore_transform     = TRUE,
        n_proteins_clustered = nrow(zscore_mat),
        n_clusters           = k,
        clustering_method    = "fuzzy c-means (e1071::cmeans)",
        rng_seed             = .tcl_seed,
        fuzzifier            = round(m_est, 3),
        condition_names      = mat_cols,
        batch_corrected      = tcl_batch_corrected,
        n_proteins_batch_not_estimable = tcl_n_bc_not_est,
        n_proteins_single_plex = tcl_n_single_plex,
        n_plexes             = tcl_n_plexes,
        software_used        = "e1071"
      )
      incProgress(0.1)
    })
    logger(paste0("TrendClusters [TMT]: Done  -  ", nrow(tcl$zscore_wide), " proteins, k=", k))
  }

  # ===========================================================================
  # DIA Trend Clusters (unified route  -" called when rv$data_source == "dia")
  # ===========================================================================
  .run_tcl_dia <- function() {
    if (is.null(rv$dia_master_cache)) {
      showNotification("Please complete the DIA Upload Wizard first.", type = "warning")
      return()
    }
    grps <- rv$dia_cond_groups
    if (length(grps) < 3L) {
      showNotification("Trend Clusters requires at least 3 conditions.", type = "warning")
      return()
    }

    k       <- max(2L, min(20L, as.integer(isolate(input$tcl_n_clusters) %||% 6L)))
    tcl$n_clusters      <- k
    tcl$membership_thresh <- isolate(input$tcl_membership) %||% 0.5
    tcl$method          <- "fcm"
    # The missing-value mode control is implemented for the TMT route only.
    # Say so rather than silently ignoring it.
    .tcl_mm <- isolate(input$tcl_missing_mode) %||% "omit"
    if (!identical(.tcl_mm, "omit")) {
      logger(paste0("TrendClusters [DIA]: missing-value mode '", .tcl_mm,
                    "' is not implemented on the DIA route and was ignored; ",
                    "missingness is handled by group-aware filtering."))
      showNotification(
        paste0("Missing-value mode '", .tcl_mm, "' does not apply to DIA data. ",
               "DIA trend clustering uses group-aware filtering instead."),
        type = "warning", duration = 8)
    }
    mat         <- rv$dia_master_cache$protein_matrix
    acc_to_gene <- rv$dia_master_cache$acc_to_gene %||% character(0)

    withProgress(message = "Running DIA Trend Clusters...", value = 0.1, {
      # Build condition means matrix (respect custom order if set)
      default_labels <- names(grps)
      cond_labels <- if (!is.null(tcl$cond_order) && length(tcl$cond_order) >= 3) {
        intersect(tcl$cond_order, default_labels)
      } else {
        default_labels
      }
      if (length(cond_labels) < 3L) cond_labels <- default_labels
      all_prots   <- rownames(mat)
      cond_mat    <- matrix(NA_real_, nrow = length(all_prots), ncol = length(cond_labels),
                            dimnames = list(all_prots, cond_labels))
      # A condition mean is only formed when at least 50% of that condition's
      # replicates are observed. Below that the mean is set NA rather than being
      # built from one surviving replicate and treated as equal evidence.
      .tcl_min_rep_frac <- 0.5
      cond_nobs <- matrix(0L, nrow = length(all_prots), ncol = length(cond_labels),
                          dimnames = list(all_prots, cond_labels))
      for (cond in cond_labels) {
        samps <- intersect(grps[[cond]], colnames(mat))
        if (length(samps) > 0) {
          sub_c <- mat[, samps, drop = FALSE]
          n_ok  <- rowSums(!is.na(sub_c))
          mu_c  <- rowMeans(sub_c, na.rm = TRUE)
          mu_c[(n_ok / length(samps)) < .tcl_min_rep_frac] <- NA_real_
          cond_mat[, cond]  <- mu_c
          cond_nobs[, cond] <- n_ok
        }
      }
      incProgress(0.2)

      n_valid   <- apply(cond_mat, 1, function(r) sum(!is.na(r)))
      cond_mat_f <- cond_mat[n_valid >= 2L, , drop = FALSE]
      if (nrow(cond_mat_f) == 0) {
        showNotification("No proteins with >=2 valid conditions.", type = "error"); return()
      }
      incProgress(0.15)

      zres       <- rowwise_zscore(cond_mat_f, cond_labels)
      zscore_mat <- zres$z
      cond_mat_f <- cond_mat_f[zres$kept_idx, , drop = FALSE]
      if (nrow(zscore_mat) == 0) {
        showNotification("No proteins remain after z-scoring.", type = "error"); return()
      }
      incProgress(0.15)

      # cmeans requires a fully complete matrix  -" drop proteins with residual NAs
      # (rowwise_zscore retains rows with ≥2 non-NA z-scores, but partial NAs
      #  remain for proteins missing in some conditions)
      complete_r <- rowSums(is.na(zscore_mat)) == 0
      if (any(!complete_r)) {
        logger(paste0("TrendClusters [DIA]: dropping ", sum(!complete_r),
                      " proteins with residual NAs before cmeans"))
        zscore_mat <- zscore_mat[complete_r, , drop = FALSE]
        cond_mat_f <- cond_mat_f[complete_r, , drop = FALSE]
      }
      if (nrow(zscore_mat) == 0) {
        showNotification("No complete-case proteins for clustering.", type = "error"); return()
      }

      # Proteins excluded from clustering because they lack a supported mean in
      # every condition. These are not imputed into the clustering: on/off
      # proteins would otherwise be forced to an extreme z-score and pile into
      # one cluster purely from the substitution. Reported descriptively instead.
      .tcl_excluded <- setdiff(rownames(cond_mat), rownames(zscore_mat))
      if (length(.tcl_excluded) > 0) {
        .exc_nobs <- cond_nobs[.tcl_excluded, , drop = FALSE]
        .exc_reps <- vapply(cond_labels, function(cd)
          length(intersect(grps[[cd]], colnames(mat))), integer(1))
        .exc_gene <- acc_to_gene[.tcl_excluded]
        .exc_gene[is.na(.exc_gene)] <- .tcl_excluded[is.na(.exc_gene)]
        .exc_df <- data.frame(Accession = .tcl_excluded, Gene = .exc_gene,
                              stringsAsFactors = FALSE)
        for (ci in seq_along(cond_labels))
          .exc_df[[paste0("Detected_", cond_labels[ci])]] <-
            paste0(.exc_nobs[, ci], "/", .exc_reps[ci])
        .exc_df$n_conditions_supported <-
          rowSums(!is.na(cond_mat[.tcl_excluded, , drop = FALSE]))
        rownames(.exc_df) <- NULL
        tcl$excluded_df <- .exc_df
      } else {
        tcl$excluded_df <- NULL
      }
      logger(paste0("TrendClusters [DIA]: ", length(.tcl_excluded),
                    " proteins routed out of clustering (no supported mean in ",
                    "every condition); ", nrow(zscore_mat), " clustered."))

      # Optional variability prefilter on the condition-mean matrix (log2 units).
      # Defaults to 0 (no filtering) until a UI control supplies a threshold.
      .tcl_min_sd <- suppressWarnings(as.numeric(isolate(input$tcl_min_sd) %||% 0))
      if (!is.finite(.tcl_min_sd)) .tcl_min_sd <- 0
      .tcl_n_lowsd <- 0L
      if (.tcl_min_sd > 0) {
        .sd_row <- apply(cond_mat_f, 1, sd, na.rm = TRUE)
        .sd_keep <- is.finite(.sd_row) & .sd_row >= .tcl_min_sd
        .tcl_n_lowsd <- sum(!.sd_keep)
        zscore_mat <- zscore_mat[.sd_keep, , drop = FALSE]
        cond_mat_f <- cond_mat_f[.sd_keep, , drop = FALSE]
        logger(paste0("TrendClusters [DIA]: SD prefilter (>= ", .tcl_min_sd,
                      " log2) removed ", .tcl_n_lowsd, " proteins."))
        if (nrow(zscore_mat) < 2) {
          showNotification("SD prefilter left fewer than 2 proteins.", type = "error")
          return()
        }
      }

      if (nrow(zscore_mat) <= k) k <- max(2L, nrow(zscore_mat) - 1L)
      # Schwaemmle & Jensen fuzzifier estimator, matching the TMT, LF-cache and
      # Single-Shot routes. Replaces a route-specific ad hoc formula.
      N     <- nrow(zscore_mat); D <- ncol(zscore_mat)
      m_est <- 1 + (1418/N + 22.05) * D^(-2) +
        (12.33/N + 0.243) * D^(-0.0406 * log(N) - 0.1134)

      # cmeans picks initial centres at random from the data rows. Seed it so a
      # given dataset always yields the same partition and the same cluster IDs.
      .tcl_seed <- 42L
      set.seed(.tcl_seed)
      cl <- tryCatch(
        e1071::cmeans(zscore_mat, centers = k, iter.max = 200,
                      m = m_est, method = "cmeans", dist = "euclidean"),
        error = function(e) {
          showNotification(paste0("Clustering failed: ", e$message), type = "error"); NULL
        }
      )
      if (is.null(cl)) return()
      incProgress(0.2)

      # Xie-Beni validity scan. Reported only: k is not changed.
      .tcl_xb_fun <- function(x, cen, u, m) {
        d2 <- outer(rowSums(x^2), rowSums(cen^2), "+") - 2 * tcrossprod(x, cen)
        d2[d2 < 0] <- 0
        cd <- as.matrix(dist(cen))^2
        diag(cd) <- Inf
        sum((u^m) * d2) / (nrow(x) * min(cd))
      }
      .tcl_k_scan <- NULL
      .tcl_k_best <- NA_integer_
      if (nrow(zscore_mat) <= 20000L) {
        setProgress(detail = "Scanning k (Xie-Beni)...")
        .k_range <- 2:10
        .xb <- vapply(.k_range, function(kk) {
          set.seed(.tcl_seed)
          fit <- tryCatch(
            e1071::cmeans(zscore_mat, centers = kk, iter.max = 100,
                          m = m_est, method = "cmeans", dist = "euclidean"),
            error = function(e) NULL)
          if (is.null(fit)) NA_real_
          else .tcl_xb_fun(zscore_mat, fit$centers, fit$membership, m_est)
        }, numeric(1))
        .tcl_k_scan <- setNames(round(.xb, 4), paste0("k", .k_range))
        if (any(is.finite(.xb))) .tcl_k_best <- .k_range[which.min(.xb)]
        logger(paste0("TrendClusters [DIA]: Xie-Beni  -  ",
                      paste(paste0(names(.tcl_k_scan), "=", .tcl_k_scan),
                            collapse = ", "),
                      "; minimum at k=", .tcl_k_best,
                      " (k used = ", k, ", unchanged)"))
      }

      accessions    <- rownames(zscore_mat)
      gene_names_v  <- acc_to_gene[accessions]
      gene_names_v[is.na(gene_names_v)] <- accessions[is.na(gene_names_v)]

      result <- build_cluster_long_df(zscore_mat, accessions, gene_names_v,
                                       cond_labels, cl$cluster, cl$membership)
      zscore_wide <- as.data.frame(zscore_mat)
      colnames(zscore_wide) <- paste0("ZScore_", cond_labels)
      zscore_wide$Accession     <- accessions
      zscore_wide$GeneName      <- gene_names_v
      zscore_wide$Cluster       <- result$cluster_map[as.character(cl$cluster)]
      zscore_wide$MaxMembership <- apply(cl$membership, 1, max)

      avg_wide <- as.data.frame(cond_mat_f)
      colnames(avg_wide) <- paste0("Log2_", cond_labels)
      prot_df_v <- data.frame(Accession = accessions, Gene = gene_names_v,
                               stringsAsFactors = FALSE)
      for (ci in seq_along(cond_labels))
        prot_df_v[[paste0("Mean_", cond_labels[ci])]] <- round(cond_mat_f[, ci], 3)

      tcl$zscore_df     <- result$long_df
      tcl$cluster_means <- result$cluster_means
      tcl$cluster_info  <- result$cluster_info
      tcl$protein_df    <- prot_df_v
      tcl$col_labels    <- cond_labels
      tcl$zscore_wide   <- zscore_wide
      tcl$avg_wide      <- avg_wide
      tcl$memberships   <- result$memberships
      tcl$acc_to_gene   <- setNames(gene_names_v, accessions)
      tcl$pipeline      <- "dia"
      tcl$render_count  <- isolate(tcl$render_count) + 1L

      # Membership diagnostics. Recorded only: nothing is filtered or refitted
      # on the basis of these numbers.
      .tcl_maxmem   <- apply(cl$membership, 1, max)
      .tcl_mem_med  <- round(median(.tcl_maxmem), 3)
      .tcl_mem_frac <- round(mean(.tcl_maxmem >= 0.5), 3)
      logger(paste0("TrendClusters [DIA]: median max-membership = ", .tcl_mem_med,
                    "; fraction >= 0.5 = ", .tcl_mem_frac,
                    "; fuzzifier = ", round(m_est, 3)))

      rv$provenance$trend_clusters <- list(
        source               = "dia",
        n_proteins_input     = nrow(cond_mat),
        n_proteins_clustered = nrow(zscore_mat),
        min_replicate_frac   = .tcl_min_rep_frac,
        sd_prefilter         = .tcl_min_sd,
        n_removed_low_sd     = .tcl_n_lowsd,
        n_proteins_routed_out = length(.tcl_excluded),
        inclusion_rule       = paste0(
          "a condition mean requires >=", round(.tcl_min_rep_frac * 100),
          "% valid replicates; clustering requires a supported mean in every condition"),
        n_clusters           = k,
        clustering_method    = "fuzzy c-means (e1071::cmeans)",
        rng_seed             = .tcl_seed,
        fuzzifier            = round(m_est, 3),
        zscore_transform     = TRUE,
        k_validity_index     = "Xie-Beni",
        k_validity_scan      = .tcl_k_scan,
        k_validity_min_at    = .tcl_k_best,
        membership_threshold  = tcl$membership_thresh,
        membership_max_median = .tcl_mem_med,
        membership_frac_above_half = .tcl_mem_frac,
        condition_names      = cond_labels,
        software_used        = "e1071"
      )
      incProgress(0.1)
      logger(paste0("TrendClusters [DIA]: Done  -  ", nrow(zscore_mat),
                    " proteins, k=", k))
    })
  }

  # ===========================================================================
  # GENERATE  -" MSstats peptide-level pipeline → protein summaries → clustering
  # ===========================================================================
  observeEvent(input$tcl_plot_btn, {
    # Route to TMT or DIA path if appropriate data source is active
    src <- rv$data_source
    if (identical(src, "tmt")) {
      .run_tcl_tmt()
      return()
    } else if (identical(src, "dia")) {
      .run_tcl_dia()
      return()
    }

    # LF DDA path — requires wizard cache
    if (is.null(rv$lf_master_cache)) {
      showNotification(
        "Please complete the Upload Wizard (Summarize Data step) before running Trend Clusters.",
        type = "warning", duration = 8)
      return()
    }

    grps <- rv$lf_cond_groups
    if (length(grps) < 2) {
      showNotification("Define at least 2 condition groups first.", type = "warning")
      return()
    }
    if (length(grps) < 3) {
      showNotification(
        paste0(
          "Trend Clusters requires at least 3 conditions to detect meaningful ",
          "expression trajectories. With only 2 groups, use the DE module instead."
        ),
        type = "warning", duration = 10)
      return()
    }

    k <- as.integer(if (!is.null(input$tcl_n_clusters)) input$tcl_n_clusters else 6L)
    if (k < 2) k <- 2L
    if (k > 20) k <- 20L
    tcl$n_clusters <- k

    # Require e1071 for fuzzy c-means (no K-means fallback)
    if (!requireNamespace("e1071", quietly = TRUE)) {
      showNotification(
        "e1071 package is required for fuzzy c-means clustering. Install with: install.packages('e1071')",
        type = "error", duration = 10)
      return()
    }

    pipeline <- if (!is.null(rv$lf_pipeline)) rv$lf_pipeline else "ip"
    tcl$pipeline <- pipeline
    tcl$membership_thresh <- if (!is.null(input$tcl_membership))
      input$tcl_membership else 0.5

    # Reset viewing state
    tcl$highlight_gene <- NULL
    updateTextInput(session, "tcl_gene_search", value = "")

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
            "running Trend Clusters. This may take a few minutes."
          ),
          tags$p(class = "text-muted small mb-0",
            "Pipeline: ", tags$b(pipeline),
            if (identical(pipeline, "ip")) paste0(" (", isolate(rv$lf_ip_sub), ")") else "")
        ),
        size = "m", easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("tcl_confirm_resummarize", "Continue",
                       class = "btn-primary", icon = icon("play"))
        )
      ))
      return()
    }

    # Use cache if available and matches current pipeline
    if (!is.null(cache) && !is.null(meta) &&
        identical(meta$pipeline, pipeline) &&
        (pipeline != "ip" || identical(meta$ip_sub, isolate(rv$lf_ip_sub)))) {
      logger(paste0("TCL [", pipeline, "]: Using cached summarization"))
      .run_tcl_from_cache(cache, k)
      return()
    }

    # No matching cache — direct user to wizard
    showNotification(
      "Please complete the Upload Wizard (Summarize Data step) before running Trend Clusters.",
      type = "warning", duration = 8)
  })

  observeEvent(input$tcl_confirm_resummarize, {
    removeModal()
    rv$lf_cache_stale <- FALSE
    # Re-trigger the run
    shinyjs::click("tcl_plot_btn")
  })

  # ===========================================================================
  # PIPELINE CACHE: Run Trend Clusters from pre-built lf_master_cache
  # ===========================================================================
  .run_tcl_from_cache <- function(cache, k) {
    processed   <- cache$processed
    acc_to_gene <- cache$acc_to_gene
    pld         <- processed$ProteinLevelData
    grps        <- rv$lf_cond_groups
    cond_names  <- names(grps)
    n_conds     <- length(cond_names)

    if (n_conds < 3) {
      showNotification(
        paste0(
          "Trend Clusters requires at least 3 conditions to detect meaningful ",
          "expression trajectories. With only 2 groups, use the DE module instead."
        ),
        type = "warning", duration = 10)
      return()
    }

    pipeline    <- isolate(rv$lf_pipeline)
    pipeline_label <- cache$pipeline_label %||% pipeline

    withProgress(message = paste0("Trend Clusters \u2014 From Cache (k=", k, ")"), value = 0, {
      incProgress(0.15, detail = "Pivoting protein-level data...")

      pld_wide <- tryCatch(
        tidyr::pivot_wider(pld[, c("Protein", "RUN", "LogIntensities")],
                           names_from = "RUN", values_from = "LogIntensities"),
        error = function(e) NULL
      )
      if (is.null(pld_wide)) {
        showNotification("Cache pivot failed  -  try re-summarizing.", type = "error"); return()
      }

      proteins_v <- as.character(pld_wide$Protein)
      mat        <- as.matrix(pld_wide[, -1]); rownames(mat) <- proteins_v

      incProgress(0.30, detail = "Computing condition means...")
      # Build condition means matrix (one row per protein, one col per condition)
      # Pivot by RUN gives a true protein × sample matrix; map RUN → GROUP for aggregation
      run_to_cond_map <- setNames(as.character(pld$GROUP), as.character(pld$RUN))
      run_to_cond_map <- run_to_cond_map[!duplicated(names(run_to_cond_map))]
      cn_mat          <- colnames(mat)
      cond_of_col     <- run_to_cond_map[cn_mat]
      cond_of_col[is.na(cond_of_col)] <- "Unknown"

      cond_mean_list <- lapply(cond_names, function(cname) {
        cols <- cn_mat[cond_of_col == cname]
        if (length(cols) == 0) return(rep(NA_real_, nrow(mat)))
        rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
      })
      cond_mat <- do.call(cbind, cond_mean_list)
      colnames(cond_mat) <- cond_names; rownames(cond_mat) <- proteins_v

      # Filter rows with all-NA or insufficient values
      valid_r <- rowSums(!is.na(cond_mat) & is.finite(cond_mat)) >= 2
      cond_mat_f <- cond_mat[valid_r, , drop = FALSE]
      if (nrow(cond_mat_f) < 5) {
        showNotification("Too few proteins with valid values across conditions.", type = "error")
        return()
      }

      incProgress(0.45, detail = "Z-score normalization...")
      zs_result  <- rowwise_zscore(cond_mat_f, cond_names)
      zscore_mat <- zs_result$z
      kept_idx   <- zs_result$kept_idx
      prot_accessions <- rownames(cond_mat_f)[kept_idx]
      gene_names_v    <- acc_to_gene[prot_accessions]
      gene_names_v[is.na(gene_names_v) | gene_names_v == ""] <-
        prot_accessions[is.na(gene_names_v) | gene_names_v == ""]
      cond_mat_f <- cond_mat_f[kept_idx, , drop = FALSE]

      if (nrow(zscore_mat) < 2) {
        showNotification("Fewer than 2 proteins with non-zero variance  -  cannot cluster.",
                         type = "error"); return()
      }
      if (nrow(zscore_mat) <= k) {
        k <- max(2L, nrow(zscore_mat) - 1L)
        showNotification(
          paste0("Too few proteins after variance filter. Auto-reduced to k=", k, "."),
          type = "warning", duration = 6)
      }

      incProgress(0.55, detail = "Fuzzy c-means clustering...")
      rownames(zscore_mat) <- prot_accessions

      # cmeans requires a fully complete matrix  -" drop any rows with residual NAs
      # (rowwise_zscore retains rows with ≥2 non-NA z-scores, but partial NAs
      #  remain for proteins missing in some conditions)
      complete_r      <- rowSums(is.na(zscore_mat)) == 0
      n_incomplete    <- sum(!complete_r)
      if (n_incomplete > 0) {
        logger(paste0("TrendClusters [Cache]: dropping ", n_incomplete,
                      " proteins with residual NAs before cmeans"))
        zscore_mat      <- zscore_mat[complete_r, , drop = FALSE]
        prot_accessions <- prot_accessions[complete_r]
        gene_names_v    <- gene_names_v[complete_r]
        cond_mat_f      <- cond_mat_f[complete_r, , drop = FALSE]
      }

      if (nrow(zscore_mat) < 2) {
        showNotification("Fewer than 2 complete-case proteins for clustering.", type = "error")
        return()
      }
      if (nrow(zscore_mat) <= k) {
        k <- max(2L, nrow(zscore_mat) - 1L)
        showNotification(paste0("Auto-reduced k=", k, " (too few complete proteins)."),
                         type = "warning", duration = 6)
      }

      N <- nrow(zscore_mat); D <- ncol(zscore_mat)
      m_est <- 1 + (1418/N + 22.05) * D^(-2) +
        (12.33/N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
      cl <- tryCatch({
        e1071::cmeans(zscore_mat, centers = k, m = m_est,
                      method = "cmeans", dist = "euclidean")
      }, error = function(e) {
        logger(paste0("TrendClusters [Cache]: cmeans error: ", e$message))
        showNotification(paste("Clustering error:", e$message), type = "error")
        NULL
      })
      if (is.null(cl)) return()
      clusters    <- cl$cluster
      memberships <- cl$membership

      incProgress(0.70, detail = "Sorting clusters...")
      result <- build_cluster_long_df(zscore_mat, prot_accessions, gene_names_v,
                                       cond_names, clusters, memberships)
      long_df       <- result$long_df
      cluster_means <- result$cluster_means
      cluster_info  <- result$cluster_info
      cluster_map   <- result$cluster_map
      memberships   <- result$memberships

      if (is.null(long_df) || nrow(long_df) == 0) {
        showNotification("Clustering produced no results.", type = "error"); return()
      }

      incProgress(0.80, detail = "Building output tables...")

      # Wide-form Z-score for CSV export
      zscore_wide <- as.data.frame(zscore_mat)
      colnames(zscore_wide) <- paste0("ZScore_", cond_names)
      zscore_wide$Accession <- prot_accessions
      zscore_wide$GeneName  <- gene_names_v
      zscore_wide$Cluster   <- clusters
      zscore_wide$Cluster   <- cluster_map[as.character(zscore_wide$Cluster)]
      zscore_wide$MaxMembership <- apply(memberships, 1, max)

      avg_wide <- as.data.frame(cond_mat_f)
      colnames(avg_wide) <- paste0("Log2_", cond_names)

      # Mean intensities per condition per protein
      prot_df_v <- data.frame(
        Accession = prot_accessions,
        Gene      = gene_names_v,
        stringsAsFactors = FALSE
      )
      for (ci in seq_along(cond_names)) {
        prot_df_v[[paste0("Mean_", cond_names[ci])]] <- round(cond_mat_f[, ci], 3)
      }

      tcl$zscore_df     <- long_df
      tcl$cluster_means <- cluster_means
      tcl$cluster_info  <- cluster_info
      tcl$protein_df    <- prot_df_v
      tcl$col_labels    <- cond_names
      tcl$zscore_wide   <- zscore_wide
      tcl$avg_wide      <- avg_wide
      tcl$memberships   <- memberships
      tcl$acc_to_gene   <- acc_to_gene

      # Provenance
      rv$provenance$trend_clusters <- list(
        pipeline           = pipeline, ip_sub = isolate(rv$lf_ip_sub),
        pipeline_label     = pipeline_label,
        n_clusters         = k,
        n_proteins         = nrow(zscore_mat),
        n_proteins_input   = nrow(cond_mat),
        n_conditions       = n_conds,
        condition_names    = cond_names,
        zscore_transform   = TRUE,
        ebn_ref_proteins   = cache$ebn_ref_vector %||% character(0),
        lod_n_imputed      = cache$lod_n_imputed %||% 0L,
        igg_conditions     = rv$lf_igg_conditions %||% character(0),
        membership_thresh  = tcl$membership_thresh,
        clustering_method  = "fuzzy c-means (e1071::cmeans)",
        fuzzifier          = round(m_est, 3),
        software_used      = "e1071",
        from_cache         = TRUE
      )

      incProgress(1, detail = "Done!")
      logger(paste0("TCL [Cache]: ", nrow(zscore_mat), " proteins in ", k,
                    " clusters across ", n_conds, " conditions"))
    })
  }


  # ===========================================================================
  # Highlight control
  # ===========================================================================
  observeEvent(input$tcl_gene_search_btn, {
    gene <- trimws(input$tcl_gene_search)
    tcl$highlight_gene <- if (nchar(gene) == 0) NULL else gene
    tcl$render_count <- isolate(tcl$render_count) + 1L
  })

  # ===========================================================================
  # POST-ANALYSIS: Membership threshold live filter
  # When the user moves the slider AFTER analysis, re-render with filter.
  # ===========================================================================
  observeEvent(input$tcl_membership, {
    tcl$membership_thresh <- input$tcl_membership
    # The slider re-renders without re-clustering, so the recorded threshold
    # would otherwise go stale. Applies to whichever route produced the plot
    # currently on screen; the diagnostics themselves are not recomputed.
    if (!is.null(rv$provenance$trend_clusters))
      rv$provenance$trend_clusters$membership_threshold <- input$tcl_membership
    # Re-render only if analysis already ran
    if (!is.null(tcl$zscore_df)) {
      tcl$render_count <- isolate(tcl$render_count) + 1L
    }
  }, ignoreInit = TRUE)

  # ===========================================================================
  # Profile Plot  -" renderPlotly (delegates to shared render_trend_subplots)
  # ===========================================================================
  output$tcl_profile_plotly <- renderPlotly({
    req(tcl$zscore_df, tcl$cluster_means)
    tcl$render_count   # force-refresh dependency

    pipeline_label <- switch(
      tcl$pipeline %||% "global",
      ip     = if (identical(isolate(rv$lf_ip_sub), "invariant")) "Invariant-protein median (ProDA)" else "IP (SAINTexpress input, no normalization)",
      dia    = "DIA (log2, median-centered)",
      tmt    = "TMT",
      global = "Global (eqMedians + AFT)",
      tcl$pipeline
    )

    long_df       <- tcl$zscore_df
    cluster_means <- tcl$cluster_means
    col_labels    <- tcl$col_labels

    # Apply custom condition order if the user has set one
    cond_order <- tcl$cond_order
    if (!is.null(cond_order) && length(cond_order) > 0 && !is.null(col_labels)) {
      new_order <- intersect(cond_order, col_labels)
      if (length(new_order) == length(col_labels)) {
        old_positions <- match(new_order, col_labels)
        old_to_new    <- setNames(seq_along(new_order), old_positions)
        long_df$CondIdx       <- old_to_new[as.character(long_df$CondIdx)]
        cluster_means$CondIdx <- old_to_new[as.character(cluster_means$CondIdx)]
        long_df$Condition       <- factor(as.character(long_df$Condition),       levels = new_order)
        cluster_means$Condition <- factor(as.character(cluster_means$Condition), levels = new_order)
        col_labels <- new_order
      }
    }

    render_trend_subplots(
      long_df       = long_df,
      cluster_means = cluster_means,
      cluster_info  = tcl$cluster_info,
      col_labels    = col_labels,
      k             = tcl$n_clusters,
      highlight     = tcl$highlight_gene,
      exact_match   = isTRUE(input$tcl_exact_match),
      title         = paste0("Replicated Trend Profiles \u2014 ", pipeline_label,
                              " (k = ", tcl$n_clusters, ")"),
      method        = "fcm",
      mem_thresh    = tcl$membership_thresh
    )
  })

  # ===========================================================================
  # Cluster membership table  -" renderDT
  # ===========================================================================
  output$tcl_cluster_table <- renderDT({
    req(tcl$zscore_df)
    long_df <- tcl$zscore_df
    prot_df <- long_df[!duplicated(long_df$Accession),
                       c("Accession", "GeneName", "Cluster", "Membership")]
    prot_df <- prot_df[order(prot_df$Cluster, prot_df$GeneName), ]
    rownames(prot_df) <- NULL
    prot_df$Membership <- round(prot_df$Membership, 3)
    colnames(prot_df) <- c("Accession", "Gene", "Cluster", "Membership")

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
  # Proteins routed out of clustering (insufficient condition coverage)
  #
  # PARITY NOTE (retained deliberately — not dead code to remove):
  # This table has no uiOutput slot in ui.R, so it is never displayed. It is kept
  # because the PTM Trends tab exposes the equivalent view (ptm_trends_excluded_table,
  # ui.R ~3026) together with a download button (ptm_trends_download_excluded,
  # ui.R ~1598). Wiring this one up is a candidate for OSMSS v4 so the label-free
  # Trend Clusters tab reaches parity with PTM Trends.
  # ===========================================================================
  output$tcl_excluded_table <- renderDT({
    req(tcl$excluded_df)
    datatable(tcl$excluded_df,
      options = list(pageLength = 25, scrollY = "400px", dom = "frtip"),
      filter = "top", rownames = FALSE, class = "compact stripe hover")
  })

  # ===========================================================================
  # CSV Download
  # ===========================================================================
  output$tcl_download_csv <- downloadHandler(
    filename = function() {
      paste0("TrendClusters_k", tcl$n_clusters, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(tcl$zscore_wide, tcl$protein_df, tcl$avg_wide)

      zw  <- tcl$zscore_wide
      pdf <- tcl$protein_df
      aw  <- tcl$avg_wide
      grps <- rv$lf_cond_groups

      # Build export: Cluster | Accession | Gene | Membership | Log2 per condition | Z-score per condition
      export_df <- data.frame(
        Cluster   = zw$Cluster,
        Accession = zw$Accession,
        Gene      = zw$GeneName,
        stringsAsFactors = FALSE
      )

      # Membership score (always present  -" fuzzy c-means only)
      if ("MaxMembership" %in% colnames(zw)) {
        export_df$MaxMembership <- round(zw$MaxMembership, 4)
      }

      # Log2 intensity per condition (MSstats-summarized)
      for (col in colnames(aw)) {
        export_df[[col]] <- aw[[col]]
      }

      # Raw abundance columns from protein data (match by accession)
      all_rep_cols <- unlist(grps, use.names = FALSE)
      if ("Accession" %in% colnames(pdf)) {
        row_idx <- match(zw$Accession, pdf[["Accession"]])
        for (col in all_rep_cols) {
          if (col %in% colnames(pdf)) {
            export_df[[col]] <- as.numeric(pdf[[col]][row_idx])
          }
        }
      }

      # Z-score columns
      zcols <- grep("^ZScore_", colnames(zw), value = TRUE)
      for (col in zcols) {
        export_df[[col]] <- zw[[col]]
      }

      # Sort by cluster then gene
      export_df <- export_df[order(export_df$Cluster, export_df$Gene), ]

      write.csv(export_df, file, row.names = FALSE)
      logger(paste0("TrendClusters: CSV exported  -  ", nrow(export_df), " proteins"))
    }
  )

  # ===========================================================================
  # Status summary
  # ===========================================================================
  output$tcl_status <- renderPrint({
    if (is.null(tcl$zscore_df)) {
      cat("Load data, define condition groups, and generate the analysis to see results.\n")
    } else {
      info <- tcl$cluster_info
      info <- info[order(info$Cluster), ]
      n_total <- sum(info$n_proteins)
      pipeline_label <- switch(
        tcl$pipeline %||% "global",
        ip     = if (identical(isolate(rv$lf_ip_sub), "invariant")) "Invariant-protein median (ProDA)" else "IP (SAINTexpress input, no normalization)",
        dia    = "DIA (log2, median-centered)",
        tmt    = "TMT",
        global = "Global (eqMedians + AFT)",
        tcl$pipeline
      )
      cat(paste0("Pipeline: ", pipeline_label, "\n"))
      cat(paste0("Clustering: Fuzzy c-means\n"))
      cat(paste0("Total proteins clustered: ", n_total, "\n"))
      cat(paste0("Number of clusters (k): ", tcl$n_clusters, "\n"))
      cat(paste0("Conditions: ", length(tcl$col_labels), "\n"))
      cat(paste0("  (", paste(tcl$col_labels, collapse = " -> "), ")\n\n"))

      # Show group composition
      grps <- rv$lf_cond_groups
      cat("Condition groups:\n")
      for (i in seq_along(grps)) {
        cat(paste0("  ", names(grps)[i], ": ", length(grps[[i]]), " replicate(s)\n"))
      }
      cat("\n")

      for (i in seq_len(nrow(info))) {
        cat(paste0("  Cluster ", info$Cluster[i], ": ",
                    info$n_proteins[i], " proteins\n"))
      }
      hl <- tcl$highlight_gene
      if (!is.null(hl) && nchar(hl) > 0) {
        long_df <- tcl$zscore_df
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
