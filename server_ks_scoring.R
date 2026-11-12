# server_ks_scoring.R
# ==============================================================================
# Kinase-Substrate (KS) Scoring — Shiny Module (OSMSS v4)
#
# Implements the PhosR kinaseSubstrateScore() algorithm in base R (no PhosR
# package required).  The module accepts a PhosR-format CSV upload and produces:
#   • KS Site Heatmap     — combined + motif-fallback scores (pheatmap)
#   • Kinase Activity     — median substrate Z-score per kinase (pheatmap)
#
# Modes
#   Pre-filter ON  : user selects sites (and optionally kinases) from dropdowns
#   Auto (DE-guided): sites split by padj<0.05 & log2FC>0 / <0 → two heatmaps
#   Fallback        : no DE columns present → all sites scored together
# ==============================================================================

server_ks_scoring <- function(input, output, session, rv, logger) {

  # ---------------------------------------------------------------------------
  # Reference data — reactive on species selector (human / mouse)
  # ---------------------------------------------------------------------------
  # Built on first use from the installed PhosR package (see ensure_phosr_reference() in global.R)
  data_dir <- PHOSR_REF_DIR

  .ks_ref_ready <- function() {
    if (!all(file.exists(phosr_ref_files())))
      showNotification("Preparing PhosR reference data (first use only)...", type = "message", duration = 5)
    ok <- ensure_phosr_reference(logger)
    if (!ok) showNotification(
      tagList("Kinase-Site Analysis needs the PhosR package. Install it in R with ",
              tags$code("BiocManager::install('PhosR')"), ", then restart OSMSS."),
      type = "error", duration = NULL)
    ok
  }

  .kinase_family_cache <- NULL
  .kinase_family <- function() {
    if (is.null(.kinase_family_cache)) {
      req(.ks_ref_ready())
      .kinase_family_cache <<- readRDS(phosr_ref_files()[["kinase_family"]])
    }
    .kinase_family_cache
  }

  .ks_ref <- reactive({
    req(.ks_ref_ready())
    sp  <- input$ks_species %||% "human"
    ml  <- readRDS(file.path(data_dir, paste0("phosr_motif_",    sp, ".rds")))
    sl  <- readRDS(file.path(data_dir, paste0("phosr_substrate_", sp, ".rds")))
    # Mouse motif weights: human NumInputSeq (default) or PhosR's native mouse weights
    if (sp == "mouse" && identical(input$ks_mouse_weights, "mouse"))
      ml$NumInputSeq <- readRDS(phosr_ref_files()[["weights_mouse_native"]])
    logger(paste0("[KS] Reference: ", sp,
                  if (sp == "mouse") paste0(" | motif weights: ",
                                            if (identical(input$ks_mouse_weights, "mouse")) "mouse (PhosR native)" else "human (default)")))
    list(motif_list     = ml,
         substrate_list = sl,
         num_input_seq  = ml$NumInputSeq,
         all_kinase_ids = names(ml$NumInputSeq))
  })

  # ---------------------------------------------------------------------------
  # Module state
  # ---------------------------------------------------------------------------
  ks_rv <- reactiveValues(
    df                   = NULL,     # parsed CSV (data.frame)
    mat_raw              = NULL,     # sites × samples matrix (after dedup)
    seqs_all             = NULL,     # named: PhosR_ID -> 15-mer flanking sequence
    has_de_cols          = FALSE,    # TRUE when log2FC / adj.pvalue columns present
    plex_labels          = NULL,     # named char: sample -> Mixture (TMT only; NULL for DIA)
    condition_labels     = NULL,     # named char: sample -> Condition (TMT only; NULL for DIA)
    mode                 = NULL,     # "prefilter" | "auto" | "all"
    result_up            = NULL,     # list from .run_pipeline() for up / pre-filter
    result_down          = NULL,     # list from .run_pipeline() for down (auto only)
    result_activity_full = NULL,     # full-dataset activity result (multi-condition rerun)
    activity_mode        = "de",    # "de" | "full" — which result the activity tab shows
    status               = "Upload a PhosR CSV to begin.",
    running              = FALSE
  )

  # ---------------------------------------------------------------------------
  # Helpers — identical to run_ks_scoring.R
  # ---------------------------------------------------------------------------
  .score_motif <- function(seqwin_vec, pwm) {
    aa_order <- rownames(pwm)
    vapply(seqwin_vec, function(seq15) {
      if (nchar(seq15) < 15L) return(0)
      score <- 0
      for (pos in seq_len(15L)) {
        aa <- substr(seq15, pos, pos)
        if (aa %in% aa_order) score <- score + pwm[aa, pos]
      }
      score
    }, numeric(1L))
  }

  .minmax <- function(mat) {
    result <- apply(mat, 2L, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (rng[2L] == rng[1L]) return(rep(0, length(x)))
      (x - rng[1L]) / (rng[2L] - rng[1L])
    })
    # apply on a 1-row matrix returns a named vector (dims dropped); restore shape
    if (!is.matrix(result))
      matrix(result, nrow = nrow(mat), dimnames = dimnames(mat))
    else
      result
  }

  # ---------------------------------------------------------------------------
  # Core scoring pipeline
  #   mat_filt  : pre-filtered (and deduplicated) sites × samples matrix
  #   seqs_filt : named character PhosR_ID -> 15-mer
  # Returns a named list with all score matrices and kinase summary.
  # ---------------------------------------------------------------------------
  .run_pipeline <- function(mat_filt, seqs_filt, ref) {
    motif_list     <- ref$motif_list
    substrate_list <- ref$substrate_list
    num_input_seq  <- ref$num_input_seq
    NUM_MOTIF <- 5L
    NUM_SUB   <- 1L

    # 1. Impute remaining NAs (column mean)
    mat_imp <- mat_filt
    for (j in seq_len(ncol(mat_imp))) {
      col    <- mat_imp[, j]
      na_idx <- which(is.na(col))
      if (length(na_idx) > 0L)
        mat_imp[na_idx, j] <- mean(col, na.rm = TRUE)
    }

    # 2. Row-wise Z-score
    row_means <- rowMeans(mat_imp)
    row_sds   <- apply(mat_imp, 1L, sd)
    row_sds[row_sds == 0 | is.na(row_sds)] <- 1
    mat_std   <- sweep(sweep(mat_imp, 1L, row_means, "-"), 1L, row_sds, "/")

    # 3. Extract 15-mer windows (centre ±7)
    seq_win <- vapply(seqs_filt[rownames(mat_std)], function(s) {
      if (is.na(s) || nchar(s) < 15L)
        return(paste(rep("_", 15L), collapse = ""))
      mid <- (nchar(s) + 1L) %/% 2L
      substr(s, max(1L, mid - 7L), min(nchar(s), mid + 7L))
    }, character(1L))

    # 4. Motif scoring
    motif_kinases  <- names(num_input_seq)[num_input_seq >= NUM_MOTIF]
    motif_filtered <- motif_list[motif_kinases]
    motif_mat      <- matrix(NA_real_, nrow(mat_std), length(motif_kinases),
                             dimnames = list(rownames(mat_std), motif_kinases))
    for (i in seq_along(motif_kinases)) {
      k   <- motif_kinases[i]
      pwm <- motif_filtered[[k]]
      if (!is.matrix(pwm) || ncol(pwm) < 15L) next
      motif_mat[, i] <- .score_motif(seq_win, pwm)
    }
    motif_mat <- .minmax(motif_mat)

    # 5. Profile scoring
    n_subs_detected <- vapply(substrate_list, function(subs)
      sum(subs %in% rownames(mat_std)), integer(1L))

    ks_profiles <- lapply(substrate_list, function(subs) {
      present <- intersect(subs, rownames(mat_std))
      if (length(present) == 0L) return(NA_real_)
      if (length(present) == 1L) return(mat_std[present, ])
      apply(mat_std[present, , drop = FALSE], 2L, median)
    })

    profile_kinases  <- names(n_subs_detected)[n_subs_detected >= NUM_SUB]
    profile_profiles <- ks_profiles[profile_kinases]
    valid_profiles   <- !vapply(profile_profiles,
                                function(x) all(is.na(x)), logical(1L))
    profile_kinases  <- profile_kinases[valid_profiles]
    profile_profiles <- profile_profiles[profile_kinases]

    if (length(profile_kinases) == 0L) {
      # No PSP substrates in this site set — skip profile scoring entirely.
      # Combined scoring will fall back to motif-only via the supplement step.
      profile_mat_ref <- matrix(NA_real_, 0L, ncol(mat_std),
                                dimnames = list(NULL, colnames(mat_std)))
      profile_mat     <- matrix(NA_real_, nrow(mat_std), 0L,
                                dimnames = list(rownames(mat_std), character(0)))
    } else {
      profile_mat_ref <- do.call(rbind, profile_profiles)
      rownames(profile_mat_ref) <- profile_kinases
      profile_mat <- matrix(NA_real_, nrow(mat_std), length(profile_kinases),
                            dimnames = list(rownames(mat_std), profile_kinases))
      for (i in seq_len(nrow(mat_std))) {
        site_vec <- mat_std[i, ]
        cors     <- apply(profile_mat_ref, 1L,
                          function(kin_vec) suppressWarnings(cor(site_vec, kin_vec)))
        cors[is.na(cors)] <- 0
        profile_mat[i, ]  <- (cors + 1) / 2
      }
    }

    # 6. Combined score (log-rank weighted average)
    shared_kinases <- intersect(colnames(motif_mat), colnames(profile_mat))
    w1 <- log(rank(num_input_seq[shared_kinases]) + 1)
    w2 <- log(rank(n_subs_detected[shared_kinases]) + 1)
    w3 <- w1 + w2
    combined_mat <- matrix(NA_real_, nrow(mat_std), length(shared_kinases),
                           dimnames = list(rownames(mat_std), shared_kinases))
    for (i in seq_along(shared_kinases)) {
      k <- shared_kinases[i]
      combined_mat[, i] <- (w1[i] / w3[i]) * motif_mat[, k] +
                           (w2[i] / w3[i]) * profile_mat[, k]
    }

    # 7. Motif fallback: supplement combined with motif-only kinases
    motif_supplement <- setdiff(colnames(motif_mat), colnames(combined_mat))
    if (length(motif_supplement) > 0L) {
      extra       <- motif_mat[rownames(combined_mat), motif_supplement, drop = FALSE]
      score_mat   <- cbind(combined_mat, extra)
      score_label <- "KS Score (combined + motif fallback)"
    } else {
      score_mat   <- combined_mat
      score_label <- "KS Combined Score"
    }

    # 8. Kinase activity matrix: kinases × samples (median substrate Z-score)
    #    For kinases with PSP substrates in data: use true median profile.
    #    For all other scored kinases (motif-fallback): compute proxy profile =
    #    mean Z-score of the top-3 highest-scoring sites. This lets the activity
    #    heatmap show all scored kinases, not just the PSP-confirmed few.
    TOP_PROXY_SITES <- 3L
    all_scored_kinases <- colnames(score_mat)
    psp_kinases  <- intersect(shared_kinases, rownames(profile_mat_ref))
    proxy_needed <- setdiff(all_scored_kinases, psp_kinases)

    # Build base: PSP-confirmed rows
    psp_mat <- if (length(psp_kinases) > 0L)
      profile_mat_ref[psp_kinases, , drop = FALSE]
    else
      matrix(NA_real_, 0L, ncol(mat_std),
             dimnames = list(character(0), colnames(mat_std)))

    # Build proxy rows for motif-fallback kinases
    proxy_mat <- matrix(NA_real_, length(proxy_needed), ncol(mat_std),
                        dimnames = list(proxy_needed, colnames(mat_std)))
    if (length(proxy_needed) > 0L && nrow(mat_std) > 0L) {
      for (k in proxy_needed) {
        if (!(k %in% colnames(score_mat))) next
        kin_scores  <- score_mat[, k]
        top_sorted  <- sort(kin_scores, decreasing = TRUE, na.last = NA)
        n_top       <- min(TOP_PROXY_SITES, length(top_sorted))
        if (n_top == 0L) next                  # all-NA column (e.g. mouse kinase missing PWM)
        top_sites   <- names(top_sorted[seq_len(n_top)])
        proxy_mat[k, ] <- colMeans(mat_std[top_sites, , drop = FALSE], na.rm = TRUE)
      }
    }

    ks_activity_mat  <- rbind(psp_mat, proxy_mat)
    ks_activity_type <- c(rep("PSP confirmed",  nrow(psp_mat)),
                          rep("Proxy (top sites)", nrow(proxy_mat)))
    names(ks_activity_type) <- rownames(ks_activity_mat)

    # 9. Kinase summary — covers ALL scored kinases (combined + motif fallback)
    #    so that the activity heatmap renderer can intersect against it.
    all_scored_for_summary <- colnames(score_mat)
    mean_all <- colMeans(score_mat[, all_scored_for_summary, drop = FALSE], na.rm = TRUE)
    kin_summary   <- data.frame(
      Kinase        = all_scored_for_summary,
      N_Motif_Seqs  = num_input_seq[all_scored_for_summary],
      N_Subs_Det    = n_subs_detected[all_scored_for_summary],
      Mean_Combined = mean_all,
      stringsAsFactors = FALSE
    )
    kin_summary <- kin_summary[order(kin_summary$Mean_Combined, decreasing = TRUE), ]

    list(
      score_mat        = score_mat,
      motif_mat        = motif_mat,
      combined_mat     = combined_mat,
      profile_mat      = profile_mat,
      ks_activity_mat  = ks_activity_mat,
      ks_activity_type = ks_activity_type,
      kin_summary      = kin_summary,
      score_label      = score_label,
      n_subs_detected  = n_subs_detected,
      shared_kinases   = shared_kinases
    )
  }

  # ---------------------------------------------------------------------------
  # CSV upload: parse, deduplicate, update dropdowns
  # ---------------------------------------------------------------------------
  observeEvent(input$ks_csv_upload, {
    req(input$ks_csv_upload)
    tryCatch({
      df <- read.csv(input$ks_csv_upload$datapath,
                     stringsAsFactors = FALSE, check.names = FALSE)

      required <- c("PhosR_ID", "Flanking_Sequence_15mer")
      if (!all(required %in% colnames(df))) {
        ks_rv$status <- paste0(
          "Missing required columns: ",
          paste(setdiff(required, colnames(df)), collapse = ", "))
        return()
      }

      meta_cols   <- c("PhosR_ID", "Gene", "Protein", "Residue",
                       "Flanking_Sequence_15mer", "log2FC", "adj.pvalue")
      sample_cols <- setdiff(colnames(df), meta_cols)

      if (length(sample_cols) == 0L) {
        ks_rv$status <- "No sample columns found. Expected columns beyond metadata."
        return()
      }

      # Detect and extract TMT plex metadata row (PhosR_ID == "__plex__")
      # and condition metadata row (PhosR_ID == "__condition__"),
      # both written by server_ptm_de.R for TMT exports.
      plex_row_idx <- which(df$PhosR_ID == "__plex__")
      cond_row_idx <- which(df$PhosR_ID == "__condition__")
      meta_row_idx <- c(plex_row_idx, cond_row_idx)

      if (length(plex_row_idx) > 0L) {
        plex_row        <- df[plex_row_idx[1L], sample_cols, drop = FALSE]
        plex_labels_vec <- setNames(as.character(unlist(plex_row[1L, ])), sample_cols)
        plex_labels_vec <- plex_labels_vec[!is.na(plex_labels_vec)]
        ks_rv$plex_labels <- plex_labels_vec
      } else {
        ks_rv$plex_labels <- NULL
      }

      if (length(cond_row_idx) > 0L) {
        cond_row        <- df[cond_row_idx[1L], sample_cols, drop = FALSE]
        cond_labels_vec <- setNames(as.character(unlist(cond_row[1L, ])), sample_cols)
        cond_labels_vec <- cond_labels_vec[!is.na(cond_labels_vec)]
        ks_rv$condition_labels <- cond_labels_vec
      } else {
        ks_rv$condition_labels <- NULL
      }

      if (length(meta_row_idx) > 0L) {
        df <- df[-meta_row_idx, , drop = FALSE]
        # Coerce sample columns back to numeric — read.csv inferred character
        # because the metadata rows contained strings (e.g. "F1", "TPA 0 hr").
        df[, sample_cols] <- lapply(df[, sample_cols],
                                    function(x) suppressWarnings(as.numeric(x)))
      }

      mat_raw        <- as.matrix(df[, sample_cols, drop = FALSE])
      rownames(mat_raw) <- df$PhosR_ID

      # Deduplicate by PhosR_ID (median per sample)
      if (anyDuplicated(rownames(mat_raw))) {
        uid     <- unique(rownames(mat_raw))
        mat_raw <- do.call(rbind, lapply(uid, function(pid) {
          rows <- mat_raw[rownames(mat_raw) == pid, , drop = FALSE]
          if (nrow(rows) == 1L) return(rows)
          matrix(
            apply(rows, 2L, function(x) {
              if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
            }),
            nrow = 1L, dimnames = list(pid, colnames(rows))
          )
        }))
        df <- df[!duplicated(df$PhosR_ID), , drop = FALSE]
      }

      seqs_all <- setNames(df$Flanking_Sequence_15mer, df$PhosR_ID)

      ks_rv$df          <- df
      ks_rv$mat_raw              <- mat_raw
      ks_rv$seqs_all             <- seqs_all
      ks_rv$has_de_cols          <- all(c("log2FC", "adj.pvalue") %in% colnames(df))
      ks_rv$result_up            <- NULL
      ks_rv$result_down          <- NULL
      ks_rv$result_activity_full <- NULL
      ks_rv$activity_mode        <- "de"
      ks_rv$mode                 <- NULL
      plex_note <- if (!is.null(ks_rv$plex_labels))
        paste0(" | TMT: ", length(unique(ks_rv$plex_labels)), " plex(es) detected — RUV correction will be applied")
      else ""
      ks_rv$status               <- paste0(
        "Loaded: ", nrow(mat_raw), " sites, ", length(sample_cols), " samples.", plex_note)

      # Populate dropdowns (server-side for performance with large site lists)
      updateSelectizeInput(session, "ks_site_filter",
                           choices = df$PhosR_ID, server = TRUE)
      updateSelectizeInput(session, "ks_kinase_filter",
                           choices = .ks_ref()$all_kinase_ids, server = TRUE)
      # Populate sample-order picker with current column order (client-side;
      # list is short enough and drag_drop plugin requires client-side mode)
      updateSelectizeInput(session, "ks_sample_order",
                           choices  = sample_cols,
                           selected = sample_cols,
                           server   = FALSE)
    }, error = function(e) {
      ks_rv$status <- paste("Error reading CSV:", conditionMessage(e))
    })
  })

  # ---------------------------------------------------------------------------
  # Run button — orchestrates the scoring pipeline
  # ---------------------------------------------------------------------------
  observeEvent(input$ks_run_btn, {
    if (is.null(ks_rv$df) || is.null(ks_rv$mat_raw)) {
      ks_rv$status <- "Please upload a PhosR CSV in Step 1 first."
      return()
    }

    ks_rv$running <- TRUE
    ks_rv$status  <- "Running KS scoring... (this may take 15-60 seconds)"

    showNotification("KS scoring running — this may take 15–60 seconds...",
                     id = "ks_running_notif", type = "message", duration = NULL)

    withProgress(message = "KS Scoring", value = 0, {

    tryCatch({
      df       <- ks_rv$df
      mat_raw  <- ks_rv$mat_raw
      seqs_all <- ks_rv$seqs_all

      incProgress(0.05, detail = "Filtering sites by missingness...")
      # Missingness filter: keep sites with >=50% non-NA
      pct_valid <- rowMeans(!is.na(mat_raw))
      mat_filt  <- mat_raw[pct_valid >= 0.5, , drop = FALSE]
      seqs_filt <- seqs_all[rownames(mat_filt)]

      # ---------------------------------------------------------------------------
      # RUV (Remove Unwanted Variation) plex correction — TMT only
      # Applied when the CSV contains a "__plex__" metadata row (written by the
      # OSMSS TMT PhosR export).  Uses PhosR::RUVphospho() with k=1 (one factor
      # of unwanted variation, appropriate for a single plex/batch covariate).
      # mat_filt is replaced with the corrected matrix; NA structure is preserved.
      # ---------------------------------------------------------------------------
      plex_labels      <- ks_rv$plex_labels
      condition_labels <- ks_rv$condition_labels
      ruv_applied <- FALSE
      logger(paste0("[KS RUV] plex_labels present: ", !is.null(plex_labels),
                    " | ruv available: ", requireNamespace("ruv", quietly=TRUE)))
      if (!is.null(plex_labels))
        logger(paste0("[KS RUV] plex_labels n=", length(plex_labels),
                      ", unique plexes: ", paste(unique(plex_labels), collapse=",")))
      if (!is.null(plex_labels) && length(plex_labels) >= 2L &&
          requireNamespace("ruv", quietly = TRUE)) {
        tryCatch({
          incProgress(0, detail = "Applying RUV plex correction...")
          logger(paste0("[KS RUV] condition_labels: ",
                        if (!is.null(condition_labels))
                          paste0("n=", length(condition_labels), ", unique: ",
                                 paste(unique(condition_labels), collapse=","))
                        else "NULL"))
          logger(paste0("[KS RUV] mat_filt dims: ", nrow(mat_filt), " rows x ", ncol(mat_filt), " cols"))
          logger(paste0("[KS RUV] plex_labels names (first 3): ",
                        paste(head(names(plex_labels), 3), collapse=" | ")))
          logger(paste0("[KS RUV] mat_filt colnames (first 3): ",
                        paste(head(colnames(mat_filt), 3), collapse=" | ")))
          common_cols  <- intersect(colnames(mat_filt), names(plex_labels))
          logger(paste0("[KS RUV] common_cols matched: ", length(common_cols),
                        " / ", ncol(mat_filt)))
          if (length(common_cols) >= 2L) {
            mat_for_ruv  <- mat_filt[, common_cols, drop = FALSE]

            na_mask     <- is.na(mat_for_ruv)
            mat_imputed <- mat_for_ruv
            for (j in seq_len(ncol(mat_imputed))) {
              col <- mat_imputed[, j]
              mat_imputed[is.na(col), j] <- mean(col, na.rm = TRUE)
            }
            keep_rows   <- rowSums(!is.na(mat_imputed)) > 0L
            mat_imputed <- mat_imputed[keep_rows, , drop = FALSE]
            na_mask_sub <- na_mask[keep_rows, , drop = FALSE]

            # M = design matrix of variation to PRESERVE (biological signal).
            # If condition labels are available, encode them so RUV doesn't
            # accidentally absorb biological differences between conditions.
            # Otherwise fall back to intercept-only (removes strongest factor).
            cond_vec <- condition_labels[common_cols]
            if (!is.null(cond_vec) && length(unique(cond_vec[!is.na(cond_vec)])) >= 2L) {
              cond_fac  <- factor(cond_vec)
              M_design  <- model.matrix(~ cond_fac - 1)   # one-hot: required by RUVIII
              logger("[KS RUV] using condition design matrix (one-hot)")
            } else {
              M_design  <- matrix(1, nrow = ncol(mat_imputed), ncol = 1)
              logger("[KS RUV] using intercept-only matrix")
            }
            logger(paste0("[KS RUV] mat_imputed dims before RUV: ",
                          nrow(mat_imputed), " rows x ", ncol(mat_imputed), " cols"))
            # Call RUVIII directly (from the 'ruv' package — no PhosR required).
            # Y must be samples x sites; ctl = all sites as pseudo-controls
            # (standard practice for phosphoproteomics without stable reference sites).
            ctl_idx       <- seq_len(nrow(mat_imputed))
            mat_corrected <- t(ruv::RUVIII(
              Y   = t(mat_imputed),
              M   = M_design,
              ctl = ctl_idx,
              k   = 1L
            ))
            logger(paste0("[KS RUV] mat_corrected dims after RUV: ",
                          nrow(mat_corrected), " rows x ", ncol(mat_corrected), " cols"))
            mat_corrected[na_mask_sub] <- NA_real_
            mat_filt[rownames(mat_corrected), common_cols] <- mat_corrected
            ruv_applied <- TRUE
            logger("[KS Scoring] RUVphospho plex correction applied (k=1)")
          }
        }, error = function(e) {
          logger(paste0("[KS Scoring] RUVphospho failed, proceeding without correction: ",
                        e$message))
        })
      }

      ref                 <- .ks_ref()
      pre_filter_on       <- isTRUE(input$ks_prefilter_on)
      show_all_kinases    <- isTRUE(input$ks_show_all_kinases)
      show_all_sig_sites  <- isTRUE(input$ks_show_all_sig_sites)
      TOP_K               <- if (show_all_kinases) Inf else 60L
      TOP_N_SITES         <- if (show_all_sig_sites) Inf else 50L
      selected_sites_in   <- input$ks_site_filter      # character (may be NULL)
      selected_kinases_in <- input$ks_kinase_filter     # character (may be NULL)

      # ----- Branch: pre-filter mode (explicit site list) -----
      # Kinase-only filter (no sites selected) falls through to the auto/fallback
      # branches below, which apply the kinase filter there.
      if (pre_filter_on && length(selected_sites_in) > 0L) {
        found <- selected_sites_in[selected_sites_in %in% rownames(mat_filt)]
        if (length(found) == 0L) {
          ks_rv$status  <- paste0(
            "None of the ", length(selected_sites_in),
            " selected sites passed the 50% missingness filter.")
          ks_rv$running <- FALSE
          return()
        }
        mat_use  <- mat_filt[found, , drop = FALSE]
        seqs_use <- seqs_filt[found]
        incProgress(0.5, detail = paste0("Scoring ", length(found), " site(s)..."))
        res      <- .run_pipeline(mat_use, seqs_use, ref)
        incProgress(0.4, detail = "Ranking kinases...")

        # Kinase selection: user-specified or auto top-K by mean score
        # NOTE: max is wrong here — after per-kinase minmax normalisation every
        # kinase's max == 1.0, so sort() falls back to alphabetical order.
        # Mean across sites correctly ranks kinases that score broadly high.
        if (length(selected_kinases_in) > 0L) {
          found_kin          <- selected_kinases_in[
            selected_kinases_in %in% colnames(res$score_mat)]
          user_kinase_filter <- TRUE
        } else {
          kin_mean_sorted <- sort(colMeans(res$score_mat, na.rm = TRUE),
                                  decreasing = TRUE, na.last = NA)
          n_keep    <- min(if (is.finite(TOP_K)) TOP_K else length(kin_mean_sorted),
                          length(kin_mean_sorted))
          found_kin <- names(kin_mean_sorted[seq_len(n_keep)])
          user_kinase_filter <- FALSE
        }

        missing_sites        <- setdiff(selected_sites_in, rownames(mat_filt))
        res$selected_sites   <- found
        res$selected_kinases <- found_kin
        res$missing_sites    <- missing_sites
        res$mode             <- "prefilter"
        res$cluster_rows     <- FALSE
        res$cluster_cols     <- !user_kinase_filter

        ks_rv$result_up   <- res
        ks_rv$result_down <- NULL
        ks_rv$mode        <- "prefilter"

        n_miss <- length(missing_sites)
        ks_rv$status <- paste0(
          "Done. ", length(found), " site(s) scored",
          if (n_miss > 0L) paste0(" (", n_miss, " not in dataset)") else "",
          ", ", length(found_kin), " kinase(s) shown.")

      # ----- Branch: auto DE-guided (two directional heatmaps) -----
      } else if (ks_rv$has_de_cols) {
        df_sig <- df[!is.na(df$adj.pvalue) & !is.na(df$log2FC) &
                     df$adj.pvalue < 0.05 &
                     df$PhosR_ID %in% rownames(mat_filt), ]

        df_up   <- df_sig[df_sig$log2FC > 0, ]
        df_down <- df_sig[df_sig$log2FC < 0, ]
        df_up   <- df_up[order(abs(df_up$log2FC),   decreasing = TRUE), ]
        df_down <- df_down[order(abs(df_down$log2FC), decreasing = TRUE), ]

        TOP_N   <- 50L
        ids_up   <- df_up$PhosR_ID[seq_len(min(
          if (is.finite(TOP_N_SITES)) TOP_N_SITES else nrow(df_up), nrow(df_up)))]
        ids_down <- df_down$PhosR_ID[seq_len(min(
          if (is.finite(TOP_N_SITES)) TOP_N_SITES else nrow(df_down), nrow(df_down)))]

        .score_directional <- function(ids, dir_mode) {
          if (length(ids) == 0L) return(NULL)
          m   <- mat_filt[ids, , drop = FALSE]
          s   <- seqs_filt[ids]
          res <- .run_pipeline(m, s, ref)
          # Kinase filter: use user-specified list when pre-filter is on with kinases
          # selected, otherwise auto top-K by mean score.
          if (pre_filter_on && length(selected_kinases_in) > 0L) {
            found_kin  <- selected_kinases_in[
              selected_kinases_in %in% colnames(res$score_mat)]
            clust_cols <- FALSE
          } else {
            kin_mean_sorted <- sort(colMeans(res$score_mat, na.rm = TRUE),
                                    decreasing = TRUE, na.last = NA)
            n_keep     <- min(if (is.finite(TOP_K)) TOP_K else length(kin_mean_sorted),
                             length(kin_mean_sorted))
            found_kin  <- names(kin_mean_sorted[seq_len(n_keep)])
            clust_cols <- TRUE
          }
          res$selected_sites   <- ids
          res$selected_kinases <- found_kin
          res$missing_sites    <- character(0)
          res$mode             <- dir_mode
          res$cluster_rows     <- TRUE
          res$cluster_cols     <- clust_cols
          res
        }

        incProgress(0.40, detail = paste0("Scoring ", length(ids_up), " up-regulated site(s)..."))
        ks_rv$result_up   <- .score_directional(ids_up,   "auto_up")
        incProgress(0.40, detail = paste0("Scoring ", length(ids_down), " down-regulated site(s)..."))
        ks_rv$result_down <- .score_directional(ids_down, "auto_down")
        ks_rv$mode        <- "auto"

        ks_rv$status <- paste0(
          "Done. Up: ", length(ids_up), " site(s), Down: ",
          length(ids_down), " site(s) (padj<0.05",
          if (show_all_sig_sites) ", all sig sites" else ", top 50 by |log2FC|",
          ").")  

      # ----- Fallback: no DE columns, score all sites -----
      } else {
        incProgress(0.5, detail = paste0("Scoring all ", nrow(mat_filt), " site(s)..."))
        res <- .run_pipeline(mat_filt, seqs_filt, ref)
        incProgress(0.4, detail = "Ranking kinases...")
        kin_mean_sorted <- sort(colMeans(res$score_mat, na.rm = TRUE),
                                decreasing = TRUE, na.last = NA)
        n_keep   <- min(if (is.finite(TOP_K)) TOP_K else length(kin_mean_sorted),
                       length(kin_mean_sorted))

        # Cap displayed sites by TOP_N_SITES (default 50) ranked by row variance.
        # No log2FC is available here, so variance is the best proxy for
        # "most interesting" sites. show_all_sig_sites overrides this cap.
        all_site_ids <- rownames(mat_filt)
        if (is.finite(TOP_N_SITES) && length(all_site_ids) > TOP_N_SITES) {
          row_vars      <- apply(mat_filt, 1L, var, na.rm = TRUE)
          row_vars[is.na(row_vars)] <- 0
          top_by_var    <- names(sort(row_vars, decreasing = TRUE)[seq_len(TOP_N_SITES)])
          display_sites <- top_by_var
          site_note     <- paste0("top ", TOP_N_SITES, " by variance")
        } else {
          display_sites <- all_site_ids
          site_note     <- "all sites"
        }

        # Kinase filter: use user-specified list when pre-filter is on with kinases
        # selected, otherwise auto top-K by mean score.
        if (pre_filter_on && length(selected_kinases_in) > 0L) {
          display_kinases <- selected_kinases_in[
            selected_kinases_in %in% colnames(res$score_mat)]
          clust_cols <- FALSE
        } else {
          display_kinases <- names(kin_mean_sorted[seq_len(n_keep)])
          clust_cols <- TRUE
        }
        res$selected_sites   <- display_sites
        res$selected_kinases <- display_kinases
        res$missing_sites    <- character(0)
        res$mode             <- "all"
        res$cluster_rows     <- TRUE
        res$cluster_cols     <- clust_cols

        ks_rv$result_up   <- res
        ks_rv$result_down <- NULL
        ks_rv$mode        <- "all"

        ks_rv$status <- paste0(
          "Done. ", length(display_sites), " site(s) shown (", site_note,
          ", no DE columns in CSV — run PTM DE before exporting to enable DE-guided mode).")
      }

      ks_rv$running <- FALSE
      removeNotification("ks_running_notif")
      showNotification(ks_rv$status, type = "message", duration = 8)

    }, error = function(e) {
      ks_rv$status  <- paste("Scoring error:", conditionMessage(e))
      ks_rv$running <- FALSE
      removeNotification("ks_running_notif")
      showNotification(paste("KS Scoring failed:", conditionMessage(e)),
                       type = "error", duration = 15)
    })

    }) # end withProgress
  })

  # ---------------------------------------------------------------------------
  # Full-dataset activity rerun (multi-condition button)
  # Runs .run_pipeline() on ALL quantified sites (ignores DE filter).
  # Only the activity matrix result is stored; up/down heatmaps are unchanged.
  # ---------------------------------------------------------------------------
  observeEvent(input$ks_run_full_activity, {
    req(!is.null(ks_rv$mat_raw))

    ks_rv$running <- TRUE
    showNotification(
      "Recomputing kinase activity with all PTMs — this may take 30–90 seconds...",
      id = "ks_full_notif", type = "message", duration = NULL)

    withProgress(message = "KS Activity (all PTMs)", value = 0, {
    tryCatch({
      mat_raw  <- ks_rv$mat_raw
      seqs_all <- ks_rv$seqs_all
      ref      <- .ks_ref()

      incProgress(0.1, detail = "Filtering sites by missingness...")
      pct_valid <- rowMeans(!is.na(mat_raw))
      mat_filt  <- mat_raw[pct_valid >= 0.5, , drop = FALSE]
      seqs_filt <- seqs_all[rownames(mat_filt)]

      # RUV plex correction (same logic as main run handler)
      plex_labels      <- ks_rv$plex_labels
      condition_labels <- ks_rv$condition_labels
      logger(paste0("[KS Full Activity RUV] plex_labels present: ", !is.null(plex_labels),
                    " | ruv available: ", requireNamespace("ruv", quietly=TRUE)))
      if (!is.null(plex_labels))
        logger(paste0("[KS Full Activity RUV] plex_labels n=", length(plex_labels),
                      ", unique plexes: ", paste(unique(plex_labels), collapse=",")))
      if (!is.null(plex_labels) && length(plex_labels) >= 2L &&
          requireNamespace("ruv", quietly = TRUE)) {
        tryCatch({
          logger(paste0("[KS Full Activity RUV] condition_labels: ",
                        if (!is.null(condition_labels))
                          paste0("n=", length(condition_labels), ", unique: ",
                                 paste(unique(condition_labels), collapse=","))
                        else "NULL"))
          logger(paste0("[KS Full Activity RUV] mat_filt dims: ",
                        nrow(mat_filt), " rows x ", ncol(mat_filt), " cols"))
          logger(paste0("[KS Full Activity RUV] plex_labels names (first 3): ",
                        paste(head(names(plex_labels), 3), collapse=" | ")))
          logger(paste0("[KS Full Activity RUV] mat_filt colnames (first 3): ",
                        paste(head(colnames(mat_filt), 3), collapse=" | ")))
          common_cols <- intersect(colnames(mat_filt), names(plex_labels))
          logger(paste0("[KS Full Activity RUV] common_cols matched: ",
                        length(common_cols), " / ", ncol(mat_filt)))
          if (length(common_cols) >= 2L) {
            mat_for_ruv <- mat_filt[, common_cols, drop = FALSE]
            na_mask     <- is.na(mat_for_ruv)
            mat_imputed <- mat_for_ruv
            for (j in seq_len(ncol(mat_imputed))) {
              col <- mat_imputed[, j]
              mat_imputed[is.na(col), j] <- mean(col, na.rm = TRUE)
            }
            keep_rows      <- rowSums(!is.na(mat_imputed)) > 0L
            mat_imputed    <- mat_imputed[keep_rows, , drop = FALSE]
            na_mask_sub    <- na_mask[keep_rows, , drop = FALSE]
            cond_vec <- condition_labels[common_cols]
            if (!is.null(cond_vec) && length(unique(cond_vec[!is.na(cond_vec)])) >= 2L) {
              cond_fac  <- factor(cond_vec)
              M_design  <- model.matrix(~ cond_fac - 1)   # one-hot: required by RUVIII
            } else {
              M_design  <- matrix(1, nrow = ncol(mat_imputed), ncol = 1)
            }
            logger(paste0("[KS Full Activity RUV] mat_imputed dims before RUV: ",
                          nrow(mat_imputed), " rows x ", ncol(mat_imputed), " cols"))
            ctl_idx       <- seq_len(nrow(mat_imputed))
            mat_corrected <- t(ruv::RUVIII(
              Y   = t(mat_imputed),
              M   = M_design,
              ctl = ctl_idx,
              k   = 1L
            ))
            logger(paste0("[KS Full Activity RUV] mat_corrected dims after RUV: ",
                          nrow(mat_corrected), " rows x ", ncol(mat_corrected), " cols"))
            mat_corrected[na_mask_sub] <- NA_real_
            mat_filt[rownames(mat_corrected), common_cols] <- mat_corrected
            logger("[KS Full Activity] RUVphospho plex correction applied (k=1)")
          }
        }, error = function(e) {
          logger(paste0("[KS Full Activity] RUVphospho failed: ", e$message))
        })
      }

      incProgress(0.2, detail = paste0("Scoring all ", nrow(mat_filt), " sites..."))
      res_full <- .run_pipeline(mat_filt, seqs_filt, ref)

      # Select kinases: honour show_all_kinases toggle
      show_all_kinases <- isTRUE(input$ks_show_all_kinases)
      TOP_K <- if (show_all_kinases) Inf else 60L
      kin_mean_sorted <- sort(colMeans(res_full$score_mat, na.rm = TRUE),
                              decreasing = TRUE, na.last = NA)
      n_keep   <- min(if (is.finite(TOP_K)) TOP_K else length(kin_mean_sorted),
                     length(kin_mean_sorted))
      res_full$selected_sites   <- rownames(mat_filt)
      res_full$selected_kinases <- names(kin_mean_sorted[seq_len(n_keep)])
      res_full$missing_sites    <- character(0)
      res_full$mode             <- "all"
      res_full$cluster_rows     <- TRUE
      res_full$cluster_cols     <- TRUE

      incProgress(0.7, detail = "Done.")
      ks_rv$result_activity_full <- res_full
      ks_rv$activity_mode        <- "full"
      ks_rv$running              <- FALSE
      removeNotification("ks_full_notif")
      showNotification(
        paste0("Full-dataset activity computed: ", nrow(mat_filt),
               " sites, ", length(res_full$selected_kinases), " kinases."),
        type = "message", duration = 8)

    }, error = function(e) {
      ks_rv$running <- FALSE
      removeNotification("ks_full_notif")
      showNotification(paste("Full-activity rerun failed:", conditionMessage(e)),
                       type = "error", duration = 15)
    })
    }) # end withProgress
  })

  # ---------------------------------------------------------------------------
  # pheatmap render helpers
  # ---------------------------------------------------------------------------

  # Build row labels from PhosR_IDs ("GENE RESIDUE" format)
  .row_labels <- function(ids, df_local) {
    gene_map <- setNames(df_local$Gene,    df_local$PhosR_ID)
    res_map  <- setNames(df_local$Residue, df_local$PhosR_ID)
    vapply(ids, function(pid) {
      gene <- gene_map[pid]
      r    <- res_map[pid]
      if (!is.na(gene) && nchar(trimws(gene)) > 0L)
        paste0(trimws(gene), " ", r)
      else
        pid
    }, character(1L))
  }

  # Render KS site heatmap (pheatmap); called inside renderPlot
  .render_ks_heatmap <- function(res) {
    req(!is.null(res))

    score_mat        <- res$score_mat
    selected_sites   <- res$selected_sites
    selected_kinases <- res$selected_kinases
    missing_sites    <- res$missing_sites

    req(length(selected_sites) > 0L, length(selected_kinases) > 0L)

    heat_mat <- score_mat[selected_sites, selected_kinases, drop = FALSE]

    # Append N/A rows for sites not found / below quantification cutoff
    if (length(missing_sites) > 0L) {
      na_rows   <- matrix(NA_real_,
                          nrow     = length(missing_sites),
                          ncol     = ncol(heat_mat),
                          dimnames = list(missing_sites, colnames(heat_mat)))
      heat_disp <- rbind(heat_mat, na_rows)
    } else {
      heat_disp <- heat_mat
    }

    # Row labels
    found_labs <- .row_labels(selected_sites, ks_rv$df)
    if (length(missing_sites) > 0L) {
      miss_labs  <- setNames(
        paste0("[N/A] ", sub(";$", "", missing_sites)),
        missing_sites)
      row_labels <- c(found_labs, miss_labs)
    } else {
      row_labels <- found_labs
    }

    # Kinase family annotation (columns)
    kfam_df <- as.data.frame(.kinase_family(), stringsAsFactors = FALSE)
    ann_col <- data.frame(
      Group  = kfam_df[match(selected_kinases, rownames(kfam_df)), "kinase_group"],
      Family = kfam_df[match(selected_kinases, rownames(kfam_df)), "kinase_family"],
      stringsAsFactors = TRUE
    )
    rownames(ann_col) <- selected_kinases
    ann_col[is.na(ann_col$Group),  "Group"]  <- "Unknown"
    ann_col[is.na(ann_col$Family), "Family"] <- "Unknown"

    # PSP / status row annotation
    psp_known_sites <- unique(unlist(.ks_ref()$substrate_list, use.names = FALSE))
    psp_flag <- ifelse(
      rownames(heat_disp) %in% psp_known_sites, "PSP known",
      ifelse(rownames(heat_disp) %in% missing_sites, "not in dataset",
             "motif only"))

    psp_lev <- c("PSP known", "motif only", "not in dataset")
    ann_row <- data.frame(
      PSP = factor(psp_flag, levels = psp_lev),
      row.names = rownames(heat_disp)
    )
    ann_colors <- list(
      PSP = c("PSP known"      = "#E69F00",
              "motif only"     = "#CCCCCC",
              "not in dataset" = "#EEEEEE")
    )
    if (length(missing_sites) > 0L) {
      status_vec     <- c(rep("found",              nrow(heat_mat)),
                          rep("below cutoff / N/A", length(missing_sites)))
      ann_row$Status <- factor(status_vec, levels = c("found", "below cutoff / N/A"))
      ann_colors$Status <- c("found" = "#F5F5F5", "below cutoff / N/A" = "#CC3333")
    }

    # Title
    dir_suffix <- switch(res$mode,
      auto_up   = " — Upregulated (padj<0.05, log2FC>0)",
      auto_down = " — Downregulated (padj<0.05, log2FC<0)",
      prefilter = " — Selected Sites",
      all       = " — All Sites",
      "")
    missing_note <- if (length(missing_sites) > 0L)
      paste0("\n[", length(missing_sites), " site(s) not in dataset: ",
             paste(head(sub(";$", "", missing_sites), 4L), collapse = ", "),
             if (length(missing_sites) > 4L) ", ..." else "", "]")
    else ""

    title_str <- paste0(
      res$score_label, dir_suffix, "\n",
      nrow(heat_mat), " sites \u00d7 ", length(selected_kinases), " kinases",
      missing_note)

    pheatmap::pheatmap(
      mat               = heat_disp,
      labels_row        = row_labels,
      annotation_col    = ann_col,
      annotation_row    = ann_row,
      annotation_colors = ann_colors,
      cluster_rows      = res$cluster_rows && nrow(heat_mat) > 1L,
      cluster_cols      = res$cluster_cols && ncol(heat_mat) > 1L,
      clustering_method = "ward.D2",
      na_col            = "#EEEEEE",
      cellwidth         = 9,
      cellheight        = 9,
      border_color      = "grey85",
      fontsize_row      = 6,
      fontsize_col      = 6,
      main              = title_str
    )
  }

  # Render kinase activity heatmap
  .render_activity_heatmap <- function(res, show_all = FALSE, sample_order = NULL) {
    req(!is.null(res))

    kin_summary      <- res$kin_summary
    ks_activity_mat  <- res$ks_activity_mat
    ks_activity_type <- res$ks_activity_type  # named char: kinase -> "PSP confirmed" / "Proxy"

    # All scored kinases are now in ks_activity_mat; fall back gracefully
    # if for some reason ks_activity_type is missing (old result object)
    if (is.null(ks_activity_type))
      ks_activity_type <- setNames(
        rep("PSP confirmed", nrow(ks_activity_mat)), rownames(ks_activity_mat))

    all_in_mat <- intersect(kin_summary$Kinase, rownames(ks_activity_mat))

    # When a user kinase filter was applied (cluster_cols == FALSE signals this),
    # restrict the activity heatmap to the same kinase set used in the KS heatmap.
    user_kin_filter <- !isTRUE(res$cluster_cols)
    if (user_kin_filter) {
      act_kinases <- intersect(res$selected_kinases, rownames(ks_activity_mat))
    } else if (show_all) {
      act_kinases <- all_in_mat
    } else {
      act_kinases <- head(all_in_mat, 40L)
    }
    req(length(act_kinases) > 0L)

    act_mat   <- ks_activity_mat[act_kinases, , drop = FALSE]
    act_types <- ks_activity_type[act_kinases]

    # Drop all-NA rows (kinases whose proxy could not be computed, e.g. mouse
    # kinases with no PWM entry). pheatmap's hclust cannot handle all-NA rows.
    valid_rows <- rowSums(!is.na(act_mat)) > 0L
    act_mat    <- act_mat[valid_rows,  , drop = FALSE]
    act_types  <- act_types[valid_rows]
    act_kinases <- rownames(act_mat)

    # Apply custom sample order when provided
    if (!is.null(sample_order) && length(sample_order) > 0L) {
      valid_cols   <- sample_order[sample_order %in% colnames(act_mat)]
      leftover     <- setdiff(colnames(act_mat), valid_cols)   # any not in list
      act_mat      <- act_mat[, c(valid_cols, leftover), drop = FALSE]
      cluster_cols_val <- FALSE
    } else {
      cluster_cols_val <- ncol(act_mat) > 1L
    }

    n_psp   <- sum(act_types == "PSP confirmed",   na.rm = TRUE)
    n_proxy <- sum(act_types == "Proxy (top sites)", na.rm = TRUE)
    title_str <- paste0(
      "Kinase Activity \u2014 Z-score profile (",
      if (user_kin_filter) paste0(nrow(act_mat), " selected kinases")
      else if (show_all) paste0(nrow(act_mat), " kinases") else "top 40",
      ": ", n_psp, " PSP-confirmed, ", n_proxy, " proxy",
      if (!is.null(sample_order) && length(sample_order) > 0L) ", custom sample order" else "",
      ")")

    kfam_df <- as.data.frame(.kinase_family(), stringsAsFactors = FALSE)
    ann_row <- data.frame(
      Source = factor(act_types,
                      levels = c("PSP confirmed", "Proxy (top sites)")),
      Group  = kfam_df[match(act_kinases, rownames(kfam_df)), "kinase_group"],
      Family = kfam_df[match(act_kinases, rownames(kfam_df)), "kinase_family"],
      stringsAsFactors = TRUE
    )
    rownames(ann_row) <- act_kinases
    ann_row[is.na(ann_row$Group),  "Group"]  <- "Unknown"
    ann_row[is.na(ann_row$Family), "Family"] <- "Unknown"
    ann_colors <- list(
      Source = c("PSP confirmed"    = "#E69F00",
                 "Proxy (top sites)" = "#BBBBBB")
    )

    pheatmap::pheatmap(
      mat               = act_mat,
      annotation_row    = ann_row,
      annotation_colors = ann_colors,
      cluster_rows      = nrow(act_mat) > 1L,
      cluster_cols      = cluster_cols_val,
      clustering_method = "ward.D2",
      cellwidth         = 12,
      cellheight        = 9,
      border_color      = "grey85",
      fontsize_row      = 7,
      fontsize_col      = 7,
      main              = title_str
    )
  }

  # ---------------------------------------------------------------------------
  # Dynamic height helper
  # ---------------------------------------------------------------------------
  .plot_height <- function(n_sites, base = 120, per_row = 10, min_h = 400, max_h = 2000) {
    max(min_h, min(max_h, base + n_sites * per_row))
  }

  # ---------------------------------------------------------------------------
  # renderPlot outputs
  # ---------------------------------------------------------------------------
  output$ks_heatmap_up <- renderPlot({
    req(!is.null(ks_rv$result_up))
    .render_ks_heatmap(ks_rv$result_up)
  }, height = function() {
    res <- ks_rv$result_up
    if (is.null(res)) return(600)
    .plot_height(length(res$selected_sites) + length(res$missing_sites))
  })

  output$ks_heatmap_down <- renderPlot({
    req(!is.null(ks_rv$result_down))
    .render_ks_heatmap(ks_rv$result_down)
  }, height = function() {
    res <- ks_rv$result_down
    if (is.null(res)) return(600)
    .plot_height(length(res$selected_sites))
  })

  # Helper: pick the activity result to display based on current mode
  .activity_res <- function() {
    if (ks_rv$activity_mode == "full" && !is.null(ks_rv$result_activity_full))
      ks_rv$result_activity_full
    else if (!is.null(ks_rv$result_up))
      ks_rv$result_up
    else
      ks_rv$result_down
  }

  output$ks_activity_heatmap <- renderPlot({
    res <- .activity_res()
    req(!is.null(res))
    custom_order <- if (isTRUE(input$ks_custom_sample_order)) input$ks_sample_order else NULL
    .render_activity_heatmap(res,
                             show_all     = isTRUE(input$ks_show_all_kinases),
                             sample_order = custom_order)
  }, height = function() {
    res <- .activity_res()
    if (is.null(res)) return(540)
    n_kin <- if (isTRUE(input$ks_show_all_kinases))
      nrow(res$ks_activity_mat)
    else
      min(40L, nrow(res$ks_activity_mat))
    max(400, min(2000, 120 + n_kin * 14))
  })

  # ---------------------------------------------------------------------------
  # Status text
  # ---------------------------------------------------------------------------
  output$ks_status_text <- renderText({ ks_rv$status })

  # ---------------------------------------------------------------------------
  # Visibility flags (used by conditionalPanel in ui.R)
  # NOTE: must use renderText(), not reactive() — only render functions send
  # values to the JS client for conditionalPanel evaluation.
  # ---------------------------------------------------------------------------
  output$ks_has_result_up   <- renderText({ if (!is.null(ks_rv$result_up))   "yes" else "" })
  output$ks_has_result_down <- renderText({ if (!is.null(ks_rv$result_down))  "yes" else "" })
  output$ks_activity_mode   <- renderText({ ks_rv$activity_mode })
  output$ks_has_plex        <- renderText({ if (!is.null(ks_rv$plex_labels))  "yes" else "" })

  outputOptions(output, "ks_has_result_up",   suspendWhenHidden = FALSE)
  outputOptions(output, "ks_has_result_down",  suspendWhenHidden = FALSE)
  outputOptions(output, "ks_activity_mode",    suspendWhenHidden = FALSE)
  outputOptions(output, "ks_has_plex",         suspendWhenHidden = FALSE)

  # ---------------------------------------------------------------------------
  # Download handlers — PDF with dynamic dimensions
  # ---------------------------------------------------------------------------
  # Inches per cell for pheatmap output (cellwidth/cellheight = 9pt = 9/72 in)
  .pdf_dims_heatmap <- function(n_sites, n_kinases) {
    cell_in <- 9 / 72
    w <- max(10, n_kinases * cell_in + 5)   # +5 for row labels + annotations
    h <- max(8,  n_sites   * cell_in + 4)   # +4 for col labels + annotations
    list(width = w, height = h)
  }
  .pdf_dims_activity <- function(n_kinases, n_samples) {
    cell_in_r <- 9  / 72   # cellheight for kinase rows
    cell_in_c <- 12 / 72   # cellwidth  for sample cols
    w <- max(8, n_samples  * cell_in_c + 4)
    h <- max(8, n_kinases  * cell_in_r + 3)
    list(width = w, height = h)
  }

  .pdf_download <- function(render_fn, name_stem, dims_fn = NULL) {
    downloadHandler(
      filename = function() {
        paste0(name_stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        tmp <- tempfile(fileext = ".pdf")
        dims <- if (!is.null(dims_fn)) dims_fn() else list(width = 18, height = 10)
        tryCatch({
          grDevices::pdf(tmp, width = dims$width, height = dims$height)
          render_fn()
          grDevices::dev.off()
          file.copy(tmp, file, overwrite = TRUE)
        }, error = function(e) {
          try(grDevices::dev.off(), silent = TRUE)
          writeLines(paste("PDF generation failed:", conditionMessage(e)), file)
        })
      },
      contentType = "application/pdf"
    )
  }

  output$ks_pdf_up <- .pdf_download(
    function() .render_ks_heatmap(ks_rv$result_up),
    "ks_heatmap_up",
    dims_fn = function() {
      res <- ks_rv$result_up
      if (is.null(res)) return(list(width = 18, height = 10))
      .pdf_dims_heatmap(
        length(res$selected_sites) + length(res$missing_sites),
        length(res$selected_kinases)
      )
    }
  )

  output$ks_pdf_down <- .pdf_download(
    function() .render_ks_heatmap(ks_rv$result_down),
    "ks_heatmap_down",
    dims_fn = function() {
      res <- ks_rv$result_down
      if (is.null(res)) return(list(width = 18, height = 10))
      .pdf_dims_heatmap(
        length(res$selected_sites),
        length(res$selected_kinases)
      )
    }
  )

  output$ks_pdf_activity <- .pdf_download(
    function() {
      res          <- .activity_res()
      custom_order <- if (isTRUE(input$ks_custom_sample_order)) input$ks_sample_order else NULL
      .render_activity_heatmap(res,
                               show_all     = isTRUE(input$ks_show_all_kinases),
                               sample_order = custom_order)
    },
    "ks_kinase_activity",
    dims_fn = function() {
      res <- .activity_res()
      if (is.null(res)) return(list(width = 14, height = 10))
      n_kin <- if (isTRUE(input$ks_show_all_kinases))
        nrow(res$ks_activity_mat)
      else
        min(40L, nrow(res$ks_activity_mat))
      .pdf_dims_activity(n_kin, ncol(res$ks_activity_mat))
    }
  )

  # ---------------------------------------------------------------------------
  # Kinase ranking plot (ggplot2 lollipop) — all kinases from score_mat,
  # ranked by mean score, colored by kinase group.  Works at any scale.
  # ---------------------------------------------------------------------------
  .render_kinase_ranking <- function(res) {
    req(!is.null(res))
    score_mat <- res$score_mat
    req(ncol(score_mat) > 0L)

    kin_means <- sort(colMeans(score_mat, na.rm = TRUE), decreasing = TRUE,
                      na.last = NA)
    kfam_df   <- as.data.frame(.kinase_family(), stringsAsFactors = FALSE)

    df_plot <- data.frame(
      Kinase     = names(kin_means),
      Mean_Score = as.numeric(kin_means),
      stringsAsFactors = FALSE
    )
    df_plot$Group <- kfam_df[match(df_plot$Kinase, rownames(kfam_df)), "kinase_group"]
    df_plot$Group[is.na(df_plot$Group)] <- "Other"
    df_plot$Kinase <- factor(df_plot$Kinase, levels = rev(df_plot$Kinase))

    # Group colour palette (up to ~10 groups)
    grp_lvls <- unique(as.character(df_plot$Group))
    pal <- setNames(
      c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
        "#D55E00","#CC79A7","#999999","#44AA99","#882255","#DDDDDD"),
      c(grp_lvls, rep("Other", max(0L, 11L - length(grp_lvls))))[seq_along(grp_lvls)]
    )

    ggplot2::ggplot(df_plot,
        ggplot2::aes(x = Mean_Score, y = Kinase, colour = Group)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = Mean_Score, y = Kinase, yend = Kinase),
        linewidth = 0.4, alpha = 0.6) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = pal, name = "Kinase Group") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.04))) +
      ggplot2::labs(
        x     = "Mean KS Score",
        y     = NULL,
        title = paste0("Kinase Ranking \u2014 ",
                       nrow(df_plot), " kinases, ",
                       length(res$selected_sites), " sites")
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        axis.text.y     = ggplot2::element_text(size = 6),
        legend.position = "right",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank()
      )
  }

  output$ks_ranking_plot <- renderPlot({
    res <- if (!is.null(ks_rv$result_up)) ks_rv$result_up else ks_rv$result_down
    req(!is.null(res))
    .render_kinase_ranking(res)
  }, height = function() {
    res <- if (!is.null(ks_rv$result_up)) ks_rv$result_up else ks_rv$result_down
    if (is.null(res)) return(600)
    n_scored <- sum(!is.nan(colMeans(res$score_mat, na.rm = TRUE)))
    max(500, min(4000, n_scored * 9 + 120))
  })

  output$ks_pdf_ranking <- .pdf_download(
    function() {
      res <- if (!is.null(ks_rv$result_up)) ks_rv$result_up else ks_rv$result_down
      print(.render_kinase_ranking(res))
    },
    "ks_kinase_ranking",
    dims_fn = function() {
      res <- if (!is.null(ks_rv$result_up)) ks_rv$result_up else ks_rv$result_down
      if (is.null(res)) return(list(width = 8, height = 12))
      n_kin <- sum(!is.nan(colMeans(res$score_mat, na.rm = TRUE)))
      list(width = 8, height = max(6, n_kin * 0.14 + 1.5))
    }
  )
}
