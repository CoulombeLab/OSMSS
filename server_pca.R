# server_pca.R ==================================================================
# Server logic for Proteome Discoverer > PCA Analysis module.
# Two processing pipelines:
#   1. Global (Whole Tissue)    - equalizeMedians → TMP → AFT imputation → PCA
#   2. IP (Immunoprecipitation) - invariant-protein median normalization (ProDA) → PCA
#
# Both use peptide-level summarization driven by the Upload Wizard's condition
# groups; PCA itself consumes rv$lf_master_cache.
#
# Output: PC1 vs PC2 scatter plot + scree plot + CV analysis.
# Called from server.R as: server_pca(input, output, session, rv, logger)
# ==============================================================================

server_pca <- function(input, output, session, rv, logger) {

  # =========================================================================
  # Reactive storage
  # =========================================================================
  pca_rv <- reactiveValues(
    # Pipeline mode
    pipeline         = "global",  # "global" | "ip"
    # Condition groups (Global / IP modes)
    cond_groups      = list(),    # named list: cond_name -> c(col1, col2, ...)
    # Results
    pca_result       = NULL,      # prcomp() output
    pca_scores       = NULL,      # data.frame of PC scores + Sample + Group
    pca_mat          = NULL,      # transposed matrix used for PCA (rows=samples, cols=proteins)
    cv_groups        = NULL,      # named list: group_name -> vector of sample row names
    # Reference-panel state
    pending          = NULL,      # list(...) saved state for the IgG modal
    ebn_refs         = NULL       # character vector of invariant reference protein accessions
  )
  register_experiment_reset(rv, pca_rv)
  methodsTextServer("pca_methods", provenance = rv$provenance, section = "pca")

  # =========================================================================
  # SUMMARIZATION CACHE — Summarize / Export / Import
  # =========================================================================
  output$pca_export_summarization <- downloadHandler(
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
      logger("PCA: Summarization exported")
    }
  )

  observe({
    if (is.null(rv$lf_master_cache)) {
      shinyjs::disable("pca_export_summarization")
    } else {
      shinyjs::enable("pca_export_summarization")
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

  # =========================================================================
  # Reset state when new CSV is loaded
  # =========================================================================
  observeEvent(rv$pd_protein_headers, {
    req(rv$pd_protein_headers)
    pca_rv$pca_result  <- NULL
    pca_rv$pca_scores  <- NULL
    pca_rv$pca_mat     <- NULL
    pca_rv$cv_groups   <- NULL
    pca_rv$pending     <- NULL
    pca_rv$ebn_refs    <- NULL
    logger("PCA: Shared CSV data updated, state reset")
  })

  # =========================================================================
  # CV Group Assignment (all modes)
  # =========================================================================
  observeEvent(input$pca_cv_assign_groups, {
    req(pca_rv$pca_scores)
    scores_df <- pca_rv$pca_scores

    # If groups already assigned via MSstats pipeline, pre-fill from those
    current_groups <- pca_rv$cv_groups
    sample_names <- scores_df$Sample

    col_inputs <- lapply(seq_along(sample_names), function(i) {
      sname <- sample_names[i]
      # Default: use Group column if available, else derive from sample name
      default_group <- if (!is.null(scores_df$Group) && !is.na(scores_df$Group[i])) {
        scores_df$Group[i]
      } else if (!is.null(current_groups)) {
        grp <- names(which(sapply(current_groups, function(g) sname %in% g)))
        if (length(grp) > 0) grp[1] else paste0("Group", ((i - 1) %% 2) + 1)
      } else {
        paste0("Group", ((i - 1) %% 2) + 1)
      }

      fluidRow(
        column(7, tags$p(class = "small mb-0 mt-2", tags$code(sname))),
        column(5, textInput(paste0("pca_cv_grp_", i), label = NULL,
                            value = default_group, width = "100%"))
      )
    })

    showModal(modalDialog(
      title = tags$span(icon("layer-group"), " Assign Samples to Groups"),
      tagList(
        helpText(class = "small",
          "Enter a group name for each sample. ",
          "CV will be calculated within each group."),
        hr(),
        tagList(col_inputs)
      ),
      size = "m", easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("pca_cv_grp_confirm", "Confirm Groups", class = "btn-success")
      )
    ))
  })

  observeEvent(input$pca_cv_grp_confirm, {
    req(pca_rv$pca_scores)
    sample_names <- pca_rv$pca_scores$Sample

    group_map <- list()
    for (i in seq_along(sample_names)) {
      grp_name <- trimws(input[[paste0("pca_cv_grp_", i)]])
      if (nchar(grp_name) == 0) grp_name <- paste0("Group", i)
      if (is.null(group_map[[grp_name]])) group_map[[grp_name]] <- c()
      group_map[[grp_name]] <- c(group_map[[grp_name]], sample_names[i])
    }

    pca_rv$cv_groups <- group_map
    logger(paste0("PCA CV: Assigned ", length(group_map), " groups: ",
                  paste(names(group_map), collapse = ", ")))
    removeModal()
  })

  # CV group assignment label
  output$pca_cv_group_label <- renderUI({
    grps <- pca_rv$cv_groups
    if (is.null(grps) || length(grps) == 0) {
      tags$p(class = "text-muted small", "No groups assigned yet.")
    } else {
      tagList(
        lapply(names(grps), function(g) {
          tags$p(class = "small mb-1",
                 tags$span(class = "badge bg-secondary", g),
                 paste0(" (", length(grps[[g]]), " samples)"))
        })
      )
    }
  })

  # =========================================================================
  # RUN PCA - Dispatcher
  # =========================================================================
  .pca_dispatch <- function() {
    pipeline <- isolate(rv$lf_pipeline)
    if (is.null(pipeline)) pipeline <- "global"
    pca_rv$pipeline <- pipeline

    # Reset previous results
    pca_rv$pca_result <- NULL
    pca_rv$pca_scores <- NULL
    pca_rv$pca_mat    <- NULL
    pca_rv$cv_groups  <- NULL
    pca_rv$pending    <- NULL

    cache <- rv$lf_master_cache
    meta  <- rv$lf_master_cache_meta
    if (!is.null(cache) && !is.null(meta) &&
        identical(meta$pipeline, pipeline) &&
        (pipeline != "ip" || identical(meta$ip_sub, isolate(rv$lf_ip_sub)))) {
      logger(paste0("PCA [", pipeline, "]: Using cached summarization"))
      .run_pca_from_cache(cache)
    } else {
      showNotification(
        "No summarized data found. Please run the Upload Wizard to summarize your data before running PCA.",
        type = "warning", duration = 8
      )
      return(invisible(NULL))
    }
  }

  # =========================================================================
  # TMT PCA (unified route — called when rv$data_source == "tmt")
  # =========================================================================
  .run_pca_tmt <- function() {
    is_limma_route <- isTRUE(rv$tmt_analysis_route %in% c("limma_protein", "limma_qf"))
    if (is_limma_route && is.null(rv$tmt_limma_cache)) {
      showNotification("Please run TMT Summarize Data first (Upload Data tab).",
                       type = "warning")
      return()
    }
    if (!is_limma_route && is.null(rv$tmt_master_cache)) {
      showNotification("Please run TMT Summarize Data first (Upload Data tab).",
                       type = "warning")
      return()
    }

    # ── Limma-route PCA (protein_matrix directly from cache) ──────────────
    if (is_limma_route) {
      lc  <- rv$tmt_limma_cache
      ann <- lc$annotation_map
      conditions <- lc$conditions
      if (is.null(conditions) || length(conditions) < 2) {
        showNotification("Need at least 2 conditions.", type = "warning"); return()
      }
      shinyjs::disable("pca_run_btn")
      on.exit(shinyjs::enable("pca_run_btn"), add = TRUE)
      withProgress(message = "Running TMT PCA (Limma)...", value = 0, {
        incProgress(0.3, detail = "Building matrix...")
        pm <- lc$protein_matrix   # proteins × channels
        # exclude Norm (bridge) channels
        use_chans <- ann$ColumnName[!ann$Condition %in% "Norm"]
        use_chans <- intersect(use_chans, colnames(pm))
        mat <- t(pm[, use_chans, drop = FALSE])  # samples × proteins
        sample_groups <- ann$Condition[match(rownames(mat), ann$ColumnName)]
        valid_cols <- colSums(is.na(mat)) == 0
        mat <- mat[, valid_cols, drop = FALSE]
        if (ncol(mat) < 2 || nrow(mat) < 3) {
          showNotification("Too few valid proteins/samples.", type = "error"); return()
        }
        incProgress(0.5, detail = "Computing PCA...")

        # ── Multi-plex batch correction for visualization ────────────────────
        # For designs with ≥2 plexes, save an uncorrected PCA (→ pca_result_de /
        # pca_scores_de shown in the "Uncorrected PCA" tab), then apply
        # limma::removeBatchEffect with Mixture as batch factor for the main tab.
        mixture_vec <- ann$Mixture[match(rownames(mat), ann$ColumnName)]
        n_plexes    <- length(unique(mixture_vec[!is.na(mixture_vec)]))
        did_bc_tmt  <- FALSE
        # See tmt_plex_correctable(): with the Norm channels dropped, a plex holding
        # one sample cannot have its offset separated from that sample's own value.
        plex_ok_pca <- tmt_plex_correctable(mixture_vec)
        if (n_plexes >= 2 && !plex_ok_pca)
          logger(paste0("PCA [TMT limma]: plex batch correction skipped  -  at least one ",
                        "plex holds a single non-reference channel (correction would set ",
                        "each such sample to its condition mean)."))

        if (n_plexes >= 2 && plex_ok_pca && requireNamespace("limma", quietly = TRUE)) {
          # 1. Save uncorrected PCA
          pca_raw          <- prcomp(mat, center = TRUE, scale. = FALSE)
          scores_raw       <- as.data.frame(pca_raw$x)
          scores_raw$Sample <- rownames(scores_raw)
          scores_raw$Group  <- sample_groups
          pca_rv$pca_result_de <- pca_raw
          pca_rv$pca_scores_de <- scores_raw

          # 2. Batch-correct (proteins × samples, then transpose back)
          # The flag is set on the SUCCESS path only. It used to be set
          # unconditionally after the tryCatch, so a caught failure still recorded
          # batch_corrected = TRUE in provenance and the methods text asserted a
          # correction that never happened.
          mat_t_bc <- tryCatch({
            .pca_bc <- limma::removeBatchEffect(
              t(mat),
              batch  = factor(mixture_vec),
              design = model.matrix(~ factor(sample_groups))
            )
            did_bc_tmt <- TRUE
            .pca_bc
          },
            error = function(e) {
              logger(paste0("PCA [TMT limma]: removeBatchEffect failed: ", e$message,
                            "  -  PCA shown WITHOUT plex batch correction."))
              did_bc_tmt <<- FALSE
              t(mat)
            }
          )
          mat        <- t(mat_t_bc)   # back to samples × proteins
          if (isTRUE(did_bc_tmt))
            logger(paste0("PCA [TMT limma]: plex batch correction applied (",
                          n_plexes, " plexes)."))
        } else {
          pca_rv$pca_result_de <- NULL
          pca_rv$pca_scores_de <- NULL
        }

        pca_out   <- prcomp(mat, center = TRUE, scale. = FALSE)
        scores_df <- as.data.frame(pca_out$x)
        scores_df$Sample <- rownames(scores_df)
        scores_df$Group  <- sample_groups
        pca_rv$pca_result <- pca_out
        pca_rv$pca_scores <- scores_df
        pca_rv$pca_mat    <- mat
        pca_rv$cv_groups  <- split(rownames(mat), sample_groups)
        .pca_var <- summary(pca_out)$importance[2, ] * 100
        rv$provenance$pca <- list(
          source          = "tmt",
          pipeline        = "tmt",
          pipeline_label  = "TMT quantitative proteomics pipeline",
          route           = rv$tmt_analysis_route,
          n_samples       = as.numeric(nrow(mat)),
          n_proteins      = as.numeric(ncol(mat)),
          centering       = TRUE,
          scaling         = FALSE,
          batch_corrected = did_bc_tmt,
          pc1_var         = as.numeric(.pca_var[1]),
          pc2_var         = as.numeric(.pca_var[min(2, length(.pca_var))]),
          condition_names = as.character(unique(scores_df$Group)),
          software_used   = "prcomp (base R)"
        )
        incProgress(0.2)
      })
      logger(paste0("PCA [TMT limma]: Done. ", nrow(mat), " samples, ", ncol(mat), " proteins."))
      return()
    }

    conditions <- rv$tmt_conditions
    if (is.null(conditions) || length(conditions) < 2) {
      showNotification("Need at least 2 conditions.", type = "warning")
      return()
    }

    logger("PCA [TMT]: Starting pipeline (scale.=FALSE, using cached summarization)")

    shinyjs::disable("pca_run_btn")
    on.exit(shinyjs::enable("pca_run_btn"), add = TRUE)

    processed <- rv$tmt_master_cache

    withProgress(message = "Running TMT PCA...", value = 0, {
      incProgress(0.50, detail = "Building protein \u00d7 sample matrix...")
      pld <- processed$ProteinLevelData
      if (is.null(pld) || nrow(pld) == 0) {
        showNotification("No protein data in TMT cache.", type = "error"); return()
      }
      pld <- pld[!is.na(pld$Abundance), ]
      pld$Protein      <- as.character(pld$Protein)
      pld$Condition    <- as.character(pld$Condition)
      pld$BioReplicate <- as.character(pld$BioReplicate)
      pld <- pld[pld$Condition != "Norm", ]
      if (nrow(pld) == 0) { showNotification("No data after removing Norm.", type = "error"); return() }

      wide <- reshape(
        pld[, c("Protein", "Abundance", "BioReplicate")],
        idvar = "BioReplicate", timevar = "Protein", direction = "wide"
      )
      colnames(wide) <- sub("^Abundance\\.", "", colnames(wide))
      rownames(wide) <- wide$BioReplicate
      wide$BioReplicate <- NULL

      bio_cond_map <- unique(pld[, c("BioReplicate", "Condition")])
      rownames(bio_cond_map) <- bio_cond_map$BioReplicate
      sample_groups <- bio_cond_map[rownames(wide), "Condition"]

      mat <- as.matrix(wide)
      valid_cols <- colSums(is.na(mat)) == 0
      mat <- mat[, valid_cols, drop = FALSE]

      if (ncol(mat) < 2 || nrow(mat) < 3) {
        showNotification("Too few valid proteins/samples.", type = "error"); return()
      }

      # Batch correction (unbridged multi-plex)
      design_opt_pca <- rv$tmt_master_cache_design
      ann_pca        <- rv$tmt_master_cache_ann
      use_ref_norm_pca <- identical(as.integer(design_opt_pca), 2L)
      did_batch_pca <- FALSE
      if (!use_ref_norm_pca && !is.null(ann_pca)) {
        bio_ann_pca  <- ann_pca[!ann_pca$Condition %in% "Norm", , drop = FALSE]
        n_plexes_pca <- length(unique(bio_ann_pca$Run))
        if (n_plexes_pca > 1) {
          bio_run_map  <- unique(pld[, c("BioReplicate", "Run"), drop = FALSE])
          sample_runs  <- bio_run_map$Run[match(rownames(mat), bio_run_map$BioReplicate)]
          if (!any(is.na(sample_runs)) && length(unique(sample_runs)) > 1 &&
              !tmt_plex_correctable(sample_runs)) {
            logger(paste0("PCA [TMT]: plex batch correction skipped  -  at least one ",
                          "plex holds a single non-reference channel (correction would ",
                          "set each such sample to its condition mean)."))
          } else if (!any(is.na(sample_runs)) && length(unique(sample_runs)) > 1) {
            mat_t <- t(mat)
            mat_t <- limma::removeBatchEffect(
              mat_t,
              batch  = factor(sample_runs),
              design = model.matrix(~ factor(sample_groups))
            )
            mat <- t(mat_t)
            did_batch_pca <- TRUE
          }
        }
      }

      incProgress(0.70, detail = "Computing PCA...")
      pca_out   <- prcomp(mat, center = TRUE, scale. = FALSE)
      scores_df <- as.data.frame(pca_out$x)
      scores_df$Sample <- rownames(scores_df)
      scores_df$Group  <- sample_groups

      pca_rv$pca_result <- pca_out
      pca_rv$pca_scores <- scores_df
      pca_rv$pca_mat    <- mat
      pca_rv$cv_groups  <- split(rownames(mat), sample_groups)

      .pca_var <- summary(pca_out)$importance[2, ] * 100
      rv$provenance$pca <- list(
        source           = "tmt",
        pipeline         = "tmt",
        pipeline_label   = "TMT quantitative proteomics pipeline",
        n_samples        = as.numeric(nrow(mat)),
        n_proteins       = as.numeric(ncol(mat)),
        design_option    = as.numeric(design_opt_pca),
        global_norm      = TRUE,
        reference_norm   = use_ref_norm_pca,
        mb_impute        = FALSE,
        batch_corrected  = did_batch_pca,
        centering        = TRUE,
        scaling          = FALSE,
        pc1_var          = as.numeric(.pca_var[1]),
        pc2_var          = as.numeric(.pca_var[2]),
        condition_names  = as.character(unique(scores_df$Group)),
        software_used    = "prcomp (base R)"
      )
      incProgress(1, detail = "Done!")
      logger(paste0("PCA [TMT]: Done! ", nrow(scores_df), " samples, ",
                    ncol(mat), " proteins"))
    })
  }

  # =========================================================================
  # DIA PCA (unified route — called when rv$data_source == "dia")
  # =========================================================================
  .run_pca_dia <- function() {
    if (is.null(rv$dia_master_cache)) {
      showNotification("Please complete the DIA Upload Wizard first.",
                       type = "warning")
      return()
    }

    shinyjs::disable("pca_run_btn")
    on.exit(shinyjs::enable("pca_run_btn"), add = TRUE)

    mat  <- rv$dia_master_cache$protein_matrix
    grps <- rv$dia_cond_groups
    cond_map <- setNames(
      rep(names(grps), vapply(grps, length, integer(1))),
      unlist(grps)
    )
    use_samps <- intersect(colnames(mat), unlist(grps))
    mat <- mat[, use_samps, drop = FALSE]

    withProgress(message = "Running DIA PCA...", value = 0.2, {
      # Row-median impute NAs
      mat <- t(apply(mat, 1, function(r) {
        med <- median(r, na.rm = TRUE)
        r[is.na(r)] <- med
        r
      }))
      mat <- mat[complete.cases(mat), , drop = FALSE]
      n_prots <- nrow(mat)

      # Optional batch correction
      if (isTRUE(rv$dia_batch_corrected) && !is.null(rv$dia_batch_groups) &&
          requireNamespace("limma", quietly = TRUE)) {
        bat_map <- setNames(
          rep(names(rv$dia_batch_groups),
              vapply(rv$dia_batch_groups, length, integer(1))),
          unlist(rv$dia_batch_groups)
        )
        batch_vec <- bat_map[colnames(mat)]
        mat <- tryCatch(
          limma::removeBatchEffect(mat, batch = batch_vec),
          error = function(e) mat
        )
      }
      incProgress(0.5)

      pca_res <- prcomp(t(mat), center = TRUE, scale. = FALSE)
      var_exp <- summary(pca_res)$importance[2, ] * 100

      pca_df <- as.data.frame(pca_res$x[, 1:min(2, ncol(pca_res$x))])
      pca_df$Sample <- rownames(pca_df)
      pca_df$Group  <- cond_map[pca_df$Sample]   # rename Condition -> Group for LF DDA compat

      # t(mat): rows=samples, cols=proteins — expected format for pca_rv$pca_mat
      mat_sp <- t(mat)

      pca_rv$pca_result <- pca_res
      pca_rv$pca_scores <- pca_df
      pca_rv$pca_mat    <- mat_sp
      pca_rv$cv_groups  <- lapply(grps, function(s) intersect(s, colnames(mat)))

      rv$provenance$pca <- list(
        source          = "dia",
        pipeline        = "dia",
        pipeline_label  = "DIA proteomics pipeline",
        n_proteins      = n_prots,
        n_samples       = ncol(mat),
        batch_corrected = isTRUE(rv$dia_batch_corrected),
        centering       = TRUE,
        scaling         = FALSE,
        pc1_var         = round(var_exp[1], 1),
        pc2_var         = round(var_exp[min(2, length(var_exp))], 1),
        condition_names = names(grps),
        software_used   = "prcomp (base R)"
      )
      incProgress(0.3)
    })
    logger(paste0("[PCA DIA] Completed: ", n_prots, " proteins, ",
                  length(use_samps), " samples."))
  }

  # =========================================================================
  # RUN PCA button handler — source router + LF DDA path
  # =========================================================================
  observeEvent(input$pca_run_btn, {
    # Route to TMT or DIA path if appropriate data source is active
    src <- rv$data_source
    if (identical(src, "tmt")) {
      .run_pca_tmt()
      return()
    } else if (identical(src, "dia")) {
      .run_pca_dia()
      return()
    }

    # LF DDA path (unchanged)
    if (is.null(rv$lf_master_cache) && is.null(rv$pd_raw_lines)) {
      showNotification("Please run the Upload Wizard first to summarize data.", type = "warning")
      return()
    }

    shinyjs::disable("pca_run_btn")
    on.exit(shinyjs::enable("pca_run_btn"), add = TRUE)

    pipeline <- isolate(rv$lf_pipeline)
    if (is.null(pipeline)) pipeline <- "global"

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
            "running PCA. This may take a few minutes."
          ),
          tags$p(class = "text-muted small mb-0",
            "Pipeline: ", tags$b(pipeline),
            if (identical(pipeline, "ip")) paste0(" (", isolate(rv$lf_ip_sub), ")") else "")
        ),
        size = "m", easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("pca_confirm_resummarize", "Continue",
                       class = "btn-primary", icon = icon("play"))
        )
      ))
      return()
    }

    .pca_dispatch()
  })

  observeEvent(input$pca_confirm_resummarize, {
    removeModal()
    rv$lf_cache_stale <- FALSE
    .pca_dispatch()
  })

  # =========================================================================
  #  PIPELINE CACHE: Run PCA from pre-built lf_master_cache
  # =========================================================================
  .run_pca_from_cache <- function(cache) {
    # Single-trunk architecture: always use cache$processed.
    # When batch correction is active, apply protein-level removeBatchEffect
    # for PCA visualization only. DE results are unaffected.
    processed_vis <- cache$processed
    acc_to_gene   <- cache$acc_to_gene
    pld           <- processed_vis$ProteinLevelData
    grps          <- rv$lf_cond_groups

    withProgress(message = "PCA \u2014 From Cache...", value = 0, {
      incProgress(0.2, detail = "Pivoting protein-level data...")

      # Pivot ProteinLevelData to wide matrix (rows = proteins, cols = runs)
      pld_wide <- tryCatch(
        tidyr::pivot_wider(pld[, c("Protein", "RUN", "LogIntensities")],
                           names_from = "RUN", values_from = "LogIntensities"),
        error = function(e) NULL
      )
      if (is.null(pld_wide)) {
        showNotification("Cache pivot failed — try re-summarizing.", type = "error")
        return()
      }

      proteins_v <- pld_wide$Protein
      mat <- as.matrix(pld_wide[, -1])
      rownames(mat) <- proteins_v
      # Row count as pivoted, before complete-case selection and before the
      # column-median fill below. provenance$n_prots_input used to be read off
      # mat at the end, which is a post-filter count under an "input" name.
      .n_prots_raw <- nrow(mat)

      # Build run-to-label map from ProteinLevelData metadata
      run_label_map <- unique(pld[, c("RUN", "GROUP", "SUBJECT")])
      run_to_label <- setNames(as.character(run_label_map$SUBJECT),
                               as.character(run_label_map$RUN))
      run_to_group <- setNames(as.character(run_label_map$GROUP),
                               as.character(run_label_map$RUN))

      incProgress(0.4, detail = "Filtering & imputing...")
      # IP/AP-MS: presence/absence IS the result, so select complete cases
      # instead of filling gaps. Gated on the experiment type (rv$lf_pipeline),
      # not on SAINT identifiers: the ProDA IP route sets a ProDA pipeline_label
      # and de_engine "proda", so the old string test never
      # fired for it and the route silently fell through to imputation. The
      # original SAINT test is OR-ed in so that route cannot regress. Other
      # routes (global LF) reach the block unchanged.
      .is_ip_route_pca <- identical(isolate(rv$lf_pipeline) %||% "", "ip") ||
        identical(cache$pipeline_label, "IP (SAINTexpress)")
      .n_prots_complete_case <- NA_integer_
      if (.is_ip_route_pca) {
        mat <- mat[stats::complete.cases(mat), , drop = FALSE]
        .n_prots_complete_case <- nrow(mat)
        logger(paste0("PCA [Cache]: IP complete-case selection \u2014 ",
                      nrow(mat), " proteins (no imputation)"))
      }
      # Column-median (per-sample) imputation for remaining NAs. This loop is
      # UNCONDITIONAL: it runs on every route. On IP routes it is a no-op only
      # because complete.cases() above already removed every NA. It is not
      # skipped, and it fills whatever any earlier missing-value handling left.
      # Cells fabricated by the fill below. Unlike the heatmap path, valid_r /
      # mat_clean RETAINS these rows, so they reach prcomp as well as the batch
      # model. Counted before the loop mutates mat.
      .n_prebatch_imputed_pca <- as.integer(sum(is.na(mat)))
      for (ci in seq_len(ncol(mat))) {
        na_m  <- is.na(mat[, ci])
        if (any(na_m)) mat[na_m, ci] <- median(mat[!na_m, ci], na.rm = TRUE)
      }
      if (.n_prebatch_imputed_pca > 0L)
        logger(paste0("PCA [Cache]: column-median fill replaced ",
                      .n_prebatch_imputed_pca, " missing cell(s); these values",
                      " reach prcomp and, when enabled, the batch model"))
      valid_r <- rowSums(is.na(mat)) == 0 & rowSums(is.finite(mat)) == ncol(mat)
      mat_clean <- mat[valid_r, , drop = FALSE]

      if (nrow(mat_clean) < 3) {
        showNotification("Too few valid proteins after filtering.", type = "error")
        return()
      }

      # Batch correction for visualization (limma::removeBatchEffect)
      batch_applied <- FALSE
      bg <- rv$lf_batch_groups
      if (isTRUE(rv$lf_batch_corrected) && !is.null(bg) && length(bg) >= 2) {
        require_pkg("limma")
        cn_clean <- colnames(mat_clean)

        # Remap batch groups: the user assigned original PD column names, but the
        # ProDA PLD uses synthesised sample_name keys (e.g. "Cond1_Rep1").
        # Use cache$design_df (sample → sample_name) to translate before matching.
        bg_remapped <- bg
        design_map  <- cache$design_df
        if (!is.null(design_map) &&
            all(c("sample", "sample_name") %in% colnames(design_map))) {
          col_to_run <- setNames(design_map$sample_name, design_map$sample)
          bg_remapped <- lapply(bg, function(runs) {
            mapped <- col_to_run[runs]
            c(runs[is.na(mapped)], mapped[!is.na(mapped)])
          })
        }

        # Map run IDs → batch labels
        batch_vec <- rep(NA_character_, length(cn_clean))
        for (bl in names(bg_remapped)) {
          batch_vec[cn_clean %in% bg_remapped[[bl]]] <- bl
        }
        # Only correct if all runs have a batch and ≥2 batches present
        if (!any(is.na(batch_vec)) && length(unique(batch_vec)) >= 2) {
          # Save uncorrected PCA BEFORE applying correction
          mat_t_raw   <- t(mat_clean)
          pca_raw     <- prcomp(mat_t_raw, center = TRUE, scale. = FALSE)
          labels_raw  <- run_to_label[colnames(mat_clean)]
          labels_raw[is.na(labels_raw)] <- colnames(mat_clean)[is.na(labels_raw)]
          groups_raw  <- run_to_group[colnames(mat_clean)]
          groups_raw[is.na(groups_raw)] <- "Unknown"
          scores_raw  <- as.data.frame(pca_raw$x)
          scores_raw$Sample <- labels_raw
          scores_raw$Group  <- groups_raw
          pca_rv$pca_result_de <- pca_raw
          pca_rv$pca_scores_de <- scores_raw

          mat_clean <- limma::removeBatchEffect(mat_clean, batch = batch_vec)
          batch_applied <- TRUE
          logger(paste0("PCA [Cache]: removeBatchEffect applied (",
                        length(unique(batch_vec)), " batches)"))
        } else {
          pca_rv$pca_result_de <- NULL
          pca_rv$pca_scores_de <- NULL
          logger(paste0("PCA [Cache]: batch correction skipped (",
                        sum(is.na(batch_vec)), " runs unmatched to any batch)"))
        }
      } else {
        pca_rv$pca_result_de <- NULL
        pca_rv$pca_scores_de <- NULL
      }

      # Match run IDs in column names to labels
      cn         <- colnames(mat_clean)
      labels     <- run_to_label[cn]
      labels[is.na(labels)] <- cn[is.na(labels)]
      groups_vec <- run_to_group[cn]
      groups_vec[is.na(groups_vec)] <- "Unknown"

      incProgress(0.6, detail = "Running PCA...")

      mat_t   <- t(mat_clean)

      pca_out <- prcomp(mat_t, center = TRUE, scale. = FALSE)

      pca_rv$pca_result <- pca_out
      pca_rv$pca_mat    <- mat_t
      pca_rv$pipeline   <- isolate(rv$lf_pipeline)

      pca_sdev <- pca_out$sdev
      pca_var  <- (pca_sdev^2) / sum(pca_sdev^2) * 100

      scores_df <- as.data.frame(pca_out$x)
      scores_df$Sample <- labels
      scores_df$Group  <- groups_vec
      pca_rv$pca_scores <- scores_df
      pca_rv$cv_groups  <- split(seq_len(nrow(scores_df)), groups_vec)

      # Provenance
      pipeline <- isolate(rv$lf_pipeline); ip_sub <- isolate(rv$lf_ip_sub)
      # SAINT cache uses raw log2 intensities — no panel normalisation or LOD
      # imputation was applied.  Flag this so the methods-text branch fires
      # even before a SAINT DE result is available.
      .is_saint_pca <- identical(cache$pipeline_label, "IP (SAINTexpress)")
      # ProDA cache carries normalization_method; MSstats cache does not.
      .is_proda_pca <- !.is_saint_pca && !is.null(cache$normalization_method)
      rv$provenance$pca <- list(
        pipeline         = pipeline, ip_sub = ip_sub,
        skip_impute      = isTRUE(isolate(rv$lf_skip_impute)),
        pipeline_label   = cache$pipeline_label %||% pipeline,
        is_saint_pca     = .is_saint_pca,
        is_proda_pca     = .is_proda_pca,
        ip_complete_case = isTRUE(.is_ip_route_pca),
        n_prots_complete_case = .n_prots_complete_case,
        n_samples        = as.numeric(nrow(mat_t)),
        n_proteins       = as.numeric(ncol(mat_t)),
        n_prots_input        = as.numeric(.n_prots_raw),
        n_prots_after_filter = as.numeric(nrow(mat)),
        prebatch_imputation_method = "column_median",
        prebatch_n_imputed         = .n_prebatch_imputed_pca,
        # Engine before pipeline label: a global (non-IP) run analysed with the
        # ProDA engine is the default non-IP path, and testing pipeline first
        # reported "equalizeMedians" from the label regardless of what the
        # cache said -- including when normalization threw and the cache says
        # *_FAILED_not_applied, which the switch default correctly maps to
        # "none". The SAINT test stays first.
        normalization    = if (.is_saint_pca) "none"
                           else if (.is_proda_pca) switch(cache$normalization_method,
                             global_median              = "equalizeMedians",
                             invariant_median                 = "invariantMedian",
                             invariant_low_confidence         = "invariantMedian",
                             invariant_fallback_global_median = "equalizeMedians",
                             "none"
                           )
                           else if (pipeline == "global") "equalizeMedians"
                           else "FALSE",
        summary_method   = if (.is_proda_pca) NA_character_ else "TMP",
        mb_impute        = (pipeline == "global"),
        lod_imputation   = if (.is_saint_pca || .is_proda_pca) FALSE else (pipeline == "ip"),
        ebn_ref_proteins = cache$ebn_ref_vector %||% character(0),
        igg_conditions   = rv$lf_igg_conditions %||% character(0),
        batch_corrected  = batch_applied,
        centering        = TRUE, scaling = FALSE,
        pc1_var          = round(pca_var[1], 1),
        pc2_var          = round(pca_var[2], 1),
        condition_names  = names(grps),
        software_used    = "prcomp (base R)",
        from_cache       = TRUE
      )

      incProgress(1, detail = "Complete!")
      logger(paste0("PCA [Cache]: Complete. ", nrow(scores_df), " samples, ",
                    ncol(mat_t), " proteins."))
    })
  }

  # =========================================================================
  # PCA Scatter Plot (PC1 vs PC2) - group coloring
  # =========================================================================
  output$pca_scatter <- renderPlotly({
    req(pca_rv$pca_result, pca_rv$pca_scores)

    pca_out   <- pca_rv$pca_result
    scores_df <- pca_rv$pca_scores

    var_explained <- summary(pca_out)$importance[2, ] * 100
    pc1_var <- round(var_explained[1], 1)
    pc2_var <- round(var_explained[2], 1)

    scores_df$hover <- paste0(
      "<b>Sample:</b> ", scores_df$Sample,
      "<br><b>Group:</b> ", scores_df$Group,
      "<br><b>PC1:</b> ", round(scores_df$PC1, 2),
      "<br><b>PC2:</b> ", round(scores_df$PC2, 2)
    )

    # Color by group
    n_groups <- length(unique(scores_df$Group))
    if (n_groups <= 1) {
      # Single group - single color
      p <- plot_ly(scores_df, x = ~PC1, y = ~PC2, type = "scatter",
                   mode = "markers+text",
                   text = ~Sample, textposition = "top center",
                   marker = list(size = 12, opacity = 0.8,
                                 color = "#3498db",
                                 line = list(width = 1, color = "#2c3e50")),
                   hoverinfo = "text", hovertext = ~hover)
    } else {
      # Multiple groups - color by group
      group_colors <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12",
                        "#9b59b6", "#1abc9c", "#e67e22", "#34495e",
                        "#16a085", "#c0392b", "#2980b9", "#27ae60")
      p <- plot_ly(scores_df, x = ~PC1, y = ~PC2, type = "scatter",
                   mode = "markers+text",
                   text = ~Sample, textposition = "top center",
                   color = ~Group,
                   colors = group_colors[seq_len(min(n_groups, length(group_colors)))],
                   marker = list(size = 12, opacity = 0.8,
                                 line = list(width = 1, color = "#2c3e50")),
                   hoverinfo = "text", hovertext = ~hover)
    }

    pipeline <- pca_rv$pipeline
    pipeline_title <- switch(pipeline,
      "global"    = "Global (eqMedians + AFT)",
      "ip"        = if (identical(isolate(rv$lf_ip_sub), "invariant")) "Invariant-protein median (ProDA)" else "IP (SAINTexpress input, no normalization)",
      "PCA"
    )

    p %>% layout(
      title = list(text = paste0("PCA \u2014 ", pipeline_title),
                   font = list(size = 16)),
      xaxis = list(title = paste0("PC1 (", pc1_var, "% variance)"),
                   zeroline = TRUE, zerolinecolor = "#ccc"),
      yaxis = list(title = paste0("PC2 (", pc2_var, "% variance)"),
                   zeroline = TRUE, zerolinecolor = "#ccc"),
      hovermode = "closest"
    )
  })

  # =========================================================================
  # PCA Plot 2: Uncorrected (raw TMP output for batch comparison)
  # When batch correction is active, the main PCA shows protein-level
  # removeBatchEffect. This second tab shows the raw TMP output so users
  # can compare before/after correction.
  # =========================================================================
  output$pca_scatter_de <- renderPlotly({
    req(pca_rv$pca_result_de, pca_rv$pca_scores_de)

    pca_out   <- pca_rv$pca_result_de
    scores_df <- pca_rv$pca_scores_de

    var_explained <- summary(pca_out)$importance[2, ] * 100
    pc1_var <- round(var_explained[1], 1)
    pc2_var <- round(var_explained[2], 1)

    scores_df$hover <- paste0(
      "<b>Sample:</b> ", scores_df$Sample,
      "<br><b>Group:</b> ", scores_df$Group,
      "<br><b>PC1:</b> ", round(scores_df$PC1, 2),
      "<br><b>PC2:</b> ", round(scores_df$PC2, 2)
    )

    n_groups <- length(unique(scores_df$Group))
    group_colors <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12",
                      "#9b59b6", "#1abc9c", "#e67e22", "#34495e",
                      "#16a085", "#c0392b", "#2980b9", "#27ae60")
    if (n_groups <= 1) {
      p <- plot_ly(scores_df, x = ~PC1, y = ~PC2, type = "scatter",
                   mode = "markers+text",
                   text = ~Sample, textposition = "top center",
                   marker = list(size = 12, opacity = 0.8,
                                 color = "#3498db",
                                 line = list(width = 1, color = "#2c3e50")),
                   hoverinfo = "text", hovertext = ~hover)
    } else {
      p <- plot_ly(scores_df, x = ~PC1, y = ~PC2, type = "scatter",
                   mode = "markers+text",
                   text = ~Sample, textposition = "top center",
                   color = ~Group,
                   colors = group_colors[seq_len(min(n_groups, length(group_colors)))],
                   marker = list(size = 12, opacity = 0.8,
                                 line = list(width = 1, color = "#2c3e50")),
                   hoverinfo = "text", hovertext = ~hover)
    }

    p %>% layout(
      title = list(
        text = "PCA \u2014 Uncorrected (Raw TMP)",
        font = list(size = 16)),
      xaxis = list(title = paste0("PC1 (", pc1_var, "% variance)"),
                   zeroline = TRUE, zerolinecolor = "#ccc"),
      yaxis = list(title = paste0("PC2 (", pc2_var, "% variance)"),
                   zeroline = TRUE, zerolinecolor = "#ccc"),
      hovermode = "closest"
    )
  })

  # =========================================================================
  # Scree Plot
  # =========================================================================
  output$pca_scree <- renderPlotly({
    req(pca_rv$pca_result)

    pca_out <- pca_rv$pca_result
    var_explained <- summary(pca_out)$importance[2, ] * 100
    cum_var       <- summary(pca_out)$importance[3, ] * 100
    n_pcs <- length(var_explained)

    scree_df <- data.frame(
      PC         = factor(paste0("PC", seq_len(n_pcs)),
                          levels = paste0("PC", seq_len(n_pcs))),
      Variance   = var_explained,
      Cumulative = cum_var
    )

    plot_ly(scree_df) %>%
      add_bars(x = ~PC, y = ~Variance, name = "Individual",
               marker = list(color = "#3498db", opacity = 0.8)) %>%
      add_trace(x = ~PC, y = ~Cumulative, name = "Cumulative",
                type = "scatter", mode = "lines+markers",
                line = list(color = "#e74c3c", width = 2),
                marker = list(color = "#e74c3c", size = 8)) %>%
      layout(
        title = list(text = "Scree Plot \u2014 Variance Explained",
                     font = list(size = 16)),
        xaxis = list(title = "Principal Component"),
        yaxis = list(title = "% Variance Explained", range = c(0, 105)),
        legend = list(orientation = "h", y = -0.2),
        barmode = "overlay"
      )
  })

  # =========================================================================
  # CV Analysis Plot (density of CV per group)
  # =========================================================================
  output$pca_cv_plot <- renderPlotly({
    req(pca_rv$cv_groups, pca_rv$pca_mat)

    groups <- pca_rv$cv_groups
    mat    <- pca_rv$pca_mat   # rows = samples, cols = proteins

    if (length(groups) == 0) return(plotly_empty())

    # Calculate CV for each protein within each group
    cv_data <- list()
    for (grp_name in names(groups)) {
      grp_samples <- groups[[grp_name]]
      grp_idx     <- which(rownames(mat) %in% grp_samples |
                           pca_rv$pca_scores$Sample %in% grp_samples)
      if (length(grp_idx) < 2) next

      grp_mat <- mat[grp_idx, , drop = FALSE]

      # Convert from log2 to linear for CV calculation
      grp_linear <- 2^grp_mat

      col_means <- colMeans(grp_linear, na.rm = TRUE)
      col_sds   <- apply(grp_linear, 2, sd, na.rm = TRUE)
      cv_vals   <- (col_sds / col_means) * 100

      cv_vals <- cv_vals[is.finite(cv_vals) & !is.na(cv_vals)]

      if (length(cv_vals) > 0) {
        cv_data[[grp_name]] <- data.frame(
          CV    = cv_vals,
          Group = grp_name,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(cv_data) == 0) {
      return(plotly_empty() %>% layout(title = "No valid CV data"))
    }

    all_cv <- do.call(rbind, cv_data)

    group_colors <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12",
                      "#9b59b6", "#1abc9c", "#e67e22", "#34495e")

    p <- plot_ly()
    for (i in seq_along(names(groups))) {
      grp <- names(groups)[i]
      grp_cv <- all_cv[all_cv$Group == grp, ]
      if (nrow(grp_cv) == 0) next

      dens <- density(grp_cv$CV, from = 0,
                      to = max(all_cv$CV, na.rm = TRUE) * 1.1,
                      n = 256, na.rm = TRUE)

      col <- group_colors[((i - 1) %% length(group_colors)) + 1]

      p <- p %>% add_trace(
        x = dens$x, y = dens$y,
        type = "scatter", mode = "lines",
        name = paste0(grp, " (median=", round(median(grp_cv$CV), 1), "%)"),
        fill = "tozeroy",
        fillcolor = paste0(col, "33"),
        line = list(color = col, width = 2)
      )

      med_cv <- median(grp_cv$CV)
      med_y  <- approx(dens$x, dens$y, xout = med_cv)$y
      p <- p %>% add_segments(
        x = med_cv, xend = med_cv, y = 0, yend = med_y,
        line = list(color = col, width = 2, dash = "dash"),
        showlegend = FALSE, hoverinfo = "text",
        text = paste0(grp, " median: ", round(med_cv, 1), "%")
      )
    }

    p %>% layout(
      title = list(text = "Coefficient of Variation \u2014 Density by Group",
                   font = list(size = 16)),
      xaxis = list(title = "CV (%)", rangemode = "tozero"),
      yaxis = list(title = "Density"),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "x unified"
    )
  })

  # =========================================================================
  # Status
  # =========================================================================
  output$pca_status <- renderPrint({
    if (is.null(pca_rv$pca_result)) {
      cat("Upload a CSV, configure the pipeline, and run PCA.\n")
    } else {
      pca_out <- pca_rv$pca_result
      var_ex  <- summary(pca_out)$importance[2, ] * 100
      pipeline <- pca_rv$pipeline
      # Check if DE was run via SAINT — override label if so
      de_prov <- tryCatch(isolate(rv$provenance$de), error = function(e) NULL)
      is_saint_de <- !is.null(de_prov) &&
        identical(de_prov$analysis_mode %||% de_prov$de_engine, "saint")
      pipeline_label <- if (is_saint_de) {
        bait <- de_prov$bait_protein %||% "unknown"
        paste0("IP (SAINTexpress) — Bait: ", bait)
      } else {
        switch(pipeline,
          "global"    = "Global (eqMedians + AFT)",
          "ip"        = if (identical(isolate(rv$lf_ip_sub), "invariant")) "Invariant-protein median (ProDA)" else "IP (SAINTexpress input, no normalization)",
          "PCA"
        )
      }
      cat(paste0("Pipeline: ", pipeline_label, "\n"))
      cat(paste0("Samples: ", nrow(pca_rv$pca_scores), "\n"))
      cat(paste0("Proteins: ", ncol(pca_rv$pca_mat), "\n"))
      cat(paste0("PC1: ", round(var_ex[1], 1), "% | PC2: ",
                 round(var_ex[2], 1), "%\n"))
      cat(paste0("Total PCs: ", length(var_ex), "\n"))

      # Show groups
      scores_df <- pca_rv$pca_scores
      if (!is.null(scores_df$Group)) {
        grp_tbl <- table(scores_df$Group)
        cat("\nSample groups:\n")
        for (g in names(grp_tbl)) {
          cat(paste0("  ", g, ": ", grp_tbl[g], " sample(s)\n"))
        }
      }

      # Show invariant reference proteins
      ebn_refs <- pca_rv$ebn_refs
      if (!is.null(ebn_refs) && length(ebn_refs) > 0) {
        cat(paste0("\nInvariant Reference Proteins (", length(ebn_refs), "):\n"))
        for (ref in ebn_refs) {
          cat(paste0("  ", ref, "\n"))
        }
      }
    }
  })
}
