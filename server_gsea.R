# server_gsea.R ================================================================
# Gene Set Enrichment Analysis (GSEA) module.
# Supports two input modes:
#   (a) Pre-Ranked (.rnk file): signed ranking metric → fgsea via clusterProfiler
#   (b) Expression Matrix (.gct + .cls): abundance matrix → compute ranking → GSEA
#
# Uses msigdbr for gene set databases and clusterProfiler::GSEA() for analysis.
# Called from server.R as: server_gsea(input, output, session, rv, logger)
# ==============================================================================

server_gsea <- function(input, output, session, rv, logger) {

  # ===========================================================================
  # Reactive storage
  # ===========================================================================
  gsea_rv <- reactiveValues(
    mode          = NULL,   # "ranked" | "expression"
    ranked_vec    = NULL,   # named numeric vector (gene → score)
    gct_matrix    = NULL,   # expression matrix (genes × samples)
    gct_genes     = NULL,   # gene name vector
    cls_labels    = NULL,   # phenotype labels (factor)
    gene_sets     = NULL,   # term2gene data.frame for clusterProfiler
    gsea_result   = NULL,   # clusterProfiler GSEA result object
    n_genes_input = 0L,
    status_msg    = ""
  )
  register_experiment_reset(rv, gsea_rv)

  # ===========================================================================
  # MSigDB cache helpers
  # ===========================================================================
  gsea_cache_dir <- file.path("data", "GSEA")
  if (!dir.exists(gsea_cache_dir)) dir.create(gsea_cache_dir, recursive = TRUE)

  gsea_cache_path <- function(species, db_species, category, subcategory) {
    sp_str   <- gsub(" ", "", species)
    sub_str  <- if (is.null(subcategory) || subcategory == "") "ALL" else gsub(":", "-", subcategory)
    fname    <- paste0("msigdb_", sp_str, "_db", db_species, "_", category, "_", sub_str, ".rds")
    file.path(gsea_cache_dir, fname)
  }

  gsea_cache_info_for <- function(species, db_species, category, subcategory) {
    cp <- gsea_cache_path(species, db_species, category, subcategory)
    if (file.exists(cp)) {
      cached <- tryCatch(readRDS(cp), error = function(e) NULL)
      if (!is.null(cached)) {
        list(exists = TRUE, date = cached$download_date, path = cp)
      } else {
        list(exists = FALSE, date = NULL, path = cp)
      }
    } else {
      list(exists = FALSE, date = NULL, path = cp)
    }
  }

  # ===========================================================================
  # Cache status display (reactive to current selections)
  # ===========================================================================
  output$gsea_cache_info <- renderUI({
    species   <- input$gsea_organism
    category  <- input$gsea_category
    subcategory <- input$gsea_subcategory
    req(species, category)

    db_species <- "HS"
    if (species == "Mus musculus") {
      mt <- input$gsea_mouse_geneset_type
      if (!is.null(mt) && mt == "native") db_species <- "MM"
    }
    if (is.null(subcategory) || subcategory == "") subcategory <- NULL

    info <- gsea_cache_info_for(species, db_species, category, subcategory)
    if (info$exists) {
      tagList(
        tags$p(class = "small text-success mt-2",
               icon("database"), " Cached gene sets found",
               tags$br(),
               tags$span(class = "text-muted",
                         paste0("Last downloaded: ", info$date))),
        checkboxInput("gsea_redownload",
                      "Re-download gene sets (may take time)",
                      value = FALSE)
      )
    } else {
      tags$p(class = "small text-muted mt-2",
             icon("cloud-arrow-down"),
             " Gene sets will be downloaded on first run.")
    }
  })

  # ===========================================================================
  # Parse uploaded .rnk file
  # ===========================================================================
  observeEvent(input$gsea_rnk_file, {
    req(input$gsea_rnk_file)
    fpath <- input$gsea_rnk_file$datapath
    fname <- input$gsea_rnk_file$name

    tryCatch({
      rnk <- read.delim(fpath, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("Gene", "Score"))
      rnk <- rnk[!is.na(rnk$Gene) & nchar(trimws(rnk$Gene)) > 0, ]
      rnk <- rnk[is.finite(rnk$Score), ]   # also drops NA

      # Deduplicate: keep highest absolute score
      rnk <- rnk[order(-abs(rnk$Score)), ]
      rnk <- rnk[!duplicated(rnk$Gene), ]

      ranked_vec <- setNames(rnk$Score, rnk$Gene)
      ranked_vec <- sort(ranked_vec, decreasing = TRUE)

      gsea_rv$ranked_vec    <- ranked_vec
      gsea_rv$n_genes_input <- length(ranked_vec)
      gsea_rv$mode          <- "ranked"
      gsea_rv$gsea_result   <- NULL

      logger(paste0("GSEA: Loaded .rnk '", fname, "' with ",
                    length(ranked_vec), " genes"))
      showNotification(
        paste0("Loaded ", length(ranked_vec), " genes from .rnk file"),
        type = "message", duration = 5)
    }, error = function(e) {
      logger(paste0("!!! [GSEA] .rnk parse error: ", e$message))
      showNotification(paste0("Failed to parse .rnk: ", e$message), type = "error")
    })
  })

  # ===========================================================================
  # Parse uploaded .gct file
  # ===========================================================================
  observeEvent(input$gsea_gct_file, {
    req(input$gsea_gct_file)
    fpath <- input$gsea_gct_file$datapath
    fname <- input$gsea_gct_file$name

    tryCatch({
      all_lines <- readLines(fpath, warn = FALSE)

      # GCT v1.2: line 1 = "#1.2", line 2 = "nGenes\tnSamples", line 3+ = data
      if (!grepl("^#1\\.2", all_lines[1])) {
        showNotification("Invalid GCT file: missing #1.2 header.", type = "error")
        return()
      }

      dims <- as.integer(strsplit(all_lines[2], "\t")[[1]])
      data_text <- paste(all_lines[3:length(all_lines)], collapse = "\n")
      gct_df <- read.delim(text = data_text, header = TRUE, stringsAsFactors = FALSE,
                           check.names = FALSE)

      gene_names <- gct_df$NAME
      sample_cols <- setdiff(colnames(gct_df), c("NAME", "Description"))
      expr_mat <- as.matrix(gct_df[, sample_cols, drop = FALSE])
      rownames(expr_mat) <- gene_names
      mode(expr_mat) <- "numeric"

      gsea_rv$gct_matrix    <- expr_mat
      gsea_rv$gct_genes     <- gene_names
      gsea_rv$n_genes_input <- length(gene_names)
      gsea_rv$gsea_result   <- NULL

      logger(paste0("GSEA: Loaded .gct '", fname, "' with ",
                    nrow(expr_mat), " genes, ", ncol(expr_mat), " samples"))
      showNotification(
        paste0("Loaded ", nrow(expr_mat), " genes × ", ncol(expr_mat), " samples"),
        type = "message", duration = 5)
    }, error = function(e) {
      logger(paste0("!!! [GSEA] .gct parse error: ", e$message))
      showNotification(paste0("Failed to parse .gct: ", e$message), type = "error")
    })
  })

  # ===========================================================================
  # Parse uploaded .cls file
  # ===========================================================================
  observeEvent(input$gsea_cls_file, {
    req(input$gsea_cls_file)
    fpath <- input$gsea_cls_file$datapath
    fname <- input$gsea_cls_file$name

    tryCatch({
      cls_lines <- readLines(fpath, warn = FALSE)
      cls_lines <- cls_lines[nchar(trimws(cls_lines)) > 0]

      if (length(cls_lines) < 3) {
        showNotification("Invalid .cls file: expected 3 lines.", type = "error")
        return()
      }

      # Line 2: # ClassName1 ClassName2
      class_names <- strsplit(trimws(cls_lines[2]), "\\s+")[[1]]
      class_names <- class_names[class_names != "#"]

      # Line 3: 0 0 0 1 1 1
      label_indices <- as.integer(strsplit(trimws(cls_lines[3]), "\\s+")[[1]])
      labels <- factor(class_names[label_indices + 1], levels = class_names)

      gsea_rv$cls_labels  <- labels
      gsea_rv$gsea_result <- NULL

      logger(paste0("GSEA: Loaded .cls '", fname, "' — ",
                    paste(table(labels), collapse = "/"), " samples (",
                    paste(class_names, collapse = " vs "), ")"))
      showNotification(
        paste0("Phenotype: ", paste(class_names, collapse = " vs "),
               " (", paste(table(labels), collapse = "/"), " samples)"),
        type = "message", duration = 5)
    }, error = function(e) {
      logger(paste0("!!! [GSEA] .cls parse error: ", e$message))
      showNotification(paste0("Failed to parse .cls: ", e$message), type = "error")
    })
  })

  # ===========================================================================
  # File status indicators
  # ===========================================================================
  output$gsea_rnk_status <- renderUI({
    if (is.null(gsea_rv$ranked_vec))
      tags$p(class = "text-muted small", icon("circle-xmark"), " No .rnk loaded")
    else
      tags$p(class = "small text-success",
             icon("circle-check"), " ",
             gsea_rv$n_genes_input, " genes loaded")
  })

  output$gsea_gct_status <- renderUI({
    if (is.null(gsea_rv$gct_matrix))
      tags$p(class = "text-muted small", icon("circle-xmark"), " No .gct loaded")
    else
      tags$p(class = "small text-success",
             icon("circle-check"), " ",
             nrow(gsea_rv$gct_matrix), " genes × ",
             ncol(gsea_rv$gct_matrix), " samples")
  })

  output$gsea_cls_status <- renderUI({
    if (is.null(gsea_rv$cls_labels))
      tags$p(class = "text-muted small", icon("circle-xmark"), " No .cls loaded")
    else {
      tbl <- table(gsea_rv$cls_labels)
      tags$p(class = "small text-success",
             icon("circle-check"), " ",
             paste(names(tbl), "=", tbl, collapse = ", "))
    }
  })

  # ===========================================================================
  # MSigDB category choices (reactive on subcategory)
  # ===========================================================================
  msigdb_categories <- list(
    "H \u2014 Hallmark"                    = "H",
    "C1 \u2014 Positional"                 = "C1",
    "C2 \u2014 Curated Pathways"           = "C2",
    "C3 \u2014 Regulatory Target"          = "C3",
    "C4 \u2014 Computational"              = "C4",
    "C5 \u2014 Gene Ontology"              = "C5",
    "C6 \u2014 Oncogenic Signatures"       = "C6",
    "C7 \u2014 Immunologic Signatures"     = "C7",
    "C8 \u2014 Cell Type Signatures"       = "C8"
  )

  # Mapping: human collection IDs → mouse-native collection IDs
  # Mouse MSigDB uses M-prefixed names. C4/C6 have no mouse equivalent.
  human_to_mouse_collection <- c(
    "H"  = "MH",
    "C1" = "M1",
    "C2" = "M2",
    "C3" = "M3",
    "C4" = NA_character_,   # not available in mouse-native
    "C5" = "M5",
    "C6" = NA_character_,   # not available in mouse-native
    "C7" = "M7",
    "C8" = "M8"
  )

  # Mouse-native subcategories (different from human)
  mouse_native_subcats <- list(
    "H"  = c("(none)" = ""),
    "C1" = c("(none)" = ""),
    "C2" = c("All" = "", "CP:BIOCARTA" = "CP:BIOCARTA",
             "CP:REACTOME" = "CP:REACTOME",
             "CP:WIKIPATHWAYS" = "CP:WIKIPATHWAYS"),
    "C3" = c("All" = "", "GTRD" = "GTRD", "MIRDB" = "MIRDB"),
    "C5" = c("All" = "", "GO:BP" = "GO:BP", "GO:CC" = "GO:CC",
             "GO:MF" = "GO:MF", "MPT" = "MPT"),
    "C7" = c("(none)" = ""),
    "C8" = c("(none)" = "")
  )

  # Human subcategories
  human_subcats <- list(
    "H"  = c("(none)" = ""),
    "C1" = c("(none)" = ""),
    "C2" = c("All" = "", "CP:BIOCARTA" = "CP:BIOCARTA",
             "CP:KEGG_MEDICUS" = "CP:KEGG_MEDICUS",
             "CP:PID" = "CP:PID", "CP:REACTOME" = "CP:REACTOME",
             "CP:WIKIPATHWAYS" = "CP:WIKIPATHWAYS"),
    "C3" = c("All" = "", "TFT:GTRD" = "TFT:GTRD", "TFT:TFT_Legacy" = "TFT:TFT_Legacy",
             "MIR:MIRDB" = "MIR:MIRDB", "MIR:MIR_Legacy" = "MIR:MIR_Legacy"),
    "C4" = c("(none)" = ""),
    "C5" = c("All" = "", "GO:BP" = "GO:BP", "GO:CC" = "GO:CC", "GO:MF" = "GO:MF",
             "HPO" = "HPO"),
    "C6" = c("(none)" = ""),
    "C7" = c("All" = "", "IMMUNESIGDB" = "IMMUNESIGDB", "VAX" = "VAX"),
    "C8" = c("(none)" = "")
  )

  # Helper: is current mode mouse-native?
  is_mouse_native <- reactive({
    input$gsea_organism == "Mus musculus" &&
      !is.null(input$gsea_mouse_geneset_type) &&
      input$gsea_mouse_geneset_type == "native"
  })

  # Update subcategory choices when category OR mouse-type changes
  observe({
    cat     <- input$gsea_category
    req(cat)
    native  <- is_mouse_native()

    if (native) {
      # C4 and C6 not available in mouse-native
      if (cat %in% c("C4", "C6")) {
        choices <- c("(not available for mouse-native)" = "")
      } else {
        choices <- if (cat %in% names(mouse_native_subcats))
          mouse_native_subcats[[cat]] else c("(none)" = "")
      }
    } else {
      choices <- if (cat %in% names(human_subcats))
        human_subcats[[cat]] else c("(none)" = "")
    }
    updateSelectInput(session, "gsea_subcategory", choices = choices)
  })

  # ===========================================================================
  # RUN GSEA
  # ===========================================================================
  observeEvent(input$gsea_run_btn, {

    # ── Immediately show a persistent notification ──
    showNotification(
      "Preparing GSEA analysis...(v2.1 Downsampling active)",
      id = "gsea_progress", type = "message",
      duration = NULL, closeButton = FALSE
    )

    mode <- input$gsea_input_mode

    # --- Validate inputs ---
    if (mode == "ranked") {
      if (is.null(gsea_rv$ranked_vec)) {
        removeNotification("gsea_progress")
        showNotification("Please upload a .rnk file first.", type = "warning")
        return()
      }
    } else {
      if (is.null(gsea_rv$gct_matrix)) {
        removeNotification("gsea_progress")
        showNotification("Please upload a .gct file first.", type = "warning")
        return()
      }
      if (is.null(gsea_rv$cls_labels)) {
        removeNotification("gsea_progress")
        showNotification("Please upload a .cls file first.", type = "warning")
        return()
      }
      if (ncol(gsea_rv$gct_matrix) != length(gsea_rv$cls_labels)) {
        removeNotification("gsea_progress")
        showNotification(
          paste0("Sample count mismatch: .gct has ", ncol(gsea_rv$gct_matrix),
                 " samples but .cls has ", length(gsea_rv$cls_labels), " labels."),
          type = "error")
        return()
      }
    }

    # --- Load dependencies ---
    require_pkg("msigdbr")
    require_pkg("clusterProfiler")

    shinyjs::disable("gsea_run_btn")
    on.exit({
      shinyjs::enable("gsea_run_btn")
      removeNotification("gsea_progress")
    }, add = TRUE)

    gsea_rv$gsea_result <- NULL

    # ── Determine species & db_species ──
    species    <- input$gsea_organism
    db_species <- "HS"              # default: human MSigDB
    if (species == "Mus musculus") {
      mt <- input$gsea_mouse_geneset_type
      if (!is.null(mt) && mt == "native") db_species <- "MM"
    }

    category    <- input$gsea_category
    subcategory <- input$gsea_subcategory
    if (is.null(subcategory) || subcategory == "") subcategory <- NULL

    # ── Map collection names for mouse-native ──
    api_collection    <- category
    api_subcollection <- subcategory
    if (db_species == "MM") {
      mapped <- human_to_mouse_collection[[category]]
      if (is.na(mapped)) {
        showNotification(
          paste0("Collection '", category,
                 "' is not available in mouse-native MSigDB. ",
                 "Try using human orthologs instead."),
          type = "error")
        return()
      }
      api_collection <- mapped
      # Note: mouse subcollection names may differ
    }

    logger(paste0("GSEA: Species=", species,
                  " db_species=", db_species,
                  " Collection=", api_collection,
                  " Subcollection=", ifelse(is.null(api_subcollection), "none", api_subcollection)))

    # ═══════════════════════════════════════════════════════════════════════
    # Step 1 — Load / cache gene sets
    # ═══════════════════════════════════════════════════════════════════════
    cache_path <- gsea_cache_path(species, db_species, api_collection, api_subcollection)
    use_cache  <- file.exists(cache_path) &&
                  (is.null(input$gsea_redownload) || !input$gsea_redownload)

    if (use_cache) {
      showNotification(
        "Loading cached gene sets...",
        id = "gsea_progress", type = "message",
        duration = NULL, closeButton = FALSE
      )
      cached <- tryCatch(readRDS(cache_path), error = function(e) {
        logger(paste0("GSEA: Cache read failed (", e$message, ") — will re-download"))
        NULL
      })
      use_cache <- !is.null(cached)
    }

    if (use_cache) {
      gene_sets <- cached$data
      logger(paste0("GSEA: Using cached gene sets from ", cached$download_date))
    } else {
      showNotification(
        "Downloading gene sets from MSigDB (this may take a moment)...",
        id = "gsea_progress", type = "message",
        duration = NULL, closeButton = FALSE
      )

      gene_sets <- tryCatch({
        args <- list(species = species)
        if (db_species != "HS") args$db_species <- db_species
        args$collection <- api_collection
        if (!is.null(api_subcollection)) args$subcollection <- api_subcollection
        suppressWarnings(do.call(msigdbr::msigdbr, args))
      }, error = function(e) {
        logger(paste0("!!! [GSEA] msigdbr error: ", e$message))
        NULL
      })

      if (is.null(gene_sets) || nrow(gene_sets) == 0) {
        showNotification("No gene sets found for this species/category.",
                         type = "error")
        return()
      }

      # ── Save to local cache ──
      tryCatch({
        if (!dir.exists(gsea_cache_dir)) dir.create(gsea_cache_dir, recursive = TRUE)
        saveRDS(list(data = gene_sets, download_date = Sys.Date()), cache_path)
        logger(paste0("GSEA: Cached gene sets to ", cache_path))
      }, error = function(e) {
        logger(paste0("GSEA: Could not cache gene sets: ", e$message))
      })
    }

    # Reset the re-download checkbox after use
    if (!is.null(input$gsea_redownload) && input$gsea_redownload) {
      updateCheckboxInput(session, "gsea_redownload", value = FALSE)
    }

    term2gene <- gene_sets[, c("gs_name", "gene_symbol")]
    n_sets    <- length(unique(term2gene$gs_name))
    logger(paste0("GSEA: Loaded ", n_sets, " gene sets (",
                  nrow(term2gene), " gene-set pairs)"))

    showNotification(
      paste0("Loaded ", n_sets, " gene sets. Building ranked list..."),
      id = "gsea_progress", type = "message",
      duration = NULL, closeButton = FALSE
    )

    # ═══════════════════════════════════════════════════════════════════════
    # Step 2 — Build ranked list
    # ═══════════════════════════════════════════════════════════════════════
    if (mode == "ranked") {
      gene_list <- gsea_rv$ranked_vec
    } else {
      showNotification(
        "Computing ranking from expression data (limma)...",
        id = "gsea_progress", type = "message",
        duration = NULL, closeButton = FALSE
      )
      require_pkg("limma")

      expr_mat <- gsea_rv$gct_matrix
      labels   <- gsea_rv$cls_labels

      # Sanitize factor levels to syntactically valid R names
      # (e.g. "48h" → "X48h") — required by limma::makeContrasts
      safe_levels <- make.names(levels(factor(labels)))
      group  <- factor(make.names(as.character(labels)), levels = safe_levels)
      design <- model.matrix(~ 0 + group)
      colnames(design) <- levels(group)

      contrast_str <- paste0(levels(group)[2], " - ", levels(group)[1])
      contrast_mat <- limma::makeContrasts(contrasts = contrast_str,
                                           levels = design)

      fit  <- limma::lmFit(expr_mat, design)
      fit2 <- limma::contrasts.fit(fit, contrast_mat)
      fit2 <- limma::eBayes(fit2)
      tt   <- limma::topTable(fit2, number = Inf, sort.by = "none")

      ranking   <- sign(tt$logFC) * (-log10(pmax(tt$P.Value, 1e-300)))
      gene_list <- setNames(ranking, rownames(expr_mat))

      # Drop non-finite values (zero-variance rows produce NaN in limma)
      n_nonfinite <- sum(!is.finite(gene_list))
      if (n_nonfinite > 0) {
        gene_list <- gene_list[is.finite(gene_list)]
        logger(paste0("GSEA: Dropped ", n_nonfinite,
                      " genes with non-finite ranking (likely zero-variance rows)"))
      }

      gene_list <- sort(gene_list, decreasing = TRUE)
      gene_list <- gene_list[!duplicated(names(gene_list))]

      logger(paste0("GSEA: Computed ranking for ", length(gene_list), " genes"))
    }

    # ═══════════════════════════════════════════════════════════════════════
    # Step 3 — Run GSEA
    # ═══════════════════════════════════════════════════════════════════════
    showNotification(
      paste0("Running GSEA on ", length(gene_list), " genes × ",
             n_sets, " gene sets (10 000 permutations)..."),
      id = "gsea_progress", type = "message",
      duration = NULL, closeButton = FALSE
    )

    set.seed(42)
    gsea_res <- tryCatch({
      suppressWarnings(clusterProfiler::GSEA(
        geneList      = gene_list,
        TERM2GENE     = term2gene,
        minGSSize     = 15,
        maxGSSize     = 500,
        pvalueCutoff  = 1,
        pAdjustMethod = "BH",
        seed          = TRUE,
        nPermSimple   = 10000,
        verbose       = FALSE
      ))
    }, error = function(e) {
      logger(paste0("!!! [GSEA] clusterProfiler error: ", e$message))
      NULL
    })

    if (is.null(gsea_res)) {
      showNotification("GSEA analysis failed. Check log for details.",
                       type = "error")
      return()
    }

    n_sig   <- sum(gsea_res@result$p.adjust < 0.05)
    n_total <- nrow(gsea_res@result)

    gsea_rv$gsea_result <- gsea_res
    gsea_rv$mode        <- mode
    gsea_rv$status_msg  <- paste0(
      "GSEA complete: ", n_total, " gene sets tested, ",
      n_sig, " significant (FDR < 0.05)")

    logger(paste0("GSEA: Complete — ", n_total, " sets tested, ",
                  n_sig, " significant (FDR < 0.05)"))

    # Remove the persistent progress notification (on.exit also does this)
    removeNotification("gsea_progress")

    # Final success notification (separate ID so on.exit doesn't kill it)
    showNotification(
      paste0(n_sig, " / ", n_total, " gene sets significant (FDR < 0.05)"),
      id = "gsea_done", type = "message", duration = 8)
  })

  # ===========================================================================
  # GSEA Results — Summary Table
  # ===========================================================================
  output$gsea_result_table <- DT::renderDT({
    req(gsea_rv$gsea_result)
    res <- gsea_rv$gsea_result@result

    show_df <- data.frame(
      Pathway       = res$Description,
      SetSize       = res$setSize,
      ES            = round(res$enrichmentScore, 4),
      NES           = round(res$NES, 4),
      PValue        = signif(res$pvalue, 4),
      FDR           = signif(res$p.adjust, 4),
      stringsAsFactors = FALSE
    )

    show_df <- show_df[order(show_df$FDR, -abs(show_df$NES)), ]

    DT::datatable(show_df, rownames = FALSE, filter = "top",
      selection = "single",
      options = list(
        pageLength = 20, scrollX = TRUE,
        dom = "Bfrtip", buttons = c("csv", "excel"),
        order = list(list(5, "asc"))
      ),
      extensions = "Buttons", class = "compact stripe hover"
    ) |>
      DT::formatStyle("NES",
        color = DT::styleInterval(0, c("#3498db", "#e74c3c")),
        fontWeight = "bold"
      ) |>
      DT::formatStyle("FDR",
        backgroundColor = DT::styleInterval(
          c(0.01, 0.05, 0.25),
          c("#d4edda", "#d1ecf1", "#fff3cd", "#f8f9fa")
        )
      )
  })

  # ===========================================================================
  # GSEA Enrichment Plot — interactive (plotly)
  # ===========================================================================
  # ===========================================================================
  # Helper: resolve selected pathway ID from table click
  # ===========================================================================
  .gsea_selected_pathway <- reactive({
    req(gsea_rv$gsea_result)
    res <- gsea_rv$gsea_result@result
    sel_row <- input$gsea_result_table_rows_selected
    if (!is.null(sel_row) && length(sel_row) == 1) {
      sorted_res <- res[order(res$p.adjust, -abs(res$NES)), ]
      sorted_res$ID[sel_row]
    } else {
      res$ID[which.max(abs(res$NES))]
    }
  })

  # ===========================================================================
  # Helper: compute enrichment data for a pathway (shared by plot + table)
  # ===========================================================================
  .gsea_pathway_data <- reactive({
    pathway_id <- .gsea_selected_pathway()
    gsea_obj <- gsea_rv$gsea_result
    res <- gsea_obj@result

    gene_list <- gsea_obj@geneList
    t2g <- gsea_obj@geneSets[[pathway_id]]
    ranked_names <- names(gene_list)
    n <- length(ranked_names)
    hits <- ranked_names %in% t2g

    n_hit <- sum(hits)
    n_miss <- n - n_hit
    if (n_hit == 0) return(NULL)

    hit_scores <- abs(gene_list[hits])
    nr <- sum(hit_scores)

    running_es <- numeric(n)
    es_running <- 0
    for (i in seq_len(n)) {
      if (hits[i]) {
        es_running <- es_running + abs(gene_list[i]) / nr
      } else {
        es_running <- es_running - 1 / n_miss
      }
      running_es[i] <- es_running
    }

    core_genes <- strsplit(res$core_enrichment[res$ID == pathway_id], "/")[[1]]

    list(
      pathway_id   = pathway_id,
      pathway_name = res$Description[res$ID == pathway_id],
      nes_val      = round(res$NES[res$ID == pathway_id], 3),
      fdr_val      = signif(res$p.adjust[res$ID == pathway_id], 4),
      es_val       = round(res$enrichmentScore[res$ID == pathway_id], 4),
      gene_list    = gene_list,
      ranked_names = ranked_names,
      hits         = hits,
      running_es   = running_es,
      core_genes   = core_genes,
      n            = n
    )
  })

  # ===========================================================================
  # Enrichment Plot — classic 3-panel (Running ES, barcode, ranked metric)
  # ===========================================================================
  output$gsea_enrichment_plot <- renderPlotly({
    pd <- .gsea_pathway_data()
    if (is.null(pd)) {
      return(plotly_empty() %>% layout(title = "No hits for this pathway"))
    }

    n           <- pd$n
    rank_pos    <- seq_len(n)
    running_es  <- pd$running_es
    hits        <- pd$hits
    gene_list   <- pd$gene_list
    ranked_names <- pd$ranked_names
    hit_positions <- which(hits)
    peak_idx    <- which.max(abs(running_es))

    # --- Panel 1: Running ES curve ---
    p1 <- plot_ly() %>%
      add_trace(x = rank_pos, y = running_es,
                type = "scatter", mode = "lines",
                name = "Running ES",
                line = list(color = "#18bc9c", width = 2),
                hoverinfo = "text",
                text = paste0("Rank: ", rank_pos,
                              "<br>Running ES: ", round(running_es, 4),
                              "<br>Gene: ", ranked_names)) %>%
      add_trace(x = c(1, n), y = c(0, 0),
                type = "scatter", mode = "lines",
                showlegend = FALSE,
                line = list(color = "#999", width = 1, dash = "dot"),
                hoverinfo = "none") %>%
      layout(
        yaxis = list(title = "Enrichment Score (ES)"),
        xaxis = list(title = "", showticklabels = FALSE),
        shapes = list(
          list(type = "line",
               x0 = peak_idx, x1 = peak_idx,
               y0 = 0, y1 = running_es[peak_idx],
               line = list(color = "#e74c3c", width = 1, dash = "dash"))
        ),
        annotations = list(
          list(x = peak_idx, y = running_es[peak_idx],
               text = paste0("ES = ", round(running_es[peak_idx], 3)),
               showarrow = TRUE, arrowhead = 2, arrowsize = 0.8,
               font = list(size = 10, color = "#e74c3c"))
        )
      )

    # --- Panel 2: Barcode hits ---
    p2 <- plot_ly() %>%
      add_trace(x = hit_positions,
                y = rep(0, length(hit_positions)),
                type = "scatter", mode = "markers",
                name = "Gene Set Hits",
                marker = list(symbol = "line-ns", size = 12,
                              color = "#000", opacity = 0.7,
                              line = list(width = 1, color = "#000")),
                hoverinfo = "text",
                text = paste0("Gene: ", ranked_names[hit_positions],
                              "<br>Rank: ", hit_positions)) %>%
      layout(
        xaxis = list(title = "", showticklabels = FALSE, range = c(1, n)),
        yaxis = list(title = "", showticklabels = FALSE,
                     fixedrange = TRUE, range = c(-0.5, 0.5)),
        margin = list(t = 0, b = 0)
      )

    # --- Panel 3: Ranked metric (gene scores) ---
# --- Downsample the ranked metric for the background trace ---
    # We keep all 'hits' to ensure accuracy, but aggressively sample the empty space.
    n_total <- length(gene_list)

    if (n_total > 2000) {
      # Always keep the first, last, peak, and all hit positions
      keep_idx <- unique(c(1, n, peak_idx, hit_positions))
      
      # Sample the remaining background to max out around 1500 points
      bg_idx <- setdiff(seq_len(n_total), keep_idx)
      sample_size <- min(1500, length(bg_idx))
      set.seed(42)
      bg_sampled <- sample(bg_idx, sample_size)
      
      # Recombine and sort
      final_idx <- sort(unique(c(keep_idx, bg_sampled)))
      
      plot_x <- rank_pos[final_idx]
      plot_y <- as.numeric(gene_list)[final_idx]
      plot_names <- ranked_names[final_idx]
    } else {
      plot_x <- rank_pos
      plot_y <- as.numeric(gene_list)
      plot_names <- ranked_names
    }

    # --- Panel 3: Ranked metric (gene scores) ---
    p3 <- plot_ly() %>%
      add_trace(x = plot_x, y = plot_y,
                type = "scatter", mode = "lines",
                name = "Ranking Metric",
                line = list(color = "#aaa", width = 1),
                fill = "tozeroy",
                fillcolor = "rgba(180,180,180,0.3)",
                hoverinfo = "text",
                text = paste0("Rank: ", plot_x,
                              "<br>Score: ", round(plot_y, 4),
                              "<br>Gene: ", plot_names)) %>%
      add_trace(x = c(1, n), y = c(0, 0),
                type = "scatter", mode = "lines",
                showlegend = FALSE,
                line = list(color = "#999", width = 1, dash = "dot"),
                hoverinfo = "none") %>%
      layout(
        xaxis = list(title = "Rank in Ordered Gene List"),
        yaxis = list(title = "Ranked Metric")
      )

    # --- Combine as 3-panel subplot ---
    p <- subplot(p1, p2, p3, nrows = 3, shareX = TRUE,
                 heights = c(0.55, 0.08, 0.37)) %>%
      layout(
        title = list(
          text = paste0(pd$pathway_name,
                        "<br><sup>NES = ", pd$nes_val,
                        " | FDR = ", pd$fdr_val,
                        " | ES = ", pd$es_val, "</sup>"),
          font = list(size = 14)
        ),
        showlegend = FALSE,
        hovermode = "x unified",
        margin = list(t = 70)
      )
    p
  })

  # ===========================================================================
  # Leading Edge Detail Table — per-gene data for selected pathway
  # ===========================================================================
  output$gsea_leading_edge_table <- DT::renderDT({
    pd <- .gsea_pathway_data()
    req(pd)

    gene_list    <- pd$gene_list
    ranked_names <- pd$ranked_names
    hits         <- pd$hits
    running_es   <- pd$running_es
    core_genes   <- pd$core_genes

    hit_positions <- which(hits)

    detail_df <- data.frame(
      Gene              = ranked_names[hit_positions],
      Rank              = hit_positions,
      `Rank Metric`     = round(gene_list[hit_positions], 4),
      `Running ES`      = round(running_es[hit_positions], 4),
      `Core Enrichment` = ifelse(ranked_names[hit_positions] %in% core_genes,
                                 "Yes", "No"),
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    )

    DT::datatable(detail_df, rownames = FALSE,
      selection = "none",
      options = list(
        pageLength = 25, scrollX = TRUE,
        dom = "Bfrtip", buttons = c("csv", "excel"),
        order = list(list(1, "asc"))
      ),
      extensions = "Buttons",
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("Core Enrichment",
        backgroundColor = DT::styleEqual(
          c("Yes", "No"), c("#d4edda", "#f8f9fa")
        ),
        fontWeight = DT::styleEqual(c("Yes", "No"), c("bold", "normal"))
      ) |>
      DT::formatStyle("Rank Metric",
        color = DT::styleInterval(0, c("#3498db", "#e74c3c"))
      )
  })

  output$gsea_leading_edge_header <- renderUI({
    pd <- .gsea_pathway_data()
    req(pd)
    n_core <- sum(pd$ranked_names[which(pd$hits)] %in% pd$core_genes)
    n_hits <- sum(pd$hits)
    tags$div(
      tags$h6(class = "mb-1",
        icon("list"), " ", pd$pathway_name
      ),
      tags$p(class = "small text-muted mb-0",
        paste0(n_hits, " genes in set, ", n_core,
               " in leading edge (core enrichment). ",
               "NES = ", pd$nes_val, ", FDR = ", pd$fdr_val))
    )
  })

  # ===========================================================================
  # Dot plot — top pathways overview
  # ===========================================================================
  output$gsea_dotplot <- renderPlotly({
    req(gsea_rv$gsea_result)
    res <- gsea_rv$gsea_result@result

    # Top 20 by FDR, split by direction
    res_up   <- res[res$NES > 0, ]
    res_down <- res[res$NES < 0, ]
    res_up   <- head(res_up[order(res_up$p.adjust), ], 10)
    res_down <- head(res_down[order(res_down$p.adjust), ], 10)
    top <- rbind(res_up, res_down)

    if (nrow(top) == 0) {
      return(plotly_empty() %>% layout(title = "No results to display"))
    }

    top$Direction <- ifelse(top$NES > 0, "Upregulated", "Downregulated")
    top$NegLog10FDR <- -log10(pmax(top$p.adjust, 1e-300))

    # Truncate long names
    top$ShortName <- ifelse(nchar(top$Description) > 50,
                            paste0(substr(top$Description, 1, 47), "..."),
                            top$Description)

    # Order by NES
    top <- top[order(top$NES), ]
    top$ShortName <- factor(top$ShortName, levels = top$ShortName)

    p <- plot_ly(data = top,
      y = ~ShortName, x = ~NES,
      type = "scatter", mode = "markers",
      color = ~Direction,
      colors = c("Upregulated" = "#e74c3c", "Downregulated" = "#3498db"),
      size = ~setSize, sizes = c(8, 25),
      marker = list(opacity = 0.8,
                    line = list(width = 0.5, color = "#333")),
      hoverinfo = "text",
      text = ~paste0(
        "<b>", Description, "</b>",
        "<br>NES: ", round(NES, 3),
        "<br>FDR: ", signif(p.adjust, 4),
        "<br>Set Size: ", setSize,
        "<br>ES: ", round(enrichmentScore, 4)
      )
    )

    p <- layout(p,
      title = list(text = "Top Enriched Gene Sets", font = list(size = 16)),
      xaxis = list(title = "Normalized Enrichment Score (NES)",
                   zeroline = TRUE, zerolinecolor = "#999"),
      yaxis = list(title = "", tickfont = list(size = 10)),
      margin = list(l = 250),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "closest"
    )
    p
  })

  # ===========================================================================
  # NES bar chart
  # ===========================================================================
  output$gsea_nes_barplot <- renderPlotly({
    req(gsea_rv$gsea_result)
    res <- gsea_rv$gsea_result@result

    sig <- res[res$p.adjust < 0.25, ]
    if (nrow(sig) == 0) sig <- head(res[order(res$p.adjust), ], 20)
    sig <- head(sig[order(-abs(sig$NES)), ], 30)
    sig <- sig[order(sig$NES), ]

    sig$ShortName <- ifelse(nchar(sig$Description) > 45,
                            paste0(substr(sig$Description, 1, 42), "..."),
                            sig$Description)
    sig$ShortName <- factor(sig$ShortName, levels = sig$ShortName)
    sig$Color <- ifelse(sig$NES > 0, "#e74c3c", "#3498db")

    p <- plot_ly(data = sig,
      y = ~ShortName, x = ~NES,
      type = "bar", orientation = "h",
      marker = list(color = sig$Color,
                    line = list(width = 0.5, color = "#333")),
      hoverinfo = "text",
      text = ~paste0(
        "<b>", Description, "</b>",
        "<br>NES: ", round(NES, 3),
        "<br>FDR: ", signif(p.adjust, 4)
      )
    )

    p <- layout(p,
      title = list(text = "NES Bar Chart (FDR < 0.25)", font = list(size = 16)),
      xaxis = list(title = "Normalized Enrichment Score (NES)"),
      yaxis = list(title = "", tickfont = list(size = 10)),
      margin = list(l = 280),
      showlegend = FALSE
    )
    p
  })

  # ===========================================================================
  # Status
  # ===========================================================================
  output$gsea_status <- renderPrint({
    if (!is.null(gsea_rv$gsea_result)) {
      res <- gsea_rv$gsea_result@result
      n_up   <- sum(res$NES > 0 & res$p.adjust < 0.05)
      n_down <- sum(res$NES < 0 & res$p.adjust < 0.05)
      cat(gsea_rv$status_msg, "\n")
      cat(paste0("  Upregulated:   ", n_up, " (FDR < 0.05)\n"))
      cat(paste0("  Downregulated: ", n_down, " (FDR < 0.05)\n"))
    } else {
      cat("Upload data and run GSEA to see results.\n")
    }
  })

  # ===========================================================================
  # Download results CSV
  # ===========================================================================
  output$gsea_download_csv <- downloadHandler(
    filename = function() {
      paste0("GSEA_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(gsea_rv$gsea_result)
      write.csv(gsea_rv$gsea_result@result, file, row.names = FALSE)
      logger("GSEA: Results CSV downloaded")
    }
  )
}
