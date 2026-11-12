# server.R (Development_v2) - Thin dispatcher
# Defines shared reactive values, tab sync, and delegates to module servers.
# ==============================================================================

## Ensure module server functions are available (needed when launched via runApp)

if (!exists("server_residue_prediction", mode = "function")) {
    source("server_residue_prediction.R", local = .GlobalEnv)
}
if (!exists("server_local_ms", mode = "function")) {
    source("server_local_ms.R", local = .GlobalEnv)
}
if (!exists("server_trend_clusters", mode = "function")) {
    source("server_trend_clusters.R", local = .GlobalEnv)
}
if (!exists("server_de", mode = "function")) {
  source("server_de.R", local = .GlobalEnv)
}
if (!exists("server_pca", mode = "function")) {
  source("server_pca.R", local = .GlobalEnv)
}
if (!exists("server_heatmap", mode = "function")) {
  source("server_heatmap.R", local = .GlobalEnv)
}
if (!exists("server_protein_coverage", mode = "function")) {
  source("server_protein_coverage.R", local = .GlobalEnv)
}
if (!exists("server_fragpipe_db", mode = "function")) {
  source("server_fragpipe_db.R", local = .GlobalEnv)
}
if (!exists("server_rank_fc", mode = "function")) {
  source("server_rank_fc.R", local = .GlobalEnv)
}
if (!exists("server_trends", mode = "function")) {
  source("server_trends.R", local = .GlobalEnv)
}
if (!exists("server_tmt_upload", mode = "function")) {
  source("server_tmt_upload.R", local = .GlobalEnv)
}
if (!exists("server_gsea", mode = "function")) {
  source("server_gsea.R", local = .GlobalEnv)
}
if (!exists("server_tmt_ptm_upload", mode = "function")) {
  source("server_tmt_ptm_upload.R", local = .GlobalEnv)
}
if (!exists("server_ptm_de", mode = "function")) {
  source("server_ptm_de.R", local = .GlobalEnv)
}
if (!exists("server_ptm_trends", mode = "function")) {
  source("server_ptm_trends.R", local = .GlobalEnv)
}
if (!exists("server_ptm_single", mode = "function")) {
  source("server_ptm_single.R", local = .GlobalEnv)
}
if (!exists("server_dia_upload", mode = "function")) {
  source("server_dia_upload.R", local = .GlobalEnv)
}
if (!exists("server_dia_ptm_upload", mode = "function")) {
  source("server_dia_ptm_upload.R", local = .GlobalEnv)
}
if (!exists("server_prerun", mode = "function")) {
  source("server_prerun.R", local = .GlobalEnv)
}
if (!exists("server_ptm_pca", mode = "function")) {
  source("server_ptm_pca.R", local = .GlobalEnv)
}
if (!exists("server_ks_scoring", mode = "function")) {
  source("server_ks_scoring.R", local = .GlobalEnv)
}
if (!exists("server_cross_exp", mode = "function")) {
  source("server_cross_exp.R", local = .GlobalEnv)
}
if (!exists("server_ptm_cross_exp", mode = "function")) {
  source("server_ptm_cross_exp.R", local = .GlobalEnv)
}

server <- function(input, output, session) {
  message("[DBG] *** server() called ***")
  # Short session identifier for log tracing (first 8 chars of token)
  .sid <- substr(session$token, 1, 8)

  osmss_file_log("Session started", level = "INFO", module = "Server",
                 session_id = .sid)
  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] DEBUG: Server Session Started (OSMSS) [", .sid, "]"))

  # ===========================================================================
  # Session inactivity timeout (3 hours)
  # Client-side JS tracks idle time; warns at 2.5 hr, disconnects at 3 hr.
  # Any click, keypress, or mousemove resets the timer.
  # ===========================================================================
  insertUI(
    selector = "body",
    where    = "beforeEnd",
    immediate = TRUE,
    ui = tags$script(HTML("
      (function() {
        var WARN_MS  = 150 * 60 * 1000;  // 2.5 hours (150 minutes)
        var KILL_MS  = 180 * 60 * 1000;  // 3 hours   (180 minutes)
        var warnTimer, killTimer, warned = false;

        function resetTimers() {
          clearTimeout(warnTimer);
          clearTimeout(killTimer);
          if (warned) {
            // Dismiss warning modal if user returned
            Shiny.setInputValue('session_timeout_reset', Math.random());
            warned = false;
          }
          warnTimer = setTimeout(function() {
            warned = true;
            Shiny.setInputValue('session_timeout_warn', Math.random());
          }, WARN_MS);
          killTimer = setTimeout(function() {
            Shiny.setInputValue('session_timeout_kill', Math.random());
          }, KILL_MS);
        }

        ['click', 'keypress', 'mousemove', 'scroll', 'touchstart'].forEach(function(evt) {
          document.addEventListener(evt, resetTimers, {passive: true});
        });
        resetTimers();
      })();
    "))
  )

  # Show warning modal at 2.5 hours of inactivity
  observeEvent(input$session_timeout_warn, {
    showModal(modalDialog(
      title = tags$span(icon("clock"), " Session Timeout Warning"),
      tags$div(class = "alert alert-warning mb-0",
        tags$p(icon("triangle-exclamation"),
               tags$strong(" Your session has been idle for 2.5 hours.")),
        tags$p("It will disconnect automatically in 30 minutes to free server resources."),
        tags$p("Click anywhere or press any key to stay connected.")
      ),
      size = "m",
      easyClose = TRUE,
      footer = modalButton("Stay Connected")
    ))
  }, ignoreInit = TRUE)

  # Dismiss warning when user comes back
  observeEvent(input$session_timeout_reset, {
    removeModal()
  }, ignoreInit = TRUE)

  # Disconnect at 3 hours of inactivity
  observeEvent(input$session_timeout_kill, {
    showModal(modalDialog(
      title = tags$span(icon("power-off"), " Session Expired"),
      tags$div(class = "alert alert-danger mb-0",
        tags$p(icon("circle-xmark"),
               tags$strong(" Your session has been disconnected due to inactivity.")),
        tags$p("Please reload the page to start a new session.")
      ),
      size = "m",
      easyClose = FALSE,
      footer = NULL
    ))
    session$close()
  }, ignoreInit = TRUE)
  
  # Tab sync: sidebar nav -> hidden content panels
  observeEvent(input$suite_nav, { 
    updateTabsetPanel(session, "main_content_area", selected = input$suite_nav) 
  })
  
  # Shared reactive values (used by both modules)
  rv <- reactiveValues(
    log_text    = c("System Ready..."),
    data_source = NULL,   # NULL | "lf_dda" | "tmt" | "dia" — set by whichever wizard last completed
    # --- Async task cancellation relay ---
    # Set to TRUE by the global cancel button observer; each module server resets
    # it to FALSE after handling the cancellation for its own task.
    async_cancel_requested = FALSE,
    # --- Shared Proteome Discoverer CSV data ---
    pd_protein_headers  = NULL,   # character vector of protein-level column names
    pd_raw_lines        = NULL,   # all lines from the CSV
    pd_protein_line_idx = NULL,   # indices of protein data rows
    pd_ratio_cols       = NULL,   # recommended FC columns
    pd_pval_cols        = NULL,   # recommended P-value columns
    pd_other_cols       = NULL,   # remaining columns
    pd_filename         = NULL,   # original filename for display
    pd_gene_col         = NULL,   # user-selected gene ID column (when GN= missing)
    pd_has_gn           = TRUE,   # whether Description contains GN= tags
    # --- Single-Shot LFQ state (independent from LF DDA wizard) ---
    ss_pd_protein_headers  = NULL,
    ss_pd_raw_lines        = NULL,
    ss_pd_protein_line_idx = NULL,
    ss_pd_has_gn           = TRUE,
    ss_pd_gene_col         = NULL,
    ss_pd_filename         = NULL,
    # --- Label-Free input format ---
    lf_input_format     = "pd_nested",  # "pd_nested" | "maxquant"
    # --- Label-Free shared pipeline state (synced across DE/PCA/HM/TCL) ---
    lf_pipeline         = "global",      # "global" | "ip"
    lf_ip_sub           = "invariant",   # "invariant" | "no_norm" | "saint"
    lf_bait_protein     = NULL,          # gene symbol for bait protein normalization
    lf_bait_groups      = NULL,          # character vector of group names that contain the bait (NULL = all groups)
    lf_igg_conditions   = NULL,          # character vector of IgG/empty-ctrl condition names
    lf_factorial_design = NULL,          # NULL or list(enabled,f1_name,f1_levels,f2_name,f2_levels,assignments)
    lf_skip_impute      = FALSE,         # TRUE = skip AFT (Global) or LOD (IP)
    lf_cond_groups      = list(),        # named list: group_name -> c(col1, col2, ...) shared across DE/PCA/HM/TCL
    # --- Batch correction (PCA/Heatmap visualization + ProDA model when valid) ---
    lf_batch_groups     = NULL,           # named list: batch_label -> c(sample_names) or NULL
    lf_batch_corrected  = FALSE,          # TRUE when batch groups are assigned
    lf_batch_for_model  = FALSE,          # TRUE when batch is valid (non-confounded) for ProDA model
    lf_batch_filter_enabled = FALSE,      # TRUE when user opts in to batch-spanning protein filter
    # --- Label-Free summarization cache (shared across all label-free modules) ---
    lf_master_cache     = NULL,          # full ProDA summarization output list
    lf_master_cache_meta = NULL,         # list: pipeline/ip_sub/skip_impute/bait/igg/groups snapshot
    lf_cache_stale      = FALSE,         # TRUE when cache invalidated by config change (alerts user)
    # Shared-peptide handling for the peptide-aggregating route (B).
    # Always one of the three accepted settings: this is the input to the matrix
    # build, which rejects anything else. Routes that aggregate no peptides
    # ignore it and record NA_character_ in their OUTPUT instead.
    lf_shared_peptides  = "unique_only",
    # FALSE only when the input was normalized upstream: FragPipe with IonQuant
    # normalization on, per the required fragpipe.workflow. PD and MaxQuant: always
    # TRUE (OSMSS reads their raw columns). Reset with each new upload; effective on
    # the ProDA global route only (lf_apply_normalization_effective).
    lf_apply_normalization = TRUE,
    lf_fp_workflow          = NULL,     # fp_read_workflow() result for the FragPipe upload
    lf_fp_workflow_filename = NULL,
    lf_shared_pep_stats = NULL,          # counts from the matrix build, carried onto the cache
    lf_import_in_progress = FALSE,       # TRUE during summarization import (suppresses invalidation)
    # --- Label-Free input mode ---
    lf_input_mode          = "peptide",   # "peptide" (PD, default) | "maxquant"
    uniprot_lookup_pending    = FALSE,    # flag: UniProt gene name lookup offered after summarization
    # --- MaxQuant input mode state ---
    mq_evidence_path          = NULL,     # stable copy of uploaded evidence.txt
    mq_proteingroups_path     = NULL,     # stable copy of uploaded proteinGroups.txt
    mq_evidence_filename      = NULL,     # original filename for display
    mq_proteingroups_filename = NULL,     # original filename for display
    mq_raw_files              = NULL,     # unique Raw.file values from evidence.txt (for Step 4 group assignment)
    mq_pg_has_intensity       = FALSE,    # TRUE if proteinGroups.txt has raw 'Intensity <experiment>' columns (used)
    mq_mbr_on                 = NULL,     # TRUE/FALSE: match between runs detected in evidence.txt; NULL = unknown
    mq_n_match_features       = NULL,     # evidence.txt features transferred by matching
    mq_n_evidence_features    = NULL,     # evidence.txt features (rows)
    mq_remove_contaminants    = FALSE,    # TRUE = drop MaxQuant 'Potential contaminant' groups (default: keep)
    lf_build_meta             = NULL,     # aggregation route + route stats from the matrix build, carried onto the cache
    mq_annotation             = NULL,     # data.frame: Raw.file, Condition, BioReplicate, Run, IsotopeLabelType
    # --- MNAR routing + pairing state ---
    lf_expect_mnar            = FALSE,    # TRUE when user answered YES to MNAR question
    lf_wizard_q1_ip           = NULL,     # "yes" | "no" | NULL — Q1: did you do a pulldown/IP?
    lf_wizard_q2a_ctrl        = NULL,     # "trap" | "depleted" | "cond" | NULL — Q2a: control type (IP path)
    lf_wizard_q2b_abs         = NULL,     # "absent" | "quant" | NULL — Q2b: absent proteins? (non-IP path)
    lf_subject_pairing        = NULL,     # named char vec: sample_name -> subject_id, or NULL
    mq_rawfile_to_experiment  = NULL,     # named char vec: Raw.file -> Experiment (from evidence.txt)
    # --- Upload Wizard state ---
    wizard_completed      = FALSE,         # TRUE once wizard completed (gates post-wizard sidebar)
    wizard_modal_open     = FALSE,         # TRUE while wizard modal is on screen
    wizard_search_program = "pd",          # "pd" | "maxquant" — set by upload auto-detection
    wizard_data_types     = character(0),  # level quantified from: "protein" (Route A) | "peptide" (Route B) | "maxquant"
    wizard_step              = 3L,            # current wizard step (3-8; steps 1-2 retired, upload is auto-detected)
    wizard_reachable_steps   = c(3L, 4L, 6L, 7L),  # ordered sequence of internal step numbers (5=factorial skipped by default)
    pd_levels_available      = character(0),  # levels found in the PD CSV: subset of c("protein","peptide")
    mq_pg_headers            = NULL,          # proteinGroups.txt column names, make.names form (evidence cross-check)
    mq_pg_experiments        = NULL,          # experiment names from proteinGroups.txt quantity columns
    lf_wiz_detect_msg        = NULL,          # message for the last rejected upload (shown in step 3)
    wizard_has_batch_step    = FALSE,          # TRUE when IP+ProDA path (adds step 8)
    lf_de_engine          = "proda",        # "proda" | "saint" — derived on wizard completion (Phase 2)
    lf_pd_selected_fc     = NULL,          # FC column chosen in wizard Step 5 (PD Direct)
    lf_pd_selected_pval   = NULL,          # P-value column chosen in wizard Step 5 (PD Direct)
    # --- Shared TMT data (files shared between Volcano & Trends modules) ---
    tmt_files          = list(),   # parsed TMT file structs
    tmt_n_files        = 0L,
    tmt_conditions     = NULL,     # unique condition names (sorted)
    tmt_files_source   = NULL,     # "volcano" | "trends" — who last uploaded
    # Residue prediction state
    bind_data = NULL, bind_report = "", 
    pred_data = NULL, known_data = NULL, struct_df = NULL, 
    cons_df = NULL, acc_df = NULL, final_plot = NULL, 
    multisite_data = NULL, heatmap_data = NULL, pdb_path = NULL, 
    kin_fam_db = NULL, selected_res = NULL, render_trigger = 0,
    plotly_ready = FALSE,
    last_render_module = NULL, last_render_args = NULL,
    # Configurable paths
    path_kin_fam = PATH_KINASE_FAMILIES,
    path_st_mat = DEFAULT_ST_MATRIX,
    path_y_mat = DEFAULT_Y_MATRIX,
    path_st_fav = DEFAULT_ST_FAV,
    # Run-level cache
    cache_key = NULL, cache_seq = NULL, cache_cons = NULL, 
    cache_acc = NULL, cache_struct = NULL,
    # ── Master TMT cache (shared across Volcano / Heatmap / PCA / Trends) ──
    tmt_master_cache        = NULL,   # proteinSummarization() output object
    tmt_master_cache_design = NULL,   # design_option used to build cache
    tmt_master_cache_ann    = NULL,   # annotation_map used to build cache
    # ── TMT Wizard metadata (shared across TMT modules) ──
    tmt_wizard_completed      = FALSE,  # mirrors tmt_rv$wizard_completed
    tmt_analysis_route        = NULL,   # "limma_protein" | "limma_qf"
    tmt_wizard_search_program = NULL,   # "pd" | "other"
    tmt_wizard_data_level     = NULL,   # "protein" | "peptide"
    tmt_limma_cache           = NULL,   # Limma-route summarization cache
    # ── FragPipe DB session cache ──
    fpdb_experiment_db      = NULL,   # data.frame from MS_Experiment_DB (or NULL)
    fpdb_fragpipe_db        = NULL,   # data.frame from FragpipeDB (or NULL)
    fpdb_db_loaded          = FALSE,  # TRUE once user imports or creates a DB
    # ── DIA Label-Free state ──
    dia_wizard_step          = 1L,
    dia_wizard_modal_open    = FALSE,
    dia_wizard_completed     = FALSE,
    dia_software             = "spectronaut", # detected: "spectronaut" | "perseus" | "diann"
    dia_diann_log            = NULL,      # diann_read_log() result for DIA-NN's optional report.log.txt
    dia_diann_log_filename   = NULL,
    dia_raw_data             = NULL,      # list from parse_spectronaut_pivot()
    dia_protein_matrix       = NULL,      # numeric matrix: raw intensities (proteins x samples)
    dia_sample_names         = NULL,      # character vector: sample column names
    dia_gene_names           = NULL,      # named char: accession -> gene
    dia_accessions           = NULL,      # character vector: protein accessions
    dia_filename             = NULL,      # character: uploaded filename
    dia_cond_groups          = list(),    # named list: condition_name -> c(sample_names)
    dia_subject_pairing      = NULL,      # named char: sample_name -> subject_id, or NULL
    dia_factorial_design     = NULL,      # list or NULL
    dia_master_cache         = NULL,      # list: full summarization result
    dia_master_cache_meta    = NULL,      # list: pipeline metadata snapshot
    dia_cache_stale          = FALSE,
    dia_batch_groups         = NULL,      # named list: "Batch 1" -> c(sample_names), or NULL
    dia_batch_corrected      = FALSE,
    dia_batch_for_model      = FALSE,     # TRUE = non-confounded, add to design matrix
    # ── DIA PTM (Phospho) state ──
    dia_ptm_wizard_step          = 1L,
    dia_ptm_wizard_modal_open    = FALSE,
    dia_ptm_wizard_completed     = FALSE,
    dia_ptm_phospho_raw          = NULL,  # list from parse_spectronaut_phospho()
    dia_ptm_global_raw           = NULL,  # list from parse_spectronaut_pivot()
    dia_ptm_phospho_filename     = NULL,
    dia_ptm_global_filename      = NULL,
    dia_ptm_cond_groups          = list(),
    dia_ptm_subject_pairing      = NULL,
    dia_ptm_master_cache         = NULL,  # list from dia_ptm_summarize()
    dia_ptm_de_results           = NULL,  # named list of dia_ptm_limma_de() results,
                                          # keyed "<test> vs <ref>"; read by PTM Single
                                          # and the PTM Trends site gate.  The dia_
                                          # prefix is load-bearing: it is what makes
                                          # this field participate in
                                          # rv$reset_experiment_state().
    dia_ptm_software             = "spectronaut",  # "spectronaut" | "maxquant"
    dia_ptm_acquisition          = "dia",           # "dia" | "dda"
    dia_ptm_mode                 = "with_global",   # "with_global" | "ptm_only" (wizard step 1)

    # ── Prerun Experiments load triggers ──
    prerun_load_lf_bundle      = NULL,
    prerun_load_dia_bundle     = NULL,
    prerun_load_tmt_bundle     = NULL,
    prerun_load_dia_ptm_bundle = NULL,
    prerun_load_tmt_ptm_bundle = NULL,
    # ── TMT PTM wizard completion flag ──
    tmt_ptm_wizard_completed   = FALSE
  )

# ===========================================================================
  # 0. EXPERIMENT STATE RESET MANAGER (Snapshot Pattern)
  # ===========================================================================
  # Captures the clean initial state of `rv`. Provides a helper to dynamically
  # restore ONLY the analytical pipeline variables back to factory defaults.
  #
  # This is the single "a new dataset begins" reset. It is called by a saved
  # experiment import (server_prerun.R) and at the start of each wizard file
  # upload (LF, TMT, DIA), before the new data is written.
  #   keep - regexes of rv keys to leave alone; uploads pass the wizard's own
  #            in-progress answers plus RESET_KEEP_ON_UPLOAD (global.R).
  # Modules whose own reactiveValues hold results derived from the dataset
  # register with register_experiment_reset() (global.R); their hooks run at
  # the end of every reset, so no module can show the previous dataset.
  local({
    init_state <- isolate(reactiveValuesToList(rv))
    reset_hooks <- list()

    # Identify prefixes for all analysis/pipeline variables we want to wipe
    exp_prefixes <- c("^lf_", "^pd_", "^ss_pd_", "^mq_",
                      "^wizard_", "^tmt_", "^dia_")

    rv$on_experiment_reset <- function(fn) {
      reset_hooks[[length(reset_hooks) + 1L]] <<- fn
    }

    rv$reset_experiment_state <- function(keep = character(0)) {
      current_keys <- isolate(names(rv))

      for (k in current_keys) {
        # If the variable belongs to an experiment module. data_source has no
        # prefix but names the active dataset, so it is reset with them.
        if (identical(k, "data_source") ||
            any(vapply(exp_prefixes, function(p) grepl(p, k), logical(1)))) {

          # CRITICAL: Do NOT reset the trigger signals that are actively loading the new data!
          if (grepl("^prerun_load_", k)) next
          if (any(vapply(keep, function(p) grepl(p, k), logical(1)))) next

          if (k %in% names(init_state)) {
            # If it existed at startup (e.g., FALSE), restore it
            rv[[k]] <- init_state[[k]]
          } else {
            # If it didn't exist at startup (it was NULL), force it back to NULL
            rv[[k]] <- NULL
          }
        }
      }
      for (fn in reset_hooks) fn()
    }
  })

  # ── Shared PTM data store (single copy for all three PTM modules) ──
  rv$ptm_data <- reactiveValues(
    ptm_long        = NULL,
    protein_long    = NULL,
    gene_lookup     = NULL,
    parsed_files    = NULL,
    ptm_parsed_files  = NULL,
    prot_parsed_files = NULL,
    protein_headers = character(0),
    protein_df      = NULL,
    conditions      = NULL,
    annotation_map  = NULL,
    design_option   = NULL,
    selected_ptms   = NULL,
    n_files         = NULL,
    source          = NULL,
    prob_threshold  = NULL,
    data_format     = NULL,           # "pd_nested" | "flat_peptide"
    cond_order      = NULL            # user-set condition order, shared by the
                                      # PTM Trends and PTM Single reorder dialogs
  )

  # ===========================================================================
  # DIA/TMT PTM DE store lifecycle
  # ===========================================================================
  # rv$dia_ptm_de_results holds limma results keyed by contrast label
  # ("<test> vs <ref>").  That key is not unique across datasets: two experiments
  # with the same condition names produce the same key, so a result must never
  # outlive the summarization cache it was computed on.
  #
  # This is the single invalidation point.  It is here, rather than at each site
  # that writes the cache, because there are already six such sites across
  # server_dia_ptm_upload.R (three software-switch resets, Summarize, the prerun
  # bundle loader) and server_ptm_de.R (the k-NN memo write) - and a seventh added
  # later would silently miss a per-site clear.  Watching the cache covers them all.
  #
  # The comparison is on ptm_cache_fingerprint(), not on the event: the k-NN memo
  # write mutates the cache without changing the dataset, and reacting to that
  # would wipe the store immediately after the DE run that filled it.
  #
  # priority is raised so that if a cache change and a consumer read ever land in
  # the same flush, the clear happens first.
  local({
    last_fp <- NULL
    observeEvent(rv$dia_ptm_master_cache, priority = 1000L, {
      fp <- ptm_cache_fingerprint(rv$dia_ptm_master_cache)
      if (!identical(fp, last_fp)) {
        if (!is.null(rv$dia_ptm_de_results)) {
          logger(paste0("[PTM] Summarization cache changed - cleared ",
                        length(rv$dia_ptm_de_results),
                        " stored DE result(s); they were computed on the previous cache."))
          rv$dia_ptm_de_results <- NULL
        }
        last_fp <<- fp
      }
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
  })

  # ===========================================================================
  # Provenance Ledger — centralized state-based pipeline tracking (v1.0)
  # Used to auto-generate publication-ready methods text.
  # Analytical modules overwrite their dedicated node atomically.
  # ===========================================================================
  rv$provenance <- reactiveValues(
    schema_version = "1.0",

    # 1. Environment Metadata (Populated immediately on app startup)
    env = list(
      osmss_version      = OSMSS_VERSION,
      r_version          = as.character(getRversion()),
      msstats_version    = tryCatch(as.character(packageVersion("MSstats")), error = function(e) NA_character_),
      limma_version      = tryCatch(as.character(packageVersion("limma")),   error = function(e) NA_character_),
      proda_version      = tryCatch(as.character(packageVersion("proDA")),   error = function(e) NA_character_),
      qfeatures_version  = tryCatch(as.character(packageVersion("QFeatures")), error = function(e) NA_character_),
      msigdbr_version    = tryCatch(as.character(packageVersion("msigdbr")), error = function(e) NA_character_)
    ),

    # 2. Global Input Data (Populated upon successful file upload/parsing)
    global = list(
      input_filename         = character(),
      total_initial_proteins = numeric(),
      input_software_source  = character(),
      search_program         = character(),  # "Proteome Discoverer" | "Other"
      data_types             = character(),  # e.g. c("protein","peptide")
      wizard_version         = "1.0"
    ),

    # 3. Differential Expression State (Overwritten by server_de.R)
    de = list(
      analysis_mode        = character(),
      starting_proteins    = numeric(),
      dropped_proteins     = numeric(),
      drop_reason          = character(),   # human-readable gate description
      imputation_method    = character(),
      fdr_threshold        = numeric(),
      logfc_threshold      = numeric(),
      software_used        = character(),
      lod_imputation       = FALSE,
      # MSstats-specific parameters
      normalization        = character(),   # "FALSE" | "equalizeMedians" | "invariantMedian"
      summary_method       = character(),   # e.g. "TMP"
      censored_int         = character(),   # e.g. "NA"
      mb_impute            = logical(),     # model-based imputation flag
      experiment_type      = character(),   # "ip" | "ip_bait" | "ip_invariant" | "total"
      # Imputation-specific parameters
      norm_strategy        = character(),   # "within-condition median" | "global median"
      pvalue_adjustment    = character(),   # e.g. "BH"
      # Phase 2 — ProDA engine fields
      de_engine            = character(),   # "proda" | "saint"
      normalization_method = character(),   # "invariant_median" | "none" | "equalizeMedians"
      aggregation_method   = character(),   # "robustSummary" | "MSstats_TMP"
      proda_version        = character(),   # packageVersion("proDA")
      qfeatures_version    = character()    # packageVersion("QFeatures")
    ),

    # 4. PCA State (Overwritten by server_pca.R)
    pca = list(
      pipeline           = character(),  # "global" | "ip"
      n_samples          = numeric(),
      n_proteins         = numeric(),    # proteins after NA filtering
      n_proteins_input   = numeric(),    # proteins before filtering
      n_proteins_dropped = numeric(),
      normalization      = character(),  # "equalizeMedians" | "invariantMedian" | "FALSE" | "pd_normalized"
      summary_method     = character(),  # "TMP" or character()
      mb_impute          = logical(),
      lod_imputation     = logical(),
      ebn_ref_proteins   = character(),
      centering          = logical(),
      scaling            = logical(),
      pc1_var            = numeric(),
      pc2_var            = numeric(),
      condition_names    = character(),
      software_used      = character()
    ),

    # 5. Trend Clusters State (Overwritten by server_trend_clusters.R)
    trend_clusters = list(
      pipeline             = character(),  # "global" | "ip"
      n_proteins_input     = numeric(),
      n_proteins_clustered = numeric(),
      n_clusters           = numeric(),
      clustering_method    = character(),  # "fuzzy c-means (e1071::cmeans)"
      distance_metric      = character(),  # "euclidean"
      fuzzifier            = numeric(),    # computed m_est
      normalization        = character(),
      summary_method       = character(),
      mb_impute            = logical(),
      lod_imputation       = logical(),
      ebn_ref_proteins     = character(),
      zscore_transform     = logical(),
      condition_names      = character(),
      software_used        = character()
    ),

    # 6. Heatmap State (Overwritten by server_heatmap.R)
    heatmap = list(
      pipeline             = character(),  # "ip" | "global" | "pd_norm"
      n_proteins_input     = numeric(),
      n_proteins_displayed = numeric(),
      selection_method     = character(),  # "top_variance" | "gene_filter"
      top_n                = numeric(),
      gene_filter          = character(),
      avg_by_condition     = logical(),
      normalization        = character(),
      summary_method       = character(),
      mb_impute            = logical(),
      lod_imputation       = logical(),
      ebn_ref_proteins     = character(),
      row_scaling          = character(),  # "z-score"
      row_clustering       = character(),  # "hierarchical (ward.D2, euclidean)"
      col_clustering       = character(),  # "hierarchical..." | "grouped" | "manual"
      condition_names      = character(),
      software_used        = character()
    ),

    # 7. Rank Fold Change State (Overwritten by server_rank_fc.R)
    rank_fc = list(
      n_proteins_input    = numeric(),   # nrow(protein_df) before NA/Inf filtering
      n_proteins_plotted  = numeric(),   # nrow(plot_df) after valid-row filtering
      n_imputed           = numeric(),   # sum(was_imputed) — 0 when impute OFF
      condition_a         = character(), # column name selected as Condition A
      condition_b         = character(), # column name selected as Condition B
      fc_method           = character(), # always "log2(B/A)"
      imputation_enabled  = logical(),   # TRUE/FALSE from input$rfc_impute_toggle
      imputation_method   = character(), # "Perseus-style random normal" or character()
      imp_downshift       = numeric(),   # downshift parameter (default 1.8), or numeric()
      imp_width           = numeric(),   # width parameter (default 0.3), or numeric()
      software_used       = character()  # "plotly"
    ),

    # 8. Trends Analysis State (Overwritten by server_trends.R)
    trends = list(
      n_proteins_input     = numeric(),   # nrow(protein_df) before finite-row filtering
      n_proteins_clustered = numeric(),   # nrow(zscore_mat) after filtering
      n_clusters           = numeric(),   # k
      clustering_method    = character(), # "fuzzy c-means (e1071::cmeans)"
      distance_metric      = character(), # "euclidean"
      fuzzifier            = numeric(),   # round(m_est, 3)
      zscore_transform     = logical(),   # always TRUE
      condition_names      = character(), # col_labels (the short condition labels)
      software_used        = character()  # "e1071"
    ),

    # 9. TMT Pairwise Volcano State (Overwritten by server_tmt_volcano.R)
    tmt_volcano = list(
      analysis_mode              = character(),   # "hybrid_limma" | "limma_protein" | "limma_qf"
      design_option              = numeric(),     # 1 = single-plex, 2 = bridge, 3 = balanced
      comparison                 = character(),   # e.g. "GroupA vs GroupB"
      n_proteins                 = numeric(),
      n_significant              = numeric(),     # sum(adj.pvalue < 0.05)
      global_norm                = logical(),     # always TRUE
      reference_norm             = logical(),     # TRUE only for design_option == 2
      mb_impute                  = logical(),     # always FALSE
      censored_quantile          = numeric(),     # 0.999
      pvalue_adjustment          = character(),   # "BH"
      software_used              = character(),   # "limma" | "limma + QFeatures"
      # hybrid_limma fields
      n_proteins_before_filter   = numeric(),
      n_proteins_filtered        = numeric(),
      # Limma-specific fields
      imputation_method          = character(),   # "k-NN (impute::impute.knn)" | character()
      imputation_k               = numeric(),     # adaptive k value
      missingness_threshold      = numeric(),     # 0.50
      n_proteins_dropped_missing = numeric(),     # proteins dropped for <50% coverage
      interference_filtered      = logical(),     # TRUE if interference filter applied
      interference_threshold     = numeric(),     # 30 (percent)
      aggregation_method         = character(),   # "robustSummary" | character()
      n_peptides_input           = numeric(),     # for limma_qf
      ebayes_trend               = logical(),     # TRUE
      data_level                 = character()    # "protein" | "peptide" | "psm"
    ),

    # 10. TMT Trend Clusters State (Overwritten by server_tmt_trends.R)
    tmt_trends = list(
      n_proteins_input     = numeric(),
      n_proteins_clustered = numeric(),
      n_clusters           = numeric(),
      clustering_method    = character(),  # "fuzzy c-means (e1071::cmeans)"
      distance_metric      = character(),  # "euclidean"
      fuzzifier            = numeric(),
      missing_mode         = character(),  # "lod" | "omit"
      n_omitted            = numeric(),
      n_imputed            = numeric(),
      lod_value            = numeric(),
      batch_corrected      = logical(),
      global_norm          = logical(),
      reference_norm       = logical(),
      zscore_transform     = logical(),    # always TRUE
      condition_names      = character(),
      software_used        = character()   # "e1071"
    ),

    # 11. TMT PCA State (Overwritten by server_tmt_pca.R)
    tmt_pca = list(
      n_samples            = numeric(),
      n_proteins           = numeric(),
      design_option        = numeric(),
      global_norm          = logical(),
      reference_norm       = logical(),
      mb_impute            = logical(),    # always FALSE
      batch_corrected      = logical(),
      centering            = logical(),    # always TRUE
      scaling              = logical(),    # always FALSE
      pc1_var              = numeric(),
      pc2_var              = numeric(),
      condition_names      = character(),
      software_used        = character()   # "prcomp (base R)"
    ),

    # 12. TMT Heatmap State (Overwritten by server_tmt_heatmap.R)
    tmt_heatmap = list(
      n_proteins_input     = numeric(),
      n_proteins_displayed = numeric(),
      missing_mode         = character(),  # "lod" | "omit"
      n_omitted            = numeric(),
      n_imputed            = numeric(),
      lod_value            = numeric(),
      selection_method     = character(),  # "top_variance" | "gene_filter"
      top_n                = numeric(),
      gene_filter          = character(),
      avg_by_condition     = logical(),
      global_norm          = logical(),
      reference_norm       = logical(),
      row_scaling          = character(),  # "z-score"
      row_clustering       = character(),  # "hierarchical (ward.D2, euclidean)"
      col_clustering       = character(),
      condition_names      = character(),
      software_used        = character()   # "plotly"
    ),

    # 13. PTM Differential Expression State (Overwritten by server_ptm_de.R)
    ptm_de = list(
      analysis_mode            = character(),  # "hybrid_limma"
      design_option            = numeric(),    # 1 | 2 | 3
      comparison               = character(),  # "GroupA vs GroupB"
      n_ptm_sites              = numeric(),    # total PTM sites in result
      n_significant            = numeric(),    # adj.pvalue < 0.05
      n_ptm_sites_unadjusted   = numeric(),    # rows in PTM.Model
      n_significant_unadjusted = numeric(),    # adj.p < 0.05 in PTM.Model
      n_proteins               = numeric(),    # rows in PROTEIN.Model
      prob_threshold           = numeric(),    # localization probability % (e.g., 75)
      selected_ptms            = character(),  # e.g., "Phospho"
      global_norm              = logical(),    # always TRUE
      reference_norm           = logical(),    # TRUE if design_option == 2
      mb_impute                = logical(),    # always FALSE
      pvalue_adjustment        = character(),  # "BH"
      batch_corrected          = logical(),    # TRUE if hybrid_limma
      n_plexes                 = numeric(),    # number of TMT plexes
      condition_names          = character(),  # all condition names
      software_used            = character()   # "limma"
    ),

    # 14. PTM Trend Clusters State (Overwritten by server_ptm_trends.R)
    ptm_trends = list(
      n_ptm_sites_input        = numeric(),    # before filtering
      n_ptm_sites_clustered    = numeric(),    # after filtering + variance purge
      n_clusters               = numeric(),    # k value
      clustering_method        = character(),  # "fuzzy c-means (e1071::cmeans)"
      distance_metric          = character(),  # "euclidean"
      fuzzifier                = numeric(),    # m value
      missing_mode             = character(),  # "omit" | "lod"
      n_omitted                = numeric(),
      n_imputed                = numeric(),
      lod_value                = numeric(),    # 1st percentile (LOD mode only)
      batch_corrected          = logical(),
      global_norm              = logical(),    # always TRUE
      reference_norm           = logical(),    # TRUE if design_option == 2
      zscore_transform         = logical(),    # always TRUE
      prob_threshold           = numeric(),    # localization probability %
      selected_ptms            = character(),  # e.g., "Phospho"
      n_plexes                 = numeric(),
      condition_names          = character(),
      software_used            = character()   # "e1071"
    ),

    # 15. PTM Single Inspector State (Overwritten by server_ptm_single.R)
    ptm_single = list(
      ptm_id                   = character(),  # e.g., "Q9Y6K9_S473_Phospho"
      gene                     = character(),  # resolved gene name
      site                     = character(),  # e.g., "S473"
      modification             = character(),  # e.g., "Phospho"
      reference_condition      = character(),  # user-selected reference
      n_conditions             = numeric(),    # total conditions
      n_contrasts              = numeric(),    # number of non-ref conditions compared
      missing_conds            = character(),  # conditions where PTM had no data
      has_protein              = logical(),    # whether parent protein data exists
      design_option            = numeric(),    # 1 | 2 | 3
      prob_threshold           = numeric(),
      selected_ptms            = character(),
      global_norm              = logical(),    # always TRUE
      reference_norm           = logical(),
      batch_corrected          = logical(),
      n_plexes                 = numeric(),
      condition_names          = character(),
      software_used            = character()   # "limma"
    ),

    # 16. DIA Global Input State (Populated upon wizard completion)
    dia_global = list(
      input_file           = character(),
      software             = character(),   # "Spectronaut"
      data_level           = character(),   # "protein"
      n_proteins_input     = numeric(),
      n_proteins_processed = numeric(),
      n_samples            = numeric(),
      n_conditions         = numeric(),
      condition_names      = character(),
      normalization        = character(),   # "column-wise median centering (log2)"
      engine               = character(),   # "limma (eBayes: trend=TRUE, robust=TRUE)"
      pairing              = logical(),
      factorial            = logical(),
      timestamp            = character()
    ),

    # 17. DIA Differential Expression State (Overwritten by server_dia_de.R)
    dia_de = list(
      comparison           = character(),
      engine               = character(),   # "limma"
      ebayes_trend         = logical(),     # TRUE
      ebayes_robust        = logical(),     # TRUE
      n_tested             = numeric(),
      n_significant        = numeric(),
      fc_threshold         = numeric(),
      p_threshold          = numeric(),
      imputation           = character(),   # "k-NN"
      missingness_filter   = character(),   # ">=2 non-NA in >=1 group"
      batch_in_model       = logical(),
      pairing_method       = character(),
      factorial            = logical(),
      timestamp            = character()
    ),

    # 18. DIA PCA State (Overwritten by server_dia_pca.R)
    dia_pca = list(
      n_proteins           = numeric(),
      n_samples            = numeric(),
      batch_corrected      = logical(),
      centering            = logical(),     # always TRUE
      scaling              = logical(),     # always FALSE
      pc1_var              = numeric(),
      pc2_var              = numeric(),
      condition_names      = character(),
      software_used        = character()    # "prcomp (base R)"
    ),

    # 19. DIA Heatmap State (Overwritten by server_dia_heatmap.R)
    dia_heatmap = list(
      n_proteins_input     = numeric(),
      n_proteins_displayed = numeric(),
      missing_mode         = character(),   # "lod" | "omit"
      row_scaling          = character(),   # "z-score" | "abundance"
      batch_corrected      = logical(),
      condition_names      = character(),
      software_used        = character()    # "plotly"
    ),

    # 20. DIA Trend Clusters State (Overwritten by server_dia_trends.R)
    dia_trends = list(
      n_proteins_input     = numeric(),
      n_proteins_clustered = numeric(),
      n_clusters           = numeric(),
      clustering_method    = character(),   # "fuzzy c-means (e1071::cmeans)"
      distance_metric      = character(),   # "euclidean"
      fuzzifier            = numeric(),
      zscore_transform     = logical(),     # always TRUE
      condition_names      = character(),
      software_used        = character()    # "e1071"
    ),

    # 21. DIA PTM State (Overwritten by server_dia_ptm_upload.R / PTM modules)
    dia_ptm = list(
      phospho_file        = character(),
      global_file         = character(),
      software            = character(),
      data_level          = character(),
      n_sites_input       = numeric(),
      n_sites_output      = numeric(),
      n_sites_adjusted    = numeric(),
      n_shared_expanded   = numeric(),
      prob_threshold      = numeric(),
      n_samples           = numeric(),
      n_conditions        = numeric(),
      condition_names     = character(),
      normalization       = character(),
      stoichiometry       = character(),
      engine              = character(),
      pairing             = logical(),
      timestamp           = character()
    )
  )
  
  # Shared logger (capped at 500 lines) — also persists to daily log file
  # NOTE: isolate() prevents reactive dependency when called inside render contexts
  # Accepts optional level/module overrides; auto-detects from legacy prefixes.
  logger <- function(txt, level = NULL, module = NULL) { 
    timestamp <- format(Sys.time(), "%H:%M:%S")
    rv$log_text <- c(isolate(rv$log_text), paste0("[", timestamp, "] ", txt))
    if(length(rv$log_text) > 500) {
      rv$log_text <- tail(rv$log_text, 500)
    }
    # Persist to file — auto-parse level/module from message if not supplied
    parsed <- osmss_parse_log_msg(txt)
    osmss_file_log(
      msg        = parsed$msg,
      level      = if (!is.null(level))  level  else parsed$level,
      module     = if (!is.null(module)) module else parsed$module,
      session_id = .sid
    )
  }

  # --- Session lifecycle logging ---
  session$onSessionEnded(function() {
    osmss_file_log("Session ended (client disconnected)", level = "INFO",
                   module = "Server", session_id = .sid)
  })

  # Global error handler — catches unhandled errors inside this session
  options(shiny.error = function() {
    err_msg <- geterrmessage()
    osmss_file_log(paste0("UNHANDLED ERROR: ", err_msg), level = "ERROR",
                   module = "Server", session_id = .sid)
  })
  
  # ===========================================================================
  # Shared CSV parser for Proteome Discoverer modules
  # Both fileInputs (de_csv_file, pca_csv_file, etc.) trigger the same parse
  # ===========================================================================

  parse_pd_csv <- function(file_info) {
    csv_path <- file_info$datapath
    orig_name <- file_info$name
    if (!file.exists(csv_path)) {
      showNotification("Uploaded file could not be read.", type = "error")
      return()
    }
    
    logger(paste0("PD CSV: Loading '", orig_name, "'..."))
    
    withProgress(message = "Parsing Discoverer CSV...", value = 0, {
      all_lines <- readLines(csv_path, warn = FALSE)
      # PD wraps GO/Pfam annotations across lines inside quoted fields. Rejoin
      # them first, or the protein rows (picked one physical line each) leave a
      # quote open and every column after it - Gene Symbol included - is lost.
      # Same step the TMT uploader applies (server_tmt_upload.R).
      all_lines <- collapse_multiline_csv(all_lines)
      incProgress(0.3, detail = "Identifying protein rows...")
      
      header_line <- all_lines[1]
      protein_headers <- as.character(
        read.csv(text = header_line, header = FALSE,
                 stringsAsFactors = FALSE, check.names = FALSE)
      )
      protein_headers <- protein_headers[!is.na(protein_headers) &
                                           nchar(trimws(protein_headers)) > 0]
      
      # Protein rows: FDR Confidence column holds High/Medium/Low, fields counted
      # quote-aware (pd_protein_row_idx in global.R, shared with the TMT parser)
      protein_idx <- pd_protein_row_idx(all_lines, header_line)
      
      incProgress(0.6, detail = "Building column selections...")
      
      ratio_cols <- grep("^Abundance Ratio:", protein_headers, value = TRUE)
      ratio_cols <- ratio_cols[!grepl("Variability|P-Value", ratio_cols)]
      pval_cols  <- grep("Abundance Ratio Adj\\. P-Value", protein_headers, value = TRUE)
      other_cols <- setdiff(protein_headers, c(ratio_cols, pval_cols))
      
      # Store in shared rv
      rv$pd_protein_headers  <- protein_headers
      rv$pd_raw_lines        <- all_lines
      rv$pd_protein_line_idx <- protein_idx
      rv$pd_ratio_cols       <- ratio_cols
      rv$pd_pval_cols        <- pval_cols
      rv$pd_other_cols       <- other_cols
      rv$pd_filename         <- orig_name
      rv$pd_gene_col         <- NULL   # reset on new upload
      rv$lf_cond_groups      <- list() # reset groups on new upload
      rv$lf_bait_groups      <- NULL   # reset bait groups on new upload
      rv$lf_batch_groups     <- NULL   # reset batch on new upload
      rv$lf_batch_corrected  <- FALSE
      rv$lf_batch_for_model  <- FALSE  # reset batch-in-model flag on new upload
      rv$lf_factorial_design <- NULL   # reset factorial design on new upload

      # --- Provenance: populate global input node ---
      # Detect input software source from column naming patterns
      .prov_source <- if (any(grepl("^Abundance Ratio:", protein_headers))) {
        "Proteome Discoverer"
      } else if (any(grepl("^MaxLFQ", protein_headers, ignore.case = TRUE))) {
        "FragPipe"
      } else {
        "Unknown"
      }
      rv$provenance$global <- list(
        input_filename         = orig_name,
        total_initial_proteins = length(protein_idx),
        input_software_source  = .prov_source
      )
      
      # --- Auto-detect gene ID source via shared helper ---
      incProgress(0.8, detail = "Checking gene ID format...")
      protein_lines_sample <- all_lines[protein_idx[seq_len(min(50, length(protein_idx)))]]
      sample_block <- paste(c(header_line, protein_lines_sample), collapse = "\n")
      sample_df <- tryCatch(
        read.csv(text = sample_block, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, na.strings = c("", "NA")),
        error = function(e) { logger(paste0("!!! [Server] CSV sample parse error: ", e$message)); NULL }
      )
      
      gd <- detect_gene_column(protein_headers, sample_df = sample_df,
                               rv = rv, logger = logger)
      if (gd$source == "column") {
        rv$pd_has_gn   <- FALSE
        rv$pd_gene_col <- gd$column
      } else if (gd$source == "gn") {
        rv$pd_has_gn <- TRUE
      } else {
        rv$pd_has_gn <- FALSE
      }
      
      incProgress(1, detail = "Done!")
    })
    
    logger(paste0("PD CSV: Found ", length(protein_idx), " proteins, ",
                  length(protein_headers), " columns from '", orig_name, "'"))
    showNotification(
      paste0("Loaded ", length(protein_idx), " proteins from ", orig_name),
      type = "message")
    
    # If neither Gene Symbol column nor GN= tags found, prompt user
    if (!rv$pd_has_gn && is.null(rv$pd_gene_col)) {
      showModal(build_gene_col_modal("pd_gene_col_radio",
                                     "pd_gene_col_confirm",
                                     rv$pd_protein_headers))
    }
  }
  






  observeEvent(input$rfc_csv_file, {
    req(input$rfc_csv_file)
    parse_pd_csv(input$rfc_csv_file)
    # Populate Single-Shot shared state from the freshly parsed rv$pd_* fields
    rv$ss_pd_protein_headers  <- rv$pd_protein_headers
    rv$ss_pd_raw_lines        <- rv$pd_raw_lines
    rv$ss_pd_protein_line_idx <- rv$pd_protein_line_idx
    rv$ss_pd_has_gn           <- rv$pd_has_gn
    rv$ss_pd_gene_col         <- rv$pd_gene_col
    rv$ss_pd_filename         <- rv$pd_filename
  })

  observeEvent(input$trd_csv_file, {
    req(input$trd_csv_file)
    parse_pd_csv(input$trd_csv_file)
    # Populate Single-Shot shared state from the freshly parsed rv$pd_* fields
    rv$ss_pd_protein_headers  <- rv$pd_protein_headers
    rv$ss_pd_raw_lines        <- rv$pd_raw_lines
    rv$ss_pd_protein_line_idx <- rv$pd_protein_line_idx
    rv$ss_pd_has_gn           <- rv$pd_has_gn
    rv$ss_pd_gene_col         <- rv$pd_gene_col
    rv$ss_pd_filename         <- rv$pd_filename
  })


  
  # ===========================================================================
  # Cache invalidation — any change to pipeline config or data clears the
  # summarization cache so the user must re-summarize before running.
  # ===========================================================================
  observeEvent({
    rv$lf_pipeline; rv$lf_ip_sub; rv$lf_bait_protein
    rv$lf_igg_conditions; rv$lf_skip_impute
    rv$lf_cond_groups; rv$pd_raw_lines
    rv$lf_shared_peptides; rv$lf_apply_normalization
  }, {
    if (isTRUE(rv$lf_import_in_progress)) return()
    # Skip if the change was triggered by a TMT or PTM upload (pd_raw_lines is
    # shared between DDA and TMT; TMT writes it during its own file upload).
    active_tab <- isolate(input$suite_nav)
    if (!is.null(active_tab) && grepl("^TMT|^PTM", active_tab)) return()
    if (!is.null(rv$lf_master_cache)) {
      rv$lf_master_cache      <- NULL
      rv$lf_master_cache_meta <- NULL
      rv$lf_cache_stale       <- TRUE
      # The wizard's "Done! N proteins processed." line described the cleared cache
      output$wiz_summarize_status <- renderUI({ NULL })
      logger("LF Cache: invalidated (pipeline config or data changed)")
      # If the wizard modal is currently open, use a notification instead of
      # showModal — otherwise the wizard gets replaced and the user is forced
      # to re-open it.
      if (isTRUE(rv$wizard_modal_open)) {
        showNotification(
          paste0("Pipeline settings changed \u2014 previous summarization cache cleared. ",
                 "You will need to re-summarize."),
          type = "warning", duration = 6
        )
      } else {
        showModal(modalDialog(
          title = tags$span(icon("triangle-exclamation"), " Pipeline Settings Changed"),
          tags$div(class = "alert alert-warning mb-0",
            tags$p(
              tags$b("Your analysis pipeline or data has changed."),
              tags$br(),
              "The previous summarization cache has been cleared."
            ),
            tags$p(
              "You can either:",
              tags$ul(
                tags$li(tags$b("Re-summarize"), " using the Summarize button in any module before running."),
                tags$li(tags$b("Run directly"), " \u2014 the module will automatically re-summarize with the new settings",
                       " (this may take a few minutes).")
              )
            )
          ),
          size = "m", easyClose = FALSE,
          footer = modalButton("Got it")
        ))
      }
    }
  }, ignoreInit = TRUE)


  # ── Data Status Banners (for simplified module sidebars) ──
  .make_status_banner <- function() {
    src <- rv$data_source
    # DIA path
    if (identical(src, "dia")) {
      if (!is.null(rv$dia_master_cache)) {
        m <- rv$dia_master_cache_meta
        return(tags$div(class = "alert alert-success py-2 px-3 small mb-2",
          icon("circle-check"),
          tags$b(" Summarized (DIA). "),
          if (!is.null(m$n_proteins)) paste0(m$n_proteins, " proteins")
        ))
      } else {
        return(tags$div(class = "alert alert-info py-2 px-3 small mb-2",
          icon("circle-info"),
          tags$b(" DIA data loading."), " Go to ", tags$b("Upload Data"), " to summarize."
        ))
      }
    }
    # TMT path
    if (identical(src, "tmt")) {
      if (!is.null(rv$tmt_master_cache)) {
        return(tags$div(class = "alert alert-success py-2 px-3 small mb-2",
          icon("circle-check"),
          tags$b(" Summarized (TMT). "),
          if (!is.null(rv$tmt_master_cache$ProteinLevelData))
            paste0(length(unique(rv$tmt_master_cache$ProteinLevelData$Protein)), " proteins")
        ))
      } else if (!is.null(rv$tmt_limma_cache)) {
        n_prots   <- nrow(rv$tmt_limma_cache$protein_matrix)
        route_lbl <- switch(rv$tmt_limma_cache$route %||% "limma_protein",
          limma_protein = "Limma (Protein-Level)",
          limma_qf      = "QFeatures + Limma",
          "Limma"
        )
        return(tags$div(class = "alert alert-success py-2 px-3 small mb-2",
          icon("circle-check"),
          tags$b(" Summarized (TMT). "),
          paste0(n_prots, " proteins (", route_lbl, ")")
        ))
      } else {
        return(tags$div(class = "alert alert-info py-2 px-3 small mb-2",
          icon("circle-info"),
          tags$b(" TMT data loading."), " Go to ", tags$b("Upload Data"), " to summarize."
        ))
      }
    }
    # LF DDA path (original + Prerun Experiments import)
    # Show "No data loaded" only when neither a raw file nor a prerun cache exists.
    if (is.null(rv$pd_raw_lines) && is.null(rv$mq_proteingroups_path) &&
        is.null(rv$lf_master_cache)) {
      tags$div(class = "alert alert-warning py-2 px-3 small mb-2",
        icon("triangle-exclamation"),
        tags$b(" No data loaded."),
        " Go to ", tags$b("Upload Data"), " to upload and summarize."
      )
    } else if (is.null(rv$lf_master_cache)) {
      tags$div(class = "alert alert-info py-2 px-3 small mb-2",
        icon("circle-info"),
        tags$b(" Data loaded."),
        " Go to ", tags$b("Upload Data"), " to summarize before running analysis."
      )
    } else {
      m <- rv$lf_master_cache_meta
      tags$div(class = "alert alert-success py-2 px-3 small mb-2",
        icon("circle-check"),
        tags$b(" Summarized. "),
        if (!is.null(m$n_proteins)) paste0(m$n_proteins, " proteins, "),
        if (!is.null(m$pipeline)) paste0("Pipeline: ", m$pipeline)
      )
    }
  }
  output$de_data_status_banner  <- renderUI(.make_status_banner())
  output$pca_data_status_banner <- renderUI(.make_status_banner())
  output$hm_data_status_banner  <- renderUI(.make_status_banner())
  output$tcl_data_status_banner <- renderUI(.make_status_banner())
  output$cov_data_status_banner <- renderUI(.make_status_banner())

  output$tmt_upload_main_status <- renderUI({
    has_msstats <- !is.null(rv$tmt_master_cache)
    has_limma   <- !is.null(rv$tmt_limma_cache)
    if (!has_msstats && !has_limma) {
      tags$div(class = "text-muted",
        tags$p("Use the sidebar to upload TMT files, configure your experimental ",
               "design, and summarize your data.")
      )
    } else if (has_limma) {
      n_prots <- nrow(rv$tmt_limma_cache$protein_matrix)
      tags$div(class = "alert alert-success mt-3",
        icon("circle-check"),
        tags$b(" Summarization complete! "),
        paste0(n_prots, " proteins ready for analysis. ",
               "Navigate to any TMT analysis module to begin.")
      )
    } else {
      n_prots <- length(unique(rv$tmt_master_cache$ProteinLevelData$Protein))
      tags$div(class = "alert alert-success mt-3",
        icon("circle-check"),
        tags$b(" Summarization complete! "),
        paste0(n_prots, " proteins ready for analysis. ",
               "Navigate to any TMT analysis module to begin.")
      )
    }
  })

  # ── PTM Data Status Banners ──
  .make_ptm_status_banner <- function() {
    if (is.null(rv$dia_ptm_master_cache)) {
      tags$div(class = "alert alert-warning py-2 px-3 small mb-2",
        icon("triangle-exclamation"),
        tags$b(" No PTM data loaded."),
        " Complete the PTM Upload Wizard to begin."
      )
    } else {
      n_ptms <- rv$dia_ptm_master_cache$n_sites_output %||% "?"
      tags$div(class = "alert alert-success py-2 px-3 small mb-2",
        icon("circle-check"),
        tags$b(" PTM data ready. "),
        paste0(n_ptms, " PTM sites summarized.")
      )
    }
  }
  output$ptm_de_data_status_banner  <- renderUI(.make_ptm_status_banner())
  output$ptrnd_data_status_banner   <- renderUI(.make_ptm_status_banner())
  output$psg_data_status_banner     <- renderUI(.make_ptm_status_banner())

  output$upload_status_display <- renderUI({
    if (is.null(rv$pd_raw_lines) && is.null(rv$mq_proteingroups_path)) {
      tags$div(class = "text-muted text-center py-5",
        icon("cloud-arrow-up", class = "fa-3x mb-3"),
        tags$h5("No data uploaded yet"),
        tags$p("Use the sidebar to upload your data file.")
      )
    } else if (is.null(rv$lf_master_cache)) {
      tags$div(class = "alert alert-info",
        icon("circle-info"),
        tags$b(" Data loaded. "),
        "Use the sidebar to configure pipeline and summarize."
      )
    } else {
      m <- rv$lf_master_cache_meta
      tags$div(class = "alert alert-success",
        icon("circle-check"),
        tags$b(" Summarization complete! "),
        if (!is.null(m$n_proteins)) paste0(m$n_proteins, " proteins ready for analysis.")
      )
    }
  })

  output$upload_data_preview <- renderUI({
    if (!is.null(rv$pd_protein_headers)) {
      tags$div(
        tags$h6("Detected Columns:"),
        tags$p(class = "small text-muted",
          paste(head(rv$pd_protein_headers, 20), collapse = ", "),
          if (length(rv$pd_protein_headers) > 20) "..." else ""
        )
      )
    }
  })

  # Handle Gene ID column selection (when GN= is missing from Description)
  observeEvent(input$pd_gene_col_confirm, {
    sel <- input$pd_gene_col_radio
    if (!is.null(sel) && sel != "") {
      rv$pd_gene_col <- sel
      logger(paste0("PD CSV: Gene ID column set to '", sel, "'"))
      showNotification(paste0("Gene ID column set to: ", sel), type = "message")
    }
    removeModal()
  })
  
  # Clone for Rank FC sidebar
  output$rfc_file_status <- renderUI({
    if (is.null(rv$pd_filename)) {
      tags$p(class = "text-muted small", icon("circle-xmark"), " No file loaded")
    } else {
      tags$p(class = "small text-success",
             icon("circle-check"), " ", rv$pd_filename,
             tags$br(),
             tags$span(class = "text-muted",
                       length(rv$pd_protein_line_idx), " proteins"))
    }
  })

  # Clone for Trends sidebar
  output$trd_file_status <- renderUI({
    if (is.null(rv$pd_filename)) {
      tags$p(class = "text-muted small", icon("circle-xmark"), " No file loaded")
    } else {
      tags$p(class = "small text-success",
             icon("circle-check"), " ", rv$pd_filename,
             tags$br(),
             tags$span(class = "text-muted",
                       length(rv$pd_protein_line_idx), " proteins"))
    }
  })

  # ===========================================================================
  # UPLOAD WIZARD — Modal-based guided setup replacing sidebar accordion
  # Writes to existing rv$ fields on completion. (Phase 1: existing pipelines)
  # ===========================================================================

  # ── Reactive outputs for conditionalPanel conditions and sidebar ──
  output$wiz_completed_flag <- reactive({ isTRUE(rv$wizard_completed) })
  outputOptions(output, "wiz_completed_flag", suspendWhenHidden = FALSE)

  output$wiz_export_btn <- renderUI({
    if (!isTRUE(rv$wizard_completed)) return(NULL)
    if (rv$lf_pipeline %in% c("global", "ip")) {
      downloadButton("upload_export_summarization",
        tagList(icon("download"), " Export Summarization (.rds)"),
        class = "btn-outline-info btn-sm w-100 mb-2")
    } else {
      NULL
    }
  })

  # ── Build wizard modal shell ──
  .wiz_build_modal <- function() {
    modalDialog(
      title     = uiOutput("wiz_modal_title"),
      uiOutput("wiz_step_content"),
      footer    = uiOutput("wiz_footer"),
      size      = "l",
      easyClose = FALSE
    )
  }

  # ── New-dataset reset for the LF wizard ──
  # An accepted primary upload (PD CSV or MaxQuant proteinGroups.txt) replaces
  # the whole experiment (all data types, every module's results) before any of
  # it is written. evidence.txt is the second half of a MaxQuant upload and
  # never resets. A rejected file resets nothing.
  .lf_wiz_begin_new_data <- function() {
    rv$reset_experiment_state(keep = c(RESET_KEEP_ON_UPLOAD,
      "^wizard_step$", "^wizard_modal_open$",
      "^lf_shared_peptides$", "^mq_remove_contaminants$"))
    logger("LF Wizard: New file uploaded - previous experiment state reset")
  }

  observeEvent(input$wiz_open, {
    rv$wizard_step <- 3L
    rv$wizard_modal_open <- TRUE
    showModal(.wiz_build_modal())
  })

  observeEvent(input$wiz_reopen, {
    rv$wizard_step <- 3L
    rv$wizard_modal_open <- TRUE
    showModal(.wiz_build_modal())
    # Reset inputs so they re-bind properly in the new modal session
    updateTextInput(session, "wiz_new_group_name", value = "")
  })

  # ── Modal title with step progress indicator ──
  output$wiz_modal_title <- renderUI({
    step      <- rv$wizard_step
    reachable <- rv$wizard_reachable_steps %||% c(3L, 4L, 6L, 7L)
    n_steps   <- length(reachable)
    eff_step  <- which(reachable == step)
    eff_step  <- if (length(eff_step) == 0L) n_steps else eff_step
    lbls_all  <- c("Search Program", "Data Types", "File Upload",
                   "Condition Groups", "Factorial Design", "Experiment Design",
                   "Confirm & Summarize", "Batch Correction")
    lbl <- if (step >= 1L && step <= 8L) lbls_all[[step]] else ""
    tagList(
      tags$div(class = "d-flex justify-content-between align-items-center",
        tags$span(icon("wand-magic-sparkles"), " Upload Wizard"),
        tags$span(class = "text-muted small",
                  paste0("Step ", eff_step, " / ", n_steps, " \u2014 ", lbl))
      ),
      tags$div(class = "progress mt-2", style = "height: 5px;",
        tags$div(class = "progress-bar bg-primary", role = "progressbar",
                 style = paste0("width:", round(eff_step / n_steps * 100), "%;"))
      )
    )
  })

  # ── Step content router ──
  output$wiz_step_content <- renderUI({
    switch(as.character(rv$wizard_step),
      "3" = .wiz_step3_ui(),
      "4" = .wiz_step4_ui(),
      "5" = .wiz_step6_factorial_ui(),
      "6" = .wiz_step5_ui(),
      "7" = .wiz_step7_summarize_ui(),
      "8" = .wiz_step8_batch_ui(),
      tags$p("Unknown step.")
    )
  })

  # ── Footer: Cancel (left) + Back / Next or Done (right) ──
  output$wiz_footer <- renderUI({
    step          <- rv$wizard_step
    has_batch     <- isTRUE(rv$wizard_has_batch_step)
    # show Next for steps 3-6, and step 7 when batch step follows
    show_next     <- step < 7L || (step == 7L && has_batch)
    cancel_btn <- actionButton("wiz_cancel",
      tagList(icon("xmark"), " Cancel"),
      class = "btn-outline-secondary btn-sm")
    back_btn <- if (step > 3L)
      actionButton("wiz_back", tagList(icon("arrow-left"), " Back"),
                   class = "btn-outline-secondary ms-1")
    else NULL
    right_btn <- if (show_next) {
      # Step 7 → batch: block Next until summarization is done
      cache_ready   <- !is.null(rv$lf_master_cache)
      step7_blocked <- (step == 7L && has_batch && !cache_ready)
      if (step7_blocked)
        tags$button(tagList(" Next ", icon("arrow-right")),
          class = "btn btn-secondary ms-1", type = "button", disabled = NA)
      else
        actionButton("wiz_next",
          tagList(" Next ", icon("arrow-right")), class = "btn-primary ms-1")
    } else {
      # Step 7 (no batch) or step 8: Done turns green when cache is ready
      cache_ready  <- !is.null(rv$lf_master_cache)
      if (cache_ready)
        actionButton("wiz_close", tagList(icon("check"), " Done"),
                     class = "btn-success ms-1")
      else
        tags$button(tagList(icon("check"), " Done"),
          class = "btn btn-secondary ms-1", type = "button", disabled = NA)
    }
    tags$div(
      class = "d-flex justify-content-between align-items-center w-100",
      cancel_btn,
      tags$div(back_btn, right_btn)
    )
  })

  observeEvent(input$wiz_back, {
    req(rv$wizard_step > 3L)
    step      <- isolate(rv$wizard_step)
    reachable <- isolate(rv$wizard_reachable_steps) %||% c(3L, 4L, 6L, 7L)
    if (step == 8L) {
      rv$wizard_step <- 7L
    } else if (step == 7L) {
      rv$wizard_step <- 6L
    } else if (step == 6L) {
      # Go to factorial (step 5) if reachable, else to step 4
      rv$wizard_step <- if (5L %in% reachable) 5L else 4L
    } else {
      rv$wizard_step <- step - 1L
    }
  })

  observeEvent(input$wiz_cancel, {
    removeModal()
    rv$wizard_modal_open <- FALSE
  })

  # ── Step 3: File Upload (auto-detected) ──
  # One file box. Whether the file is a Proteome Discoverer protein export or
  # MaxQuant proteinGroups.txt is read from its header row
  # (detect_lf_upload_type, global.R); that sets the search program and the
  # quantification route the retired steps 1-2 used to ask for. MaxQuant
  # evidence.txt is requested once proteinGroups.txt is accepted. The step body
  # reads no reactive value, so an upload never re-renders (and clears) its own
  # file box; what was loaded is shown by wiz_loaded_ui.
  .wiz_step3_ui <- function() {
    # Blue question mark with a hover hint (bslib tooltip, as used in ui.R)
    .hint <- function(text)
      bslib::tooltip(
        tags$span(class = "ms-1", style = "cursor: help; color: #3498db;",
                  icon("circle-question")),
        text, placement = "right")
    tagList(
      tags$h5("Upload Search Results"),
      tags$div(class = "border rounded p-2 mb-3 small",
        tags$div(class = "fw-bold mb-1", "OSMSS supports the following files:"),
        tags$ul(class = "mb-0 ps-3",
          tags$li("Proteome Discoverer files (protein, peptides or both, saved as a .CSV)",
            .hint(paste0("Not sure? Look for columns that say \u201cAbundance:\u201d. ",
                         "Make sure you save the export as CSV!"))),
          tags$li("MaxQuant files (proteinGroups.txt)",
            .hint(paste0("Not sure? Look for files named evidence.txt and ",
                         "proteinGroups.txt in the MaxQuant output folder ",
                         "(combined/txt). Upload proteinGroups.txt first; ",
                         "evidence.txt is asked for next."))),
          tags$li("FragPipe files (combined_protein.tsv + fragpipe.workflow)",
            .hint(paste0("Not sure? Both files are in the FragPipe output folder. Upload ",
                         "combined_protein.tsv first; fragpipe.workflow is asked for next and is ",
                         "required: it records whether IonQuant normalized the intensities.")))
        ),
        tags$div(class = "text-muted mt-1",
          "The file type is detected automatically from its column headers.")
      ),
      fileInput("wiz_file_any", "Choose file:",
                accept = c(".csv", ".txt", ".tsv"), width = "100%"),
      uiOutput("wiz_detect_msg_ui"),
      uiOutput("wiz_loaded_ui")
    )
  }

  # Last rejected upload. Cleared by the next accepted file.
  output$wiz_detect_msg_ui <- renderUI({
    msg <- rv$lf_wiz_detect_msg
    if (is.null(msg)) return(NULL)
    tags$div(class = "alert alert-danger py-2 px-3 small mt-2",
      icon("circle-xmark"), tags$b(paste0(" ", msg$file, ": ")), msg$text)
  })

  # Shared-peptide choice. Shown only where a peptide-aggregating route will run
  # (PD peptides re-summarized here, Route B). Reads the committed setting
  # reactively, so returning to step 3 shows what Next saved.
  .wiz_shared_pep_radio <- function() {
    shared_prev <- rv$lf_shared_peptides
    if (is.null(shared_prev) || is.na(shared_prev)) shared_prev <- "unique_only"
    tagList(
      tags$h6(class = "mt-2", "Peptides matching more than one protein"),
      radioButtons("wiz_shared_peptides", NULL,
        choiceNames = list(
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("filter"), " Discard shared peptides"),
            tags$span(class = "badge bg-success ms-2", "Recommended"),
            tags$br(),
            tags$small(class = "text-muted",
              "Peptides whose protein accession list names more than one ",
              "protein are removed before protein-level summarization. This is ",
              "the MSstats default (useUniquePeptide = TRUE).")),
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("triangle-exclamation"), " Keep all peptides"),
            tags$br(),
            tags$small(class = "text-muted",
              "Shared peptides are kept and attributed to the ",
              tags$strong("first accession listed"), " in the export. That is ",
              "position in a string, not evidence: it is not razor assignment, ",
              "and it can flatten the fold change of the protein that truly ",
              "contributed the signal.")),
          tags$div(class = "p-2 border rounded mb-1",
            tags$strong(icon("scale-balanced"), " Assign by parsimony (razor)"),
            tags$br(),
            tags$small(class = "text-muted",
              "Each shared peptide is assigned to the single candidate protein ",
              "with the most unique-peptide evidence across the dataset. ",
              "Peptides whose candidates have no unique support at all are ",
              "removed. This assignment is derived here from the exported ",
              "peptide table; it is not taken from the search engine."))
        ),
        choiceValues = c("unique_only", "keep_all", "razor"),
        selected = shared_prev
      )
    )
  }

  # What was loaded, and the choices that depend on it. Re-rendered when a
  # primary file (PD CSV or proteinGroups.txt) is accepted and, for PD, when Next
  # commits the route (so Back shows the saved choice) - never by evidence.txt,
  # whose file box lives in the MaxQuant branch.
  output$wiz_loaded_ui <- renderUI({
    sp <- rv$wizard_search_program
    if (identical(sp, "maxquant") && !is.null(rv$mq_proteingroups_filename)) {
      n_exp <- length(isolate(rv$mq_pg_experiments))
      tagList(
        tags$div(class = "alert alert-success py-2 px-3 small mt-2",
          icon("circle-check"), tags$b(" MaxQuant proteinGroups.txt: "),
          paste0(rv$mq_proteingroups_filename, " (", n_exp, " experiment",
                 if (n_exp != 1L) "s", ").")),
        fileInput("wiz_file_mq_evidence",
                  "Now upload evidence.txt from the same MaxQuant search:",
                  accept = c(".txt", ".tsv"), width = "100%"),
        helpText(class = "small text-muted mt-0",
          "evidence.txt supplies the run list used to build condition groups and shows ",
          "whether match between runs was on. Protein intensities are read from ",
          "proteinGroups.txt. Each run must have a matching quantity column there."),
        uiOutput("wiz_mq_upload_status"),
        tags$hr(),
        checkboxInput("wiz_mq_remove_contaminants", "Remove potential contaminant groups",
                      value = isTRUE(isolate(rv$mq_remove_contaminants))),
        helpText(class = "small text-muted mt-0",
          "MaxQuant flags a group as a potential contaminant when it contains any ",
          "contaminant-database entry (e.g. keratins, albumin, trypsin), even if another ",
          "protein leads the group. Off (default): these groups are kept and their leading ",
          "IDs carry MaxQuant's CON__ prefix. On: they are removed before analysis. ",
          "Reverse (decoy) and only-identified-by-site groups are always removed."),
        uiOutput("wiz_mq_protein_id_ui")
      )
    } else if (identical(sp, "pd") && !is.null(rv$pd_filename)) {
      lv     <- rv$pd_levels_available %||% character(0)
      n_prot <- length(isolate(rv$pd_protein_line_idx))
      tagList(
        tags$div(class = "alert alert-success py-2 px-3 small mt-2",
          icon("circle-check"), tags$b(" Proteome Discoverer export: "),
          paste0(rv$pd_filename, " (", format(n_prot, big.mark = ","),
                 " proteins; levels found: ", paste(lv, collapse = " + "), ").")),
        if (all(c("protein", "peptide") %in% lv)) {
          tagList(
            checkboxInput("wiz_pd_from_peptides",
              "Re-summarize proteins from the nested peptides",
              value = identical(rv$wizard_data_types, "peptide")),
            helpText(class = "small text-muted mt-0",
              "Off (default): Proteome Discoverer's protein abundances are used directly. ",
              "On: OSMSS summarizes each protein from its peptides (QFeatures ",
              "robustSummary), with the shared-peptide handling chosen below."),
            conditionalPanel(condition = "input.wiz_pd_from_peptides",
                             .wiz_shared_pep_radio())
          )
        } else if (identical(lv, "peptide")) {
          tagList(
            tags$div(class = "alert alert-info py-2 px-3 small mt-2",
              icon("circle-info"),
              " No protein rows were found, so proteins are summarized from the ",
              "nested peptides."),
            .wiz_shared_pep_radio()
          )
        } else {
          tags$div(class = "alert alert-info py-2 px-3 small mt-2",
            icon("circle-info"),
            " No nested peptide rows were found: Proteome Discoverer's protein ",
            "abundances are used directly.")
        }
      )
    } else if (identical(sp, "fragpipe") && !is.null(rv$pd_filename)) {
      # FragPipe: combined_protein.tsv was converted to PD layout on upload; the
      # required fragpipe.workflow is asked for here (its status in its own output,
      # so this panel and its file box are not re-rendered by that upload).
      hdr <- isolate(rv$pd_protein_headers)
      tagList(
        tags$div(class = "alert alert-success py-2 px-3 small mt-2",
          icon("circle-check"), tags$b(" FragPipe combined_protein.tsv: "),
          paste0(rv$pd_filename, " (", format(length(isolate(rv$pd_protein_line_idx)), big.mark = ","),
                 " proteins, ", sum(startsWith(hdr, "Abundance: ")), " samples; per-sample ",
                 "\u2018Intensity\u2019 columns are used, not \u2018MaxLFQ Intensity\u2019).")),
        fileInput("wiz_file_fp_workflow",
                  "Now upload fragpipe.workflow from the same FragPipe run (required):",
                  accept = c(".workflow", ".txt"), width = "100%"),
        helpText(class = "small text-muted mt-0",
          "fragpipe.workflow is saved in the FragPipe output folder, next to combined_protein.tsv. It ",
          "records whether IonQuant normalized the intensities across runs, which decides whether OSMSS ",
          "normalizes them again, and the quantification settings reported in the methods text."),
        uiOutput("wiz_fp_workflow_status")
      )
    } else {
      tags$div(class = "alert alert-secondary py-2 px-3 small mt-2",
        icon("circle-info"), " Awaiting file upload.")
    }
  })

  output$wiz_fp_workflow_status <- renderUI({
    wf <- rv$lf_fp_workflow
    if (is.null(wf))
      return(tags$div(class = "alert alert-secondary py-2 px-3 small mt-2",
        icon("circle-info"), " Waiting for fragpipe.workflow (required)."))
    ions <- if (isTRUE(is.finite(wf$top_n)) && wf$top_n > 0L) paste0("top ", wf$top_n, " ions") else "all ions"
    tags$div(class = "alert alert-success py-2 px-3 small mt-2",
      icon("circle-check"), tags$b(paste0(" ", rv$lf_fp_workflow_filename,
                                          if (!is.na(wf$version)) paste0(" (FragPipe ", wf$version, ")") else "", ": ")),
      if (isTRUE(wf$normalization))
        paste0("IonQuant normalization was ON, so the intensities are already normalized across runs and ",
               "OSMSS will not normalize them again (in OSMSS\u2019s DDA benchmark, normalizing FragPipe ",
               "intensities a second time increased false positives).")
      else paste0("IonQuant normalization was OFF, so OSMSS will apply its global median normalization."),
      tags$br(),
      tags$span(class = "text-muted", paste0(
        "Intensities from ", ions,
        if (!is.na(wf$uniqueness)) paste0(", ", wf$uniqueness, " peptides") else "",
        "; match between runs ", if (isTRUE(wf$mbr)) "on" else "off", ".")))
  })

  # Required fragpipe.workflow: its ionquant.normalization decides whether OSMSS normalizes
  observeEvent(input$wiz_file_fp_workflow, {
    req(input$wiz_file_fp_workflow, identical(rv$wizard_search_program, "fragpipe"))
    f  <- input$wiz_file_fp_workflow
    wf <- fp_read_workflow(f$datapath)
    if (!is.null(wf$error)) {
      rv$lf_wiz_detect_msg <- list(file = f$name, text = wf$error)
      showNotification(wf$error, type = "error", duration = 10)
      return()
    }
    rv$lf_fp_workflow          <- wf
    rv$lf_fp_workflow_filename <- f$name
    rv$lf_apply_normalization  <- !isTRUE(wf$normalization)
    rv$lf_wiz_detect_msg       <- NULL
    logger(paste0("LF Wizard: fragpipe.workflow ", f$name, " (FragPipe ", wf$version %||% "?",
                  "): IonQuant normalization ", if (isTRUE(wf$normalization)) "on - OSMSS normalization off"
                  else "off - OSMSS normalization on"))
  })

  # ── Primary upload: detect, then load as PD or MaxQuant proteinGroups ──
  observeEvent(input$wiz_file_any, {
    req(input$wiz_file_any)
    f   <- input$wiz_file_any
    det <- detect_lf_upload_type(f$datapath)
    logger(paste0("LF Wizard: '", f$name, "' detected as ", det$kind))
    .reject <- function(text) {
      rv$lf_wiz_detect_msg <- list(file = f$name, text = text)
      showNotification(text, type = "error", duration = 10)
    }

    if (identical(det$kind, "pd_nested")) {
      .lf_wiz_begin_new_data()
      rv$wizard_search_program <- "pd"
      rv$lf_input_mode         <- "peptide"
      rv$lf_input_format       <- "pd_nested"
      ok <- tryCatch({ parse_pd_csv(f); TRUE }, error = function(e) {
        logger(paste0("!!! LF Wizard: PD CSV parse failed: ", e$message))
        FALSE
      })
      lines    <- rv$pd_raw_lines
      has_prot <- length(rv$pd_protein_line_idx) > 0L
      # Same nested-peptide test as Route B (build_proda_protein_matrix), which
      # needs at least 10 peptide rows
      pep_hdr  <- grep("^,.*Annotated Sequence", lines)
      has_pep  <- length(pep_hdr) > 0L &&
                  length(setdiff(grep("^,", lines), pep_hdr)) >= 10L
      if (!ok || (!has_prot && !has_pep)) {
        .lf_wiz_begin_new_data()   # drop the partly loaded file
        .reject(if (!ok) "The file could not be parsed as a Proteome Discoverer export."
                else paste0("No protein rows (FDR confidence High/Medium/Low) and no nested ",
                            "peptide rows were found. Export the Proteins table with ",
                            "peptides nested under each protein."))
        return()
      }
      rv$pd_levels_available <- c(if (has_prot) "protein", if (has_pep) "peptide")
      # Default route: PD protein abundances (Route A) when present. The step 3
      # switch can move it to peptide re-summarization (Route B) on Next.
      rv$wizard_data_types   <- if (has_prot) "protein" else "peptide"
      rv$lf_wiz_detect_msg   <- NULL

    } else if (identical(det$kind, "mq_proteingroups")) {
      hdrs <- make.names(det$headers, unique = TRUE)
      exps <- mq_pg_experiments(det$headers)
      if (length(exps) == 0L) {
        .reject(paste0("proteinGroups.txt has no 'Intensity <experiment>' or 'LFQ intensity ",
                       "<experiment>' columns, so it carries no protein quantities."))
        return()
      }
      # Namespace the stable copy per session: a bare filename collides between
      # concurrent sessions sharing tempdir() on a multi-user server.
      stable <- file.path(tempdir(), paste0("wiz_mq_pg_", substr(session$token, 1, 8), "_",
                                            basename(f$name)))
      file.copy(f$datapath, stable, overwrite = TRUE)
      .lf_wiz_begin_new_data()
      rv$wizard_search_program     <- "maxquant"
      rv$wizard_data_types         <- "maxquant"
      rv$lf_input_mode             <- "maxquant"
      rv$lf_input_format           <- "maxquant"
      rv$mq_pg_headers             <- hdrs
      rv$mq_pg_experiments         <- exps
      rv$mq_pg_has_intensity       <- any(grepl("^Intensity\\.", hdrs))
      rv$mq_proteingroups_path     <- stable
      rv$mq_proteingroups_filename <- f$name
      rv$lf_wiz_detect_msg         <- NULL

    } else if (identical(det$kind, "fp_combined_protein")) {
      cv <- fp_combined_protein_to_pd(f$datapath)
      if (!is.null(cv$error)) { .reject(cv$error); return() }
      # FragPipe takes the PD protein-direct route (Route A): the table is held in
      # PD layout ("Abundance: <sample>" columns, Accession = Protein ID).
      .lf_wiz_begin_new_data()
      rv$wizard_search_program <- "fragpipe"
      rv$wizard_data_types     <- "protein"
      rv$lf_input_mode         <- "peptide"
      rv$lf_input_format       <- "pd_nested"
      rv$pd_levels_available   <- "protein"
      rv$pd_raw_lines          <- cv$lines
      rv$pd_protein_line_idx   <- cv$protein_idx
      rv$pd_protein_headers    <- cv$headers
      rv$pd_other_cols         <- cv$headers
      rv$pd_ratio_cols         <- character(0)
      rv$pd_pval_cols          <- character(0)
      rv$pd_has_gn             <- FALSE
      rv$pd_gene_col           <- if (isTRUE(cv$has_gene)) "Gene Symbol" else NULL
      rv$pd_filename           <- f$name
      rv$provenance$global     <- list(input_filename = f$name, total_initial_proteins = cv$n_proteins,
                                       input_software_source = "FragPipe")
      rv$lf_wiz_detect_msg     <- NULL
      logger(paste0("LF Wizard: FragPipe combined_protein.tsv: ", cv$n_proteins, " proteins, ",
                    length(cv$samples), " samples", if (cv$n_decoys_removed > 0L)
                      paste0(", ", cv$n_decoys_removed, " decoy rows removed") else ""))

    } else {
      .reject(det$message)
    }
  })

  # ── MaxQuant evidence.txt (asked for after proteinGroups.txt) ──
  observeEvent(input$wiz_file_mq_evidence, {
    req(input$wiz_file_mq_evidence, rv$mq_proteingroups_path)
    f <- input$wiz_file_mq_evidence
    .reject <- function(text) {
      rv$lf_wiz_detect_msg <- list(file = f$name, text = text)
      showNotification(text, type = "error", duration = 10)
    }
    tryCatch({
      stable <- file.path(tempdir(), paste0("wiz_mq_evidence_", substr(session$token, 1, 8), "_",
                                            basename(f$name)))
      file.copy(f$datapath, stable, overwrite = TRUE)
      # data.table::fread preserves original column names (spaces, etc.); normalise to R dot-notation
      hdrs_orig <- colnames(data.table::fread(stable, nrows = 0, data.table = FALSE))
      hdrs      <- make.names(hdrs_orig, unique = TRUE)
      required  <- c("Raw.file", "Sequence", "Proteins", "Intensity")
      missing   <- required[!required %in% hdrs]
      if (length(missing) > 0) {
        .reject(paste0("This is not a MaxQuant evidence.txt: missing column(s) ",
                       paste(missing, collapse = ", "), "."))
        return()
      }
      raw_col_idx <- which(hdrs == "Raw.file")
      raw_dt      <- data.table::fread(stable, select = raw_col_idx, data.table = FALSE)
      raw_files   <- unique(as.character(raw_dt[[1]]))
      # Raw.file -> Experiment mapping (used by Route MQ for quantity column names)
      raw_to_exp <- NULL
      if ("Experiment" %in% hdrs) {
        exp_col_idx <- which(hdrs == "Experiment")
        exp_dt <- data.table::fread(stable, select = c(raw_col_idx, exp_col_idx),
                                    data.table = FALSE)
        raw_exp_pairs <- unique(exp_dt)
        raw_to_exp <- setNames(as.character(raw_exp_pairs[[2]]),
                               as.character(raw_exp_pairs[[1]]))
      }
      # Same-search check: every run must resolve to a proteinGroups.txt quantity
      # column by the rule the matrix build uses (mq_resolve_quant_col).
      q_cols    <- vapply(raw_files, mq_resolve_quant_col, character(1L),
                          raw_to_exp = raw_to_exp, cn_norm = rv$mq_pg_headers)
      unmatched <- raw_files[is.na(q_cols)]
      if (length(unmatched) > 0) {
        .reject(paste0(length(unmatched), " of ", length(raw_files), " runs in this evidence.txt ",
                       "have no matching quantity column in proteinGroups.txt (e.g. ",
                       paste(head(unmatched, 3), collapse = ", "), "). Both files must come ",
                       "from the same MaxQuant search."))
        return()
      }
      rv$mq_evidence_path         <- stable
      rv$mq_evidence_filename     <- f$name
      rv$mq_raw_files             <- raw_files
      rv$mq_rawfile_to_experiment <- raw_to_exp
      if (!is.null(raw_to_exp))
        logger(paste0("[MQ] Raw.file->Experiment map: ", length(raw_to_exp), " entries"))
      # Match between runs: features transferred by matching have a Type containing
      # "MATCH" (MULTI-MATCH, MULTI-MATCH-MSMS, TIMS-MULTI-MATCH) and a filled
      # "Match score" / "Match time difference". Either marker counts, so a new
      # MaxQuant type label does not hide it. Off is MaxQuant's default.
      .mbr_cols <- intersect(c("Type", "Match.score", "Match.time.difference"), hdrs)
      if (length(.mbr_cols) > 0) {
        .mb <- data.table::fread(stable, select = match(.mbr_cols, hdrs), data.table = FALSE)
        names(.mb) <- hdrs[match(.mbr_cols, hdrs)]
        .is_match <- rep(FALSE, nrow(.mb))
        if ("Type" %in% names(.mb)) .is_match <- .is_match | grepl("MATCH", .mb$Type, fixed = TRUE)
        for (.mc in intersect(c("Match.score", "Match.time.difference"), names(.mb)))
          .is_match <- .is_match | (!is.na(.mb[[.mc]]) & nzchar(as.character(.mb[[.mc]])))
        rv$mq_mbr_on              <- any(.is_match)
        rv$mq_n_match_features    <- sum(.is_match)
        rv$mq_n_evidence_features <- nrow(.mb)
        logger(paste0("[MQ] Match between runs: ", if (any(.is_match)) "on" else "off", " (",
                      sum(.is_match), " of ", nrow(.mb), " evidence features transferred by matching)"))
      } else {
        rv$mq_mbr_on <- NULL; rv$mq_n_match_features <- NULL; rv$mq_n_evidence_features <- NULL
      }
      rv$lf_wiz_detect_msg <- NULL
    }, error = function(e) {
      showNotification(paste0("evidence.txt upload error: ", e$message),
                       type = "error", duration = 8)
    })
  })

  output$wiz_mq_upload_status <- renderUI({
    ev_ok <- !is.null(rv$mq_evidence_filename)
    pg_ok <- !is.null(rv$mq_proteingroups_filename)
    if (ev_ok && pg_ok) {
      tagList(
        tags$div(class = "alert alert-success py-2 px-3 small mt-2",
          icon("circle-check"),
          paste0(" ", rv$mq_evidence_filename,
                 " (", length(rv$mq_raw_files), " runs) + ",
                 rv$mq_proteingroups_filename, " loaded.")),
        # Match between runs: information only; the analysis is the same either way
        if (!is.null(rv$mq_mbr_on))
          tags$div(class = "alert alert-secondary py-2 px-3 small mt-1",
            icon("circle-info"),
            if (isTRUE(rv$mq_mbr_on))
              paste0(" Match between runs was on in MaxQuant (",
                     format(rv$mq_n_match_features, big.mark = ","), " of ",
                     format(rv$mq_n_evidence_features, big.mark = ","),
                     " features transferred by matching).")
            else
              paste0(" Match between runs was off in MaxQuant (MaxQuant's default): proteins are ",
                     "quantified only in samples where one of their peptides was identified by ",
                     "MS/MS, so more values are missing. The analysis is unchanged."))
      )
    } else {
      tags$div(class = "alert alert-secondary py-2 px-3 small mt-2",
        icon("circle-info"), " Waiting for evidence.txt (required).")
    }
  })

  output$wiz_mq_protein_id_ui <- renderUI({
    req(rv$mq_proteingroups_filename)
    # proteinGroups.txt intensities are protein-level values MaxQuant already
    # computed, including its razor/unique peptide assignment, so OSMSS offers
    # no protein-ID or razor choice on this route.
    quant_box <- tags$div(class = "alert alert-info py-2 px-3 small mt-2",
      icon("circle-info"),
      tags$b(" Protein quantification: "),
      "raw protein-group \u2018Intensity\u2019 values are read directly from proteinGroups.txt, one row per ",
      "protein group, labelled by its leading protein, and normalized once by OSMSS. ",
      "MaxQuant\u2019s \u2018LFQ intensity\u2019 values are not used: they are already normalized by ",
      "MaxQuant, and normalizing them again would apply a second normalization. Shared-peptide (razor) assignment ",
      "was already made by MaxQuant, following its \u2018Peptides used for protein ",
      "quantification\u2019 setting (see parameters.txt; default: razor + unique). ",
      "Reverse hits and only-identified-by-site groups are removed; potential contaminant ",
      "groups are ", if (isTRUE(rv$mq_remove_contaminants)) "removed" else "kept",
      " (option above).")
    lfq_box <- if (!isTRUE(rv$mq_pg_has_intensity))
      tags$div(class = "alert alert-warning py-2 px-3 small mt-2",
        icon("triangle-exclamation"),
        " No raw \u2018Intensity\u2019 columns found \u2014 \u2018LFQ intensity\u2019 columns will ",
        "be used. These are already normalized by MaxQuant, and OSMSS will normalize them again; ",
        "re-export proteinGroups.txt with Intensity columns if possible.")
    tagList(quant_box, lfq_box)
  })

  # MaxQuant contaminant option -> rv (read by build_proda_protein_matrix, MQ route)
  observeEvent(input$wiz_mq_remove_contaminants, {
    rv$mq_remove_contaminants <- isTRUE(input$wiz_mq_remove_contaminants)
  })

  # ── Step 4: Condition Groups ──
  .wiz_step4_ui <- function() {
    use_mq  <- identical(rv$wizard_search_program, "maxquant")
    grps    <- rv$lf_cond_groups
    n_grps  <- length(grps)

    # Compute available items
    if (use_mq) {
      avail <- rv$mq_raw_files %||% character(0)
    } else {
      all_c <- rv$pd_other_cols %||% character(0)
      # Raw PD abundances are "Abundance: <anything>". Derived PD columns put a
      # qualifier before the first colon -- "Abundances (Grouped):",
      # "Abundances (Scaled):", "Abundances (Normalized):", "Abundances Count:",
      # "Abundance Ratio:" -- so classify on that prefix only (keywords in the
      # sample label, e.g. "Normal", must not drop a raw column).
      abund_prefix <- trimws(sub(":.*$", "", all_c))
      abund <- all_c[grepl(":", all_c) &
                     grepl("^Abundance$", abund_prefix, ignore.case = TRUE)]
      if (length(abund) == 0) {
        # Fallback for exports without the "Abundance:" form
        abund <- all_c[grepl("^Abundance", all_c, ignore.case = TRUE)]
        abund <- abund[!grepl("\\(|ratio|p-val|variab|count|norm|grouped|scaled",
                              abund, ignore.case = TRUE)]
      }
      avail <- if (length(abund) > 0) abund else all_c
    }
    assigned  <- unlist(grps, use.names = FALSE)
    avail_now <- setdiff(avail, assigned)
    # For MQ: prefer Experiment labels over Raw.file names in the UI
    raw_to_exp    <- if (use_mq) rv$mq_rawfile_to_experiment else NULL
    .item_labels  <- function(x) if (!is.null(raw_to_exp) && length(raw_to_exp) > 0)
      ifelse(x %in% names(raw_to_exp), raw_to_exp[x], x) else x
    avail_choices <- if (length(avail_now) > 0) setNames(avail_now, .item_labels(avail_now)) else avail_now

    tagList(
      tags$h5("Define Condition Groups"),
      tags$p(class = "small text-muted",
        "Assign ", if (use_mq) "runs" else "sample columns",
        " to named condition groups. At least 2 groups required for DE."),

      # ── Existing groups ──
      if (n_grps > 0)
        tags$div(
          tags$label(class = "form-label fw-bold small mb-1",
                     paste0("Defined groups (", n_grps, "):")),
          tags$div(style = "max-height: 200px; overflow-y: auto; margin-bottom: 8px;",
            lapply(seq_along(grps), function(i) {
              gname <- names(grps)[[i]]
              items <- grps[[i]]
              tags$div(class = "card mb-1",
                tags$div(class = "card-body py-1 px-2",
                  tags$div(class = "d-flex align-items-center gap-2",
                    tags$strong(gname, class = "me-auto"),
                    tags$small(class = "text-muted badge bg-secondary",
                               paste0(length(items), " item", if (length(items) != 1) "s")),
                    actionButton(paste0("wiz_del_grp_", i), "",
                                 icon = icon("trash"),
                                 class = "btn-sm btn-outline-danger py-0 px-1")
                  ),
                  tags$div(class = "text-muted small mt-1",
                    style = "font-size: 0.78em; word-break: break-all;",
                    paste(.item_labels(items), collapse = ", "))
                )
              )
            })
          ),
          actionButton("wiz_clear_all_groups", tagList(icon("trash"), " Clear All"),
                       class = "btn-outline-danger btn-sm mb-2")
        )
      else
        tags$div(class = "alert alert-warning py-2 px-3 small mb-2",
          icon("triangle-exclamation"), " No groups defined yet."),

      tags$hr(class = "my-2"),

      # ── Add new group ──
      tags$label(class = "form-label fw-bold small",
                 if (use_mq) "Add new group from runs:" else "Add new group from columns:"),
      textInput("wiz_new_group_name", NULL,
                placeholder = "Group name (e.g. Control)", value = ""),
      if (length(avail_now) > 0) {
        tags$div(
          tags$div(class = "alert alert-info py-1 px-2 small mb-1",
            icon("circle-info"),
            paste0(" ", length(avail_now), " unassigned item",
                   if (length(avail_now) != 1) "s", " \u2014 check boxes to select:")),
          checkbox_filter_controls("wiz_new_group_search", "wiz_new_group_cols",
                                   placeholder = "Search (contains)..."),
          tags$div(
            style = paste0(
              "max-height: 220px; overflow-y: auto; border: 1px solid #dee2e6;",
              " border-radius: 4px; padding: 6px 10px; background: #fff;"
            ),
            checkboxGroupInput(
              "wiz_new_group_cols", NULL,
              choices  = avail_choices,
              selected = character(0),
              width    = "100%"
            )
          )
        )
      } else if (n_grps >= 2) {
        tags$div(class = "alert alert-success py-2 px-3 small mb-1",
          icon("circle-check"),
          paste0(" All items assigned. ", n_grps, " groups defined. Click Next to continue."))
      } else {
        tags$div(class = "alert alert-warning py-2 px-3 small mb-1",
          icon("triangle-exclamation"),
          " All available items are already assigned to groups.")
      },
      actionButton("wiz_add_group_btn", tagList(icon("plus"), " Add Group"),
                   class = "btn-success btn-sm w-100 mt-2"),
      # ── Pairing section ──
      tags$hr(class = "my-3"),
      checkboxInput("wiz_enable_pairing",
        tagList(
          icon("link", class = "me-1"),
          tags$strong("Samples are paired/matched across conditions"),
          tags$span(class = "text-muted small",
            " (e.g., same patient pre and post treatment)")
        ),
        value = !is.null(rv$lf_subject_pairing)),
      conditionalPanel(
        condition = "input.wiz_enable_pairing == true",
        tags$div(class = "alert alert-info py-2 px-3 small mt-2",
          tags$strong(icon("circle-info", class = "me-1"), "Paired Design"),
          tags$br(),
          "Assign a Subject ID to each sample below. ",
          "Samples sharing the same Subject ID will be treated as ",
          "repeated measures from the same biological unit.",
          tags$br(), tags$br(),
          "For ProDA: Subject IDs are used to account for repeated measures in the ",
          "differential expression model.",
          tags$br(),
          "For ProDA: Subject ID is added as a fixed covariate in the design formula."
        ),
        uiOutput("wiz_pairing_table")
      )
    )
  }

  # ── Pairing table: one row per sample, Subject ID selectize ──
  output$wiz_pairing_table <- renderUI({
    grps <- rv$lf_cond_groups
    if (length(grps) == 0) return(tags$p(class = "text-muted small", "Define condition groups first."))
    all_samples <- unlist(lapply(names(grps), function(g) {
      lapply(grps[[g]], function(s) list(sample = s, condition = g))
    }), recursive = FALSE)
    n_samps     <- length(all_samples)
    n_conds     <- length(grps)
    n_subjects  <- ceiling(n_samps / n_conds)
    subject_suggestions <- paste0("Subject_", seq_len(max(n_subjects, 1L)))
    saved_pairing <- rv$lf_subject_pairing  # named char vec or NULL
    rows <- lapply(all_samples, function(x) {
      sname   <- x$sample
      input_id <- paste0("wiz_pair_", make.names(sname))
      saved_val <- if (!is.null(saved_pairing) && sname %in% names(saved_pairing))
        saved_pairing[[sname]] else ""
      tags$tr(
        tags$td(class = "small", sname),
        tags$td(class = "small text-muted", x$condition),
        tags$td(selectizeInput(input_id, NULL,
          choices  = c(setNames("", "Select or type..."), subject_suggestions),
          selected = saved_val,
          options  = list(create = TRUE, placeholder = "Select or type..."),
          width    = "100%"))
      )
    })
    tags$table(class = "table table-sm table-bordered mt-2",
      tags$thead(tags$tr(
        tags$th(class = "small", "Sample"),
        tags$th(class = "small", "Condition"),
        tags$th(class = "small", "Subject ID")
      )),
      tags$tbody(rows)
    )
  })

  observeEvent(input$wiz_add_group_btn, {
    gname <- trimws(input$wiz_new_group_name %||% "")
    gcols <- input$wiz_new_group_cols
    if (!nzchar(gname)) {
      showNotification("Enter a group name.", type = "warning"); return()
    }
    if (length(gcols) == 0) {
      showNotification("Check at least one item to include in the group.",
                       type = "warning"); return()
    }
    cur <- rv$lf_cond_groups
    cur[[gname]] <- gcols
    rv$lf_cond_groups <- cur
    updateTextInput(session, "wiz_new_group_search", value = "")
    showNotification(
      paste0("Group \u2018", gname, "\u2019 saved (", length(gcols), " item",
             if (length(gcols) != 1) "s", ")."),
      type = "message", duration = 3
    )
  })

  observeEvent(input$wiz_clear_all_groups, {
    rv$lf_cond_groups <- list()
  })

  # Delete individual group buttons (up to 20 groups)
  lapply(1:20, function(i) {
    observeEvent(input[[paste0("wiz_del_grp_", i)]], {
      cur <- rv$lf_cond_groups
      if (i <= length(cur)) {
        deleted_name <- names(cur)[[i]]
        cur[[i]] <- NULL
        rv$lf_cond_groups <- cur
        showNotification(paste0("Group \u2018", deleted_name, "\u2019 removed."),
                         type = "message", duration = 2)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
  })

  # ── Step 6 content: Experimental Design (card-button Q&A) ──
  # Two-question flow:
  #   Q1  wiz_q1_ip      — "yes" (did a pulldown) | "no" (lysate comparison)
  #   Q2a wiz_q2a_ctrl   — "trap" | "depleted" | "cond"  (IP path)
  #   Q2b wiz_q2b_abs    — "absent" | "quant"             (non-IP path)
  .wiz_step5_ui <- function() {
    q1_sel  <- rv$lf_wizard_q1_ip    %||% character(0)
    q2a_sel <- rv$lf_wizard_q2a_ctrl %||% character(0)
    q2b_sel <- rv$lf_wizard_q2b_abs  %||% character(0)
    fd      <- rv$lf_factorial_design
    is_fact <- !is.null(fd) && isTRUE(fd$enabled)

    # Reactive "no-IP / quant" terminal label (changes with factorial state)
    no_engine_desc <- if (is_fact) {
      tagList(tags$strong("ProDA + Global Median (factorial)"),
              " \u2014 additive/interaction formula for multi-factor designs.")
    } else {
      tagList(tags$strong("ProDA + Global Median"),
              " \u2014 standard differential expression for whole-proteome comparisons.")
    }

    # ── Card helper ──────────────────────────────────────────────────────────
    .card <- function(input_id, value, ico, title, subtitle, sel) {
      is_sel <- length(sel) > 0L && identical(sel, value)
      tags$div(
        class            = paste0("wiz-exp-card", if (isTRUE(is_sel)) " wiz-exp-card--selected"),
        `data-cardgroup` = input_id,
        onclick          = sprintf(
          "document.querySelectorAll('.wiz-exp-card[data-cardgroup=\"%s\"]').forEach(function(e){e.classList.remove('wiz-exp-card--selected');}); this.classList.add('wiz-exp-card--selected'); Shiny.setInputValue('%s','%s',{priority:'event'});",
          input_id, input_id, value),
        icon(ico, class = "wiz-card-icon fa-fw"),
        tags$div(
          tags$div(class = "fw-semibold lh-sm", title),
          tags$div(class = "small text-muted mt-1", subtitle)
        )
      )
    }

    # ── Terminal result helper (expandable <details>) ─────────────────────────
    .terminal <- function(chip_label, algo, norm, best_for) {
      tags$div(
        class = "wiz-terminal-wrap mt-2",
        tags$details(
          class = "wiz-result-details",
          tags$summary(
            tags$span(class = "wiz-result-chip",
              icon("circle-check", class = "me-1"),
              chip_label,
              tags$small(class = "ms-2 wiz-chip-hint", "\u25b8 details")
            )
          ),
          tags$div(
            class = "wiz-result-body small mt-2",
            tags$div(tags$span(class = "fw-semibold", "Algorithm: "), algo),
            tags$div(class = "mt-1", tags$span(class = "fw-semibold", "Normalization: "), norm),
            tags$div(class = "mt-1", tags$span(class = "fw-semibold", "Best for: "), best_for)
          )
        )
      )
    }

    tagList(
      tags$style(HTML("
        .wiz-card-row { display:flex; gap:.75rem; margin-bottom:.25rem; }
        .wiz-card-col { display:flex; flex-direction:column; gap:.5rem; margin-bottom:.25rem; }
        .wiz-exp-card {
          flex:1; display:flex; align-items:flex-start; gap:.75rem; cursor:pointer;
          border:2px solid #dee2e6; border-radius:.5rem; padding:.75rem; background:#fff;
          transition:border-color .15s, background .15s; user-select:none;
        }
        .wiz-exp-card:hover { border-color:#86b7fe; background:#f8f9ff; }
        .wiz-exp-card--selected { border-color:#0d6efd !important; background:#e7f0ff !important; }
        .wiz-card-icon { font-size:1.2rem; padding-top:.1rem; color:#0d6efd; flex-shrink:0; }
        .wiz-q-label { font-size:.875rem; font-weight:600; margin-bottom:.5rem; margin-top:.75rem; }
        .wiz-result-chip {
          display:inline-flex; align-items:center; gap:.3rem;
          background:#d1e7dd; color:#0a3622; border-radius:1rem;
          padding:.25rem .75rem; font-size:.82rem; font-weight:600;
          cursor:pointer;
        }
        .wiz-chip-hint { opacity:.65; }
        .wiz-result-details summary { list-style:none; }
        .wiz-result-details summary::-webkit-details-marker { display:none; }
        .wiz-result-body {
          border-left:3px solid #a3cfbb; padding-left:.75rem; margin-top:.4rem; color:#444;
        }
      ")),

      tags$h5("Experimental Design"),
      tags$p(class = "small text-muted mb-3",
        "Answer two quick questions so OSMSS can select the right statistical pipeline."),

      # ── Q1 ───────────────────────────────────────────────────────────────────
      tags$div(class = "wiz-q-label mt-0",
        "Did both your samples go through a pulldown or immunoprecipitation (IP) step?"),
      tags$div(class = "wiz-card-row",
        .card("wiz_q1_ip", "yes", "filter",
              "Yes \u2014 I did a pulldown / IP",
              "FLAG/HA/GFP pull-down, co-IP, BioID, antibody capture",
              q1_sel),
        .card("wiz_q1_ip", "no", "vials",
              "No \u2014 comparing cell or tissue lysates",
              "Whole-proteome, total cell extract, CRISPR KO, drug treatment",
              q1_sel)
      ),

      # ── Q1 = yes → Q2a ───────────────────────────────────────────────────────
      conditionalPanel(
        condition = "input.wiz_q1_ip == 'yes'",
        uiOutput("wiz_apms_group_filter_ui"),
        tags$div(class = "wiz-q-label",
          "What is in your control sample?"),
        tags$div(class = "wiz-card-col",
          .card("wiz_q2a_ctrl", "trap", "circle-xmark",
                "Non-specific capture reagent",
                "IgG, GFP-alone, FLAG-alone, His-tag, empty vector \u2014 same IP protocol, no bait",
                q2a_sel),
          .card("wiz_q2a_ctrl", "depleted", "layer-group",
                "Different cellular fraction or compartment",
                "e.g. cytoplasmic lysate as control for a nuclear pulldown",
                q2a_sel),
          .card("wiz_q2a_ctrl", "cond", "code-compare",
                "Same bait pulled down from a different biological condition",
                "e.g. IP of protein X in WT vs. KO, or treated vs. untreated",
                q2a_sel)
        ),
        # Q2a = trap → SAINTexpress
        conditionalPanel(
          condition = "input.wiz_q2a_ctrl == 'trap'",
          .terminal(
            "SAINTexpress (intensity mode)",
            "SAINTexpress v3.6.3 \u2014 Bayesian interaction scoring",
            "None (raw intensities passed directly)",
            "FLAG/HA/GFP pull-downs, co-IP, proximity labeling (BioID/APEX)"
          ),
          uiOutput("wiz_bait_gene_ui")
        ),
        # Q2a = depleted → ProDA + invariant-protein median normalization (fractionated)
        conditionalPanel(
          condition = "input.wiz_q2a_ctrl == 'depleted'",
          .terminal(
            "ProDA + invariant-protein median normalization (fractionated AP-MS)",
            "ProDA probabilistic dropout model",
            paste0("Invariant-protein median normalization (median over a data-driven ",
                   "set of unchanged proteins; robust to large one-directional shifts)"),
            "Nuclear, chromatin, or organelle-localized pulldowns where the control fraction is depleted"
          )
        ),
        # Q2a = cond → ProDA + invariant-protein median normalization (quantitative AP-MS)
        conditionalPanel(
          condition = "input.wiz_q2a_ctrl == 'cond'",
          .terminal(
            "ProDA + invariant-protein median normalization (quantitative AP-MS)",
            "ProDA probabilistic dropout model",
            paste0("Invariant-protein median normalization (median over a data-driven ",
                   "set of unchanged proteins; robust to large one-directional shifts)"),
            "Comparing interactomes across conditions (WT vs. KO, \u00b1 drug, treated vs. untreated)"
          )
        )
      ),

      # ── Q1 = no → Q2b ────────────────────────────────────────────────────────
      conditionalPanel(
        condition = "input.wiz_q1_ip == 'no'",
        tags$div(class = "wiz-q-label",
          "Are any proteins expected to be completely absent (not just lower) in one condition?"),
        tags$div(class = "wiz-card-row",
          .card("wiz_q2b_abs", "absent", "scissors",
                "Yes \u2014 proteins are entirely absent",
                "CRISPR KO, siRNA knockdown, gene deletion",
                q2b_sel),
          .card("wiz_q2b_abs", "quant", "chart-bar",
                "No \u2014 quantitative differences only",
                "Drug treatment, time course, cell line comparison",
                q2b_sel)
        ),
        # Q2b = absent → ProDA + invariant-protein median normalization
        conditionalPanel(
          condition = "input.wiz_q2b_abs == 'absent'",
          .terminal(
            "ProDA + invariant-protein median normalization (knockout / knockdown)",
            "ProDA probabilistic dropout model",
            paste0("Invariant-protein median normalization (median over a data-driven ",
                   "set of unchanged proteins; robust to large one-directional shifts)"),
            "CRISPR knockouts, gene deletions, siRNA knockdowns"
          )
        ),
        # Q2b = quant → ProDA + Global Median
        conditionalPanel(
          condition = "input.wiz_q2b_abs == 'quant'",
          tags$div(
            class = "wiz-terminal-wrap mt-2",
            tags$details(
              class = "wiz-result-details",
              tags$summary(
                tags$span(class = "wiz-result-chip",
                  icon("circle-check", class = "me-1"),
                  "Quantitative analysis",
                  tags$small(class = "ms-2 wiz-chip-hint", "\u25b8 details")
                )
              ),
              tags$div(class = "wiz-result-body small mt-2", no_engine_desc)
            )
          )
        )
      )
    )
  }

  # ── Step 6: Factorial Design (optional; only shown when IP pipeline + >2 groups) ──
  .wiz_step6_factorial_ui <- function() {
    fd <- rv$lf_factorial_design
    tagList(
      tags$h5("Factorial Design"),
      tags$div(class = "alert alert-secondary py-2 px-3 small mb-3",
        icon("circle-info", class = "me-1"),
        tags$strong("This step is optional."),
        " If your experiment has two experimental variables \u2014 for example, ",
        tags$strong("different treatments"), " AND ",
        tags$strong("different time points or genotypes"), " \u2014 you can define a ",
        tags$strong("two-factor (factorial) design"), ". OSMSS will fit a statistical model ",
        "that accounts for both variables simultaneously, improving statistical power by ",
        "separating the effects of each factor.",
        tags$br(), tags$br(),
        "If you have only one variable, leave this unchecked and ",
        "press ", tags$strong("Next"), " to continue with a standard one-factor model."
      ),
      checkboxInput("wiz_factorial_mode",
        tagList(icon("table-cells-large", class = "me-1"),
                tags$strong("Enable two-factor (factorial) design")),
        value = isTRUE(fd$enabled)),
      conditionalPanel(
        condition = "input.wiz_factorial_mode",
        tags$div(class = "border rounded p-3 mt-2",
          tags$small(class = "text-muted d-block mb-3",
            "Define the two experimental factors below. For each group you defined in ",
            "Step 4, assign which level of each factor it belongs to."),
          fluidRow(
            column(6, textInput("wiz_f1_name", "Factor 1 name",
              value = fd$f1_name %||% "ip_type",
              placeholder = "e.g. ip_type")),
            column(6, textInput("wiz_f1_levels",
              HTML("Factor 1 levels <small class='text-muted'>(control first)</small>"),
              value = paste(fd$f1_levels %||% c("IgG", "Bait"), collapse = ", "),
              placeholder = "e.g. IgG, Bait"))
          ),
          fluidRow(
            column(6, textInput("wiz_f2_name", "Factor 2 name",
              value = fd$f2_name %||% "time",
              placeholder = "e.g. time")),
            column(6, textInput("wiz_f2_levels", "Factor 2 levels",
              value = paste(fd$f2_levels %||% character(0), collapse = ", "),
              placeholder = "e.g. D0, D10"))
          ),
          tags$label(class = "form-label small fw-bold mt-2 mb-1",
                     "Assign each group to its factor levels:"),
          uiOutput("wiz_factorial_assignments"),
          tags$hr(class = "my-3"),
          tags$div(class = "alert alert-secondary py-2 px-3 small mt-2 mb-0",
            icon("circle-check", class = "me-1"),
            tags$strong("Interaction model only:"), tags$code(" ~ f1 * f2"),
            tags$br(),
            "OSMSS always fits a full interaction model for factorial designs. ",
            "This means condition-specific comparisons are estimated independently ",
            "\u2014 e.g. \u201cBait D0 vs Bait D10\u201d and \u201cIgG D0 vs IgG D10\u201d are separate contrasts. ",
            "An additive model (pooling effects across groups) is not offered because it ",
            "would obscure biologically distinct responses in IP-MS experiments."
          )
        )
      )
    )
  }

  # ── Per-group factor assignment table ──
  output$wiz_factorial_assignments <- renderUI({
    if (!isTRUE(input$wiz_factorial_mode)) return(NULL)
    groups <- names(rv$lf_cond_groups)
    if (length(groups) == 0)
      return(tags$p(class = "text-muted small", "No groups defined yet."))
    f1_name   <- input$wiz_f1_name   %||% "Factor 1"
    f2_name   <- input$wiz_f2_name   %||% "Factor 2"
    f1_levels <- trimws(strsplit(input$wiz_f1_levels %||% "", ",")[[1]])
    f1_levels <- f1_levels[nzchar(f1_levels)]
    f2_levels <- trimws(strsplit(input$wiz_f2_levels %||% "", ",")[[1]])
    f2_levels <- f2_levels[nzchar(f2_levels)]
    if (length(f1_levels) == 0) f1_levels <- "(define levels above)"
    if (length(f2_levels) == 0) f2_levels <- "(define levels above)"
    fd <- rv$lf_factorial_design
    rows <- lapply(groups, function(g) {
      prev_f1 <- fd$assignments[[g]]$f1
      prev_f2 <- fd$assignments[[g]]$f2
      tags$div(class = "d-flex align-items-center gap-2 mb-1",
        tags$span(class = "text-truncate fw-bold small",
                  style = "width:130px; min-width:130px", g),
        selectInput(paste0("wiz_fact_f1_", make.names(g)),
                    label = NULL,
                    choices  = c("(select)" = "", f1_levels),
                    selected = prev_f1 %||% "",
                    width    = "145px"),
        selectInput(paste0("wiz_fact_f2_", make.names(g)),
                    label = NULL,
                    choices  = c("(select)" = "", f2_levels),
                    selected = prev_f2 %||% "",
                    width    = "145px")
      )
    })
    # Validation hint: check for duplicates
    tagList(
      tags$div(class = "d-flex align-items-center gap-2 mb-1 mt-2",
        tags$span(class = "fw-bold small text-muted",
                  style = "width:130px; min-width:130px", "Group"),
        tags$span(class = "fw-bold small text-muted", style = "width:145px",
                  f1_name, tags$small(class="text-danger ms-1", "(ctrl first)")),
        tags$span(class = "fw-bold small text-muted", style = "width:145px", f2_name)
      ),
      tagList(rows)
    )
  })

  # ── Step 7: Confirm & Summarize ──
  .wiz_step7_summarize_ui <- function() {
    pipeline <- rv$lf_pipeline
    ip_sub   <- rv$lf_ip_sub
    # NULL unless proteins are summarized from PD peptides (Route B): on the
    # protein-direct routes (PD protein abundances, MaxQuant) the row is omitted
    # rather than reporting a filter that could not have run. The stored setting
    # can be left over from an earlier PD file, so the route decides, not the value.
    .sp_lf <- if (identical(rv$wizard_data_types, "peptide")) rv$lf_shared_peptides else NULL
    shared_pep_label <- if (is.null(.sp_lf) || is.na(.sp_lf)) NULL else switch(
      .sp_lf,
      unique_only = "Discarded (peptides matching >1 protein removed)",
      keep_all    = "Kept, attributed to first listed accession",
      razor       = "Assigned by parsimony (most unique-peptide evidence)",
      as.character(.sp_lf))
    tagList(
      tags$h5("Review & Confirm"),
      tags$div(class = "card mb-3",
        tags$div(class = "card-body py-2",
          tags$h6(class = "card-title", "Configuration Summary"),
          tags$table(class = "table table-sm mb-0",
            tags$tbody(
              tags$tr(tags$td(class = "fw-bold", "Source:"),
                      tags$td(switch(rv$wizard_search_program,
                        "pd"       = "Proteome Discoverer",
                        "maxquant" = "MaxQuant",
                        "fragpipe" = "FragPipe",
                        rv$wizard_search_program))),
              tags$tr(tags$td(class = "fw-bold", "Quantified from:"),
                      tags$td(if (identical(rv$wizard_search_program, "fragpipe"))
                        "FragPipe combined_protein.tsv intensities (IonQuant)"
                      else switch(paste(rv$wizard_data_types, collapse = "+"),
                        "protein"  = "Proteome Discoverer protein abundances",
                        "peptide"  = "Peptides, summarized to proteins by OSMSS",
                        "maxquant" = "MaxQuant proteinGroups.txt intensities",
                        paste(rv$wizard_data_types, collapse = " + ")))),
              # Omitted on protein-direct input, where no peptide filter runs.
              if (!is.null(shared_pep_label))
                tags$tr(tags$td(class = "fw-bold", "Shared peptides:"),
                        tags$td(shared_pep_label)),
              tags$tr(tags$td(class = "fw-bold", "Pipeline:"), tags$td(pipeline)),
              tags$tr(tags$td(class = "fw-bold", "Normalization:"),
                      tags$td(
                        if (identical(rv$lf_de_engine, "saint")) "SAINTexpress (internal)"
                        else if (identical(ip_sub, "invariant")) "Invariant-protein median normalization"
                        else if (!lf_apply_normalization_effective(rv))
                          "None in OSMSS (already normalized by IonQuant, per fragpipe.workflow)"
                        else "Global median (OSMSS)")),
              if (!is.null(ip_sub) && identical(pipeline, "ip") && !identical(ip_sub, "invariant"))
                tags$tr(tags$td(class = "fw-bold", "IP Method:"), tags$td(ip_sub)),
              tags$tr(tags$td(class = "fw-bold", "Groups:"),
                      tags$td(paste0(length(rv$lf_cond_groups), ": ",
                        paste(names(rv$lf_cond_groups), collapse = ", "))))
            )
          )
        )
      ),
      tagList(
          tags$p(class = "small text-muted",
            "Click Summarize to process your data. This may take 1\u20133 minutes. ",
            "Once complete, the ", tags$strong("Done"), " button will turn green."),
          actionButton("wiz_summarize_btn",
            tagList(icon("gears"), " Summarize Data"),
            class = "btn-primary w-100 mb-2"),
          uiOutput("wiz_summarize_status")
        )
    )
  }

  output$wiz_summarize_status <- renderUI({ NULL })

  # ── Bait gene picker — populated from gene names already parsed from the data ──
  # ── AP-MS multi-condition group filter (step 6, shown when AP-MS + >2 groups) ──
  output$wiz_apms_group_filter_ui <- renderUI({
    if (rv$wizard_step != 6L) return(NULL)
    grps <- names(rv$lf_cond_groups)
    if (length(grps) <= 2L) return(NULL)
    tagList(
      tags$hr(class = "my-2"),
      tags$div(class = "alert alert-danger py-2 px-3 small mt-2",
        icon("triangle-exclamation", class = "me-2"),
        tags$strong("Statistical Warning: Multiple Conditions in a Bait vs. Control Design"),
        tags$br(), tags$br(),
        "Analyzing multiple conditions simultaneously in a Bait vs. Control design will ",
        "compromise your statistics. In this specific setup, proDA artificially pools variance ",
        "across fundamentally incompatible backgrounds, and SAINT fails to resolve complex ",
        "noise models.",
        tags$br(), tags$br(),
        tags$strong("Action Recommended:"),
        " Please select exactly one Bait and one Control group to analyze for this session."
      ),
      tags$div(class = "border rounded p-3 mt-2 mb-1",
        tags$p(class = "small fw-semibold mb-2",
          icon("filter", class = "me-1"),
          "Select the groups to include in this analysis:"),
        checkboxGroupInput(
          "wiz_apms_keep_groups",
          label    = NULL,
          choices  = grps,
          selected = grps,
          width    = "100%"
        ),
        tags$p(class = "small text-muted mt-1 mb-0",
          icon("circle-info", class = "me-1"),
          "All groups are pre-selected (no change). Deselect unwanted groups to remove them ",
          "before analysis. If you keep all groups, the analysis proceeds as currently configured."
        )
      )
    )
  })

  output$wiz_bait_gene_ui <- renderUI({
    # Extract unique gene names / accessions from whatever is already parsed.
    # Called at Step 6 render time; build_proda_protein_matrix() has NOT run yet.
    gene_choices <- tryCatch({
      choices <- character(0)

      # --- PD nested protein data (most common path) ---
      if (!is.null(rv$pd_raw_lines) && length(rv$pd_protein_line_idx) > 0) {
        pd_raw <- rv$pd_raw_lines
        pd_idx <- rv$pd_protein_line_idx
        block  <- paste(c(pd_raw[1], pd_raw[pd_idx]), collapse = "\n")
        prot_df <- tryCatch(
          read.csv(text = block, header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = c("", "NA")),
          error = function(e) NULL
        )
        if (!is.null(prot_df) && nrow(prot_df) > 0) {
          accessions <- as.character(prot_df[["Accession"]])
          gene_names <- if (!isTRUE(rv$pd_has_gn) && !is.null(rv$pd_gene_col) &&
                             rv$pd_gene_col %in% colnames(prot_df)) {
            as.character(prot_df[[rv$pd_gene_col]])
          } else if ("Gene Symbol" %in% colnames(prot_df)) {
            as.character(prot_df[["Gene Symbol"]])
          } else if ("Gene.Symbol" %in% colnames(prot_df)) {
            as.character(prot_df[["Gene.Symbol"]])
          } else {
            stringr::str_replace(
              stringr::str_extract(prot_df[["Description"]], "GN=([^ ]+)"),
              "^GN=", "")
          }
          gene_names[is.na(gene_names) | gene_names == ""] <-
            accessions[is.na(gene_names) | gene_names == ""]
          choices <- sort(unique(gene_names[nzchar(gene_names)]))
        }
      }

      # --- MaxQuant proteinGroups path ---
      if (length(choices) == 0 && !is.null(rv$mq_proteingroups_path)) {
        mq <- tryCatch(
          utils::read.delim(rv$mq_proteingroups_path, stringsAsFactors = FALSE,
                            nrows = 2000, check.names = FALSE),
          error = function(e) NULL
        )
        if (!is.null(mq)) {
          gn_col <- intersect(c("Gene names", "Gene Names", "Gene.names"), colnames(mq))[1]
          if (!is.na(gn_col)) {
            raw_gn <- unlist(strsplit(mq[[gn_col]], ";"))
            choices <- sort(unique(raw_gn[nzchar(trimws(raw_gn))]))
          } else if ("Fasta headers" %in% colnames(mq)) {
            # No parsed gene column: fall back to GN= tags, as the matrix build does
            raw_gn <- stringr::str_replace(
              unlist(stringr::str_extract_all(as.character(mq[["Fasta headers"]]), "GN=[^ ;]+")),
              "^GN=", "")
            choices <- sort(unique(raw_gn[nzchar(raw_gn)]))
          }
        }
      }

      choices
    }, error = function(e) character(0))

    if (length(gene_choices) > 0) {
      # Offer a searchable dropdown pre-populated with gene names from the data
      selectizeInput(
        "wiz_bait_gene",
        label = tags$span(icon("tag", class = "me-1"), "Bait Gene Name:"),
        choices  = gene_choices,
        selected = {
          cur <- isolate(input$wiz_bait_gene) %||% ""
          if (nzchar(cur) && cur %in% gene_choices) cur else NULL
        },
        options  = list(
          placeholder = "Search or type gene name\u2026",
          maxOptions  = 500L
        ),
        width = "100%"
      )
    } else {
      # Fallback to free-text if gene names cannot be parsed yet
      textInput(
        "wiz_bait_gene",
        label = tags$span(icon("tag", class = "me-1"), "Bait Gene Name:"),
        value       = isolate(input$wiz_bait_gene) %||% "",
        placeholder = "e.g., BRCA1, ALOX5AP, FLAG-MYC",
        width       = "100%"
      )
    }
  })

  # ── Step 8: Batch Correction ──
  .wiz_step8_batch_ui <- function() {
    tagList(
      tags$h5("Batch Correction (Optional)"),
      # PCA chart auto-loaded
      tags$div(class = "alert alert-warning py-2 px-3 mb-2",
        icon("circle-info", class = "me-1"),
        tags$strong("Pre-Modeling QC:"), " Use this PCA to explore your data\u2019s variance.",
        tags$br(),
        tags$strong("Rule of Thumb:"), " Always enable batch correction if samples were physically ",
        "processed in batches, even if they already appear to overlap perfectly here."
      ),
      tags$p(class = "fw-semibold mb-1",
        icon("braille"), " Uncorrected PCA \u2014 colored by condition:"),
      plotlyOutput("wiz_batch_pca", height = "240px"),
      tags$hr(),
      checkboxInput("wiz_use_batch",
        tagList(icon("layer-group"), " Enable batch correction for this dataset"),
        value = FALSE),
      conditionalPanel("input.wiz_use_batch == true",
        tags$div(class = "alert alert-warning py-2 px-3 small mb-3",
          icon("triangle-exclamation", class = "me-1"),
          tags$strong("Complex Designs & Missing Data"), tags$br(),
          "Combining a two-factor design with batch correction is resource-intensive. ",
          "If your dataset has too many missing values across specific batches, ",
          "the statistical model will fail. ",
          tags$strong("If your model fails, you have two fallback options:"),
          tags$br(), tags$br(),
          tags$strong("Option 1: Drop Batch Correction"),
          tags$span(class = "badge bg-success ms-1 small", "Recommended"), tags$br(),
          "Uncheck the box above. You lose the technical correction, but preserve your full ",
          "dataset and protect your biological questions.",
          tags$br(), tags$br(),
          tags$strong("Option 2: Strict Filtering"), tags$br(),
          "If batch effects are massive and must be removed, use the evaluation tool below. ",
          "This stabilizes the math by permanently deleting proteins with sparse coverage.",
          tags$br(),
          actionButton("wiz_show_batch_filter", "Open Option 2 \u2014 filter by batch",
            class = "btn-sm mt-2",
            style = "background-color: transparent; border: 1px solid #000; color: #000;",
            icon = icon("filter"))
        ),
        # ── Filter panel (hidden until user clicks "Try option 3") ──
        conditionalPanel("input.wiz_show_batch_filter % 2 == 1",
          tags$div(class = "border border-warning rounded p-3 bg-light mb-3",
            tags$div(class = "d-flex justify-content-between align-items-start mb-2",
              tags$div(
                tags$strong(icon("filter", class = "me-1"), "Batch-Spanning Filter"),
                tags$p(class = "small text-muted mb-0 mt-1",
                  "Proteins entirely absent in any batch \u00d7 condition cell have inestimable ",
                  "batch coefficients — this causes proDA\u2019s weighted matrix to collapse. ",
                  "The filter requires \u22651 valid value in ", tags$em("every"), " cell. ",
                  "Assign your batches below first, then click ", tags$strong("Preview"), "."
                )
              )
            ),
            actionButton("wiz_batch_filter_preview", "Preview filter impact",
              class = "btn-outline-secondary btn-sm mb-2", icon = icon("magnifying-glass")),
            uiOutput("wiz_batch_filter_result"),
            uiOutput("wiz_batch_filter_checkbox")
          )
        ),
        tags$div(class = "mt-2",
          tags$label(class = "form-label fw-semibold", "Number of batches:"),
          numericInput("wiz_batch_n_batches", NULL,
            value = 2L, min = 2L, max = 20L, step = 1L, width = "120px"),
          tags$p(class = "small text-muted mt-1 mb-2",
            "Assign each sample to a batch below."),
          uiOutput("wiz_batch_assign_rows"),
          tags$div(class = "mt-2",
            uiOutput("wiz_batch_confound_alert")
          )
        )
      ),
    )
  }

  # ── Navigation: Next with per-step validation and state commit ──
  observeEvent(input$wiz_next, {
    step      <- isolate(rv$wizard_step)
    next_step <- step + 1L
    valid <- switch(as.character(step),
      "3" = {
        if (identical(rv$wizard_search_program, "maxquant") &&
            !is.null(rv$mq_proteingroups_path)) {
          if (is.null(rv$mq_evidence_path)) {
            showNotification("Please upload evidence.txt from the same MaxQuant search.",
                             type = "warning"); FALSE
          } else TRUE
        } else if (identical(rv$wizard_search_program, "pd") && !is.null(rv$pd_filename)) {
          TRUE
        } else if (identical(rv$wizard_search_program, "fragpipe") && !is.null(rv$pd_filename)) {
          if (is.null(rv$lf_fp_workflow)) {
            showNotification("Please upload fragpipe.workflow from the same FragPipe run (required).",
                             type = "warning"); FALSE
          } else TRUE
        } else {
          showNotification("Please upload a file first.", type = "warning"); FALSE
        }
      },
      "4" = {
        if (length(rv$lf_cond_groups) < 2) {
          showNotification("Define at least 2 condition groups.", type = "warning")
          FALSE
        } else TRUE
      },
      "5" = TRUE,
      "6" = {
        q1_val <- input$wiz_q1_ip %||% character(0)
        if (!length(q1_val)) {
          showNotification("Please select whether you performed a pulldown / IP.", type = "warning")
          FALSE
        } else if (identical(q1_val, "yes")) {
          q2a_val <- input$wiz_q2a_ctrl %||% character(0)
          if (!length(q2a_val)) {
            showNotification("Please describe the control sample for your pulldown.", type = "warning")
            FALSE
          } else if (identical(q2a_val, "trap")) {
            bait_gene_val <- trimws(input$wiz_bait_gene %||% "")
            if (!nzchar(bait_gene_val)) {
              showNotification("Please enter a bait gene name for SAINT analysis.", type = "warning")
              FALSE
            } else TRUE
          } else TRUE
        } else {
          q2b_val <- input$wiz_q2b_abs %||% character(0)
          if (!length(q2b_val)) {
            showNotification("Please indicate whether proteins are entirely absent in one condition.", type = "warning")
            FALSE
          } else TRUE
        }
      },
      TRUE
    )
    if (!valid) return()

    # Commit state for the completed step (suppress cache invalidation during writes)
    rv$lf_import_in_progress <- TRUE
    on.exit(rv$lf_import_in_progress <- FALSE, add = TRUE)

    switch(as.character(step),
      "3" = {
        # Search program and input mode were set when the file was detected.
        # PD: commit the quantification route. Route A (PD protein abundances)
        # unless the file has no protein rows or the user switched to
        # re-summarizing from the nested peptides (Route B).
        if (identical(rv$wizard_search_program, "pd")) {
          lv       <- rv$pd_levels_available %||% character(0)
          from_pep <- "peptide" %in% lv &&
                      (!("protein" %in% lv) || isTRUE(input$wiz_pd_from_peptides))
          dtypes   <- if (from_pep) "peptide" else "protein"
          rv$wizard_data_types <- dtypes
          rv$lf_input_mode     <- "peptide"
          rv$lf_input_format   <- "pd_nested"
          # Shared-peptide handling. This is the INPUT to the matrix build, which
          # rejects any value outside the three accepted settings, so it must stay
          # valid even on routes that never consult it: hiding a control means it
          # is not editable, not that it has no value. The separate PROVENANCE
          # record is what carries NA on non-peptide routes - the routes
          # themselves set shared_peptide_handling = NA_character_ in their output,
          # which is what keeps the methods text silent there.
          rv$lf_shared_peptides <- if (from_pep)
            (input$wiz_shared_peptides %||% "unique_only") else "unique_only"
        }
      },
      "4" = {
        # ── Collect pairing assignments ──
        if (isTRUE(input$wiz_enable_pairing)) {
          grps_all <- rv$lf_cond_groups
          all_samps <- unlist(grps_all, use.names = FALSE)
          pairing_vec <- setNames(
            vapply(all_samps, function(s) {
              v <- input[[paste0("wiz_pair_", make.names(s))]] %||% ""
              trimws(as.character(v))
            }, character(1L)),
            all_samps
          )
          # Validate: each non-empty subject ID must appear in ≥2 different conditions
          non_empty <- pairing_vec[nzchar(pairing_vec)]
          if (length(non_empty) < length(all_samps)) {
            showNotification(
              "Pairing: some samples have no Subject ID assigned — pairing will be ignored.",
              type = "warning")
            rv$lf_subject_pairing <- NULL
          } else {
            # Check that each subject ID spans multiple conditions
            samp_to_cond <- setNames(
              rep(names(grps_all), lengths(grps_all)), unlist(grps_all, use.names=FALSE))
            subj_conds <- tapply(samp_to_cond[names(non_empty)], non_empty, unique, simplify=FALSE)
            confounded <- any(vapply(subj_conds, length, integer(1)) < 2L)
            if (confounded) {
              showNotification(
                "Pairing: all samples for some subjects are in the same condition — pairing removed.",
                type = "warning")
              rv$lf_subject_pairing <- NULL
            } else {
              rv$lf_subject_pairing <- non_empty
              logger(paste0("[Wizard] Pairing enabled: ", length(unique(non_empty)), " subjects, ",
                            length(non_empty), " samples"))
            }
          }
        } else {
          rv$lf_subject_pairing <- NULL
        }
        # ── Build ProDA annotation tables ──
        grps     <- rv$lf_cond_groups
        pairing  <- rv$lf_subject_pairing
        .biorep <- function(sample_name, group_name, run_idx) {
          if (!is.null(pairing) && sample_name %in% names(pairing) && nzchar(pairing[[sample_name]]))
            pairing[[sample_name]]
          else
            paste0(group_name, "_", run_idx)
        }
        if (identical(rv$wizard_search_program, "maxquant")) {
          ann_list <- list()
          for (gname in names(grps)) {
            for (ri in seq_along(grps[[gname]])) {
              sname <- grps[[gname]][[ri]]
              ann_list[[length(ann_list) + 1]] <- data.frame(
                Raw.file         = sname,
                Condition        = gname,
                BioReplicate     = .biorep(sname, gname, ri),
                Run              = sname,
                IsotopeLabelType = "L",
                stringsAsFactors = FALSE
              )
            }
          }
          rv$mq_annotation <- do.call(rbind, ann_list)
        }
        # ── Factorial step skip logic ──
        # Factorial (step 5) is only applicable when >2 condition groups
        factorial_applicable <- length(rv$lf_cond_groups) > 2L
        reachable <- c(
          3L, 4L,
          if (factorial_applicable) 5L,
          6L, 7L, 8L
        )
        rv$wizard_reachable_steps <- reachable
        rv$wizard_has_batch_step  <- TRUE
        rv$wizard_step <- if (factorial_applicable) 5L else 6L
        return()
      },
      "5" = {
        # ── Step 5: Factorial Design commit ──
        if (isTRUE(input$wiz_factorial_mode)) {
          groups  <- names(rv$lf_cond_groups)
          f1_raw  <- trimws(strsplit(input$wiz_f1_levels %||% "", ",")[[1]])
          f1_levs <- f1_raw[nzchar(f1_raw)]
          f2_raw  <- trimws(strsplit(input$wiz_f2_levels %||% "", ",")[[1]])
          f2_levs <- f2_raw[nzchar(f2_raw)]
          if (length(f1_levs) >= 2 && length(f2_levs) >= 2) {
            assigns <- setNames(
              lapply(groups, function(g) list(
                f1 = input[[paste0("wiz_fact_f1_", make.names(g))]] %||% "",
                f2 = input[[paste0("wiz_fact_f2_", make.names(g))]] %||% ""
              )),
              groups
            )
            rv$lf_factorial_design <- list(
              enabled     = TRUE,
              interaction = TRUE,
              f1_name     = input$wiz_f1_name %||% "ip_type",
              f1_levels   = f1_levs,
              f2_name     = input$wiz_f2_name %||% "time",
              f2_levels   = f2_levs,
              assignments = assigns
            )
          } else {
            rv$lf_factorial_design <- NULL
            showNotification(
              "Factorial design: define at least 2 levels for each factor, or uncheck to skip.",
              type = "warning")
          }
        } else {
          rv$lf_factorial_design <- NULL
        }
        rv$wizard_step <- 6L
        return()
      },
      "6" = {
        # ── Step 6: Experimental Design commit ──
        q1_ip    <- input$wiz_q1_ip    %||% "no"
        q2a_ctrl <- input$wiz_q2a_ctrl %||% character(0)
        q2b_abs  <- input$wiz_q2b_abs  %||% character(0)
        is_ip    <- identical(q1_ip, "yes")
        mnar     <- is_ip || identical(q2b_abs, "absent")
        rv$lf_expect_mnar    <- mnar
        rv$lf_wizard_q1_ip   <- q1_ip
        rv$lf_wizard_q2a_ctrl <- q2a_ctrl
        rv$lf_wizard_q2b_abs  <- q2b_abs
        # ── AP-MS group filter: apply user's selection before the rest of step 6 commits ──
        if (is_ip && length(rv$lf_cond_groups) > 2L) {
          keep_grps <- input$wiz_apms_keep_groups %||% names(rv$lf_cond_groups)
          keep_grps <- intersect(keep_grps, names(rv$lf_cond_groups))
          if (length(keep_grps) >= 2L && length(keep_grps) < length(rv$lf_cond_groups)) {
            removed <- setdiff(names(rv$lf_cond_groups), keep_grps)
            rv$lf_cond_groups      <- rv$lf_cond_groups[keep_grps]
            rv$lf_factorial_design <- NULL   # group composition changed; clear factorial
            logger(paste0("[Wizard] AP-MS group filter: keeping [",
                          paste(keep_grps, collapse = ", "), "]; removed: [",
                          paste(removed, collapse = ", "), "]"))
            showNotification(
              tagList(icon("filter", class = "me-1"),
                      paste0("Analysis restricted to ", length(keep_grps), " groups: ",
                             paste(keep_grps, collapse = " + "))),
              type = "message", duration = 6L
            )
          }
        }
        rv$lf_skip_impute    <- FALSE
        rv$lf_bait_protein   <- NULL
        rv$lf_bait_groups    <- NULL
        if (is_ip) {
          if (identical(q2a_ctrl, "trap")) {
            rv$lf_pipeline     <- "ip"
            rv$lf_ip_sub       <- "saint"
            # ── Factorial gate: SAINTexpress has no covariate support → redirect to ProDA ──
            has_factorial <- !is.null(rv$lf_factorial_design) &&
                             isTRUE(rv$lf_factorial_design$enabled)
            if (has_factorial && isTRUE(.PRODA_AVAILABLE)) {
              rv$lf_de_engine    <- "proda"
              rv$lf_ip_sub       <- "invariant"
              rv$lf_bait_protein <- trimws(input$wiz_bait_gene %||% "")
              logger("[Wizard] AP-MS trap-ctrl + factorial: SAINTexpress skipped, routing to ProDA + invariant-protein median normalization")
              showNotification(
                tagList(
                  tags$strong("Factorial design detected: SAINTexpress bypassed."),
                  tags$br(),
                  "SAINTexpress models a single bait vs. control comparison and ",
                  "cannot accommodate a multi-factor design. ",
                  "Your data will be analysed with ", tags$b("ProDA + invariant-protein median normalization"),
                  " instead, which supports factorial contrast matrices."
                ),
                type = "warning", duration = 12L
              )
            } else {
              rv$lf_de_engine    <- "saint"
              rv$lf_bait_protein <- trimws(input$wiz_bait_gene %||% "")
              logger(paste0("[Wizard] AP-MS trap-ctrl: SAINTexpress pathway. Bait: ",
                            rv$lf_bait_protein))
            }
          } else {
            # depleted or cond → ProDA + invariant-protein median normalization
            rv$lf_pipeline     <- "ip"
            rv$lf_ip_sub       <- "invariant"
            rv$lf_de_engine    <- "proda"
            logger(paste0("[Wizard] AP-MS ctrl=", q2a_ctrl, ": routing to ProDA + invariant-protein median normalization"))
          }
        } else if (identical(q2b_abs, "absent")) {
          # non-IP, proteins entirely absent → invariant-protein median normalization
          # (MNAR-aware ProDA route; the panel's invariance gate is centred on the
          # bulk of complete proteins, see select_invariant_protein_panel)
          rv$lf_pipeline       <- "ip"
          rv$lf_ip_sub         <- "invariant"
          rv$lf_de_engine      <- "proda"
          logger("[Wizard] Non-IP, absent proteins: routing to ProDA + invariant-protein median normalization")
        } else {
          # non-IP, quantitative differences only → Global Median
          rv$lf_pipeline       <- "global"
          rv$lf_ip_sub         <- "no_norm"
          rv$lf_de_engine      <- "proda"
          logger("[Wizard] Non-IP, quantitative: routing to ProDA + Global Median")
        }
        rv$lf_igg_conditions <- character(0)
        # Sync hidden legacy inputs so existing conditionalPanels work
        updateRadioButtons(session,  "upload_pipeline",    selected = rv$lf_pipeline)
        updateRadioButtons(session,  "upload_ip_sub",      selected = rv$lf_ip_sub)
        updateCheckboxInput(session, "upload_skip_impute", value    = FALSE)
        rv$wizard_step <- 7L
        return()
      },
      "7" = {
        # Move to batch correction step
        rv$wizard_step <- 8L
        return()
      }
    )

    rv$wizard_step <- next_step
  })

  # ── Summarize button in wizard step 7 ──
  observeEvent(input$wiz_summarize_btn, {
    use_proda <- identical(rv$lf_de_engine, "proda") && isTRUE(.PRODA_AVAILABLE)
    use_saint <- identical(rv$lf_de_engine, "saint") && isTRUE(.SAINT_AVAILABLE)

    if (identical(lf_wiz_task$status(), "running")) {
      showNotification("Summarization is already running.", type = "warning")
      return()
    }

    output$wiz_summarize_status <- renderUI({
      tags$div(class = "alert alert-info py-2 px-3 small mt-2",
        icon("spinner", class = "fa-spin"), " Summarizing\u2026 Please wait.")
    })

    # --- Determine mode and pre-build protein matrix synchronously (ProDA/SAINT) ---
    mode <- if (use_saint) "saint" else if (use_proda) "proda" else "maxquant"

    # All modes: build_proda_protein_matrix runs synchronously here.
    mat_res       <- NULL
    genes_from_data <- FALSE
    aggregation_method <- NULL
    if (mode %in% c("proda", "saint", "maxquant")) {
      req(rv$lf_cond_groups)
      # %||% only substitutes for NULL. This value can also arrive as
      # NA_character_ (a prerun bundle exported from a non-peptide route stores
      # NA, and %||% passes it through on import), as character(0), or empty -
      # none of which the matrix build accepts. Coerce here rather than
      # loosening validation, which would let a genuinely wrong value pass
      # silently. A valid setting, including "keep_all" restored from a legacy
      # bundle, is passed through untouched.
      mat_res <- build_proda_protein_matrix(
        rv, rv$lf_cond_groups,
        shared_peptides = { .sp <- rv$lf_shared_peptides
                            if (length(.sp) != 1L || is.na(.sp) ||
                                !nzchar(.sp)) "unique_only" else .sp },
        logger = logger)
      if (!is.null(mat_res$error)) {
        output$wiz_summarize_status <- renderUI({
          tags$div(class = "alert alert-danger py-2 px-3 small mt-2",
            icon("circle-xmark"), paste0(" Summarization failed: ", mat_res$error))
        })
        return()
      }
      genes_from_data    <- isTRUE(mat_res$genes_from_data)
      # Same rule as the DE-module build (server_de.R): the protein-direct
      # routes name themselves via input_route; the peptide route
      # aggregates with robustSummary.
      aggregation_method <- mat_res$aggregation_method %||% mat_res$input_route %||% "robustSummary"
      # Carried onto the cache in the result observer; the async task does not see it.
      rv$lf_build_meta <- list(aggregation_method = aggregation_method,
                               mq_stats           = mat_res$mq_stats,
                               # search software and FragPipe workflow settings, for provenance / methods
                               norm_software      = rv$wizard_search_program,
                               fragpipe_settings  = rv$lf_fp_workflow)
      # Shared-peptide provenance. NULL on the routes that aggregate no
      # peptides, so the methods renderer stays silent there.
      rv$lf_shared_pep_stats <- if (is.null(mat_res$shared_peptide_handling) ||
                                    is.na(mat_res$shared_peptide_handling)) NULL else list(
        shared_peptide_handling     = mat_res$shared_peptide_handling,
        n_peptides_total            = mat_res$n_peptides_total,
        n_peptides_shared           = mat_res$n_peptides_shared,
        n_peptides_removed_shared   = mat_res$n_peptides_removed_shared,
        shared_detection_method     = mat_res$shared_detection_method,
        shared_detection_applicable = mat_res$shared_detection_applicable,
        n_peptides_razor_assigned             = mat_res$n_peptides_razor_assigned,
        n_peptides_dropped_no_unique_evidence = mat_res$n_peptides_dropped_no_unique_evidence,
        razor_basis                           = mat_res$razor_basis)
    }

    lf_wiz_cancelled(FALSE)
    show_async_overlay(session,
                       "Summarizing data\u2026",
                       paste0("Mode: ", mode))

    lf_wiz_task$invoke(
      mode               = mode,
      protein_matrix     = if (!is.null(mat_res)) mat_res$protein_matrix else NULL,
      design_df          = if (!is.null(mat_res)) mat_res$design_df      else NULL,
      acc_to_gene        = if (!is.null(mat_res)) mat_res$acc_to_gene    else NULL,
      aggregation_method = aggregation_method,
      genes_from_data    = genes_from_data,
      ip_sub             = rv$lf_ip_sub,
      igg_conditions     = rv$lf_igg_conditions %||% character(0),
      skip_impute        = FALSE,
      apply_normalization = lf_apply_normalization_effective(rv)
    )
  })

  # ── UniProt Gene Name Lookup: Decline ──
  observeEvent(input$uniprot_lookup_decline, {
    removeModal()
    rv$uniprot_lookup_pending <- FALSE
    # Re-open the wizard at the current step so user can proceed
    showModal(.wiz_build_modal())
  })

  # ── UniProt Gene Name Lookup (after summarization) ──
  observeEvent(input$uniprot_lookup_confirm, {
    removeModal()
    if (!isTRUE(rv$uniprot_lookup_pending)) return()
    rv$uniprot_lookup_pending <- FALSE

    cache <- rv$lf_master_cache
    if (is.null(cache) || is.null(cache$acc_to_gene)) {
      showNotification("No protein data available for lookup.", type = "warning")
      return()
    }

    accessions <- names(cache$acc_to_gene)
    showNotification(
      paste0("Looking up gene names for ", length(accessions), " proteins from UniProt..."),
      type = "message", duration = NULL, id = "uniprot_lookup_progress"
    )

    new_a2g <- tryCatch(
      lookup_gene_names_uniprot(accessions, logger = logger),
      error = function(e) {
        showNotification(paste0("UniProt lookup failed: ", e$message), type = "error", duration = 10)
        NULL
      }
    )

    removeNotification("uniprot_lookup_progress")

    if (is.null(new_a2g)) return()

    # Update acc_to_gene in master cache
    rv$lf_master_cache$acc_to_gene <- new_a2g

    # Also update in master_cache_meta if stored there
    if (!is.null(rv$lf_master_cache_meta$acc_to_gene))
      rv$lf_master_cache_meta$acc_to_gene <- new_a2g

    # Update ProteinLevelData gene annotations if present (for downstream modules)
    if (!is.null(cache$processed$ProteinLevelData)) {
      # Some modules use a "Gene" column if present
      pld <- rv$lf_master_cache$processed$ProteinLevelData
      pld$Gene <- new_a2g[pld$Protein]
      pld$Gene[is.na(pld$Gene)] <- pld$Protein[is.na(pld$Gene)]
      rv$lf_master_cache$processed$ProteinLevelData <- pld
    }

    n_resolved <- sum(new_a2g != names(new_a2g))
    showNotification(
      paste0("Gene names resolved: ", n_resolved, " / ", length(new_a2g), " proteins"),
      type = if (n_resolved > 0) "message" else "warning",
      duration = 8
    )
    logger(paste0("UniProt lookup complete: ", n_resolved, " / ", length(new_a2g), " gene names resolved"))

    # Re-open the wizard so the user can continue from the current step
    showModal(.wiz_build_modal())
  })

  # ── Close wizard ──
  # ── Batch step outputs ──

  # Auto-compute uncorrected PCA from cached protein matrix when step 8 loads.
  # Uses cache$processed$ProteinLevelData so the result is identical to the
  # Run PCA module (same data source, same imputation, same prcomp settings).
  # NAs are filled column-wise for prcomp only — ProDA DE never sees these values.
  output$wiz_batch_pca <- renderPlotly({
    req(rv$wizard_step == 8L, rv$lf_master_cache)
    cache <- rv$lf_master_cache

    pld <- tryCatch(cache$processed$ProteinLevelData, error = function(e) NULL)
    req(!is.null(pld))

    pld_wide <- tryCatch(
      tidyr::pivot_wider(pld[, c("Protein", "RUN", "LogIntensities")],
                         names_from = "RUN", values_from = "LogIntensities"),
      error = function(e) NULL
    )
    req(!is.null(pld_wide))

    mat <- as.matrix(pld_wide[, -1L])
    rownames(mat) <- pld_wide$Protein

    # Column-wise median imputation — visualization only, matches main PCA module
    for (ci in seq_len(ncol(mat))) {
      na_m <- is.na(mat[, ci])
      if (any(na_m)) mat[na_m, ci] <- median(mat[!na_m, ci], na.rm = TRUE)
    }
    valid_r   <- rowSums(is.na(mat)) == 0L & rowSums(is.finite(mat)) == ncol(mat)
    mat_clean <- mat[valid_r, , drop = FALSE]
    req(nrow(mat_clean) >= 3L)

    # PCA — center=TRUE, scale.=FALSE, matching main PCA module exactly
    pca_res <- prcomp(t(mat_clean), center = TRUE, scale. = FALSE)
    scores  <- as.data.frame(pca_res$x[, 1:min(2L, ncol(pca_res$x)), drop = FALSE])
    colnames(scores) <- c("PC1", "PC2")[seq_len(ncol(scores))]
    scores$sample    <- rownames(scores)

    # Condition labels from PLD GROUP column (matches main PCA module)
    run_grp          <- unique(pld[, c("RUN", "GROUP")])
    cond_map         <- setNames(as.character(run_grp$GROUP), as.character(run_grp$RUN))
    scores$condition <- cond_map[scores$sample] %||% scores$sample

    # Hover label: original sample name (as in the batch rows below), not "<Group>_Rep<i>"
    ddf <- cache$design_df
    scores$label <- scores$sample
    if (!is.null(ddf) && all(c("sample", "sample_name") %in% colnames(ddf))) {
      orig <- setNames(as.character(ddf$sample), ddf$sample_name)[scores$sample]
      ok   <- !is.na(orig) & nzchar(orig)
      scores$label[ok] <- orig[ok]
    }

    pca_sdev <- pca_res$sdev
    pct      <- round((pca_sdev^2 / sum(pca_sdev^2) * 100)[1:min(2L, length(pca_sdev))], 1)

    plot_ly(scores, x = ~PC1, y = ~PC2, color = ~condition,
            text = ~label, hovertemplate = "%{text}<extra>%{marker.color}</extra>",
            type = "scatter", mode = "markers",
            marker = list(size = 10, opacity = 0.85)) |>
      layout(
        xaxis = list(title = paste0("PC1 (", pct[1], "%)")),
        yaxis = list(title = paste0("PC2 (", if (length(pct) >= 2L) pct[2L] else 0, "%)")),
        legend = list(orientation = "v", x = 1.02, y = 0.5),
        margin = list(l = 50, r = 10, t = 10, b = 40)
      ) |>
      plotly::config(displayModeBar = FALSE)
  })

  # Per-sample batch assignment dropdowns
  output$wiz_batch_assign_rows <- renderUI({
    req(isTRUE(input$wiz_use_batch), rv$lf_master_cache)
    n_bat     <- max(2L, as.integer(input$wiz_batch_n_batches %||% 2L))
    bat_choices <- paste0("Batch ", seq_len(n_bat))
    ddf       <- rv$lf_master_cache$design_df
    req(!is.null(ddf))
    # Show the name the sample had when it was put in a group (design_df$sample);
    # sample_name is the internal "<Group>_Rep<i>" label, whose numbering follows
    # the order of the group's checkbox list. Inputs stay keyed on sample_name.
    orig_of   <- if ("sample" %in% colnames(ddf))
                   setNames(as.character(ddf$sample), ddf$sample_name)
                 else setNames(ddf$sample_name, ddf$sample_name)
    fd        <- rv$lf_factorial_design
    is_fact   <- isTRUE(fd$enabled)
    # Previous assignments (persist across re-renders)
    prev_batch <- isolate(rv$lf_batch_groups)
    prev_map   <- if (!is.null(prev_batch)) {
      unlist(lapply(names(prev_batch), function(b) setNames(rep(b, length(prev_batch[[b]])), prev_batch[[b]])))
    } else character(0)
    header_row <- tags$div(
      class = "d-flex gap-2 align-items-center fw-semibold small border-bottom pb-1 mb-2",
      tags$div(style = "width:160px;", "Sample"),
      tags$div(style = "width:120px;", "Condition"),
      if (is_fact) tags$div(style = "width:90px;", fd$f1_name %||% "Factor 1"),
      if (is_fact) tags$div(style = "width:90px;", fd$f2_name %||% "Factor 2"),
      tags$div(style = "width:110px;", "Batch")
    )
    rows <- lapply(ddf$sample_name, function(sn) {
      cond     <- ddf$condition[ddf$sample_name == sn]
      prev_sel <- prev_map[sn] %||% bat_choices[1L]
      f1_val   <- if (is_fact) fd$assignments[[cond]]$f1 %||% "" else NULL
      f2_val   <- if (is_fact) fd$assignments[[cond]]$f2 %||% "" else NULL
      disp     <- unname(orig_of[sn])
      if (is.na(disp) || !nzchar(disp)) disp <- sn
      tags$div(
        class = "d-flex gap-2 align-items-center mb-1",
        tags$div(class = "small text-truncate", style = "width:160px;",
                 title = paste0(disp, "  (internal: ", sn, ")"), disp),
        tags$div(class = "small badge bg-secondary text-truncate", style = "width:120px; font-size:0.72em;", cond),
        if (is_fact) tags$div(class = "small text-muted", style = "width:90px;", f1_val),
        if (is_fact) tags$div(class = "small text-muted", style = "width:90px;", f2_val),
        tags$div(style = "width:110px;",
          selectInput(
            inputId  = paste0("wiz_batch_sel_", make.names(sn)),
            label    = NULL,
            choices  = bat_choices,
            selected = prev_sel,
            width    = "110px"
          )
        )
      )
    })
    tagList(header_row, rows)
  })

  # Confounding check alert
  output$wiz_batch_confound_alert <- renderUI({
    req(isTRUE(input$wiz_use_batch), rv$lf_master_cache)
    ddf       <- rv$lf_master_cache$design_df
    req(!is.null(ddf))
    n_bat     <- max(2L, as.integer(input$wiz_batch_n_batches %||% 2L))
    bat_choices <- paste0("Batch ", seq_len(n_bat))
    bat_assign <- vapply(ddf$sample_name, function(sn) {
      input[[paste0("wiz_batch_sel_", make.names(sn))]] %||% bat_choices[1L]
    }, character(1L))
    cond_vec <- setNames(ddf$condition, ddf$sample_name)[names(bat_assign)]
    tbl      <- table(bat_assign, cond_vec)
    confounded <- any(colSums(tbl > 0L) < 2L)
    if (confounded) {
      confound_conds <- colnames(tbl)[colSums(tbl > 0L) < 2L]
      tags$div(class = "alert alert-warning py-2 px-3 small",
        icon("triangle-exclamation"), " ",
        tags$strong("Confounding detected: "),
        paste(confound_conds, collapse = ", "),
        " appear in only one batch. Batch will be used for ", tags$em("visualization only"),
        " (not added to the ProDA model)."
      )
    } else {
      tags$div(class = "alert alert-success py-2 px-3 small",
        icon("circle-check"), " ",
        "No confounding detected. Batch will be added to the ProDA model."
      )
    }
  })

  # ── Batch-spanning filter: preview and opt-in ──
  # Helper: compute which proteins survive the batch-spanning filter.
  # Returns list(keep = logical[n_prots], n_total, n_keep, n_drop, cell_table)
  .compute_batch_filter <- function() {
    cache <- rv$lf_master_cache
    if (is.null(cache) || is.null(cache$protein_matrix) || is.null(cache$design_df))
      return(NULL)
    mat <- cache$protein_matrix
    ddf <- cache$design_df
    n_bat <- max(2L, as.integer(input$wiz_batch_n_batches %||% 2L))
    bat_choices <- paste0("Batch ", seq_len(n_bat))
    # Read current batch assignments from the UI
    bat_assign <- vapply(ddf$sample_name, function(sn) {
      input[[paste0("wiz_batch_sel_", make.names(sn))]] %||% NA_character_
    }, character(1L))
    if (any(is.na(bat_assign))) return(NULL)  # assignments not yet set
    # Build batch × condition cell labels per sample
    cond_vec <- setNames(ddf$condition, ddf$sample_name)[names(bat_assign)]
    cell_labels <- paste0(bat_assign, ":", cond_vec)
    unique_cells <- sort(unique(cell_labels))
    n_cells <- length(unique_cells)
    # For each protein: does it have ≥1 valid value in EVERY cell?
    keep <- apply(mat, 1L, function(row_vals) {
      valid_cells <- cell_labels[!is.na(row_vals[names(bat_assign)])]
      length(unique(valid_cells)) >= n_cells
    })
    # Cell table for display
    cell_tbl <- table(bat_assign, cond_vec)
    list(keep = keep, n_total = length(keep), n_keep = sum(keep),
         n_drop = sum(!keep), cell_table = cell_tbl, n_cells = n_cells)
  }

  # reactiveValues to hold filter preview (so checkbox can render without re-clicking)
  .batch_filter_preview <- reactiveVal(NULL)

  observeEvent(input$wiz_batch_filter_preview, {
    result <- .compute_batch_filter()
    .batch_filter_preview(result)
  })

  output$wiz_batch_filter_result <- renderUI({
    result <- .batch_filter_preview()
    if (is.null(result)) return(NULL)
    pct_drop <- round(100 * result$n_drop / max(result$n_total, 1L), 1)
    pct_keep <- round(100 * result$n_keep / max(result$n_total, 1L), 1)
    severity <- if (pct_drop > 50) "alert-danger" else if (pct_drop > 20) "alert-warning" else "alert-info"
    tags$div(class = paste("alert", severity, "py-2 px-3 small mt-2 mb-1"),
      tags$strong(paste0(result$n_total, " proteins total.")),
      tags$br(),
      icon(if (pct_drop > 50) "triangle-exclamation" else "circle-info", class = "me-1"),
      tags$strong(paste0(result$n_drop, " proteins (", pct_drop, "%) would be removed")),
      " \u2014 missing in at least one batch \u00d7 condition cell.",
      tags$br(),
      tags$strong(paste0(result$n_keep, " proteins (", pct_keep, "%) would remain")),
      " with complete batch \u00d7 condition coverage ",
      paste0("(", result$n_cells, " cells).")
    )
  })

  output$wiz_batch_filter_checkbox <- renderUI({
    result <- .batch_filter_preview()
    if (is.null(result)) return(NULL)
    prev_val <- isolate(rv$lf_batch_filter_enabled)
    tagList(
      checkboxInput("wiz_batch_filter_apply",
        tagList(icon("filter", class = "me-1"),
          tags$strong(paste0("Apply filter (keep ", result$n_keep, " of ", result$n_total, " proteins)"))),
        value = isTRUE(prev_val)),
      if (result$n_drop > 0 && result$n_keep > 0)
        tags$p(class = "small text-muted mt-0 mb-0",
          "Filtered proteins will be excluded from the batch-corrected model only. ",
          "If you later disable batch correction, the full matrix is restored.")
      else if (result$n_keep == 0)
        tags$p(class = "small text-danger mt-0 mb-0 fw-bold",
          "No proteins survive this filter. Batch correction is not feasible with this data.")
    )
  })

  observeEvent(input$wiz_close, {
    step <- isolate(rv$wizard_step)
    # Commit batch assignments when closing from step 8
    if (step == 8L && isTRUE(isolate(input$wiz_use_batch))) {
      ddf    <- isolate(rv$lf_master_cache$design_df)
      if (!is.null(ddf)) {
        n_bat      <- max(2L, as.integer(isolate(input$wiz_batch_n_batches) %||% 2L))
        bat_choices <- paste0("Batch ", seq_len(n_bat))
        bat_assign <- vapply(ddf$sample_name, function(sn) {
          isolate(input[[paste0("wiz_batch_sel_", make.names(sn))]]) %||% bat_choices[1L]
        }, character(1L))
        # Build named list: batch_label -> c(sample_names)
        batch_groups <- split(names(bat_assign), bat_assign)
        # Check confounding
        cond_vec   <- setNames(ddf$condition, ddf$sample_name)[names(bat_assign)]
        tbl        <- table(bat_assign, cond_vec)
        confounded <- any(colSums(tbl > 0L) < 2L)
        rv$lf_batch_groups    <- batch_groups
        rv$lf_batch_corrected <- TRUE
        rv$lf_batch_for_model <- !confounded
        rv$lf_batch_filter_enabled <- isTRUE(isolate(input$wiz_batch_filter_apply))
        if (confounded) {
          showNotification(
            "Batch groups are confounded with condition. Batch used for visualization only (not in ProDA model).",
            type = "warning", duration = 8L
          )
        }
      }
    } else if (step == 8L) {
      # Step 8 but Use Batch unchecked
      rv$lf_batch_groups    <- NULL
      rv$lf_batch_corrected <- FALSE
      rv$lf_batch_for_model <- FALSE
      rv$lf_batch_filter_enabled <- FALSE
    }
    removeModal()
    rv$wizard_modal_open <- FALSE
    rv$wizard_completed  <- TRUE

    # ── Set active data source (LF DDA wizard completed last) ──
    # The previous experiment (TMT, DIA, older LF) was already cleared when the
    # file was uploaded (.lf_wiz_begin_new_data -> reset_experiment_state).
    rv$data_source <- "lf_dda"
  })

  # ── Sidebar badges (post-wizard) ──
  output$wiz_sidebar_badges <- renderUI({
    pipeline <- rv$lf_pipeline
    ip_sub   <- rv$lf_ip_sub %||% ""
    n_grps   <- length(rv$lf_cond_groups)
    sp_lab   <- switch(rv$wizard_search_program %||% "",
                       "pd"       = "Proteome Discoverer",
                       "maxquant" = "MaxQuant",
                       "fragpipe" = "FragPipe",
                       "\u2014")
    dt_lab   <- if (length(rv$wizard_data_types) > 0)
      paste(rv$wizard_data_types, collapse = " + ") else "\u2014"
    pl_lab   <- switch(paste0(pipeline, "_", ip_sub),
      "global_"    = "Global (ProDA)",
      "ip_no_norm" = "IP \u2014 Global Median",
      "ip_invariant" = "Invariant-protein median (ProDA)",
      pipeline
    )
    badge <- function(ico, val)
      tags$div(class = "d-flex align-items-center gap-2 py-1 border-bottom",
        tags$span(icon(ico), class = "text-secondary"),
        tags$span(class = "small", val)
      )
    tags$div(class = "card mb-2",
      tags$div(class = "card-body py-2 px-3 small",
        badge("database",    sp_lab),
        badge("file-csv",    dt_lab),
        badge("microscope",  pl_lab),
        badge("layer-group", paste0(n_grps, " group(s)")),
        if (!is.null(rv$lf_master_cache))
          badge("circle-check",
            paste0(rv$lf_master_cache_meta$n_proteins %||% "?", " proteins cached"))
        else
          badge("triangle-exclamation", "Not yet summarized")
      )
    )
  })

  # ── LF Wizard summarization async task ──────────────────────────────────────
  lf_wiz_task <- ExtendedTask$new(function(
      mode, protein_matrix, design_df, acc_to_gene, aggregation_method,
      genes_from_data, ip_sub, igg_conditions, skip_impute, apply_normalization = TRUE) {
    promises::future_promise({
      compute_lf_summarization(
        mode               = mode,
        protein_matrix     = protein_matrix,
        design_df          = design_df,
        acc_to_gene        = acc_to_gene,
        aggregation_method = aggregation_method,
        genes_from_data    = genes_from_data,
        ip_sub             = ip_sub,
        igg_conditions     = igg_conditions,
        skip_impute        = skip_impute,
        apply_normalization = apply_normalization
      )
    }, seed = TRUE)
  })
  lf_wiz_cancelled <- reactiveVal(FALSE)

  # LF Wizard summarization result observer ------------------------------------
  observeEvent(lf_wiz_task$result(), {
    hide_async_overlay(session)
    if (isTRUE(lf_wiz_cancelled())) {
      lf_wiz_cancelled(FALSE)
      output$wiz_summarize_status <- renderUI({ NULL })
      return()
    }
    res <- lf_wiz_task$result()
    for (msg in res$messages %||% character(0)) {
      showNotification(msg, type = "warning", duration = 8)
    }
    if (!isTRUE(res$success)) {
      output$wiz_summarize_status <- renderUI({
        tags$div(class = "alert alert-danger py-2 px-3 small mt-2",
          icon("circle-xmark"), paste0(" Summarization failed: ", res$error))
      })
      return()
    }
    result <- res$result
    n_p    <- result$n_prots_processed

    use_proda <- identical(rv$lf_de_engine, "proda") && isTRUE(.PRODA_AVAILABLE)
    use_saint <- identical(rv$lf_de_engine, "saint") && isTRUE(.SAINT_AVAILABLE)

    rv$lf_master_cache <- result
    # Carried from the synchronous matrix build; the async task does not see it.
    rv$lf_master_cache$shared_peptide_stats <- rv$lf_shared_pep_stats
    # Route provenance, matching what the DE-module build attaches (server_de.R).
    # Only filled when absent, so a value the summarization set itself wins.
    for (.bm in names(rv$lf_build_meta)) {
      if (is.null(rv$lf_master_cache[[.bm]]) && !is.null(rv$lf_build_meta[[.bm]]))
        rv$lf_master_cache[[.bm]] <- rv$lf_build_meta[[.bm]]
    }
    rv$lf_master_cache_meta <- list(
      pipeline       = rv$lf_pipeline,
      ip_sub         = rv$lf_ip_sub,
      skip_impute    = FALSE,
      n_proteins     = n_p,
      de_engine      = if (use_saint) "saint" else "proda",
      format         = if (use_saint) "saint"  else "proda",
      expect_mnar    = isTRUE(rv$lf_expect_mnar),
      pipeline_label = result$pipeline_label,
      bait           = rv$lf_bait_protein %||% NULL
    )

    rv$lf_rep_quality_flag <- res$rep_quality_flag
    rv$lf_master_cache_meta$rep_quality_flag <- res$rep_quality_flag
    rv$lf_cache_stale <- FALSE

    .is_inv_result <- grepl("^invariant_", result$normalization_method %||% "")
    if (isTRUE(result$ebn_fallback)) {
      showNotification(
        tagList(
          tags$strong("Invariant-Protein Normalization Notice: "),
          "No panel of ≥", as.character(.INVARIANT_MIN_REFS),
          " invariant reference proteins could be identified. ",
          "The analysis automatically fell back to global median normalization."
        ),
        type = "warning", duration = 10L
      )
    } else if (.is_inv_result) {
      # Panel applied. Report low confidence and any disagreement between the
      # mode and median centres of the invariance gate (the mode is kept).
      .inv_msgs <- character(0)
      if (identical(result$normalization_method, "invariant_low_confidence"))
        .inv_msgs <- c(.inv_msgs, paste0(
          "The reference panel has ", length(result$ebn_ref_vector %||% character(0)),
          " proteins or was drawn at a relaxed step (ideal ≥",
          as.character(.INVARIANT_MIN_REFS_IDEAL), " at the strict step)."))
      if (isTRUE(result$invariant_mode_ambiguous))
        .inv_msgs <- c(.inv_msgs, paste0(
          "The between-condition differences of complete proteins show competing peaks; ",
          "the invariance gate was centred on the peak nearest zero instead of the tallest peak."))
      if (isTRUE(result$invariant_centre_disagree))
        .inv_msgs <- c(.inv_msgs, paste0(
          "Mode and median of the between-condition differences disagree by ",
          sprintf("%.2f", as.numeric(result$invariant_centre_gap %||% NA)),
          " log₂, which indicates a large one-directional shift. ",
          "The mode was used."))
      if (length(.inv_msgs))
        showNotification(
          tagList(tags$strong("Invariant-Protein Normalization: "),
                  paste(.inv_msgs, collapse = " "),
                  " Results should be interpreted with this in mind."),
          type = "warning", duration = 12L)
    }

    rv$provenance$global <- list(
      input_filename         = if (identical(rv$wizard_search_program, "maxquant"))
        paste0(rv$mq_evidence_filename %||% "evidence.txt", " + ",
               rv$mq_proteingroups_filename %||% "proteinGroups.txt")
      else if (identical(rv$wizard_search_program, "fragpipe"))
        paste0(rv$pd_filename, " + ", rv$lf_fp_workflow_filename %||% "fragpipe.workflow")
      else rv$pd_filename,
      total_initial_proteins = result$n_prots_input,
      input_software_source  = switch(rv$wizard_search_program,
        "pd"       = "Proteome Discoverer",
        "maxquant" = "MaxQuant",
        "fragpipe" = "FragPipe",
        rv$wizard_search_program),
      search_program         = switch(rv$wizard_search_program,
        "pd"       = "Proteome Discoverer",
        "maxquant" = "MaxQuant",
        "fragpipe" = "FragPipe",
        rv$wizard_search_program),
      data_types             = rv$wizard_data_types,
      expect_mnar            = isTRUE(rv$lf_expect_mnar),
      wizard_version         = "2.0"
    )

    flag <- res$rep_quality_flag
    output$wiz_summarize_status <- renderUI({
      flag_ui <- if (!is.null(flag) && !identical(flag$severity, "ok")) {
        pct   <- round(flag$worst_rate * 100)
        samp  <- flag$worst_samp %||% "unknown"
        eff_n <- flag$eff_n
        if (identical(flag$severity, "block")) {
          tags$div(
            class = "alert alert-danger py-2 px-3 small mt-2",
            icon("circle-xmark"), tags$b(" Replicate Quality Issue \u2014 Statistical Analysis Degraded"),
            tags$br(),
            paste0("Sample \u2018", samp, "\u2019 has ", pct, "% missing values \u2014 ",
                   "substantially higher than other replicates. ",
                   "This prevents ProDA and limma from estimating a stable variance prior across proteins, ",
                   "causing all p-values to collapse to 1."),
            tags$br(), tags$br(),
            tags$b("Recommendation: "), "Consider excluding this replicate before re-running. ",
            "As a contingency, if your data includes nested peptide-level rows, ",
            "DE analysis will automatically use a ",
            tags$b("Peptide-Level Background Z-Test (Iterative 3\u03c3 Clipping)"),
            " \u2014 the same architectural approach as Proteome Discoverer\u2019s ",
            "\u2018t-test (Background Based)\u2019 method. If only protein-level data is available, you will be prompted to either ",
            "remove the low-quality replicate or re-export with peptide groups included.",
            tags$br(),
            tags$span(class = "text-muted",
              "Note: This contingency mode strictly requires valid paired peptide observations ",
              "and does not use imputation, making p-values more conservative than Proteome Discoverer\u2019s.")
          )
        } else {
          tags$div(
            class = "alert alert-warning py-2 px-3 small mt-2",
            icon("triangle-exclamation"), tags$b(" Replicate Quality Warning"),
            tags$br(),
            paste0("Sample \u2018", samp, "\u2019 has ", pct, "% missing values. ",
                   "With effective n \u2248 ", eff_n, " per group, variance-based methods may have reduced power. ",
                   "Consider excluding this replicate if results appear non-significant.")
          )
        }
      } else NULL
      tagList(
        tags$div(class = "alert alert-success py-2 px-3 small mt-2",
          icon("circle-check"),
          paste0(" Done! ", n_p, " proteins processed.")),
        flag_ui
      )
    })

    if (!isTRUE(res$genes_from_data)) {
      rv$uniprot_lookup_pending <- TRUE
      showModal(modalDialog(
        title = tags$span(icon("dna"), " Gene Name Lookup"),
        tags$p("No gene names were found in the uploaded data. ",
               "Protein accession numbers will be used as labels."),
        tags$p("Would you like to look up gene names from UniProt? ",
               tags$b("(Requires internet connection)")),
        tags$p(class = "text-muted small",
               "This queries the UniProt ID Mapping API directly. ",
               "Only the protein accession numbers from your data are sent."),
        footer = tagList(
          actionButton("uniprot_lookup_decline", "No, use accessions",
                       icon = icon("xmark")),
          actionButton("uniprot_lookup_confirm", "Yes, look up gene names",
                       class = "btn-primary", icon = icon("magnifying-glass"))
        ),
        size = "m",
        easyClose = FALSE
      ))
    }
  })

  # LF Wizard cancel observer --------------------------------------------------
  observeEvent(rv$async_cancel_requested, {
    if (isTRUE(rv$async_cancel_requested) &&
        identical(lf_wiz_task$status(), "running")) {
      lf_wiz_cancelled(TRUE)
      rv$async_cancel_requested <- FALSE
    }
  })

  # ── Async task cancel relay ──────────────────────────────────────────────────
  # The global cancel button (in the async overlay) sets rv$async_cancel_requested.
  # Each module server observes this flag and handles cancellation for its own task.
  observeEvent(input$async_cancel_btn, {
    rv$async_cancel_requested <- TRUE
    hide_async_overlay(session)
    showNotification("Task cancelled.", type = "warning", duration = 4)
  })

  # Delegate to module servers
  message("[DBG] >> server_residue_prediction")
  server_residue_prediction(input, output, session, rv, logger)
  message("[DBG] >> server_local_ms")
  server_local_ms(input, output, session, rv, logger)
  message("[DBG] >> server_de")
  server_de(input, output, session, rv, logger)
  message("[DBG] >> server_pca")
  server_pca(input, output, session, rv, logger)
  message("[DBG] >> server_heatmap")
  server_heatmap(input, output, session, rv, logger)
  message("[DBG] >> server_protein_coverage")
  server_protein_coverage(input, output, session, rv, logger)
  message("[DBG] >> server_rank_fc")
  server_rank_fc(input, output, session, rv, logger)
  message("[DBG] >> server_trends")
  server_trends(input, output, session, rv, logger)
  message("[DBG] >> server_trend_clusters")
  server_trend_clusters(input, output, session, rv, logger)
  message("[DBG] >> server_fragpipe_db")
  server_fragpipe_db(input, output, session, rv, logger)
  message("[DBG] >> server_tmt_upload")
  server_tmt_upload(input, output, session, rv, logger)
  message("[DBG] >> server_tmt_ptm_upload")
  server_tmt_ptm_upload(input, output, session, rv, logger)
  message("[DBG] >> server_ptm_de")
  server_ptm_de(input, output, session, rv, logger)
  message("[DBG] >> server_ptm_trends")
  server_ptm_trends(input, output, session, rv, logger)
  message("[DBG] >> server_ptm_single")
  server_ptm_single(input, output, session, rv, logger)
  message("[DBG] >> server_ptm_pca")
  server_ptm_pca(input, output, session, rv, logger)
  message("[DBG] >> server_ks_scoring")
  server_ks_scoring(input, output, session, rv, logger)
  message("[DBG] >> server_cross_exp")
  server_cross_exp(input, output, session, rv, logger)
  message("[DBG] >> server_ptm_cross_exp")
  server_ptm_cross_exp(input, output, session, rv, logger)
  message("[DBG] >> server_gsea")
  server_gsea(input, output, session, rv, logger)
  message("[DBG] >> server_dia_upload")
  server_dia_upload(input, output, session, rv, logger)
  message("[DBG] >> server_dia_ptm_upload")
  server_dia_ptm_upload(input, output, session, rv, logger)
  message("[DBG] >> server_prerun")
  server_prerun(input, output, session, rv, logger)
  message("[DBG] >> shinyjs::hide overlay")
  shinyjs::hide("osmss_loading_overlay", anim = TRUE, animType = "fade", time = 0.5)

  # ═══════════════════════════════════════════════════════════════════════
  # Sync Step 0 (Input Type) within module groups
  # ═══════════════════════════════════════════════════════════════════════

  # ── Label Free DDA ──
  local({
    ids <- c("de_input_type", "cov_input_type", "pca_input_type",
             "tcl_input_type", "hm_input_type")
    for (src in ids) {
      local({
        my_src <- src
        others <- setdiff(ids, my_src)
        observeEvent(input[[my_src]], {
          val <- input[[my_src]]
          for (tgt in others)
            updateRadioButtons(session, tgt, selected = val)
        }, ignoreInit = TRUE)
      })
    }
  })

  # ── TMT Analysis ── (input type now centralized in tmt_upload_input_type)

  # ── TMT-PTM ── (input type now centralized in ptm_upload_input_type)

}
