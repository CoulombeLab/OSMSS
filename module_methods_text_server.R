# module_methods_text_server.R ==================================================
# Reusable Shiny server module: renders publication-ready methods text from the
# provenance ledger using glue templates.  Handles NULL values gracefully.
# ==============================================================================

methodsTextServer <- function(id, provenance, section = "de") {
  moduleServer(id, function(input, output, session) {

    # Helper: safe check for a non-empty, non-NULL scalar string
    .has <- function(x) length(x) > 0 && !all(is.na(x)) && any(nzchar(as.character(x)))

    # =========================================================================
    # Shared fragment helpers
    # =========================================================================

    # Data-source sentence for the MaxQuant protein-direct route
    # (aggregation_method "maxquant_lfq_direct"). mq_stats may be absent on
    # older provenance; the counts are then omitted.
    .mq_direct_sentence <- function(mq) {
      .n <- function(x) if (length(x) == 1L && !is.na(x)) format(as.integer(x), big.mark = ",") else NULL
      qt <- if (identical(mq$quant_type, "Intensity")) "raw Intensity" else "LFQ intensity"
      # contaminants_removed is absent on older provenance, where they were always removed
      contam_kept <- isFALSE(mq$contaminants_removed)
      rm_parts <- c(
        if (!contam_kept && !is.null(.n(mq$n_removed_contaminant))) paste0(.n(mq$n_removed_contaminant), " potential contaminant"),
        if (!is.null(.n(mq$n_removed_reverse)))     paste0(.n(mq$n_removed_reverse), " reverse (decoy)"),
        if (!is.null(.n(mq$n_removed_only_by_site))) paste0(.n(mq$n_removed_only_by_site), " only-identified-by-site"))
      kept_sent <- if (contam_kept)
        paste0("Potential contaminant groups (MaxQuant 'Potential contaminant' = +",
               if (!is.null(.n(mq$n_flagged_contaminant))) paste0("; ", .n(mq$n_flagged_contaminant), " group(s)") else "",
               ") were retained; their leading IDs carry the CON__ prefix. ") else ""
      # Match between runs (detected from evidence.txt at upload); silent when unknown
      mbr_sent <- if (isTRUE(mq$mbr_on))
        paste0("Match between runs had been enabled in MaxQuant",
               if (!is.null(.n(mq$n_match_features)) && !is.null(.n(mq$n_evidence_features)))
                 paste0(" (", .n(mq$n_match_features), " of ", .n(mq$n_evidence_features),
                        " evidence features were transferred by matching)") else "",
               ". ")
      else if (isFALSE(mq$mbr_on))
        paste0("Match between runs had not been enabled in MaxQuant (no evidence features were ",
               "transferred by matching), so proteins were quantified only in samples in which at ",
               "least one of their peptides was identified by MS/MS. ")
      else ""
      norm_note <- if (identical(mq$quant_type, "Intensity"))
        paste0("Raw Intensity values (summed peptide intensities, not normalized by MaxQuant) were used ",
               "rather than MaxQuant's LFQ intensities, which are already normalized by the MaxLFQ ",
               "algorithm, so that the data were normalized once, by OSMSS, as described below. ")
      else
        paste0("LFQ intensity values, which are already normalized across runs by MaxQuant's MaxLFQ ",
               "algorithm, were used because the export contained no raw Intensity columns; the OSMSS ",
               "normalization described below was applied in addition. ")
      paste0(
        "Protein-group ", qt, " values were taken directly from the MaxQuant ",
        "proteinGroups.txt output, with each protein group labelled by its leading ",
        "protein; peptide-to-protein assignment, including razor-peptide handling, was ",
        "as configured in MaxQuant and no additional aggregation was performed within OSMSS. ",
        norm_note, mbr_sent,
        if (length(rm_parts) > 0 && !is.null(.n(mq$n_groups_input)))
          paste0("Of ", .n(mq$n_groups_input), " protein groups, ",
                 paste(rm_parts, collapse = ", "), " group(s) were removed. ")
        else if (contam_kept) "Reverse (decoy) hits and only-identified-by-site groups were removed. "
        else "Potential contaminants, reverse (decoy) hits and only-identified-by-site groups were removed. ",
        kept_sent,
        "Non-positive values were treated as missing and intensities were log₂-transformed",
        if (!is.null(.n(mq$n_proteins_output)))
          paste0(", yielding ", .n(mq$n_proteins_output), " protein groups with at least one quantified value")
        else "",
        ".")
    }

    # IP sub-method human label
    .ip_sub_label <- function(ip_sub) {
      switch(as.character(ip_sub %||% "invariant"),
        no_norm   = "Global Median Normalization",
        invariant = "Invariant-Protein Median Normalization",
        saint     = "SAINTexpress (no normalization)",
        as.character(ip_sub)
      )
    }

    # Full pipeline human description
    .pipeline_label_full <- function(pipeline, ip_sub = NULL, skip_imp = FALSE,
                                      de_engine = NULL, expect_mnar = NULL) {
      if (identical(de_engine, "proda") && identical(pipeline, "global")) {
        "quantitative proteome pipeline (ProDA + global median normalization)"
      } else if (identical(pipeline, "global")) {
        paste0("global proteome pipeline (ProDA + median normalization",
               if (isTRUE(skip_imp)) ", no imputation" else "", ")")
      } else if (identical(pipeline, "ip")) {
        paste0("MNAR-aware enrichment pipeline \u2014 ", .ip_sub_label(ip_sub))
      } else {
        pipeline %||% "unspecified pipeline"
      }
    }

    .norm_label_shared <- function(norm_str, ebn_proteins = character()) {
      has_ebn <- length(ebn_proteins) > 0 && !all(is.na(ebn_proteins))
      switch(as.character(norm_str),
        equalizeMedians = "median equalization (equalizeMedians)",
        invariantMedian = {
          base <- paste0(
            "invariant-protein median normalization, wherein invariant reference ",
            "proteins were algorithmically selected based on complete detection ",
            "across all runs, high abundance, low within-condition SD (log₂) and ",
            "low between-condition amplitude evaluated on condition means centred ",
            "on the bulk of completely detected proteins"
          )
          if (has_ebn && length(ebn_proteins) <= 5) {
            paste0(base, " (", paste(ebn_proteins, collapse = ", "), ")")
          } else {
            base
          }
        },
        `FALSE`           = "no normalization",
        `none (IP/Co-IP)` = "no normalization (IP/Co-IP enrichment design)",
        pd_normalized     = "normalization applied within Proteome Discoverer",
        norm_str
      )
    }

    .lod_sentence_shared <- function() {
      paste0(
        "Feature-level LOD imputation was applied before protein-level ",
        "summarization: for each run, a LOD value was calculated as 50% of the ",
        "minimum observed intensity in that run. For proteins with partial ",
        "within-condition dropout (detected in some but not all runs of a condition), ",
        "representative peptide features at the LOD intensity were injected into ",
        "missing runs. This allows ProDA protein-level summarization to propagate realistic ",
        "peptide-level variance into protein estimates."
      )
    }

    .ebn_list_widget <- function(ebn_proteins, ip_sub = NULL) {
      if (is.null(ebn_proteins) || length(ebn_proteins) == 0) return(NULL)
      list(
        tags$hr(style = "margin: 8px 0;"),
        tags$p(tags$strong("Invariant-Protein Median Normalization \u2014 Reference Panel"),
               style = "margin-bottom: 4px; font-size: 0.88em;"),
        tags$p(paste(ebn_proteins, collapse = ", "),
               style = "font-family: monospace; font-size: 0.84em; word-break: break-all;")
      )
    }

    # =========================================================================
    # Template dictionary
    # =========================================================================

    render_de_methods <- function(prov) {
      env <- prov$env
      de  <- prov$de
      # Defined per-branch below; initialised here so every source (dia, tmt,
      # lf, ip...) reaches the shared return with a valid value.
      ebn_list_below <- NULL

      # Guard: nothing to render yet — accept either old analysis_mode or new pipeline field
      if (is.null(de) || (!.has(de$analysis_mode) && !.has(de$pipeline))) return(NULL)

      sentences <- character()

      # Normalise: prefer new pipeline fields, fall back to old
      .de_pipeline  <- de$pipeline  %||% "global"
      .de_ip_sub    <- de$ip_sub    %||% NULL
      .de_skip_imp  <- isTRUE(de$skip_impute)
      .de_analysis  <- de$analysis_mode %||%
                       (de$de_engine %||% "proda")
      .de_am_proda  <- identical(.de_analysis, "proda") || identical(de$de_engine, "proda")
      .de_am_saint  <- identical(.de_analysis, "saint") || identical(de$de_engine, "saint")
      .de_am_welch  <- identical(.de_analysis, "welch") || identical(de$de_engine, "welch")

      # =========================================================
      # PEPTIDE-LEVEL BACKGROUND Z-TEST branch — early return
      # Triggered automatically when a replicate has >30% missing values.
      # Implements PD's "t-test (Background Based)" architecture using
      # peptide-level variance and iterative 3σ clipping.
      # =========================================================
      if (.de_am_welch) {
        rq        <- de$rep_quality_flag %||% NULL
        worst_pct  <- if (!is.null(rq)) round(rq$worst_rate * 100) else "high"
        worst_samp <- if (!is.null(rq)) (rq$worst_samp %||% "one replicate") else "one replicate"
        n_prots    <- de$n_prots_output  %||% de$n_prots_processed %||% "N"
        n_pval_na  <- de$n_prots_pval_na %||% NULL
        n_pep      <- de$n_peptides_used %||% NULL
        bg_sd      <- de$background_sd   %||% NULL
        ctrl_lbl   <- de$ctrl_group %||% "control"
        cond_lbl   <- de$cond_group %||% "condition"

        pep_note <- if (!is.null(n_pep) && !is.null(bg_sd)) {
          paste0(" The global background standard deviation, estimated from ",
                 format(as.integer(n_pep), big.mark = ","),
                 " paired peptide-level log\u2082 fold-changes after iterative 3\u03c3 outlier ",
                 "rejection, was \u03c3 = ", bg_sd, " log\u2082 units.")
        } else ""

        na_note <- if (!is.null(n_pval_na) && as.integer(n_pval_na) > 0L) {
          paste0(" ", format(as.integer(n_pval_na), big.mark = ","),
                 " protein(s) with no valid paired peptide observations (absent in at least ",
                 "one group) were excluded from p-value calculation; fold-changes from the ",
                 if (identical(.de_ip_sub, "invariant")) "invariant-protein-normalised"
                 else "normalised",
                 " protein matrix were still reported for these proteins.")
        } else ""

        text_out <- paste0(
          "Due to elevated missing values in sample \u2018", worst_samp, "\u2019 (",
          worst_pct, "% of proteins absent), global variance-estimation methods ",
          "(proDA, limma empirical Bayes) would produce degenerate statistical results. ",
          "Differential protein abundance was therefore assessed using a ",
          "Peptide-Level Background Z-Test with iterative 3\u03c3 clipping, ",
          "replicating the architecture of Proteome Discoverer\u2019s \u2018t-test (Background Based)\u2019 ",
          "method. ",

          "Briefly, per-peptide log\u2082 fold-changes were computed for all peptides ",
          "with valid observations in both the \u2018", cond_lbl, "\u2019 and \u2018", ctrl_lbl,
          "\u2019 groups (strict requirement: no imputation). ",
          "The global background standard deviation was estimated from these paired ",
          "peptide fold-changes using iterative 3\u03c3 clipping to remove biological outliers ",
          "while preserving the natural technical noise distribution. ",
          "For each protein, a z-statistic was calculated as: ",
          "z = mean(peptide log\u2082FC) \u00d7 \u221a(N\u2095\u2091\u209a) / \u03c3\u2095\u2091\u209a, ",
          "where N\u2095\u2091\u209a is the number of valid paired peptides for that protein and ",
          "\u03c3\u2095\u2091\u209a is the global background standard deviation. ",
          "Two-tailed p-values were derived from the standard normal distribution and ",
          "adjusted for multiple testing using the Benjamini\u2013Hochberg method.",
          pep_note, na_note, " ",
          "This contingency analysis was applied to ", n_prots, " proteins. ",
          "Proteins with no paired peptide evidence receive NA p-values; ",
          "fold-changes are still reported using the ",
          if (identical(.de_ip_sub, "invariant")) "invariant-protein-normalised"
          else "normalised",
          " protein matrix. ",
          "Results should be interpreted with caution given the reduced replicate quality. ",
          "This method requires paired peptide data; protein-only exports cannot provide ",
          "the peptide-level variance estimates needed for this test."
        )
        return(list(text = text_out, ebn_list = NULL))
      }

      # =========================================================
      # SAINT branch — early return with SAINTexpress-specific text
      # =========================================================
      if (.de_am_saint) {
        pipeline_desc <- de$pipeline_label %||% "IP (SAINTexpress)"
        bait_prot     <- de$bait_protein   %||% "unknown"

        sentences <- character()

        # Pipeline opening sentence
        sentences <- c(sentences, glue::glue(
          "Protein-protein interaction scoring was performed using SAINTexpress ",
          "(v3.6.3, intensity module), ",
          "following the {pipeline_desc} pipeline."
        ))

        # Bait and comparison
        sentences <- c(sentences, glue::glue(
          "The bait protein was {bait_prot}. ",
          "The comparison tested was: {de$cond_group} (bait pull-down) vs. ",
          "{de$ctrl_group} (negative control)."
        ))

        # Aggregation / data source
        agg_method <- de$aggregation_method %||% "pd_protein_direct"
        if (identical(agg_method, "pd_protein_direct") && identical(de$norm_software, "fragpipe")) {
          # FragPipe input is converted to PD protein layout and takes the same route
          sentences <- c(sentences, fp_methods_sentence(de$fragpipe_settings))
        } else if (identical(agg_method, "pd_protein_direct")) {
          sentences <- c(sentences, paste0(
            "Protein-level abundances were taken directly from the Proteome ",
            "Discoverer export; peptide- and PSM-level filtering was as configured ",
            "in the Proteome Discoverer consensus workflow. Non-positive values were ",
            "treated as missing and abundances were log\u2082-transformed."))
        } else if (identical(agg_method, "maxquant_lfq_direct")) {
          sentences <- c(sentences, .mq_direct_sentence(de$mq_stats))
        } else {
          sentences <- c(sentences, paste0(
            "Protein-level abundance values were derived from peptide-level data ",
            "by robust summarization (", agg_method, ")."
          ))
        }

        # Pre-SAINT entry filter (fields from server_de.R; absent on older caches)
        if (!is.null(de$n_pre_saint_input)) {
          .n_psi <- as.integer(de$n_pre_saint_input   %||% 0L)
          .n_psd <- as.integer(de$n_dropped_pre_saint %||% 0L)
          sentences <- c(sentences, paste0(
            "Prior to scoring, proteins with fewer than two quantified values ",
            "across all bait and control runs were excluded (",
            format(.n_psd, big.mark = ","), " of ",
            format(.n_psi, big.mark = ","),
            "). No other filtering was applied by OSMSS."))
        }

        # No OSMSS-side normalization or imputation prior to scoring
        sentences <- c(sentences, paste0(
          "No normalization or imputation was applied by OSMSS prior to scoring. ",
          "Interactions were scored with SAINTexpress-int using all control ",
          "purifications as virtual controls (-L set to the number of control ",
          "purifications; -R left at the program default), and are reported as ",
          "SaintScore and Bayesian FDR (BFDR); no Benjamini\u2013Hochberg or ",
          "permutation correction was applied, as SAINTexpress controls the error ",
          "rate internally. SAINTexpress applies its own internal normalization ",
          "across purifications during scoring. No values were imputed for the ",
          "purpose of scoring."))

        # Statistical model
        sentences <- c(sentences, paste0(
          "SAINTexpress models the distribution of prey protein intensities across ",
          "bait and control purifications using a mixture model with a true-interaction ",
          "component and a background contaminant component. ",
          "For each prey protein, SAINTexpress computes a SAINT probability score (SP) ",
          "reflecting the posterior probability of a true interaction, and a Bayesian ",
          "False Discovery Rate (BFDR) calibrated across all scored preys. ",
          "Interactions with BFDR \u2264\u20090.05 were considered statistically significant."
        ))

        # Fold change is reported separately from SAINT scoring (fields from
        # compute_saint_de() provenance; absent on caches predating them).
        if (!is.null(de$fc_floor_source)) {
          .fc_src <- de$fc_floor_source
          .fc_val <- if (is.numeric(de$fc_floor_value))
            round(unname(de$fc_floor_value), 2) else de$fc_floor_value
          .n_ca <- as.integer(de$lod_n_imputed      %||% 0L)
          .n_ba <- as.integer(de$n_bait_absent      %||% 0L)
          .n_sr <- as.integer(de$n_single_replicate %||% 0L)
          sentences <- c(sentences, paste0(
            "Log\u2082 fold changes are reported separately from SAINT scoring and ",
            "were computed from the un-normalized log\u2082 abundances as the ",
            "difference of arm means; they are therefore not on the same basis as ",
            "the SAINT score. For preys detected in the bait arm only, the control ",
            "mean cannot be measured and a detection-limit substitute was used (",
            .fc_src, ", = ", .fc_val, " log\u2082 units); these are substitutions, ",
            "not measurements, are flagged FC_Origin = \"control-absent\", and are ",
            "marked with triangular symbols in the volcano plot (n = ",
            format(.n_ca, big.mark = ","), "). Fold changes are not reported for ",
            "preys absent from the bait arm (n = ", format(.n_ba, big.mark = ","),
            "), because negative enrichment is not interpretable in AP-MS. ",
            "Proteins detected in fewer than two replicates of the bait arm are ",
            "flagged (SingleReplicate, n = ", format(.n_sr, big.mark = ","),
            ") and marked with star symbols; they were not excluded."))
        }

        # Display description (volcano + SaintScore tiers)
        graph_paragraph <- paste(
          "Interaction scoring results are displayed as a volcano plot with Log2 fold",
          "change on the x-axis and \u2212log\u2081\u2080(BFDR) on the y-axis.",
          "Proteins are assigned confidence tiers based on the SAINT probability score:",
          "High (SP\u2009\u22650.80), Medium (SP\u2009\u22650.50), or Low (SP\u2009<\u20090.50).",
          "Proteins detected exclusively in the bait condition and absent from all",
          "controls are labeled Exclusive (SP\u2009\u22650.80), High (SP\u2009\u22650.50\u2013<0.80),",
          "or Medium (SP\u2009<\u20090.50).",
          "PiScore and the confidence tiers are display fields provided to help navigate",
          "the results interactively; they should not be used as reported statistics.",
          "For publication, use the SAINTexpress SaintScore and BFDR together with the",
          "reported log\u2082 fold change."
        )

        # Protein count trajectory
        n_in  <- as.integer(de$n_prots_input    %||% NA)
        n_out <- as.integer(de$n_prots_output   %||% NA)
        if (!is.na(n_in) && !is.na(n_out)) {
          sentences <- c(sentences, paste0(
            "Of the ", format(n_in, big.mark = ","),
            " proteins in the abundance matrix, ",
            format(n_out, big.mark = ","),
            " were returned by SAINTexpress after interaction scoring."
          ))
        }

        # R version
        if (.has(env$r_version)) {
          sentences <- c(sentences, glue::glue(
            "All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."
          ))
        }

        return(list(
          text      = paste(paste(sentences, collapse = " "), graph_paragraph, sep = "\n\n"),
          ebn_list  = NULL,
          empty_msg = "Run a SAINT Differential Expression analysis to generate methods text."
        ))
      }
      # =========================================================
      # End SAINT branch
      # =========================================================

      # =========================================================
      # TMT branch — delegate to the TMT volcano renderer
      # =========================================================
      # compute_tmt_de() writes BOTH prov$de and prov$tmt_volcano. Everything
      # currently routes through the DE tab, so without this branch TMT falls
      # through to the label-free/ProDA default and prints a paragraph
      # describing a pipeline that never ran. Delegating rather than
      # duplicating keeps one source of truth for the TMT text.
      #
      # Deliberately does NOT test prov$tmt_volcano for existence: that ledger
      # entry persists for the whole session, so a DIA run following a TMT run
      # would be misrouted. Only fields written by THIS run are used.
      .is_tmt <- !identical(de$source, "dia") &&
                 (identical(de$source, "tmt") ||
                  identical(de$analysis_mode, "limma_protein") ||
                  identical(de$analysis_mode, "limma_qf"))
      if (.is_tmt) {
        tmt_out <- render_tmt_volcano_methods(prov)
        if (!is.null(tmt_out) && .has(tmt_out$text)) {
          return(list(text      = tmt_out$text,
                      ebn_list  = NULL,
                      empty_msg = "Run a TMT analysis to generate methods text."))
        }
        return(list(
          text      = paste0("TMT analysis detected, but the provenance required to ",
                             "describe it is missing. Re-run the Upload Wizard rather ",
                             "than loading a saved bundle."),
          ebn_list  = NULL,
          empty_msg = "Run a TMT analysis to generate methods text."))
      }
      # =========================================================
      # End TMT branch
      # =========================================================

      # =========================================================
      # DIA branch — limma on protein-level DIA matrix
      # =========================================================
      if (identical(de$source, "dia")) {
        cond_t    <- de$cond_group  %||% "test"
        cond_r    <- de$ctrl_group  %||% "reference"
        n_tested  <- de$n_tested    %||% NA
        n_sig     <- de$n_sig       %||% NA
        n_nmar    <- de$n_nmar      %||% 0L
        is_lod    <- isTRUE(de$perseus_lod)
        is_paired <- isTRUE(de$subject_paired)
        is_batch  <- isTRUE(de$batch_corrected)

        # ── Upload/summarization-stage provenance (written by server_dia_upload.R) ──
        g          <- prov$dia_global
        g_software <- if (!is.null(g) && .has(g$software)) g$software else "Spectronaut"
        g_sw_ver   <- if (!is.null(g) && .has(g$software_version)) paste(g_software, g$software_version) else g_software
        # normalization_source: "user", or the file it was read from (report.log.txt / fragpipe.workflow)
        g_from_log <- !is.null(g) && .has(g$normalization_source) && !identical(g$normalization_source, "user")
        g_log_file <- if (g_from_log && identical(g$normalization_source, "fragpipe.workflow"))
                        "FragPipe’s fragpipe.workflow" else "DIA-NN’s log, report.log.txt"
        g_export   <- if (identical(g_software, "DIA-NN"))
                        paste0(" protein-group matrix (report.pg_matrix.tsv",
                               if (!is.null(g) && .has(g$run_within)) paste0("; DIA-NN run within ", g$run_within),
                               ")")
                      else " protein-level export"
        g_via      <- if (!is.null(g) && .has(g$exported_via))
                        paste0(" (exported through ", g$exported_via, ")") else ""
        g_log_sent <- if (!is.null(g) && isTRUE(g$input_log2))
          "The export held log2 intensities, which were used without a second log transformation. "
        else "Zero values were set to missing and intensities were log2-transformed. "
        g_filter_sent <- if (identical(g_software, "DIA-NN"))
          " DIA-NN filters this matrix at 1% protein-group FDR."
        else if (!is.null(g) && isTRUE(g$qvalue_filtered))
          " The export was reported by the user as having protein-group Q-value filtering (< 0.01) applied upstream."
        else ""
        g_in       <- if (!is.null(g) && .has(g$n_proteins_input))
                        as.integer(g$n_proteins_input) else NA_integer_
        g_out      <- if (!is.null(g) && .has(g$n_proteins_processed))
                        as.integer(g$n_proteins_processed) else NA_integer_

        # Normalization check reported with either branch (dia_normalization_check(), server_dia_upload.R)
        chk_sent <- if (!is.null(g) && .has(g$norm_check_within_sd) && !is.na(g$norm_check_within_sd)) {
          paste0(" Before summarization, per-run medians of the ", as.integer(g$norm_check_n_complete),
                 " proteins quantified in every run were compared between replicate runs of each ",
                 "condition (SD ", formatC(as.numeric(g$norm_check_within_sd), format = "f", digits = 3),
                 " log2; values of 0.02 or less are consistent with cross-run normalized input).")
        } else ""
        norm_sent <- if (!is.null(g) && isFALSE(g$normalization_applied)) paste0(
          "Protein group quantities were taken from the ", g_sw_ver,
          g_export, g_via, ", which had been cross-run normalized in ", g_software,
          if (g_from_log) paste0(" (", g$upstream_norm_mode, " normalization, as recorded in ", g_log_file, "). ")
          else " (the software's default, as indicated in the upload wizard). ",
          g_log_sent,
          "No additional normalization was applied, because a ",
          "second, column-wise median centering of already normalized data can shift the unchanged ",
          "background when a large part of the proteome differs between conditions.",
          chk_sent,
          if (!g_from_log && (identical(g$norm_check_verdict, "not_normalized") ||
                              identical(g$norm_check_verdict, "borderline")))
            " The user retained this setting although the check did not clearly indicate normalized input."
          else "",
          g_filter_sent
        ) else paste0(
          "Protein group quantities were taken from the ", g_sw_ver,
          g_export, g_via, ". ", g_log_sent,
          if (g_from_log)
            paste0("According to ", g_log_file, ", DIA-NN cross-run normalization was switched off (--no-norm).", chk_sent, " ")
          else if (!is.null(g) && isFALSE(g$input_normalized_upstream))
            paste0("The user reported that cross-run normalization had not been applied in ", g_software, ".", chk_sent, " ")
          else "",
          "Column-wise median centering was then applied to ",
          "align global signal levels across runs. Median centering assumes that the ",
          "majority of proteins are unchanged between samples; for comparisons between ",
          "deliberately distinct proteomes (e.g. subcellular or biochemical fractions), ",
          "absolute fold-change magnitudes should be interpreted with this assumption in mind.",
          if (!is.null(g$normalization_used_cc) && !is.na(g$normalization_used_cc)) {
            if (isTRUE(g$normalization_used_cc))
              paste0(" Centering constants were derived from the ",
                     as.integer(g$normalization_n_complete),
                     " proteins quantified in every sample, so that the correction ",
                     "is not biased by sample-specific missingness.")
            else
              paste0(" Too few proteins were quantified in every sample to derive ",
                     "centering constants from complete cases, so constants were ",
                     "computed across all observed values. Samples with more ",
                     "missingness therefore contribute medians drawn from a ",
                     "higher-abundance subset, and the correction retains some ",
                     "dependence on the missingness pattern.")
          } else "",
          g_filter_sent
        )

        # Fix B: the summarization-stage filter runs before any comparison exists
        # and was previously unreported.
        prefilter_sent <- if (!is.na(g_in) && !is.na(g_out)) {
          paste0(
            "The input file contained ", g_in, " protein groups; ", g_out,
            " were retained after requiring at least 50% valid values in at least one ",
            "condition group across the full experiment. Proteins failing this global ",
            "filter were excluded from all downstream analysis, including the ",
            "presence/absence analysis described below.",
            if (!is.null(g$n_dup_removed) && g$n_dup_removed > 0)
              paste0(" A further ", as.integer(g$n_dup_removed),
                     " protein groups sharing a leading accession with an earlier entry ",
                     "were collapsed during parsing.")
            else ""
          )
        } else {
          "[Methods generator: protein-group counts unavailable - summarization provenance missing. Re-run the Upload Wizard rather than loading a saved bundle, or add counts manually.]"
        }

        lod_impute_sent <- if (is_lod) {
          paste0(
            "Perseus-style LOD imputation was applied prior to statistical testing: ",
            "for each sample, missing values were replaced by random draws from a ",
            "downshifted normal distribution (mean \u2212 1.8 \u00d7 SD, width 0.3 \u00d7 SD per sample column). ",
            "Statistical significance was assessed using an s0-adjusted Welch t-test ",
            "(s0\u2009=\u20090.1, following Perseus defaults). ",
            "Significance was determined by permutation-based FDR (250 label permutations; ",
            "the expected false-positive count at each threshold was taken as the median ",
            "null count across permutations), rather than by analytical p-value adjustment."
          )
        } else {
          .k          <- de$knn_k %||% NA
          .mvf        <- round((de$min_valid_frac %||% 0.5) * 100)
          .imax       <- de$max_missing_per_cond %||% NA
          .n_imp      <- as.integer(de$n_imputed  %||% 0L)
          .n_rows     <- as.integer(de$n_knn_rows %||% 0L)
          .n_na       <- as.integer(de$n_left_na  %||% 0L)

          # A fraction resolves to a different integer at every replicate count
          # (>=2 of 3 = 67%, >=2 of 4 = 50%). Report the rule that was actually
          # applied alongside the fraction, but only when every count is known.
          .mn_t       <- de$min_n_test     %||% NA
          .mn_r       <- de$min_n_ref      %||% NA
          .ns_t       <- de$n_test_samples %||% NA
          .ns_r       <- de$n_ref_samples  %||% NA
          .resolved   <- if (!any(is.na(c(.mn_t, .mn_r, .ns_t, .ns_r))))
            paste0(" (\u2265", as.integer(.mn_t), " of ", as.integer(.ns_t),
                   " and \u2265", as.integer(.mn_r), " of ", as.integer(.ns_r), ")")
          else ""

          paste0(
            "Proteins were retained for testing if at least ", .mvf,
            "% of replicates in ",
            if (isTRUE(de$require_both)) "both compared conditions" else "at least one compared condition",
            .resolved,
            " had a valid (non-missing) measurement. ",
            if (!is.na(.imax) && .imax > 0)
              paste0("Missing values were imputed by k-nearest neighbours (k\u2009=\u2009", .k,
                     ", impute::impute.knn) only for proteins with at most ",
                     as.integer(.imax), " missing value",
                     if (as.integer(.imax) == 1L) "" else "s",
                     " within either condition (",
                     .n_imp, " values across ", .n_rows, " proteins). ")
            else "",
            if (.n_na > 0L)
              paste0("The remaining ", .n_na, " missing values were not imputed; ",
                     "limma fits each protein on its available observations. ")
            else "",
            "Differential abundance was assessed using limma",
            if (.has(env$limma_version)) paste0(" (v", env$limma_version, ")") else "",
            " with empirical Bayes ",
            "variance moderation (eBayes with trend\u2009=\u2009TRUE, robust\u2009=\u2009TRUE) and ",
            "Benjamini\u2013Hochberg p-value adjustment. The trend\u2009=\u2009TRUE option explicitly models the ",
            "mean-intensity variance relationship common in label-free mass spectrometry data.",
            # Model scope (dia_fit_scope_guard(); absent in results from before 2026-10-04)
            if (.has(de$fit_scope_reason) && !is.na(de$fit_scope_reason)) {
              .ratio <- if (.has(de$fit_scope_ratio) && !is.na(de$fit_scope_ratio))
                          formatC(as.numeric(de$fit_scope_ratio), format = "f", digits = 2) else NULL
              .thr   <- de$fit_scope_threshold %||% 1.5
              switch(de$fit_scope_reason,
                homogeneous = paste0(" A single linear model was fitted across all ", de$n_conditions_fitted,
                  " conditions (", de$n_samples_fitted, " samples) and the comparison was extracted as a ",
                  "contrast, so the residual variance was pooled across conditions. This was done because ",
                  "the typical within-condition standard deviation differed between conditions by a factor of ",
                  .ratio, ", within the threshold of ", .thr, "."),
                heterogeneous = paste0(" Only the samples of the two compared conditions were fitted, because ",
                  "the typical within-condition standard deviation differed between conditions by a factor of ",
                  .ratio, ", above the threshold of ", .thr, " for pooling the residual variance across all conditions."),
                not_assessable = paste0(" Only the samples of the two compared conditions were fitted, because ",
                  "within-condition variability could not be compared across conditions (fewer than two ",
                  "conditions with at least three samples)."),
                "")
            } else "",
            if (is_paired) " A within-subject paired design was accommodated using limma::duplicateCorrelation." else "",
            if (is_batch)  " Batch effects were included as a fixed covariate in the limma design matrix." else ""
          )
        }

        nmar_sent <- if (!is_lod && as.integer(n_nmar %||% 0L) > 0L) {
          paste0(
            as.integer(n_nmar), " proteins that failed the 50% detection filter were ",
            "excluded from limma and subjected to a presence/absence analysis using ",
            "Fisher\u2019s exact test (Benjamini\u2013Hochberg corrected) to determine whether ",
            "their absence pattern was non-random (NMAR); ",
            "Fisher p-values were corrected within this excluded set only. This correction ",
            "is independent of the correction applied to the limma-tested set, so no joint ",
            "false discovery rate is claimed across the two sets. ",
            if (length(de$fisher_has_power) > 0 && !is.na(de$fisher_has_power)) {
              if (isFALSE(de$fisher_has_power))
                paste0("At this replicate count the smallest attainable two-sided Fisher ",
                       "p-value",
                       if (length(de$fisher_min_achievable_p) > 0 &&
                           !is.na(de$fisher_min_achievable_p))
                         paste0(" (", signif(as.numeric(de$fisher_min_achievable_p), 3), ")")
                       else "",
                       " exceeds 0.05, so no protein in this set can reach significance; ",
                       "the p-values are reported as NA and the set is presented as a ranked ",
                       "descriptive list of detection frequencies and substitution-based ",
                       "fold-change magnitudes. ")
              else ""
            } else {
              paste0("Where the replicate count makes the smallest attainable two-sided ",
                     "Fisher p-value exceed 0.05, the p-values are not interpretable and are ",
                     "reported as NA; the excluded set is then presented as a ranked ",
                     "descriptive list of detection frequencies and substitution-based ",
                     "fold-change magnitudes. ")
            },
            if (!is.na(g_out) && !is.na(as.integer(n_tested %||% NA))) {
              .unacc <- g_out - as.integer(n_tested) - as.integer(n_nmar %||% 0L)
              if (.unacc > 0L)
                paste0("A further ", .unacc, " proteins had no detections in either of the two ",
                       "compared conditions and were not analysed for this contrast.")
              else ""
            } else ""
          )
        } else ""

        disclaimer_sent <- if (!is_lod) {
          paste0(
            "Note: limma applied to DIA protein-level matrices does not explicitly model ",
            "DIA-specific chromatographic uncertainty from precursor scoring and cross-run ",
            "alignment. This trade-off is widely accepted in the field; the majority of ",
            "published DIA proteomics studies apply limma or equivalent linear models ",
            "directly to Spectronaut/DIA-NN protein-level output. ",
            "For workflows requiring precursor-level uncertainty propagation, MSstatsDIA ",
            "provides a fully DIA-aware alternative."
          )
        } else ""

        .n_unt <- as.integer(de$n_untestable %||% 0L)
        n_sent <- if (!is.na(as.integer(n_tested %||% NA))) {
          paste0(as.integer(n_tested), " proteins were tested; ",
                 as.integer(n_sig %||% 0L), " were significant (adj.p\u00a0<\u00a00.05).",
                 if (.n_unt > 0L)
                   paste0(" A further ", .n_unt, " proteins entered the model but had ",
                          "insufficient observations to yield an estimate and are reported ",
                          "without a p-value.")
                 else "")
        } else ""

        text_out <- paste(
          paste0("DIA differential abundance analysis was performed on the comparison: ",
                 cond_t, " vs. ", cond_r, "."),
          norm_sent,
          if (nzchar(prefilter_sent)) prefilter_sent else NULL,
          lod_impute_sent,
          if (nzchar(n_sent)) n_sent else NULL,
          if (nzchar(nmar_sent)) nmar_sent else NULL,
          if (nzchar(disclaimer_sent)) disclaimer_sent else NULL,
          sep = " "
        )

        if (.has(env$r_version))
          text_out <- paste(text_out,
            paste0("All analyses were conducted", if (.has(env$osmss_version)) paste0(" with OSMSS (v", env$osmss_version, ")") else "", " in R (v", env$r_version, ")."))

        return(list(
          text      = trimws(text_out),
          ebn_list  = NULL,
          empty_msg = "Run a DIA Differential Expression analysis to generate methods text."
        ))
      }
      # =========================================================
      # End DIA branch
      # =========================================================

      # =========================================================
      # ProDA branch — early return with ProDA-specific text
      # =========================================================
      if (.de_am_proda) {
        proda_v <- de$proda_version %||% (env$proda_version %||% NULL)
        qf_v    <- de$qfeatures_version %||% (env$qfeatures_version %||% NULL)

        pipeline_desc <- de$pipeline_label %||% "IP (ProDA)"

        pl_sent <- if (!is.null(proda_v) && !is.na(proda_v)) {
          glue::glue(
            "Differential expression analysis was performed using proDA (v{proda_v}), ",
            "a probabilistic model for label-free differential expression analysis, ",
            "following the {pipeline_desc} pipeline."
          )
        } else {
          glue::glue(
            "Differential expression analysis was performed using proDA, ",
            "a probabilistic model for label-free differential expression analysis, ",
            "following the {pipeline_desc} pipeline."
          )
        }
        # Aggregation
        agg_method <- de$aggregation_method  # NULL = unknown / old provenance
        if (identical(agg_method, "robustSummary")) {
          agg_tool <- if (!is.null(qf_v) && !is.na(qf_v)) {
            glue::glue("QFeatures (v{qf_v}) / MsCoreUtils::robustSummary (MASS::rlm)")
          } else {
            "QFeatures / MsCoreUtils::robustSummary (MASS::rlm)"
          }
          agg_n_prot <- if (!is.na(de$n_prots_input %||% NA_integer_)) as.integer(de$n_prots_input) else NULL
          # Shared-peptide handling on the peptide-aggregating route (PD nested
          # peptide). Silent when the slot is absent (older provenance), on the
          # protein-direct routes, and when the
          # mapped protein column carried no candidate lists to act on.
          sp_st <- de$shared_peptide_stats
          shared_pep_str <- if (is.null(sp_st) || !isTRUE(sp_st$shared_detection_applicable) ||
                                length(sp_st$shared_peptide_handling) == 0 ||
                                is.na(sp_st$shared_peptide_handling)) {
            ""
          } else {
            .sp_tot <- if (length(sp_st$n_peptides_total) > 0 &&
                           !is.na(sp_st$n_peptides_total))
                         as.integer(sp_st$n_peptides_total) else NA_integer_
            .sp_rm  <- if (length(sp_st$n_peptides_removed_shared) > 0 &&
                           !is.na(sp_st$n_peptides_removed_shared))
                         as.integer(sp_st$n_peptides_removed_shared) else NA_integer_
            .sp_sh  <- if (length(sp_st$n_peptides_shared) > 0 &&
                           !is.na(sp_st$n_peptides_shared))
                         as.integer(sp_st$n_peptides_shared) else NA_integer_
            .sp_rz  <- if (length(sp_st$n_peptides_razor_assigned) > 0 &&
                           !is.na(sp_st$n_peptides_razor_assigned))
                         as.integer(sp_st$n_peptides_razor_assigned) else NA_integer_
            .sp_dr  <- if (length(sp_st$n_peptides_dropped_no_unique_evidence) > 0 &&
                           !is.na(sp_st$n_peptides_dropped_no_unique_evidence))
                         as.integer(sp_st$n_peptides_dropped_no_unique_evidence) else NA_integer_
            switch(as.character(sp_st$shared_peptide_handling),
              unique_only = paste0(
                " Peptides matching more than one protein were removed before ",
                "protein-level summarization, following the MSstats default ",
                "(useUniquePeptide = TRUE)",
                if (!is.na(.sp_rm) && !is.na(.sp_tot))
                  paste0("; ", format(.sp_rm, big.mark = ","), " of ",
                         format(.sp_tot, big.mark = ","), " peptide row",
                         if (.sp_tot != 1) "s" else "", " ",
                         if (.sp_rm != 1) "were" else "was", " removed")
                else "",
                ". Shared peptides were identified from the list of protein ",
                "accessions reported for each peptide."),
              keep_all = paste0(
                " Peptides matching more than one protein were retained and ",
                "attributed to the first accession listed in the export",
                if (!is.na(.sp_sh) && !is.na(.sp_tot))
                  paste0(" (", format(.sp_sh, big.mark = ","), " of ",
                         format(.sp_tot, big.mark = ","), " peptide row",
                         if (.sp_tot != 1) "s" else "", ")")
                else "",
                ". That attribution is positional and is not evidence-based: it ",
                "is not a razor assignment, and the shared signal is credited to ",
                "one candidate protein rather than apportioned among them."),
              razor = paste0(
                " Peptides matching more than one protein were assigned to a ",
                "single protein by parsimony: each was attributed to the candidate ",
                "protein supported by the greatest number of unique peptide ",
                "sequences across the dataset, with ties resolved alphabetically",
                if (!is.na(.sp_rz))
                  paste0(" (", format(.sp_rz, big.mark = ","), " peptide row",
                         if (.sp_rz != 1) "s" else "", " ",
                         if (.sp_rz != 1) "were" else "was", " reassigned)")
                else "",
                ". Peptides whose candidate proteins had no unique-peptide support ",
                "at all could not be assigned and were removed",
                if (!is.na(.sp_dr))
                  paste0(" (", format(.sp_dr, big.mark = ","), " peptide row",
                         if (.sp_dr != 1) "s" else "", ")")
                else "",
                ". This parsimony assignment was derived from the exported peptide ",
                "table by this pipeline; a peptide-level export carries no razor ",
                "annotation, so the assignment was not made by the search engine."),
              "")
          }
          outcome_note <- if (!is.null(agg_n_prot)) {
            paste0(" Quantitative data were obtained for ",
                   format(agg_n_prot, big.mark = ","), " proteins.")
          } else ""

          sentences <- c(sentences, paste0(
            "Protein abundances were quantified via a peptide\u2192protein summarization ",
            "pipeline using ", agg_tool, ". This approach uses iteratively reweighted ",
            "least squares to down-weight outlier peptides and returns a protein-level NA ",
            "when ALL peptides for a given protein are missing in a sample, thereby ",
            "preserving genuine signal dropout rather than imputing arbitrary values.",
            shared_pep_str, outcome_note
          ))
        } else if (identical(agg_method, "pd_protein_direct") && identical(de$norm_software, "fragpipe")) {
          sentences <- c(sentences, paste(fp_methods_sentence(de$fragpipe_settings),
            "No additional peptide-to-protein aggregation was performed within OSMSS."))
        } else if (identical(agg_method, "pd_protein_direct")) {
          sentences <- c(sentences,
            "Protein-level abundance values were used directly from Proteome Discoverer ",
            "output; no additional peptide-to-protein aggregation was performed within OSMSS."
          )
        } else if (identical(agg_method, "maxquant_lfq_direct")) {
          sentences <- c(sentences, .mq_direct_sentence(de$mq_stats))
        }

        # DE intro and comparison (after data provenance)
        sentences <- c(sentences, pl_sent)

        if (.has(de$ctrl_group) && .has(de$cond_group)) {
          sentences <- c(sentences, glue::glue(
            "The comparison tested was: {de$cond_group} vs. {de$ctrl_group} (control)."
          ))
        }

        # Normalization
        norm_method <- de$normalization_method %||% (de$normalization %||% "none")
        ebn_prots   <- de$ebn_ref_proteins %||% character(0)
        has_ebn     <- length(ebn_prots) > 0 && !all(is.na(ebn_prots))

        norm_sent <- switch(as.character(norm_method),
          "invariant_median" = ,
          "invariant_low_confidence" = {
            .n_inv <- length(ebn_prots[!is.na(ebn_prots)])
            .inv_detail <- if (has_ebn && length(ebn_prots) <= 5) {
              paste0(" (", paste(ebn_prots, collapse = ", "), ")")
            } else if (has_ebn) {
              paste0(" (n = ", length(ebn_prots), " algorithmically selected proteins)")
            } else ""
            .ctr_m <- as.character(de$invariant_centre_method %||% "mode")
            .gap   <- suppressWarnings(as.numeric(de$invariant_centre_gap %||% NA))
            .centre_sent <- paste0(
              "Invariance was judged on condition means centred, per condition, on the ",
              if (identical(.ctr_m, "peak_nearest_zero")) "density peak nearest zero" else "density mode",
              " of the between-condition deviations of all completely detected proteins, ",
              "so that a shift common to most proteins between conditions (for example ",
              "after a knockout, from bait enrichment, or from equal loading of samples of ",
              "different composition) does not exclude truly unchanged proteins from the ",
              "panel. All conditions, including any negative-control conditions, were ",
              "included in the reference-protein search.",
              if (isTRUE(de$invariant_mode_ambiguous))
                paste0(" The distribution of these deviations showed competing peaks, so the ",
                       "peak nearest zero (the smallest implied loading correction) rather ",
                       "than the tallest peak was used as the centre.") else "",
              if (isTRUE(de$invariant_centre_disagree) && !is.na(.gap))
                paste0(" The mode and median of these deviations differed by ",
                       signif(.gap, 2), " log₂, indicating a large one-directional ",
                       "shift; the mode was used.") else "")
            .inv_note <- if (identical(as.character(norm_method), "invariant_low_confidence")) {
              .few   <- .n_inv < .INVARIANT_MIN_REFS_IDEAL
              .loose <- identical(de$ebn_strict, FALSE)
              .stp   <- suppressWarnings(as.integer(de$ebn_step %||% NA))
              .nb <- paste(c(
                if (.few) paste0("fewer than ", .INVARIANT_MIN_REFS_IDEAL,
                                 " reference proteins were identified"),
                if (.loose) paste0("the panel was drawn at relaxation step ",
                                   if (!is.na(.stp)) .stp else "n",
                                   ", after the between-condition amplitude criterion ",
                                   "had been relaxed")),
                collapse = ", and ")
              if (!nzchar(.nb)) .nb <- "the reference panel was flagged as low confidence"
              paste0(" NOTE: ", .nb, ", indicating reduced normalization confidence. ",
                     "Results should be interpreted with caution.")
            } else ""
            paste0(
              "Invariant-protein median normalization was performed using ",
              "proDA::median_normalization with ", .n_inv,
              " invariant reference proteins as spike-in references", .inv_detail, ".",
              .inv_note,
              " These proteins were selected from the data for complete detection across ",
              "all runs, high abundance, low within-condition SD (log₂) and low ",
              "between-condition amplitude. ", .centre_sent
            )
          },
          "invariant_fallback_global_median" = paste0(
            "Invariant-protein median normalization was attempted but fewer than ",
            .INVARIANT_MIN_REFS, " invariant reference proteins were identified. The pipeline ",
            "automatically fell back to global column-median normalization using ",
            "proDA::median_normalization to align per-sample intensity baselines."
          ),
          "global_median" = paste0(
            "Prior to model fitting, global column-median normalization was applied ",
            "using proDA::median_normalization to align per-sample intensity baselines. ",
            "This corrects for systematic inter-sample technical variation without ",
            "requiring a set of spike-in reference proteins."
          ),
          "none_input_normalized" = paste0(
            "No normalization was applied in OSMSS: ",
            if (identical(de$norm_software, "fragpipe"))
              paste0("IonQuant had already normalized the intensities across runs (fragpipe.workflow: ",
                     "ionquant.normalization = 1)")
            else "the input had been normalized across runs upstream",
            "."
          ),
          "none" = paste0(
            "No normalization was applied to the protein-level abundance matrix."
          ),
          paste0("Normalization method: ", norm_method, ".")
        )
        sentences <- c(sentences, norm_sent)

        # Which relaxation step of the reference-panel search was reached, and
        # the thresholds in force at that step. Stated as fact: this is what
        # "invariant" meant on this run.
        .ebn_step_n  <- suppressWarnings(as.integer(de$ebn_step %||% NA))
        if (!is.na(.ebn_step_n) && .ebn_step_n > 0) {
          .ebn_steps_tot <- suppressWarnings(as.integer(de$ebn_n_steps %||% NA))
          .ebn_var_lbl   <- if (!is.null(de$ebn_threshold_sd))
            paste0("within-condition SD \u2264 ",
                   signif(as.numeric(de$ebn_threshold_sd), 3), " log\u2082")
          else if (!is.null(de$ebn_threshold_cv))
            paste0("within-condition |SD/mean| \u2264 ",
                   signif(as.numeric(de$ebn_threshold_cv), 3), "%")
          else NULL
          .ebn_amp_lbl <- if (!is.null(de$ebn_threshold_amplitude) &&
                              !is.na(de$ebn_threshold_amplitude))
            paste0("between-condition amplitude \u2264 ",
                   signif(as.numeric(de$ebn_threshold_amplitude), 3), " log\u2082")
          else NULL
          .ebn_abd_lbl <- if (!is.null(de$ebn_threshold_abundance) &&
                              !is.na(de$ebn_threshold_abundance))
            paste0("global median abundance \u2265 ",
                   signif(as.numeric(de$ebn_threshold_abundance), 4), " log\u2082")
          else NULL
          .ebn_crit <- paste(c(.ebn_var_lbl, .ebn_amp_lbl, .ebn_abd_lbl),
                             collapse = ", ")
          .ebn_n_sel <- suppressWarnings(as.integer(de$ebn_n_selected %||% NA))
          sentences <- c(sentences, paste0(
            "The reference panel was drawn at step ", .ebn_step_n,
            if (!is.na(.ebn_steps_tot)) paste0(" of ", .ebn_steps_tot) else "",
            " of the threshold-relaxation schedule",
            if (nzchar(.ebn_crit)) paste0(", where the criteria in force were ",
                                          .ebn_crit) else "",
            if (!is.na(.ebn_n_sel)) paste0("; ", .ebn_n_sel,
                                           " proteins met them") else "",
            "."))
        }

        # Statistical model
        stat_model_sent <- paste0(
          "ProDA models protein abundance using a linear probabilistic dropout model ",
          "with location and variance moderation across all proteins. Missing values are ",
          "treated as informative non-random events: the probability of non-detection is ",
          "modeled as a sigmoid function of the latent protein abundance, enabling valid ",
          "inference even for proteins with sparse detection. Statistical testing used the ",
          "Wald test on the fitted linear model coefficients, with p-values adjusted for ",
          "multiple testing using the Benjamini-Hochberg (BH) method."
        )

        # Adaptive fit strategy annotation
        fit_strategy <- de$proda_fit_strategy %||% "full"
        if (identical(fit_strategy, "early_stop")) {
          max_iter_used <- de$proda_max_iter_used %||% NA
          iter_part <- if (!is.na(max_iter_used)) paste0(" (", max_iter_used, " iterations)") else ""
          stat_model_sent <- paste0(
            stat_model_sent,
            " The EM parameter estimation was terminated early", iter_part,
            " because numerical instability in the dropout-curve scale parameters ",
            "was detected at a later iteration; the last stable parameter estimates were retained."
          )
        } else if (identical(fit_strategy, "unmoderated")) {
          stat_model_sent <- paste0(
            stat_model_sent,
            " Note: location and variance moderation were disabled as a fallback ",
            "because the moderated model encountered numerical instability ",
            "during EM parameter estimation."
          )
        }

        # Empirical Bayes df correction annotation
        if (isTRUE(de$proda_df_corrected)) {
          stat_model_sent <- paste0(
            stat_model_sent,
            " To ensure strict adherence to Empirical Bayes methodology, ",
            "a structural correction was applied to the proDA Wald test output. While the ",
            "native test_diff function correctly utilizes the posterior moderated variance to ",
            "calculate the standard error, it evaluates the t-statistic using the unmoderated ",
            "residual degrees of freedom, leading to severe statistical incoherence in low-n ",
            "designs. To correct this, Wald test p-values were explicitly recalculated using ",
            "the feature-specific moderated degrees of freedom extracted from the fitted model ",
            "prior to Benjamini-Hochberg FDR control."
          )
        }

        # Determine design formula string used in ProDA (for both batch and factorial sentences)
        .proda_formula <- if (isTRUE(de$factorial_enabled)) {
          if (isTRUE(de$interaction_model))
            if (isTRUE(de$batch_corrected)) "~ Factor 1 * Factor 2 + batch"
            else                            "~ Factor 1 * Factor 2"
          else
            if (isTRUE(de$batch_corrected)) "~ Factor 1 + Factor 2 + batch"
            else                            "~ Factor 1 + Factor 2"
        } else {
          # Include subject if present
          subj_in <- isTRUE(de$subject_pairing_used)
          if (isTRUE(de$batch_corrected) && subj_in) "~ batch + subject + condition"
          else if (isTRUE(de$batch_corrected))        "~ batch + condition"
          else if (subj_in)                           "~ subject + condition"
          else                                        "~ condition"
        }

        # Batch correction in ProDA model
        if (isTRUE(de$batch_corrected)) {
          n_batches <- if (!is.null(de$batch_groups)) length(de$batch_groups) else 2L
          batch_sent <- paste0(
            "Batch effects were modeled as a fixed covariate in the ProDA linear model ",
            "(design formula: ", .proda_formula, "), treating batch as an additive offset. ",
            "This parameterisation allows the EM algorithm to estimate and partial out ",
            "systematic between-batch shifts in the mean dropout curve, yielding ",
            "batch-adjusted fold-change estimates and test statistics."
          )
          if (!is.null(de$batch_groups) && length(de$batch_groups) >= 2) {
            # Members are the original sample names, as supplied when groups were defined
            members <- vapply(names(de$batch_groups), function(b)
              paste0(b, " (", paste(de$batch_groups[[b]], collapse = ", "), ")"), character(1))
            batch_sent <- paste0(batch_sent,
              " A total of ", n_batches, " batches were defined, with samples assigned as ",
              "follows: ", paste(members, collapse = "; "), ".")
          }
          stat_model_sent <- paste0(stat_model_sent, " ", batch_sent)
        }

        # Subject/paired design notation
        if (isTRUE(de$subject_pairing_used) && !isTRUE(de$factorial_enabled)) {
          n_subj <- de$subject_n %||% NA_integer_
          subj_sent <- paste0(
            "A paired (repeated-measures) experimental design was used. Subject identity ",
            "was included as a fixed covariate in the ProDA design formula (",
            .proda_formula, ") to account for within-subject correlation across conditions",
            if (!is.na(n_subj)) paste0(" (", n_subj, " subjects)") else "",
            "."
          )
          sentences <- c(sentences, subj_sent)
        }

        # Factorial design notation
        if (isTRUE(de$factorial_enabled)) {
          f1 <- de$factorial_f1 %||% "Factor 1"
          f2 <- de$factorial_f2 %||% "Factor 2"
          model_type <- if (isTRUE(de$interaction_model)) "interaction" else "additive"
          effect_desc <- if (isTRUE(de$interaction_model))
            "including main effects and their cross-product term to capture condition-specific responses across levels of both factors"
          else
            "estimating independent main effects for each factor"
          sentences <- c(sentences, paste0(
            "A two-factor (factorial) experimental design was specified, with factors ",
            tags$strong(f1), " and ", tags$strong(f2), ". ",
            "ProDA was fitted using a ", model_type, " model (", .proda_formula, "), ",
            effect_desc, "."
          ))
        }

        sentences <- c(sentences, stat_model_sent)

        # The variance prior is estimated from whatever samples entered the fit.
        if (identical(de$variance_prior_scope %||% "", "pairwise_subset")) {
          .n_cd <- suppressWarnings(as.integer(
            de$n_conditions_dropped_from_model %||% NA))
          sentences <- c(sentences, paste0(
            "Each pairwise comparison was fitted on the samples belonging to that 
comparison only",
            if (!is.na(.n_cd) && .n_cd > 0)
              paste0(" (", .n_cd, " other condition",
                     if (.n_cd != 1) "s were" else " was", " excluded from the model)")
            else "",
            ". The empirical Bayes variance prior is therefore specific to each 
contrast, and effect sizes, moderated variances and adjusted p-values from 
different contrasts of the same dataset are not directly comparable with one 
another."))
        }

        # Contrast-level minimum evidence.
        .n_ce <- suppressWarnings(as.integer(de$n_contrast_empty_removed %||% NA))
        .n_le <- suppressWarnings(as.integer(de$n_low_evidence %||% NA))
        if (!is.na(.n_ce)) {
          sentences <- c(sentences, paste0(
            "Before fitting, proteins with no observation in either compared 
condition were removed, since they carry no information about the contrast",
            if (.n_ce > 0) paste0(" (", format(.n_ce, big.mark = ","),
                                  " protein", if (.n_ce != 1) "s" else "",
                                  " removed)") else " (none on this dataset)",
            ". Proteins with fewer than two observations in at least one of
the two compared arms were retained, because proDA models incomplete rows
explicitly, and are flagged in the results table with the LowEvidence column",
            if (!is.na(.n_le)) paste0(" (", format(.n_le, big.mark = ","),
                                      " so flagged)") else "",
            ". The flag is descriptive only and does not alter the model fit,
fold changes or p-values; for these proteins the estimate for the sparse arm
relies largely on the proDA dropout model and should be interpreted with
caution."))
        }

        # Protein-count trajectory: QFeatures → ProDA core filter → pre-filter
        n_prots_qf        <- as.integer(de$n_prots_input    %||% de$starting_proteins %||% NA)
        n_prots_processed <- as.integer(de$n_prots_processed %||% NA)
        traj_sent <- character()

        # Step 1: QFeatures → ProDA minimum-detection filter
        if (!is.na(n_prots_qf) && !is.na(n_prots_processed) && n_prots_processed < n_prots_qf) {
          n_core_removed <- n_prots_qf - n_prots_processed
          traj_sent <- paste0(
            "Of the ", format(n_prots_qf, big.mark = ","),
            " proteins in the quantitative matrix, ",
            format(n_prots_processed, big.mark = ","),
            " passed proDA\u2019s minimum-detection threshold (",
            format(n_core_removed, big.mark = ","), " excluded)."
          )
        } else if (!is.na(n_prots_qf)) {
          traj_sent <- paste0(
            "A total of ", format(n_prots_qf, big.mark = ","),
            " quantified proteins were submitted to the analysis pipeline."
          )
        }

        if (nchar(traj_sent) > 0)
          sentences <- c(sentences, traj_sent)

        # R version
        if (.has(env$r_version)) {
          sentences <- c(sentences, glue::glue(
            "All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."
          ))
        }

        graph_paragraph <- paste(
          "Differential expression results are displayed as a volcano plot of log2 fold",
          "change against \u2212log10(BH-adjusted p-value): proteins with adjusted p < 0.05",
          "and log2 fold change > 1 or < \u22121 are shown in red (up) or blue (down),",
          "other proteins in grey, and LowEvidence proteins (fewer than two",
          "observations in at least one arm) as grey squares regardless of their",
          "fold change or p-value."
        )

        return(list(
          text      = paste(paste(sentences, collapse = " "), graph_paragraph, sep = "\n\n"),
          ebn_list  = if (has_ebn && length(ebn_prots) > 5) ebn_prots else NULL,
          ebn_ip_sub = .de_ip_sub,
          empty_msg = "Run a Differential Expression analysis to generate methods text."
        ))
      }
      # =========================================================
      # End ProDA branch — fall through to PD-Direct
      # =========================================================

      # --- 1. Pipeline description ---
      pipeline_desc <- if (!is.null(de$pipeline_label) && nzchar(de$pipeline_label)) {
        de$pipeline_label
      } else {
        .pipeline_label_full(.de_pipeline, .de_ip_sub, .de_skip_imp)
      }

      # --- 2. Software & version ---
      sw <- if (.has(de$software_used)) de$software_used else "an unspecified tool"
      sw_ver <- ""
      if (identical(sw, "limma") && .has(env$limma_version)) {
        sw_ver <- glue::glue(" (v{env$limma_version})")
      }
      sentences <- c(sentences, glue::glue(
        "Differential expression analysis was performed using {sw}{sw_ver} with the {pipeline_desc}."
      ))

      # --- Ctrl / Cond group mention ---
      if (.has(de$ctrl_group) && .has(de$cond_group)) {
        sentences <- c(sentences, glue::glue(
          "The comparison tested was: {de$cond_group} vs. {de$ctrl_group} (control)."
        ))
      }


      # --- 3b. Batch correction (single-trunk, dual-routing architecture) ---
      if (isTRUE(de$batch_corrected)) {
        n_batches <- if (!is.null(de$batch_groups)) length(de$batch_groups) else NULL
        bc_count  <- if (!is.null(n_batches)) paste0(" (", n_batches, " batches)") else ""
        sentences <- c(sentences, paste0(
          "Batch correction was performed at the protein level", bc_count,
          ". All features were first summarized to protein-level intensities by ProDA ",
          "in a single trunk, without any feature-level batch correction."
        ))
        sentences <- c(sentences, paste0(
          "For differential expression testing, batch effects were modeled as a fixed ",
          "covariate in a limma empirical Bayes moderated t-test (eBayes with robust = TRUE ",
          "and trend = TRUE) with a design matrix of the form ~ 0 + Condition + Batch. The ",
          "robust option protects LOD-imputed zero-variance proteins from aggressive variance ",
          "shrinkage, and the trend option models the mean-variance relationship commonly ",
          "observed in mass spectrometry data."
        ))
        sentences <- c(sentences, paste0(
          "For visualization (PCA plots, heatmaps), protein-level batch correction was ",
          "applied using limma::removeBatchEffect on the TMP-summarized protein matrix, ",
          "preserving biological condition differences via the design matrix."
        ))
      }

      # --- 4. Starting proteins ---
      if (length(de$starting_proteins) > 0 && !is.na(de$starting_proteins)) {
        sentences <- c(sentences, glue::glue(
          "A total of {as.integer(de$starting_proteins)} quantified proteins were submitted to the analysis pipeline."
        ))
      }

      # --- 5. Dropped proteins with gate reason ---
      if (length(de$dropped_proteins) > 0 && !is.na(de$dropped_proteins) &&
          de$dropped_proteins > 0) {
        if (.has(de$drop_reason)) {
          sentences <- c(sentences, glue::glue(
            "Prior to statistical testing, {as.integer(de$dropped_proteins)} proteins were excluded: {de$drop_reason}."
          ))
        } else {
          sentences <- c(sentences, glue::glue(
            "Prior to statistical testing, {as.integer(de$dropped_proteins)} proteins were removed due to insufficient evidence."
          ))
        }
      }

      # --- 6. Imputation-mode details (limma pipeline) ---
      if (identical(de$analysis_mode, "imputation")) {
        imp_parts <- character()
        if (.has(de$norm_strategy)) {
          imp_parts <- c(imp_parts, glue::glue(
            "Abundance values were log2-transformed and normalized using {de$norm_strategy} normalization"
          ))
        }
        if (.has(de$imputation_method)) {
          imp_parts <- c(imp_parts, glue::glue(
            "remaining missing values were imputed using a {de$imputation_method} strategy"
          ))
        }
        if (length(imp_parts) > 0) {
          sentences <- c(sentences, paste0(
            paste(imp_parts, collapse = "; "), "."
          ))
        }
        if (.has(de$pvalue_adjustment)) {
          sentences <- c(sentences, glue::glue(
            "Statistical significance was evaluated using limma's empirical Bayes moderated t-test (eBayes), with p-values adjusted for multiple testing using the Benjamini-Hochberg (BH) method."
          ))
        }
      }

      # --- 8. Significance thresholds ---
      has_fdr <- length(de$fdr_threshold) > 0 && !is.na(de$fdr_threshold)
      has_fc  <- length(de$logfc_threshold) > 0 && !is.na(de$logfc_threshold)
      if (has_fdr && has_fc) {
        sentences <- c(sentences, glue::glue(
          "Significance was assessed at an adjusted p-value threshold of {de$fdr_threshold} and a log2 fold-change cutoff of \u00b1{de$logfc_threshold}."
        ))
      } else if (has_fdr) {
        sentences <- c(sentences, glue::glue(
          "Significance was assessed at an adjusted p-value threshold of {de$fdr_threshold}."
        ))
      }

      # --- 9. R version ---
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue(
          "All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."
        ))
      }

      # --- 10. Graph description (static) ---
      graph_paragraph <- paste(
        "Differential expression results are displayed as a standard volcano plot with \u2212log10(adjusted p-value) on the Y-axis and log2 fold-change on the X-axis."
      )

      list(text = paste(paste(sentences, collapse = " "), graph_paragraph, sep = "\n\n"),
           ebn_list = ebn_list_below,
           ebn_ip_sub = .de_ip_sub,
           empty_msg = "Run a Differential Expression analysis to generate methods text.")
    }

    # =========================================================================
    # PCA methods template
    # =========================================================================

    render_pca_methods <- function(prov) {
      env <- prov$env
      pc  <- prov$pca
      de  <- prov$de

      if (is.null(pc) || !.has(pc$pipeline)) return(NULL)

      # ── SAINT branch: PCA data came from the SAINTexpress cache ──────────────
      # Primary check: the explicit is_saint_pca flag written by server_pca.R.
      # Fallback: consult DE provenance only when is_saint_pca is NULL (legacy
      # cache that predates the flag), NOT when it is explicitly FALSE — stale
      # SAINT DE provenance from a prior run must not override a fresh ProDA PCA.
      .de_am_saint_pca <- isTRUE(pc$is_saint_pca) ||
        (is.null(pc$is_saint_pca) && !is.null(de) &&
           identical(de$analysis_mode %||% de$de_engine, "saint"))
      if (.de_am_saint_pca) {
        # Complete-case selection arrived with the is_saint_pca flag. Caches
        # predating it ran column-median imputation, so only claim complete
        # cases when the flag is explicitly present.
        .pca_complete_case <- isTRUE(pc$is_saint_pca)
        bait_prot   <- de$bait_protein   %||% NULL
        ctrl_group  <- de$ctrl_group     %||% NULL
        cond_group  <- de$cond_group     %||% NULL
        sw          <- if (.has(pc$software_used)) pc$software_used else "prcomp (base R)"
        n_prots     <- pc$n_proteins %||% NA
        n_samps     <- pc$n_samples  %||% NA
        pc1_v       <- pc$pc1_var    %||% NA
        pc2_v       <- pc$pc2_var    %||% NA
        sentences <- character()
        sentences <- c(sentences, paste0(
          "Principal component analysis (PCA) was performed using ", sw,
          " on the log\u2082-transformed protein intensity matrix from the ",
          "IP (SAINTexpress) pipeline. ",
          "The analysis included ", if (!is.na(n_samps)) n_samps else "n",
          " samples and ",
          if (!is.na(n_prots)) n_prots else "n",
          if (.pca_complete_case)
            paste0(" proteins retained after complete-case selection, i.e. proteins ",
                   "quantified in every sample; no values were imputed for the PCA.")
          else
            paste0(" proteins retained after filtering for at least 50% quantification ",
                   "across samples; remaining missing values were replaced with the ",
                   "per-sample median prior to PCA.")
        ))
        # Include bait / group detail only when DE provenance is available
        if (!is.null(bait_prot) && !is.null(cond_group) && !is.null(ctrl_group)) {
          sentences <- c(sentences, paste0(
            "Samples comprised ", cond_group, " (bait) and ", ctrl_group,
            " (control), with bait protein ", bait_prot, "."
          ))
        }
        sentences <- c(sentences,
          "No normalisation was applied to the protein matrix prior to PCA; ",
          "values reflect the log\u2082-transformed intensities used for SAINT scoring ",
          "(SAINTexpress incorporates its own internal normalisation model)."
        )
        sentences <- c(sentences, paste0(
          "Data were submitted to prcomp with centering (mean subtracted) and no ",
          "unit-variance scaling",
          if (!is.na(pc1_v) && !is.na(pc2_v))
            paste0(". PC1 explained ", pc1_v, "% and PC2 explained ", pc2_v,
                   "% of the total variance")
          else "",
          "."
        ))
        if (.has(env$r_version))
          sentences <- c(sentences,
            paste0("All analyses were conducted", if (.has(env$osmss_version)) paste0(" with OSMSS (v", env$osmss_version, ")") else "", " in R (v", env$r_version, ").")
          )
        return(list(
          text      = paste(sentences, collapse = " "),
          ebn_list  = NULL,
          empty_msg = "Run PCA to generate methods text."
        ))
      }
      # ── End SAINT branch ───────────────────────────────────────────────────

      sentences      <- character()
      ebn_list_below <- NULL

      .pc_skip_imp <- isTRUE(pc$skip_impute)

      # 1. Software + pipeline description
      sw <- if (.has(pc$software_used)) pc$software_used else "prcomp (base R)"
      pipeline_desc <- if (!is.null(pc$pipeline_label) && nzchar(pc$pipeline_label %||% "")) {
        pc$pipeline_label
      } else {
        .pipeline_label_full(pc$pipeline, pc$ip_sub %||% NULL, .pc_skip_imp)
      }
      sentences <- c(sentences, glue::glue(
        "Principal component analysis (PCA) was performed using {sw} following the {pipeline_desc}."
      ))

      # 2. Matrix dimensions
      has_n  <- length(pc$n_proteins) > 0 && !is.na(pc$n_proteins)
      has_s  <- length(pc$n_samples)  > 0 && !is.na(pc$n_samples)
      if (has_n && has_s) {
        sentences <- c(sentences, glue::glue(
          "The analysis included {as.integer(pc$n_samples)} samples and {as.integer(pc$n_proteins)} proteins."
        ))
      }

      # 3. Dropped proteins
      n_drop_field <- pc$n_proteins_dropped %||% (
        if (!is.null(pc$n_prots_input) && !is.null(pc$n_proteins))
          as.numeric(pc$n_prots_input) - as.numeric(pc$n_proteins)
        else NULL
      )
      has_drop <- !is.null(n_drop_field) && !is.na(n_drop_field) && as.numeric(n_drop_field) > 0
      if (has_drop) {
        sentences <- c(sentences, glue::glue(
          "{as.integer(n_drop_field)} proteins were removed prior to PCA due to missing values across all samples."
        ))
      }

      # 4. Preprocessing pipeline
      if (.has(pc$normalization)) {
        ebn_proteins <- if (!is.null(pc$ebn_ref_proteins)) pc$ebn_ref_proteins else character()
        has_ebn      <- length(ebn_proteins) > 0 && !all(is.na(ebn_proteins))
        norm_label   <- .norm_label_shared(pc$normalization, ebn_proteins)

        ip_sub_str <- if (identical(pc$pipeline, "ip") && !is.null(pc$ip_sub %||% NULL) && nzchar(pc$ip_sub %||% "")) {
          paste0(" (", .ip_sub_label(pc$ip_sub), ")")
        } else ""
        pipeline_display <- switch(pc$pipeline,
          ip     = paste0("IP/Co-IP (pull-down) enrichment", ip_sub_str),
          global = "global proteome",
          pc$pipeline
        )
        sm_str <- if (.has(pc$summary_method)) {
          sm <- if (identical(pc$summary_method, "TMP")) "Tukey's Median Polish (TMP)" else pc$summary_method
          glue::glue("protein summarization by {sm}")
        }
        parts <- c(
          glue::glue("{pipeline_display} pipeline"),
          if (!is.null(sm_str)) sm_str,
          norm_label,
          if (isTRUE(pc$mb_impute)) "model-based imputation (MBimpute)"
        )
        sentences <- c(sentences, paste0(
          "Data were processed through the ",
          paste(parts, collapse = ", "), "."
        ))
        if (has_ebn && length(ebn_proteins) > 5) ebn_list_below <- ebn_proteins
      }

      # 4a. IgG exclusion from calibration (PCA)
      pc_igg_conds <- pc$igg_conditions %||% character(0)
      if (length(pc_igg_conds) > 0 && !all(is.na(pc_igg_conds))) {
        sentences <- c(sentences, paste0(
          "IgG/empty negative control condition(s) (",
          paste(pc_igg_conds, collapse = ", "),
          ") were excluded from reference-panel calibration to prevent corruption of the ",
          "normalization baseline by fundamentally different protein profiles."
        ))
      }

      # 5. LOD imputation / skip
      if (.pc_skip_imp) {
        sentences <- c(sentences, "Imputation was disabled for this analysis.")
      } else if (isTRUE(pc$lod_imputation)) {
        sentences <- c(sentences, .lod_sentence_shared())
      }

      # 5b. Complete-case selection on IP/AP-MS routes. Absence is signal in
      # these experiments, so QC PCA is run on proteins quantified in every
      # sample rather than on an imputed matrix.
      if (isTRUE(pc$ip_complete_case)) {
        .pc_ncc <- suppressWarnings(as.integer(pc$n_prots_complete_case %||% NA))
        sentences <- c(sentences, paste0(
          "Because this is an immunoprecipitation/AP-MS experiment, in which 
non-detection carries information, the QC PCA was computed on complete cases 
only\u2014proteins quantified in every sample",
          if (!is.na(.pc_ncc)) paste0(" (", format(.pc_ncc, big.mark = ","),
                                      " proteins)") else "",
          ". No values were imputed for this analysis."))
      }

      # 5c. Column-median fill on the PCA path. These cells reach prcomp
      # directly -- the row filter downstream retains them -- and, when batch
      # correction is enabled, the batch model as well.
      .pc_pbn <- suppressWarnings(as.integer(pc$prebatch_n_imputed %||% NA))
      if (!is.na(.pc_pbn) && .pc_pbn > 0) {
        sentences <- c(sentences, paste0(
          "PCA requires a complete matrix, so the ",
          format(.pc_pbn, big.mark = ","), " value",
          if (.pc_pbn != 1) "s" else "",
          " still missing at that stage ",
          if (.pc_pbn != 1) "were" else "was",
          " replaced with the median of the observed values in the same sample ",
          "(per-sample median) before the principal component decomposition",
          if (isTRUE(pc$batch_corrected))
            ", and before batch correction was applied" else "",
          ". These substituted values are retained in the matrix submitted to ",
          "prcomp and contribute to the component loadings."))
      }

      # 6. PCA settings
      settings <- c(
        if (length(pc$centering) > 0 && !is.na(pc$centering))
          if (isTRUE(pc$centering)) "centering (mean subtracted)" else "no centering",
        if (length(pc$scaling) > 0 && !is.na(pc$scaling))
          if (isTRUE(pc$scaling)) "unit-variance scaling" else "no unit-variance scaling"
      )
      if (length(settings) > 0) {
        sentences <- c(sentences, glue::glue(
          "Data were submitted to prcomp with {paste(settings, collapse = ' and ')}."
        ))
      }

      # 7. Variance explained
      has_pc1 <- length(pc$pc1_var) > 0 && !is.na(pc$pc1_var)
      has_pc2 <- length(pc$pc2_var) > 0 && !is.na(pc$pc2_var)
      if (has_pc1 && has_pc2) {
        sentences <- c(sentences, glue::glue(
          "PC1 explained {round(pc$pc1_var, 1)}% and PC2 explained {round(pc$pc2_var, 1)}% of the total variance."
        ))
      }

      # 8. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = ebn_list_below,
           ebn_ip_sub = pc$ip_sub %||% NULL,
           empty_msg = "Run a PCA analysis to generate methods text.")
    }

    # =========================================================================
    # Trend Clusters methods template
    # =========================================================================

    render_trend_methods <- function(prov) {
      env <- prov$env
      tc  <- prov$trend_clusters

      # DIA runs write source = "dia" and no $pipeline field, so the guard below
      # returned NULL and no methods text was produced for DIA trend clusters.
      # TMT is routed the same way on the line below: the Trend Clusters tab
      # instantiates this module with section = "trend_clusters", so the
      # dedicated TMT renderer was otherwise unreachable. Both branches key on
      # source only, so LF and IP/SAINT are unaffected.
      if (identical(tc$source %||% "", "dia")) return(render_dia_trends_methods(prov))
      if (identical(tc$source %||% "", "tmt")) return(render_tmt_trends_methods(prov))

      if (is.null(tc) || !.has(tc$pipeline)) return(NULL)

      sentences      <- character()
      ebn_list_below <- NULL
      .tc_skip_imp   <- isTRUE(tc$skip_impute)

      # 1. Algorithm + pipeline context
      algo <- if (.has(tc$clustering_method)) tc$clustering_method else "fuzzy c-means (e1071::cmeans)"
      pipeline_desc <- if (!is.null(tc$pipeline_label) && nzchar(tc$pipeline_label %||% "")) {
        tc$pipeline_label
      } else {
        .pipeline_label_full(tc$pipeline, tc$ip_sub %||% NULL, .tc_skip_imp)
      }
      sentences <- c(sentences, glue::glue(
        "Protein abundance trends were analyzed using fuzzy c-means clustering ({algo}) as implemented in the R package e1071, following the {pipeline_desc}."
      ))

      # 2. Cluster count & fuzzifier
      has_k <- length(tc$n_clusters) > 0 && !is.na(tc$n_clusters)
      has_m <- length(tc$fuzzifier)  > 0 && !is.na(tc$fuzzifier)
      if (has_k) {
        m_str <- if (has_m) round(tc$fuzzifier, 3) else "auto"
        sentences <- c(sentences, glue::glue(
          "A total of {as.integer(tc$n_clusters)} clusters were specified, with a fuzzifier exponent m\u00a0=\u00a0{m_str}."
        ))
        sentences <- c(sentences, paste0(
          "Clusters are numbered by their mean direction at the last condition, 
so the numbering is display order only: it depends on the condition order 
chosen in the interface and is not stable between runs or between two 
orderings of the same conditions. Cluster numbers should not be cited as 
identifiers."))
      }

      # 3. Protein counts
      has_in <- length(tc$n_proteins_input)     > 0 && !is.na(tc$n_proteins_input)
      has_np <- length(tc$n_proteins)           > 0 && !is.na(tc$n_proteins)
      has_cl <- length(tc$n_proteins_clustered) > 0 && !is.na(tc$n_proteins_clustered)
      n_out  <- (tc$n_proteins_clustered %||% tc$n_proteins) %||% NULL
      if (has_in && !is.null(n_out)) {
        if (as.integer(n_out) != as.integer(tc$n_proteins_input)) {
          sentences <- c(sentences, glue::glue(
            "{as.integer(tc$n_proteins_input)} proteins were submitted; {as.integer(n_out)} were retained after removing those with incomplete condition coverage."
          ))
        } else {
          sentences <- c(sentences, glue::glue(
            "{as.integer(n_out)} proteins were submitted for clustering."
          ))
        }
      }

      # 4. Preprocessing
      if (.has(tc$normalization)) {
        ebn_proteins <- if (!is.null(tc$ebn_ref_proteins)) tc$ebn_ref_proteins else character()
        has_ebn      <- length(ebn_proteins) > 0 && !all(is.na(ebn_proteins))
        norm_label   <- .norm_label_shared(tc$normalization, ebn_proteins)

        ip_sub_str <- if (identical(tc$pipeline, "ip") && !is.null(tc$ip_sub %||% NULL)) {
          paste0(" \u2014 ", .ip_sub_label(tc$ip_sub))
        } else ""
        pipeline_display <- if (grepl("^ip$", tc$pipeline, ignore.case = TRUE)) {
          paste0("IP/Co-IP", ip_sub_str)
        } else {
          "global proteome"
        }
        sm_str <- if (.has(tc$summary_method)) {
          sm <- if (identical(tc$summary_method, "TMP")) "Tukey's Median Polish (TMP)" else tc$summary_method
          glue::glue("protein summarization by {sm}")
        }
        parts <- c(
          glue::glue("{pipeline_display} pipeline"),
          if (!is.null(sm_str)) sm_str,
          norm_label,
          if (isTRUE(tc$mb_impute)) "model-based imputation (MBimpute)"
        )
        sentences <- c(sentences, paste0(
          "Prior to clustering, protein abundances were quantified using ProDA with the ",
          paste(parts, collapse = ", "), "."
        ))
        if (has_ebn && length(ebn_proteins) > 5) ebn_list_below <- ebn_proteins
      }

      # 4a. IgG exclusion from calibration
      tc_igg_conds <- tc$igg_conditions %||% character(0)
      if (length(tc_igg_conds) > 0 && !all(is.na(tc_igg_conds))) {
        sentences <- c(sentences, paste0(
          "IgG/empty negative control condition(s) (",
          paste(tc_igg_conds, collapse = ", "),
          ") were excluded from reference-panel calibration to prevent corruption of the ",
          "normalization baseline by fundamentally different protein profiles."
        ))
      }

      # 5. LOD imputation / skip
      if (.tc_skip_imp) {
        sentences <- c(sentences, "Imputation was disabled for this analysis.")
      } else if (isTRUE(tc$lod_imputation)) {
        n_lod <- tc$lod_n_imputed %||% NULL
        lod_count_str <- if (!is.null(n_lod) && !is.na(n_lod) && as.integer(n_lod) > 0)
          paste0("; ", as.integer(n_lod), " missing value(s) imputed") else ""
        sentences <- c(sentences, paste0(
          sub("\\.", paste0(lod_count_str, "."), .lod_sentence_shared(), fixed = TRUE)
        ))
      }

      # 6. Z-score
      if (isTRUE(tc$zscore_transform)) {
        sentences <- c(sentences,
          "Protein abundance profiles were z-score normalized across conditions prior to clustering to emphasize relative abundance patterns."
        )
      }

      # 7. Distance metric
      if (.has(tc$distance_metric)) {
        sentences <- c(sentences, glue::glue("Clustering used {tc$distance_metric} distance."))
      }

      # 8. Conditions
      if (length(tc$condition_names) > 0 && any(nzchar(tc$condition_names))) {
        sentences <- c(sentences, glue::glue(
          "Conditions analyzed: {paste(tc$condition_names, collapse = ', ')}."
        ))
      }

      # 9. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = ebn_list_below,
           ebn_ip_sub = tc$ip_sub %||% NULL,
           empty_msg = "Run a Trend Clusters analysis to generate methods text.")
    }

      # =========================================================================
      # Heatmap methods template
      # =========================================================================

      render_heatmap_methods <- function(prov) {
        env <- prov$env
        hm  <- prov$heatmap
        de  <- prov$de

      # DIA runs write source = "dia" into this same ledger entry but carry no
      # $pipeline field, so without this the guard below returned NULL and the
      # DIA heatmap produced no methods text at all. TMT is routed the same way
      # on the line below: the Heatmap tab instantiates this module with
      # section = "heatmap", so the dedicated TMT renderer was otherwise
      # unreachable. Both branches key on source only, so LF and IP/SAINT are
      # unaffected.
      if (identical(hm$source %||% "", "dia")) return(render_dia_heatmap_methods(prov))
      if (identical(hm$source %||% "", "tmt")) return(render_tmt_heatmap_methods(prov))

      if (is.null(hm) || !.has(hm$pipeline)) return(NULL)

      # ── SAINT branch: DE was run via SAINTexpress ──────────────────────────
      .de_am_saint_hm <- !is.null(de) &&
        identical(de$analysis_mode %||% de$de_engine, "saint")
      if (.de_am_saint_hm) {
        bait_prot  <- de$bait_protein %||% "unknown"
        n_shown    <- hm$n_proteins_displayed %||% NA
        n_total    <- hm$n_proteins           %||% NA
        row_sc     <- hm$row_scaling          %||% "z-score"
        col_cl     <- hm$col_clustering       %||% "hierarchical"
        conds      <- if (length(hm$condition_names %||% character()) > 0)
                        paste(hm$condition_names, collapse = ", ") else "all conditions"
        sentences  <- character()
        sentences <- c(sentences, paste0(
          "Protein abundance data from the IP (SAINTexpress) pipeline ",
          "(bait protein: ", bait_prot, ") were visualized as a hierarchically ",
          "clustered heatmap. ",
          "Abundance values are log\u2082-transformed intensities from the same ",
          "per-sample total-intensity normalized matrix used for SAINT scoring."
        ))
        if (!is.na(n_shown)) {
          sel_desc <- if (identical(hm$selection_method, "gene_filter"))
            paste0(n_shown, " proteins were selected by manual gene filter.")
          else
            paste0("The top ", hm$top_n %||% n_shown,
                   " most variable proteins were selected; ",
                   n_shown, " proteins are shown.")
          sentences <- c(sentences, sel_desc)
        }
        sentences <- c(sentences, paste0(
          "Rows were scaled by ", row_sc, ". ",
          "Columns were ordered by ", col_cl, " clustering. ",
          "Conditions displayed: ", conds, "."
        ))
        if (.has(env$r_version))
          sentences <- c(sentences,
            paste0("All analyses were conducted", if (.has(env$osmss_version)) paste0(" with OSMSS (v", env$osmss_version, ")") else "", " in R (v", env$r_version, ").")
          )
        return(list(
          text      = paste(sentences, collapse = " "),
          ebn_list  = NULL,
          empty_msg = "Run SAINT DE analysis first to generate heatmap methods text."
        ))
      }
      # ── End SAINT branch ───────────────────────────────────────────────────

      sentences      <- character()
      ebn_list_below <- NULL
      .hm_skip_imp   <- isTRUE(hm$skip_impute)

      # 1. Overview with pipeline description
      pipeline_desc <- if (!is.null(hm$pipeline_label) && nzchar(hm$pipeline_label %||% "")) {
        hm$pipeline_label
      } else {
        .pipeline_label_full(hm$pipeline, hm$ip_sub %||% NULL, .hm_skip_imp)
      }
      sentences <- c(sentences, glue::glue(
        "Protein abundance data were visualized as a hierarchically clustered heatmap following the {pipeline_desc}."
      ))

      # 2. Protein selection
      has_sel <- .has(hm$selection_method)
      has_n   <- length(hm$n_proteins_displayed) > 0 && !is.na(hm$n_proteins_displayed)
      if (has_sel) {
        if (identical(hm$selection_method, "gene_filter")) {
          n_filt <- length(hm$gene_filter)
          sentences <- c(sentences, glue::glue(
            "{as.integer(hm$n_proteins_displayed)} protein(s) were displayed following manual gene selection ({n_filt} gene(s) specified)."
          ))
        } else {
          sentences <- c(sentences, glue::glue(
            "The top {as.integer(hm$top_n)} most variable proteins (by row variance) were selected for display, resulting in {as.integer(hm$n_proteins_displayed)} proteins shown."
          ))
        }
      } else if (has_n) {
        sentences <- c(sentences, glue::glue(
          "{as.integer(hm$n_proteins_displayed)} proteins were displayed."
        ))
      }

      # 3. Condition averaging
      has_avg <- length(hm$avg_by_condition) > 0 && !is.na(hm$avg_by_condition)
      if (has_avg && isTRUE(hm$avg_by_condition)) {
        sentences <- c(sentences,
          "Abundances were averaged within each experimental condition for display."
        )
      }

      # 4. Preprocessing
      if (identical(hm$pipeline, "pd_norm")) {
        sentences <- c(sentences,
          "Protein abundances were pre-normalized within Proteome Discoverer prior to import."
        )
      } else if (.has(hm$normalization)) {
        ebn_proteins <- if (!is.null(hm$ebn_ref_proteins)) hm$ebn_ref_proteins else character()
        has_ebn      <- length(ebn_proteins) > 0 && !all(is.na(ebn_proteins))
        norm_label   <- .norm_label_shared(hm$normalization, ebn_proteins)

        ip_sub_str <- if (identical(hm$pipeline, "ip") && !is.null(hm$ip_sub %||% NULL) && nzchar(hm$ip_sub %||% "")) {
          paste0(" \u2014 ", .ip_sub_label(hm$ip_sub))
        } else ""
        pipeline_display <- switch(as.character(hm$pipeline),
          ip     = paste0("IP/Co-IP enrichment", ip_sub_str),
          global = "global proteome",
          hm$pipeline
        )
        sm_str <- if (.has(hm$summary_method)) {
          sm <- if (identical(hm$summary_method, "TMP")) "Tukey's Median Polish (TMP)" else hm$summary_method
          glue::glue("protein summarization by {sm}")
        }
        parts <- c(
          glue::glue("{pipeline_display} pipeline"),
          if (!is.null(sm_str)) sm_str,
          norm_label,
          if (isTRUE(hm$mb_impute)) "model-based imputation (MBimpute)"
        )
        sentences <- c(sentences, paste0(
          "Protein abundances were quantified using ProDA with the ",
          paste(parts, collapse = ", "), "."
        ))
        if (has_ebn && length(ebn_proteins) > 5) ebn_list_below <- ebn_proteins
      }

      # 4a. IgG exclusion from calibration (heatmap)
      hm_igg_conds <- hm$igg_conditions %||% character(0)
      if (length(hm_igg_conds) > 0 && !all(is.na(hm_igg_conds))) {
        sentences <- c(sentences, paste0(
          "IgG/empty negative control condition(s) (",
          paste(hm_igg_conds, collapse = ", "),
          ") were excluded from reference-panel calibration to prevent corruption of the ",
          "normalization baseline by fundamentally different protein profiles."
        ))
      }

      # 5. LOD imputation / skip
      if (.hm_skip_imp && !identical(hm$pipeline, "pd_norm")) {
        sentences <- c(sentences, "Imputation was disabled for this analysis.")
      } else if (isTRUE(hm$lod_imputation)) {
        n_lod <- hm$lod_n_imputed %||% NULL
        lod_count_str <- if (!is.null(n_lod) && !is.na(n_lod) && as.integer(n_lod) > 0)
          paste0("; ", as.integer(n_lod), " missing value(s) imputed") else ""
        sentences <- c(sentences, paste0(
          sub("\\.", paste0(lod_count_str, "."), .lod_sentence_shared(), fixed = TRUE)
        ))
      }

      # 5c. Display fill. Since Pass E this runs AFTER batch correction, so the
      # fabricated cells no longer enter the model that estimates the batch
      # coefficients. The position is recorded in provenance so old and new
      # runs can be told apart.
      .hm_pbn <- suppressWarnings(as.integer(hm$prebatch_n_imputed %||% NA))
      .hm_pbp <- hm$prebatch_imputation_position %||% NA_character_
      if (!is.na(.hm_pbn) && .hm_pbn > 0) {
        sentences <- c(sentences, if (identical(.hm_pbp,
                                               "after_batch_correction")) {
          paste0(
            "Batch correction was computed from the observed values only; ",
            "limma::removeBatchEffect estimates each protein from the samples ",
            "in which it was quantified, so no values were substituted ",
            "beforehand. Afterwards, for display purposes only, the ",
            format(.hm_pbn, big.mark = ","), " remaining missing value",
            if (.hm_pbn != 1) "s were" else " was",
            " replaced with the median of the observed values in the same ",
            "sample (per-sample median). These substituted values do not ",
            "influence the estimated batch coefficients. Proteins with too ",
            "few observations to estimate a batch coefficient are shown ",
            "uncorrected rather than removed.")
        } else {
          paste0(
            "Batch correction requires a complete matrix, so the ",
            format(.hm_pbn, big.mark = ","), " value",
            if (.hm_pbn != 1) "s" else "",
            " still missing at that stage ",
            if (.hm_pbn != 1) "were" else "was",
            " replaced with the median of the observed values in the same sample ",
            "(per-sample median) before limma::removeBatchEffect was applied. ",
            "These substituted values contribute to the estimated batch ",
            "coefficients, which are then applied to every value in the matrix.")
        })
      }

      # 6. Row scaling
      if (.has(hm$row_scaling)) {
        sentences <- c(sentences, glue::glue(
          "Rows (proteins) were {hm$row_scaling}-scaled prior to clustering."
        ))
      }

      # Missing-cell display on the ProDA route.
      if (!is.null(de) &&
          identical(de$analysis_mode %||% de$de_engine %||% "", "proda")) {
        sentences <- c(sentences, paste0(
          "Cells with no measurement are drawn as gaps rather than filled: they 
are not imputed for display, so a blank cell denotes a protein not quantified in 
that sample and is distinct from a protein measured at its own row mean. 
Missing cells are still substituted with the row mean for the purpose of the 
hierarchical clustering distance only."))
      }

      # 7. Clustering
      has_rc <- .has(hm$row_clustering)
      has_cc <- .has(hm$col_clustering)
      if (has_rc && has_cc) {
        if (identical(hm$row_clustering, hm$col_clustering) &&
            !hm$col_clustering %in% c("manual", "grouped")) {
          sentences <- c(sentences, glue::glue(
            "Both rows and columns were clustered using {hm$row_clustering}."
          ))
        } else {
          col_str <- switch(hm$col_clustering,
            manual  = "manually ordered by the user",
            grouped = "grouped by experimental condition",
            hm$col_clustering
          )
          sentences <- c(sentences, glue::glue(
            "Rows were clustered using {hm$row_clustering}; columns were {col_str}."
          ))
        }
      }

      # 8. Conditions
      if (length(hm$condition_names) > 0 && any(nzchar(hm$condition_names %||% character()))) {
        sentences <- c(sentences, glue::glue(
          "Conditions included: {paste(hm$condition_names, collapse = ', ')}."
        ))
      }

      # 9. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = ebn_list_below,
           ebn_ip_sub = hm$ip_sub %||% NULL,
           empty_msg = "Run a Heatmap analysis to generate methods text.")
    }

    # =========================================================================
    # Rank Fold Change methods template
    # =========================================================================

    render_rank_fc_methods <- function(prov) {
      env <- prov$env
      rfc <- prov$rank_fc

      if (is.null(rfc) || !.has(rfc$condition_a)) return(NULL)

      sentences <- character()

      # 1. Overview
      sentences <- c(sentences,
        "Single-shot rank fold-change analysis was performed on Proteome Discoverer abundance data."
      )

      # 2. Conditions & FC method
      if (.has(rfc$condition_a) && .has(rfc$condition_b)) {
        fc_m <- if (.has(rfc$fc_method)) rfc$fc_method else "log2(B/A)"
        sentences <- c(sentences, glue::glue(
          "Fold changes were computed as {fc_m} between Condition B (\"{rfc$condition_b}\") and Condition A (\"{rfc$condition_a}\")."
        ))
      }

      # 3. Protein counts
      has_input   <- length(rfc$n_proteins_input)   > 0 && !is.na(rfc$n_proteins_input)
      has_plotted <- length(rfc$n_proteins_plotted) > 0 && !is.na(rfc$n_proteins_plotted)
      if (has_input && has_plotted) {
        n_dropped <- as.integer(rfc$n_proteins_input) - as.integer(rfc$n_proteins_plotted)
        if (n_dropped > 0) {
          sentences <- c(sentences, glue::glue(
            "{as.integer(rfc$n_proteins_input)} proteins were parsed; {as.integer(n_dropped)} were excluded due to missing or non-finite abundance values in one or both conditions, leaving {as.integer(rfc$n_proteins_plotted)} proteins for analysis."
          ))
        } else {
          sentences <- c(sentences, glue::glue(
            "{as.integer(rfc$n_proteins_plotted)} proteins were included in the analysis."
          ))
        }
      }

      # 4. Imputation
      if (length(rfc$imputation_enabled) > 0 && !is.na(rfc$imputation_enabled)) {
        if (isTRUE(rfc$imputation_enabled) && .has(rfc$imputation_method)) {
          imp_parts <- rfc$imputation_method
          has_ds <- length(rfc$imp_downshift) > 0 && !is.na(rfc$imp_downshift)
          has_w  <- length(rfc$imp_width)     > 0 && !is.na(rfc$imp_width)
          if (has_ds && has_w) {
            imp_parts <- glue::glue(
              "{rfc$imputation_method} (downshift\u00a0=\u00a0{rfc$imp_downshift}, width\u00a0=\u00a0{rfc$imp_width})"
            )
          }
          n_imp <- if (length(rfc$n_imputed) > 0 && !is.na(rfc$n_imputed)) as.integer(rfc$n_imputed) else 0L
          sentences <- c(sentences, glue::glue(
            "Missing values were imputed using {imp_parts}; {n_imp} protein(s) received at least one imputed value."
          ))
        } else {
          sentences <- c(sentences,
            "No imputation was applied; proteins with missing values in either condition were excluded."
          )
        }
      }

      # 5. Ranking
      sentences <- c(sentences,
        "Proteins were ranked by descending absolute log2 fold-change magnitude."
      )

      # 6. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a Rank Fold Change analysis to generate methods text.")
    }

    # =========================================================================
    # Trends Analysis methods template
    # =========================================================================

    render_trends_methods <- function(prov) {
      env <- prov$env
      tr  <- prov$trends

      if (is.null(tr) || length(tr$n_clusters) == 0 || is.na(tr$n_clusters)) return(NULL)

      sentences <- character()

      # 1. Algorithm
      algo <- if (.has(tr$clustering_method)) tr$clustering_method else "fuzzy c-means (e1071::cmeans)"
      sentences <- c(sentences, glue::glue(
        "Single-shot protein abundance trends were analyzed using {algo} as implemented in the R package e1071."
      ))

      # 2. Cluster count & fuzzifier
      has_k <- length(tr$n_clusters) > 0 && !is.na(tr$n_clusters)
      has_m <- length(tr$fuzzifier)  > 0 && !is.na(tr$fuzzifier)
      if (has_k) {
        m_str <- if (has_m) round(tr$fuzzifier, 3) else "auto"
        sentences <- c(sentences, glue::glue(
          "Data were partitioned into {as.integer(tr$n_clusters)} clusters, with a fuzzifier exponent m\u00a0=\u00a0{m_str}."
        ))
      }

      # 3. Protein counts
      has_in <- length(tr$n_proteins_input)     > 0 && !is.na(tr$n_proteins_input)
      has_cl <- length(tr$n_proteins_clustered) > 0 && !is.na(tr$n_proteins_clustered)
      if (has_in && has_cl) {
        n_dropped <- as.integer(tr$n_proteins_input) - as.integer(tr$n_proteins_clustered)
        if (n_dropped > 0) {
          sentences <- c(sentences, glue::glue(
            "{as.integer(tr$n_proteins_input)} proteins were parsed from the Proteome Discoverer export; {as.integer(n_dropped)} were removed due to non-finite abundance values, leaving {as.integer(tr$n_proteins_clustered)} proteins for clustering."
          ))
        } else {
          sentences <- c(sentences, glue::glue(
            "{as.integer(tr$n_proteins_clustered)} proteins were submitted for clustering."
          ))
        }
      }

      # 4. Preprocessing
      sentences <- c(sentences,
        "Abundance values were log2-transformed and row-wise z-score normalized (subtract row mean, divide by row standard deviation) to emphasize relative abundance patterns across conditions."
      )

      # 5. Distance metric
      if (.has(tr$distance_metric)) {
        sentences <- c(sentences, glue::glue(
          "Clustering used {tr$distance_metric} distance."
        ))
      }

      # 6. Conditions
      if (length(tr$condition_names) > 0 && any(nzchar(tr$condition_names))) {
        sentences <- c(sentences, glue::glue(
          "Conditions analyzed: {paste(tr$condition_names, collapse = ', ')}."
        ))
      }

      # 7. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a Trends Analysis to generate methods text.")
    }

    # =========================================================================
    # TMT Pairwise Volcano methods template
    # =========================================================================

    render_tmt_volcano_methods <- function(prov) {
      env <- prov$env
      tv  <- prov$tmt_volcano

      if (is.null(tv) || !.has(tv$analysis_mode)) return(NULL)

      # All current routes produce limma results. The analysis_mode will be
      # "limma_protein", "limma_qf", or "limma" (fallback).
      is_limma <- tv$analysis_mode %in% c("limma_protein", "limma_qf", "limma")
      if (!is_limma) {
        return(list(
          text      = paste0("Methods text is not available for the analysis mode '",
                             tv$analysis_mode,
                             "'. Re-run the TMT analysis with the current pipeline."),
          ebn_list  = NULL,
          empty_msg = "Run a TMT Volcano analysis to generate methods text."
        ))
      }

      sentences <- character()
      design_opt <- as.integer(tv$design_option %||% 1L)
      r_v        <- if (.has(env$r_version)) env$r_version else "?"
      n_drop  <- if (length(tv$n_proteins_dropped_missing) > 0 &&
                     !is.na(tv$n_proteins_dropped_missing))
                   as.integer(tv$n_proteins_dropped_missing) else 0L
      k_val   <- if (length(tv$imputation_k) > 0 && !is.na(tv$imputation_k))
                   as.integer(tv$imputation_k) else "adaptive"
      has_n   <- length(tv$n_proteins)    > 0 && !is.na(tv$n_proteins)
      has_ns  <- length(tv$n_significant) > 0 && !is.na(tv$n_significant)

      # Pooled reference channels per plex (min, max over plexes). Absent on
      # caches made before multiple Norm channels were allowed: one per plex.
      .nbr <- if (length(tv$n_bridge_per_plex) == 2L && !anyNA(tv$n_bridge_per_plex))
                as.integer(tv$n_bridge_per_plex) else c(1L, 1L)
      .multi_bridge <- .nbr[2] > 1L
      .nbr_txt <- if (.nbr[1] == .nbr[2]) as.character(.nbr[1])
                  else paste0(.nbr[1], " to ", .nbr[2])

      # Experiment design description
      design_desc <- switch(as.character(design_opt),
        "1" = "a single-plex TMT experiment",
        "2" = if (.multi_bridge)
                paste0("a multi-plex TMT experiment with internal reference (bridge) channels (",
                       .nbr_txt, " pooled reference channels per plex)")
              else "a multi-plex TMT experiment with an internal reference (bridge) channel",
        # Balance is a claim about the experiment, so only make it when the check
        # actually ran and passed. isFALSE() is deliberate: it is TRUE only when
        # the flag is present and FALSE. On legacy bundles predating the design
        # assertions the flag is absent and neither claim is made.
        "3" = paste0("a multi-plex TMT experiment without a bridge channel",
                     if (isFALSE(tv$design_warn_batch_confounded))
                       ", with a balanced design across runs"
                     else if (isTRUE(tv$design_warn_batch_confounded))
                       paste0(", in which at least one condition was not represented in ",
                              "two or more plexes, so the plex covariate is partially ",
                              "confounded with condition for those conditions")
                     else "",
                     if (length(tv$n_plexes) > 0 && !is.na(tv$n_plexes))
                       paste0(" (", as.integer(tv$n_plexes), " plexes)") else ""),
        "a TMT experiment"
      )

      # Batch handling on the no-bridge route, stated once and used in both
      # route branches so the two cannot drift. Standard practice is to say
      # explicitly which matrix carries the correction: the covariate for
      # inference, removeBatchEffect for display only.
      plex_covariate_str <- paste0(
        "A run (plex) covariate term was included in the limma design matrix to ",
        "account for systematic between-run variation.",
        if (identical(tv$cross_plex_normalization, "model_covariate") ||
            is.null(tv$cross_plex_normalization))
          paste0(" No cross-plex rescaling was applied to the analysis matrix; the ",
                 "covariate estimates and removes the plex offset per protein within ",
                 "the model.")
        else "",
        " For visualization only (principal component analysis, heatmaps and trend ",
        "clusters), plex effects were removed with limma::removeBatchEffect, with ",
        "condition membership preserved as a covariate. Those corrected matrices ",
        "were used for display only and never for statistical testing.")

      # Channel-normalization check (single plex; tmt_normalization_check(), upload
      # wizard step 6), reported with either answer. Absent for multi-plex designs
      # and for caches made before the check existed.
      .tmt_chk_sent <- if (.has(tv$norm_check_within_sd) && !is.na(tv$norm_check_within_sd))
        paste0(" Before summarization, per-channel medians of the ", as.integer(tv$norm_check_n_complete),
               " proteins quantified in every channel were compared between replicate channels of each ",
               "condition (SD ", formatC(as.numeric(tv$norm_check_within_sd), format = "f", digits = 3),
               " log2; values of 0.02 or less are consistent with input already normalized across channels).")
      else ""

      # Channels marked "Exclude" in the design grid (never in the matrix).
      excl_sent <- tmt_excluded_sentence(tv$excluded_channels %||% character(0))

      # Normalization description (depends on route and design)
      norm_desc <- if (identical(as.integer(design_opt), 2L)) {
        # The level at which centering happens differs by route after the
        # protein-level pass was removed from limma_qf: peptide level before
        # aggregation there, protein level on limma_protein. IRS follows in both.
        paste0(if (identical(tv$analysis_mode, "limma_qf"))
                 paste0("Column-wise median centering was applied once, at the peptide ",
                        "level, before aggregation to proteins, with no further centering ",
                        "after aggregation. This was ")
               else
                 paste0("Column-wise median centering was applied at the protein level to ",
                        "align per-channel loading. This was "),
               "followed by Internal Reference Scaling (IRS): for each protein, per-plex ",
               "bridge-channel means were used to estimate and correct multiplicative scaling ",
               "differences between plexes, aligning intensities to a global bridge reference.",
               if (.multi_bridge)
                 paste0(" Where a plex carried more than one bridge channel, the per-plex ",
                        "reference was the mean of its bridge channels (log2), computed from ",
                        "those observed for the protein.")
               else "")
      } else if (isFALSE(tv$normalization_applied)) {
        # Single plex, input reported as already normalized across channels (wizard step 6)
        paste0("No column-wise median centering was applied: in the upload wizard, the user reported ",
               "the reporter ion abundances as already normalized across channels before upload, and ",
               "they were used as they were (log2).", .tmt_chk_sent,
               if (identical(tv$norm_check_verdict, "not_normalized") ||
                   identical(tv$norm_check_verdict, "borderline"))
                 " The user kept this setting although the check did not clearly indicate normalized input."
               else "")
      } else {
        paste0(if (identical(tv$analysis_mode, "limma_qf"))
                 paste0("Column-wise median centering was applied once, at the peptide level, ",
                        "before aggregation to proteins. No further centering was applied after ",
                        "aggregation.")
               else
                 "Column-wise median centering was applied to align global signal levels across TMT channels.",
               .tmt_chk_sent)
      }

      # ── Shared sentences (identical in both routes) ──────────────────────────
      .kv <- if (length(tv$knn_k) > 0 && !is.na(tv$knn_k)) as.integer(tv$knn_k) else 10L

      # Eligibility cap applied within a plex (T4). Reported numerically so the
      # rule that was actually applied is on the record.
      .pmof <- if (length(tv$plex_min_obs_frac) > 0 && !is.na(tv$plex_min_obs_frac))
                 as.numeric(tv$plex_min_obs_frac) else 0.5
      .n_wp <- if (length(tv$n_values_imputed_within_plex) > 0 &&
                   !is.na(tv$n_values_imputed_within_plex))
                 as.integer(tv$n_values_imputed_within_plex) else NA_integer_

      # A channel over the per-channel missingness limit makes impute.knn error;
      # the caller catches it and NOTHING is imputed. Without this the methods
      # text would still describe k-NN imputation as performed. Appended to both
      # imputation branches.
      .n_cmax <- if (length(tv$n_channels_over_colmax) > 0 &&
                     !is.na(tv$n_channels_over_colmax))
                   as.integer(tv$n_channels_over_colmax) else 0L
      # Only reachable when values WERE filled - the fills-zero case is handled by
      # its own branch below. Stated factually so it cannot contradict the
      # imputed-value count reported in the same sentence.
      colmax_str <- if (.n_cmax > 0L)
        paste0(" ", .n_cmax, " channel", if (.n_cmax != 1L) "s" else "",
               " exceeded the per-channel missingness limit (colmax\u2009=\u20090.8), which ",
               "aborts the k-nearest-neighbour computation for any matrix affected; ",
               "the count above is the number of values actually filled.")
      else ""
      # Within-plex route: the event that leaves values unimputed is a per-plex
      # k-NN failure, counted by n_plex_knn_failed. colmax is deliberately NOT
      # mentioned here - that pre-check runs only on the global route, so on this
      # route n_channels_over_colmax is always 0 and would never fire.
      .n_pfail <- if (length(tv$n_plex_knn_failed) > 0 && !is.na(tv$n_plex_knn_failed))
                    as.integer(tv$n_plex_knn_failed) else 0L
      # Only reachable when values WERE filled in at least one plex; the
      # nothing-filled case is handled by its own branch below.
      plexfail_str <- if (.n_pfail > 0L)
        paste0(" The k-nearest-neighbour computation failed for ", .n_pfail, " plex",
               if (.n_pfail != 1L) "es" else "",
               ", so no values were imputed for the affected plex",
               if (.n_pfail != 1L) "es" else "",
               " and those missing values were retained for model fitting; the count ",
               "above is the number of values actually filled.")
      else ""
      # Effective fill count for the route actually taken. When a channel blew the
      # colmax limit AND nothing was filled, the text must not simultaneously
      # describe k-NN imputation as performed - it was attempted and aborted.
      .n_imp_eff <- if (isTRUE(tv$within_plex_knn)) {
        if (length(tv$n_values_imputed_within_plex) > 0 &&
            !is.na(tv$n_values_imputed_within_plex))
          as.integer(tv$n_values_imputed_within_plex) else NA_integer_
      } else {
        if (length(tv$n_values_imputed_global) > 0 &&
            !is.na(tv$n_values_imputed_global))
          as.integer(tv$n_values_imputed_global) else NA_integer_
      }

      # Each route has its own failure mechanism: per-plex k-NN failure within
      # plexes, a channel over the colmax limit globally.
      .route_failed <- if (isTRUE(tv$within_plex_knn)) .n_pfail > 0L else .n_cmax > 0L
      impute_sent <- if (.route_failed && !is.na(.n_imp_eff) && .n_imp_eff == 0L) {
        if (isTRUE(tv$within_plex_knn))
          paste0("Imputation by k-nearest neighbours was attempted but not performed: ",
                 "the computation failed for ", .n_pfail, " plex",
                 if (.n_pfail != 1L) "es" else "",
                 ". No values were imputed and all missing values were retained for ",
                 "model fitting, so every protein was fitted on its observed values only.")
        else
          paste0("Imputation by k-nearest neighbours was attempted but not performed: ",
                 .n_cmax, " channel", if (.n_cmax != 1L) "s" else "",
                 " exceeded the per-channel missingness limit (colmax\u2009=\u20090.8), which ",
                 "aborts the k-nearest-neighbour computation for the whole matrix. No ",
                 "values were imputed and all missing values were retained for model ",
                 "fitting, so every protein was fitted on its observed values only.")
      } else if (isTRUE(tv$within_plex_knn)) {
        paste0("Missing values were imputed by k-nearest neighbours applied independently ",
               "within each plex (k\u2009=\u2009", .kv, ", impute::impute.knn), so that neighbours ",
               "were drawn from proteins measured in the same run",
               if (!is.na(.n_wp)) paste0(" (", .n_wp, " values imputed)") else "",
               ". Imputation was applied only after the coverage filter, never before it. ",
               "A protein was eligible for imputation within a plex only if at least ",
               round(.pmof * 100), "% of that plex's channels were observed for it; this ",
               "threshold matches the ", round(.pmof * 100), "% coverage rule used by the ",
               "protein-level filter, so a single coverage rule governs the pipeline. ",
               "Proteins below it retained their missing values, which limma fits on the ",
               "observations available. The k-nearest-neighbour search was computed over the ",
               "full plex submatrix so that donors were drawn from every protein in the plex, ",
               "and imputed values were written back only to eligible proteins. No sample-mean ",
               "fallback was permitted (rowmax\u2009=\u20091), the clustering approximation was ",
               "disabled (maxp\u2009=\u2009Inf), the per-channel missingness limit was left at its ",
               "default (colmax\u2009=\u20090.8), and a fixed random seed (42) was used. Reporter ion ",
               "intensities recorded as zero or negative were treated as confidently absent ",
               "rather than missing and were not filled. No imputation was performed across plexes.",
               plexfail_str)
      } else if (isTRUE(tv$cross_plex_impute_disabled)) {
        paste0("No imputation was performed. Missing values were retained and limma fitted ",
               "each protein on its available observations; filling values for a plex in which ",
               "a protein was never identified would fabricate an entire block of measurements.")
      } else if (identical(tv$impute_pkg_available, FALSE)) {
        paste0("No imputation was performed: the 'impute' package was unavailable at run ",
               "time. Every protein was fitted on its observed values only, so proteins with ",
               "missing values contributed fewer observations to the model and results differ ",
               "from a run with imputation enabled.")
      } else {
        .n_imp  <- if (length(tv$n_values_imputed_global) > 0 &&
                       !is.na(tv$n_values_imputed_global))
                     as.integer(tv$n_values_imputed_global) else NA_integer_
        .n_zero <- if (length(tv$n_zero_origin_held) > 0 &&
                       !is.na(tv$n_zero_origin_held))
                     as.integer(tv$n_zero_origin_held) else 0L
        paste0("Missing values were imputed by k-nearest neighbours (k\u2009=\u2009", .kv,
               ", impute::impute.knn) in a single pass over the entire experiment, carried ",
               "out before any pairwise comparison was selected, so that imputed values do ",
               "not depend on the contrast being tested. A protein was eligible only if it ",
               "had at most one missing value in every condition in the design, not only in ",
               "the two conditions being compared",
               if (!is.na(.n_imp)) paste0(" (", .n_imp, " values imputed)") else "",
               ". Imputation was applied only after the coverage filter, never before it. ",
               "Neighbours were drawn from all annotated channels: the k-nearest-neighbour ",
               "search was computed over the full filtered protein matrix, so that donors ",
               "were not restricted to near-complete proteins, and imputed values were ",
               "written back only to eligible proteins. No sample-mean fallback ",
               "was permitted (rowmax\u2009=\u20091), the clustering approximation was disabled ",
               "(maxp\u2009=\u2009Inf), the per-channel missingness limit was left at its ",
               "default (colmax\u2009=\u20090.8), and a fixed random seed (42) was used. ",
               "Reporter ion intensities ",
               "recorded as zero or negative were treated as confidently absent rather than ",
               "missing and were not filled",
               if (.n_zero > 0L) paste0(" (", .n_zero, " such values retained as missing)")
               else "",
               ". All remaining missing values were retained for model fitting.",
               colmax_str)
      }

      .n_gf <- if (length(tv$n_proteins_global_filter) > 0 &&
                   !is.na(tv$n_proteins_global_filter))
                 as.integer(tv$n_proteins_global_filter) else NA_integer_
      .n_oa <- if (length(tv$n_proteins_one_arm_excluded) > 0 &&
                   !is.na(tv$n_proteins_one_arm_excluded))
                 as.integer(tv$n_proteins_one_arm_excluded) else NA_integer_

      # Plex-coverage filter at summarization (protein route, multi-plex) and the
      # table of proteins not tested. Both absent on bundles made before they were
      # recorded, and then omitted.
      .n_pc <- if (.has(tv$n_proteins_plex_coverage_dropped)) as.integer(tv$n_proteins_plex_coverage_dropped) else 0L
      .n_nt <- if (.has(tv$n_proteins_not_tested_listed)) as.integer(tv$n_proteins_not_tested_listed) else 0L
      .plexcov_sent <- if (.n_pc > 0L && .has(tv$plex_coverage_min_plexes) && .has(tv$plex_coverage_n_plexes))
        paste0("At summarization, proteins quantified in fewer than ", as.integer(tv$plex_coverage_min_plexes),
               " of the ", as.integer(tv$plex_coverage_n_plexes), " plexes were removed (", .n_pc, " proteins). ")
        else ""
      .listed_sent <- if (.n_nt > 0L)
        paste0(" The ", .n_nt, " proteins removed by the coverage filters were not tested; they are listed in a ",
               "separate table with their observed channels per condition", if (.n_pc > 0L || isTRUE(design_opt > 1L))
               " and per plex" else "", " and the reason for removal.")
        else ""

      filter_sent <- if (isTRUE(tv$global_filter_applied)) {
        paste0(.plexcov_sent, "Missingness filtering was applied in two stages. First, before any ",
               "comparison was selected, proteins were required to have at least 50% valid ",
               "measurements in at least one condition",
               if (!is.na(.n_gf)) paste0(" (", .n_gf, " proteins retained)") else "",
               "; this fixes a single protein set for every comparison, so that the ",
               "empirical Bayes prior, the imputed values and the fitted protein set do ",
               "not depend on which contrast is tested. No further filtering was applied ",
               "for the comparison reported here: every protein in that set entered the ",
               "model, and the requested contrast was extracted from it. Proteins with at ",
               "least 50% valid measurements in only one of the two compared conditions ",
               "were flagged as low coverage but tested normally",
               if (!is.na(.n_oa)) paste0(" (", .n_oa, " proteins)") else "",
               "; their limited evidence is reflected in a larger standard error and in ",
               "per-protein residual degrees of freedom rather than in a separate ",
               "exclusion rule. Proteins with no observations at all in one of the ",
               "compared conditions have no estimable contrast, return a missing ",
               "p-value and are therefore absent from the false-discovery-rate correction",
               if (length(tv$n_proteins_unestimable) > 0 &&
                   !is.na(tv$n_proteins_unestimable))
                 paste0(" (", as.integer(tv$n_proteins_unestimable), " proteins)")
               else "",
               ".", .listed_sent)
      #
      # DO NOT DELETE THIS BRANCH.
      # It describes the old per-pair-only missingness filter. Since the
      # contrast-independent global filter was introduced, global_filter_applied is
      # TRUE on every route, so every CURRENT run takes the branch above and this
      # one is unreachable for new analyses. It is retained solely so that prerun
      # bundles generated BEFORE that change still render methods text that matches
      # how they were actually processed.
      } else {
        paste0("Proteins with fewer than 50% valid measurements in either of the two ",
               "compared conditions were excluded from testing (", n_drop,
               " proteins removed); coverage in the remaining conditions was not required, ",
               "as those channels contribute only to variance estimation.")
      }

      model_sent <- paste0(
        "Statistical testing was performed using limma. A single linear model was fitted ",
        "across all ", if (length(tv$n_conditions_model) > 0 && !is.na(tv$n_conditions_model))
                          paste0(as.integer(tv$n_conditions_model), " ") else "",
        "conditions in the experiment and the requested contrast extracted via contrasts.fit, ",
        "rather than refitting on each pairwise subset. This yields ",
        if (isTRUE(tv$global_filter_applied))
          "one consistent empirical Bayes prior across comparisons and substantially more residual degrees of freedom"
        else
          "substantially more residual degrees of freedom",
        if (length(tv$df_residual) > 0 && !is.na(tv$df_residual))
          paste0(" (", as.integer(tv$df_residual), " residual df per protein)") else "",
        ", and assumes a common residual variance across conditions. ",
        if (isTRUE(tv$array_weights_applied))
          paste0("That assumption was not left unchecked: sample weights were estimated ",
                 "per condition with limma::arrayWeights (var.group = condition) and ",
                 "applied to the fit, so a condition with higher residual variance ",
                 "contributes proportionally less to the pooled estimate",
                 if (length(tv$array_weights) > 0 && all(!is.na(tv$array_weights)))
                   paste0(" (weights: ",
                          paste(round(as.numeric(tv$array_weights), 2), collapse = ", "),
                          ")") else "",
                 if (length(tv$residual_var_ratio) > 0 && !is.na(tv$residual_var_ratio))
                   paste0(". The ratio of the largest to the smallest per-condition ",
                          "median residual standard deviation was ",
                          signif(as.numeric(tv$residual_var_ratio), 3), ". ")
                 else ". ")
        else "",
        "Variance was moderated ",
        "by empirical Bayes with a mean-intensity trend and robust estimation (eBayes: ",
        "trend\u2009=\u2009TRUE, robust\u2009=\u2009TRUE), followed by Benjamini\u2013Hochberg FDR correction ",
        "within this comparison.",
        # Design caveats recorded by the pre-model assertions. Batch confounding is
        # deliberately NOT repeated here - it is already stated in design_desc.
        # Absent flags (legacy bundles) emit nothing.
        if (isTRUE(tv$design_warn_low_replication) ||
            isTRUE(tv$design_warn_low_residual_df)) {
          paste0(" Note on the design: ",
                 if (isTRUE(tv$design_warn_low_replication))
                   paste0("at least one condition had fewer than two replicates, ",
                          "contributing a coefficient to the model but no residual ",
                          "degrees of freedom, so its variance estimate rests entirely ",
                          "on the empirical Bayes prior")
                 else "",
                 if (isTRUE(tv$design_warn_low_replication) &&
                     isTRUE(tv$design_warn_low_residual_df)) "; " else "",
                 if (isTRUE(tv$design_warn_low_residual_df))
                   paste0("the number of model parameters approached the number of ",
                          "samples, leaving few residual degrees of freedom")
                 else "",
                 ".")
        } else "")

      if (identical(tv$analysis_mode, "limma_protein")) {
        sentences <- c(sentences,
          glue::glue("Differential expression analysis was performed on protein-level TMT reporter ion intensities from {design_desc}."),
          excl_sent,
          "Raw abundances were log2-transformed (zero and non-positive values were replaced with NA prior to transformation).",
          norm_desc,
          filter_sent,
          impute_sent,
          if (identical(as.integer(design_opt), 3L))
            plex_covariate_str
          else character(0),
          model_sent
        )
      } else {
        # limma_qf (peptide/PSM aggregation route)
        data_lvl <- if (.has(tv$data_level) && tv$data_level %in% c("protein", "peptide", "psm", "limma_protein", "limma_qf")) {
          switch(tv$data_level,
            "limma_protein" = "protein",
            "limma_qf"      = "peptide",
            tv$data_level)
        } else "peptide"
        data_lvl_cap <- paste0(toupper(substring(data_lvl, 1, 1)), substring(data_lvl, 2))
        interf_str <- if (isTRUE(tv$interference_filtered)) {
          paste0(" Precursor rows with co-isolation interference exceeding ",
                 if (!is.na(tv$interference_threshold)) tv$interference_threshold else 50,
                 "% were excluded prior to aggregation; rows with no reported interference ",
                 "value were retained",
                 if (length(tv$n_interf_kept_na) > 0 && !is.na(tv$n_interf_kept_na) &&
                     length(tv$n_interf_kept_pass) > 0 && !is.na(tv$n_interf_kept_pass))
                   paste0(" (", format(as.integer(tv$n_interf_kept_pass), big.mark = ","),
                          " rows retained on passing the threshold, ",
                          format(as.integer(tv$n_interf_kept_na), big.mark = ","),
                          " retained because no interference value was reported)")
                 else "",
                 ".")
        } else {
          " NOTE: Isolation Interference column was not found \u2014 co-isolation interference was not filtered."
        }
        # Peptide-level confidence is inherited, not re-applied here. Stating it
        # explicitly stops a reader assuming this pipeline performed its own
        # peptide-level FDR control.
        pep_fdr_str <- if (identical(tv$peptide_fdr_source, "qvalue_filtered_here")) {
          " Peptide-level identifications were additionally filtered on q-value during summarization."
        } else {
          paste0(" Peptide-level false discovery rate control was inherited from the ",
                 "search engine export (e.g. Proteome Discoverer) and was not ",
                 "re-applied downstream; no additional peptide- or PSM-level q-value ",
                 "filter was applied by this pipeline.")
        }
        n_pep    <- if (length(tv$n_peptides_input) > 0 && !is.na(tv$n_peptides_input))
                      as.integer(tv$n_peptides_input) else NA
        pep_str  <- if (!is.na(n_pep)) glue::glue(" ({n_pep} {data_lvl} rows input)") else ""
        interf_str <- paste0(interf_str, pep_fdr_str)
        .n_sing <- if (length(tv$n_agg_singleton) > 0 && !is.na(tv$n_agg_singleton))
                     as.integer(tv$n_agg_singleton) else 0L
        .n_cmea <- if (length(tv$n_agg_fallback_colmeans) > 0 &&
                       !is.na(tv$n_agg_fallback_colmeans))
                     as.integer(tv$n_agg_fallback_colmeans) else 0L
        agg_fallback_str <- if (.n_sing > 0L || .n_cmea > 0L) {
          paste0("Robust regression is undefined for a single observation and can fail to ",
                 "converge, so a fallback was used where required: ",
                 if (.n_sing > 0L)
                   paste0("for ", format(.n_sing, big.mark = ","),
                          " groups represented by a single feature, that feature's values ",
                          "were carried through directly")
                 else "",
                 if (.n_sing > 0L && .n_cmea > 0L) "; " else "",
                 if (.n_cmea > 0L)
                   paste0("for ", format(.n_cmea, big.mark = ","),
                          " groups the robust fit did not converge and the column mean ",
                          "of the contributing features was used instead")
                 else "",
                 ".")
        } else ""
        # Shared-peptide handling. The filter precedes aggregation, so this sits
        # between the data-level sentence and the summarization sentence rather
        # than being concatenated onto interf_str, which already absorbs
        # pep_fdr_str. Silent only when the setting is present but unknown (NA).
        .sp_mode <- if (length(tv$shared_peptide_handling) > 0)
                      tv$shared_peptide_handling else "keep_all"
        .sp_tot  <- if (length(tv$n_peptides_total) > 0 && !is.na(tv$n_peptides_total))
                      as.integer(tv$n_peptides_total) else NA_integer_
        .sp_rm   <- if (length(tv$n_peptides_removed_shared) > 0 &&
                        !is.na(tv$n_peptides_removed_shared))
                      as.integer(tv$n_peptides_removed_shared) else NA_integer_
        .sp_how  <- if (length(tv$shared_detection_method) > 0 &&
                        !is.na(tv$shared_detection_method))
                      as.character(tv$shared_detection_method) else "master_accession_count"
        .sp_how_txt <- switch(.sp_how,
          quan_info         = "the search engine's own unique/razor annotation",
          n_proteins_column = "the search engine's per-peptide protein-count column",
          master_accession_count = "the list of master protein accessions reported for each peptide",
          mixed             = "a mixture of the above, differing between input files",
          "the list of master protein accessions reported for each peptide")
        shared_pep_str <- if (is.na(.sp_mode)) {
          ""
        } else if (identical(.sp_mode, "unique_only")) {
          paste0("Peptides matching more than one master protein were removed before ",
                 "protein-level summarization, following the MSstatsTMT default ",
                 "(useUniquePeptide = TRUE)",
                 if (!is.na(.sp_rm) && !is.na(.sp_tot))
                   paste0("; ", format(.sp_rm, big.mark = ","), " of ",
                          format(.sp_tot, big.mark = ","), " peptide row",
                          if (.sp_tot != 1) "s" else "", " ",
                          if (.sp_rm != 1) "were" else "was", " removed")
                 else "",
                 ". Shared peptides were identified from ", .sp_how_txt, ".")
        } else if (identical(.sp_mode, "razor")) {
          .rz_a <- if (length(tv$n_peptides_razor_assigned) > 0 &&
                       !is.na(tv$n_peptides_razor_assigned))
                     as.integer(tv$n_peptides_razor_assigned) else NA_integer_
          .rz_d <- if (length(tv$n_peptides_dropped_no_unique_evidence) > 0 &&
                       !is.na(tv$n_peptides_dropped_no_unique_evidence))
                     as.integer(tv$n_peptides_dropped_no_unique_evidence) else NA_integer_
          paste0("Peptides matching more than one master protein were assigned to a ",
                 "single protein by parsimony: each was attributed to the candidate ",
                 "protein supported by the greatest number of unique peptide ",
                 "sequences across the dataset, with ties resolved alphabetically",
                 if (!is.na(.rz_a))
                   paste0(" (", format(.rz_a, big.mark = ","), " peptide row",
                          if (.rz_a != 1) "s" else "", " ",
                          if (.rz_a != 1) "were" else "was", " reassigned)")
                 else "",
                 ". Peptides whose candidate proteins had no unique-peptide support ",
                 "at all could not be assigned and were removed",
                 if (!is.na(.rz_d))
                   paste0(" (", format(.rz_d, big.mark = ","), " peptide row",
                          if (.rz_d != 1) "s" else "", ")")
                 else "",
                 ". This parsimony assignment was derived from the exported peptide ",
                 "table by this pipeline; a peptide-level export carries no razor ",
                 "annotation, so the assignment was not made by the search engine. ",
                 "Shared peptides were identified from ", .sp_how_txt, ".")
        } else {
          paste0("Peptides matching more than one master protein were retained and ",
                 "attributed to the first accession listed in the export",
                 if (!is.na(.sp_tot) && !is.na(tv$n_peptides_shared %||% NA) &&
                     length(tv$n_peptides_shared) > 0 && !is.na(tv$n_peptides_shared))
                   paste0(" (", format(as.integer(tv$n_peptides_shared), big.mark = ","),
                          " of ", format(.sp_tot, big.mark = ","), " peptide row",
                          if (.sp_tot != 1) "s" else "", ")")
                 else "",
                 ". That attribution is positional and is not evidence-based: it is ",
                 "not a razor assignment, and the shared signal is credited to one ",
                 "candidate protein rather than apportioned among them.")
        }
        sentences <- c(sentences,
          glue::glue("Differential expression analysis was performed on TMT reporter ion intensities from {design_desc}."),
          excl_sent,
          glue::glue("{data_lvl_cap}-level data were used as input, taken from raw reporter ion ",
                     "abundance columns; Proteome Discoverer's pre-scaled abundance columns were ",
                     "not used, as they are already normalised.{interf_str}"),
          shared_pep_str,
          # The centering clause lives in norm_desc, which is emitted on the very
          # next line and is the more complete statement. Stating it here too gave
          # the reader the same fact twice in consecutive sentences.
          glue::glue("{data_lvl_cap}-to-protein summarization was performed using QFeatures with robust regression summarization (MsCoreUtils::robustSummary, MASS::rlm), ",
                     "producing protein-level log2 abundances{pep_str}."),
          # Up to three estimators can coexist in one output matrix. Claiming
          # robustSummary for every protein is an overclaim, so say how many took
          # each fallback. Nothing is emitted when both counters are zero or absent.
          agg_fallback_str,
          norm_desc,
          filter_sent,
          impute_sent,
          if (identical(as.integer(design_opt), 3L))
            plex_covariate_str
          else character(0),
          model_sent
        )
      }

      # Remove any empty strings from optional sentences
      sentences <- sentences[nchar(sentences) > 0]

      if (has_n && has_ns) {
        # Single-peptide quantification is one of the first things a proteomics
        # reviewer asks about. These proteins are tested and FDR-corrected like
        # any other and are flagged in the results table; reporting the count is
        # a sensitivity statement, not an exclusion. Silent when the keys are absent.
        .sp_n <- if (length(tv$n_tested_single_peptide) > 0 &&
                     !is.na(tv$n_tested_single_peptide))
                   as.integer(tv$n_tested_single_peptide) else NA_integer_
        .sp_s <- if (length(tv$n_significant_single_peptide) > 0 &&
                     !is.na(tv$n_significant_single_peptide))
                   as.integer(tv$n_significant_single_peptide) else NA_integer_
        .sp_str <- if (!is.na(.sp_n) && !is.na(.sp_s) && .sp_n > 0L)
          glue::glue(" Of the proteins analysed, {format(.sp_n, big.mark = \",\")} ",
                     "were quantified from a single peptide, of which ",
                     "{format(.sp_s, big.mark = \",\")} were significant; these were ",
                     "retained in the model and in the FDR correction and are flagged ",
                     "in the results table.")
        else ""
        sentences <- c(sentences, glue::glue(
          "{as.integer(tv$n_proteins)} proteins were included in the final analysis; ",
          "{as.integer(tv$n_significant)} were statistically significant after Benjamini-Hochberg FDR correction (adj. p-value < 0.05).{.sp_str}"
        ))
      }
      if (.has(env$r_version)) {
        sentences <- c(sentences,
          glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a TMT Volcano analysis to generate methods text.")
    }

    # =========================================================================
    # TMT Trend Clusters methods template
    # =========================================================================

    render_tmt_trends_methods <- function(prov) {
      env <- prov$env
      # The unified trend-clusters module writes rv$provenance$trend_clusters
      # with source = "tmt"; prov$tmt_trends is the older per-source ledger
      # entry, which no live route populates. Accept whichever is present so
      # both wirings produce text. Mirrors render_dia_trends_methods.
      tt  <- prov$tmt_trends
      if (!.has(tt$n_clusters) &&
          identical(prov$trend_clusters$source %||% "", "tmt"))
        tt <- prov$trend_clusters

      if (is.null(tt) || length(tt$n_clusters) == 0 || is.na(tt$n_clusters)) return(NULL)

      sentences <- character()

      # 1. Algorithm
      algo <- if (.has(tt$clustering_method)) tt$clustering_method else "fuzzy c-means (e1071::cmeans)"
      sentences <- c(sentences, glue::glue(
        "TMT protein abundance trends were analyzed using {algo} as implemented in the R package e1071."
      ))

      # 2. Summarization pipeline
      sentences <- c(sentences,
        "Protein-level abundances were obtained by limma-based summarization (QFeatures robustSummary with peptide-level median centering, or protein-level summarization with protein-level median centering). The summarized matrix retains its missing values; any imputation described here applies to this visualization only and is independent of the imputation used for differential expression testing."
      )

      # 3. Missing data handling
      if (.has(tt$missing_mode)) {
        if (identical(tt$missing_mode, "lod")) {
          has_lod <- length(tt$lod_value) > 0 && !is.na(tt$lod_value)
          # The floor became per-sample on the TMT routes; lod_scope records
          # which construction a given run used. Runs predating that change (an
          # imported prerun bundle) carry no lod_scope and did use a single
          # matrix-wide value, so their original wording is kept for them.
          .lod_per_sample <- identical(tt$lod_scope %||% "", "per_sample_1st_percentile")
          lod_str <- if (!has_lod) "" else if (.lod_per_sample) {
            glue::glue(", computed separately for each sample as the 1st percentile of that sample's observed values (minimum-detection substitution, equivalent to imputeLCMD::impute.MinDet with q = 0.01); the median per-sample floor was {round(tt$lod_value, 4)} log2 abundance units")
          } else {
            glue::glue(" (LOD = {round(tt$lod_value, 4)} log2 abundance units, 1st percentile of the protein-level matrix)")
          }
          n_imp   <- if (length(tt$n_imputed) > 0 && !is.na(tt$n_imputed)) as.integer(tt$n_imputed) else 0L
          sentences <- c(sentences, glue::glue(
            "Missing protein abundances were imputed using a limit-of-detection (LOD) strategy{lod_str}; {n_imp} missing values across all conditions were replaced."
          ))
        } else {
          n_om <- if (length(tt$n_omitted) > 0 && !is.na(tt$n_omitted)) as.integer(tt$n_omitted) else 0L
          sentences <- c(sentences, glue::glue(
            "Proteins with any missing abundance values across the selected conditions were excluded from clustering ({n_om} proteins removed)."
          ))
        }
      }

      # 4. Batch correction
      if (length(tt$batch_corrected) > 0 && isTRUE(tt$batch_corrected)) {
        sentences <- c(sentences,
          "Plex-level batch effects were removed prior to clustering using limma::removeBatchEffect, with condition membership preserved as a covariate."
        )
        # Proteins seen in one plex only have no between-plex contrast, so they
        # pass through the correction unadjusted. Disclosed rather than left to
        # imply that every protein was corrected.
        .tt_sp <- if (length(tt$n_proteins_single_plex) > 0 &&
                      !is.na(tt$n_proteins_single_plex))
                    as.integer(tt$n_proteins_single_plex) else 0L
        if (.tt_sp > 0) {
          sentences <- c(sentences, glue::glue(
            "{.tt_sp} proteins were observed in a single plex; these carry no plex adjustment, as the batch term is not identifiable from a single plex."
          ))
        }
      }

      # 5. Protein counts
      has_in <- length(tt$n_proteins_input)     > 0 && !is.na(tt$n_proteins_input)
      has_cl <- length(tt$n_proteins_clustered) > 0 && !is.na(tt$n_proteins_clustered)
      if (has_in && has_cl) {
        sentences <- c(sentences, glue::glue(
          "{as.integer(tt$n_proteins_clustered)} proteins were submitted for clustering (of {as.integer(tt$n_proteins_input)} input proteins after missing-data filtering and zero-variance removal)."
        ))
      }

      # 6. Z-score transform
      sentences <- c(sentences,
        "Abundance values were row-wise z-score normalized (per protein: subtract mean, divide by standard deviation) to emphasize relative abundance patterns across conditions."
      )

      # 7. Cluster count & fuzzifier
      has_k <- length(tt$n_clusters) > 0 && !is.na(tt$n_clusters)
      has_m <- length(tt$fuzzifier)  > 0 && !is.na(tt$fuzzifier)
      if (has_k) {
        m_str <- if (has_m) round(tt$fuzzifier, 3) else "auto"
        sentences <- c(sentences, glue::glue(
          "Data were partitioned into {as.integer(tt$n_clusters)} fuzzy clusters with fuzzifier exponent m\u00a0=\u00a0{m_str} and euclidean distance."
        ))
      }

      # 8. Conditions
      if (length(tt$condition_names) > 0 && any(nzchar(tt$condition_names))) {
        sentences <- c(sentences, glue::glue(
          "Conditions analyzed: {paste(tt$condition_names, collapse = ', ')}."
        ))
      }

      # 9. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a TMT Trend Clusters analysis to generate methods text.")
    }

    # =========================================================================
    # TMT PCA methods template
    # =========================================================================

    render_tmt_pca_methods <- function(prov) {
      env <- prov$env
      tp  <- prov$tmt_pca

      if (is.null(tp) || length(tp$n_samples) == 0 || is.na(tp$n_samples)) return(NULL)

      sentences <- character()

      # 1. Introduction
      sentences <- c(sentences,
        "TMT principal component analysis (PCA) was performed to assess sample-level variation and clustering."
      )

      # 2. Summarization pipeline
      sentences <- c(sentences,
        "Protein-level abundances were obtained by limma-based summarization (QFeatures robustSummary with peptide-level median centering, or protein-level summarization with protein-level median centering). The summarized matrix retains its missing values; any imputation described here applies to this visualization only and is independent of the imputation used for differential expression testing."
      )

      # 3. Design option
      if (length(tp$design_option) > 0 && !is.na(tp$design_option)) {
        design_str <- switch(as.character(as.integer(tp$design_option)),
          "1" = "single-plex design",
          "2" = "multi-plex design with a bridge/reference channel",
          "3" = "balanced multi-plex design",
          paste0("design option ", as.integer(tp$design_option))
        )
        sentences <- c(sentences, glue::glue(
          "The experiment used a {design_str}."
        ))
      }
      # 4. Batch correction
      if (length(tp$batch_corrected) > 0 && isTRUE(tp$batch_corrected)) {
        sentences <- c(sentences,
          "Inter-plex batch effects were removed prior to PCA using limma::removeBatchEffect, with condition membership as a design covariate."
        )
      }

      # 5. PCA settings
      sentences <- c(sentences,
        "PCA was computed using the base R prcomp() function with centering enabled (center = TRUE) and without variance scaling (scale. = FALSE)."
      )

      # 6. Matrix dimensions
      has_ns <- length(tp$n_samples)  > 0 && !is.na(tp$n_samples)
      has_np <- length(tp$n_proteins) > 0 && !is.na(tp$n_proteins)
      if (has_ns && has_np) {
        sentences <- c(sentences, glue::glue(
          "The input matrix comprised {as.integer(tp$n_samples)} samples and {as.integer(tp$n_proteins)} proteins."
        ))
      }

      # 7. Variance explained
      has_pc1 <- length(tp$pc1_var) > 0 && !is.na(tp$pc1_var)
      has_pc2 <- length(tp$pc2_var) > 0 && !is.na(tp$pc2_var)
      if (has_pc1 && has_pc2) {
        sentences <- c(sentences, glue::glue(
          "PC1 and PC2 explained {round(tp$pc1_var, 1)}% and {round(tp$pc2_var, 1)}% of the total variance, respectively."
        ))
      }

      # 8. Conditions
      if (length(tp$condition_names) > 0 && any(nzchar(tp$condition_names))) {
        sentences <- c(sentences, glue::glue(
          "Conditions analyzed: {paste(tp$condition_names, collapse = ', ')}."
        ))
      }

      # 9. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a TMT PCA analysis to generate methods text.")
    }

    # =========================================================================
    # TMT Heatmap methods template
    # =========================================================================

    render_tmt_heatmap_methods <- function(prov) {
      env <- prov$env
      # See render_tmt_trends_methods: the unified heatmap module writes
      # rv$provenance$heatmap with source = "tmt"; prov$tmt_heatmap is the older
      # per-source ledger entry, which no live route populates.
      th  <- prov$tmt_heatmap
      if (!.has(th$missing_mode) &&
          identical(prov$heatmap$source %||% "", "tmt"))
        th <- prov$heatmap

      if (is.null(th) || !.has(th$missing_mode)) return(NULL)

      sentences <- character()

      # 1. Introduction
      sentences <- c(sentences,
        "TMT protein abundance data were visualized as a hierarchically clustered heatmap."
      )

      # 2. Summarization pipeline
      sentences <- c(sentences,
        "Protein-level abundances were obtained by limma-based summarization (QFeatures robustSummary with peptide-level median centering, or protein-level summarization with protein-level median centering). The summarized matrix retains its missing values; any imputation described here applies to this visualization only and is independent of the imputation used for differential expression testing."
      )

      # 2b. Batch correction
      if (isTRUE(th$batch_corrected)) {
        sentences <- c(sentences,
          "Inter-plex batch effects were removed using limma::removeBatchEffect, with plex (Run) as the batch factor and biological condition retained in the design matrix."
        )
        # See render_tmt_trends_methods.
        .th_sp <- if (length(th$n_proteins_single_plex) > 0 &&
                      !is.na(th$n_proteins_single_plex))
                    as.integer(th$n_proteins_single_plex) else 0L
        if (.th_sp > 0) {
          sentences <- c(sentences, glue::glue(
            "{.th_sp} proteins were observed in a single plex; these carry no plex adjustment, as the batch term is not identifiable from a single plex."
          ))
        }
      }

      # 3. Missing data handling
      if (.has(th$missing_mode)) {
        if (identical(th$missing_mode, "lod")) {
          has_lod <- length(th$lod_value) > 0 && !is.na(th$lod_value)
          # See render_tmt_trends_methods: wording is driven by lod_scope, and
          # runs without it predate the per-sample change and kept a single
          # matrix-wide floor.
          .lod_per_sample <- identical(th$lod_scope %||% "", "per_sample_1st_percentile")
          lod_str <- if (!has_lod) "" else if (.lod_per_sample) {
            glue::glue(", computed separately for each sample as the 1st percentile of that sample's observed values (minimum-detection substitution, equivalent to imputeLCMD::impute.MinDet with q = 0.01); the median per-sample floor was {round(th$lod_value, 4)} log2 units")
          } else {
            glue::glue(" (LOD = {round(th$lod_value, 4)} log2 units, 1st percentile of the matrix)")
          }
          n_om  <- if (length(th$n_omitted) > 0 && !is.na(th$n_omitted)) as.integer(th$n_omitted) else 0L
          n_imp <- if (length(th$n_imputed) > 0 && !is.na(th$n_imputed)) as.integer(th$n_imputed) else 0L
          sentences <- c(sentences, glue::glue(
            "Proteins with fewer than two valid replicates in all conditions were excluded ({n_om} proteins removed); remaining missing values were imputed using a limit-of-detection strategy{lod_str} ({n_imp} values imputed)."
          ))
        } else {
          n_om <- if (length(th$n_omitted) > 0 && !is.na(th$n_omitted)) as.integer(th$n_omitted) else 0L
          sentences <- c(sentences, glue::glue(
            "Proteins with any missing abundance values were excluded prior to visualization ({n_om} proteins removed)."
          ))
        }
      }

      # 4. Protein selection
      if (.has(th$selection_method)) {
        if (identical(th$selection_method, "gene_filter")) {
          gf <- if (length(th$gene_filter) > 0 && any(nzchar(th$gene_filter))) {
            paste(th$gene_filter, collapse = ", ")
          } else "user-specified genes"
          sentences <- c(sentences, glue::glue(
            "The following proteins were highlighted: {gf}."
          ))
        } else {
          n_top_val <- if (length(th$top_n) > 0 && !is.na(th$top_n)) as.integer(th$top_n) else 50L
          sentences <- c(sentences, glue::glue(
            "The top {n_top_val} proteins by within-condition variance were selected for display."
          ))
        }
      }

      # 5. Condition averaging
      if (length(th$avg_by_condition) > 0 && !is.na(th$avg_by_condition)) {
        if (isTRUE(th$avg_by_condition)) {
          sentences <- c(sentences,
            "Abundance values were averaged across biological replicates within each condition prior to display."
          )
        }
      }

      # 6. Row scaling
      if (.has(th$row_scaling) && identical(th$row_scaling, "normalized-abundance")) {
        sentences <- c(sentences,
          "Heatmap cells display limma-summarized log2 abundance values directly (no row-wise scaling), reflecting the absolute magnitude of protein expression."
        )
      } else {
        sentences <- c(sentences,
          "Each protein row was z-score scaled (subtract row mean, divide by row standard deviation) to normalize expression levels for visual comparison."
        )
      }

      # 7. Clustering
      row_clust_str <- if (.has(th$row_clustering)) th$row_clustering else "hierarchical (ward.D2, euclidean)"
      col_clust_str <- if (.has(th$col_clustering)) th$col_clustering else "hierarchical (ward.D2, euclidean)"
      sentences <- c(sentences, glue::glue(
        "Row (protein) clustering used {row_clust_str}; column (sample/condition) ordering used {col_clust_str}."
      ))

      # 8. Summary counts
      has_np <- length(th$n_proteins_displayed) > 0 && !is.na(th$n_proteins_displayed)
      has_ni <- length(th$n_proteins_input)     > 0 && !is.na(th$n_proteins_input)
      if (has_np && has_ni) {
        sentences <- c(sentences, glue::glue(
          "{as.integer(th$n_proteins_displayed)} proteins are displayed (of {as.integer(th$n_proteins_input)} proteins in the input matrix, before missing-data filtering)."
        ))
      }

      # 9. Conditions
      if (length(th$condition_names) > 0 && any(nzchar(th$condition_names))) {
        sentences <- c(sentences, glue::glue(
          "Conditions included: {paste(th$condition_names, collapse = ', ')}."
        ))
      }

      # 10. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a TMT Heatmap analysis to generate methods text.")
    }

    # Spectronaut PTM uploads: sentence for matrices NOT re-centred because the export was
    # cross-run normalized upstream (server_dia_ptm_upload.R step 6). Returns NULL when
    # both matrices were centred (earlier behaviour; old bundles carry no flags). In
    # PTM-only mode there is no global matrix, and only the PTM answer is reported.
    .ptm_upstream_norm_sentence <- function(t) {
      if (is.null(t) || !(isFALSE(t$ptm_centering_applied) || isFALSE(t$global_centering_applied)))
        return(NULL)
      .sd <- function(x) if (.has(x) && !is.na(x)) formatC(as.numeric(x), format = "f", digits = 3) else NULL
      chk <- c(if (!is.null(.sd(t$norm_check_ptm_within_sd)))    paste0("PTM export SD ", .sd(t$norm_check_ptm_within_sd)),
               if (!is.null(.sd(t$norm_check_global_within_sd))) paste0("global export SD ", .sd(t$norm_check_global_within_sd)))
      chk_sent <- if (length(chk))
        paste0(" Before summarization, per-run medians of the features quantified in every run were compared ",
               "between replicate runs of each condition (", paste(chk, collapse = "; "),
               " log2; values of 0.02 or less are consistent with cross-run normalized input).")
      else ""
      if (identical(t$ptm_mode, "ptm_only")) return(paste0(
        "The phospho-enriched export had been cross-run normalized in Spectronaut (as indicated in the ",
        "upload wizard) and was not median-centred again.", chk_sent))
      both <- isFALSE(t$ptm_centering_applied) && isFALSE(t$global_centering_applied)
      if (both) paste0(
        "The phospho-enriched and global protein exports had both been cross-run normalized in Spectronaut ",
        "(the software's default, as indicated in the upload wizard), so neither matrix was median-centred ",
        "again, and the stoichiometry-adjusted values were not re-centred: with both experiments normalized ",
        "upstream, their per-run offsets cancel in the subtraction up to a constant.", chk_sent)
      else paste0(
        "The ", if (isFALSE(t$ptm_centering_applied)) "phospho-enriched" else "global protein",
        " export had been cross-run normalized in Spectronaut (as indicated in the upload wizard) and was not ",
        "median-centred again; the ", if (isFALSE(t$ptm_centering_applied)) "global protein" else "phosphosite",
        " matrix was column-median-centred (every column shifted to the grand median of the column medians).",
        chk_sent)
    }

    # =========================================================================
    # PTM Differential Expression Methods
    # =========================================================================
    render_ptm_de_methods <- function(prov) {
      env <- prov$env
      pd  <- prov$ptm_de

      if (is.null(pd) || !.has(pd$analysis_mode)) return(NULL)

      sentences <- character()

      # 1. PTM type & software introduction
      ptms_str <- if (length(pd$selected_ptms) > 1) {
        paste(pd$selected_ptms, collapse = " and ")
      } else if (.has(pd$selected_ptms)) {
        pd$selected_ptms
      } else {
        "post-translational modification"
      }
      sentences <- c(sentences, glue::glue(
        "Post-translational modification (PTM) differential expression analysis targeting ",
        "{ptms_str} modifications was performed using {pd$software_used}."
      ))

      # 2. Comparison statement
      if (.has(pd$comparison)) {
        sentences <- c(sentences, glue::glue("The comparison tested was: {pd$comparison}."))
      }

      # 3. Analysis mode detail
      mode_detail <- switch(pd$analysis_mode,
        tmt_ptm_limma = {
          .do       <- as.integer(pd$design_option %||% 1L)
          .norm_str <- if (identical(.do, 2L)) {
            "IRS normalization via bridge channel(s) was applied independently to the phospho-enriched and global protein matrices."
          } else if (identical(.do, 3L)) {
            paste0("Column-wise median centering was applied independently to the phospho-enriched and global protein matrices. ",
                   "Plex membership was included as a run covariate in the limma design matrix to account for batch effects across multiple TMT experiments.")
          } else {
            "Column-wise median centering was applied independently to the phospho-enriched and global protein matrices."
          }
          .prob_str <- if (length(pd$prob_threshold) > 0 && !is.na(pd$prob_threshold)) {
            paste0("ptmRS/phosphoRS site localization probabilities were filtered at the observation level (\u2265",
                   as.numeric(pd$prob_threshold),
                   "%); individual peptide-channel rows below the threshold were discarded prior to site aggregation.")
          } else {
            "ptmRS/phosphoRS site localization probabilities were filtered at a user-defined threshold."
          }
          paste0(
            "TMT phosphoproteomics data from Proteome Discoverer were analyzed using a limma-based pipeline. ",
            .prob_str, " ",
            "For each PTM site and TMT channel, the highest-intensity qualifying peptide was selected ",
            "(following the MaxQuant Phospho(STY)Sites convention). ",
            "Intensities were log2-transformed (zero values set to NA). ",
            .norm_str, " ",
            "Stoichiometric adjustment was performed per sample as: adj = log\u2082(phospho site) \u2212 log\u2082(global protein), ",
            "isolating modification-specific fold-change from global protein abundance changes. ",
            "Missing values were not imputed; lmFit accommodates unbalanced missingness natively (MNAR assumption). ",
            "Differential expression was tested using limma: lmFit \u2192 eBayes(robust = TRUE) \u2192 topTable, ",
            "with an intensity-variance trend (trend = TRUE) for the absolute PTM and global protein models and ",
            "trend = FALSE for the stoichiometry-adjusted model, whose row means are log-ratios rather than ",
            "abundances. Benjamini-Hochberg FDR correction was applied independently within each model."
          )
        },
        dia_limma = {
          .prob_str <- if (length(pd$prob_threshold) > 0 && !is.na(pd$prob_threshold)) {
            paste0("Site localization probabilities were assessed by a consensus approach across all samples ",
                   "(maximum probability per candidate position; weakest-link criterion for multi-phosphopeptides) ",
                   "and filtered at \u2265", as.numeric(pd$prob_threshold), "%.")
          } else {
            "Site localization probabilities were filtered at a user-defined threshold."
          }
          paste0(
            "DIA phosphoproteomics data from Spectronaut were analyzed using a limma-based pipeline. ",
            .prob_str, " ",
            "Precursor intensities across all charge states mapping to the same phosphosite were summed per sample ",
            "(linear scale) before log2 transformation. ",
            if (!is.null(.ptm_upstream_norm_sentence(prov$dia_ptm)))
              paste0(.ptm_upstream_norm_sentence(prov$dia_ptm), " ")
            else paste0("Both the phospho-enriched and global protein matrices were independently subjected to ",
                        "column-wise median centering (per-sample median subtracted). "),
            "Sites with fewer than 2 valid values in all conditions were excluded. ",
            "Missing values in the absolute PTM and global protein matrices were imputed via k-nearest neighbors ",
            "(k = max(2, floor(n/2))); the stoichiometric adjusted matrix was left as-is ",
            "(lmFit fits on available observations). ",
            "Stoichiometric adjustment was performed per sample as: adj = log\u2082(phospho) \u2212 log\u2082(protein). ",
            "Differential expression was tested using limma: lmFit \u2192 eBayes(robust = TRUE) \u2192 topTable, ",
            "with an intensity-variance trend (trend = TRUE) for the absolute PTM and global protein models and ",
            "trend = FALSE for the stoichiometry-adjusted model, whose row means are log-ratios rather than ",
            "abundances. Benjamini-Hochberg FDR correction was applied independently within each model."
          )
        },
        dia_pivot_limma = {
          .prob_str <- if (length(pd$prob_threshold) > 0 && !is.na(pd$prob_threshold)) {
            paste0("Site localization probabilities were filtered at \u2265", as.numeric(pd$prob_threshold), "%.")
          } else {
            "Site localization probabilities were filtered at a user-defined threshold."
          }
          paste0(
            "DIA phosphoproteomics data from Spectronaut (PTM Pivot format) were analyzed using a limma-based pipeline. ",
            .prob_str, " ",
            "Phosphosite intensities were summed per site per sample before log2 transformation. ",
            if (!is.null(.ptm_upstream_norm_sentence(prov$dia_ptm)))
              paste0(.ptm_upstream_norm_sentence(prov$dia_ptm), " ")
            else "Both the phospho-enriched and global protein matrices were independently subjected to column-wise median centering. ",
            "Sites with fewer than 2 valid values in all conditions were excluded. ",
            "Missing values in the absolute PTM and global protein matrices were imputed via k-nearest neighbors; ",
            "the stoichiometric adjusted matrix was handled natively by lmFit. ",
            "Stoichiometric adjustment was performed per sample as: adj = log\u2082(phospho) \u2212 log\u2082(protein). ",
            "Differential expression was tested using limma eBayes(robust = TRUE), with an intensity-variance ",
            "trend (trend = TRUE) for the absolute PTM and global protein models and trend = FALSE for the ",
            "stoichiometry-adjusted model, whose row means are log-ratios rather than abundances. ",
            "Benjamini-Hochberg FDR correction was applied independently within each model."
          )
        },
        lf_dda_limma = {
          .prob_str <- if (length(pd$prob_threshold) > 0 && !is.na(pd$prob_threshold)) {
            paste0("Site localization probabilities were filtered at \u2265", as.numeric(pd$prob_threshold), "%.")
          } else {
            "Site localization probabilities were filtered at a user-defined threshold."
          }
          paste0(
            "Label-free DDA phosphoproteomics data from MaxQuant were analyzed using a limma-based pipeline. ",
            .prob_str, " ",
            "Phosphopeptide intensities were summed per site per sample before log2 transformation. ",
            "Both the phospho-enriched and global protein matrices were independently subjected to column-wise median centering. ",
            "Sites with fewer than 2 valid values in all conditions were excluded. ",
            "Missing values were imputed using minimum-probability (MinProb) imputation: values were drawn from a ",
            "per-column downshifted Gaussian (center: column mean \u2212 1.8\u00d7SD; width: 0.3\u00d7SD), following the Perseus convention. ",
            "Stoichiometric adjustment was performed per sample as: adj = log\u2082(phospho) \u2212 log\u2082(protein). ",
            "Differential expression was tested using limma eBayes(trend=TRUE, robust=TRUE) with Benjamini-Hochberg FDR correction."
          )
        },
        glue::glue("Statistical testing was performed using the {pd$analysis_mode} PTM differential expression pipeline.")
      )
      sentences <- c(sentences, mode_detail)

      # 7. Results counts
      has_n   <- length(pd$n_ptm_sites) > 0 && !is.na(pd$n_ptm_sites)
      has_ns  <- length(pd$n_significant) > 0 && !is.na(pd$n_significant)
      if (has_n && has_ns) {
        sentences <- c(sentences, glue::glue(
          "The adjusted (stoichiometric) model included {as.integer(pd$n_ptm_sites)} PTM sites, ",
          "of which {as.integer(pd$n_significant)} were statistically significant ",
          "after Benjamini-Hochberg FDR correction (adj. p-value < 0.05)."
        ))
      }
      has_nu  <- length(pd$n_ptm_sites_unadjusted) > 0 && !is.na(pd$n_ptm_sites_unadjusted)
      has_nsu <- length(pd$n_significant_unadjusted) > 0 && !is.na(pd$n_significant_unadjusted)
      if (has_nu && has_nsu) {
        sentences <- c(sentences, glue::glue(
          "The absolute PTM model contained {as.integer(pd$n_ptm_sites_unadjusted)} sites ",
          "({as.integer(pd$n_significant_unadjusted)} significant)."
        ))
      }
      has_np <- length(pd$n_proteins) > 0 && !is.na(pd$n_proteins)
      if (has_np) {
        sentences <- c(sentences, glue::glue(
          "The global protein model contained {as.integer(pd$n_proteins)} proteins."
        ))
      }

      # 8. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a PTM Differential Expression analysis to generate methods text.")
    }

    # =========================================================================
    # PTM Trend Clusters Methods
    # =========================================================================
    render_ptm_trends_methods <- function(prov) {
      env <- prov$env
      pt  <- prov$ptm_trends

      if (is.null(pt) || length(pt$n_clusters) == 0 || is.na(pt$n_clusters)) return(NULL)

      sentences <- character()

      # 1. PTM type & algorithm introduction
      ptms_str <- if (length(pt$selected_ptms) > 1) {
        paste(pt$selected_ptms, collapse = " and ")
      } else if (.has(pt$selected_ptms)) {
        pt$selected_ptms
      } else {
        "post-translational modification"
      }
      sentences <- c(sentences, glue::glue(
        "Post-translational modification (PTM) trend analysis targeting {ptms_str} modifications ",
        "was performed using stoichiometric (protein-adjusted) abundance profiles."
      ))
      sentences <- c(sentences,
        "Fuzzy c-means clustering was applied to identify co-regulated PTM site groups across experimental conditions."
      )

      # 3. Stoichiometric adjustment
      sentences <- c(sentences,
        "For each PTM site, the adjusted abundance was calculated as PTM log2-abundance minus the corresponding global protein log2-abundance, isolating the modification-specific signal from overall protein expression changes."
      )

      # 4. Missing data handling
      if (.has(pt$missing_mode)) {
        has_nom  <- length(pt$n_omitted) > 0 && !is.na(pt$n_omitted)
        has_nimp <- length(pt$n_imputed) > 0 && !is.na(pt$n_imputed)
        has_lod  <- length(pt$lod_value) > 0 && !is.na(pt$lod_value)
        if (identical(pt$missing_mode, "omit")) {
          if (has_nom) {
            sentences <- c(sentences, glue::glue(
              "{as.integer(pt$n_omitted)} PTM sites with incomplete condition coverage were excluded."
            ))
          }
        } else {
          lod_str <- if (has_lod) glue::glue(" (LOD = {round(pt$lod_value, 4)})") else ""
          imp_str <- if (has_nimp) glue::glue("{as.integer(pt$n_imputed)} values were imputed. ") else ""
          omit_str <- if (has_nom) glue::glue("{as.integer(pt$n_omitted)} PTM sites that failed the existence filter (\\u22652 valid replicates in \\u22651 condition) were excluded prior to imputation. ") else ""
          sentences <- c(sentences, glue::glue(
            "Missing values were imputed at the limit of detection (LOD), defined as the 1st percentile ",
            "of the adjusted abundance matrix{lod_str}. {imp_str}{omit_str}"
          ))
        }
      }

      # 5. Batch correction
      if (isTRUE(pt$batch_corrected)) {
        n_pl <- if (length(pt$n_plexes) > 0 && !is.na(pt$n_plexes)) as.integer(pt$n_plexes) else "multiple"
        sentences <- c(sentences, glue::glue(
          "Per-plex batch effects were removed using limma::removeBatchEffect with condition membership ",
          "as a covariate ({n_pl} plexes)."
        ))
      }

      # 6. Localization probability
      if (length(pt$prob_threshold) > 0 && !is.na(pt$prob_threshold)) {
        sentences <- c(sentences, glue::glue(
          "PTM site localization was filtered at a {as.numeric(pt$prob_threshold)}% probability threshold."
        ))
      }

      # 7. Z-score transform
      sentences <- c(sentences,
        "Abundance values were row-wise z-score normalized (centered to mean, scaled by standard deviation) prior to clustering."
      )

      # 8. Clustering parameters
      has_k <- length(pt$n_clusters) > 0 && !is.na(pt$n_clusters)
      has_m <- length(pt$fuzzifier) > 0 && !is.na(pt$fuzzifier)
      if (has_k && has_m) {
        sentences <- c(sentences, glue::glue(
          "Data were partitioned into {as.integer(pt$n_clusters)} fuzzy clusters using the fuzzy c-means ",
          "algorithm (e1071::cmeans) with fuzzifier exponent m = {pt$fuzzifier} and Euclidean distance."
        ))
      }

      # 9. PTM site counts
      has_nc <- length(pt$n_ptm_sites_clustered) > 0 && !is.na(pt$n_ptm_sites_clustered)
      has_ni <- length(pt$n_ptm_sites_input) > 0 && !is.na(pt$n_ptm_sites_input)
      if (has_nc && has_ni) {
        sentences <- c(sentences, glue::glue(
          "{as.integer(pt$n_ptm_sites_clustered)} PTM sites were included in the final clustering ",
          "(of {as.integer(pt$n_ptm_sites_input)} initial sites)."
        ))
      }

      # 10. Conditions
      if (length(pt$condition_names) > 0) {
        cond_str <- paste(pt$condition_names, collapse = " \u2192 ")
        sentences <- c(sentences, glue::glue("Conditions analyzed: {cond_str}."))
      }

      # 11. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run a PTM Trend Clusters analysis to generate methods text.")
    }

    # =========================================================================
    # PTM Single Inspector Methods
    # =========================================================================
    render_ptm_single_methods <- function(prov) {
      env <- prov$env
      ps  <- prov$ptm_single

      if (is.null(ps) || !.has(ps$ptm_id)) return(NULL)

      sentences <- character()

      # 1. Introduction
      gene_site <- if (.has(ps$gene) && .has(ps$site)) {
        glue::glue("{ps$gene}_{ps$site}")
      } else if (.has(ps$ptm_id)) {
        ps$ptm_id
      } else {
        "the selected PTM"
      }
      mod_str <- if (.has(ps$modification)) glue::glue(" ({ps$modification})") else ""
      sentences <- c(sentences, glue::glue(
        "Single PTM inspection was performed for {gene_site}{mod_str}."
      ))

      # 2. Reference condition & contrasts
      has_ref  <- .has(ps$reference_condition)
      has_nc   <- length(ps$n_contrasts) > 0 && !is.na(ps$n_contrasts)
      if (has_ref && has_nc) {
        sentences <- c(sentences, glue::glue(
          "All non-reference conditions were compared against the reference condition ",
          "via limma-based moderated t-tests (adj.method = BH),",
          " yielding {as.integer(ps$n_contrasts)} pairwise contrast(s)."
        ))
      }

      # 3. Design option
      if (length(ps$design_option) > 0 && !is.na(ps$design_option)) {
        design_str <- switch(as.character(as.integer(ps$design_option)),
          "1" = "single-plex design (one TMT experiment)",
          "2" = "multi-plex design with a bridge/reference channel (reference normalization enabled)",
          "3" = "balanced multi-plex design without a bridge channel",
          paste0("design option ", as.integer(ps$design_option))
        )
        sentences <- c(sentences, glue::glue("The experiment used a {design_str}."))
      }
      if (length(ps$n_plexes) > 0 && !is.na(ps$n_plexes) && ps$n_plexes > 0) {
        sentences <- c(sentences, glue::glue(
          "The experiment included {as.integer(ps$n_plexes)} TMT plex(es)."
        ))
      }

      # 5. Localization probability
      if (length(ps$prob_threshold) > 0 && !is.na(ps$prob_threshold)) {
        sentences <- c(sentences, glue::glue(
          "PTM site localization was filtered at a {as.numeric(ps$prob_threshold)}% probability threshold."
        ))
      }

      # 6. Three-panel output
      sentences <- c(sentences,
        "Three complementary views were generated: (A) Summarized Abundance \u2014 condition-level boxplots showing raw PTM and global protein log2-abundances without imputation; (B) Stoichiometric Shape \u2014 condition-mean line profiles for protein, adjusted PTM (stoichiometric), and total PTM; (C) Mechanistic Fold Changes \u2014 grouped fold-change bars (adjusted, unadjusted, and protein-level log2FC vs. reference)."
      )

      # 7. Data completeness warnings
      if (isFALSE(ps$has_protein)) {
        sentences <- c(sentences,
          "Note: No global protein data was available for the parent protein; adjusted (stoichiometric) results may be absent."
        )
      }
      if (length(ps$missing_conds) > 0) {
        mc_str <- paste(ps$missing_conds, collapse = ", ")
        sentences <- c(sentences, glue::glue(
          "The following conditions lacked sufficient PTM data and were excluded from fold-change calculations: {mc_str}."
        ))
      }

      # 8. Batch correction
      if (isTRUE(ps$batch_corrected)) {
        n_pl <- if (length(ps$n_plexes) > 0 && !is.na(ps$n_plexes)) as.integer(ps$n_plexes) else "multiple"
        sentences <- c(sentences, glue::glue(
          "Per-plex batch effects were removed using limma::removeBatchEffect with condition membership ",
          "as a covariate ({n_pl} plexes)."
        ))
      }

      # 9. Conditions
      if (length(ps$condition_names) > 0) {
        cond_str <- paste(ps$condition_names, collapse = ", ")
        sentences <- c(sentences, glue::glue("Conditions analyzed: {cond_str}."))
      }

      # 10. R version
      if (.has(env$r_version)) {
        sentences <- c(sentences, glue::glue("All analyses were conducted{if (.has(env$osmss_version)) paste0(' with OSMSS (v', env$osmss_version, ')') else ''} in R (v{env$r_version})."))
      }

      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Select a PTM site and run the inspector to generate methods text.")
    }

    # =========================================================================
    # Render output
    # =========================================================================

    # ── DIA PCA methods ──
    render_dia_pca_methods <- function(prov) {
      g   <- prov$dia_global
      p   <- prov$dia_pca
      env <- prov$env

      if (is.null(p) || !.has(p$n_proteins))
        return(list(text = NULL, ebn_list = NULL,
                    empty_msg = "Run DIA PCA to generate methods text."))

      rv_v  <- if (.has(env$r_version)) env$r_version else "?"
      n_p   <- as.integer(p$n_proteins)
      n_s   <- as.integer(p$n_samples)
      pc1v  <- if (.has(p$pc1_var)) round(p$pc1_var, 1) else "?"
      pc2v  <- if (.has(p$pc2_var)) round(p$pc2_var, 1) else "?"
      batch_str <- if (isTRUE(p$batch_corrected))
        " Batch effects were removed for visualization using limma::removeBatchEffect()." else ""

      sentences <- c(
        glue::glue(
          "Principal Component Analysis was performed on DIA protein-level log2 {if (!is.null(g) && isFALSE(g$normalization_applied)) 'intensities (cross-run normalized upstream, not re-centered)' else 'median-centered intensities'} ({n_p} proteins x {n_s} samples)."
        ),
        "Missing values were imputed per row using row-median values; proteins with no valid measurements were excluded.",
        "PCA was computed using prcomp() with centering (center = TRUE) and no variance scaling (scale. = FALSE).",
        if (nchar(batch_str) > 0) batch_str,
        glue::glue("PC1 explained {pc1v}% and PC2 explained {pc2v}% of total variance."),
        if (.has(rv_v)) glue::glue("All analyses were performed in R (v{rv_v}).")
      )
      sentences <- sentences[vapply(sentences, function(s) !is.null(s) && nzchar(s), logical(1))]
      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run DIA PCA to generate methods text.")
    }

    # ── DIA Heatmap methods ──
    render_dia_heatmap_methods <- function(prov) {
      # The unified heatmap module writes rv$provenance$heatmap with
      # source = "dia"; prov$dia_heatmap is the older per-source ledger entry.
      # Accept whichever is present so both wirings produce text.
      h   <- prov$dia_heatmap
      if (!.has(h$n_proteins_displayed) &&
          identical(prov$heatmap$source %||% "", "dia")) h <- prov$heatmap
      env <- prov$env

      if (is.null(h) || !(.has(h$n_proteins_displayed) ||
                          identical(h$source %||% "", "dia")))
        return(list(text = NULL, ebn_list = NULL,
                    empty_msg = "Run DIA Heatmap to generate methods text."))

      rv_v   <- if (.has(env$r_version)) env$r_version else "?"
      n_disp <- if (.has(h$n_proteins_displayed))
        as.integer(h$n_proteins_displayed) else NA_integer_
      scaling <- if (.has(h$row_scaling)) h$row_scaling else "z-score"
      miss_m  <- if (.has(h$missing_mode)) h$missing_mode else "omit"
      batch_str <- if (isTRUE(h$batch_corrected))
        paste0("Batch effects were removed using limma::removeBatchEffect() prior to ",
               "visualization, with experimental condition supplied as the design matrix so ",
               "that condition differences were preserved.") else ""

      miss_str <- if (miss_m == "omit") {
        paste0("Proteins with any missing value across the displayed samples were excluded",
               if (.has(h$n_omitted) && as.integer(h$n_omitted) > 0)
                 paste0(" (", as.integer(h$n_omitted), " proteins removed)") else "",
               ".")
      } else {
        paste0(
          "Missing values were replaced with the 1st percentile of the corresponding ",
          "sample's own observed log2 intensity distribution, a per-sample low-intensity ",
          "substitution",
          if (.has(h$n_imputed) && as.integer(h$n_imputed) > 0)
            paste0("; ", as.integer(h$n_imputed), " values were substituted") else "",
          ". This substitution is applied for display only; no statistical test is ",
          "performed on the heatmap matrix."
        )
      }

      sentences <- c(
        paste0("Protein abundance data from the DIA summarization (log2-transformed, ",
               if (isFALSE(prov$dia_global$normalization_applied))
                 "protein-group intensities cross-run normalized upstream and not re-centered"
               else "column-median-centered protein-group intensities",
               ") were visualized as a clustered heatmap."),
        if (identical(h$selection_method %||% "", "gene_filter")) {
          .gf <- h$gene_filter
          paste0(
            "Proteins were not ranked; ",
            if (!is.na(n_disp)) paste0(n_disp, " protein(s) were ") else "proteins were ",
            "selected manually by gene name",
            if (length(.gf) > 0 && any(nzchar(as.character(.gf))))
              paste0(" (", paste(as.character(.gf), collapse = ", "), ")") else "",
            "."
          )
        } else {
          paste0(
            "A heatmap was generated for ",
            if (!is.na(n_disp)) paste0("the top ", n_disp, " proteins") else
              "the most variable proteins",
            " ranked by across-sample variance. Variance was computed on observed values ",
            "only; substituted values were excluded from the calculation, and proteins with ",
            "fewer than two observed values were not eligible for ranking."
          )
        },
        miss_str,
        if (isTRUE(h$avg_by_condition))
          paste0("Replicate samples were averaged within each condition before display, so ",
                 "each column is a condition mean rather than an individual sample.")
        else NULL,
        if (nchar(batch_str) > 0) batch_str,
        glue::glue("Rows were scaled using {scaling} transformation."),
        if (.has(h$z_score_order))
          glue::glue("Scaling was applied as: {h$z_score_order}.") else NULL,
        if (identical(scaling, "normalized-abundance"))
          paste0("Cells with no measurement are rendered as gaps rather than being assigned ",
                 "a numeric value.") else NULL,
        {
          .row_clustered <- !identical(h$row_clustering %||%
                                       "hierarchical (ward.D2, euclidean)", "none")
          .row_clustered <- .row_clustered &&
            !grepl("^none", h$row_clustering %||% "hierarchical (ward.D2, euclidean)")
          .col_manual <- identical(h$col_clustering %||% "", "manual")
          if (.row_clustered && .col_manual)
            paste0("Rows were ordered by hierarchical clustering (Ward.D2 linkage on ",
                   "Euclidean distance); column order was set manually by the user.")
          else if (.row_clustered)
            paste0("Hierarchical clustering (Ward.D2 linkage on Euclidean distance) was ",
                   "applied to both rows and columns. Clustering was used for display ",
                   "ordering only: no dendrogram was cut and no cluster assignments were ",
                   "derived.")
          else if (.col_manual)
            paste0("Too few rows were displayed for row clustering, so rows are shown in ",
                   "their original order; column order was set manually by the user.")
          else
            paste0("Too few rows were displayed for row clustering, so rows are shown in ",
                   "their original order. Columns were ordered by hierarchical clustering ",
                   "(Ward.D2 linkage on Euclidean distance), for display ordering only.")
        },
        if (length(h$condition_names) > 0 && any(nzchar(as.character(h$condition_names))))
          glue::glue("Conditions included: {paste(h$condition_names, collapse = ', ')}.")
        else NULL,
        if (.has(rv_v)) glue::glue("All analyses were performed in R (v{rv_v}).")
      )
      sentences <- sentences[vapply(sentences, function(s) !is.null(s) && nzchar(s), logical(1))]
      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run DIA Heatmap to generate methods text.")
    }

    # ── DIA Trend Clusters methods ──
    render_dia_trends_methods <- function(prov) {
      # See render_dia_heatmap_methods: the unified trend-clusters module writes
      # rv$provenance$trend_clusters with source = "dia".
      t   <- prov$dia_trends
      if (!.has(t$n_proteins_clustered) &&
          identical(prov$trend_clusters$source %||% "", "dia"))
        t <- prov$trend_clusters
      env <- prov$env

      if (is.null(t) || !.has(t$n_proteins_clustered))
        return(list(text = NULL, ebn_list = NULL,
                    empty_msg = "Run DIA Trend Clusters to generate methods text."))

      rv_v  <- if (.has(env$r_version)) env$r_version else "?"
      n_cl  <- as.integer(t$n_clusters)
      n_p   <- as.integer(t$n_proteins_clustered)
      m_val <- if (.has(t$fuzzifier)) round(t$fuzzifier, 2) else "?"
      conds <- if (.has(t$condition_names)) paste(t$condition_names, collapse = ", ") else "?"

      sentences <- c(
        glue::glue(
          "Protein abundance trend analysis was performed using fuzzy c-means clustering (e1071::cmeans)."
        ),
        glue::glue(
          "Condition-mean log2 intensities were computed from the DIA summarization ({n_p} proteins across {conds})."
        ),
        if (.has(t$min_replicate_frac))
          glue::glue(
            "A condition mean was formed only where at least ",
            "{round(as.numeric(t$min_replicate_frac) * 100)}% of that condition's replicates ",
            "were observed; below that threshold the condition mean was treated as missing."
          ) else NULL,
        if (.has(t$n_proteins_routed_out) && as.integer(t$n_proteins_routed_out) > 0)
          glue::glue(
            "Clustering required a supported mean in every condition. ",
            "{as.integer(t$n_proteins_routed_out)} proteins did not meet this requirement and ",
            "were reported separately as a descriptive detection table rather than being ",
            "imputed into the clustering, so that proteins undetectable in one condition were ",
            "not assigned a fabricated profile."
          ) else NULL,
        if (.has(t$sd_prefilter) && as.numeric(t$sd_prefilter) > 0)
          glue::glue(
            "Proteins with a standard deviation below {as.numeric(t$sd_prefilter)} log2 units ",
            "across condition means were excluded before clustering ",
            "({as.integer(t$n_removed_low_sd %||% 0L)} proteins)."
          ) else NULL,
        "Each protein's expression profile was z-score normalized (row-wise) prior to clustering.",
        glue::glue(
          "Proteins were partitioned into k = {n_cl} clusters using fuzzy c-means with ",
          "Euclidean distance and fuzzifier m = {m_val}, estimated from the number of ",
          "features and the number of conditions using the method of Schw\u00e4mmle and Jensen."
        ),
        if (.has(t$rng_seed))
          glue::glue(
            "Cluster centres are initialised at random; the random number generator was ",
            "seeded (seed = {as.integer(t$rng_seed)}) so that repeated runs on the same data ",
            "return identical partitions and identical cluster numbering."
          ) else NULL,
        if (.has(t$k_validity_min_at))
          glue::glue(
            "The number of clusters was specified by the user. For reference, a Xie-Beni ",
            "validity scan over k = 2 to 10, using the same fuzzifier and random seed, ",
            "reached its minimum at k = {as.integer(t$k_validity_min_at)}; this was reported ",
            "only and did not alter the k used."
          ) else NULL,
        if (.has(t$membership_max_median))
          glue::glue(
            "Fuzzy c-means assigns each protein a membership value in every cluster; ",
            "proteins were assigned to the cluster of highest membership. The median ",
            "maximum membership was {as.numeric(t$membership_max_median)}",
            if (.has(t$membership_frac_above_half))
              glue::glue(", and {round(as.numeric(t$membership_frac_above_half) * 100, 1)}% ",
                         "of clustered proteins exceeded a maximum membership of 0.5")
            else "",
            "."
          ) else NULL,
        if (.has(t$membership_threshold) && as.numeric(t$membership_threshold) > 0)
          glue::glue(
            "In the profile plots, individual protein traces are drawn only for proteins ",
            "whose maximum membership reaches {as.numeric(t$membership_threshold)}; the ",
            "cluster average line is computed over all members of the cluster, including ",
            "those not drawn. The panel header reports the number drawn against the ",
            "cluster total."
          ) else NULL,
        if (.has(rv_v)) glue::glue("All analyses were performed in R (v{rv_v}).")
      )
      sentences <- sentences[vapply(sentences, function(s) !is.null(s) && nzchar(s), logical(1))]
      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Run DIA Trend Clusters to generate methods text.")
    }

    # ── DIA PTM render functions ──────────────────────────────────────────────

    render_dia_ptm_de_methods <- function(prov) {
      t   <- prov$dia_ptm
      pd  <- prov$ptm_de     # DE-stage provenance: filter, imputation, counts
      env <- prov$env
      if (is.null(t) || !.has(t$n_sites_output))
        return(list(text = NULL, ebn_list = NULL,
                    empty_msg = "Complete the PTM Upload Wizard and run DE to generate methods text."))

      rv_v   <- if (.has(env$r_version))    env$r_version   else "?"
      lm_v   <- if (.has(env$limma_version)) env$limma_version else "?"
      n_si   <- if (.has(t$n_sites_input))   as.integer(t$n_sites_input)   else "?"
      n_so   <- if (.has(t$n_sites_passing_prob)) as.integer(t$n_sites_passing_prob)
                else if (.has(t$n_sites_output))  as.integer(t$n_sites_output)  else "?"
      n_adj  <- if (.has(t$n_sites_adjusted)) as.integer(t$n_sites_adjusted) else "?"
      prob   <- if (.has(t$prob_threshold))   suppressWarnings(as.numeric(t$prob_threshold)) else NA_real_
      conds  <- if (.has(t$condition_names))  paste(t$condition_names, collapse = ", ") else "?"
      paired <- isTRUE(t$pairing)
      is_dda    <- identical(t$acquisition_mode, "DDA")
      is_tmt_pd <- identical(t$acquisition_mode, "TMT_PD")

      if (is_dda) {
        # ── DDA MaxQuant methods text ──
        sentences <- c(
          glue::glue(
            "Phosphoproteomics data were acquired by data-dependent acquisition (DDA) and ",
            "processed with MaxQuant / Perseus. Site-level PTM intensities and matched ",
            "whole-proteome protein-group intensities were exported as tab-separated matrices."
          ),
          glue::glue(
            "A total of {n_si} phosphosite entries were parsed from the PTM matrix. ",
            if (.has(t$has_probability) && isTRUE(t$has_probability))
              glue::glue("Sites with localization probability < {prob}% were excluded, retaining {n_so} high-confidence phosphosites. ")
            else
              glue::glue("No localization probability column was present; all {n_so} sites were retained. "),
            "Non-phospho modifications (e.g., Carbamidomethyl) were removed prior to analysis."
          ),
          glue::glue(
            "Missing values in both the PTM and global protein matrices were imputed ",
            "independently using the Perseus MinProb method: for each sample column, values ",
            "were drawn from a Gaussian distribution with mean = column_mean \u2212 1.8 \u00d7 column_SD ",
            "and SD = 0.3 \u00d7 column_SD (random seed 42). This models left-censored (MNAR) ",
            "missingness expected from DDA acquisition."
          ),
          glue::glue(
            "Stoichiometric adjustment was performed by subtracting the matched global protein ",
            "log2 abundance from each phosphosite log2 abundance (PTM_log2 \u2212 global_log2), ",
            "yielding {n_adj} adjusted site rows. Protein identifiers were resolved through ",
            "the PTM.CollapseKey column; shared-peptide sites were expanded per protein group."
          ),
          if (paired)
            "A paired experimental design was accounted for using limma::duplicateCorrelation()."
          else NULL,
          glue::glue(
            "Differential analysis was performed with limma (v{lm_v}) using empirical Bayes ",
            "moderation with robust = TRUE throughout. An intensity-variance trend ",
            "(trend = TRUE) was fitted for the absolute PTM and global protein models. The ",
            "stoichiometry-adjusted model was moderated without a trend (trend = FALSE), ",
            "because its row means are log-ratios rather than abundances and therefore carry ",
            "no intensity information for a trend to be fitted against. Multiple testing ",
            "correction was applied using the Benjamini-Hochberg procedure, independently ",
            "within each model. All analyses were performed in R (v{rv_v})."
          )
        )
      } else if (is_tmt_pd) {
        # ── TMT Proteome Discoverer methods text ──
        n_plx    <- if (.has(t$n_plexes))      as.integer(t$n_plexes)      else 1L
        dopt     <- if (.has(t$design_option)) as.integer(t$design_option) else 1L
        ptms     <- if (.has(t$selected_ptms)) paste(t$selected_ptms, collapse = " and ") else "Phospho"
        plex_str <- if (n_plx == 1L) "a single TMT plex" else paste0(n_plx, " TMT plexes")

        # Every number below is read from provenance. Do not hard-code a count, a
        # threshold or a parameter value that provenance already carries: that is
        # how this branch came to name a column the parser never reads, and to
        # describe a shared-peptide expansion this route does not perform. Where a
        # field is absent the sentence depending on it is dropped rather than
        # asserted with a default - a wrong number in a methods section is worse
        # than a missing one.
        .abund_src  <- if (.has(t$abundance_source)) as.character(t$abundance_source) else NA_character_
        .n_pep_in   <- if (.has(t$n_pep_groups_input))     as.integer(t$n_pep_groups_input)     else NULL
        .n_pep_kp   <- if (.has(t$n_pep_groups_conf_kept)) as.integer(t$n_pep_groups_conf_kept) else NULL
        .pep_basis  <- if (.has(t$pep_confidence_basis))   as.character(t$pep_confidence_basis) else NA_character_
        .pepq       <- if (.has(t$pep_qvalue_max))  as.numeric(t$pep_qvalue_max)  else NULL
        .protq      <- if (.has(t$prot_qvalue_max)) as.numeric(t$prot_qvalue_max) else NULL
        .interf_on  <- isTRUE(t$interference_filtered)
        .interf_th  <- if (.has(t$interference_threshold)) as.numeric(t$interference_threshold) else NULL
        .n_int_rm   <- if (.has(t$n_interf_removed))   as.integer(t$n_interf_removed)   else NULL
        .n_int_na   <- if (.has(t$n_interf_kept_na))   as.integer(t$n_interf_kept_na)   else NULL
        .n_obs_coll <- if (.has(t$n_obs_collected))   as.integer(t$n_obs_collected)   else NULL
        .n_obs_low  <- if (.has(t$n_obs_below_own_threshold)) as.integer(t$n_obs_below_own_threshold) else NULL
        .n_obs_np   <- if (.has(t$n_obs_no_own_prob)) as.integer(t$n_obs_no_own_prob) else NULL
        .n_razor    <- if (.has(t$n_razor_reassigned)) as.integer(t$n_razor_reassigned) else NULL
        .n_shared   <- if (.has(t$n_shared_peptide_sites)) as.integer(t$n_shared_peptide_sites) else NULL
        .n_multipf  <- if (.has(t$n_sites_multi_peptidoform)) as.integer(t$n_sites_multi_peptidoform) else NULL
        .min_valid  <- if (.has(t$min_valid_per_cond)) as.integer(t$min_valid_per_cond) else NULL
        .n_aftval   <- if (.has(t$n_sites_after_valid_filter)) as.integer(t$n_sites_after_valid_filter) else NULL
        .n_plexdet  <- if (.has(t$n_plexes_detected)) as.integer(t$n_plexes_detected) else NULL
        .n_dropplex <- if (.has(t$n_sites_dropped_plex)) as.integer(t$n_sites_dropped_plex) else NULL
        .site_frac  <- if (.has(t$site_min_plex_frac)) as.numeric(t$site_min_plex_frac) else NULL
        .prot_frac  <- if (.has(t$prot_min_plex_frac)) as.numeric(t$prot_min_plex_frac) else NULL
        .site_min_n <- if (.has(t$site_min_plex_n))    as.integer(t$site_min_plex_n)    else NULL
        .n_post_irs   <- if (.has(t$n_sites_dropped_post_irs)) as.integer(t$n_sites_dropped_post_irs) else NULL
        .irs_policy   <- if (.has(t$irs_policy)) as.character(t$irs_policy) else NA_character_
        .irs_blanked  <- if (.has(t$irs_n_blanked)) as.integer(t$irs_n_blanked) else NULL
        .irs_nobridge <- if (.has(t$irs_n_plex_no_bridge)) as.integer(t$irs_n_plex_no_bridge) else NULL
        .norm_basis <- if (.has(t$normalization_basis)) as.character(t$normalization_basis) else NULL
        .n_nopar    <- if (.has(t$n_no_parent_protein)) as.integer(t$n_no_parent_protein) else NULL
        .ctr_sprd   <- if (.has(t$adj_centering_delta))
                         round(diff(range(as.numeric(t$adj_centering_delta), na.rm = TRUE)), 3) else NULL
        # DE-stage fields, from prov$ptm_de.
        .de_min     <- if (.has(pd$de_min_valid_per_cond)) as.integer(pd$de_min_valid_per_cond) else NULL
        .de_eff     <- if (.has(pd$de_min_valid_per_cond_effective)) pd$de_min_valid_per_cond_effective else NULL
        .de_eff_str <- if (!is.null(.de_eff))
                         paste0(names(.de_eff), " ", as.integer(.de_eff), collapse = "; ") else NULL
        .n_test     <- if (.has(pd$n_ptm_sites))           as.integer(pd$n_ptm_sites)           else NULL
        .n_est      <- if (.has(pd$n_ptm_sites_estimable)) as.integer(pd$n_ptm_sites_estimable) else NULL
        .imp_scp    <- if (.has(pd$imputation_scope)) as.character(pd$imputation_scope) else NA_character_
        .plex_cov   <- isTRUE(pd$plex_covariate_used)
        .n_batches  <- if (.has(pd$n_plex_batches)) as.integer(pd$n_plex_batches) else NULL
        .n_nmar     <- if (.has(pd$n_nmar)) as.integer(pd$n_nmar) else NULL
        .nmar_pow   <- pd$nmar_has_power
        .nmar_minp  <- if (.has(pd$nmar_min_achievable_p)) signif(as.numeric(pd$nmar_min_achievable_p), 3) else NULL

        sentences <- c(
          glue::glue(
            "Phosphoproteomics data were acquired using isobaric tandem mass tag (TMT) labeling ",
            "and quantified with Proteome Discoverer. Paired phospho-enriched and global ",
            "(unenriched) fractions were exported for {plex_str} as nested Proteome Discoverer ",
            "CSV reports, each carrying a protein-level section and a nested Peptide Groups ",
            "section. Quantification was taken from the per-channel reporter-ion columns of the ",
            "Peptide Groups section (columns named 'Abundance: F<plex>: <tag>, ...'), which are ",
            "the raw, unnormalized reporter intensities",
            if (identical(.abund_src, "scaled_only"))
              "; no raw abundance columns were present in this export and Proteome Discoverer's 'Abundances (Scaled)' columns were used instead, so the normalization described below was applied on top of Proteome Discoverer's own"
            else "",
            "."
          ),
          tmt_excluded_sentence(t$excluded_channels %||% character(0)),
          if (!is.null(.n_pep_in) && identical(.pep_basis, "qvalue"))
            glue::glue(
              "Peptide groups were filtered on identification confidence using the Proteome ",
              "Discoverer Qvality q-value at q < {.pepq}, retaining {.n_pep_kp} of {.n_pep_in} ",
              "peptide groups across all phospho-enriched files."
            )
          else if (!is.null(.n_pep_in) && identical(.pep_basis, "confidence"))
            glue::glue(
              "No Qvality q-value column was present in the export, so peptide groups were ",
              "filtered on the Proteome Discoverer Confidence tier, retaining only ",
              "high-confidence groups ({.n_pep_kp} of {.n_pep_in} across all phospho-enriched files)."
            )
          else if (!is.null(.n_pep_in))
            glue::glue(
              "No peptide-level identification-confidence column was present in the export, so no ",
              "peptide confidence filter was applied downstream; all {.n_pep_in} exported peptide ",
              "groups were carried forward with the confidence already imposed by the Proteome ",
              "Discoverer workflow."
            )
          else NULL,
          if (.interf_on && !is.null(.interf_th))
            glue::glue(
              "Co-isolation interference was filtered at the peptide-group level: groups with a ",
              "reported isolation interference above {.interf_th}% were removed",
              if (!is.null(.n_int_rm))
                glue::glue(" ({.n_int_rm} peptide group{if (.n_int_rm == 1L) '' else 's'})")
              else "",
              ". Groups carrying no interference value were retained rather than discarded",
              if (!is.null(.n_int_na))
                glue::glue(" ({.n_int_na} peptide group{if (.n_int_na == 1L) '' else 's'})")
              else "",
              ", so the filter removes measured co-isolation and does not act as a completeness ",
              "requirement."
            )
          else
            glue::glue(
              "No isolation-interference column was found in the export, so no co-isolation ",
              "filter was applied."
            ),
          glue::glue(
            "Peptide-to-site assignment was resolved by parsing the 'Modifications in Master ",
            "Proteins' column, which encodes each modification with its residue position, ",
            "modification type and ptmRS (or IMP-ptmRS) site-localization probability. ",
            "{n_si} distinct candidate {ptms} site{if (identical(n_si, 1L)) '' else 's'} ",
            "{if (identical(n_si, 1L)) 'was' else 'were'} parsed",
            if (!is.null(.n_obs_coll))
              glue::glue(" from {.n_obs_coll} peptide-group-by-channel observations")
            else "",
            "."
          ),
          glue::glue(
            "Localization was decided per site rather than per observation, following the class-I ",
            "convention: a site was retained when its maximum ptmRS probability across all ",
            "observations reached {prob}%, and its intensities were then taken from every channel ",
            "in which it was measured. {n_so} of {n_si} candidate sites met that criterion",
            if (!is.null(.n_obs_low) && !is.null(.n_obs_np))
              glue::glue(". Because the criterion is applied per site, the retained set includes ",
                         "individual observations whose own recorded probability was below the ",
                         "threshold ({.n_obs_low}) or absent ({.n_obs_np})")
            else "",
            "."
          ),
          glue::glue(
            "Peptides matching more than one master protein were assigned to a single protein by ",
            "the razor rule: the master protein with the greatest total PSM support, summed ",
            "across every phospho-enriched file so that the assignment is identical in all plexes",
            if (!is.null(.n_razor))
              glue::glue(" ({.n_razor} site{if (.n_razor == 1L) '' else 's'} ",
                         "{if (.n_razor == 1L) 'was' else 'were'} assigned to a ",
                         "protein other than the first listed in the export)")
            else "",
            ". Each such site is therefore tested once, under its razor protein; the remaining ",
            "master proteins are retained and reported in the Also_Maps_To column of the results ",
            "table, and sites resting on a shared peptide are flagged Is_Shared = TRUE",
            if (!is.null(.n_shared))
              glue::glue(" ({.n_shared} site{if (.n_shared == 1L) '' else 's'})")
            else "",
            ". This follows the MaxQuant Phospho (STY)Sites convention of one row per site ",
            "carrying a designated leading protein alongside the full list of matching proteins. ",
            "There is at present no single agreed standard for shared-peptide handling in PTM ",
            "quantification: other tools assign the measurement to every matching protein, or ",
            "merge the matching proteins into one combined site identifier, and under either of ",
            "those policies the same measurement would appear as more than one tested row. The ",
            "choice determines which protein a site is reported under and the size of the tested ",
            "set; it does not change the measured intensities."
          ),
          glue::glue(
            "For each retained site, peptide-group reporter intensities mapping to that site ",
            "within a channel were summed on the linear scale and the sum was log2-transformed. ",
            "Distinct peptidoforms of the same site are summed together rather than one being ",
            "selected over another, and none is discarded",
            if (!is.null(.n_multipf))
              glue::glue(" ({.n_multipf} of the retained sites carry more than one distinct peptidoform)")
            else "",
            "."
          ),
          if (!is.null(.min_valid))
            glue::glue(
              "Sites were retained at summarization when quantified in at least {.min_valid} ",
              "channel{if (.min_valid == 1L) '' else 's'} of at least one condition",
              if (!is.null(.n_aftval))
                glue::glue(", leaving {.n_aftval} site{if (.n_aftval == 1L) '' else 's'}")
              else "",
              ". This filter is deliberately looser than the per-comparison filter applied at the ",
              "differential-testing stage."
            )
          else NULL,
          if (!is.null(.n_plexdet) && .n_plexdet >= 2L) {
            if (!is.null(.site_frac) && .site_frac > 0)
              glue::glue(
                "Sites quantified in too few plexes were removed on a completeness fraction of ",
                "{.site_frac}",
                if (!is.null(.site_min_n))
                  glue::glue(" (at least {.site_min_n} of {.n_plexdet} plexes)")
                else "",
                if (!is.null(.n_dropplex))
                  glue::glue(": {.n_dropplex} site{if (.n_dropplex == 1L) '' else 's'} ",
                             "{if (.n_dropplex == 1L) 'was' else 'were'} dropped")
                else "",
                if (!is.null(.prot_frac))
                  glue::glue("; the matched global protein matrix used a completeness fraction of {.prot_frac}")
                else "",
                "."
              )
            # A bundle predating the site_min_plex_frac provenance field tells us
            # nothing about which rule ran, so assert nothing. Claiming no filter
            # was applied to a run that may well have applied one is the same
            # failure this branch was rewritten to remove. Both inner conditions
            # false leaves the block NULL, which c() drops.
            else if (!is.null(.site_frac))
              glue::glue(
                "No plex-coverage filter was applied: a site was not required to be quantified ",
                "in any minimum number of plexes. Retention was governed by the per-condition ",
                "valid-value rule stated above and, at the testing stage, by the per-comparison ",
                "rule applied to each contrast. This follows the reference implementation for ",
                "multi-plex isobaric data ",
                "(MSstatsTMT), which applies no run-coverage filter and instead reduces the ",
                "model for a feature measured in a single run. A coverage requirement would ",
                "preferentially remove sites whose compared conditions sit within a single ",
                "plex, and those are the sites least exposed to between-plex variation, since ",
                "both arms share a batch."
              )
          } else NULL,
          if (dopt == 2L && !is.null(.irs_blanked) && identical(.irs_policy, "na"))
            glue::glue(
              "Internal Reference Scaling derives a per-site correction from that site's ",
              "bridge-channel value in each plex, so a site with no bridge value in a plex ",
              "cannot be aligned there. Rather than removing such sites, the affected ",
              "site-by-plex blocks were set missing and each site was retained on the plexes ",
              "in which it could be corrected ({.irs_blanked} site-by-plex ",
              "block{if (.irs_blanked == 1L) '' else 's'} set missing)",
              if (!is.null(.n_post_irs) && .n_post_irs > 0L)
                glue::glue("; {.n_post_irs} site{if (.n_post_irs == 1L) '' else 's'} then fell ",
                           "below the per-condition valid-value rule and ",
                           "{if (.n_post_irs == 1L) 'was' else 'were'} removed")
              else "",
              ". This mirrors the reference-channel normalization in MSstatsTMT, which drops the ",
              "affected run for a feature rather than the feature itself.",
              if (!is.null(.irs_nobridge) && .irs_nobridge > 0L)
                glue::glue(" {.irs_nobridge} plex carried no bridge channel at all; its values ",
                           "could not be aligned to the other plexes and were set missing ",
                           "throughout.")
              else ""
            )
          else NULL,
          glue::glue(
            "Channel loading was corrected by median centering on the log2 scale: a per-channel ",
            "shift was applied so that every channel reaches a common reference level",
            if (!is.null(.norm_basis))
              glue::glue(", with the shifts derived from {.norm_basis}")
            else "",
            ". Only biological channels set the reference level; bridge channels were ",
            "brought to the same level without contributing to it.",
            if (dopt == 2L)
              " Plex-to-plex differences were then removed by Internal Reference Scaling: for each site, the geometric mean of its bridge-channel intensities across the plexes in which it carried one was taken as the reference, and each of those plexes was scaled to it. Bridge channels were excluded from the output matrices after correction."
            else if (dopt == 3L)
              " Plex-to-plex differences were not corrected arithmetically; plex was instead carried as a fixed covariate in the linear model, so that inter-plex loading differences are absorbed by the model rather than by a per-feature correction."
            else ""
          ),
          if (!is.null(.protq))
            glue::glue(
              "Matched global protein abundances were summarized from the unenriched fractions of ",
              "the same plexes, with protein-level rows filtered at q < {.protq}, normalized on ",
              "the same biological channels, and subject to the same normalization route as the ",
              "site matrix."
            )
          else NULL,
          glue::glue(
            "Stoichiometric adjustment was performed per sample by subtracting the matched parent ",
            "protein log2 abundance from each site log2 abundance (site_log2 \u2212 protein_log2), ",
            "yielding {n_adj} adjusted rows. Each measured site yields exactly one adjusted row: ",
            "sites mapping to more than one master protein are not expanded into separate ",
            "hypothesis rows, so the adjusted row count equals the number of measured sites",
            if (!is.null(.n_nopar) && .n_nopar > 0L)
              glue::glue(". {.n_nopar} site{if (.n_nopar == 1L) '' else 's'} had no quantified ",
                         "parent protein in the global export and {if (.n_nopar == 1L) 'is' else 'are'} ",
                         "all-missing in the adjusted matrix")
            else "",
            "."
          ),
          if (!is.null(.ctr_sprd))
            glue::glue(
              "Because the site and protein matrices were centred independently, their per-channel ",
              "offsets do not cancel under subtraction, so the adjusted matrix was ",
              "column-median-centred again after the subtraction to remove that residual offset ",
              "(range of the applied shifts: {.ctr_sprd} log2 units)."
            )
          else NULL,
          if (!is.null(.de_min))
            glue::glue(
              "For each pairwise comparison, a site was carried into the linear model only when ",
              "BOTH compared conditions carried at least {.de_min} valid ",
              "value{if (.de_min == 1L) '' else 's'}",
              if (!is.null(.de_eff_str))
                glue::glue(" (effective per-condition requirement: {.de_eff_str})")
              else "",
              ". The threshold is a flat count and is not scaled by the replicate number. A fold ",
              "change requires both arms; sites meeting the threshold in only one condition are ",
              "excluded from the model and assessed instead for a non-random detection pattern by Fisher's exact test on detected / not-detected counts."
            )
          else NULL,
          if (identical(.imp_scp, "none"))
            glue::glue(
              "No missing-value imputation was performed on any of the three matrices. All models ",
              "were fitted on the observed values only, so a site contributes to a coefficient ",
              "only in the samples where it was actually quantified."
            )
          else NULL,
          if (.plex_cov)
            glue::glue(
              "Plex was included in the design matrix as a fixed batch covariate",
              if (!is.null(.n_batches)) glue::glue(" ({.n_batches} plexes)") else "",
              ", so inter-plex loading differences are absorbed by the model rather than inflating ",
              "within-condition variance. Sites quantified in a single plex are rank-deficient in ",
              "the plex coefficients; limma sets those coefficients to NA, and the condition ",
              "contrast, which places zero weight on them, is still estimated."
            )
          else NULL,
          if (!is.null(.n_nmar) && .n_nmar > 0L)
            glue::glue(
              "Sites excluded by the both-arms filter cannot be fitted and were instead assessed ",
              "for a non-random detection pattern by Fisher's exact test on detected / ",
              "not-detected counts, corrected by the Benjamini-Hochberg procedure within that ",
              "excluded set only ({.n_nmar} site{if (.n_nmar == 1L) '' else 's'}). This correction ",
              "is independent of the correction applied to the modelled set, so no joint false ",
              "discovery rate is claimed across the two",
              if (identical(.nmar_pow, FALSE) && !is.null(.nmar_minp))
                glue::glue(". At this replicate count the smallest attainable two-sided Fisher ",
                           "p-value is {.nmar_minp}, so no member of this set can reach ",
                           "significance; p-values are reported as NA and the set is presented as ",
                           "a ranked descriptive list of detection frequencies")
              else "",
              "."
            )
          else NULL,
          if (!is.null(.n_test) && !is.null(.n_est) && .n_est < .n_test)
            glue::glue(
              "Of the {.n_test} adjusted rows carried into testing, {.n_est} returned an estimable ",
              "coefficient; the remainder appear in the exported table without statistics."
            )
          else NULL,
          if (paired)
            "A paired experimental design was accounted for using limma::duplicateCorrelation()."
          else NULL,
          glue::glue(
            "Differential analysis was performed with limma (v{lm_v}). Linear models were fitted ",
            "by lmFit and moderated by empirical Bayes with robust = TRUE throughout, at limma's ",
            "default hyperparameters (proportion = 0.01, stdev.coef.lim = c(0.1, 4), ",
            "winsor.tail.p = c(0.05, 0.1)). An intensity-variance trend (trend = TRUE) was fitted ",
            "for the absolute site and global protein models; the stoichiometry-adjusted model was ",
            "moderated without a trend (trend = FALSE), because its row means are log-ratios ",
            "rather than abundances and carry no intensity information for a trend to be fitted ",
            "against. Statistics were extracted with topTable(coef = 1, number = Inf, ",
            "sort.by = \"none\", adjust.method = \"BH\")."
          ),
          glue::glue(
            "Three models were fitted on the same samples: the stoichiometry-adjusted sites, the ",
            "absolute site abundances, and the global protein abundances. They are built from the ",
            "same measurements and are not independent evidence of one another. The ",
            "Benjamini-Hochberg procedure was applied independently within each result table; no ",
            "correction was applied across the three tables or across comparisons."
          ),
          glue::glue(
            "Scope. Quantification is read at the peptide-group level as exported by Proteome ",
            "Discoverer; PSM-level records are not imported, so PSM-level filtering, ",
            "re-quantification and per-PSM interference handling are those applied by the ",
            "Proteome Discoverer workflow and were not repeated here. No contaminant or decoy ",
            "filtering was applied downstream of the export, so entries flagged as contaminants ",
            "by the search are retained unless they were removed upstream. Proteome Discoverer's ",
            "own 'Abundances (Grouped)', 'Abundances (Scaled)', 'Abundance Ratio', 'Abundance ",
            "Ratio (log2)' and 'Abundance Ratio Adj. P-Value' columns were deliberately not used: ",
            "all normalization, adjustment and testing reported here were performed downstream ",
            "from the raw per-channel abundances."
          ),
          if (.has(rv_v)) glue::glue("All analyses were performed in R (v{rv_v}).") else NULL
        )
      } else {
        # ── DIA Spectronaut methods text ──
        # Every number below is read from provenance. Do not hard-code a count,
        # a parameter value or a threshold that provenance already carries: that
        # is how this branch came to describe a median polish it never ran.
        .de_min   <- if (.has(pd$de_min_valid_per_cond)) as.integer(pd$de_min_valid_per_cond) else NULL
        .knn_cap  <- if (.has(pd$knn_max_missing_per_cond)) as.integer(pd$knn_max_missing_per_cond) else NULL
        .imp_scp  <- if (.has(pd$imputation_scope)) as.character(pd$imputation_scope) else NA_character_
        .n_test   <- if (.has(pd$n_ptm_sites)) as.integer(pd$n_ptm_sites) else NULL
        .n_est    <- if (.has(pd$n_ptm_sites_estimable)) as.integer(pd$n_ptm_sites_estimable) else NULL
        .n_nopar  <- if (.has(t$n_no_parent_protein)) as.integer(t$n_no_parent_protein) else NULL
        .n_distin <- if (.has(pd$n_distinct_measured_sites)) as.integer(pd$n_distinct_measured_sites) else NULL
        .n_multipf<- if (.has(t$n_sites_multi_peptidoform_retained))
                       as.integer(t$n_sites_multi_peptidoform_retained) else NULL
        .n_lowloc <- if (.has(t$n_obs_below_own_threshold)) as.integer(t$n_obs_below_own_threshold) else NULL
        .n_obs_coll   <- if (.has(t$n_obs_collected))   as.integer(t$n_obs_collected)   else NULL
        .n_obs_noprob <- if (.has(t$n_obs_no_own_prob)) as.integer(t$n_obs_no_own_prob) else NULL
        .norm_basis   <- if (.has(t$normalization_basis)) as.character(t$normalization_basis) else NULL
        .norm_delta   <- if (.has(t$normalization_delta_log2))
                           round(as.numeric(t$normalization_delta_log2), 3) else NULL
        .n_fill_ptm   <- if (.has(pd$n_knn_values_cache_ptm))     as.integer(pd$n_knn_values_cache_ptm)     else NULL
        .n_fill_prot  <- if (.has(pd$n_knn_values_cache_protein)) as.integer(pd$n_knn_values_cache_protein) else NULL
        .n_mpfb       <- if (.has(t$n_medpolish_fallback))        as.integer(t$n_medpolish_fallback)        else NULL
        .ctr_sprd <- if (.has(t$adj_centering_delta))
                       round(diff(range(as.numeric(t$adj_centering_delta), na.rm = TRUE)), 3) else NULL
        .n_skip_ptm  <- if (.has(pd$n_knn_rows_skipped_cache_ptm))
                          as.integer(pd$n_knn_rows_skipped_cache_ptm) else NULL
        .n_skip_prot <- if (.has(pd$n_knn_rows_skipped_cache_protein))
                          as.integer(pd$n_knn_rows_skipped_cache_protein) else NULL
        # Effective per-condition missingness thresholds. Named vector; the names
        # are the condition labels, so both are rendered together.
        .de_eff     <- if (.has(pd$de_min_valid_per_cond_effective))
                         pd$de_min_valid_per_cond_effective else NULL
        .de_eff_str <- if (!is.null(.de_eff))
                         paste0(names(.de_eff), " ", as.integer(.de_eff),
                                collapse = "; ") else NULL
        # Site summarization method. Where absent, neither method is asserted.
        .site_summ  <- if (.has(t$site_summarization))
                         as.character(t$site_summarization) else NA_character_
        # PTM-only mode (no global proteome): one absolute PTM model, no stoichiometry.
        # Provenance written before the mode existed has no ptm_mode and had a global file.
        .ptm_only   <- identical(t$ptm_mode, "ptm_only")

        # This branch renders for BOTH the Spectronaut precursor route and the
        # Spectronaut PTM pivot route, which run different summarizers. Pick the
        # phrasing from the route discriminator in provenance. data_level is the
        # precise field ("precursor" vs "site (pivot)"); software is the fallback.
        # When neither is present, "unknown" selects the phrasing that asserts the
        # least - never the precursor wording.
        .route <- if (.has(t$data_level)) {
                    if (grepl("pivot",     as.character(t$data_level), ignore.case = TRUE)) "pivot"
                    else if (grepl("precursor", as.character(t$data_level), ignore.case = TRUE)) "precursor"
                    else "unknown"
                  } else if (.has(t$software)) {
                    if (grepl("pivot",       as.character(t$software), ignore.case = TRUE)) "pivot"
                    else if (grepl("spectronaut", as.character(t$software), ignore.case = TRUE)) "precursor"
                    else "unknown"
                  } else "unknown"
        .has_prob <- !(.has(t$has_probability) && !isTRUE(t$has_probability))

        sentences <- c(
          if (identical(.route, "precursor"))
            glue::glue(
              "Phosphoproteomics data were acquired using Spectronaut (DIA) and exported as a ",
              "precursor-level report containing per-sample EG.TotalQuantity, ",
              "EG.PTMPositions, EG.PTMProbabilities, and EG.PTMSites columns for Phospho (STY)."
            )
          else if (identical(.route, "pivot"))
            glue::glue(
              "Phosphoproteomics data were acquired using Spectronaut (DIA) and exported as a ",
              "site-level PTM pivot report, with one row per phosphosite collapse key and one ",
              "quantity column per sample."
            )
          else
            glue::glue(
              "Phosphoproteomics data were acquired by data-independent acquisition (DIA) and ",
              "exported for site-level analysis."
            ),
          glue::glue(
            if (identical(.route, "precursor"))
              "A total of {n_si} precursor-level observations were parsed. "
            else if (identical(.route, "pivot"))
              "A total of {n_si} site-level entries were parsed. "
            else
              "A total of {n_si} entries were parsed. ",
            if (.has_prob)
              "Phosphosites with maximum localization probability < {prob}% were excluded, retaining {n_so} high-confidence phosphosites."
            else
              "No localization probability column was present in the export, so all {n_so} sites were retained."
          ),
          if (identical(.route, "precursor"))
            glue::glue(
              "Phosphosite localization was determined by consensus across all biological replicates: ",
              "for each precursor row, the maximum localization probability per candidate position ",
              "was computed across all samples, and the top-N positions (by peptide multiplicity) ",
              "selected from these consensus probabilities were used as the site identifier for every ",
              "sample in that row. Ties were broken deterministically in the N-to-C terminal direction. ",
              "Because the criterion is applied per site across all samples rather than per observation, ",
              "once a site passes its intensity is taken from every matched sample, including samples ",
              "whose own localization probability for that row was absent or below the threshold",
              if (!is.null(.n_lowloc) && !is.null(.n_obs_coll))
                glue::glue(" ({.n_lowloc} of {.n_obs_coll} retained precursor-by-sample ",
                           "observations had a recorded probability below the threshold; ",
                           "{.n_obs_noprob} had none recorded)")
              else if (!is.null(.n_lowloc))
                glue::glue(" ({.n_lowloc} of the retained precursor-by-sample observations)")
              else "",
              "."
            )
          else if (identical(.route, "pivot"))
            glue::glue(
              "Phosphosite localization was resolved by Spectronaut before export: site identities ",
              "were taken directly from the collapse-key column and no consensus was computed across ",
              "precursor rows. ",
              if (.has_prob)
                "Where several exported rows shared a collapse key, the maximum localization probability across those rows was used as that site's probability for the threshold above."
              else
                "The export carried no per-site localization probability, so no probability threshold was applied."
            )
          else NULL,
          if (identical(.route, "pivot"))
            glue::glue(
              "For each site, all exported intensities sharing that site's collapse key within a ",
              "sample were summed on the linear scale and the sum was then log2-transformed. ",
              "Summation is keyed on the site identifier, so rows differing only in their reported ",
              "parent protein are summed together."
            )
          else if (identical(.site_summ, "medpolish"))
            glue::glue(
              "For each site, precursor intensities were summed on the linear scale across charge ",
              "states within each peptidoform and log2-transformed. Site-level abundance per sample ",
              "was then estimated by Tukey median polish (stats::medpolish) over the peptidoform by ",
              "sample matrix, taking the overall term plus the per-sample column effect. A ",
              "peptidoform detected in some samples but not others therefore loads onto the ",
              "peptidoform term rather than shifting the site estimate, so the site value is ",
              "insensitive to which peptidoforms were detected where",
              if (!is.null(.n_multipf))
                glue::glue(" - the case for {.n_multipf} of the retained sites, which carry more than one distinct modified sequence")
              else "",
              ". Sites with a single peptidoform reduce to the log2 of the summed intensities.",
              if (!is.null(.n_mpfb) && .n_mpfb > 0L)
                glue::glue(" Median polish did not return a finite estimate for {.n_mpfb} site",
                           "{if (.n_mpfb == 1L) '' else 's'}, which fell back to summation.")
              else ""
            )
          else if (identical(.site_summ, "sum"))
            glue::glue(
              "For each site, all precursor intensities mapping to that site within a sample were ",
              "summed on the linear scale and the sum was then log2-transformed. Summation is keyed ",
              "on the site identifier rather than on the peptidoform, so distinct modified sequences ",
              "reporting the same site are summed together",
              if (!is.null(.n_multipf))
                glue::glue(" ({.n_multipf} of the retained sites carry more than one distinct modified sequence)")
              else "",
              "."
            )
          else
            glue::glue(
              "For each site, precursor intensities were aggregated to a single log2 abundance per ",
              "sample."
            ),
          if (!is.null(.ptm_upstream_norm_sentence(t)))
            .ptm_upstream_norm_sentence(t)
          else if (.ptm_only)
            "The phosphosite matrix was column-median-centred: every column was shifted to the grand median of the column medians."
          else glue::glue(
            "The phosphosite and global protein matrices were each column-median-centred ",
            "independently: every column was shifted to the grand median of the column medians. ",
            if (!is.null(.norm_basis))
              "The global protein matrix was centred on {.norm_basis}"
            else
              "The global protein matrix was centred on a complete-case basis where enough complete rows were available and on an all-observed basis otherwise",
            if (!is.null(.norm_basis) && !is.null(.norm_delta))
              "; the maximum per-column discrepancy between the complete-case and all-observed centring bases was {.norm_delta} log2 units"
            else "",
            "."
          ),
          # adj_centering_delta is written only by the summarizer that performs the
          # post-subtraction centring, so its presence is what makes this sentence true.
          if (!is.null(.ctr_sprd))
            glue::glue(
              "Because the two matrices are centred independently, their per-column offsets do not ",
              "cancel under subtraction, so the stoichiometry-adjusted matrix was column-median-centred ",
              "again after the subtraction to remove that residual offset (range of the applied ",
              "shifts: {.ctr_sprd} log2 units)."
            )
          else NULL,
          if (.ptm_only)
            glue::glue(
              "No matched global proteome was supplied, so no stoichiometric adjustment was ",
              "performed: phosphosite abundances were analysed as absolute values, and a change in ",
              "a site therefore includes any change in the abundance of its parent protein. ",
              "Peptides with multiple protein identifications were not expanded into per-protein ",
              "rows; such a site was tested once, under all of its proteins, and flagged as shared."
            )
          else glue::glue(
            "Stoichiometric adjustment was performed by subtracting the corresponding global ",
            "protein log2 abundance (from a matched whole-proteome Spectronaut export) from each ",
            "phosphosite abundance, yielding {n_adj} adjusted site rows at summarization",
            if (!is.null(.n_nopar))
              glue::glue(", of which {.n_nopar} row{if (.n_nopar == 1L) '' else 's'} ",
                         "{if (.n_nopar == 1L) 'was' else 'were'} all-missing because the parent ",
                         "protein was not quantified in the global export, and ",
                         "{if (.n_nopar == 1L) 'was' else 'were'} dropped at the missingness filter")
            else "",
            ".",
            if (identical(t$data_level, "precursor"))
              glue::glue(
                " Peptides with multiple protein identifications were expanded into per-protein ",
                "hypothesis rows (suffixed '_SharedInference'), so the tested row count exceeds the ",
                "number of independently measured sites",
                if (!is.null(.n_distin))
                  glue::glue(" ({.n_distin} distinct measured sites among the rows tested)")
                else "",
                ". There is at present no single agreed standard for shared-peptide handling in ",
                "PTM quantification. Expansion is one of several policies in use; others assign ",
                "the measurement to a single protein by a razor rule, or merge the matching ",
                "proteins into one combined site identifier. Under expansion the rows sharing a ",
                "measurement are not independent hypotheses, so the multiple-testing correction ",
                "is applied over a row set larger than the number of distinct measurements. This ",
                "is disclosed rather than corrected: the policy affects the size of the tested ",
                "set and which protein a site is reported under, not the measured intensities."
              )
            else if (identical(t$data_level, "site (pivot)"))
              glue::glue(
                " Each phosphosite yields a single adjusted row; sites mapping to more than one ",
                "parent protein are not expanded into separate hypothesis rows, so the tested row ",
                "count equals the number of measured sites."
              )
            else ""
          ),
          if (!is.null(.de_min))
            glue::glue(
              "Sites were retained for testing where at least one condition carried a valid ",
              "value in at least half of that condition's samples, subject to a floor of ",
              "{.de_min}",
              if (!is.null(.de_eff))
                glue::glue(" (effective per-condition requirement: {.de_eff_str})")
              else "",
              "."
            )
          else NULL,
          glue::glue(
            if (.ptm_only) "Missing values in the absolute phosphosite matrix were imputed by "
            else "Missing values in the absolute phosphosite and global protein matrices were imputed by ",
            "k-nearest neighbours (impute::impute.knn, k = 10, colmax = 0.8, maxp = Inf, seed 42)",
            if (!is.null(.knn_cap))
              glue::glue(", restricted to rows with at most {.knn_cap} missing ",
                         "{if (.knn_cap == 1L) 'value' else 'values'} in every condition; rows above ",
                         "that cap retained their missing values and were fitted on the observations that exist")
            else "",
            if (identical(.imp_scp, "cache"))
              ", with donor neighbours drawn from every sample in the dataset rather than from the two conditions being compared, so imputed values do not depend on which comparison is run"
            else if (identical(.imp_scp, "contrast"))
              ", with donor neighbours drawn from the two conditions being compared"
            else "",
            if (!is.null(.n_skip_ptm))
              glue::glue("; {.n_skip_ptm} phosphosite row{if (.n_skip_ptm == 1L) '' else 's'} ",
                         "{if (.n_skip_ptm == 1L) 'was' else 'were'} excluded from imputation by that cap",
                         if (!.ptm_only && !is.null(.n_skip_prot))
                           glue::glue(", and {.n_skip_prot} protein row{if (.n_skip_prot == 1L) '' else 's'}")
                         else "",
                         if (.ptm_only && !is.null(.n_fill_ptm))
                           glue::glue("; k-NN filled {.n_fill_ptm} phosphosite values in total")
                         else if (!is.null(.n_fill_ptm) && !is.null(.n_fill_prot))
                           glue::glue("; k-NN filled {.n_fill_ptm} phosphosite and ",
                                      "{.n_fill_prot} protein values in total")
                         else "")
            else "",
            if (.ptm_only) "."
            else ". The stoichiometry-adjusted matrix was not imputed and was fitted on observed values only."
          ),
          if (isTRUE(pd$knn_cap_exhaustive))
            glue::glue(
              "Note: at this replicate count and missingness level the eligibility cap excluded ",
              "essentially every row from imputation, so the k-nearest-neighbour step filled little ",
              "or nothing and the models were fitted on observed values throughout."
            )
          else NULL,
          if (isTRUE(pd$knn_cap_degenerate))
            glue::glue(
              "Note: at this replicate count the imputation eligibility cap admits every row that ",
              "passed the missingness filter, so it did not further restrict imputation."
            )
          else NULL,
          if (.ptm_only && !is.null(.n_test) && !is.null(.n_est))
            glue::glue("Of the {.n_test} sites carried into testing after the missingness ",
                       "filter, {.n_est} returned an estimable coefficient.")
          else if (!is.null(.n_test) && !is.null(.n_est) && .n_est < .n_test)
            glue::glue(
              "Because the adjusted matrix is not imputed, sites with no observations in one of the ",
              "two contrasted conditions yield no estimate: of the {.n_test} adjusted rows carried ",
              "into testing after the missingness filter, {.n_est} returned an estimable ",
              "coefficient. The remainder appear in the exported table without statistics."
            )
          else if (!is.null(.n_test) && !is.null(.n_est))
            glue::glue("Of the {.n_test} adjusted rows carried into testing after the missingness ",
                       "filter, {.n_est} returned an estimable coefficient.")
          else NULL,
          if (.ptm_only)
            glue::glue(
              "Sites for which more than 75% of values in at least one contrasted condition were ",
              "missing prior to imputation are marked with a distinct symbol in the volcano display ",
              "and can be hidden interactively. Under the imputation rule above such rows are never ",
              "imputed, so the marker indicates heavy missingness rather than imputation. ",
              "A single model was fitted, on the absolute phosphosite abundances. No correction ",
              "was applied across comparisons."
            )
          else glue::glue(
            "Sites for which more than 75% of values in at least one contrasted condition were ",
            "missing prior to imputation are marked with a distinct symbol in the volcano display ",
            "and can be hidden interactively. Under the imputation rule above such rows are never ",
            "imputed, so the marker indicates heavy missingness rather than imputation; for ",
            "the stoichiometry-adjusted table, where no imputation is performed, it is computed ",
            "from the absolute phosphosite matrix. ",
            "Three models were fitted on the same samples: the stoichiometry-adjusted sites, the ",
            "absolute phosphosite abundances, and the global protein abundances. They are built ",
            "from the same measurements and are not independent evidence of one another. No ",
            "correction was applied across the three result tables or across comparisons."
          ),
          if (paired)
            "A paired experimental design was accounted for using limma::duplicateCorrelation()."
          else NULL,
          if (.ptm_only)
            glue::glue(
              "Differential analysis was performed with limma (v{lm_v}) using empirical Bayes ",
              "moderation with robust = TRUE and an intensity-variance trend (trend = TRUE). ",
              "Multiple testing correction was applied using the Benjamini-Hochberg procedure. ",
              "All analyses were performed in R (v{rv_v})."
            )
          else glue::glue(
            "Differential analysis was performed with limma (v{lm_v}) using empirical Bayes ",
            "moderation with robust = TRUE throughout. An intensity-variance trend ",
            "(trend = TRUE) was fitted for the absolute PTM and global protein models. The ",
            "stoichiometry-adjusted model was moderated without a trend (trend = FALSE), ",
            "because its row means are log-ratios rather than abundances and therefore carry ",
            "no intensity information for a trend to be fitted against. Multiple testing ",
            "correction was applied using the Benjamini-Hochberg procedure, independently ",
            "within each model. All analyses were performed in R (v{rv_v})."
          )
        )
      }

      sentences <- sentences[vapply(sentences, function(s) !is.null(s) && nzchar(s), logical(1))]
      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Complete PTM DE to generate methods text.")
    }

    render_dia_ptm_trends_methods <- function(prov) {
      t   <- prov$dia_ptm
      env <- prov$env
      if (is.null(t) || !.has(t$n_sites_output))
        return(list(text = NULL, ebn_list = NULL,
                    empty_msg = "Complete the PTM Upload Wizard and run PTM Trends to generate methods text."))

      rv_v   <- if (.has(env$r_version)) env$r_version else "?"
      n_so   <- if (.has(t$n_sites_output)) as.integer(t$n_sites_output) else "?"
      conds  <- if (.has(t$condition_names)) paste(t$condition_names, collapse = ", ") else "?"
      is_dda <- identical(t$acquisition_mode, "DDA")
      # The acquisition mode has three values, not two. Reading it as a DDA/DIA
      # binary described every TMT Proteome Discoverer dataset as Spectronaut DIA.
      is_tmt_pd_tr <- identical(t$acquisition_mode, "TMT_PD")
      src_lbl <- if (is_dda) "DDA MaxQuant / Perseus"
                 else if (is_tmt_pd_tr) "TMT Proteome Discoverer"
                 else "DIA Spectronaut"
      imp_note <- if (is_dda)
        " Missing values were imputed using MinProb (Perseus-style downshifted Gaussian) prior to stoichiometric adjustment."
      else ""
      # PTM-only mode: the absolute PTM matrix is clustered (no global proteome)
      ptm_only <- identical(t$ptm_mode, "ptm_only")
      val_lbl  <- if (ptm_only) "absolute" else "adjusted"

      r <- prov$dia_ptm_trends_run
      rule_sent <- if (!is.null(r)) {
        .thr <- if (!is.null(r$min_valid_per_cond))
          paste0(" (realised thresholds: ",
                 paste0(names(r$min_valid_per_cond), " ",
                        as.integer(r$min_valid_per_cond), collapse = "; "), ")")
        else ""
        glue::glue(
          "No missing values were imputed. A PTM site entered the clustering only if it ",
          "was quantified in at least ",
          if (identical(r$stringency, "strict")) "all replicates"
          else "half of the replicates (minimum two, capped at the replicate count)",
          " of every condition{.thr}; {as.integer(r$n_sites_excluded)} of ",
          "{as.integer(r$n_sites_input)} sites failed this rule and are reported in the ",
          "excluded-site table rather than being filled in or silently dropped. ",
          "{as.integer(r$n_sites_clustered)} sites were clustered."
        )
      } else ""

      param_sent <- if (!is.null(r)) {
        # Rounded: an estimated fuzzifier is a full-precision double and would
        # otherwise render as 1.38724193847 in a methods paragraph.
        .fz     <- if (.has(r$fuzzifier)) signif(as.numeric(r$fuzzifier), 4) else "?"
        .fz_src <- if (.has(r$fuzzifier_source)) as.character(r$fuzzifier_source) else NA_character_
        .fz_txt <- if (identical(.fz_src, "estimated (Schwaemmle-Jensen)"))
          glue::glue(", fuzzifier m = {.fz}, estimated from the dimensions of the ",
                     "clustered matrix by the method of Schwaemmle and Jensen ",
                     "(Bioinformatics 2010) rather than fixed at a default")
        else if (identical(.fz_src, "user-specified"))
          glue::glue(", fuzzifier m = {.fz}, set by the analyst")
        else
          glue::glue(", fuzzifier m = {.fz}")
        .zv     <- if (.has(r$n_sites_zero_variance)) as.integer(r$n_sites_zero_variance) else NULL
        glue::glue(
          "Per-condition mean {val_lbl} log2 abundances were z-score normalized row-wise and ",
          "partitioned by fuzzy c-means (e1071::cmeans, Euclidean distance, k = ",
          "{as.integer(r$n_clusters)}{.fz_txt}, iter.max = 200, ",
          "random seed {as.integer(r$seed)}).",
          if (!is.null(.zv) && .zv > 0L)
            glue::glue(" {.zv} site{if (.zv == 1L) '' else 's'} with no between-condition ",
                       "variance {if (.zv == 1L) 'was' else 'were'} excluded before ",
                       "clustering, because a flat profile z-scores to zero throughout and ",
                       "would anchor a centroid while carrying no trend.")
          else "",
          if (.has(r$k_selection))
            " The number of clusters was chosen by the analyst on inspection of the resulting profiles, not by a cluster-validity criterion; it should be read as a display parameter rather than an estimate of the number of underlying response classes."
          else "",
          " Clusters are numbered by descending mean ",
          "z-score at the last condition, so cluster indices depend on the condition order ",
          "shown. Individual site traces are displayed for sites whose maximum membership ",
          "is at least {r$membership_thresh}; this is a display threshold and does not ",
          "affect cluster assignment or the reported cluster sizes."
        )
      } else glue::glue(
        "Per-condition mean {val_lbl} log2 abundances were computed and z-score normalized ",
        "(row-wise) prior to fuzzy c-means clustering (e1071::cmeans, Euclidean distance)."
      )

      sentences <- c(
        if (ptm_only)
          glue::glue(
            "PTM trend analysis was performed on absolute PTM site abundances ",
            "({n_so} sites across conditions: {conds}) derived from {src_lbl} data. No matched ",
            "global proteome was supplied, so the abundances are not adjusted for protein ",
            "abundance and a site's profile includes any change in its parent protein."
          )
        else glue::glue(
          "PTM trend analysis was performed on stoichiometry-adjusted PTM site abundances ",
          "({n_so} sites across conditions: {conds}) derived from {src_lbl} data.{imp_note}"
        ),
        if (nzchar(rule_sent)) rule_sent else NULL,
        if (!is.null(r) && identical(r$site_gate, "de_sig") &&
            !is.null(r$site_gate_applied) && grepl("^de_sig", r$site_gate_applied))
          glue::glue(
            "Clustering was restricted to PTM sites reaching an adjusted p-value below ",
            "0.05 in the ",
            if (ptm_only) "absolute PTM" else "stoichiometry-adjusted",
            " differential expression model of at least ",
            "one tested contrast ({r$site_gate_applied}); sites with no significant change ",
            "were excluded, because row-wise z-scoring rescales every profile to unit ",
            "amplitude and would otherwise give unchanged sites the same apparent response ",
            "as responsive ones."
          )
        else NULL,
        param_sent,
        if (.has(rv_v)) glue::glue("All analyses were performed in R (v{rv_v}).") else NULL
      )
      sentences <- sentences[vapply(sentences, function(s) !is.null(s) && nzchar(s), logical(1))]
      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Complete PTM Trends to generate methods text.")
    }

    render_dia_ptm_single_methods <- function(prov) {
      t   <- prov$dia_ptm
      env <- prov$env
      if (is.null(t) || !.has(t$n_sites_output))
        return(list(text = NULL, ebn_list = NULL,
                    empty_msg = "Complete the PTM Upload Wizard and inspect a site to generate methods text."))

      rv_v   <- if (.has(env$r_version)) env$r_version else "?"
      prob   <- if (.has(t$prob_threshold)) suppressWarnings(as.numeric(t$prob_threshold)) else NA_real_
      # The acquisition mode has three values, not two. Read as a DDA/DIA binary,
      # this described every TMT Proteome Discoverer dataset as Spectronaut DIA and
      # claimed its quantities were TMP-summarized precursors. They are neither:
      # TMT reads peptide-group reporter-ion intensities and sums them per site.
      is_dda    <- identical(t$acquisition_mode, "DDA")
      is_tmt_pd <- identical(t$acquisition_mode, "TMT_PD")
      src_lbl <- if (is_dda) "DDA MaxQuant / Perseus dataset"
                 else if (is_tmt_pd) "TMT Proteome Discoverer dataset"
                 else "DIA Spectronaut dataset"
      quant_lbl <- if (is_dda)
        "log2 site intensities (per-sample sum, MinProb-imputed)"
      else if (is_tmt_pd)
        "log2 site abundances (peptide-group reporter-ion intensities summed within each site and channel)"
      else
        "log2-scale absolute precursor-level quantities (TMP-summarized)"

      s <- prov$dia_ptm_single_run
      stat_sent <- if (is.null(s)) NULL
        else if (!isTRUE(s$de_available) && isTRUE(s$nmar_available))
          glue::glue(
            "Abundances are plotted as measured values only; no value shown in this view is ",
            "imputed. This site was not fitted by the linear model: it did not meet the ",
            "valid-value requirement in both compared conditions, so it was assessed for a ",
            "non-random detection pattern by Fisher's exact test on detected / not-detected ",
            "counts, corrected within that excluded set only. No fold change or moderated ",
            "p-value is reported for it, and no joint false discovery rate is claimed across ",
            "the modelled and presence/absence sets. The reference condition was ",
            "{s$reference_condition}."
          )
        else if (isTRUE(s$de_available))
          glue::glue(
            "Abundances are plotted as measured values only; no value shown in this view is ",
            "imputed. Fold changes and adjusted p-values are not recomputed here - they are ",
            "read from the limma fit produced for the same contrast by the PTM differential ",
            "expression analysis (empirical Bayes moderation, Benjamini-Hochberg correction ",
            "across all tested sites), so a site's statistics are identical in both views. ",
            "The reference condition was {s$reference_condition} ({as.integer(s$n_contrasts)} ",
            "contrast(s))."
          )
        else
          glue::glue(
            "Abundances are plotted as measured values only; no value shown in this view is ",
            "imputed. No differential expression result was stored for the selected reference ",
            "condition ({s$reference_condition}), so no fold change or p-value is reported for ",
            "this site; the measured abundances and per-condition detection counts are shown ",
            "without inference."
          )

      # Modification label: written by the Single run since the module stopped
      # hard-coding "Phospho". Absent on older provenance, in which case the
      # generic wording is used rather than asserting a modification.
      .mod_txt <- if (.has(s$modification)) as.character(s$modification) else NA_character_
      # Plex labelling is a TMT-only display behaviour and only engages when the
      # run actually spans two or more plexes; n_plexes comes from the wizard.
      .n_plx   <- if (.has(t$n_plexes)) as.integer(t$n_plexes) else NA_integer_

      sentences <- c(
        glue::glue(
          if (!is.na(.mod_txt) && !identical(.mod_txt, "Phospho"))
            "Individual {.mod_txt} site abundances were visualized using {quant_lbl} "
          else
            "Individual phosphosite abundances were visualized using {quant_lbl} ",
          "from the {src_lbl}. ",
          if (is_tmt_pd)
            glue::glue("Sites were filtered on ptmRS site-localization probability at a site ",
                       "maximum of {prob}% or greater.")
          else if (!is_dda)
            glue::glue("Sites were filtered to localization probability >= {prob}%.")
          else ""
        ),
        if (is_tmt_pd && !is.na(.n_plx) && .n_plx >= 2L)
          glue::glue(
            "Replicate labels in the abundance panel carry the TMT plex each channel ",
            "belongs to. The values shown are the measured, summarized abundances, so ",
            "they still contain any between-plex offset; that offset is removed at the ",
            "model stage (by internal reference scaling, or by carrying plex as a ",
            "covariate) and is therefore absent from the fold changes reported below, ",
            "but not from the points above."
          )
        else NULL,
        if (is_tmt_pd)
          glue::glue(
            "Where the site rests on a peptide shared by more than one master protein, ",
            "the protein shown is the razor assignment and the alternatives are listed ",
            "alongside it; the measured quantities are identical under either assignment."
          )
        else NULL,
        stat_sent,
        if (identical(t$ptm_mode, "ptm_only"))
          paste0("No matched global proteome was supplied, so only absolute site abundances and ",
                 "fold changes are shown; they are not adjusted for the abundance of the parent protein.")
        else NULL,
        if (.has(rv_v)) glue::glue("Visualization was performed in R (v{rv_v}).") else NULL
      )
      sentences <- sentences[vapply(sentences, function(s) !is.null(s) && nzchar(s), logical(1))]
      list(text = paste(sentences, collapse = " "), ebn_list = NULL,
           empty_msg = "Inspect a PTM site to generate methods text.")
    }

    output$methods_md <- renderUI({
      prov <- reactiveValuesToList(provenance)
      # Use is.reactive() first (handles modern Shiny R6 reactives that may not
      # pass is.function()), then fall back to is.function() for plain lambdas,
      # then treat as a literal string.
      sec  <- if (shiny::is.reactive(section) || is.function(section)) section() else section

      result <- switch(sec,
        de             = render_de_methods(prov),
        pca            = render_pca_methods(prov),
        trend_clusters = render_trend_methods(prov),
        heatmap        = render_heatmap_methods(prov),
        rank_fc        = render_rank_fc_methods(prov),
        trends         = render_trends_methods(prov),
        tmt_volcano    = render_tmt_volcano_methods(prov),
        tmt_trends     = render_tmt_trends_methods(prov),
        tmt_pca        = render_tmt_pca_methods(prov),
        tmt_heatmap    = render_tmt_heatmap_methods(prov),
        ptm_de         = render_ptm_de_methods(prov),
        ptm_trends     = render_ptm_trends_methods(prov),
        ptm_single     = render_ptm_single_methods(prov),
        dia_pca        = render_dia_pca_methods(prov),
        dia_heatmap    = render_dia_heatmap_methods(prov),
        dia_trends     = render_dia_trends_methods(prov),
        dia_ptm_de     = render_dia_ptm_de_methods(prov),
        dia_ptm_trends = render_dia_ptm_trends_methods(prov),
        dia_ptm_single = render_dia_ptm_single_methods(prov),
        render_de_methods(prov)
      )

      # Determine the empty message AFTER obtaining the result, so that
      # render functions returning list(text=NULL, empty_msg="...") can
      # surface their specific guidance text.
      empty_msg <- if (!is.null(result) && !is.null(result$empty_msg))
        result$empty_msg
      else
        "Run an analysis to generate methods text."

      # Show placeholder when: result is NULL, OR result$text is NULL/blank
      # (the latter happens for DIA paths before the wizard is completed).
      if (is.null(result) || is.null(result$text) || !nzchar(trimws(result$text))) {
        return(tags$p(class = "text-muted fst-italic", empty_msg))
      }

      main_text <- result$text
      ebn_list  <- result$ebn_list

      # For DE: prepend global context sentence
      if (sec == "de") {
        global_text <- NULL
        g <- prov$global
        if (!is.null(g) && length(g$input_filename) > 0 && nzchar(g$input_filename)) {
          global_text <- glue::glue(
            "Protein-level quantification data was exported from {g$input_software_source} ",
            "(\"{g$input_filename}\", {as.integer(g$total_initial_proteins)} proteins)."
          )
        }
        main_text <- paste(c(global_text, main_text), collapse = " ")
      }

      # Build UI elements
      ui_elements <- list(tags$p(main_text))

      # Append the invariant reference panel list if >5 proteins
      if (!is.null(ebn_list) && length(ebn_list) > 0) {
        ui_elements <- c(ui_elements, .ebn_list_widget(ebn_list, result$ebn_ip_sub))
      }

      tags$div(
        style = "background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 12px 16px; font-size: 0.92em; line-height: 1.6;",
        tagList(ui_elements)
      )
    })
  })
}
