# server_residue_prediction.R ==================================================
# All server logic for the Residue PTM Prediction module set:
#   - Kinase Prediction, 14-3-3, ScanSite, Multi-Site, Protein Kinome
# Called from server.R as: server_residue_prediction(input, output, session, rv, logger)
# =============================================================================

server_residue_prediction <- function(input, output, session, rv, logger) {
  
  # LOG OUTPUTS
  output$console_log_kin      <- renderText({ paste(rv$log_text, collapse = "\n") })
  output$console_log_1433     <- renderText({ paste(rv$log_text, collapse = "\n") })
  output$console_log_bind     <- renderText({ paste(rv$log_text, collapse = "\n") })
  output$console_log_multi    <- renderText({ paste(rv$log_text, collapse = "\n") })
  output$console_log_map      <- renderText({ paste(rv$log_text, collapse = "\n") })
  
  # INITIALIZATION
  observe({
    isolate({
      if(file.exists(rv$path_kin_fam)) {
        db <- load_kinase_families(rv$path_kin_fam)
        if(!is.null(db)) rv$kin_fam_db <- db
      }
      # TKL kinase list — ensure_global_resources() triggers lazy-load if needed
      k_list <- fetch_tkl_kinase_list()
      if(length(k_list) > 0) updateSelectizeInput(session, "tkl_kinase", choices = k_list, selected = "CDK2", server = TRUE)
    })
  })
  
  # PURGE cache when user changes target inputs (debounced to avoid per-keystroke firing)
  .input_key_raw <- reactive({
    req(input$uniprot_id, input$gene_name)
    paste0(input$uniprot_id, "_", input$gene_name)
  })
  .input_key <- debounce(.input_key_raw, 1000)  # wait 1s after last keystroke

  observe({
    current_input_key <- .input_key()
    if (!is.null(rv$cache_key) && rv$cache_key != current_input_key) {
      logger(paste0("!!! INPUT CHANGE DETECTED: Purging Cache for ", rv$cache_key))
      rv$cache_key <- NULL; rv$cache_seq <- NULL; rv$cache_cons <- NULL; rv$cache_acc <- NULL; rv$cache_struct <- NULL
    }
  })
  
  # RUN ANALYSIS LOGIC
  observeEvent(input$run_btn, {
    req(input$uniprot_id, input$gene_name)
    
    # Only handle tabs that belong to this module
    residue_tabs <- c("Kinase Prediction", "14-3-3 Prediction",
                       "ScanSite - Upstream/Downstream",
                       "Multi-Site Analysis", "Protein Kinome")
    if (!(input$suite_nav %in% residue_tabs)) return()
    
    rv$log_text <- c(paste0(">>> INITIALIZING ", input$gene_name))

    # Reference data each user downloads (Kinase Library, PhosphoSitePlus): stop and explain if missing
    need_tkl <- input$suite_nav %in% c("Multi-Site Analysis", "Protein Kinome") ||
                (input$suite_nav == "Kinase Prediction" && grepl("TKL", input$mode_select))
    need_psp <- input$suite_nav == "Kinase Prediction" && grepl("PSP|Combined", input$mode_select)
    miss_tkl <- need_tkl && !tkl_data_available()
    miss_psp <- need_psp && !psp_data_available()
    if (miss_tkl || miss_psp) {
      logger("!!! Reference data missing (Kinase Library / PhosphoSitePlus). See the download instructions.")
      showModal(missing_data_modal(need_tkl = miss_tkl, need_psp = miss_psp))
      return()
    }

    withProgress(message = 'Initializing...', value = 0, {
      progress_wrapper <- function(msg) { incProgress(0, detail = msg) }
      
      # Ensure run-level cache is filled
      current_key <- paste0(input$uniprot_id, "_", input$gene_name)
      if (is.null(rv$cache_key) || rv$cache_key != current_key) {
        logger(paste0(">>> CACHE MISS: Fetching heavy resources for ", input$gene_name))
        
        incProgress(0.05, detail = "Downloading UniProt sequence...")
        rv$cache_seq <- fetch_uniprot_seq(input$uniprot_id)
        
        incProgress(0.10, detail = "Running evolutionary conservation (MSA)...")
        rv$cache_cons <- get_evolutionary_profile_shiny(input$uniprot_id, input$gene_name, logger, progress_wrapper)
        
        incProgress(0.15, detail = "Calculating solvent accessibility (SASA)...")
        rv$cache_acc <- get_accessibility_shiny(input$uniprot_id, logger, progress_wrapper)
        
        incProgress(0.10, detail = "Fetching UniProt domain structure...")
        rv$cache_struct <- tryCatch({ require_pkg("drawProteins"); drawProteins::feature_to_dataframe(drawProteins::get_features(input$uniprot_id)) }, error=function(e) { logger(paste0("!!! [ResPred] drawProteins fetch error: ", e$message)); NULL })
        
        rv$cache_key <- current_key
        logger(">>> CACHE WARMED.")
      } else {
        incProgress(0.40, detail = "Using cached data...")
        logger(">>> CACHE HIT: Reusing cached sequence, MSA, SASA, and structure.")
      }
      
      # --- SCANSITE LOGIC ---
      if (input$suite_nav == "ScanSite - Upstream/Downstream") {
        req(input$bind_sites)
        s_vec <- as.numeric(unlist(strsplit(input$bind_sites, ",")))
        s_vec <- s_vec[!is.na(s_vec)]
        
        incProgress(0.2, "Fetching Sequence...")
        seq_str <- rv$cache_seq
        
        if(!is.null(seq_str)) {
          incProgress(0.5, "Querying Scansite...")
          raw_data <- get_scansite_predictions(input$uniprot_id, seq_str, logger)
          
          if(!is.null(raw_data)) {
            rv$bind_report <- report_site_status(s_vec, seq_str, raw_data)
            rv$bind_data <- raw_data %>% filter(Position %in% s_vec)
            logger(paste0(">>> Scansite Analysis Complete. Found ", nrow(rv$bind_data), " hits."))
          }
        } else { logger("!!! UniProt Fetch Failed.") }
        
      # --- MULTI-SITE LOGIC ---
      } else if (input$suite_nav == "Multi-Site Analysis") {
        req(input$target_sites)
        sites <- parse_sites(input$target_sites)
        if (length(sites) > 0) {
          raw_multi <- get_tkl_shiny(input$uniprot_id, input$gene_name, rv$cache_acc, rv$cache_cons, logger, function(x){}, specific_positions = sites)
          if (!is.null(raw_multi)) {
            all_kinases <- unique(raw_multi$Kinase)
            grid <- expand.grid(Kinase = all_kinases, Position = sites, stringsAsFactors = FALSE)
            rv$multisite_data <- grid %>% left_join(raw_multi, by = c("Kinase", "Position")) %>% mutate(Score = replace_na(Score, 0))
          }
        }
        
      # --- KINOME LOGIC ---
      } else if (input$suite_nav == "Protein Kinome") {
        rv$selected_res <- NULL
        incProgress(0.1, detail = "Fetching PDB...")
        rv$pdb_path <- download_pdb_r_global(input$uniprot_id)
        incProgress(0.3, detail = "Calculating Conservation...")
        rv$cons_df <- rv$cache_cons
        incProgress(0.6, detail = "Running TKL...")
        raw_tkl <- get_tkl_shiny(input$uniprot_id, input$gene_name, NULL, NULL, logger, function(x){})
        
        if (!is.null(raw_tkl)) {
          incProgress(0.9, detail = "Mapping...")
          raw_tkl$Family <- assign_kinase_family(raw_tkl$Kinase, rv$kin_fam_db)
          raw_tkl$Position <- as.integer(raw_tkl$Position)
          if(!is.null(rv$cons_df)) {
            rv$cons_df$Conservation[is.na(rv$cons_df$Conservation)] <- 0
            raw_tkl <- raw_tkl %>% left_join(rv$cons_df, by="Position")
          } else { raw_tkl$Conservation <- 0 }
          
          rv$heatmap_data <- raw_tkl %>% mutate(Site_Label = paste0(Residue, Position)) %>% select(Site_Label, Position, Kinase, Family, Score, Conservation)
          rv$render_trigger <- 1
        }
        
      # --- STANDARD LOGIC (Kinase Prediction / 14-3-3) ---
      } else {
        rv$acc_df <- rv$cache_acc
        rv$cons_df <- rv$cache_cons
        rv$struct_df <- rv$cache_struct
        if(is.null(rv$struct_df)) rv$struct_df <- data.frame(type="CHAIN", begin=1, end=2000)
        
        if(input$suite_nav == "Kinase Prediction") {
          if(grepl("TKL", input$mode_select)) {
            incProgress(0.10, detail = "Running kinase predictions...")
            rv$pred_data <- get_tkl_shiny(input$uniprot_id, input$gene_name, rv$acc_df, rv$cons_df, logger, progress_wrapper)
          } else { rv$pred_data <- NULL }
          
          if(grepl("PSP", input$mode_select) || grepl("Combined", input$mode_select)) {
            incProgress(0.05, detail = "Querying PhosphoSitePlus...")
            rv$known_data <- get_psp_shiny(input$uniprot_id, input$gene_name, rv$acc_df, rv$cons_df, logger, progress_wrapper)
          } else { rv$known_data <- NULL }
          
        } else { # 14-3-3 Mode
          incProgress(0.10, detail = "Querying 14-3-3 prediction server...")
          raw_1433 <- get_1433_data(input$uniprot_id, logger)
          if(!is.null(raw_1433)) {
            raw_1433$DSSP <- NA_real_; raw_1433$Conservation <- NA_real_; raw_1433$Kinase <- "14-3-3"
            if(!is.null(rv$acc_df)) raw_1433$DSSP <- rv$acc_df$DSSP[match(raw_1433$Position, rv$acc_df$Position)]
            if(!is.null(rv$cons_df)) raw_1433$Conservation <- rv$cons_df$Conservation[match(raw_1433$Position, rv$cons_df$Position)]
            rv$pred_data <- raw_1433
          } else { rv$pred_data <- NULL }
          rv$known_data <- NULL
        }
        
        domains <- subset(rv$struct_df, type %in% c("DOMAIN", "REGION", "MOTIF", "REPEAT"))
        if(nrow(domains) > 0) domains <- domains[order(domains$begin), ]
        rv$domain_choices <- if(nrow(domains) > 0) unique(paste0(domains$description, " [", domains$begin, "-", domains$end, "]")) else NULL
        if(!is.null(rv$domain_choices)) rv$domain_choices <- unique(rv$domain_choices)
        
        show_psp <- (input$suite_nav == "Kinase Prediction" && (grepl("PSP", input$mode_select) || grepl("Combined", input$mode_select)))
        
        showModal(modalDialog(
          title = "Visualization Settings",
          h4("Select Domains"),
          div(class = "d-flex gap-2 mb-2", actionButton("sel_all_dom", "Select All", class = "btn-sm btn-outline-primary"), actionButton("sel_none_dom", "Deselect All", class = "btn-sm btn-outline-secondary")),
          if(!is.null(rv$domain_choices)) checkboxGroupInput("selected_domains", NULL, choices = rv$domain_choices, selected = rv$domain_choices) else p("No domains."),
          hr(),
          if(show_psp) tagList(h4("PSP Evidence"), checkboxGroupInput("psp_evidence", NULL, choices = c("LT_LIT", "MS_LIT", "MS_CST"), selected = c("LT_LIT", "MS_LIT", "MS_CST"))) else hidden(checkboxGroupInput("psp_evidence", "", choices = "LT_LIT", selected = "LT_LIT")),
          footer = tagList(modalButton("Cancel"), actionButton("render_final_btn", "Generate Plot", class = "btn-success"))
        ))
      }
      incProgress(0.25)
      logger(">>> DATA ACQUISITION COMPLETE.")
    })
  })
  
  observeEvent(input$sel_all_dom, { updateCheckboxGroupInput(session, "selected_domains", selected = rv$domain_choices) })
  observeEvent(input$sel_none_dom, { updateCheckboxGroupInput(session, "selected_domains", selected = character(0)) })
  
  observeEvent(input$render_final_btn, {
    removeModal()
    logger("Rendering Plot...")
    target_k <- NULL; mode_flag <- NULL
    
    if (input$suite_nav == "Kinase Prediction") {
      final_psp <- NULL
      if(!is.null(rv$known_data)) final_psp <- get_psp_shiny(input$uniprot_id, input$gene_name, rv$acc_df, rv$cons_df, logger, function(x){}, input$psp_evidence)
      target_k <- input$tkl_kinase
      mode_flag <- input$mode_select
      rv$last_render_module <- "Kinase Prediction"
      rv$last_render_args <- list(pred_data = rv$pred_data, known_data = final_psp,
                                   struct_df = rv$struct_df, uniprot_id = input$uniprot_id,
                                   query_val = target_k, mode_type = mode_flag,
                                   selected_domains = input$selected_domains)
      rv$final_plot <- plot_master_shiny(rv$pred_data, final_psp, rv$struct_df, input$uniprot_id, target_k, mode_flag, input$selected_domains,
                                          filter_score_pct = input$kin_filter_score %||% 0,
                                          filter_cons_pct  = input$kin_filter_cons  %||% 0,
                                          filter_acc_pct   = input$kin_filter_acc   %||% 0)
    } else {
      rv$last_render_module <- "14-3-3 Prediction"
      rv$last_render_args <- list(pred_data = rv$pred_data, known_data = NULL,
                                   struct_df = rv$struct_df, uniprot_id = input$uniprot_id,
                                   query_val = "14-3-3", mode_type = "1433",
                                   selected_domains = input$selected_domains)
      rv$final_plot <- plot_master_shiny(rv$pred_data, NULL, rv$struct_df, input$uniprot_id, "14-3-3", "1433", input$selected_domains,
                                          filter_score_pct = input$f1433_filter_score %||% 0,
                                          filter_cons_pct  = input$f1433_filter_cons  %||% 0,
                                          filter_acc_pct   = input$f1433_filter_acc   %||% 0)
    }
    
    if(is.null(rv$final_plot)) { logger("!!! Plot Failed: No data."); showNotification("Plot Failed: No data found.", type="error") }
    else { logger(">>> PLOT SUCCESS.") }
  })
  
  # --- FILTER APPLY HANDLERS ---
  observeEvent(input$kin_apply_filters, {
    req(rv$last_render_args)
    if (!identical(rv$last_render_module, "Kinase Prediction")) return()
    logger(paste0(">>> Applying Kinase filters: Score>=", input$kin_filter_score,
                   "%, Cons>=", input$kin_filter_cons, "%, Acc>=", input$kin_filter_acc, "%"))
    args <- rv$last_render_args
    rv$final_plot <- plot_master_shiny(args$pred_data, args$known_data, args$struct_df,
                                        args$uniprot_id, args$query_val, args$mode_type,
                                        args$selected_domains,
                                        filter_score_pct = input$kin_filter_score,
                                        filter_cons_pct  = input$kin_filter_cons,
                                        filter_acc_pct   = input$kin_filter_acc)
    if (is.null(rv$final_plot)) showNotification("No data remaining after filtering.", type = "warning")
  })
  
  observeEvent(input$f1433_apply_filters, {
    req(rv$last_render_args)
    if (!identical(rv$last_render_module, "14-3-3 Prediction")) return()
    logger(paste0(">>> Applying 14-3-3 filters: Score>=", input$f1433_filter_score,
                   "%, Cons>=", input$f1433_filter_cons, "%, Acc>=", input$f1433_filter_acc, "%"))
    args <- rv$last_render_args
    rv$final_plot <- plot_master_shiny(args$pred_data, args$known_data, args$struct_df,
                                        args$uniprot_id, args$query_val, args$mode_type,
                                        args$selected_domains,
                                        filter_score_pct = input$f1433_filter_score,
                                        filter_cons_pct  = input$f1433_filter_cons,
                                        filter_acc_pct   = input$f1433_filter_acc)
    if (is.null(rv$final_plot)) showNotification("No data remaining after filtering.", type = "warning")
  })

  # Main plot outputs
  output$main_plot_kin  <- renderPlotly({ req(rv$final_plot); rv$final_plot })
  output$main_plot_1433 <- renderPlotly({ req(rv$final_plot); rv$final_plot })
  
  # ScanSite outputs
  output$bind_plot_down <- renderPlotly({
    req(rv$bind_data)
    plot_data <- rv$bind_data %>% filter(Category == "Downstream Binder") %>% mutate(Strength = -log10(Percentile/100))
    if(nrow(plot_data) == 0) {
      return(plot_ly() %>% layout(annotations = list(list(
        text = "No downstream predictions found", showarrow = FALSE,
        xref = "paper", yref = "paper", x = 0.5, y = 0.5, font = list(size = 14)
      ))))
    }
    plot_ly(data = plot_data, x = ~Strength, y = ~reorder(Motif_Name, Strength),
            type = 'bar', color = ~as.factor(Position), orientation = 'h') %>%
      layout(barmode = 'group',
             title = list(text = paste0("Downstream Binding Domains<br><sup>", ATTR_SCANSITE, "</sup>")),
             margin = list(t = 80),
             xaxis = list(title = "Strength (-log10 Percentile)"),
             yaxis = list(title = "Binding Partner"),
             legend = list(title = list(text = "Residue")))
  })
  
  output$bind_plot_up <- renderPlotly({
    req(rv$bind_data)
    plot_data <- rv$bind_data %>% filter(Category == "Upstream Kinase") %>% mutate(Strength = -log10(Percentile/100))
    if(nrow(plot_data) == 0) {
      return(plot_ly() %>% layout(annotations = list(list(
        text = "No upstream kinase predictions found", showarrow = FALSE,
        xref = "paper", yref = "paper", x = 0.5, y = 0.5, font = list(size = 14)
      ))))
    }
    plot_ly(data = plot_data, x = ~Strength, y = ~reorder(Motif_Name, Strength),
            type = 'bar', color = ~as.factor(Position), orientation = 'h') %>%
      layout(barmode = 'group',
             title = list(text = paste0("Upstream Kinases<br><sup>", ATTR_SCANSITE, "</sup>")),
             margin = list(t = 80),
             xaxis = list(title = "Strength (-log10 Percentile)"),
             yaxis = list(title = "Kinase"),
             legend = list(title = list(text = "Residue")))
  })
  
  output$bind_report    <- renderText({ rv$bind_report })
  output$table_raw_bind <- DT::renderDT({ req(rv$bind_data); rv$bind_data })
  
  # Multi-Site outputs
  output$multisite_table <- DT::renderDT({
    req(rv$multisite_data)
    wide_df <- rv$multisite_data %>% select(Kinase, Position, Score) %>% pivot_wider(names_from = Position, values_from = Score, names_prefix = "Site_")
    site_cols <- grep("Site_", names(wide_df))
    mat_vals <- as.matrix(wide_df[, site_cols])
    wide_df$Average <- rowMeans(mat_vals, na.rm = TRUE)
    wide_df$SD <- apply(mat_vals, 1, sd, na.rm = TRUE)
    wide_df %>% arrange(desc(Average), SD) %>% mutate(across(where(is.numeric), \(x) round(x, 3)))
  }, options = list(pageLength = 15))
  
  output$multisite_plot <- renderPlotly({
    req(rv$multisite_data)
    agg_df <- rv$multisite_data %>% group_by(Kinase) %>% summarise(Avg = mean(Score), SD = sd(Score), SE = sd(Score) / sqrt(n()), .groups = "drop") %>% arrange(desc(Avg), SD) %>% slice_head(n = 20)
    ordered_kinases <- rev(agg_df$Kinase)
    agg_df$Kinase <- factor(agg_df$Kinase, levels = ordered_kinases)
    points_df <- rv$multisite_data %>% filter(Kinase %in% ordered_kinases) %>% mutate(Kinase = factor(Kinase, levels = ordered_kinases))
    plot_ly(data = agg_df, x = ~Avg, y = ~Kinase, type = 'bar', orientation = 'h',
            marker = list(color = 'steelblue', opacity = 0.5),
            error_x = list(type = 'data', array = agg_df$SE, color = 'black'),
            name = 'Average', showlegend = FALSE) %>%
      add_markers(data = points_df, x = ~Score, y = ~Kinase,
                  color = ~as.factor(Position), opacity = 0.8,
                  marker = list(size = 8),
                  name = ~as.factor(Position)) %>%
      layout(title = list(text = paste0("Top 20 Predicted Kinases for ", input$gene_name,
                                        "<br><sup>Ranked by Highest Average & Lowest Standard Deviation</sup>")),
             xaxis = list(title = "Prediction Score"),
             yaxis = list(title = "Kinase", categoryorder = "array", categoryarray = ordered_kinases),
             legend = list(title = list(text = "Residue"), orientation = "h", y = -0.15),
             barmode = 'overlay')
  })
  
  # Kinome heatmap
  output$heatmap_plotly <- renderPlotly({
    req(rv$heatmap_data)
    rv$plotly_ready <- FALSE
    min_score <- input$heatmap_score_filter
    min_cons  <- input$cons_threshold
    plot_data <- rv$heatmap_data %>% filter(Score >= min_score) %>% filter(Conservation >= min_cons)
    if(nrow(plot_data) < 2) return(NULL)
    
    if (input$heatmap_sort == "struct" && !is.null(rv$pdb_path)) {
      unique_pos <- unique(plot_data$Position)
      sorted_pos <- get_3d_sort_order(rv$pdb_path, unique_pos)
      if (!is.null(sorted_pos)) {
        plot_data_uniq <- plot_data %>% select(Position, Site_Label) %>% distinct()
        ordered_labels <- plot_data_uniq$Site_Label[match(sorted_pos, plot_data_uniq$Position)]
        ordered_labels <- ordered_labels[!is.na(ordered_labels)]
        plot_data$Site_Label <- factor(plot_data$Site_Label, levels = ordered_labels)
      } else {
        plot_data$Site_Label <- reorder(plot_data$Site_Label, plot_data$Position)
      }
    } else {
      plot_data$Site_Label <- reorder(plot_data$Site_Label, plot_data$Position)
    }
    
    n_rows <- length(unique(plot_data$Site_Label))
    calc_height <- max(600, n_rows * 15)
    ordered_kinases <- plot_data %>% select(Family, Kinase) %>% distinct() %>% arrange(Family, Kinase)
    gap_size <- 2; current_x <- 1; last_fam <- ordered_kinases$Family[1]
    x_map <- data.frame(Kinase = character(), X = numeric(), stringsAsFactors = FALSE)
    for(i in 1:nrow(ordered_kinases)) {
      fam <- ordered_kinases$Family[i]
      if(fam != last_fam) { current_x <- current_x + gap_size; last_fam <- fam }
      x_map <- rbind(x_map, data.frame(Kinase = ordered_kinases$Kinase[i], X = current_x))
      current_x <- current_x + 1
    }
    plot_data <- plot_data %>% left_join(x_map, by="Kinase")
    fam_centers <- plot_data %>% group_by(Family) %>% summarise(Center = mean(X))
    y_max_rows <- length(unique(plot_data$Site_Label))
    
    p <- ggplot(plot_data, aes(x = X, y = Site_Label, fill = Score, customdata = Position, text = paste("Family:", Family, "<br>Kinase:", Kinase, "<br>Site:", Site_Label))) +
      geom_tile(color = "white", linewidth = 0.1) +
      annotate("text", x = fam_centers$Center, y = y_max_rows + 2, label = fam_centers$Family, size = 3.5, fontface = "bold", vjust = 0) +
      scale_fill_gradientn(colors = c("#f7fcf5", "#41b6c4", "#081d58"), values = scales::rescale(c(0, 0.5, 1))) +
      scale_x_continuous(breaks = x_map$X, labels = x_map$Kinase) +
      theme_minimal() + labs(y = NULL, x = NULL) + coord_cartesian(clip = "off") + expand_limits(y = y_max_rows + 5) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), axis.text.y = element_text(size = 9), panel.grid = element_blank(), plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"))
    
    plotly_obj <- ggplotly(p, height = calc_height, tooltip = "text", source = "A") %>% event_register("plotly_click") %>% layout(autosize = TRUE, margin = list(t = 180, b = 100))
    isolate(rv$plotly_ready <- TRUE)
    plotly_obj
  })
  
  observeEvent({
    req(isTRUE(rv$plotly_ready))
    event_data("plotly_click", source = "A")
  }, {
    click_data <- event_data("plotly_click", source = "A")
    if(!is.null(click_data)) {
      res_num <- as.character(click_data$customdata)
      if(length(res_num) > 0 && res_num != "") {
        rv$selected_res <- res_num
        rv$render_trigger <- rv$render_trigger + 1
        logger(paste("Click Detected. Re-rendering for residue:", res_num))
      }
    }
  })
  
  observeEvent(input$highlight_btn, {
    req(input$manual_residues)
    res_vec <- unlist(strsplit(input$manual_residues, ","))
    cleaned <- as.numeric(trimws(res_vec))
    cleaned <- cleaned[!is.na(cleaned)]
    if(length(cleaned) > 0) {
      sel_str <- paste(cleaned, collapse = " or ")
      rv$selected_res <- sel_str
      rv$render_trigger <- rv$render_trigger + 1
      logger(paste("Manual Input. Re-rendering for:", sel_str))
    }
  })
  
  output$structure <- NGLVieweR::renderNGLVieweR({
    req(rv$pdb_path, rv$render_trigger)
    require_pkg("NGLVieweR")
    viewer <- NGLVieweR::NGLVieweR(rv$pdb_path) %>%
      NGLVieweR::addRepresentation("cartoon", param = list(name = "protein", color = "lightgray")) %>%
      NGLVieweR::stageParameters(backgroundColor = "white") %>%
      NGLVieweR::setQuality("high")
    if (!is.null(rv$selected_res)) {
      viewer <- viewer %>% NGLVieweR::addRepresentation("spacefill", param = list(name = "highlight", sele = rv$selected_res, color = "red", radiusScale = 2))
    }
    viewer
  })
  
  # Download handlers
  # output$download_plot removed: Plotly's built-in ModeBar camera icon handles export

  output$download_multisite_csv <- downloadHandler(
    filename = function() { paste0(input$gene_name, "_MultiSite_Analysis.csv") },
    content = function(file) {
      req(rv$multisite_data)
      wide_df <- rv$multisite_data %>% select(Kinase, Position, Score) %>% pivot_wider(names_from = Position, values_from = Score, names_prefix = "Pos_")
      site_cols <- grep("Pos_", names(wide_df))
      mat_vals <- as.matrix(wide_df[, site_cols])
      wide_df$Average <- rowMeans(mat_vals, na.rm = TRUE)
      wide_df$SD <- apply(mat_vals, 1, sd, na.rm = TRUE)
      write.csv(wide_df %>% arrange(desc(Average), SD), file, row.names = FALSE)
    }
  )
}
