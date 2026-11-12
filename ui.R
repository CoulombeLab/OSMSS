# ui.R (OSMSS) - Three-Panel Layout
# [Nav Sidebar] | [Input Sidebar] | [Content Area]
# ==============================================================================

ui <- page_sidebar(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  useShinyjs(),
  # Force uppercase on UniProt ID fields so case mismatches never corrupt lookups
  tags$script(HTML(
    "$(document).on('input', '#uniprot_id, #ms_uniprot_id', function() {",
    "  var pos = this.selectionStart;",
    "  this.value = this.value.toUpperCase();",
    "  this.setSelectionRange(pos, pos);",
    "});"
  )),
  #### SPINNING LOADING Screen ####
  tags$div(
    id = "osmss_loading_overlay",
    tags$div(class = "text-center",
      tags$div(class = "spinner-border text-primary mb-3", role = "status"),
      tags$h5(class = "text-muted", "Initializing OSMSS..."),
      tags$p(class = "small text-muted", "Loading modules. This may take a moment.")
    )
  ),
  # EXTENDED TASK: Async task blocking overlay — shown while ExtendedTask runs in background
  tags$div(
    id = "async_task_overlay",
    style = "display:none;",
    class = "async-overlay",
    tags$div(
      class = "async-overlay-content",
      tags$div(class = "spinner-border text-primary mb-3", role = "status",
               tags$span(class = "visually-hidden", "Running...")),
      tags$h4(id = "async_task_label", "Processing..."),
      tags$p(id = "async_task_detail", class = "text-muted", ""),
      actionButton("async_cancel_btn", "Cancel",
                   class = "btn btn-outline-danger mt-3",
                   icon = icon("xmark"))
    )
  ),
  tags$head(tags$style(HTML("
    #console_log { background-color: #1e1e1e; color: #00ff00; font-family: monospace;
                   font-size: 12px; height: 500px; overflow-y: scroll; padding: 10px; }
    .nav-stacked > li > a { padding: 6px 12px; font-size: 0.9em; }
    /* Custom collapsible sidebar menu */
    .nav-menu-item { display: block; padding: 7px 14px; color: #2c3e50; text-decoration: none;
                     font-size: 0.9em; border-radius: 4px; cursor: pointer; }
    .nav-menu-item:hover { background-color: #e8f6f3; color: #18bc9c; text-decoration: none; }
    .nav-menu-item.active { background-color: #18bc9c; color: white; }
    .nav-menu-item.disabled { color: #adb5bd; cursor: not-allowed; opacity: 0.6; }
    .nav-menu-item.disabled:hover { background-color: transparent; color: #adb5bd; }
    /* Featured nav item — Saved Experiments */
    .nav-menu-featured { background-color: #e5f2f7; border-left: 3px solid #18bc9c;
                         padding-left: 11px; font-weight: 600; font-size: 1.2em; margin-bottom: 15px;
                         color: #2c3e50; border: 1px solid; border-color: #e0e0e0;}
    .nav-menu-featured.active { background-color: #143d4c; color: white; border-left-color: #12a083; }
    /* Group item — Collapsible */
    .nav-group-btn { display: block; width: 100%; padding: 7px 14px; color: #2c3e50;
                     background: #e5f2f7; border: 1px solid #ddd; border-radius: 4px;
                     font-size: 1.2em; font-weight: 600; cursor: pointer; text-align: left;
                     margin-bottom: 15px; }
    .nav-group-btn:hover { background-color: #d5f5e3; }
    .nav-group-btn .caret-icon { float: right; transition: transform 0.2s; }
    .nav-group-btn.open .caret-icon { transform: rotate(90deg); }
    .nav-group-btn.open { margin-bottom: 2px;}
    .nav-submenu { display: none; padding-left: 12px; margin-bottom: 4px; }
    .nav-submenu.open { display: block; }
    .nav-submenu .nav-submenu { padding-left: 10px; }
    .nav-submenu .nav-group-btn { font-size: 0.85em; padding: 5px 10px; font-weight: 500; margin-bottom: 2px;
                                   background: #f5f7f8; border-color: #e0e0e0; }                                                         
    /* External Resources modal overlay */
    .ext-modal-backdrop { display: none; position: fixed; inset: 0; background: rgba(0,0,0,0.35);
                          z-index: 1060; justify-content: center; align-items: center; }
    .ext-modal-backdrop.open { display: flex; }
    .ext-modal-panel { background: #fff; border-radius: 10px; box-shadow: 0 8px 32px rgba(0,0,0,0.2);
                       padding: 20px 24px 16px; min-width: 280px; max-width: 360px; position: relative; }
    .ext-modal-panel h5 { margin: 0 0 12px; font-size: 1em; font-weight: 700; color: #2c3e50; }
    .ext-modal-panel .ext-close { position: absolute; top: 10px; right: 14px;
                                  background: none; border: none; font-size: 1.2em;
                                  color: #adb5bd; cursor: pointer; line-height: 1; }
    .ext-modal-panel .ext-close:hover { color: #e74c3c; }
    .ext-modal-panel a { display: block; padding: 7px 10px; color: #2c3e50;
                         text-decoration: none; font-size: 0.9em; border-radius: 4px; }
    .ext-modal-panel a:hover { background: #e8f6f3; color: #18bc9c; }
    .ext-modal-panel a i { width: 20px; text-align: center; margin-right: 8px; color: #6c757d; }
    /* Collapsible step sections */
    .step-section { border: 1px solid #dee2e6; border-radius: 6px; margin-bottom: 8px;
                    background: #fff; }
    .step-section > summary.step-header { padding: 8px 12px; font-weight: 600; font-size: 0.9em;
                                          cursor: pointer; list-style: none; display: flex;
                                          align-items: center; gap: 8px; background: #f8f9fa;
                                          border-radius: 6px; user-select: none; }
    .step-section[open] > summary.step-header { border-bottom: 1px solid #dee2e6;
                                                 border-radius: 6px 6px 0 0; }
    .step-section > summary.step-header::-webkit-details-marker { display: none; }
    .step-section > summary.step-header::before { content: ''; display: inline-block;
                                                   width: 0; height: 0; flex-shrink: 0;
                                                   border-left: 6px solid #6c757d;
                                                   border-top: 5px solid transparent;
                                                   border-bottom: 5px solid transparent;
                                                   transition: transform 0.2s; }
    .step-section[open] > summary.step-header::before { transform: rotate(90deg); }
    .step-badge { display: inline-flex; align-items: center; justify-content: center;
                  width: 22px; height: 22px; border-radius: 50%; background: #18bc9c;
                  color: white; font-size: 0.75em; font-weight: 700; flex-shrink: 0; }
    .step-badge-yellow { display: inline-flex; align-items: center; justify-content: center;
                  width: 22px; height: 22px; border-radius: 50%; background: #f39c12;
                  color: white; font-size: 0.65em; font-weight: 700; flex-shrink: 0; }
    .step-body { padding: 10px 12px 8px; }
    /* Ensure input sidebar scrolls instead of clipping */
    .bslib-page-sidebar > .bslib-sidebar-layout > .sidebar > .sidebar-content,
    .bslib-sidebar-layout > .sidebar > .sidebar-content {
      overflow-y: auto !important;
      display: block !important;
    }
    .bslib-sidebar-layout > .sidebar {
      overflow-y: auto !important;
    }
    /* Compact import-summarization file-input to match export button */
    .import-summ-compact .form-group { margin-bottom: 4px; }
    .import-summ-compact .input-group .form-control { display: none; }
    .import-summ-compact .input-group-btn,
    .import-summ-compact .input-group-prepend { width: 100%; display: block; }
    .import-summ-compact .progress { display: none; }
    .import-summ-compact .btn-file {
      width: 100%; border-radius: 0.25rem !important;
      font-size: 0.875rem; padding: 0.25rem 0.5rem;
      background: transparent; color: #3498db;
      border: 1px solid #3498db; font-weight: 400; text-align: center;
    }
    .import-summ-compact .btn-file:hover,
    .import-summ-compact .btn-file:focus {
      background: #3498db; color: #fff;
    }
    /* Full-page loading overlay — shown on startup, hidden by server after init */
    #osmss_loading_overlay {
      position: fixed; top: 0; left: 0; width: 100%; height: 100%;
      background: #ffffff; z-index: 9999;
      display: flex; align-items: center; justify-content: center;
      flex-direction: column;
    }
    #osmss_loading_overlay .spinner-border {
      width: 3rem; height: 3rem;
    }
    /* Async ExtendedTask blocking overlay */
    .async-overlay {
      position: fixed;
      top: 0; left: 0; right: 0; bottom: 0;
      background: rgba(255,255,255,0.85);
      z-index: 10000;
      display: flex;
      align-items: center;
      justify-content: center;
      text-align: center;
    }
    .async-overlay-content {
      padding: 2rem;
      border-radius: 12px;
      background: white;
      box-shadow: 0 4px 24px rgba(0,0,0,0.15);
      min-width: 280px;
    }
  "))),
  tags$head(tags$script(HTML("
    $(document).on('click', '.nav-group-btn:not(.disabled)', function() {
      var target = $(this).data('target');
      $(this).toggleClass('open');
      $('#' + target).toggleClass('open');
    });
    $(document).on('click', '.nav-menu-item:not(.disabled)', function() {
      var tab = $(this).data('tab');
      if (tab) {
        $('.nav-menu-item').removeClass('active');
        $(this).addClass('active');
        Shiny.setInputValue('suite_nav', tab, {priority: 'event'});
      }
    });
    $(document).on('shiny:connected', function() {
      Shiny.setInputValue('suite_nav', 'Home', {priority: 'event'});
    });
    // External Resources modal
    $(document).on('click', '#ext_resources_btn', function(e) {
      e.preventDefault();
      $('#ext_resources_modal').addClass('open');
    });
    $(document).on('click', '.ext-modal-backdrop, .ext-close', function(e) {
      if (e.target === this) $('#ext_resources_modal').removeClass('open');
    });
    // Fix fileInput binding inside BS5 modals
    $(document).on('shown.bs.modal', '.modal', function() {
      Shiny.bindAll($(this));
    });
    // Wizard group pickers: search + Select All (see checkbox_filter_controls)
    function cbxFilterRows(wrap) {
      return $('#' + $(wrap).attr('data-target') + ' .checkbox');
    }
    function cbxFilterApply(wrap) {
      var term = ($(wrap).find('input[type=text]').val() || '').toLowerCase();
      cbxFilterRows(wrap).each(function() {
        var lbl = $(this).find('label span').text().toLowerCase();
        this.style.display = (term === '' || lbl.indexOf(term) !== -1) ? '' : 'none';
      });
    }
    $(document).on('input change', '.cbx-filter input[type=text]', function() {
      cbxFilterApply($(this).closest('.cbx-filter'));
    });
    // A re-rendered checklist starts unfiltered; re-apply the current term.
    $(document).on('shiny:bound', function(e) {
      if (!e.target || !e.target.id) return;
      $('.cbx-filter').filter(function() {
        return $(this).attr('data-target') === e.target.id;
      }).each(function() { cbxFilterApply(this); });
    });
    $(document).on('click', '.cbx-filter-select-all', function() {
      var boxes = cbxFilterRows($(this).closest('.cbx-filter'))
        .filter(function() { return this.style.display !== 'none'; })
        .find('input[type=checkbox]');
      boxes.prop('checked', true);
      boxes.first().trigger('change');
    });
  "))),

title = tags$a(
    "One Stop MS Shop (OSMSS)", 
    href = "#", 
    
    # --- The Simplified JavaScript Action ---
    onclick = "
      Shiny.setInputValue('suite_nav', 'Home');  
      document.querySelectorAll('.nav-menu-featured.active').forEach(function(el) {
        el.classList.remove('active');
      });
    ", 
    
    # --- The Visual Styling ---
    style = "color: inherit; text-decoration: none; font-size: 1.3em; font-weight: 600;" 
  ),

  # =========================================================================
  # LEFT SIDEBAR: Navigation only (narrow)
  # =========================================================================
  sidebar = sidebar(
    width = 275, open = "desktop",
    
    # Saved Experiments
    tags$a(class = "nav-menu-item nav-menu-featured", `data-tab` = "Saved Experiments",
           icon("folder-open"), " Saved Experiments"),

    # --- Analyze MS Data ---
    tags$button(class = "nav-group-btn", `data-target` = "submenu_discoverer",
                icon("microscope"), " Analyze MS Data",
                tags$span(class = "caret-icon", icon("caret-right"))),
    tags$div(id = "submenu_discoverer", class = "nav-submenu",

      # -- Protein sub-group --
      tags$button(class = "nav-group-btn", `data-target` = "submenu_protein",
                  icon("circle-nodes"), " Protein",
                  tags$span(class = "caret-icon", icon("caret-right"))),
      tags$div(id = "submenu_protein", class = "nav-submenu",

        # Upload Data
        tags$button(class = "nav-group-btn", `data-target` = "submenu_upload_data",
                    icon("upload"), " Upload Data",
                    tags$span(class = "caret-icon", icon("caret-right"))),
        tags$div(id = "submenu_upload_data", class = "nav-submenu",
          tags$a(class = "nav-menu-item", `data-tab` = "Upload Data",
                 icon("layer-group"), " Label Free DDA"),
          tags$a(class = "nav-menu-item", `data-tab` = "DIA Upload Data",
                 icon("atom"), " Label Free DIA"),
          tags$a(class = "nav-menu-item", `data-tab` = "TMT Upload Data",
                 icon("flask-vial"), " TMT")
        ),

        # Explore Data
        tags$button(class = "nav-group-btn", `data-target` = "submenu_explore_data",
                    icon("magnifying-glass-chart"), " Explore Data",
                    tags$span(class = "caret-icon", icon("caret-right"))),
        tags$div(id = "submenu_explore_data", class = "nav-submenu",
          tags$a(class = "nav-menu-item", `data-tab` = "PCA Analysis",
                 icon("braille"), " PCA Analysis"),
          tags$a(class = "nav-menu-item", `data-tab` = "Differential Expression",
                 icon("chart-area"), " Differential Expression"),
          tags$a(class = "nav-menu-item", `data-tab` = "Trend Clusters",
                 icon("chart-line"), " Trend Clusters"),
          tags$a(class = "nav-menu-item", `data-tab` = "Heatmap",
                 icon("grip"), " Heatmap"),
          tags$a(class = "nav-menu-item", `data-tab` = "GSEA Analysis",
                 icon("chart-line"), " GSEA Analysis"),
          tags$a(class = "nav-menu-item", `data-tab` = "Protein Coverage",
                 icon("map"), " Protein Coverage"),
          tags$a(class = "nav-menu-item", `data-tab` = "Cross-Experiment Comparison",
                 icon("circle-nodes"), " Cross-Experiment")
        )
      ),

      # -- PTM Analysis sub-group --
      tags$button(class = "nav-group-btn", `data-target` = "submenu_pd_tmt_ptm",
                  icon("star-of-life"), " PTM",
                  tags$span(class = "caret-icon", icon("caret-right"))),
      tags$div(id = "submenu_pd_tmt_ptm", class = "nav-submenu",

        # Upload sub-group
        tags$button(class = "nav-group-btn", `data-target` = "submenu_ptm_upload",
                    icon("upload"), " Upload",
                    tags$span(class = "caret-icon", icon("caret-right"))),
        tags$div(id = "submenu_ptm_upload", class = "nav-submenu",
          tags$a(class = "nav-menu-item", `data-tab` = "PTM Upload Data",
                 icon("flask-vial"), " TMT-PTM"),
          tags$a(class = "nav-menu-item", `data-tab` = "Label Free PTM Upload",
                 icon("flask"), " Label Free PTM")
        ),

        # Explore sub-group
        tags$button(class = "nav-group-btn", `data-target` = "submenu_ptm_explore",
                    icon("magnifying-glass-chart"), " Explore",
                    tags$span(class = "caret-icon", icon("caret-right"))),
        tags$div(id = "submenu_ptm_explore", class = "nav-submenu",
          tags$a(class = "nav-menu-item", `data-tab` = "PTM PCA",
                 icon("braille"), " PTM PCA"),
          tags$a(class = "nav-menu-item", `data-tab` = "PTM DE",
                 icon("volcano"), " PTM DE"),
          tags$a(class = "nav-menu-item", `data-tab` = "PTM Trends",
                 icon("chart-line"), " PTM Trends"),
          tags$a(class = "nav-menu-item", `data-tab` = "PTM Single",
                 icon("bullseye"), " Single PTM"),
          tags$a(class = "nav-menu-item", `data-tab` = "KS Analysis",
                 icon("diagram-project"), " Kinase-Site Analysis"),
          tags$a(class = "nav-menu-item", `data-tab` = "PTM Cross-Experiment Comparison",
                 icon("circle-nodes"), " Cross-Experiment")
        )
      ),

      # -- Single-Shot LFQ sub-group --
      tags$button(class = "nav-group-btn", `data-target` = "submenu_pd_singleshot",
                  icon("crosshairs"), " N=1 Label Free Experiments",
                  tags$span(class = "caret-icon", icon("caret-right"))),
      tags$div(id = "submenu_pd_singleshot", class = "nav-submenu",
        tags$a(class = "nav-menu-item", `data-tab` = "Rank Fold Change",
               icon("arrow-down-wide-short"), " Rank Fold Change"),
        tags$a(class = "nav-menu-item", `data-tab` = "Trends Analysis",
               icon("chart-line"), " Trends Analysis"),
        tags$a(class = "nav-menu-item", `data-tab` = "Protein Coverage",
               icon("map"), " Protein Coverage"),
        tags$a(class = "nav-menu-item", `data-tab` = "GSEA Analysis",
               icon("chart-line"), " GSEA Analysis")
      )
    ),
    
    # --- PTM Prediction Tools ---    
    tags$button(class = "nav-group-btn", `data-target` = "submenu_kinase",
                icon("chart-line"), " Kinase Prediction",
                tags$span(class = "caret-icon", icon("caret-right"))),
    tags$div(id = "submenu_kinase", class = "nav-submenu",
      tags$a(class = "nav-menu-item", `data-tab` = "Kinase Prediction",
             icon("chart-line"), " Kinase Prediction"),
      tags$a(class = "nav-menu-item", `data-tab` = "Multi-Site Analysis",
             icon("table-cells"), " Multi-Site Analysis"),
      tags$a(class = "nav-menu-item", `data-tab` = "Protein Kinome",
             icon("network-wired"), " Kinome Analysis"),
      tags$a(class = "nav-menu-item", `data-tab` = "ScanSite - Upstream/Downstream",
             icon("magnet"), " Scansite Up/Downstream"),
      tags$a(class = "nav-menu-item", `data-tab` = "14-3-3 Prediction",
             icon("dna"), " 14-3-3 Prediction")
    ),

         # --- Fragpipe PTM Database ---      
    tags$button(class = "nav-group-btn", `data-target` = "submenu_fragpipe",
                icon("flask"), " PTM Database",
                tags$span(class = "caret-icon", icon("caret-right"))),
    tags$div(id = "submenu_fragpipe", class = "nav-submenu",
      tags$a(class = "nav-menu-item", `data-tab` = "Load Database",
             icon("upload"), " Load Database"),
      tags$a(class = "nav-menu-item", `data-tab` = "Manage Database",
             icon("database"), " Manage Database"),
      tags$a(class = "nav-menu-item", `data-tab` = "FragPipe PTM Analysis",
             icon("chart-bar"), " Explore PTMs (Graph)"),
      tags$a(class = "nav-menu-item", `data-tab` = "FragPipe Data Table",
             icon("table"), " Explore PTMs (Table)")
    ),
    
    # --- External Resources
    tags$a(class = "nav-menu-item", id = "ext_resources_btn",
           icon("arrow-up-right-from-square"), " External Resources"),
  ),
  # External Resources modal (rendered outside the sidebar so it floats freely)
  tags$div(id = "ext_resources_modal", class = "ext-modal-backdrop",
    tags$div(class = "ext-modal-panel",
      tags$button(class = "ext-close", HTML("&times;")),
      tags$h5(icon("arrow-up-right-from-square"), " External Resources"),
      tags$a(href = "https://www.uniprot.org/", target = "_blank",
             icon("database"), "UniProt"),
      tags$a(href = "https://www.phosphosite.org/", target = "_blank",
             icon("atom"), "PhosphoSite Plus"),
      tags$a(href = "https://kinase-library.mit.edu/", target = "_blank",
             icon("book"), "The Kinase Library"),
      tags$a(href = "https://scansite4.mit.edu/", target = "_blank",
             icon("magnifying-glass"), "Scansite 4.0"),
      tags$a(href = "https://fragpipe.nesvilab.org/", target = "_blank",
             icon("flask"), "FragPipe"),
      tags$a(href = "https://www.compbio.dundee.ac.uk/1433pred/", target = "_blank",
             icon("dna"), "14-3-3-Pred"),
      tags$a(href = "https://string-db.org/", target = "_blank",
             icon("circle-nodes"), "STRING-db"),
      tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4102138/", target = "_blank",
        icon("flask"), "SAINT Express"),
      tags$a(href = "https://www.bioconductor.org/packages/release/bioc/html/proDA.html", target = "_blank",
        icon("circle-nodes"), "ProDA"),
      tags$a(href = "https://www.proteomexchange.org/", target = "_blank",
             icon("globe"), "ProteomeXchange"),
      tags$a(href = "https://www.coulombelab.com/", target = "_blank",
             icon("flask"), "Coulombe lab"),
      tags$a(href = "mailto:erezco@umich.edu?subject=OSMSS%3A%20Request", target = "_blank",
             icon("envelope"), "Contact Us")             
    )
  ),
  
  # =========================================================================
  # MAIN AREA: Input Sidebar + Content Panels
  # =========================================================================
  layout_sidebar(
    fillable = TRUE, fill = TRUE,
    
    # --- INPUT SIDEBAR (context-sensitive per module) ---
    sidebar = sidebar(
      width = 300, open = "desktop",
      
      # HOME
      conditionalPanel(
        condition = "input.suite_nav == 'Home'",
        h4(icon("house"), " Welcome"),
        p("Select a module from the navigation menu."),
        hr(),
      ),
      
      # =====================================================================
      # ADJACENT MENU: Saved Experiments Controls
      # =====================================================================
      conditionalPanel(
        condition = "input.suite_nav == 'Saved Experiments'",
        h4(icon("folder-open"), " Saved Experiments"),
        
        # Visible styled folder button (Stacked, full width)
        tags$button(
          id    = "prerun_folder_btn_js",
          class = "btn btn-primary w-100 mb-2",
          icon("folder-open"), " Select Folder"
        ),
        
        # Refresh Button (Stacked, full width)
        actionButton("prerun_refresh_btn",
                     tagList(icon("rotate"), " Refresh"),
                     class = "btn-outline-secondary w-100 mb-3",
                     title = "Re-process the last uploaded batch"),
        
        # Status Output
        uiOutput("prerun_folder_status", class = "mb-3"),
        
        # Instructional Hint
        tags$div(class = "alert alert-info py-2 px-3 small mb-2",
                 icon("circle-info"), " Select the folder where you extracted ",
                 "your experiment ZIPs. The app scans all sub-folders for ",
                 tags$code(".preview.rds"), " files and shows them as tiles. ",
                 "Clicking ", tags$b("Load"), " uploads only that experiment’s ",
                 "full ", tags$code(".rds"), " file."
        ),
        
        # --- HIDDEN LOGIC ELEMENTS ---
        tags$input(
          id               = "prerun_folder_trigger",
          type             = "file",
          webkitdirectory  = "",
          multiple         = "",
          style            = "display:none;"
        ),
        tags$div(
          id = "prerun_file_input_wrapper",
          fileInput("prerun_files", label = NULL, multiple = TRUE, accept = ".rds")
        ),
        tags$div(
          id = "prerun_auto_rds_wrapper",
          fileInput("prerun_auto_rds_input", label = NULL, multiple = FALSE, accept = ".rds")
        )
      ),

      # FRAGPIPE DATABASE LOAD DATABASE
      conditionalPanel(
        condition = "input.suite_nav == 'Load Database'",
        h4(icon("upload"), " Load Database"),
        helpText(class = "small text-muted",
          "Import an existing database or start a new one.",
          " Use Export in Manage Database to save your work."
        )
      ),
      
      # FRAGPIPE DATABASE MANAGE DATABASE
      conditionalPanel(
        condition = "input.suite_nav == 'Manage Database'",
        h4(icon("database"), " Manage FragPipe DB"),
        actionButton("fpdb_add_exp_btn", "Add Experiment",
                     class = "btn-success w-100 mb-3", icon = icon("plus-circle")),
        hr(),
        helpText(class = "small text-muted",
          "Browse experiments in Card view or Table view.",
          " Click a card to see full details, edit metadata, or upload data."
        ),
        hr(),
        downloadButton("fpdb_export_btn", "Export Database",
                       class = "btn-outline-primary w-100 mb-3",
                       icon = icon("file-export")),
        hr(),
        tags$label(class = "fw-bold small", "Database Summary"),
        verbatimTextOutput("fpdb_db_summary")
      ),
      
      # FRAGPIPE DATABASE PTM ANALYSIS
      conditionalPanel(
        condition = "input.suite_nav == 'FragPipe PTM Analysis'",
        h4(icon("flask"), " Explore PTMs"),

        # Step 1: Select Protein
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Select Protein"),
          tags$div(class = "step-body",
            textInput("ms_gene_name", "Gene Name:", value = "KRT18"),
            textInput("ms_uniprot_id", "UniProt ID:", value = "P05783"),
            tags$div(style = "display:none;",  # Human only in this release
            radioButtons("ms_organism", NULL, inline = TRUE,
              choiceNames = list(
                "Human",
                tags$span(class = "text-muted", "Mouse ", tags$em("(coming soon)"))
              ),
              choiceValues = c("Human", "Mouse"),
              selected = "Human"
            ),
            tags$script(HTML(
              "$(document).ready(function(){
                 $('input[name=ms_organism][value=Mouse]').prop('disabled',true)
                   .closest('label').css('opacity','0.45').css('cursor','not-allowed');
               });"
            )))
          )
        ),

        # FRAGPIPE DATABASE Step 2: View Mode
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "View Mode"),
          tags$div(class = "step-body",
            radioButtons("ms_view_mode", NULL,
                          choices = c("Aggregate Graph" = "count",
                                      "Confidence Graph" = "confidence"),
                          selected = "count", inline = TRUE),
            helpText(class = "small text-muted",
              "Aggregate: Y = experiment count per site. ",
              "Confidence: Y = best confidence score per site.")
          )
        ),

        # FRAGPIPE DATABASE: Generate Button
        actionButton("ms_plot_btn", "Generate Plot",
                     class = "btn-success w-100 my-2", icon = icon("chart-bar")),

        # FRAGPIPE DATABASE: Step 3: Protein Domains
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Protein Domains"),
          tags$div(class = "step-body",
            helpText("Select domains to overlay on the plot:"),
            div(class = "d-flex gap-2 mb-2",
              actionButton("ms_dom_all", "Select All", class = "btn-sm btn-outline-primary"),
              actionButton("ms_dom_none", "Deselect All", class = "btn-sm btn-outline-secondary")
            ),
            checkboxGroupInput("ms_selected_domains", NULL, choices = character(0))
          )
        ),

        # FRAGPIPE DATABASE: Confidence Score Info
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "text-muted", icon("circle-info")), "Confidence Score"),
          tags$div(class = "step-body",
            tags$p(style = "font-size: 0.85em; color: #555;",
              tags$b("Score (0\u2013100)"), " = ",
              "Mod Type Prior (25) + Localization (20) + Score Delta (20) + ",
              "PeptideProphet (15) + Uniqueness (10) + Mass Accuracy (4) + ",
              "Missed Cleavages (3) + Peptide Confirmation (3)"
            ),
            tags$ul(style = "font-size: 0.85em; color: #555; padding-left: 18px;",
              tags$li(tags$span(style="color:green;", "\u25CF"), " High: \u226570"),
              tags$li(tags$span(style="color:orange;", "\u25CF"), " Mixed: 45\u201369"),
              tags$li(tags$span(style="color:red;", "\u25CF"), " Low: <45")
            )
          )
        )
      ),
      
      # FRAGPIPE DATABASE TABLE
      conditionalPanel(
        condition = "input.suite_nav == 'FragPipe Data Table'",
        h4(icon("table"), " Filter Data"),
        textInput("ms_table_gene", "Gene Name:", value = "KRT18"),
        tags$div(style = "display:none;",  # Human only in this release
        radioButtons("ms_table_organism", NULL, inline = TRUE,
          choiceNames = list(
            "Human",
            tags$span(class = "text-muted", "Mouse ", tags$em("(coming soon)"))
          ),
          choiceValues = c("Human", "Mouse"),
          selected = "Human"
        ),
        tags$script(HTML(
          "$(document).ready(function(){
             $('input[name=ms_table_organism][value=Mouse]').prop('disabled',true)
               .closest('label').css('opacity','0.45').css('cursor','not-allowed');
           });"
        ))),
        actionButton("ms_table_btn", "Load Table",
                     class = "btn-info w-100", icon = icon("magnifying-glass")),
        br(), br(),
        helpText("Shows all raw data from FragpipeDB for the selected gene.")
      ),

      # =====================================================================
      # UPLOAD DATA (Centralized upload/pipeline/groups/summarize for LFQ)
      # =====================================================================
      conditionalPanel(
        condition = "input.suite_nav == 'Upload Data'",
        h4(icon("upload"), " Upload Data"),

        #Hidden inputs (kept for conditionalPanel compat in DE/PCA/HM/TCL)
        tags$div(style = "display:none;",
          radioButtons("upload_pipeline", NULL,
            choices  = c("global", "ip"), selected = "global"),
          radioButtons("upload_ip_sub", NULL,
            choices  = c("invariant", "no_norm"), selected = "invariant"),
          checkboxInput("upload_skip_impute", NULL, value = FALSE)
        ),

        # Wizard Trigger (before first completion)
        conditionalPanel(
          condition = "!output.wiz_completed_flag",
          tags$div(class = "alert alert-info py-2 px-3 small mb-2",
            icon("wand-magic-sparkles"),
            " Use the wizard to configure your upload settings."
          ),
          actionButton("wiz_open",
            tagList(icon("wand-magic-sparkles"), " Start Setup Wizard"),
            class = "btn-primary w-100 mb-2")
        ),

        # Post-Wizard Sidebar (after completion)
        conditionalPanel(
          condition = "output.wiz_completed_flag",
          uiOutput("wiz_sidebar_badges"),
          tags$hr(class = "my-2"),
          actionButton("wiz_reopen",
            tagList(icon("rotate"), " Re-open Wizard"),
            class = "btn-outline-secondary btn-sm w-100 mb-2")
        ),

        # Batch Correction status (sidebar, shown after wizard assigns batches)
        conditionalPanel(
          condition = "output.upload_batch_active",
          tags$details(class = "step-section",
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", icon("layer-group")),
              "Batch Correction"),
            tags$div(class = "step-body",
              uiOutput("upload_batch_status"),
              actionButton("upload_batch_clear", "Clear Batches",
                           class = "btn-outline-danger btn-sm w-100 mt-1",
                           icon = icon("xmark"))
            )
          )
        ),
        #Export Summarization (server-rendered based on cache state)
        uiOutput("wiz_export_btn"),
      ),

      # DIFFERENTIAL EXPRESSION
      conditionalPanel(
        condition = "input.suite_nav == 'Differential Expression'",
        h4(icon("chart-area"), " Differential Expression"),

        # Data Status Banner
        uiOutput("de_data_status_banner"),

        # Control / Condition group pickers
        uiOutput("de_ctrl_group_picker"),
        uiOutput("de_cond_group_picker"),

        # Factorial model type (only shown when factorial design is active)
        uiOutput("de_factorial_model_ui"),

        # Step 1: Run Analysis
        shinyjs::disabled(
          actionButton("de_run_btn", "Run DE Analysis",
                       class = "btn-success w-100 my-2", icon = icon("play"))
        ),
        # Perseus LOD toggle -- hidden in this release (not described in the paper); to re-enable,
        # set the condition back to "output.de_is_dia"
        conditionalPanel(
          condition = "false",
          actionButton("de_lod_toggle",
                       "Activate Perseus LOD Mode",
                       class = "btn-secondary w-100 mb-2",
                       icon = icon("wand-magic-sparkles")),
          conditionalPanel(
            condition = "output.de_lod_active",
            tags$div(class = "mt-1 mb-2",
              checkboxInput("de_lod_ack", label = tagList(
                tags$span(class = "small",
                  "I understand that this mode uses a full ",
                  tags$b('"Perseus-equivalent pipeline"'),
                  ", including LOD imputation normally used for DDA data. ",
                  "For new analysis of DIA data \u2014 standard limma mode is recommended. ",
                  "Activating this mode will run ",
                  tags$b("(1)"), " log\u2082 transform + column-wise median centering, ",
                  tags$b("(2)"), " LOD imputation (per-column Gaussian: mean \u2212 1.8\u00d7SD, width 0.3\u00d7SD), ",
                  tags$b("(3)"), " Welch t-test with s0=0.1 fudge factor, ",
                  tags$b("(4)"), " permutation FDR (250 perms, median null count). ",
                  "Covariates and batch correction will be ",
                  tags$b("disabled.")
                )
              ), value = FALSE)
            )
          )
        ),

        # Step 2: Plot Options (Classic DE volcano)
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Plot Options"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted mb-2",
              "Volcano plot: \u2212log10(adj. P-value) on the Y-axis, log2 fold change on the X-axis."),
            div(class = "d-flex align-items-center gap-2 mb-2",
              checkboxInput("de_show_capped", "Show capped values", value = TRUE),
              tooltip(
                tags$span(class = "text-muted", style = "cursor: help; font-size: 0.9em;",
                          icon("circle-info")),
                paste0("PD caps FC at 100 (or 0.01) and P-values at a floor. ",
                       "Capped values create artificial clusters. ",
                       "An AdjPvalue of 0 is treated as capped."),
                placement = "right"
              )
            )
          )
        ),

        # Step 3: Highlight & Filters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Highlight & Filters"),
          tags$div(class = "step-body",
            # Gene search
            textInput("de_gene_search", "Highlight Gene:",
                      placeholder = "e.g. Krt8"),
            checkboxInput("de_exact_match", "Exact match only", value = FALSE),
            actionButton("de_gene_search_btn", "Highlight",
                         class = "btn-warning w-100 mb-2", icon = icon("search")),
            hr(),
            # Highlight Rules
            div(class = "d-flex align-items-center justify-content-between mb-2",
              h6(class = "mb-0", icon("highlighter"), " Highlight Rules"),
              tooltip(
                actionButton("de_add_highlight", NULL,
                             class = "btn-outline-warning btn-sm",
                             icon = icon("plus")),
                "Add highlight rule (column + operator + value)",
                placement = "right"
              )
            ),
            helpText(class = "small text-muted",
              "Highlight proteins matching rules. Multiple rules are OR'd."),
            uiOutput("de_active_highlights"),
            actionButton("de_clear_highlights", "Clear Highlights",
                         class = "btn-outline-warning btn-sm w-100",
                         icon = icon("eraser")),
            hr(),
            # Data Filters
            div(class = "d-flex align-items-center justify-content-between mb-2",
              h6(class = "mb-0", icon("filter"), " Data Filters"),
              tooltip(
                actionButton("de_open_filters", NULL,
                             class = "btn-outline-primary btn-sm",
                             icon = icon("sliders")),
                "Open filter editor",
                placement = "right"
              )
            ),
            helpText(class = "small text-muted",
              "Filter (remove) proteins by any column. ",
              "Conditions within a group use AND or OR; groups are AND'd."),
            uiOutput("de_active_filters"),
            actionButton("de_clear_filters", "Clear All Filters",
                         class = "btn-outline-danger btn-sm w-100",
                         icon = icon("trash"))
          )
        ),

        # Step 4: Export
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Export"),
          tags$div(class = "step-body",
            downloadButton("de_download_csv", "Download Results CSV",
                           class = "btn-primary btn-sm w-100 mb-2"),
            # DIA only: proteins excluded from limma and tested by Fisher (NMAR tab)
            conditionalPanel(
              condition = "output.de_nmar_download_ready",
              downloadButton("de_download_nmar", "Download NMAR Analysis",
                             class = "btn-primary btn-sm w-100 mb-2")
            )
          )
        ),

        # Step 5: GSEA Export
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "5"), "GSEA Export"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted mb-2",
              "Export analysis results in GSEA-compatible formats. ",
              "Available after running DE analysis."),

            # .rnk
            conditionalPanel(
              condition = "output.de_gsea_rnk_ready",
              tags$label(class = "form-label fw-bold small",
                         icon("arrow-down-wide-short"), " Ranked List (.rnk)"),
              helpText(class = "small text-muted",
                "Signed \u03c0-score: sign(Log2FC) \u00d7 |Log2FC| \u00d7 \u2212log10(adj. p) ",
                "(BFDR for SAINT). Use with GSEA Pre-Ranked."),
              downloadButton("de_download_rnk",
                             "Download .rnk",
                             class = "btn-outline-info btn-sm w-100 mb-3"),
            ),

            # .gct + .cls
            conditionalPanel(
              condition = "output.de_gsea_gct_ready",
              tags$label(class = "form-label fw-bold small mt-2",
                         icon("table-cells"), " Expression Matrix (.gct + .cls)"),
              helpText(class = "small text-muted",
                "Post-normalization log2 abundances per replicate. ",
                "Use with standard GSEA (requires both files)."),
              downloadButton("de_download_gct",
                             "Download .gct",
                             class = "btn-outline-info btn-sm w-100 mb-2"),
              downloadButton("de_download_cls",
                             "Download .cls",
                             class = "btn-outline-info btn-sm w-100 mb-2")
            )
          )
        )
      ),  # end DE conditionalPanel
      
      # RANK FOLD CHANGE (Single-Shot)
      conditionalPanel(
        condition = "input.suite_nav == 'Rank Fold Change'",
        h4(icon("arrow-down-wide-short"), " Rank Fold Change"),
        helpText(class = "small text-muted mb-2",
          "Single-shot analysis for experiments without replicates. ",
          "Compare two conditions by abundance-based fold change."),

        # Step 1: Upload File
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Upload File"),
          tags$div(class = "step-body",
            fileInput("rfc_csv_file", "Upload Discoverer CSV",
                      accept = c(".csv"), width = "100%"),
            helpText(class = "small text-muted",
                     "Proteome Discoverer export with individual Abundance columns."),
            uiOutput("rfc_file_status")
          )
        ),

        # Step 2: Choose Abundance Columns
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Choose Abundance Columns"),
          tags$div(class = "step-body",
            tags$div(class = "alert alert-info py-2 px-3 small mb-2",
              icon("circle-info"),
              " Select one Abundance column per condition. The module will ",
              "Log2-transform, compute fold change, and Z-score the data."),
            tags$label(class = "form-label", "Condition A (e.g. Control):"),
            div(class = "d-flex align-items-center gap-2 mb-2",
              actionButton("rfc_pick_col_a", "Select Column...",
                           class = "btn-outline-secondary btn-sm w-100",
                           icon = icon("columns"))
            ),
            uiOutput("rfc_col_a_label"),
            tags$label(class = "form-label mt-2", "Condition B (e.g. Treatment):"),
            div(class = "d-flex align-items-center gap-2 mb-2",
              actionButton("rfc_pick_col_b", "Select Column...",
                           class = "btn-outline-secondary btn-sm w-100",
                           icon = icon("columns"))
            ),
            uiOutput("rfc_col_b_label")
          )
        ),

        # Step 2b: Imputation Option
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", icon("wand-magic-sparkles")),
            "Imputation (Optional)"),
          tags$div(class = "step-body",
            checkboxInput("rfc_impute_toggle",
                          "Enable Perseus-style imputation for missing values",
                          value = FALSE),
            conditionalPanel(
              condition = "input.rfc_impute_toggle == true",
              tags$div(class = "alert alert-info py-2 px-3 small mb-2",
                style = "border-left: 4px solid #3498db;",
                icon("circle-info"),
                " Missing (NA) abundance values will be replaced with ",
                "low-intensity draws sampled from a downshifted normal ",
                "distribution (mean \u2212 1.8 SD, width 0.3 SD), following the ",
                "Perseus imputation method. This reveals proteins detected ",
                "in only one condition as extreme fold changes. ",
                tags$b("No statistical testing is performed."),
                tags$br(), tags$br(),
                tags$span(class = "text-muted",
                  icon("triangle-exclamation"),
                  " Imputed fold changes are estimates. ",
                  "Imputed points will be marked with a distinct shape (triangle) ",
                  "on the plots.")
              ),
              tags$div(class = "small mb-2",
                tags$details(
                  tags$summary(class = "text-primary fw-semibold",
                               style = "cursor: pointer;",
                               icon("list-ol"), " Imputation Pipeline"),
                  tags$ol(class = "mb-1 mt-1", style = "padding-left: 18px;",
                    tags$li("Replace Abundance \u2264 0 with NA"),
                    tags$li("Log2 transform"),
                    tags$li("Per-column Perseus imputation (downshift 1.8 SD)"),
                    tags$li("Compute Log2 FC from imputed values"),
                    tags$li("Flag imputed proteins with distinct visual marker")
                  )
                )
              ),
              numericInput("rfc_impute_downshift", "Downshift (SDs):",
                           value = 1.8, min = 0.5, max = 5, step = 0.1),
              numericInput("rfc_impute_width", "Width (fraction of SD):",
                           value = 0.3, min = 0.05, max = 1, step = 0.05),
              helpText(class = "small text-muted",
                "Defaults match Perseus: downshift = 1.8 SD, width = 0.3 SD. ",
                "Adjust only if you have a specific reason.")
            )
          )
        ),

        # Generate Button
        actionButton("rfc_plot_btn", "Generate Rank FC Plot",
                     class = "btn-success w-100 my-2",
                     icon = icon("arrow-down-wide-short")),

        # Step 3: Viewing Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Viewing Parameters"),
          tags$div(class = "step-body",
            # Highlight
            h6(icon("highlighter"), " Highlight"),
            textInput("rfc_gene_search", "Highlight Gene:",
                      placeholder = "e.g. Krt8"),
            checkboxInput("rfc_exact_match", "Exact match only", value = FALSE),
            actionButton("rfc_gene_search_btn", "Highlight",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("search")),
            hr(),
            # Label Top N
            h6(icon("tags"), " Label Top Fold Changes"),
            numericInput("rfc_top_n", "Number of top proteins:",
                         value = 20, min = 1, max = 200, step = 1),
            actionButton("rfc_label_top_btn", "Toggle Top Labels",
                         class = "btn-info w-100 mb-2",
                         icon = icon("tags")),
            hr(),
            # Gutter width
            h6(icon("arrows-left-right"), " FC Threshold (gutter)"),
            numericInput("rfc_gutter_width", "Log2 FC threshold:",
                         value = 1, min = 0, max = 10, step = 0.25),
            helpText(class = "small text-muted",
                     "Points outside \u00b1 this threshold are colored."),
            hr(),
            # Data Filters
            div(class = "d-flex align-items-center justify-content-between mb-2",
              h6(class = "mb-0", icon("filter"), " Data Filters"),
              tooltip(
                actionButton("rfc_open_filters", NULL,
                             class = "btn-outline-primary btn-sm",
                             icon = icon("sliders")),
                "Open filter editor",
                placement = "right"
              )
            ),
            uiOutput("rfc_active_filters"),
            actionButton("rfc_clear_filters", "Clear All Filters",
                         class = "btn-outline-danger btn-sm w-100",
                         icon = icon("trash"))
          )
        ),

        # Step 4: Download Results
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Download Results"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted mb-2",
              "Download the Rank FC result table as CSV. ",
              "Includes proteins, Log2FC, rank, and imputation flags."),
            conditionalPanel(
              condition = "output.rfc_csv_ready",
              downloadButton("rfc_download_csv",
                             "Download Result CSV",
                             class = "btn-outline-success btn-sm w-100 mb-2")
            )
          )
        ),

        # Step 5: GSEA Export
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "5"), "GSEA Export"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted mb-2",
              "Export ranked gene list for GSEA Pre-Ranked analysis. ",
              "Available after generating the Rank FC plot."),
            conditionalPanel(
              condition = "output.rfc_gsea_rnk_ready",
              tags$label(class = "form-label fw-bold small",
                         icon("arrow-down-wide-short"), " Ranked List (.rnk)"),
              helpText(class = "small text-muted",
                "Genes ranked by Log2 fold change (B vs A). ",
                "Use with GSEA Pre-Ranked."),
              downloadButton("rfc_download_rnk",
                             "Download .rnk",
                             class = "btn-outline-info btn-sm w-100 mb-2")
            )
          )
        )
      ),

      # TRENDS ANALYSIS (Single-Shot)
      conditionalPanel(
        condition = "input.suite_nav == 'Trends Analysis'",
        h4(icon("chart-line"), " Trends Analysis"),
        helpText(class = "small text-muted mb-2",
          "Cluster proteins by expression pattern across conditions. ",
          "Identifies groups with similar trends (up, down, flat, etc.)."),

        # Step 1: Upload File
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Upload File"),
          tags$div(class = "step-body",
            fileInput("trd_csv_file", "Upload Discoverer CSV",
                      accept = c(".csv"), width = "100%"),
            helpText(class = "small text-muted",
                     "Proteome Discoverer export with individual Abundance columns."),
            uiOutput("trd_file_status")
          )
        ),

        # Step 2: Choose Abundance Columns
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Choose Abundance Columns"),
          tags$div(class = "step-body",
            tags$div(class = "alert alert-info py-2 px-3 small mb-2",
              icon("circle-info"),
              " Select all Abundance columns representing your conditions or ",
              "timepoints. Each column = one point on the X-axis."),
            actionButton("trd_pick_cols", "Select Columns...",
                         class = "btn-outline-secondary btn-sm w-100",
                         icon = icon("list-check")),
            uiOutput("trd_col_label")
          )
        ),

        # Step 3: Clustering Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Clustering Parameters"),
          tags$div(class = "step-body",
            numericInput("trd_n_clusters", "Number of clusters (k):",
                         value = 6, min = 2, max = 20, step = 1),
            helpText(class = "small text-muted",
                     "k=6 is typically a good starting point for biological data.")
          )
        ),

        # Generate Button
        actionButton("trd_plot_btn", "Generate Trends Analysis",
                     class = "btn-success w-100 my-2",
                     icon = icon("chart-line")),

        # Step 4: Viewing Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Viewing Parameters"),
          tags$div(class = "step-body",
            h6(icon("highlighter"), " Highlight Gene"),
            textInput("trd_gene_search", "Gene name:",
                      placeholder = "e.g. Krt8"),
            checkboxInput("trd_exact_match", "Exact match only", value = FALSE),
            actionButton("trd_gene_search_btn", "Highlight",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("search")),
            helpText(class = "small text-muted",
                     "Highlights the protein's trend line in orange across all clusters.")
          )
        )
      ),

      # TREND CLUSTERS (Replicated Quantitation)
      conditionalPanel(
        condition = "input.suite_nav == 'Trend Clusters'",
        h4(icon("chart-line"), " Trend Clusters"),

        uiOutput("tcl_data_status_banner"),

        # Step 1: Condition Order
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Condition Order"),
          tags$div(class = "step-body",
            tags$div(class = "alert alert-info py-2 px-3 small mb-2",
              icon("circle-info"),
              " Set the X-axis order (left to right). Typically chronological ",
              "or dose-increasing."),
            uiOutput("tcl_cond_order_label"),
            actionButton("tcl_reorder_btn", "Reorder Conditions",
                         class = "btn-outline-secondary btn-sm w-100",
                         icon = icon("arrows-up-down"))
          )
        ),

        # Step 2: Missing Data Handling
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Missing Data Handling"),
          tags$div(class = "step-body",
            radioButtons("tcl_missing_mode", NULL,
              choices  = c("Omit incomplete proteins" = "omit",
                           "LOD imputation"           = "lod"),
              selected = "omit", inline = FALSE),
            helpText(class = "small text-muted mb-2",
              tags$b("Omit:"), " Drops any protein missing a condition (strictest). ",
              tags$br(),
              tags$b("LOD:"), " Imputes missing values at the 1st-percentile of ",
              "the log2 protein matrix.")
          )
        ),

        # Step 3: Clustering Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Clustering Parameters"),
          tags$div(class = "step-body",
            numericInput("tcl_n_clusters", "Number of clusters (k):",
                         value = 6, min = 2, max = 20, step = 1),
            helpText(class = "small text-muted",
                     "k=6 is typically a good starting point for biological data.")
          )
        ),

        # Generate Button
        actionButton("tcl_plot_btn", "Generate Trend Clusters",
                     class = "btn-success w-100 my-2",
                     icon = icon("chart-line")),

        # Step 4: Viewing Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Viewing Parameters"),
          tags$div(class = "step-body",
            sliderInput("tcl_membership", "Membership threshold:",
                        min = 0, max = 1, value = 0.5, step = 0.05),
            helpText(class = "small text-muted",
                     "Filter proteins by fuzzy membership score. ",
                     "0 = show all, 1 = only core members. ",
                     "Changes apply immediately after analysis."),
            h6(icon("highlighter"), " Highlight Gene"),
            textInput("tcl_gene_search", "Gene name:",
                      placeholder = "e.g. Krt8"),
            checkboxInput("tcl_exact_match", "Exact match only", value = FALSE),
            actionButton("tcl_gene_search_btn", "Highlight",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("search")),
            helpText(class = "small text-muted",
                     "Highlights the protein's trend line in orange across all clusters.")
          )
        )
      ),

      # ═══════════════════════════════════════════════════════════
      # TMT UPLOAD DATA
      # ═══════════════════════════════════════════════════════════
      conditionalPanel(
        condition = "input.suite_nav == 'TMT Upload Data'",
        h4(icon("upload"), " TMT Upload Data"),

        # -- Pre-wizard: info prompt + Start Wizard button --
        conditionalPanel(
          condition = "!output.tmt_wiz_completed_flag",
          tags$div(class = "alert alert-info py-2 px-3 small mb-2",
            icon("circle-info"),
            " Use the Setup Wizard to configure your TMT experiment. ",
            "Select your software, data type, upload files, and define your ",
            "experimental design."),
          actionButton("tmt_wiz_open", "Start Setup Wizard",
            class = "btn-primary w-100 mb-2",
            icon = icon("hat-wizard"))
        ),

        # -- Post-wizard: summary badges + Re-open + Export --
        conditionalPanel(
          condition = "output.tmt_wiz_completed_flag",
          uiOutput("tmt_wiz_sidebar_badges"),
          actionButton("tmt_wiz_reopen", "Re-open Wizard",
            class = "btn-outline-secondary btn-sm w-100 mb-2",
            icon = icon("pen-to-square")),
          downloadButton("tmt_upload_export_summarization",
            "Export Summarization (.rds)",
            class = "btn-outline-info btn-sm w-100 mb-2")
        )
      ),

      # ═══════════════════════════════════════════════════════════
      # PTM UPLOAD DATA
      # ─────────────────────────────────────────────────────────────────────
      conditionalPanel(
        condition = "input.suite_nav == 'PTM Upload Data'",
        h4(icon("flask-vial"), " TMT PTM Upload"),

        # Pre-wizard: show Open Wizard button
        conditionalPanel(
          condition = "!output.tmt_ptm_wiz_completed_flag",
          tags$div(class = "alert alert-info py-2 px-3 small mb-2",
            icon("circle-info"),
            " Upload and configure your TMT phosphoproteomics data using the wizard."),
          actionButton("tmt_ptm_wiz_open", "Open PTM Upload Wizard",
                       class = "btn-primary w-100 mt-1",
                       icon = icon("wand-magic-sparkles"))
        ),

        # Post-wizard: badges + re-open + export
        conditionalPanel(
          condition = "output.tmt_ptm_wiz_completed_flag",
          uiOutput("tmt_ptm_wiz_sidebar_badges"),
          actionButton("tmt_ptm_wiz_reopen", "Re-open Wizard",
                       class = "btn-outline-primary w-100 mt-1 mb-1",
                       icon = icon("rotate")),
          downloadButton("tmt_ptm_upload_export_btn",
                         "Export Summarization (.zip)",
                         class = "btn-outline-success w-100 mb-2",
                         icon = icon("download"))
        )
      ),

      # ─────────────────────────────────────────────────────────────────────
      # PTM PCA
      conditionalPanel(
        condition = "input.suite_nav == 'PTM PCA'",
        h4(icon("braille"), " PTM PCA"),

        uiOutput("ptm_pca_data_status"),

        # ── Step 1: Filter Settings ──
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Filter Settings"),
          tags$div(class = "step-body",
            sliderInput("ptm_pca_valid_pct",
                        "Min. % valid values per site:",
                        min = 0, max = 100, value = 50, step = 5),
            helpText(class = "small text-muted",
              "PTM sites (or proteins) with fewer valid values than this ",
              "threshold are excluded before PCA. Remaining NAs in surviving ",
              "rows are replaced with the row-median.")
          )
        ),

        # ── Run ──
        shinyjs::disabled(
          actionButton("ptm_pca_run_btn", "Run PTM PCA",
                       class = "btn-secondary w-100 my-2",
                       icon = icon("play"))
        ),

        hr(),
        verbatimTextOutput("ptm_pca_status")
      ),

      conditionalPanel(
        condition = "input.suite_nav == 'PTM DE'",
        h4(icon("star-of-life"), " PTM Analysis - DE"),

        uiOutput("ptm_de_data_status_banner"),

        # Step 1: Select Comparison
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Select Comparison"),
          tags$div(class = "step-body",
            selectInput("ptm_group_b", "Control Group:",
                        choices = NULL, width = "100%"),
            selectInput("ptm_group_a", "Condition Group:",
                        choices = NULL, width = "100%")                        
          )
        ),

        # Run Button
        shinyjs::disabled(
          actionButton("ptm_run_btn", "Run PTM DE",
                       class = "btn-secondary w-100 my-2", icon = icon("play"))
        ),

        # Step 2: Viewing Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Viewing Parameters"),
          tags$div(class = "step-body",
            h6(icon("layer-group"), " Analysis View"),
            radioButtons("ptm_view_mode", NULL,
              choices = c(
                "Adjusted to Total Protein" = "adjusted",
                "Absolute PTM"              = "unadjusted",
                "Protein View"              = "protein"
              ),
              selected = "adjusted", inline = FALSE
            ),
            # Views offered and their help follow the loaded cache (PTM-only: Absolute PTM alone)
            uiOutput("ptm_view_mode_help"),
            hr(),
            numericInput("ptm_fc_thresh", "Log2 FC threshold:",
                         value = 1, min = 0, max = 10, step = 0.25),
            numericInput("ptm_p_thresh", "Adj. p-value threshold:",
                         value = 0.05, min = 0.001, max = 1, step = 0.01),
            numericInput("ptm_n_labels", "Label top N significant:",
                         value = 15, min = 0, max = 100, step = 1),
            checkboxInput("ptm_show_heavy_imputed",
                          tags$span(icon("triangle-exclamation"),
                                    " Show heavily imputed sites (LOD, triangles)"),
                          value = TRUE),
            helpText(class = "small text-muted",
                     "Sites with >75% missing values in any group are drawn with values imputed near the ",
                     "detection limit (Perseus-style), as triangles; colored by significance direction. ",
                     "Their fold changes rest on imputed values and should be read with care."),
            hr(),
            h6(icon("highlighter"), " Highlight"),
            textInput("ptm_gene_search", "Gene name:",
                      placeholder = "e.g. VIM"),
            actionButton("ptm_gene_search_btn", "Highlight Gene",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("search")),
            hr(),
            textInput("ptm_ptm_search", "PTM (Gene_Site):",
                      placeholder = "e.g. Krt16_S44"),
            actionButton("ptm_ptm_search_btn", "Highlight PTM",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("crosshairs")),
            checkboxInput("ptm_exact_match", "Exact match only", value = FALSE),
            hr(),
            h6(icon("highlighter"), " Highlight Rules"),
            actionButton("ptm_de_add_highlight", "Add Highlight Rule",
                         class = "btn-outline-warning w-100 mb-2",
                         icon = icon("plus")),
            actionButton("ptm_de_clear_highlights", "Clear All Highlights",
                         class = "btn-outline-secondary btn-sm w-100 mb-2",
                         icon = icon("eraser")),
            uiOutput("ptm_de_active_highlights"),
            hr(),
            h6(icon("filter"), " Data Filters"),
            actionButton("ptm_de_open_filters", "Add Data Filter",
                         class = "btn-outline-info w-100 mb-2",
                         icon = icon("plus")),
            actionButton("ptm_de_clear_filters", "Clear All Filters",
                         class = "btn-outline-secondary btn-sm w-100 mb-2",
                         icon = icon("eraser")),
            uiOutput("ptm_de_active_filters")
          )
        ),

        # Step 3: Export & Save
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Export & Save"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted",
              "Download PTM DE results or save/load full analysis."),
            h6(icon("database"), " Save / Load Analysis"),
            fileInput("ptm_load_rds", NULL,
                      accept = ".rds", width = "100%",
                      buttonLabel = list(icon("upload"), " Load Analysis"),
                      placeholder = "No file selected")
          )
        ),
        downloadButton("ptm_download_csv", "Download Results CSV",
                       class = "btn-primary w-100 mb-2"),
        conditionalPanel(
          condition = "output.ptm_fc_data_ready",
          downloadButton("ptm_download_fc_csv", "Download FC x FC CSV",
                         class = "btn-outline-primary w-100 mb-2")
        ),
        downloadButton("ptm_save_rds", "Save Analysis (.rds)",
                       class = "btn-outline-success btn-sm w-100"),

        # -- Step 4: Export Data for Post-Analysis --
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Export Data for Post-Analysis"),
          tags$div(class = "step-body",
            tags$h6(icon("flask"), " PhosR Export"),
            helpText(class = "small text-muted mb-1",
              "Per-sample quantification matrix with flanking sequences for PhosR ",
              "kinase-substrate scoring. Multi-phospho peptides are expanded to ",
              "individual site rows. Requires a PTM upload (DIA or TMT)."),
            uiOutput("ptm_phosr_export_ui"),
            hr(),
            tags$h6(icon("chart-bar"), " KSEA Export"),
            helpText(class = "small text-muted mb-1",
              "DE results formatted for KSEA: linear FC (2^log2FC), ",
              "semicolon-separated residues, flanking peptide from local proteome. ",
              "Requires PTM DE results."),
            uiOutput("ptm_ksea_export_ui")
          )
        )
      ),

      # PTM TRENDS
      conditionalPanel(
        condition = "input.suite_nav == 'PTM Trends'",
        h4(icon("chart-line"), " PTM Trend Clusters"),

        uiOutput("ptrnd_data_status_banner"),

        # Step 1: Site Eligibility
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Site Eligibility"),
          tags$div(class = "step-body",
            radioButtons("ptrnd_stringency", label = NULL, inline = FALSE,
              choiceNames = list(
                tags$span("Standard ",
                  tags$span(class = "text-muted small",
                    "(\u226550% of replicates valid in every condition)")),
                tags$span("Strict ",
                  tags$span(class = "text-muted small",
                    "(all replicates valid in every condition)"))
              ),
              choiceValues = list("standard", "strict"),
              selected = "standard"
            ),
            helpText(class = "small text-muted",
              "Nothing is imputed. Clustering needs a complete condition profile, so a ",
              "site must meet the chosen threshold in ", tags$b("every"), " condition.",
              tags$br(),
              tags$b("Standard:"), " at least half that condition's replicates (minimum 2, ",
              "capped at the replicate count).", tags$br(),
              tags$b("Strict:"), " every replicate.", tags$br(),
              "Sites that fail are listed in the Excluded Sites tab, not silently dropped.")
          )
        ),

        # Step 2: Condition Order
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Condition Order"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted",
              "Set the X-axis order for trend profiles (e.g. chronological)."),
            uiOutput("ptm_trends_cond_label"),
            actionButton("ptm_trends_reorder_btn", "Reorder Conditions",
                         class = "btn-outline-secondary btn-sm w-100 mb-2",
                         icon = icon("arrows-up-down"))
          )
        ),

        # Step 3: Clustering Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Clustering Parameters"),
          tags$div(class = "step-body",
            numericInput("ptm_trends_n_clusters", "Number of clusters (k):",
                         value = 6, min = 2, max = 20, step = 1),
            helpText(class = "small text-muted",
                     "k=6 is typically a good starting point."),
            checkboxInput("ptrnd_fuzz_auto",
                          "Estimate fuzziness from the data (recommended)",
                          value = TRUE),
            sliderInput("ptrnd_fuzziness", "Fuzziness (m), manual:",
                        min = 1.1, max = 5.0, value = 2.0, step = 0.1),
            helpText(class = "small text-muted",
                     "Controls cluster overlap. m = 2 is a legacy default from ",
                     "low-dimensional clustering, not a standard for omics data: ",
                     "with many sites and few conditions it drives memberships ",
                     "toward 1/k and under-resolves structure. With estimation on, ",
                     "m is computed from the shape of the clustered matrix by the ",
                     "method of Schwaemmle and Jensen (Bioinformatics 2010), as ",
                     "implemented in Mfuzz::mestimate; the value actually used is ",
                     "written to the log and to the methods text. Untick to force ",
                     "the manual value above. Higher values = softer boundaries."),
            hr(),
            radioButtons("ptrnd_site_gate", "Sites to cluster:", inline = FALSE,
              choiceNames = list(
                tags$span("All eligible sites"),
                tags$span("DE-significant only ",
                  tags$span(class = "text-muted small", "(adj.p < 0.05)"))
              ),
              choiceValues = list("none", "de_sig"),
              selected = "none"
            ),
            helpText(class = "small text-muted",
                     "Profiles are z-scored per site, so an unchanged site is drawn with ",
                     "the same amplitude as a strongly responding one. Restricting to ",
                     "sites significant in at least one PTM DE contrast removes that. ",
                     "Requires a completed PTM DE run on this dataset; if none exists, or ",
                     "if fewer sites qualify than there are clusters, all eligible sites ",
                     "are clustered and a notice is shown.")
          )
        ),

        # Run Button
        shinyjs::disabled(
          actionButton("ptm_trends_run_btn", "Run PTM Trend Analysis",
                       class = "btn-secondary w-100 my-2",
                       icon = icon("chart-line"))
        ),

        # Step 4: Viewing Parameters
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Viewing Parameters"),
          tags$div(class = "step-body",
            sliderInput("ptm_trends_membership", "Membership threshold:",
                        min = 0, max = 1, value = 0.7, step = 0.05),
            helpText(class = "small text-muted",
                     "Filter PTM sites by fuzzy membership score. ",
                     "0 = show all, 1 = only core members."),
            h6(icon("highlighter"), " Highlight"),
            textInput("ptm_trends_gene_search", "Gene name:",
                      placeholder = "e.g. Krt16"),
            actionButton("ptm_trends_gene_search_btn", "Highlight Gene",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("search")),
            hr(),
            textInput("ptm_trends_ptm_search", "PTM (Gene_Site):",
                      placeholder = "e.g. Krt16_S44"),
            actionButton("ptm_trends_ptm_search_btn", "Highlight PTM",
                         class = "btn-warning w-100 mb-2",
                         icon = icon("crosshairs")),
            checkboxInput("ptm_trends_exact_match", "Exact match only", value = FALSE),
            helpText(class = "small text-muted",
                     "Highlights the PTM trend line(s) in orange across all clusters.")
          )
        ),

        # Step 5: Export & Save
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "5"), "Export & Save"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted",
              "Download results or save/load the full analysis."),
            h6(icon("database"), " Save / Load Analysis"),
            fileInput("ptm_trends_load_rds", NULL,
                      accept = ".rds", width = "100%",
                      buttonLabel = list(icon("upload"), " Load Analysis"),
                      placeholder = "No file selected")
          )
        ),
        # Main / Total PTM CSV buttons; labels follow the loaded cache (PTM-only: absolute, no Total)
        uiOutput("ptm_trends_download_ui"),
        downloadButton("ptm_trends_download_excluded", "Download Excluded Sites CSV",
                       class = "btn-outline-secondary w-100 mb-2"),
        downloadButton("ptm_trends_save_rds", "Save Analysis (.rds)",
                       class = "btn-outline-success btn-sm w-100")
      ),

      # SINGLE PTM
      conditionalPanel(
        condition = "input.suite_nav == 'PTM Single'",
        h4(icon("bullseye"), " Single PTM Inspector"),

        uiOutput("psg_data_status_banner"),

        # Step 1: Condition Order
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Condition Order"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted",
              "Set the X-axis order for plots (e.g. chronological)."),
            uiOutput("psg_cond_label"),
            actionButton("psg_reorder_btn", "Reorder Conditions",
                         class = "btn-outline-secondary btn-sm w-100 mb-2",
                         icon = icon("arrows-up-down"))
          )
        ),

        # Step 2: Select PTM & Plot
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Select PTM & Plot"),
          tags$div(class = "step-body",
            conditionalPanel(
              condition = "output.psg_cache_ready",
              selectizeInput("psg_ptm_select", "PTM Site:",
                             choices = NULL, options = list(
                               placeholder = "Select a PTM site...",
                               maxOptions = 5000
                             )),
              selectInput("psg_reference", "Reference Condition:",
                          choices = NULL),
              helpText(class = "small text-muted",
                "Fold changes are computed relative to the reference condition."),
              actionButton("psg_plot_btn", "Plot PTM",
                           class = "btn-success w-100 my-2",
                           icon = icon("chart-bar"))
            ),
            conditionalPanel(
              condition = "!output.psg_cache_ready",
              tags$p(class = "text-muted small fst-italic",
                     icon("triangle-exclamation"),
                     " Run Summarization in PTM Upload Data first to enable PTM selection.")
            )
          )
        ),

        # Step 3: Export
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Export"),
          tags$div(class = "step-body",
            downloadButton("psg_download_csv", "Download Data (CSV)",
                           class = "btn-primary w-100 mb-2",
                           icon = icon("file-csv")),
            downloadButton("psg_save_rds", "Save Analysis (.rds)",
                           class = "btn-outline-success btn-sm w-100")
          )
        ),

        verbatimTextOutput("psg_status")
      ),

      # KS ANALYSIS
      conditionalPanel(
        condition = "input.suite_nav == 'KS Analysis'",
        h4(icon("diagram-project"), " Kinase-Site Analysis"),

        # ── Step 1: Upload ──
        tags$details(class = "step-section", open = TRUE,
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Upload CSV"),
          tags$div(class = "step-body",
            helpText(class = "small text-muted",
              "Upload a PhosR-format CSV exported from the PTM DE module.",
              tags$br(),
              "Required columns: ", tags$code("PhosR_ID"), ", ",
              tags$code("Flanking_Sequence_15mer"), ".",
              tags$br(),
              "Optional DE columns: ", tags$code("log2FC"), ", ",
              tags$code("adj.pvalue"), "."),
            radioButtons("ks_species", label = "Species:",
                         choices  = c("Human" = "human", "Mouse" = "mouse"),
                         selected = "human", inline = TRUE),
            conditionalPanel(
              condition = "input.ks_species == 'mouse'",
              radioButtons("ks_mouse_weights", label = "Motif weights (mouse):",
                           choices  = c("Human (default)" = "human",
                                        "Mouse (PhosR native)" = "mouse"),
                           selected = "human", inline = TRUE),
              p(class = "text-muted", style = "font-size: 0.8em;",
                "Number of input sequences per kinase motif, used to rank and weight motif scores.")
            ),
            fileInput("ks_csv_upload", NULL,
                      accept = c(".csv", "text/csv"),
                      buttonLabel = "Choose CSV...")
          )
        ),

        # ── Step 2: Pre-filter (optional) ──
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Pre-filter (optional)"),
          tags$div(class = "step-body",
            checkboxInput("ks_prefilter_on", "Enable site pre-filter", value = FALSE),
            conditionalPanel(
              condition = "input.ks_prefilter_on",
              helpText(class = "small text-muted mb-1",
                "Select sites to score. Profile scores reflect co-regulation",
                "among ONLY these sites — reproduces the PhosR manuscript approach."),
              selectizeInput("ks_site_filter", "Sites (PhosR IDs):",
                             choices  = NULL,
                             multiple = TRUE,
                             options  = list(
                               placeholder = "Select sites...",
                               maxOptions  = 500,
                               plugins     = list("remove_button")
                             )),
              helpText(class = "small text-muted mb-1 mt-2",
                "Optionally restrict the kinase columns shown."),
              selectizeInput("ks_kinase_filter", "Kinases:",
                             choices  = NULL,
                             multiple = TRUE,
                             options  = list(
                               placeholder = "Select kinases (default: top 60)...",
                               maxOptions  = 400,
                               plugins     = list("remove_button")
                             ))
            ),
            conditionalPanel(
              condition = "!input.ks_prefilter_on",
              tags$div(class = "alert alert-info py-2 px-3 small mb-0",
                icon("circle-info"),
                " When pre-filter is off and log2FC / adj.pvalue columns are",
                " present, sites are split into upregulated and downregulated",
                " (padj\u202f<\u202f0.05) and scored separately.")
            )
          )
        ),

        # ── Step 3: Run ──
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Run"),
          tags$div(class = "step-body",
            checkboxInput("ks_show_all_kinases",
                          "Show all kinases (no top-N limit)",
                          value = FALSE),
            checkboxInput("ks_show_all_sig_sites",
                          "Show all significant sites (no top-50 limit)",
                          value = FALSE),
            helpText(class = "small text-muted mb-2",
              "By default, top 60 kinases and top 50 sites (by |log2FC|) are shown. ",
              "Both options apply to KS and Activity heatmaps."),
            actionButton("ks_run_btn", "Run PhosR Predictions",
                         class = "btn-success w-100 my-2",
                         icon  = icon("play")),
            conditionalPanel(
              condition = "output.ks_has_plex",
              tags$div(class = "d-flex align-items-center gap-2 mt-1",
                tags$span(class = "badge bg-info text-dark",
                  icon("shield-halved"), " TMT detected"),
                tags$span(class = "small text-muted",
                  "Plex batch correction (RUVphospho) will be applied automatically")
              )
            )
          )
        ),

        # ── Step 4: Sample Order (Kinase Activity only) ──
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Sample Order"),
          tags$div(class = "step-body",
            checkboxInput("ks_custom_sample_order",
                          "Custom sample order (Activity chart only)",
                          value = FALSE),
            conditionalPanel(
              condition = "input.ks_custom_sample_order",
              helpText(class = "small text-muted mb-2",
                icon("circle-info"),
                " Drag to reorder samples. Order applies to the Kinase Activity",
                " heatmap only (clustering is disabled). Re-run to apply."),
              selectizeInput("ks_sample_order", NULL,
                             choices  = NULL,
                             multiple = TRUE,
                             options  = list(
                               plugins     = list("drag_drop", "remove_button"),
                               placeholder = "Upload a CSV first..."
                             ))
            )
          )
        ),

        # ── Step 5: Download ──
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "5"), "Download PDFs"),
          tags$div(class = "step-body",
            conditionalPanel(
              condition = "output.ks_has_result_up",
              downloadButton("ks_pdf_up", "KS Heatmap (Up / Pre-filter)",
                             class = "btn-primary w-100 mb-2",
                             icon  = icon("file-pdf"))
            ),
            conditionalPanel(
              condition = "output.ks_has_result_down",
              downloadButton("ks_pdf_down", "KS Heatmap (Down)",
                             class = "btn-primary w-100 mb-2",
                             icon  = icon("file-pdf"))
            ),
            conditionalPanel(
              condition = "output.ks_has_result_up || output.ks_has_result_down",
              downloadButton("ks_pdf_activity", "Kinase Activity Heatmap",
                             class = "btn-outline-primary w-100 mb-1",
                             icon  = icon("file-pdf")),
              tags$p(class = "small text-muted fst-italic mb-2",
                style = "font-size: 0.75rem;",
                icon("circle-info", class = "me-1"),
                "For multi-condition experiments, see the note in the ",
                tags$strong("\u201cKinase Activity\u201d"), " results tab before downloading.")
            ),
            conditionalPanel(
              condition = "output.ks_has_result_up || output.ks_has_result_down",
              downloadButton("ks_pdf_ranking", "Kinase Ranking Chart",
                             class = "btn-outline-primary w-100 mb-2",
                             icon  = icon("file-pdf"))
            ),
            conditionalPanel(
              condition = "!output.ks_has_result_up && !output.ks_has_result_down",
              helpText(class = "small text-muted fst-italic",
                icon("triangle-exclamation"),
                " Run analysis first to enable downloads.")
            )
          )
        ),

        # Status — always visible, outside any collapsible section
        verbatimTextOutput("ks_status_text")
      ),

      # PCA ANALYSIS
      conditionalPanel(
        condition = "input.suite_nav == 'PCA Analysis'",
        h4(icon("braille"), " PCA Analysis"),

        uiOutput("pca_data_status_banner"),


        # Step 1: CV Group Assignment
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "CV Group Assignment"),
          tags$div(class = "step-body",
            helpText(class = "small",
              "Assign samples to groups for the CV analysis tab. ",
              "Groups are set automatically from the condition definitions; ",
              "use this to override them."),
            actionButton("pca_cv_assign_groups", "Assign / Edit Groups...",
                         class = "btn-outline-info btn-sm w-100",
                         icon = icon("users")),
            uiOutput("pca_cv_group_label")
          )
        ),


        # Run Button
        actionButton("pca_run_btn", "Run PCA",
                     class = "btn-success w-100 my-2", icon = icon("play"))
      ),

      # HEATMAP
      conditionalPanel(
        condition = "input.suite_nav == 'Heatmap'",
        h4(icon("grip"), " Heatmap"),

        uiOutput("hm_data_status_banner"),

        # Missing Data
        tags$details(class = "step-section", open = TRUE,
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", icon("circle-question")), " Missing Data"),
          tags$div(class = "step-body",
            radioButtons("hm_missing_mode", label = NULL,
              choices = c(
                "Omit proteins with any missing value" = "omit",
                "LOD imputation (1st percentile per column)" = "lod"
              ),
              selected = "omit"
            )
          )
        ),


        # Generate Button
        actionButton("hm_generate", "Generate Heatmap",
                     class = "btn-success w-100 my-2", icon = icon("fire")),


        # Heatmap Scale
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", icon("ruler-combined")), " Heatmap Scale"),
          tags$div(class = "step-body",
            radioButtons("hm_scale_mode", NULL,
              choices = c(
                "Row Z-Score (relative pattern)" = "zscore",
                "Normalized Abundance (log2)"     = "abundance"
              ),
              selected = "zscore", inline = FALSE),
            helpText(class = "small text-muted mb-1",
              tags$b("Z-Score:"), " Each protein is centered to its own mean and scaled ",
              "by its standard deviation. Colors show how far each sample deviates from ",
              "that protein\u2019s average across the experiment \u2014 emphasizes relative patterns. ",
              "A value of +2 means \u20182 SD above this protein\u2019s mean,\u2019 not \u2018abundance = 2.\u2019",
              tags$br(),
              tags$b("Normalized Abundance:"), " Plots summarized log2 values ",
              "directly. Colors reflect the actual magnitude of protein expression. ",
              "Useful for comparing absolute levels between proteins or conditions."
            )
          )
        ),


        # Step 1: Display Options
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Display Options"),
          tags$div(class = "step-body",
            numericInput("hm_n_top", "Number of Top Genes:",
                         value = 50, min = 5, max = 500, step = 5),
            selectizeInput("hm_gene_filter", "Select Specific Genes:",
                           choices = NULL, multiple = TRUE,
                           options = list(placeholder = "Search genes...",
                                          maxOptions = 10000,
                                          plugins = list("remove_button"))),
            tags$hr(),
            checkboxInput("hm_avg_by_cond", "Average by Condition",
                          value = FALSE),
            helpText(class = "small text-muted mt-n1",
              "Collapse replicates to per-condition means. Variance and top-gene",
              "ranking recomputed on averaged groups. (IP/Global only)"),
            tags$hr(),
            tags$label(class = "form-label fw-bold small",
              icon("arrows-up-down"), " Column Order"),
            uiOutput("hm_col_order_label"),
            actionButton("hm_reorder_cols", "Reorder Columns...",
                         class = "btn-outline-secondary btn-sm w-100 mb-2",
                         icon = icon("arrows-up-down")),
            hr(),
            actionButton("hm_refresh", "Refresh View",
                         class = "btn-outline-info btn-sm w-100",
                         icon = icon("sync"))
          )
        ),


        # Step 2: Save / Load Analysis
        # Download button
        downloadButton("hm_download_csv", "Download Results CSV",
                       class = "btn-primary w-100 mb-2")
      ),

      # PROTEIN COVERAGE
      conditionalPanel(
        condition = "input.suite_nav == 'Protein Coverage'",
        h4(icon("map"), " Protein Coverage"),

        uiOutput("cov_data_status_banner"),


        # Step 1: Configure Protein
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Configure Protein"),
          tags$div(class = "step-body",
            textInput("cov_accession", "Protein Accession:",
                      placeholder = "e.g. P60709"),
            tags$label(class = "form-label mt-2", "Positions Column (peptide-level):"),
            selectInput("cov_pos_col", label = NULL,
                        choices = character(0)),
            helpText(class = "small text-muted",
              "Usually called \"Positions in Master Proteins\". ",
              "Auto-detected from the peptide header row.")
          )
        ),


        # Generate Button
        actionButton("cov_generate_btn", "Generate Coverage Map",
                     class = "btn-success w-100 my-2", icon = icon("map")),


        # Step 2: Domain Overlay
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Domain Overlay"),
          tags$div(class = "step-body",
            uiOutput("cov_domain_selector")
          )
        )
      ),


      # DIA UPLOAD DATA
      conditionalPanel(
        condition = "input.suite_nav == 'DIA Upload Data'",
        h4(icon("upload"), " DIA Upload Data"),

        # Pre-wizard: show open button
        conditionalPanel(
          condition = "!output.dia_wiz_completed_flag",
          tags$div(class = "alert alert-info py-2 px-3 small",
            icon("circle-info"), " Upload and configure your DIA data using the wizard."),
          actionButton("dia_wiz_open", "Open Upload Wizard",
                       class = "btn-primary w-100 mt-2",
                       icon = icon("wand-magic-sparkles"))
        ),

        # Post-wizard: badges + re-open + export
        conditionalPanel(
          condition = "output.dia_wiz_completed_flag",
          uiOutput("dia_wiz_sidebar_badges"),
          actionButton("dia_wiz_reopen", "Re-open Wizard",
                       class = "btn-outline-primary w-100 mt-1 mb-1",
                       icon = icon("rotate")),
          downloadButton("dia_upload_export_summ_btn", "Export Summarization (.rds)",
                         class = "btn-outline-success w-100 mb-2")
        ),

        # Batch correction status
        conditionalPanel(
          condition = "output.dia_batch_active",
          tags$div(class = "alert alert-warning py-2 px-3 small mt-2",
            icon("layer-group"), " Batch correction active.")
        )
      ),


      # LABEL FREE PTM UPLOAD
      conditionalPanel(
        condition = "input.suite_nav == 'Label Free PTM Upload'",
        h4(icon("flask"), " Label Free PTM Upload"),

        # Pre-wizard: show open button
        conditionalPanel(
          condition = "!output.dia_ptm_wiz_completed_flag",
          tags$div(class = "alert alert-info py-2 px-3 small",
            icon("circle-info"),
            " Upload and configure your DIA phosphoproteomics data using the wizard."),
          actionButton("dia_ptm_wiz_open", "Open PTM Upload Wizard",
                       class = "btn-primary w-100 mt-2",
                       icon = icon("wand-magic-sparkles"))
        ),

        # Post-wizard: badges + re-open + export
        conditionalPanel(
          condition = "output.dia_ptm_wiz_completed_flag",
          uiOutput("dia_ptm_wiz_sidebar_badges"),
          actionButton("dia_ptm_wiz_reopen", "Re-open Wizard",
                       class = "btn-outline-primary w-100 mt-1 mb-1",
                       icon = icon("rotate")),
          downloadButton("dia_ptm_upload_export_summ_btn",
                         "Export Summarization (.zip)",
                         class = "btn-outline-success w-100 mb-2")
        )
      ),

      # ================================================================
      # RESIDUE PREDICTION shared inputs (all 6 modules)
      # ================================================================
      conditionalPanel(
        condition = paste0(
          "['Kinase Prediction','14-3-3 Prediction',",
          "'ScanSite - Upstream/Downstream','Multi-Site Analysis','Protein Kinome']",
          ".indexOf(input.suite_nav) !== -1"
        ),
        h4(icon("dna"), " Target Protein"),
        textInput("uniprot_id", "UniProt ID:", value = "Q04695"), 
        textInput("gene_name", "Gene Symbol:", value = "KRT17"), 
        hr(),
        
        # ScanSite-specific
        conditionalPanel(
          condition = "input.suite_nav == 'ScanSite - Upstream/Downstream'",
          textAreaInput("bind_sites", "Target Sites (UniProt Numbering):", 
                        value = "13, 30", rows = 2, placeholder = "e.g. 13, 30"),
          p(class="text-muted", style="font-size: 0.8em;",
            "Note: Use UniProt numbering (includes signal peptide).")
        ),
        
        # Kinase-specific
        conditionalPanel(
          condition = "input.suite_nav == 'Kinase Prediction'",
          selectInput("mode_select", "Type:",
                      choices = c(
                        "PhosphoPlus" = "Single_PSP",
                        "Kinase Library" = "Single_TKL",
                        "PhosphoPlus+Kinase Library Overlay (Log2)" = "Combined_Log2_TKL",
                        "PhosphoPlus+Kinase Library Overlay (Linear)" = "Combined_Linear_TKL"
                      ), selected = "Single_PSP", selectize = FALSE),
          conditionalPanel(
            "input.mode_select.includes('TKL')",
            selectizeInput("tkl_kinase", "Target Kinase (TKL):",
                           choices = NULL,
                           options = list(placeholder = 'Loading TKL...'))
          )
        ),
        
        # Multi-Site-specific
        conditionalPanel(
          condition = "input.suite_nav == 'Multi-Site Analysis'",
          textAreaInput("target_sites", "Target Sites (for Kinase Pred):",
                        value = "18, 23, 44", rows = 2)
        ),
        
        # Kinome-specific
        conditionalPanel(
          condition = "input.suite_nav == 'Protein Kinome'",
          sliderInput("cons_threshold", "Conservation Filter:",
                      min = 0, max = 1, value = 0, step = 0.1),
          radioButtons("heatmap_sort", "Sort Sites By:",
                       choices = c("Sequence Position (N->C)" = "seq",
                                   "3D Structure Clustering" = "struct")),
          sliderInput("heatmap_score_filter", "Min Prediction Score:",
                      min = 0, max = 5, value = 0.3, step = 0.1),
          hr(),
          textInput("manual_residues", "Highlight Sites (e.g. 18, 23):", value = ""),
          actionButton("highlight_btn", "Show on Structure", class="btn-info w-100")
        ),
        
        actionButton("run_btn", "Initialize Analysis",
                     class = "btn-primary w-100", icon = icon("rocket")),
        br(), br(),
        
        # Filters (Kinase Prediction)
        conditionalPanel(
          condition = "input.suite_nav == 'Kinase Prediction'",
          tags$details(class = "step-section",
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", icon("filter")), "Filters"),
            tags$div(class = "step-body",
              tags$p(class = "small text-muted mb-2",
                tags$b("0%"), " = show all | ",
                tags$b("99%"), " = 99th percentile only"),
              sliderInput("kin_filter_score", "Prediction Strength:",
                          min = 0, max = 99, value = 0, step = 1, post = "%"),
              sliderInput("kin_filter_cons", "Conservation:",
                          min = 0, max = 99, value = 0, step = 1, post = "%"),
              sliderInput("kin_filter_acc", "Accessibility:",
                          min = 0, max = 99, value = 0, step = 1, post = "%"),
              actionButton("kin_apply_filters", "Apply Filters",
                           class = "btn-primary btn-sm w-100", icon = icon("filter"))
            )
          )
        ),
        
        # Filters (14-3-3 Prediction)
        conditionalPanel(
          condition = "input.suite_nav == '14-3-3 Prediction'",
          tags$details(class = "step-section",
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", icon("filter")), "Filters"),
            tags$div(class = "step-body",
              tags$p(class = "small text-muted mb-2",
                tags$b("0%"), " = show all | ",
                tags$b("99%"), " = 99th percentile only"),
              sliderInput("f1433_filter_score", "Prediction Strength:",
                          min = 0, max = 99, value = 0, step = 1, post = "%"),
              sliderInput("f1433_filter_cons", "Conservation:",
                          min = 0, max = 99, value = 0, step = 1, post = "%"),
              sliderInput("f1433_filter_acc", "Accessibility:",
                          min = 0, max = 99, value = 0, step = 1, post = "%"),
              actionButton("f1433_apply_filters", "Apply Filters",
                           class = "btn-primary btn-sm w-100", icon = icon("filter"))
            )
          )
        ),
        
        # Downloads
        conditionalPanel(
          "input.suite_nav == 'Multi-Site Analysis'",
          downloadButton("download_multisite_csv", "Download CSV",
                         class = "btn-success w-100")
        ),
        # download_plot button removed: Plotly ModeBar provides built-in export
      ),
      
      # UNDER CONSTRUCTION - no inputs (DIA modules have their own panels above)

      conditionalPanel(
        condition = "input.suite_nav == 'Under Construction'",
        h4(icon("triangle-exclamation"), " Under Construction"),
        p(class = "text-muted", "This module is not yet available.")
      ),

      # GSEA ANALYSIS
      conditionalPanel(
        condition = "input.suite_nav == 'GSEA Analysis'",
        h4(icon("chart-line"), " GSEA Analysis"),

        # Step 1: Input Mode
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Input Data"),
          tags$div(class = "step-body",
            radioButtons("gsea_input_mode", "Analysis Type:",
              choices = c(
                "Pre-Ranked (.rnk)"          = "ranked",
                "Expression Matrix (.gct + .cls)" = "expression"
              ),
              selected = "ranked"
            ),

            # .rnk upload
            conditionalPanel(
              condition = "input.gsea_input_mode == 'ranked'",
              fileInput("gsea_rnk_file", "Upload .rnk File",
                        accept = c(".rnk", ".txt", ".tsv"), width = "100%"),
              helpText(class = "small text-muted",
                "Tab-delimited: Gene \t Score. ",
                "Export from DE module (Step 8) or supply your own."),
              uiOutput("gsea_rnk_status")
            ),

            # .gct + .cls upload
            conditionalPanel(
              condition = "input.gsea_input_mode == 'expression'",
              fileInput("gsea_gct_file", "Upload .gct File",
                        accept = c(".gct"), width = "100%"),
              uiOutput("gsea_gct_status"),
              fileInput("gsea_cls_file", "Upload .cls File",
                        accept = c(".cls"), width = "100%"),
              uiOutput("gsea_cls_status"),
              helpText(class = "small text-muted",
                "GCT v1.2 expression matrix + CLS phenotype file. ",
                "Export both from DE module (Step 8).")
            )
          )
        ),

        # Step 2: Organism
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "2"), "Organism"),
          tags$div(class = "step-body",
            radioButtons("gsea_organism", NULL,
              choices = c(
                "Human"  = "Homo sapiens",
                "Mouse"  = "Mus musculus"
              ),
              selected = "Homo sapiens",
              inline = TRUE
            ),
            # Mouse-specific: ortholog vs native gene sets
            conditionalPanel(
              condition = "input.gsea_organism == 'Mus musculus'",
              radioButtons("gsea_mouse_geneset_type", "Gene Set Source:",
                choices = c(
                  "Human orthologs (recommended)" = "ortholog",
                  "Mouse-native gene sets"        = "native"
                ),
                selected = "ortholog"
              ),
              helpText(class = "small text-muted",
                "Human orthologs maps MSigDB gene sets to mouse via ortholog ",
                "mapping. Mouse-native uses gene sets curated directly for mouse.")
            )
          )
        ),

        # Step 3: Gene Set Collection
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "3"), "Gene Set Collection"),
          tags$div(class = "step-body",
            selectInput("gsea_category", "MSigDB Category:",
              choices = c(
                "H \u2014 Hallmark"                    = "H",
                "C1 \u2014 Positional"                 = "C1",
                "C2 \u2014 Curated Pathways"           = "C2",
                "C3 \u2014 Regulatory Target"          = "C3",
                "C4 \u2014 Computational"              = "C4",
                "C5 \u2014 Gene Ontology"              = "C5",
                "C6 \u2014 Oncogenic Signatures"       = "C6",
                "C7 \u2014 Immunologic Signatures"     = "C7",
                "C8 \u2014 Cell Type Signatures"       = "C8"
              ),
              selected = "H"
            ),
            selectInput("gsea_subcategory", "Subcategory (if applicable):",
              choices = c("(none)" = ""),
              selected = ""
            ),
            # Cache status (rendered by server_gsea.R)
            uiOutput("gsea_cache_info")
          )
        ),

        # Run Button
        actionButton("gsea_run_btn", "Run GSEA",
                     class = "btn-success w-100 my-2", icon = icon("play")),

        # Download
        tags$details(class = "step-section",
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "4"), "Export Results"),
          tags$div(class = "step-body",
            downloadButton("gsea_download_csv", "Download Results CSV",
                           class = "btn-primary w-100 mb-2")
          )
        )
      ),

      # CROSS-EXPERIMENT COMPARISON
      conditionalPanel(
        condition = "input.suite_nav == 'Cross-Experiment Comparison'",

        helpText(class = "small text-muted",
          "Upload 2 or more DE Result CSV files from ",
          tags$b("Protein \u2192 Differential Expression \u2192 Download Results CSV"),
          "."),

        # -- Step 1: Upload files --
        tags$details(class = "step-section", open = TRUE,
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Upload DE CSV Files"),
          tags$div(class = "step-body",
            fileInput("ce_files",
                      label    = NULL,
                      multiple = TRUE,
                      accept   = ".csv",
                      width    = "100%",
                      buttonLabel   = tags$span(icon("upload"), " Browse"),
                      placeholder   = "No files selected"),
            helpText(class = "small text-muted",
              "Accepts: MSstats, SAINT, ProDA, and DIA-limma outputs."),
            uiOutput("ce_status_msg")
          )
        ),

        # -- Step 2: Experiment labels (only when >= 2 files loaded) --
        conditionalPanel(
          condition = "output.ce_ready",
          tags$details(class = "step-section", open = TRUE,
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", "2"), "Experiment Labels"),
            tags$div(class = "step-body",
              helpText(class = "small text-muted",
                "Edit labels used in the Venn diagrams and tables:"),
              uiOutput("ce_exp_name_inputs")
            )
          )
        ),

        # -- Step 3: Highlight gene (press to render) --
        conditionalPanel(
          condition = "output.ce_ready",
          tags$details(class = "step-section",
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", "3"), "Highlight Gene"),
            tags$div(class = "step-body",
              helpText(class = "small text-muted",
                "Type a gene name and press Highlight to mark it in the Venn diagrams."),
              textInput("ce_search_input", NULL,
                        placeholder = "e.g. EGFR, TP53",
                        width = "100%"),
              div(class = "d-flex gap-2 mb-1",
                actionButton("ce_highlight_btn", "Highlight",
                             class = "btn-warning btn-sm flex-grow-1",
                             icon = icon("crosshairs")),
                actionButton("ce_clear_btn", "Clear",
                             class = "btn-outline-secondary btn-sm",
                             icon = icon("xmark"))
              ),
              uiOutput("ce_search_result_ui")
            )
          )
        )
      ),

      # PTM CROSS-EXPERIMENT COMPARISON
      conditionalPanel(
        condition = "input.suite_nav == 'PTM Cross-Experiment Comparison'",
        h4(icon("atom"), " PTM Cross-Experiment"),

        helpText(class = "small text-muted",
          "Upload 2 or more ", tags$b("PTM DE Master CSV"), " files from ",
          tags$b("PTM \u2192 Pairwise DE \u2192 Download CSV"), "."),

        # -- Step 1: Upload files --
        tags$details(class = "step-section", open = TRUE,
          tags$summary(class = "step-header",
            tags$span(class = "step-badge", "1"), "Upload PTM DE CSV Files"),
          tags$div(class = "step-body",
            fileInput("cptm_files",
                      label         = NULL,
                      multiple      = TRUE,
                      accept        = ".csv",
                      width         = "100%",
                      buttonLabel   = tags$span(icon("upload"), " Browse"),
                      placeholder   = "No files selected"),
            helpText(class = "small text-muted",
              "Accepts PTM_DE_Master_*.csv outputs from the PTM DE routes."),
            uiOutput("cptm_status_msg")
          )
        ),

        # -- Step 2: Experiment labels --
        conditionalPanel(
          condition = "output.cptm_ready",
          tags$details(class = "step-section", open = TRUE,
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", "2"), "Experiment Labels"),
            tags$div(class = "step-body",
              helpText(class = "small text-muted",
                "Edit labels for the Venn diagrams and tables:"),
              uiOutput("cptm_exp_name_inputs")
            )
          )
        ),

        # -- Step 3: Highlight gene / PTM site (press to render) --
        conditionalPanel(
          condition = "output.cptm_ready",
          tags$details(class = "step-section",
            tags$summary(class = "step-header",
              tags$span(class = "step-badge", "3"), "Highlight Gene / PTM"),
            tags$div(class = "step-body",
              helpText(class = "small text-muted",
                "Type a gene or PTM site and press Highlight to mark it in the Venn diagrams."),
              textInput("cptm_search_input", NULL,
                        placeholder = "e.g. EGFR, EGFR_Y1068_Phospho",
                        width = "100%"),
              div(class = "d-flex gap-2 mb-1",
                actionButton("cptm_highlight_btn", "Highlight",
                             class = "btn-warning btn-sm flex-grow-1",
                             icon = icon("crosshairs")),
                actionButton("cptm_clear_btn", "Clear",
                             class = "btn-outline-secondary btn-sm",
                             icon = icon("xmark"))
              ),
              uiOutput("cptm_search_result_ui")
            )
          )
        )
      ),
    ),
    
    # =========================================================================
    # CONTENT PANELS (hidden tabset, switched by suite_nav)
    # =========================================================================
    tabsetPanel(
      id = "main_content_area", type = "hidden",
      
      # --- HOME ---
      tabPanel("Home",
        card(
          card_header("Welcome to OSMSS"),
          card_body(
            h3("One Stop MS Shop (OSMSS)"),
            p(paste0("Version ", OSMSS_VERSION, " (July 8th, 2026)"), br(), br(), "A bioinformatics suite for analysis of protein-level MS data and phosphoproteomics.",
              br(), "Select a module from the navigation to begin."),
            # News & updates: collapsed by default, newest first
            tags$details(class = "alert alert-info py-2 px-3 mb-3",
              tags$summary(style = "cursor: pointer;",
                icon("bullhorn"), tags$strong(" NEWS & UPDATES")),
              tags$ul(class = "mb-0 mt-2 ps-3",
                tags$li(tags$strong("10/16/26"),
                  " TMT channel assignment can now exclude channels from the analysis ",
                  "(e.g. leaving an IgG control out of a multi-condition comparison, ",
                  "which improves the statistics)!"),
                tags$li(tags$strong("10/14/26"),
                  " DIA upload now identifies your file automatically and can work with ",
                  "Spectronaut, DIA-NN, FragPipe and Perseus data!"),
                tags$li(tags$strong("10/13/26"),
                  " DDA upload now identifies your file automatically and can work with ",
                  "Proteome Discoverer, MaxQuant and FragPipe data!"))),
            p(class = "text-muted", "System Ready.")
          )
        )
      ),
      
      # --- LOCAL MS PANELS ---
      tabPanel("Load Database",
        fluidRow(
          column(6,
            card(
              card_header(icon("file-import"), " Import Existing Database"),
              card_body(
                tags$p(class = "text-muted small",
                  "Upload your FragpipeDB.csv and MS_Experiment_DB.csv files to load them into this session."),
                fileInput("fpdb_import_fragpipe",
                          "PTM Dataset (FragpipeDB.csv)",
                          accept = ".csv"),
                fileInput("fpdb_import_experiments",
                          "Experiment Dataset (MS_Experiment_DB.csv)",
                          accept = ".csv"),
                actionButton("fpdb_import_btn", "Import Database",
                             class = "btn-success w-100 mt-2",
                             icon = icon("file-import"))
              )
            )
          ),
          column(6,
            card(
              card_header(icon("plus-circle"), " Create New Database"),
              card_body(
                tags$p(class = "text-muted small",
                  "Start with an empty database. Add experiments and data through the Manage Database page.",
                  " Use Export Database when done to save your work."),
                actionButton("fpdb_new_db_btn", "Create Empty Database",
                             class = "btn-primary w-100 mt-2",
                             icon = icon("plus-circle"))
              )
            )
          )
        ),
        br(),
        uiOutput("fpdb_load_status")
      ),
      tabPanel("Manage Database",
        tags$style(HTML("
          .fpdb-list-item:hover { background-color: #f0faf7; }
          .fpdb-list-active { background-color: #e8f5f0 !important; border-left-color: #18bc9c !important; }
          .fpdb-list-item { transition: background-color 0.15s; }
          .nav-tabs .nav-link { font-size: 0.9em; padding: 6px 14px; }
          #fpdb_card_grid.recalculating,
          #fpdb_detail_panel.recalculating,
          #fpdb_exp_table.recalculating { min-height: 200px; position: relative; opacity: 0.35; }
          #fpdb_card_grid.recalculating::before,
          #fpdb_exp_table.recalculating::before {
            content: 'Loading database...'; position: absolute; top: 50%; left: 50%;
            transform: translate(-50%, -50%); background: rgba(255,255,255,0.95);
            padding: 12px 24px; border-radius: 8px; font-size: 0.95em; color: #666;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1); z-index: 100; white-space: nowrap;
          }
        ")),
        uiOutput("fpdb_no_db_banner"),
        navset_tab(
          nav_panel(
            title = tags$span(icon("id-card"), " Card View"),
            card(
              card_body(
                class = "p-3",
                style = "position: relative; min-height: 400px;",
                # Loading overlay — visible by default, hidden by server once data loads
                tags$div(
                  id = "fpdb_init_loading",
                  style = "position: absolute; top: 0; left: 0; right: 0; bottom: 0;
                           background: #fff; z-index: 999; display: flex;
                           align-items: center; justify-content: center; flex-direction: column;",
                  tags$div(class = "text-center",
                    tags$div(class = "spinner-border text-primary mb-3",
                             role = "status",
                             style = "width: 3rem; height: 3rem;"),
                    tags$h5(class = "text-muted", "Initializing Database..."),
                    tags$p(class = "small text-muted",
                           "Loading experiment data and statistics. This may take a moment.")
                  )
                ),
                layout_columns(
                  col_widths = c(5, 7),
                  # Left: Card grid
                  tags$div(
                    style = "max-height: 80vh; overflow-y: auto;",
                    uiOutput("fpdb_card_grid")
                  ),
                  # Right: Detail panel
                  tags$div(
                    class = "border-start ps-3",
                    style = "max-height: 80vh; overflow-y: auto;",
                    uiOutput("fpdb_detail_panel")
                  )
                )
              )
            )
          ),
          nav_panel(
            title = tags$span(icon("table"), " Table View"),
            card(
              card_header("Experiment Database"),
              card_body(
                DTOutput("fpdb_exp_table"),
                br(),
                verbatimTextOutput("fpdb_db_summary")
              )
            )
          )
        )
      ),
      
      tabPanel("FragPipe PTM Analysis",
        card(
          card_header("Modification Viewer"),
          card_body(
            uiOutput("ms_plot_container"),
            br(),
            verbatimTextOutput("ms_plot_status")
          )
        )
      ),
      
      tabPanel("FragPipe Data Table",
        card(
          card_header("FragPipe Data Explorer"),
          card_body(
            DTOutput("ms_gene_table"),
            br(),
            verbatimTextOutput("ms_table_status")
          )
        )
      ),
      
      # --- UPLOAD DATA PANEL ---
      tabPanel("Upload Data",
        card(
          card_header(icon("upload"), " Upload Data"),
          card_body(
            uiOutput("upload_status_display"),
            tags$hr(),
            uiOutput("upload_data_preview")
          )
        )
      ),

      # --- PROTEOME DISCOVERER PANELS ---
      tabPanel("Differential Expression",
        navset_card_underline(
          title = "Differential Expression Analysis",
          nav_panel("DE Plot",
            plotlyOutput("de_plot", height = "650px"),
            br(),
            verbatimTextOutput("de_status")
          ),
          nav_panel("Data Table",
            DT::DTOutput("de_data_table")
          ),
          nav_panel("Methods",
            methodsTextUI("de_methods")
          ),
          # NMAR tab -- DIA only; content managed server-side via de_nmar_content_ui
          nav_panel(
            title = tagList(icon("microscope"), " NMAR Analysis"),
            value = "de_nmar_tab",
            uiOutput("de_nmar_content_ui")
          ),
          # Not Tested & Low Coverage -- TMT only; content managed server-side via
          # de_tmt_onearm_content_ui. Descriptive table, no inference.
          nav_panel(
            title = tagList(icon("circle-half-stroke"), " Not Tested & Low Coverage"),
            value = "de_tmt_onearm_tab",
            uiOutput("de_tmt_onearm_content_ui")
          )
        )
      ),

      # --- SINGLE-SHOT PANELS ---
      tabPanel("Rank Fold Change",
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            "Single-Shot — Rank Fold Change"
          ),
          card_body(
            navset_tab(
              nav_panel(
                title = tags$span(icon("chart-area"), " MA Plot"),
                plotlyOutput("rfc_ma_plot", height = "650px")
              ),
              nav_panel(
                title = tags$span(icon("arrows-up-down"), " Diagonal Scatter"),
                plotlyOutput("rfc_diag_plot", height = "650px")
              ),
              nav_panel(
                title = tags$span(icon("file-lines"), " Methods"),
                methodsTextUI("rfc_methods")
              )
            ),
            br(),
            verbatimTextOutput("rfc_status"),
            conditionalPanel(
              condition = "input.rfc_impute_toggle == true",
              tags$div(class = "alert alert-warning py-2 px-3 small mt-2",
                icon("triangle-exclamation"),
                " Imputation is active. Triangle markers (\u25b2) indicate proteins ",
                "with imputed abundance values. These fold changes are estimates ",
                "and should be interpreted with caution."
              )
            )
          )
        )
      ),

      tabPanel("Trends Analysis",
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            tags$span("Single-Shot — Trends Analysis"),
            downloadButton("trd_download_csv", "Download CSV",
                           class = "btn-outline-secondary btn-sm",
                           icon = icon("download"))
          ),
          card_body(
            navset_tab(
              nav_panel(
                title = tags$span(icon("chart-line"), " Profile Plot"),
                plotlyOutput("trd_profile_plotly", height = "700px")
              ),
              nav_panel(
                title = tags$span(icon("table"), " Cluster Table"),
                DTOutput("trd_cluster_table")
              ),
              nav_panel(
                title = tags$span(icon("file-lines"), " Methods"),
                methodsTextUI("trd_methods")
              )
            ),
            br(),
            verbatimTextOutput("trd_status")
          )
        )
      ),

      tabPanel("Trend Clusters",
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            tags$span("Replicated Quantitation — Trend Clusters"),
            downloadButton("tcl_download_csv", "Download CSV",
                           class = "btn-outline-secondary btn-sm",
                           icon = icon("download"))
          ),
          card_body(
            navset_tab(
              nav_panel(
                title = tags$span(icon("chart-line"), " Profile Plot"),
                plotlyOutput("tcl_profile_plotly", height = "700px")
              ),
              nav_panel(
                title = tags$span(icon("table"), " Cluster Table"),
                DTOutput("tcl_cluster_table")
              ),
              nav_panel(
                title = tags$span(icon("file-lines"), " Methods"),
                methodsTextUI("tcl_methods")
              )
            ),
            br(),
            verbatimTextOutput("tcl_status")
          )
        )
      ),

      tabPanel("PCA Analysis",
        navset_card_underline(
          title = "PCA Output",
          nav_panel("PCA Plot",
            plotlyOutput("pca_scatter", height = "600px"),
            br(),
            verbatimTextOutput("pca_status")
          ),
          nav_panel("Uncorrected PCA",
            plotlyOutput("pca_scatter_de", height = "600px")
          ),
          nav_panel("Scree Plot",
            plotlyOutput("pca_scree", height = "500px")
          ),
          nav_panel("CV Analysis",
            plotlyOutput("pca_cv_plot", height = "500px")
          ),
          nav_panel("Methods",
            methodsTextUI("pca_methods")
          )
        )
      ),

      tabPanel("Heatmap",
        navset_card_underline(
          title = "Heatmap Analysis",
          nav_panel("Heatmap",
            uiOutput("hm_plot_container"),
            uiOutput("hm_ebn_info"),
            br(),
            verbatimTextOutput("hm_status")
          ),
          nav_panel("Statistics",
            DT::DTOutput("hm_stats_table")
          ),
          nav_panel("Methods",
            methodsTextUI("hm_methods")
          )
        )
      ),

      tabPanel("Protein Coverage",
        card(
          card_header("Proteome Discoverer - Protein Coverage Map"),
          card_body(
            plotlyOutput("cov_plot", height = "650px"),
            br(),
            verbatimTextOutput("cov_status")
          )
        )
      ),

      # --- CROSS-EXPERIMENT COMPARISON ---
      tabPanel("Cross-Experiment Comparison",

        # Info / instructions banner (hidden once files are loaded)
        uiOutput("ce_panel_header"),

        # Venn diagrams — shown only when >= 2 files are ready
        conditionalPanel(
          condition = "output.ce_ready",

          # Two Venn diagrams side by side
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header(
                class = "d-flex align-items-center gap-2",
                tags$span(style = "color:#dc3545;", icon("arrow-up")),
                " Upregulated Genes"
              ),
              card_body(
                helpText(class = "small text-muted",
                         "Log\u2082FC > 1 \u2227 FDR / BFDR < 0.05"),
                plotOutput("ce_venn_up", height = "400px")
              )
            ),
            card(
              card_header(
                class = "d-flex align-items-center gap-2",
                tags$span(style = "color:#0d6efd;", icon("arrow-down")),
                " Downregulated Genes"
              ),
              card_body(
                helpText(class = "small text-muted",
                         "Log\u2082FC < \u22121 \u2227 FDR / BFDR < 0.05"),
                plotOutput("ce_venn_down", height = "400px")
              )
            )
          ),

          # Per-experiment summary stats
          card(
            card_header(
              class = "d-flex align-items-center gap-2",
              icon("table"), " Per-Experiment Summary"
            ),
            card_body(
              helpText(class = "small text-muted mb-2",
                tags$b("\u2020 Unique:"),
                " gene regulated in this experiment but absent (in the same direction) from every other uploaded experiment."
              ),
              DT::DTOutput("ce_summary_table")
            )
          ),

          # Gene-level overlap table
          card(
            card_header(
              class = "d-flex justify-content-between align-items-center",
              tags$span(class = "d-flex align-items-center gap-2",
                        icon("table-list"), " Regulated Gene Overlap"),
              downloadButton("ce_dl_overlap", "Download CSV",
                             class = "btn-outline-secondary btn-sm py-0")
            ),
            card_body(
              helpText(class = "small text-muted mb-2",
                "All regulated genes across experiments. ",
                tags$span(style = "background:#ffe0e0; padding: 1px 5px; border-radius:3px;
                           color:#9b0000; font-weight:bold;", "Up"),
                " / ",
                tags$span(style = "background:#dde8ff; padding: 1px 5px; border-radius:3px;
                           color:#003580; font-weight:bold;", "Down"),
                " per experiment. Sort or filter to find shared / unique genes."
              ),
              DT::DTOutput("ce_overlap_table")
            )
          )
        )
      ),

      tabPanel("TMT Upload Data",
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            tags$span("TMT Upload Data"),
            tags$span(class = "badge bg-info", "Data Management")
          ),
          card_body(
            class = "p-3",
            tags$div(class = "text-center py-4",
              icon("upload", class = "fa-3x text-muted mb-3"), tags$br(),
              tags$h5("Upload and configure your TMT data here."),
              tags$p(class = "text-muted",
                "Use the sidebar to upload files, set up your experimental design, ",
                "choose your analysis mode, and summarize your data. ",
                "Once summarized, navigate to any TMT analysis module to begin."),
              uiOutput("tmt_upload_main_status")
            )
          )
        )
      ),

      tabPanel("PTM Upload Data",
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            tags$span("PTM Upload Data"),
            tags$span(class = "badge bg-info", "Data Management")
          ),
          card_body(
            class = "p-3",
            tags$div(class = "py-3",
              uiOutput("tmt_ptm_upload_main_status")
            )
          )
        )
      ),

      tabPanel("PTM PCA",
        card(
          card_header("PTM Principal Component Analysis"),
          card_body(
            navset_tab(
              id = "ptm_pca_tabs",
              nav_panel(
                title = tags$span(icon("braille"), " PTM PCA (Unadjusted)"),
                value = "ptm_pca",
                plotlyOutput("ptm_pca_ptm_scatter", height = "500px"),
                br(),
                plotlyOutput("ptm_pca_ptm_scree", height = "280px")
              ),
              nav_panel(
                title = tags$span(icon("braille"), " Global Protein PCA"),
                value = "prot_pca",
                plotlyOutput("ptm_pca_prot_scatter", height = "500px"),
                br(),
                plotlyOutput("ptm_pca_prot_scree", height = "280px")
              ),
              nav_panel(
                title = tags$span(icon("chart-bar"), " CV Analysis"),
                value = "cv",
                plotlyOutput("ptm_pca_cv_plot", height = "450px")
              )
            )
          )
        )
      ),

      tabPanel("PTM DE",
        card(
          card_header("PTM Pairwise Differential Expression"),
          card_body(
            navset_tab(
              id = "ptm_de_tabs",
              nav_panel(
                title = tags$span(icon("volcano"), " Volcano Plot"),
                value = "volcano",
                plotlyOutput("ptm_volcano_plotly", height = "650px")
              ),
              nav_panel(
                title = tags$span(icon("braille"), " FC vs Global"),
                value = "fc_scatter",
                plotlyOutput("ptm_fc_scatter_plotly", height = "650px")
              ),
              nav_panel(
                title = tags$span(icon("table"), " DE Table"),
                value = "de_table",
                DT::DTOutput("ptm_de_table")
              ),
              nav_panel(
                title = tags$span(icon("eye-slash"), " Detection Test"),
                value = "nmar",
                tags$div(class = "alert alert-info py-2 px-3 small mt-2",
                  tags$b("Detection concerns, not fold-change results."),
                  " This table lists every site in the comparison whose detection ",
                  "pattern needs looking at, in two groups shown in the Assessment ",
                  "column. Sites ", tags$b("excluded from the model"),
                  " had too few valid values for limma to fit, so no fold change ",
                  "exists for them; they are assessed instead by Fisher's exact test ",
                  "on detected / not-detected counts, corrected within that excluded ",
                  "set only. Sites ", tags$b("modelled by limma"),
                  " passed the valid-value filter but have a condition with zero ",
                  "detections or with more than 75% of samples missing, so their fold ",
                  "change rests on very little data. The two carry separate statistic ",
                  "columns and no joint false discovery rate is claimed across them — ",
                  "nor across this table and the DE table. Where the replicate count ",
                  "makes the smallest attainable Fisher p-value exceed 0.05, those ",
                  "p-values are reported as NA and that part of the table is a ranked ",
                  "descriptive list."),
                DT::DTOutput("ptm_nmar_table")
              ),
              nav_panel(
                title = tags$span(icon("file-lines"), " Methods"),
                value = "methods",
                methodsTextUI("ptm_de_methods")
              )
            ),
            br(),
            verbatimTextOutput("ptm_status")
          )
        )
      ),

      tabPanel("PTM Trends",
        card(
          card_header(
            class = "d-flex justify-content-between align-items-center",
            tags$span("PTM Trend Clusters (Fuzzy c-means Clustering)")
          ),
          card_body(
            navset_tab(
              id = "ptm_trends_tabs",
              nav_panel(
                # "Adjusted PTM Trends", or "Absolute PTM Trends" for a PTM-only cache
                title = tags$span(icon("chart-line"), " ",
                                  textOutput("ptm_trends_main_tab_title", inline = TRUE)),
                value = "main",
                plotlyOutput("ptm_trends_profile_plotly", height = "700px"),
                br(),
                DT::DTOutput("ptm_trends_cluster_table")
              ),
              nav_panel(
                title = tags$span(icon("chart-line"), " Total PTM Trends"),
                value = "total",
                plotlyOutput("ptm_total_profile_plotly", height = "700px"),
                br(),
                DT::DTOutput("ptm_total_cluster_table")
              ),
              nav_panel(
                title = tags$span(icon("filter-circle-xmark"), " Excluded Sites"),
                tags$p(class = "small text-muted mt-2",
                  "Sites that did not meet the site-eligibility threshold in every ",
                  "condition. They are excluded from the cluster figures because ",
                  "clustering needs a complete profile and nothing is imputed — they ",
                  "are not missing from the dataset, and their measured values remain ",
                  "viewable in the Single PTM tab."),
                DT::DTOutput("ptm_trends_excluded_table")
              ),
              nav_panel(
                title = tags$span(icon("file-lines"), " Methods"),
                methodsTextUI("ptm_trnd_methods")
              )
            ),
            br(),
            verbatimTextOutput("ptm_trends_status")
          )
        )
      ),

      tabPanel("PTM Single",
        card(
          card_header("Single PTM Inspector"),
          card_body(
            uiOutput("psg_na_warning"),
            navset_tab(
              nav_panel(
                title = tags$span(icon("chart-line"), " Mechanistic Fold Changes"),
                plotlyOutput("psg_fc_plot", height = "500px")
              ),
              nav_panel(
                title = tags$span(icon("chart-column"), " Summarized Abundance"),
                plotlyOutput("psg_abundance_plot", height = "500px")
              ),
              nav_panel(
                title = tags$span(icon("file-lines"), " Methods"),
                methodsTextUI("ptm_sg_methods")
              )
            ),
            br(),
            verbatimTextOutput("psg_plot_status")
          )
        )
      ),

      # --- KS ANALYSIS PANEL ---
      tabPanel("KS Analysis",
        card(
          full_screen = TRUE,
          style = "height: calc(100vh - 80px); min-height: 500px;",
          card_header(
            class = "d-flex justify-content-between align-items-center",
            tags$span(icon("diagram-project"), " Kinase-Site Analysis"),
            tags$span(class = "badge bg-secondary", "PhosR algorithm")
          ),
          card_body(
            class = "p-2 overflow-auto",
            style = "height: 100%;",
            navset_tab(
              id = "ks_result_tabs",

              # Tab 1: Upregulated / pre-filter heatmap
              nav_panel(
                title = tags$span(icon("arrow-up"), " KS Heatmap (Up / Selected)"),
                conditionalPanel(
                  condition = "output.ks_has_result_up",
                  div(style = "overflow-x: auto;",
                    plotOutput("ks_heatmap_up", width = "100%"))
                ),
                conditionalPanel(
                  condition = "!output.ks_has_result_up",
                  div(class = "text-center py-5 text-muted",
                    icon("diagram-project", class = "fa-3x mb-3"), tags$br(),
                    tags$h5("No results yet"),
                    tags$p("Upload a CSV and click ",
                      tags$b("Run PhosR Predictions"), " in the sidebar.")
                  )
                )
              ),

              # Tab 2: Downregulated heatmap (auto DE mode only)
              nav_panel(
                title = tags$span(icon("arrow-down"), " KS Heatmap (Down)"),
                conditionalPanel(
                  condition = "output.ks_has_result_down",
                  div(style = "overflow-x: auto;",
                    plotOutput("ks_heatmap_down", width = "100%"))
                ),
                conditionalPanel(
                  condition = "!output.ks_has_result_down",
                  div(class = "text-center py-5 text-muted",
                    icon("circle-info", class = "fa-2x mb-3"), tags$br(),
                    tags$p(class = "fst-italic",
                      "Downregulated heatmap is produced when the CSV contains ",
                      tags$code("adj.pvalue"), " and ", tags$code("log2FC"),
                      " columns and pre-filter is OFF.")
                  )
                )
              ),

              # Tab 3: Kinase activity heatmap
              nav_panel(
                title = tags$span(icon("gauge-high"), " Kinase Activity"),
                conditionalPanel(
                  condition = "output.ks_has_result_up || output.ks_has_result_down",
                  div(class = "alert alert-warning mx-2 mt-2 mb-1 py-2 px-3",
                    style = "font-size: 0.82rem; border-left: 4px solid #f0ad4e;",
                    tags$div(class = "d-flex gap-2 align-items-start",
                      icon("triangle-exclamation", class = "mt-1 flex-shrink-0"),
                      tags$div(
                        tags$strong("Kinase Activity — scope depends on your experimental design:"),
                        tags$ol(class = "mb-0 mt-1 ps-3",
                          tags$li(
                            tags$strong("PTMs used: "),
                            "By default, only PTMs significant in the comparison used to",
                            " generate this file are included. For a ", tags$strong("2-condition"),
                            " experiment (e.g. lesion vs. control) this is ideal — it focuses",
                            " on sites that genuinely changed, preventing non-changing sites",
                            " from diluting kinase activity scores. For a ",
                            tags$strong("multi-condition"), " time course (e.g. 0, 12, 24, 48 h),",
                            " this can be limiting: sites that only change at 48 h won't appear",
                            " at all if the input was derived from a 0 h vs. 12 h comparison,",
                            " causing kinase activity at later time points to be underestimated."
                          ),
                          tags$li(class = "mt-1",
                            tags$strong("Kinases shown: "),
                            "Unless \u201cShow all kinases\u201d is enabled, only kinases",
                            " whose substrates were significant in that specific comparison",
                            " are displayed. Enable it to see the full kinase set."
                          )
                        ),
                        conditionalPanel(
                          condition = "output.ks_has_plex",
                          tags$div(class = "mt-2 d-flex align-items-center gap-2",
                            tags$span(class = "badge bg-info text-dark flex-shrink-0",
                              icon("shield-halved"), " TMT"),
                            tags$span(class = "small",
                              "Plex/batch effects are automatically corrected using ",
                              tags$strong("RUVphospho"), " (k\u00a0=\u00a01) before scoring.",
                              " All three result tabs benefit from this correction.")
                          )
                        )
                      )
                    )
                  ),
                  conditionalPanel(
                    condition = "output.ks_has_result_up || output.ks_has_result_down",
                    div(class = "d-flex align-items-center gap-3 mx-2 mb-2",
                      actionButton("ks_run_full_activity",
                                   "Recompute Activity \u2014 All PTMs (multi-condition)",
                                   class = "btn-success btn-sm",
                                   icon  = icon("rotate")),
                      conditionalPanel(
                        condition = "output.ks_activity_mode == 'full'",
                        tags$span(class = "badge bg-success",
                          icon("check"), " Showing full-dataset activity")
                      )
                    )
                  ),
                  div(style = "overflow-x: auto;",
                    plotOutput("ks_activity_heatmap", width = "100%",
                               height = "540px"))
                ),
                conditionalPanel(
                  condition = "!output.ks_has_result_up && !output.ks_has_result_down",
                  div(class = "text-center py-5 text-muted",
                    icon("gauge-high", class = "fa-3x mb-3"), tags$br(),
                    tags$h5("No results yet"),
                    tags$p("Run the analysis to see kinase activity profiles.")
                  )
                )
              ),

              # Tab 4: Kinase ranking (lollipop chart — all kinases, compact)
              nav_panel(
                title = tags$span(icon("ranking-star"), " Kinase Ranking"),
                conditionalPanel(
                  condition = "output.ks_has_result_up || output.ks_has_result_down",
                  tags$p(class = "text-muted small px-2 pt-2 mb-1",
                    icon("circle-info"),
                    " All scored kinases ranked by mean KS score, coloured by kinase group.",
                    " Scroll to explore. Download PDF for a publication-ready version."),
                  div(style = "overflow-y: auto;",
                    plotOutput("ks_ranking_plot", width = "100%"))
                ),
                conditionalPanel(
                  condition = "!output.ks_has_result_up && !output.ks_has_result_down",
                  div(class = "text-center py-5 text-muted",
                    icon("ranking-star", class = "fa-3x mb-3"), tags$br(),
                    tags$h5("No results yet"),
                    tags$p("Run the analysis to see the kinase ranking.")
                  )
                )
              )
            )
          )
        )
      ),

      # --- PTM CROSS-EXPERIMENT COMPARISON ---
      tabPanel("PTM Cross-Experiment Comparison",

        # Info banner (shown before files are loaded)
        uiOutput("cptm_panel_header"),

        conditionalPanel(
          condition = "output.cptm_ready",

          navset_card_underline(
            title = tags$span(icon("atom"), " PTM Cross-Experiment Comparison"),

            # ── Tab 1: Adjusted PTM (site-level) ─────────────────────────────
            nav_panel(
              title = tags$span(icon("star-of-life"), " Adjusted PTM"),
              value = "cptm_tab_adj",
              layout_columns(
                col_widths = c(6, 6),
                card(
                  card_header(tags$span(style = "color:#dc3545;", icon("arrow-up")),
                               " Up-regulated PTM Sites"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Log\u2082FC > 1 \u2227 Adjusted adj.p < 0.05"),
                    plotOutput("cptm_venn_adj_up", height = "380px")
                  )
                ),
                card(
                  card_header(tags$span(style = "color:#0d6efd;", icon("arrow-down")),
                               " Down-regulated PTM Sites"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Log\u2082FC < \u22121 \u2227 Adjusted adj.p < 0.05"),
                    plotOutput("cptm_venn_adj_down", height = "380px")
                  )
                )
              ),
              card(
                card_header(icon("table"), " Per-Experiment Summary (Adjusted PTM)"),
                card_body(
                  helpText(class = "small text-muted mb-2",
                    tags$b("\u2020 Unique:"), " site present in this experiment but absent ",
                    "(same direction) from every other experiment."),
                  DT::DTOutput("cptm_summary_adj")
                )
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  tags$span(icon("table-list"), " PTM Site Overlap (Adjusted)"),
                  downloadButton("cptm_dl_overlap_adj", "CSV",
                                 class = "btn-outline-secondary btn-sm py-0")
                ),
                card_body(
                  helpText(class = "small text-muted mb-2",
                    "All regulated PTM sites. Filter or sort to find shared / unique sites. ",
                    tags$span(style = "background:#fff3cd; padding:1px 5px; border-radius:3px;",
                              "Discordant"),
                    " = site up in one experiment, down in another."
                  ),
                  DT::DTOutput("cptm_overlap_adj")
                )
              )
            ),

            # ── Tab 2: Absolute PTM (site-level) ─────────────────────────────
            nav_panel(
              title = tags$span(icon("wave-square"), " Absolute PTM"),
              value = "cptm_tab_abs",
              layout_columns(
                col_widths = c(6, 6),
                card(
                  card_header(tags$span(style = "color:#dc3545;", icon("arrow-up")),
                               " Up-regulated PTM Sites"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Log\u2082FC > 1 \u2227 Absolute adj.p < 0.05"),
                    plotOutput("cptm_venn_abs_up", height = "380px")
                  )
                ),
                card(
                  card_header(tags$span(style = "color:#0d6efd;", icon("arrow-down")),
                               " Down-regulated PTM Sites"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Log\u2082FC < \u22121 \u2227 Absolute adj.p < 0.05"),
                    plotOutput("cptm_venn_abs_down", height = "380px")
                  )
                )
              ),
              card(
                card_header(icon("table"), " Per-Experiment Summary (Absolute PTM)"),
                card_body(DT::DTOutput("cptm_summary_abs"))
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  tags$span(icon("table-list"), " PTM Site Overlap (Absolute)"),
                  downloadButton("cptm_dl_overlap_abs", "CSV",
                                 class = "btn-outline-secondary btn-sm py-0")
                ),
                card_body(DT::DTOutput("cptm_overlap_abs"))
              )
            ),

            # ── Tab 3: Protein Only ───────────────────────────────────────────
            nav_panel(
              title = tags$span(icon("circle-nodes"), " Protein Only"),
              value = "cptm_tab_prot",
              layout_columns(
                col_widths = c(6, 6),
                card(
                  card_header(tags$span(style = "color:#dc3545;", icon("arrow-up")),
                               " Up-regulated Proteins"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Protein-level model: Log\u2082FC > 1 \u2227 adj.p < 0.05"),
                    plotOutput("cptm_venn_prot_up", height = "380px")
                  )
                ),
                card(
                  card_header(tags$span(style = "color:#0d6efd;", icon("arrow-down")),
                               " Down-regulated Proteins"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Protein-level model: Log\u2082FC < \u22121 \u2227 adj.p < 0.05"),
                    plotOutput("cptm_venn_prot_down", height = "380px")
                  )
                )
              ),
              card(
                card_header(icon("table"), " Per-Experiment Summary (Protein Only)"),
                card_body(DT::DTOutput("cptm_summary_prot"))
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  tags$span(icon("table-list"), " Protein Overlap"),
                  downloadButton("cptm_dl_overlap_prot", "CSV",
                                 class = "btn-outline-secondary btn-sm py-0")
                ),
                card_body(DT::DTOutput("cptm_overlap_prot"))
              )
            ),

            # ── Tab 4: Adjusted — Collapsed by Protein ────────────────────────
            nav_panel(
              title = tags$span(icon("layer-group"),
                                " Adjusted \u2014 Collapsed by Protein"),
              value = "cptm_tab_collprot",
              tags$div(
                class = "alert alert-secondary py-2 px-3 small mt-2 mb-3",
                icon("circle-info"),
                " Uses the ", tags$b("Adjusted PTM"), " model but collapses all sites ",
                "per gene. A gene is counted as regulated if ",
                tags$b("any"), " of its sites meets the threshold."
              ),
              layout_columns(
                col_widths = c(6, 6),
                card(
                  card_header(tags$span(style = "color:#dc3545;", icon("arrow-up")),
                               " Up-regulated Proteins"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Any site: Adjusted Log\u2082FC > 1 \u2227 adj.p < 0.05"),
                    plotOutput("cptm_venn_collprot_up", height = "380px")
                  )
                ),
                card(
                  card_header(tags$span(style = "color:#0d6efd;", icon("arrow-down")),
                               " Down-regulated Proteins"),
                  card_body(
                    helpText(class = "small text-muted",
                             "Any site: Adjusted Log\u2082FC < \u22121 \u2227 adj.p < 0.05"),
                    plotOutput("cptm_venn_collprot_down", height = "380px")
                  )
                )
              ),
              card(
                card_header(icon("table"),
                             " Per-Experiment Summary (Adjusted, Collapsed by Protein)"),
                card_body(DT::DTOutput("cptm_summary_collprot"))
              ),
              card(
                card_header(
                  class = "d-flex justify-content-between align-items-center",
                  tags$span(icon("table-list"), " Protein Overlap (Collapsed)"),
                  downloadButton("cptm_dl_overlap_collprot", "CSV",
                                 class = "btn-outline-secondary btn-sm py-0")
                ),
                card_body(
                  helpText(class = "small text-muted mb-2",
                    tags$span(style = "background:#fff3cd; padding:1px 5px; border-radius:3px;",
                              "Discordant"),
                    " = protein has at least one up-regulated AND one down-regulated site."
                  ),
                  DT::DTOutput("cptm_overlap_collprot")
                )
              )
            )

          ) # end navset_card_underline
        )   # end conditionalPanel
      ),

      # --- GSEA ANALYSIS PANEL ---
      tabPanel("GSEA Analysis",
        # Status banner
        card(
          card_body(
            class = "py-2",
            verbatimTextOutput("gsea_status")
          )
        ),

        # Results table
        card(
          card_header(icon("table"), " GSEA Results"),
          card_body(
            class = "p-2",
            helpText(class = "small text-muted mb-2",
              "Click a row to view its enrichment plot below. ",
              "Colour legend: ",
              tags$span(style = "color:#e74c3c; font-weight:bold;", "Red = positive NES"),
              " / ",
              tags$span(style = "color:#3498db; font-weight:bold;", "Blue = negative NES"),
              ". Background: ",
              tags$span(style = "background:#d4edda; padding:0 4px;", "FDR < 0.01"), " ",
              tags$span(style = "background:#d1ecf1; padding:0 4px;", "< 0.05"), " ",
              tags$span(style = "background:#fff3cd; padding:0 4px;", "< 0.25")
            ),
            DT::DTOutput("gsea_result_table", width = "100%")
          )
        ),

        # Enrichment plot + Leading Edge + Overview tabs
        navset_card_tab(
          id = "gsea_detail_tabs",

          # Tab 1: Enrichment Plot (3-panel)
          nav_panel(
            title = tagList(icon("chart-area"), " Enrichment Plot"),
            card_body(
              class = "p-2",
              helpText(class = "small text-muted mb-2",
                "Classic 3-panel enrichment plot: running enrichment score, ",
                "gene set hits (barcode), and ranked metric distribution. ",
                "Click a pathway in the table above to update."),
              plotlyOutput("gsea_enrichment_plot", height = "550px")
            )
          ),

          # Tab 2: Leading Edge Table
          nav_panel(
            title = tagList(icon("list"), " Leading Edge"),
            card_body(
              class = "p-2",
              uiOutput("gsea_leading_edge_header"),
              hr(class = "my-2"),
              helpText(class = "small text-muted mb-2",
                "Gene-level details for the selected pathway. ",
                tags$span(style = "background:#d4edda; padding:0 4px; font-weight:bold;",
                          "Green"),
                " = core enrichment genes (leading edge)."),
              DT::DTOutput("gsea_leading_edge_table", width = "100%")
            )
          ),

          # Tab 3: Overview (Dot plot + NES)
          nav_panel(
            title = tagList(icon("chart-bar"), " Overview Plots"),
            card_body(
              class = "p-2",
              layout_columns(
                col_widths = c(6, 6),
                card(
                  card_header(icon("circle"), " Dot Plot"),
                  card_body(
                    class = "p-2",
                    plotlyOutput("gsea_dotplot", height = "500px")
                  )
                ),
                card(
                  card_header(icon("chart-bar"), " NES Bar Chart"),
                  card_body(
                    class = "p-2",
                    plotlyOutput("gsea_nes_barplot", height = "500px")
                  )
                )
              )
            )
          )
        )
      ),

      tabPanel("DIA Upload Data",
        card(
          card_header(
            "DIA Upload Data",
            tags$span(class = "badge bg-info", "Data Management")
          ),
          card_body(class = "p-3",
            tags$div(class = "text-center py-4",
              icon("upload", class = "fa-3x text-muted mb-3"), tags$br(),
              tags$h5("Upload and configure your DIA data here."),
              tags$p(class = "text-muted",
                "Use the sidebar wizard to upload, configure, and summarize."),
              uiOutput("dia_upload_main_status")
            )
          )
        )
      ),

      tabPanel("Label Free PTM Upload",
        card(
          card_header(
            "Label Free Phosphoproteomics Upload",
            tags$span(class = "badge bg-success", "Label Free PTM")
          ),
          card_body(class = "p-3",
            tags$div(class = "text-center py-4",
              icon("flask", class = "fa-3x text-muted mb-3"), tags$br(),
              tags$h5("Upload and configure your label-free phosphoproteomics data here."),
              tags$p(class = "text-muted",
                "Use the sidebar wizard to upload the phospho-enriched and global ",
                "proteomics files, assign conditions, and summarize. ",
                "Results will be available in PTM DE, PTM Trends, and PTM Single."),
              uiOutput("dia_ptm_upload_main_status")
            )
          )
        )
      ),

      tabPanel("Under Construction",
        card(
          card_header(
            class = "bg-warning text-dark",
            icon("triangle-exclamation"), " Module Under Construction"
          ),
          card_body(
            class = "text-center py-5",
            tags$div(
              style = "font-size: 5em; color: #f0ad4e; margin-bottom: 20px;",
              icon("helmet-safety")
            ),
            h2("Coming Soon"),
            p(class = "lead text-muted",
              "This module is currently under active development."),
            p(class = "text-muted",
              "Check back in a future release for full functionality."),
            hr(),
            p(class = "small text-muted",
              "OSMSS ", tags$b("v5"), " — Coulombe Lab, University of Michigan")
          )
        )
      ),

      # --- RESIDUE PREDICTION PANELS ---
      tabPanel("Kinase Prediction",
        navset_card_underline(
          title = "Kinase Output",
          nav_panel("Visualization", plotlyOutput("main_plot_kin", height = "600px")),
          nav_panel("Log", verbatimTextOutput("console_log_kin"))
        )
      ),
      
      tabPanel("14-3-3 Prediction",
        navset_card_underline(
          title = "14-3-3 Output",
          nav_panel("Visualization",
                    p(class="text-muted", ATTR_1433),
                    plotlyOutput("main_plot_1433", height = "600px")),
          nav_panel("Log", verbatimTextOutput("console_log_1433"))
        )
      ),
      
      tabPanel("ScanSite - Upstream/Downstream",
        navset_card_underline(
          title = "Downstream Binding Analysis",
          nav_panel("Downstream Binders", 
                    p(class="text-muted", ATTR_SCANSITE),
                    p(class="text-muted",
                      "Green = Predicted Downstream Binding Domains (SH2, 14-3-3, etc)."),
                    plotlyOutput("bind_plot_down", height = "600px")),
          nav_panel("Upstream Kinases",
                    p(class="text-muted", ATTR_SCANSITE),
                    p(class="text-muted",
                      "Red = Predicted Upstream Kinase Motifs."),
                    plotlyOutput("bind_plot_up", height = "600px")),
          nav_panel("Site Context Report", verbatimTextOutput("bind_report")),
          nav_panel("Full Data", DT::dataTableOutput("table_raw_bind")),
          nav_panel("Log", verbatimTextOutput("console_log_bind"))
        )
      ),
      
      tabPanel("Multi-Site Analysis",
        navset_card_underline(
          title = "Comparative Analysis",
          nav_panel("Graph View", plotlyOutput("multisite_plot", height = "700px")),
          nav_panel("Table View", DT::DTOutput("multisite_table")),
          nav_panel("Log", verbatimTextOutput("console_log_multi"))
        )
      ),
      
      tabPanel("Protein Kinome",
        navset_card_underline(
          title = "Hypothesis Generator",
          nav_panel("Interactive Heatmap",
                    plotlyOutput("heatmap_plotly", height = "auto"),
                    hr(),
                    h4("3D Structural Context"),
                    NGLVieweR::NGLVieweROutput("structure", height = "500px")),
          nav_panel("Log", verbatimTextOutput("console_log_map"))
        )
      ),

# =====================================================================
      # SAVED EXPERIMENTS (Main Content Window)
      # =====================================================================
      tabPanel("Saved Experiments",
        tags$style(HTML("
          .prerun-tile-grid     { min-height:300px; padding-top: 15px; }
          /* Hide Shiny fileInputs driven programmatically by JS */
          #prerun_file_input_wrapper,
          #prerun_auto_rds_wrapper { display:none !important; }
        ")),
        
        tags$script(HTML("
          // ── Folder-level state (populated on folder select, cleared on re-select) ──
          window.prerunFullRdsMap      = {};
          window.prerunPreviewRelPaths = [];

          // ── Open folder picker ─────────────────────────────────────────────────
          $(document).on('click', '#prerun_folder_btn_js', function() {
            document.getElementById('prerun_folder_trigger').click();
          });

          // ── Folder selected ────────────────────────────────────────────────────
          $(document).on('change', '#prerun_folder_trigger', function(e) {
            var allFiles = Array.from(e.target.files);

            var previews = allFiles.filter(function(f) {
              return f.name.toLowerCase().endsWith('.preview.rds');
            });
            var fullRds = allFiles.filter(function(f) {
              var n = f.name.toLowerCase();
              return n.endsWith('.rds') && !n.endsWith('.preview.rds');
            });

            window.prerunFullRdsMap = {};
            fullRds.forEach(function(f) {
              window.prerunFullRdsMap[f.name] = f;
            });

            if (previews.length === 0) {
              window.prerunPreviewRelPaths = [];
              Shiny.setInputValue('prerun_no_previews', Math.random(), {priority: 'event'});
              e.target.value = '';
              return;
            }

            window.prerunPreviewRelPaths = previews.map(function(f) {
              return f.webkitRelativePath || f.name;
            });

            Shiny.setInputValue('prerun_preview_paths', window.prerunPreviewRelPaths, {priority: 'event'});

            var dt = new DataTransfer();
            previews.forEach(function(f) { dt.items.add(f); });
            var inp = document.querySelector('#prerun_file_input_wrapper input[type=\"file\"]');
            inp.files = dt.files;
            $(inp).trigger('change');

            var rel     = previews[0].webkitRelativePath || '';
            var folder  = rel ? rel.split('/')[0] : '(files)';
            var folders = {};
            previews.forEach(function(f) {
              var parts = (f.webkitRelativePath || f.name).split('/');
              if (parts.length > 1) folders[parts.slice(0, -1).join('/')] = 1;
            });

            Shiny.setInputValue('prerun_folder_info', {
              folder:   folder,
              n:        previews.length,
              nFolders: Object.keys(folders).length
            }, {priority: 'event'});
          });

          // ── Load button click (plain button, class prerun-load-btn) ───────────
          $(document).on('click', '.prerun-load-btn', function() {
            var idx          = parseInt($(this).data('idx'), 10);
            var relPath      = window.prerunPreviewRelPaths[idx - 1] || '';
            var previewName  = relPath.split('/').pop() || relPath;
            var expectedName = previewName.replace(/\\.preview\\.rds$/i, '.rds');
            var fullFile     = window.prerunFullRdsMap[expectedName];

            if (fullFile) {
              Shiny.setInputValue('prerun_auto_load_request',
                {idx: idx, name: expectedName, stamp: Math.random()},
                {priority: 'event'});

              var autoInp = document.querySelector('#prerun_auto_rds_wrapper input[type=\"file\"]');
              var dt2 = new DataTransfer();
              dt2.items.add(fullFile);
              autoInp.files = dt2.files;
              $(autoInp).trigger('change');
            } else {
              Shiny.setInputValue('prerun_manual_load_request',
                {idx: idx, stamp: Math.random()},
                {priority: 'event'});
            }
          });
        ")),
        
        tags$div(class = "prerun-tile-grid",
          uiOutput("prerun_grid")
        )
      )

    )
  )
)
