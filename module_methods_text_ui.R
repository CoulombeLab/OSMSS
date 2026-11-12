# module_methods_text_ui.R =====================================================
# Reusable Shiny UI module: renders publication-ready methods text from the
# provenance ledger.  Includes a Markdown rendering box and a Copy-to-Clipboard
# button implemented via the browser Clipboard API (no extra R packages).
# ==============================================================================

methodsTextUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      class = "card mt-3",
      tags$div(
        class = "card-header d-flex justify-content-between align-items-center",
        tags$span(icon("file-lines"), " Auto-Generated Methods Text"),
        tags$button(
          id = ns("copy_btn"),
          class = "btn btn-outline-secondary btn-sm",
          onclick = paste0(
            "var el = document.getElementById('", ns("methods_md"), "');",
            "if (el) { navigator.clipboard.writeText(el.innerText).then(",
            "function() { ",
            "  var btn = document.getElementById('", ns("copy_btn"), "');",
            "  btn.innerHTML = '<i class=\"fa fa-check\"></i> Copied!';",
            "  setTimeout(function() { btn.innerHTML = '<i class=\"fa fa-copy\"></i> Copy to Clipboard'; }, 2000);",
            "},",
            "function() { alert('Clipboard copy failed. Please select and copy manually.'); }",
            "); }"
          ),
          icon("copy"), " Copy to Clipboard"
        )
      ),
      tags$div(
        class = "card-body",
        uiOutput(ns("methods_md"))
      )
    )
  )
}
