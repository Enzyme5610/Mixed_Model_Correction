library(shiny)
library(bslib)
library(lme4)
library(car)

# ---- Helpers ---------------------------------------------------------------

# Split a pasted list on commas, semicolons, tabs, spaces or newlines, then
# convert names the same way read.csv() does (e.g. "HLA-A" -> "HLA.A").
parse_genes <- function(txt) {
  g <- trimws(unlist(strsplit(txt, "[,;\\s]+", perl = TRUE)))
  g <- g[nzchar(g)]
  unique(make.names(g))
}

# Pick the column whose name matches `want` (case-insensitive), else the first.
guess_col <- function(cols, want) {
  hit <- cols[tolower(cols) == tolower(want)]
  if (length(hit)) hit[1] else cols[1]
}

# Same model and test as the original script:
#   f <- reformulate(c("Tx", "(1|Line/Batch)"), var)
#   MM_Form <- lmer(f, data = datos)
#   Anova(MM_Form, type = "II", test = "F")
fit_one <- function(var, datos, tx, line, batch) {
  out <- list(gene = var, anova = NULL, note = "")
  if (!var %in% names(datos)) {
    out$note <- "Column not found in file"
    return(out)
  }
  if (!is.numeric(datos[[var]])) {
    out$note <- "Column is not numeric"
    return(out)
  }
  f <- reformulate(c(tx, sprintf("(1|%s/%s)", line, batch)), var)
  warn <- character()
  MM_Form <- tryCatch(
    withCallingHandlers(
      suppressMessages(lmer(f, data = datos)),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(MM_Form, "error")) {
    out$note <- paste("Model failed:", conditionMessage(MM_Form))
    return(out)
  }
  out$anova <- tryCatch(
    Anova(MM_Form, type = "II", test = "F"),
    error = function(e) e
  )
  notes <- c(if (isSingular(MM_Form)) "Singular fit (a random-effect variance is ~0)", warn)
  if (inherits(out$anova, "error")) {
    notes <- c(notes, paste("ANOVA failed:", conditionMessage(out$anova)))
    out$anova <- NULL
  }
  out$note <- paste(unique(notes), collapse = "; ")
  out
}

summarise_results <- function(res, tx) {
  rows <- lapply(res, function(r) {
    a <- r$anova
    if (is.null(a) || !tx %in% rownames(a)) {
      return(data.frame(Gene = r$gene, F = NA, NumDF = NA, DenDF = NA,
                        p = NA, Note = r$note))
    }
    data.frame(Gene = r$gene, F = a[tx, "F"], NumDF = a[tx, "Df"],
               DenDF = a[tx, "Df.res"], p = a[tx, "Pr(>F)"], Note = r$note)
  })
  do.call(rbind, rows)
}

# ---- UI --------------------------------------------------------------------

ui <- page_sidebar(
  title = "Mixed Model Correction",
  sidebar = sidebar(
    width = 360,
    fileInput("file", "1. Upload data (.csv)", accept = c(".csv", "text/csv")),
    downloadLink("example", "Download an example file"),
    hr(),
    textAreaInput("genes", "2. Genes / parameters to test",
                  placeholder = "GRIA1, GRIA2, GRIN1, GAD1 ...", rows = 5),
    actionLink("all_numeric", "Fill with all numeric columns"),
    hr(),
    tags$b("3. Columns"),
    selectInput("tx", "Treatment (fixed effect)", choices = NULL),
    selectInput("line", "Line (random effect)", choices = NULL),
    selectInput("batch", "Batch (random effect, nested in Line)", choices = NULL),
    actionButton("run", "Run models", class = "btn-primary")
  ),
  navset_card_tab(
    nav_panel("Summary",
      uiOutput("formula_txt"),
      tableOutput("summary"),
      downloadButton("dl_summary", "Download summary (.csv)")
    ),
    nav_panel("Full ANOVA output", verbatimTextOutput("full")),
    nav_panel("Data preview", tableOutput("preview")),
    nav_panel("About",
      markdown("
Each gene/parameter is fit with a linear mixed model

`value ~ Tx + (1 | Line/Batch)`

using **lme4**, and Treatment is tested with
`car::Anova(type = 'II', test = 'F')` (Type II Wald F test with
Kenward-Roger degrees of freedom).

**Data format:** one row per sample, with columns for Treatment, Line and
Batch plus one numeric column per gene/parameter. Column names are converted
the same way R's `read.csv()` does, so `HLA-A` becomes `HLA.A`.
")
    )
  )
)

# ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {

  dat <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })

  observeEvent(dat(), {
    cols <- names(dat())
    updateSelectInput(session, "tx", choices = cols, selected = guess_col(cols, "Tx"))
    updateSelectInput(session, "line", choices = cols, selected = guess_col(cols, "Line"))
    updateSelectInput(session, "batch", choices = cols, selected = guess_col(cols, "Batch"))
  })

  observeEvent(input$all_numeric, {
    d <- dat()
    num <- names(d)[vapply(d, is.numeric, logical(1))]
    num <- setdiff(num, c(input$tx, input$line, input$batch))
    updateTextAreaInput(session, "genes", value = paste(num, collapse = ", "))
  })

  results <- eventReactive(input$run, {
    datos <- dat()
    MM_Vars <- parse_genes(input$genes)
    validate(need(length(MM_Vars) > 0, "Enter at least one gene/parameter."))
    factors <- c(input$tx, input$line, input$batch)
    validate(need(!anyDuplicated(factors),
                  "Treatment, Line and Batch must be different columns."))
    for (f in factors) datos[[f]] <- as.factor(datos[[f]])

    res <- withProgress(message = "Fitting models", value = 0, {
      lapply(seq_along(MM_Vars), function(i) {
        incProgress(1 / length(MM_Vars), detail = MM_Vars[i])
        fit_one(MM_Vars[i], datos, input$tx, input$line, input$batch)
      })
    })
    list(res = res, tx = input$tx,
         formula = sprintf("gene ~ %s + (1|%s/%s)", input$tx, input$line, input$batch))
  })

  summary_df <- reactive(summarise_results(results()$res, results()$tx))

  output$formula_txt <- renderUI({
    tags$p("Model: ", tags$code(results()$formula))
  })

  output$summary <- renderTable(summary_df(), digits = 4, na = "",
                                display = c("s", "s", "f", "d", "g", "g", "s"))

  output$full <- renderPrint({
    for (r in results()$res) {
      cat("\n")
      print(r$gene)
      if (!is.null(r$anova)) print(r$anova)
      if (nzchar(r$note)) cat("Note:", r$note, "\n")
      cat("\n")
    }
  })

  output$preview <- renderTable(head(dat(), 50))

  output$dl_summary <- downloadHandler(
    filename = function() paste0("mixed_model_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(summary_df(), file, row.names = FALSE)
  )

  output$example <- downloadHandler(
    filename = "example_data.csv",
    content = function(file) file.copy("example_data.csv", file)
  )
}

shinyApp(ui, server)
