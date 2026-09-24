library(shiny)
library(bslib)
library(lme4)
library(car)
library(emmeans)

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

random_terms <- function(choice, line, batch) {
  switch(choice,
    nested  = sprintf("(1|%s/%s)", line, batch),
    crossed = c(sprintf("(1|%s)", line), sprintf("(1|%s)", batch)),
    line    = sprintf("(1|%s)", line),
    batch   = sprintf("(1|%s)", batch)
  )
}

fit_one <- function(gene, dat, tx, rand) {
  out <- list(gene = gene, fit = NULL, anova = NULL, note = "")
  if (!gene %in% names(dat)) {
    out$note <- "Column not found in file"
    return(out)
  }
  if (!is.numeric(dat[[gene]])) {
    out$note <- "Column is not numeric"
    return(out)
  }
  f <- reformulate(c(tx, rand), gene)
  warn <- character()
  fit <- tryCatch(
    withCallingHandlers(
      suppressMessages(lmer(f, data = dat)),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out$note <- paste("Model failed:", conditionMessage(fit))
    return(out)
  }
  out$fit <- fit
  # Type II Wald F test with Kenward-Roger df (same as the original script).
  out$anova <- tryCatch(
    Anova(fit, type = "II", test.statistic = "F"),
    error = function(e) e
  )
  notes <- c(if (isSingular(fit)) "Singular fit (a random-effect variance is ~0)", warn)
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
  out <- do.call(rbind, rows)
  out$p_BH <- p.adjust(out$p, method = "BH")
  out[, c("Gene", "F", "NumDF", "DenDF", "p", "p_BH", "Note")]
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
    tags$b("3. Model"),
    selectInput("tx", "Treatment (fixed effect)", choices = NULL),
    selectInput("line", "Line (random effect)", choices = NULL),
    selectInput("batch", "Batch (random effect)", choices = NULL),
    radioButtons("rand", "Random-effect structure", choices = c(
      "Batch nested in Line: (1|Line/Batch)" = "nested",
      "Line and Batch crossed: (1|Line) + (1|Batch)" = "crossed",
      "Line only: (1|Line)" = "line",
      "Batch only: (1|Batch)" = "batch"
    )),
    checkboxInput("pairwise", "Pairwise treatment comparisons (emmeans, Tukey)", FALSE),
    actionButton("run", "Run models", class = "btn-primary")
  ),
  navset_card_tab(
    nav_panel("Summary",
      uiOutput("formula_txt"),
      tableOutput("summary"),
      downloadButton("dl_summary", "Download summary (.csv)")
    ),
    nav_panel("Full ANOVA output", verbatimTextOutput("full")),
    nav_panel("Pairwise comparisons", verbatimTextOutput("pairs")),
    nav_panel("Data preview", tableOutput("preview")),
    nav_panel("About",
      markdown("
Each gene/parameter is fit with a linear mixed model

`value ~ Treatment + (1 | Line/Batch)`

using **lme4**. Treatment is tested with a Type II Wald F test with
Kenward-Roger degrees of freedom (`car::Anova(type = 'II', test.statistic = 'F')`).
`p_BH` is the Benjamini-Hochberg adjusted p-value across all genes tested in the run.

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
    d <- dat()
    genes <- parse_genes(input$genes)
    validate(need(length(genes) > 0, "Enter at least one gene/parameter."))
    factors <- c(input$tx, input$line, input$batch)
    validate(need(!anyDuplicated(factors),
                  "Treatment, Line and Batch must be different columns."))
    for (f in factors) d[[f]] <- as.factor(d[[f]])
    rand <- random_terms(input$rand, input$line, input$batch)

    res <- withProgress(message = "Fitting models", value = 0, {
      lapply(seq_along(genes), function(i) {
        incProgress(1 / length(genes), detail = genes[i])
        fit_one(genes[i], d, input$tx, rand)
      })
    })
    list(res = res, tx = input$tx, pairwise = input$pairwise,
         formula = paste("gene ~", paste(c(input$tx, rand), collapse = " + ")))
  })

  summary_df <- reactive(summarise_results(results()$res, results()$tx))

  output$formula_txt <- renderUI({
    tags$p("Model: ", tags$code(results()$formula))
  })

  output$summary <- renderTable(summary_df(), digits = 4, na = "",
                                display = c("s", "s", "f", "d", "g", "g", "g", "s"))

  output$full <- renderPrint({
    for (r in results()$res) {
      cat("\n[1]", dQuote(r$gene, FALSE), "\n")
      if (!is.null(r$anova)) print(r$anova)
      if (nzchar(r$note)) cat("Note:", r$note, "\n")
    }
  })

  output$pairs <- renderPrint({
    rr <- results()
    if (!rr$pairwise) {
      cat("Tick 'Pairwise treatment comparisons' and re-run to see these.\n")
      return(invisible())
    }
    for (r in rr$res) {
      cat("\n====", r$gene, "====\n")
      if (is.null(r$fit)) { cat(r$note, "\n"); next }
      em <- tryCatch(emmeans(r$fit, rr$tx), error = function(e) e)
      if (inherits(em, "error")) { cat("emmeans failed:", conditionMessage(em), "\n"); next }
      print(em)
      print(pairs(em, adjust = "tukey"))
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
