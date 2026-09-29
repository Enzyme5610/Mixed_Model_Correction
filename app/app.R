library(shiny)
library(bslib)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)

# ---- Model (colleague's Kenward-Roger script) ------------------------------

# Reads and checks the CSV exactly as the script does. Returns the data with
# formula-safe names plus the original parameter labels, or stops with the
# script's error message.
prepare_data <- function(path) {
  datos <- read.csv(path, header = TRUE, check.names = FALSE)
  if (ncol(datos) < 4L) stop("The CSV must contain at least four columns.")
  original_names <- names(datos)
  parameter_columns <- seq.int(4L, ncol(datos))
  parameter_labels <- original_names[parameter_columns]
  if (any(!nzchar(trimws(parameter_labels)))) {
    stop("Each parameter column (column 4 onward) must have a header.")
  }
  if (!setequal(original_names[1:3], c("Tx", "Line", "Batch"))) {
    stop("The first three columns must be named Tx, Line, and Batch (in any order).")
  }
  # Make headers safe for formulas while retaining original labels in output.
  names(datos) <- make.names(original_names, unique = TRUE)
  datos$Tx <- as.factor(datos$Tx)
  datos$Line <- as.factor(datos$Line)
  datos$Batch <- as.factor(datos$Batch)
  list(datos = datos, MM_Vars = names(datos)[parameter_columns],
       parameter_labels = parameter_labels)
}

run_models <- function(prep, adjust, progress = function(n, label) NULL) {
  datos <- prep$datos
  n <- length(prep$MM_Vars)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    var <- prep$MM_Vars[i]
    label <- prep$parameter_labels[i]
    progress(n, label)
    f <- reformulate(c("Tx", "(1|Line/Batch)"), var)
    warn <- character()
    res <- withCallingHandlers(
      tryCatch({
        MM_Form <- lmerTest::lmer(f, data = datos, REML = TRUE)
        # Type II F tests using Kenward-Roger degrees of freedom.
        anova_result <- stats::anova(MM_Form, type = "II", ddf = "Kenward-Roger")

        # Pairwise comparisons between treatments from the same model.
        em <- emmeans(MM_Form, specs = "Tx", lmer.df = "kenward-roger")
        means <- as.data.frame(summary(em))
        pw <- as.data.frame(summary(pairs(em, adjust = adjust)))
        list(MM_Form = MM_Form, anova = anova_result, means = means, pairs = pw)
      }, error = function(e) {
        stop("Model failed for '", label, "': ", conditionMessage(e), call. = FALSE)
      }),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    result_table <- as.data.frame(res$anova)
    res$table <- data.frame(
      Parameter = label,
      Term = rownames(result_table),
      DF_method = "Kenward-Roger",
      result_table,
      row.names = NULL,
      check.names = FALSE
    )
    res$pairs_table <- data.frame(
      Parameter = label,
      res$pairs,
      DF_method = "Kenward-Roger",
      P_adjust = adjust_label(adjust, nlevels(datos$Tx)),
      check.names = FALSE
    )
    res$var <- var
    res$label <- label
    res$warnings <- unique(warn)
    out[[i]] <- res
  }
  out
}

adjust_label <- function(adjust, k) {
  if (k <= 2) return("none (one comparison)")
  c(tukey = "Tukey", bonferroni = "Bonferroni")[[adjust]]
}

# ---- Plot --------------------------------------------------------------------

format_p <- function(p) {
  ifelse(p < 0.0001, "p < 0.0001", paste("p =", signif(p, 3)))
}

# Raw values per sample (colored by Line) with the model's estimated mean and
# 95% CI for each treatment, and pairwise p-values as brackets.
draw_plot <- function(datos, res, adjust) {
  y <- datos[[res$var]]
  lev <- levels(datos$Tx)
  k <- length(lev)
  lines <- levels(datos$Line)
  cols <- hcl.colors(length(lines), "Dark 3")
  means <- res$means
  pw <- res$pairs

  rng <- range(c(y, means$lower.CL, means$upper.CL), na.rm = TRUE)
  h <- diff(rng)
  if (h == 0) h <- 1
  n_pairs <- nrow(pw)
  top <- rng[2] + h * (0.08 + 0.1 * n_pairs)

  op <- par(mar = c(4.5, 4.5, 5, 8.5))
  on.exit(par(op))
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(rng[1] - 0.05 * h, top),
       xaxt = "n", xlab = "Tx", ylab = res$label, las = 1)
  title(main = res$label, line = 2.8)
  axis(1, at = seq_len(k), labels = lev)
  mtext("Tx + (1|Line/Batch), Kenward-Roger df", side = 3, line = 1.2, cex = 0.8)
  mtext(paste("Pairwise adjustment:", adjust_label(adjust, k)),
        side = 3, line = 0.3, cex = 0.8)

  # Raw data, jittered, left of each group center.
  set.seed(1)
  xj <- as.integer(datos$Tx) - 0.12 + runif(length(y), -0.08, 0.08)
  points(xj, y, pch = 19, col = adjustcolor(cols[as.integer(datos$Line)], 0.8))

  # Estimated marginal means with 95% CI, right of each group center.
  xe <- match(as.character(means$Tx), lev) + 0.15
  arrows(xe, means$lower.CL, xe, means$upper.CL,
         angle = 90, code = 3, length = 0.05, lwd = 2)
  points(xe, means$emmean, pch = 23, bg = "white", cex = 1.4, lwd = 2)

  # Pairwise brackets (emmeans orders pairs like combn: 1-2, 1-3, ..., 2-3, ...).
  prs <- combn(k, 2)
  for (i in seq_len(n_pairs)) {
    a <- prs[1, i]; b <- prs[2, i]
    yb <- rng[2] + h * (0.06 + 0.1 * (i - 1))
    tick <- h * 0.02
    segments(c(a, a, b), c(yb - tick, yb, yb), c(a, b, b), c(yb, yb, yb - tick))
    text((a + b) / 2, yb, format_p(pw$p.value[i]), pos = 3, cex = 0.8, offset = 0.2)
  }

  usr <- par("usr")
  legend(usr[2] + 0.02 * diff(usr[1:2]), usr[4], legend = lines, title = "Line",
         col = cols, pch = 19, bty = "n", xpd = TRUE, cex = 0.9)
  legend(usr[2] + 0.02 * diff(usr[1:2]), usr[3] + 0.25 * diff(usr[3:4]),
         legend = c("Sample", "Model mean\n± 95% CI"), pch = c(19, 23),
         col = c("grey40", "black"), pt.bg = "white", bty = "n", xpd = TRUE,
         cex = 0.8, y.intersp = 1.4)
}

# ---- UI --------------------------------------------------------------------

ui <- page_sidebar(
  title = "Mixed Model Correction",
  sidebar = sidebar(
    width = 320,
    fileInput("file", "1. Upload data (.csv)", accept = c(".csv", "text/csv")),
    helpText("First three columns: Tx, Line, Batch (any order).",
             "Every column from the 4th onward is analyzed as a parameter."),
    downloadLink("example", "Download an example file"),
    hr(),
    radioButtons("adjust", "2. Pairwise p-value adjustment",
                 choices = c("Tukey" = "tukey", "Bonferroni" = "bonferroni")),
    helpText("With only two treatment groups there is a single comparison,",
             "so both give the same p-value."),
    actionButton("run", "3. Run models", class = "btn-primary")
  ),
  navset_card_tab(
    nav_panel("ANOVA results",
      tableOutput("anova_table"),
      downloadButton("dl_anova", "Download results (.csv)")
    ),
    nav_panel("ANOVA output", verbatimTextOutput("anova_print")),
    nav_panel("Pairwise",
      tableOutput("pairs_table"),
      downloadButton("dl_pairs", "Download pairwise (.csv)")
    ),
    nav_panel("Plots",
      selectInput("plot_param", "Parameter", choices = NULL),
      plotOutput("plot", height = "480px"),
      div(
        downloadButton("dl_png", "PNG (this parameter)"),
        downloadButton("dl_pdf", "PDF (this parameter)"),
        downloadButton("dl_pdf_all", "PDF (all parameters)")
      )
    ),
    nav_panel("Data preview", tableOutput("preview")),
    nav_panel("About",
      markdown("
Each parameter is fit with a linear mixed model

`parameter ~ Tx + (1 | Line/Batch)`

using `lmerTest::lmer(REML = TRUE)`, and Tx is tested with
`anova(type = 'II', ddf = 'Kenward-Roger')`.

**Pairwise comparisons** between treatments use estimated marginal means
from the same model (`emmeans`, Kenward-Roger df), with Tukey or Bonferroni
adjustment as selected before running.

**Plots** show each sample (colored by Line) and the model's estimated mean
with 95% confidence interval for each treatment, with pairwise p-values.
")
    )
  )
)

# ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {

  prep <- reactive({
    req(input$file)
    tryCatch(prepare_data(input$file$datapath),
             error = function(e) validate(conditionMessage(e)))
  })

  results <- eventReactive(input$run, {
    p <- prep()
    res <- withProgress(message = "Fitting models", value = 0,
      tryCatch(run_models(p, input$adjust,
                          function(n, label) incProgress(1 / n, detail = label)),
               error = function(e) validate(conditionMessage(e))))
    list(res = res, prep = p, adjust = input$adjust,
         base = tools::file_path_sans_ext(input$file$name))
  })

  observeEvent(results(), {
    labels <- vapply(results()$res, `[[`, "", "label")
    updateSelectInput(session, "plot_param", choices = labels)
  })

  anova_df <- reactive(do.call(rbind, lapply(results()$res, `[[`, "table")))
  pairs_df <- reactive(do.call(rbind, lapply(results()$res, `[[`, "pairs_table")))

  # Show p-values to 3 significant figures on screen; the CSVs keep full precision.
  show_p <- function(df) {
    for (col in intersect(c("Pr(>F)", "p.value"), names(df))) {
      df[[col]] <- as.character(signif(df[[col]], 3))
    }
    df
  }
  output$anova_table <- renderTable(show_p(anova_df()), digits = 4)
  output$pairs_table <- renderTable(show_p(pairs_df()), digits = 4)

  output$anova_print <- renderPrint({
    r <- results()
    cat("Parameters: ", paste(r$prep$parameter_labels, collapse = ", "), "\n", sep = "")
    for (x in r$res) {
      cat("\n", x$label, "\n", sep = "")
      print(x$anova)
      if (length(x$warnings)) cat("Warning:", x$warnings, sep = "\n  ")
      cat("\n")
    }
  })

  current <- reactive({
    r <- results()
    req(input$plot_param)
    r$res[[match(input$plot_param, vapply(r$res, `[[`, "", "label"))]]
  })

  output$plot <- renderPlot(draw_plot(results()$prep$datos, current(), results()$adjust),
                            res = 110)

  output$preview <- renderTable(head(prep()$datos, 50))

  output$dl_anova <- downloadHandler(
    filename = function() paste0(results()$base, "_MM_KR_results.csv"),
    content = function(file) write.csv(anova_df(), file, row.names = FALSE)
  )

  output$dl_pairs <- downloadHandler(
    filename = function() paste0(results()$base, "_MM_KR_pairwise.csv"),
    content = function(file) write.csv(pairs_df(), file, row.names = FALSE)
  )

  plot_name <- function(ext) {
    paste0(results()$base, "_", make.names(input$plot_param), ".", ext)
  }

  output$dl_png <- downloadHandler(
    filename = function() plot_name("png"),
    content = function(file) {
      plotPNG(function() draw_plot(results()$prep$datos, current(), results()$adjust),
              filename = file, width = 2100, height = 1500, res = 300)
    }
  )

  output$dl_pdf <- downloadHandler(
    filename = function() plot_name("pdf"),
    content = function(file) {
      pdf(file, width = 7, height = 5)
      draw_plot(results()$prep$datos, current(), results()$adjust)
      dev.off()
    }
  )

  output$dl_pdf_all <- downloadHandler(
    filename = function() paste0(results()$base, "_MM_KR_plots.pdf"),
    content = function(file) {
      pdf(file, width = 7, height = 5)
      for (x in results()$res) draw_plot(results()$prep$datos, x, results()$adjust)
      dev.off()
    }
  )

  output$example <- downloadHandler(
    filename = "example_data.csv",
    content = function(file) file.copy("example_data.csv", file)
  )
}

shinyApp(ui, server)
