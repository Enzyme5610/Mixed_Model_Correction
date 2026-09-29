library(shiny)
library(bslib)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)

# ---- Model (Kenward-Roger) ------------------------------

prepare_data <- function(path) {
  datos <- read.csv(path, header = TRUE, check.names = FALSE)
  original_names <- names(datos)
  if (!"Tx" %in% original_names) stop("The table must include a column named Tx.")
  names(datos) <- make.names(original_names, unique = TRUE)
  datos$Tx <- as.factor(datos$Tx)
  for (v in intersect(c("Line", "Batch"), names(datos))) datos[[v]] <- as.factor(datos[[v]])
  cols <- which(!original_names %in% c("Tx", "Line", "Batch") &
                nzchar(trimws(original_names)))  # skip blank headers
  list(datos = datos, MM_Vars = names(datos)[cols],
       parameter_labels = original_names[cols],
       numeric = vapply(datos[cols], is.numeric, logical(1)))
}

# Random effects by number of Lines/Batches, as in the script
random_term <- function(datos) {
  nl <- nlevels(datos$Line)
  nb <- nlevels(datos$Batch)
  if (nl > 1 && nb > 1) "(1|Line/Batch)"
  else if (nl > 1) "(1|Line)"
  else if (nb > 1) "(1|Batch)"
  else stop("Something is wrong with the number of Lines or Batches. Please check your data.")
}

run_models <- function(prep, vars, adjust, progress = function(n, label) NULL) {
  datos <- prep$datos
  rand <- random_term(datos)
  n <- length(vars)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    var <- vars[i]
    label <- prep$parameter_labels[match(var, prep$MM_Vars)]
    progress(n, label)
    f <- reformulate(c("Tx", rand), var)
    warn <- character()
    res <- withCallingHandlers(
      tryCatch({
        MM_Form <- lmerTest::lmer(f, data = datos, REML = TRUE)
        anova_result <- stats::anova(MM_Form, type = "II", ddf = "Kenward-Roger")
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
    res$rand <- rand
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

plot_defaults <- list(layout = "side", color = TRUE, size = 1, brackets = TRUE)

draw_plot <- function(datos, res, adjust, opt = plot_defaults) {
  y <- datos[[res$var]]
  lev <- levels(datos$Tx)
  k <- length(lev)
  lines <- levels(datos$Line)
  by_line <- opt$color && length(lines) > 0
  cols <- if (by_line) hcl.colors(length(lines), "Dark 3") else "grey45"
  grp <- if (by_line) as.integer(datos$Line) else 1
  means <- res$means
  pw <- res$pairs

  rng <- range(c(y, means$lower.CL, means$upper.CL), na.rm = TRUE)
  h <- diff(rng)
  if (h == 0) h <- 1
  n_pairs <- if (opt$brackets) nrow(pw) else 0
  top <- rng[2] + h * (0.05 + 0.1 * n_pairs)

  # Shrink text and margins below 5 in
  op <- par(cex = min(1, min(dev.size("in")) / 5),
            mar = c(3, 4.5, if (k > 2) 4.5 else 3.8, 7.5))
  on.exit(par(op))
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(rng[1] - 0.05 * h, top),
       xaxt = "n", xlab = "", ylab = res$label, las = 1)
  axis(1, at = seq_len(k), labels = lev)
  if (k > 2) {
    title(main = res$label, line = 2.6)
    mtext("Linear mixed model, Kenward-Roger", side = 3, line = 1.2, cex = 0.8 * par("cex"))
    mtext(paste0("Pairwise p-values: ", adjust_label(adjust, k), "-adjusted"),
          side = 3, line = 0.3, cex = 0.8 * par("cex"))
  } else {
    title(main = res$label, line = 1.6)
    mtext("Linear mixed model, Kenward-Roger", side = 3, line = 0.4, cex = 0.8 * par("cex"))
  }

  x <- as.integer(datos$Tx)
  xe <- match(as.character(means$Tx), lev)
  set.seed(1)
  if (opt$layout == "side") {
    xj <- x - 0.12 + runif(length(y), -0.08, 0.08)
    xe <- xe + 0.15
  } else {
    xj <- x + runif(length(y), -0.15, 0.15)
  }

  # Samples
  points(xj, y, pch = 19, cex = opt$size,
         col = adjustcolor(cols[grp], 0.75))

  # Model means ± 95% CI
  arrows(xe, means$lower.CL, xe, means$upper.CL,
         angle = 90, code = 3, length = 0.05, lwd = 2)
  points(xe, means$emmean, pch = 23, bg = "white", cex = 1.4, lwd = 2)

  # Brackets; emmeans pair order matches combn
  prs <- combn(k, 2)
  for (i in seq_len(n_pairs)) {
    a <- prs[1, i]; b <- prs[2, i]
    yb <- rng[2] + h * (0.06 + 0.1 * (i - 1))
    tick <- h * 0.02
    segments(c(a, a, b), c(yb - tick, yb, yb), c(a, b, b), c(yb, yb, yb - tick))
    text((a + b) / 2, yb, format_p(pw$p.value[i]), pos = 3, cex = 0.8, offset = 0.2)
  }

  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2])
  if (by_line) {
    legend(lx, usr[4], legend = lines, title = "Line", col = cols, pch = 19,
           bty = "n", xpd = TRUE, cex = 0.85)
  }
  legend(lx, usr[3] + 0.25 * diff(usr[3:4]),
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
    helpText("Needs columns named Tx, Line and Batch (exact spelling), in any position."),
    downloadLink("example", "Download an example file"),
    hr(),
    selectizeInput("params", "2. Parameters to analyze", choices = NULL,
                   multiple = TRUE, options = list(plugins = list("remove_button"))),
    helpText("Numeric columns are preselected."),
    hr(),
    radioButtons("adjust", "3. Pairwise p-value adjustment",
                 choices = c("Tukey" = "tukey", "Bonferroni" = "bonferroni")),
    helpText("With only two treatment groups there is a single comparison,",
             "so both give the same p-value."),
    actionButton("run", "4. Run models", class = "btn-primary")
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
    nav_panel("Plots", div(  # plain div: no fill layout
      layout_columns(
        col_widths = c(4, 4, 4),
        div(
          selectInput("plot_param", "Parameter", choices = NULL),
          radioButtons("layout", "Means", inline = TRUE,
                       choices = c("Beside dots" = "side", "Over dots" = "overlay"))
        ),
        div(
          checkboxInput("color_line", "Color dots by Line", TRUE),
          checkboxInput("brackets", "Show p-values", TRUE),
          sliderInput("pt_size", "Dot size", min = 0.4, max = 2, value = 1, step = 0.1)
        ),
        div(
          numericInput("w", "Width (in)", value = 5, min = 3, max = 12, step = 0.5),
          checkboxInput("square", "Square", TRUE),
          conditionalPanel("!input.square",
            numericInput("h", "Height (in)", value = 5, min = 3, max = 12, step = 0.5))
        )
      ),
      plotOutput("plot", width = "auto", height = "auto", fill = FALSE),
      div(
        downloadButton("dl_png", "PNG (this parameter)"),
        downloadButton("dl_pdf", "PDF (this parameter)"),
        downloadButton("dl_pdf_all", "PDF (all parameters)")
      )
    )),
    nav_panel("Data preview", tableOutput("preview")),
    nav_panel("About",
      markdown("
Each parameter is fit with a linear mixed model

`parameter ~ Tx + (1 | Line/Batch)`

using `lmerTest::lmer(REML = TRUE)`, and Tx is tested with
`anova(type = 'II', ddf = 'Kenward-Roger')`. With only one Line the model
uses `(1 | Batch)`; with only one Batch, `(1 | Line)`.

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

  observeEvent(prep(), {
    p <- prep()
    updateSelectizeInput(session, "params",
                         choices = setNames(p$MM_Vars, p$parameter_labels),
                         selected = p$MM_Vars[p$numeric])
  })

  results <- eventReactive(input$run, {
    p <- prep()
    validate(need(length(input$params) > 0, "Select at least one parameter."))
    res <- withProgress(message = "Fitting models", value = 0,
      tryCatch(run_models(p, input$params, input$adjust,
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

  # 3 sig. figs on screen; CSVs keep full precision
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
    labels <- vapply(r$res, `[[`, "", "label")
    cat("Parameters: ", paste(labels, collapse = ", "), "\n", sep = "")
    cat("Model: parameter ~ Tx + ", r$res[[1]]$rand, "\n", sep = "")
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

  opt <- reactive(list(layout = input$layout, color = input$color_line,
                       size = input$pt_size, brackets = input$brackets))

  # Size in inches, clamped to 3-12
  dims <- reactive({
    fit <- function(x) if (is.numeric(x) && !is.na(x)) min(max(x, 3), 12) else 5
    w <- fit(input$w)
    c(w = w, h = if (isTRUE(input$square)) w else fit(input$h))
  })

  output$plot <- renderPlot(
    draw_plot(results()$prep$datos, current(), results()$adjust, opt()),
    width = function() dims()[["w"]] * 96, height = function() dims()[["h"]] * 96,
    res = 96)

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
      d <- dims()
      plotPNG(function() draw_plot(results()$prep$datos, current(), results()$adjust, opt()),
              filename = file, width = d[["w"]] * 300, height = d[["h"]] * 300, res = 300)
    }
  )

  output$dl_pdf <- downloadHandler(
    filename = function() plot_name("pdf"),
    content = function(file) {
      pdf(file, width = dims()[["w"]], height = dims()[["h"]])
      draw_plot(results()$prep$datos, current(), results()$adjust, opt())
      dev.off()
    }
  )

  output$dl_pdf_all <- downloadHandler(
    filename = function() paste0(results()$base, "_MM_KR_plots.pdf"),
    content = function(file) {
      pdf(file, width = dims()[["w"]], height = dims()[["h"]])
      for (x in results()$res) draw_plot(results()$prep$datos, x, results()$adjust, opt())
      dev.off()
    }
  )

  output$example <- downloadHandler(
    filename = "example_data.csv",
    content = function(file) file.copy("example_data.csv", file)
  )
}

shinyApp(ui, server)
