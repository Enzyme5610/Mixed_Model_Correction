library(shiny)
library(bslib)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)

# ---- Model (Kenward-Roger) ------------------------------

prepare_data <- function(path) {
  datos <- read.csv(path, header = TRUE, check.names = FALSE, na.strings = c("", "NA"))
  original_names <- names(datos)
  if (!"Tx" %in% original_names) stop("The table must include a column named Tx.")
  names(datos) <- make.names(original_names, unique = TRUE)
  datos <- datos[!is.na(datos$Tx), ]  # drop empty rows
  datos$Tx <- as.factor(datos$Tx)
  for (v in intersect(c("Line", "Batch"), names(datos))) datos[[v]] <- as.factor(datos[[v]])
  cols <- which(!original_names %in% c("Tx", "Line", "Batch") &
                nzchar(trimws(original_names)))  # skip blank headers
  list(datos = datos, MM_Vars = names(datos)[cols],
       parameter_labels = original_names[cols],
       numeric = vapply(datos[cols], is.numeric, logical(1)))
}

# Random effects by design and number of Lines/Batches, as in the script
random_term <- function(datos, design = "nested") {
  nl <- nlevels(datos$Line)
  nb <- nlevels(datos$Batch)
  if (nl > 1 && nb > 1) {
    if (design == "crossed") return("(1|Line) + (1|Batch)")
    # One row per line x batch: nested can't be fit, reduces to (1|Line)
    one_each <- all(table(interaction(datos$Line, datos$Batch, drop = TRUE)) == 1)
    if (one_each) "(1|Line)" else "(1|Line/Batch)"
  }
  else if (nl > 1) "(1|Line)"
  else if (nb > 1) "(1|Batch)"
  else stop("Something is wrong with the number of Lines or Batches. Please check your data.")
}

run_models <- function(prep, vars, adjust, progress = function(n, label) NULL,
                       design = "nested") {
  datos <- prep$datos
  rand <- random_term(datos, design)
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
      Comparison = paste(levels(datos$Tx), collapse = " vs "),
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

plot_defaults <- list(type = "dots", layout = "side", dots = TRUE, color = TRUE,
                      size = 1, brackets = TRUE, scale = "raw", ref = NULL, ylab = "",
                      err = "ci", dir = "both", caps = TRUE, center = "diamond")

draw_plot <- function(datos, res, adjust, opt = plot_defaults) {
  opt <- modifyList(plot_defaults, Filter(Negate(is.null), opt))
  lev <- levels(datos$Tx)
  k <- length(lev)
  lines <- levels(datos$Line)
  by_line <- opt$color && length(lines) > 0
  cols <- if (by_line) hcl.colors(length(lines), "Dark 3") else "grey45"
  grp <- if (by_line) as.integer(datos$Line) else 1
  means <- res$means
  pw <- res$pairs

  # Reference group plotted first
  ref <- if (isTRUE(opt$ref %in% lev)) opt$ref else lev[1]
  ord <- c(ref, setdiff(lev, ref))
  pos <- match(lev, ord)
  xs <- pos[as.integer(datos$Tx)]
  xm <- pos[match(as.character(means$Tx), lev)]

  # Display scale; stats stay on entered values
  ref_mean <- means$emmean[as.character(means$Tx) == ref]
  scale <- if (opt$scale == "ratio" && ref_mean == 0) "raw" else opt$scale
  tf <- switch(scale, raw = identity,
               ratio = function(v) v / ref_mean,
               fc = function(v) 2^-(v - ref_mean))
  ylab <- switch(scale, raw = res$label,
                 ratio = paste("Relative to", ref),
                 fc = bquote("Fold change vs" ~ .(ref) ~ (2^{-Delta*Delta*Ct})))
  if (nzchar(trimws(opt$ylab))) ylab <- opt$ylab  # user label wins
  y <- tf(datos[[res$var]])

  # Center and error bars: model-based (CI, SE) or raw (SEM, SD)
  g <- as.character(means$Tx)
  if (opt$err %in% c("sem", "sd")) {
    m <- tapply(y, datos$Tx, mean, na.rm = TRUE)[g]
    s <- tapply(y, datos$Tx, sd, na.rm = TRUE)[g]
    if (opt$err == "sem") s <- s / sqrt(tapply(!is.na(y), datos$Tx, sum)[g])
    lo <- m - s; hi <- m + s
  } else {
    m <- tf(means$emmean)
    a <- if (opt$err == "se") means$emmean - means$SE else means$lower.CL
    b <- if (opt$err == "se") means$emmean + means$SE else means$upper.CL
    lo <- pmin(tf(a), tf(b)); hi <- pmax(tf(a), tf(b))
  }
  if (opt$dir == "up") lo <- m

  rng <- range(c(y, lo, hi, if (opt$type == "bar") 0), na.rm = TRUE)
  h <- diff(rng)
  if (h == 0) h <- 1
  n_pairs <- if (opt$brackets) nrow(pw) else 0
  top <- rng[2] + h * (0.05 + 0.1 * n_pairs)

  # Shrink text and margins below 5 in
  op <- par(cex = min(1, min(dev.size("in")) / 5),
            mar = c(3, 4.5, if (k > 2) 4.5 else 3.8, 7.5))
  on.exit(par(op))
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(rng[1] - 0.05 * h, top),
       xaxt = "n", xlab = "", ylab = ylab, las = 1)
  axis(1, at = seq_len(k), labels = ord)
  if (scale != "raw") abline(h = 1, lty = 3, col = "grey60")
  if (k > 2) {
    title(main = res$label, line = 2.6)
    mtext("Linear mixed model, Kenward-Roger", side = 3, line = 1.2, cex = 0.8 * par("cex"))
    mtext(paste0("Pairwise p-values: ", adjust_label(adjust, k), "-adjusted"),
          side = 3, line = 0.3, cex = 0.8 * par("cex"))
  } else {
    title(main = res$label, line = 1.6)
    mtext("Linear mixed model, Kenward-Roger", side = 3, line = 0.4, cex = 0.8 * par("cex"))
  }

  set.seed(1)
  side <- opt$type == "dots" && opt$layout == "side"
  xj <- xs + if (side) -0.12 + runif(length(y), -0.08, 0.08) else runif(length(y), -0.15, 0.15)
  xe <- xm + if (side) 0.15 else 0

  if (opt$type == "bar") {
    rect(xm - 0.3, 0, xm + 0.3, m, col = "grey88", border = "grey30")
  }
  if (opt$type == "violin") {
    for (g in seq_len(k)) {
      v <- y[xs == g & !is.na(y)]
      if (length(v) < 2 || diff(range(v)) == 0) next
      d <- density(v, from = min(v), to = max(v))
      w <- d$y / max(d$y) * 0.35
      polygon(c(g - w, rev(g + w)), c(d$x, rev(d$x)), col = "grey92", border = "grey45")
    }
  }

  # Samples
  show_dots <- opt$type == "dots" || opt$dots
  if (show_dots) {
    points(xj, y, pch = 19, cex = opt$size, col = adjustcolor(cols[grp], 0.75))
  }

  # Error bars and mean
  if (opt$caps) {
    suppressWarnings(arrows(xe, lo, xe, hi, angle = 90, length = 0.05, lwd = 2,
                            code = if (opt$dir == "up") 2 else 3))
  } else segments(xe, lo, xe, hi, lwd = 2)
  line_mark <- opt$type != "bar" && opt$center == "line"
  if (line_mark) segments(xe - 0.15, m, xe + 0.15, m, lwd = 3)
  else if (opt$type != "bar") points(xe, m, pch = 23, bg = "white", cex = 1.4, lwd = 2)

  # Brackets; emmeans pair order matches combn on model levels
  prs <- combn(k, 2)
  for (i in seq_len(n_pairs)) {
    a <- min(pos[prs[, i]]); b <- max(pos[prs[, i]])
    yb <- rng[2] + h * (0.06 + 0.1 * (i - 1))
    tick <- h * 0.02
    segments(c(a, a, b), c(yb - tick, yb, yb), c(a, b, b), c(yb, yb, yb - tick))
    text((a + b) / 2, yb, format_p(pw$p.value[i]), pos = 3, cex = 0.8, offset = 0.2)
  }

  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2])
  if (by_line && show_dots) {
    legend(lx, usr[4], legend = lines, title = "Line", col = cols, pch = 19,
           bty = "n", xpd = TRUE, cex = 0.85)
  }
  lab <- switch(opt$err, ci = "Model mean\n± 95% CI", se = "Model mean\n± SE",
                sem = "Mean ± SEM", sd = "Mean ± SD")
  key <- data.frame(lab = lab, pch = if (opt$type == "bar") 22 else if (line_mark) NA else 23,
                    bg = if (opt$type == "bar") "grey88" else "white", col = "black",
                    lty = if (line_mark) 1 else NA)
  if (opt$type == "violin") key <- rbind(data.frame(lab = "Distribution", pch = 22, bg = "grey92", col = "grey45", lty = NA), key)
  if (show_dots) key <- rbind(data.frame(lab = "Sample", pch = 19, bg = NA, col = "grey40", lty = NA), key)
  legend(lx, usr[3] + 0.3 * diff(usr[3:4]), legend = key$lab, pch = key$pch, lty = key$lty,
         lwd = 3, col = key$col, pt.bg = key$bg, pt.lwd = 1, bty = "n", xpd = TRUE,
         cex = 0.8, y.intersp = 1.4)
}

# ---- UI --------------------------------------------------------------------

ui <- page_sidebar(
  title = "Mixed Model Correction",
  sidebar = sidebar(
    width = 320,
    # Save files in the browser (Shinylive download links are unreliable)
    tags$script(HTML("
      Shiny.addCustomMessageHandler('save_file', function(m) {
        const bytes = Uint8Array.from(atob(m.data), c => c.charCodeAt(0));
        const url = URL.createObjectURL(new Blob([bytes], {type: m.type}));
        const a = document.createElement('a');
        a.href = url; a.download = m.name;
        document.body.appendChild(a); a.click(); a.remove();
        setTimeout(() => URL.revokeObjectURL(url), 5000);
      });")),
    fileInput("file", "1. Upload data (.csv)", accept = c(".csv", "text/csv")),
    helpText("Needs columns named Tx, Line and Batch (exact spelling), in any position."),
    actionLink("example", "Download an example file"),
    hr(),
    selectizeInput("params", "2. Parameters to analyze", choices = NULL,
                   multiple = TRUE, options = list(plugins = list("remove_button"))),
    helpText("Numeric columns are preselected."),
    hr(),
    radioButtons("design", "3. Do lines share batches?", choices = c(
      "No, each line has its own batches" = "nested",
      "Yes, all lines were run together in each batch" = "crossed")),
    helpText("Example of \"Yes\": L1-L4 all on the same qPCR plate labeled B1."),
    hr(),
    radioButtons("adjust", "4. Pairwise p-value adjustment",
                 choices = c("Tukey" = "tukey", "Bonferroni" = "bonferroni")),
    helpText("With only two treatment groups there is a single comparison,",
             "so both give the same p-value."),
    actionButton("run", "5. Run models", class = "btn-primary"),
    div(class = "small text-muted mt-3",
        "Original script: Dr. Luis Gustavo Hernandez Carballo", br(),
        "Shiny app and visualizations: Prachetas Jai Patel")
  ),
  navset_card_tab(
    nav_panel("ANOVA results",
      tableOutput("anova_table"),
      actionButton("dl_anova", "Download results (.csv)", icon = icon("download"))
    ),
    nav_panel("ANOVA output", verbatimTextOutput("anova_print")),
    nav_panel("Pairwise",
      tableOutput("pairs_table"),
      actionButton("dl_pairs", "Download pairwise (.csv)", icon = icon("download"))
    ),
    nav_panel("Plots", div(  # plain div: no fill layout
      layout_columns(
        col_widths = c(4, 4, 4),
        div(
          selectInput("plot_param", "Parameter", choices = NULL),
          radioButtons("type", "Plot type", inline = TRUE,
                       choices = c("Dots" = "dots", "Bar" = "bar", "Violin" = "violin")),
          conditionalPanel("input.type == 'dots'",
            radioButtons("layout", "Means", inline = TRUE,
                         choices = c("Beside dots" = "side", "Over dots" = "overlay"))),
          conditionalPanel("input.type != 'dots'",
            checkboxInput("dots", "Show dots", TRUE)),
          selectInput("err", "Error bars", choices = c(
            "95% CI (model)" = "ci", "SE (model)" = "se", "SEM" = "sem", "SD" = "sd")),
          conditionalPanel("input.err == 'sem' || input.err == 'sd'",
            helpText("SEM and SD use the raw values and ignore Line and Batch.")),
          radioButtons("dir", NULL, inline = TRUE,
                       choices = c("Both directions" = "both", "Above only" = "up")),
          checkboxInput("caps", "Caps", TRUE),
          conditionalPanel("input.type != 'bar'",
            radioButtons("center", "Mean marker", inline = TRUE,
                         choices = c("Diamond" = "diamond", "Line" = "line")))
        ),
        div(
          selectInput("scale", "Y axis", choices = c(
            "Values as entered" = "raw",
            "Relative to reference (linear data)" = "ratio",
            "Fold change 2^-ΔΔCt (ΔCt data)" = "fc")),
          textInput("ylab", "Y-axis label (optional)", placeholder = "Name (units)"),
          selectInput("ref", "Reference group", choices = NULL),
          checkboxInput("color_line", "Color dots by Line", TRUE),
          checkboxInput("brackets", "Show p-values", TRUE)
        ),
        div(
          sliderInput("pt_size", "Dot size", min = 0.4, max = 2, value = 1, step = 0.1),
          numericInput("w", "Width (in)", value = 5, min = 3, max = 12, step = 0.5),
          checkboxInput("square", "Square", TRUE),
          conditionalPanel("!input.square",
            numericInput("h", "Height (in)", value = 5, min = 3, max = 12, step = 0.5))
        )
      ),
      plotOutput("plot", width = "auto", height = "auto", fill = FALSE),
      div(
        actionButton("dl_png", "PNG (this parameter)", icon = icon("download")),
        actionButton("dl_pdf", "PDF (this parameter)", icon = icon("download")),
        actionButton("dl_pdf_all", "PDF (all parameters)", icon = icon("download"))
      )
    )),
    nav_panel("Data preview", tableOutput("preview")),
    nav_panel("About",
      markdown("
Each parameter is fit with a linear mixed model

`parameter ~ Tx + (1 | Line/Batch)`

using `lmerTest::lmer(REML = TRUE)`, and Tx is tested with
`anova(type = 'II', ddf = 'Kenward-Roger')`. With only one Line the model
uses `(1 | Batch)`; with only one Batch, `(1 | Line)`. With one row per
line and batch (no replicates), the nested term can't be estimated and the
model reduces to `(1 | Line)`.

If lines share batches (all lines run together in each batch, e.g. on the
same qPCR plate), answer \"Yes\" in step 3 to use the crossed model
`parameter ~ Tx + (1 | Line) + (1 | Batch)`, which removes run-wide shifts
shared by all lines.

**Pairwise comparisons** between treatments use estimated marginal means
from the same model (`emmeans`, Kenward-Roger df), with Tukey or Bonferroni
adjustment as selected before running.

**Plots** show each sample (colored by Line) and the model's estimated mean
with 95% confidence interval for each treatment, with pairwise p-values, as
dots, bars or violins. The Y axis can show values relative to a reference
group: as a ratio for linear data (e.g. electrophysiology), or as fold change
2^-ΔΔCt when the values are ΔCt (qPCR). Statistics always use the values as
entered.

**Error bars** can show the model's 95% CI or SE (matching the statistics),
or the SEM or SD of the raw values (which ignore Line and Batch).

---

**Credits**

- Original R script: **Dr. Luis Gustavo Hernandez Carballo**
- Shiny app and visualizations: **Prachetas Jai Patel**

If you use this tool in a publication, poster or presentation, please
acknowledge both authors.
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
                          function(n, label) incProgress(1 / n, detail = label),
                          design = input$design),
               error = function(e) validate(conditionMessage(e))))
    list(res = res, prep = p, adjust = input$adjust,
         base = tools::file_path_sans_ext(input$file$name))
  })

  observeEvent(results(), {
    labels <- vapply(results()$res, `[[`, "", "label")
    updateSelectInput(session, "plot_param", choices = labels)
    lev <- levels(results()$prep$datos$Tx)
    ctrl <- grep("^(control|ctrl|gfp)", lev, ignore.case = TRUE, value = TRUE)
    updateSelectInput(session, "ref", choices = lev, selected = c(ctrl, lev)[1])
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

  opt <- reactive(list(type = input$type, layout = input$layout, dots = input$dots,
                       color = input$color_line, size = input$pt_size,
                       brackets = input$brackets, scale = input$scale, ref = input$ref,
                       ylab = input$ylab, err = input$err, dir = input$dir,
                       caps = input$caps, center = input$center))

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

  save_file <- function(name, type, write) {
    tmp <- tempfile()
    write(tmp)
    session$sendCustomMessage("save_file", list(
      name = name, type = type,
      data = jsonlite::base64_enc(readBin(tmp, "raw", file.size(tmp)))))
  }

  # Results or NULL (with a hint) if models haven't run
  ready <- function() {
    r <- tryCatch(if (input$run > 0) results(), error = function(e) NULL)
    if (is.null(r)) showNotification("Run the models first.", type = "warning")
    r
  }

  observeEvent(input$dl_anova, {
    r <- ready(); req(r)
    save_file(paste0(r$base, "_MM_KR_results.csv"), "text/csv",
              function(f) write.csv(anova_df(), f, row.names = FALSE))
  })

  observeEvent(input$dl_pairs, {
    r <- ready(); req(r)
    save_file(paste0(r$base, "_MM_KR_pairwise.csv"), "text/csv",
              function(f) write.csv(pairs_df(), f, row.names = FALSE))
  })

  plot_name <- function(r, ext) paste0(r$base, "_", make.names(input$plot_param), ".", ext)

  observeEvent(input$dl_png, {
    r <- ready(); req(r)
    d <- dims()
    save_file(plot_name(r, "png"), "image/png", function(f) {
      plotPNG(function() draw_plot(r$prep$datos, current(), r$adjust, opt()),
              filename = f, width = d[["w"]] * 300, height = d[["h"]] * 300, res = 300)
    })
  })

  observeEvent(input$dl_pdf, {
    r <- ready(); req(r)
    save_file(plot_name(r, "pdf"), "application/pdf", function(f) {
      pdf(f, width = dims()[["w"]], height = dims()[["h"]])
      draw_plot(r$prep$datos, current(), r$adjust, opt())
      dev.off()
    })
  })

  observeEvent(input$dl_pdf_all, {
    r <- ready(); req(r)
    save_file(paste0(r$base, "_MM_KR_plots.pdf"), "application/pdf", function(f) {
      pdf(f, width = dims()[["w"]], height = dims()[["h"]])
      for (x in r$res) draw_plot(r$prep$datos, x, r$adjust, opt())
      dev.off()
    })
  })

  observeEvent(input$example, {
    save_file("example_data.csv", "text/csv",
              function(f) file.copy("example_data.csv", f))
  })
}

shinyApp(ui, server)
