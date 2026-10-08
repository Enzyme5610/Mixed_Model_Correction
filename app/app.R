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
  # Levels in CSV order
  for (v in intersect(c("Tx", "Line", "Batch"), names(datos)))
    datos[[v]] <- factor(datos[[v]], levels = unique(datos[[v]]))
  cols <- which(!original_names %in% c("Tx", "Line", "Batch") &
                nzchar(trimws(original_names)))  # skip blank headers
  # Metadata (dates, IDs, wells) and constant columns aren't preselected
  meta <- grepl("date|day|well|passage|(^|[^a-z])id($|[^a-z])", original_names[cols], ignore.case = TRUE) |
    vapply(datos[cols], function(x) length(unique(na.omit(x))) < 2, TRUE)
  list(datos = datos, MM_Vars = names(datos)[cols],
       parameter_labels = original_names[cols],
       numeric = vapply(datos[cols], is.numeric, logical(1)) & !meta)
}

# Random effects by design and number of Lines/Batches, as in the script
random_term <- function(datos, design = "nested") {
  nl <- nlevels(datos$Line)
  nb <- nlevels(datos$Batch)
  if (nl > 1 && nb > 1) {
    # Crossed needs a batch shared by lines
    if (design == "crossed" && length(shared_batches(datos))) return("(1|Line) + (1|Batch)")
    # Nested needs a line with repeat batches and replicates within them
    repeats <- any(tapply(datos$Batch, datos$Line, function(b) length(unique(b))) > 1, na.rm = TRUE)
    one_each <- all(table(interaction(datos$Line, datos$Batch, drop = TRUE)) == 1)
    if (repeats && !one_each) "(1|Line/Batch)" else "(1|Line)"
  }
  else if (nl > 1) "(1|Line)"
  else if (nb > 1) "(1|Batch)"
  else stop("Something is wrong with the number of Lines or Batches. Please check your data.")
}

run_models <- function(prep, vars, adjust, progress = function(n, label) NULL,
                       design = "nested") {
  datos <- prep$datos
  n <- length(vars)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    var <- vars[i]
    label <- prep$parameter_labels[match(var, prep$MM_Vars)]
    progress(n, label)
    # Lines/Batches counted where this parameter was measured
    rand <- tryCatch(random_term(droplevels(datos[!is.na(datos[[var]]), ]), design),
                     error = function(e) stop("'", label, "': ", conditionMessage(e), call. = FALSE))
    f <- reformulate(c("Tx", rand), var)
    warn <- character()
    res <- withCallingHandlers(
      tryCatch({
        MM_Form <- lmerTest::lmer(f, data = datos, REML = TRUE)
        anova_result <- stats::anova(MM_Form, type = "II", ddf = "Kenward-Roger")
        em <- emmeans(MM_Form, specs = "Tx", lmer.df = "kenward-roger")
        means <- as.data.frame(summary(em))
        pw <- as.data.frame(summary(pairs(em, adjust = adjust), infer = TRUE))
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
    # Possible outliers: scaled model residual beyond 3
    r <- residuals(res$MM_Form, type = "pearson", scaled = TRUE)
    res$out <- data.frame(i = match(names(r), rownames(datos)), resid = unname(r))
    res$out <- res$out[abs(res$out$resid) > 3, ]
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

# Bracket label: p-value or GraphPad Prism-style stars
sig_label <- function(p, opt) {
  if (opt$labels != "stars") return(format_p(p))
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                   c("****", "***", "**", "*", "ns"), right = FALSE))
}
# Bracket label; asterisks sit high in the line, so lower them to match "ns"
bracket_text <- function(x, y, p, opt, cex) {
  lab <- sig_label(p, opt)
  star <- opt$labels == "stars" && lab != "ns"
  text(x, y - if (star) 0.35 * strheight("*", cex = cex * 1.3) else 0, lab,
       pos = 3, offset = 0.2, cex = if (star) cex * 1.3 else cex)
}

# Star key as a legend block from y
star_key <- function(opt, x, y) {
  if (opt$labels != "stars" || !opt$brackets) return(y)
  y <- y - strheight("M", cex = 0.75)  # gap below the key above
  g <- legend(x, y, legend = c("p ≥ 0.05", "p < 0.05", "p < 0.01", "p < 0.001", "p < 0.0001"),
              title = "Significance", title.adj = 0, pch = NA, x.intersp = 3.2,
              bty = "n", xpd = TRUE, cex = 0.75)
  text(g$rect$left + 0.1 * g$rect$w, g$text$y, c("ns", "*", "**", "***", "****"),
       adj = c(0, 0.5), xpd = TRUE, cex = 0.75)
  g$rect$top - g$rect$h
}

# Group order: user's order if valid, else reference first then CSV order
tx_order <- function(datos, opt) {
  lev <- levels(datos$Tx)
  ord <- unlist(opt$order)
  if (length(ord) == length(lev) && setequal(ord, lev)) return(ord)
  ref <- if (isTRUE(opt$ref %in% lev)) opt$ref else lev[1]
  c(ref, setdiff(lev, ref))
}

shape_set <- c("● Circle" = 19, "■ Square" = 15, "▲ Triangle" = 17, "◆ Diamond" = 18,
               "○ Open circle" = 1, "□ Open square" = 0, "△ Open triangle" = 2,
               "◇ Open diamond" = 5, "▽ Open down triangle" = 6, "× Cross" = 4,
               "+ Plus" = 3, "✱ Star" = 8)

# Per-level color or shape by Line, Batch or Group; user picks override defaults
dot_style <- function(datos, opt, by, what = "col") {
  f <- switch(by, line = datos$Line, batch = datos$Batch, group = datos$Tx, NULL)
  if (is.null(f) || nlevels(f) == 0) return(NULL)
  lev <- if (by == "group") tx_order(datos, list(ref = opt$ref)) else levels(f)
  n <- length(lev)
  def <- if (what == "pch") rep_len(shape_set, n)
         else if (by == "group") c("#666666", hcl.colors(max(n - 1, 1), "Dark 3"))[seq_len(n)]
         else hcl.colors(n, "Dark 3")
  picks <- if (what == "pch") opt$shapes else opt$cols
  val <- vapply(seq_len(n), function(i) {
    p <- picks[[paste0(by, ":", lev[i])]]
    as.character(if (is.null(p)) def[i] else p)
  }, "")
  if (what == "pch") val <- as.integer(val)
  list(lev = lev, val = val, idx = match(as.character(f), lev),
       title = c(line = "Line", batch = "Batch", group = "Group")[[by]])
}

# Opaque tint toward white (EMF has no transparency)
tint <- function(col, a) rgb(t(1 - a * (1 - col2rgb(col) / 255)))

# Dot color and shape legends, stacked from y
dot_legend <- function(x, y, cs, ss) {
  same <- identical(cs$title, ss$title)
  keys <- list(
    if (!is.null(cs)) list(s = cs, col = cs$val, pch = if (same) ss$val else 19),
    if (!is.null(ss) && !same) list(s = ss, col = "grey30", pch = ss$val))
  for (k in Filter(Negate(is.null), keys)) {
    g <- legend(x, y, legend = k$s$lev, title = k$s$title, title.adj = 0, col = k$col,
                pch = k$pch, bty = "n", xpd = TRUE, cex = 0.85)
    y <- g$rect$top - g$rect$h
  }
  y
}

plot_defaults <- list(type = "bar", layout = "side", dots = TRUE, color_by = "line",
                      shape_by = "none", cols = list(), shapes = list(), order = NULL,
                      labels = "stars",
                      size = 1, brackets = TRUE, scale = "raw", ref = NULL, ylab = "",
                      err = "sem", center = "diamond",
                      rot = "auto", adj_batch = FALSE, outliers = TRUE)

# Each row's estimated batch shift from the model (0 if no batch term)
batch_shift <- function(datos, res) {
  re <- lme4::ranef(res$MM_Form)
  key <- if ("Batch" %in% names(re)) list(re$Batch, as.character(datos$Batch))
         else if ("Batch:Line" %in% names(re))
           list(re[["Batch:Line"]], paste(datos$Batch, datos$Line, sep = ":"))
  if (is.null(key)) return(NULL)
  key[[1]][key[[2]], 1]
}

# Display values, center and error bars for one parameter
plot_stats <- function(datos, res, opt) {
  lev <- levels(datos$Tx)
  means <- res$means
  ref <- if (isTRUE(opt$ref %in% lev)) opt$ref else lev[1]  # plotted first
  ref_mean <- means$emmean[as.character(means$Tx) == ref]
  # Ratios need positive values; otherwise direction can flip (e.g. negative dCt)
  ratio_ok <- ref_mean > 0 && all(datos[[res$var]] > 0, na.rm = TRUE)
  scale <- if (opt$scale == "ratio" && !ratio_ok) "raw" else opt$scale
  # Display scale; stats stay on entered values
  tf <- switch(scale, raw = identity,
               ratio = function(v) v / ref_mean,
               fc = function(v) 2^-(v - ref_mean))
  raw <- datos[[res$var]]
  shift <- if (opt$adj_batch) batch_shift(datos, res)
  y <- tf(if (is.null(shift)) raw else raw - shift)  # display only
  g <- as.character(means$Tx)
  # Model-based (CI, SE) or raw (SEM, SD)
  if (opt$err %in% c("sem", "sd", "none")) {
    m <- tapply(y, datos$Tx, mean, na.rm = TRUE)[g]
    s <- tapply(y, datos$Tx, sd, na.rm = TRUE)[g]
    if (opt$err == "sem") s <- s / sqrt(tapply(!is.na(y), datos$Tx, sum)[g])
    if (opt$err == "none") s <- 0
    lo <- m - s; hi <- m + s
  } else {
    m <- tf(means$emmean)
    a <- if (opt$err == "se") means$emmean - means$SE else means$lower.CL
    b <- if (opt$err == "se") means$emmean + means$SE else means$upper.CL
    lo <- pmin(tf(a), tf(b)); hi <- pmax(tf(a), tf(b))
  }
  if (opt$type == "bar") lo <- m  # bars: error bar above only
  list(ref = ref, ord = tx_order(datos, opt), scale = scale, y = y,
       g = g, m = m, lo = lo, hi = hi, adjusted = !is.null(shift))
}

model_note <- function(opt, adjusted) {
  paste0("Linear mixed model, Kenward-Roger",
         if (opt$adj_batch) if (adjusted) "; batch-adjusted values"
                            else "; batch adjustment not possible (no batch term)")
}

scale_note <- function(opt, scale) {
  if (opt$scale == "ratio" && scale == "raw") {
    mtext("Relative scale needs positive values; showing values as entered (use fold change for ΔCt)",
          side = 1, line = par("mar")[1] - 1, adj = 0, cex = 0.6 * par("cex"))
  }
}

scale_label <- function(st, raw, opt) {
  if (nzchar(trimws(opt$ylab))) return(opt$ylab)  # user label wins
  switch(st$scale, raw = raw,
         ratio = paste("Relative to", st$ref),
         fc = bquote("Fold change vs" ~ .(st$ref) ~ (2^{-Delta*Delta*Ct})))
}

draw_err <- function(xe, st, opt, half = 0.15) {
  if (opt$err == "none") return()
  suppressWarnings(arrows(xe, st$lo, xe, st$hi, angle = 90, length = 0.05, lwd = 2,
                          code = if (opt$type == "bar") 2 else 3))
  if (opt$type == "bar") return()
  if (opt$center == "line") segments(xe - half, st$m, xe + half, st$m, lwd = 3)
  else points(xe, st$m, pch = 23, bg = "white", cex = 1.4, lwd = 2)
}

# X labels at 0, 45 or 90 degrees; "auto" picks 90 when crowded
label_angle <- function(labels, opt) {
  if (opt$rot != "auto") return(as.numeric(opt$rot))
  # Estimated label width vs space per label (inches); R drops overlapping labels
  cex <- min(1, min(dev.size("in")) / 5)
  slot <- (dev.size("in")[1] - 12 * par("cin")[2] * cex) / length(labels)
  if (max(nchar(labels)) * 0.75 * par("cin")[1] * cex > 0.9 * slot) 90 else 0
}
bottom_mar <- function(labels, angle) {
  if (angle == 0) 3 else 1.5 + max(nchar(labels)) * if (angle == 90) 0.65 else 0.5
}
x_labels <- function(at, labels, angle) {
  if (angle == 45) {
    axis(1, at = at, labels = FALSE)
    text(at, par("usr")[3] - 0.04 * diff(par("usr")[3:4]), labels,
         srt = 45, adj = 1, xpd = TRUE)
  } else axis(1, at = at, labels = labels, las = if (angle == 90) 2 else 1)
}

# Red ring around flagged cells
ring_outliers <- function(x, y, idx, opt) {
  if (opt$outliers && length(idx)) points(x[idx], y[idx], pch = 1, cex = opt$size * 2,
                                          col = "red", lwd = 1.5)
}
outlier_row <- data.frame(lab = "Possible outlier", pch = 1, lty = NA, bg = NA, col = "red")

# Legend rows for the summary marks: box = box, line = line
type_key <- function(opt) {
  row <- function(lab, pch = NA, lty = NA, bg = NA, col = "black")
    data.frame(lab = lab, pch = pch, lty = lty, bg = bg, col = col)
  line_mark <- opt$center == "line"
  err <- opt$err != "none"
  mean_row <- if (err) row(err_label(opt), pch = if (line_mark) NA else 23,
                           lty = if (line_mark) 1 else NA, bg = "white")
  switch(opt$type,
    bar = rbind(row(if (opt$err %in% c("ci", "se")) "Model mean" else "Mean", pch = 22, bg = "grey88"),
                if (err) row(sub("^.*(±)", "\\1", err_label(opt)), lty = 1)),
    box = rbind(row("Median", lty = 1, col = "grey30"), row("IQR", pch = 22, bg = "grey88", col = "grey30"), mean_row),
    violin = rbind(row("Distribution", pch = 22, bg = "grey92", col = "grey45"), mean_row),
    mean_row)
}

err_label <- function(opt) {
  switch(opt$err, ci = "Model mean\n± 95% CI", se = "Model mean\n± SE",
         sem = "Mean ± SEM", sd = "Mean ± SD", none = "Mean")
}

draw_plot <- function(datos, res, adjust, opt = plot_defaults) {
  opt <- modifyList(plot_defaults, Filter(Negate(is.null), opt))
  lev <- levels(datos$Tx)
  k <- length(lev)
  cs <- dot_style(datos, opt, opt$color_by)
  ss <- dot_style(datos, opt, opt$shape_by, "pch")
  gs <- dot_style(datos, opt, "group")
  pw <- res$pairs

  st <- plot_stats(datos, res, opt)
  ord <- st$ord
  pos <- match(lev, ord)
  xs <- pos[as.integer(datos$Tx)]
  xm <- pos[match(st$g, lev)]
  scale <- st$scale
  ylab <- scale_label(st, res$label, opt)
  y <- st$y; m <- st$m; lo <- st$lo; hi <- st$hi

  rng <- range(c(y, lo, hi, if (opt$type == "bar") 0), na.rm = TRUE)
  h <- diff(rng)
  if (h == 0) h <- 1
  n_pairs <- if (opt$brackets) nrow(pw) else 0
  top <- rng[2] + h * (0.05 + 0.1 * n_pairs)

  # Shrink text and margins below 5 in
  ang <- label_angle(ord, opt)
  op <- par(cex = min(1, min(dev.size("in")) / 5), tcl = -0.25, mgp = c(3, 0.6, 0),
            mar = c(bottom_mar(ord, ang), 4.5, if (k > 2) 4.5 else 3.8, 7.5))
  on.exit(par(op))
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(rng[1] - 0.05 * h, top),
       xaxt = "n", xlab = "", ylab = ylab, las = 1)
  x_labels(seq_len(k), ord, ang)
  if (scale != "raw") abline(h = 1, lty = 3, col = "grey60")
  if (k > 2) {
    title(main = res$label, line = 2.6)
    mtext(model_note(opt, st$adjusted), side = 3, line = 1.2, cex = 0.8 * par("cex"))
    mtext(paste0("Pairwise p-values: ", adjust_label(adjust, k), "-adjusted"),
          side = 3, line = 0.3, cex = 0.8 * par("cex"))
  } else {
    title(main = res$label, line = 1.6)
    mtext(model_note(opt, st$adjusted), side = 3, line = 0.4, cex = 0.8 * par("cex"))
  }

  set.seed(1)
  side <- opt$type == "dots" && opt$layout == "side"
  xj <- xs + if (side) -0.12 + runif(length(y), -0.08, 0.08) else runif(length(y), -0.15, 0.15)
  xe <- xm + if (side) 0.15 else if (opt$type == "box") 0.3 else 0  # box: mean beside

  if (opt$type == "bar") {
    rect(xm - 0.3, 0, xm + 0.3, m, border = "grey30",
         col = tint(gs$val[match(st$g, gs$lev)], 0.3))
  }
  if (opt$type == "violin") {
    for (g in seq_len(k)) {
      v <- y[xs == g & !is.na(y)]
      if (length(v) < 2 || diff(range(v)) == 0) next
      d <- density(v, from = min(v), to = max(v))
      w <- d$y / max(d$y) * 0.35
      polygon(c(g - w, rev(g + w)), c(d$x, rev(d$x)), border = "grey45",
              col = tint(gs$val[match(ord[g], gs$lev)], 0.2))
    }
  }
  if (opt$type == "box") {
    boxplot(split(y, factor(xs, levels = seq_len(k))), at = seq_len(k), add = TRUE,
            axes = FALSE, outline = FALSE, boxwex = 0.4, border = "grey30",
            col = tint(gs$val[match(ord, gs$lev)], 0.3))
  }

  # Samples
  show_dots <- opt$type == "dots" || opt$dots
  if (show_dots) {
    col <- if (is.null(cs)) "grey45" else cs$val[cs$idx]
    points(xj, y, pch = if (is.null(ss)) 19 else ss$val[ss$idx], cex = opt$size,
           col = tint(col, 0.75))
    ring_outliers(xj, y, res$out$i, opt)
  }

  draw_err(xe, st, opt)

  # Brackets; emmeans pair order matches combn on model levels
  prs <- combn(k, 2)
  for (i in seq_len(n_pairs)) {
    a <- min(pos[prs[, i]]); b <- max(pos[prs[, i]])
    yb <- rng[2] + h * (0.06 + 0.1 * (i - 1))
    tick <- h * 0.02
    segments(c(a, a, b), c(yb - tick, yb, yb), c(a, b, b), c(yb, yb, yb - tick))
    bracket_text((a + b) / 2, yb, pw$p.value[i], opt, 0.8)
  }
  scale_note(opt, scale)

  # Right-hand column: dots, then mean/error key, then stars
  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2])
  ky <- if (show_dots) dot_legend(lx, usr[4], cs, ss) else usr[4]
  key <- type_key(opt)
  if (show_dots) key <- rbind(data.frame(lab = "Sample", pch = 19, lty = NA, bg = NA, col = "grey40"), key)
  if (show_dots && opt$outliers && nrow(res$out)) key <- rbind(key, outlier_row)
  if (NROW(key)) {
    g <- legend(lx, ky, legend = key$lab, pch = key$pch, lty = key$lty,
                lwd = 3, col = key$col, pt.bg = key$bg, pt.lwd = 1, bty = "n", xpd = TRUE,
                cex = 0.8, y.intersp = 1.4)
    ky <- g$rect$top - g$rect$h
  }
  star_key(opt, lx, ky)
}

# All parameters in one figure: parameters on x, groups side by side
draw_multi <- function(datos, results, adjust, opt = plot_defaults) {
  opt <- modifyList(plot_defaults, Filter(Negate(is.null), opt))
  n <- length(results)
  sts <- lapply(results, function(r) plot_stats(datos, r, opt))
  # Shared axis: if any parameter can't use the relative scale, none do
  if (opt$scale == "ratio" && any(vapply(sts, function(s) s$scale == "raw", TRUE)))
    sts <- lapply(results, function(r) plot_stats(datos, r, modifyList(opt, list(scale = "raw"))))
  ord <- sts[[1]]$ord
  k <- length(ord)
  w <- 0.8 / k
  off <- (seq_len(k) - (k + 1) / 2) * w
  gs <- dot_style(datos, opt, "group")
  cols <- gs$val[match(ord, gs$lev)]
  gi <- match(as.character(datos$Tx), ord)
  cs <- if (opt$color_by != "group") dot_style(datos, opt, opt$color_by)  # groups have a legend
  ss <- dot_style(datos, opt, opt$shape_by, "pch")
  labels <- vapply(results, `[[`, "", "label")

  all_v <- unlist(lapply(sts, function(s) c(s$y, s$lo, s$hi)))
  rng <- range(c(all_v, if (opt$type == "bar") 0), na.rm = TRUE)
  h <- diff(rng)
  if (h == 0) h <- 1
  n_pairs <- if (opt$brackets) nrow(results[[1]]$pairs) else 0
  top <- rng[2] + h * (0.05 + 0.1 * n_pairs)

  ang <- label_angle(labels, opt)
  op <- par(cex = min(1, min(dev.size("in")) / 5), tcl = -0.25, mgp = c(3, 0.6, 0),
            mar = c(bottom_mar(labels, ang), 4.5, if (k > 2) 4.5 else 3.8, 7.5))
  on.exit(par(op))
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(rng[1] - 0.05 * h, top), xaxt = "n",
       xlab = "", ylab = scale_label(sts[[1]], "Value", opt), las = 1)
  x_labels(seq_len(n), labels, ang)
  if (sts[[1]]$scale != "raw") abline(h = 1, lty = 3, col = "grey60")
  title(main = paste(ord, collapse = " vs "), line = if (k > 2) 2.6 else 1.6)
  mtext(model_note(opt, all(vapply(sts, `[[`, TRUE, "adjusted"))), side = 3, line = if (k > 2) 1.2 else 0.4,
        cex = 0.8 * par("cex"))
  if (k > 2) mtext(paste0("Pairwise p-values: ", adjust_label(adjust, k), "-adjusted"),
                   side = 3, line = 0.3, cex = 0.8 * par("cex"))

  set.seed(1)
  show_dots <- opt$type == "dots" || opt$dots
  for (i in seq_len(n)) {
    st <- sts[[i]]
    xm <- i + off[match(st$g, ord)]
    xs <- i + off[gi]
    if (opt$type == "bar") {
      rect(xm - w * 0.4, 0, xm + w * 0.4, st$m, border = "grey30",
           col = tint(cols[match(st$g, ord)], 0.3))
    }
    if (opt$type == "violin") {
      for (j in seq_len(k)) {
        v <- st$y[gi == j & !is.na(st$y)]
        if (length(v) < 2 || diff(range(v)) == 0) next
        d <- density(v, from = min(v), to = max(v))
        hw <- d$y / max(d$y) * w * 0.45
        x0 <- i + off[j]
        polygon(c(x0 - hw, rev(x0 + hw)), c(d$x, rev(d$x)), border = "grey45",
                col = tint(cols[j], 0.2))
      }
    }
    if (opt$type == "box") {
      boxplot(split(st$y, factor(gi, levels = seq_len(k))), at = i + off, add = TRUE,
              axes = FALSE, outline = FALSE, boxwex = w * 0.5, border = "grey30",
              col = tint(cols, 0.3))
    }
    if (show_dots) {
      col <- if (opt$color_by == "group") cols[gi] else if (is.null(cs)) "grey45" else cs$val[cs$idx]
      xd <- xs + runif(length(xs), -w * 0.25, w * 0.25)
      points(xd, st$y, pch = if (is.null(ss)) 19 else ss$val[ss$idx],
             cex = opt$size, col = tint(col, 0.75))
      ring_outliers(xd, st$y, results[[i]]$out$i, opt)
    }
    draw_err(xm + if (opt$type == "box") w * 0.35 else 0, st, opt, half = w * 0.3)

    # Brackets at one shared height; emmeans pair order matches combn on model levels
    pw <- results[[i]]$pairs
    prs <- combn(levels(datos$Tx), 2)
    for (p in seq_len(n_pairs)) {
      xa <- sort(i + off[match(prs[, p], ord)])
      yb <- rng[2] + h * (0.04 + 0.08 * (p - 1))
      segments(c(xa[1], xa[1], xa[2]), c(yb - h * 0.015, yb, yb),
               c(xa[1], xa[2], xa[2]), c(yb, yb, yb - h * 0.015))
      bracket_text(mean(xa), yb, pw$p.value[p], opt, 0.65)
    }
  }
  scale_note(opt, sts[[1]]$scale)

  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2])
  g <- legend(lx, usr[4], legend = ord, title = "Group", title.adj = 0, pch = 22, pt.cex = 1.6,
              pt.bg = tint(cols, 0.5), col = cols, bty = "n", xpd = TRUE, cex = 0.85)
  ky <- g$rect$top - g$rect$h
  if (show_dots) ky <- dot_legend(lx, ky, cs, ss)
  key <- type_key(opt)
  if (show_dots && opt$outliers && any(vapply(results, function(r) nrow(r$out) > 0, TRUE)))
    key <- rbind(key, outlier_row)
  if (NROW(key)) {
    g <- legend(lx, ky, legend = key$lab, pch = key$pch, lty = key$lty, lwd = 3,
                col = key$col, pt.bg = key$bg, pt.lwd = 1, bty = "n", xpd = TRUE,
                cex = 0.8, y.intersp = 1.4)
    ky <- g$rect$top - g$rect$h
  }
  star_key(opt, lx, ky)
}

# Batch labels used by more than one line (crossed needs one)
shared_batches <- function(d) {
  if (nlevels(d$Line) < 2 || nlevels(d$Batch) < 2) return(list())
  Filter(function(x) length(x) > 1, lapply(split(as.character(d$Line), d$Batch), unique))
}

# Rows per Line x Batch, split by Tx
design_table <- function(d) {
  n <- table(d$Line, d$Batch, d$Tx)
  cell <- apply(n, 1:2, function(x) paste(paste(names(x)[x > 0], x[x > 0]), collapse = ", "))
  data.frame(Line = rownames(cell), cell, check.names = FALSE, row.names = NULL)
}

# Where each parameter was measured
coverage <- function(p, vars) {
  do.call(rbind, lapply(vars, function(v) {
    s <- droplevels(p$datos[!is.na(p$datos[[v]]), ])
    data.frame(Parameter = p$parameter_labels[match(v, p$MM_Vars)], n = nrow(s),
               Lines = paste(levels(s$Line), collapse = ", "),
               Batches = paste(levels(s$Batch), collapse = ", "),
               Testable = nlevels(s$Line) > 1 || nlevels(s$Batch) > 1, check.names = FALSE)
  }))
}

# Mini spreadsheet; color = one run
batch_sheet <- function(rows, caption, formula) {
  cols <- c(g = "#0F6E56", p = "#534AB7", o = "#D85A30", b = "#185FA5",
            y = "#854F0B", k = "#993556")
  tagList(
    tags$table(class = "mmc-sheet",
      tags$tr(tags$th("Tx"), tags$th("Line"), tags$th("Batch")),
      lapply(rows, function(r) tags$tr(style = paste0("background:", cols[[r[4]]]),
        tags$td(r[1]), tags$td(r[2]), tags$td(r[3])))),
    tags$small(class = "text-muted d-block", caption),
    tags$small(class = "text-muted font-monospace d-block", formula))
}

# Welcome tab: steps, batches, designs and fallbacks
welcome_page <- function() {
  step <- function(n, title, sub) div(class = "mmc-step",
    strong(n, " ", title), tags$small(class = "text-muted d-block", sub))
  arrow <- span(class = "mmc-arrow", "→")
  fb <- function(data, model) tags$tr(tags$td(data), tags$td(model))
  div(class = "p-2", style = "max-width: 900px;",
    h3("Welcome to Mixed Model Correction!"),
    p("Tests treatment effects on cells grouped in lines and batches. Replaces the",
      "t-test and one-way ANOVA for this kind of data. Runs in your browser; your",
      "data stays on your computer."),
    h5(class = "mt-4", "How to use"),
    div(class = "mmc-steps",
      step("1", "Upload CSV", "Tx, Line, Batch + parameters"), arrow,
      step("2", "Parameters", "Numeric columns preselected"), arrow,
      step("3", "Batch design", "Nested unless lines shared a batch"), arrow,
      step("4", "Adjustment", "Tukey or Bonferroni; pick your control"), arrow,
      step("5", "Run models", "Results and plots open in the tabs")),
    h5(class = "mt-4", "Batch design (step 3)"),
    p("A batch is one run (one color below): one ephys recording day or one qPCR plate.",
      "Pick the design by how the cells were run. One file = one experiment."),
    div(class = "row g-3",
      div(class = "col-md-6", div(class = "border rounded p-2 h-100",
        strong("Nested"), p(class = "small mb-1",
          "Each line run on its own. Number each line's runs B1, B2…"),
        batch_sheet(list(
          c("Control", "L1", "B1", "g"), c("Control", "L1", "B2", "p"),
          c("AD", "L2", "B1", "o"), c("AD", "L2", "B2", "b"),
          c("Control", "L3", "B1", "y"), c("AD", "L4", "B1", "k")),
          "L2's B1 is not L1's B1. L3 and L4 have one batch each.",
          "(1|Line/Batch)"))),
      div(class = "col-md-6", div(class = "border rounded p-2 h-100",
        strong("Crossed"), p(class = "small mb-1",
          "Lines run together. Number the shared runs B1, B2…"),
        batch_sheet(list(
          c("Control", "L1", "B1", "g"), c("AD", "L2", "B1", "g"),
          c("Control", "L3", "B1", "g"), c("Control", "L1", "B2", "p"),
          c("AD", "L2", "B2", "p"), c("AD", "L4", "B2", "p")),
          "B1 is one day or plate for every line.",
          "(1|Line) + (1|Batch)")))),
    p(class = "small text-muted", "Same line in both groups (KD, OE, drug)? Use the same line",
      "label; either design works."),
    h5(class = "mt-4", "Which model runs"),
    tags$table(class = "table table-sm w-auto",
      tags$tr(tags$th("Your data"), tags$th("Model used")),
      fb("Crossed, lines share a batch", tags$code("(1|Line) + (1|Batch)")),
      fb("Nested, some lines with repeat batches", tags$code("(1|Line/Batch)")),
      fb("Several lines, one batch each", tags$code("(1|Line)")),
      fb("One line, several batches", tagList(tags$code("(1|Batch)"), " (that line only)")),
      fb("One line, one batch", "Can't be tested")),
    h5(class = "mt-4", "Not covered"),
    tags$ul(class = "small",
      tags$li("Data: small counts, percentages near 0 or 100%, scores, omics, skewed data with zeros"),
      tags$li("Designs: two batch sources, extra levels (e.g. coverslip), genotype × treatment,",
              "repeated measures per cell, matched pairs")),
    p(class = "small text-muted", "Full details: ",
      tags$a(href = "https://github.com/Enzyme5610/Mixed_Model_Correction", "README")))
}

# Plot download formats (raster at 300 dpi)
raster <- function(dev, ...) function(f, w, h) dev(f, width = w, height = h, units = "in", res = 300, ...)
plot_formats <- list(
  png  = list(label = "PNG", ext = "png", type = "image/png", open = raster(png)),
  tiff = list(label = "TIFF", ext = "tiff", type = "image/tiff",
              open = raster(tiff, compression = "lzw")),
  jpeg = list(label = "JPEG", ext = "jpg", type = "image/jpeg", open = raster(jpeg, quality = 95)),
  pdf  = list(label = "PDF", ext = "pdf", type = "application/pdf", open = pdf),
  svg  = list(label = "SVG", ext = "svg", type = "image/svg+xml", open = svg),
  # Plain EMF (no EMF+) so PowerPoint can ungroup it into shapes
  emf  = list(label = "EMF", ext = "emf", type = "image/emf",
              open = function(f, w, h) devEMF::emf(f, w, h, emfPlus = FALSE, family = "Aptos")),
  pdf_all = list(type = "application/pdf", open = pdf))

dl_item <- function(f, label) {
  tags$li(tags$a(class = "dropdown-item mmc-dl", href = "#", `data-fmt` = f, label))
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
      });
      // Color and shape pickers, kept per session
      const picks = {cols: {}, shapes: {}};
      document.addEventListener('change', e => {
        const el = e.target.closest('.mmc-pick'); if (!el) return;
        picks[el.dataset.input][el.dataset.key] = el.value;
        Shiny.setInputValue(el.dataset.input, picks[el.dataset.input]);
      });
      // External links open in a new tab
      document.addEventListener('click', e => {
        const a = e.target.closest('a[href^=http]');
        if (a) { a.target = '_blank'; a.rel = 'noopener'; }
      });
      // Plot download menu
      document.addEventListener('click', e => {
        const a = e.target.closest('.mmc-dl'); if (!a) return;
        e.preventDefault();
        Shiny.setInputValue('dl_plot', a.dataset.fmt, {priority: 'event'});
      });
      // Group order: drag (mouse or touch), or arrow buttons
      let drag = null;
      const sendOrder = ul => Shiny.setInputValue('tx_order',
        [...ul.children].map(li => li.dataset.lev), {priority: 'event'});
      document.addEventListener('pointerdown', e => {
        const li = e.target.closest('.mmc-order li');
        if (!li || e.target.closest('input, button')) return;
        drag = li; li.classList.add('active'); e.preventDefault();
      });
      document.addEventListener('pointermove', e => {
        const li = drag && document.elementFromPoint(e.clientX, e.clientY)?.closest('.mmc-order li');
        if (!li || li === drag || li.parentNode !== drag.parentNode) return;
        const r = li.getBoundingClientRect();
        li.parentNode.insertBefore(drag, e.clientY > r.top + r.height / 2 ? li.nextSibling : li);
      });
      document.addEventListener('pointerup', () => {
        if (!drag) return;
        drag.classList.remove('active'); sendOrder(drag.parentNode); drag = null;
      });
      document.addEventListener('click', e => {
        const b = e.target.closest('.mmc-move'); if (!b) return;
        const li = b.closest('li'), ul = li.parentNode;
        if (b.dataset.dir === 'up' && li.previousElementSibling) ul.insertBefore(li, li.previousElementSibling);
        if (b.dataset.dir === 'down' && li.nextElementSibling) ul.insertBefore(li.nextElementSibling, li);
        sendOrder(ul);
      });")),
    tags$style(".mmc-order li { cursor: grab; touch-action: none; } .mmc-pick[type=color] { width: 2em; height: 1.6em; padding: 0; border: 0; }
      .mmc-sheet { border-collapse: collapse; table-layout: fixed; width: 100%; margin: 4px 0 2px; font: .75rem monospace; }
      .mmc-sheet th { background: #f1f1f1; color: #555; border: 1px solid #ccc; padding: 1px 6px; font-weight: 400; }
      .mmc-sheet td { color: #fff; font-weight: 600; border: 1px solid #fff; padding: 1px 6px; }
      .mmc-steps { display: flex; flex-wrap: wrap; gap: .4rem; align-items: stretch; }
      .mmc-step { border: 1px solid #dee2e6; border-radius: .5rem; padding: .4rem .6rem; flex: 1 1 110px; max-width: 160px; }
      .mmc-arrow { align-self: center; font-size: 1.4rem; color: #888; }
      #design .radio { margin-bottom: .6rem; }"),
    fileInput("file", "1. Upload data (.csv)", accept = c(".csv", "text/csv")),
    helpText("Needs columns named Tx, Line and Batch (exact spelling), in any position."),
    actionLink("example", "Download an example file"),
    hr(),
    selectizeInput("params", "2. Parameters to analyze", choices = NULL,
                   multiple = TRUE, options = list(plugins = list("remove_button"))),
    helpText("Numeric columns are preselected."),
    uiOutput("param_warning"),
    hr(),
    radioButtons("design", "3. Batch design",
      choices = c("Nested" = "nested", "Crossed" = "crossed")),
    uiOutput("design_note"),
    hr(),
    radioButtons("adjust", "4. Pairwise p-value adjustment",
                 choices = c("Tukey" = "tukey", "Bonferroni" = "bonferroni")),
    helpText("With only two treatment groups there is a single comparison,",
             "so both give the same p-value."),
    selectInput("ref", "Reference (control) group", choices = NULL),
    helpText("Sets comparison direction in tables and plots; p-values don't change."),
    actionButton("run", "5. Run models", class = "btn-primary"),
    div(class = "small text-muted mt-3",
        "Original script: Dr. Luis Gustavo Hernandez Carballo", br(),
        "Shiny app and visualizations: Prachetas Jai Patel")
  ),
  navset_card_tab(
    id = "tabs",
    nav_panel("Welcome", welcome_page()),
    nav_panel("ANOVA results",
      helpText("Fold-change columns appear when the plot Y axis is set to fold change",
               "(ΔCt data); confidence intervals are in the Pairwise tab."),
      tableOutput("anova_table"),
      actionButton("dl_anova", "Download results (.csv)", icon = icon("download")),
      h6(class = "mt-4", "Possible outliers"),
      helpText("Cells more than 3 SD from the model's prediction for their group, line and",
               "batch. Nothing is removed; check these cells and edit the CSV if needed."),
      tableOutput("outlier_table"),
      checkboxInput("show_resid", "Show residual checks", FALSE),
      conditionalPanel("input.show_resid",
        selectInput("resid_param", "Parameter", choices = NULL),
        plotOutput("resid_plot", width = "640px", height = "320px", fill = FALSE),
        helpText("Grey band: 95% range expected if the model fits. Left: points inside",
                 "the band, residuals roughly normal. Right: even spread around 0, equal",
                 "variance. Red: possible outlier."))
    ),
    nav_panel("ANOVA output", verbatimTextOutput("anova_print")),
    nav_panel("Pairwise",
      helpText("Fold-change columns (2^-estimate) appear when the plot Y axis is set",
               "to fold change (ΔCt data)."),
      tableOutput("pairs_table"),
      actionButton("dl_pairs", "Download pairwise (.csv)", icon = icon("download"))
    ),
    nav_panel("Plots", layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(width = 300, accordion(
        open = "1. What to plot",
        accordion_panel("1. What to plot",
          radioButtons("fig", NULL, choices = c(
            "One parameter" = "one", "All parameters in one figure" = "all")),
          conditionalPanel("input.fig == 'one'",
            selectInput("plot_param", "Parameter", choices = NULL)),
          conditionalPanel("input.fig == 'all'",
            selectizeInput("multi", "Parameters", choices = NULL, multiple = TRUE,
                           options = list(plugins = list("remove_button")))),
          radioButtons("type", "Plot type", inline = TRUE, selected = "bar",
                       choices = c("Dots" = "dots", "Bar" = "bar", "Box" = "box",
                                   "Violin" = "violin")),
          conditionalPanel("input.type == 'dots'",
            radioButtons("layout", "Mean and error bar", inline = TRUE,
                         choices = c("Beside dots" = "side", "Over dots" = "overlay"))),
          conditionalPanel("input.type != 'dots'",
            checkboxInput("dots", "Show dots", TRUE)),
          checkboxInput("outliers", "Circle possible outliers", TRUE),
          checkboxInput("adj_batch", "Batch-adjusted values", FALSE),
          conditionalPanel("input.adj_batch",
            helpText("Dots minus each batch's estimated shift. Statistics are",
                     "unchanged; don't re-test adjusted values."))),
        accordion_panel("2. Group order & colors",
          helpText("Drag a group (or use the arrows) to reorder. Click a swatch to",
                   "change its color."),
          uiOutput("group_ui")),
        accordion_panel("3. Dot colors & shapes",
          selectInput("color_by", "Color dots by", choices = c(
            "Line" = "line", "Batch" = "batch", "Group (set in 2)" = "group", "None" = "none")),
          selectInput("shape_by", "Shape dots by", choices = c(
            "None" = "none", "Line" = "line", "Batch" = "batch")),
          uiOutput("dot_style_ui"),
          sliderInput("pt_size", "Dot size", min = 0.4, max = 2, value = 1, step = 0.1)),
        accordion_panel("4. Mean & error bars",
          selectInput("err", "Error bars", selected = "sem", choices = c(
            "95% CI (model)" = "ci", "SE (model)" = "se", "SEM" = "sem", "SD" = "sd",
            "None" = "none")),
          conditionalPanel("input.err == 'sem' || input.err == 'sd'",
            helpText("SEM and SD use the raw values and ignore Line and Batch.")),
          conditionalPanel("input.type != 'bar' && input.err != 'none'",
            radioButtons("center", "Mean marker", inline = TRUE,
                         choices = c("Diamond" = "diamond", "Line" = "line")))),
        accordion_panel("5. Axes & significance",
          selectInput("scale", "Y axis", choices = c(
            "Values as entered" = "raw",
            "Relative to reference (linear data)" = "ratio",
            "Fold change 2^-ΔΔCt (ΔCt data)" = "fc")),
          helpText("Relative: for positive measurements (e.g. amplitude). For ΔCt, use fold change."),
          textInput("ylab", "Y-axis label (optional)", placeholder = "Name (units)"),
          selectInput("rot", "X label angle", choices = c(
            "Auto" = "auto", "Horizontal" = "0", "45°" = "45", "Vertical" = "90")),
          radioButtons("sig", "Significance", inline = TRUE, selected = "stars",
                       choices = c("p-values" = "p", "Stars" = "stars", "Hide" = "none"))),
        accordion_panel("6. Figure size",
          numericInput("w", "Width (in)", value = 5, min = 3, max = 12, step = 0.5),
          checkboxInput("square", "Square", TRUE),
          conditionalPanel("!input.square",
            numericInput("h", "Height (in)", value = 5, min = 3, max = 12, step = 0.5)))
      )),
      plotOutput("plot", width = "auto", height = "auto", fill = FALSE),
      div(class = "dropdown",
        tags$button(type = "button", class = "btn btn-default dropdown-toggle",
                    `data-bs-toggle` = "dropdown", icon("download"), "Download plot"),
        tags$ul(class = "dropdown-menu", lapply(names(plot_formats), function(f) {
          if (f == "pdf_all") return(tagList(tags$li(tags$hr(class = "dropdown-divider")),
                                              dl_item(f, "PDF, all parameters")))
          dl_item(f, plot_formats[[f]]$label)
        })))
    )),
    nav_panel("Data preview",
      h6("Design: rows per line and batch"), tableOutput("design_tbl"),
      h6("Where each parameter was measured"), tableOutput("coverage_tbl"),
      h6("First 50 rows"), tableOutput("preview")),
    nav_panel("About",
      markdown("
**Model.** Each parameter: `parameter ~ Tx + (1 | Line/Batch)`, fit with
`lmerTest::lmer` (REML); Tx tested by Type II F test, Kenward-Roger df.
Crossed (step 3, lines shared a batch): `(1 | Line) + (1 | Batch)`.
One Line: `(1 | Batch)`. One Batch, or no line with repeat batches:
`(1 | Line)`.

**FDR.** q-values (Benjamini-Hochberg) for panels of many similar
parameters (e.g. gene panels). Adjusts Tx p-values across all parameters in
a run.

**Pairwise.** `emmeans` from the same model, Tukey or Bonferroni adjusted.
The reference group sets direction only; p-values don't change.

**Plots.** Display only; statistics always use the values as entered.
Model CI/SE match the statistics; SEM/SD ignore Line and Batch.
Batch-adjusted values subtract each batch's estimated shift; don't re-test
them.

---

**Credits**

- Original R script: **Dr. Luis Gustavo Hernandez Carballo**
- Shiny app and visualizations: **Prachetas Jai Patel**

Please acknowledge both authors if you use this tool. Loading animation:
[loading-bar](https://github.com/loadingio/loading-bar) (MIT). References:
[README](https://github.com/Enzyme5610/Mixed_Model_Correction#references).
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
    lev <- levels(p$datos$Tx)
    ctrl <- grep("^(control|ctrl|gfp)", lev, ignore.case = TRUE, value = TRUE)
    updateSelectInput(session, "ref", choices = lev, selected = c(ctrl, lev)[1])
  })

  # Box plots default to no mean/error bar; restore SEM when leaving box
  prev_type <- reactiveVal("bar")
  observeEvent(input$type, {
    if (input$type == "box") updateSelectInput(session, "err", selected = "none")
    else if (prev_type() == "box" && identical(input$err, "none"))
      updateSelectInput(session, "err", selected = "sem")
    prev_type(input$type)
  }, ignoreInit = TRUE)

  # Model that will run, with fallbacks
  output$design_note <- renderUI({
    if (is.null(input$file)) return(NULL)
    d <- prep()$datos
    f <- tryCatch(random_term(d, input$design), error = function(e) NULL)
    div(class = "alert alert-info small py-1 px-2 mb-2",
      sprintf("%d line(s), %d batch(es). Model: ", nlevels(d$Line), nlevels(d$Batch)),
      if (is.null(f)) "can't be tested." else tags$code(f),
      if (identical(input$design, "crossed") && !length(shared_batches(d)))
        " No batch holds several lines, so crossed isn't possible.",
      if (nlevels(d$Line) == 1) " Results apply to this line only.")
  })

  output$param_warning <- renderUI({
    req(input$params)
    cv <- coverage(prep(), input$params)
    bad <- cv$Parameter[!cv$Testable]
    if (length(bad)) div(class = "alert alert-warning small py-1 px-2",
      paste0("One line and one batch only, can't be tested: ", paste(bad, collapse = ", "),
             ". Unselect them."))
  })

  output$design_tbl <- renderTable(design_table(prep()$datos))
  output$coverage_tbl <- renderTable({
    p <- prep()
    cv <- coverage(p, p$MM_Vars[p$numeric])
    cv$Testable <- ifelse(cv$Testable, "yes", "no")
    cv
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
    nav_select("tabs", "ANOVA results")
    labels <- vapply(results()$res, `[[`, "", "label")
    updateSelectInput(session, "plot_param", choices = labels)
    updateSelectizeInput(session, "multi", choices = labels, selected = labels)
  })

  # Group order: reset to reference first on new data or reference
  tx_ord <- reactiveVal(NULL)
  observeEvent(list(prep(), input$ref), tx_ord(NULL))
  observeEvent(input$tx_order, tx_ord(unlist(input$tx_order)))

  picker <- function(key, val, type = "cols") {
    tags$input(type = "color", class = "mmc-pick", `data-input` = type,
               `data-key` = key, value = val, title = "Pick a color")
  }

  output$group_ui <- renderUI({
    d <- prep()$datos
    o <- list(order = tx_ord(), ref = input$ref, cols = isolate(input$cols))
    gs <- dot_style(d, o, "group")
    arrow <- function(dir, g) tags$button(type = "button", class = "btn btn-sm btn-link p-0 mmc-move",
      `data-dir` = dir, `aria-label` = paste("Move", g, dir), if (dir == "up") "▲" else "▼")
    tags$ul(class = "list-group mmc-order", lapply(tx_order(d, o), function(g) tags$li(
      class = "list-group-item d-flex align-items-center gap-2 py-1",
      `data-lev` = g,
      span("⠿", class = "text-muted"), picker(paste0("group:", g), gs$val[gs$lev == g]),
      span(g, class = "flex-grow-1"), arrow("up", g), arrow("down", g))))
  })

  # Color and/or shape pickers for each Line or Batch level
  output$dot_style_ui <- renderUI({
    d <- prep()$datos
    o <- list(cols = isolate(input$cols), shapes = isolate(input$shapes))
    by <- unique(setdiff(c(input$color_by, input$shape_by), c("none", "group")))
    lapply(by, function(b) {
      cs <- if (input$color_by == b) dot_style(d, o, b)
      ss <- if (input$shape_by == b) dot_style(d, o, b, "pch")
      s <- if (is.null(cs)) ss else cs
      tags$table(class = "table table-sm align-middle mb-2",
        tags$caption(class = "caption-top", s$title),
        lapply(seq_along(s$lev), function(i) {
          key <- paste0(b, ":", s$lev[i])
          tags$tr(tags$td(s$lev[i]),
            if (!is.null(cs)) tags$td(picker(key, cs$val[i])),
            if (!is.null(ss)) tags$td(tags$select(
              class = "mmc-pick form-select form-select-sm", `data-input` = "shapes",
              `data-key` = key, `aria-label` = paste("Shape for", s$lev[i]),
              lapply(seq_along(shape_set), function(j) tags$option(
                value = shape_set[[j]], selected = if (shape_set[[j]] == ss$val[i]) NA,
                names(shape_set)[j])))))
        }))
    })
  })

  # Tables follow the reference group and group order; pairs as "other - Ref"
  ref_order <- reactive({
    lev <- levels(results()$prep$datos$Tx)
    ref <- if (isTRUE(input$ref %in% lev)) input$ref else lev[1]
    list(lev = lev, ref = ref,
         ord = tx_order(results()$prep$datos, list(order = tx_ord(), ref = ref)))
  })

  anova_df <- reactive({
    df <- do.call(rbind, lapply(results()$res, `[[`, "table"))
    ro <- ref_order()
    df <- cbind(df[1], Comparison = paste(ro$ord, collapse = " vs "), df[-1])
    # Benjamini-Hochberg across parameters in this run
    i <- match("Pr(>F)", names(df))
    df <- cbind(df[seq_len(i)], q_FDR = p.adjust(df[[i]], "BH"), df[-seq_len(i)])
    # dCt data: each group's fold change vs the reference (model means)
    if (identical(input$scale, "fc")) {
      for (g in setdiff(ro$ord, ro$ref)) {
        df[[paste("FC", g, "vs", ro$ref)]] <- vapply(results()$res, function(r) {
          m <- setNames(r$means$emmean, r$means$Tx)
          2^-(m[[g]] - m[[ro$ref]])
        }, 0)
      }
    }
    # Plain column names
    lab <- c("Sum Sq" = "Sum of squares", "Mean Sq" = "Mean square", NumDF = "df (Tx)",
             DenDF = "df (error)", "F value" = "F", "Pr(>F)" = "p", q_FDR = "q (FDR)")
    i <- names(df) %in% names(lab)
    names(df)[i] <- lab[names(df)[i]]
    df
  })

  pairs_df <- reactive({
    df <- do.call(rbind, lapply(results()$res, `[[`, "pairs_table"))
    ro <- ref_order()
    prs <- combn(ro$lev, 2)  # emmeans pair order
    i <- rep(seq_len(ncol(prs)), length.out = nrow(df))
    flip <- prs[1, i] == ro$ref
    a <- ifelse(flip, prs[2, i], prs[1, i]); b <- ifelse(flip, prs[1, i], prs[2, i])
    df$contrast <- paste(a, "-", b)
    df$estimate <- ifelse(flip, -df$estimate, df$estimate)
    df$t.ratio <- ifelse(flip, -df$t.ratio, df$t.ratio)
    lo <- df$lower.CL
    df$lower.CL <- ifelse(flip, -df$upper.CL, lo)
    df$upper.CL <- ifelse(flip, -lo, df$upper.CL)
    # dCt data: fold change of the first group vs the second (2^-estimate)
    if (identical(input$scale, "fc")) {
      df$Fold_change <- 2^-df$estimate
      df$FC_lower <- 2^-df$upper.CL
      df$FC_upper <- 2^-df$lower.CL
    }
    df
  })

  # 3 sig. figs on screen; CSVs keep full precision
  show_p <- function(df) {
    for (col in intersect(c("p", "q (FDR)", "p.value"), names(df))) {
      df[[col]] <- as.character(signif(df[[col]], 3))
    }
    df
  }
  output$anova_table <- renderTable(show_p(anova_df()), digits = 4)

  # Residual checks, drawn only when shown
  observeEvent(results(), updateSelectInput(session, "resid_param",
    choices = vapply(results()$res, `[[`, "", "label")))
  output$resid_plot <- renderPlot({
    r <- results()
    x <- r$res[[match(input$resid_param, vapply(r$res, `[[`, "", "label"))]]
    req(x)
    rs <- residuals(x$MM_Form, type = "pearson", scaled = TRUE)
    col <- ifelse(abs(rs) > 3, "red", "grey40")
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), tcl = -0.25, mgp = c(2.5, 0.6, 0))
    # Q-Q with pointwise 95% band around the quartile line
    q <- qqnorm(rs, plot.it = FALSE)
    z <- sort(q$x); pp <- pnorm(z); n <- length(rs)
    b <- diff(quantile(rs, c(.25, .75), names = FALSE)) / diff(qnorm(c(.25, .75)))
    a <- quantile(rs, .25, names = FALSE) - b * qnorm(.25)
    se <- b * sqrt(pp * (1 - pp) / n) / dnorm(z)
    lo <- a + b * z - 1.96 * se; hi <- a + b * z + 1.96 * se
    plot(q, type = "n", las = 1, main = "Normal Q-Q", ylim = range(rs, lo, hi),
         xlab = "Theoretical quantiles", ylab = "Scaled residual")
    polygon(c(z, rev(z)), c(lo, rev(hi)), col = "grey90", border = NA)
    abline(a, b)
    points(q, pch = 19, col = col)
    # Residuals vs fitted; band = 95% expected range
    fv <- fitted(x$MM_Form)
    plot(fv, rs, type = "n", las = 1, main = "Residuals vs fitted", ylim = range(rs, -3.2, 3.2),
         xlab = "Fitted value", ylab = "Scaled residual")
    rect(par("usr")[1], -1.96, par("usr")[2], 1.96, col = "grey90", border = NA)
    abline(h = c(-3, 0, 3), lty = c(3, 1, 3), col = "grey60")
    points(fv, rs, pch = 19, col = col)
  }, width = 640, height = 320, res = 96)
  output$outlier_table <- renderTable({
    r <- results(); d <- r$prep$datos
    o <- do.call(rbind, lapply(r$res, function(x) if (nrow(x$out)) data.frame(
      Parameter = x$label, "CSV row" = as.integer(rownames(d)[x$out$i]) + 1L,
      Tx = d$Tx[x$out$i], Line = d$Line[x$out$i], Batch = d$Batch[x$out$i],
      Value = d[[x$var]][x$out$i], "Scaled residual" = x$out$resid, check.names = FALSE)))
    if (is.null(o)) data.frame(Result = "None flagged") else o
  }, digits = 3)
  output$pairs_table <- renderTable(show_p(pairs_df()), digits = 4)

  output$anova_print <- renderPrint({
    r <- results()
    labels <- vapply(r$res, `[[`, "", "label")
    cat("Parameters: ", paste(labels, collapse = ", "), "\n", sep = "")
    for (x in r$res) {
      cat("\n", x$label, " ~ Tx + ", x$rand, "\n", sep = "")
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
                       color_by = input$color_by, shape_by = input$shape_by,
                       cols = input$cols, shapes = input$shapes, order = tx_ord(),
                       labels = input$sig, size = input$pt_size,
                       brackets = !identical(input$sig, "none"), scale = input$scale, ref = input$ref,
                       ylab = input$ylab, err = input$err,
                       center = input$center, rot = input$rot,
                       adj_batch = input$adj_batch, outliers = input$outliers))

  # Size in inches, clamped to 3-12
  dims <- reactive({
    fit <- function(x) if (is.numeric(x) && !is.na(x)) min(max(x, 3), 12) else 5
    w <- fit(input$w)
    c(w = w, h = if (isTRUE(input$square)) w else fit(input$h))
  })

  # Current figure: one parameter or all selected in one figure
  draw_current <- function(r) {
    if (identical(input$fig, "all")) {
      sel <- r$res[vapply(r$res, `[[`, "", "label") %in% input$multi]
      validate(need(length(sel) > 0, "Select at least one parameter."))
      draw_multi(r$prep$datos, sel, r$adjust, opt())
    } else draw_plot(r$prep$datos, current(), r$adjust, opt())
  }

  output$plot <- renderPlot(
    draw_current(results()),
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

  plot_name <- function(r, ext) {
    what <- if (identical(input$fig, "all")) "all_parameters" else make.names(input$plot_param)
    paste0(r$base, "_", what, ".", ext)
  }

  observeEvent(input$dl_plot, {
    r <- ready(); req(r)
    f <- plot_formats[[input$dl_plot]]
    d <- dims()
    name <- if (input$dl_plot == "pdf_all") paste0(r$base, "_MM_KR_plots.pdf")
            else plot_name(r, f$ext)
    save_file(name, f$type, function(file) {
      f$open(file, d[["w"]], d[["h"]])
      on.exit(dev.off())
      if (input$dl_plot == "pdf_all") for (x in r$res) draw_plot(r$prep$datos, x, r$adjust, opt())
      else draw_current(r)
    })
  })

  observeEvent(input$example, {
    save_file("example_data.csv", "text/csv",
              function(f) file.copy("example_data.csv", f))
  })
}

shinyApp(ui, server)
