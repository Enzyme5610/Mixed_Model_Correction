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

# Bracket label: p-value or GraphPad Prism-style stars
sig_label <- function(p, opt) {
  if (opt$labels != "stars") return(format_p(p))
  as.character(cut(p, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                   c("****", "***", "**", "*", "ns"), right = FALSE))
}
star_key <- function(opt) {
  if (opt$labels == "stars" && opt$brackets) {
    mtext("ns p ≥ 0.05   * p < 0.05   ** p < 0.01   *** p < 0.001   **** p < 0.0001",
          side = 1, line = par("mar")[1] - 1, adj = 1, cex = 0.6 * par("cex"))
  }
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
    g <- legend(x, y, legend = k$s$lev, title = k$s$title, col = k$col, pch = k$pch,
                bty = "n", xpd = TRUE, cex = 0.85)
    y <- g$rect$top - g$rect$h
  }
  y
}

plot_defaults <- list(type = "dots", layout = "side", dots = TRUE, color_by = "line",
                      shape_by = "none", cols = list(), shapes = list(), order = NULL,
                      labels = "p",
                      size = 1, brackets = TRUE, scale = "raw", ref = NULL, ylab = "",
                      err = "ci", dir = "both", caps = TRUE, center = "diamond",
                      rot = "auto", adj_batch = FALSE)

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
  if (opt$caps) {
    suppressWarnings(arrows(xe, st$lo, xe, st$hi, angle = 90, length = 0.05, lwd = 2,
                            code = if (opt$dir == "up") 2 else 3))
  } else segments(xe, st$lo, xe, st$hi, lwd = 2)
  if (opt$type == "bar") return()
  if (opt$center == "line") segments(xe - half, st$m, xe + half, st$m, lwd = 3)
  else points(xe, st$m, pch = 23, bg = "white", cex = 1.4, lwd = 2)
}

# X labels at 0, 45 or 90 degrees; "auto" picks 90 when crowded
label_angle <- function(labels, opt, crowded = FALSE) {
  if (opt$rot == "auto") if (crowded) 90 else 0 else as.numeric(opt$rot)
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

err_label <- function(opt) {
  switch(opt$err, ci = "Model mean\n± 95% CI", se = "Model mean\n± SE",
         sem = "Mean ± SEM", sd = "Mean ± SD")
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
  op <- par(cex = min(1, min(dev.size("in")) / 5),
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
  xe <- xm + if (side) 0.15 else 0

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

  # Samples
  show_dots <- opt$type == "dots" || opt$dots
  if (show_dots) {
    col <- if (is.null(cs)) "grey45" else cs$val[cs$idx]
    points(xj, y, pch = if (is.null(ss)) 19 else ss$val[ss$idx], cex = opt$size,
           col = tint(col, 0.75))
  }

  draw_err(xe, st, opt)
  line_mark <- opt$type != "bar" && opt$center == "line"

  # Brackets; emmeans pair order matches combn on model levels
  prs <- combn(k, 2)
  for (i in seq_len(n_pairs)) {
    a <- min(pos[prs[, i]]); b <- max(pos[prs[, i]])
    yb <- rng[2] + h * (0.06 + 0.1 * (i - 1))
    tick <- h * 0.02
    segments(c(a, a, b), c(yb - tick, yb, yb), c(a, b, b), c(yb, yb, yb - tick))
    text((a + b) / 2, yb, sig_label(pw$p.value[i], opt), pos = 3, cex = 0.8, offset = 0.2)
  }
  star_key(opt)
  scale_note(opt, scale)

  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2])
  ky <- usr[3] + 0.3 * diff(usr[3:4])
  if (show_dots) ky <- min(ky, dot_legend(lx, usr[4], cs, ss))
  key <- data.frame(lab = err_label(opt), pch = if (opt$type == "bar") 22 else if (line_mark) NA else 23,
                    bg = if (opt$type == "bar") "grey88" else "white", col = "black",
                    lty = if (line_mark) 1 else NA)
  if (opt$type == "violin") key <- rbind(data.frame(lab = "Distribution", pch = 22, bg = "grey92", col = "grey45", lty = NA), key)
  if (show_dots) key <- rbind(data.frame(lab = "Sample", pch = 19, bg = NA, col = "grey40", lty = NA), key)
  legend(lx, ky, legend = key$lab, pch = key$pch, lty = key$lty,
         lwd = 3, col = key$col, pt.bg = key$bg, pt.lwd = 1, bty = "n", xpd = TRUE,
         cex = 0.8, y.intersp = 1.4)
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

  ang <- label_angle(labels, opt, crowded = n > 4 || max(nchar(labels)) > 10)
  op <- par(cex = min(1, min(dev.size("in")) / 5),
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
    if (show_dots) {
      col <- if (opt$color_by == "group") cols[gi] else if (is.null(cs)) "grey45" else cs$val[cs$idx]
      points(xs + runif(length(xs), -w * 0.25, w * 0.25), st$y,
             pch = if (is.null(ss)) 19 else ss$val[ss$idx],
             cex = opt$size, col = tint(col, 0.75))
    }
    draw_err(xm, st, opt, half = w * 0.3)

    # Brackets at one shared height; emmeans pair order matches combn on model levels
    pw <- results[[i]]$pairs
    prs <- combn(levels(datos$Tx), 2)
    for (p in seq_len(n_pairs)) {
      xa <- sort(i + off[match(prs[, p], ord)])
      yb <- rng[2] + h * (0.04 + 0.08 * (p - 1))
      segments(c(xa[1], xa[1], xa[2]), c(yb - h * 0.015, yb, yb),
               c(xa[1], xa[2], xa[2]), c(yb, yb, yb - h * 0.015))
      text(mean(xa), yb, sig_label(pw$p.value[p], opt), pos = 3, cex = 0.65, offset = 0.15)
    }
  }
  star_key(opt)
  scale_note(opt, sts[[1]]$scale)

  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2])
  g <- legend(lx, usr[4], legend = ord, title = "Group", pch = 22, pt.cex = 1.6,
              pt.bg = tint(cols, 0.5), col = cols, bty = "n", xpd = TRUE, cex = 0.85)
  ky <- usr[3] + 0.2 * diff(usr[3:4])
  if (show_dots) ky <- min(ky, dot_legend(lx, g$rect$top - g$rect$h, cs, ss))
  line_mark <- opt$type != "bar" && opt$center == "line"
  legend(lx, ky, legend = err_label(opt),
         pch = if (opt$type == "bar") NA else if (line_mark) NA else 23,
         lty = if (line_mark || opt$type == "bar") 1 else NA, lwd = 2, pt.bg = "white",
         bty = "n", xpd = TRUE, cex = 0.8)
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
              open = function(f, w, h) devEMF::emf(f, w, h, emfPlus = FALSE)),
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
    tags$style(".mmc-order li { cursor: grab; touch-action: none; } .mmc-pick[type=color] { width: 2em; height: 1.6em; padding: 0; border: 0; }"),
    fileInput("file", "1. Upload data (.csv)", accept = c(".csv", "text/csv")),
    helpText("Needs columns named Tx, Line and Batch (exact spelling), in any position."),
    actionLink("example", "Download an example file"),
    hr(),
    selectizeInput("params", "2. Parameters to analyze", choices = NULL,
                   multiple = TRUE, options = list(plugins = list("remove_button"))),
    helpText("Numeric columns are preselected."),
    hr(),
    radioButtons("design", "3. Do lines share batches?",
      choiceNames = list(
        tagList("No, each batch has one line", br(),
                tags$small(class = "text-muted font-monospace", "(1|Line/Batch)")),
        tagList("Yes, multiple lines were run per batch", br(),
                tags$small(class = "text-muted font-monospace", "(1|Line) + (1|Batch)"))),
      choiceValues = c("nested", "crossed")),
    helpText("Example of \"Yes\": L1-L4 all on the same qPCR plate labeled B1."),
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
    nav_panel("ANOVA results",
      helpText("Fold-change columns appear when the plot Y axis is set to fold change",
               "(ΔCt data); confidence intervals are in the Pairwise tab."),
      tableOutput("anova_table"),
      actionButton("dl_anova", "Download results (.csv)", icon = icon("download"))
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
          radioButtons("type", "Plot type", inline = TRUE,
                       choices = c("Dots" = "dots", "Bar" = "bar", "Violin" = "violin")),
          conditionalPanel("input.type == 'dots'",
            radioButtons("layout", "Mean and error bar", inline = TRUE,
                         choices = c("Beside dots" = "side", "Over dots" = "overlay"))),
          conditionalPanel("input.type != 'dots'",
            checkboxInput("dots", "Show dots", TRUE)),
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
          selectInput("err", "Error bars", choices = c(
            "95% CI (model)" = "ci", "SE (model)" = "se", "SEM" = "sem", "SD" = "sd")),
          conditionalPanel("input.err == 'sem' || input.err == 'sd'",
            helpText("SEM and SD use the raw values and ignore Line and Batch.")),
          radioButtons("dir", NULL, inline = TRUE,
                       choices = c("Both directions" = "both", "Above only" = "up")),
          checkboxInput("caps", "Caps", TRUE),
          conditionalPanel("input.type != 'bar'",
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
          radioButtons("sig", "Significance", inline = TRUE,
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
    nav_panel("Data preview", tableOutput("preview")),
    nav_panel("About",
      markdown("
**Model.** Each parameter: `parameter ~ Tx + (1 | Line/Batch)`, fit with
`lmerTest::lmer` (REML); Tx tested by Type II F test, Kenward-Roger df.
If lines share batches (step 3 \"Yes\"): `(1 | Line) + (1 | Batch)`.
One Line: `(1 | Batch)`. One Batch, or no replicates per line and batch:
`(1 | Line)`.

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
    df$Comparison <- paste(ro$ord, collapse = " vs ")
    # dCt data: each group's fold change vs the reference (model means)
    if (identical(input$scale, "fc")) {
      for (g in setdiff(ro$ord, ro$ref)) {
        df[[paste("FC", g, "vs", ro$ref)]] <- vapply(results()$res, function(r) {
          m <- setNames(r$means$emmean, r$means$Tx)
          2^-(m[[g]] - m[[ro$ref]])
        }, 0)
      }
    }
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
                       color_by = input$color_by, shape_by = input$shape_by,
                       cols = input$cols, shapes = input$shapes, order = tx_ord(),
                       labels = input$sig, size = input$pt_size,
                       brackets = !identical(input$sig, "none"), scale = input$scale, ref = input$ref,
                       ylab = input$ylab, err = input$err, dir = input$dir,
                       caps = input$caps, center = input$center, rot = input$rot,
                       adj_batch = input$adj_batch))

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
