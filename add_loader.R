# Adds the loading screen to the built site (run by the deploy workflow).
site <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "site"
w <- file.path(site, "shinylive/webr/webr-worker.js")

# Install the download reporter only if webR's worker is a classic script
ok <- file.exists(w) && !any(grepl("^\\s*(import|export) ", readLines(w, warn = FALSE)))
if (ok) invisible(file.copy("mmc-worker.js", file.path(site, "shinylive/webr/mmc-worker.js"), overwrite = TRUE))
pkgs <- length(list.files(file.path(site, "shinylive/webr/packages"), "\\.tgz$", recursive = TRUE))

loader <- readLines("loader.html", encoding = "UTF-8")
loader <- sub("__MMC_SHIM__", tolower(ok), loader, fixed = TRUE)
loader <- sub("__MMC_PKGS__", max(pkgs, 1), loader, fixed = TRUE)
f <- file.path(site, "index.html")
x <- readLines(f, encoding = "UTF-8")
i <- grep("</body>", x, fixed = TRUE)
writeLines(append(x, loader, i - 1), f, useBytes = TRUE)
cat("Loading screen added; download tracking:", ok, "; packages:", pkgs, "\n")
