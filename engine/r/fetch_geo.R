#!/usr/bin/env Rscript

# fetch_geo.R — pull metadata (and attempt raw counts) for a GEO series.
#
# Usage:
#   Rscript fetch_geo.R --accession GSE12345 --outdir data/GSE12345
#
# Writes:
#   metadata.csv        one row per sample: geo_accession, title, + parsed
#                       characteristics_ch1 fields (candidate condition/batch cols)
#   supp_files.json     list of downloaded supplementary files (likely counts)
#   fetch_report.json   what was found + guidance flags
#
# NOTE: GEO stores RNA-seq counts in heterogeneous supplementary files. This
# script reliably retrieves *metadata*; it downloads supplementary files and
# reports them, but the user/app must confirm which file is the counts matrix.

suppressPackageStartupMessages({
  library(GEOquery)
  library(jsonlite)
})

parse_args <- function(a) { o <- list(); i <- 1
  while (i <= length(a)) { o[[sub("^--", "", a[i])]] <- a[i+1]; i <- i + 2 }; o }
args <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot(!is.null(args$accession), !is.null(args$outdir))
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

acc <- toupper(trimws(args$accession))
if (!grepl("^GSE[0-9]+$", acc))
  stop("Invalid GEO series accession: ", acc, " (expected GSE<digits>)")
# collect into an environment (reference semantics) so tryCatch handlers can
# update it from any scope without the pitfalls of <<- at top level
R <- new.env()
R$accession <- acc; R$metadata_ok <- FALSE; R$supp_listed <- FALSE
R$notes <- character(0)
add_note <- function(msg) R$notes <- c(R$notes, msg)

# ---- metadata --------------------------------------------------------------
tryCatch({
  gse <- getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)
  pdata <- Biobase::pData(gse[[1]])
  # keep identifiers + characteristics fields (these hold condition/batch info)
  keep_cols <- grep("geo_accession|title|source_name|characteristics|:ch1$",
                    names(pdata), value = TRUE, ignore.case = TRUE)
  keep_cols <- unique(c("geo_accession", "title", keep_cols))
  keep_cols <- keep_cols[keep_cols %in% names(pdata)]
  md <- pdata[, keep_cols, drop = FALSE]
  names(md) <- sub(":ch1$", "", names(md))   # tidy ":ch1" suffixes
  write.csv(md, file.path(args$outdir, "metadata.csv"), row.names = FALSE)
  R$metadata_ok <- TRUE
  R$n_samples <- nrow(md)
  R$candidate_columns <- setdiff(names(md), c("geo_accession", "title"))
}, error = function(e) add_note(paste("metadata:", conditionMessage(e))))

# ---- supplementary files: LIST only, do NOT download (can be huge) ----------
tryCatch({
  sf <- getGEOSuppFiles(acc, makeDirectory = FALSE, fetch_files = FALSE)
  files <- if (!is.null(sf) && "fname" %in% colnames(sf)) as.character(sf$fname)
           else if (!is.null(sf) && nrow(sf) > 0) basename(rownames(sf))
           else character(0)
  urls  <- if (!is.null(sf) && "url" %in% colnames(sf)) as.character(sf$url)
           else if (!is.null(sf) && nrow(sf) > 0) rownames(sf)
           else character(0)
  looks_counts <- grepl("count|raw|matrix|expr|fpkm|tpm|rsem|featurecounts",
                        files, ignore.case = TRUE)
  R$supp_listed <- length(files) > 0
  R$supp_files <- as.list(files)
  R$likely_counts_files <- as.list(files[looks_counts])
  R$likely_counts_urls <- as.list(urls[looks_counts])   # direct GEO download URLs
  if (!any(looks_counts) && length(files) > 0)
    add_note("No author supplementary file name obviously matches 'counts'.")
}, error = function(e) add_note(paste("supp:", conditionMessage(e))))

if (!R$metadata_ok)
  add_note("Metadata fetch failed; fall back to manual upload.")

report <- as.list(R)
write_json(report, file.path(args$outdir, "fetch_report.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat(toJSON(report, auto_unbox = TRUE, pretty = TRUE))
