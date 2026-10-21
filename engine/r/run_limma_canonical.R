#!/usr/bin/env Rscript

# run_limma_canonical.R — faithful, runnable version of Dr. Zhang lab's
# DE_analysis.qmd (limma-voom, AD vs Control for GSE125583-style data).
# Same logic: TMM, filter rowSums(log2cpm > 1) >= min_group_n, design ~Condition,
# significance abs(logFC) > lfc & P.Value < pval (RAW p-value, ~1.5-fold default).
#
# Usage:
#   Rscript run_limma_canonical.R \
#     --counts   expression_matrix_SYMBOL.csv   (col 1 = gene symbol, rest = GSMs) \
#     --meta     metadata.csv                   (has geo_accession + diagnosis) \
#     --outdir   results/GSE125583 \
#     [--diagnosis_col diagnosis] [--accession_col geo_accession] \
#     [--lfc 0.58] [--pval 0.05]

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(edgeR); library(limma)
  library(ggplot2); library(reshape2); library(jsonlite)
})

parse_args <- function(a){o<-list();i<-1;while(i<=length(a)){o[[sub("^--","",a[i])]]<-a[i+1];i<-i+2};o}
A <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot(!is.null(A$counts), !is.null(A$meta), !is.null(A$outdir))
diag_col <- if (is.null(A$diagnosis_col)) "diagnosis" else A$diagnosis_col
acc_col  <- if (is.null(A$accession_col)) "geo_accession" else A$accession_col
lfc_cut  <- as.numeric(if (is.null(A$lfc)) 0.58 else A$lfc)
pval_cut <- as.numeric(if (is.null(A$pval)) 0.05 else A$pval)
dir.create(A$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load counts (first column = gene symbol) ------------------------------
counts_raw <- read.csv(A$counts, check.names = FALSE, stringsAsFactors = FALSE)
colnames(counts_raw)[1] <- "Gene"
counts_raw$Gene <- make.unique(as.character(counts_raw$Gene))
rownames(counts_raw) <- counts_raw$Gene
counts_mat <- as.matrix(sapply(counts_raw[, -1, drop = FALSE], as.numeric))
rownames(counts_mat) <- counts_raw$Gene

# ---- Load metadata, derive Condition (Control vs AD) -----------------------
meta <- read.csv(A$meta, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(acc_col %in% names(meta), diag_col %in% names(meta))
meta[[acc_col]]  <- as.character(meta[[acc_col]])
meta[[diag_col]] <- as.character(meta[[diag_col]])

sample_names <- colnames(counts_mat)
meta2 <- meta %>%
  filter(.data[[acc_col]] %in% sample_names) %>%
  mutate(Condition = case_when(
    tolower(.data[[diag_col]]) %in% c("control") ~ "Control",
    grepl("alzheimer", tolower(.data[[diag_col]])) ~ "AD",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(Condition))

meta2 <- meta2[match(sample_names, meta2[[acc_col]]), , drop = FALSE]
keep_samps <- !is.na(meta2$Condition)
metadata <- data.frame(
  Sample    = meta2[[acc_col]][keep_samps],
  Condition = factor(meta2$Condition[keep_samps], levels = c("Control", "AD")))
counts_mat <- counts_mat[, metadata$Sample, drop = FALSE]

cat("Samples per group:\n"); print(table(metadata$Condition))

# ---- DGEList, TMM, filter (canonical) --------------------------------------
dge <- DGEList(counts = counts_mat, samples = metadata, group = metadata$Condition)
dge <- calcNormFactors(dge, method = "TMM")
cat(sprintf("Genes before filtering: %d\n", nrow(dge)))
log2cpm_raw <- cpm(dge, log = TRUE, prior.count = 1)
min_n <- min(table(metadata$Condition))
keep  <- rowSums(log2cpm_raw > 1) >= min_n
dge_filt <- dge[keep, , keep.lib.sizes = FALSE]
dge_filt <- calcNormFactors(dge_filt, method = "TMM")
cat(sprintf("Genes after filtering: %d\n", nrow(dge_filt)))

# ---- voom + limma (canonical) ----------------------------------------------
design <- model.matrix(~ Condition, data = metadata)
rownames(design) <- colnames(dge_filt)
png(file.path(A$outdir, "voom_mean_var.png"), width = 900, height = 600)
v <- voom(dge_filt, design, plot = TRUE); dev.off()
fit <- eBayes(lmFit(v, design))
deg <- topTable(fit, coef = "ConditionAD", number = Inf, sort.by = "P")

deg_sig   <- subset(deg, abs(logFC) > lfc_cut & P.Value < pval_cut)
up_genes   <- rownames(deg_sig[deg_sig$logFC > 0, , drop = FALSE])
down_genes <- rownames(deg_sig[deg_sig$logFC < 0, , drop = FALSE])

# ---- save ------------------------------------------------------------------
write.csv(deg, file.path(A$outdir, "DEGs_AD_vs_Control_voom.csv"), row.names = TRUE)
write.csv(cpm(dge_filt, log = TRUE, prior.count = 1),
          file.path(A$outdir, "GSE125583_log2cpm_values.csv"))
writeLines(up_genes,   file.path(A$outdir, "Upregulated_AD_vs_Control.txt"))
writeLines(down_genes, file.path(A$outdir, "Downregulated_AD_vs_Control.txt"))
write_json(list(
  n_control = sum(metadata$Condition == "Control"),
  n_AD      = sum(metadata$Condition == "AD"),
  genes_before = nrow(dge), genes_after = nrow(dge_filt),
  lfc_cut = lfc_cut, pval_cut = pval_cut,
  n_sig = nrow(deg_sig), n_up = length(up_genes), n_down = length(down_genes)),
  file.path(A$outdir, "summary.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("\nDEGs (|logFC|>%.2f & P<%.2f): %d  (up %d / down %d)\n",
            lfc_cut, pval_cut, nrow(deg_sig), length(up_genes), length(down_genes)))
cat("Top 10 DEGs:\n"); print(head(deg[, c("logFC","P.Value","adj.P.Val")], 10))
cat("Results in:", A$outdir, "\n")
