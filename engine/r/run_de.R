#!/usr/bin/env Rscript

# run_de.R — unified bulk RNA-seq differential expression engine
# Runs DESeq2, edgeR, and limma-voom on the same counts + design and emits
# standardized result tables, a 3-method comparison, batch-effect diagnostics,
# and plots (volcano per method + heatmap of top genes).
#
# Usage:
#   Rscript run_de.R \
#     --counts   path/to/counts.csv      (genes x samples; first col = gene id) \
#     --samples  path/to/samples.csv      (sample, condition, [batch]) \
#     --contrast "treatment,control"      (test vs reference level of `condition`) \
#     --outdir   path/to/results \
#     [--batch   batch]                   (column in samples.csv to model/correct) \
#     [--padj 0.05] [--lfc 1] [--topn 50]
#
# Outputs in --outdir:
#   results_<method>.csv        gene, log2FC, pvalue, padj, method, direction
#   results_combined.csv        all methods stacked
#   comparison.csv              per-gene significance across methods
#   comparison_summary.json     counts of up/down/overlap per method
#   volcano_<method>.png        volcano plot per method
#   heatmap_top.png             z-scored heatmap of top genes
#   pca.png                     PCA colored by condition (and batch if given)
#   batch_report.json           PC-vs-metadata association (batch-effect signal)
#   run_meta.json               parameters used

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# Versions of R and every R package the portal's analysis uses (DE engine, GEO
# fetch, enrichment, ID mapping), recorded with each run for methods sections.
software_versions <- function() {
  pkgs <- c("DESeq2", "edgeR", "limma", "GEOquery", "clusterProfiler", "ReactomePA",
            "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "ggplot2", "pheatmap", "jsonlite")
  v <- vapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)),
                                         error = function(e) "not installed"), "")
  bioc <- tryCatch(as.character(BiocManager::version()), error = function(e) "unknown")
  c(list(R = paste(R.version$major, R.version$minor, sep = "."), Bioconductor = bioc), as.list(v))
}

# ---- tiny arg parser -------------------------------------------------------
parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    val <- args[i + 1]
    out[[key]] <- val
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot(!is.null(args$counts), !is.null(args$samples),
          !is.null(args$contrast), !is.null(args$outdir))

sig_cut    <- as.numeric(args$padj %||% 0.05)   # threshold value (applied to sig_metric)
lfc_cut    <- as.numeric(args$lfc  %||% 1)
topn       <- as.integer(args$topn %||% 50)
batch_col  <- args$batch                         # may be NULL
sig_metric <- args$sig_metric %||% "padj"        # "padj" or "pvalue" (lab canonical uses pvalue)
filter_mode <- args$filter %||% "canonical"      # "canonical" (TMM log2cpm) or "count"
groups_mode <- args$groups_mode %||% "subset"    # "subset" (only test+ref) or "all" groups
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- load data -------------------------------------------------------------
read_table_auto <- function(path) {
  sep <- if (grepl("\\.tsv$|\\.txt$", path)) "\t" else ","
  # read WITHOUT row.names so duplicate gene symbols don't error; make.unique them
  df <- read.delim(path, sep = sep, header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
  genes <- make.unique(as.character(df[[1]]))
  mat <- as.matrix(sapply(df[, -1, drop = FALSE], as.numeric))
  rownames(mat) <- genes
  mat
}

counts <- read_table_auto(args$counts)
counts[is.na(counts)] <- 0
storage.mode(counts) <- "integer"

samples <- read.delim(args$samples, sep = if (grepl("\\.tsv$|\\.txt$", args$samples)) "\t" else ",",
                      header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# normalize: require a 'sample' and 'condition' column
stopifnot("sample" %in% names(samples), "condition" %in% names(samples))
rownames(samples) <- samples$sample

# align counts columns to samples
common <- intersect(colnames(counts), samples$sample)
if (length(common) < 2) stop("Fewer than 2 samples match between counts columns and sample sheet.")
counts  <- counts[, common, drop = FALSE]
samples <- samples[common, , drop = FALSE]
samples$condition <- factor(samples$condition)

contrast <- strsplit(args$contrast, ",")[[1]]
test_lvl <- trimws(contrast[1]); ref_lvl <- trimws(contrast[2])
samples$condition <- relevel(samples$condition, ref = ref_lvl)
if (!all(c(test_lvl, ref_lvl) %in% levels(samples$condition)))
  stop(sprintf("Contrast levels '%s' / '%s' not found in condition column.", test_lvl, ref_lvl))
# For a pairwise comparison, optionally drop all other groups so normalization,
# the voom/variance trend, the gene filter and the model use ONLY the two groups
# being compared (matches the typical standalone two-group analysis).
if (groups_mode == "subset") {
  keep_s  <- as.character(samples$condition) %in% c(test_lvl, ref_lvl)
  samples <- samples[keep_s, , drop = FALSE]
  counts  <- counts[, rownames(samples), drop = FALSE]
  samples$condition <- relevel(droplevels(samples$condition), ref = ref_lvl)
  message(sprintf("groups=subset: modeling only %s vs %s (%d samples).",
                  test_lvl, ref_lvl, ncol(counts)))
}

# coefficient that edgeR/limma must test: the test level (ref is the baseline).
# With >2 groups this is NOT simply the first condition column, so name it explicitly.
test_coef <- paste0("condition", test_lvl)

has_batch <- !is.null(batch_col) && batch_col %in% names(samples)
if (has_batch) samples[[batch_col]] <- factor(samples[[batch_col]])

# design formula
design_formula <- if (has_batch)
  as.formula(sprintf("~ %s + condition", batch_col)) else ~ condition

# filter genes ONCE so all three methods see the same gene set (fair comparison).
min_n <- min(table(samples$condition))
if (filter_mode == "canonical") {
  # lab canonical: TMM normalize, keep genes with log2CPM > 1 in >= min_group_n samples
  dge0 <- calcNormFactors(DGEList(counts = counts), method = "TMM")
  log2cpm0 <- cpm(dge0, log = TRUE, prior.count = 1)
  keep <- rowSums(log2cpm0 > 1) >= min_n
} else {
  keep <- rowSums(counts >= 10) >= max(2, floor(min_n / 2))
}
counts <- counts[keep, , drop = FALSE]
message(sprintf("Filter (%s): kept %d genes across %d samples.",
                filter_mode, nrow(counts), ncol(counts)))

standardize <- function(df, method) {
  df$method <- method
  sig_val <- if (sig_metric == "pvalue") df$pvalue else df$padj
  df$direction <- ifelse(sig_val < sig_cut & df$log2FC >= lfc_cut, "up",
                  ifelse(sig_val < sig_cut & df$log2FC <= -lfc_cut, "down", "ns"))
  df$direction[is.na(df$direction)] <- "ns"
  df[, c("gene", "log2FC", "pvalue", "padj", "method", "direction")]
}

# ---- DESeq2 ----------------------------------------------------------------
# DESeq2 prints a note when condition levels contain spaces/punctuation (e.g. GEO's
# "Alzheimer's disease"). Harmless, so hide only that note; other messages pass.
hide_level_note <- function(expr) withCallingHandlers(expr, message = function(m) {
  if (grepl("levels of factors in the design", conditionMessage(m), fixed = TRUE))
    invokeRestart("muffleMessage")
})

run_deseq2 <- function() {
  dds <- DESeqDataSetFromMatrix(counts, colData = samples, design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", test_lvl, ref_lvl))
  data.frame(gene = rownames(res), log2FC = res$log2FoldChange,
             pvalue = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
}

# ---- edgeR (glmQLF) --------------------------------------------------------
run_edger <- function() {
  d <- DGEList(counts = counts, group = samples$condition)
  d <- calcNormFactors(d, method = "TMM")
  mm <- model.matrix(design_formula, data = samples)
  d <- estimateDisp(d, mm)
  fit <- glmQLFit(d, mm)
  stopifnot(test_coef %in% colnames(mm))   # test the user's test-vs-ref, not the 1st group
  qlf <- glmQLFTest(fit, coef = test_coef)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tt), log2FC = tt$logFC,
             pvalue = tt$PValue, padj = tt$FDR, stringsAsFactors = FALSE)
}

# ---- limma-voom ------------------------------------------------------------
run_limma <- function() {
  d <- DGEList(counts = counts)
  d <- calcNormFactors(d, method = "TMM")
  mm <- model.matrix(design_formula, data = samples)
  v <- voom(d, mm)
  fit <- lmFit(v, mm)
  fit <- eBayes(fit)
  stopifnot(test_coef %in% colnames(mm))
  tt <- topTable(fit, coef = test_coef, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), log2FC = tt$logFC,
             pvalue = tt$P.Value, padj = tt$adj.P.Val, stringsAsFactors = FALSE)
}

message("Running DESeq2 ..."); r_deseq <- standardize(hide_level_note(run_deseq2()), "DESeq2")
message("Running edgeR ...");  r_edger <- standardize(run_edger(),  "edgeR")
message("Running limma-voom ..."); r_limma <- standardize(run_limma(), "limma-voom")

methods <- list(DESeq2 = r_deseq, edgeR = r_edger, `limma-voom` = r_limma)
for (m in names(methods)) {
  write.csv(methods[[m]], file.path(args$outdir, sprintf("results_%s.csv", m)), row.names = FALSE)
}
combined <- do.call(rbind, methods)
write.csv(combined, file.path(args$outdir, "results_combined.csv"), row.names = FALSE)

# ---- 3-method comparison ---------------------------------------------------
genes <- sort(unique(combined$gene))
comp <- data.frame(gene = genes, stringsAsFactors = FALSE)
for (m in names(methods)) {
  d <- methods[[m]]; rownames(d) <- d$gene
  comp[[paste0(m, "_log2FC")]] <- d[genes, "log2FC"]
  comp[[paste0(m, "_padj")]]   <- d[genes, "padj"]
  comp[[paste0(m, "_dir")]]    <- d[genes, "direction"]
}
sig_mat <- sapply(names(methods), function(m) comp[[paste0(m, "_dir")]] %in% c("up", "down"))
comp$n_methods_sig <- rowSums(sig_mat)
write.csv(comp, file.path(args$outdir, "comparison.csv"), row.names = FALSE)

summary_list <- lapply(names(methods), function(m) {
  d <- methods[[m]]
  list(up = sum(d$direction == "up"), down = sum(d$direction == "down"),
       total_sig = sum(d$direction %in% c("up", "down")))
})
names(summary_list) <- names(methods)
overlap_all <- sum(comp$n_methods_sig == length(methods))
write_json(list(methods = summary_list,
                significant_in_all = overlap_all,
                significant_in_any = sum(comp$n_methods_sig >= 1),
                sig_metric = sig_metric, sig_cut = sig_cut, lfc_cut = lfc_cut,
                filter_mode = filter_mode, n_genes_tested = nrow(counts)),
           file.path(args$outdir, "comparison_summary.json"), auto_unbox = TRUE, pretty = TRUE)

# ---- volcano per method ----------------------------------------------------
volcano <- function(df, method) {
  metric_vals <- if (sig_metric == "pvalue") df$pvalue else df$padj
  df <- df[is.finite(metric_vals) & !is.na(metric_vals), ]
  metric_vals <- if (sig_metric == "pvalue") df$pvalue else df$padj
  df$neglogp <- -log10(metric_vals + 1e-300)
  ylab <- if (sig_metric == "pvalue") "-log10(p-value)" else "-log10(adjusted p)"
  p <- ggplot(df, aes(log2FC, neglogp, color = direction)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c(up = "#d62728", down = "#1f77b4", ns = "grey70")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(sig_cut), linetype = "dashed", color = "grey40") +
    labs(title = paste("Volcano:", method), x = "log2 fold change", y = ylab) +
    theme_bw()
  ggsave(file.path(args$outdir, sprintf("volcano_%s.png", method)), p,
         width = 6, height = 5, dpi = 120)
}
for (m in names(methods)) volcano(methods[[m]], m)

# ---- normalized expression (for heatmap + batch PCA + gene profiles) -------
dds_norm <- hide_level_note(DESeqDataSetFromMatrix(counts, colData = samples, design = ~ 1))
vsd <- tryCatch(hide_level_note(assay(vst(dds_norm, blind = TRUE))),
                error = function(e) log2(edgeR::cpm(counts) + 1))
write.csv(data.frame(gene = rownames(vsd), vsd, check.names = FALSE),
          file.path(args$outdir, "normalized_expression.csv"), row.names = FALSE)

# per-method normalized expression (each method's native space) for gene profiles
dge_tmm <- calcNormFactors(DGEList(counts = counts), method = "TMM")
log2cpm_all <- cpm(dge_tmm, log = TRUE, prior.count = 1)
v_all <- voom(dge_tmm, model.matrix(design_formula, data = samples))
norm_by_method <- list(DESeq2 = vsd, edgeR = log2cpm_all, `limma-voom` = v_all$E)
for (m in names(norm_by_method)) {
  mat <- norm_by_method[[m]]
  write.csv(data.frame(gene = rownames(mat), mat, check.names = FALSE),
            file.path(args$outdir, sprintf("norm_%s.csv", m)), row.names = FALSE)
}

# ---- heatmap of top genes (union of top by DESeq2 padj) --------------------
top_genes <- head(r_deseq$gene[order(r_deseq$padj)], topn)
top_genes <- top_genes[top_genes %in% rownames(vsd)]
if (length(top_genes) >= 2) {
  mat <- vsd[top_genes, , drop = FALSE]
  mat <- t(scale(t(mat)))
  ann <- samples[, c("condition", if (has_batch) batch_col), drop = FALSE]
  png(file.path(args$outdir, "heatmap_top.png"), width = 900, height = 1100, res = 120)
  pheatmap(mat, annotation_col = ann, show_rownames = length(top_genes) <= 60,
           main = sprintf("Top %d genes (z-scored)", length(top_genes)))
  dev.off()
}

# ---- batch-effect diagnostics: PCA + PC~metadata association ---------------
pca <- prcomp(t(vsd), scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pcdf <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   condition = samples$condition)
if (has_batch) pcdf$batch <- samples[[batch_col]]
p_pca <- ggplot(pcdf, aes(PC1, PC2, color = condition,
                          shape = if (has_batch) batch else NULL)) +
  geom_point(size = 3) +
  labs(title = "PCA (VST)",
       x = sprintf("PC1 (%.1f%%)", 100 * var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * var_exp[2])) + theme_bw()
ggsave(file.path(args$outdir, "pca.png"), p_pca, width = 6, height = 5, dpi = 120)

# association of top PCs with metadata columns via ANOVA R^2
assoc <- list()
meta_cols <- c("condition", if (has_batch) batch_col)
for (pc in 1:min(3, ncol(pca$x))) {
  for (mc in meta_cols) {
    fit <- summary(lm(pca$x[, pc] ~ samples[[mc]]))
    assoc[[length(assoc) + 1]] <- list(PC = paste0("PC", pc), variable = mc,
                                       var_explained = round(100 * var_exp[pc], 2),
                                       r2 = round(fit$r.squared, 3),
                                       p = signif(stats::pf(fit$fstatistic[1],
                                            fit$fstatistic[2], fit$fstatistic[3],
                                            lower.tail = FALSE), 3))
  }
}
batch_flag <- FALSE
if (has_batch) {
  # flag if batch explains >=25% of a top-2 PC
  for (a in assoc) if (a$variable == batch_col && a$PC %in% c("PC1","PC2") && a$r2 >= 0.25) batch_flag <- TRUE
}
write_json(list(has_batch_column = has_batch,
                batch_effect_suspected = batch_flag,
                associations = assoc),
           file.path(args$outdir, "batch_report.json"), auto_unbox = TRUE, pretty = TRUE)

write_json(list(counts = args$counts, samples = args$samples,
                contrast = list(test = test_lvl, reference = ref_lvl),
                batch = batch_col %||% NA, sig_metric = sig_metric,
                sig_cut = sig_cut, lfc = lfc_cut, filter_mode = filter_mode,
                n_genes_tested = nrow(counts), n_samples = ncol(counts),
                software = software_versions()),
           file.path(args$outdir, "run_meta.json"), auto_unbox = TRUE, pretty = TRUE)

message("Done. Results in: ", args$outdir)
