#!/usr/bin/env Rscript

# enrich.R — functional over-representation analysis (ORA) with clusterProfiler.
# Gold-standard Bioconductor enrichment: hypergeometric test against GO / KEGG /
# Reactome, with a proper background universe (the genes actually tested in DE).
#
# Usage:
#   Rscript enrich.R \
#     --genes     genes.txt        (gene SYMBOLS, one per line — the DEGs) \
#     --universe  universe.txt     (all tested gene SYMBOLS; optional but recommended) \
#     --organism  human|mouse \
#     --databases GO_BP,GO_MF,GO_CC,KEGG,Reactome \
#     --outdir    results/enrichment \
#     [--padj 0.05] [--qval 0.2]
#
# Writes, per database: <db>.csv with columns
#   ID, Term, Pvalue, AdjPvalue, qvalue, Count, GeneRatio, BgRatio, Genes
# and enrich_index.json summarizing which databases returned terms.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(jsonlite)
})

parse_args <- function(a){o<-list();i<-1;while(i<=length(a)){o[[sub("^--","",a[i])]]<-a[i+1];i<-i+2};o}
A <- parse_args(commandArgs(trailingOnly = TRUE))
stopifnot(!is.null(A$genes), !is.null(A$organism), !is.null(A$outdir))
padj_cut <- as.numeric(if (is.null(A$padj)) 0.05 else A$padj)
qval_cut <- as.numeric(if (is.null(A$qval)) 0.2 else A$qval)
dbs <- strsplit(if (is.null(A$databases)) "GO_BP,KEGG,Reactome" else A$databases, ",")[[1]]
dir.create(A$outdir, showWarnings = FALSE, recursive = TRUE)

org_pkg  <- if (A$organism == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
kegg_org <- if (A$organism == "mouse") "mmu" else "hsa"
reactome_org <- if (A$organism == "mouse") "mouse" else "human"
suppressPackageStartupMessages(library(org_pkg, character.only = TRUE))
OrgDb <- get(org_pkg)

read_genes <- function(p) unique(trimws(readLines(p, warn = FALSE)))
genes <- read_genes(A$genes)
universe <- if (!is.null(A$universe) && file.exists(A$universe)) read_genes(A$universe) else NULL

# SYMBOL -> ENTREZID (needed for KEGG/Reactome; GO can use SYMBOL directly)
to_entrez <- function(syms) {
  m <- suppressWarnings(bitr(syms, "SYMBOL", "ENTREZID", OrgDb = OrgDb))
  m$ENTREZID[!is.na(m$ENTREZID)]
}
entrez_genes <- if (any(dbs %in% c("KEGG", "Reactome"))) to_entrez(genes) else NULL
entrez_univ  <- if (!is.null(universe) && any(dbs %in% c("KEGG", "Reactome"))) to_entrez(universe) else NULL

standardize <- function(er) {
  if (is.null(er) || nrow(as.data.frame(er)) == 0) return(NULL)
  df <- as.data.frame(er)
  data.frame(
    ID = df$ID, Term = df$Description, Pvalue = df$pvalue,
    AdjPvalue = df$p.adjust, qvalue = df$qvalue, Count = df$Count,
    GeneRatio = df$GeneRatio, BgRatio = df$BgRatio,
    Genes = df$geneID, stringsAsFactors = FALSE)
}

run_one <- function(db) {
  tryCatch({
    if (grepl("^GO_", db)) {
      ont <- sub("^GO_", "", db)  # BP / MF / CC
      er <- enrichGO(gene = genes, OrgDb = OrgDb, keyType = "SYMBOL", ont = ont,
                     universe = universe, pvalueCutoff = padj_cut,
                     qvalueCutoff = qval_cut, pAdjustMethod = "BH", readable = FALSE)
    } else if (db == "KEGG") {
      er <- enrichKEGG(gene = entrez_genes, organism = kegg_org,
                       universe = entrez_univ, pvalueCutoff = padj_cut,
                       qvalueCutoff = qval_cut, pAdjustMethod = "BH")
      if (!is.null(er)) er <- setReadable(er, OrgDb, keyType = "ENTREZID")
    } else if (db == "Reactome") {
      if (!requireNamespace("ReactomePA", quietly = TRUE)) stop("ReactomePA not installed")
      er <- ReactomePA::enrichPathway(gene = entrez_genes, organism = reactome_org,
                       universe = entrez_univ, pvalueCutoff = padj_cut,
                       qvalueCutoff = qval_cut, pAdjustMethod = "BH", readable = TRUE)
    } else {
      stop(paste("unknown database:", db))
    }
    standardize(er)
  }, error = function(e) { message(sprintf("%s failed: %s", db, conditionMessage(e))); NULL })
}

index <- list()
for (db in dbs) {
  res <- run_one(db)
  if (!is.null(res)) {
    write.csv(res, file.path(A$outdir, sprintf("%s.csv", db)), row.names = FALSE)
    index[[db]] <- nrow(res)
  } else {
    index[[db]] <- 0
  }
}

write_json(list(databases = index, n_genes = length(genes),
                n_universe = if (is.null(universe)) NA else length(universe),
                organism = A$organism, padj = padj_cut),
           file.path(A$outdir, "enrich_index.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Enrichment done:\n"); print(unlist(index))
