#!/usr/bin/env bash
# One-time environment setup for labmates running this app natively (no Docker).
# Installs the R/Bioconductor packages the engine scripts need, then the Python deps.
#
# Usage: ./setup.sh
set -e
cd "$(dirname "$0")"

echo "== Checking R =="
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found on PATH. Install R >= 4.4 first: https://cran.r-project.org/"
  exit 1
fi
Rscript -e 'if (getRversion() < "4.4") stop("R >= 4.4 required, found ", getRversion())'
Rscript --version

echo
echo "== Installing R packages (Bioconductor + CRAN) =="
echo "This installs DESeq2/edgeR/limma (RNA) and"
echo "clusterProfiler/org.Hs.eg.db/org.Mm.eg.db/ReactomePA (enrichment). First run"
echo "can take 15-30+ minutes — Bioconductor source/binary downloads are large."
Rscript - <<'RSCRIPT'
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "DESeq2", "edgeR", "limma", "GEOquery",              # RNA-seq engine
  "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "ReactomePA"  # enrichment
)
cran_pkgs <- c("dplyr", "ggplot2", "jsonlite", "pheatmap", "readr", "reshape2")

missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")

still_missing <- c(bioc_pkgs, cran_pkgs)[!sapply(c(bioc_pkgs, cran_pkgs), requireNamespace, quietly = TRUE)]
if (length(still_missing)) {
  cat("\nFAILED to install:", paste(still_missing, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("\nAll R packages installed OK.\n")
}
RSCRIPT

echo
echo "== Checking Python =="
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 not found on PATH. Install Python >= 3.10 first."
  exit 1
fi
python3 --version

echo
echo "== Installing Python packages =="
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt

echo
echo "== Done =="
echo "Run the app with: ./run.sh"
echo "(Optional) The AI literature-synthesis feature needs the 'claude' CLI logged"
echo "in locally (claude.ai/code) — the app works fully without it, that panel just"
echo "won't be available."
