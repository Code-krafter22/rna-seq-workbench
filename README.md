---
title: RNA-seq Workbench
emoji: 🧬
colorFrom: green
colorTo: gray
sdk: docker
app_port: 8501
pinned: false
short_description: Bulk RNA-seq DE with DESeq2, edgeR and limma-voom
---

# Bulk RNA-seq Differential Expression App (transcriptomics only)

A Streamlit web app for bulk RNA-seq differential expression. Takes a raw counts
matrix (from NCBI GEO or manual upload) and runs **three** DE pipelines side by side
(**DESeq2**, **edgeR**, **limma-voom**), checks for **batch effects**, and reports
up/down-regulated genes with **volcano plots**, a **heatmap**, a **cross-method
comparison**, and GO/KEGG/Reactome enrichment. You can then inspect a **gene of
interest** and run a **multi-database literature search** on it.

## Architecture

```
app/streamlit_app.py         Wizard UI + orchestration (Python / Streamlit)
engine/r/run_de.R            RNA DE engine: DESeq2 + edgeR + limma-voom, batch PCA, plots
engine/r/fetch_geo.R         GEO metadata + supplementary-counts fetch (GEOquery)
engine/r/enrich.R            GO/KEGG/Reactome ORA (clusterProfiler + ReactomePA)
engine/python/litreview.py   Multi-database literature search (PubMed + PubMed Central +
                             Europe PMC + OpenAlex)
engine/python/claude_llm.py  AI evidence synthesis via the Claude CLI (`claude -p`)
engine/python/pubmed.py      Standalone PubMed E-utilities helper
engine/python/geo.py         GEO counts-matrix detection/download
data/                        Inputs (uploaded / fetched)
results/                     Engine outputs (tables, plots, JSON)
```

## Wizard flow

1. **Input** — GEO accession (auto-fetch) or manual upload.
2. **Design** — condition column, contrast (test vs reference), optional batch column.
3. **Run** — DESeq2 + edgeR + limma-voom.
4. **Results** — per-method DE counts, batch-effect verdict + PCA, volcanoes, heatmap,
   downloadable comparison table, GO/KEGG/Reactome enrichment.
5. **Gene of interest** — cross-method stats + expression profile by group.
6. **Literature** — PubMed + PubMed Central + Europe PMC + OpenAlex search, plus optional AI
   synthesis via the Claude CLI.

## Setup

**Native (no Docker)** — for running locally or having labmates run their own copy:

```bash
./setup.sh   # installs R/Bioconductor packages + Python deps (first run: 15-30+ min)
./run.sh     # or: streamlit run app/streamlit_app.py
```

Then open http://localhost:8501. Needs **R ≥ 4.4** and **Python ≥ 3.10** already
installed; `setup.sh` handles everything else (`DESeq2`/`edgeR`/`limma`/`GEOquery` for
RNA, `clusterProfiler`/`org.Hs.eg.db`/
`org.Mm.eg.db`/`ReactomePA` for enrichment, plus `requirements.txt`).

**Docker** — avoids installing R/Bioconductor locally at all, and guarantees everyone
runs the identical pinned environment:

```bash
docker build -t rna-seq-app .
docker run -p 8501:8501 -v "$(pwd)/data:/app/data" -v "$(pwd)/results:/app/results" rna-seq-app
```

Then open http://localhost:8501. The volume mounts persist uploaded/fetched data and
results outside the container. Pass `-e APP_PASSWORD=yourpassword` if you're exposing
this beyond localhost (see `app/streamlit_app.py`'s `_require_auth`).

**Note**: the AI literature-synthesis panel shells out to a local `claude` CLI using
*your* logged-in Claude subscription — it only works where that CLI is authenticated
(your own machine), not inside Docker or a bare server. Everything else works fully
either way.

### Sharing with labmates

For a handful of labmates, running locally (native or Docker, per-person) is simpler
and cheaper than hosting a shared instance: this app needs R + Bioconductor alongside
Python, which free hosting tiers (Streamlit Community Cloud, Render/Railway/Fly.io
free tiers) generally can't run reliably — Bioconductor installs are slow/heavy and
free-tier RAM (512MB-1GB) is tight for DESeq2/limma on real datasets. If you do
want one shared instance later, a small paid VM running the Dockerfile above is the
straightforward path.

## Try it with synthetic data

`data/counts.csv` + `data/samples.csv` are a small synthetic RNA dataset (30 up, 30
down genes, mild batch). In the Input step choose **Manual upload**, upload those two
files, set condition = `condition`, contrast `treatment` vs `control`, batch = `batch`.

## Roadmap

- [ ] Wire the Python engine (PyDESeq2 + `inmoose` for edgeR/limma + pycombat).
- [ ] UpSet/Venn plot for method overlap.
- [ ] `ComBat_seq` batch *correction* option (currently batch is modeled in design).
- [ ] Async job queue + progress for large datasets.
- [ ] GSEA on the DE gene sets.
