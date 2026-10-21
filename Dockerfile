# Bulk RNA-seq DE app — R/Bioconductor + Python/Streamlit in one image.
#
# Build:  docker build -t rna-seq-app .
# Run:    docker run -p 8501:8501 -v "$(pwd)/data:/app/data" -v "$(pwd)/results:/app/results" rna-seq-app
# Then open http://localhost:8501
#
# Optional: pass -e APP_PASSWORD=yourpassword to gate the app behind a shared
# password if you're exposing it beyond localhost (see app/streamlit_app.py).
#
# NOTE: the AI literature-synthesis feature shells out to a local `claude` CLI
# on the HOST machine's subscription login — it does not work inside this
# container (no login state, no network policy for it here). Everything else
# (RNA DE, enrichment, PubMed/Europe PMC/OpenAlex search)
# works fully in-container.

FROM bioconductor/bioconductor_docker:latest

WORKDIR /app

# ---- System deps: bioconductor_docker is R-focused, don't assume Python/curl
# are present ------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-venv python3-pip curl \
    && rm -rf /var/lib/apt/lists/*

# ---- R / Bioconductor packages (own layer: rarely changes, slow to build) --
RUN Rscript -e ' \
    bioc_pkgs <- c("DESeq2","edgeR","limma","GEOquery", \
                   "clusterProfiler","org.Hs.eg.db","org.Mm.eg.db","ReactomePA"); \
    cran_pkgs <- c("dplyr","ggplot2","jsonlite","pheatmap","readr","reshape2"); \
    BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE); \
    install.packages(cran_pkgs, repos = "https://cloud.r-project.org"); \
    missing <- c(bioc_pkgs, cran_pkgs)[!sapply(c(bioc_pkgs, cran_pkgs), requireNamespace, quietly = TRUE)]; \
    if (length(missing)) stop("FAILED to install: ", paste(missing, collapse = ", "))'

# ---- Python deps in a venv (own layer: changes more often than R packages) -
# A venv sidesteps PEP 668 "externally managed environment" pip restrictions
# on newer Debian/Ubuntu base images, regardless of what this base ships with.
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# ---- App code ----------------------------------------------------------------
COPY . /app

RUN mkdir -p /app/data /app/results

EXPOSE 8501
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

ENTRYPOINT ["streamlit", "run", "app/streamlit_app.py", \
           "--server.address=0.0.0.0", "--server.port=8501"]
