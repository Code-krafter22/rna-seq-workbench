"""Check for and fetch NCBI-generated standardized RNA-seq raw counts for a GSE.

NCBI runs a uniform RNA-seq pipeline for many human/mouse GEO series and exposes
a clean genes x samples raw-counts matrix (rows = Entrez GeneID). This module:

  * check_counts(gse)   -> is that standardized counts matrix available? (no download)
  * download_counts(...) -> download it, map Entrez GeneID -> Symbol, write a CSV
                            compatible with the DE engine (first col = gene symbol)

Detection rule (validated): the counts download URL returns
content-type=application/octet-stream (a gzip) when available, or text/html when
the series has no standardized counts.
"""
from __future__ import annotations

import gzip
import io
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests

BASE = "https://www.ncbi.nlm.nih.gov/geo/download/"
# (organism key, genome build, Bioconductor annotation package for Entrez -> symbol)
BUILDS = [("human", "GRCh38.p13", "org.Hs.eg.db"), ("mouse", "GRCm39", "org.Mm.eg.db")]
ORG_DB = {org: db for org, _, db in BUILDS}

# GEO series accessions are always GSE + digits; reject anything else so a
# malicious "accession" can never be turned into a filesystem path segment
# (e.g. "../../etc") by callers building `DATA / gse`.
ACCESSION_RE = re.compile(r"^GSE\d+$")

# Hosts GEO is allowed to redirect/point us to for supplementary-file downloads.
ALLOWED_DOWNLOAD_HOSTS = {"www.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov"}


def is_valid_accession(gse: str) -> bool:
    return bool(ACCESSION_RE.match((gse or "").strip().upper()))


def _assert_safe_download_url(url: str) -> None:
    """Only allow https downloads from known NCBI/GEO hosts (defense against
    SSRF if a supplementary-file URL is ever manipulated upstream)."""
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in ALLOWED_DOWNLOAD_HOSTS:
        raise ValueError(f"Refusing to download from disallowed URL: {url}")


@dataclass
class CountsInfo:
    available: bool
    gse: str
    organism: str = ""
    build: str = ""
    url: str = ""
    filename: str = ""
    size_mb: float = 0.0


def _counts_url(gse: str, build: str) -> tuple[str, str]:
    fn = f"{gse}_raw_counts_{build}_NCBI.tsv.gz"
    return f"{BASE}?type=rnaseq_counts&acc={gse}&format=file&file={fn}", fn


def entrez_to_symbol(ids, organism: str) -> dict[str, str]:
    """Map Entrez GeneIDs -> HGNC/MGI symbols with the local org.*.eg.db package
    (already required for enrichment). NCBI's annot.tsv.gz download now sits
    behind a reCAPTCHA page, so it can't be fetched from a script."""
    r_code = (f"suppressMessages(library({ORG_DB[organism]})); "
              "ids <- readLines(file('stdin')); "
              f"s <- AnnotationDbi::mapIds({ORG_DB[organism]}, ids, 'SYMBOL', 'ENTREZID', "
              "multiVals = 'first'); "
              "writeLines(paste(ids, ifelse(is.na(s), '', s), sep = '\\t'))")
    p = subprocess.run(["Rscript", "-e", r_code], input="\n".join(map(str, ids)),
                       capture_output=True, text=True, timeout=300, check=True)
    pairs = (line.split("\t") for line in p.stdout.splitlines())
    return {k: v for k, v in pairs if v}


ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def ncbi_gene_names(ids: list[str], timeout: int = 60) -> dict[str, str]:
    """Symbols for Entrez IDs that org.*.eg.db dropped because NCBI retired or
    merged them (e.g. 100131096 = TNRC6C-AS1, now merged into TNRC6C). The
    E-utilities record keeps the symbol the ID had when the counts were made."""
    out = {}
    for i in range(0, len(ids), 400):  # ponytail: no API key, so ~3 req/s cap; fine for a few thousand ids
        r = requests.post(ESUMMARY, data={"db": "gene", "retmode": "json",
                                          "id": ",".join(ids[i:i + 400])}, timeout=timeout)
        r.raise_for_status()
        res = r.json()["result"]
        out.update({u: res[u]["name"] for u in res.get("uids", []) if res[u].get("name")})
    return out


def check_counts(gse: str, timeout: int = 30) -> CountsInfo:
    """Probe (human then mouse) for the standardized counts matrix without
    downloading it. Returns availability + metadata about the file."""
    gse = gse.strip().upper()
    if not is_valid_accession(gse):
        raise ValueError(f"Invalid GEO series accession: {gse!r} (expected GSE<digits>)")
    for organism, build, _ in BUILDS:
        url, fn = _counts_url(gse, build)
        try:
            r = requests.get(url, stream=True, timeout=timeout)
        except requests.RequestException:
            continue
        ctype = r.headers.get("content-type", "")
        cdisp = r.headers.get("content-disposition", "")
        clen = int(r.headers.get("content-length", 0) or 0)
        r.close()
        if "octet-stream" in ctype or "attachment" in cdisp:
            return CountsInfo(True, gse, organism, build, url, fn,
                              round(clen / 1e6, 1))
    return CountsInfo(False, gse)


def _read_gz_tsv(content: bytes, **kw) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(gzip.decompress(content)), sep="\t", **kw)


def download_counts(info: CountsInfo, outdir: str, timeout: int = 180,
                    map_symbols: bool = True) -> dict:
    """Download the standardized counts, map Entrez GeneID -> Symbol (best effort),
    and write `<gse>_counts_symbol.csv` (first column = gene). Returns a summary."""
    if not is_valid_accession(info.gse):
        raise ValueError(f"Invalid GEO series accession: {info.gse!r}")
    _assert_safe_download_url(info.url)
    r = requests.get(info.url, timeout=timeout)
    r.raise_for_status()
    counts = _read_gz_tsv(r.content)
    counts = counts.rename(columns={counts.columns[0]: "GeneID"})
    n_genes, n_samples = counts.shape[0], counts.shape[1] - 1

    symbols_mapped = False
    counts = counts.rename(columns={"GeneID": "Gene"})
    counts["Gene"] = counts["Gene"].astype(str)
    if map_symbols and info.organism in ORG_DB:
        try:
            sym_map = entrez_to_symbol(counts["Gene"], info.organism)
            symbols_mapped = bool(sym_map)
        except (OSError, subprocess.SubprocessError):
            sym_map = {}  # keep Entrez IDs; UI reports symbols_mapped=False
        missing = [g for g in counts["Gene"] if g not in sym_map]
        if missing:
            try:
                sym_map.update(ncbi_gene_names(missing))
            except (requests.RequestException, ValueError, KeyError):
                pass  # offline / NCBI hiccup: those rows keep their Entrez ID
        counts["Gene"] = counts["Gene"].map(sym_map).fillna(counts["Gene"])

    out = Path(outdir) / f"{info.gse}_counts_symbol.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out, index=False)
    return {"path": str(out), "n_genes": n_genes, "n_samples": n_samples,
            "symbols_mapped": symbols_mapped, "organism": info.organism,
            "build": info.build, "preview": counts.head(5)}


def download_supp(url: str, gse: str, outdir: str, timeout: int = 240) -> dict:
    """Download an author supplementary file (e.g. GSE..._counts.tsv.gz) directly
    from GEO, decompress if gzipped, parse it (delimiter sniffed), and write a
    counts CSV (first column = gene id). Returns shape + a small preview.

    Author count files are heterogeneous, so the caller should let the user
    eyeball the preview (first col = gene, remaining cols = samples)."""
    if not is_valid_accession(gse):
        raise ValueError(f"Invalid GEO series accession: {gse!r}")
    _assert_safe_download_url(url)
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    raw = r.content
    fn = url.rstrip("/").split("/")[-1]
    if fn.endswith(".gz"):
        raw = gzip.decompress(raw)
        fn = fn[:-3]
    # delimiter: trust extension, else let pandas sniff
    if fn.endswith((".tsv", ".tab", ".txt")):
        df = pd.read_csv(io.BytesIO(raw), sep="\t")
    elif fn.endswith(".csv"):
        df = pd.read_csv(io.BytesIO(raw), sep=",")
    else:
        df = pd.read_csv(io.BytesIO(raw), sep=None, engine="python")
    df = df.rename(columns={df.columns[0]: "Gene"})

    # are the values raw integer counts? (DESeq2/edgeR need counts, not TPM/FPKM)
    num = df.iloc[:, 1:].select_dtypes("number")
    looks_integer = bool(not num.empty
                         and (num.fillna(0).to_numpy() % 1 == 0).all())

    out = Path(outdir) / f"{gse}_supp_counts.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    return {"path": str(out), "filename": fn, "n_rows": int(df.shape[0]),
            "n_cols": int(df.shape[1] - 1), "columns": list(df.columns[:10]),
            "looks_integer": looks_integer, "preview": df.head(5)}


if __name__ == "__main__":
    import sys
    gse = sys.argv[1] if len(sys.argv) > 1 else "GSE125583"
    info = check_counts(gse)
    print(f"{gse}: available={info.available} organism={info.organism} "
          f"build={info.build} size={info.size_mb}MB")
    print("file:", info.filename)
