"""Bulk RNA-seq (transcriptomics) differential expression web app (Streamlit MVP).

GEO/upload -> DESeq2 + edgeR + limma-voom -> results, gene of interest, literature.

The heavy lifting lives in engine/r/*.R (called via subprocess) and
engine/python/*.py. This file is orchestration + UI only.
"""
from __future__ import annotations

import hmac
import json
import os
import platform
import re
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
import requests
import streamlit as st

ROOT = Path(__file__).resolve().parent.parent
ENGINE_R = ROOT / "engine" / "r"
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DATA.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "engine" / "python"))
from geneplot import build_gene_figure, figure_bytes, PALETTES  # noqa: E402
from litreview import search_literature  # noqa: E402
from claude_llm import available as llm_available, synthesize_literature  # noqa: E402
from geo import check_counts, download_counts, download_supp, is_valid_accession  # noqa: E402
from enrichment import LIBRARIES as ENRICHR_LIBS, enrich as enrichr_enrich  # noqa: E402
from enrichment import top_terms_figure as enr_figure, figure_bytes as enr_bytes  # noqa: E402
sys.path.insert(0, str(ROOT / "app"))
import ui  # noqa: E402

METHODS = ["DESeq2", "edgeR", "limma-voom"]
STEPS = ["Input", "Design", "Run", "Results", "Gene of interest", "Literature"]

st.set_page_config(page_title="RNA-seq Workbench", layout="wide")
ui.inject_css()



def _require_auth() -> None:
    """Gate the whole app behind a shared password when APP_PASSWORD is set.

    This app has no per-user accounts and runs Rscript subprocesses + arbitrary
    GEO downloads on the host it's deployed on: fine for local, single-user use
    (the default, with APP_PASSWORD unset), but if it's ever bound to a network
    interface or tunneled, set APP_PASSWORD so it isn't wide open.
    """
    pw = os.environ.get("APP_PASSWORD")
    if not pw:
        return
    if st.session_state.get("_authed"):
        return
    ui.auth_header()
    entry = st.text_input("Password", type="password")
    if st.button("Log in", type="primary") or entry:
        if hmac.compare_digest(entry, pw):
            st.session_state["_authed"] = True
            st.rerun()
        elif entry:
            st.error("Incorrect password.")
    if not st.session_state.get("_authed"):
        st.stop()


_require_auth()
ss = st.session_state
ss.setdefault("mode", None)
ss.setdefault("lit_result", None)
ss.setdefault("lit_synthesis", None)
ss.setdefault("gene_of_interest", None)
ss.setdefault("counts_path", None)
ss.setdefault("samples_df", None)
ss.setdefault("run_dir", None)
ss.setdefault("geo", None)
ss.setdefault("geo_counts_loaded", None)
ss.setdefault("geo_err", None)
ss.setdefault("supp_loaded", None)
ss.setdefault("sample_col_auto", None)
ss.setdefault("samples_std", None)
ss.setdefault("contrast", None)
ss.setdefault("batch_col", None)
ss.setdefault("enrich", None)
ss.setdefault("enrich_meta", None)



def run_cmd(cmd: list[str], timeout: int = 1200) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def _guess(cols: list[str], wanted: list[str]) -> int:
    """Return index of the first column whose name matches a wanted keyword."""
    low = [c.lower() for c in cols]
    for w in wanted:
        for i, c in enumerate(low):
            if w in c:
                return i
    return 0


UPLOAD_TYPES = ["csv", "tsv", "txt", "xlsx", "xls"]


def _read_uploaded_table(upload, index_col=None) -> pd.DataFrame:
    """Read an uploaded CSV/TSV/TXT/XLSX/XLS file into a DataFrame."""
    name = upload.name.lower()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(upload, index_col=index_col)
    sep = "\t" if name.endswith((".tsv", ".txt")) else ","
    return pd.read_csv(upload, sep=sep, index_col=index_col)



def validate_clean_matrix(df: pd.DataFrame) -> list[str]:
    """Light sanity checks for an uploaded matrix (first column = gene id, rest =
    numeric sample columns). Catches the most common real mistake: extra
    title/description rows above the real header."""
    errors = []
    if df.shape[1] < 2:
        errors.append(f"Only {df.shape[1]} sample column(s) found after the first (gene "
                      "id) column. Need at least 2 samples for a differential comparison.")
    numeric = df.apply(pd.to_numeric, errors="coerce")
    if numeric.notna().sum().sum() / max(1, df.size) < 0.5:
        errors.append("Most values in this file aren't numeric. Check that the first row is "
                      "the real header (sample names) and the first column is the gene id. "
                      "a common cause is extra title/description rows above the actual data "
                      "(e.g. exported straight from a publication's supplementary table).")
    return errors


def _counts_columns(path) -> set | None:
    """Sample column names of a loaded matrix (all but the first/gene col)."""
    if not path:
        return None
    try:
        cols = pd.read_csv(path, nrows=0).columns.tolist()
        return {str(c) for c in cols[1:]}
    except Exception:  # noqa: BLE001
        return None


def _best_sample_col(md, cols: list[str], counts_cols: set) -> int:
    """Index of the metadata column whose values best overlap the counts headers."""
    best_i, best_n = 0, -1
    for i, c in enumerate(cols):
        n = len(set(md[c].astype(str)) & counts_cols)
        if n > best_n:
            best_i, best_n = i, n
    return best_i


# words that mark a level as the control / baseline group (-> reference side)
CONTROL_KEYWORDS = {"control", "healthy", "normal", "ctrl", "ctl", "baseline",
                    "wildtype", "wt", "vehicle", "untreated", "mock", "sham",
                    "hc", "naive", "placebo", "uninfected", "ctr"}


def _control_level_index(levels: list[str]) -> int:
    """Index of the level that looks like the control/baseline group, so it can
    default to the REFERENCE side. Token match first, then substring; else 0."""
    for i, lv in enumerate(levels):
        tokens = {t for t in re.split(r"[^a-z0-9]+", str(lv).lower()) if t}
        if tokens & CONTROL_KEYWORDS:
            return i
    for i, lv in enumerate(levels):
        if any(k in str(lv).lower() for k in
               ("control", "healthy", "normal", "baseline", "untreated", "vehicle")):
            return i
    return 0


def counts_match_note():
    """At load time: report which metadata column matches the matrix sample IDs,
    and stash it so the Design step defaults to it."""
    cc = _counts_columns(ss.get("counts_path"))
    sdf = ss.get("samples_df")
    if not cc or sdf is None:
        return
    cols = list(sdf.columns)
    best = cols[_best_sample_col(sdf, cols, cc)]
    n = len(set(sdf[best].astype(str)) & cc)
    if n >= 2:
        ss["sample_col_auto"] = best
        st.success(f"Sample IDs match metadata column `{best}` "
                   f"({n}/{len(cc)}). Design step will use it automatically.")
    else:
        st.warning("Sample IDs don't clearly match any metadata column "
                   f"(best `{best}`: {n}/{len(cc)}). Set the Sample ID column "
                   "manually in the Design step.")


# ============================================================================
# STEP RENDERERS
# ============================================================================
def render_input_step() -> None:
    ui.step_header(
        1, len(STEPS), "Input data",
        "Bring in a counts matrix from GEO or a local file. Raw integer counts "
        "work best; normalized values are flagged before analysis.")
    src = st.radio("Source", ["NCBI GEO accession", "Manual upload"],
                   horizontal=True, key="input_src")

    if src == "NCBI GEO accession":
        acc = st.text_input("GEO series accession", placeholder="e.g. GSE125583",
                            key="geo_acc")

        if st.button("Check GEO", disabled=not acc, key="check_geo"):
            gse = acc.strip().upper()
            if not is_valid_accession(gse):
                st.error(f"Invalid GEO series accession `{gse}`. Expected a format "
                         "like `GSE125583`.")
                st.stop()
            outdir = DATA / gse
            with st.spinner(f"Checking {gse}: raw-counts availability + metadata…"):
                info = check_counts(gse)
                p = run_cmd(["Rscript", str(ENGINE_R / "fetch_geo.R"),
                             "--accession", gse, "--outdir", str(outdir)])
            report, meta_df = {}, None
            rf = outdir / "fetch_report.json"
            if rf.exists():
                report = json.loads(rf.read_text())
                if report.get("metadata_ok") and (outdir / "metadata.csv").exists():
                    meta_df = pd.read_csv(outdir / "metadata.csv")
            if meta_df is not None:
                ss["samples_df"] = meta_df
            ss["geo"] = {"gse": gse, "outdir": str(outdir), "counts": info,
                                  "report": report, "meta_loaded": meta_df is not None}
            ss["geo_counts_loaded"] = None
            ss["geo_err"] = None if (meta_df is not None or info.available) \
                else (p.stderr or "")[-1500:]

        g = ss.get("geo")
        if g:
            info, report = g["counts"], g["report"]
            if g["meta_loaded"]:
                st.success(f"Metadata for {report.get('n_samples', '?')} samples loaded.")
                st.write("**Candidate columns for condition/batch:**",
                         report.get("candidate_columns"))
                st.dataframe(ss["samples_df"], width="stretch", height=240)
            else:
                st.warning("Metadata not fetched. You can still upload a sample sheet "
                           "in Manual upload mode.")

            st.subheader("Raw counts matrix")
            if info.available:
                st.success("Standardized raw counts matrix **available**: NCBI "
                           f"{info.organism} ({info.build}), {info.size_mb} MB · "
                           f"`{info.filename}`")
                choice = st.radio("Do you want to use this counts matrix?",
                                  ["Proceed with it", "Keep / don't use it"],
                                  horizontal=True, key="geo_choice")
                if choice == "Proceed with it":
                    if st.button("Load NCBI counts", key="load_ncbi", icon=":material/download:"):
                        with st.spinner("Downloading + mapping Entrez IDs → symbols…"):
                            try:
                                res = download_counts(info, g["outdir"])
                                ss["counts_path"] = res["path"]
                                ss["geo_counts_loaded"] = res
                            except Exception as e:  # noqa: BLE001
                                st.error(f"Download failed: {e}")
                    if ss.get("geo_counts_loaded"):
                        res = ss["geo_counts_loaded"]
                        ids = "gene symbols" if res["symbols_mapped"] else "Entrez IDs"
                        st.success(f"Counts loaded: {res['n_genes']:,} genes × "
                                   f"{res['n_samples']} samples ({ids}).")
                        st.caption("Raw counts preview (genes × samples):")
                        st.dataframe(res["preview"], width="stretch")
                        counts_match_note()
                else:
                    st.caption("Not using the NCBI counts. Upload your own below.")
            else:
                st.warning("No NCBI standardized raw counts matrix for this GSE "
                           "(typical for microarray or older RNA-seq series).")
                likely = report.get("likely_counts_files") or []
                urls = report.get("likely_counts_urls") or []
                supp = report.get("supp_files") or []
                if isinstance(likely, str):          # auto_unbox safety
                    likely, urls = [likely], [urls]
                if likely:
                    st.info("Author supplementary file(s) that look like counts. "
                            "download & load directly (no manual step):")
                    for name, url in zip(likely, urls):
                        if st.button(f"Download & load {name}", key=f"supp_{name}", icon=":material/download:"):
                            with st.spinner(f"Downloading {name} from GEO…"):
                                try:
                                    res = download_supp(url, g["gse"], g["outdir"])
                                    ss["counts_path"] = res["path"]
                                    ss["supp_loaded"] = res
                                except Exception as e:  # noqa: BLE001
                                    st.error(f"Download failed: {e}")
                    if ss.get("supp_loaded"):
                        res = ss["supp_loaded"]
                        st.success(f"Loaded `{res['filename']}`: {res['n_rows']:,} "
                                   f"rows × {res['n_cols']} sample columns.")
                        if not res["looks_integer"]:
                            st.warning("Values are **not integers**. This file looks "
                                       "like normalized expression (TPM/FPKM/CPM), not raw "
                                       "counts. DESeq2/edgeR need raw counts; limma-voom or "
                                       "rounding may work, but prefer a raw-count file.")
                        st.caption("Raw counts preview (first column = gene id, "
                                   "remaining columns = samples):")
                        st.dataframe(res["preview"], width="stretch")
                        counts_match_note()
                elif supp:
                    st.caption(f"Author supplementary files (not obviously counts): {supp}")

            with st.expander("Raw check report"):
                st.json({"counts_available": info.available, "organism": info.organism,
                         "build": info.build, "file": info.filename,
                         "size_mb": info.size_mb, **report})

        if ss.get("geo_err"):
            with st.expander("GEO metadata fetch: details"):
                st.code(ss["geo_err"])

        st.divider()
        st.subheader("Or upload a counts matrix")
        st.caption("If NCBI counts aren't available or you prefer your own "
                   "(genes × samples, first column = gene id).")

    up = st.file_uploader("Counts matrix (CSV/TSV/Excel: genes × samples)", type=UPLOAD_TYPES, key="matrix_up")
    if up is not None:
        df = _read_uploaded_table(up, index_col=0)
        errors = validate_clean_matrix(df)
        if errors:
            st.error("This file doesn't look like a valid gene × sample matrix:\n\n"
                     + "\n".join(f"- {e}" for e in errors))
        else:
            cp = DATA / "matrix_uploaded.csv"
            df.to_csv(cp)
            ss["counts_path"] = str(cp)
            st.success(f"Matrix loaded: {df.shape[0]} rows × {df.shape[1]} samples")
            st.caption("Preview (rows × samples):")
            st.dataframe(df.head(), width="stretch")
            counts_match_note()

    if src == "Manual upload":
        sup = st.file_uploader("Sample sheet (sample, condition, [batch])",
                               type=UPLOAD_TYPES, key="samples_up")
        if sup is not None:
            sdf = _read_uploaded_table(sup)
            ss["samples_df"] = sdf
            st.success(f"Sample sheet loaded: {len(sdf)} samples")
            st.dataframe(sdf, width="stretch")


def render_design_step() -> None:
    ui.step_header(
        2, len(STEPS), "Experimental design",
        "Match the sample sheet to the count matrix, choose the comparison, "
        "and put batch in the model if it belongs there.")
    sdf = ss.get("samples_df")
    if sdf is None:
        st.warning("Load a sample sheet or GEO metadata in the Input step first.")
        return

    md = sdf.copy()
    st.subheader("Sample sheet")
    st.caption(f"{len(md)} samples · {len(md.columns)} metadata columns")
    st.dataframe(md, width="stretch")
    cols = list(md.columns)

    # default the Sample ID column to whichever metadata column actually
    # matches the loaded matrix headers (avoids "fewer than 2 samples match")
    counts_cols = _counts_columns(ss.get("counts_path"))
    auto = ss.get("sample_col_auto")
    if auto in cols:                                # detected at load time
        sample_idx = cols.index(auto)
    elif counts_cols:
        sample_idx = _best_sample_col(md, cols, counts_cols)
    else:
        sample_idx = _guess(cols, ["sample", "geo_accession", "title"])
    st.subheader("Match samples to the count matrix")
    sample_col = st.selectbox(
        "Sample ID column", cols, index=sample_idx,
        help="Must contain the SAME ids as your matrix's column headers.",
        key="sample_col")

    if counts_cols:
        overlap = len(set(md[sample_col].astype(str)) & counts_cols)
        if overlap >= 2:
            st.success(f"{overlap}/{len(counts_cols)} matrix columns match "
                       f"`{sample_col}`. Good to go.")
        else:
            best = cols[_best_sample_col(md, cols, counts_cols)]
            st.error(f"Only {overlap}/{len(counts_cols)} matrix columns match "
                     f"`{sample_col}`. The engine needs the column whose values are "
                     f"the matrix headers. Try **`{best}`**.")
    st.subheader("Comparison")
    cond_col = st.selectbox(
        "Condition / group column", cols,
        index=_guess(cols, ["condition", "diagnosis", "disease", "status",
                            "group", "treatment", "genotype", "phenotype",
                            "source_name"]),
        help="The variable that defines your comparison groups (e.g. diagnosis).",
        key="cond_col")
    use_batch = st.checkbox("Model a batch variable (recommended if batches exist)",
                            key="use_batch")
    batch_col = (st.selectbox("Batch column", [c for c in cols if c != cond_col],
                              key="batch_sel")
                 if use_batch else None)

    levels = sorted(md[cond_col].dropna().astype(str).unique().tolist())
    if len(levels) < 2:
        st.error(f"Condition column '{cond_col}' has <2 groups.")
        return

    # default the control/healthy group to the REFERENCE side, disease to test
    ref_idx = _control_level_index(levels)
    test_idx = next((i for i in range(len(levels)) if i != ref_idx), 0)
    with st.container(border=True):
        c1, c2 = st.columns(2)
        test_lvl = c1.selectbox("Test group (numerator) · disease / treatment",
                                levels, index=test_idx, key="test_lvl")
        ref_lvl = c2.selectbox("Reference group (denominator) · control / baseline",
                               levels, index=ref_idx, key="ref_lvl")
        st.caption(f"Comparing **{test_lvl}** vs **{ref_lvl}** → positive log2FC = "
                   f"higher in **{test_lvl}**. Keep control/healthy as the reference.")

    if test_lvl == ref_lvl:
        st.error("Test and reference must be different groups.")
        st.stop()

    # build standardized sample sheet for the engine
    std = pd.DataFrame({"sample": md[sample_col].astype(str),
                        "condition": md[cond_col].astype(str)})
    if batch_col:
        std["batch"] = md[batch_col].astype(str)
    std_path = DATA / "samples_std.csv"
    std.to_csv(std_path, index=False)
    ss["samples_std"] = str(std_path)
    ss["contrast"] = f"{test_lvl},{ref_lvl}"
    ss["batch_col"] = "batch" if batch_col else None
    st.success(f"Design set: {test_lvl} vs {ref_lvl}"
               + (f", batch = {batch_col}" if batch_col else ""))
    st.caption("Standardized sample sheet sent to the engine")
    st.dataframe(std, width="stretch")


def render_rna_run_step() -> None:
    ui.step_header(
        3, len(STEPS), "Run the analysis",
        "Choose the significance convention and gene filter, then run all three "
        "methods against the same design.")
    engine = st.radio("Analysis engine", ["R (Bioconductor)", "Python (experimental)"],
                      horizontal=True, key="rna_engine")

    with st.container(border=True):
        metric_label = st.radio(
            "Significance metric",
            ["Raw p-value (lab canonical)", "Adjusted p / FDR"],
            horizontal=True, key="rna_metric",
            help="Your DE_analysis.qmd uses raw P.Value < 0.05 & |logFC| > 0.58.")
        sig_metric = "pvalue" if metric_label.startswith("Raw") else "padj"

        filt_label = st.radio(
            "Gene filter (shared by all 3 methods)",
            ["Canonical: TMM, log2CPM > 1 in ≥ min-group-n", "Count-based: counts ≥ 10"],
            horizontal=True, key="rna_filter")
        filter_mode = "canonical" if filt_label.startswith("Canonical") else "count"

        groups_label = st.radio(
            "Samples to model (for the pairwise contrast)",
            ["Only the two contrasted groups (recommended)",
             "All groups in the condition column"],
            horizontal=True, key="rna_groups",
            help="Subsetting to the two groups matches a typical standalone two-group "
                 "analysis. Modeling all groups shares variance across all samples.")
        groups_mode = "subset" if groups_label.startswith("Only") else "all"

        c1, c2, c3 = st.columns(3)
        default_lfc = 0.58 if sig_metric == "pvalue" else 1.0
        padj = c1.number_input("Significance cutoff", 0.0, 1.0, 0.05, 0.01, key="rna_padj")
        lfc = c2.number_input("|log2FC| cutoff", 0.0, 10.0, default_lfc, 0.01, key="rna_lfc")
        topn = c3.number_input("Top genes (heatmap)", 10, 200, 50, 10, key="rna_topn")

    ready = ss.counts_path and ss.get("samples_std") and ss.get("contrast")
    if not ready:
        st.warning("Finish Input and Design first: counts matrix, sample sheet "
                   "and a contrast are required.")
    if engine.startswith("Python"):
        st.info("The Python engine (PyDESeq2 + inmoose) is stubbed for this MVP; "
                "the R engine is fully wired. Use R for now.")

    if st.button("Run 3-method DE", disabled=not ready or engine.startswith("Python"),
                key="rna_run_btn", type="primary", icon=":material/play_arrow:"):
        run_dir = RESULTS / "run"
        cmd = ["Rscript", str(ENGINE_R / "run_de.R"),
               "--counts", ss.counts_path, "--samples", ss["samples_std"],
               "--contrast", ss["contrast"], "--outdir", str(run_dir),
               "--padj", str(padj), "--lfc", str(lfc), "--topn", str(int(topn)),
               "--sig_metric", sig_metric, "--filter", filter_mode,
               "--groups_mode", groups_mode]
        if ss.get("batch_col"):
            cmd += ["--batch", ss["batch_col"]]
        with st.spinner("Running DESeq2 · edgeR · limma-voom…"):
            p = run_cmd(cmd)
        if p.returncode == 0:
            ss.run_dir = str(run_dir)
            st.success("Analysis complete. Open Results to review the three methods.")
            st.code(p.stderr[-1500:] or "ok")
        else:
            st.error("Engine failed.")
            st.code(p.stderr[-3000:])


def render_results_step() -> None:
    ui.step_header(
        4, len(STEPS), "Results",
        "Compare hit counts across methods, check the batch structure, then read "
        "the volcanoes, heatmap and enrichment.")
    methods = METHODS
    run_dir = ss.get("run_dir")
    if not run_dir:
        st.warning("Run the analysis in the Run step first.")
        return

    rd = Path(run_dir)
    summary = json.loads((rd / "comparison_summary.json").read_text())
    batch = json.loads((rd / "batch_report.json").read_text())
    st.subheader("Method comparison" if len(methods) > 1 else "Summary")
    metric_name = "raw p" if summary.get("sig_metric") == "pvalue" else "FDR"
    st.caption(f"Significance: {metric_name} < {summary.get('sig_cut')} & "
               f"|log2FC| > {summary.get('lfc_cut')}  ·  "
               f"filter: {summary.get('filter_mode')}  ·  "
               f"{summary.get('n_genes_tested')} genes tested")
    cols = st.columns(len(methods))
    for i, m in enumerate(methods):
        s = summary["methods"][m]
        cols[i].metric(
            m, f"{s['total_sig']:,} DE",
            f"{s['up']:,} up · {s['down']:,} down",
            delta_color="off", border=True)
    if len(methods) > 1:
        st.caption(f"Significant in all {len(methods)} methods: "
                   f"{summary['significant_in_all']:,} · in any: "
                   f"{summary['significant_in_any']:,}")

    st.subheader("Batch-effect check")
    if batch["batch_effect_suspected"]:
        st.error("Batch effect suspected: batch explains a top principal "
                 "component. It is being modeled in the design.")
    elif batch["has_batch_column"]:
        st.success("Batch modeled; no dominant batch effect on top PCs.")
    else:
        st.info("No batch column provided.")
    with st.container(border=True):
        st.image(str(rd / "pca.png"), width="stretch",
                 caption="Principal components of normalized expression")
    with st.expander("PC ↔ metadata associations"):
        st.dataframe(pd.DataFrame(batch["associations"]))

    st.subheader("Volcano plot" + ("s" if len(methods) > 1 else ""))
    vc = st.columns(len(methods))
    for i, m in enumerate(methods):
        with vc[i].container(border=True):
            st.image(str(rd / f"volcano_{m}.png"), width="stretch", caption=m)

    st.subheader("Heatmap (top genes)")
    hm = rd / "heatmap_top.png"
    if hm.exists():
        with st.container(border=True):
            st.image(str(hm), width="stretch")

    st.subheader("Cross-method table" if len(methods) > 1 else "Results table")
    comp = pd.read_csv(rd / "comparison.csv")
    only_sig = st.checkbox("Show only genes significant in ≥1 method",
                           value=True, key="only_sig")
    view = comp[comp["n_methods_sig"] >= 1] if only_sig else comp
    st.dataframe(view.sort_values("n_methods_sig", ascending=False),
                 width="stretch", height=400)
    st.download_button("Download comparison.csv", comp.to_csv(index=False),
                       "comparison.csv", key="dl_comparison",
                       icon=":material/download:")

    # ---- Functional enrichment (clusterProfiler: GO / KEGG / Reactome) ---
    st.subheader("Functional enrichment")
    st.caption("ORA via **clusterProfiler** (hypergeometric test), with the "
               "tested genes as background. GO uses org.Hs/Mm.eg.db; KEGG/Reactome "
               "query online.")
    dir_cols = [f"{m}_dir" for m in methods]
    ENRICH_DBS = {"GO Biological Process": "GO_BP",
                  "GO Molecular Function": "GO_MF",
                  "GO Cellular Component": "GO_CC",
                  "KEGG": "KEGG", "Reactome": "Reactome"}

    def _gene_set(source: str, direction: str) -> list[str]:
        if source == "Consensus (≥2 methods)":
            up = (comp[dir_cols] == "up").sum(axis=1) >= 2
            dn = (comp[dir_cols] == "down").sum(axis=1) >= 2
        else:
            up = comp[f"{source}_dir"] == "up"
            dn = comp[f"{source}_dir"] == "down"
        mask = {"Up": up, "Down": dn, "Both": up | dn}[direction]
        return comp.loc[mask, "gene"].astype(str).tolist()

    geo_org = ""
    if isinstance(ss.get("geo"), dict):
        geo_org = getattr(ss["geo"].get("counts"), "organism", "")
    source_opts = (["Consensus (≥2 methods)"] + methods) if len(methods) > 1 else methods
    e1, e2, e3 = st.columns(3)
    e_source = e1.selectbox("Gene set", source_opts, key="enr_source")
    e_dir = e2.radio("Direction", ["Down", "Up", "Both"], horizontal=True, key="enr_dir")
    e_org = e3.selectbox("Organism", ["human", "mouse"],
                         index=1 if geo_org == "mouse" else 0, key="enr_org")
    e_libs = st.multiselect("Databases", list(ENRICH_DBS),
                            default=["GO Biological Process", "KEGG", "Reactome"],
                            key="enr_libs")
    e_genes = _gene_set(e_source, e_dir)
    st.caption(f"{len(e_genes)} genes · {e_dir.lower()} · {e_source}  "
               f"(background: {len(comp):,} tested)")

    use_enrichr = not has_r_package("clusterProfiler")
    if use_enrichr:
        st.info("clusterProfiler isn't installed on this server, so enrichment uses the "
                "Enrichr web service (same GO / KEGG / Reactome libraries, no custom background).")
    if st.button("Run enrichment", disabled=not e_genes or not e_libs, key="enr_run") and use_enrichr:
        with st.spinner("Querying Enrichr…"):
            libs = {d: ENRICHR_LIBS[d] for d in e_libs}
            try:
                res = enrichr_enrich(e_genes, libs)
                ss["enrich"] = {d: df[df["Adj P-value"] < 0.05] if "Adj P-value" in df else pd.DataFrame()
                                for d, df in res.items()}
                ss["enrich_meta"] = {"n": len(e_genes), "dir": e_dir, "source": e_source,
                                     "org": e_org, "engine": "Enrichr"}
            except requests.RequestException as e:
                st.error(f"Enrichr request failed: {e}")
    elif ss.get("enr_run"):
        edir = rd / "enrichment"
        edir.mkdir(exist_ok=True)
        (edir / "genes.txt").write_text("\n".join(e_genes))
        (edir / "universe.txt").write_text("\n".join(comp["gene"].astype(str)))
        cmd = ["Rscript", str(ENGINE_R / "enrich.R"),
               "--genes", str(edir / "genes.txt"),
               "--universe", str(edir / "universe.txt"),
               "--organism", e_org, "--outdir", str(edir),
               "--databases", ",".join(ENRICH_DBS[d] for d in e_libs)]
        with st.spinner("Running clusterProfiler (GO / KEGG / Reactome)…"):
            p = run_cmd(cmd, timeout=600)
        if p.returncode == 0:
            results = {}
            for disp in e_libs:
                f = edir / f"{ENRICH_DBS[disp]}.csv"
                results[disp] = pd.read_csv(f) if f.exists() else pd.DataFrame()
            ss["enrich"] = results
            ss["enrich_meta"] = {"n": len(e_genes), "dir": e_dir,
                                          "source": e_source, "org": e_org}
        else:
            st.error("Enrichment failed.")
            st.code(p.stderr[-2500:])

    er = ss.get("enrich")
    if er:
        meta = ss.get("enrich_meta", {})
        st.caption(f"{meta.get('engine', 'clusterProfiler')} ORA: {meta.get('n')} "
                   f"{str(meta.get('dir', '')).lower()} genes "
                   f"({meta.get('source')}, {meta.get('org')}).")
        for tab, (lib, edf) in zip(st.tabs(list(er.keys())), er.items()):
            with tab:
                if edf is None or edf.empty:
                    st.warning("No enriched terms at FDR < 0.05.")
                else:
                    plot_df = edf.rename(columns={"AdjPvalue": "Adj P-value"})
                    fig = enr_figure(plot_df, lib, top_n=12)
                    st.pyplot(fig, width="stretch")
                    st.dataframe(edf, width="stretch", height=300,
                                 hide_index=True)
                    d1, d2 = st.columns(2)
                    d1.download_button("Download CSV", edf.to_csv(index=False),
                                       f"enrichment_{lib}.csv", "text/csv",
                                       key=f"enr_csv_{lib}",
                                       icon=":material/download:")
                    d2.download_button("Download PNG", enr_bytes(fig, "png"),
                                       f"enrichment_{lib}.png", "image/png",
                                       key=f"enr_png_{lib}",
                                       icon=":material/download:")

    render_software_versions(rd)


# What each package does in the portal, for the versions table / methods section.
R_SOFTWARE = {
    "R": "Language runtime for the analysis engine",
    "Bioconductor": "Release the Bioconductor packages come from",
    "DESeq2": "Differential expression (negative binomial GLM)",
    "edgeR": "Differential expression (quasi-likelihood); TMM normalization, filtering",
    "limma": "Differential expression (limma-voom)",
    "GEOquery": "GEO sample metadata and supplementary files",
    "clusterProfiler": "GO / KEGG over-representation analysis",
    "ReactomePA": "Reactome pathway enrichment",
    "AnnotationDbi": "Gene ID lookups",
    "org.Hs.eg.db": "Human gene annotation (Entrez ID to symbol, GO)",
    "org.Mm.eg.db": "Mouse gene annotation (Entrez ID to symbol, GO)",
    "ggplot2": "Volcano and PCA plots",
    "pheatmap": "Heatmap of top genes",
    "jsonlite": "Engine result files",
}
PY_SOFTWARE = {
    "streamlit": "Web portal",
    "pandas": "Tables, file parsing",
    "numpy": "Numerics",
    "matplotlib": "Gene-of-interest and enrichment plots",
    "requests": "GEO, NCBI and literature database queries",
    "openpyxl": "Excel (.xlsx) uploads",
}


def _py_version(pkg: str) -> str:
    try:
        return version(pkg)
    except PackageNotFoundError:
        return "not installed"


@st.cache_resource
def has_r_package(pkg: str) -> bool:
    """True when Rscript can load `pkg` (Streamlit Cloud's Debian R lacks some)."""
    try:
        return subprocess.run(["Rscript", "-e", f"quit(status = !requireNamespace('{pkg}', quietly = TRUE))"],
                              capture_output=True, timeout=60).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def render_software_versions(rd: Path) -> None:
    st.subheader("Software versions")
    st.caption("Exact versions used for this run, for reproducibility and the methods "
               "section of a paper. R versions are recorded by the engine at run time.")
    meta_path = rd / "run_meta.json"
    r_ver = json.loads(meta_path.read_text()).get("software", {}) if meta_path.exists() else {}
    rows = [{"Software": name, "Language": "R", "Version": r_ver.get(name, "re-run to record"),
             "Used for": why} for name, why in R_SOFTWARE.items()]
    rows.append({"Software": "Python", "Language": "Python",
                 "Version": platform.python_version(), "Used for": "Language runtime for the portal"})
    rows += [{"Software": pkg, "Language": "Python", "Version": _py_version(pkg), "Used for": why}
             for pkg, why in PY_SOFTWARE.items()]
    df = pd.DataFrame(rows)
    st.dataframe(df, width="stretch", hide_index=True, height=35 * (len(df) + 1) + 3)
    st.download_button("Download versions (CSV)", df.to_csv(index=False),
                       "software_versions.csv", "text/csv", key="dl_versions",
                       icon=":material/download:")


def render_gene_of_interest_step() -> None:
    ui.step_header(
        5, len(STEPS), "Gene of interest",
        "Pull one gene across all three methods and inspect its normalized "
        "expression by group.")
    if not ss.get("run_dir"):
        st.warning("Run the analysis in the Run step first.")
        return
    methods = METHODS

    rd = Path(ss["run_dir"])
    comp = pd.read_csv(rd / "comparison.csv")
    expr = pd.read_csv(rd / "normalized_expression.csv", index_col=0)
    samples = pd.read_csv(ss["samples_std"])

    gene = st.selectbox("Gene of interest", comp["gene"].tolist(), key="goi_gene")
    ss["gene_of_interest"] = gene
    row = comp[comp["gene"] == gene].iloc[0]

    st.subheader(f"{gene}: across methods")
    tbl = pd.DataFrame({
        "method": methods,
        "log2FC": [row[f"{m}_log2FC"] for m in methods],
        "padj": [row[f"{m}_padj"] for m in methods],
        "direction": [row[f"{m}_dir"] for m in methods],
    })
    st.dataframe(tbl, width="stretch")

    st.subheader(f"{gene}: expression profile")

    # condition order: reference first, then test (from the chosen contrast)
    test_lvl, ref_lvl = (ss.get("contrast", "B,A").split(",") + ["", ""])[:2]
    conditions = [c for c in (ref_lvl, test_lvl) if c in set(samples["condition"])]
    if not conditions:
        conditions = sorted(samples["condition"].dropna().unique().tolist())

    c1, c2 = st.columns([2, 1])
    scheme = c1.selectbox("Color scheme", list(PALETTES.keys()), key="goi_scheme")
    font_size = c2.slider("Font size", 10, 28, 16, 1, key="goi_font")

    # load each method's native normalized expression for this gene
    per_method, missing = {}, []
    for m in methods:
        f = rd / f"norm_{m}.csv"
        mat = pd.read_csv(f, index_col=0) if f.exists() else expr
        if gene in mat.index:
            long = mat.loc[gene].reset_index()
            long.columns = ["sample", "expression"]
            per_method[m] = long.merge(samples, on="sample", how="left")
        else:
            missing.append(m)

    if not per_method:
        st.info("Gene not in the normalized expression matrices (filtered out).")
        return

    fig = build_gene_figure(gene, per_method, conditions, scheme, font_size)
    with st.container(border=True):
        st.pyplot(fig, width="stretch")
    if missing:
        st.caption(f"Not available for: {', '.join(missing)} (filtered out).")

    long_all = pd.concat(
        [d.assign(method=m) for m, d in per_method.items()], ignore_index=True)
    d1, d2, d3, d4 = st.columns(4)
    d1.download_button("Data (CSV)", long_all.to_csv(index=False),
                       f"{gene}_expression.csv", "text/csv", key="goi_dl_csv",
                       icon=":material/download:")
    d2.download_button("PNG", figure_bytes(fig, "png"),
                       f"{gene}_expression.png", "image/png", key="goi_dl_png",
                       icon=":material/download:")
    d3.download_button("SVG", figure_bytes(fig, "svg"),
                       f"{gene}_expression.svg", "image/svg+xml", key="goi_dl_svg",
                       icon=":material/download:")
    d4.download_button("JPEG", figure_bytes(fig, "jpeg"),
                       f"{gene}_expression.jpg", "image/jpeg", key="goi_dl_jpg",
                       icon=":material/download:")


def render_literature_step() -> None:
    ui.step_header(
        6, len(STEPS), "Literature",
        "Search PubMed, PubMed Central, Europe PMC and OpenAlex for the selected gene, "
        "then optionally synthesize the retrieved evidence.")
    st.caption("Multi-database search → de-duplicate → rank by citations → flag "
               "Tier-1 venues and classify impact.")
    default_gene = ss.get("gene_of_interest", "") or ""
    with st.container(border=True):
        c1, c2, c3, c4 = st.columns([3, 3, 1, 1])
        gene = c1.text_input("Gene", value=default_gene,
                             help="Must appear in each paper: in the title/abstract, "
                                  "or anywhere in the text (see scope below).")
        context = c2.text_input("Context (optional: disease / tissue)",
                                placeholder="e.g. Alzheimer")
        year_from = c3.number_input("From year", 1990, 2026, 2015, 1)
        limit = c4.number_input("Max papers", 10, 50, 25, 5)
        scope = st.radio(
            "Where to look for the gene",
            ["Full text (title, abstract and body)", "Title and abstract only"],
            horizontal=True,
            help="Full text also searches inside article bodies (PubMed Central, "
                 "Europe PMC open-access papers, OpenAlex). Finds many more papers, "
                 "including ones that only mention the gene in passing. Title and "
                 "abstract only keeps papers that are mainly about the gene.")
    full_text = scope.startswith("Full")

    if st.button("Search literature", disabled=not gene,
                 icon=":material/search:"):
        with st.spinner("Searching PubMed, PubMed Central, Europe PMC and OpenAlex…"):
            try:
                ss["lit_result"] = search_literature(
                    gene, context, limit=int(limit), year_from=int(year_from),
                    full_text=full_text)
                ss["lit_synthesis"] = None
            except Exception as e:  # noqa: BLE001
                ss["lit_result"] = None
                st.error(f"Literature query failed: {e}")

    res = ss.get("lit_result")
    if res:
        for note in res.notes:
            st.warning(note)

        if res.reported:
            hits = "  ·  ".join(f"{db}: {n:,}" for db, n in res.total_by_db.items() if n)
            st.success(f"Reported in literature: {hits}")
        elif not res.databases:
            st.error("No database responded. Check your internet connection and "
                     "try again (see warnings above).")
        else:
            st.warning("No matching papers; possibly novel in this context.")

        m1, m2, m3 = st.columns(3)
        m1.metric("Databases searched", len(res.databases), border=True)
        m2.metric("Papers retrieved", len(res.papers), border=True)
        m3.metric("High-impact", res.n_high_impact, border=True)
        st.caption(f"Query: `{res.query}`  ·  databases: {', '.join(res.databases)}")
        # Google Scholar has no API (and blocks automated queries), so link out to it.
        scholar_q = f'"{res.gene}" {res.context}'.strip()
        st.link_button("Open this search in Google Scholar",
                       "https://scholar.google.com/scholar?"
                       + urlencode({"q": scholar_q, "as_ylo": int(year_from)}),
                       icon=":material/open_in_new:")

        if res.papers:
            df = pd.DataFrame([{
                "Citations": p.citations, "Year": p.year, "Title": p.title,
                "Venue": p.venue, "Tier-1": "★" if p.tier1 else "",
                "Impact": p.classification,
                "Preprint": "✓" if p.is_preprint else "",
                "Source": p.source, "URL": p.url,
            } for p in res.papers])
            st.dataframe(
                df, width="stretch", height=440, hide_index=True,
                column_config={"URL": st.column_config.LinkColumn("Link",
                               display_text="open")})
            st.download_button("Download citations (CSV)", df.to_csv(index=False),
                               f"{res.gene}_literature.csv", "text/csv",
                               icon=":material/download:")

        st.divider()
        st.markdown("**AI evidence synthesis**: summarizes the retrieved papers "
                    "using your Claude subscription via `claude -p` (no API key).")
        if not res.papers:
            st.caption("No papers to synthesize; the gene appears unstudied in this "
                       "context (try removing the context or lowering the year).")
        elif not llm_available():
            st.caption("`claude` CLI not found on PATH; install/login to enable.")
        elif st.button("Synthesize evidence"):
            with st.spinner("Asking Claude to synthesize the evidence…"):
                try:
                    ss["lit_synthesis"] = synthesize_literature(
                        res.gene, res.context, res.papers)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Synthesis failed: {e}")
        if ss.get("lit_synthesis"):
            st.info(ss["lit_synthesis"])


# ============================================================================
# STEPS + sidebar + dispatch
# ============================================================================
def _go(target):
    """Navigate to a step. Runs as a button on_click callback (i.e. BEFORE the
    next rerun instantiates the sidebar radio), so writing its key is allowed."""
    ss["nav"] = target


if "nav" not in ss or ss["nav"] not in STEPS:
    ss["nav"] = STEPS[0]

ui.sidebar_brand()
step = st.sidebar.radio("Workflow", STEPS, key="nav",
                        label_visibility="collapsed")
ui.sidebar_status(
    counts_label="Loaded" if ss.counts_path else "Not loaded",
    samples_label=(f"{len(ss.samples_df)} samples"
                   if ss.samples_df is not None else "Not loaded"),
    results_label="Ready" if ss.run_dir else "Not run",
    contrast_label=ui.contrast_label(ss.get("contrast")),
    batch_label="Modeled" if ss.get("batch_col") else None,
)

ui.institution_band()
ui.hero(
    sample_count=len(ss.samples_df) if ss.samples_df is not None else None,
    counts_loaded=bool(ss.counts_path),
    run_ready=bool(ss.run_dir),
    contrast=ss.get("contrast"),
)


def nav_buttons():
    """Back / Next buttons rendered at the bottom of every step."""
    i = STEPS.index(step)
    ui.nav_progress(i, len(STEPS), STEPS[i])
    with st.container(key="nav_row"):
        cols = st.columns([1, 1.5])
        if i > 0:
            cols[0].button("Back", on_click=_go, args=(STEPS[i - 1],),
                           width="stretch", key=f"back_{i}",
                           icon=":material/arrow_back:")
        if i < len(STEPS) - 1:
            cols[1].button(f"Next: {STEPS[i + 1]}", type="primary",
                           on_click=_go, args=(STEPS[i + 1],),
                           width="stretch", key=f"next_{i}",
                           icon=":material/arrow_forward:", icon_position="right")


i = STEPS.index(step)

if i == 0:
    render_input_step()
elif i == 1:
    render_design_step()
elif i == 2:
    render_rna_run_step()
elif i == 3:
    render_results_step()
elif i == 4:
    render_gene_of_interest_step()
elif i == 5:
    render_literature_step()


# Back / Next navigation at the bottom of every step
nav_buttons()
ui.footer()

# Streamlit keeps the scroll position across reruns; jump to the top whenever the
# step changes (Next/Back buttons or the sidebar).
if ss.get("_shown_step") not in (None, step):
    ui.scroll_to_top(step)
ss["_shown_step"] = step
