"""Functional enrichment (over-representation) via the Enrichr API.

Takes a list of gene SYMBOLS (e.g. the up- or down-regulated DEGs) and tests them
against curated gene-set libraries (GO, KEGG, Reactome, WikiPathways). Enrichr is
symbol-based and web-hosted, so no local annotation databases (org.Hs.eg.db) or
Entrez mapping are needed. Human/mouse symbols supported.

Ref: Chen et al. 2013 (BMC Bioinformatics); https://maayanlab.cloud/Enrichr/
"""
from __future__ import annotations

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

ENRICHR = "https://maayanlab.cloud/Enrichr"

# display name -> Enrichr library id
LIBRARIES = {
    "GO Biological Process": "GO_Biological_Process_2023",
    "KEGG": "KEGG_2021_Human",
    "Reactome": "Reactome_2022",
    "GO Molecular Function": "GO_Molecular_Function_2023",
    "WikiPathways": "WikiPathway_2023_Human",
    "GO Cellular Component": "GO_Cellular_Component_2023",
}


def _add_list(genes: list[str], timeout: int) -> int:
    payload = "\n".join(sorted({g.strip() for g in genes if g and g.strip()}))
    r = requests.post(
        f"{ENRICHR}/addList",
        files={"list": (None, payload), "description": (None, "rna-seq-app")},
        timeout=timeout)
    r.raise_for_status()
    return r.json()["userListId"]


def _enrich_library(user_list_id: int, library: str, timeout: int) -> pd.DataFrame:
    r = requests.get(f"{ENRICHR}/enrich",
                     params={"userListId": user_list_id, "backgroundType": library},
                     timeout=timeout)
    r.raise_for_status()
    rows = []
    for e in r.json().get(library, []):
        # e = [rank, term, pvalue, zscore, combined, [genes], adj_p, old_p, old_adj_p]
        rows.append({"Term": e[1], "P-value": e[2], "Adj P-value": e[6],
                     "Combined score": round(e[4], 1),
                     "N genes": len(e[5]), "Overlap genes": ", ".join(e[5])})
    df = pd.DataFrame(rows)
    return df.sort_values("Adj P-value").reset_index(drop=True) if not df.empty else df


def enrich(genes: list[str], libraries: dict[str, str], timeout: int = 30) -> dict:
    """Run Enrichr ORA for `genes` against each library in `libraries`
    ({display_name: library_id}). Returns {display_name: DataFrame}."""
    if not genes:
        return {}
    uid = _add_list(genes, timeout)
    out = {}
    for name, lib in libraries.items():
        try:
            out[name] = _enrich_library(uid, lib, timeout)
        except Exception as e:  # noqa: BLE001
            out[name] = pd.DataFrame({"error": [str(e)]})
    return out


def top_terms_figure(df: pd.DataFrame, title: str, top_n: int = 12,
                     color: str = "#1f77b4", font_size: int = 12):
    """Horizontal bar of the top terms by -log10(adjusted p)."""
    d = df.head(top_n).iloc[::-1]
    y = np.arange(len(d))
    vals = -np.log10(d["Adj P-value"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.5 * len(d) + 1)))
    ax.barh(y, vals, color=color, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([t[:60] for t in d["Term"]], fontsize=font_size)
    ax.set_xlabel("-log10(adjusted p)", fontsize=font_size)
    ax.set_title(title, fontsize=font_size + 2, fontweight="bold")
    ax.axvline(-np.log10(0.05), ls="--", color="grey", lw=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


def figure_bytes(fig, fmt: str = "png") -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=200, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.getvalue()


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path:
        genes = [l.strip() for l in open(path) if l.strip()]
    else:
        genes = ["SST", "BDNF", "VGF", "ADCYAP1", "SNAP25", "SYT1", "GAD1",
                 "PVALB", "CALB1", "NEFL", "SYP", "GABRA1", "GRIN1", "DLG4"]
    res = enrich(genes, {k: LIBRARIES[k] for k in ["GO Biological Process", "KEGG"]})
    for lib, df in res.items():
        print(f"\n=== {lib} ({len(df)} terms) ===")
        if not df.empty and "error" not in df.columns:
            print(df.head(6)[["Term", "Adj P-value", "N genes"]].to_string(index=False))
        else:
            print(df)
