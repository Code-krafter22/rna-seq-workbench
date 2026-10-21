"""PubMed literature search for a gene of interest via NCBI E-utilities.

Used by the Streamlit app to answer: "Is this gene reported in the literature
(optionally in a given context/disease), and what are the top papers?"

No API key required for low volume; set NCBI_API_KEY env var to raise rate limits.
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass, field

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY")


@dataclass
class Paper:
    pmid: str
    title: str
    journal: str
    year: str

    @property
    def url(self) -> str:
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"


@dataclass
class LitResult:
    gene: str
    context: str
    query: str
    total: int
    papers: list[Paper] = field(default_factory=list)

    @property
    def reported(self) -> bool:
        return self.total > 0


def _params(extra: dict) -> dict:
    p = {"db": "pubmed", "retmode": "json", **extra}
    if API_KEY:
        p["api_key"] = API_KEY
    return p


def search_gene(gene: str, context: str = "", retmax: int = 10,
                timeout: int = 20) -> LitResult:
    """Search PubMed for a gene, optionally narrowed by a context term
    (e.g. a disease or tissue). Returns total hit count + top papers."""
    gene = gene.strip()
    terms = [f'"{gene}"[Title/Abstract] OR "{gene}"[Gene]']
    if context.strip():
        terms.append(f'"{context.strip()}"[Title/Abstract]')
    query = " AND ".join(f"({t})" for t in terms)

    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params=_params({"term": query, "retmax": retmax,
                                     "sort": "relevance"}), timeout=timeout)
    r.raise_for_status()
    data = r.json().get("esearchresult", {})
    total = int(data.get("count", 0))
    ids = data.get("idlist", [])

    papers: list[Paper] = []
    if ids:
        time.sleep(0.34 if not API_KEY else 0.1)  # respect rate limits
        s = requests.get(f"{EUTILS}/esummary.fcgi",
                         params=_params({"id": ",".join(ids)}), timeout=timeout)
        s.raise_for_status()
        result = s.json().get("result", {})
        for pmid in result.get("uids", []):
            doc = result.get(pmid, {})
            papers.append(Paper(
                pmid=pmid,
                title=doc.get("title", "(no title)"),
                journal=doc.get("fulljournalname") or doc.get("source", ""),
                year=(doc.get("pubdate", "") or "").split(" ")[0],
            ))
    return LitResult(gene=gene, context=context, query=query,
                     total=total, papers=papers)


if __name__ == "__main__":
    import sys
    g = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    ctx = sys.argv[2] if len(sys.argv) > 2 else ""
    res = search_gene(g, ctx)
    print(f"Query: {res.query}")
    print(f"Reported: {res.reported}  (total hits: {res.total})")
    for p in res.papers:
        print(f"  [{p.year}] {p.title}  — {p.journal}\n    {p.url}")
