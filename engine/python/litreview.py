"""Multi-database literature search for a gene of interest.

Follows the methodology of the `literature-review` skill, but runs live inside
the app via free, key-less APIs:

  * Europe PMC  — indexes PubMed + PMC + preprints (bioRxiv/medRxiv), returns
                  citation counts (citedByCount) and DOIs.
  * OpenAlex — open index of ~250M works (journals, preprints, conference papers)
                  with citation counts; no key needed.

Results are merged, de-duplicated (by DOI, then normalized title), ranked by
citation count, classified by the skill's age/citation thresholds, and flagged
when published in a Tier-1 venue. This mirrors the skill's "prioritize
high-impact papers" + "multi-database + dedupe + citation ranking" steps in a
form a deployed app can call without an agent.
"""
from __future__ import annotations

import datetime
import os
import re
import time
from dataclasses import dataclass, field

import requests

CURRENT_YEAR = datetime.date.today().year

# Optional NCBI key raises PubMed E-utilities rate limits (not required).
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")

# A shared session with a descriptive User-Agent. Some APIs throttle or block the
# default "python-requests" UA; a real UA avoids intermittent blocks.
_SESSION = requests.Session()
_SESSION.headers.update({
    "User-Agent": "rna-seq-app/1.0 (bulk RNA-seq DE; gene literature search)",
    "Accept": "application/json",
})


def _get(url: str, params: dict, timeout: int, retries: int = 3, headers: dict | None = None):
    """GET with small exponential backoff so a transient blip doesn't wipe results."""
    last_exc = None
    for attempt in range(retries):
        try:
            return _SESSION.get(url, params=params, timeout=timeout, headers=headers)
        except requests.RequestException as e:
            last_exc = e
            time.sleep(0.6 * (attempt + 1))
    raise last_exc

# Tier-1 venues from the literature-review skill's prioritization guidance.
TIER1 = [
    "nature", "science", "cell", "new england journal", "lancet", "jama",
    "pnas", "proceedings of the national academy", "nature medicine",
    "nature biotechnology",
]


def is_tier1(venue: str) -> bool:
    v = (venue or "").lower()
    return any(t in v for t in TIER1)


def classify(citations: int | None, year: int | None) -> str:
    """Skill's impact classification by paper age + citation count."""
    if not citations or not year:
        return ""
    age = CURRENT_YEAR - year
    if age <= 3:
        if citations >= 100:
            return "Highly Influential"
        if citations >= 20:
            return "Noteworthy"
    elif age <= 7:
        if citations >= 500:
            return "Landmark"
        if citations >= 100:
            return "Significant"
    else:
        if citations >= 1000:
            return "Foundational"
        if citations >= 500:
            return "Seminal"
    return ""


@dataclass
class Paper:
    title: str
    year: int | None
    venue: str
    citations: int
    doi: str = ""
    pmid: str = ""
    source: str = ""           # which database(s) surfaced it
    is_preprint: bool = False

    @property
    def url(self) -> str:
        if self.doi:
            return f"https://doi.org/{self.doi}"
        if self.pmid:
            return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"
        return ""

    @property
    def tier1(self) -> bool:
        return is_tier1(self.venue)

    @property
    def classification(self) -> str:
        return classify(self.citations, self.year)


@dataclass
class LitResult:
    gene: str
    context: str
    query: str
    total_by_db: dict = field(default_factory=dict)  # {"PubMed": n, "Europe PMC": n}
    databases: list[str] = field(default_factory=list)
    papers: list[Paper] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)   # per-source warnings/errors

    @property
    def total(self) -> int:
        return max(self.total_by_db.values(), default=0)

    @property
    def reported(self) -> bool:
        return any(v > 0 for v in self.total_by_db.values()) or bool(self.papers)

    @property
    def n_high_impact(self) -> int:
        return sum(1 for p in self.papers if p.tier1 or p.classification)


def _norm_title(t: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def _pubmed_term(gene: str, context: str, year_from: int, field: str = "Title/Abstract") -> str:
    """Gene REQUIRED in `field` (relevance gate) + context + year filter.
    field="body" on db=pmc searches the article full text."""
    term = f'"{gene.strip()}"[{field}]'
    if context.strip():
        term += f' AND "{context.strip()}"'
    term += f" AND {year_from}:3000[pdat]"
    return term


def _epmc_query(gene: str, context: str, year_from: int, full_text: bool = False) -> str:
    g = gene.strip()
    # unfielded = title, abstract AND full text of open-access articles
    term = f'"{g}"' if full_text else f'(TITLE:"{g}" OR ABSTRACT:"{g}")'
    if context.strip():
        term += f' AND "{context.strip()}"'
    term += f" AND (PUB_YEAR:[{year_from} TO 3000])"
    return term


def _display_query(gene: str, context: str, year_from: int, full_text: bool = False) -> str:
    ctx = f' AND "{context.strip()}"' if context.strip() else ""
    where = "[Full text]" if full_text else "[Title/Abstract]"
    return f'"{gene.strip()}"{where}{ctx} AND year>={year_from}'


def _search_pubmed(gene: str, context: str, year_from: int,
                   limit: int, timeout: int, db: str = "pubmed") -> tuple[int, list[Paper]]:
    """NCBI PubMed via E-utilities: total hit count + top records.
    (Citation counts aren't available from E-utilities, so those are 0 here and
    get filled in by Europe PMC / OpenAlex during the merge.)"""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    common = {"db": db, "retmode": "json"}
    if NCBI_API_KEY:
        common["api_key"] = NCBI_API_KEY
    term = _pubmed_term(gene, context, year_from, "body" if db == "pmc" else "Title/Abstract")
    r = _get(f"{base}/esearch.fcgi",
             {**common, "term": term, "retmax": limit, "sort": "relevance"}, timeout)
    r.raise_for_status()
    d = r.json().get("esearchresult", {})
    total = int(d.get("count", 0))
    ids = d.get("idlist", [])
    papers: list[Paper] = []
    if ids:
        time.sleep(0.34 if not NCBI_API_KEY else 0.1)
        s = _get(f"{base}/esummary.fcgi", {**common, "id": ",".join(ids)}, timeout)
        s.raise_for_status()
        result = s.json().get("result", {})
        for uid in result.get("uids", []):
            doc = result.get(uid, {})
            ids = {a.get("idtype"): a.get("value", "") for a in doc.get("articleids", [])}
            doi = ids.get("doi", "").lower()
            pmid = uid if db == "pubmed" else ids.get("pmid", "")
            yr = (doc.get("pubdate", "") or "").split(" ")[0]
            papers.append(Paper(
                title=doc.get("title", "(no title)"),
                year=int(yr) if yr.isdigit() else None,
                venue=doc.get("fulljournalname") or doc.get("source", "") or "",
                citations=0, doi=doi, pmid=pmid,
                source="PubMed" if db == "pubmed" else "PubMed Central",
                is_preprint=False))
    return total, papers


def _search_europepmc(gene: str, context: str, year_from: int,
                      limit: int, timeout: int, full_text: bool = False) -> tuple[int, list[Paper]]:
    r = _get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        params={"query": _epmc_query(gene, context, year_from, full_text), "format": "json",
                "pageSize": limit, "sort": "CITED desc", "resultType": "core"},
        timeout=timeout)
    r.raise_for_status()
    data = r.json()
    total = int(data.get("hitCount", 0))
    papers = []
    for it in data.get("resultList", {}).get("result", []):
        journal = (it.get("journalInfo", {}) or {}).get("journal", {}).get("title") \
            or it.get("journalTitle", "") or it.get("bookOrReportDetails", {}).get("publisher", "")
        papers.append(Paper(
            title=it.get("title", "(no title)"),
            year=int(it["pubYear"]) if it.get("pubYear", "").isdigit() else None,
            venue=journal or "",
            citations=int(it.get("citedByCount", 0) or 0),
            doi=(it.get("doi") or "").lower(),
            pmid=it.get("pmid", "") or "",
            source="Europe PMC",
            is_preprint=(it.get("source") == "PPR"),
        ))
    return total, papers


def _search_openalex(gene: str, context: str, year_from: int,
                     limit: int, timeout: int, full_text: bool = False) -> tuple[int, list[Paper]]:
    """OpenAlex works search: gene required in title/abstract (or anywhere in the
    indexed full text), most-cited first."""
    field = "fulltext.search" if full_text else "title_and_abstract.search"
    flt = f'{field}:"{gene.strip()}",from_publication_date:{year_from}-01-01'
    params = {"filter": flt, "sort": "cited_by_count:desc", "per-page": limit,
              "select": "title,publication_year,primary_location,cited_by_count,doi,ids,type"}
    if context.strip():
        params["search"] = context.strip()
    r = _get("https://api.openalex.org/works", params, timeout)
    r.raise_for_status()
    data = r.json()
    papers = []
    for it in data.get("results", []):
        src = ((it.get("primary_location") or {}).get("source") or {})
        papers.append(Paper(
            title=it.get("title") or "(no title)",
            year=it.get("publication_year"),
            venue=src.get("display_name", "") or "",
            citations=int(it.get("cited_by_count", 0) or 0),
            doi=(it.get("doi") or "").lower().removeprefix("https://doi.org/"),
            pmid=((it.get("ids") or {}).get("pmid") or "").rstrip("/").rsplit("/", 1)[-1],
            source="OpenAlex",
            is_preprint=(it.get("type") == "preprint"),
        ))
    return int((data.get("meta") or {}).get("count", 0)), papers


def _merge(papers: list[Paper]) -> list[Paper]:
    """De-duplicate by DOI (then normalized title); keep max citation count and
    record which databases surfaced each paper."""
    by_key: dict[str, Paper] = {}
    for p in papers:
        key = p.doi or _norm_title(p.title)
        if not key:
            continue
        if key in by_key:
            cur = by_key[key]
            cur.citations = max(cur.citations, p.citations)
            srcs = set(cur.source.split(" + ")) | {p.source}
            cur.source = " + ".join(sorted(srcs))
            cur.doi = cur.doi or p.doi
            cur.pmid = cur.pmid or p.pmid
            cur.venue = cur.venue or p.venue
            cur.year = cur.year or p.year
        else:
            by_key[key] = p
    return sorted(by_key.values(), key=lambda x: x.citations, reverse=True)


def search_literature(gene: str, context: str = "", limit: int = 25,
                      year_from: int = 2015, timeout: int = 25,
                      full_text: bool = False) -> LitResult:
    """full_text=True also matches the gene inside article bodies (PubMed Central,
    Europe PMC open-access full text, OpenAlex full-text index): far more hits,
    including papers that only mention the gene in passing."""
    query = _display_query(gene, context, year_from, full_text)
    databases, all_papers, total_by_db, notes = [], [], {}, []

    try:
        cnt, pm = _search_pubmed(gene, context, year_from, limit, timeout)
        all_papers += pm
        databases.append("PubMed")
        total_by_db["PubMed"] = cnt
    except Exception as e:  # noqa: BLE001
        notes.append(f"PubMed unavailable: {type(e).__name__}: {e}")

    if full_text:  # PubMed itself only indexes title/abstract/MeSH; PMC has bodies
        try:
            time.sleep(0.34 if not NCBI_API_KEY else 0.1)
            cnt, pmc = _search_pubmed(gene, context, year_from, limit, timeout, db="pmc")
            all_papers += pmc
            databases.append("PubMed Central")
            total_by_db["PubMed Central"] = cnt
        except Exception as e:  # noqa: BLE001
            notes.append(f"PubMed Central unavailable: {type(e).__name__}: {e}")

    try:
        cnt, epmc = _search_europepmc(gene, context, year_from, limit, timeout, full_text)
        all_papers += epmc
        databases.append("Europe PMC")
        total_by_db["Europe PMC"] = cnt
    except Exception as e:  # noqa: BLE001
        notes.append(f"Europe PMC unavailable: {type(e).__name__}: {e}")

    try:
        cnt, oa = _search_openalex(gene, context, year_from, limit, timeout, full_text)
        all_papers += oa
        databases.append("OpenAlex")
        total_by_db["OpenAlex"] = cnt
    except Exception as e:  # noqa: BLE001
        notes.append(f"OpenAlex unavailable: {type(e).__name__}: {e}")

    # final safety net: enforce year filter on every merged paper
    merged = [p for p in _merge(all_papers) if p.year is None or p.year >= year_from]
    return LitResult(gene=gene, context=context, query=query,
                     total_by_db=total_by_db, databases=databases,
                     papers=merged, notes=notes)


if __name__ == "__main__":
    import sys
    g = sys.argv[1] if len(sys.argv) > 1 else "BDNF"
    ctx = sys.argv[2] if len(sys.argv) > 2 else "Alzheimer"
    res = search_literature(g, ctx)
    print(f"Query: {res.query}  | databases: {res.databases}")
    print(f"Reported: {res.reported}  total(EuropePMC): {res.total}  "
          f"high-impact: {res.n_high_impact}")
    for p in res.papers[:10]:
        tag = "★T1" if p.tier1 else ""
        cls = f"[{p.classification}]" if p.classification else ""
        pp = "(preprint)" if p.is_preprint else ""
        print(f"  {p.citations:>5} cites | {p.year} | {p.venue[:40]:40} {tag}{cls}{pp}")
        print(f"        {p.title[:90]}")
