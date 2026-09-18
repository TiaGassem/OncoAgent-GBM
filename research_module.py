"""OncoAgent-GBM: PubMed literature research and bibliography module."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Optional
from datetime import datetime

import requests


PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PUBMED_ELINKS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"


@dataclass
class PubMedArticle:
    """Represents a single PubMed article."""
    pmid: str = ""
    title: str = ""
    authors: list[str] = field(default_factory=list)
    journal: str = ""
    pub_date: str = ""
    abstract: str = ""
    doi: str = ""
    keywords: list[str] = field(default_factory=list)
    pmc_id: str = ""
    url: str = ""


def search_pubmed(
    query: str,
    max_results: int = 20,
    sort: str = "relevance",
    use_history: bool = False,
) -> list[PubMedArticle]:
    """Search PubMed and return article metadata."""
    search_params = {
        "db": "pubmed",
        "term": query,
        "retmax": str(max_results),
        "retmode": "xml",
        "sort": sort,
    }
    if use_history:
        search_params["usehistory"] = "y"

    try:
        resp = requests.get(PUBMED_ESEARCH, params=search_params, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        return []

    root = ET.fromstring(resp.content)
    id_list = [id_elem.text for id_elem in root.findall(".//Id") if id_elem.text]

    if not id_list:
        return []

    return fetch_articles_by_id(id_list)


def fetch_articles_by_id(pmid_list: list[str]) -> list[PubMedArticle]:
    """Fetch detailed article information by PMID list."""
    if not pmid_list:
        return []

    params = {
        "db": "pubmed",
        "id": ",".join(pmid_list),
        "retmode": "xml",
        "rettype": "abstract",
    }

    try:
        resp = requests.get(PUBMED_EFETCH, params=params, timeout=30)
        resp.raise_for_status()
    except requests.RequestException:
        return []

    return _parse_pubmed_xml(resp.text)


def _parse_pubmed_xml(xml_text: str) -> list[PubMedArticle]:
    """Parse PubMed efetch XML response into article objects."""
    articles = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return articles

    for article_elem in root.findall(".//PubmedArticle"):
        article = PubMedArticle()

        pmid_elem = article_elem.find(".//PMID")
        if pmid_elem is not None and pmid_elem.text:
            article.pmid = pmid_elem.text.strip()

        title_elem = article_elem.find(".//ArticleTitle")
        if title_elem is not None and title_elem.text:
            article.title = title_elem.text.strip()

        journal_elem = article_elem.find(".//Journal/Title")
        if journal_elem is not None and journal_elem.text:
            article.journal = journal_elem.text.strip()

        year_elem = article_elem.find(".//PubDate/Year")
        month_elem = article_elem.find(".//PubDate/Month")
        day_elem = article_elem.find(".//PubDate/Day")
        if year_elem is not None and year_elem.text:
            parts = [year_elem.text.strip()]
            if month_elem is not None and month_elem.text:
                parts.append(month_elem.text.strip())
            if day_elem is not None and day_elem.text:
                parts.append(day_elem.text.strip())
            article.pub_date = " ".join(parts)

        authors = []
        for author_elem in article_elem.findall(".//Author"):
            last = author_elem.find("LastName")
            first = author_elem.find("ForeName")
            if last is not None and last.text:
                name = last.text.strip()
                if first is not None and first.text:
                    name += f" {first.text.strip()}"
                authors.append(name)
        article.authors = authors

        abstract_parts = []
        for abs_text in article_elem.findall(".//AbstractText"):
            label = abs_text.get("Label", "")
            text = abs_text.text or ""
            text = "".join(abs_text.itertext()).strip()
            if label:
                abstract_parts.append(f"{label}: {text}")
            else:
                abstract_parts.append(text)
        article.abstract = " ".join(abstract_parts)

        article.doi = _extract_doi(article_elem)
        article.pmc_id = _extract_pmc(article_elem)

        article.url = f"https://pubmed.ncbi.nlm.nih.gov/{article.pmid}/"

        articles.append(article)

    return articles


def _extract_doi(article_elem: ET.Element) -> str:
    """Extract DOI from article element."""
    for eid in article_elem.findall(".//ArticleId"):
        if eid.get("IdType") == "doi" and eid.text:
            return eid.text.strip()
    for eid in article_elem.findall(".//ELocationID"):
        if eid.get("EIdType") == "doi" and eid.text:
            return eid.text.strip()
    return ""


def _extract_pmc(article_elem: ET.Element) -> str:
    """Extract PMC ID from article element."""
    for eid in article_elem.findall(".//ArticleId"):
        if eid.get("IdType") == "pmc" and eid.text:
            return eid.text.strip()
    return ""


def format_citation_apa(article: PubMedArticle) -> str:
    """Format an article as APA citation."""
    if not article.authors:
        author_str = "Unknown Author"
    elif len(article.authors) == 1:
        author_str = article.authors[0]
    elif len(article.authors) <= 20:
        author_str = ", ".join(article.authors[:-1]) + f", & {article.authors[-1]}"
    else:
        author_str = ", ".join(article.authors[:19]) + f", ... {article.authors[-1]}"

    year = article.pub_date.split()[0] if article.pub_date else "n.d."

    citation = f"{author_str} ({year}). {article.title}. {article.journal}"
    if article.doi:
        citation += f". https://doi.org/{article.doi}"
    else:
        citation += f". {article.url}"

    return citation


def format_citation_bibtex(article: PubMedArticle) -> str:
    """Format an article as BibTeX entry."""
    if article.authors:
        first_author_last = article.authors[0].split()[-1] if article.authors else "Unknown"
    else:
        first_author_last = "Unknown"

    year = article.pub_date.split()[0] if article.pub_date else "0000"
    key = f"{first_author_last}{year}pmid{article.pmid}"

    authors_bibtex = " and ".join(article.authors) if article.authors else "Unknown"

    entry = (
        f"@article{{{key},\n"
        f"  author = {{{authors_bibtex}}},\n"
        f"  title = {{{article.title}}},\n"
        f"  journal = {{{article.journal}}},\n"
        f"  year = {{{year}}},\n"
        f"  pmid = {{{article.pmid}}},\n"
    )

    if article.doi:
        entry += f"  doi = {{{article.doi}}},\n"
    if article.pmc_id:
        entry += f"  pmcid = {{{article.pmc_id}}},\n"

    entry += f"  url = {{{article.url}}}\n}}"
    return entry


# Pre-built GBM search queries for the research module
GBM_SEARCH_QUERIES = {
    "glioblastoma_general": (
        "(glioblastoma OR GBM) AND drug discovery AND (treatment OR therapy)"
    ),
    "phosphatases_gbm": (
        "(glioblastoma OR GBM) AND (protein tyrosine phosphatase OR PTP OR phosphatase)"
        "AND (therapy OR drug target OR inhibitor)"
    ),
    "bbb_penetration": (
        "blood-brain barrier AND (glioblastoma OR brain tumor) "
        "AND (drug delivery OR CNS penetration)"
    ),
    "ptpn_targets": (
        "(PTPN1 OR PTPN2 OR PTPN11 OR SHP-2 OR PTP1B) AND glioblastoma"
    ),
    "egfrviii_therapy": (
        "(EGFRvIII OR EGFR variant III) AND glioblastoma AND (therapy OR immunotherapy OR CAR-T)"
    ),
    "mgmt_methylation": (
        "(MGMT OR O6-methylguanine-DNA methyltransferase) AND glioblastoma "
        "AND (temozolomide OR methylation OR biomarker)"
    ),
    "idh1_mutation": (
        "(IDH1 OR isocitrate dehydrogenase) AND glioblastoma "
        "AND (mutation OR inhibitor OR 2-hydroxyglutarate)"
    ),
    "dual_specificity_phosphatases": (
        "(DUSP OR dual-specificity phosphatase) AND glioblastoma AND (inhibitor OR ERK)"
    ),
    "molecular_docking_gbm": (
        "(molecular docking OR virtual screening) AND glioblastoma AND (drug design)"
    ),
}


def get_gbm_queries() -> dict[str, str]:
    """Return pre-built GBM literature search queries."""
    return GBM_SEARCH_QUERIES


GBM_REFERENCE_GUIDE = {
    "title": "OncoAgent-GBM Reference Guide: GBM Oncogenes and Phosphatase Targets",
    "version": "1.0",
    "last_updated": datetime.now().isoformat(),
    "sections": {
        "EGFRvIII": {
            "description": (
                "Epidermal Growth Factor Receptor variant III (EGFRvIII) is a constitutively "
                "active mutant lacking exons 2-7, found in approximately 50% of GBM tumors. "
                "It drives proliferation through RAS/MAPK, PI3K/AKT/mTOR, and JAK/STAT pathways."
            ),
            "key_papers": [
                "Huang HS et al. (1997) J Biol Chem. EGFRvIII in gliomas.",
                "Nishikawa R et al. (1994) Nat Med. EGFRvIII is tumor-specific.",
                "Sugawa N et al. (1997) Cancer Res. EGFRvIII in GBM cell lines.",
            ],
            "therapeutic_approaches": [
                "Monoclonal antibodies (cetuximab, nimotuzumab)",
                "Small molecule TKIs (erlotinib, gefitinib, lapatinib)",
                "EGFRvIII-targeted CAR-T cell therapy",
                "Bispecific antibodies (EGFRvIII/CD3)",
                "ADC (antibody-drug conjugates)",
            ],
        },
        "MGMT": {
            "description": (
                "O-6-methylguanine-DNA methyltransferase (MGMT) is a DNA repair enzyme. "
                "Its promoter methylation silences expression, predicting improved response "
                "to temozolomide (TMZ). MGMT methylation status is the most important "
                "biomarker for TMZ sensitivity in GBM."
            ),
            "key_papers": [
                "Esteller M et al. (2000) NEJM. MGMT silencing predicts TMZ response.",
                "Hegi ME et al. (2005) NEJM. MGMT and TMZ benefit.",
                "Stupp R et al. (2005) NEJM. Stupp protocol (RT + TMZ).",
            ],
            "clinical_significance": (
                "MGMT promoter methylation = ~6 month survival benefit with TMZ. "
                "Unmethylated MGMT tumors have limited TMZ benefit; alternative strategies "
                "include MGMT inhibitor O6-benzylguanine (in clinical trials)."
            ),
        },
        "IDH1": {
            "description": (
                "Isocitrate Dehydrogenase 1 (IDH1) R132H mutation produces the "
                "oncometabolite 2-hydroxyglutarate (2-HG), causing epigenetic dysregulation "
                "and HIF-1alpha stabilization. Found in ~8-15% of GBM (higher in secondary GBM)."
            ),
            "key_papers": [
                "Parsons DW et al. (2008) Science. IDH1/2 mutations in glioma.",
                "Yan H et al. (2009) NEJM. IDH1/2 mutations in glioma classification.",
                "Dang L et al. (2009) Nature. IDH mutations produce 2-HG.",
            ],
            "therapeutic_approaches": [
                "Ivosidenib (AG-120) - approved IDH1 inhibitor (AML, being studied in glioma)",
                "Enasidenib (AG-221) - IDH2 inhibitor",
                "Targeting 2-HG downstream epigenetic effects",
                "Combination with hypomethylating agents",
            ],
        },
        "PTPN_phosphatases": {
            "description": (
                "Protein Tyrosine Phosphatases (PTPs) are critical regulators of GBM signaling. "
                "PTP1B (PTPN1) overexpression drives insulin/IGF-1 resistance. SHP-2 (PTPN11) "
                "is an oncogenic driver promoting RAS/MAPK and PI3K/AKT. TC-PTP (PTPN2) "
                "regulates JAK/STAT. DUSP6/MKP-3 is a key ERK regulator."
            ),
            "key_papers": [
                "Bharadwaj G et al. (2014) Cell Signal. SHP-2 in GBM.",
                "Hoekstra E et al. (2016) Oncotarget. PTP1B in cancer.",
                "Rush J et al. (1995) Biochem J. Phosphoproteomics of signaling.",
                "Julien SG et al. (2014) Nat Rev Cancer. PTPs in cancer.",
            ],
            "therapeutic_strategies": [
                "PTP1B inhibitors (MS407, trodusquemine) - preclinical anti-GBM",
                "SHP-099 (allosteric SHP-2 inhibitor) - Phase I trials",
                "SHP2-D1 inhibitors targeting the auto-inhibitory interface",
                "DUSP6 inhibitors synergizing with MEK inhibitors",
                "Substrate-trapping mutants for target validation",
            ],
        },
        "multi_phosphatase_targeting": {
            "description": (
                "Multi-phosphatase targeting aims to simultaneously inhibit multiple "
                "phosphatases to overcome resistance mechanisms in GBM. The catalytic "
                "domains share structural homology but have distinct selectivity pockets."
            ),
            "design_strategies": [
                "Fragment-based drug design targeting conserved catalytic cysteine",
                "Allosteric inhibitors exploiting non-conserved regulatory sites",
                "Proteolysis-targeting chimeras (PROTACs) for phosphatase degradation",
                "Dual-occupancy docking in multi-target pharmacophore models",
            ],
        },
    },
}


def get_gbm_reference_guide() -> dict:
    """Retrieve the full GBM reference guide."""
    return GBM_REFERENCE_GUIDE


def format_multiple_citations(
    articles: list[PubMedArticle], format_type: str = "apa"
) -> list[str]:
    """Format multiple articles in the specified citation format."""
    formatters = {
        "apa": format_citation_apa,
        "bibtex": format_citation_bibtex,
    }
    formatter = formatters.get(format_type, format_citation_apa)
    return [formatter(a) for a in articles]
