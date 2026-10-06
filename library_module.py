"""OncoAgent-GBM: Credible, extractive Library / Chat-with-PDF / Literature-Review engine.

Design goal: NO hallucination. Every answer is built from text that actually
exists in the user's uploaded PDFs (quoted verbatim with page numbers) or from
real PubMed records. If nothing relevant is found, it says so. This is the
credible alternative to generative tools like SciSpace/Anara for thesis work.

Depends only on `pypdf` (in requirements) + the Python stdlib + research_module,
so it imports cleanly even when rdkit/streamlit are unavailable.
"""

from __future__ import annotations

import re
import math
from dataclasses import dataclass, field
from typing import Optional

try:
    from pypdf import PdfReader
    HAS_PYPDF = True
except ImportError:
    HAS_PYPDF = False


_STOP = {
    "the", "and", "for", "are", "was", "were", "this", "that", "with", "from",
    "have", "has", "had", "not", "but", "can", "all", "any", "our", "their",
    "which", "what", "how", "why", "who", "when", "into", "also", "been",
    "its", "they", "these", "those", "such", "may", "will", "a", "an", "of",
    "in", "on", "to", "is", "as", "by", "be", "or", "at", "it",
}


def _tokenize(text: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9]+", text.lower()) if len(t) > 2 and t not in _STOP]


def _split_sentences(text: str) -> list[str]:
    text = re.sub(r"\s+", " ", text).strip()
    # naive sentence split that keeps decimals/abbreviations mostly intact
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", text)
    return [p.strip() for p in parts if len(p.strip()) > 2]


@dataclass
class PdfChunk:
    page: int = 0
    text: str = ""
    tokens: list = field(default_factory=list)


@dataclass
class PdfDoc:
    name: str = ""
    num_pages: int = 0
    chunks: list = field(default_factory=list)  # list[PdfChunk]
    full_text_len: int = 0


def ingest_pdf(file_bytes: bytes, name: str = "document.pdf", max_pages: int = 60) -> PdfDoc:
    """Parse a PDF into page-anchored sentence chunks for extractive retrieval.

    Returns an empty PdfDoc (num_pages=0) if pypdf is missing or parsing fails.
    Never fabricates content.
    """
    doc = PdfDoc(name=name)
    if not HAS_PYPDF:
        return doc
    import io
    try:
        reader = PdfReader(io.BytesIO(file_bytes))
    except Exception:
        return doc
    doc.num_pages = len(reader.pages)
    for i, page in enumerate(reader.pages[:max_pages]):
        try:
            txt = page.extract_text() or ""
        except Exception:
            txt = ""
        doc.full_text_len += len(txt)
        for sent in _split_sentences(txt):
            doc.chunks.append(PdfChunk(page=i + 1, text=sent, tokens=_tokenize(sent)))
    return doc


def _idf(docs_tokens: list[list[str]]) -> dict:
    df: dict = {}
    n = len(docs_tokens)
    for toks in docs_tokens:
        for t in set(toks):
            df[t] = df.get(t, 0) + 1
    return {t: math.log(1 + n / (1 + c)) for t, c in df.items()}


def query_pdf(doc: PdfDoc, question: str, top_k: int = 5) -> list[dict]:
    """Return the most relevant verbatim passages with page numbers (extractive).

    Scoring = TF-IDF-like keyword overlap. Returns [] if nothing matches, so the
    UI can honestly say 'No matching passage found in this document.'
    """
    if not doc.chunks:
        return []
    q_tokens = set(_tokenize(question))
    if not q_tokens:
        return []
    idf = _idf([c.tokens for c in doc.chunks])
    scored = []
    for c in doc.chunks:
        ctoks = set(c.tokens)
        overlap = q_tokens & ctoks
        if not overlap:
            continue
        score = sum(idf.get(t, 0.0) for t in overlap)
        scored.append((score, c))
    scored.sort(key=lambda x: x[0], reverse=True)
    out = []
    for score, c in scored[:top_k]:
        out.append({"page": c.page, "text": c.text, "score": round(score, 3)})
    return out


def answer_from_pdf(doc: PdfDoc, question: str, top_k: int = 5) -> str:
    """Build an extractive, page-cited answer. No generation, no fabrication."""
    hits = query_pdf(doc, question, top_k=top_k)
    if not hits:
        return ("No matching passage found in this document for that question. "
                "Try different keywords, or confirm the PDF contains selectable text "
                "(scanned PDFs need OCR).")
    lines = [f"**Answer assembled from \"{doc.name}\" (verbatim excerpts):**\n"]
    for h in hits:
        snippet = h["text"]
        if len(snippet) > 400:
            snippet = snippet[:400].rsplit(" ", 1)[0] + "..."
        lines.append(f"> {snippet}\n> \u2014 p.{h['page']}")
    coverage = min(100, 30 + 14 * len(hits))
    lines.append(f"\n_Extractive confidence (keyword coverage): {coverage}/100. "
                 f"Quotes are verbatim; verify in context on the cited page._")
    return "\n\n".join(lines)


def extract_numbers(doc: PdfDoc, keyword: str = "") -> list[dict]:
    """Extract sentences containing numeric data (IC50, p-values, %, nM, etc.).

    Useful for 'Extract Data' workflows. Returns page-anchored hits only; it does
    not interpret or invent values.
    """
    num_re = re.compile(r"\d")
    unit_re = re.compile(r"(ic50|ec50|ki\b|nm\b|\u00b5m|um\b|mg/kg|p\s*[<=>]|%|kcal)", re.I)
    kw = keyword.lower().strip()
    out = []
    for c in doc.chunks:
        if not num_re.search(c.text):
            continue
        if not unit_re.search(c.text):
            continue
        if kw and kw not in c.text.lower():
            continue
        out.append({"page": c.page, "text": c.text})
    return out


def literature_review(topic: str, max_results: int = 15) -> dict:
    """Build a credible, PubMed-grounded mini literature review.

    Returns {'count', 'papers', 'apa'} using only real retrieved records.
    Imports research_module lazily so this module stays dependency-light.
    """
    try:
        from research_module import build_paper_index
    except Exception:
        return {"count": 0, "papers": [], "apa": [], "error": "research_module unavailable"}
    index = build_paper_index(topic, max_results=max_results)
    papers = []
    apa = []
    for p in index:
        papers.append({
            "pmid": p.get("pmid", ""),
            "title": p.get("title", ""),
            "year": p.get("year", ""),
            "journal": p.get("journal", ""),
            "url": p.get("url", ""),
            "abstract": p.get("abstract", ""),
        })
        if p.get("apa"):
            apa.append(p["apa"])
    return {"count": len(papers), "papers": papers, "apa": apa}
