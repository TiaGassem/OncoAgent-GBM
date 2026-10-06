"""OncoAgent-GBM: Extractive 'chat with your PDF / library' engine.

Credibility-first alternative to SciSpace/Anara-style PDF chat.

Key design decision: this engine is EXTRACTIVE, not generative. It never
writes new sentences about the papers. It locates and returns the author's
OWN sentences, verbatim, each tagged with the source filename and page
number. Because no free text is synthesized, there is nothing to
hallucinate -- every statement is a direct, verifiable quote the user can
check against the original PDF page. This is deliberately more defensible
in a thesis defense than fluent-but-unverifiable AI summaries.

Dependencies: pypdf only (already in requirements.txt). Ranking is a
self-contained TF-IDF so no heavy ML stack is needed.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

try:
    from pypdf import PdfReader
    HAS_PYPDF = True
except ImportError:
    HAS_PYPDF = False


# ------------------------------------------------------------------
# Data structures
# ------------------------------------------------------------------
@dataclass
class Passage:
    """A verbatim passage extracted from a source document."""
    source: str = ""          # filename
    page: int = 0             # 1-based page number
    text: str = ""            # verbatim passage text
    tokens: Counter = field(default_factory=Counter)


@dataclass
class Answer:
    """An extractive answer: ranked verbatim quotes with provenance."""
    question: str = ""
    quotes: list = field(default_factory=list)   # list of dicts
    confidence: int = 0                          # 0-100
    note: str = ""


# ------------------------------------------------------------------
# Tokenisation + light stopword removal
# ------------------------------------------------------------------
_STOP = {
    "the", "a", "an", "and", "or", "of", "to", "in", "on", "for", "is",
    "are", "was", "were", "be", "been", "by", "with", "as", "at", "that",
    "this", "these", "those", "it", "its", "from", "we", "our", "their",
    "which", "has", "have", "had", "not", "can", "may", "also", "than",
    "what", "how", "why", "does", "do", "did", "into", "between", "both",
}


def _tokenize(text: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9]+", text.lower())
            if len(t) > 2 and t not in _STOP]


def _split_sentences(text: str) -> list[str]:
    """Conservative sentence splitter that keeps scientific notation intact."""
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return []
    # Split on sentence end followed by capital/space; avoid splitting e.g. 'et al.'
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", text)
    return [p.strip() for p in parts if len(p.strip()) > 25]


# ------------------------------------------------------------------
# Domain synonym expansion (GBM / CDC25)
# ------------------------------------------------------------------
# Purpose: catch passages that use DIFFERENT wording than the question
# (e.g. query 'cell cycle arrest' should also surface 'G2/M checkpoint').
# This expands the QUERY only -- it never changes the quoted passage, so the
# engine stays fully extractive and verifiable. Expanded (synonym) terms are
# down-weighted so exact author wording always ranks first.
GBM_SYNONYMS: dict[str, list[str]] = {
    "cdc25": ["cdc25a", "cdc25b", "cdc25c", "phosphatase", "dual-specificity"],
    "glioblastoma": ["gbm", "glioma", "astrocytoma", "glial", "neuro-oncology"],
    "overexpressed": ["overexpression", "upregulated", "elevated", "increased", "higher"],
    "prognosis": ["survival", "outcome", "mortality", "malignancy", "aggressiveness"],
    "proliferation": ["growth", "viability", "division", "anti-proliferative", "cytotoxic"],
    "apoptosis": ["cell death", "programmed", "caspase", "pro-apoptotic"],
    "cycle": ["g2/m", "g2", "mitosis", "mitotic", "checkpoint", "arrest", "cyclin", "cdk"],
    "arrest": ["checkpoint", "g2/m", "block", "halt", "cycle"],
    "inhibitor": ["inhibition", "antagonist", "blocker", "ic50", "ki", "potency"],
    "nsc": ["nsc-95397", "nsc95397", "quinone", "naphthoquinone"],
    "mechanism": ["pathway", "mode", "action", "target", "binding"],
    "migration": ["invasion", "motility", "wound", "scratch", "metastasis"],
    "dose": ["concentration", "dose-response", "micromolar", "nanomolar"],
    "expression": ["mrna", "transcript", "protein", "level"],
    "temozolomide": ["tmz", "alkylating", "chemotherapy"],
    "resistance": ["resistant", "refractory", "recurrence", "relapse"],
}

# build reverse lookup so any synonym also points back to its group
_SYN_LOOKUP: dict[str, set] = {}
for _k, _vs in GBM_SYNONYMS.items():
    group = set([_k] + _vs)
    for _w in group:
        _SYN_LOOKUP.setdefault(_w, set()).update(group - {_w})

_SYNONYM_WEIGHT = 0.5  # expanded terms count half as much as exact matches


def _expand_terms(q_terms: Counter) -> dict[str, float]:
    """Return {term: weight}. Exact query terms = full weight; synonyms = half."""
    weighted: dict[str, float] = {}
    for term, tf in q_terms.items():
        weighted[term] = max(weighted.get(term, 0.0), float(tf))
        for syn in _SYN_LOOKUP.get(term, ()):  # tokenised single words only
            for sw in re.split(r"[^a-z0-9]+", syn):
                if len(sw) > 1 and sw not in weighted:
                    weighted[sw] = _SYNONYM_WEIGHT * tf
    return weighted



# ------------------------------------------------------------------
# PDF extraction
# ------------------------------------------------------------------
def extract_pdf_pages(file_bytes: bytes) -> list[tuple[int, str]]:
    """Return [(page_number_1based, page_text), ...] from a PDF byte stream."""
    if not HAS_PYPDF:
        raise RuntimeError("pypdf is not installed. Add 'pypdf>=4.0' to requirements.")
    import io
    reader = PdfReader(io.BytesIO(file_bytes))
    pages = []
    for i, page in enumerate(reader.pages, start=1):
        try:
            txt = page.extract_text() or ""
        except Exception:
            txt = ""
        pages.append((i, txt))
    return pages


# ------------------------------------------------------------------
# Index building
# ------------------------------------------------------------------
def build_passage_index(
    docs: list[tuple[str, bytes]],
    passages_per_page_mode: str = "sentence",
) -> list[Passage]:
    """Build a passage index from one or more (filename, pdf_bytes) docs.

    passages_per_page_mode:
      - 'sentence': each sentence is a passage (precise, best for quoting)
      - 'page': whole page is a passage (coarser)
    """
    index: list[Passage] = []
    for filename, data in docs:
        try:
            pages = extract_pdf_pages(data)
        except Exception:
            continue
        for page_num, page_text in pages:
            if not page_text.strip():
                continue
            if passages_per_page_mode == "page":
                units = [page_text]
            else:
                units = _split_sentences(page_text)
            for unit in units:
                toks = Counter(_tokenize(unit))
                if not toks:
                    continue
                index.append(Passage(source=filename, page=page_num,
                                     text=unit.strip(), tokens=toks))
    return index


def _compute_idf(index: list[Passage]) -> dict[str, float]:
    n = len(index)
    df: Counter = Counter()
    for p in index:
        for term in p.tokens:
            df[term] += 1
    return {t: math.log((n + 1) / (c + 1)) + 1.0 for t, c in df.items()}


# ------------------------------------------------------------------
# Extractive query
# ------------------------------------------------------------------
def answer_extractive(
    index: list[Passage],
    question: str,
    top_k: int = 5,
    idf: Optional[dict] = None,
) -> Answer:
    """Return the top_k verbatim passages most relevant to the question.

    Scoring = sum over query terms of (query_tf * passage_tf * idf). No text
    is generated; quotes are returned exactly as written in the PDF.
    """
    q_terms = Counter(_tokenize(question))
    ans = Answer(question=question)
    if not index:
        ans.note = "No documents indexed. Upload a PDF first."
        return ans
    if not q_terms:
        ans.note = "Query had no searchable terms."
        return ans

    idf = idf or _compute_idf(index)

    # expand query with domain synonyms (query-side only; passages stay verbatim)
    weighted_terms = _expand_terms(q_terms)
    exact_terms = set(q_terms)

    scored = []
    for p in index:
        score = 0.0
        matched = 0
        for term, w in weighted_terms.items():
            if term in p.tokens:
                score += w * p.tokens[term] * idf.get(term, 1.0)
                if term in exact_terms:
                    matched += 1
        if score > 0:
            # length normalisation to avoid favouring very long passages
            norm = score / math.sqrt(sum(p.tokens.values()))
            scored.append((norm, matched, p))

    if not scored:
        ans.note = "No validated passage found in the uploaded document(s) for this query."
        ans.confidence = 0
        return ans

    scored.sort(key=lambda x: x[0], reverse=True)
    top = scored[:top_k]

    # Confidence: coverage of EXACT query terms by the best hit + depth of support
    best_matched = max((m for _, m, _ in top), default=0)
    coverage = best_matched / max(1, len(exact_terms))
    depth = min(1.0, len(top) / top_k)
    ans.confidence = int(round(min(100, 100 * (0.7 * coverage + 0.3 * depth))))

    for norm, matched, p in top:
        ans.quotes.append({
            "source": p.source,
            "page": p.page,
            "quote": p.text,
            "terms_matched": matched,
            "score": round(norm, 4),
        })
    return ans


def format_answer_markdown(ans: Answer) -> str:
    """Render an extractive answer as source-cited Markdown."""
    if not ans.quotes:
        return (f"**Question:** {ans.question}\n\n"
                f"{ans.note or 'No validated evidence found.'}\n\n"
                f"**Credibility score:** {ans.confidence}/100")
    lines = [f"**Question:** {ans.question}", "",
             "**Verbatim evidence from your document(s)** "
             "(nothing is paraphrased or generated):", ""]
    for i, q in enumerate(ans.quotes, start=1):
        lines.append(f"{i}. \u201c{q['quote']}\u201d  \n"
                     f"   \u2014 *{q['source']}, p.{q['page']}* "
                     f"(terms matched: {q['terms_matched']})")
    lines.append("")
    lines.append(f"**Credibility score:** {ans.confidence}/100 "
                 f"(extractive: every line above is a direct quote you can verify on the cited page)")
    return "\n".join(lines)
