"""Provenance-carrying adapters for public neuro-oncology data sources.

Design rules
------------
* No scientific facts are typed into this app. Facts come from the source at query
  time and every record carries a ``Provenance`` (source, endpoint, retrieval time,
  licence note).
* If a source cannot be reached or returns an unexpected shape, the adapter returns
  ``Result(ok=False, error=...)``. It never falls back to remembered values.
* Pure Python + ``requests``; no Streamlit import, so it can be used as a library/CLI
  and unit-tested with mocked HTTP.

STATUS: written against the documented public REST APIs but NOT exercised against the
live services in the build sandbox (no outbound access). Run ``python -m data_sources``
on a networked machine to smoke-test each adapter before relying on it.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Any

import requests

USER_AGENT = "OncoAgent-GBM/2.0 (research-education; open-source)"
TIMEOUT = 20

CT_API = "https://clinicaltrials.gov/api/v2/studies"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
EPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CBIO_API = "https://www.cbioportal.org/api"


@dataclass
class Provenance:
    source: str
    endpoint: str
    retrieved_utc: str
    license_note: str
    query: dict = field(default_factory=dict)


@dataclass
class Result:
    ok: bool
    data: Any = None
    provenance: Provenance | None = None
    error: str = ""
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


_CACHE: dict[tuple, tuple[float, Any]] = {}
CACHE_TTL = 3600


def _get_json(url: str, params: dict | None = None) -> Any:
    key = (url, tuple(sorted((params or {}).items())))
    hit = _CACHE.get(key)
    if hit and time.time() - hit[0] < CACHE_TTL:
        return hit[1]
    last_exc: Exception | None = None
    for attempt in range(2):  # one retry on transient failure
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT,
                             headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            if r.status_code == 429:
                time.sleep(1.5)
                continue
            r.raise_for_status()
            data = r.json()
            _CACHE[key] = (time.time(), data)
            return data
        except (requests.RequestException, ValueError) as exc:
            last_exc = exc
            time.sleep(0.5)
    raise RuntimeError(f"{type(last_exc).__name__}: {last_exc}")


def _fail(source: str, endpoint: str, query: dict, license_note: str, exc: Exception) -> Result:
    return Result(False, None, Provenance(source, endpoint, _now(), license_note, query),
                  error=f"{source} unavailable or returned unexpected data ({exc}). "
                        "No fallback values are shown.")


# ---------------------------------------------------------------- ClinicalTrials.gov
def clinicaltrials_search(condition: str = "glioblastoma", term: str = "", status: str = "",
                          page_size: int = 20) -> Result:
    """List registered studies. Registry listing only: NOT eligibility or trial matching."""
    q = {"query.cond": condition, "pageSize": max(1, min(int(page_size), 50)), "format": "json"}
    if term:
        q["query.term"] = term
    if status:
        q["filter.overallStatus"] = status
    lic = "ClinicalTrials.gov data are public (US NLM); cite the NCT ID."
    try:
        raw = _get_json(CT_API, q)
        rows = []
        for s in raw.get("studies", []):
            p = s.get("protocolSection", {})
            ident = p.get("identificationModule", {})
            rows.append({
                "nct_id": ident.get("nctId"),
                "title": ident.get("briefTitle"),
                "status": p.get("statusModule", {}).get("overallStatus"),
                "phases": ", ".join(p.get("designModule", {}).get("phases", []) or []),
                "interventions": "; ".join(i.get("name", "") for i in
                                           p.get("armsInterventionsModule", {}).get("interventions", []) or []),
            })
        return Result(True, rows, Provenance("ClinicalTrials.gov API v2", CT_API, _now(), lic, q))
    except Exception as exc:  # noqa: BLE001
        return _fail("ClinicalTrials.gov", CT_API, q, lic, exc)


# ---------------------------------------------------------------- ChEMBL
def chembl_search(name: str, limit: int = 5) -> Result:
    """Look up a molecule by name; returns ChEMBL ID, max clinical phase, SMILES."""
    ep = f"{CHEMBL_API}/molecule/search.json"
    q = {"q": name, "limit": max(1, min(int(limit), 20))}
    lic = "ChEMBL is CC BY-SA 3.0 (EMBL-EBI); cite ChEMBL and the release."
    try:
        raw = _get_json(ep, q)
        rows = []
        for m in raw.get("molecules", []):
            st = m.get("molecule_structures") or {}
            rows.append({
                "chembl_id": m.get("molecule_chembl_id"),
                "pref_name": m.get("pref_name"),
                "max_phase": m.get("max_phase"),
                "canonical_smiles": st.get("canonical_smiles"),
            })
        return Result(True, rows, Provenance("ChEMBL REST API", ep, _now(), lic, q))
    except Exception as exc:  # noqa: BLE001
        return _fail("ChEMBL", ep, q, lic, exc)


def chembl_mechanism(chembl_id: str) -> Result:
    ep = f"{CHEMBL_API}/mechanism.json"
    q = {"molecule_chembl_id": chembl_id}
    lic = "ChEMBL is CC BY-SA 3.0 (EMBL-EBI)."
    try:
        raw = _get_json(ep, q)
        rows = [{"mechanism_of_action": m.get("mechanism_of_action"),
                 "action_type": m.get("action_type"),
                 "target_chembl_id": m.get("target_chembl_id")}
                for m in raw.get("mechanisms", [])]
        return Result(True, rows, Provenance("ChEMBL REST API", ep, _now(), lic, q))
    except Exception as exc:  # noqa: BLE001
        return _fail("ChEMBL", ep, q, lic, exc)


# ---------------------------------------------------------------- Europe PMC
def europepmc_search(query: str, page_size: int = 10) -> Result:
    q = {"query": query, "format": "json", "resultType": "lite",
         "pageSize": max(1, min(int(page_size), 50))}
    lic = "Europe PMC metadata; abstracts/full text carry their own licences."
    try:
        raw = _get_json(EPMC_API, q)
        rows = [{"pmid": r.get("pmid"), "doi": r.get("doi"), "title": r.get("title"),
                 "authors": r.get("authorString"), "journal": r.get("journalTitle"),
                 "year": r.get("pubYear"), "source": r.get("source"), "id": r.get("id")}
                for r in raw.get("resultList", {}).get("result", [])]
        return Result(True, rows, Provenance("Europe PMC REST", EPMC_API, _now(), lic, q))
    except Exception as exc:  # noqa: BLE001
        return _fail("Europe PMC", EPMC_API, q, lic, exc)


# ---------------------------------------------------------------- cBioPortal
def cbioportal_studies(keyword: str = "glioma") -> Result:
    """Public cBioPortal studies matching a keyword, with sample counts and PMIDs."""
    ep = f"{CBIO_API}/studies"
    q = {"keyword": keyword, "projection": "SUMMARY"}
    lic = "cBioPortal public studies; cite the original study and cBioPortal."
    try:
        raw = _get_json(ep, q)
        rows = [{"study_id": s.get("studyId"), "name": s.get("name"),
                 "samples": s.get("allSampleCount"), "pmid": s.get("pmid"),
                 "citation": s.get("citation")} for s in raw]
        return Result(True, rows, Provenance("cBioPortal API", ep, _now(), lic, q))
    except Exception as exc:  # noqa: BLE001
        return _fail("cBioPortal", ep, q, lic, exc)


def cbioportal_mutation_fraction(study_id: str, hugo_symbol: str) -> Result:
    """Fraction of *sequenced* samples in a study with >=1 mutation in a gene.

    Association/frequency in one cohort only. Not a risk estimate for any individual.
    """
    lic = "cBioPortal public studies; cite the original study and cBioPortal."
    q = {"study_id": study_id, "gene": hugo_symbol}
    ep = f"{CBIO_API}/molecular-profiles/{study_id}_mutations/mutations"
    try:
        gene = _get_json(f"{CBIO_API}/genes/{hugo_symbol}")
        entrez = gene.get("entrezGeneId")
        if entrez is None:
            raise ValueError("gene not found")
        sl = _get_json(f"{CBIO_API}/sample-lists/{study_id}_sequenced")
        denom = sl.get("sampleCount")
        if not denom:
            raise ValueError("no sequenced sample list for this study")
        muts = _get_json(ep, {"sampleListId": f"{study_id}_sequenced",
                              "entrezGeneId": entrez, "projection": "SUMMARY"})
        mutated = len({m.get("uniqueSampleKey") or m.get("sampleId") for m in muts})
        data = {"study_id": study_id, "gene": hugo_symbol, "mutated_samples": mutated,
                "sequenced_samples": denom, "fraction": mutated / denom}
        warns = ["Single-cohort frequency; cohort composition, assay and filtering vary. "
                 "Not applicable to any individual."]
        return Result(True, data, Provenance("cBioPortal API", ep, _now(), lic, q), warnings=warns)
    except Exception as exc:  # noqa: BLE001
        return _fail("cBioPortal", ep, q, lic, exc)


# ---------------------------------------------------------------- link-outs (no data shown)
def cell_line_links(name: str) -> dict[str, str]:
    """Authoritative places to check a cell line's identity/mutations. Link-outs only."""
    from urllib.parse import quote_plus
    n = quote_plus(name.strip())
    return {
        "Cellosaurus (identity, misidentification, STR)": f"https://www.cellosaurus.org/search?query={n}",
        "DepMap portal (CCLE mutations, dependencies)": "https://depmap.org/portal/",
        "Cell Model Passports": "https://cellmodelpassports.sanger.ac.uk/",
    }


if __name__ == "__main__":  # smoke test on a networked machine
    for label, res in [("ct", clinicaltrials_search(page_size=2)),
                       ("chembl", chembl_search("temozolomide", 1)),
                       ("epmc", europepmc_search("glioblastoma CDC25", 2)),
                       ("cbio", cbioportal_studies("glioma"))]:
        print(label, "OK" if res.ok else f"FAIL: {res.error}")
