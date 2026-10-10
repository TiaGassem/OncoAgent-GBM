"""Streamlit UI for live public-data lookups and the Limitations page."""
from __future__ import annotations

import pandas as pd
import streamlit as st

import data_sources as ds
from runcard import build_runcard, to_json

APP_VERSION = "2.0.0"
DISCLAIMER = "Research and education only. Not a medical device. Do not enter patient data."


@st.cache_data(ttl=3600, show_spinner=False)
def _cached(fn_name: str, **kwargs):
    return getattr(ds, fn_name)(**kwargs).as_dict()


def _provenance_caption(res: dict) -> None:
    p = res.get("provenance")
    if p:
        st.caption(f"Source: {p['source']} | retrieved {p['retrieved_utc']} UTC | {p['license_note']}")
    for w in res.get("warnings") or []:
        st.warning(w)


def _download_card(kind: str, res: dict, params: dict) -> None:
    card = build_runcard(kind, inputs=params, parameters={},
                         results_summary={"ok": res["ok"], "n_rows": len(res["data"]) if isinstance(res.get("data"), list) else None},
                         app_version=APP_VERSION,
                         data_sources=[res["provenance"]] if res.get("provenance") else [])
    st.download_button("Download run card (JSON)", to_json(card), file_name=f"runcard_{kind}.json",
                       mime="application/json", key=f"rc_{kind}")


def tab_data_explorer() -> None:
    st.markdown('<div class="section-header">Live Data Explorer (public sources)</div>', unsafe_allow_html=True)
    st.error(DISCLAIMER + " Results are population-level or registry information, "
             "never advice for an individual.")
    st.caption("Nothing on this page is typed into the app: every table is fetched from the source "
               "when you press the button, with source and retrieval time shown. If a source is "
               "unreachable you get an error, not remembered values.")

    t_mut, t_cell, t_drug, t_trial, t_lit = st.tabs(
        ["Mutations (cBioPortal)", "Cell lines", "Drugs (ChEMBL)", "Trial registry", "Literature"])

    with t_mut:
        kw = st.text_input("Study keyword", "glioma", key="cb_kw")
        if st.button("List public studies", key="cb_go"):
            r = _cached("cbioportal_studies", keyword=kw)
            if r["ok"]:
                st.dataframe(pd.DataFrame(r["data"]), use_container_width=True)
                st.session_state["cb_studies"] = [d["study_id"] for d in r["data"]]
            else:
                st.error(r["error"])
            _provenance_caption(r)
        studies = st.session_state.get("cb_studies", [])
        if studies:
            c1, c2 = st.columns(2)
            sid = c1.selectbox("Study", studies, key="cb_sid")
            gene = c2.text_input("Gene (HGNC symbol)", "IDH1", key="cb_gene").strip().upper()
            if st.button("Mutated fraction in this cohort", key="cb_frac"):
                r = _cached("cbioportal_mutation_fraction", study_id=sid, hugo_symbol=gene)
                if r["ok"]:
                    d = r["data"]
                    st.metric(f"{gene} mutated / sequenced", f"{d['mutated_samples']} / {d['sequenced_samples']}",
                              f"{100 * d['fraction']:.1f}%")
                    _download_card("cbioportal_fraction", r, {"study_id": sid, "gene": gene})
                else:
                    st.error(r["error"])
                _provenance_caption(r)

    with t_cell:
        st.info("Cell-line identity, misidentification and mutation status are NOT stored here. "
                "Check the authoritative sources below before using any line. Authentication "
                "(STR profiling) and provenance problems are well documented for some glioma lines "
                "(e.g. the ATCC U-87 MG origin issue; verify in Cellosaurus).")
        name = st.text_input("Cell line", "U-87 MG", key="cl_name")
        for label, url in ds.cell_line_links(name).items():
            st.markdown(f"- [{label}]({url})")
        st.caption("Cell-line results do not transfer automatically to patients.")

    with t_drug:
        dn = st.text_input("Drug name", "temozolomide", key="dr_name")
        if st.button("Search ChEMBL", key="dr_go"):
            r = _cached("chembl_search", name=dn, limit=5)
            if r["ok"] and r["data"]:
                st.dataframe(pd.DataFrame(r["data"]), use_container_width=True)
                st.session_state["chembl_ids"] = [d["chembl_id"] for d in r["data"] if d.get("chembl_id")]
            elif r["ok"]:
                st.info("No ChEMBL molecule found for that name.")
            else:
                st.error(r["error"])
            _provenance_caption(r)
        ids = st.session_state.get("chembl_ids", [])
        if ids:
            cid = st.selectbox("Show mechanism for", ids, key="dr_cid")
            if st.button("Get mechanism", key="dr_mech"):
                r = _cached("chembl_mechanism", chembl_id=cid)
                if r["ok"] and r["data"]:
                    st.dataframe(pd.DataFrame(r["data"]), use_container_width=True)
                elif r["ok"]:
                    st.info("No mechanism recorded in ChEMBL for this molecule.")
                else:
                    st.error(r["error"])
                _provenance_caption(r)
        st.caption("Brain penetration is not reported here. BBB evidence must come from experimental "
                   "or clinical publications, not from this app.")

    with t_trial:
        c1, c2, c3 = st.columns(3)
        cond = c1.text_input("Condition", "glioblastoma", key="ct_cond")
        term = c2.text_input("Extra term (optional)", "", key="ct_term")
        status = c3.selectbox("Status", ["", "RECRUITING", "ACTIVE_NOT_RECRUITING", "COMPLETED", "TERMINATED"], key="ct_status")
        if st.button("Search registry", key="ct_go"):
            r = _cached("clinicaltrials_search", condition=cond, term=term, status=status, page_size=20)
            if r["ok"]:
                df = pd.DataFrame(r["data"])
                if not df.empty:
                    df["link"] = "https://clinicaltrials.gov/study/" + df["nct_id"].astype(str)
                st.dataframe(df, use_container_width=True,
                             column_config={"link": st.column_config.LinkColumn("Registry page")})
                _download_card("trial_registry", r, {"condition": cond, "term": term, "status": status})
            else:
                st.error(r["error"])
            _provenance_caption(r)
        st.warning("This is a registry listing only. The app does not assess eligibility or match "
                   "anyone to a trial. Eligibility is decided by trial investigators; consult a "
                   "qualified clinician.")

    with t_lit:
        q = st.text_input("Europe PMC query", "glioblastoma AND CDC25", key="ep_q")
        if st.button("Search Europe PMC", key="ep_go"):
            r = _cached("europepmc_search", query=q, page_size=15)
            if r["ok"]:
                st.dataframe(pd.DataFrame(r["data"]), use_container_width=True)
            else:
                st.error(r["error"])
            _provenance_caption(r)


def tab_limitations() -> None:
    st.markdown('<div class="section-header">Limitations & Honesty</div>', unsafe_allow_html=True)
    st.error(DISCLAIMER)
    st.markdown("""
**What this tool is:** an open research and teaching aid for computational drug-discovery workflows in
glioblastoma / neuro-oncology.

**What it is not:** a diagnostic, prognostic, treatment-selection or trial-matching tool. It makes no
statement about any individual.

**Computation**
- Docking (AutoDock Vina) is **non-covalent**, rigid-receptor, simplified preparation. Scores are
  relative hypotheses, not binding affinities or efficacy. CDC25 inhibitors are often covalent/redox-active
  quinones; Vina cannot model that chemistry.
- Validate any setup by redocking a co-crystallised ligand (heavy-atom RMSD, target <= 2.0 Å) before trusting
  a new score.
- BBB and ADMET outputs are rule-based heuristics, not measurements.
- Redox-cycling quinones and other PAINS-type compounds can look active for artefactual reasons.

**Data**
- Frequencies, studies, trials and drug records are fetched live and carry source + retrieval time.
- Cohort associations (TCGA, CGGA, cBioPortal) are not causal and not predictive for individuals.
- Cell-line findings do not transfer automatically to patients; check line identity in Cellosaurus.

**Known gaps (not yet implemented)**: CGGA/TCGA survival analysis with multiple-testing correction, DepMap
dependency, GDSC/PRISM sensitivity, Cellosaurus API integration, Arabic right-to-left layout, formal expert review.

**Privacy**: do not enter patient data or personal genomic files. Lab-notebook entries stay in your session.

**Errors**: please report them via the GitHub issue tracker; fixes are logged in CHANGELOG.md.
""")
