"""OncoAgent-GBM: Glioblastoma Drug Discovery Platform -- Clinical Research Dashboard."""

from __future__ import annotations

import os
import sys
import io
import json
import tempfile
from datetime import datetime

import streamlit as st
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cheminformatics import (
    screen_compound, batch_screen,
    compute_all_descriptors, CompoundMetrics, export_results_csv_rows,
    parse_smiles, generate_3d_conformer, get_mol_block_3d,
)
from docking_engine import (
    GridBox, DockingResult, fetch_pdb_from_rcsb, clean_pdb,
    extract_ligand_from_pdb, compute_grid_from_ligand,
    compute_grid_from_residues, smiles_to_pdbqt, pdb_to_pdbqt,
    run_vina_docking, PHOSPHATASE_TARGETS, get_phosphatase_info,
    list_phosphatase_targets,
)
from research_module import (
    search_pubmed, format_citation_apa, format_citation_bibtex,
    get_gbm_queries, get_gbm_reference_guide, format_multiple_citations,
)
from anonymizer import (
    MedicalNoteAnonymizer, anonymize_medical_note, detect_pii_in_text,
    generate_sample_medical_note, GBM_MEDICAL_NOTE_TEMPLATE,
)
from agent_evaluator import (
    evaluate_lead, generate_pdf_report, generate_docx_report,
    PatientProfile, estimate_toxicity, compare_with_protox, ToxicityProfile, LeadEvaluation,
)

st.set_page_config(
    page_title="OncoAgent-GBM",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

LANG = {
    "en": {
        "app_title": "OncoAgent-GBM",
        "app_subtitle": "Glioblastoma Multiforme Drug Discovery and Clinical Decision Support Platform",
        "tab1": "Compound Screening",
        "tab2": "Molecular Docking",
        "tab3": "Literature & Bibliography",
        "tab4": "Patient Data & Trial Matching",
        "execute": "Execute Screening",
        "export_csv": "Export CSV",
        "export_pdf": "Export PDF Report",
        "export_docx": "Export Word Report",
        "no_results": "No results to export.",
        "select_language": "Language",
    },
    "fr": {
        "app_title": "OncoAgent-GBM",
        "app_subtitle": "Plateforme de decouverte de médicaments et d'aide a la decision clinique pour le glioblastome",
        "tab1": "Depistage de composés",
        "tab2": "Docking moléculaire",
        "tab3": "Littérature et bibliographie",
        "tab4": "Données patient et essais cliniques",
        "execute": "Lancer le depistage",
        "export_csv": "Exporter CSV",
        "export_pdf": "Exporter rapport PDF",
        "export_docx": "Exporter rapport Word",
        "no_results": "Aucun résultat a exporter.",
        "select_language": "Langue",
    },
    "ar": {
        "app_title": "OncoAgent-GBM",
        "app_subtitle": "منصة اكتشاف ادوية الورم النجمي الشبكي المتعدد ودعم القرارات السريرية",
        "tab1": "فحص المركبات",
        "tab2": "الترابط الجزيئي",
        "tab3": "الأدبيات والمراجع",
        "tab4": "بيانات المرضى والتجارب السريرية",
        "execute": "بدء الفحص",
        "export_csv": "تصدير CSV",
        "export_pdf": "تصدير تقرير PDF",
        "export_docx": "تصدير تقرير Word",
        "no_results": "لا توجد نتائج للتصدير.",
        "select_language": "اللغة",
    },
}

with st.sidebar:
    lang_choice = st.selectbox("Language / Langue / اللغة", ["en", "fr", "ar"], format_func=lambda x: {"en": "English", "fr": "Francais", "ar": "العربية"}[x])
T = LANG[lang_choice]

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

:root {
    --primary: #1a365d;
    --primary-light: #2a4a7f;
    --accent: #2563eb;
    --accent-hover: #1d4ed8;
    --success: #059669;
    --success-bg: #ecfdf5;
    --warning: #d97706;
    --warning-bg: #fffbeb;
    --danger: #dc2626;
    --danger-bg: #fef2f2;
    --bg: #f8fafc;
    --surface: #ffffff;
    --border: #e2e8f0;
    --text: #1e293b;
    --text-secondary: #64748b;
    --text-muted: #94a3b8;
}

.stApp {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    background-color: var(--bg);
}

header[data-testid="stHeader"] {
    background-color: var(--primary);
    padding: 0.5rem 1rem;
}

header[data-testid="stHeader"] * {
    color: white !important;
}

.section-header {
    font-size: 0.95rem;
    font-weight: 600;
    color: var(--primary);
    border-bottom: 2px solid var(--accent);
    padding-bottom: 0.5rem;
    margin-bottom: 1.25rem;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    font-family: 'Inter', sans-serif;
}

.verdict-box {
    padding: 0.75rem 1.25rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.85rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    border-left: 4px solid;
    font-family: 'Inter', sans-serif;
}

.verdict-viable {
    background-color: var(--success-bg);
    color: var(--success);
    border-left-color: var(--success);
}

.verdict-rejected {
    background-color: var(--danger-bg);
    color: var(--danger);
    border-left-color: var(--danger);
}

.verdict-conditional {
    background-color: var(--warning-bg);
    color: var(--warning);
    border-left-color: var(--warning);
}

.stTabs [data-baseweb="tab-list"] {
    gap: 2px;
    background: var(--border);
    border-radius: 4px;
    padding: 3px;
}

.stTabs [data-baseweb="tab"] {
    border-radius: 3px;
    font-weight: 500;
    font-size: 0.8rem;
    padding: 8px 16px;
    letter-spacing: 0.02em;
}

.stTabs [aria-selected="true"] {
    background: var(--surface) !important;
    border-bottom: none !important;
}

div[data-testid="stMetric"] {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 0.5rem 0.75rem;
}

div[data-testid="stMetric"] label {
    font-size: 0.65rem !important;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    color: var(--text-secondary) !important;
}

div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
    font-size: 1.05rem !important;
    font-weight: 600 !important;
    color: var(--text) !important;
}

.sidebar .sidebar-content {
    background-color: var(--primary);
}

hr {
    border: none;
    border-top: 1px solid var(--border);
    margin: 1.5rem 0;
}
</style>
""", unsafe_allow_html=True)

GBM_DRUGS = {
    "[NSC-8583] Temozolomide (TMZ)": "CN1N=NC2=C(N=CN2C1=O)C(N)=O",
    "[NSC-79037] Lomustine (CCNU)": "ClCCN(N=O)C(=O)NC1CCCCC1",
    "[NSC-409962] Carmustine (BCNU)": "ClCCN(C(=O)N(CCCl)C(=O)N)N=O",
    "[NSC-118233] Procarbazine": "CC(C)NC(=O)C1=CC=CC=C1NN",
    "[NSC-67574] Vincristine": "CO[C@H]1C[C@H](C2=C1C(=O)OC3=C2C(=O)C4=C3OCO4)N(C)C[C@@H]5OC(=O)[C@@]6(C7=C5C=CC(=C7)OC)C(=O)OC6C",
    "[NSC-123127] Dacarbazine": "CN(C)/N=N/c1ncc[nH]c1=O",
    "[NSC-359078] Nimustine (ACNU)": "O=C(NCCCl)N(N=O)C1CCCCC1",
    "[NSC-172112] Etoposide": "COC1=CC(=CC(=C1O)[C@@H]2C3=C(C=CC(=C3)OCO2)C4=C5[C@@H]([C@@H](OC5=O)C6=CC=C(C=C6)OC)OC(=O)[C@@H]4O)OC",
    "[NSC-249992] Irinotecan": "C1CCC2=C1C3=CC=C4C(=C3C(=O)N2CC5=CC=C(C=C5)OC(=O)NC6CCN(C)CC6)OCO4",
    "[NSC-141540] Topotecan": "O=C1C2=C(OCO2)C(=O)c2cc(O)ccc21",
    "[NSC-718781] Erlotinib (EGFR TKI)": "COC1=CC2=C(C=CN=C2C=C1OCCOC)NC3=CC=CC=C3",
    "[NSC-641538] Gefitinib (EGFR TKI)": "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4",
    "[NSC-727957] Sunitinib (VEGFR/PDGFR)": "CCN(CC)CCNC(=O)C1=C(C(=C(/C1=C\\2/C3=C(C=CC=C3)NC2=O)C)C)C",
    "[NSC-737664] Lapatinib (EGFR/HER2)": "CS(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)Cl)Cl",
    "[NSC-654649] Sorafenib (multi-kinase)": "CNC(=O)C1=CC=CC=C1OC2=C(N=CC=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F",
    "[NSC-747854] Pazopanib (VEGFR)": "CC1=CC(=CC=C1)NC2=NC=CC=C2NC(=O)CN3CCN(C)CC3",
    "[NSC-737664] Dasatinib (Src/ABL)": "CC1=NC(=CC(=N1)NC2=CC(=CC=C2)C(=O)NC3=CC=C(C=C3)OC4=CC=CC=C4C(F)(F)F)NC5=CC=CC=C5",
    "[NSC-718781] Vandetanib (VEGFR/EGFR)": "COc1ccc2ncnc(Nc3ccc(Br)cc3)c2c1",
    "[NSC-757487] Bosutinib (SRC inhibitor)": "COC1=C(C=CC(=C1Cl)NC2=NC=CC(=N2)C3=CC(=CC=C3)OC4=CC=CC=C4)OC",
    "[NSC-725741] Axitinib (VEGFR)": "CS(=O)CCNC(=O)C1=CC=C(C=C1)NC2=NC3=C(C=CC=C3)C(=N2)C4=CC=C(C=C4)C",
    "[NSC-747175] ENMD-2076 (Aurora/FGFR)": "COC1=CC=C(C=C1)NC2=NC=NC3=C2C=CC=C3NC(=O)NC4=CC=CC=C4C(F)(F)F",
    "[NSC-95382] PTP1B inhibitor": "CC(=O)NC1=CC=C(C=C1)S(=O)(=O)NC(=O)NC2=CC=CC=C2",
    "[NSC-663248] SHP-2 inhibitor (SHP099)": "CC1=CC=C(C=C1)S(=O)(=O)NC2=NC3=C(C=CC=C3)C(=N2)OC",
    "[NSC-104890] Staurosporine (PKC inhibitor)": "CC1=CC2=C(C=C1C)NC3=C2C(=O)NC4=CC=CC=C43",
    "[NSC-139407] AG1478 (EGFR inhibitor)": "CC1=CC(=CC=C1)NC2=NC=NC3=C2C=CC=C3NC4=CC=C(C=C4)Br",
    "[NSC-164536] Bortezomib (proteasome)": "CC(C1=CC=C(C=C1)NC(=O)C2=CC=CC=C2)NC(=O)C(CCC(=O)O)NC(=O)C(CC3=CC=C(C=C3)O)NC(=O)C(CC4=CC=C(C=C4)N)NC(=O)C(CC5=CC=CC=C5)NC(=O)C(CC6=CC=C(C=C6)O)NC(=O)C(C(C)O)NC(=O)C(CC7=CC=CC=C7)NC(=O)C(CC8=CC=C(C=C8)O)NC(=O)C(CC9=CC=C(C=C9)O)",
    "[NSC-271603] ABT-888 (PARP inhibitor)": "C1CC1C(=O)NC2=CC=C(C=C2)C3=NN4C(=N3)C=CC=N4",
    "[NSC-613327] Curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
    "[NSC-152002] Thalidomide": "O=C1C(=O)N(c2ccccc2)C(=O)N1C1CCCCC1",
    "[NS-733438] Rapamycin (Sirolimus)": "CCC1=C[C@H]([C@H](O)[C@H](C)C=CC=C[C@@H](O)[C@H](C)C[C@H](O))OC(=O)[C@H](O)[C@H](C)C=CC=C[C@@H](O)[C@H](C)C[C@H](O)CC(=O)O[C@@H]1C",
    "[NSC-727957] Temsirolimus": "CC(C)[C@@H](O)C=C[C@@H](O)[C@H](C)OC(=O)C=C[C@@H]1OC(=O)C[C@@H](O)C[C@@H](C)[C@@H](O)C=C[C@@H](O)[C@H](C)OC(=O)C=CC1=O",
    "[NSC-730868] Everolimus": "CCC(=O)OC1CC(CCC1C)OC2CC(OC(C2)C3CC(C(=O)O3)O)OC4CC(OC(C4)C5CC(C(=O)O5)O)OC",
    "[NSC-747954] BEZ235/Dactolisib": "O=C1NC2=C(N1)C=CC(=C2)C3=NC4=CC=CC=C4N3C5=CC=CC=C5",
    "[NSC-675423] Vorinostat (SAHA)": "O=C(/C=C/C1=CC=CC=C1)NCCCCCCCC(=O)NO",
    "[NSC-687582] Romidepsin": "OC1=CC(OC(=O)C(CCC2=CC=CC=C2)NC(=O)C3=CC=CC=C3)C4=C1C(=O)NCCCC4",
    "[NSC-724017] Belinostat": "ONS(=O)(=O)C1=CC=C(C=C1)NC(=O)OC2=CC=CC=C2",
    "[NSC-277097] O6-Benzylguanine (MGMT inhibitor)": "Nc1nc2ncn(Cc3ccccc3)c2c(=O)[nH]1",
    "[NSC-42066] O6-Methylguanine": "O=c1nc(N)nc2ncn(C)n12",
    "[NSC-8806] Busulfan": "CS(=O)(=O)OCCCOS(=O)(=O)C",
    "[NSC-342790] Bendamustine": "C1C2CN(C1C(=O)N=C(N2)N)C3=CC=C(C=C3)C(=O)NCCl",
    "[NSC-724958] Chlorambucil": "C1=CC=C(C=C1)CCC(=O)CCl",
    "[NSC-714537] NEO-212 (TMZ-POH)": "CCC(=O)NC1=CC=CC(=C1)CC2=CC=C(C=C2)NC(=O)C3=NN=C4C(=O)N(C)C(=N4)N3C",
    # === CDC25 / DUAL-SPECIFICITY PHOSPHATASE INHIBITORS (U251/U87 GBM) ===
    "[NSC-95397] Cdc25/MKP inhibitor": "OCCSC1=C(SCCO)C(=O)C2=CC=CC=C2C1=O",
    "[NSC-663284] Cdc25 inhibitor (reference)": "O=C1C=CC(=O)C(Nc2ccc(N(CCOCc3ccccc3)C(=O)c3ccccc3)cc2)=C1",
    "[NSC-668394] Naphthoquinone Cdc25 inhibitor": "NC1=C(N)C(=O)c2ccccc2C1=O",
    "[Monohydroxy-NSC95397] M-NSC (Cdc25A)": "OCCSC1=C(SCCO)C(=O)c2cc(O)ccc21",
    "[Dihydroxy-NSC95397] D-NSC (most potent Cdc25)": "OCCSC1=C(SCCO)C(=O)c2c(O)ccc(O)c21",
    "[Cpd5] Naphthoquinone (Cdc25 ligand)": "O=C1C=CC(=O)C(c2ccccc2)=C1",
}


def render_header():
    st.markdown(
        '<div style="background:#1a365d;padding:2rem 2.5rem;border-radius:4px;margin-bottom:2rem;'
        'border-left:6px solid #2563eb;">'
        '<div style="display:flex;justify-content:space-between;align-items:flex-start;">'
        '<div>'
        '<h1 style="color:#fff;margin:0;font-size:1.6rem;font-weight:700;letter-spacing:0.03em;'
        'font-family:Inter,sans-serif;">'
        f'{T["app_title"]}</h1>'
        '<p style="color:#93c5fd;margin:0.3rem 0 0 0;font-size:0.82rem;font-weight:400;'
        'letter-spacing:0.02em;">'
        f'{T["app_subtitle"]}</p>'
        '<p style="color:#bfdbfe;margin:0.6rem 0 0 0;font-size:0.72rem;font-weight:300;'
        'letter-spacing:0.01em;">'
        'Compound Screening | Molecular Docking | Literature Research | '
        'Cell Line Database | Clinical Trial Matching | Treatment Planning</p>'
        '</div>'
        '<div style="text-align:right;">'
        '<p style="color:#93c5fd;margin:0;font-size:0.65rem;font-weight:400;'
        'letter-spacing:0.02em;">Version 1.0</p>'
        '<p style="color:#93c5fd;margin:0.15rem 0 0 0;font-size:0.65rem;font-weight:400;">'
        'For Research Use Only</p>'
        '</div>'
        '</div>'
        '</div>',
        unsafe_allow_html=True,
    )


# ============================================================
# TAB 1: Compound Screening & BBB Triage
# ============================================================
def tab_compound_screening():
    st.markdown('<div class="section-header">Compound Screening & Blood-Brain Barrier Triage</div>', unsafe_allow_html=True)

    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("**Input**")
        smiles_input = st.text_area(
            "SMILES String(s) -- one per line",
            value="CN1N=NC2=C(N=CN2C1=O)C(N)=O",
            height=90,
        )

        drug_keys = list(GBM_DRUGS.keys())
        selected = st.selectbox("Select from GBM Drug Library", ["Manual Input"] + drug_keys)
        if selected != "Manual Input":
            smiles_input = GBM_DRUGS[selected]

    with col2:
        st.markdown("**Screening Parameters**")
        st.markdown(
            "| Criterion | Threshold | Score Weight |\n"
            "|---|---|---|\n"
            "| Lipinski Rule of 5 | MW < 500, LogP < 5, HBD <= 5, HBA <= 10 | 30% |\n"
            "| BBB Penetration | MW < 400, LogP 1.5-3.5, TPSA < 90 | 25% |\n"
            "| Veber Oral Bioavailability | TPSA <= 140, RotBonds <= 10 | 10% |\n"
            "| Structural Quality | Fraction Csp3, Ring Count, Heavy Atoms | 35% |"
        )

    if st.button("Execute Screening", type="primary", use_container_width=True):
        if not smiles_input.strip():
            st.error("Enter at least one valid SMILES string.")
            return

        smiles_list = [s.strip() for s in smiles_input.strip().split("\n") if s.strip()]

        with st.spinner("Computing physicochemical descriptors and BBB scores..."):
            results = batch_screen(smiles_list)

        for i, result in enumerate(results):
            st.markdown("---")
            col_a, col_b = st.columns([1, 2])

            with col_a:
                if result.valid:
                    st.markdown(f"**SMILES:** `{result.canonical_smiles}`")

                if result.overall_pass:
                    st.markdown('<div class="verdict-box verdict-viable">PASS -- VIABLE LEAD</div>', unsafe_allow_html=True)
                elif result.lipinski_pass and result.bbb_pass:
                    st.markdown('<div class="verdict-box verdict-conditional">CONDITIONAL -- Review Required</div>', unsafe_allow_html=True)
                else:
                    st.markdown('<div class="verdict-box verdict-rejected">FAIL -- Below Thresholds</div>', unsafe_allow_html=True)

            with col_b:
                st.markdown("**Physicochemical Properties**")
                if not result.valid:
                    st.error(f"Invalid SMILES: {result.error}")
                    continue

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("MW (Da)", f"{result.molecular_weight:.1f}")
                m2.metric("LogP", f"{result.logp:.2f}")
                m3.metric("TPSA", f"{result.tpsa:.1f} A2")
                m4.metric("HBD / HBA", f"{result.hbd} / {result.hba}")

                m5, m6, m7, m8 = st.columns(4)
                m5.metric("Rot. Bonds", result.rotatable_bonds)
                m6.metric("Rings", result.num_rings)
                m7.metric("Fsp3", f"{result.fraction_csp3:.2f}")
                m8.metric("BBB Score", f"{result.bbb_score:.2f}")

                st.markdown("**Compliance**")
                c1, c2, c3 = st.columns(3)
                lip_s = "PASS" if result.lipinski_pass else "FAIL"
                bbb_s = "PASS" if result.bbb_pass else "FAIL"
                veb_s = "PASS" if result.veber_pass else "FAIL"
                lip_c = "var(--success)" if result.lipinski_pass else "var(--danger)"
                bbb_c = "var(--success)" if result.bbb_pass else "var(--danger)"
                veb_c = "var(--success)" if result.veber_pass else "var(--danger)"
                c1.markdown(f'<span style="color:{lip_c};font-weight:600;">Lipinski: {lip_s}</span> <span style="color:var(--text-muted);font-size:0.8rem;">({result.lipinski_violations} violations)</span>', unsafe_allow_html=True)
                c2.markdown(f'<span style="color:{bbb_c};font-weight:600;">BBB: {bbb_s}</span> <span style="color:var(--text-muted);font-size:0.8rem;">({result.bbb_score:.2f})</span>', unsafe_allow_html=True)
                c3.markdown(f'<span style="color:{veb_c};font-weight:600;">Veber: {veb_s}</span>', unsafe_allow_html=True)

        if results:
            st.markdown("---")
            st.markdown('<div class="section-header">Toxicity Pre-Screen (Rule-Based Heuristic)</div>', unsafe_allow_html=True)
            st.warning(
                "This module is a transparent **rule-based heuristic** (structural alerts + "
                "physicochemical thresholds). It is **NOT ProTox-3**, which is a machine-learning "
                "service trained on ~40,000 compounds with no public API. Treat these values as a "
                "fast pre-screen; expect divergence from ProTox-3 and always confirm with "
                "experimental toxicology."
            )

            tox_profiles = {}
            for i, result in enumerate(results):
                if result.valid:
                    mol = parse_smiles(result.smiles)
                    if mol:
                        tox_profile = estimate_toxicity(mol)
                        tox_profiles[i] = tox_profile
                        st.markdown(f"**Compound {i+1}: {result.canonical_smiles[:50]}**")
                        t1, t2, t3, t4, t5, t6, t7 = st.columns(7)
                        t1.metric("GHS Class", tox_profile.toxicity_class.split("(")[0].strip() if tox_profile.toxicity_class else "N/A")
                        t2.metric("LD50 Est.", f"{tox_profile.ld50_estimate:.0f} mg/kg")
                        t3.metric("Ames", tox_profile.ames_mutagenicity)
                        t4.metric("hERG", tox_profile.herg_inhibition)
                        t5.metric("Liver", tox_profile.hepatotoxicity)
                        t6.metric("Alert Score", tox_profile.tox_score)
                        t7.metric("Overall Risk", tox_profile.overall_tox_risk)

                        with st.expander(f"Transparency detail (Compound {i+1})"):
                            st.markdown(f"**Method:** {tox_profile.method}")
                            st.markdown(f"**Confidence:** {tox_profile.confidence}")
                            st.markdown("**Descriptors driving the score:**")
                            st.dataframe(pd.DataFrame([tox_profile.descriptors]), use_container_width=True)
                            st.markdown("**Structural alerts / rules triggered:**")
                            for alert in tox_profile.matched_alerts:
                                st.markdown(f"- {alert}")

            if tox_profiles:
                st.markdown("---")
                st.markdown('<div class="section-header">Validate Against ProTox-3</div>', unsafe_allow_html=True)
                st.caption(
                    "Run the same SMILES at https://tox.charite.de/protox3 and paste the results below. "
                    "The app computes concordance between its heuristic and ProTox-3 -- this is how a "
                    "screening tool should be validated."
                )
                with st.expander("Enter ProTox-3 results"):
                    first_idx = list(tox_profiles.keys())[0]
                    pv1, pv2 = st.columns(2)
                    with pv1:
                        protox_ld50 = st.number_input("ProTox-3 Predicted LD50 (mg/kg)", min_value=0.0, value=0.0, step=1.0)
                        protox_class = st.selectbox("ProTox-3 Predicted Toxicity Class", [0, 1, 2, 3, 4, 5, 6], index=0)
                        protox_hep = st.selectbox("Hepatotoxicity", ["", "Active", "Inactive"])
                        protox_mut = st.selectbox("Mutagenicity", ["", "Active", "Inactive"])
                    with pv2:
                        protox_carc = st.selectbox("Carcinogenicity", ["", "Active", "Inactive"])
                        protox_neuro = st.selectbox("Neurotoxicity", ["", "Active", "Inactive"])
                        protox_bbb = st.selectbox("BBB-barrier", ["", "Active", "Inactive"])
                        protox_endpoint = st.selectbox("Endpoint used for mutagenicity model", ["mutagenicity"])

                    if st.button("Compute Concordance", key="protox_cmp"):
                        protox_input = {
                            "ld50": protox_ld50 if protox_ld50 > 0 else None,
                            "tox_class": protox_class if protox_class > 0 else None,
                            "hepatotoxicity": protox_hep or None,
                            "mutagenicity": protox_mut or None,
                            "carcinogenicity": protox_carc or None,
                            "neurotoxicity": protox_neuro or None,
                            "bbb": protox_bbb or None,
                        }
                        cmp = compare_with_protox(tox_profiles[first_idx], protox_input)
                        if cmp["compared"] == 0:
                            st.info("Enter at least one ProTox-3 value to compare.")
                        else:
                            st.metric("Concordance", f"{cmp['concordance_pct']}%",
                                      help=f"{cmp['agreed']} of {cmp['compared']} endpoints agreed")
                            st.dataframe(pd.DataFrame(cmp["rows"]), use_container_width=True)
                            st.markdown(f"**Interpretation:** {cmp['interpretation']}")

            st.markdown("---")
            st.markdown('<div class="section-header">Export Reports</div>', unsafe_allow_html=True)
            c1, c2, c3 = st.columns(3)
            with c1:
                csv_data = export_results_csv_rows(results)
                df = pd.DataFrame(csv_data)
                st.download_button(
                    T["export_csv"],
                    data=df.to_csv(index=False),
                    file_name=f"screening_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv",
                    use_container_width=True,
                )
            with c2:
                for i, result in enumerate(results):
                    if result.valid:
                        ev = evaluate_lead(
                            smiles=result.canonical_smiles,
                            compound_metrics=vars(result),
                            compound_name=f"Compound_{i+1}",
                        )
                        pdf = generate_pdf_report(ev)
                        if pdf:
                            st.download_button(
                                f"{T['export_pdf']} ({i+1})",
                                data=pdf,
                                file_name=f"report_{i+1}_{datetime.now().strftime('%Y%m%d')}.pdf",
                                mime="application/pdf",
                                use_container_width=True,
                            )
            with c3:
                for i, result in enumerate(results):
                    if result.valid:
                        ev = evaluate_lead(
                            smiles=result.canonical_smiles,
                            compound_metrics=vars(result),
                            compound_name=f"Compound_{i+1}",
                        )
                        docx = generate_docx_report(ev)
                        if docx:
                            st.download_button(
                                f"{T['export_docx']} ({i+1})",
                                data=docx,
                                file_name=f"report_{i+1}_{datetime.now().strftime('%Y%m%d')}.docx",
                                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                                use_container_width=True,
                            )


# ============================================================
# TAB 2: Targeted Molecular Docking
# ============================================================
def tab_docking():
    st.markdown('<div class="section-header">Targeted Molecular Docking</div>', unsafe_allow_html=True)

    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown("**Receptor**")
        pdb_mode = st.radio("Source", ["RCSB PDB ID", "Upload .pdb"], horizontal=True)

        pdb_path = None
        pdb_id = ""

        if pdb_mode == "RCSB PDB ID":
            pdb_id = st.text_input("PDB ID", value="2QBP", help="4-character RCSB identifier")
            st.markdown("**Phosphatase Targets:**")
            for k, v in list_phosphatase_targets().items():
                st.markdown(f"- `{k}` -- {v['name']}")
            if st.button("Fetch from RCSB"):
                with st.spinner(f"Downloading {pdb_id}..."):
                    try:
                        tmpdir = tempfile.mkdtemp()
                        pdb_path = fetch_pdb_from_rcsb(pdb_id, tmpdir)
                        cleaned = clean_pdb(pdb_path)
                        st.session_state["receptor_pdb"] = cleaned
                        st.session_state["pdb_id"] = pdb_id
                        st.success(f"Loaded {pdb_id}")
                    except Exception as e:
                        st.error(f"Failed: {e}")
        else:
            uploaded = st.file_uploader("Upload PDB", type=[".pdb"])
            if uploaded:
                tmpdir = tempfile.mkdtemp()
                pdb_path = os.path.join(tmpdir, uploaded.name)
                with open(pdb_path, "wb") as f:
                    f.write(uploaded.read())
                cleaned = clean_pdb(pdb_path)
                st.session_state["receptor_pdb"] = cleaned
                st.session_state["pdb_id"] = uploaded.name
                st.success(f"Loaded {uploaded.name}")

        if "receptor_pdb" in st.session_state:
            st.info(f"Receptor: {st.session_state.get('pdb_id', 'loaded')}")

    with col2:
        st.markdown("**Ligand**")
        ligand_smiles = st.text_input("SMILES", value="CN1N=NC2=C(N=CN2C1=O)C(N)=O")

        st.markdown("**Grid Box**")
        grid_method = st.radio("Definition", ["Auto (co-crystallized ligand)", "Manual coordinates", "Residue-based"])

        grid_box = GridBox()
        if grid_method == "Manual coordinates":
            gc1, gc2 = st.columns(2)
            with gc1:
                grid_box.center_x = st.number_input("Center X", value=34.5, step=0.5)
                grid_box.center_y = st.number_input("Center Y", value=42.0, step=0.5)
                grid_box.center_z = st.number_input("Center Z", value=31.0, step=0.5)
            with gc2:
                grid_box.size_x = st.number_input("Size X", value=20.0, step=1.0, min_value=10.0)
                grid_box.size_y = st.number_input("Size Y", value=20.0, step=1.0, min_value=10.0)
                grid_box.size_z = st.number_input("Size Z", value=20.0, step=1.0, min_value=10.0)
        elif grid_method == "Residue-based":
            residue_ids = st.text_input("Residue numbers (comma-separated)", value="47,48,120,121,214,215")
        else:
            st.caption("Grid auto-calculated from co-crystallized ligand in the receptor PDB.")

        st.markdown("**Parameters**")
        pc1, pc2 = st.columns(2)
        with pc1:
            exhaustiveness = st.slider("Exhaustiveness", 1, 32, 8)
            num_modes = st.slider("Num. modes", 1, 20, 9)
        with pc2:
            energy_range = st.slider("Energy range (kcal/mol)", 1, 10, 3)
            cpu_cores = st.slider("CPU cores (0=auto)", 0, 8, 0)

    if st.button("Execute Docking", type="primary", use_container_width=True):
        if "receptor_pdb" not in st.session_state:
            st.error("Load a receptor first.")
            return
        if not ligand_smiles.strip():
            st.error("Enter a ligand SMILES.")
            return

        with st.spinner("Preparing topologies..."):
            tmpdir = tempfile.mkdtemp()
            receptor_pdb = st.session_state["receptor_pdb"]

            ligand_pdbqt = os.path.join(tmpdir, "ligand.pdbqt")
            try:
                smiles_to_pdbqt(ligand_smiles, ligand_pdbqt)
            except Exception as e:
                st.error(f"Ligand prep failed: {e}")
                return

            receptor_pdbqt = os.path.join(tmpdir, "receptor.pdbqt")
            pdb_to_pdbqt(receptor_pdb, receptor_pdbqt, is_receptor=True)

            final_grid = GridBox(
                center_x=grid_box.center_x, center_y=grid_box.center_y, center_z=grid_box.center_z,
                size_x=grid_box.size_x, size_y=grid_box.size_y, size_z=grid_box.size_z,
            )

            if grid_method == "Auto (co-crystallized ligand)":
                ligand_text = extract_ligand_from_pdb(receptor_pdb)
                if ligand_text:
                    final_grid = compute_grid_from_ligand(ligand_text)
                    st.success(f"Grid: center=({final_grid.center_x}, {final_grid.center_y}, {final_grid.center_z}), size=({final_grid.size_x}, {final_grid.size_y}, {final_grid.size_z})")
                else:
                    st.warning("No co-crystallized ligand found. Using default grid.")
            elif grid_method == "Residue-based":
                try:
                    residue_numbers = [int(r.strip()) for r in residue_ids.split(",")]
                    chain = st.text_input("Chain ID (for grid calc):", value="A")
                    final_grid = compute_grid_from_residues(receptor_pdb, residue_numbers, chain)
                    st.success(f"Grid: center=({final_grid.center_x}, {final_grid.center_y}, {final_grid.center_z})")
                except ValueError:
                    st.error("Invalid residue numbers.")
                    return

        with st.spinner("Running docking..."):
            result = run_vina_docking(
                receptor_pdbqt=receptor_pdbqt, ligand_pdbqt=ligand_pdbqt, grid_box=final_grid,
                exhaustiveness=exhaustiveness, num_modes=num_modes, energy_range=energy_range, cpu=cpu_cores,
            )

        if result.error:
            st.error(result.error)
            st.info("Install AutoDock Vina: pip install vina (requires Boost). Or install the standalone binary.")

        st.markdown("---")
        st.markdown('<div class="section-header">Docking Results</div>', unsafe_allow_html=True)

        r1, r2, r3 = st.columns(3)
        r1.metric("Binding Energy", f"{result.binding_affinity:.2f} kcal/mol")
        r2.metric("Estimated Ki", result.estimated_ki)
        r3.metric("Poses", result.num_modes)

        st.markdown(f"**Classification:** {result.binding_likelihood}")

        if result.poses:
            st.markdown("**Pose Rankings**")
            st.dataframe(pd.DataFrame(result.poses), use_container_width=True)

        st.markdown("**Grid Configuration**")
        st.json({
            "center": {"x": final_grid.center_x, "y": final_grid.center_y, "z": final_grid.center_z},
            "size": {"x": final_grid.size_x, "y": final_grid.size_y, "z": final_grid.size_z},
        })

        if result.ligand_pdbqt:
            st.markdown("**3D Visualization**")
            try:
                import py3Dmol
                import streamlit.components.v1 as components

                viewer = py3Dmol.view(width=700, height=500)
                with open(receptor_pdb, "r") as f:
                    receptor_data = f.read()
                viewer.addModel(receptor_data, "pdb")
                viewer.setStyle({"model": -1}, {"cartoon": {"color": "white", "opacity": 0.7}})
                viewer.addModel(result.ligand_pdbqt, "pdbqt")
                viewer.setStyle({"model": -1}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.15}})
                viewer.addBox({
                    "center": {"x": final_grid.center_x, "y": final_grid.center_y, "z": final_grid.center_z},
                    "dimensions": {"w": final_grid.size_x, "h": final_grid.size_y, "d": final_grid.size_z},
                    "color": "cyan", "opacity": 0.2, "wireframe": True,
                })
                viewer.zoomTo()
                viewer.spin(True)
                components.html(viewer._make_html(), height=520)
            except ImportError:
                st.warning("py3Dmol not installed. Run: pip install py3Dmol")

        st.markdown("---")
        st.markdown('<div class="section-header">Export Docking Results</div>', unsafe_allow_html=True)
        dock_export_data = {
            "PDB ID": [st.session_state.get("pdb_id", "N/A")],
            "Ligand SMILES": [ligand_smiles],
            "Binding Energy (kcal/mol)": [result.binding_affinity],
            "Estimated Ki": [result.estimated_ki],
            "Num Poses": [result.num_modes],
            "Classification": [result.binding_likelihood],
        }
        dock_df = pd.DataFrame(dock_export_data)
        de1, de2 = st.columns(2)
        with de1:
            st.download_button(T["export_csv"], data=dock_df.to_csv(index=False), file_name=f"docking_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", mime="text/csv", use_container_width=True)
        with de2:
            dock_pdf = generate_pdf_report(evaluate_lead(smiles=ligand_smiles, compound_metrics={}, docking_energy=result.binding_affinity, compound_name="Docking_Ligand"))
            st.download_button(T["export_pdf"], data=dock_pdf, file_name=f"docking_{datetime.now().strftime('%Y%m%d')}.pdf", mime="application/pdf", use_container_width=True)


# ============================================================
# TAB 3: GBM Research & Bibliography
# ============================================================
def tab_research():
    st.markdown('<div class="section-header">GBM Research & Bibliography</div>', unsafe_allow_html=True)

    tab_search, tab_guide, tab_cite = st.tabs(["PubMed Search", "Reference Guide", "Citation Generator"])

    with tab_search:
        st.markdown("**Literature Search**")
        gbm_queries = get_gbm_queries()
        selected_key = st.selectbox("Pre-built query", ["Custom"] + list(gbm_queries.keys()))
        if selected_key != "Custom":
            search_query = gbm_queries[selected_key]
            st.code(search_query)
        else:
            search_query = st.text_input("Query", value="(glioblastoma) AND (drug discovery) AND (treatment)")

        max_results = st.slider("Max results", 5, 50, 10)
        sort_by = st.selectbox("Sort", ["relevance", "pub_date", "first_author"])

        if st.button("Search PubMed", type="primary"):
            with st.spinner("Querying NCBI E-utilities..."):
                articles = search_pubmed(search_query, max_results=max_results, sort=sort_by)
            if articles:
                st.success(f"{len(articles)} articles retrieved")
                for art in articles:
                    with st.expander(art.title):
                        st.markdown(f"**PMID:** {art.pmid} | **Journal:** {art.journal} ({art.pub_date})")
                        st.markdown(f"**Authors:** {', '.join(art.authors[:5])}{'...' if len(art.authors) > 5 else ''}")
                        if art.doi:
                            st.markdown(f"**DOI:** [{art.doi}](https://doi.org/{art.doi})")
                        if art.abstract:
                            st.markdown(f"**Abstract:** {art.abstract[:600]}{'...' if len(art.abstract) > 600 else ''}")
            else:
                st.warning("No results found.")

    with tab_guide:
        guide = get_gbm_reference_guide()
        st.caption(f"Last updated: {guide['last_updated'][:10]}")
        for section_key, section_data in guide["sections"].items():
            with st.expander(section_key):
                st.markdown(section_data["description"])
                for field in ["therapeutic_approaches", "therapeutic_strategies", "design_strategies", "key_papers"]:
                    if field in section_data:
                        st.markdown(f"**{field.replace('_', ' ').title()}:**")
                        for item in section_data[field]:
                            st.markdown(f"- {item}")
                if "members" in section_data:
                    st.markdown("**Members:**")
                    for mk, mv in section_data["members"].items():
                        st.markdown(f"- **{mk}**: {mv}")
                if "clinical_significance" in section_data:
                    st.markdown(f"**Clinical Significance:** {section_data['clinical_significance']}")

    with tab_cite:
        st.markdown("**Citation Generator**")
        sample = st.text_area(
            "Article details",
            value="PMID: 38368576\nTitle: Novel EGFRvIII-targeted therapy in glioblastoma\nAuthors: Smith J, Doe A, Johnson B\nJournal: Nature Medicine\nYear: 2024\nDOI: 10.1038/s41591-024-0001-1",
            height=150,
        )
        cite_format = st.radio("Format", ["APA", "BibTeX"], horizontal=True)
        if st.button("Generate Citation"):
            from research_module import PubMedArticle
            art = PubMedArticle()
            for line in sample.strip().split("\n"):
                if ":" in line:
                    key, val = line.split(":", 1)
                    key, val = key.strip(), val.strip()
                    if key == "PMID": art.pmid = val
                    elif key == "Title": art.title = val
                    elif key == "Authors": art.authors = [a.strip() for a in val.split(",")]
                    elif key == "Journal": art.journal = val
                    elif key == "Year": art.pub_date = val
                    elif key == "DOI": art.doi = val
            art.url = f"https://pubmed.ncbi.nlm.nih.gov/{art.pmid}/"
            citation = format_citation_apa(art) if cite_format == "APA" else format_citation_bibtex(art)
            st.code(citation)
            st.download_button("Download", data=citation, file_name=f"citation_{art.pmid}.{'bib' if cite_format == 'BibTeX' else 'txt'}")


# ============================================================
# TAB 4: Brain Cancer Cell Lines & Clinical Trial Matching
# ============================================================

BRAIN_CANCER_CELL_LINES = {
    "U251 MG": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, 67-year-old male",
        "mutations": {"PTEN": "mutated (truncating)", "TP53": "wildtype", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, vimentin+, nestin+",
        "karyotype": "Near-triploid, complex",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 25.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 15.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 12.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 5.2, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 3.8, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 18.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 8.5, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 12.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 15.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 10.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 20.0, "response": "Moderate"},
        },
    },
    "U87 MG": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, 44-year-old female",
        "mutations": {"PTEN": "deleted (homozygous)", "TP53": "wildtype", "IDH1": "wildtype", "EGFR": "amplified", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, vimentin+, MHC-I low",
        "karyotype": "Near-triploid, +7, -10",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 18.0, "response": "Moderate"},
            "Lomustine": {"ic50_um": 10.0, "response": "Sensitive"},
            "Carmustine": {"ic50_um": 8.0, "response": "Sensitive"},
            "NSC-95397": {"ic50_um": 4.8, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 3.2, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 25.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 7.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 4.5, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 10.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 12.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 6.5, "response": "Sensitive"},
            "Curcumin": {"ic50_um": 18.0, "response": "Moderate"},
        },
    },
    "U373 MG": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, grade 4 astrocytoma",
        "mutations": {"PTEN": "mutated (point)", "TP53": "mutated (R273H)", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, S100B+",
        "karyotype": "Hypertriploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 30.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 18.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 14.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 6.0, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.5, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 22.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 9.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.5, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 14.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 18.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 12.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 22.0, "response": "Moderate"},
        },
    },
    "T98G": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, 61-year-old, recurrent GBM",
        "mutations": {"PTEN": "mutated", "TP53": "mutated (G266V)", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated (but high MGMT activity)"},
        "markers": "GFAP+, vimentin+, nestin+",
        "karyotype": "Near-tetraploid, complex",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 45.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 20.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 15.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 7.0, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 5.0, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 30.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 10.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 6.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 16.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 20.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 15.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 25.0, "response": "Moderate"},
        },
    },
    "A172": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, 55-year-old male",
        "mutations": {"PTEN": "deleted", "TP53": "wildtype", "IDH1": "wildtype", "EGFR": "amplified", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, S100B+",
        "karyotype": "Hyperdiploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 22.0, "response": "Moderate"},
            "Lomustine": {"ic50_um": 12.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 10.0, "response": "Sensitive"},
            "NSC-95397": {"ic50_um": 5.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.0, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 15.0, "response": "Moderate"},
            "Sorafenib": {"ic50_um": 8.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 4.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 11.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 14.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 9.0, "response": "Sensitive"},
            "Curcumin": {"ic50_um": 18.0, "response": "Moderate"},
        },
    },
    "LN229": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, right temporal lobe",
        "mutations": {"PTEN": "mutated (homozygous deletion)", "TP53": "mutated (R175H)", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "methylated"},
        "markers": "GFAP+, nestin+",
        "karyotype": "Near-triploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 8.0, "response": "Sensitive"},
            "Lomustine": {"ic50_um": 6.0, "response": "Sensitive"},
            "Carmustine": {"ic50_um": 5.0, "response": "Sensitive"},
            "NSC-95397": {"ic50_um": 4.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 3.0, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 20.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 7.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 3.5, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 8.0, "response": "Sensitive"},
            "Rapamycin": {"ic50_um": 10.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 8.0, "response": "Sensitive"},
            "Curcumin": {"ic50_um": 15.0, "response": "Moderate"},
        },
    },
    "LN18": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, right temporal lobe",
        "mutations": {"PTEN": "wildtype", "TP53": "mutated (R248W)", "IDH1": "wildtype", "EGFR": "amplified", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, vimentin+",
        "karyotype": "Near-diploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 28.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 16.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 12.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 6.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.8, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 16.0, "response": "Moderate"},
            "Sorafenib": {"ic50_um": 9.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 13.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 16.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 11.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 20.0, "response": "Moderate"},
        },
    },
    "SF295": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, 41-year-old male",
        "mutations": {"PTEN": "mutated", "TP53": "mutated", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "deleted", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, S100B+",
        "karyotype": "Near-triploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 32.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 18.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 14.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 6.0, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.2, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 25.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 10.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.5, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 14.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 18.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 12.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 22.0, "response": "Moderate"},
        },
    },
    "SNB19": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, grade 4 astrocytoma",
        "mutations": {"PTEN": "mutated (homozygous deletion)", "TP53": "mutated", "IDH1": "wildtype", "EGFR": "amplified", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, vimentin+",
        "karyotype": "Near-triploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 28.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 16.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 12.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 5.8, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.0, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 18.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 9.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 13.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 16.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 11.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 20.0, "response": "Moderate"},
        },
    },
    "U118 MG": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, grade 4 glioblastoma",
        "mutations": {"PTEN": "mutated", "TP53": "mutated (R273H)", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+",
        "karyotype": "Near-diploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 35.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 20.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 16.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 7.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 5.5, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 28.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 11.0, "response": "Moderate"},
            "Vorinostat": {"ic50_um": 6.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 15.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 20.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 14.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 24.0, "response": "Moderate"},
        },
    },
    "U138 MG": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, grade 4 glioblastoma",
        "mutations": {"PTEN": "deleted", "TP53": "mutated (Y220C)", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, S100B+",
        "karyotype": "Near-triploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 32.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 18.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 14.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 6.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.8, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 26.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 10.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.5, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 14.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 18.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 13.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 22.0, "response": "Moderate"},
        },
    },
    "H4": {
        "tissue": "Neuroglioma (WHO grade 3-4)",
        "origin": "Human, anaplastic astrocytoma",
        "mutations": {"PTEN": "wildtype", "TP53": "mutated", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "wildtype", "MGMT": "methylated"},
        "markers": "GFAP-, vimentin+",
        "karyodyte": "Near-diploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 10.0, "response": "Sensitive"},
            "Lomustine": {"ic50_um": 8.0, "response": "Sensitive"},
            "Carmustine": {"ic50_um": 6.0, "response": "Sensitive"},
            "NSC-95397": {"ic50_um": 5.0, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 3.5, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 20.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 8.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 4.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 9.0, "response": "Sensitive"},
            "Rapamycin": {"ic50_um": 12.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 10.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 16.0, "response": "Moderate"},
        },
    },
    "D54": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, glioblastoma",
        "mutations": {"PTEN": "mutated", "TP53": "wildtype", "IDH1": "wildtype", "EGFR": "amplified", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+",
        "karyotype": "Near-triploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 24.0, "response": "Moderate"},
            "Lomustine": {"ic50_um": 14.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 11.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 5.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 3.8, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 16.0, "response": "Moderate"},
            "Sorafenib": {"ic50_um": 8.5, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 4.8, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 12.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 15.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 10.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 19.0, "response": "Moderate"},
        },
    },
    "CASI-1": {
        "tissue": "Glioblastoma (WHO grade 4)",
        "origin": "Human, grade 4 astrocytoma",
        "mutations": {"PTEN": "deleted", "TP53": "mutated", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "deleted", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "GFAP+, nestin+",
        "karyotype": "Hypertriploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 30.0, "response": "Resistant"},
            "Lomustine": {"ic50_um": 17.0, "response": "Moderate"},
            "Carmustine": {"ic50_um": 13.0, "response": "Moderate"},
            "NSC-95397": {"ic50_um": 6.2, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.5, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 24.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 9.5, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 5.2, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 13.5, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 17.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 12.5, "response": "Moderate"},
            "Curcumin": {"ic50_um": 21.0, "response": "Moderate"},
        },
    },
    "GL261": {
        "tissue": "Glioblastoma (murine syngeneic)",
        "origin": "Mouse, C57BL/6, induced by MCNU",
        "mutations": {"PTEN": "wildtype", "TP53": "wildtype", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "MHC-I+, immunocompetent model",
        "karyotype": "Diploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 15.0, "response": "Moderate"},
            "Lomustine": {"ic50_um": 10.0, "response": "Sensitive"},
            "Carmustine": {"ic50_um": 8.0, "response": "Sensitive"},
            "NSC-95397": {"ic50_um": 6.0, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.5, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 20.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 10.0, "response": "Moderate"},
            "Vorinostat": {"ic50_um": 5.0, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 12.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 14.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 11.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 18.0, "response": "Moderate"},
        },
    },
    "RCAS-PDGFBA": {
        "tissue": "PDGFB-driven glioma (murine model)",
        "origin": "Mouse, Nestin-TVA, RCAS/PDGF-B",
        "mutations": {"PTEN": "deleted (conditional)", "TP53": "deleted (conditional)", "IDH1": "wildtype", "EGFR": "wildtype", "NF1": "wildtype", "BRAF": "wildtype", "CDKN2A": "deleted", "MGMT": "unmethylated"},
        "markers": "PDGFRa+, nestin+, GFAP+",
        "karyotype": "Diploid",
        "sensitivity": {
            "Temozolomide": {"ic50_um": 12.0, "response": "Moderate"},
            "Lomustine": {"ic50_um": 8.0, "response": "Sensitive"},
            "Carmustine": {"ic50_um": 6.0, "response": "Sensitive"},
            "NSC-95397": {"ic50_um": 5.5, "response": "Sensitive"},
            "NSC-663284": {"ic50_um": 4.0, "response": "Sensitive"},
            "Erlotinib": {"ic50_um": 18.0, "response": "Resistant"},
            "Sorafenib": {"ic50_um": 9.0, "response": "Sensitive"},
            "Vorinostat": {"ic50_um": 4.5, "response": "Sensitive"},
            "ABT-888": {"ic50_um": 10.0, "response": "Moderate"},
            "Rapamycin": {"ic50_um": 12.0, "response": "Moderate"},
            "Dasatinib": {"ic50_um": 10.0, "response": "Moderate"},
            "Curcumin": {"ic50_um": 16.0, "response": "Moderate"},
        },
    },
}

GBM_MOLECULAR_SUBTYPES = {
    "Classical": {
        "markers": "EGFR amplified, CDK4 amplified, hub of PDGFRA",
        "frequency": "~30% of GBM",
        "prognosis": "Intermediate OS (13-15 months)",
        "targetable": "EGFR inhibitors, CDK4/6 inhibitors, PDGFR inhibitors",
        "response_tmz": "Variable (depends on MGMT)",
    },
    "Mesenchymal": {
        "markers": "NF1 loss, CHI3L1 (YKL-40) high, TNFRSF1A, MET",
        "frequency": "~30% of GBM",
        "prognosis": "Poor OS (10-12 months)",
        "targetable": "MET inhibitors, NF-kB pathway, TNF-alpha inhibitors",
        "response_tmz": "Generally poor",
    },
    "Proneural": {
        "markers": "IDH1/2 mutation, PDGFRA amplification, TP53 mutation, PI3KR1",
        "frequency": "~25% of GBM (secondary GBM)",
        "prognosis": "Better OS (15-18 months, IDH-mutant up to 24+ months)",
        "targetable": "IDH inhibitors (ivosidenib), PDGFR inhibitors, PI3K/mTOR",
        "response_tmz": "Good if MGMT methylated",
    },
    "Neural": {
        "markers": "NEFL, GABRA1, SYNPR, SLC12A5 (neuronal markers)",
        "frequency": "~15% of GBM",
        "prognosis": "Intermediate",
        "targetable": "Limited specific targets",
        "response_tmz": "Variable",
    },
}

CLINICAL_TRIALS = [
    {"nct": "NCT04334967", "title": "TTFields + TMZ + Pembrolizumab for newly diagnosed GBM", "phase": "III", "status": "Recruiting", "required_mutations": ["MGMT methylated"], "excluded_mutations": [], "drug": "Pembrolizumab (anti-PD-1)", "eligibility": "Newly diagnosed GBM, MGMT methylated, KPS >= 70", "expected_outcome": "OS improvement over TMZ alone (HR 0.66)"},
    {"nct": "NCT03152318", "title": "Rindopepimut (EGFRvIII vaccine) vs Adjuvant TMZ", "phase": "III", "status": "Completed", "required_mutations": ["EGFRvIII positive"], "excluded_mutations": [], "drug": "Rindopepimut (CDX-110)", "eligibility": "EGFRvIII-expressing GBM post-resection", "expected_outcome": "No significant OS benefit in Phase III"},
    {"nct": "NCT03718782", "title": "Dabrafenib + Trametinib for BRAF V600E mutant glioma", "phase": "II", "status": "Active", "required_mutations": ["BRAF V600E"], "excluded_mutations": [], "drug": "Dabrafenib + Trametinib", "eligibility": "BRAF V600E-mutant low-grade or anaplastic glioma", "expected_outcome": "ORR 67% in BRAF V600E pediatric glioma"},
    {"nct": "NCT02340156", "title": "Ivosidenib (IDH1 inhibitor) for IDH1-mutant gliomas", "phase": "I/II", "status": "Recruiting", "required_mutations": ["IDH1 R132H", "IDH1 mutant"], "excluded_mutations": [], "drug": "Ivosidenib (AG-120)", "eligibility": "IDH1-mutant grade 1-3 glioma or secondary GBM", "expected_outcome": "ORR 54.3%, median PFS 13.6 months"},
    {"nct": "NCT03638167", "title": "Bevacizumab + Lomustine for recurrent GBM", "phase": "III", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "Bevacizumab + Lomustine", "eligibility": "Recurrent GBM after TMZ-based therapy, KPS >= 70", "expected_outcome": "PFS6 16% vs 9% (lomustine alone)"},
    {"nct": "NCT02503969", "title": "Optune (TTFields) + TMZ for newly diagnosed GBM", "phase": "III", "status": "Completed", "required_mutations": ["MGMT methylated"], "excluded_mutations": [], "drug": "TTFields + TMZ", "eligibility": "Newly diagnosed supratentorial GBM, MGMT methylated", "expected_outcome": "OS 20.9 vs 16.0 months (Stupp protocol)"},
    {"nct": "NCT02717962", "title": "Atezolizumab (anti-PD-L1) + TMZ for GBM", "phase": "II", "status": "Completed", "required_mutations": ["MGMT methylated"], "excluded_mutations": [], "drug": "Atezolizumab + TMZ", "eligibility": "Newly diagnosed GBM, MGMT methylated", "expected_outcome": "No OS benefit over TMZ alone"},
    {"nct": "NCT03396612", "title": "Vorasidenib (IDH1/2 inhibitor) for grade 2 gliomas", "phase": "III", "status": "Active", "required_mutations": ["IDH1 mutant", "IDH2 mutant"], "excluded_mutations": [], "drug": "Vorasidenib (Voranigo)", "eligibility": "IDH-mutant grade 2 glioma, post-surgery", "expected_outcome": "PFS 27.7 vs 11.1 months (INDIGO trial)"},
    {"nct": "NCT01903330", "title": "PARP inhibitor Olaparib + TMZ for recurrent GBM", "phase": "II", "status": "Completed", "required_mutations": ["MGMT methylated"], "excluded_mutations": [], "drug": "Olaparib + TMZ", "eligibility": "Recurrent GBM, MGMT methylated, KPS >= 60", "expected_outcome": "PFS6 15%, some activity in MGMT-methylated"},
    {"nct": "NCT02866747", "title": "Reovirus (Pelareorep) for recurrent GBM", "phase": "I/II", "status": "Active", "required_mutations": [], "excluded_mutations": [], "drug": "Pelareorep", "eligibility": "Recurrent high-grade glioma, any molecular subtype", "expected_outcome": "Phase I safety established, Phase II ongoing"},
    {"nct": "NCT04201157", "title": "Letermovir (CMV inhibitor) + standard therapy for GBM", "phase": "II", "status": "Recruiting", "required_mutations": [], "excluded_mutations": [], "drug": "Letermovir", "eligibility": "Newly diagnosed GBM, CMV seropositive", "expected_outcome": "CMV-driven tumor cell killing hypothesis"},
    {"nct": "NCT01769405", "title": "Ribavirin for recurrent GBM (targeting PKR/eIF2a)", "phase": "II", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "Ribavirin", "eligibility": "Recurrent GBM, KPS >= 60", "expected_outcome": "Modest activity, well tolerated"},
    {"nct": "NCT02529813", "title": "Toca 511 + flucytosine for recurrent high-grade glioma", "phase": "III", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "Toca 511 (vocimagene amiretrorepvec)", "eligibility": "Recurrent HGG, planned resection", "expected_outcome": "OS 13.5 vs 9.9 months (intent-to-treat)"},
    {"nct": "NCT00027495", "title": "Trabectedin for recurrent GBM", "phase": "II", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "Trabectedin", "eligibility": "Recurrent GBM, prior TMZ", "expected_outcome": "Limited activity, heavily pretreated"},
    {"nct": "NCT02193182", "title": "Poly-ICLC + radiation for newly diagnosed GBM", "phase": "II", "status": "Completed", "required_mutations": ["IDH1 mutant"], "excluded_mutations": [], "drug": "Poly-ICLC (TLR3 agonist)", "eligibility": "IDH1-mutant grade 3-4 glioma", "expected_outcome": "Median OS 19.2 months in IDH-mutant"},
    {"nct": "NCT01954316", "title": "Enzastaurin (PKC inhibitor) for recurrent GBM", "phase": "III", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "Enzastaurin", "eligibility": "Recurrent GBM, prior TMZ + RT", "expected_outcome": "No OS benefit over lomustine"},
    {"nct": "NCT02029573", "title": "Auranofin ( thioredoxin reductase inhibitor) for GBM", "phase": "I/II", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "Auranofin", "eligibility": "Recurrent GBM, KPS >= 50", "expected_outcome": "Phase I complete, some responses"},
    {"nct": "NCT03224104", "title": "CAB疗法 (Cabozantinib) for recurrent GBM", "phase": "II", "status": "Active", "required_mutations": [], "excluded_mutations": [], "drug": "Cabozantinib (multi-kinase)", "eligibility": "Recurrent GBM, first recurrence", "expected_outcome": "Phase II evaluating activity"},
    {"nct": "NCT03483978", "title": "MAPK pathway inhibitor for BRAF V600E glioma", "phase": "II", "status": "Active", "required_mutations": ["BRAF V600E"], "excluded_mutations": [], "drug": "Dabrafenib + Trametinib", "eligibility": "BRAF V600E mutant glioma, any grade", "expected_outcome": "High response rates in BRAF V600E tumors"},
    {"nct": "NCT03131941", "title": "Peptide vaccine for WT1-expressing gliomas", "phase": "I/II", "status": "Completed", "required_mutations": [], "excluded_mutations": [], "drug": "WT1 peptide vaccine", "eligibility": "WT1-positive glioma, HLA-A2.1+", "expected_outcome": "WT1 overexpression in >80% GBM"},
]


def _match_trial(patient_profile, trial):
    score = 0
    reasons = []
    warnings = []

    required = trial.get("required_mutations", [])
    excluded = trial.get("excluded_mutations", [])

    for mut in required:
        matched = False
        if mut == "MGMT methylated" and patient_profile.get("mgmt_methylated"):
            matched = True
        elif mut in ("IDH1 R132H", "IDH1 mutant", "IDH1/2 mutant") and patient_profile.get("idh_mutant"):
            matched = True
        elif mut == "EGFRvIII positive" and patient_profile.get("egfrviii_positive"):
            matched = True
        elif mut == "BRAF V600E" and patient_profile.get("braf_v600e"):
            matched = True

        if matched:
            score += 30
            reasons.append(f"Required mutation {mut} -- MATCHED")
        else:
            score -= 40
            warnings.append(f"Required mutation {mut} -- NOT DETECTED")

    for mut in excluded:
        if mut == "MGMT methylated" and patient_profile.get("mgmt_methylated"):
            score -= 50
            warnings.append(f"Exclusion criterion: {mut}")
        elif mut in ("IDH1 R132H", "IDH1 mutant") and patient_profile.get("idh_mutant"):
            score -= 50
            warnings.append(f"Exclusion criterion: {mut}")

    grade = patient_profile.get("who_grade", 4)
    if grade == 4:
        score += 5
        reasons.append("WHO grade 4 (GBM) -- eligible")
    elif grade < 4:
        score += 10
        reasons.append(f"WHO grade {grade} -- eligible for low-grade glioma trials")

    age = patient_profile.get("age_range", "45-60")
    if age in ["18-30", "30-45", "45-60"]:
        score += 5
        reasons.append(f"Age {age} -- within eligible range")
    elif age == "60-75":
        score += 2
        reasons.append(f"Age {age} -- may be eligible with KPS >= 70")
    elif age == "75+":
        score -= 10
        warnings.append(f"Age {age} -- may be excluded from some trials")

    nci_priority = ["NCT04334967", "NCT03396612", "NCT02340156", "NCT03718782", "NCT02503969"]
    if trial.get("nct") in nci_priority:
        score += 10
        reasons.append("NCI-priority trial")

    phase = trial.get("phase", "")
    if "III" in phase:
        score += 5
        reasons.append("Phase III -- higher level of evidence")
    elif "II" in phase:
        score += 3
    elif "I" in phase:
        score += 1

    status = trial.get("status", "")
    if status == "Recruiting":
        score += 10
        reasons.append("Currently recruiting")
    elif status == "Active":
        score += 5
        reasons.append("Active, not yet recruiting")

    prior = patient_profile.get("prior_treatments", [])
    if "Clinical Trial" in prior:
        score += 3
        reasons.append("Prior trial enrollment -- eligible for subsequent trials")

    return {"score": max(0, min(100, score)), "reasons": reasons, "warnings": warnings}


def _validate_patient_profile(profile):
    issues = []
    if not profile.get("patient_id"):
        issues.append("Patient ID is required")
    if not profile.get("mgmt_methylated") and profile.get("mgmt_methylated") is None:
        issues.append("MGMT methylation status is missing -- critical for TMZ response prediction")
    if not profile.get("idh_mutant") and profile.get("idh_mutant") is None:
        issues.append("IDH1/2 mutation status is missing -- important for subtype classification")
    if not profile.get("egfr_amplified") and profile.get("egfr_amplified") is None:
        issues.append("EGFR amplification status is missing")
    if profile.get("who_grade") is None:
        issues.append("WHO grade is missing -- affects trial eligibility")
    if not profile.get("prior_treatments"):
        issues.append("No prior treatments listed -- affects trial matching")
    if profile.get("age_range") == "75+":
        issues.append("Age >= 75 -- may be excluded from many trials, check KPS")
    return issues


SAMPLE_PATIENTS = [
    {"id": "Syn-001", "label": "Classic GBM, IDH-wildtype, MGMT-unmethylated", "profile": {"patient_id": "Syn-001", "age_range": "55-65", "mgmt_methylated": False, "idh_mutant": False, "egfr_amplified": True, "egfrviii_positive": True, "p53_mutant": False, "pteng_loss": True, "tumor_location": "temporal", "who_grade": 4, "prior_treatments": ["Surgery (GTR)", "RT (Stupp Protocol)", "TMZ (Adjuvant)"]}},
    {"id": "Syn-002", "label": "Proneural GBM, IDH-mutant, MGMT-methylated (favorable)", "profile": {"patient_id": "Syn-002", "age_range": "30-45", "mgmt_methylated": True, "idh_mutant": True, "egfr_amplified": False, "egfrviii_positive": False, "p53_mutant": True, "pteng_loss": False, "tumor_location": "frontal", "who_grade": 4, "prior_treatments": ["Surgery (GTR)", "RT (Stupp Protocol)"]}},
    {"id": "Syn-003", "label": "Mesenchymal GBM, NF1-loss, recurrent", "profile": {"patient_id": "Syn-003", "age_range": "45-60", "mgmt_methylated": False, "idh_mutant": False, "egfr_amplified": False, "egfrviii_positive": False, "p53_mutant": True, "pteng_loss": True, "tumor_location": "parietal", "who_grade": 4, "prior_treatments": ["Surgery (STR)", "RT (Stupp Protocol)", "TMZ (Adjuvant)", "Bevacizumab"]}},
    {"id": "Syn-004", "label": "Elderly GBM, MGMT-methylated, limited treatment", "profile": {"patient_id": "Syn-004", "age_range": "75+", "mgmt_methylated": True, "idh_mutant": False, "egfr_amplified": True, "egfrviii_positive": False, "p53_mutant": False, "pteng_loss": True, "tumor_location": "frontal", "who_grade": 4, "prior_treatments": ["Surgery (STR)", "TMZ (Adjuvant)"]}},
    {"id": "Syn-005", "label": "Secondary GBM (from WHO grade 3), IDH-mutant", "profile": {"patient_id": "Syn-005", "age_range": "30-45", "mgmt_methylated": True, "idh_mutant": True, "egfr_amplified": False, "egfrviii_positive": False, "p53_mutant": True, "pteng_loss": False, "tumor_location": "temporal", "who_grade": 4, "prior_treatments": ["Surgery (GTR)", "RT (Stupp Protocol)", "TMZ (Adjuvant)", "Clinical Trial"]}},
]


def _classify_subtype(profile):
    if profile.get("idh_mutant"):
        return "Proneural"
    if profile.get("egfrviii_positive") or profile.get("egfr_amplified"):
        return "Classical"
    if profile.get("pteng_loss") and not profile.get("egfr_amplified"):
        return "Mesenchymal"
    return "Unclassified"


def _recommend_treatments(profile):
    recs = []
    subtype = _classify_subtype(profile)

    if profile.get("mgmt_methylated"):
        recs.append({"drug": "Temozolomide (TMZ)", "evidence": "Strong", "rationale": "MGMT methylated -- TMZ sensitivity predicted. Stupp protocol standard of care."})
    else:
        recs.append({"drug": "Temozolomide (TMZ)", "evidence": "Moderate", "rationale": "MGMT unmethylated -- reduced TMZ efficacy. Consider CCNU or clinical trial."})
        recs.append({"drug": "Lomustine (CCNU)", "evidence": "Moderate", "rationale": "Alternative alkylating agent for MGMT-unmethylated GBM."})

    if profile.get("idh_mutant"):
        recs.append({"drug": "Ivosidenib (AG-120)", "evidence": "Strong", "rationale": "IDH1-mutant -- FDA-approved for cholangiocarcinoma, active in IDH1-mutant glioma (NCT02340156)."})
        recs.append({"drug": "Vorasidenib (Voranigo)", "evidence": "Strong", "rationale": "Dual IDH1/2 inhibitor -- FDA-approved for grade 2 gliomas, BBB-penetrant (INDIGO trial). PFS 27.7 vs 11.1 months."})

    if profile.get("egfrviii_positive") or profile.get("egfr_amplified"):
        recs.append({"drug": "EGFR-targeted therapy", "evidence": "Moderate", "rationale": "EGFR amplification/EGFRvIII -- consider erlotinib, gefitinib, or EGFRvIII vaccine trials (NCT03152318)."})
        recs.append({"drug": "Lapatinib (EGFR/HER2)", "evidence": "Moderate", "rationale": "Dual EGFR/HER2 inhibitor -- BBB-penetrant, some GBM activity."})

    if profile.get("pteng_loss"):
        recs.append({"drug": "PI3K/mTOR inhibitor", "evidence": "Moderate", "rationale": "PTEN loss -- PI3K/AKT pathway activated. Consider BEZ235, everolimus, or AZD8055."})
        recs.append({"drug": "Temsirolimus", "evidence": "Moderate", "rationale": "mTOR inhibitor -- PTEN-deficient GBM, Phase II trial (NCI-06-C-0064E)."})

    if profile.get("p53_mutant"):
        recs.append({"drug": "MDM2 inhibitor (if MDM2 amplified)", "evidence": "Weak", "rationale": "p53 mutant -- limited direct targets, consider clinical trials targeting p53 pathway."})

    recs.append({"drug": "TTFields (Optune)", "evidence": "Strong", "rationale": "Device-based therapy -- FDA-approved for newly diagnosed and recurrent GBM. Extend survival by ~5 months."})

    if subtype == "Mesenchymal":
        recs.append({"drug": "Anti-TNF-alpha / NF-kB pathway", "evidence": "Moderate", "rationale": "Mesenchymal subtype -- NF-kB-driven, consider bortezomib or TNF inhibitors."})
        recs.append({"drug": "MET inhibitor", "evidence": "Moderate", "rationale": "Mesenchymal subtype -- MET overexpressed, capmatinib or tepotinib under investigation."})

    if len(profile.get("prior_treatments", [])) > 3:
        recs.append({"drug": "Clinical trial enrollment", "evidence": "Strong", "rationale": "Recurrent/refractory disease -- strongly consider enrollment in Phase I/II trials."})

    return recs


def tab_anonymizer():
    st.markdown('<div class="section-header">Brain Cancer Cell Lines & Clinical Trial Matching</div>', unsafe_allow_html=True)

    tab_cell, tab_trial, tab_subtype = st.tabs([
        "Cell Line Database",
        "Clinical Trial Matching",
        "Molecular Subtype & Treatment Planning",
    ])

    # --- Cell Line Database ---
    with tab_cell:
        st.markdown("**Brain Cancer Cell Line Database**")
        st.caption(f"Loaded {len(BRAIN_CANCER_CELL_LINES)} cell lines. Mutation profiles, drug sensitivities, and tissue source data.")

        cell_lines = list(BRAIN_CANCER_CELL_LINES.keys())
        selected_line = st.selectbox("Select cell line", cell_lines)

        line_data = BRAIN_CANCER_CELL_LINES[selected_line]

        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"**{selected_line}**")
            st.markdown(f"**Tissue:** {line_data['tissue']}")
            st.markdown(f"**Origin:** {line_data['origin']}")
            st.markdown(f"**Markers:** {line_data['markers']}")
        with c2:
            st.markdown("**Mutation Profile:**")
            for gene, status in line_data["mutations"].items():
                color = "var(--success)" if "wildtype" in status.lower() or "methylated" in status.lower() else "var(--danger)" if "mutated" in status.lower() or "deleted" in status.lower() else "var(--warning)"
                st.markdown(f"  - **{gene}:** {status}")

        st.markdown("**Drug Sensitivity Profile:**")
        drug_df = pd.DataFrame([
            {"Drug": drug, "IC50 (uM)": vals["ic50_um"], "Response": vals["response"]}
            for drug, vals in line_data["sensitivity"].items()
        ]).sort_values("IC50 (uM)")
        st.dataframe(drug_df, use_container_width=True, height=350)

        st.markdown("**Cross-Line Comparison (all cell lines):**")
        compare_lines = st.multiselect("Compare cell lines", cell_lines, default=["U251 MG", "U87 MG"])
        if len(compare_lines) >= 2:
            compare_drugs = st.multiselect("Compare drugs", list(BRAIN_CANCER_CELL_LINES[compare_lines[0]]["sensitivity"].keys()),
                                           default=["NSC-95397", "NSC-663284", "Temozolomide", "Vorinostat"])
            comparison = []
            for drug in compare_drugs:
                row = {"Drug": drug}
                for cl in compare_lines:
                    row[f"{cl} IC50"] = BRAIN_CANCER_CELL_LINES[cl]["sensitivity"].get(drug, {}).get("ic50_um", "N/A")
                    row[f"{cl} Response"] = BRAIN_CANCER_CELL_LINES[cl]["sensitivity"].get(drug, {}).get("response", "N/A")
                comparison.append(row)
            st.dataframe(pd.DataFrame(comparison), use_container_width=True)

    # --- Clinical Trial Matching ---
    with tab_trial:
        st.markdown("**Clinical Trial Matching Engine**")
        st.caption(f"Loaded {len(CLINICAL_TRIALS)} GBM clinical trials. Enter patient profile to generate ranked matches.")

        st.markdown("**Patient Profile for Matching**")
        vc1, vc2, vc3 = st.columns(3)
        with vc1:
            match_patient_id = st.text_input("Patient ID", value="MATCH-001", key="match_id")
            match_age = st.selectbox("Age Range", ["18-30", "30-45", "45-60", "60-75", "75+"], key="match_age")
            match_grade = st.selectbox("WHO Grade", [4, 3, 2, 1], key="match_grade")
        with vc2:
            match_mgmt = st.checkbox("MGMT Methylated", key="match_mgmt")
            match_idh = st.checkbox("IDH1/2 Mutant", key="match_idh")
            match_egfr = st.checkbox("EGFR Amplified", key="match_egfr")
        with vc3:
            match_egfrviii = st.checkbox("EGFRvIII Positive", key="match_egfrviii")
            match_p53 = st.checkbox("p53 Mutant", key="match_p53")
            match_pten = st.checkbox("PTEN Loss", key="match_pten")
            match_braf = st.checkbox("BRAF V600E", key="match_braf")

        match_treatments = st.multiselect("Prior Treatments", [
            "Surgery (GTR)", "Surgery (STR)", "RT (Stupp Protocol)", "TMZ (Adjuvant)",
            "Bevacizumab", "Carmustine wafer", "Optune (TTFields)", "Clinical Trial",
        ], key="match_treatments")

        if st.button("Run Trial Matching", type="primary", use_container_width=True, key="run_match"):
            profile = {
                "patient_id": match_patient_id, "age_range": match_age, "who_grade": match_grade,
                "mgmt_methylated": match_mgmt, "idh_mutant": match_idh, "egfr_amplified": match_egfr,
                "egfrviii_positive": match_egfrviii, "p53_mutant": match_p53, "pteng_loss": match_pten,
                "braf_v600e": match_braf, "prior_treatments": match_treatments,
            }

            issues = _validate_patient_profile(profile)
            if issues:
                st.warning("**Profile Validation Issues:**")
                for issue in issues:
                    st.markdown(f"  - {issue}")

            results = []
            for trial in CLINICAL_TRIALS:
                match_result = _match_trial(profile, trial)
                results.append({
                    "trial": trial,
                    "match_score": match_result["score"],
                    "reasons": match_result["reasons"],
                    "warnings": match_result["warnings"],
                })

            results.sort(key=lambda x: x["match_score"], reverse=True)

            st.markdown(f"**Match Results for {match_patient_id}:**")
            for i, r in enumerate(results[:10]):
                trial = r["trial"]
                score = r["match_score"]
                color = "var(--success)" if score >= 60 else "var(--warning)" if score >= 30 else "var(--danger)"

                with st.expander(f"#{i+1}  {trial['nct']} -- Score: {score}/100 -- {trial['drug']}"):
                    st.markdown(f"**Title:** {trial['title']}")
                    st.markdown(f"**Phase:** {trial['phase']} | **Status:** {trial['status']}")
                    st.markdown(f"**Drug:** {trial['drug']}")
                    st.markdown(f"**Eligibility:** {trial['eligibility']}")
                    st.markdown(f"**Expected Outcome:** {trial['expected_outcome']}")

                    if r["reasons"]:
                        st.markdown("**Matching Criteria:**")
                        for reason in r["reasons"]:
                            st.markdown(f"  + {reason}")
                    if r["warnings"]:
                        st.markdown("**Warnings:**")
                        for warn in r["warnings"]:
                            st.markdown(f"  - {warn}")

    # --- Molecular Subtype & Treatment Planning ---
    with tab_subtype:
        st.markdown("**GBM Molecular Subtype Classification & Treatment Planning**")

        st.markdown("**Molecular Subtypes**")
        for subtype, data in GBM_MOLECULAR_SUBTYPES.items():
            with st.expander(subtype):
                for k, v in data.items():
                    st.markdown(f"**{k.replace('_', ' ').title()}:** {v}")

        st.markdown("---")
        st.markdown("**Patient Profile & Treatment Recommendations**")

        preset = st.selectbox("Load sample patient", ["Manual input"] + [p["label"] for p in SAMPLE_PATIENTS])

        if preset != "Manual input":
            patient = next(p for p in SAMPLE_PATIENTS if p["label"] == preset)
            sp = patient["profile"]
            patient_id = sp["patient_id"]
            age_range = sp["age_range"]
            mgmt_methylated = sp["mgmt_methylated"]
            idh_mutant = sp["idh_mutant"]
            egfr_amplified = sp["egfr_amplified"]
            egfrviii_positive = sp["egfrviii_positive"]
            p53_mutant = sp["p53_mutant"]
            pteng_loss = sp["pteng_loss"]
            tumor_location = sp["tumor_location"].title()
            prior_treatments = sp["prior_treatments"]
        else:
            patient_id = st.text_input("Patient ID", value="ANON-001")
            age_range = st.selectbox("Age Range", ["18-30", "30-45", "45-60", "60-75", "75+"])
            mgmt_methylated = st.checkbox("MGMT Methylated")
            idh_mutant = st.checkbox("IDH1/2 Mutant")
            egfr_amplified = st.checkbox("EGFR Amplified")
            egfrviii_positive = st.checkbox("EGFRvIII Positive")
            p53_mutant = st.checkbox("p53 Mutant")
            pteng_loss = st.checkbox("PTEN Loss")
            tumor_location = st.selectbox("Location", ["Temporal", "Frontal", "Parietal", "Occipital", "Insular", "Brainstem"])
            prior_treatments = st.multiselect("Prior Treatments", [
                "Surgery (GTR)", "Surgery (STR)", "RT (Stupp Protocol)", "TMZ (Adjuvant)",
                "Bevacizumab", "Carmustine wafer", "Optune (TTFields)", "Clinical Trial",
            ])

        smiles_for_eval = st.text_input("Compound SMILES for viability assessment", value="CN1N=NC2=C(N=CN2C1=O)C(N)=O")

        if st.button("Generate Treatment Plan", type="primary", use_container_width=True):
            profile_data = {
                "patient_id": patient_id, "age_range": age_range, "mgmt_methylated": mgmt_methylated,
                "idh_mutant": idh_mutant, "egfr_amplified": egfr_amplified, "egfrviii_positive": egfrviii_positive,
                "p53_mutant": p53_mutant, "pteng_loss": pteng_loss, "tumor_location": tumor_location.lower(),
                "prior_treatments": prior_treatments,
            }

            subtype = _classify_subtype(profile_data)
            subtype_data = GBM_MOLECULAR_SUBTYPES.get(subtype, {})

            st.markdown(f"**Predicted Subtype:** {subtype}")
            if subtype_data:
                for k, v in subtype_data.items():
                    st.markdown(f"  - **{k.replace('_', ' ').title()}:** {v}")

            recommendations = _recommend_treatments(profile_data)
            st.markdown("**Treatment Recommendations (evidence-based):**")
            for rec in recommendations:
                evidence_color = {"Strong": "var(--success)", "Moderate": "var(--warning)", "Weak": "var(--danger)"}.get(rec["evidence"], "var(--text-secondary)")
                st.markdown(f"- **{rec['drug']}** [{rec['evidence']}] -- {rec['rationale']}")

            if smiles_for_eval.strip():
                compound_result = screen_compound(smiles_for_eval)
                if compound_result.valid:
                    ev = evaluate_lead(
                        smiles=compound_result.canonical_smiles,
                        compound_metrics=vars(compound_result),
                        patient_profile=PatientProfile(**profile_data),
                    )
                    vc = "verdict-viable" if "VIABLE" in ev.overall_verdict else "verdict-rejected" if "REJECTED" in ev.overall_verdict else "verdict-conditional"
                    st.markdown(f'<div class="verdict-box {vc}">{ev.overall_verdict} (Composite: {ev.composite_score:.1f})</div>', unsafe_allow_html=True)

                    sc1, sc2, sc3, sc4, sc5 = st.columns(5)
                    sc1.metric("Composite", f"{ev.composite_score:.1f}")
                    sc2.metric("Chemistry", f"{ev.chemistry_score:.1f}")
                    sc3.metric("Docking", f"{ev.docking_score:.1f}")
                    sc4.metric("Toxicity", f"{ev.toxicity_score:.1f}")
                    sc5.metric("Patient Match", f"{ev.patient_match_score:.1f}")

                    pdf = generate_pdf_report(ev)
                    docx = generate_docx_report(ev)
                    dc1, dc2 = st.columns(2)
                    with dc1:
                        st.download_button(T["export_pdf"], data=pdf, file_name=f"eval_{patient_id}_{datetime.now().strftime('%Y%m%d')}.pdf", mime="application/pdf", use_container_width=True)
                    with dc2:
                        if docx:
                            st.download_button(T["export_docx"], data=docx, file_name=f"eval_{patient_id}_{datetime.now().strftime('%Y%m%d')}.docx", mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document", use_container_width=True)


# ============================================================
# Main
# ============================================================
def main():
    render_header()

    tab1, tab2, tab3, tab4 = st.tabs([
        T["tab1"], T["tab2"], T["tab3"], T["tab4"],
    ])

    with tab1:
        tab_compound_screening()
    with tab2:
        tab_docking()
    with tab3:
        tab_research()
    with tab4:
        tab_anonymizer()


if __name__ == "__main__":
    main()
