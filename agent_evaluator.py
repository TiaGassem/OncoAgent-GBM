"""OncoAgent-GBM: Lead compound evaluator with multi-criteria viability scoring, PDF and DOCX audit."""

from __future__ import annotations

import os
import io
import json
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

from fpdf import FPDF

try:
    from docx import Document
    from docx.shared import Inches, Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT
    HAS_DOCX = True
except ImportError:
    HAS_DOCX = False

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@dataclass
class ToxicityProfile:
    """ProTox-3 inspired toxicity metrics (computed locally)."""
    ames_mutagenicity: str = "unknown"
    herg_inhibition: str = "unknown"
    hepatotoxicity: str = "unknown"
    skin_sensitization: str = "unknown"
    oral_toxicity: str = "unknown"
    carcinogenicity: str = "unknown"
    ld50_estimate: float = 0.0  # mg/kg
    toxicity_class: str = "unknown"
    overall_tox_risk: str = "unknown"


@dataclass
class PatientProfile:
    """Anonymized patient mutation profile."""
    patient_id: str = "ANON-000"
    age_range: str = "50-65"
    mgmt_methylated: bool = False
    idh_mutant: bool = False
    egfr_amplified: bool = False
    egfrviii_positive: bool = False
    p53_mutant: bool = False
    pteng_loss: bool = False
    tumor_location: str = "temporal"
    who_grade: int = 4
    prior_treatments: list[str] = field(default_factory=list)


@dataclass
class LeadEvaluation:
    """Complete lead compound evaluation result."""
    smiles: str = ""
    compound_name: str = ""
    overall_verdict: str = "PENDING"
    composite_score: float = 0.0
    chemistry_score: float = 0.0
    docking_score: float = 0.0
    toxicity_score: float = 0.0
    patient_match_score: float = 0.0
    sub_scores: dict = field(default_factory=dict)
    chemistry_details: dict = field(default_factory=dict)
    docking_details: dict = field(default_factory=dict)
    toxicity_details: dict = field(default_factory=dict)
    patient_match_details: dict = field(default_factory=dict)
    rationale: str = ""
    recommendations: list[str] = field(default_factory=list)
    timestamp: str = ""
    pdf_bytes: bytes = b""


# --- Tox Estimation (ProTox-3 inspired, rule-based) ---

LIVER_TOX_PATTERNS = [
    "c=O", "O=C(c)", "c(=O)c", "NC(=S)", "c1ccccc1O",
    "O=C(O)", "C(=O)N", "c1ccc(N)cc1",
]

HERG_PATTERNS = [
    "C1CCCCN1", "c1ccc2[nH]ccc2c1", "n1ccnc1",
]

MUTAGENIC_PATTERNS = [
    "N=N", "N=Nc1", "c1nnn", "N(=O)=O", "[N+](=O)[O-]",
    "C=NN", "NNC=O",
]


def _count_substructure_matches(mol, smarts_list):
    if mol is None:
        return 0
    count = 0
    for smarts in smarts_list:
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat:
                count += len(mol.GetSubstructMatches(pat))
        except Exception:
            continue
    return count


def estimate_toxicity(mol) -> ToxicityProfile:
    """Estimate toxicity using structural alerts and physicochemical rules."""
    profile = ToxicityProfile()
    if mol is None:
        return profile

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    aromatic_rings = len([r for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)])

    tox_score = 0  # higher = more toxic

    mut_hits = _count_substructure_matches(mol, MUTAGENIC_PATTERNS)
    if mut_hits > 0:
        profile.ames_mutagenicity = "likely_positive"
        tox_score += 2
    else:
        profile.ames_mutagenicity = "likely_negative"

    herg_hits = _count_substructure_matches(mol, HERG_PATTERNS)
    if herg_hits > 1 or (logp > 3 and aromatic_rings >= 2):
        profile.herg_inhibition = "risk"
        tox_score += 2
    else:
        profile.herg_inhibition = "low_risk"

    liver_hits = _count_substructure_matches(mol, LIVER_TOX_PATTERNS)
    if liver_hits > 2 or logp > 4:
        profile.hepatotoxicity = "risk"
        tox_score += 1
    else:
        profile.hepatotoxicity = "low_risk"

    if mw > 600 or logp > 5:
        profile.skin_sensitization = "moderate_risk"
        tox_score += 1
    else:
        profile.skin_sensitization = "low_risk"

    if logp > 3.5 and mw > 400:
        profile.oral_toxicity = "Class 4 (harmful)"
    elif logp > 5:
        profile.oral_toxicity = "Class 3 (toxic)"
    else:
        profile.oral_toxicity = "Class 5 (slightly hazardous)"

    if mut_hits > 1 or (logp > 4 and mw > 500):
        profile.carcinogenicity = "alert"
        tox_score += 2
    else:
        profile.carcinogenicity = "low_risk"

    # Estimate LD50 using Haggin (logP-based heuristic)
    if logp < 1:
        profile.ld50_estimate = 5000
    elif logp < 2:
        profile.ld50_estimate = 2000
    elif logp < 3:
        profile.ld50_estimate = 500
    elif logp < 4:
        profile.ld50_estimate = 300
    else:
        profile.ld50_estimate = 100

    if tox_score >= 6:
        profile.toxicity_class = "Class 2 (dangerous)"
        profile.overall_tox_risk = "HIGH"
    elif tox_score >= 4:
        profile.toxicity_class = "Class 3 (harmful)"
        profile.overall_tox_risk = "MODERATE"
    elif tox_score >= 2:
        profile.toxicity_class = "Class 4 (caution)"
        profile.overall_tox_risk = "LOW-MODERATE"
    else:
        profile.toxicity_class = "Class 5 (low risk)"
        profile.overall_tox_risk = "LOW"

    return profile


# --- Scoring Functions ---

def score_chemistry(metrics: dict) -> tuple[float, dict]:
    """Score compound chemistry (0-100 scale)."""
    score = 100.0
    details = {}

    mw = metrics.get("molecular_weight", 0)
    logp = metrics.get("logp", 0)
    tpsa = metrics.get("tpsa", 0)
    hbd = metrics.get("hbd", 0)
    hba = metrics.get("hba", 0)
    violations = metrics.get("lipinski_violations", 0)
    bbb_score = metrics.get("bbb_score", 0)
    veber = metrics.get("veber_pass", False)
    rotatable = metrics.get("rotatable_bonds", 0)

    # Lipinski compliance (0-30 pts)
    lipinski_pts = max(0, 30 - violations * 10)
    score -= (30 - lipinski_pts)
    details["lipinski_pts"] = lipinski_pts

    # BBB compliance (0-25 pts)
    bbb_pts = bbb_score * 25
    score -= (25 - bbb_pts)
    details["bbb_pts"] = round(bbb_pts, 1)

    # TPSA (0-15 pts)
    if tpsa < 60:
        tpsa_pts = 15
    elif tpsa < 90:
        tpsa_pts = 10
    elif tpsa < 120:
        tpsa_pts = 5
    else:
        tpsa_pts = 0
    score -= (15 - tpsa_pts)
    details["tpsa_pts"] = tpsa_pts

    # Veber (0-10 pts)
    veber_pts = 10 if veber else 0
    score -= (10 - veber_pts)
    details["veber_pts"] = veber_pts

    # Rotatable bonds (0-10 pts)
    if rotatable <= 5:
        rot_pts = 10
    elif rotatable <= 10:
        rot_pts = 5
    else:
        rot_pts = 0
    score -= (10 - rot_pts)
    details["rot_pts"] = rot_pts

    # MW penalty
    if mw > 500:
        score -= 10
        details["mw_penalty"] = -10
    elif mw > 450:
        score -= 5
        details["mw_penalty"] = -5
    else:
        details["mw_penalty"] = 0

    score = max(0, min(100, score))
    return round(score, 1), details


def score_docking(binding_energy: float) -> tuple[float, dict]:
    """Score docking result (0-100 scale)."""
    details = {}
    if binding_energy >= 0:
        score = 0
    elif binding_energy <= -12:
        score = 100
    elif binding_energy <= -10:
        score = 90
    elif binding_energy <= -8:
        score = 75
    elif binding_energy <= -6:
        score = 55
    elif binding_energy <= -4:
        score = 30
    else:
        score = 10

    details["binding_energy"] = binding_energy
    details["binding_category"] = (
        "EXCELLENT" if score >= 85 else
        "GOOD" if score >= 70 else
        "MODERATE" if score >= 50 else
        "WEAK" if score >= 30 else
        "VERY_WEAK"
    )
    return score, details


def score_toxicity(profile: ToxicityProfile) -> tuple[float, dict]:
    """Score toxicity profile (0-100, higher = safer)."""
    score = 100.0
    details = {}

    tox_risk = profile.overall_tox_risk
    if tox_risk == "HIGH":
        score -= 50
    elif tox_risk == "MODERATE":
        score -= 30
    elif tox_risk == "LOW-MODERATE":
        score -= 15
    else:
        score -= 0

    if profile.ames_mutagenicity == "likely_positive":
        score -= 15
        details["ames_penalty"] = -15

    if profile.herg_inhibition == "risk":
        score -= 15
        details["herg_penalty"] = -15

    if profile.hepatotoxicity == "risk":
        score -= 10
        details["liver_penalty"] = -10

    if profile.carcinogenicity == "alert":
        score -= 10
        details["carcinogenicity_penalty"] = -10

    score = max(0, min(100, score))
    details["overall_tox_risk"] = tox_risk
    return round(score, 1), details


def match_patient(profile: PatientProfile, compound_metrics: dict) -> tuple[float, dict]:
    """Score patient-compound matching based on molecular profile."""
    score = 50.0
    details = {}

    if profile.mgmt_methylated:
        score += 10
        details["mgmt_match"] = "favorable"
    else:
        score -= 5
        details["mgmt_match"] = "unfavorable"
        details["mgmt_note"] = "Unmethylated MGMT - TMZ resistance expected"

    if profile.idh_mutant:
        score += 5
        details["idh_note"] = "IDH1 mutant - consider IDH1 inhibitor"

    if profile.egfrviii_positive:
        score += 10
        details["egfrviii_note"] = "EGFRvIII+ - consider targeted therapy"
    elif profile.egfr_amplified:
        score += 5
        details["egfr_note"] = "EGFR amplified - EGFR-targeted therapy candidate"

    if profile.pteng_loss:
        score -= 5
        details["pten_note"] = "PTEN loss - PI3K pathway activation"

    if profile.p53_mutant:
        score -= 5
        details["p53_note"] = "p53 mutant - increased genomic instability"

    bbb = compound_metrics.get("bbb_score", 0)
    if bbb > 0.6:
        score += 10
        details["bbb_match"] = "good_brain_penetration"
    elif bbb > 0.4:
        score += 5
        details["bbb_match"] = "moderate_brain_penetration"
    else:
        score -= 10
        details["bbb_match"] = "poor_brain_penetration"

    score = max(0, min(100, score))
    return round(score, 1), details


# --- Main Evaluation ---

def evaluate_lead(
    smiles: str,
    compound_metrics: Optional[dict] = None,
    docking_energy: Optional[float] = None,
    patient_profile: Optional[PatientProfile] = None,
    compound_name: str = "",
) -> LeadEvaluation:
    """Run full lead compound evaluation."""
    result = LeadEvaluation(smiles=smiles, compound_name=compound_name)
    result.timestamp = datetime.now().isoformat()

    if compound_metrics is None:
        compound_metrics = {}

    # Chemistry scoring
    result.chemistry_score, result.chemistry_details = score_chemistry(compound_metrics)

    # Docking scoring
    if docking_energy is not None:
        result.docking_score, result.docking_details = score_docking(docking_energy)
    else:
        result.docking_score = 50.0
        result.docking_details = {"note": "No docking data provided; score set to 50 (neutral)"}

    # Toxicity scoring
    if HAS_RDKIT and smiles:
        mol = Chem.MolFromSmiles(smiles)
        tox_profile = estimate_toxicity(mol)
    else:
        tox_profile = ToxicityProfile()

    result.toxicity_score, result.toxicity_details = score_toxicity(tox_profile)
    result.toxicity_details["profile"] = {
        "ames": tox_profile.ames_mutagenicity,
        "herg": tox_profile.herg_inhibition,
        "hepatotoxicity": tox_profile.hepatotoxicity,
        "oral_tox": tox_profile.oral_toxicity,
        "carcinogenicity": tox_profile.carcinogenicity,
        "ld50_mg_kg": tox_profile.ld50_estimate,
        "tox_class": tox_profile.toxicity_class,
    }

    # Patient matching
    if patient_profile:
        result.patient_match_score, result.patient_match_details = match_patient(
            patient_profile, compound_metrics
        )
    else:
        result.patient_match_score = 50.0
        result.patient_match_details = {"note": "No patient profile provided"}

    # Composite score (weighted average)
    weights = {"chemistry": 0.30, "docking": 0.30, "toxicity": 0.25, "patient": 0.15}
    result.composite_score = round(
        weights["chemistry"] * result.chemistry_score
        + weights["docking"] * result.docking_score
        + weights["toxicity"] * result.toxicity_score
        + weights["patient"] * result.patient_match_score,
        1,
    )

    result.sub_scores = {
        "chemistry": result.chemistry_score,
        "docking": result.docking_score,
        "toxicity": result.toxicity_score,
        "patient_match": result.patient_match_score,
        "composite": result.composite_score,
    }

    # Verdict
    if result.composite_score >= 70 and result.chemistry_score >= 60:
        result.overall_verdict = "VIABLE LEAD"
    elif result.composite_score >= 50:
        result.overall_verdict = "CONDITIONAL - Requires optimization"
    else:
        result.overall_verdict = "REJECTED"

    # Generate rationale
    result.rationale = _generate_rationale(result)
    result.recommendations = _generate_recommendations(result)

    return result


def _generate_rationale(result: LeadEvaluation) -> str:
    """Generate human-readable rationale for the evaluation."""
    parts = [f"Composite Score: {result.composite_score}/100"]

    if result.chemistry_score >= 70:
        parts.append("Strong drug-like properties with good BBB penetration potential")
    elif result.chemistry_score >= 50:
        parts.append("Moderate chemistry profile; structural optimization recommended")
    else:
        parts.append("Poor drug-like properties; significant redesign needed")

    if result.docking_score >= 70:
        parts.append("Favorable predicted binding to target")
    elif result.docking_score >= 50:
        parts.append("Moderate binding affinity; further optimization needed")
    else:
        parts.append("Weak predicted binding; unlikely to be effective")

    if result.toxicity_score >= 70:
        parts.append("Low toxicity risk profile")
    elif result.toxicity_score >= 50:
        parts.append("Moderate toxicity concerns; in vitro validation recommended")
    else:
        parts.append("Significant toxicity flags; structural modification required")

    return "; ".join(parts)


def _generate_recommendations(result: LeadEvaluation) -> list[str]:
    """Generate actionable recommendations."""
    recs = []

    if result.chemistry_score < 60:
        recs.append("Optimize molecular properties for BBB penetration (reduce MW, adjust LogP)")

    if result.docking_score < 60:
        recs.append("Improve target binding through structure-activity relationship (SAR) studies")

    if result.toxicity_score < 60:
        recs.append("Address toxicity concerns: evaluate mutagenic alerts and hERG liability")

    tox = result.toxicity_details.get("profile", {})
    if tox.get("ames") == "likely_positive":
        recs.append("Remove or modify mutagenic structural alerts (e.g., nitro groups, azo bonds)")

    if tox.get("herg") == "risk":
        recs.append("Reduce hERG inhibition risk: decrease lipophilicity or add polar groups")

    if result.composite_score >= 70:
        recs.append("Compound is viable for advancing to in vitro enzymatic and cell-based assays")

    if result.patient_match_score > 60:
        recs.append("Patient-specific molecular profile is favorable for this compound")

    if not recs:
        recs.append("No specific recommendations at this time")

    return recs


# --- PDF Report Generation ---

class AuditPDF(FPDF):
    """Custom PDF class for the lead evaluation audit report."""

    def header(self):
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, "OncoAgent-GBM Lead Evaluation Report", align="C", new_x="LMARGIN", new_y="NEXT")
        self.set_font("Helvetica", "", 9)
        self.cell(0, 5, "Automated Lead Viability Audit - Glioblastoma Drug Discovery", align="C", new_x="LMARGIN", new_y="NEXT")
        self.ln(3)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}} | OncoAgent-GBM Platform | Confidential", align="C")

    def section_title(self, title: str):
        self.set_font("Helvetica", "B", 12)
        self.set_fill_color(41, 128, 185)
        self.set_text_color(255, 255, 255)
        self.cell(0, 8, f"  {title}", fill=True, new_x="LMARGIN", new_y="NEXT")
        self.set_text_color(0, 0, 0)
        self.ln(2)

    def key_value(self, key: str, value: str, bold_value: bool = False):
        self.set_font("Helvetica", "B", 10)
        self.cell(55, 6, f"{key}:", new_x="END")
        style = "B" if bold_value else ""
        self.set_font("Helvetica", style, 10)
        self.cell(0, 6, str(value), new_x="LMARGIN", new_y="NEXT")

    def score_bar(self, label: str, score: float, max_score: float = 100):
        self.set_font("Helvetica", "", 10)
        self.cell(55, 6, f"{label}:", new_x="END")
        bar_width = 80
        bar_height = 5
        x = self.get_x()
        y = self.get_y()

        self.set_fill_color(220, 220, 220)
        self.rect(x, y, bar_width, bar_height, "F")

        fill = min(score / max_score, 1.0)
        if fill >= 0.7:
            self.set_fill_color(46, 204, 113)
        elif fill >= 0.5:
            self.set_fill_color(241, 196, 15)
        else:
            self.set_fill_color(231, 76, 60)

        self.rect(x, y, bar_width * fill, bar_height, "F")
        self.set_xy(x + bar_width + 3, y)
        self.set_font("Helvetica", "B", 10)
        self.cell(20, 6, f"{score:.1f}", new_x="LMARGIN", new_y="NEXT")

    def add_bullet(self, text: str):
        self.set_font("Helvetica", "", 9)
        self.cell(8, 5, chr(8226), new_x="END")
        self.multi_cell(0, 5, text)


def generate_pdf_report(eval_result: LeadEvaluation) -> bytes:
    """Generate a comprehensive PDF audit report for the lead evaluation."""
    pdf = AuditPDF()
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_page()

    # Verdict Banner
    if "VIABLE" in eval_result.overall_verdict:
        pdf.set_fill_color(46, 204, 113)
    elif "REJECTED" in eval_result.overall_verdict:
        pdf.set_fill_color(231, 76, 60)
    else:
        pdf.set_fill_color(241, 196, 15)

    pdf.set_font("Helvetica", "B", 14)
    pdf.set_text_color(255, 255, 255)
    pdf.cell(0, 12, f"VERDICT: {eval_result.overall_verdict}", align="C", fill=True,
             new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)
    pdf.ln(5)

    # Compound Information
    pdf.section_title("1. Compound Information")
    pdf.key_value("SMILES", eval_result.smiles)
    pdf.key_value("Compound Name", eval_result.compound_name or "Unnamed")
    pdf.key_value("Evaluation Date", eval_result.timestamp[:19])
    pdf.ln(3)

    # Score Summary
    pdf.section_title("2. Evaluation Scores (Weighted Composite)")
    pdf.score_bar("Chemistry (30%)", eval_result.chemistry_score)
    pdf.score_bar("Docking (30%)", eval_result.docking_score)
    pdf.score_bar("Toxicity Safety (25%)", eval_result.toxicity_score)
    pdf.score_bar("Patient Match (15%)", eval_result.patient_match_score)
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(55, 7, "COMPOSITE SCORE:", new_x="END")
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 7, f"{eval_result.composite_score:.1f} / 100", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3)

    # Chemistry Details
    pdf.section_title("3. Chemistry Profile")
    chem = eval_result.chemistry_details
    for k, v in chem.items():
        pdf.key_value(k.replace("_", " ").title(), str(v))
    pdf.ln(2)

    # Docking Details
    pdf.section_title("4. Docking Analysis")
    dock = eval_result.docking_details
    for k, v in dock.items():
        pdf.key_value(k.replace("_", " ").title(), str(v))
    pdf.ln(2)

    # Toxicity Details
    pdf.section_title("5. Toxicity Assessment (ProTox-3 Inspired)")
    tox = eval_result.toxicity_details.get("profile", {})
    if tox:
        pdf.key_value("Ames Mutagenicity", tox.get("ames", "unknown"))
        pdf.key_value("hERG Inhibition Risk", tox.get("herg", "unknown"))
        pdf.key_value("Hepatotoxicity", tox.get("hepatotoxicity", "unknown"))
        pdf.key_value("Oral Toxicity", tox.get("oral_tox", "unknown"))
        pdf.key_value("Carcinogenicity", tox.get("carcinogenicity", "unknown"))
        pdf.key_value("Est. LD50", f"{tox.get('ld50_mg_kg', 0)} mg/kg")
        pdf.key_value("Toxicity Class", tox.get("tox_class", "unknown"))
    for k, v in eval_result.toxicity_details.items():
        if k != "profile":
            pdf.key_value(k.replace("_", " ").title(), str(v))
    pdf.ln(2)

    # Patient Match
    if eval_result.patient_match_details:
        pdf.section_title("6. Patient Profile Match")
        for k, v in eval_result.patient_match_details.items():
            pdf.key_value(k.replace("_", " ").title(), str(v))
        pdf.ln(2)

    # Rationale
    pdf.section_title("7. Evaluation Rationale")
    pdf.set_font("Helvetica", "", 10)
    pdf.multi_cell(0, 5, eval_result.rationale)
    pdf.ln(3)

    # Recommendations
    pdf.section_title("8. Recommendations")
    for rec in eval_result.recommendations:
        pdf.add_bullet(rec)
    pdf.ln(3)

    # Disclaimer
    pdf.section_title("Disclaimer")
    pdf.set_font("Helvetica", "I", 8)
    pdf.multi_cell(0, 4,
        "This report is generated by the OncoAgent-GBM automated evaluation system for "
        "research purposes only. All predictions are computational estimates and must be "
        "validated experimentally. This does not constitute medical advice or clinical "
        "decision-making guidance. ProTox-3 metrics are estimated using structural alerts "
        "and physicochemical rules, not full ML models."
    )

    return pdf.output()


def generate_docx_report(eval_result: LeadEvaluation) -> bytes:
    """Generate a comprehensive Word document audit report."""
    if not HAS_DOCX:
        return b""

    doc = Document()

    style = doc.styles['Normal']
    font = style.font
    font.name = 'Calibri'
    font.size = Pt(10)

    title = doc.add_heading('OncoAgent-GBM Lead Evaluation Report', level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER

    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = subtitle.add_run('Automated Lead Viability Audit - Glioblastoma Drug Discovery')
    run.font.size = Pt(11)
    run.font.color.rgb = RGBColor(100, 100, 100)

    doc.add_paragraph(f'Generated: {eval_result.timestamp[:19]}')

    if "VIABLE" in eval_result.overall_verdict:
        color = RGBColor(46, 204, 113)
    elif "REJECTED" in eval_result.overall_verdict:
        color = RGBColor(231, 76, 60)
    else:
        color = RGBColor(241, 196, 15)

    verdict_para = doc.add_paragraph()
    verdict_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = verdict_para.add_run(f'VERDICT: {eval_result.overall_verdict}')
    run.font.size = Pt(14)
    run.bold = True
    run.font.color.rgb = color

    doc.add_heading('1. Compound Information', level=1)
    doc.add_paragraph(f'SMILES: {eval_result.smiles}')
    doc.add_paragraph(f'Compound Name: {eval_result.compound_name or "Unnamed"}')

    doc.add_heading('2. Evaluation Scores', level=1)
    table = doc.add_table(rows=6, cols=2)
    table.style = 'Light Grid Accent 1'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    scores = [
        ('Chemistry (30%)', f'{eval_result.chemistry_score:.1f}'),
        ('Docking (30%)', f'{eval_result.docking_score:.1f}'),
        ('Toxicity Safety (25%)', f'{eval_result.toxicity_score:.1f}'),
        ('Patient Match (15%)', f'{eval_result.patient_match_score:.1f}'),
        ('COMPOSITE SCORE', f'{eval_result.composite_score:.1f} / 100'),
        ('Verdict', eval_result.overall_verdict),
    ]
    for i, (label, value) in enumerate(scores):
        table.rows[i].cells[0].text = label
        table.rows[i].cells[1].text = value

    doc.add_heading('3. Chemistry Profile', level=1)
    for k, v in eval_result.chemistry_details.items():
        doc.add_paragraph(f'{k.replace("_", " ").title()}: {v}')

    doc.add_heading('4. Docking Analysis', level=1)
    for k, v in eval_result.docking_details.items():
        doc.add_paragraph(f'{k.replace("_", " ").title()}: {v}')

    doc.add_heading('5. Toxicity Assessment (ProTox-3 Inspired)', level=1)
    tox = eval_result.toxicity_details.get("profile", {})
    if tox:
        tox_table = doc.add_table(rows=len(tox), cols=2)
        tox_table.style = 'Light Grid Accent 1'
        tox_table.alignment = WD_TABLE_ALIGNMENT.CENTER
        tox_labels = {
            'ames': 'Ames Mutagenicity', 'herg': 'hERG Inhibition Risk',
            'hepatotoxicity': 'Hepatotoxicity', 'oral_tox': 'Oral Toxicity',
            'carcinogenicity': 'Carcinogenicity', 'ld50_mg_kg': 'Est. LD50 (mg/kg)',
            'tox_class': 'Toxicity Class',
        }
        for i, (k, v) in enumerate(tox.items()):
            tox_table.rows[i].cells[0].text = tox_labels.get(k, k)
            tox_table.rows[i].cells[1].text = str(v)

    if eval_result.patient_match_details:
        doc.add_heading('6. Patient Profile Match', level=1)
        for k, v in eval_result.patient_match_details.items():
            doc.add_paragraph(f'{k.replace("_", " ").title()}: {v}')

    doc.add_heading('7. Evaluation Rationale', level=1)
    doc.add_paragraph(eval_result.rationale)

    doc.add_heading('8. Recommendations', level=1)
    for rec in eval_result.recommendations:
        doc.add_paragraph(rec, style='List Bullet')

    doc.add_heading('Disclaimer', level=1)
    disclaimer = doc.add_paragraph()
    run = disclaimer.add_run(
        'This report is generated by the OncoAgent-GBM automated evaluation system for '
        'research purposes only. All predictions are computational estimates and must be '
        'validated experimentally. This does not constitute medical advice or clinical '
        'decision-making guidance. ProTox-3 metrics are estimated using structural alerts '
        'and physicochemical rules, not full ML models.'
    )
    run.font.size = Pt(8)
    run.italic = True

    buffer = io.BytesIO()
    doc.save(buffer)
    buffer.seek(0)
    return buffer.read()
