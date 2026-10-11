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


# ------------------------------------------------------------------
# CREDIBILITY LAYER -- validated source links.
# Every evaluation/report must be traceable to these real references.
# No clinical predictions, no diagnosis; research use only.
# ------------------------------------------------------------------
SOURCE_LINKS = {
    "Clinical trials": "https://clinicaltrials.gov/",
    "Cell lines (Cellosaurus)": "https://cellosaurus.org/",
    "Cell lines (NCI DTP)": "https://dtp.cancer.gov/",
    "Toxicity / structural alerts (PubChem)": "https://pubchem.ncbi.nlm.nih.gov/",
    "Docking (SwissDock)": "https://www.swissdock.ch/",
    "Docking (CB-Dock2)": "https://cadd.labshare.cn/cb-dock2/",
    "Docking (AutoDock Vina)": "https://autodock.scripps.edu/",
    "Receptor structures (RCSB PDB)": "https://www.rcsb.org/",
    "Literature (PubMed)": "https://pubmed.ncbi.nlm.nih.gov/",
}

RESEARCH_USE_DISCLAIMER = (
    "For Research Use Only. Not medical advice and not a diagnosis. "
    "All outputs are computational pre-screens that must be verified "
    "experimentally and against the cited primary sources."
)


def get_source_links() -> dict:
    """Return the validated source links used across OncoAgent-GBM reports."""
    return dict(SOURCE_LINKS)


@dataclass
class ToxicityProfile:
    """Rule-based (heuristic) toxicity screen.

    NOTE: This is NOT ProTox-3. ProTox-3 is a machine-learning service
    (Charite, ~40,000 compound training set) with no public API. This
    module applies structural alerts and physicochemical thresholds only,
    and is intended as a fast, transparent pre-screen -- not a replacement
    for ProTox-3 or experimental toxicology.
    """
    ames_mutagenicity: str = "unknown"
    herg_inhibition: str = "unknown"
    hepatotoxicity: str = "unknown"
    skin_sensitization: str = "unknown"
    oral_toxicity: str = "unknown"
    carcinogenicity: str = "unknown"
    ld50_estimate: float = 0.0  # mg/kg
    toxicity_class: str = "unknown"
    overall_tox_risk: str = "unknown"
    # --- transparency fields ---
    tox_score: int = 0
    matched_alerts: list = field(default_factory=list)
    descriptors: dict = field(default_factory=dict)
    confidence: str = "low (rule-based)"
    method: str = "Structural alerts + LogP/MW thresholds"


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
    sub_scores: dict = field(default_factory=dict)
    chemistry_details: dict = field(default_factory=dict)
    docking_details: dict = field(default_factory=dict)
    toxicity_details: dict = field(default_factory=dict)
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
    "N=N", "N=Nc1ccccc1", "c1nn[nH]n1", "[N+](=O)[O-]",
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
    """Rule-based (heuristic) toxicity screen.

    This is a transparent pre-screen, NOT a machine-learning model and NOT
    ProTox-3. Every output is traceable to a structural alert or a
    physicochemical threshold, and all driving factors are recorded in the
    profile so results can be audited and compared against external tools.
    """
    profile = ToxicityProfile()
    if mol is None:
        profile.method = "N/A (no structure)"
        return profile

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)
    aromatic_rings = len([r for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)])

    profile.descriptors = {
        "MW (Da)": round(mw, 1),
        "LogP": round(logp, 2),
        "TPSA (A2)": round(tpsa, 1),
        "HBA": hba,
        "HBD": hbd,
        "Aromatic rings": aromatic_rings,
    }

    tox_score = 0  # higher = more toxic
    alerts = []

    # --- Mutagenicity (Ames) ---
    mut_hits = _count_substructure_matches(mol, MUTAGENIC_PATTERNS)
    if mut_hits > 0:
        profile.ames_mutagenicity = "likely_positive"
        tox_score += 2
        alerts.append(f"Mutagenicity: {mut_hits} structural alert match(es) among azo/nitro/nitroso SMARTS")
    else:
        profile.ames_mutagenicity = "likely_negative"

    # --- hERG inhibition ---
    herg_hits = _count_substructure_matches(mol, HERG_PATTERNS)
    if herg_hits > 1 or (logp > 3 and aromatic_rings >= 2):
        profile.herg_inhibition = "risk"
        tox_score += 2
        reason = f"{herg_hits} basic-amine alert(s)" if herg_hits > 1 else f"LogP {logp:.2f} > 3 with {aromatic_rings} aromatic rings"
        alerts.append(f"hERG risk: {reason}")
    else:
        profile.herg_inhibition = "low_risk"

    # --- Hepatotoxicity ---
    liver_hits = _count_substructure_matches(mol, LIVER_TOX_PATTERNS)
    if liver_hits > 2 or logp > 4:
        profile.hepatotoxicity = "risk"
        tox_score += 1
        reason = f"{liver_hits} reactive-group alert(s)" if liver_hits > 2 else f"LogP {logp:.2f} > 4"
        alerts.append(f"Hepatotoxicity risk: {reason}")
    else:
        profile.hepatotoxicity = "low_risk"

    # --- Skin sensitization ---
    if mw > 600 or logp > 5:
        profile.skin_sensitization = "moderate_risk"
        tox_score += 1
        alerts.append(f"Skin sensitization: MW {mw:.0f} > 600 or LogP {logp:.2f} > 5")

    # --- Carcinogenicity ---
    if mut_hits > 1 or (logp > 4 and mw > 500):
        profile.carcinogenicity = "alert"
        tox_score += 2
        alerts.append("Carcinogenicity: multiple mutagenic alerts or high LogP + high MW")
    else:
        profile.carcinogenicity = "low_risk"

    # LD50 / GHS class are NOT estimated: no validated model is implemented here.
    # Use ProTox-3 or experimental data (see the ProTox-3 comparison input).
    profile.ld50_estimate = 0.0
    profile.toxicity_class = "Not estimated (use ProTox-3 or experimental data)"
    profile.oral_toxicity = "Not estimated"

    # Overall flag = structural-alert burden only (a screening heuristic, not a toxicity prediction)
    if tox_score >= 6:
        profile.overall_tox_risk = "HIGH"
    elif tox_score >= 3:
        profile.overall_tox_risk = "MODERATE"
    elif tox_score >= 1:
        profile.overall_tox_risk = "LOW-MODERATE"
    else:
        profile.overall_tox_risk = "LOW"

    profile.tox_score = tox_score
    profile.matched_alerts = alerts
    if not alerts:
        profile.matched_alerts = ["No structural alerts matched"]

    return profile


def compare_with_protox(profile: ToxicityProfile, protox: dict) -> dict:
    """Compare the local heuristic screen against user-supplied ProTox-3 results.

    Parameters
    ----------
    profile : ToxicityProfile  -- output of estimate_toxicity()
    protox  : dict             -- values copied from the ProTox-3 web report.
        Recognised keys (all optional):
          ld50 (mg/kg, float), tox_class (1-6, int),
          hepatotoxicity ("Active"/"Inactive"), mutagenicity,
          carcinogenicity, neurotoxicity, bbb ("Active"/"Inactive")

    Returns
    -------
    dict with a per-endpoint agreement table and an overall concordance score.
    """
    rows = []
    agreed = 0
    compared = 0

    # LD50
    if protox.get("ld50"):
        try:
            px = float(protox["ld50"])
            rows.append({
                "Endpoint": "LD50 (mg/kg)",
                "Heuristic": "not estimated",
                "ProTox-3": f"{px:.0f}",
                "Agreement": "n/a (no local LD50 model)",
            })
        except (TypeError, ValueError):
            pass

    # Toxicity class
    if protox.get("tox_class"):
        try:
            pc = int(protox["tox_class"])
            our_c = int(profile.toxicity_class.split()[1].strip("("))
            diff = abs(pc - our_c)
            if diff == 0:
                verdict = "Exact match"
            elif diff == 1:
                verdict = "Adjacent class"
            else:
                verdict = f"Differ by {diff} classes"
            rows.append({
                "Endpoint": "Acute toxicity class",
                "Heuristic": f"Class {our_c}",
                "ProTox-3": f"Class {pc}",
                "Agreement": verdict,
            })
            compared += 1
            agreed += 1 if diff <= 1 else 0
        except (TypeError, ValueError, IndexError):
            pass

    # Binary endpoints
    binary_map = [
        ("hepatotoxicity", "Hepatotoxicity", profile.hepatotoxicity, "risk"),
        ("mutagenicity", "Mutagenicity", profile.ames_mutagenicity, "likely_positive"),
        ("carcinogenicity", "Carcinogenicity", profile.carcinogenicity, "alert"),
    ]
    for key, label, our_val, our_pos in binary_map:
        if key in protox and protox[key] not in (None, ""):
            p_active = str(protox[key]).strip().lower().startswith("active")
            our_active = (our_val == our_pos)
            match = (p_active == our_active)
            rows.append({
                "Endpoint": label,
                "Heuristic": "Positive" if our_active else "Negative",
                "ProTox-3": "Active" if p_active else "Inactive",
                "Agreement": "Yes" if match else "No",
            })
            compared += 1
            agreed += 1 if match else 0

    concordance = round(100.0 * agreed / compared, 1) if compared else 0.0

    if compared == 0:
        interpretation = "No ProTox-3 values provided for comparison."
    elif concordance >= 75:
        interpretation = "Good concordance with ProTox-3 on the supplied endpoints."
    elif concordance >= 50:
        interpretation = "Partial concordance. Divergence expected: the local screen is rule-based, ProTox-3 is machine-learning."
    else:
        interpretation = "Low concordance. Treat the heuristic as a rough pre-screen only; defer to ProTox-3 and experimental data."

    return {
        "rows": rows,
        "compared": compared,
        "agreed": agreed,
        "concordance_pct": concordance,
        "interpretation": interpretation,
    }


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


# --- Main Evaluation ---

def evaluate_lead(
    smiles: str,
    compound_metrics: Optional[dict] = None,
    docking_energy: Optional[float] = None,
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
    docking_provided = docking_energy is not None
    if docking_provided:
        result.docking_score, result.docking_details = score_docking(docking_energy)
    else:
        result.docking_score = 0.0
        result.docking_details = {"note": "No docking data provided; docking is excluded from the composite (not scored)."}

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

    # Composite = weighted mean over the components that were actually provided.
    # (Screening aid only: weights are author conventions, not validated.)
    weights = {"chemistry": 0.40, "toxicity": 0.30}
    scores = {"chemistry": result.chemistry_score, "toxicity": result.toxicity_score}
    if docking_provided:
        weights["docking"] = 0.30
        scores["docking"] = result.docking_score
    wsum = sum(weights.values())
    result.composite_score = round(sum(weights[k] * scores[k] for k in weights) / wsum, 1)

    result.sub_scores = {
        "chemistry": result.chemistry_score,
        "docking": result.docking_score if docking_provided else None,
        "toxicity": result.toxicity_score,
        "composite": result.composite_score,
    }

    # Verdict
    if result.composite_score >= 70 and result.chemistry_score >= 60:
        result.overall_verdict = "MEETS SCREENING THRESHOLDS"
    elif result.composite_score >= 50:
        result.overall_verdict = "PARTIAL - some thresholds not met"
    else:
        result.overall_verdict = "BELOW SCREENING THRESHOLDS"

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
        recs.append("Passes the in-silico screening thresholds; confirm with experiments (enzymatic and cell-based assays) before drawing conclusions")

    if not recs:
        recs.append("No specific recommendations at this time")

    return recs


# --- PDF Report Generation ---

def _latin1_safe(text) -> str:
    """Convert text to a latin-1 safe string for FPDF core fonts."""
    if text is None:
        return ""
    text = str(text)
    replacements = {
        "\u2022": "-", "\u2013": "-", "\u2014": "-", "\u2018": "'", "\u2019": "'",
        "\u201c": '"', "\u201d": '"', "\u2026": "...", "\u00a0": " ",
        "\u2192": "->", "\u2265": ">=", "\u2264": "<=", "\u00b5": "u", "\u03b1": "alpha",
        "\u03b2": "beta", "\u03b3": "gamma", "\u0394": "delta", "\u03bc": "u",
    }
    for k, v in replacements.items():
        text = text.replace(k, v)
    return text.encode("latin-1", "replace").decode("latin-1")


class AuditPDF(FPDF):
    def multi_cell(self, w, h=None, text="", *args, **kwargs):
        """Always start a wrapped block at the left margin (fixes 'Not enough
        horizontal space' when the cursor was left at the right edge by cell(new_x='END'))."""
        self.set_x(self.l_margin)
        kwargs.setdefault("new_x", "LMARGIN")
        kwargs.setdefault("new_y", "NEXT")
        return super().multi_cell(w, h, text, *args, **kwargs)

    """Custom PDF class for the lead evaluation audit report."""

    def header(self):
        self.set_font("Helvetica", "B", 16)
        self.cell(0, 10, _latin1_safe("OncoAgent-GBM Lead Evaluation Report"), align="C", new_x="LMARGIN", new_y="NEXT")
        self.set_font("Helvetica", "", 9)
        self.cell(0, 5, _latin1_safe("In-silico screening report - research use only, not a viability or clinical assessment"), align="C", new_x="LMARGIN", new_y="NEXT")
        self.ln(3)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.cell(0, 10, _latin1_safe(f"Page {self.page_no()}/{{nb}} | OncoAgent-GBM Platform | Confidential"), align="C")

    def section_title(self, title: str):
        self.set_font("Helvetica", "B", 12)
        self.set_fill_color(41, 128, 185)
        self.set_text_color(255, 255, 255)
        self.cell(0, 8, _latin1_safe(f"  {title}"), fill=True, new_x="LMARGIN", new_y="NEXT")
        self.set_text_color(0, 0, 0)
        self.ln(2)

    def key_value(self, key: str, value: str, bold_value: bool = False):
        self.set_font("Helvetica", "B", 10)
        self.cell(55, 6, _latin1_safe(f"{key}:"), new_x="END")
        style = "B" if bold_value else ""
        self.set_font("Helvetica", style, 10)
        self.cell(0, 6, _latin1_safe(value), new_x="LMARGIN", new_y="NEXT")

    def score_bar(self, label: str, score: float, max_score: float = 100):
        self.set_font("Helvetica", "", 10)
        self.cell(55, 6, _latin1_safe(f"{label}:"), new_x="END")
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
        self.cell(8, 5, "-", new_x="END")
        self.multi_cell(0, 5, _latin1_safe(text))


def _build_pdf_report(eval_result: LeadEvaluation) -> bytes:
    """Generate a comprehensive PDF audit report for the lead evaluation."""
    pdf = AuditPDF()
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_page()

    # Verdict Banner
    if "MEETS" in eval_result.overall_verdict:
        pdf.set_fill_color(46, 204, 113)
    elif "BELOW" in eval_result.overall_verdict:
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
    pdf.score_bar("Chemistry (40%)", eval_result.chemistry_score)
    pdf.score_bar("Docking (30%, if provided)", eval_result.docking_score)
    pdf.score_bar("Structural-alert screen (30%)", eval_result.toxicity_score)
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
    pdf.section_title("5. Structural-alert screen (heuristic; not a toxicity prediction)")
    tox = eval_result.toxicity_details.get("profile", {})
    if tox:
        pdf.key_value("Ames Mutagenicity", tox.get("ames", "unknown"))
        pdf.key_value("hERG Inhibition Risk", tox.get("herg", "unknown"))
        pdf.key_value("Hepatotoxicity", tox.get("hepatotoxicity", "unknown"))
        pdf.key_value("Oral Toxicity", tox.get("oral_tox", "unknown"))
        pdf.key_value("Carcinogenicity", tox.get("carcinogenicity", "unknown"))
        pdf.key_value("LD50", "not estimated (use ProTox-3 / experiment)")
        pdf.key_value("Toxicity Class", tox.get("tox_class", "unknown"))
    for k, v in eval_result.toxicity_details.items():
        if k != "profile":
            pdf.key_value(k.replace("_", " ").title(), str(v))
    pdf.ln(2)

    # Rationale
    pdf.section_title("7. Evaluation Rationale")
    pdf.set_font("Helvetica", "", 10)
    pdf.multi_cell(0, 5, _latin1_safe(eval_result.rationale))
    pdf.ln(3)

    # Recommendations
    pdf.section_title("8. Recommendations")
    for rec in eval_result.recommendations:
        pdf.add_bullet(rec)
    pdf.ln(3)

    # Validated sources (credibility layer)
    pdf.section_title("9. Validated Sources")
    pdf.set_font("Helvetica", "", 9)
    for label, url in SOURCE_LINKS.items():
        pdf.multi_cell(0, 5, _latin1_safe(f"- {label}: {url}"))
    pdf.ln(3)

    # Disclaimer
    pdf.section_title("Disclaimer")
    pdf.set_font("Helvetica", "I", 8)
    pdf.multi_cell(0, 4, _latin1_safe(
        "This report is generated by the OncoAgent-GBM automated evaluation system for "
        "research purposes only. All predictions are computational estimates and must be "
        "validated experimentally. This does not constitute medical advice or clinical "
        "decision-making guidance. The toxicity module is a transparent rule-based "
        "heuristic (structural alerts + physicochemical thresholds) and is NOT ProTox-3; "
        "results may differ from machine-learning toxicity services."
    ))

    return bytes(pdf.output())


def generate_pdf_report(eval_result: LeadEvaluation) -> bytes:
    """Safe wrapper: never raises, returns a minimal PDF on any failure."""
    try:
        return _build_pdf_report(eval_result)
    except Exception as exc:
        try:
            fallback = FPDF()
            fallback.set_auto_page_break(auto=True, margin=20)
            fallback.add_page()
            fallback.set_font("Helvetica", "B", 14)
            fallback.cell(0, 12, _latin1_safe("OncoAgent-GBM Lead Evaluation Report"),
                          align="C", new_x="LMARGIN", new_y="NEXT")
            fallback.ln(4)
            fallback.set_font("Helvetica", "", 10)
            fallback.multi_cell(0, 5, _latin1_safe(f"Verdict: {eval_result.overall_verdict}"))
            fallback.multi_cell(0, 5, _latin1_safe(f"Composite Score: {eval_result.composite_score}/100"))
            fallback.multi_cell(0, 5, _latin1_safe(f"SMILES: {eval_result.smiles}"))
            fallback.ln(2)
            fallback.multi_cell(0, 5, _latin1_safe("Note: full report rendering encountered an issue "
                                                   "and this is a reduced summary. Detail: " + str(exc)))
            return bytes(fallback.output())
        except Exception:
            return b""


def _build_docx_report(eval_result: LeadEvaluation) -> bytes:
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
    run = subtitle.add_run('In-silico screening report - research use only, not a viability or clinical assessment')
    run.font.size = Pt(11)
    run.font.color.rgb = RGBColor(100, 100, 100)

    doc.add_paragraph(f'Generated: {eval_result.timestamp[:19]}')

    if "MEETS" in eval_result.overall_verdict:
        color = RGBColor(46, 204, 113)
    elif "BELOW" in eval_result.overall_verdict:
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
        ('Chemistry (40%)', f'{eval_result.chemistry_score:.1f}'),
        ('Docking (30%, if provided)', f'{eval_result.docking_score:.1f}'),
        ('Structural-alert screen (30%)', f'{eval_result.toxicity_score:.1f}'),
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

    doc.add_heading('5. Structural-alert screen (heuristic; not a toxicity prediction)', level=1)
    tox = eval_result.toxicity_details.get("profile", {})
    if tox:
        tox_table = doc.add_table(rows=len(tox), cols=2)
        tox_table.style = 'Light Grid Accent 1'
        tox_table.alignment = WD_TABLE_ALIGNMENT.CENTER
        tox_labels = {
            'ames': 'Ames Mutagenicity', 'herg': 'hERG Inhibition Risk',
            'hepatotoxicity': 'Hepatotoxicity', 'oral_tox': 'Oral Toxicity',
            'carcinogenicity': 'Carcinogenicity', 'ld50_mg_kg': 'LD50 (mg/kg)',
            'tox_class': 'Toxicity Class',
        }
        for i, (k, v) in enumerate(tox.items()):
            tox_table.rows[i].cells[0].text = tox_labels.get(k, k)
            tox_table.rows[i].cells[1].text = str(v)

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
        'decision-making guidance. The toxicity module is a transparent rule-based '
        'heuristic and is NOT ProTox-3; results may differ from machine-learning services.'
    )
    run.font.size = Pt(8)
    run.italic = True

    # Validated sources (credibility layer)
    doc.add_heading("Validated Sources", level=1)
    for label, url in SOURCE_LINKS.items():
        doc.add_paragraph(f"{label}: {url}", style="List Bullet")

    buffer = io.BytesIO()
    doc.save(buffer)
    buffer.seek(0)
    return buffer.read()


def generate_docx_report(eval_result: LeadEvaluation) -> bytes:
    """Safe wrapper around DOCX generation; returns b'' on failure."""
    try:
        return _build_docx_report(eval_result)
    except Exception:
        return b""
