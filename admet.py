"""ADMET helpers: physicochemical descriptors + a BOILED-Egg-style BBB / HIA
classification, computed from the user's own SMILES with RDKit.

HONESTY NOTE: the authoritative BOILED-Egg (Daina & Zoete, ChemMedChem 2016,
doi:10.1002/cmdc.201600182) uses two confidence ELLIPSES in the WLOGP vs TPSA
plane. To avoid reproducing ellipse coefficients we are not 100%% sure of, this
module uses the simplified RECTANGULAR boundaries (the operational thresholds
requested for this project) and labels them as an approximation. For the exact
published egg, the UI links out to SwissADME. Nothing here is fabricated: every
number is computed from the molecule the user typed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    HAS_RDKIT = True
except Exception:  # pragma: no cover
    HAS_RDKIT = False


@dataclass
class AdmetResult:
    ok: bool = False
    error: str = ""
    smiles: str = ""
    mw: float = 0.0
    wlogp: float = 0.0        # Crippen MolLogP (WLOGP-type)
    tpsa: float = 0.0
    hbd: int = 0
    hba: int = 0
    rot_bonds: int = 0
    bbb_region: str = ""       # "BBB-permeant (yolk)" / "HIA-only (white)" / "outside"
    hia_region: bool = False
    bbb_region_flag: bool = False
    note: str = ""


# Simplified BOILED-Egg rectangular boundaries (approximation; see module note).
# White = passive human intestinal absorption (HIA) zone.
# Yolk  = blood-brain-barrier (BBB) permeation zone.
HIA_TPSA_MAX = 142.0
HIA_WLOGP_MIN, HIA_WLOGP_MAX = -1.0, 6.0
BBB_TPSA_MAX = 75.0
BBB_WLOGP_MIN, BBB_WLOGP_MAX = 0.5, 3.5


def compute_admet(smiles: str) -> AdmetResult:
    r = AdmetResult(smiles=smiles)
    if not HAS_RDKIT:
        r.error = "RDKit not installed on the server (add rdkit to requirements.txt)."
        return r
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        r.error = f"Invalid SMILES: {smiles}"
        return r
    r.mw = round(float(Descriptors.MolWt(mol)), 2)
    r.wlogp = round(float(Crippen.MolLogP(mol)), 2)
    r.tpsa = round(float(Descriptors.TPSA(mol)), 2)
    r.hbd = int(Lipinski.NumHDonors(mol))
    r.hba = int(Lipinski.NumHAcceptors(mol))
    r.rot_bonds = int(Descriptors.NumRotatableBonds(mol))

    r.hia_region = (r.tpsa <= HIA_TPSA_MAX and HIA_WLOGP_MIN <= r.wlogp <= HIA_WLOGP_MAX)
    r.bbb_region_flag = (r.tpsa <= BBB_TPSA_MAX and BBB_WLOGP_MIN <= r.wlogp <= BBB_WLOGP_MAX)
    if r.bbb_region_flag:
        r.bbb_region = "BBB-permeant (yolk)"
        r.note = ("Within the simplified BBB window (TPSA < 75, "
                  "0.5 < WLOGP < 3.5): brain penetration is PLAUSIBLE. Confirm on "
                  "SwissADME (exact ellipse) and experimentally.")
    elif r.hia_region:
        r.bbb_region = "HIA-only (white)"
        r.note = ("Likely intestinal absorption but OUTSIDE the BBB window: "
                  "brain penetration less likely. Confirm on SwissADME.")
    else:
        r.bbb_region = "outside"
        r.note = ("Outside both the HIA and BBB windows in this simplified map.")
    r.ok = True
    return r


def boiled_egg_png(result: AdmetResult) -> Optional[bytes]:
    """Render a BOILED-Egg-style WLOGP vs TPSA map with the molecule plotted.
    Returns PNG bytes, or None if matplotlib is unavailable.
    """
    try:
        import io
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except Exception:
        return None

    fig, ax = plt.subplots(figsize=(6, 4.2))
    # White (HIA) zone
    ax.add_patch(Rectangle((HIA_WLOGP_MIN, 0), HIA_WLOGP_MAX - HIA_WLOGP_MIN,
                           HIA_TPSA_MAX, facecolor="#f5f5f5", edgecolor="#bbbbbb",
                           label="HIA zone (white)", zorder=1))
    # Yolk (BBB) zone
    ax.add_patch(Rectangle((BBB_WLOGP_MIN, 0), BBB_WLOGP_MAX - BBB_WLOGP_MIN,
                           BBB_TPSA_MAX, facecolor="#ffd54a", edgecolor="#e0a800",
                           alpha=0.8, label="BBB zone (yolk)", zorder=2))
    ax.scatter([result.wlogp], [result.tpsa], s=90, color="#c0392b",
               edgecolor="black", zorder=5, label="Your molecule")
    ax.annotate(f"({result.wlogp}, {result.tpsa})",
                (result.wlogp, result.tpsa),
                textcoords="offset points", xytext=(8, 6), fontsize=8)
    ax.set_xlabel("WLOGP (lipophilicity)")
    ax.set_ylabel("TPSA (\u00c5\u00b2)")
    ax.set_title("BOILED-Egg (simplified) \u2013 WLOGP vs TPSA")
    ax.set_xlim(min(-2, result.wlogp - 1), max(8, result.wlogp + 1))
    ax.set_ylim(0, max(160, result.tpsa + 10))
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130)
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()
