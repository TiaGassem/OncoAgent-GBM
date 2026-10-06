"""OncoAgent-GBM: docking fallback for environments WITHOUT AutoDock Vina.

Why this exists
---------------
On Streamlit Cloud the Vina binary / `vina` python package is usually NOT
installed, so `run_vina_docking` returns 0.00 kcal/mol and 0 poses (exactly
the empty card the user saw). This module gives a usable result in that case
WITHOUT fabricating fake docking numbers.

It offers TWO honest sources, in priority order:

1. STORED real results -- the user's OWN previously-computed docking values
   (AutoDock4, AutoDock Vina, PLANTS) for the CDC25A/B/C + NSC-95397 system.
   These are real numbers from real runs, shown with full provenance
   (program, number of seeds/runs, and any [VERIFY] flag). This is a results
   database, not a prediction.

2. A transparent HEURISTIC affinity ESTIMATE for arbitrary compounds where no
   real run exists. It is a reproducible empirical score from basic molecular
   descriptors -- it is explicitly NOT AutoDock Vina, has NO 3D pose and NO
   receptor, and is labelled for triage/ranking only. It must be confirmed by
   a real docking run before any scientific claim.

Nothing here invents a Vina score and calls it Vina. Every output states its
source and its limits.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ------------------------------------------------------------------
# 1. STORED real docking results (user's own runs) -- provenance kept
# ------------------------------------------------------------------
# Values previously computed and verified by the user. kcal/mol for AD4 and
# Vina; PLANTS reports ChemPLP (a different, unitless score -- NEVER mix the
# two on one axis). Flags preserve the known data-quality caveats.
STORED_DOCKING = {
    "NSC-95397": {
        "CDC25A": {
            "autodock4_kcal_mol": -5.29,
            "autodock4_seeds": 10,
            "autodock4_sd": 0.00,  # suspected copy error across 10 seeds
            "vina_kcal_mol": -4.36,
            "vina_seeds": 3,
            "plants_chemplp": None,
            "flags": ["[VERIFY] AD4 SD=0.00 across 10 seeds -- check raw .dlg logs"],
        },
        "CDC25B": {
            "autodock4_kcal_mol": -6.41,
            "autodock4_seeds": 10,
            "autodock4_sd": None,
            "vina_kcal_mol": -5.99,
            "vina_seeds": 3,
            "plants_chemplp": None,
            "flags": [],
        },
        "CDC25C": {
            "autodock4_kcal_mol": -7.30,
            "autodock4_seeds": 10,
            "autodock4_sd": None,
            "vina_kcal_mol": -6.03,
            "vina_seeds": 3,
            "plants_chemplp": -69.21,
            "plants_runs": 1,
            "flags": ["PLANTS ChemPLP is unitless -- do not compare to kcal/mol"],
        },
    },
}

# Ranking observed across all three programs (computational observation only;
# ranking is NOT the same as experimental selectivity).
STORED_RANKING_NOTE = (
    "Observed computational ranking (all 3 programs): CDC25C > CDC25B > CDC25A. "
    "This is a docking-score ranking only -- it is NOT proof of biochemical "
    "selectivity, which requires enzyme-inhibition assays."
)


def normalize_target(name: str) -> str:
    """Map loose target text to CDC25A/B/C if possible."""
    t = (name or "").upper().replace(" ", "").replace("-", "")
    for iso in ("CDC25A", "CDC25B", "CDC25C"):
        if iso in t:
            return iso
    return name


def normalize_ligand(name: str) -> str:
    t = (name or "").upper().replace(" ", "").replace("_", "-")
    if "95397" in t:
        return "NSC-95397"
    return name


def lookup_stored(ligand: str, target: str) -> Optional[dict]:
    """Return the stored real result dict for a ligand/target, or None."""
    lig = normalize_ligand(ligand)
    tgt = normalize_target(target)
    return STORED_DOCKING.get(lig, {}).get(tgt)


def reference_rows(ligand: str = "NSC-95397") -> list:
    """Flat, EDITABLE seed rows for a ligand (one row per isoform x program).

    IMPORTANT: these are convenience defaults transcribed from the user's own
    thesis records. They are NOT computed by the app, and real docking scores
    depend on the exact receptor, grid, ligand preparation, and number of
    seeds -- so the user must verify/edit them against their own .dlg/.log
    files. The app only stores and organises what the user confirms.
    """
    lig = normalize_ligand(ligand)
    data = STORED_DOCKING.get(lig, {})
    rows = []
    for iso in ("CDC25A", "CDC25B", "CDC25C"):
        d = data.get(iso, {})
        if d.get("autodock4_kcal_mol") is not None:
            rows.append({"Isoform": iso, "Program": "AutoDock4",
                         "Score": d["autodock4_kcal_mol"], "Unit": "kcal/mol",
                         "Runs/Seeds": d.get("autodock4_seeds", 0)})
        if d.get("vina_kcal_mol") is not None:
            rows.append({"Isoform": iso, "Program": "AutoDock Vina",
                         "Score": d["vina_kcal_mol"], "Unit": "kcal/mol",
                         "Runs/Seeds": d.get("vina_seeds", 0)})
        if d.get("plants_chemplp") is not None:
            rows.append({"Isoform": iso, "Program": "PLANTS",
                         "Score": d["plants_chemplp"], "Unit": "ChemPLP (unitless)",
                         "Runs/Seeds": d.get("plants_runs", 0)})
    return rows


def reference_flags(ligand: str = "NSC-95397") -> list:
    """Collect all [VERIFY]/caveat flags for a ligand's stored values."""
    lig = normalize_ligand(ligand)
    data = STORED_DOCKING.get(lig, {})
    flags = []
    for iso, d in data.items():
        for f in d.get("flags", []):
            flags.append(f"{iso}: {f}")
    return flags


# ------------------------------------------------------------------
# 2. Transparent heuristic affinity ESTIMATE (NOT docking)
# ------------------------------------------------------------------
@dataclass
class HeuristicEstimate:
    estimated_affinity_kcal_mol: float = 0.0
    estimated_ki: str = ""
    rationale: str = ""
    descriptors_used: dict = field(default_factory=dict)
    is_docking: bool = False  # always False -- this is NOT a docking result
    disclaimer: str = ""


def _ki_from_energy(dg_kcal_mol: float, temp_k: float = 298.15) -> str:
    """Convert a free energy to an *order-of-magnitude* Ki string.

    Uses dG = R*T*ln(Ki). This is the standard relation; the Ki here is only
    as good as the (estimated) energy, so it is reported as approximate.
    """
    import math
    R = 1.987204259e-3  # kcal/mol/K
    try:
        ki_molar = math.exp(dg_kcal_mol / (R * temp_k))
    except OverflowError:
        return "n/a"
    # format with sensible units
    if ki_molar >= 1e-3:
        return f"~{ki_molar * 1e3:.1f} mM"
    if ki_molar >= 1e-6:
        return f"~{ki_molar * 1e6:.1f} uM"
    if ki_molar >= 1e-9:
        return f"~{ki_molar * 1e9:.1f} nM"
    return f"~{ki_molar * 1e12:.1f} pM"


def heuristic_affinity_estimate(descriptors: dict) -> HeuristicEstimate:
    """Reproducible empirical affinity ESTIMATE from molecular descriptors.

    This is a transparent, deterministic heuristic -- NOT a trained model and
    NOT a docking engine. It encodes well-known medicinal-chemistry trends
    (favourable: moderate lipophilicity and aromatic contact area; penalised:
    excessive size, polarity, and flexibility) into a simple additive score,
    anchored so typical small-molecule binders land in the -4 to -9 kcal/mol
    range that real Vina/AD4 produce. It gives a defensible *ranking* signal
    when no docking binary is available; it does not produce a 3D pose.

    `descriptors` expects keys (any missing default to neutral values):
      mw, logp, hbd, hba, tpsa, rot_bonds, aromatic_rings
    """
    d = {
        "mw": float(descriptors.get("mw", 350.0)),
        "logp": float(descriptors.get("logp", 2.5)),
        "hbd": float(descriptors.get("hbd", 1)),
        "hba": float(descriptors.get("hba", 4)),
        "tpsa": float(descriptors.get("tpsa", 60.0)),
        "rot_bonds": float(descriptors.get("rot_bonds", 4)),
        "aromatic_rings": float(descriptors.get("aromatic_rings", 2)),
    }

    # Additive empirical terms (kcal/mol). Signs follow standard SAR intuition.
    base = -4.0
    # lipophilic contact: helps up to logP~4, then plateaus
    lip = -0.45 * min(d["logp"], 4.0)
    # aromatic ring stacking / shape complementarity
    arom = -0.40 * min(d["aromatic_rings"], 4.0)
    # size: mild favourable up to ~400 Da, penalise beyond ~500 Da
    if d["mw"] <= 400:
        size = -0.0015 * d["mw"]
    else:
        size = -0.60 + 0.004 * (d["mw"] - 400)  # penalty grows past 400
    # desolvation penalty from high polarity
    polar = 0.012 * max(0.0, d["tpsa"] - 60.0)
    # entropy penalty from rotatable bonds
    flex = 0.15 * max(0.0, d["rot_bonds"] - 3)
    # excess H-bond donors slightly unfavourable for membrane/pocket fit
    hbd_pen = 0.10 * max(0.0, d["hbd"] - 3)

    dg = base + lip + arom + size + polar + flex + hbd_pen
    # clamp to the physically plausible small-molecule window
    dg = max(-12.0, min(-2.0, dg))
    dg = round(dg, 2)

    est = HeuristicEstimate(
        estimated_affinity_kcal_mol=dg,
        estimated_ki=_ki_from_energy(dg),
        descriptors_used=d,
        is_docking=False,
        rationale=(
            f"base {base:+.2f} | lipophilicity {lip:+.2f} | aromatic {arom:+.2f} "
            f"| size {size:+.2f} | polarity {polar:+.2f} | flexibility {flex:+.2f} "
            f"| HBD {hbd_pen:+.2f}  ==>  {dg:+.2f} kcal/mol"
        ),
        disclaimer=(
            "HEURISTIC ESTIMATE -- this is NOT an AutoDock Vina result. It uses "
            "molecular descriptors only: there is NO receptor, NO 3D pose, and "
            "NO scoring of actual protein-ligand contacts. Use it only to RANK "
            "candidates for triage. Any binding claim must be confirmed by a "
            "real docking run (Vina / AutoDock4 / PLANTS) and, ultimately, by "
            "an experimental enzyme-inhibition assay."
        ),
    )
    return est


def vina_unavailable_message() -> str:
    return (
        "AutoDock Vina is not installed in this environment, so a live docking "
        "run cannot be executed here. Showing stored real results where "
        "available, otherwise a clearly-labelled heuristic estimate. To run "
        "real docking, use the generated Vina command on a machine with Vina "
        "installed (see the 'Vina command' expander)."
    )
