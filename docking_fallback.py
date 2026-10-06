"""OncoAgent-GBM: docking fallback for environments WITHOUT a live engine.

Why this exists
---------------
If the AutoDock Vina binary is not on PATH (e.g. `packages.txt` missing
`autodock-vina`), or a docking run fails/times out, `run_vina_docking`
returns no usable pose. Instead of showing a blank card OR fabricating a
Vina score, this module provides ONE honest, input-dependent fallback:

  A transparent HEURISTIC affinity ESTIMATE computed only from the molecular
  descriptors of the ligand the user typed. It is explicitly NOT AutoDock
  Vina: there is NO receptor, NO 3D pose, and NO scoring of real contacts.
  It is labelled for triage/ranking only and must be confirmed by a real
  docking run before any scientific claim.

There are NO stored, hardcoded, or user-specific docking numbers anywhere in
this module. Every value shown is derived live from the current user's own
input, so nobody's results are baked in.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# ------------------------------------------------------------------
# Transparent heuristic affinity ESTIMATE (NOT docking)
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
    lip = -0.45 * min(d["logp"], 4.0)
    arom = -0.40 * min(d["aromatic_rings"], 4.0)
    if d["mw"] <= 400:
        size = -0.0015 * d["mw"]
    else:
        size = -0.60 + 0.004 * (d["mw"] - 400)
    polar = 0.012 * max(0.0, d["tpsa"] - 60.0)
    flex = 0.15 * max(0.0, d["rot_bonds"] - 3)
    hbd_pen = 0.10 * max(0.0, d["hbd"] - 3)

    dg = base + lip + arom + size + polar + flex + hbd_pen
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
        "run cannot be executed here. Showing a clearly-labelled heuristic "
        "estimate derived from your ligand only. To run real docking, add "
        "`autodock-vina` + `openbabel` to packages.txt, or use the generated "
        "Vina command on a machine with Vina installed."
    )
