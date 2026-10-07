"""OncoAgent-GBM: docking validation utilities (redocking RMSD).

The single most important credibility test for any docking setup is REDOCKING:
take a protein that was solved WITH its ligand (a co-crystal structure), remove
the ligand, dock it back, and measure how far the top predicted pose is from the
real crystallographic pose (heavy-atom RMSD, in Angstrom).

Accepted community convention:
    RMSD <= 2.0 A  -> the docking protocol REPRODUCES the native pose (good).
    2.0 < RMSD <= 3.0 A -> borderline.
    RMSD > 3.0 A  -> the protocol does NOT reproduce the native pose for this case.

This module only MEASURES. It never fabricates a pose or a score. If atom counts
do not match (different ligand prep / hydrogens), it says so instead of guessing.

NOTE on correctness: this computes a direct heavy-atom RMSD assuming the atom
order is preserved between reference and docked pose (true for redocking the SAME
molecule prepared the same way). It does NOT do symmetry-corrected / Hungarian
atom matching, so for highly symmetric ligands the value can be an upper bound.
This limitation is stated to the user, not hidden.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple
import math


@dataclass
class RMSDResult:
    ok: bool = False
    error: str = ""
    rmsd: float = 0.0
    n_atoms_ref: int = 0
    n_atoms_pose: int = 0
    verdict: str = ""
    note: str = ""


def _parse_heavy_atom_xyz(text: str) -> List[Tuple[str, float, float, float]]:
    """Parse ATOM/HETATM coordinates from a PDB or PDBQT block.

    Returns a list of (element_or_name, x, y, z) for HEAVY atoms only
    (hydrogens are skipped so prep differences in H do not corrupt the RMSD).
    Only the FIRST model/pose is read (stops at the first ENDMDL after atoms).
    """
    atoms: List[Tuple[str, float, float, float]] = []
    seen_atom = False
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            seen_atom = True
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            name = line[12:16].strip()
            element = line[76:78].strip() if len(line) >= 78 else ""
            tag = (element or name).upper()
            if tag.startswith("H") and tag not in ("HG", "HF", "HO", "HE", "HS"):
                # skip hydrogens (names like H, H1, HA...). Keep rare metals.
                if element == "H" or (not element and name[:1] == "H"):
                    continue
            atoms.append((tag, x, y, z))
        elif rec == "ENDMDL" and seen_atom:
            break
    return atoms


def redock_rmsd(reference_block: str, docked_pose_block: str) -> RMSDResult:
    """Compute heavy-atom RMSD (A) between a reference (native) ligand pose and a
    docked pose. Both inputs are PDB or PDBQT text."""
    r = RMSDResult()
    ref = _parse_heavy_atom_xyz(reference_block)
    pose = _parse_heavy_atom_xyz(docked_pose_block)
    r.n_atoms_ref = len(ref)
    r.n_atoms_pose = len(pose)
    if not ref or not pose:
        r.error = ("Could not read heavy-atom coordinates from one of the inputs. "
                   "Make sure both are valid PDB/PDBQT ligand blocks.")
        return r
    if len(ref) != len(pose):
        r.error = (f"Heavy-atom counts differ (reference {len(ref)} vs docked "
                   f"{len(pose)}). Redocking RMSD requires the SAME molecule "
                   "prepared the same way. Re-extract the native ligand and dock "
                   "that exact molecule; do not compare two different ligands.")
        return r
    sq = 0.0
    for (_, ax, ay, az), (_, bx, by, bz) in zip(ref, pose):
        sq += (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
    r.rmsd = round(math.sqrt(sq / len(ref)), 3)
    if r.rmsd <= 2.0:
        r.verdict = "PASS \u2014 reproduces the native pose (RMSD \u2264 2.0 \u00c5)"
    elif r.rmsd <= 3.0:
        r.verdict = "BORDERLINE (2.0\u20133.0 \u00c5) \u2014 interpret with caution"
    else:
        r.verdict = "FAIL \u2014 does NOT reproduce the native pose (RMSD > 3.0 \u00c5)"
    r.note = ("Direct heavy-atom RMSD assuming preserved atom order (no symmetry "
              "correction). For symmetric ligands treat this as an upper bound.")
    r.ok = True
    return r
