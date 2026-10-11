"""OncoAgent-GBM: simple structure analysis on user-supplied PDB text (NumPy only).
 - ligand-residue contacts and H-bond *candidates* (distance only; no angle check)
 - AlphaFold pLDDT summary (B-factor column of AlphaFold models)"""
from __future__ import annotations
import numpy as np


def parse_pdb(text: str) -> list:
    atoms = []
    for ln in text.splitlines():
        if ln[:6].strip() in ("ATOM", "HETATM") and len(ln) >= 54:
            try:
                atoms.append({
                    "rec": ln[:6].strip(), "name": ln[12:16].strip(), "resname": ln[17:20].strip(),
                    "chain": ln[21:22].strip(), "resseq": int(ln[22:26]),
                    "xyz": (float(ln[30:38]), float(ln[38:46]), float(ln[46:54])),
                    "bfac": float(ln[60:66]) if len(ln) >= 66 and ln[60:66].strip() else float("nan"),
                    "elem": (ln[76:78].strip() if len(ln) >= 78 else ln[12:16].strip()[:1]).upper(),
                })
            except ValueError:
                continue
    return atoms


def _is_h(a):
    return a["elem"] in ("H", "D") or (not a["elem"] and a["name"][:1] == "H")


def ligand_contacts(protein_pdb: str, ligand_pdb: str, cutoff: float = 4.0, hbond_cutoff: float = 3.5):
    """Residues with any heavy atom within `cutoff` A of any ligand heavy atom, plus N/O pairs
    within `hbond_cutoff` (candidates only)."""
    prot = [a for a in parse_pdb(protein_pdb) if a["rec"] == "ATOM" and not _is_h(a)]
    lig = [a for a in parse_pdb(ligand_pdb) if not _is_h(a)]
    if not prot or not lig:
        raise ValueError("could not read heavy atoms from protein and/or ligand")
    P = np.array([a["xyz"] for a in prot]); L = np.array([a["xyz"] for a in lig])
    D = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2)
    contacts, hb = {}, []
    for i, a in enumerate(prot):
        j = int(D[i].argmin())
        if D[i, j] <= cutoff:
            key = (a["chain"], a["resseq"], a["resname"])
            contacts[key] = min(contacts.get(key, 99.0), float(D[i, j]))
        for j2 in np.where(D[i] <= hbond_cutoff)[0]:
            if a["elem"] in ("N", "O") and lig[j2]["elem"] in ("N", "O"):
                hb.append({"residue": f"{a['resname']}{a['resseq']}{a['chain']}", "protein_atom": a["name"],
                           "ligand_atom": lig[j2]["name"] or lig[j2]["elem"], "distance_A": round(float(D[i, j2]), 2)})
    res = [{"chain": c, "resseq": r, "resname": n, "min_distance_A": round(d, 2)}
           for (c, r, n), d in sorted(contacts.items(), key=lambda kv: (kv[0][0], kv[0][1]))]
    return res, hb


def plddt_summary(pdb_text: str) -> dict:
    """AlphaFold models store pLDDT in the B-factor column. Bands: >90 very high, 70-90 confident,
    50-70 low, <50 very low. Only meaningful for AlphaFold-style models, NOT experimental structures."""
    ca = [a["bfac"] for a in parse_pdb(pdb_text) if a["name"] == "CA" and a["rec"] == "ATOM"]
    ca = np.array([b for b in ca if np.isfinite(b)])
    if len(ca) == 0:
        raise ValueError("no CA atoms with B-factor values")
    return {"n_residues": int(len(ca)), "mean_plddt": round(float(ca.mean()), 1),
            "very_high_>90": round(float((ca > 90).mean()), 3),
            "confident_70_90": round(float(((ca > 70) & (ca <= 90)).mean()), 3),
            "low_50_70": round(float(((ca > 50) & (ca <= 70)).mean()), 3),
            "very_low_<=50": round(float((ca <= 50).mean()), 3)}
