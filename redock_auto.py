"""Automatic redocking validation (the standard credibility test for a docking set-up).

Input: a PDB file of a protein crystallised WITH a ligand. The app (1) extracts the ligand,
(2) removes it from the receptor, (3) builds the Vina box around the crystal ligand,
(4) docks the ligand back with several random seeds, (5) measures heavy-atom RMSD to the crystal pose.
Pass criterion (community convention): top-ranked pose RMSD <= 2.0 A.
Limits: no symmetry correction; ligand bond orders are perceived by OpenBabel from 3D geometry;
a pass on one structure does not validate the protocol for a different target.
"""
from __future__ import annotations
import os, shutil, subprocess, tempfile
import numpy as np
import docking_engine as de
from validation import redock_rmsd
import structure_tools as stt

IGNORE = {"HOH", "WAT", "DOD", "NA", "CL", "K", "MG", "ZN", "CA", "MN", "FE", "CU", "CO", "NI", "CD",
          "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "ACT", "DMS", "MPD", "FMT", "TRS", "BME", "NO3", "IOD", "BR"}


def list_ligands(pdb_text: str, min_heavy: int = 6) -> list:
    groups = {}
    for ln in pdb_text.splitlines():
        if ln.startswith("HETATM") and len(ln) >= 54:
            rn = ln[17:20].strip()
            if rn in IGNORE:
                continue
            el = (ln[76:78].strip() if len(ln) >= 78 else ln[12:16].strip()[:1]).upper()
            if el in ("H", "D"):
                continue
            key = (rn, ln[21:22].strip(), int(ln[22:26]))
            groups.setdefault(key, []).append(ln)
    return [{"resname": k[0], "chain": k[1], "resseq": k[2], "heavy_atoms": len(v)}
            for k, v in sorted(groups.items()) if len(v) >= min_heavy]


def _extract(pdb_text, resname, chain, resseq):
    lig, prot = [], []
    for ln in pdb_text.splitlines():
        if ln.startswith("ATOM") and len(ln) >= 54:
            prot.append(ln)
        elif ln.startswith("HETATM") and ln[17:20].strip() == resname and ln[21:22].strip() == chain \
                and int(ln[22:26]) == resseq and (ln[76:78].strip().upper() if len(ln) >= 78 else "C") not in ("H", "D"):
            lig.append(ln)
    if not lig:
        raise ValueError(f"ligand {resname} {chain}{resseq} not found")
    return lig, prot


def run_redock(pdb_text: str, resname: str, chain: str, resseq: int, seeds=(42, 7, 11),
               exhaustiveness: int = 8, pad: float = 8.0, num_modes: int = 9) -> dict:
    obabel = shutil.which("obabel")
    if not obabel:
        raise RuntimeError("OpenBabel ('obabel') not found on this server; real docking cannot run.")
    if not de.engine_status()["can_dock"]:
        raise RuntimeError("AutoDock Vina not available on this server.")
    lig, prot = _extract(pdb_text, resname, chain, resseq)
    with tempfile.TemporaryDirectory() as td:
        rec_pdb, lig_pdb = f"{td}/rec.pdb", f"{td}/lig.pdb"
        open(rec_pdb, "w").write("\n".join(prot) + "\nEND\n")
        open(lig_pdb, "w").write("\n".join(lig) + "\nEND\n")
        rec_clean = de.clean_pdb(rec_pdb, remove_hetero=True)
        de.pdb_to_pdbqt(rec_clean, f"{td}/receptor.pdbqt", is_receptor=True)
        p = subprocess.run([obabel, lig_pdb, "-O", f"{td}/lig.pdbqt", "-h", "--partialcharge", "gasteiger"],
                           capture_output=True, text=True, timeout=120)
        if not os.path.exists(f"{td}/lig.pdbqt") or os.path.getsize(f"{td}/lig.pdbqt") == 0:
            raise RuntimeError(f"ligand preparation failed: {p.stderr[-200:]}")
        ref_pdbqt = open(f"{td}/lig.pdbqt").read()
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lig])
        c, ext = xyz.mean(0), xyz.max(0) - xyz.min(0)
        size = np.clip(ext + 2 * pad, 20.0, 40.0)
        grid = de.GridBox(center_x=float(c[0]), center_y=float(c[1]), center_z=float(c[2]),
                          size_x=float(size[0]), size_y=float(size[1]), size_z=float(size[2]))
        runs, first_poses = [], None
        for sd in seeds:
            r = de.run_vina_docking(f"{td}/receptor.pdbqt", f"{td}/lig.pdbqt", grid,
                                    exhaustiveness=exhaustiveness, num_modes=num_modes, seed=sd)
            if r.error:
                raise RuntimeError(f"docking failed (seed {sd}): {r.error}")
            blocks = [b + "ENDMDL" for b in r.ligand_pdbqt.split("ENDMDL") if "ATOM" in b]
            rm = [redock_rmsd(ref_pdbqt, b) for b in blocks]
            if any(not x.ok for x in rm):
                raise RuntimeError(rm[[not x.ok for x in rm].index(True)].error)
            vals = [x.rmsd for x in rm]
            runs.append({"seed": sd, "top_pose_energy": r.binding_affinity, "top_pose_rmsd": vals[0],
                         "best_pose_rmsd": min(vals), "best_pose_rank": int(np.argmin(vals)) + 1})
            if first_poses is None:
                first_poses = r.ligand_pdbqt
        rec_text = open(rec_clean).read()
    top = [x["top_pose_rmsd"] for x in runs]
    ok = sum(v <= 2.0 for v in top)
    return {"ligand": f"{resname} {chain}{resseq}", "n_heavy": len(lig), "runs": runs,
            "n_pass": ok, "n_seeds": len(runs), "mean_top_rmsd": float(np.mean(top)),
            "grid": {"center": {"x": grid.center_x, "y": grid.center_y, "z": grid.center_z},
                     "size": {"x": grid.size_x, "y": grid.size_y, "z": grid.size_z}},
            "receptor_pdb": rec_text, "poses_pdbqt": first_poses,
            "reference_pdb": "\n".join(lig) + "\nEND\n",
            "verdict": ("PASS" if ok == len(runs) else "PARTIAL" if ok else "FAIL")}
