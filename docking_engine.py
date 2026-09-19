"""OncoAgent-GBM: Molecular docking engine for targeted multi-phosphatase docking."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

try:
    import meeko
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    HAS_MEEKO = True
except ImportError:
    HAS_MEEKO = False

try:
    from vina import Vina
    HAS_VINA_PYTHON = True
except ImportError:
    HAS_VINA_PYTHON = False


@dataclass
class GridBox:
    """Docking grid box definition."""
    center_x: float = 0.0
    center_y: float = 0.0
    center_z: float = 0.0
    size_x: float = 20.0
    size_y: float = 20.0
    size_z: float = 20.0


@dataclass
class DockingResult:
    """Result of a single docking run."""
    binding_affinity: float = 0.0
    num_modes: int = 0
    rmsd_lower: float = 0.0
    rmsd_upper: float = 0.0
    poses: list = field(default_factory=list)
    ligand_pdbqt: str = ""
    grid_box: Optional[GridBox] = None
    error: str = ""
    estimated_ki: str = ""
    binding_likelihood: str = ""


def fetch_pdb_from_rcsb(pdb_id: str, output_dir: str) -> str:
    """Download a PDB file from RCSB PDB database."""
    pdb_id = pdb_id.strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    import requests
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(response.text)

    return output_path


def clean_pdb(pdb_path: str, remove_water: bool = True, remove_hetero: bool = False) -> str:
    """Clean a PDB file by removing water, alternate locations, and optionally heteroatoms."""
    lines = []
    seen_alt = {}  # track alternate conformer selection

    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            record = line[:6].strip() if len(line) >= 6 else ""

            if record == "HETATM" and remove_water:
                resname = line[17:20].strip()
                if resname in ("HOH", "WAT", "SO4", "PO4", "GOL"):
                    continue

            if record in ("ATOM", "HETATM"):
                alt_loc = line[16] if len(line) > 16 else " "
                if alt_loc not in (" ", "A", ""):
                    continue
                if alt_loc == "A":
                    chain_res_atom = (line[20:22], line[22:26], line[27:30], line[76:78].strip())
                    if chain_res_atom in seen_alt:
                        continue
                    seen_alt[chain_res_atom] = True

            if record == "HETATM" and remove_hetero:
                continue

            lines.append(line)

    cleaned_path = pdb_path.replace(".pdb", "_cleaned.pdb")
    with open(cleaned_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    return cleaned_path


def extract_ligand_from_pdb(
    pdb_path: str,
    ligand_resname: Optional[str] = None,
    exclude_residues: Optional[list[str]] = None,
) -> Optional[str]:
    """Extract a co-crystallized ligand from a PDB file.

    If ligand_resname is given, extract that specific residue.
    Otherwise, extract the first HETATM residue that is not standard.
    """
    exclude = set(exclude_residues or ["HOH", "WAT", "SO4", "PO4", "GOL", "DMS", "ACT"])

    ligand_atoms = []
    target_resname = ligand_resname

    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname in exclude:
                    continue

                if line.startswith("HETATM"):
                    if target_resname is None:
                        target_resname = resname
                    if resname == target_resname:
                        ligand_atoms.append(line)

    if not ligand_atoms:
        return None

    ligand_pdb = "".join(ligand_atoms)
    return ligand_pdb


def compute_grid_from_ligand(ligand_pdb_text: str, padding: float = 10.0) -> GridBox:
    """Compute grid box coordinates from a ligand's atomic positions."""
    coords = []
    for line in ligand_pdb_text.strip().split("\n"):
        if line.startswith("HETATM") or line.startswith("ATOM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue

    if not coords:
        return GridBox()

    coords = np.array(coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)

    return GridBox(
        center_x=round(float(center[0]), 3),
        center_y=round(float(center[1]), 3),
        center_z=round(float(center[2]), 3),
        size_x=round(float(span[0]) + 2 * padding, 1),
        size_y=round(float(span[1]) + 2 * padding, 1),
        size_z=round(float(span[2]) + 2 * padding, 1),
    )


def compute_grid_from_residues(
    pdb_path: str,
    residue_numbers: list[int],
    chain_id: str = "A",
    padding: float = 10.0,
) -> GridBox:
    """Compute grid box from specified residue numbers (catalytic pocket definition)."""
    coords = []
    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_seq = int(line[22:26].strip())
                    chain = line[21]
                    if res_seq in residue_numbers and chain == chain_id:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue

    if not coords:
        return GridBox()

    coords = np.array(coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)

    return GridBox(
        center_x=round(float(center[0]), 3),
        center_y=round(float(center[1]), 3),
        center_z=round(float(center[2]), 3),
        size_x=round(float(span[0]) + 2 * padding, 1),
        size_y=round(float(span[1]) + 2 * padding, 1),
        size_z=round(float(span[2]) + 2 * padding, 1),
    )


def smiles_to_pdbqt(smiles: str, output_path: str) -> str:
    """Convert SMILES to a PDBQT file for Vina docking."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    if HAS_MEEKO:
        preparator = MoleculePreparation(
            hydrate=False,
            flexible_amides=False,
            rigiditate_bonds=False,
        )
        mol_setups = preparator.prepare(mol)
        if mol_setups:
            setup = mol_setups[0]
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_path, "w", encoding="utf-8") as f:
                    f.write(pdbqt_string)
                return output_path

    mol_pdb = Chem.MolToPdbBlock(mol)
    pdb_path = output_path.replace(".pdbqt", "_tmp.pdb")
    with open(pdb_path, "w", encoding="utf-8") as f:
        f.write(mol_pdb)

    pdbqt_path = _convert_pdb_to_pdbqt(pdb_path, output_path, is_ligand=True)
    return pdbqt_path


def pdb_to_pdbqt(pdb_path: str, output_path: str, is_receptor: bool = True) -> str:
    """Convert a PDB file to PDBQT format for Vina docking."""
    return _convert_pdb_to_pdbqt(pdb_path, output_path, is_ligand=not is_receptor)


def _convert_pdb_to_pdbqt(
    pdb_path: str, output_path: str, is_ligand: bool = False
) -> str:
    """Convert PDB to PDBQT using a simple atom-typing heuristic.

    For production use, AutoDockTools or Meeko is preferred. This provides
    a basic fallback for receptor preparation.
    """
    ATOM_CHARGES = {
        "C": 0.0, "CA": 0.0, "CB": 0.0, "CG": 0.0, "CD": 0.0, "CE": 0.0,
        "CZ": 0.0, "CH2": 0.0, "NE": -0.35, "NH": -0.35, "NH1": -0.40,
        "NH2": -0.40, "ND": -0.35, "NZ": -0.25, "N": -0.35, "O": -0.40,
        "OG": -0.40, "OH": -0.40, "OD": -0.40, "OE": -0.40, "SD": -0.10,
        "SG": -0.10, "P": 0.50, "F": -0.20, "CL": -0.15, "BR": -0.15,
        "I": -0.15, "FE": 0.30, "MG": 0.40, "CA2": 0.40, "ZN": 0.30,
    }

    VINA_TYPES = {
        "C": "C", "CA": "C", "CB": "C", "CG": "C", "CD": "C",
        "N": "NA", "NE": "NA", "NH": "NA", "NH1": "NA", "NH2": "NA",
        "ND": "NA", "NZ": "NA",
        "O": "OA", "OG": "OA", "OH": "OA", "OD": "OA", "OE": "OA",
        "S": "SA", "SD": "SA", "SG": "SA",
        "P": "P", "F": "F", "CL": "CL", "BR": "BR", "I": "I",
    }

    lines_out = []
    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]

                key = atom_name
                charge = ATOM_CHARGES.get(key, ATOM_CHARGES.get(element, 0.0))
                vina_type = VINA_TYPES.get(element.upper(), "C")
                if vina_type is None:
                    vina_type = "C"

                pdbqt_line = (
                    f"ATOM  {line[6:12]}{line[12:54]}"
                    f" {charge:+.3f}     {vina_type}  "
                )
                lines_out.append(pdbqt_line + "\n")

    with open(output_path, "w", encoding="utf-8") as f:
        f.writelines(lines_out)

    return output_path


def run_vina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    grid_box: GridBox,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: int = 3,
    cpu: int = 0,
) -> DockingResult:
    """Run AutoDock Vina docking and return results."""
    result = DockingResult(grid_box=grid_box)

    if HAS_VINA_PYTHON:
        return _run_vina_python(
            receptor_pdbqt, ligand_pdbqt, grid_box,
            exhaustiveness, num_modes, energy_range, cpu, result
        )

    return _run_vina_cli(
        receptor_pdbqt, ligand_pdbqt, grid_box,
        exhaustiveness, num_modes, energy_range, cpu, result
    )


def _run_vina_python(
    receptor_pdbqt: str, ligand_pdbqt: str, grid_box: GridBox,
    exhaustiveness: int, num_modes: int, energy_range: int,
    cpu: int, result: DockingResult,
) -> DockingResult:
    """Run docking using the Python Vina binding."""
    try:
        v = Vina(sf_name="vina")
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(
            center=[grid_box.center_x, grid_box.center_y, grid_box.center_z],
            box_size=[grid_box.size_x, grid_box.size_y, grid_box.size_z],
        )
        v.dock(
            exhaustiveness=exhaustiveness,
            n_poses=num_modes,
            energy_range=energy_range,
            cpu=cpu,
        )

        energies = v.energies()
        result.binding_affinity = round(float(energies[0][0]), 2)
        result.num_modes = len(energies)
        result.rmsd_lower = round(float(energies[0][1]), 3) if len(energies[0]) > 1 else 0
        result.rmsd_upper = round(float(energies[0][2]), 3) if len(energies[0]) > 2 else 0

        for i, e in enumerate(energies):
            result.poses.append({
                "rank": i + 1,
                "binding_energy_kcal_mol": round(float(e[0]), 2),
                "rmsd_lower": round(float(e[1]), 3) if len(e) > 1 else 0,
                "rmsd_upper": round(float(e[2]), 3) if len(e) > 2 else 0,
            })

        result.ligand_pdbqt = v.poses()
        result.estimated_ki = _energy_to_ki(result.binding_affinity)
        result.binding_likelihood = _classify_binding(result.binding_affinity)

    except Exception as e:
        result.error = f"Vina Python docking failed: {str(e)}"

    return result


def _run_vina_cli(
    receptor_pdbqt: str, ligand_pdbqt: str, grid_box: GridBox,
    exhaustiveness: int, num_modes: int, energy_range: int,
    cpu: int, result: DockingResult,
) -> DockingResult:
    """Run docking using AutoDock Vina CLI."""
    vina_path = shutil.which("vina") or shutil.which("autodock_vina")
    if vina_path is None:
        result.error = (
            "AutoDock Vina not found. Install via: pip install vina, "
            "or install AutoDock Vina binary and ensure it is on PATH."
        )
        return result

    with tempfile.TemporaryDirectory() as tmpdir:
        output_pdbqt = os.path.join(tmpdir, "docked_ligand.pdbqt")

        cmd = [
            vina_path,
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", str(grid_box.center_x),
            "--center_y", str(grid_box.center_y),
            "--center_z", str(grid_box.center_z),
            "--size_x", str(grid_box.size_x),
            "--size_y", str(grid_box.size_y),
            "--size_z", str(grid_box.size_z),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range),
            "--out", output_pdbqt,
        ]

        if cpu > 0:
            cmd.extend(["--cpu", str(cpu)])

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300
            )
            stdout = proc.stdout
            stderr = proc.stderr

            if proc.returncode != 0:
                result.error = f"Vina CLI error: {stderr}"
                return result

            for line in stdout.split("\n"):
                line = line.strip()
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            result.binding_affinity = float(parts[3])
                            result.rmsd_lower = float(parts[4])
                            result.rmsd_upper = float(parts[5])
                            result.num_modes += 1
                            result.poses.append({
                                "rank": result.num_modes,
                                "binding_energy_kcal_mol": float(parts[3]),
                                "rmsd_lower": float(parts[4]),
                                "rmsd_upper": float(parts[5]),
                            })
                        except (ValueError, IndexError):
                            continue

            if os.path.exists(output_pdbqt):
                with open(output_pdbqt, "r", encoding="utf-8") as f:
                    result.ligand_pdbqt = f.read()

            result.estimated_ki = _energy_to_ki(result.binding_affinity)
            result.binding_likelihood = _classify_binding(result.binding_affinity)

        except subprocess.TimeoutExpired:
            result.error = "Vina docking timed out after 300 seconds."
        except FileNotFoundError:
            result.error = "AutoDock Vina executable not found on PATH."

    return result


def _energy_to_ki(energy_kcal: float) -> str:
    """Convert binding energy (kcal/mol) to estimated Ki (nM) using the Boltzmann relation.

    Ki = exp(dG / (R * T)) where R = 1.987 cal/(mol*K), T = 298.15 K.
    """
    import math
    R = 1.987  # cal/(mol*K)
    T = 298.15  # K
    dg_cal = energy_kcal * 1000
    ki_mol = math.exp(dg_cal / (R * T))
    ki_nM = ki_mol * 1e9
    if ki_nM < 1:
        return f"{ki_nM:.3f} nM"
    elif ki_nM < 1000:
        return f"{ki_nM:.1f} nM"
    elif ki_nM < 1e6:
        return f"{ki_nM/1000:.1f} uM"
    else:
        return f"{ki_nM/1e6:.1f} mM"


def _classify_binding(energy_kcal: float) -> str:
    """Classify binding likelihood from docking score."""
    if energy_kcal <= -10.0:
        return "EXCELLENT - Strong binding predicted"
    elif energy_kcal <= -8.0:
        return "GOOD - Favorable binding"
    elif energy_kcal <= -6.0:
        return "MODERATE - Weak binding possible"
    elif energy_kcal <= -4.0:
        return "WEAK - Unlikely effective"
    else:
        return "VERY WEAK - Not recommended"


# Phosphatase target database for GBM
PHOSPHATASE_TARGETS = {
    "PTPN1": {
        "name": "PTP1B (Protein Tyrosine Phosphatase 1B)",
        "pdb_ids": ["2QBP", "3CWE", "4I8N"],
        "relevance": "Negative regulator of insulin/IGF-1 signaling; overexpressed in GBM",
        "catalytic_residues": [47, 48, 120, 121, 128, 179, 180, 214, 215, 220, 221, 228],
    },
    "PTPN2": {
        "name": "TC-PTP (T-cell Protein Tyrosine Phosphatase)",
        "pdb_ids": ["1XEO", "2HNP", "3S97"],
        "relevance": "Regulates JAK/STAT signaling; implicated in GBM immune evasion",
        "catalytic_residues": [12, 13, 48, 49, 78, 110, 111, 115, 174, 175],
    },
    "PTPN11": {
        "name": "SHP-2 (SH2 domain-containing PTP)",
        "pdb_ids": ["2SHP", "3B07", "4HJO"],
        "relevance": "Oncogenic driver in GBM; promotes RAS/MAPK and PI3K/AKT signaling",
        "catalytic_residues": [457, 458, 460, 479, 496, 497, 502, 504],
    },
    "PTPRT": {
        "name": "RPTP-kappa (Receptor-type PTP kappa)",
        "pdb_ids": ["3F7V"],
        "relevance": "Tumor suppressor lost in GBM; regulates neuronal signaling",
        "catalytic_residues": [523, 524, 526, 544, 559, 560, 565, 567],
    },
    "DUSP6": {
        "name": "DUSP6/MKP-3 (Dual-specificity PTP 6)",
        "pdb_ids": ["1MKP", "2VSW", "3O3P"],
        "relevance": "Negative regulator of ERK1/2; dysregulated in GBM",
        "catalytic_residues": [61, 125, 172, 176, 255, 257, 258, 295, 297],
    },
}

GBM_ONCOGENE_GUIDE = {
    "EGFRvIII": {
        "full_name": "Epidermal Growth Factor Receptor variant III",
        "mutation": "Deletion of exons 2-7, resulting in constitutive activation",
        "prevalence": "~50% of GBM tumors",
        "pathways": ["RAS/MAPK", "PI3K/AKT/mTOR", "JAK/STAT"],
        "targeting_strategies": [
            "Anti-EGFR monoclonal antibodies (cetuximab, nimotuzumab)",
            "TKIs (erlotinib, gefitinib - limited efficacy in GBM)",
            "EGFRvIII-directed CAR-T cells",
            "Bispecific antibodies targeting EGFRvIII/CD3",
        ],
    },
    "MGMT": {
        "full_name": "O-6-methylguanine-DNA methyltransferase",
        "mutation": "Promoter methylation (silencing) predicts temozolomide response",
        "prevalence": "35-45% of GBM patients show promoter methylation",
        "pathways": ["DNA repair", "Alkylation damage response"],
        "clinical_significance": "MGMT methylation = better prognosis with TMZ; key biomarker",
    },
    "IDH1": {
        "full_name": "Isocitrate Dehydrogenase 1",
        "mutation": "R132H hotspot mutation produces 2-hydroxyglutarate (2-HG)",
        "prevalence": "~8-15% of GBM (higher in secondary GBM)",
        "pathways": ["Epigenetic regulation via 2-HG", "HIF-1alpha stabilization", "Oncometabolite production"],
        "targeting_strategies": [
            "IDH1 inhibitors (ivosidenib, AG-120 - approved for AML)",
            "IDH1/2 inhibitors (enasidenib for IDH2)",
            "Targeting 2-HG downstream effects",
        ],
    },
    "PTPN_phosphatases": {
        "full_name": "Protein Tyrosine Phosphatases in GBM",
        "members": {
            "PTPN1": "PTP1B - regulates insulin/IGF-1, overexpressed in GBM",
            "PTPN2": "TC-PTP - JAK/STAT regulator, immune evasion role",
            "PTPN11": "SHP-2 - oncogenic, drives RAS/MAPK and PI3K/AKT",
            "PTPRT": "RPTP-kappa - tumor suppressor lost in GBM",
            "DUSP6": "MKP-3 - ERK negative regulator, dysregulated",
        },
        "therapeutic_potential": [
            "PTPN1 inhibitors (MS407, trodusquemine) show anti-GBM activity",
            "SHP-099 (SHP-2 allosteric inhibitor) blocks RAS signaling",
            "DUSP6 inhibition synergizes with MEK inhibitors",
            "Dual-specificity phosphatases as synthetic lethal targets",
        ],
    },
}


def get_phosphatase_info(target_key: str) -> Optional[dict]:
    """Retrieve phosphatase target information."""
    return PHOSPHATASE_TARGETS.get(target_key.upper())


def list_phosphatase_targets() -> dict:
    """List all available phosphatase targets."""
    return {k: {"name": v["name"], "relevance": v["relevance"]} for k, v in PHOSPHATASE_TARGETS.items()}


def get_gbm_oncogene_guide() -> dict:
    """Retrieve the GBM oncogene reference guide."""
    return GBM_ONCOGENE_GUIDE
