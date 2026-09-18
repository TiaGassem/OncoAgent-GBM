"""OncoAgent-GBM: Cheminformatics module for compound screening and BBB triage."""

from __future__ import annotations

import io
import json
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski

try:
    from rdkit.Chem import Draw, PandasTools
except ImportError:
    Draw = None
    PandasTools = None

try:
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    rdMolDescriptors = None

try:
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:
    rdMolDraw2D = None


@dataclass
class CompoundMetrics:
    """Computed physicochemical metrics for a single compound."""
    smiles: str
    valid: bool = False
    molecular_weight: float = 0.0
    logp: float = 0.0
    tpsa: float = 0.0
    hbd: int = 0
    hba: int = 0
    rotatable_bonds: int = 0
    num_rings: int = 0
    fraction_csp3: float = 0.0
    lipinski_violations: int = 0
    lipinski_pass: bool = False
    bbb_score: float = 0.0
    bbb_pass: bool = False
    veber_pass: bool = False
    overall_pass: bool = False
    mol_block_3d: str = ""
    canonical_smiles: str = ""
    iupac_name: str = ""
    error: str = ""


def parse_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Parse a SMILES string into an RDKit Mol object."""
    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    return mol


def generate_3d_conformer(
    mol: Chem.Mol,
    num_confs: int = 1,
    max_attempts: int = 200,
    prune_rms_thresh: float = 0.1,
    random_seed: int = 42,
) -> Optional[Chem.Mol]:
    """Generate an energy-minimized 3D conformer for a molecule.

    Uses ETKDGv3 for embedding and MMFF94s for force-field optimization.
    Returns the molecule with the best conformer, or None on failure.
    """
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = 0
    params.pruneRmsThresh = prune_rms_thresh

    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    if not conf_ids:
        return None

    best_conf_id = None
    best_energy = float("inf")

    for conf_id in conf_ids:
        try:
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            if mmff_props is None:
                continue
            AllChem.MMFFOptimizeMolecule(
                mol, confId=conf_id, mmffVariant="MMFF94s", maxIters=200
            )
            energy = AllChem.MMFFGetMoleculeForceField(
                mol, mmff_props, confId=conf_id
            ).CalcEnergy()
            if energy < best_energy:
                best_energy = energy
                best_conf_id = conf_id
        except Exception:
            continue

    if best_conf_id is None and conf_ids:
        best_conf_id = conf_ids[0]
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=best_conf_id, mmffVariant="MMFF94s")
        except Exception:
            pass

    if best_conf_id is not None:
        conf = mol.GetConformer(best_conf_id)
        new_mol = Chem.RWMol(mol)
        new_mol.RemoveAllConformers()
        new_mol.AddConformer(conf, assignId=True)
        return new_mol

    return mol


def compute_lipinski(mol: Chem.Mol) -> dict:
    """Compute Lipinski's Rule of Five properties."""
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    return {
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "lipinski_violations": violations,
        "lipinski_pass": violations <= 1,
    }


def compute_bbb_score(mol: Chem.Mol) -> dict:
    """Compute Blood-Brain Barrier penetration score using medicinal chemistry rules.

    Scoring criteria (GBM CNS drug-likeness):
      - MW < 400 Da  (+2 pts)
      - MW 400-450   (+1 pt)
      - LogP 1.5-3.5 (+3 pts)
      - LogP 1.0-1.5 or 3.5-4.5 (+1 pt)
      - TPSA < 60 Å²  (+3 pts)
      - TPSA 60-90 Å² (+2 pts)
      - TPSA 90-120 Å² (+0 pts)
      - TPSA > 120 Å² (-2 pts)
      - HBD <= 2      (+2 pts)
      - HBD 3-4       (+0 pts)
      - HBA <= 5      (+1 pt)
      - Rotatable bonds <= 5 (+1 pt)
      - Heavy atoms <= 25    (+1 pt)
      - Max score = 14
    """
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotbonds = Lipinski.NumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()

    score = 0.0

    if mw < 400:
        score += 2
    elif mw < 450:
        score += 1

    if 1.5 <= logp <= 3.5:
        score += 3
    elif 1.0 <= logp < 1.5 or 3.5 < logp <= 4.5:
        score += 1

    if tpsa < 60:
        score += 3
    elif tpsa < 90:
        score += 2
    elif tpsa < 120:
        score += 0
    else:
        score -= 2

    if hbd <= 2:
        score += 2
    elif hbd <= 4:
        score += 0
    else:
        score -= 1

    if hba <= 5:
        score += 1

    if rotbonds <= 5:
        score += 1

    if heavy_atoms <= 25:
        score += 1

    max_score = 14.0
    normalized = round(max(0, score) / max_score, 3)
    pass_threshold = 0.5

    return {
        "raw_bbb_score": round(score, 2),
        "normalized_bbb_score": normalized,
        "bbb_pass": normalized >= pass_threshold,
        "details": {
            "mw": round(mw, 2),
            "logp": round(logp, 2),
            "tpsa": round(tpsa, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rotbonds,
            "heavy_atoms": heavy_atoms,
        },
    }


def compute_veber(mol: Chem.Mol) -> dict:
    """Veber's rules for oral bioavailability: TPSA <= 140 and rotatable bonds <= 10."""
    tpsa = Descriptors.TPSA(mol)
    rotbonds = Lipinski.NumRotatableBonds(mol)
    return {
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotbonds,
        "veber_pass": tpsa <= 140 and rotbonds <= 10,
    }


def generate_2d_image(mol: Chem.Mol, width: int = 350, height: int = 300) -> bytes:
    """Generate a 2D PNG image of the molecule."""
    if rdMolDraw2D is None:
        return b""
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 1.5
    opts.padding = 0.15
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def generate_2d_svg(mol: Chem.Mol, width: int = 350, height: int = 300) -> str:
    """Generate a 2D SVG string of the molecule."""
    if rdMolDraw2D is None:
        return ""
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 1.5
    opts.padding = 0.15
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def get_mol_block_3d(mol: Chem.Mol) -> str:
    """Export the 3D conformer as a MOL block string."""
    if mol.GetNumConformers() == 0:
        return ""
    return Chem.MolToMolBlock(mol)


def compute_all_descriptors(mol: Chem.Mol) -> dict:
    """Compute a comprehensive set of descriptors."""
    return {
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "num_rings": rdMolDescriptors.RingCount(mol),
        "fraction_csp3": round(Lipinski.FractionCSP3(mol), 3),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_aromatic_rings": rdMolDescriptors.NumAromaticRings(mol),
        "num_aliphatic_rings": rdMolDescriptors.NumAliphaticRings(mol),
        "num_heteroatoms": Lipinski.NumHeteroatoms(mol),
        "num_amide_bonds": rdMolDescriptors.NumAmideBonds(mol),
        "molar_refractivity": round(Descriptors.MolMR(mol), 2),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
    }


def screen_compound(smiles: str) -> CompoundMetrics:
    """Full compound screening pipeline: parse, compute descriptors, score BBB, generate 3D."""
    result = CompoundMetrics(smiles=smiles)

    mol = parse_smiles(smiles)
    if mol is None:
        result.error = "Invalid SMILES string. Could not parse molecule."
        return result

    result.valid = True
    result.canonical_smiles = Chem.MolToSmiles(mol)

    lip = compute_lipinski(mol)
    result.molecular_weight = lip["molecular_weight"]
    result.logp = lip["logp"]
    result.tpsa = round(Descriptors.TPSA(mol), 2)
    result.hbd = lip["hbd"]
    result.hba = lip["hba"]
    result.rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    result.num_rings = rdMolDescriptors.RingCount(mol)
    result.fraction_csp3 = round(Lipinski.FractionCSP3(mol), 3)
    result.lipinski_violations = lip["lipinski_violations"]
    result.lipinski_pass = lip["lipinski_pass"]

    bbb = compute_bbb_score(mol)
    result.bbb_score = bbb["normalized_bbb_score"]
    result.bbb_pass = bbb["bbb_pass"]

    veber = compute_veber(mol)
    result.veber_pass = veber["veber_pass"]

    result.overall_pass = result.lipinski_pass and result.bbb_pass and result.veber_pass

    mol_3d = generate_3d_conformer(mol)
    if mol_3d is not None:
        result.mol_block_3d = get_mol_block_3d(mol_3d)

    return result


def batch_screen(smiles_list: list[str]) -> list[CompoundMetrics]:
    """Screen a list of SMILES and return metrics for each."""
    return [screen_compound(s) for s in smiles_list]


def export_results_json(results: list[CompoundMetrics]) -> str:
    """Export screening results as a JSON string."""
    return json.dumps([asdict(r) for r in results], indent=2)


def export_results_csv_rows(results: list[CompoundMetrics]) -> list[dict]:
    """Export screening results as a list of dicts (for CSV/DataFrame)."""
    rows = []
    for r in results:
        rows.append({
            "SMILES": r.canonical_smiles or r.smiles,
            "Valid": r.valid,
            "MW": r.molecular_weight,
            "LogP": r.logp,
            "TPSA": r.tpsa,
            "HBD": r.hbd,
            "HBA": r.hba,
            "RotBonds": r.rotatable_bonds,
            "Lipinski_Violations": r.lipinski_violations,
            "Lipinski_Pass": r.lipinski_pass,
            "BBB_Score": r.bbb_score,
            "BBB_Pass": r.bbb_pass,
            "Veber_Pass": r.veber_pass,
            "Overall_Pass": r.overall_pass,
            "Error": r.error,
        })
    return rows
