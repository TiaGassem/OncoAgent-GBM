"""Formula-verified reference compounds.

Every SMILES here is checked in ``tests/test_reference_compounds.py`` against an
expected molecular formula (RDKit). A structure that does not match is a bug.
These are convenience examples for the screening tab, NOT a curated database.
For any other compound, fetch the structure live (see ``data_sources.chembl_search``)
or paste your own SMILES, and verify it yourself.
"""
from __future__ import annotations

# name -> (SMILES, expected molecular formula)
_VERIFIED: dict[str, tuple[str, str]] = {
    "Temozolomide": ("CN1N=NC2=C(N=CN2C1=O)C(N)=O", "C6H6N6O2"),
    "Lomustine (CCNU)": ("ClCCN(N=O)C(=O)NC1CCCCC1", "C9H16ClN3O2"),
    "Carmustine (BCNU)": ("ClCCNC(=O)N(CCCl)N=O", "C5H9Cl2N3O2"),
    "Dacarbazine": ("CN(C)N=Nc1[nH]cnc1C(N)=O", "C6H10N6O"),
    "Procarbazine": ("CNNCc1ccc(cc1)C(=O)NC(C)C", "C12H19N3O"),
    "Busulfan": ("CS(=O)(=O)OCCCCOS(C)(=O)=O", "C6H14O6S2"),
    "Chlorambucil": ("OC(=O)CCCc1ccc(cc1)N(CCCl)CCCl", "C14H19Cl2NO2"),
    "Bendamustine": ("Cn1c(CCCC(O)=O)nc2cc(ccc12)N(CCCl)CCCl", "C16H21Cl2N3O2"),
    "Thalidomide": ("O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1", "C13H10N2O4"),
    "Vorinostat (SAHA)": ("ONC(=O)CCCCCCC(=O)Nc1ccccc1", "C14H20N2O3"),
    "Belinostat": ("ONC(=O)/C=C/c1cccc(c1)S(=O)(=O)Nc1ccccc1", "C15H14N2O4S"),
    "Erlotinib": ("COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC", "C22H23N3O4"),
    "Gefitinib": ("COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4", "C22H24ClFN4O3"),
    "Sorafenib": ("CNC(=O)C1=CC=CC=C1OC2=C(N=CC=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F", "C21H16ClF3N4O3"),
    "Axitinib": ("CNC(=O)c1ccccc1Sc1ccc2c(/C=C/c3ccccn3)n[nH]c2c1", "C22H18N4OS"),
    "Sunitinib": ("CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\\C(=O)Nc3ccc(F)cc23)c1C", "C22H27FN4O2"),
    "Lapatinib": ("CS(=O)(=O)CCNCc1ccc(o1)-c1ccc2ncnc(Nc3ccc(OCc4cccc(F)c4)c(Cl)c3)c2c1", "C29H26ClFN4O4S"),
    "Dasatinib": ("Cc1nc(Nc2ncc(s2)C(=O)Nc2c(C)cccc2Cl)cc(n1)N1CCN(CCO)CC1", "C22H26ClN7O2S"),
    "Pazopanib": ("Cc1ccc(Nc2nccc(n2)N(C)c2ccc3c(C)n(C)nc3c2)cc1S(N)(=O)=O", "C21H23N7O2S"),
    "Bortezomib": ("CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnccn1)B(O)O", "C19H25BN4O4"),
    "Curcumin": ("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O", "C21H20O6"),
    # Quinone-type Cdc25/phosphatase inhibitor used in the literature as a tool compound.
    # NOTE: quinones are redox-cycling / assay-interference prone (PAINS-type); treat with care.
    "NSC 95397 (quinone; PAINS-prone)": ("OCCSC1=C(SCCO)C(=O)c2ccccc2C1=O", "C14H14O4S2"),
}

GBM_REFERENCE_SMILES: dict[str, str] = {k: v[0] for k, v in _VERIFIED.items()}
EXPECTED_FORMULAS: dict[str, str] = {k: v[1] for k, v in _VERIFIED.items()}
