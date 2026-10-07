# OncoAgent-GBM

**An open-source, computational drug-discovery toolkit for glioblastoma (GBM) research.**

OncoAgent-GBM is a Streamlit web application that integrates structure-based
docking, cheminformatics, pharmacokinetic prediction and literature tooling into
a single, reproducible, **no-fabrication** workflow aimed at the CDC25
phosphatase family (CDC25A/B/C) as therapeutic targets in glioblastoma.

> **Research & educational use only — NOT a clinical or diagnostic tool.**
> Every quantitative output is a *computational prediction* that must be
> confirmed experimentally. The app never invents a value: when something
> cannot be computed, it says so explicitly.

---

## Features

| Module | What it does | Engine / method |
|---|---|---|
| **Compound Screening** | Lipinski/Veber/BBB triage of SMILES | RDKit descriptors |
| **BOILED-Egg** | Gastro-intestinal absorption & BBB permeation map | RDKit TPSA/WLOGP (simplified) |
| **Molecular Docking** | Real receptor–ligand docking on YOUR inputs | AutoDock Vina + OpenBabel/Meeko |
| **Docking Validation** | Redocking RMSD vs native crystallographic pose | heavy-atom RMSD |
| **4PL IC50 Fit** | Dose–response fit on YOUR data | SciPy 4-parameter logistic |
| **Toxicity pre-screen** | Rule-based structural-alert heuristic | transparent rules (NOT ProTox) |
| **Literature** | PubMed search + APA/BibTeX citations | NCBI E-utilities |
| **PDF Chat** | Extractive Q&A over uploaded papers | passage retrieval (quotes only) |

## Scientific integrity principles

1. **No fabrication.** No hardcoded binding scores, IC50s or toxicity values are
   ever presented as results. Heuristic estimates are labelled as heuristics.
2. **Reproducibility.** Every docking run records its random seed, AutoDock Vina
   version, exact command and timestamp.
3. **Validation first.** The Docking Validation tab measures redocking RMSD so
   the docking setup can be shown to reproduce known crystallographic poses
   (RMSD ≤ 2.0 Å) before any prediction is trusted.
4. **Traceable citations.** Methods are tied to real references with DOI/PMID
   (see `oncoagent_methods.bib`).

## Known limitations (read before citing)

- Docking uses rigid-receptor Vina with simplified protein preparation
  (protonation, tautomers, waters, metal ions and flexible side chains are not
  rigorously handled). Treat scores as **relative rankings**, not affinities.
- When Meeko/OpenBabel are unavailable, ligands are docked **rigid** (no
  rotatable bonds), which changes results; this is stated in the output.
- Redocking RMSD is a direct heavy-atom RMSD without symmetry correction; for
  symmetric ligands it is an upper bound.
- BOILED-Egg boundaries are a rectangular approximation of the published egg
  ellipses — confirm on SwissADME.
- The toxicity module is a rule-based heuristic, **not** ProTox-3.
- No molecular dynamics is run; the MD section is a methods template only.

## Deployment

Deployed on Streamlit Community Cloud **from this GitHub repository** (not from a
zip). For real docking the repo must contain a `packages.txt` with exactly:

```
autodock-vina
openbabel
```

Python dependencies are in `requirements.txt` (RDKit, SciPy, Meeko, py3Dmol, …).

```bash
# local run
pip install -r requirements.txt
streamlit run app.py
```

## How to validate the docking before you report any number

1. Pick a PDB with a co-crystallised ligand relevant to your target.
2. Extract the native ligand; note its coordinates.
3. Dock that same ligand back with the same grid box.
4. Open **Docking Validation**, paste native + docked pose, read the RMSD.
5. Report RMSD for 3–5 cases. Only then interpret new predictions.

## Citing

If you connect this repository to Zenodo you obtain a citable DOI. Please also
cite AutoDock Vina, RDKit and the method references listed in
`oncoagent_methods.bib`, and verify every DOI/PMID yourself.

## License & disclaimer

Provided for research and education. No warranty. Not for clinical, diagnostic
or therapeutic use.
