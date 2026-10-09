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
| **Pocket detection** | Finds the binding cavity without a defined site | fpocket, or built-in LIGSITE-style geometric detector (NumPy/SciPy) |
| **Docking Validation** | Redocking RMSD vs native crystallographic pose | heavy-atom RMSD |
| **4PL IC50 Fit** | Dose–response fit on YOUR data | SciPy 4-parameter logistic |
| **Toxicity pre-screen** | Rule-based structural-alert heuristic | transparent rules (NOT ProTox) |
| **Literature** | PubMed search + APA/BibTeX citations | NCBI E-utilities |
| **PDF Chat** | Extractive Q&A over uploaded papers | passage retrieval (quotes only) |
| **AI Chat Assistant** | Source-cited GBM Q&A, with optional bring-your-own-key fluent mode | grounded retrieval + optional user LLM (rephrase only) |
| **Lab protocol library** | 16 open-access GBM/oncology assay templates (spheroid, viability/MTT, scratch, clonogenic, Western blot, qRT-PCR, flow cytometry cell-cycle & apoptosis, caspase-3/7, IF/ICC, transwell, comet, EdU/BrdU, docking+MD, in-silico ADMET) | templates with `[verify]` placeholders + live PubMed/PMC lookup |
| **My Lab Notebook** | Write your own notes/protocols, keep them private to your session, export to Word/PDF/Markdown, encrypted backup | python-docx + fpdf2 + Fernet (AES) |

### Privacy of the Lab Notebook
Notebook entries live **only in your browser session** — nothing is written to the
server, no shared database, no logs. Keep your notes by downloading them (Word /
PDF / Markdown) or an AES-encrypted passphrase backup you can re-import. The
passphrase is never seen or stored; lose it and the backup is unrecoverable.

### Use it on phone / tablet / PC
It is a responsive web app: open the link on any device. To get an app-like icon,
"Add to Home Screen" (iOS Safari / Android Chrome) or "Install this site as an app"
(Chrome/Edge on Windows/Mac). No app store, research use only.

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
- When no binding site is defined, the box is placed on a detected cavity
  (fpocket if present, else a built-in LIGSITE-style geometric detector). A
  whole-protein (blind) box is only a flagged, low-confidence last resort.
  For the most credible results, define the site via a co-crystallised ligand
  or catalytic residues.
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
