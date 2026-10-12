# Audit of OncoAgent-GBM_final.zip (performed 2026-10-10)

## Method
Static read of every file, `py_compile`, import-chain check, RDKit formula check of every hard-coded SMILES,
headless Streamlit load (`streamlit.testing.AppTest`), pytest. Live network calls could NOT be run in the audit sandbox.

## Findings fixed in v2.0.0
| # | Finding | Severity | Action |
|---|---|---|---|
| 1 | 27 of 31 checkable drug SMILES had a formula different from the named drug (e.g. vincristine, topotecan, bortezomib, thalidomide) | Critical | Replaced by 21 formula-verified structures + test |
| 2 | Hard-coded trial list with outcome claims and a patient "matching" engine; sample patients; "Treatment Planning" | Critical (clinical-claim risk) | Removed; registry listing only |
| 3 | Hard-coded cell-line mutations and placeholder IC50 values (with duplicate dict keys) | High | Removed; link-outs to Cellosaurus/DepMap |
| 4 | `autodock-vina`, `openbabel` in requirements.txt (not pip-installable) | High (deploy) | Moved to system packages |
| 5 | Title/prompt wording claimed "clinical decision support" | High | Reworded |
| 6 | `config.toml` at repo root was ignored; XSRF protection disabled | Medium | Moved to `.streamlit/`, enabled |
| 7 | Dockerfile healthcheck used curl not installed | Low | Fixed |
| 8 | README cited `oncoagent_methods.bib` that did not exist | Medium | Added (verify DOIs) |

## Status after v2.1.0 (2026-10-11)
Verified by running: 11 pytest tests pass, ruff (syntax/undefined-name) clean, app loads headless without exceptions,
PDF and DOCX report export produce valid files, redocking RMSD handles PDBQT polar hydrogens.

## NOT done / still open (state these honestly in the thesis and in applications)
- Docking was validated only on imatinib/Abl (1IEP, a standard test case, ligand taken from the crystal conformer). It has NOT been validated on CDC25; do that with a CDC25 structure that has a bound ligand, if one exists in the PDB.
- Live adapters (cBioPortal, ChEMBL, ClinicalTrials.gov, Europe PMC): the audit sandbox has no access to these hosts. Confirmed only that they fail safely with no fake data. Run `python scripts/smoke_live.py` on your computer.
- `protocols.py`, `library_module.py`, `assistant_llm.py`, the reference guide: typed-in text. Each claim still needs a PMID/DOI or an "author convention" label; only the MGMT figure was verified.
- Structural-alert "toxicity" screen and the composite weights are author conventions, not validated models.
- Survival and enrichment statistics exist for user-uploaded CSVs only; automatic TCGA/CGGA download, DepMap, GDSC and Cellosaurus are NOT implemented.
- Arabic/French strings in the new Statistics tab were written by an AI and need review by a native speaker (you).
- No expert review (neuro-oncologist / cell biologist) and no legal review of data-source licences.
- Redocking RMSD has no symmetry correction.

## Scope notes for v2.3.0
- Omics tools are exploratory and assume log-scale expression; they do NOT replace DESeq2/edgeR/limma on raw counts.
- H-bond 'candidates' use distance only (no angle/donor-acceptor typing). PAINS/Brenk flags are warnings, not proof.
- Not hosted (needs GPU/servers): structure prediction, protein design, generative models. Resources tab links to the originals.
- FR/AR text is limited to tab labels and the Statistics tab; it needs native-speaker review.

## v2.4.0 notes
- 3D viewer: JavaScript syntax-checked and exercised against a stubbed 3Dmol (correct labels/lines/box calls for real 1IEP data); NOT yet seen rendering in a browser by the author of this audit.
- Redocking: validated on 1IEP only, with ligand coordinates from a public example. Needs RCSB download (internet) and OpenBabel + Vina on the server.
- Omics: statistical core cross-checked against SciPy; loaders tested on synthetic messy layouts, not yet on real TCGA/CGGA/GEO exports. Send failing files' first rows + the error text so the loader can be fixed.
