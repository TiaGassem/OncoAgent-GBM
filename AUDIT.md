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
- Never run end-to-end with a real Vina + OpenBabel docking job in this audit; test it yourself on a co-crystal case (redock, RMSD <= 2 A).
- Live adapters (cBioPortal, ChEMBL, ClinicalTrials.gov, Europe PMC) were not exercised against the real services here.
- `protocols.py`, `library_module.py`, `assistant_llm.py`, the reference guide: typed-in text. Each claim still needs a PMID/DOI or an "author convention" label; only the MGMT figure was verified.
- Structural-alert "toxicity" screen and the composite weights are author conventions, not validated models.
- Not implemented: TCGA/CGGA survival statistics, DepMap, GDSC, Cellosaurus API, Arabic RTL, docking decoy enrichment.
- No expert review (neuro-oncologist / cell biologist) and no legal review of data-source licences.
- Redocking RMSD has no symmetry correction.
