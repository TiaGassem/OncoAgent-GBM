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

## NOT done / still open (be honest in your thesis and applications)
- Live adapters are untested against the real services; run `python -m data_sources` and fix any shape mismatch.
- `agent_evaluator.py` still contains a hard-coded toxicity/lead-evaluation heuristic and `PatientProfile` class; the
  app no longer exposes patient input, but review or remove it. `anonymizer.py` (medical-note PII) is unused by the UI.
- `assistant_llm.py`, `library_module.py` and `protocols.py` contain typed-in text; each protocol/claim still needs a
  PMID/DOI or an "author convention" label. Not reviewed line by line.
- CGGA/TCGA survival statistics, DepMap, GDSC, Cellosaurus API, Arabic RTL layout, docking decoy enrichment: not implemented.
- No expert (neuro-oncologist / cell biologist) review yet. Licences of data sources not legally reviewed.
- Redocking RMSD has no symmetry correction (stated in README).
