# Changelog

## 2.0.0 - 2026-10-10
### Removed (credibility)
- Hard-coded "clinical trial matching", sample patients, "treatment planning" and molecular-subtype recommendations.
- Hard-coded cell-line mutation tables and placeholder IC50 values.
- 27 drug SMILES whose molecular formula did not match the named drug (found by automated check).
### Added
- `data_sources.py`: live, provenance-tagged adapters (cBioPortal, ChEMBL, ClinicalTrials.gov registry, Europe PMC).
- Live Data Explorer tab and Limitations & Honesty tab.
- `runcard.py` + JSON run-card download for docking and lookups.
- `reference_compounds.py` (21 formula-verified structures) and pytest suite (28 tests), CI workflow.
- CITATION.cff, CONTRIBUTING, SECURITY, code of conduct, issue template, AUDIT.md.
### Fixed
- `requirements.txt` listed `autodock-vina`/`openbabel` as pip packages (not installable); now system packages.
- Streamlit config moved to `.streamlit/` (it was ignored at repo root); XSRF protection enabled.
- Dockerfile healthcheck needed `curl`; Vina/OpenBabel now installed.
- Reworded app title/subtitle/prompts to remove clinical-decision-support language.
## 1.0.0
- Initial release.
