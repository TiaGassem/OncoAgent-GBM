# Changelog

## 2.1.0 - 2026-10-11
### Removed (invented numbers / unsupported claims)
- LogP-binned "LD50 estimate" and GHS class (no validated model). LD50 is now "not estimated"; ProTox-3 values can still be entered.
- Patient profile / matching engine (`PatientProfile`, `match_patient`) and its 15% weight in the composite score.
- "VIABLE LEAD / REJECTED" verdicts. Now: MEETS / PARTIAL / BELOW SCREENING THRESHOLDS (screening aid, author-convention weights).
- Chat "Credibility 0-100" badge (was `50 + 12n + 18`). Now shows plain evidence counts. PDF-chat score renamed "query-term coverage".
- Descriptor-based "estimated kcal/mol / Ki" shown when docking is unavailable.
- Approximate PDBQT atom-typing fallback: receptor prep now requires OpenBabel and refuses otherwise.
- Unsourced prevalence numbers in the oncogene guide (EGFRvIII "~50%", IDH1 range). MGMT figure now cites Hegi 2005 (PMID 15758010, verified on PubMed).
- Fake example PMID/DOI/authors pre-filled in the citation tools.
### Fixed
- PDF export of the lead report crashed ("Not enough horizontal space").
- Docking used a receptor that still contained the co-crystal ligand; now protein-only.
- Redocking RMSD miscounted PDBQT polar hydrogens (type HD/HS), so docked poses always failed the atom-count check.
- Missing imports (`re`, `fetch_articles_by_id`) crashed Explore Paper / Citation paths.
- Composite score no longer treats missing docking as a neutral 50.
- `config.toml` moved to `.streamlit/` (it was ignored), XSRF protection on.
### Added
- `tests/` (11 tests), GitHub Actions CI, `.streamlit/config.toml`. These folders were missing from the GitHub repo.

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
