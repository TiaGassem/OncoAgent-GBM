# Changelog

## 2.4.1 - 2026-10-13
### Fixed (app stuck on "Your app is in the oven" on Streamlit Cloud)
- `.streamlit/config.toml` (active since 2.2.0, previously ignored at repo root) set `server.port`, `server.address` and `headless`; Streamlit Cloud manages these itself and overriding them can leave the app stuck on startup. Removed (the Dockerfile still passes them on its command line).
- Removed the extra `vina` pip dependency (could fail to install on some Python versions); Vina comes from `packages.txt` (`autodock-vina`), and the app reports clearly if it is missing.
### Verified
- Server starts and answers its health check in a clean venv built from requirements.txt.
- Redocking passes through the command-line Vina path as well (0.90 A, same as the Python path), tested with Vina 1.2.5. Debian's apt package may be an older Vina; not tested.

## 2.4.0 - 2026-10-12
Rebuilt after user testing (3D view unclear, omics files rejected, redocking confusing).
### Changed
- **3D viewer rewritten** (`viewer3d.py`, 3Dmol.js pinned): protein name/PDB ID, pose selector with Vina energies, pocket residues labelled by name+number, dashed distance lines with Å labels (N/O pairs <= 3.5 A), grid box, optional crystal ligand overlay. Spin is off by default. Tables of contacts and distances under the viewer.
- **Redocking is now automatic** (`redock_auto.py`): upload a PDB with a bound ligand (or an RCSB ID), pick the ligand, press run. Extracts and removes the ligand, builds the box, docks with 1-3 seeds, reports top-pose and best-pose RMSD, shows crystal-vs-docked overlay. Manual PDBQT pasting moved to an "Advanced" expander.
- **Omics inputs are forgiving** (`omics_io.py`): CSV/TSV/TXT/XLSX, auto-detects genes-in-rows vs samples-in-rows, averages duplicate genes, drops text columns, fuzzy sample-name matching, name-rule groups, and tells the user exactly what it understood and what is wrong. Templates + simulated demo data (clearly labelled) + plain-language help.
- Streamlit pinned to the tested 1.65.x; HTML views use `st.iframe` with a fallback (`components.html` is deprecated).
### Verified
- Real redocking of imatinib/Abl from a plain PDB complex through the wizard: top-pose RMSD 0.90 A with two seeds, clicked through in the running app.
- Demo omics run recovers all 4 planted genes across 8 random seeds; fresh venv from requirements.txt passes 34 tests.
### Not verified
- Viewer rendering in a real browser (JS validated with node against a stub only). Open `demo_redocking_viewer_1IEP.html` to check.

## 2.3.0 - 2026-10-11
### Added
- **Omics & Tools tab** (labels in EN/FR/AR): differential expression (Welch / Mann-Whitney + Benjamini-Hochberg FDR, volcano, PCA), gene-set over-representation from a GMT file, single-gene co-expression, on user-uploaded tables.
- **Compound tools:** PAINS and Brenk structural alerts, QED, Murcko scaffolds, Tanimoto similarity, Butina clustering (RDKit).
- **Structure tools:** ligand-residue contacts and H-bond candidates; AlphaFold pLDDT summary.
- 10 new tests cross-checked against SciPy/RDKit. On imatinib-Abl (1IEP) the contact tool recovers the literature residues (Thr315, Met318, Glu286, Asp381, His361).
- Resources tab: link-outs for what needs GPUs or specialised servers (AlphaFold DB, ESM Atlas, SwissADME, ProTox-3, Open Targets, DepMap).

## 2.2.1 - 2026-10-11
### Fixed
- **Docking crashed on every run** with the Python Vina binding (`dock()` got an unexpected `energy_range`; `cpu` belongs to the constructor). Found by a real redocking test.
- Deprecated `utcnow()`; `meeko` needs `gemmi` (added to requirements, plus `vina`).
### Verified (this release)
- Real redocking of imatinib into Abl (PDB 1IEP) through the app's own pipeline: top-pose RMSD 0.90 A (seed 42) and 0.90 A (seed 7), best energy -13.26 / -13.25 kcal/mol. Re-run with `python scripts/redock_1iep.py`.
### Added
- `scripts/redock_1iep.py`, `scripts/smoke_live.py` (live-service check; fails loudly, never shows fake values).

## 2.2.0 - 2026-10-11
### Added
- `stats_module.py` + Statistics tab: Kaplan-Meier, two-group log-rank, docking ROC AUC / enrichment factor on user-supplied CSVs; EN/FR/AR text. Log-rank reproduces the Freireich textbook values in tests.
- `PUBLISHING.md` (release + Zenodo checklist), ORCID in `CITATION.cff`, README brought up to date.

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
