---
title: OncoAgent-GBM
emoji: brain
colorFrom: blue
colorTo: indigo
sdk: streamlit
sdk_version: 1.32.0
app_file: app.py
pinned: true
license: mit
---

# OncoAgent-GBM

Glioblastoma Multiforme Drug Discovery and Clinical Decision Support Platform.

## Overview

OncoAgent-GBM is a local, privacy-first Python platform designed for Glioblastoma (GBM) drug discovery research. It integrates compound screening, molecular docking, literature research, clinical trial matching, and treatment planning into a single professional dashboard.

**Live Demo:** [https://oncoagent-gbm.streamlit.app](https://oncoagent-gbm.streamlit.app)

## Features

### 1. Compound Screening & BBB Triage
- Physicochemical analysis (MW, LogP, TPSA, HBD, HBA)
- Blood-Brain Barrier (BBB) permeability scoring
- Lipinski Rule of 5 and Veber rules compliance
- 52-compound GBM drug library including NSC-95397 and related naphthoquinones
- **Rule-based toxicity pre-screen** using structural alerts (PubChem patterns) + LogP/MW thresholds (Ames, hERG, hepatotoxicity, LD50, carcinogenicity). This is a transparent heuristic, **NOT ProTox-3** and not a machine-learning predictor.

### 2. Molecular Docking
- PDB receptor fetching from RCSB
- AutoDock Vina docking engine (when available), plus SwissDock and CB-Dock2 as web alternatives
- Validated reproducible grid (vGrid): center=(48.164, 10.08, 3.111), size=(29.5, 37.1, 26.9)
- 3D visualization with py3Dmol
- Phosphatase target database (Cdc25A/B/C, MKP-1, SHP-1/2, PTEN, PTP1B)

### 3. Cell Line Database
- 18 brain cancer cell lines (U251, U87, U373, T98G, A172, LN229, LN18, SF295, SNB19, U118, U138, H4, D54, CASI-1, GL261, RCAS-PDGFBA)
- Mutation profiles (PTEN, TP53, IDH1, EGFR, NF1, BRAF, CDKN2A, MGMT)
- Drug sensitivity data (IC50 values for 12+ compounds)
- Cross-line comparison tool

### 4. Clinical Trial Matching
- 20+ GBM clinical trials with NCT IDs
- Mutation-based eligibility matching algorithm
- Profile validation and scoring
- Evidence-based treatment recommendations

### 5. Literature & Bibliography
- PubMed API search
- APA and BibTeX citation generation
- GBM reference guide

### 6. AI Chat Assistant (source-cited)
- Context-aware GBM research Q&A with topic chips: Docking | Cell Line | Toxicity | Targets | Trials | Sources
- Literature answers backed by a real PubMed-built paper index
- No diagnosis and no clinical predictions; every reply ends with a disclaimer and validated source links

### 7. Export
- PDF audit reports (FPDF2)
- Word/DOCX reports (python-docx)
- CSV data export
- Multilingual support (English, French, Arabic)

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Frontend | Streamlit |
| Cheminformatics | RDKit |
| Docking | AutoDock Vina |
| Toxicity | Rule-based structural alerts (PubChem patterns) - NOT ProTox-3 |
| PDF | FPDF2 |
| Word | python-docx |
| PII | Presidio |
| API | PubMed E-utilities |

## Installation

```bash
git clone https://github.com/TiaGassem/OncoAgent-GBM.git
cd OncoAgent-GBM
pip install -r requirements.txt
python -m spacy download en_core_web_sm
streamlit run app.py
```

## Deployment

### Streamlit Community Cloud
1. Push to GitHub
2. Go to https://share.streamlit.io
3. Select repo and deploy

### Hugging Face Spaces
1. Create a new Space at https://huggingface.co/new-space
2. Select Streamlit SDK
3. Upload project files

## Project Structure

```
OncoAgent-GBM/
  app.py                  # Main Streamlit dashboard
  cheminformatics.py      # RDKit compound analysis
  docking_engine.py       # Molecular docking engine
  agent_evaluator.py      # Lead evaluation + PDF/DOCX reports
  research_module.py      # PubMed literature search
  anonymizer.py           # PII anonymization
  requirements.txt        # Python dependencies
  Dockerfile              # Container deployment
  .streamlit/config.toml  # Streamlit theme config
  LICENSE                 # MIT License
```

## Validated Sources (credibility layer)

Every result is traceable to a real public source:

- Clinical trials: https://clinicaltrials.gov/
- Cell lines: https://cellosaurus.org/ and https://dtp.cancer.gov/
- Toxicity / structural alerts: https://pubchem.ncbi.nlm.nih.gov/
- Docking: https://swissdock.ch/ , https://cbl-dock2.mohit.bio/ , https://autodock.scripps.edu/
- Receptor structures: https://www.rcsb.org/
- Literature: https://pubmed.ncbi.nlm.nih.gov/

## Disclaimer

This platform is for **research purposes only**. **Not medical advice.** All computational predictions must be validated experimentally. This does not constitute medical advice or clinical decision-making guidance.

## License

MIT License - see [LICENSE](LICENSE) for details.
