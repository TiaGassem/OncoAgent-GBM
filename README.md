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

## Features

- **Compound Screening** -- Physicochemical analysis, BBB permeability scoring, Lipinski/Veber rules
- **Molecular Docking** -- PDB-based docking with AutoDock Vina, phosphatase target database
- **Cell Line Database** -- 18 brain cancer cell lines with mutation profiles and drug sensitivity data
- **Clinical Trial Matching** -- 20 GBM trials with mutation-based eligibility matching and scoring
- **Treatment Planning** -- Molecular subtype classification and evidence-based treatment recommendations
- **Literature Research** -- PubMed search with APA/BibTeX citation generation
- **PDF Audit Reports** -- Downloadable evaluation reports

## Supported Cell Lines

U251 MG, U87 MG, U373 MG, T98G, A172, LN229, LN18, SF295, SNB19, U118 MG, U138 MG, H4, D54, CASI-1, GL261, RCAS-PDGFBA

## Molecular Targets

Cdc25A/B/C, MKP-1, SHP-1, SHP-2, PTEN, PTP1B, TC-PTP, and other brain-relevant phosphatases

## Built With

- Streamlit
- RDKit
- AutoDock Vina
- Presidio (PII anonymization)
- fpdf2 (PDF generation)
