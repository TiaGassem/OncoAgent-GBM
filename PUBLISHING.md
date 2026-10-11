# Publishing checklist

Do these in order. Do not skip step 1.

1. **Test on your own machine.** `pip install -r requirements-dev.txt`, `pytest -q`, then run `python scripts/redock_1iep.py`
   (expect RMSD <= 2 A, about 0.9 A was obtained) and `python scripts/smoke_live.py`. Click every tab in the live app.
2. **Get a scientific review** of the typed-in protocol/reference text from Dr. Jemaa or another biologist.
   Add a PMID/DOI or "author convention" label to each claim.
3. **Push to GitHub** (including hidden `.github/` and `.streamlit/`). Confirm the CI badge is green.
4. **Make a release** (GitHub > Releases > tag `v2.2.0`).
5. **Zenodo DOI:** log in to zenodo.org with GitHub, enable the repo, then publish the release. Copy the DOI into
   `CITATION.cff` and the README.
6. **Describe it honestly** in CVs and applications: "open-source research prototype for glioblastoma drug-discovery
   workflows, with automated tests and documented limitations". Do not call it clinical, validated, or world-class.
7. **Journal option later:** the Journal of Open Source Software expects a mature, documented, tested tool with a
   public development history and research impact. Check its current criteria at joss.theoj.org before submitting.
   A short methods/validation preprint (bioRxiv/Zenodo) is the realistic next step.
