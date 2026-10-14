"""Short, scannable in-app guides: why / what you need / steps / expect / don't / if it fails.
Text describes only behaviour that exists in the app. Research and education only."""
import streamlit as st

G = {
 "chat_pdf": dict(title="Chat with your PDFs", why="Find the exact sentence in your papers that answers a question, with file name and page number.",
   need="Text-based PDFs (not scans).", steps=["Upload one or more PDFs.", "Ask a specific question using the paper's own words (e.g. 'which cell lines were used?').", "Click the cited page in your PDF to verify."],
   expect="1-5 short verbatim quotes with file + page. It never writes new sentences, so it cannot invent facts.", dont="Don't treat the 0-100 'coverage' number as correctness, and don't expect a summary.",
   fails="No results? The PDF may be a scan (run OCR first) or your question uses different words than the paper: rephrase."),
 "screen": dict(title="Compound screening", why="Quick drug-likeness and brain-penetration (BBB) rule check before spending time on docking.",
   need="A molecule as SMILES (a text code for a structure; copy it from PubChem).", steps=["Paste the SMILES.", "Run the screen.", "Read which rule sets are met."],
   expect="Descriptors and rule-based flags (screening aid only).", dont="Don't read 'meets thresholds' as 'the drug works'.", fails="Error on SMILES? Check it on pubchem.ncbi.nlm.nih.gov or redraw it."),
 "dock": dict(title="Molecular docking", why="Predict how a molecule fits into a protein pocket and which residues it touches.",
   need="A protein structure (PDB ID or file) and a ligand (SMILES or file).", steps=["Run the Docking Validation tab first on your protein family.", "Load the protein and define the pocket (box around a known ligand if present).", "Add your ligand and run.", "Inspect the 3D view and the contact tables."],
   expect="Binding energy in kcal/mol (more negative = stronger predicted binding) and ranked poses.", dont="Don't compare energies across different proteins, and never call a docking score proof of binding.",
   fails="'Engine missing' means the server lacks Vina/OpenBabel. Ligand fails? Check the SMILES and charges."),
 "validate": dict(title="Docking validation (redocking)", why="Proves your docking set-up can reproduce a real crystal structure before you trust it on new molecules.",
   need="A PDB file that contains a protein AND its bound drug (or a 4-letter PDB ID).", steps=["Upload the PDB or type its ID.", "Pick the ligand to redock.", "Run with 2-3 seeds (1-3 minutes each)."],
   expect="RMSD <= 2.0 Å = pass (imatinib in Abl gave 0.90 Å in our test). Orange = real pose, green = docked.", dont="Don't skip it, and don't claim one pass validates other targets.",
   fails="FAIL? Check the pocket box, ligand protonation and try higher exhaustiveness. No ligand listed? That structure has no drug-like ligand."),
 "lit": dict(title="Literature, references and citations", why="Find papers on PubMed and cite them correctly.",
   need="Internet access.", steps=["Search a topic or paste a PMID/DOI.", "Open the paper and read the abstract.", "Copy the formatted citation."],
   expect="Titles, authors, PMIDs from PubMed/Crossref.", dont="Don't cite a paper you haven't opened.", fails="Empty results? Check spelling, or the service may be rate-limited: retry in a minute."),
 "live": dict(title="Live data explorer", why="Look up public data: mutations (cBioPortal), cell lines, drugs (ChEMBL), trials (ClinicalTrials.gov).",
   need="Internet access.", steps=["Choose the sub-tab.", "Enter your query (e.g. a gene or drug name).", "Check the source and date shown with each result."],
   expect="Records labelled with their source. Nothing is invented if a service is down.", dont="Don't use it for patient decisions.", fails="Error shown? The service is unreachable or changed: try later."),
 "ic50": dict(title="IC50 curve fit (4PL)", why="Estimate how much drug halves a response from dose-response data.",
   need="Concentrations and responses; at least 6 concentrations spanning the curve, ideally with replicates.", steps=["Enter concentration,response pairs.", "Run the fit.", "Check the curve shape before reading the IC50."],
   expect="An IC50 estimate and the fitted curve.", dont="Don't trust an IC50 outside your tested range or when the curve has no plateau.", fails="Fit fails? Use more points, check units and that responses are not all equal."),
 "stats": dict(title="Statistics: survival and docking enrichment", why="Test whether two patient groups survive differently, or whether your docking separates real binders from decoys.",
   need="CSV: time, event (1/0), group  |  or score, label (1 active / 0 decoy).", steps=["Upload your CSV.", "Press Run.", "Read p-value and curves (or AUC / enrichment)."],
   expect="Log-rank chi-square and p, Kaplan-Meier curves; AUC 0.5 = random, 1.0 = perfect.", dont="Don't claim significance with very few events; report n.", fails="Error? Check column names and that there are exactly 2 groups."),
 "omics": dict(title="Omics and tools", why="Find genes that differ between groups, test gene sets, check compounds and ligand contacts.",
   need="Expression table (genes x samples) + group labels, or SMILES, or a PDB + ligand.", steps=["Click 'Try simulated demo data' first to see the output.", "Download the template and fill in your data.", "Read 'What I understood from your file'.", "Use FDR < 0.05, not the raw p-value."],
   expect="Ranked genes with log2 fold change and FDR, volcano plot, PCA.", dont="Don't use raw RNA-seq counts (use DESeq2/edgeR), and don't trust results with 3 samples per group.", fails="The app tells you what it could not read. Fix the layout or send the first 5 rows of your file."),
 "assistant": dict(title="AI chat assistant", why="Ask GBM/oncology questions and get answers tied to references.",
   need="Nothing for basic use; an API key only if you enable the optional polish step.", steps=["Ask a specific question.", "Check the cited references.", "Open the PMID/DOI before using a claim."],
   expect="Answers from the app's curated notes and literature search; it says so when nothing validated is found.", dont="Don't treat it as medical advice or as a replacement for reading the paper.", fails="'No validated evidence found'? Rephrase, or use Chat with Papers on your own PDFs."),
 "notebook": dict(title="Lab notebook", why="Keep dated experiment notes and export them.",
   need="Nothing.", steps=["Write an entry.", "Download a backup before closing the page."], expect="A list of your entries you can export.",
   dont="Don't store patient data or confidential files.", fails="Entries gone? The session ended: restore from your backup file."),
}

GOALS = {
 "I have PDFs and questions": ("Chat with Papers", "Ask in the paper's own words and check the cited page."),
 "I have a molecule (SMILES)": ("Compound Screening", "Quick rule check before docking."),
 "I want to dock a molecule into a protein": ("Docking Validation, then Molecular Docking", "Validate on a crystal structure first, then dock."),
 "I have gene-expression data (CSV)": ("Omics & Tools", "Click 'Try simulated demo data' first."),
 "I have survival data": ("Statistics", "Upload time/event/group and read the log-rank result."),
 "I have dose-response data": ("4PL IC50 Fit", "Enter concentration,response pairs."),
 "I need papers or citations": ("Literature / Citations", "Search PubMed, copy formatted citations."),
 "I am lost": ("Start with the demo", "Open Omics & Tools, press 'Try simulated demo data', then Run."),
}


def render_guide(key: str):
    g = G[key]
    with st.expander(f"🧭 30-second guide: {g['title']}", expanded=False):
        st.markdown(f"**Why use it?** {g['why']}  \n**You need:** {g['need']}")
        st.markdown("**Steps**  \n" + "  \n".join(f"{i}. {s}" for i, s in enumerate(g["steps"], 1)))
        c1, c2, c3 = st.columns(3)
        c1.success(f"**Expect**\n\n{g['expect']}")
        c2.warning(f"**Don't**\n\n{g['dont']}")
        c3.info(f"**If it fails**\n\n{g['fails']}")


def start_here():
    with st.expander("👋 New here? Start here (30 seconds)", expanded=False):
        goal = st.selectbox("What do you want to do?", list(GOALS), key="start_goal")
        tab, why = GOALS[goal]
        st.success(f"**Go to the tab: {tab}.** {why}")
        st.caption("Every tab has a '🧭 30-second guide' at the top. Results are research aids, never medical advice.")
