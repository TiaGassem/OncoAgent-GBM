"""Short, scannable in-app guides: why / what you need / steps / expect / don't / if it fails.
Text describes only behaviour that exists in the app. Research and education only."""
import re
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


# ---------------- trust map ----------------
CRED = [
 ("Docking (AutoDock Vina)", "Real engine; works on ANY protein and ligand you provide. Evidence shown here: redocking test case imatinib/Abl (0.90 Å). Validate each new protein yourself (Docking Validation) when a crystal ligand exists.", "Cite Vina (Trott & Olson 2010; Eberhardt 2021). Report receptor/ligand preparation, box, exhaustiveness, seeds."),
 ("Redocking test", "Standard method (RMSD ≤ 2 Å).", "A pass supports the set-up for that protein; it does not prove scores for new ligands."),
 ("Survival / enrichment statistics", "Log-rank matches the Freireich textbook result (16.79); AUC/EF unit-tested.", "Report n and number of events."),
 ("Omics (differential expression)", "Welch test + BH-FDR cross-checked against SciPy. Exploratory.", "Not DESeq2/edgeR: use those for raw counts."),
 ("Compound alerts (PAINS/Brenk)", "RDKit's published catalogues.", "Warnings only, never proof of assay interference."),
 ("ADME / BOILED-Egg / toxicity screen", "Local approximation (RDKit descriptors, structural alerts). NOT identical to SwissADME or ProTox-3 and not validated against them.", "For papers, use SwissADME / ProTox-3 and cite them."),
 ("IC50 fit (4PL)", "Recovers known parameters on simulated data (tested).", "Needs points on both plateaus; report R² and n."),
 ("PDF chat", "Returns verbatim quotes with page numbers; cannot invent text.", "Quote the paper, not the app."),
 ("AI assistant / curated notes", "Partly unverified text; says so when nothing validated is found.", "Open the cited PMID/DOI before using a claim."),
 ("Live databases", "Adapters written, not tested at scale.", "Record source and date."),
 ("Target fishing / bioactivity prediction", "NOT built in.", "Use SwissTargetPrediction or PLATO, then verify in ChEMBL; cite them."),
]

# ---------------- sample files (so nobody is blocked by 'I have no data') ----------------
def sample_files():
    import numpy as np
    import omics_io as oio
    X, g = oio.demo_dataset()
    gr = "sample,group\n" + "\n".join(f"{s},{v}" for s, v in g.items()) + "\n"
    mp = [6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35]; me = [1,1,1,0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0]
    ct = [1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23]
    surv = "time,event,group\n" + "".join(f"{t},{e},6-MP\n" for t, e in zip(mp, me)) + "".join(f"{t},1,control\n" for t in ct)
    rng = np.random.default_rng(3)
    enr = "score,label\n" + "".join(f"{v:.2f},1\n" for v in rng.normal(-9, 1, 30)) + "".join(f"{v:.2f},0\n" for v in rng.normal(-7, 1, 300))
    return {"expression": ("SIMULATED expression table", X.to_csv(), "demo_expression.csv"), "groups": ("Group labels for it", gr, "demo_groups.csv"),
            "survival": ("Real textbook survival data (Freireich 1963, leukaemia, not GBM)", surv, "freireich_survival.csv"),
            "enrichment": ("SIMULATED docking scores (30 actives, 300 decoys)", enr, "demo_enrichment.csv")}


def _dl(keys):
    sf = sample_files()
    for k in keys:
        lab, txt, fn = sf[k]
        st.download_button(f"⬇ {lab}", txt, fn, key=f"dl_{k}_{fn}")


R = {
 "pdf": ("Chat with Papers", ["Upload your PDFs.", "Ask with the paper's own words.", "Click the cited page to verify."], None, None),
 "ic50": ("4PL IC50 Fit", ["Type concentration,response pairs (one per line).", "Use at least 6 concentrations spanning the whole curve.", "Check the curve before reading the IC50."], None, "Export from your plate reader to Excel with two columns: concentration, response."),
 "r_screen": ("Compound Screening", ["Paste the SMILES (copy it from PubChem).", "Run and read the rule checks."], None, "Local screen only. For ADME/toxicity in a paper use SwissADME and ProTox-3."),
 "r_dock": ("Docking Validation, then Molecular Docking", ["Docking Validation: choose 'RCSB PDB ID', type 1IEP, download, pick the ligand, run (shows the test passing).", "Do the same with a crystal structure of YOUR protein if one has a bound ligand.", "Then Molecular Docking: same receptor preparation, same box, 3 seeds.", "Report mean ± SD and every setting."], None,
            "Why your thesis Vina results may differ: scores depend on receptor/ligand preparation, box, exhaustiveness, Vina version and seed. We measured -13.26 vs -13.93 kcal/mol for the SAME complex with two ligand preparations. Never compare scores prepared differently."),
 "r_target": ("Live Data Explorer (ChEMBL)", ["Look up the compound's known bioactivity and mechanism.", "For target PREDICTION use SwissTargetPrediction or PLATO, then confirm in ChEMBL."], None, "Target fishing is not built into this app."),
 "r_omics": ("Omics & Tools", ["Press 'Try simulated demo data' and Run to see the output.", "Then use your own table: genes in rows, samples in columns.", "Read 'What I understood from your file'.", "Use FDR < 0.05."], ["expression", "groups"], None),
 "r_surv": ("Statistics", ["Columns: time, event (1 = event, 0 = censored), group (two groups).", "Upload and run."], ["survival"], None),
 "r_enr": ("Statistics", ["Columns: score, label (1 = known active, 0 = decoy).", "Tick 'lower is better' for Vina energies."], ["enrichment"], None),
 "r_demo": ("Omics & Tools", ["Open the 'Omics & Tools' tab.", "Press 'Try simulated demo data', then 'Run differential expression'."], None, None),
 "r_lit": ("Literature & Bibliography", ["Search your topic.", "Open each paper before citing it."], None, None),
 "r_data": ("Free public data", ["cBioPortal (cbioportal.org): pick a glioblastoma study, download clinical (survival) and expression data.", "GEO (ncbi.nlm.nih.gov/geo): search 'glioblastoma', download a series matrix.", "GlioVis and CGGA: glioma expression + survival.", "From your lab: Excel with gene names in column 1 and one column per sample.", "Always cite the dataset and check its licence."], ["expression", "survival"], None),
}
TREE = {
 "root": ("What do you have right now?", [("📄 Papers (PDF)", "pdf"), ("🧪 A molecule or a protein", "mol"), ("📊 Numbers in a table", "table"), ("🔬 Lab dose-response results", "ic50"), ("🤷 Nothing yet / I am lost", "none")]),
 "mol": ("What do you want to know?", [("Is it drug-like / can it reach the brain?", "r_screen"), ("How does it bind a protein?", "r_dock"), ("What are its known targets / bioactivity?", "r_target")]),
 "table": ("What do the numbers measure?", [("Gene expression (tumour vs normal)", "r_omics"), ("Patient survival", "r_surv"), ("Docking scores of actives vs decoys", "r_enr")]),
 "none": ("What would you like?", [("See a working demo in 10 seconds", "r_demo"), ("Find papers on my topic", "r_lit"), ("Learn how to get free public data", "r_data")]),
}
KEYWORDS = {"pdf": ["pdf", "paper", "article"], "r_screen": ["smiles", "drug-like", "lipinski", "bbb"], "r_dock": ["dock", "vina", "binding", "protein"], "r_target": ["target", "bioactiv", "plato", "chembl"],
            "r_omics": ["expression", "omics", "gene", "rna", "csv"], "r_surv": ["survival", "kaplan", "hazard"], "ic50": ["ic50", "dose", "viability"], "r_data": ["data", "dataset", "download"], "r_lit": ["pubmed", "reference", "citation"]}


def _show_result(k):
    tab, steps, files, note = R[k]
    st.success(f"**Go to: {tab}**")
    st.markdown("\n".join(f"{i}. {s}" for i, s in enumerate(steps, 1)))
    if files:
        st.caption("No data yet? Download a sample file to practise:")
        _dl(files)
    if note:
        st.warning(note)


def _guide_chat():
    import guide_bot as gb
    ss = st.session_state
    ss.setdefault("gb_hist", [("assistant", gb.OPENER)]); ss.setdefault("gb_ctx", {}); ss.setdefault("gb_sug", gb.START)

    def send(text):
        ss["gb_hist"].append(("user", text))
        md, ss["gb_ctx"], ss["gb_sug"] = gb.respond(text, ss["gb_ctx"])
        ss["gb_hist"].append(("assistant", md))
    for i, (role, md) in enumerate(ss["gb_hist"]):
        with st.chat_message(role):
            files = None
            m = re.search(r"\[\[DL:([a-z,]+)\]\]", md)
            if m:
                files, md = m.group(1).split(","), re.sub(r"\s*\[\[DL:[a-z,]+\]\]", "", md)
            st.markdown(md)
            if files and i == len(ss["gb_hist"]) - 1:
                st.caption("Sample files to practise with:")
                _dl(files)
    cols = st.columns(len(ss["gb_sug"]) if ss["gb_sug"] else 1)
    for j, sug in enumerate(ss["gb_sug"] or []):
        if cols[j].button(sug, key=f"gb_s{len(ss['gb_hist'])}_{j}"):
            send(sug); st.rerun()
    q = st.chat_input("Type your question or goal here (e.g. 'I am lost', 'what is RMSD?')", key="gb_input")
    if q:
        send(q); st.rerun()


def start_here():
    with st.expander("👋 Lost? Start here: I will guide you", expanded=False):
        _guide_chat()
        with st.expander("How much can I trust each tool? (read before citing)"):
            import pandas as pd
            st.dataframe(pd.DataFrame(CRED, columns=["Tool", "Status", "What to do / cite"]), hide_index=True)
        st.caption("Every tab has a '🧭 30-second guide'. Research and education only, never medical advice.")
