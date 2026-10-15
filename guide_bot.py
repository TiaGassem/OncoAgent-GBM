"""Conversational guide (rule-based, no LLM, nothing invented): answers questions about the agent and common
analysis concepts, asks follow-up questions, remembers context and builds a step-by-step plan.
respond(text, ctx) -> (markdown, ctx, suggestions). Research and education only."""
import re

OPENER = ("Hi! I'm your guide. Tell me in one sentence what you are studying and what you have "
          "(for example: *'I have RNA-seq data from GBM patients and want to dock temozolomide'*), or ask me anything about this app.")
START = ["I have papers (PDF)", "I have a molecule and a protein", "I have a table of numbers", "I have no data yet"]

# (keywords, tab, one-line instruction)
INTENTS = {
 "data": (["no data", "don't have data", "dont have data", "no csv", "no file", "where to get data", "public data", "download data", "dataset", "nothing yet"],
          "Free public data", "Download a sample file here to practise, or get real public data (cBioPortal, GEO, GlioVis, CGGA)."),
 "dock": (["dock", "docked", "docking", "vina", "binding", "binds", "pocket", "ligand", "inhibitor", "inhibitors", "ligands", "temozolomide"],
          "Docking Validation, then Molecular Docking", "Validate on a crystal structure of your protein family first (type a PDB ID), then dock your drug with the same preparation and 3 seeds."),
 "screen": (["smiles", "drug-like", "druglike", "lipinski", "bbb", "blood-brain", "blood brain", "brain penetr", "admet", "adme"],
            "Compound Screening", "Paste the SMILES and read the rule checks (a local screen: use SwissADME / ProTox-3 for papers)."),
 "target": (["target fishing", "target prediction", "plato", "swisstarget", "bioactivity", "off-target", "off target"],
            "Not built in", "Use PLATO or SwissTargetPrediction, then check known activity in Live Data Explorer (ChEMBL). Cite those tools."),
 "omics": (["expression", "omics", "rna-seq", "rnaseq", "transcript", "differential", "volcano", "tcga", "cgga", "geo", "microarray", "degs"],
           "Omics & Tools", "Upload genes x samples + groups (or press 'Try simulated demo data'), run, and use FDR < 0.05."),
 "surv": (["survival", "kaplan", "log-rank", "logrank", "hazard", "prognosis", "prognostic"],
          "Statistics", "Upload time, event (1/0), group and read the log-rank result and curves."),
 "ic50": (["ic50", "dose-response", "dose response", "viability", "mtt", "crystal violet", "resazurin", "cell line", "cells treated"],
          "4PL IC50 Fit", "Enter concentration,response pairs (6+ concentrations covering the whole curve)."),
 "pdf": (["pdf", "paper", "papers", "article", "articles", "read my"],
         "Chat with Papers", "Upload your PDFs, ask in the paper's own words, click the cited page to verify."),
 "lit": (["pubmed", "citation", "reference", "bibliography", "cite", "literature"],
         "Literature & Bibliography", "Search PubMed, open each paper, copy the formatted citation."),
}
ORDER = ["data", "pdf", "lit", "screen", "dock", "target", "omics", "surv", "ic50"]

# concept answers: (keywords, answer)
GLOSS = [
 (["rmsd"], "**RMSD** = average distance (Å) between the docked pose and the real crystal pose. **≤ 2.0 Å = pass** (standard convention). 0 = identical."),
 (["fdr", "p-value", "p value", "pvalue", "multiple testing"], "**FDR** is the p-value corrected for testing thousands of genes at once (Benjamini-Hochberg). Use **FDR < 0.05**, not the raw p-value."),
 (["log2fc", "fold change", "log2 fold"], "**log2FC**: +1 means 2× higher in group A than B; -1 means 2× lower; 0 = no change."),
 (["log-rank", "logrank", "kaplan"], "**Kaplan-Meier** draws survival over time; the **log-rank test** asks whether two groups' curves differ. A small p means they likely differ. Needs enough events; report n."),
 (["ic50", "4pl", "four parameter"], "**IC50** = concentration giving a 50% response. The **4PL** model fits a curve with top, bottom, slope and IC50. You need points on both plateaus; report R² and n."),
 (["kcal", "binding energy", "vina score", "affinity"], "The **Vina score** (kcal/mol) is a predicted binding score: more negative = stronger predicted binding. It is not a measured affinity. Compare scores only for the same protein prepared the same way."),
 (["smiles"], "**SMILES** is a text code for a molecule (e.g. temozolomide: `Cn1nnc2c(C(N)=O)ncn2c1=O`). Copy it from PubChem."),
 (["pdb id", "what is pdb", "pdb file", "rcsb"], "A **PDB file** is a 3D protein structure; the **PDB ID** is its 4-character code on rcsb.org (e.g. 1IEP)."),
 (["exhaustiveness", "seed"], "**Exhaustiveness** = how hard Vina searches (higher = slower, more thorough). **Seed** = random start; run 3 seeds and report mean ± SD."),
 (["grid box", "search box", "box"], "The **grid box** is the 3D region Vina searches. Centre it on the pocket (e.g. on a known ligand) and make it just big enough for the ligand."),
 (["redock", "validation", "validate"], "**Redocking** removes a crystal ligand and docks it back; if it lands ≤ 2 Å from the real pose, your set-up is credible for that protein. Use the *Docking Validation* tab (type a PDB ID)."),
 (["pains", "alert", "brenk"], "**PAINS/Brenk alerts** flag substructures that often interfere with assays. They are warnings, not proof."),
 (["pca"], "**PCA** shows whether samples cluster by group. If tumour and normal separate, the groups differ globally."),
 (["auc", "enrichment", "decoy"], "**AUC** = chance an active scores better than a decoy (0.5 random, 1.0 perfect). **EF** = how many more actives appear in the top fraction than by chance."),
 (["differs", "differ from", "different from", "different results", "thesis", "not same", "mismatch"], "Docking scores depend on receptor/ligand preparation, box, exhaustiveness, Vina version and seed (we measured -13.26 vs -13.93 for the same complex with two ligand preparations). Use identical files and settings, run 3 seeds, and report differences honestly."),
 (["can i cite", "how do i cite", "how to cite", "cite this", "cite the app", "publish", "journal"], "Cite **AutoDock Vina** (Trott & Olson 2010; Eberhardt 2021), **RDKit**, **SciPy**, and this app via `CITATION.cff`. Report every setting. For ADME/toxicity in papers use SwissADME/ProTox-3. Research use only."),
 (["error", "not working", "doesn't work", "doesnt work", "failed", "crash", "bug"], "Read the red message: it says what failed. Common fixes: wrong file layout (use the template), no internet (PDB/PubMed lookups), docking engine missing on the server. Send me the message and the first 5 rows of your file."),
 (["patient data", "privacy", "confidential", "secure"], "Do **not** upload patient data or confidential files. Use public or de-identified data only."),
 (["what is this", "what can you do", "what does this app", "what is oncoagent"], "OncoAgent-GBM is an open research toolkit for glioblastoma: docking with validation, compound screening, omics and survival statistics, IC50 fitting, cited PDF chat, literature tools and a lab notebook. It runs on **any** protein, drug or dataset you provide. Research and education only."),
 (["which tab", "where do i start", "where to start", "first"], "Quick start: press **'Try simulated demo data'** in *Omics & Tools*, then Run. Then tell me your own goal and I'll make you a plan."),
]
REFUSE = ["should i take", "treat my", "my patient", "prescribe", "what dose should", "cure", "am i going to"]
YES, NO = ["yes", "yeah", "yep", "i do", "i have", "ok"], ["no", "nope", "not yet", "nothing", "i don't", "i dont"]


def _norm(t):
    return " " + re.sub(r"\s+", " ", t.lower().strip()) + " "


def _has(t, kws):
    """Whole-word / whole-phrase matching (so 'differ' never matches 'differential')."""
    return any(re.search(r"(?<![a-z0-9])" + re.escape(k) + r"(?![a-z0-9])", t) for k in kws)


def _plan(found):
    lines, dl = [], []
    if "dock" in found and "screen" not in found:
        pass
    for k in ORDER:
        if k in found:
            _, tab, ins = INTENTS[k]
            lines.append(f"**{tab}**: {ins}")
    if "omics" in found:
        dl += ["expression", "groups"]
    if "surv" in found:
        dl.append("survival")
    if "data" in found and not dl:
        dl += ["expression", "survival"]
    return lines, dl


ONB_Q = {
 1: ("**Step 1 of 3.** What are you trying to achieve right now?", ["Understand a disease or target (read papers)", "Test a drug or molecule", "Analyse my experimental data", "Write my thesis or paper"]),
 2: ("**Step 2 of 3.** What do you already have?", ["Papers (PDF)", "A molecule (name or SMILES)", "Lab or patient-cohort data (table)", "Nothing yet"]),
 3: ("**Step 3 of 3.** How comfortable are you with this kind of analysis?", ["Beginner", "Some experience"]),
}
LOST_WORDS = ["lost", "confused", "stuck", "beginner", "don't know where", "dont know where", "help me", "where do i begin"]
STEPS = {  # key -> (tab, why, guide key)
 "pdf": ("Chat with Papers", "get exact, page-cited answers from your PDFs", "chat_pdf"),
 "lit": ("Literature & Bibliography", "find and cite papers correctly", "lit"),
 "screen": ("Compound Screening", "quick drug-likeness and brain-penetration rule check", "screen"),
 "validate": ("Docking Validation", "prove the docking set-up works (use the built-in example first)", "validate"),
 "dock": ("Molecular Docking", "dock your molecule into your protein and read the contacts", "dock"),
 "ic50": ("4PL IC50 Fit", "turn dose-response numbers into an IC50", "ic50"),
 "omics": ("Omics & Tools", "find genes that differ between groups", "omics"),
 "surv": ("Statistics", "survival curves and log-rank, or docking enrichment", "stats"),
 "notebook": ("Lab Notebook", "keep dated notes and a backup", "notebook"),
 "data": ("Free public data", "sample files + real public datasets (cBioPortal, GEO, GlioVis, CGGA)", None),
}


def _onb_goal(t):
    for k, kws in (("lit", ["understand", "disease", "target", "papers", "read", "literature"]), ("drug", ["test", "drug", "molecule", "compound"]),
                   ("data", ["analyse", "analyze", "data", "experiment"]), ("write", ["write", "thesis", "paper", "publish"])):
        if _has(t, kws):
            return k
    return None


def _onb_assets(t):
    if _has(t, ["nothing", "none", "no data", "nothing yet"]):
        return "none"
    for k, kws in (("pdf", ["papers", "pdf"]), ("mol", ["molecule", "smiles", "compound", "drug"]), ("table", ["lab", "table", "data", "cohort", "csv", "excel"])):
        if _has(t, kws):
            return k
    return None


def _roadmap(goal, assets, level):
    steps = []
    if assets == "none":
        steps.append("data")
    if goal == "lit":
        steps += ["pdf", "lit"] if assets == "pdf" else ["lit", "pdf"]
    elif goal == "drug":
        steps += ["screen", "validate", "dock"] + (["ic50"] if assets == "table" else [])
    elif goal == "data":
        steps += ["omics", "surv", "ic50"] if assets != "pdf" else ["pdf"]
    else:
        steps += ["lit", "pdf", "notebook"]
    seen, out = set(), []
    for k in steps:
        if k not in seen:
            seen.add(k); out.append(k)
    lines = []
    for i, k in enumerate(out, 1):
        tab, why, _ = STEPS[k]
        pre = "start with the example data: " if (level == "beginner" and k in ("omics", "validate", "surv")) else ""
        lines.append(f"{i}. **{tab}**: {pre}{why}.")
    md = "Here is **your personal roadmap**:\n\n" + "\n".join(lines)
    if goal == "drug":
        md += "\n\n⚠ Docking scores depend on preparation, box and seed. Validate first, report every setting, and use SwissADME/ProTox-3 for ADME/toxicity in papers."
    if goal == "write":
        md += "\n\nFor writing: read *Limitations & Honesty*, cite AutoDock Vina, RDKit, SciPy and this app (CITATION.cff), and say what is simulated or unvalidated."
    if assets == "none" or "omics" in out or "surv" in out:
        md += " [[DL:expression,groups,survival]]"
    return md, out


def _explain(n, ctx):
    keys = ctx.get("roadmap") or []
    if not (1 <= n <= len(keys)):
        return None
    gk = STEPS[keys[n - 1]][2]
    if not gk:
        return "**Free public data**: download the sample files below to practise. Real data: cBioPortal (glioblastoma study → clinical + expression), GEO (series matrix), GlioVis, CGGA. Cite the dataset. [[DL:expression,groups,survival]]"
    from guides import G
    g = G[gk]
    return (f"**{g['title']}**\n\n*Why:* {g['why']}\n\n*You need:* {g['need']}\n\n" + "\n".join(f"{i}. {x}" for i, x in enumerate(g["steps"], 1)) +
            f"\n\n*Expect:* {g['expect']}\n\n*Don't:* {g['dont']}\n\n*If it fails:* {g['fails']}")


TOUR = ("**The app in 12 lines**\n\nChat with Papers (cited quotes) · Compound Screening (drug-likeness, BBB) · Molecular Docking (any drug, any protein) · "
        "Docking Validation (redocking test) · Literature (PubMed) · Live Data Explorer (cBioPortal, ChEMBL) · 4PL IC50 Fit · Statistics (survival, enrichment) · "
        "Omics & Tools (differential genes, PCA) · Citations · AI Assistant (grounded notes) · Lab Notebook. Each tab has a 🧭 30-second guide.")


def respond(text, ctx=None):
    ctx = dict(ctx or {})
    t = _norm(text)
    m = re.search(r"(?:explain|show|open)?\s*step\s*(\d)", t)
    if m and ctx.get("roadmap"):
        e = _explain(int(m.group(1)), ctx)
        if e:
            return (e, ctx, ["Explain step 2" if m.group(1) == "1" else "Show me the tour", "Show me the tour", "Start over"])
    if _has(t, ["tour", "all the tabs", "what are the tabs"]):
        return (TOUR, ctx, ["Start over", "What is RMSD?"])
    onb = ctx.get("onb")
    if onb and len(t.split()) > 8 and any(_has(t, v[0]) for v in INTENTS.values()):
        ctx["onb"] = onb = None   # the user described their situation in free text: switch to a plan
    if onb and not _has(t, ["restart", "start over", "reset"]):
        st_ = onb["step"]
        val = _onb_goal(t) if st_ == 1 else _onb_assets(t) if st_ == 2 else ("beginner" if _has(t, ["beginner", "new", "no experience", "never"]) else "some" if _has(t, ["some", "experience", "comfortable"]) else None)
        if val is None:
            q, sug = ONB_Q[st_]
            return ("Sorry, I did not catch that. " + q, ctx, sug)
        onb[["goal", "assets", "level"][st_ - 1]] = val
        if st_ < 3:
            onb["step"] = st_ + 1
            ctx["onb"] = onb
            q, sug = ONB_Q[st_ + 1]
            return (q, ctx, sug)
        md, keys = _roadmap(onb["goal"], onb["assets"], onb["level"])
        ctx["onb"] = None; ctx["roadmap"] = keys
        return (md + "\n\nAsk me to *explain step 1*, take a *tour*, or ask any question.", ctx, ["Explain step 1", "Show me the tour", "Start over"])
    if _has(t, REFUSE):
        return ("I can't give medical advice or treatment decisions: this is a research and education tool. For patient care, please consult a clinician. "
                "I'm happy to help with a research question instead.", ctx, START)
    pend = ctx.get("pending")
    if pend == "have_table":
        ctx["pending"] = None
        if _has(t, NO):
            lines, dl = _plan({"data"})
            return ("No problem. **Start with a sample file** (below) to see how it works, then get real public data: cBioPortal (pick a glioblastoma study, download clinical + expression), "
                    "GEO (search 'glioblastoma', download a series matrix), GlioVis, CGGA. From your own lab: Excel with gene names in column 1 and one column per sample. Always cite the dataset. "
                    "[[DL:expression,groups,survival]]", ctx, ["Show me my plan again", "What is FDR?"])
        if _has(t, YES):
            return ("Great. Genes in rows, samples in columns, first column = gene names; plus a two-column group file (sample, group). "
                    "The app reports what it understood from your file. Open **Omics & Tools** or **Statistics** and upload it. [[DL:expression,groups]]", ctx, ["What is FDR?", "How do I read the result?"])
    if _has(t, ["restart", "start over", "reset"]):
        return (OPENER, {}, START)
    if _has(t, ["thank", "thanks", "merci"]):
        return ("You're welcome! Ask me anything else, or tell me your next goal.", ctx, START)
    STRONG = ["differs", "differ from", "different from", "different results", "thesis", "can i cite", "how do i cite", "how to cite", "cite this", "not working", "doesn't work", "doesnt work", "error", "crash", "patient data", "privacy", "what is this app", "what can you do"]
    for kws, ans in GLOSS:
        if any(k in STRONG and _has(t, [k]) for k in kws):
            return (ans, ctx, ["Make me a plan", "Start over"])
    # concept questions (explicit question forms win over plans)
    asking = t.strip().endswith("?") or _has(t, ["what is", "what's", "what does", "meaning", "explain", "how do i read", "why "])
    for kws, ans in GLOSS:
        if _has(t, kws) and (asking or len(t.split()) <= 5):
            return (ans, ctx, ["Make me a plan", "What is RMSD?", "Can I cite this?"])
    found = {k for k, (kws, _, _) in INTENTS.items() if _has(t, kws)}
    if _has(t, ["i have no data", "no data yet", "nothing yet", "i have no csv"]):
        found = {"data"}
    if found and not (_has(t, LOST_WORDS) and len(t.split()) <= 6):
        lines, dl = _plan(found)
        names = ", ".join(INTENTS[k][1] for k in ORDER if k in found)
        ans = "Here is your plan:\n\n" + "\n".join(f"{i}. {l}" for i, l in enumerate(lines, 1))
        if "dock" in found:
            ans += "\n\n⚠ Docking scores depend on preparation, box and seed: report all settings and validate your protein first."
        if dl:
            ctx["pending"] = "have_table"
            ans += "\n\nQuick question: **do you already have your data table (CSV/Excel)?** (yes / no)"
            if "omics" in found or "surv" in found:
                ans += " [[DL:" + ",".join(dl) + "]]"
        ctx["found"] = sorted(found)
        return (ans, ctx, ["What is FDR?", "Can I cite this?", "Start over"])
    if _has(t, LOST_WORDS + ["help", "start", "don't know", "dont know"]):
        ctx["onb"] = {"step": 1}
        q, sug = ONB_Q[1]
        return ("No problem, I'll guide you from zero. Three quick questions, then I'll build your roadmap.\n\n" + q, ctx, sug)
    for kws, ans in GLOSS:
        if _has(t, kws):
            return (ans, ctx, ["Make me a plan", "Start over"])
    if _has(t, ["plan"]) and ctx.get("found"):
        return respond(" ".join(INTENTS[k][0][0] for k in ctx["found"]), ctx)
    return ("I'm not sure I understood. Which is closest to what you need? You can also ask things like *'what is RMSD?'*, *'how do I read FDR?'*, *'can I cite this?'* or *'my docking differs from my thesis'*.", ctx, START)
