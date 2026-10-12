"""Streamlit tab: omics, compound and structure tools. Works only on data the user uploads/types."""
import io
import pandas as pd
import streamlit as st
import omics_module as om

TITLES = {  # short labels translated; long help text is English (see README)
    "en": ["Omics & Tools", "Differential expression", "Gene-set enrichment", "Co-expression",
           "Compound tools", "Structure tools", "Resources"],
    "fr": ["Omiques & outils", "Expression differentielle", "Enrichissement de jeux de genes", "Co-expression",
           "Outils composes", "Outils structure", "Ressources"],
    "ar": ["الأوميكس والأدوات", "التعبير التفاضلي", "إثراء مجموعات الجينات", "التعبير المشترك",
           "أدوات المركبات", "أدوات البنية", "مصادر"],
}
SCOPE = ("Exploratory analysis of YOUR log-scale expression table (genes x samples, first column = gene). "
         "Not a substitute for DESeq2/edgeR/limma on raw counts. Nothing is downloaded or invented.")


HELP = """
**What is this?** You give a table of gene-expression numbers (one row per gene, one column per sample) and
say which samples are in which group (e.g. *tumor* vs *normal*). The app tests every gene for a difference.

**Your expression file** (CSV, TSV, TXT or Excel): first column = gene names (e.g. `CDC25B`), other columns = samples
with numbers. Samples-in-rows files are detected and flipped automatically. Download the template below.

**Your group file**: two columns, `sample` and `group`. Or type a rule such as: names containing "tumor" = tumor.

**How to read the result.** `log2FC` = difference between groups on the log2 scale (+1 means 2x higher in group A).
`p_value` = chance of a difference this large if there were no real difference. `fdr_bh` = the p-value corrected for
testing thousands of genes at once: **use fdr_bh < 0.05, not p_value.**

**Sample size matters.** With 3 samples per group almost nothing survives correction, even if real. Aim for 6-8+ per
group, report n, and never claim significance from p_value alone.

**Limits.** Expects log-scale values (log2 of microarray/RNA-seq normalised data). For raw RNA-seq counts use DESeq2 /
edgeR / limma-voom. Exploratory analysis only; validate hits in independent cohorts and in the lab.
"""


def _show_report(lines):
    for ln in lines:
        (st.warning if ln.startswith(("WARNING", "NOTE")) else st.info)(ln)


def _get_inputs(key, need_groups=True):
    """Returns (expr, groups or None, is_demo). Always tells the user what it understood."""
    import omics_io as oio
    mode = st.radio("Data", ["Try simulated demo data", "Use my own files"], horizontal=True, key=f"{key}_mode")
    if mode.startswith("Try"):
        X, g = oio.demo_dataset()
        st.warning("SIMULATED data (not real patients): 1,200 genes, 8 tumor vs 8 normal, with CDC25B, CDC25A, EGFR up and "
                   "MGMT down on purpose. Use it to see how the tool works, then switch to your own files.")
        return X, g, True
    d1, d2 = st.columns(2)
    d1.download_button("Template: expression CSV", oio.EXPRESSION_TEMPLATE, "expression_template.csv", key=f"{key}_t1")
    d2.download_button("Template: groups CSV", oio.GROUPS_TEMPLATE, "groups_template.csv", key=f"{key}_t2")
    ef = st.file_uploader("Expression table (CSV / TSV / TXT / XLSX)", type=["csv", "tsv", "txt", "xlsx", "xls"], key=f"{key}_e")
    orient = st.selectbox("Layout", ["Auto-detect", "Genes in rows (samples in columns)", "Samples in rows (genes in columns)"], key=f"{key}_o")
    if not ef:
        st.caption("Upload your expression table to continue.")
        return None, None, False
    try:
        df = oio.read_table(ef.getvalue(), ef.name)
        expr, rep = oio.load_expression(df, None if orient.startswith("Auto") else orient.startswith("Samples"))
    except Exception as e:
        st.error(f"I could not read your expression file: {e}")
        return None, None, False
    with st.expander("What I understood from your expression file", expanded=True):
        _show_report(rep)
        st.dataframe(expr.iloc[:5, :6])
    if not need_groups:
        return expr, None, False
    how = st.radio("Groups", ["Upload group file", "Type a rule from sample names"], horizontal=True, key=f"{key}_how")
    groups = None
    if how.startswith("Upload"):
        gf = st.file_uploader("Group file (CSV: sample, group)", type=["csv", "tsv", "txt", "xlsx"], key=f"{key}_g")
        if gf:
            try:
                groups, rep2 = oio.load_groups(oio.read_table(gf.getvalue(), gf.name), expr.columns)
                with st.expander("What I understood from your group file", expanded=True):
                    _show_report(rep2)
            except Exception as e:
                st.error(f"I could not read your group file: {e}")
    else:
        st.caption(f"Your sample names look like: {', '.join(map(str, expr.columns[:5]))} ...")
        c1, c2 = st.columns(2)
        na = c1.text_input("Group A name", "tumor", key=f"{key}_na"); ka = c1.text_input("A if name contains (comma separated)", "tumor", key=f"{key}_ka")
        nb = c2.text_input("Group B name", "normal", key=f"{key}_nb"); kb = c2.text_input("B if name contains", "normal", key=f"{key}_kb")
        groups = oio.groups_from_rule(expr.columns, {na: ka, nb: kb})
        st.info(f"{(groups == na).sum()} samples in {na}, {(groups == nb).sum()} in {nb}, "
                f"{len(expr.columns) - len(groups)} unassigned (ignored).")
    return expr, groups, False


def tab_analysis(lang="en"):
    T = TITLES.get(lang, TITLES["en"])
    st.subheader(T[0])
    st.caption(SCOPE)
    with st.expander("How this works, in plain words (read first)", expanded=False):
        st.markdown(HELP)
    t_de, t_ora, t_co, t_cmp, t_str, t_res = st.tabs(T[1:])

    with t_de:
        expr, groups, demo = _get_inputs("de")
        if expr is not None and groups is not None and len(groups):
            labels = sorted(set(groups))
            if len(labels) < 2:
                st.error(f"Only {len(labels)} group found ({', '.join(labels)}). Two groups are needed.")
            else:
                a = st.selectbox("Group A", labels, key="de_a")
                b = st.selectbox("Group B (reference)", [x for x in labels if x != a], key="de_b")
                method = st.radio("Test", ["welch", "mannwhitney"], horizontal=True, key="de_m",
                                  help="Welch t-test is the usual choice; Mann-Whitney is rank-based and needs more samples.")
                lin = st.checkbox("Values are linear (TPM/FPKM): apply log2(x+1)", key="de_lin")
                if st.button("Run differential expression", key="de_run", type="primary"):
                    try:
                        gs = groups[groups.isin([a, b])]
                        res = om.differential_expression(expr[[c for c in expr.columns if c in gs.index]], gs, a, b, method, lin)
                        na_, nb_ = int((gs == a).sum()), int((gs == b).sum())
                        sig = res["fdr_bh"] < 0.05
                        st.success(f"{na_} {a} vs {nb_} {b} samples, {len(res):,} genes tested: {int(sig.sum())} with FDR < 0.05.")
                        if min(na_, nb_) < 4:
                            st.warning("Fewer than 4 samples in a group: results are very unreliable.")
                        st.dataframe(res.head(200), use_container_width=True)
                        st.download_button("Download full table (CSV)", res.to_csv(index=False), "de_results.csv")
                        import matplotlib
                        matplotlib.use("Agg")
                        import matplotlib.pyplot as plt
                        import numpy as np
                        fig, ax = plt.subplots(figsize=(6.5, 4.8))
                        ax.scatter(res["log2FC"], -np.log10(res["p_value"]), s=7, c=np.where(sig, "#dc2626", "#94a3b8"))
                        for _, r in res[sig].head(12).iterrows():
                            ax.annotate(r["gene"], (r["log2FC"], -np.log10(r["p_value"])), fontsize=8)
                        ax.set_xlabel(f"log2 fold change ({a} - {b})"); ax.set_ylabel("-log10 p-value")
                        ax.set_title("Volcano plot (red = FDR < 0.05)"); st.pyplot(fig)
                        sc, ev = om.pca(expr[[c for c in expr.columns if c in gs.index]])
                        fig2, ax2 = plt.subplots(figsize=(5, 4))
                        for lab in (a, b):
                            idx = [x for x in sc.index if gs.get(x) == lab]
                            ax2.scatter(sc.loc[idx, "PC1"], sc.loc[idx, "PC2"], label=lab)
                        ax2.set_xlabel(f"PC1 ({ev[0]:.0%})"); ax2.set_ylabel(f"PC2 ({ev[1]:.0%})"); ax2.legend()
                        ax2.set_title("PCA of samples: do the groups separate?"); st.pyplot(fig2)
                        if demo:
                            st.info("Demo check: CDC25B, CDC25A, EGFR (up) and MGMT (down) were planted, so they should appear at the top.")
                    except Exception as e:
                        st.error(f"Could not analyse: {e}")

    with t_ora:
        st.caption("Paste hit genes (one per line) and upload a gene-set file in GMT format "
                   "(e.g. MSigDB Hallmark or Reactome, downloaded from their sites; check each licence).")
        hits = st.text_area("Hit genes", height=120, key="ora_hits")
        bgf = st.file_uploader("Background genes (text, one per line) — optional but recommended", type=["txt"], key="ora_bg")
        gmt = st.file_uploader("Gene sets (.gmt)", type=["gmt", "txt"], key="ora_gmt")
        if gmt and hits.strip() and st.button("Run enrichment", key="ora_run"):
            try:
                sets = om.parse_gmt(gmt.read().decode("utf-8", "ignore"))
                hit_list = [x.strip() for x in hits.split() if x.strip()]
                bg = ([x.strip() for x in bgf.read().decode().split() if x.strip()] if bgf
                      else sorted(set().union(*sets.values())))
                if not bgf:
                    st.warning("No background given: using all genes in the gene-set file (can bias results).")
                df = om.over_representation(hit_list, bg, sets)
                st.dataframe(df.head(100))
            except Exception as e:
                st.error(f"Could not analyse: {e}")

    with t_co:
        expr, _, demo = _get_inputs("co", need_groups=False)
        if expr is not None:
            gene = st.text_input("Gene to correlate against all others (e.g. CDC25B)", "CDC25B" if demo else "", key="co_g")
            meth = st.radio("Method", ["spearman", "pearson"], horizontal=True, key="co_m")
            if gene and st.button("Correlate", key="co_run", type="primary"):
                try:
                    q = gene.strip()
                    if q not in expr.index:
                        hits = [g for g in expr.index if g.upper() == q.upper()]
                        q = hits[0] if hits else q
                    df = om.correlate_gene(expr, q, meth)
                    st.dataframe(df.head(100), use_container_width=True)
                    st.caption(f"{len(expr.columns)} samples. Correlation is not causation; with few samples p-values are unstable.")
                except Exception as e:
                    st.error(f"Could not analyse: {e}")

    with t_cmp:
        import compound_tools as ct
        smi = st.text_area("SMILES (one per line)", height=120, key="ct_smi",
                           help="Structural alerts (PAINS, Brenk) are warnings, not proof of assay interference.")
        cutoff = st.slider("Butina cluster distance cutoff", 0.1, 0.9, 0.4, 0.05, key="ct_cut")
        if st.button("Analyse compounds", key="ct_run"):
            lines = [x.strip() for x in smi.splitlines() if x.strip()]
            rows = []
            for s in lines:
                try:
                    al = ct.structural_alerts(s)
                    rows.append({"smiles": s, "QED": round(ct.qed_score(s), 3), "scaffold": ct.murcko_scaffold(s),
                                 "alerts": "; ".join(f"{c}:{d}" for c, d in al) or "none"})
                except ValueError as e:
                    rows.append({"smiles": s, "QED": None, "scaffold": "", "alerts": str(e)})
            if rows:
                st.dataframe(pd.DataFrame(rows))
            ok = [r["smiles"] for r in rows if r["QED"] is not None]
            if len(ok) >= 2:
                st.write("Butina clusters (indices):", ct.butina_clusters(ok, cutoff))

    with t_str:
        import structure_tools as stt
        mode = st.radio("Analysis", ["Ligand contacts", "AlphaFold pLDDT"], horizontal=True, key="st_mode")
        if mode == "Ligand contacts":
            pf = st.file_uploader("Protein PDB", type=["pdb"], key="st_p")
            lf = st.file_uploader("Ligand (PDB or PDBQT, same coordinate frame)", type=["pdb", "pdbqt"], key="st_l")
            cut = st.slider("Contact cutoff (Å)", 3.0, 6.0, 4.0, 0.1, key="st_c")
            if pf and lf and st.button("Find contacts", key="st_run"):
                try:
                    res, hb = stt.ligand_contacts(pf.read().decode(), lf.read().decode(), cut)
                    st.dataframe(pd.DataFrame(res))
                    st.write("H-bond candidates (N/O pairs ≤ 3.5 Å; distance only, no angle check):")
                    st.dataframe(pd.DataFrame(hb))
                except Exception as e:
                    st.error(f"Could not analyse: {e}")
        else:
            af = st.file_uploader("AlphaFold model PDB", type=["pdb"], key="st_af")
            if af and st.button("Summarise pLDDT", key="st_pl"):
                try:
                    st.json(stt.plddt_summary(af.read().decode()))
                    st.caption("pLDDT bands: >90 very high, 70-90 confident, 50-70 low, <50 very low. "
                               "Not valid for experimental structures.")
                except Exception as e:
                    st.error(f"Could not analyse: {e}")

    with t_res:
        st.markdown(
            "**Not hosted here** (need GPUs or specialised servers; use the originals):\n\n"
            "- Structure prediction: [AlphaFold DB](https://alphafold.ebi.ac.uk), [ESM Atlas](https://esmatlas.com)\n"
            "- ADME / toxicity web servers: [SwissADME](https://www.swissadme.ch), [ProTox-3](https://tox.charite.de/protox3/)\n"
            "- Target evidence: [Open Targets](https://platform.opentargets.org), [DepMap](https://depmap.org/portal/)\n"
            "- Cohorts: [cBioPortal](https://www.cbioportal.org), [CGGA](http://www.cgga.org.cn), "
            "[GlioVis](http://gliovis.bioinfo.cnio.es)\n\nCheck each site's terms before use.")
