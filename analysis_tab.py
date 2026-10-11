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


def _read_csv(f):
    return pd.read_csv(f, index_col=0)


def _expr_and_groups(key):
    c1, c2 = st.columns(2)
    ef = c1.file_uploader("Expression CSV (genes x samples)", type=["csv"], key=f"{key}_e")
    gf = c2.file_uploader("Groups CSV (columns: sample, group)", type=["csv"], key=f"{key}_g")
    if not (ef and gf):
        return None, None
    expr = _read_csv(ef)
    g = pd.read_csv(gf)
    groups = pd.Series(g["group"].astype(str).values, index=g["sample"].astype(str).values)
    expr.columns = expr.columns.astype(str)
    return expr, groups


def tab_analysis(lang="en"):
    T = TITLES.get(lang, TITLES["en"])
    st.subheader(T[0])
    st.caption(SCOPE)
    t_de, t_ora, t_co, t_cmp, t_str, t_res = st.tabs(T[1:])

    with t_de:
        expr, groups = _expr_and_groups("de")
        if expr is not None:
            labels = sorted(set(groups))
            a = st.selectbox("Group A", labels, key="de_a")
            b = st.selectbox("Group B (reference)", labels, index=min(1, len(labels) - 1), key="de_b")
            method = st.radio("Test", ["welch", "mannwhitney"], horizontal=True, key="de_m")
            lin = st.checkbox("Values are linear (TPM/FPKM): apply log2(x+1)", key="de_lin")
            if st.button("Run differential expression", key="de_run"):
                try:
                    res = om.differential_expression(expr, groups, a, b, method, lin)
                    st.dataframe(res.head(200))
                    st.download_button("Download full table (CSV)", res.to_csv(index=False), "de_results.csv")
                    import matplotlib
                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt
                    import numpy as np
                    fig, ax = plt.subplots(figsize=(6, 4.5))
                    sig = res["fdr_bh"] < 0.05
                    ax.scatter(res["log2FC"], -np.log10(res["p_value"]), s=6, c=np.where(sig, "#dc2626", "#94a3b8"))
                    ax.set_xlabel(f"log2FC ({a} - {b})"); ax.set_ylabel("-log10 p"); ax.set_title("Volcano (red: FDR < 0.05)")
                    st.pyplot(fig)
                    st.caption(f"{int(sig.sum())} genes with BH-FDR < 0.05 of {len(res)} tested. "
                               "Check sample size and batch effects before interpreting.")
                    sc, ev = om.pca(expr)
                    fig2, ax2 = plt.subplots(figsize=(5, 4))
                    for lab in labels:
                        idx = [s for s in sc.index if groups.get(s) == lab]
                        ax2.scatter(sc.loc[idx, "PC1"], sc.loc[idx, "PC2"], label=lab)
                    ax2.set_xlabel(f"PC1 ({ev[0]:.0%})"); ax2.set_ylabel(f"PC2 ({ev[1]:.0%})"); ax2.legend()
                    st.pyplot(fig2)
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
        ef = st.file_uploader("Expression CSV (genes x samples)", type=["csv"], key="co_e")
        gene = st.text_input("Gene (e.g. CDC25B)", key="co_g")
        meth = st.radio("Method", ["spearman", "pearson"], horizontal=True, key="co_m")
        if ef and gene and st.button("Correlate", key="co_run"):
            try:
                df = om.correlate_gene(_read_csv(ef), gene.strip(), meth)
                st.dataframe(df.head(100))
                st.caption("Correlation is not causation; with few samples, p-values are unstable.")
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
