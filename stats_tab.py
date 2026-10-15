"""Streamlit tab: survival (Kaplan-Meier + log-rank) and docking enrichment (AUC / EF).
Runs only on data the user uploads. UI text: English / French / Arabic."""
import io
import pandas as pd
import streamlit as st
from stats_module import kaplan_meier, logrank_test, roc_auc, enrichment_factor

TXT = {
    "en": {
        "title": "Statistics: survival and docking enrichment",
        "intro": "Upload YOUR data. Nothing is looked up or invented. Research use only.",
        "surv": "Survival (Kaplan-Meier + log-rank)",
        "surv_help": "CSV with columns: time, event (1 = event, 0 = censored), group (exactly 2 groups, e.g. high/low expression).",
        "enr": "Docking enrichment (actives vs decoys)",
        "enr_help": "CSV with columns: score, label (1 = known active, 0 = decoy). For Vina energies use 'lower is better'.",
        "upload": "Upload CSV", "run": "Run",
        "lower": "Lower score is better (e.g. Vina kcal/mol)", "frac": "Top fraction for EF",
        "chi": "Log-rank chi-square (1 d.f.)", "p": "p-value", "obs": "Observed events",
        "exp": "Expected events", "auc": "ROC AUC (0.5 = random)", "ef": "Enrichment factor",
        "warn_n": "Small samples give unreliable p-values; report the number of events.",
        "warn_dec": "Enrichment is only meaningful with property-matched decoys (e.g. DUD-E) and many actives.",
        "err": "Could not analyse the file",
    },
    "fr": {
        "title": "Statistiques : survie et enrichissement du docking",
        "intro": "Chargez VOS donnees. Rien n'est invente. Usage de recherche uniquement.",
        "surv": "Survie (Kaplan-Meier + log-rank)",
        "surv_help": "CSV avec colonnes : time, event (1 = evenement, 0 = censure), group (exactement 2 groupes).",
        "enr": "Enrichissement du docking (actifs vs leurres)",
        "enr_help": "CSV avec colonnes : score, label (1 = actif connu, 0 = leurre).",
        "upload": "Charger un CSV", "run": "Lancer",
        "lower": "Un score plus bas est meilleur (ex. Vina kcal/mol)", "frac": "Fraction superieure pour l'EF",
        "chi": "Chi-deux du log-rank (1 ddl)", "p": "valeur p", "obs": "Evenements observes",
        "exp": "Evenements attendus", "auc": "AUC ROC (0,5 = hasard)", "ef": "Facteur d'enrichissement",
        "warn_n": "Les petits echantillons donnent des valeurs p peu fiables ; indiquez le nombre d'evenements.",
        "warn_dec": "L'enrichissement n'a de sens qu'avec des leurres apparies (ex. DUD-E) et de nombreux actifs.",
        "err": "Analyse impossible",
    },
    "ar": {
        "title": "الإحصاء: البقاء وإثراء الالتحام الجزيئي",
        "intro": "ارفع بياناتك أنت. لا يتم اختلاق أي شيء. للاستخدام البحثي فقط.",
        "surv": "البقاء (كابلان-ماير + اختبار لوغ-رانك)",
        "surv_help": "ملف CSV بالأعمدة: time و event (1 = حدث، 0 = مراقَب) و group (مجموعتان بالضبط).",
        "enr": "إثراء الالتحام (النشطة مقابل الشراك)",
        "enr_help": "ملف CSV بالأعمدة: score و label (1 = نشط معروف، 0 = شرَك).",
        "upload": "رفع ملف CSV", "run": "تشغيل",
        "lower": "الدرجة الأقل أفضل (مثل Vina kcal/mol)", "frac": "الجزء الأعلى لحساب EF",
        "chi": "مربع كاي لاختبار لوغ-رانك (درجة حرية 1)", "p": "القيمة الاحتمالية p", "obs": "الأحداث المرصودة",
        "exp": "الأحداث المتوقعة", "auc": "AUC (0.5 = عشوائي)", "ef": "عامل الإثراء",
        "warn_n": "العينات الصغيرة تعطي قيم p غير موثوقة؛ اذكر عدد الأحداث.",
        "warn_dec": "الإثراء ذو معنى فقط مع شراك مطابقة الخصائص (مثل DUD-E) وعدد كبير من النشطة.",
        "err": "تعذر تحليل الملف",
    },
}


def _input(key, t, kind):
    from guides import sample_files
    mode = st.radio("Data", ["Try example data", "Use my own CSV"], horizontal=True, key=f"{key}_mode")
    if mode.startswith("Try"):
        lab, txt, _ = sample_files()[kind]
        st.warning(lab + ". For learning only.")
        return pd.read_csv(io.StringIO(txt))
    f = st.file_uploader(t["upload"], type=["csv"], key=f"{key}_csv")
    return pd.read_csv(f) if f else None


def tab_statistics(lang: str = "en"):
    t = TXT.get(lang, TXT["en"])
    st.subheader(t["title"])
    st.caption(t["intro"])
    s1, s2 = st.tabs([t["surv"], t["enr"]])

    with s1:
        st.caption(t["surv_help"])
        df0 = _input("surv", t, "survival")
        if df0 is not None and st.button(t["run"], key="surv_run"):
            try:
                df = df0.dropna(subset=["time", "event", "group"])
                lr = logrank_test(df["time"], df["event"], df["group"])
                c1, c2, c3 = st.columns(3)
                c1.metric(t["chi"], f"{lr.chi2:.2f}")
                c2.metric(t["p"], f"{lr.p_value:.2e}")
                c3.metric("n", f"{lr.n[0]} / {lr.n[1]}")
                st.write(f"{t['obs']}: {lr.observed} | {t['exp']}: ({lr.expected[0]:.2f}, {lr.expected[1]:.2f})")
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(6, 4))
                for gname, sub in df.groupby("group"):
                    km = kaplan_meier(sub["time"], sub["event"])
                    ax.step(km.times, km.survival, where="post", label=f"{gname} (n={len(sub)})")
                ax.set_ylim(0, 1.02); ax.set_xlabel("time"); ax.set_ylabel("S(t)"); ax.legend()
                st.pyplot(fig)
                st.warning(t["warn_n"])
            except Exception as e:  # shown to the user, never hidden
                st.error(f"{t['err']}: {e}")

    with s2:
        st.caption(t["enr_help"])
        lower = st.checkbox(t["lower"], value=True, key="enr_lower")
        frac = st.slider(t["frac"], 0.01, 0.5, 0.05, 0.01, key="enr_frac")
        df0 = _input("enr", t, "enrichment")
        if df0 is not None and st.button(t["run"], key="enr_run"):
            try:
                df = df0.dropna(subset=["score", "label"])
                sc = -df["score"] if lower else df["score"]
                a = sc[df["label"] == 1]; d = sc[df["label"] == 0]
                c1, c2 = st.columns(2)
                c1.metric(t["auc"], f"{roc_auc(a, d):.3f}")
                c2.metric(f"{t['ef']} @ {frac:.0%}", f"{enrichment_factor(a, d, frac):.2f}")
                st.caption(f"actives = {len(a)}, decoys = {len(d)}")
                st.warning(t["warn_dec"])
            except Exception as e:
                st.error(f"{t['err']}: {e}")
