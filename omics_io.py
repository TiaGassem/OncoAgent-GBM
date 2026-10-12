"""Forgiving loaders for expression tables and group labels, with a plain-language report of
everything that was detected or changed. Never silently guesses: every decision is reported."""
from __future__ import annotations
import io
import numpy as np
import pandas as pd

GENE_HINTS = ("gene", "symbol", "hugo", "ensembl", "id", "name", "probe", "feature")
SAMPLE_HINTS = ("sample", "sample_id", "sampleid", "id", "name", "barcode", "patient", "run")
GROUP_HINTS = ("group", "condition", "class", "status", "type", "label", "phenotype", "arm")


def read_table(data: bytes, filename: str = "") -> pd.DataFrame:
    """CSV / TSV / TXT (delimiter sniffed) or XLSX."""
    if filename.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(io.BytesIO(data))
    text = data.decode("utf-8-sig", errors="replace")
    first = text.split("\n", 1)[0]
    sep = max([",", "\t", ";", "|"], key=first.count)
    return pd.read_csv(io.StringIO(text), sep=sep)


def load_expression(df: pd.DataFrame, samples_in_rows: bool | None = None):
    """Return (expr genes x samples, report lines). samples_in_rows=None -> auto-detect."""
    rep = []
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 3:
        raise ValueError(f"Only {df.shape[1]} column(s) found. Expected one gene-name column plus at least 2 sample columns. "
                         "Check that the file is a table (CSV/TSV) and not empty.")
    first = df.columns[0]
    gene_col = None
    if not pd.api.types.is_numeric_dtype(df[first]) or first.lower().startswith("unnamed") or any(h == first.lower() for h in GENE_HINTS):
        gene_col = first
    else:
        for c in df.columns:
            if any(h == c.lower() for h in GENE_HINTS) and not pd.api.types.is_numeric_dtype(df[c]):
                gene_col = c; break
    if gene_col is None:
        raise ValueError("Could not find a column of gene names. Put gene symbols (e.g. CDC25B) in the FIRST column.")
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df = df.set_index(gene_col)
    rep.append(f"Row labels taken from column '{gene_col}'.")

    num = df.apply(pd.to_numeric, errors="coerce")
    keep = [c for c in df.columns if num[c].notna().mean() >= 0.5]
    dropped = [c for c in df.columns if c not in keep]
    if dropped:
        rep.append(f"Dropped {len(dropped)} non-numeric column(s): {', '.join(map(str, dropped[:6]))}" + (" ..." if len(dropped) > 6 else ""))
    num = num[keep]
    if num.shape[1] < 2:
        raise ValueError("Fewer than 2 numeric columns remain. Is this a table of numbers (expression values)?")

    nr, nc = num.shape
    if samples_in_rows is None:
        samples_in_rows = nr < nc
        if samples_in_rows:
            rep.append(f"Table has {nr} rows x {nc} columns: fewer rows than columns, so I assumed SAMPLES are in rows and flipped it. "
                       "If that is wrong, tick 'my file has samples in rows' accordingly.")
    if samples_in_rows:
        num = num.T
    num.index = num.index.astype(str).str.strip()
    num.columns = [str(c).strip() for c in num.columns]

    dup = int(num.index.duplicated().sum())
    if dup:
        num = num.groupby(level=0).mean()
        rep.append(f"{dup} duplicated gene name(s) were averaged.")
    n_all_nan = int(num.isna().all(axis=1).sum())
    if n_all_nan:
        num = num.dropna(how="all"); rep.append(f"Removed {n_all_nan} gene row(s) with no values.")
    miss = float(num.isna().mean().mean())
    if miss > 0:
        rep.append(f"{miss:.1%} of values are missing; tests ignore missing values per gene.")
    if (num.min().min() < 0) and num.max().max() > 50:
        rep.append("WARNING: negative AND very large values. These do not look like log-scale expression; check your data.")
    elif num.max().max() > 100:
        rep.append("NOTE: values go above 100, so this looks like raw/linear values (counts, TPM, FPKM). "
                   "Tick 'values are linear' so log2(x+1) is applied. For raw counts prefer DESeq2/edgeR.")
    rep.insert(0, f"Detected {num.shape[0]:,} genes x {num.shape[1]:,} samples.")
    return num, rep


def _pick(df, hints, fallback_idx):
    for c in df.columns:
        if str(c).strip().lower() in hints:
            return c
    return df.columns[fallback_idx]


def load_groups(df: pd.DataFrame, sample_ids):
    """Return (groups Series indexed by sample_id as in expression, report lines)."""
    if df.shape[1] < 2:
        raise ValueError("Group file needs 2 columns: sample name and group (e.g. S1,tumor).")
    scol = _pick(df, SAMPLE_HINTS, 0)
    gcol = _pick(df.drop(columns=[scol]), GROUP_HINTS, 0)
    raw = pd.Series(df[gcol].astype(str).str.strip().values, index=df[scol].astype(str).str.strip().values)
    rep = [f"Group file: sample names from '{scol}', groups from '{gcol}' ({raw.nunique()} groups: {', '.join(sorted(raw.unique())[:6])})."]
    sample_ids = [str(s).strip() for s in sample_ids]
    exact = [s for s in sample_ids if s in raw.index]
    if len(exact) < len(sample_ids):
        norm = {k.upper().replace("-", ".").replace("_", "."): v for k, v in raw.items()}
        mapped = {s: norm.get(s.upper().replace("-", ".").replace("_", ".")) for s in sample_ids}
        if sum(v is not None for v in mapped.values()) > len(exact):
            out = pd.Series({s: v for s, v in mapped.items() if v is not None})
            rep.append(f"Matched {len(out)}/{len(sample_ids)} samples after ignoring case and -/_/. differences.")
            return out, _warn_unmatched(rep, sample_ids, out.index)
    out = raw[raw.index.isin(sample_ids)]
    rep.append(f"Matched {len(out)}/{len(sample_ids)} expression samples to the group file.")
    return out, _warn_unmatched(rep, sample_ids, out.index)


def _warn_unmatched(rep, sample_ids, matched):
    miss = [s for s in sample_ids if s not in set(matched)]
    if miss:
        rep.append(f"WARNING: {len(miss)} sample(s) have no group and are ignored, e.g. {', '.join(miss[:5])}.")
    return rep


def groups_from_rule(sample_ids, rules: dict):
    """rules = {'tumor': 'tumor,gbm', 'normal': 'normal,ctrl'}; case-insensitive 'contains' on sample names."""
    out = {}
    for s in sample_ids:
        hits = [g for g, kws in rules.items() if any(k.strip() and k.strip().lower() in str(s).lower() for k in kws.split(","))]
        if len(hits) == 1:
            out[s] = hits[0]
    return pd.Series(out, dtype=object)


def demo_dataset(seed: int = 7):
    """SIMULATED data (not real patients) so the tool can be tried: 1200 genes, 8 tumour + 8 normal.
    Built-in positive controls: CDC25B +3.5, CDC25A +2.5, EGFR +3 (up in tumour), MGMT -3 (down)."""
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i:04d}" for i in range(1195)] + ["CDC25A", "CDC25B", "EGFR", "MGMT", "CDC25C"]
    cols = [f"tumor_{i}" for i in range(1, 9)] + [f"normal_{i}" for i in range(1, 9)]
    X = pd.DataFrame(rng.normal(8, 0.6, (len(genes), 16)), index=genes, columns=cols)
    t = [c for c in cols if c.startswith("tumor")]
    for gene, shift in {"CDC25B": 3.5, "CDC25A": 2.5, "EGFR": 3.0, "MGMT": -3.0}.items():
        X.loc[gene, t] += shift
    groups = pd.Series(["tumor"] * 8 + ["normal"] * 8, index=X.columns)
    return X, groups


EXPRESSION_TEMPLATE = "gene,sample1,sample2,sample3,sample4\nCDC25A,8.1,8.4,7.2,7.0\nCDC25B,9.5,9.8,7.1,7.3\nEGFR,10.2,10.9,8.0,8.2\n"
GROUPS_TEMPLATE = "sample,group\nsample1,tumor\nsample2,tumor\nsample3,normal\nsample4,normal\n"
