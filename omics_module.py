"""OncoAgent-GBM: lightweight omics statistics on USER-supplied tables (NumPy/SciPy/pandas).

Scope (be honest): exploratory analysis of log-scale expression matrices. For raw RNA-seq counts use
DESeq2/edgeR/limma-voom; this module does NOT model count data. Nothing is downloaded or invented.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN preserved)."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n); adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def differential_expression(expr: pd.DataFrame, groups: pd.Series, a: str, b: str,
                            method: str = "welch", log2_transform: bool = False) -> pd.DataFrame:
    """expr: genes x samples (log-scale unless log2_transform=True). Contrast = mean(a) - mean(b).
    method: 'welch' (t-test, unequal variance) or 'mannwhitney'."""
    if method not in ("welch", "mannwhitney"):
        raise ValueError("method must be 'welch' or 'mannwhitney'")
    sa = [s for s in expr.columns if groups.get(s) == a]
    sb = [s for s in expr.columns if groups.get(s) == b]
    if len(sa) < 2 or len(sb) < 2:
        raise ValueError(f"need >= 2 samples per group (got {len(sa)} and {len(sb)})")
    X = expr.astype(float)
    if log2_transform:
        if (X < 0).any().any():
            raise ValueError("negative values cannot be log2-transformed")
        X = np.log2(X + 1)
    A, B = X[sa].to_numpy(), X[sb].to_numpy()
    if method == "welch":
        _, p = stats.ttest_ind(A, B, axis=1, equal_var=False, nan_policy="omit")
    else:
        p = np.array([stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
                      if np.isfinite(x).sum() and np.isfinite(y).sum() else np.nan
                      for x, y in zip(A, B)])
    res = pd.DataFrame({"gene": X.index, "log2FC": A.mean(axis=1) - B.mean(axis=1),
                        "mean_a": A.mean(axis=1), "mean_b": B.mean(axis=1), "p_value": p})
    res["fdr_bh"] = bh_fdr(res["p_value"].to_numpy())
    return res.sort_values("p_value", na_position="last").reset_index(drop=True)


def pca(expr: pd.DataFrame, n_components: int = 2, scale: bool = True):
    """PCA of samples (columns). Returns (scores DataFrame samples x PCs, explained variance ratio)."""
    X = expr.astype(float).T.dropna(axis=1)
    X = X.loc[:, X.std() > 0]
    Z = (X - X.mean()) / (X.std(ddof=1) if scale else 1.0)
    U, S, _ = np.linalg.svd(Z.to_numpy(), full_matrices=False)
    k = min(n_components, len(S))
    ev = S ** 2 / np.sum(S ** 2)
    scores = pd.DataFrame(U[:, :k] * S[:k], index=X.index, columns=[f"PC{i+1}" for i in range(k)])
    return scores, ev[:k]


def parse_gmt(text: str) -> dict:
    """GMT: name<TAB>description<TAB>gene1<TAB>gene2..."""
    sets = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[0]:
            sets[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    return sets


def over_representation(hits, background, gene_sets: dict, min_size: int = 3) -> pd.DataFrame:
    """One-sided hypergeometric test (enrichment). Gene sets are intersected with the background."""
    bg = {g.upper() for g in background}
    hit = {g.upper() for g in hits} & bg
    if not hit:
        raise ValueError("none of the hit genes are in the background")
    N, n = len(bg), len(hit)
    rows = []
    for name, members in gene_sets.items():
        m = {g.upper() for g in members} & bg
        if len(m) < min_size:
            continue
        k = len(m & hit)
        p = stats.hypergeom.sf(k - 1, N, len(m), n) if k > 0 else 1.0
        rows.append({"gene_set": name, "set_size": len(m), "overlap": k,
                     "p_value": float(p), "genes": ",".join(sorted(m & hit))})
    if not rows:
        return pd.DataFrame(columns=["gene_set", "set_size", "overlap", "p_value", "fdr_bh", "genes"])
    df = pd.DataFrame(rows)
    df["fdr_bh"] = bh_fdr(df["p_value"].to_numpy())
    return df.sort_values("p_value").reset_index(drop=True)


def correlate_gene(expr: pd.DataFrame, gene: str, method: str = "spearman") -> pd.DataFrame:
    """Correlation of one gene with all others across samples, BH-adjusted."""
    if gene not in expr.index:
        raise ValueError(f"{gene} not found")
    X = expr.astype(float)
    y = X.loc[gene].to_numpy()
    fn = stats.spearmanr if method == "spearman" else stats.pearsonr
    rows = []
    for g, row in X.iterrows():
        if g == gene:
            continue
        v = row.to_numpy()
        ok = np.isfinite(v) & np.isfinite(y)
        if ok.sum() < 4 or np.std(v[ok]) == 0 or np.std(y[ok]) == 0:
            continue
        r, p = fn(y[ok], v[ok])
        rows.append((g, float(r), float(p), int(ok.sum())))
    df = pd.DataFrame(rows, columns=["gene", "r", "p_value", "n"])
    df["fdr_bh"] = bh_fdr(df["p_value"].to_numpy())
    return df.sort_values("p_value").reset_index(drop=True)
