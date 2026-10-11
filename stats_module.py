"""OncoAgent-GBM: survival and docking-enrichment statistics (pure NumPy/SciPy).

All functions compute on data the USER supplies. Nothing is looked up or invented.
Validated in tests against the classic Freireich leukaemia data (log-rank chi2 = 16.79).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.stats import chi2


@dataclass
class KMCurve:
    times: list
    survival: list
    n_at_risk: list


def kaplan_meier(time, event) -> KMCurve:
    t = np.asarray(time, float); e = np.asarray(event, int)
    if len(t) == 0 or len(t) != len(e):
        raise ValueError("time and event must be non-empty and of equal length")
    if np.any(t < 0) or not set(np.unique(e)) <= {0, 1}:
        raise ValueError("time must be >= 0 and event must be 0 (censored) or 1 (event)")
    ts, surv, risk = [0.0], [1.0], [len(t)]
    s = 1.0
    for u in np.unique(t[e == 1]):
        n = int(np.sum(t >= u)); d = int(np.sum((t == u) & (e == 1)))
        s *= 1 - d / n
        ts.append(float(u)); surv.append(s); risk.append(n)
    return KMCurve(ts, surv, risk)


@dataclass
class LogRank:
    chi2: float
    p_value: float
    observed: tuple
    expected: tuple
    n: tuple


def logrank_test(time, event, group) -> LogRank:
    """Two-group Mantel-Cox log-rank test (1 d.f.)."""
    t = np.asarray(time, float); e = np.asarray(event, int); g = np.asarray(group)
    labels = list(dict.fromkeys(g.tolist()))
    if len(labels) != 2:
        raise ValueError(f"exactly 2 groups required, got {len(labels)}")
    m = g == labels[0]
    O1 = E1 = V = 0.0
    for u in np.unique(t[e == 1]):
        at = t >= u
        n = at.sum(); n1 = (at & m).sum()
        d = ((t == u) & (e == 1)).sum(); d1 = ((t == u) & (e == 1) & m).sum()
        O1 += d1; E1 += d * n1 / n
        if n > 1:
            V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V == 0:
        raise ValueError("no variance: need events in both groups")
    stat = (O1 - E1) ** 2 / V
    O_tot = int(((e == 1)).sum())
    return LogRank(float(stat), float(chi2.sf(stat, 1)),
                   (int(O1), O_tot - int(O1)), (float(E1), O_tot - float(E1)),
                   (int(m.sum()), int((~m).sum())))


def roc_auc(scores_active, scores_decoy) -> float:
    """AUC where HIGHER score = predicted more active. (For Vina kcal/mol, pass -energy.)"""
    a = np.asarray(scores_active, float); d = np.asarray(scores_decoy, float)
    if len(a) == 0 or len(d) == 0:
        raise ValueError("need at least one active and one decoy")
    wins = (a[:, None] > d[None, :]).sum() + 0.5 * (a[:, None] == d[None, :]).sum()
    return float(wins / (len(a) * len(d)))


def enrichment_factor(scores_active, scores_decoy, fraction: float = 0.01) -> float:
    """EF at the top `fraction` of the ranked list (higher score = better)."""
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    a = np.asarray(scores_active, float); d = np.asarray(scores_decoy, float)
    allv = np.concatenate([a, d]); lab = np.concatenate([np.ones(len(a)), np.zeros(len(d))])
    order = np.argsort(-allv, kind="stable")
    k = max(1, int(np.ceil(fraction * len(allv))))
    hits = lab[order][:k].sum()
    return float((hits / k) / (len(a) / len(allv)))
