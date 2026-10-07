"""4-Parameter Logistic (4PL) dose-response / IC50 curve fitting.

This module fits the user's OWN in-vitro dose-response data. It computes
nothing unless the user provides (dose, response) pairs, and it never invents
or stores results. The model is the standard sigmoidal 4PL:

    y = d + (a - d) / (1 + (x / c) ** b)

where
    a = bottom (response at zero dose / minimum asymptote),
    d = top (maximum asymptote),
    c = inflection point = IC50 (dose at half-maximal response),
    b = Hill slope.

Returns IC50, Hill slope, the four fitted parameters, R^2 and the fitted curve
so the UI can plot it. If the fit cannot converge, it says so honestly rather
than returning a fabricated number.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except Exception:  # pragma: no cover
    HAS_SCIPY = False


def fourpl(x, a, b, c, d):
    """4PL model. x must be > 0 (dose). c is IC50."""
    x = np.asarray(x, dtype=float)
    # guard against divide-by-zero / negative base with fractional exponent
    ratio = np.where(c > 0, x / c, 0.0)
    return d + (a - d) / (1.0 + np.power(ratio, b))


@dataclass
class FitResult:
    ok: bool = False
    error: str = ""
    ic50: float = 0.0
    hill_slope: float = 0.0
    bottom: float = 0.0
    top: float = 0.0
    r_squared: float = 0.0
    ic50_in_range: bool = True
    n_points: int = 0
    x_data: List[float] = field(default_factory=list)
    y_data: List[float] = field(default_factory=list)
    x_curve: List[float] = field(default_factory=list)
    y_curve: List[float] = field(default_factory=list)


def parse_pairs(text: str) -> Tuple[List[float], List[float], str]:
    """Parse pasted dose-response data.

    Accepts one pair per line, dose and response separated by comma / tab /
    space, e.g.:
        0.1, 98
        1    82
        10   45
        100  8
    Returns (doses, responses, error_message). On any bad line it reports the
    offending line instead of guessing.
    """
    doses: List[float] = []
    resp: List[float] = []
    for i, raw in enumerate(text.strip().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p for p in line.replace(",", " ").replace("\t", " ").split() if p]
        if len(parts) < 2:
            return [], [], f"Line {i}: need two numbers (dose and response), got: '{raw}'"
        try:
            x = float(parts[0]); y = float(parts[1])
        except ValueError:
            return [], [], f"Line {i}: could not read two numbers from: '{raw}'"
        if x <= 0:
            return [], [], (f"Line {i}: dose must be > 0 for a log-scale IC50 fit "
                            f"(got {x}). Drop the zero-dose/control row or use a "
                            f"small positive dose.")
        doses.append(x); resp.append(y)
    return doses, resp, ""


def fit_4pl(doses: List[float], responses: List[float]) -> FitResult:
    """Fit the 4PL model to the user's data. No fabrication: if it cannot fit,
    `ok=False` and `error` explains why."""
    res = FitResult()
    if not HAS_SCIPY:
        res.error = "scipy is not installed on the server (add scipy to requirements.txt)."
        return res
    x = np.asarray(doses, dtype=float)
    y = np.asarray(responses, dtype=float)
    res.n_points = int(x.size)
    if x.size < 4:
        res.error = (f"Need at least 4 dose points to fit a 4-parameter model; "
                     f"you gave {x.size}.")
        return res

    # Reasonable starting guesses from the data (no hidden assumptions).
    a0 = float(np.max(y))          # top / response at low dose
    d0 = float(np.min(y))          # bottom / response at high dose
    c0 = float(np.median(x))       # IC50 guess = middle dose
    b0 = 1.0                       # Hill slope guess
    try:
        popt, _ = curve_fit(
            fourpl, x, y, p0=[a0, b0, c0, d0], maxfev=20000,
            bounds=([-np.inf, -np.inf, 1e-9, -np.inf],
                    [np.inf, np.inf, np.inf, np.inf]),
        )
    except Exception as e:
        res.error = f"Curve fit did not converge: {e}. Check your data (need a clear sigmoid)."
        return res

    a, b, c, d = [float(v) for v in popt]
    y_pred = fourpl(x, a, b, c, d)
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    res.ok = True
    res.bottom = round(d, 4)
    res.top = round(a, 4)
    res.hill_slope = round(b, 4)
    res.ic50 = round(c, 6)
    res.r_squared = round(r2, 4)
    res.ic50_in_range = bool(np.min(x) <= c <= np.max(x))
    res.x_data = [float(v) for v in x]
    res.y_data = [float(v) for v in y]
    # smooth curve across the tested dose range (log-spaced)
    xc = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 200)
    res.x_curve = [float(v) for v in xc]
    res.y_curve = [float(v) for v in fourpl(xc, a, b, c, d)]
    return res
