"""Run cards: a small JSON record that makes an analysis re-runnable and citable.

Pure Python (no Streamlit). Contains no patient data by design: callers must only
pass method parameters and public identifiers.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from importlib import metadata

SCHEMA = "oncoagent-gbm.runcard/1"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def package_versions(names=("rdkit", "numpy", "scipy", "pandas", "streamlit", "meeko", "requests")) -> dict:
    out = {}
    for n in names:
        try:
            out[n] = metadata.version(n)
        except metadata.PackageNotFoundError:
            out[n] = "not installed"
    return out


def build_runcard(kind: str, inputs: dict, parameters: dict | None = None,
                  results_summary: dict | None = None, app_version: str = "unknown",
                  data_sources: list[dict] | None = None, input_texts: dict | None = None) -> dict:
    """Assemble a run card. ``input_texts`` maps label -> raw text to be hashed (not stored)."""
    hashes = {k: sha256_text(v) for k, v in (input_texts or {}).items()}
    return {
        "schema": SCHEMA,
        "kind": kind,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "app_version": app_version,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": package_versions(),
        "inputs": inputs,
        "input_sha256": hashes,
        "parameters": parameters or {},
        "results_summary": results_summary or {},
        "data_sources": data_sources or [],
        "notice": "Research and education only. Not a medical device. "
                  "Computational outputs are hypotheses, not evidence of efficacy.",
    }


def to_json(card: dict) -> str:
    return json.dumps(card, indent=2, sort_keys=True, ensure_ascii=False)
