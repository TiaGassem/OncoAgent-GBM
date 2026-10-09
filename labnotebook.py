"""OncoAgent-GBM: personal Lab Notebook & Protocols (private, export to Word/PDF).

PRIVACY MODEL (important, and honest):
  * Your entries are kept ONLY in your own browser SESSION (Streamlit
    session_state). They are NOT written to the server disk, NOT put in any
    shared database, and NOT logged. When your session ends they are gone from
    the server — the strongest privacy is that there is no server-side copy to
    leak.
  * To keep your notes you DOWNLOAD them to your own device: Word (.docx),
    PDF (.pdf), Markdown (.md) or an encrypted backup (.json/.enc) you can
    re-import later. Your data lives where YOU keep it.
  * The optional passphrase backup encrypts the file with AES (Fernet, key
    derived from your passphrase via PBKDF2-SHA256). We never see or store the
    passphrase; lose it and the backup cannot be recovered.

NO FABRICATION: this module only stores and formats what the USER types. It
never invents protocol numbers, results, or citations.
"""

from __future__ import annotations

import base64
import io
import json
import os
from datetime import datetime

# ---- optional deps, all defensive so the app never crashes on import -------
try:
    from docx import Document
    from docx.shared import Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    HAS_DOCX = True
except Exception:
    HAS_DOCX = False

try:
    from fpdf import FPDF
    HAS_FPDF = True
except Exception:
    HAS_FPDF = False

try:
    from cryptography.fernet import Fernet
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    HAS_CRYPTO = True
except Exception:
    HAS_CRYPTO = False


CATEGORIES = ["Lab note", "Protocol", "Observation / remark", "Result log",
              "Idea / to-do"]

_DISCLAIMER = ("Research/educational use only — not a clinical or diagnostic "
               "record. Entries are user-authored; verify any protocol against "
               "your validated SOP and cite real, verified sources.")


# ============================================================
# data model helpers
# ============================================================
def new_entry(title: str = "", author: str = "", category: str = "Lab note",
              tags: str = "", body: str = "", source: str = "") -> dict:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return {
        "id": base64.urlsafe_b64encode(os.urandom(6)).decode().rstrip("="),
        "title": title.strip() or "Untitled entry",
        "author": author.strip(),
        "category": category if category in CATEGORIES else "Lab note",
        "tags": tags.strip(),
        "body": body,
        "source": source.strip(),
        "created": now,
        "updated": now,
    }


def _latin1_safe(text: str) -> str:
    """fpdf2 core fonts are latin-1 only; replace unsupported chars."""
    if text is None:
        return ""
    repl = {"\u2014": "-", "\u2013": "-", "\u2018": "'", "\u2019": "'",
            "\u201c": '"', "\u201d": '"', "\u2022": "-", "\u2026": "...",
            "\u00b1": "+/-", "\u00b7": "-", "\u2192": "->", "\u2264": "<=",
            "\u2265": ">=", "\u00b5": "u", "\u03bc": "u", "\u00b0": " deg"}
    for k, v in repl.items():
        text = text.replace(k, v)
    return text.encode("latin-1", "replace").decode("latin-1")


# ============================================================
# exporters  (each returns bytes; never raises)
# ============================================================
def entries_to_markdown(entries: list[dict], owner: str = "") -> bytes:
    lines = ["# OncoAgent-GBM — Lab Notebook", ""]
    if owner:
        lines.append(f"**Owner:** {owner}  ")
    lines.append(f"**Exported:** {datetime.now().strftime('%Y-%m-%d %H:%M')}  ")
    lines.append(f"**Entries:** {len(entries)}")
    lines.append("")
    lines.append(f"> {_DISCLAIMER}")
    lines.append("")
    for i, e in enumerate(entries, 1):
        lines.append(f"## {i}. {e.get('title','Untitled')}")
        meta = []
        if e.get("category"):
            meta.append(f"**Category:** {e['category']}")
        if e.get("author"):
            meta.append(f"**Author:** {e['author']}")
        if e.get("created"):
            meta.append(f"**Created:** {e['created']}")
        if e.get("updated") and e.get("updated") != e.get("created"):
            meta.append(f"**Updated:** {e['updated']}")
        if e.get("tags"):
            meta.append(f"**Tags:** {e['tags']}")
        if e.get("source"):
            meta.append(f"**Source/ref:** {e['source']}")
        if meta:
            lines.append("  \n".join(meta))
            lines.append("")
        lines.append(e.get("body", "") or "_(empty)_")
        lines.append("")
    return "\n".join(lines).encode("utf-8")


def entries_to_docx(entries: list[dict], owner: str = "") -> bytes:
    if not HAS_DOCX:
        return b""
    try:
        doc = Document()
        style = doc.styles["Normal"].font
        style.name = "Calibri"
        style.size = Pt(11)

        h = doc.add_heading("OncoAgent-GBM — Lab Notebook", level=0)
        h.alignment = WD_ALIGN_PARAGRAPH.CENTER
        sub = doc.add_paragraph()
        sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = sub.add_run(f"Exported {datetime.now().strftime('%Y-%m-%d %H:%M')}"
                        + (f" · Owner: {owner}" if owner else "")
                        + f" · {len(entries)} entr" + ("y" if len(entries) == 1 else "ies"))
        r.italic = True
        r.font.size = Pt(9)

        d = doc.add_paragraph()
        dr = d.add_run(_DISCLAIMER)
        dr.italic = True
        dr.font.size = Pt(8)
        dr.font.color.rgb = RGBColor(0x99, 0x1b, 0x1b)

        for i, e in enumerate(entries, 1):
            doc.add_heading(f"{i}. {e.get('title','Untitled')}", level=1)
            meta_bits = []
            if e.get("category"):
                meta_bits.append(f"Category: {e['category']}")
            if e.get("author"):
                meta_bits.append(f"Author: {e['author']}")
            if e.get("created"):
                meta_bits.append(f"Created: {e['created']}")
            if e.get("updated") and e.get("updated") != e.get("created"):
                meta_bits.append(f"Updated: {e['updated']}")
            if e.get("tags"):
                meta_bits.append(f"Tags: {e['tags']}")
            if e.get("source"):
                meta_bits.append(f"Source/ref: {e['source']}")
            if meta_bits:
                mp = doc.add_paragraph()
                mr = mp.add_run("  |  ".join(meta_bits))
                mr.italic = True
                mr.font.size = Pt(8)
                mr.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
            for para in (e.get("body", "") or "").split("\n"):
                doc.add_paragraph(para)
        buf = io.BytesIO()
        doc.save(buf)
        return buf.getvalue()
    except Exception:
        return b""


def entries_to_pdf(entries: list[dict], owner: str = "") -> bytes:
    if not HAS_FPDF:
        return b""
    try:
        pdf = FPDF()
        pdf.set_auto_page_break(auto=True, margin=18)
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 15)
        pdf.cell(0, 10, _latin1_safe("OncoAgent-GBM - Lab Notebook"),
                 align="C", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "I", 8)
        meta = f"Exported {datetime.now().strftime('%Y-%m-%d %H:%M')}"
        if owner:
            meta += f"  -  Owner: {owner}"
        meta += f"  -  {len(entries)} entries"
        pdf.cell(0, 6, _latin1_safe(meta), align="C", new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0x99, 0x1b, 0x1b)
        pdf.multi_cell(0, 4, _latin1_safe(_DISCLAIMER))
        pdf.set_text_color(0, 0, 0)
        pdf.ln(2)
        for i, e in enumerate(entries, 1):
            pdf.set_font("Helvetica", "B", 12)
            pdf.multi_cell(0, 7, _latin1_safe(f"{i}. {e.get('title','Untitled')}"))
            bits = []
            for lbl, key in (("Category", "category"), ("Author", "author"),
                             ("Created", "created"), ("Tags", "tags"),
                             ("Source/ref", "source")):
                if e.get(key):
                    bits.append(f"{lbl}: {e[key]}")
            if bits:
                pdf.set_font("Helvetica", "I", 8)
                pdf.set_text_color(0x55, 0x55, 0x55)
                pdf.multi_cell(0, 4, _latin1_safe("  |  ".join(bits)))
                pdf.set_text_color(0, 0, 0)
            pdf.set_font("Helvetica", "", 10)
            pdf.multi_cell(0, 5, _latin1_safe(e.get("body", "") or "(empty)"))
            pdf.ln(3)
        return bytes(pdf.output())
    except Exception:
        return b""


# ============================================================
# backup / restore  (JSON, optionally passphrase-encrypted)
# ============================================================
def entries_to_json(entries: list[dict], owner: str = "") -> bytes:
    payload = {"app": "OncoAgent-GBM", "kind": "labnotebook", "version": 1,
               "owner": owner,
               "exported": datetime.now().isoformat(timespec="seconds"),
               "entries": entries}
    return json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")


def _derive_key(passphrase: str, salt: bytes) -> bytes:
    kdf = PBKDF2HMAC(algorithm=hashes.SHA256(), length=32, salt=salt,
                     iterations=390000)
    return base64.urlsafe_b64encode(kdf.derive(passphrase.encode("utf-8")))


def encrypt_backup(entries: list[dict], passphrase: str, owner: str = "") -> bytes:
    """Return an encrypted .enc blob: salt(16) + Fernet token. Never raises."""
    if not HAS_CRYPTO or not passphrase:
        return b""
    try:
        salt = os.urandom(16)
        key = _derive_key(passphrase, salt)
        token = Fernet(key).encrypt(entries_to_json(entries, owner))
        return b"OAGBM1" + salt + token
    except Exception:
        return b""


def decrypt_backup(blob: bytes, passphrase: str) -> list[dict]:
    """Decrypt an .enc blob back to a list of entries. Raises on wrong key."""
    if not HAS_CRYPTO:
        raise RuntimeError("Encryption library not available in this build.")
    if not blob.startswith(b"OAGBM1"):
        raise ValueError("Not an OncoAgent-GBM encrypted backup.")
    salt = blob[6:22]
    token = blob[22:]
    key = _derive_key(passphrase, salt)
    data = Fernet(key).decrypt(token)
    payload = json.loads(data.decode("utf-8"))
    return payload.get("entries", [])


def load_json_backup(raw: bytes) -> list[dict]:
    """Load a plain (unencrypted) JSON backup. Raises on bad input."""
    payload = json.loads(raw.decode("utf-8"))
    if isinstance(payload, dict) and "entries" in payload:
        return payload["entries"]
    if isinstance(payload, list):
        return payload
    raise ValueError("Unrecognised backup format.")
