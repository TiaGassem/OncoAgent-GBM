"""OncoAgent-GBM: OPTIONAL bring-your-own-key LLM polish layer.

The chat assistant is GROUNDED by default: every answer is built from retrieved,
cited sources (live PubMed + curated validated topics) and never fabricates.

This module adds an OPTIONAL layer on top: if — and only if — the user supplies
THEIR OWN API key (OpenAI or Anthropic), it rewrites the already-retrieved,
already-cited draft into fluent prose. The key is NEVER hardcoded here; it is
passed in at call time (from a Streamlit field or st.secrets) and only used for
that one request.

HARD GROUNDING RULE sent to the model: it may ONLY rephrase/clarify the material
it is given. It must not add facts, numbers, drug doses, mechanisms or citations
that are not already present, and it must preserve every PMID/DOI. If the
material does not answer the question, it must say so. This keeps the fluent
mode consistent with the no-fabrication policy.

If no key is given, or the call fails for any reason, we return the grounded
draft unchanged — the app never breaks and never silently invents content.
"""

from __future__ import annotations

import requests


_SYSTEM = (
    "You are a glioblastoma (GBM) / oncology RESEARCH assistant. You will be "
    "given (a) a user question and (b) DRAFT material that was retrieved from "
    "cited scientific sources (PubMed abstracts, PMIDs/DOIs, validated app "
    "topics). Your ONLY job is to rewrite the draft into a clear, well-"
    "structured answer for a researcher.\n\n"
    "ABSOLUTE RULES:\n"
    "1. Do NOT add any fact, number, dose, concentration, timing, mechanism, "
    "gene, statistic or claim that is not already in the draft material.\n"
    "2. Do NOT invent, guess or add any citation, PMID, DOI or reference. "
    "Preserve every PMID/DOI exactly as given.\n"
    "3. If the draft material does not actually answer the question, say so "
    "plainly instead of filling the gap.\n"
    "4. This is research/educational support only — never give medical advice, "
    "diagnosis, or treatment recommendations for real patients.\n"
    "5. Keep it concise. Do not pad."
)


def polish_answer(question: str, grounded_markdown: str, provider: str,
                  api_key: str, model: str = "", timeout: int = 40) -> dict:
    """Optionally rephrase a grounded answer with the user's own LLM key.

    Returns {text, used_llm, error}. On ANY problem, text falls back to the
    original grounded_markdown and used_llm is False.
    """
    out = {"text": grounded_markdown, "used_llm": False, "error": ""}
    key = (api_key or "").strip()
    if not key:
        return out  # grounded-only: nothing to do

    user_block = (
        f"USER QUESTION:\n{question}\n\n"
        f"DRAFT MATERIAL (retrieved + cited — rephrase ONLY this, add nothing):\n"
        f"{grounded_markdown}"
    )

    try:
        provider = (provider or "").strip().lower()
        if provider in ("openai", "gpt"):
            mdl = model or "gpt-4o-mini"
            resp = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={"Authorization": f"Bearer {key}",
                         "Content-Type": "application/json"},
                json={
                    "model": mdl,
                    "temperature": 0.2,
                    "messages": [
                        {"role": "system", "content": _SYSTEM},
                        {"role": "user", "content": user_block},
                    ],
                },
                timeout=timeout,
            )
            if resp.status_code != 200:
                out["error"] = f"OpenAI API error {resp.status_code}: {resp.text[:200]}"
                return out
            data = resp.json()
            text = data["choices"][0]["message"]["content"].strip()
            if text:
                out["text"] = text
                out["used_llm"] = True
            return out

        if provider in ("anthropic", "claude"):
            mdl = model or "claude-3-5-haiku-latest"
            resp = requests.post(
                "https://api.anthropic.com/v1/messages",
                headers={"x-api-key": key,
                         "anthropic-version": "2023-06-01",
                         "Content-Type": "application/json"},
                json={
                    "model": mdl,
                    "max_tokens": 1024,
                    "temperature": 0.2,
                    "system": _SYSTEM,
                    "messages": [{"role": "user", "content": user_block}],
                },
                timeout=timeout,
            )
            if resp.status_code != 200:
                out["error"] = f"Anthropic API error {resp.status_code}: {resp.text[:200]}"
                return out
            data = resp.json()
            parts = data.get("content", [])
            text = "".join(p.get("text", "") for p in parts).strip()
            if text:
                out["text"] = text
                out["used_llm"] = True
            return out

        out["error"] = f"Unknown provider '{provider}'. Use 'openai' or 'anthropic'."
        return out
    except Exception as e:  # network, JSON, key errors — fall back silently
        out["error"] = f"{type(e).__name__}: {e}"
        return out
