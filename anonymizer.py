"""OncoAgent-GBM: Local PII anonymizer for medical notes using Presidio."""

from __future__ import annotations

import re
import hashlib
from dataclasses import dataclass, field
from typing import Optional

try:
    from presidio_analyzer import AnalyzerEngine, RecognizerRegistry
    from presidio_anonymizer import AnonymizerEngine
    from presidio_anonymizer.entities import OperatorConfig, OperatorResult
    HAS_PRESIDIO = True
except ImportError:
    HAS_PRESIDIO = False


@dataclass
class AnonymizerResult:
    """Result of an anonymization operation."""
    original_text: str = ""
    anonymized_text: str = ""
    entities_found: list[dict] = field(default_factory=list)
    entity_count: int = 0
    anonymization_method: str = "replace"
    error: str = ""


@dataclass
class PIIDetection:
    """A single detected PII entity."""
    entity_type: str = ""
    start: int = 0
    end: int = 0
    score: float = 0.0
    original_text: str = ""
    anonymized_text: str = ""


class MedicalNoteAnonymizer:
    """Local PII anonymizer for medical notes using Microsoft Presidio.

    All processing occurs in-memory with no external API calls.
    Supports: names, dates, phone numbers, emails, addresses, medical record numbers,
    SSNs, and other healthcare PII.
    """

    DEFAULT_ENTITIES = [
        "PERSON",
        "EMAIL_ADDRESS",
        "PHONE_NUMBER",
        "CREDIT_CARD",
        "IP_ADDRESS",
        "LOCATION",
        "MEDICAL_LICENSE",
        "US_SSN",
        "US_PASSPORT",
        "US_DRIVER_LICENSE",
        "DATE_TIME",
        "NRP",
        "NUMERIC_IDENTIFIER",
    ]

    MEDICAL_CUSTOM_ENTITIES = [
        "MRN",  # Medical Record Number
        "HEALTH_PLAN_ID",
        "ACCOUNT_NUMBER",
        "CERTIFICATE_LICENSE_NUMBER",
    ]

    def __init__(self, language: str = "en"):
        self.language = language
        self._analyzer = None
        self._anonymizer = None
        self._initialize_presidio()

    def _initialize_presidio(self):
        """Initialize Presidio analyzer and anonymizer engines."""
        if not HAS_PRESIDIO:
            return

        try:
            registry = RecognizerRegistry()
            registry.load_predefined_recognizers()
            self._analyzer = AnalyzerEngine(registry=registry)
            self._anonymizer = AnonymizerEngine()
        except Exception:
            try:
                self._analyzer = AnalyzerEngine()
                self._anonymizer = AnonymizerEngine()
            except Exception:
                self._analyzer = None
                self._anonymizer = None

    def detect_pii(
        self,
        text: str,
        entities: Optional[list[str]] = None,
        score_threshold: float = 0.3,
    ) -> list[PIIDetection]:
        """Detect PII entities in text without anonymizing."""
        if not HAS_PRESIDIO or self._analyzer is None:
            return self._regex_fallback_detect(text)

        if entities is None:
            entities = self.DEFAULT_ENTITIES

        try:
            results = self._analyzer.analyze(
                text=text,
                entities=entities,
                language=self.language,
                score_threshold=score_threshold,
            )

            detections = []
            for result in results:
                detection = PIIDetection(
                    entity_type=result.entity_type,
                    start=result.start,
                    end=result.end,
                    score=round(result.score, 3),
                    original_text=text[result.start:result.end],
                )
                detections.append(detection)

            return detections

        except Exception:
            return self._regex_fallback_detect(text)

    def anonymize(
        self,
        text: str,
        method: str = "replace",
        entities: Optional[list[str]] = None,
        score_threshold: float = 0.3,
        custom_label: Optional[str] = None,
    ) -> AnonymizerResult:
        """Anonymize PII in medical text.

        Args:
            text: Input text to anonymize.
            method: Anonymization method ('replace', 'mask', 'hash', 'redact', 'synthesize').
            entities: List of entity types to anonymize.
            score_threshold: Minimum confidence score for detection.
            custom_label: Custom replacement label (e.g., '<PATIENT>').
        """
        result = AnonymizerResult(original_text=text, anonymization_method=method)

        if not text or not text.strip():
            return result

        detections = self.detect_pii(text, entities, score_threshold)
        result.entities_found = [
            {
                "type": d.entity_type,
                "text": d.original_text,
                "score": d.score,
                "position": f"{d.start}:{d.end}",
            }
            for d in detections
        ]
        result.entity_count = len(detections)

        if not HAS_PRESIDIO or self._analyzer is None or self._anonymizer is None:
            result.anonymized_text = self._regex_fallback_anonymize(text, method)
            return result

        try:
            analyzer_results = self._analyzer.analyze(
                text=text,
                entities=entities or self.DEFAULT_ENTITIES,
                language=self.language,
                score_threshold=score_threshold,
            )

            operators = self._build_operators(method, custom_label)

            anonymized = self._anonymizer.anonymize(
                text=text,
                analyzer_results=analyzer_results,
                operators=operators,
            )

            result.anonymized_text = anonymized.text

        except Exception:
            result.anonymized_text = self._regex_fallback_anonymize(text, method)

        return result

    def _build_operators(
        self, method: str, custom_label: Optional[str]
    ) -> dict[str, OperatorConfig]:
        """Build Presidio operator configuration based on anonymization method."""
        if method == "replace":
            label = custom_label or "<PII>"
            return {
                "PERSON": OperatorConfig("replace", {"new_value": label or "<PATIENT>"}),
                "EMAIL_ADDRESS": OperatorConfig("replace", {"new_value": "<EMAIL>"}),
                "PHONE_NUMBER": OperatorConfig("replace", {"new_value": "<PHONE>"}),
                "LOCATION": OperatorConfig("replace", {"new_value": "<LOCATION>"}),
                "US_SSN": OperatorConfig("replace", {"new_value": "<SSN>"}),
                "DATE_TIME": OperatorConfig("replace", {"new_value": "<DATE>"}),
                "DEFAULT": OperatorConfig("replace", {"new_value": label or "<PII>"}),
            }
        elif method == "mask":
            return {
                "DEFAULT": OperatorConfig("mask", {
                    "masking_char": "*",
                    "chars_to_mask": 20,
                    "from_end": False,
                }),
            }
        elif method == "hash":
            return {"DEFAULT": OperatorConfig("hash", {"hash_type": "sha256"})}
        elif method == "redact":
            return {"DEFAULT": OperatorConfig("redact", {})}
        elif method == "synthesize":
            return {"DEFAULT": OperatorConfig("synthesize", {})}
        else:
            label = custom_label or "<PII>"
            return {"DEFAULT": OperatorConfig("replace", {"new_value": label})}

    def _regex_fallback_detect(self, text: str) -> list[PIIDetection]:
        """Fallback PII detection using regex patterns when Presidio is unavailable."""
        patterns = {
            "EMAIL_ADDRESS": r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b',
            "PHONE_NUMBER": r'\b(?:\+?1[-.\s]?)?(?:\(?\d{3}\)?[-.\s]?)?\d{3}[-.\s]?\d{4}\b',
            "US_SSN": r'\b\d{3}[-.\s]?\d{2}[-.\s]?\d{4}\b',
            "DATE_TIME": r'\b\d{1,2}[/-]\d{1,2}[/-]\d{2,4}\b|\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\w*\s+\d{1,2},?\s+\d{4}\b',
            "MEDICAL_RECORD": r'\b(?:MRN|MR#|Medical Record)[:\s#]*(\d{6,12})\b',
            "PERSON": r'\b(?:Dr|Mr|Mrs|Ms|Prof)\.?\s+[A-Z][a-z]+(?:\s+[A-Z][a-z]+)?\b',
        }

        detections = []
        for entity_type, pattern in patterns.items():
            for match in re.finditer(pattern, text, re.IGNORECASE):
                detections.append(PIIDetection(
                    entity_type=entity_type,
                    start=match.start(),
                    end=match.end(),
                    score=0.8,
                    original_text=match.group(),
                ))

        return sorted(detections, key=lambda d: d.start)

    def _regex_fallback_anonymize(self, text: str, method: str) -> str:
        """Fallback anonymization using regex when Presidio is unavailable."""
        detections = self._regex_fallback_detect(text)

        replacements = {}
        for det in detections:
            if method == "replace":
                type_map = {
                    "PERSON": "<PATIENT>",
                    "EMAIL_ADDRESS": "<EMAIL>",
                    "PHONE_NUMBER": "<PHONE>",
                    "US_SSN": "<SSN>",
                    "DATE_TIME": "<DATE>",
                    "MEDICAL_RECORD": "<MRN>",
                    "LOCATION": "<LOCATION>",
                }
                replacement = type_map.get(det.entity_type, "<PII>")
            elif method == "hash":
                replacement = hashlib.sha256(det.original_text.encode()).hexdigest()[:12]
            elif method == "mask":
                replacement = "*" * len(det.original_text)
            elif method == "redact":
                replacement = "[REDACTED]"
            else:
                replacement = "<PII>"

            replacements[(det.start, det.end)] = replacement

        anonymized = []
        prev_end = 0
        for (start, end), replacement in sorted(replacements.items()):
            anonymized.append(text[prev_end:start])
            anonymized.append(replacement)
            prev_end = end
        anonymized.append(text[prev_end:])

        return "".join(anonymized)


def anonymize_medical_note(
    text: str,
    method: str = "replace",
    entities: Optional[list[str]] = None,
) -> AnonymizerResult:
    """Convenience function to anonymize a medical note."""
    anonymizer = MedicalNoteAnonymizer()
    return anonymizer.anonymize(text, method=method, entities=entities)


def detect_pii_in_text(
    text: str,
    score_threshold: float = 0.3,
) -> list[PIIDetection]:
    """Convenience function to detect PII without anonymizing."""
    anonymizer = MedicalNoteAnonymizer()
    return anonymizer.detect_pii(text, score_threshold=score_threshold)


# GBM-specific anonymization patterns
GBM_MEDICAL_NOTE_TEMPLATE = """Patient: {patient_name}
DOB: {dob}
MRN: {mrn}
Diagnosis: Glioblastoma Multiforme (ICD-10: C71.x)
Date of Diagnosis: {diag_date}

Tumor Profile:
- Location: {tumor_location}
- Grade: WHO Grade IV
- Molecular markers: {molecular_markers}

Treatment History:
- Surgery: {surgery_date} ({surgery_type})
- Radiation: {rt_dates}
- Chemotherapy: TMZ {tmz_schedule}

Genetic Testing:
- MGMT status: {mgmt_status}
- IDH1/2: {idh_status}
- EGFRvIII: {egfr_status}
- 1p/19q codeletion: {codeletion_status}

Follow-up Notes:
{followup_notes}

Physician: {physician_name}
Institution: {institution}
"""


def generate_sample_medical_note() -> str:
    """Generate a sample GBM medical note for testing the anonymizer."""
    return GBM_MEDICAL_NOTE_TEMPLATE.format(
        patient_name="John Michael Smith",
        dob="03/15/1965",
        mrn="MRN-00482937",
        diag_date="01/12/2025",
        tumor_location="Right temporal lobe, extending to insular cortex",
        molecular_markers="EGFRvIII+, MGMT unmethylated, IDH1 wildtype",
        surgery_date="01/20/2025",
        surgery_type="Gross total resection (95%)",
        rt_dates="02/10/2025 - 03/14/2025 (60 Gy/30 fractions)",
        tmz_schedule="Concurrent: 75 mg/m2 daily x 42 days; Adjuvant: 150-200 mg/m2 x 5/28 days",
        mgmt_status="Unmethylated",
        idh_status="Wildtype (R132H negative)",
        egfr_status="Positive (EGFRvIII amplification detected by FISH)",
        codeletion_status="Intact (no codeletion)",
        followup_notes=(
            "Patient tolerated concurrent chemoradiation well, with grade 2 fatigue and mild nausea. "
            "MRI at 4 weeks post-RT showed expected treatment changes, no definite progression. "
            "Adjuvant TMZ initiated. Follow-up in 4 weeks. Contact patient at john.smith@email.com "
            "or (555) 234-5678. SSN: 432-10-9876."
        ),
        physician_name="Dr. Sarah Johnson",
        institution="University Medical Center Neuro-Oncology",
    )
