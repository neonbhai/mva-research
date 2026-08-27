"""Rendering layer: the submission, the dossier, and the Track 2 reports.

Layer 7 (GP-01). It reads models, the evidence ledger and the clock; nothing
imports it except the composition root. It performs no analysis of its own — a
report that could compute a score could disagree with the pipeline that produced
it — and it holds one enforcement responsibility: GP-10, no claim without
evidence, implemented in :mod:`mva.reporting.assertions`.
"""

from __future__ import annotations

from mva.reporting.assertions import (
    CONTRADICTED_MARKER,
    TIER_MARKERS,
    Assertion,
    AssertionChecker,
    weakest_tier,
)
from mva.reporting.dossier import DOSSIER_QUESTIONS, build_candidate_dossier
from mva.reporting.render import (
    DEFAULT_TEMPLATES_DIR,
    format_cell,
    markdown_table,
    render_template,
)
from mva.reporting.track1 import (
    ACCEPTED_PROBAND_ID,
    EPCR_CEILING,
    EPCR_FLOOR,
    MAX_SUBMISSION_ROWS,
    TRACK1_COLUMNS,
    SubmissionRow,
    build_submission_rows,
    composite_to_epcr,
    render_submission_csv,
    truncation_notice,
    validate_submission,
    write_submission,
)
from mva.reporting.track2 import (
    DRUG_QUESTIONS,
    NOT_MEDICAL_ADVICE,
    UNCALIBRATED_WEIGHTS,
    build_drug_report,
    build_mechanism_report,
    build_rejection_record,
    build_track2_report,
)

__all__ = [
    "ACCEPTED_PROBAND_ID",
    "CONTRADICTED_MARKER",
    "DEFAULT_TEMPLATES_DIR",
    "DOSSIER_QUESTIONS",
    "DRUG_QUESTIONS",
    "EPCR_CEILING",
    "EPCR_FLOOR",
    "MAX_SUBMISSION_ROWS",
    "NOT_MEDICAL_ADVICE",
    "TIER_MARKERS",
    "TRACK1_COLUMNS",
    "UNCALIBRATED_WEIGHTS",
    "Assertion",
    "AssertionChecker",
    "SubmissionRow",
    "build_candidate_dossier",
    "build_drug_report",
    "build_mechanism_report",
    "build_rejection_record",
    "build_submission_rows",
    "build_track2_report",
    "composite_to_epcr",
    "format_cell",
    "markdown_table",
    "render_submission_csv",
    "render_template",
    "truncation_notice",
    "validate_submission",
    "weakest_tier",
    "write_submission",
]
