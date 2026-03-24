#!/usr/bin/env python3
"""Ensure README lists currently unhandled use-case validation gaps."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
README = REPO_ROOT / "README.md"


def test_readme_lists_known_unhandled_use_case_errors_section():
    text = README.read_text(encoding="utf-8")
    assert "### Known unhandled use-case errors (current gaps)" in text

    expected_items = [
        "--workflow canoes` / `xhmm` / `clamms",
        "--workflow gcnv",
        "--workflow cnvkit",
        "--workflow dragen",
        "--workflow normalise` / `evaluate",
        "--workflow feature_extraction",
        "--workflow train",
        "--workflow full",
        "Any workflow: input files exist but are not readable",
        "Any workflow: output path exists but is not writable",
    ]
    for item in expected_items:
        assert item in text, f"README missing unhandled use-case error bullet: {item}"
