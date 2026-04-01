#!/usr/bin/env python3
"""Guard against regex syntax that breaks Nextflow/Groovy parsing."""

from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parent.parent


def test_nf_files_do_not_use_trailing_regex_i_flag():
    bad_pattern = re.compile(r"/[^/\n]*/i\b")
    for nf_path in REPO_ROOT.rglob("*.nf"):
        text = nf_path.read_text(encoding="utf-8")
        assert bad_pattern.search(text) is None, (
            f"{nf_path} contains /.../i regex syntax; use inline (?i) instead"
        )
