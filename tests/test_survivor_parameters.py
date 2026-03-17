#!/usr/bin/env python3
"""
Tests for modules/modules-survivor.nf SURVIVOR merge parameters.

Verifies that the SURVIVOR merge command uses the correct parameters as
recommended by the official SURVIVOR documentation:
  https://github.com/fritzsedlazeck/SURVIVOR/wiki/Methods-and-Parameter

Expected call:
  SURVIVOR merge <file_list> 1000 2 1 1 0 30 <output>

Parameters:
  max_dist    = 1000  Maximum distance between breakpoints (bp)
  min_callers = 2     Minimum number of callers that must support an SV
  use_type    = 1     Require same SV type (1=yes)
  use_strand  = 1     Require same strand (1=yes)
  est_dist    = 0     Do not estimate SV distance
  min_sv_size = 30    Minimum SV size (bp)
"""

import os
import re

MODULE_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'modules', 'modules-survivor.nf'
)


def _get_survivor_merge_args():
    """Parse the SURVIVOR merge command from the Nextflow module and return
    positional arguments as a list of strings (excluding the list file and
    output filename placeholders)."""
    with open(MODULE_PATH) as fh:
        content = fh.read()

    # Find the SURVIVOR merge invocation block and collect the numeric args.
    # In the .nf file, Nextflow script blocks escape backslashes as \\, so
    # line continuations appear as \\\\ in the raw file content.
    match = re.search(
        r'SURVIVOR\s+merge\s+.*?(\d+)\s*\\\\[^\n]*\n\s*(\d+)\s*\\\\[^\n]*\n\s*(\d+)\s*\\\\[^\n]*\n\s*(\d+)\s*\\\\[^\n]*\n\s*(\d+)\s*\\\\[^\n]*\n\s*(\d+)',
        content,
        re.DOTALL,
    )
    assert match, "Could not locate SURVIVOR merge numeric arguments in module file"
    return [int(g) for g in match.groups()]


def test_survivor_merge_max_dist():
    """max_dist should be 1000 bp (pairwise breakpoint distance threshold)."""
    args = _get_survivor_merge_args()
    max_dist = args[0]
    assert max_dist == 1000, f"max_dist expected 1000, got {max_dist}"


def test_survivor_merge_min_callers():
    """min_callers should be 2 (consensus: SV must be supported by >=2 callers)."""
    args = _get_survivor_merge_args()
    min_callers = args[1]
    assert min_callers == 2, (
        f"min_callers expected 2 (consensus calling per SURVIVOR docs), got {min_callers}"
    )


def test_survivor_merge_use_type():
    """use_type should be 1 (require same SV type across callers)."""
    args = _get_survivor_merge_args()
    use_type = args[2]
    assert use_type == 1, f"use_type expected 1, got {use_type}"


def test_survivor_merge_use_strand():
    """use_strand should be 1 (require same strand per SURVIVOR reference docs)."""
    args = _get_survivor_merge_args()
    use_strand = args[3]
    assert use_strand == 1, (
        f"use_strand expected 1 (per SURVIVOR reference documentation), got {use_strand}"
    )


def test_survivor_merge_est_dist():
    """est_dist should be 0 (do not estimate SV distance)."""
    args = _get_survivor_merge_args()
    est_dist = args[4]
    assert est_dist == 0, f"est_dist expected 0, got {est_dist}"


def test_survivor_merge_min_sv_size():
    """min_sv_size should be 30 bp."""
    args = _get_survivor_merge_args()
    min_sv_size = args[5]
    assert min_sv_size == 30, f"min_sv_size expected 30, got {min_sv_size}"
