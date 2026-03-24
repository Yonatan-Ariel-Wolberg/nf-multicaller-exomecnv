#!/usr/bin/env python3
"""
Tests for .github/workflows/codeql.yml structure.

Validates that:
1. The workflow has the expected jobs: analyze and analyze-scripts.
2. The analyze job uses the correct matrix languages (python).
3. The analyze-scripts job has Install Nextflow and Analyze Nextflow scripts steps.
4. The Analyze Nextflow scripts step does NOT use the deprecated -dsl2 flag
   (DSL2 is the default in Nextflow 22+; the flag was removed in modern versions).
"""

import os

import yaml
import pytest

WORKFLOW_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    ".github",
    "workflows",
    "codeql.yml",
)


@pytest.fixture(scope="module")
def workflow():
    with open(WORKFLOW_PATH) as f:
        return yaml.safe_load(f)


def test_workflow_has_analyze_job(workflow):
    assert "analyze" in workflow["jobs"]


def test_workflow_has_analyze_scripts_job(workflow):
    assert "analyze-scripts" in workflow["jobs"]


def test_analyze_job_matrix_includes_python(workflow):
    matrix_include = workflow["jobs"]["analyze"]["strategy"]["matrix"]["include"]
    languages = [entry["language"] for entry in matrix_include]
    assert "python" in languages


def test_analyze_scripts_job_has_install_nextflow_step(workflow):
    steps = workflow["jobs"]["analyze-scripts"]["steps"]
    step_names = [s.get("name", "") for s in steps]
    assert "Install Nextflow" in step_names


def test_analyze_scripts_job_has_analyze_nextflow_step(workflow):
    steps = workflow["jobs"]["analyze-scripts"]["steps"]
    step_names = [s.get("name", "") for s in steps]
    assert "Analyze Nextflow scripts" in step_names


def test_nextflow_analyze_step_does_not_use_dsl2_flag(workflow):
    """The -dsl2 flag was removed in Nextflow 22+; DSL2 is now the default."""
    steps = workflow["jobs"]["analyze-scripts"]["steps"]
    for step in steps:
        if step.get("name") == "Analyze Nextflow scripts":
            run_cmd = step.get("run", "")
            assert "-dsl2" not in run_cmd, (
                "The -dsl2 flag is not supported in Nextflow 22+. "
                "DSL2 is the default; remove the flag."
            )
            return
    pytest.fail("Analyze Nextflow scripts step not found")
