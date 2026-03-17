#!/usr/bin/env python3
"""
Tests for the normalise_cnv_caller_quality_scores.py script and its integration
into each CNV caller Nextflow module.

Validates:
  1. The script has a shebang, imports argparse and pysam, and exposes a CLI
     with --input_vcf, --output_vcf, and --caller arguments.
  2. Each caller is handled with the correct quality-score fields (FORMAT vs INFO).
  3. Each Nextflow module contains the NORMALISE_CNV_QUALITY_SCORES process with
     the correct --caller value.
  4. Each module's workflow calls the process and emits normalised_vcf /
     normalised_vcf_index.
  5. nextflow.config has a 'pysam' label entry for the container.
"""

import os
import re

import pytest

BIN_DIR  = os.path.join(os.path.dirname(__file__), '..', 'bin')
REPO_ROOT = os.path.join(os.path.dirname(__file__), '..')

SCRIPT_PATH = os.path.join(BIN_DIR, 'normalise_cnv_caller_quality_scores.py')

MODULES = {
    'CANOES':   'modules/modules-canoes.nf',
    'CLAMMS':   'modules/modules-clamms.nf',
    'XHMM':     'modules/modules-xhmm.nf',
    'INDELIBLE':'modules/modules-indelible.nf',
    'GATK':     'modules/modules-gatk-gcnv.nf',
    'CNVKIT':   'modules/modules-cnvkit.nf',
    'DRAGEN':   'modules/modules-icav2-dragen.nf',
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_file(rel_path):
    path = os.path.join(REPO_ROOT, rel_path)
    with open(path) as fh:
        return fh.read()


def _read_script():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


def _extract_process(module_text, process_name):
    match = re.search(
        rf"process {re.escape(process_name)}\s*\{{(.+?)(?=\nprocess |\nworkflow |\Z)",
        module_text,
        re.DOTALL,
    )
    return match.group(1) if match else None


def _extract_workflow(module_text, workflow_name):
    match = re.search(
        rf"workflow {re.escape(workflow_name)}\s*\{{(.+)",
        module_text,
        re.DOTALL,
    )
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def script_text():
    return _read_script()


@pytest.fixture(scope='module')
def canoes_text():
    return _read_file('modules/modules-canoes.nf')


@pytest.fixture(scope='module')
def clamms_text():
    return _read_file('modules/modules-clamms.nf')


@pytest.fixture(scope='module')
def xhmm_text():
    return _read_file('modules/modules-xhmm.nf')


@pytest.fixture(scope='module')
def indelible_text():
    return _read_file('modules/modules-indelible.nf')


@pytest.fixture(scope='module')
def gatk_text():
    return _read_file('modules/modules-gatk-gcnv.nf')


@pytest.fixture(scope='module')
def cnvkit_text():
    return _read_file('modules/modules-cnvkit.nf')


@pytest.fixture(scope='module')
def dragen_text():
    return _read_file('modules/modules-icav2-dragen.nf')


@pytest.fixture(scope='module')
def config_text():
    return _read_file('nextflow.config')


# ===========================================================================
# 1. Script structure and CLI
# ===========================================================================

class TestScriptStructure:
    """normalise_cnv_caller_quality_scores.py must be a well-formed CLI script."""

    def test_shebang_present(self, script_text):
        assert script_text.startswith('#!/usr/bin/env python3'), (
            "Script must start with #!/usr/bin/env python3 shebang"
        )

    def test_imports_argparse(self, script_text):
        assert 'import argparse' in script_text, (
            "Script must import argparse to provide CLI argument parsing"
        )

    def test_imports_pysam(self, script_text):
        assert 'import pysam' in script_text, (
            "Script must import pysam to read and write VCF files"
        )

    def test_has_main_block(self, script_text):
        assert "__name__ == '__main__'" in script_text or \
               '__name__ == "__main__"' in script_text, (
            "Script must have an if __name__ == '__main__' guard so it is "
            "usable both as a module and as a standalone CLI tool"
        )

    def test_cli_has_input_vcf(self, script_text):
        assert '--input_vcf' in script_text, (
            "CLI must accept --input_vcf argument"
        )

    def test_cli_has_output_vcf(self, script_text):
        assert '--output_vcf' in script_text, (
            "CLI must accept --output_vcf argument"
        )

    def test_cli_has_caller(self, script_text):
        assert '--caller' in script_text, (
            "CLI must accept --caller argument to identify the CNV caller"
        )

    def test_cli_choices_include_all_callers(self, script_text):
        for caller in ('CANOES', 'CLAMMS', 'XHMM', 'GATK', 'DRAGEN', 'CNVKIT', 'INDELIBLE'):
            assert caller in script_text, (
                f"CLI --caller choices must include {caller}"
            )


# ===========================================================================
# 2. Caller-specific quality-score field correctness
# ===========================================================================

class TestCallerFieldCorrectness:
    """Each caller block must read quality scores from the correct VCF location."""

    def test_canoes_reads_from_format_not_info(self, script_text):
        """CANOES stores Q_SOME in FORMAT, not INFO."""
        canoes_block = re.search(
            r'if caller == "CANOES":(.+?)(?=elif caller|$)',
            script_text, re.DOTALL
        )
        assert canoes_block, "CANOES caller block not found"
        block = canoes_block.group(1)
        # Must use sample["Q_SOME"] (FORMAT), not record.info["Q_SOME"]
        assert 'record.info' not in block, (
            "CANOES block must not access Q_SOME via record.info; "
            "Q_SOME is a FORMAT field in CANOES VCFs"
        )
        assert 'sample' in block and '"Q_SOME"' in block, (
            "CANOES block must read Q_SOME from sample FORMAT fields"
        )

    def test_clamms_reads_from_format_not_info(self, script_text):
        """CLAMMS stores Q_SOME and Q_EXACT in FORMAT, not INFO."""
        clamms_block = re.search(
            r'elif caller == "CLAMMS":(.+?)(?=elif caller|$)',
            script_text, re.DOTALL
        )
        assert clamms_block, "CLAMMS caller block not found"
        block = clamms_block.group(1)
        # Must use sample["Q_SOME"] and sample["Q_EXACT"] (FORMAT), not record.info
        assert 'record.info' not in block, (
            "CLAMMS block must not access Q_SOME/Q_EXACT via record.info; "
            "they are FORMAT fields in CLAMMS VCFs"
        )
        assert 'sample' in block and '"Q_SOME"' in block, (
            "CLAMMS block must read Q_SOME from sample FORMAT fields"
        )
        assert '"Q_EXACT"' in block, (
            "CLAMMS block must read Q_EXACT from sample FORMAT fields "
            "to gate the normalisation: Q_EXACT < 0 → QUAL_norm = 0"
        )

    def test_indelible_uses_uppercase_info_fields(self, script_text):
        """INDELIBLE INFO fields are uppercase: SR_TOTAL, AVG_MAPQ."""
        indelible_block = re.search(
            r'elif caller == "INDELIBLE":(.+?)(?=except|$)',
            script_text, re.DOTALL
        )
        assert indelible_block, "INDELIBLE caller block not found"
        block = indelible_block.group(1)
        assert '"SR_TOTAL"' in block, (
            "INDELIBLE block must use 'SR_TOTAL' (uppercase) as the INFO field key"
        )
        assert '"AVG_MAPQ"' in block, (
            "INDELIBLE block must use 'AVG_MAPQ' (uppercase) as the INFO field key"
        )
        assert '"sr_total"' not in block, (
            "INDELIBLE block must not use lowercase 'sr_total' (field is uppercase)"
        )
        assert '"avg_mapq"' not in block, (
            "INDELIBLE block must not use lowercase 'avg_mapq' (field is uppercase)"
        )

    def test_indelible_applies_mum_dad_sr_gates(self, script_text):
        """INDELIBLE block must check mum_sr < 2 and dad_sr < 2 as de-novo quality gates."""
        indelible_block = re.search(
            r'elif caller == "INDELIBLE":(.+?)(?=except|$)',
            script_text, re.DOTALL
        )
        assert indelible_block, "INDELIBLE caller block not found"
        block = indelible_block.group(1)
        assert '"MUM_SR"' in block or 'mum_sr' in block, (
            "INDELIBLE block must check MUM_SR: if mum_sr >= 2.0 → QUAL_norm = 0.0"
        )
        assert '"DAD_SR"' in block or 'dad_sr' in block, (
            "INDELIBLE block must check DAD_SR: if dad_sr >= 2.0 → QUAL_norm = 0.0"
        )

    def test_sample_accessor_uses_first_value(self, script_text):
        """sample must be obtained from record.samples.values(), not record.samples."""
        assert 'record.samples.values()' in script_text, (
            "sample must be obtained via next(iter(record.samples.values())) "
            "to correctly access the single-sample FORMAT fields in pysam"
        )
        # The old broken assignment was: sample = record.samples
        assert 'sample = record.samples\n' not in script_text and \
               'sample = record.samples  ' not in script_text, (
            "Assigning sample = record.samples is incorrect; it returns the "
            "VariantRecordSamples dict-like object, not a single sample's FORMAT fields"
        )


# ===========================================================================
# 3. NORMALISE_CNV_QUALITY_SCORES process in each module
# ===========================================================================

class TestNormaliseProcessExists:
    """Each module must contain a NORMALISE_CNV_QUALITY_SCORES process."""

    def _check_process(self, module_text, caller_name):
        body = _extract_process(module_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert body is not None, (
            f"NORMALISE_CNV_QUALITY_SCORES process not found in {caller_name} module"
        )
        return body

    def test_canoes_has_process(self, canoes_text):
        self._check_process(canoes_text, 'CANOES')

    def test_clamms_has_process(self, clamms_text):
        self._check_process(clamms_text, 'CLAMMS')

    def test_xhmm_has_process(self, xhmm_text):
        self._check_process(xhmm_text, 'XHMM')

    def test_indelible_has_process(self, indelible_text):
        self._check_process(indelible_text, 'INDELIBLE')

    def test_gatk_has_process(self, gatk_text):
        self._check_process(gatk_text, 'GATK')

    def test_cnvkit_has_process(self, cnvkit_text):
        self._check_process(cnvkit_text, 'CNVKIT')

    def test_dragen_has_process(self, dragen_text):
        self._check_process(dragen_text, 'DRAGEN')


class TestNormaliseProcessCaller:
    """NORMALISE_CNV_QUALITY_SCORES process must pass the correct --caller value."""

    def _check_caller(self, module_text, expected_caller):
        body = _extract_process(module_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert body is not None, (
            f"NORMALISE_CNV_QUALITY_SCORES process not found in module"
        )
        assert f'--caller {expected_caller}' in body, (
            f"NORMALISE_CNV_QUALITY_SCORES must pass --caller {expected_caller}"
        )

    def test_canoes_caller(self, canoes_text):
        self._check_caller(canoes_text, 'CANOES')

    def test_clamms_caller(self, clamms_text):
        self._check_caller(clamms_text, 'CLAMMS')

    def test_xhmm_caller(self, xhmm_text):
        self._check_caller(xhmm_text, 'XHMM')

    def test_indelible_caller(self, indelible_text):
        self._check_caller(indelible_text, 'INDELIBLE')

    def test_gatk_caller(self, gatk_text):
        self._check_caller(gatk_text, 'GATK')

    def test_cnvkit_caller(self, cnvkit_text):
        self._check_caller(cnvkit_text, 'CNVKIT')

    def test_dragen_caller(self, dragen_text):
        self._check_caller(dragen_text, 'DRAGEN')


class TestNormaliseProcessOutputs:
    """NORMALISE_CNV_QUALITY_SCORES must emit normalised_vcf and normalised_vcf_index."""

    def _check_outputs(self, module_text, caller_name):
        body = _extract_process(module_text, 'NORMALISE_CNV_QUALITY_SCORES')
        assert body is not None, f"NORMALISE_CNV_QUALITY_SCORES not found for {caller_name}"
        assert 'normalised_vcf' in body, (
            f"{caller_name}: process must emit normalised_vcf"
        )
        assert 'normalised_vcf_index' in body, (
            f"{caller_name}: process must emit normalised_vcf_index"
        )
        assert '.normalised.vcf.gz' in body, (
            f"{caller_name}: output files must use .normalised.vcf.gz suffix"
        )

    def test_canoes_outputs(self, canoes_text):
        self._check_outputs(canoes_text, 'CANOES')

    def test_clamms_outputs(self, clamms_text):
        self._check_outputs(clamms_text, 'CLAMMS')

    def test_xhmm_outputs(self, xhmm_text):
        self._check_outputs(xhmm_text, 'XHMM')

    def test_indelible_outputs(self, indelible_text):
        self._check_outputs(indelible_text, 'INDELIBLE')

    def test_gatk_outputs(self, gatk_text):
        self._check_outputs(gatk_text, 'GATK')

    def test_cnvkit_outputs(self, cnvkit_text):
        self._check_outputs(cnvkit_text, 'CNVKIT')

    def test_dragen_outputs(self, dragen_text):
        self._check_outputs(dragen_text, 'DRAGEN')


# ===========================================================================
# 4. Workflow integration: process called and emits added
# ===========================================================================

class TestWorkflowCallsNormalise:
    """Each workflow must invoke NORMALISE_CNV_QUALITY_SCORES after BGZIP/annotate."""

    def _check_workflow_call(self, module_text, workflow_name, caller_name):
        wf = _extract_workflow(module_text, workflow_name)
        assert wf is not None, f"{workflow_name} workflow not found"
        assert 'NORMALISE_CNV_QUALITY_SCORES' in wf, (
            f"{caller_name} workflow must call NORMALISE_CNV_QUALITY_SCORES"
        )

    def test_canoes_workflow(self, canoes_text):
        self._check_workflow_call(canoes_text, 'CANOES', 'CANOES')

    def test_clamms_workflow(self, clamms_text):
        self._check_workflow_call(clamms_text, 'CLAMMS', 'CLAMMS')

    def test_xhmm_workflow(self, xhmm_text):
        self._check_workflow_call(xhmm_text, 'XHMM', 'XHMM')

    def test_indelible_workflow(self, indelible_text):
        self._check_workflow_call(indelible_text, 'INDELIBLE', 'INDELIBLE')

    def test_gatk_workflow(self, gatk_text):
        self._check_workflow_call(gatk_text, 'GATK_GCNV', 'GATK')

    def test_cnvkit_workflow(self, cnvkit_text):
        self._check_workflow_call(cnvkit_text, 'CNVKIT', 'CNVKIT')

    def test_dragen_workflow(self, dragen_text):
        self._check_workflow_call(dragen_text, 'DRAGEN', 'DRAGEN')


class TestWorkflowEmitsNormalised:
    """Each workflow emit block must expose normalised_vcf and normalised_vcf_index."""

    def _check_emit(self, module_text, workflow_name, caller_name):
        wf = _extract_workflow(module_text, workflow_name)
        assert wf is not None, f"{workflow_name} workflow not found"
        assert 'normalised_vcf' in wf, (
            f"{caller_name} workflow must emit normalised_vcf"
        )
        assert 'normalised_vcf_index' in wf, (
            f"{caller_name} workflow must emit normalised_vcf_index"
        )

    def test_canoes_emit(self, canoes_text):
        self._check_emit(canoes_text, 'CANOES', 'CANOES')

    def test_clamms_emit(self, clamms_text):
        self._check_emit(clamms_text, 'CLAMMS', 'CLAMMS')

    def test_xhmm_emit(self, xhmm_text):
        self._check_emit(xhmm_text, 'XHMM', 'XHMM')

    def test_indelible_emit(self, indelible_text):
        self._check_emit(indelible_text, 'INDELIBLE', 'INDELIBLE')

    def test_gatk_emit(self, gatk_text):
        self._check_emit(gatk_text, 'GATK_GCNV', 'GATK')

    def test_cnvkit_emit(self, cnvkit_text):
        self._check_emit(cnvkit_text, 'CNVKIT', 'CNVKIT')

    def test_dragen_emit(self, dragen_text):
        self._check_emit(dragen_text, 'DRAGEN', 'DRAGEN')


# ===========================================================================
# 5. nextflow.config has pysam container label
# ===========================================================================

class TestNextflowConfigPysamLabel:
    """nextflow.config must declare a 'pysam' label for the container."""

    def test_pysam_label_present(self, config_text):
        assert "withLabel: 'pysam'" in config_text, (
            "nextflow.config must contain a withLabel: 'pysam' block "
            "so that NORMALISE_CNV_QUALITY_SCORES processes can resolve a container"
        )

    def test_pysam_container_uses_pysam_image(self, config_text):
        pysam_block_match = re.search(
            r"withLabel:\s*'pysam'\s*\{(.+?)(?=\})", config_text, re.DOTALL
        )
        assert pysam_block_match, "pysam label block not found in nextflow.config"
        block = pysam_block_match.group(1)
        assert 'pysam' in block.lower(), (
            "The 'pysam' container image path/URI must reference pysam"
        )
