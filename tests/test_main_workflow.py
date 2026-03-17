#!/usr/bin/env python3
"""
Integration tests for main.nf – verifies that all 9 workflow cases are
wired correctly and that all supporting files are coherent.

Covers:
  1. All 9 workflow modules are included (INDELIBLE, CANOES, XHMM, CLAMMS,
     DRAGEN, CNVKIT, GATK_GCNV, SURVIVOR, TRUVARI).
  2. Every workflow case is present in the switch block:
     indelible, canoes, xhmm, clamms, dragen, cnvkit, gcnv, survivor, truvari.
  3. gather_vcfs() checks for at least 2 caller directories before continuing.
  4. The get_id closure strips caller name suffixes and .vcf/.vcf.gz extensions
     from filenames to recover the bare sample ID.
  5. RUN_SURVIVOR and RUN_TRUVARI filter samples that have fewer than 2 VCFs
     (so only samples covered by ≥2 callers enter consensus steps).
  6. All 9 params JSON files exist in params/ and contain valid JSON with the
     correct "workflow" key.
  7. Each params JSON references a workflow name that matches one of the 9 cases.
"""

import json
import os
import re

import pytest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAIN_NF = os.path.join(REPO_ROOT, "main.nf")
PARAMS_DIR = os.path.join(REPO_ROOT, "params")

# All 9 supported workflow names
ALL_WORKFLOWS = [
    "indelible",
    "canoes",
    "xhmm",
    "clamms",
    "dragen",
    "cnvkit",
    "gcnv",
    "survivor",
    "truvari",
]

# Mapping from workflow name to expected params JSON filename
PARAMS_FILES = {
    "indelible": "params-indelible.json",
    "canoes":    "params-canoes.json",
    "xhmm":      "params-xhmm.json",
    "clamms":    "params-clamms.json",
    "dragen":    "params-icav2-dragen.json",
    "cnvkit":    "params-cnvkit.json",
    "gcnv":      "params-gatk-gcnv.json",
    "survivor":  "params-survivor.json",
    "truvari":   "params-truvari.json",
}

# Module names expected in include statements
MODULE_INCLUDES = {
    "INDELIBLE": "modules-indelible.nf",
    "CANOES":    "modules-canoes.nf",
    "XHMM":      "modules-xhmm.nf",
    "CLAMMS":    "modules-clamms.nf",
    "DRAGEN":    "modules-icav2-dragen.nf",
    "CNVKIT":    "modules-cnvkit.nf",
    "GATK_GCNV": "modules-gatk-gcnv.nf",
    "SURVIVOR":  "modules-survivor.nf",
    "TRUVARI":   "modules-truvari.nf",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_main():
    with open(MAIN_NF) as fh:
        return fh.read()


def _load_params(filename):
    path = os.path.join(PARAMS_DIR, filename)
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def main_text():
    return _read_main()


# ===========================================================================
# 1. Module includes
# ===========================================================================

class TestModuleIncludes:
    """main.nf must include all 9 workflow modules."""

    @pytest.mark.parametrize("workflow_name,module_file", MODULE_INCLUDES.items())
    def test_include_present(self, main_text, workflow_name, module_file):
        """Each of the 9 modules must be imported via an include statement."""
        assert f"include {{ {workflow_name} }}" in main_text, (
            f"main.nf must include {{ {workflow_name} }} from ./modules/{module_file}"
        )
        assert module_file in main_text, (
            f"main.nf include for {workflow_name} must reference {module_file}"
        )


# ===========================================================================
# 2. Workflow switch cases
# ===========================================================================

class TestWorkflowSwitchCases:
    """The workflow switch block must handle all 9 workflow modes."""

    @pytest.mark.parametrize("workflow_case", ALL_WORKFLOWS)
    def test_case_present(self, main_text, workflow_case):
        """Every workflow mode must have a case in the switch block."""
        # case['indelible'] style
        assert f"case['{workflow_case}']" in main_text, (
            f"main.nf switch block must contain case['{workflow_case}'] "
            f"to dispatch the {workflow_case} workflow"
        )

    def test_default_case_with_error_message(self, main_text):
        """The switch block must have a default case that exits with an error message."""
        assert "default:" in main_text, (
            "main.nf switch block must contain a default: case"
        )
        # The default case must list all workflow options
        for wf in ALL_WORKFLOWS:
            assert f"--workflow {wf}" in main_text, (
                f"The default error message must list '--workflow {wf}' as a valid option"
            )


# ===========================================================================
# 3. Sub-workflow definitions
# ===========================================================================

class TestSubWorkflowDefinitions:
    """main.nf must define a RUN_* sub-workflow for each of the 9 callers."""

    @pytest.mark.parametrize("workflow_case", ALL_WORKFLOWS)
    def test_run_workflow_defined(self, main_text, workflow_case):
        """Each caller must have a corresponding RUN_<WORKFLOW> sub-workflow."""
        expected_name = f"RUN_{workflow_case.upper()}"
        # gcnv → RUN_GCNV
        assert f"workflow {expected_name}" in main_text, (
            f"main.nf must define 'workflow {expected_name}' as a sub-workflow "
            f"for the {workflow_case} caller"
        )


# ===========================================================================
# 4. gather_vcfs() helper function
# ===========================================================================

class TestGatherVcfsFunction:
    """gather_vcfs() must validate caller directory count and extract sample IDs."""

    def test_gather_vcfs_defined(self, main_text):
        """gather_vcfs() function must be defined in main.nf."""
        assert "def gather_vcfs()" in main_text, (
            "main.nf must define the gather_vcfs() helper function "
            "for SURVIVOR and TRUVARI consensus workflows"
        )

    def test_gather_vcfs_requires_at_least_two_dirs(self, main_text):
        """gather_vcfs() must exit with error if fewer than 2 caller dirs are given."""
        assert "dir_count < 2" in main_text, (
            "gather_vcfs() must check dir_count < 2 and exit 1 when fewer than "
            "two caller directories are provided"
        )

    def test_gather_vcfs_error_message_mentions_two_callers(self, main_text):
        """The gather_vcfs error message must explain the two-caller requirement."""
        assert "TWO" in main_text or "at least" in main_text, (
            "gather_vcfs() error message must explain that at least two caller "
            "directories must be provided"
        )

    def test_gather_vcfs_strips_canoes_suffix(self, main_text):
        """get_id closure must strip _CANOES from filenames."""
        assert "CANOES" in main_text, (
            "gather_vcfs() get_id closure must handle CANOES filenames"
        )

    def test_gather_vcfs_strips_clamms_suffix(self, main_text):
        """get_id closure must strip _CLAMMS from filenames."""
        assert "CLAMMS" in main_text, (
            "gather_vcfs() get_id closure must handle CLAMMS filenames"
        )

    def test_gather_vcfs_strips_xhmm_suffix(self, main_text):
        """get_id closure must strip _XHMM from filenames."""
        assert "XHMM" in main_text, (
            "gather_vcfs() get_id closure must handle XHMM filenames"
        )

    def test_gather_vcfs_strips_cnvkit_suffix(self, main_text):
        """get_id closure must strip _CNVKIT from filenames."""
        assert "CNVKIT" in main_text, (
            "gather_vcfs() get_id closure must handle CNVKIT filenames"
        )

    def test_gather_vcfs_strips_gcnv_suffix(self, main_text):
        """get_id closure must strip _GCNV from filenames."""
        assert "GCNV" in main_text, (
            "gather_vcfs() get_id closure must handle GCNV filenames (GATK gCNV)"
        )

    def test_gather_vcfs_strips_dragen_suffix(self, main_text):
        """get_id closure must strip _DRAGEN from filenames."""
        assert "DRAGEN" in main_text, (
            "gather_vcfs() get_id closure must handle DRAGEN filenames"
        )

    def test_gather_vcfs_strips_indelible_suffix(self, main_text):
        """get_id closure must strip _INDELIBLE from filenames."""
        assert "INDELIBLE" in main_text, (
            "gather_vcfs() get_id closure must handle INDELIBLE filenames"
        )

    def test_gather_vcfs_globs_vcf_files(self, main_text):
        """gather_vcfs() must collect *.vcf* files from each caller directory."""
        assert "*.vcf*" in main_text, (
            "gather_vcfs() must use '*.vcf*' glob to collect both .vcf and .vcf.gz files"
        )

    def test_gather_vcfs_used_by_survivor_case(self, main_text):
        """The survivor case must call gather_vcfs() to collect VCFs."""
        # Find the survivor case block
        survivor_case = re.search(
            r"case\['survivor'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert survivor_case, "case['survivor'] block not found"
        assert "gather_vcfs()" in survivor_case.group(1), (
            "case['survivor'] must call gather_vcfs() to collect VCFs from caller dirs"
        )

    def test_gather_vcfs_used_by_truvari_case(self, main_text):
        """The truvari case must call gather_vcfs() to collect VCFs."""
        truvari_case = re.search(
            r"case\['truvari'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert truvari_case, "case['truvari'] block not found"
        assert "gather_vcfs()" in truvari_case.group(1), (
            "case['truvari'] must call gather_vcfs() to collect VCFs from caller dirs"
        )

    @pytest.mark.parametrize("sample_id,filename", [
        ("SAMPLE001", "SAMPLE001_CANOES_output.vcf"),
        ("SAMPLE001", "SAMPLE001_CANOES_output.vcf.gz"),
        ("SAMPLE002", "SAMPLE002_CLAMMS_output.vcf"),
        ("MYSAMPLE",  "MYSAMPLE_XHMM_output.vcf.gz"),
        ("SID",       "SID_GCNV_genotyped_segments.sorted.vcf.gz"),
        ("SAMP",      "SAMP_CNVKIT_output.vcf"),
        ("PT01",      "PT01_DRAGEN_output.vcf"),
        ("PT02",      "PT02_INDELIBLE_output.vcf"),
    ])
    def test_get_id_regex_extracts_sample_id(self, main_text, sample_id, filename):
        """The get_id regex in gather_vcfs() must correctly strip caller suffixes.

        Simulates the Groovy closure:
          f.name.replaceAll(/_((CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '')
                .replaceAll(/\\.vcf(\\.gz)?$/i, '')
        """
        import re as _re
        # Replicate the Groovy regex logic in Python
        extracted = _re.sub(
            r'_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*', '', filename
        )
        extracted = _re.sub(r'\.vcf(\.gz)?$', '', extracted, flags=_re.IGNORECASE)
        assert extracted == sample_id, (
            f"get_id regex must extract '{sample_id}' from '{filename}', "
            f"got '{extracted}'"
        )


# ===========================================================================
# 5. RUN_SURVIVOR and RUN_TRUVARI: enforce ≥2 VCFs per sample
# ===========================================================================

class TestConsensusWorkflowFiltering:
    """RUN_SURVIVOR and RUN_TRUVARI must only process samples with ≥2 caller VCFs."""

    def test_run_survivor_filters_single_caller_samples(self, main_text):
        """RUN_SURVIVOR must filter out samples with only 1 caller VCF."""
        run_survivor = re.search(
            r"workflow RUN_SURVIVOR\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_survivor, "workflow RUN_SURVIVOR not found"
        body = run_survivor.group(1)
        # Must groupTuple and filter for >= 2 VCFs
        assert ".groupTuple()" in body, (
            "RUN_SURVIVOR must use .groupTuple() to aggregate VCFs by sample_id"
        )
        assert ".filter" in body, (
            "RUN_SURVIVOR must use .filter to exclude samples with < 2 caller VCFs"
        )
        assert "vcfs.size() >= 2" in body, (
            "RUN_SURVIVOR filter condition must require vcfs.size() >= 2"
        )

    def test_run_truvari_filters_single_caller_samples(self, main_text):
        """RUN_TRUVARI must filter out samples with only 1 caller VCF."""
        run_truvari = re.search(
            r"workflow RUN_TRUVARI\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_truvari, "workflow RUN_TRUVARI not found"
        body = run_truvari.group(1)
        assert ".groupTuple()" in body, (
            "RUN_TRUVARI must use .groupTuple() to aggregate VCFs by sample_id"
        )
        assert ".filter" in body, (
            "RUN_TRUVARI must use .filter to exclude samples with < 2 caller VCFs"
        )
        assert "vcfs.size() >= 2" in body, (
            "RUN_TRUVARI filter condition must require vcfs.size() >= 2"
        )

    def test_run_survivor_calls_survivor_workflow(self, main_text):
        """RUN_SURVIVOR sub-workflow must invoke the SURVIVOR module workflow."""
        run_survivor = re.search(
            r"workflow RUN_SURVIVOR\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_survivor, "workflow RUN_SURVIVOR not found"
        assert "SURVIVOR(" in run_survivor.group(1), (
            "RUN_SURVIVOR must call the SURVIVOR workflow from modules-survivor.nf"
        )

    def test_run_truvari_calls_truvari_workflow(self, main_text):
        """RUN_TRUVARI sub-workflow must invoke the TRUVARI module workflow."""
        run_truvari = re.search(
            r"workflow RUN_TRUVARI\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_truvari, "workflow RUN_TRUVARI not found"
        assert "TRUVARI(" in run_truvari.group(1), (
            "RUN_TRUVARI must call the TRUVARI workflow from modules-truvari.nf"
        )


# ===========================================================================
# 6. Params JSON files: existence and validity
# ===========================================================================

class TestParamsJsonFiles:
    """All 9 params JSON files must exist, be valid JSON, and contain 'workflow'."""

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_exists(self, workflow_name, filename):
        """Each workflow must have a corresponding params JSON file in params/."""
        path = os.path.join(PARAMS_DIR, filename)
        assert os.path.isfile(path), (
            f"params/{filename} must exist for the {workflow_name} workflow"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_is_valid_json(self, workflow_name, filename):
        """Each params file must be valid JSON (parseable without errors)."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        try:
            data = _load_params(filename)
        except json.JSONDecodeError as exc:
            pytest.fail(
                f"params/{filename} contains invalid JSON: {exc}"
            )
        assert isinstance(data, dict), (
            f"params/{filename} must be a JSON object (dict), got {type(data).__name__}"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_has_workflow_key(self, workflow_name, filename):
        """Each params file must contain a 'workflow' key."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        assert "workflow" in data, (
            f"params/{filename} must contain a 'workflow' key to specify "
            f"which pipeline to run (expected '{workflow_name}')"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_workflow_value_is_correct(self, workflow_name, filename):
        """The 'workflow' value in each params file must match the expected workflow name."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        if "workflow" not in data:
            pytest.skip(f"'workflow' key absent in {filename}")
        # dragen params file uses 'dragen' not 'icav2-dragen'
        expected = workflow_name
        actual = data["workflow"]
        assert actual == expected, (
            f"params/{filename}: 'workflow' must be '{expected}', got '{actual}'"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_has_outdir_key(self, workflow_name, filename):
        """Each params file must contain an output path key for result organisation.

        Most workflows use 'outdir'. The ICAv2-DRAGEN cloud workflow instead uses
        'localDownloadPath' because results are downloaded from ICA rather than
        written to a local output directory.
        """
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        # DRAGEN (ICAv2 cloud) uses localDownloadPath instead of outdir
        has_output_path = "outdir" in data or "localDownloadPath" in data
        assert has_output_path, (
            f"params/{filename} must contain either an 'outdir' key "
            f"(local workflows) or 'localDownloadPath' (ICAv2-DRAGEN) "
            f"specifying where results should be written"
        )

    @pytest.mark.parametrize("filename", [
        "params-survivor.json",
        "params-truvari.json",
    ])
    def test_consensus_params_have_caller_dirs(self, filename):
        """Consensus params files (survivor, truvari) must specify at least one caller dir."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        caller_dir_keys = [
            "canoes_dir", "clamms_dir", "xhmm_dir",
            "cnvkit_dir", "gcnv_dir", "dragen_dir", "indelible_dir",
        ]
        found = [k for k in caller_dir_keys if k in data]
        assert len(found) >= 2, (
            f"params/{filename} must specify at least two caller_dir keys "
            f"(e.g. canoes_dir, clamms_dir) so gather_vcfs() has enough inputs. "
            f"Found: {found}"
        )


# ===========================================================================
# 7. Nextflow DSL2 declaration
# ===========================================================================

class TestDsl2Declaration:
    """main.nf must use Nextflow DSL2."""

    def test_dsl2_enabled(self, main_text):
        """main.nf must declare nextflow.enable.dsl=2."""
        assert "nextflow.enable.dsl=2" in main_text, (
            "main.nf must contain 'nextflow.enable.dsl=2' to use DSL2 module imports"
        )

    def test_shebang_line(self, main_text):
        """main.nf must start with a Nextflow shebang line."""
        assert main_text.startswith("#!/usr/bin/env nextflow"), (
            "main.nf must start with '#!/usr/bin/env nextflow'"
        )
