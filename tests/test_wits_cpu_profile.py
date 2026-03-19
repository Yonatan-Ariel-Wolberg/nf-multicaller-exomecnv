#!/usr/bin/env python3
"""
Tests for CPU-node execution on the ZA-Wits-Core HPC cluster.

Validates that nextflow.config is correctly configured for:
  1. Apptainer (the container runtime on the Wits cluster, successor to Singularity).
     - An `apptainer { ... }` scope must be defined and enabled.
     - Both the `singularity` and `apptainer` cacheDir settings must honour the
       NXF_SINGULARITY_CACHEDIR environment variable so that a shared image cache
       on /spaces can be used across all worker nodes.
  2. No hardcoded user-specific bind-mount paths in the default runOptions.
     - params.runOptions must default to an empty string (users customise it
       in their own config or params-file).
     - No paths like /home/<user> or /dataB/aux appear in the default
       singularity/apptainer runOptions.
  3. The `wits` profile configures the SLURM executor and provides
     Singularity/Apptainer cache paths for the cluster.
  4. No GPU-specific clusterOptions (e.g., --gres=gpu) appear in the default
     process block or the `wits` profile.  The cluster is tested CPU-only.
  5. The `icav2-dragen` process label uses parameterised values for the container
     path and ICA credentials directory instead of hardcoded user-specific paths.
"""

import os
import re

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
NEXTFLOW_CONFIG = os.path.join(REPO_ROOT, 'nextflow.config')


def _read_config():
    with open(NEXTFLOW_CONFIG) as fh:
        return fh.read()


def _global_params_section(content):
    """Return the text of the global params { } block (before 'profiles {')."""
    return content.split('profiles {')[0]


def _extract_block(content, start_pattern):
    """Return the body of the first block matching start_pattern { ... }.

    Uses brace counting to handle nested blocks correctly.
    """
    match = re.search(start_pattern + r'\s*\{', content)
    assert match, f"Block matching '{start_pattern}' not found in nextflow.config"
    start = match.end()  # position just after the opening {
    depth = 1
    pos = start
    while pos < len(content) and depth > 0:
        if content[pos] == '{':
            depth += 1
        elif content[pos] == '}':
            depth -= 1
        pos += 1
    return content[start:pos - 1]


def _wits_profile_section(content):
    """Return the text of the wits { ... } profile block."""
    # Match `wits` only inside profiles { ... }
    profiles_idx = content.find('profiles {')
    assert profiles_idx >= 0, "profiles { } block not found in nextflow.config"
    return _extract_block(content[profiles_idx:], r'\bwits')


def _process_section(content):
    """Return the text of the global process { } block (before profiles {)."""
    global_part = _global_params_section(content)
    match = re.search(r'\bprocess\s*\{(.+)', global_part, re.DOTALL)
    assert match, "process block not found before profiles section"
    return match.group(1)


def _icav2_label_section(content):
    """Return the text of the withLabel: 'icav2-dragen' { } block."""
    match = re.search(
        r"withLabel:\s*'icav2-dragen'\s*\{([^}]+)\}", content, re.DOTALL
    )
    assert match, "withLabel: 'icav2-dragen' block not found"
    return match.group(1)


# ===========================================================================
# 1. Apptainer scope
# ===========================================================================

class TestApptainerScope:
    """nextflow.config must define an apptainer { } scope."""

    def test_apptainer_scope_exists(self):
        """An apptainer { } block must be present."""
        content = _read_config()
        assert re.search(r'\bapptainer\s*\{', content), (
            "nextflow.config must define an 'apptainer { }' scope for Apptainer "
            "(the container runtime used on the Wits Rocky 9 cluster)"
        )

    def test_apptainer_enabled(self):
        """apptainer.enabled must be set to true."""
        content = _read_config()
        apptainer_block = re.search(
            r'\bapptainer\s*\{([^}]+)\}', content, re.DOTALL
        )
        assert apptainer_block, "apptainer block not found"
        assert 'enabled' in apptainer_block.group(1), (
            "apptainer block must set 'enabled = true'"
        )
        assert re.search(r'enabled\s*=\s*true', apptainer_block.group(1)), (
            "apptainer.enabled must be true"
        )

    def test_apptainer_automounts(self):
        """apptainer.autoMounts must be set to true."""
        content = _read_config()
        apptainer_block = re.search(
            r'\bapptainer\s*\{([^}]+)\}', content, re.DOTALL
        )
        assert apptainer_block, "apptainer block not found"
        assert re.search(r'autoMounts\s*=\s*true', apptainer_block.group(1)), (
            "apptainer.autoMounts must be true so that host paths accessible "
            "to the job are automatically available inside the container"
        )

    def test_apptainer_cachedir_uses_env_var(self):
        """apptainer.cacheDir must reference NXF_SINGULARITY_CACHEDIR."""
        content = _read_config()
        apptainer_block = re.search(
            r'\bapptainer\s*\{([^}]+)\}', content, re.DOTALL
        )
        assert apptainer_block, "apptainer block not found"
        assert 'NXF_SINGULARITY_CACHEDIR' in apptainer_block.group(1), (
            "apptainer.cacheDir must use NXF_SINGULARITY_CACHEDIR so that all "
            "worker nodes share a single image store on a shared filesystem "
            "(e.g. /spaces/<project>/nxf_cache)"
        )


# ===========================================================================
# 2. Singularity scope – NXF_SINGULARITY_CACHEDIR
# ===========================================================================

class TestSingularityCacheDir:
    """singularity.cacheDir must honour NXF_SINGULARITY_CACHEDIR."""

    def test_singularity_cachedir_uses_env_var(self):
        """singularity.cacheDir must reference the NXF_SINGULARITY_CACHEDIR env var."""
        content = _read_config()
        singularity_block = re.search(
            r'\bsingularity\s*\{([^}]+)\}', content, re.DOTALL
        )
        assert singularity_block, "singularity block not found"
        assert 'NXF_SINGULARITY_CACHEDIR' in singularity_block.group(1), (
            "singularity.cacheDir must use the NXF_SINGULARITY_CACHEDIR "
            "environment variable so images are cached on a shared volume "
            "accessible by all worker nodes"
        )

    def test_singularity_cachedir_has_fallback(self):
        """singularity.cacheDir must provide a fallback when env var is not set."""
        content = _read_config()
        singularity_block = re.search(
            r'\bsingularity\s*\{([^}]+)\}', content, re.DOTALL
        )
        assert singularity_block, "singularity block not found"
        # Groovy Elvis operator ?: provides the fallback
        assert '?:' in singularity_block.group(1), (
            "singularity.cacheDir must use the Groovy Elvis operator (?:) to "
            "fall back to a default path when NXF_SINGULARITY_CACHEDIR is unset"
        )


# ===========================================================================
# 3. No hardcoded user-specific paths in default runOptions
# ===========================================================================

class TestNoHardcodedUserPaths:
    """Default params.runOptions must be empty; no user-specific paths."""

    def test_run_options_default_is_empty(self):
        """params.runOptions must default to an empty string."""
        content = _read_config()
        global_section = _global_params_section(content)
        match = re.search(r"runOptions\s*=\s*'([^']*)'", global_section)
        if not match:
            match = re.search(r'runOptions\s*=\s*"([^"]*)"', global_section)
        assert match is not None, (
            "params.runOptions must be declared in the global params block"
        )
        value = match.group(1)
        assert value == '', (
            f"params.runOptions default must be '' (empty string) so users "
            f"are not forced to use Wits-specific bind mounts. Got: '{value}'"
        )

    def test_no_hardcoded_home_ywolberg(self):
        """No /home/ywolberg paths should appear in default singularity runOptions."""
        content = _read_config()
        # Check the global singularity block (before profiles)
        global_section = _global_params_section(content)
        assert '/home/ywolberg' not in global_section, (
            "nextflow.config must not contain hardcoded '/home/ywolberg' paths "
            "in the default (pre-profiles) section"
        )

    def test_no_hardcoded_datab_aux(self):
        """No /dataB/aux paths should appear in the default singularity runOptions."""
        content = _read_config()
        global_section = _global_params_section(content)
        assert '/dataB/aux' not in global_section, (
            "nextflow.config must not hardcode '/dataB/aux' in the default "
            "runOptions (it belongs in a user-specific config override)"
        )

    def test_no_hardcoded_datag_ddd(self):
        """No /dataG/ddd paths should appear in the default singularity runOptions."""
        content = _read_config()
        global_section = _global_params_section(content)
        assert '/dataG/ddd' not in global_section, (
            "nextflow.config must not hardcode '/dataG/ddd' in the default "
            "runOptions (it belongs in a user-specific config override)"
        )


# ===========================================================================
# 4. wits profile – SLURM executor and Singularity/Apptainer cache
# ===========================================================================

class TestWitsProfile:
    """The wits profile must configure the SLURM executor and image cache."""

    def test_wits_profile_exists(self):
        """A 'wits' profile must be defined in profiles { }."""
        content = _read_config()
        profiles_section = content.split('profiles {', 1)[-1]
        assert re.search(r'\bwits\s*\{', profiles_section), (
            "nextflow.config must define a 'wits' profile inside profiles { }"
        )

    def test_wits_uses_slurm_executor(self):
        """wits profile must set executor.name = 'slurm'."""
        content = _read_config()
        wits = _wits_profile_section(content)
        assert "executor.name" in wits and "slurm" in wits, (
            "wits profile must set executor.name = 'slurm' to use the "
            "SLURM workload manager on the Wits cluster"
        )

    def test_wits_sets_batch_queue(self):
        """wits profile must direct jobs to the 'batch' SLURM queue."""
        content = _read_config()
        wits = _wits_profile_section(content)
        assert "queue" in wits and "batch" in wits, (
            "wits profile must set process.queue = 'batch' for SLURM job submission"
        )

    def test_wits_singularity_cachedir_uses_env_var(self):
        """wits profile must set singularity.cacheDir using NXF_SINGULARITY_CACHEDIR."""
        content = _read_config()
        wits = _wits_profile_section(content)
        assert 'singularity.cacheDir' in wits, (
            "wits profile must set singularity.cacheDir so images are cached "
            "on a shared volume accessible to all cluster worker nodes"
        )
        assert 'NXF_SINGULARITY_CACHEDIR' in wits, (
            "wits profile's singularity.cacheDir must use NXF_SINGULARITY_CACHEDIR"
        )

    def test_wits_apptainer_cachedir_uses_env_var(self):
        """wits profile must set apptainer.cacheDir using NXF_SINGULARITY_CACHEDIR."""
        content = _read_config()
        wits = _wits_profile_section(content)
        assert 'apptainer.cacheDir' in wits, (
            "wits profile must set apptainer.cacheDir so Apptainer images are "
            "cached on a shared volume accessible to all cluster worker nodes"
        )
        assert 'NXF_SINGULARITY_CACHEDIR' in wits, (
            "wits profile's apptainer.cacheDir must use NXF_SINGULARITY_CACHEDIR"
        )

    def test_wits_no_gpu_cluster_options(self):
        """wits profile must not contain --gres=gpu (CPU-only testing)."""
        content = _read_config()
        wits = _wits_profile_section(content)
        assert '--gres=gpu' not in wits, (
            "wits profile must not request GPU resources (--gres=gpu). "
            "The current testing environment uses CPU-only worker nodes."
        )


# ===========================================================================
# 5. No GPU clusterOptions in the default process block
# ===========================================================================

class TestNoGpuInDefaultProcess:
    """The global process { } block must not contain GPU-specific clusterOptions."""

    def test_no_gres_gpu_in_default_process(self):
        """--gres=gpu must not appear in the default (non-profile) process block."""
        content = _read_config()
        proc = _process_section(content)
        assert '--gres=gpu' not in proc, (
            "The default process block must not contain '--gres=gpu'. "
            "GPU resources should only be requested in a dedicated GPU profile."
        )

    def test_no_gpu_reservation_in_default_process(self):
        """--reservation=gpu must not appear in the default process block."""
        content = _read_config()
        proc = _process_section(content)
        assert '--reservation=gpu' not in proc, (
            "The default process block must not contain '--reservation=gpu'. "
            "GPU reservations should only be requested in a dedicated GPU profile."
        )


# ===========================================================================
# 6. ICAv2-DRAGEN label – parameterised, no hardcoded user paths
# ===========================================================================

class TestIcav2DragenLabel:
    """The icav2-dragen label must use params for user-specific paths."""

    def test_no_hardcoded_home_ywolberg_in_icav2_label(self):
        """icav2-dragen label must not hardcode /home/ywolberg."""
        content = _read_config()
        label = _icav2_label_section(content)
        assert '/home/ywolberg' not in label, (
            "withLabel: 'icav2-dragen' must not hardcode '/home/ywolberg'. "
            "Use params.icav2_creds_dir or System.getProperty('user.home') instead."
        )

    def test_icav2_container_uses_param(self):
        """icav2-dragen container must reference a parameter, not a hardcoded path."""
        content = _read_config()
        label = _icav2_label_section(content)
        # The container line must reference params (not hardcode /home/ywolberg)
        container_line = [l for l in label.splitlines() if 'container' in l]
        assert container_line, "icav2-dragen label must set container"
        assert 'params' in container_line[0] or 'icav2_container' in container_line[0], (
            "icav2-dragen container must be set via params.icav2_container "
            "rather than a hardcoded /home/... path"
        )

    def test_icav2_creds_dir_uses_param_or_system(self):
        """icav2-dragen runOptions must reference params or System for creds dir."""
        content = _read_config()
        label = _icav2_label_section(content)
        assert 'icav2_creds_dir' in label or "user.home" in label, (
            "icav2-dragen runOptions must reference params.icav2_creds_dir or "
            "System.getProperty('user.home') for the ICA credentials directory "
            "instead of hardcoding a user-specific path"
        )
