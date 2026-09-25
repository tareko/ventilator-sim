"""The JavaScript port must agree with the Python model, case for case.

These are the checks that make `web/ventsim.html` the same simulator rather than
an imitation of one.  They are skipped if node is not on the PATH.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
node = pytest.mark.skipif(shutil.which("node") is None, reason="node not available")

JS_FILES = [ROOT / "web" / "ventsim.js", ROOT / "web" / "cases.json",
            ROOT / "tools" / "check_parity.js"]


@node
@pytest.mark.parametrize("path", JS_FILES, ids=lambda p: str(p.name))
def test_port_files_exist(path: Path):
    assert path.exists(), f"{path.name} is missing - regenerate the JS port"


@node
def test_javascript_matches_python_case_by_case():
    proc = subprocess.run(["node", "tools/check_parity.js"], cwd=ROOT,
                          capture_output=True, text=True)
    tail = "\n".join(proc.stdout.strip().split("\n")[-6:])
    assert proc.returncode == 0, f"parity failed:\n{tail}\n{proc.stderr[-800:]}"
    assert "PARITY OK" in proc.stdout


@node
def test_html_scenarios_mean_the_same_patients():
    dump = subprocess.run(["python3", "tools/dump_scenarios.py",
                           "web/scenarios.json"], cwd=ROOT,
                          capture_output=True, text=True)
    assert dump.returncode == 0, dump.stdout + dump.stderr
    proc = subprocess.run(["node", "tools/check_scenarios.js",
                           "web/scenarios.json"], cwd=ROOT,
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "scenarios agree" in proc.stdout


@node
def test_javascript_has_no_syntax_errors():
    for name in ["web/ventsim.js", "tools/check_parity.js", "tools/check_scenarios.js"]:
        proc = subprocess.run(["node", "--check", name], cwd=ROOT,
                              capture_output=True, text=True)
        assert proc.returncode == 0, f"{name}: {proc.stderr}"
