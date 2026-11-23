import json
import subprocess
import sys

# Modules that used to install a process-wide "ignore" filter, and the instrument APIs that import them
MODULES = [
    "drtsans",
    "drtsans.redparams",
    "drtsans.plots.api",
    "drtsans.tof.eqsans",
    "drtsans.mono.gpsans",
    "drtsans.mono.biosans",
    "drtsans.mono.biosans.reduce_spice",
]

SCRIPT = f"""
import json
import warnings
for module in {MODULES!r}:
    __import__(module)
broad = [
    category.__name__
    for action, message, category, module, _ in warnings.filters
    if action == "ignore" and message is None and module is None and category in (Warning, UserWarning)
]
print(json.dumps(broad))
"""


def test_import_installs_no_broad_ignore_filter():
    """Importing drtsans must not hide Python warnings from users for the whole process.

    The imports run in a subprocess, because pytest manages ``warnings.filters`` in the test process and the
    modules are already imported there.
    """
    result = subprocess.run([sys.executable, "-c", SCRIPT], capture_output=True, text=True, check=True)
    broad_filters = json.loads(result.stdout.strip().splitlines()[-1])
    assert broad_filters == []
