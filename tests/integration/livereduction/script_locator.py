import importlib.util
from pathlib import Path

# Add the root directory to the path
_root_dir = Path(__file__).parent.parent.parent.parent  # Go up 4 levels from test file


def _load_module(script_file: str) -> importlib.machinery.ModuleSpec:
    """
    Load a Python module dynamically from the specified file.

    Parameters
    ----------
    script_file : str
        The name of the script file to load (relative to the repository's root directory).

    Returns
    -------
    importlib.machinery.ModuleSpec
        The loaded module object.

    Notes
    -----
    - The function constructs the full path to the script file based on the project structure.
    - The module is loaded using `importlib.util` to allow dynamic imports.
    """
    script_file = Path(script_file)
    script_path = _root_dir / script_file
    spec = importlib.util.spec_from_file_location(script_file.stem, script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reduce_EQSANS_live_post_proc = _load_module("scripts/livereduction/eqsans/reduce_EQSANS_live_post_proc.py")
reduce_EQSANS_posixpath = _root_dir / "scripts" / "autoreduction" / "reduce_EQSANS.py"
