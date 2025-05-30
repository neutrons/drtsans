# flake8: noqa
from pathlib import Path

try:
    from ._version import __version__  # noqa: F401
except ImportError:
    __version__ = "unknown"

#######
# Ordered alphabetically
#######
from .api import *
from .auto_wedge import *
from .beam_finder import *
from .instruments import *
from .momentum_transfer import *
from .polarization import *
from .reductionlog import *
from .resolution import *
from .sensitivity import *
from .solid_angle import *
from .thickness_normalization import *

# FIXME the functions done as strings can't be done via __all__ because module and function have same name
__all__ = (
    ["convert_to_q", "solid_angle_correction"]
    + api.__all__
    + beam_finder.__all__
    + auto_wedge.__all__
    + instruments.__all__
    + polarization.__all__
    + reductionlog.__all__
    + thickness_normalization.__all__
    + resolution.__all__
    + sensitivity.__all__
)

# directory where to put non-source
configdir = str(Path(__file__).parent / "configuration")
scriptsdir = str(Path(__file__).parent.parent / "scripts")
