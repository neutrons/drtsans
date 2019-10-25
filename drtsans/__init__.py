from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .api import *  # noqa: F403
from .beam_finder import *  # noqa: F403
from .geometry import *  # noqa: F403
from .sensitivity import *  # noqa: F403
from .solid_angle_correction import *  # noqa: F403
from .savereductionlog import *  # noqa: F403
from .convert_to_q import *  # noqa: F403


# FIXME the functions done as strings can't be done via __all__ because module and function have same name
__all__ = (['convert_to_q', 'savereductionlog', 'solid_angle_correction']
           + api.__all__ + beam_finder.__all__ + geometry.__all__ + sensitivity.__all__)  # noqa: F405
