from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .api import *  # noqa: F403
from .sensitivity import *  # noqa: F403
from .solid_angle_correction import *  # noqa: F403
from .savereductionlog import *  # noqa: F403
from .convert_to_q import *  # noqa: F403

__all__ = (api.__all__, convert_to_q.__all__, savereductionlog.__all__, sensitivity.__all__,  # noqa: F405
           solid_angle_correction.__all__)  # noqa: F405
