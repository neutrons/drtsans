from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .api import half_polarization, subtract_background  # noqa: F401
from .sensitivity import (apply_sensitivity_correction,  # noqa: F401
                          calculate_sensitivity_correction)
from .solid_angle_correction import *  # noqa: F403, F401
from .savereductionlog import *  # noqa: F401, F403

__all__ = ['apply_sensitivity_correction', 'calculate_sensitivity_correction', 'half_polarization',  # noqa: F405
           'savereductionlog', 'solid_angle_correction', 'subtract_background']  # noqa: F405
