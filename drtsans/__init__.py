# flake8: noqa
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .api import *
from .beam_finder import *
from .instruments import *
from .sensitivity import *
from .solid_angle_correction import *
from .reductionlog import *
from .convert_to_q import *
from .resolution import *


# FIXME the functions done as strings can't be done via __all__ because module and function have same name
__all__ = ['convert_to_q', 'solid_angle_correction'] + api.__all__ + beam_finder.__all__ + instruments.__all__ +\
          reductionlog.__all__ + sensitivity.__all__
