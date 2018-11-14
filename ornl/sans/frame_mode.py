from __future__ import (absolute_import, division, print_function)

from enum import Enum


class FrameMode(Enum):
    r"""
    Selects if instrument operating in frame-skipping mode
    """
    not_skip = 0
    skip = 1
