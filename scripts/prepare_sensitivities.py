"""
    SANS sensitivities preparation script
"""
import os
import sys
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import mantid.simpleapi as msapi  # noqa E402

import drtsans  # noqa E402
from drtsans.mono import gpsans as sans  # noqa E402
from drtsans.settings import unique_workspace_dundername as uwd  # noqa E402

INSTRUMENT = 'GPSANS'

# If it is GPSANS, there could be 2 options
MOVING_DETECTORS = True


def main():
    pass
