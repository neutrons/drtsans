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

INSTRUMENT = 'GPSANS'  # From 'EQSANS', 'BIOSANS'

IPTS = 24664
FLOOD_RUNS = [1697, 1701, 1699]

DIRECT_BEAM_RUNS = [1698, 1702, 1700]
MASK_BEAM_CENTER_RADIUS = 0.5  # mm

BEAM_CENTER_MASKS = None

UNIVERSAL_MASK = 'Mask.XML'

# If it is GPSANS, there could be 2 options
MOVING_DETECTORS = True

# Normalize data
NORMALIZE_COUNTS = True  # BIO/EQ-SANS may not, But GPSANS does regardless

#
prepare_data()

calculate_sensitivities()

export()


def main():
    pass
