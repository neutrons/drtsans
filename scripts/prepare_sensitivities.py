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
FLOOD_RUNS = [1697, 1701, 1699]  # Single value integer or a list or tuple

DIRECT_BEAM_RUNS = [1698, 1702, 1700]
MASK_BEAM_CENTER_RADIUS = 0.5  # mm

BEAM_CENTER_MASKS = None

UNIVERSAL_MASK = 'Mask.XML'

# If it is GPSANS, there could be 2 options
MOVING_DETECTORS = True

# Normalize data
NORMALIZE_COUNTS = True  # BIO/EQ-SANS may not, But GPSANS does regardless

# END OF USER INPUTS

def prepare_data():
    from drtsans.mono.gpsans.api import prepare_data
    from drtsans.mono.biosans.api import prepare_data
    from drtsans.tof.eqsans.api import prepare_data

    prepare_data(mask=UNIVERSAL_MASK)




if INSTRUMENT.lower() == 'gpsans' and MOVING_DETECTORS is True:
    # Use the moving detector algorithm
    pass
else:
    # Use the sensitivity patch method
    pass
#
prepare_data()

calculate_sensitivities()

export()


def main():
    pass
