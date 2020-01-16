"""
    SANS sensitivities preparation script
"""
import sys
import warnings
from mantid.simpleapi import SaveNexusProcessed
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
warnings.simplefilter(action="ignore", category=FutureWarning)


INSTRUMENT = 'GPSANS'  # From 'EQSANS', 'BIOSANS'

# Input Flood Runs
FLOOD_RUNS = ['1697', '1701', '1699']  # Single value integer or a list or tuple

# Output
SENSITIVITY_FILE = '/HFIR/CG2/shared/sensitivity1697.nxs'

# About Masks
DIRECT_BEAM_RUNS = ['1698', '1702', '1700']
MASK_BEAM_CENTER_RADIUS = 0.5  # mm
BEAM_CENTER_MASKS = None

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
MASKED_PIXELS = '1-8,249-256'

# If it is GPSANS, there could be 2 options
MOVING_DETECTORS = True

# THRESHOLD
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 2.0


# END OF USER INPUTS


# Load data files
if INSTRUMENT == 'GPSANS':
    from drtsans.mono.gpsans.api import prepare_data
    import drtsans.mono.gpsans as mysans
elif INSTRUMENT == 'BIOSANS':
    from drtsans.mono.biosans.api import prepare_data
    import drtsans.mono.biosans as mysans
elif INSTRUMENT == 'EQSANS':
    from drtsans.tof.eqsans.api import prepare_data
    import drtsans.tof.eqsans as mysans
else:
    print('Instrument {} is not supported.  Supported are {}'
          ''.format(INSTRUMENT, 'GPSANS, EQSANS, BIOSANS'))
    sys.exit(-1)

# Process runs
if isinstance(FLOOD_RUNS, int):
    sans_runs = [FLOOD_RUNS]
else:
    sans_runs = list(FLOOD_RUNS)
if DIRECT_BEAM_RUNS is not None:
    sans_runs.extend(list(DIRECT_BEAM_RUNS))

# Extra masks
extra_mask_dict = dict()
if MASKED_PIXELS is not None:
    extra_mask_dict['Pixel'] = MASKED_PIXELS

# Load data with masking: returning to a list of workspace references
# processing includes: load, mask, normalize by monitor
workspaces = [prepare_data(data=sans_runs[i],
                           mask=UNIVERSAL_MASK, btp=extra_mask_dict,
                           flux_method='monitor') for i in range(len(sans_runs))]

# Decide algorithm to prepare sensitivities
if INSTRUMENT == 'GPSANS' and MOVING_DETECTORS is True:
    # Prepare by using moving detector algorithm
    from drtsans.mono.gpsans.prepare_sensitivity import prepare_sensitivity_correction

    # Default for GPSANS
    if MASK_BEAM_CENTER_RADIUS is None:
        MASK_BEAM_CENTER_RADIUS = 65  # mm

    # Calculate sensitivities for each file
    sens_ws = prepare_sensitivity_correction(workspaces,
                                             threshold_min=MIN_THRESHOLD,
                                             threshold_max=1.5,
                                             beam_center_radius=MASK_BEAM_CENTER_RADIUS)

else:
    # Prepare by Use the sensitivity patch method
    from drtsans.sensitivity import prepare_sensitivity_correction

    # Mask detector center
    # There is one and only 1 data workspace
    data_ws = workspaces[0]
    if BEAM_CENTER_MASKS:
        # mask beam center by Mantid XML mask file
        apply_mask(data_ws, mask=BEAM_CENTER_MASKS)
    else:
        # find beam center to mask
        xc, yc = mysans.find_beam_center(data_ws)
        # Center detector to the data workspace (change in geometry)
        mysans.center_detector(data_ws, xc, yc)
        # Mask the new beam center by 65 mm (Lisa's magic number)
        det = list(circular_mask_from_beam_center(data_ws, MASK_BEAM_CENTER_RADIUS))
        apply_mask(data_ws, mask=det)

    # working on 1 and only 1
    sens_ws = prepare_sensitivity_correction(data_ws,
                                             min_threshold=MIN_THRESHOLD,
                                             max_threshold=MAX_THRESHOLD)

# Export
SaveNexusProcessed(InputWorkspace=sens_ws, Filename=SENSITIVITY_FILE)
