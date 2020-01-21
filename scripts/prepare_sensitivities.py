"""
    SANS sensitivities preparation script

    # goal
    1. implement a universal mask_beam_center(flood_ws, beam_center_mask=None, beam_center_ws=None)
       for 3 types of mask
    2. add option for wing/main detector for BIOSANS:w


"""
import sys
import warnings
from mantid.simpleapi import SaveNexusProcessed
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
from drtsans.process_uncertainties import set_init_uncertainties
warnings.simplefilter(action="ignore", category=FutureWarning)


INSTRUMENT = 'CG2'  # From 'EQSANS', 'CG3'

# Input Flood Runs
FLOOD_RUNS = 1697, 1701, 1699  # Single value integer or a list or tuple

# Output
SENSITIVITY_FILE = '/HFIR/CG2/shared/sensitivity1697.nxs'

# About Masks
DIRECT_BEAM_RUNS = 1698, 1702, 1700
MASK_BEAM_CENTER_RADIUS = 65  # mm
BEAM_CENTER_MASKS = None

# Default mask to detector
UNIVERSAL_MASK = None  # 'Mask.XML'
MASKED_PIXELS = '1-8,249-256'

# If it is GPSANS or BIOSANS there could be 2 options to calculate detector efficiencies
MOVING_DETECTORS = True

# THRESHOLD
MIN_THRESHOLD = 0.5
MAX_THRESHOLD = 2.0

# BIO-SANS detector
WING_DETECTOR = False

# END OF USER INPUTS


# Load data files
if INSTRUMENT == 'CG2':
    from drtsans.mono.gpsans.api import prepare_data
    import drtsans.mono.gpsans as mysans
elif INSTRUMENT == 'CG3':
    from drtsans.mono.biosans.api import prepare_data
    import drtsans.mono.biosans as mysans
elif INSTRUMENT == 'EQSANS':
    from drtsans.tof.eqsans.api import prepare_data
    import drtsans.tof.eqsans as mysans
else:
    print('Instrument {} is not supported.  Supported are {}'
          ''.format(INSTRUMENT, 'CG2, EQSANS, CG3'))
    sys.exit(-1)

# Process flood runs
if isinstance(FLOOD_RUNS, int):
    sans_runs = [FLOOD_RUNS]
else:
    sans_runs = list(FLOOD_RUNS)
# Process beam center runs
if DIRECT_BEAM_RUNS is not None:
    beam_center_runs = list(DIRECT_BEAM_RUNS)
else:
    beam_center_runs = None

# Extra masks
extra_mask_dict = dict()
if MASKED_PIXELS is not None:
    extra_mask_dict['Pixel'] = MASKED_PIXELS

# Load data with masking: returning to a list of workspace references
# processing includes: load, mask, normalize by monitor
flood_workspaces = [prepare_data(data='{}_{}'.format(INSTRUMENT, sans_runs[i]),
                                 mask=UNIVERSAL_MASK, btp=extra_mask_dict,
                                 overwrite_instrument=False,
                                 flux_method='monitor') for i in range(len(sans_runs))]
if beam_center_runs is not None:
    beam_center_workspaces = [prepare_data(data='{}_{}'.format(INSTRUMENT, beam_center_runs[i]),
                                           mask=UNIVERSAL_MASK, btp=extra_mask_dict,
                                           overwrite_instrument=False,
                                           flux_method='monitor') for i in range(len(beam_center_runs))]
else:
    beam_center_workspaces = None


# TODO - After testing, moving this to mask_util
def mask_beam_center(flood_ws, beam_center_mask, beam_center_ws, beam_center_radius):
    """Mask beam center

    Mask beam center with 3 algorithms
    1. if beam center mask is present, mask by file
    2. Otherwise if beam center workspace is specified, find beam center from this workspace and mask
    3. Otherwise find beam center for flood workspace and mask itself

    Parameters
    ----------
    flood_ws : ~mantid.api.MatrixWorkspace
        Mantid workspace for flood data
    beam_center_mask : str or None
        Mantid mask file in XML format on beam center
    beam_center_ws : ~mantid.api.MatrixWorkspace or None
        Mantid workspace for beam center
    beam_center_radius : float
        beam center radius in mm

    Returns
    -------

    """
    # Calculate masking (masked file or detectors)
    if beam_center_mask is not None:
        # beam center mask XML file
        masking = beam_center_mask
    else:
        # calculate beam center mask from beam center workspace
        if beam_center_ws is None:
            # if not dedicated beam center workspace, mask itself
            beam_center_ws = flood_ws

        # Use beam center ws to find beam center
        xc, yc = mysans.find_beam_center(beam_center_ws)

        # Center detector to the data workspace (change in geometry)
        mysans.center_detector(flood_ws, xc, yc)

        # Mask the new beam center by 65 mm (Lisa's magic number)
        masking = list(circular_mask_from_beam_center(flood_ws, beam_center_radius))

    # Mask
    apply_mask(flood_ws, mask=masking)  # data_ws reference shall not be invalidated here!

    # Set uncertainties
    # output: masked are zero intensity and zero error
    masked_flood_ws = set_init_uncertainties(flood_ws)

    return masked_flood_ws


# Mask detector center
# Default for GPSANS
if MASK_BEAM_CENTER_RADIUS is None:
    MASK_BEAM_CENTER_RADIUS = 65  # mm
flood_workspaces = mask_beam_center(flood_workspaces, beam_center_runs, BEAM_CENTER_MASKS, MASK_BEAM_CENTER_RADIUS)

# Decide algorithm to prepare sensitivities
if INSTRUMENT in ['CG2', 'CG3'] and MOVING_DETECTORS is True:
    # Prepare by using moving detector algorithm
    from drtsans.mono.gpsans.prepare_sensitivity import prepare_sensitivity_correction

    # Calculate sensitivities for each file
    sens_ws = prepare_sensitivity_correction(flood_workspaces,
                                             threshold_min=MIN_THRESHOLD,
                                             threshold_max=1.5)

else:
    # Prepare by Use the sensitivity patch method
    from drtsans.sensitivity import prepare_sensitivity_correction

    # working on 1 and only 1
    sens_ws = prepare_sensitivity_correction(flood_workspaces[0],
                                             min_threshold=MIN_THRESHOLD,
                                             max_threshold=MAX_THRESHOLD)

# Export
SaveNexusProcessed(InputWorkspace=sens_ws, Filename=SENSITIVITY_FILE)
