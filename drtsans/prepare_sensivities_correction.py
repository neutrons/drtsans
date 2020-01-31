import drtsans
import drtsans.mono.gpsans
import drtsans.mono.biosans
import drtsans.tof.eqsans
from mantid.simpleapi import SaveNexusProcessed, MaskAngle
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
from drtsans.process_uncertainties import set_init_uncertainties

# Constants
CG2 = 'CG2'
CG3 = 'CG3'
EQSANS = 'EQSANS'
PIXEL = 'Pixel'
ANGLE = 'Angle'

PREPARE_DATA = {CG2: drtsans.mono.gpsans.api.prepare_data,
                CG3: drtsans.mono.biosans.api.prepare_data,
                EQSANS: drtsans.tof.eqsans.api.prepare_data}

FIND_BEAM_CENTER = {CG2: drtsans.mono.gpsans.find_beam_center,
                    CG3: drtsans.mono.biosans.find_beam_center,
                    EQSANS: drtsans.tof.eqsans.find_beam_center}

CENTER_DETECTOR = {CG2: drtsans.mono.gpsans.center_detector,
                   CG3: drtsans.mono.biosans.center_detector,
                   EQSANS: drtsans.tof.eqsans.center_detector}


class PrepareSensitivityCorrection(object):
    """Workflow class to prepare sensitivities correction file
    """
    def __init__(self, instrument, is_wing_detector=False):
        """Initialization

        Parameters
        ----------
        instrument : str
            instrument name, CG2, CG2, EQSANS
        is_wing_detector : bool
            Flag to calculate sensivitities for 'wing' detector special to BIOSANS/CG3
        """
        if instrument not in [CG2, CG2, EQSANS]:
            raise RuntimeError('Instrument {} is not supported')
        self._instrument = instrument

        # flood runs
        self._flood_runs = None
        # direct beam (center) runs
        self._direct_beam_center_runs = None
        # mask
        self._default_mask = None
        self._extra_mask_dict = dict()
        self._beam_center_radius = 65  # mm

        # BIOSANS special
        self._is_wing_detector = is_wing_detector
        self._flood_mask_angle = None

    def set_flood_runs(self, flood_runs):
        """Set flood run numbers

        Parameters
        ----------
        flood_runs : ~list or int or tuple
            list of run number as integers

        Returns
        -------
        None

        """
        # Process flood runs
        if isinstance(flood_runs, int):
            self._flood_runs = [flood_runs]
        else:
            self._flood_runs = list(flood_runs)

    def set_direct_beam_runs(self, direct_beam_runs):
        """Set direct beam runs

        Parameters
        ----------
        direct_beam_runs : ~list or int or tuple

        Returns
        -------

        """
        if isinstance(direct_beam_runs, int):
            # defined as a single value
            self._direct_beam_center_runs = [direct_beam_runs]
        else:
            # shall be a list or tuple
            self._direct_beam_center_runs = list(direct_beam_runs)

    def set_masks(self, default_mask, pixels, angles):
        """Set masks

        Parameters
        ----------
        pixels
        angles

        Returns
        -------

        """
        if default_mask is not None:
            self._default_mask = default_mask

        if pixels is not None:
            self._extra_mask_dict[PIXEL] = pixels
        if angles is not None:
            self._extra_mask_dict[ANGLE] = angles

    def set_beam_center_radius(self, radius):
        """Set beam center radius

        Parameters
        ----------
        radius : float
            radius in mm
        Returns
        -------

        """
        self._beam_center_radius = radius

    def execute(self, MOVING_DETECTORS, MIN_THRESHOLD, MAX_THRESHOLD, SENSITIVITY_FILE):

        prepare_data = PREPARE_DATA[self._instrument]

        # Prepare data
        # Load data with masking: returning to a list of workspace references
        # processing includes: load, mask, normalize by monitor
        flood_workspaces = [prepare_data(data='{}_{}'.format(self._instrument, self._flood_runs[i]),
                                         mask=self._default_mask,
                                         btp=self._extra_mask_dict,
                                         overwrite_instrument=False,
                                         flux_method='monitor')
                            for i in range(len(self._flood_runs))]

        # Calculate beam center
        beam_centers = self._calculate_beam_centers()

        # Center detectors
        self._center_detectors(flood_workspaces, beam_centers)

        # Mask beam center
        # BIOSANS - MASK_ANGLES = 1.5, 100.  # 1.5, 57.0   # None
        # self._flood_mask_angle = MASK_ANGLES[1]
        flood_workspaces = [self._mask_beam_center(flood_workspace, beam_centers, mask_angles=self._flood_mask_angle)
                            for flood_workspace in flood_workspaces]

        # Decide algorithm to prepare sensitivities
        if self._instrument in [CG2, CG3] and MOVING_DETECTORS is True:
            # Prepare by using moving detector algorithm
            from drtsans.sensitivity_correction_moving_detectors import calculate_sensitivity_correction

            # Calculate sensitivities for each file
            sens_ws = calculate_sensitivity_correction(flood_workspaces,
                                                       threshold_min=MIN_THRESHOLD,
                                                       threshold_max=1.5)

        else:
            # Prepare by Use the sensitivity patch method
            from drtsans.sensitivity_correction_patch import calculate_sensitivity_correction

            # working on 1 and only 1
            sens_ws = calculate_sensitivity_correction(flood_workspaces[0],
                                                       min_threshold=MIN_THRESHOLD,
                                                       max_threshold=MAX_THRESHOLD)

        # Export
        SaveNexusProcessed(InputWorkspace=sens_ws, Filename=SENSITIVITY_FILE)

    def _calculate_beam_centers(self, mask_angle=None):
        """Find beam centers for all flood runs

        Parameters
        ----------
        mask_angle : float or None
            MASK_ANGLE[0]

        Returns
        -------
        ~list
            List of beam centers (in tuple)

        """
        # Determine the beam center runs
        if self._direct_beam_center_runs is None:
            beam_center_runs = self._flood_runs[:]
        else:
            beam_center_runs = self._direct_beam_center_runs[:]

        # Prepare data
        # TODO FIXME - Only applied for BIOSANS with mask_angle case!!! and GPSANS moving detector
        # TODO FIXME - (cont) It is not nessary for EQSANS because data won't be modified at all!
        prepare_data = PREPARE_DATA[self._instrument]
        beam_center_workspaces = [prepare_data(data='{}_{}'.format(self._instrument, beam_center_runs[i]),
                                               mask=self._default_mask, btp=self._extra_mask_dict,
                                               overwrite_instrument=False,
                                               flux_method='monitor',
                                               output_workspace='BC_{}_{}'.format(self._instrument,
                                                                                  beam_center_runs[i]))
                                  for i in range(len(beam_center_runs))]

        # Find detector center
        find_beam_center = FIND_BEAM_CENTER[self._instrument]
        if self._instrument == CG3 and mask_angle is not None:
            # CG3: apply mask on angle
            for i in range(len(beam_center_workspaces)):
                apply_mask(beam_center_workspaces[i], Components='wing_detector')
            MaskAngle(Workspace=beam_center_runs, MinAngle=mask_angle, Angle="TwoTheta")
        # END-IF

        beam_center_list = [find_beam_center(beam_center_workspaces[i])
                            for i in range(len(beam_center_runs))]

        print('DEBUG: beam centers: {}'.format(beam_center_list))

        return beam_center_list

    def _center_detectors(self, flood_ws_list, beam_center_list):
        """

        Parameters
        ----------
        flood_ws_list
        beam_center_list

        Returns
        -------

        """
        center_detector = CENTER_DETECTOR[self._instrument]

        if self._instrument == CG3:
            # BIO SANS : 3 center value
            for i in range(len(flood_ws_list)):
                xc, yc, wc = beam_center_list[i]
                center_detector(flood_ws_list[i], xc, yc, wc)
        else:
            # GPSANS, EQSANS
            for i in range(len(flood_ws_list)):
                xc, yc = beam_center_list[i]
                center_detector(flood_ws_list[i], xc, yc)

    def _mask_beam_center(self, flood_ws, beam_center, mask_angles=None):
        """Mask beam center

        Mask beam center with 3 algorithms
        1. if beam center mask is present, mask by file
        2. Otherwise if beam center workspace is specified, find beam center from this workspace and mask
        3. Otherwise find beam center for flood workspace and mask itself

        Parameters
        ----------
        flood_ws : ~mantid.api.MatrixWorkspace
            Mantid workspace for flood data
        beam_center : tuple or str
            if tuple, beam centers (xc, yc, wc) / (xc, yc); str: beam center masks file
        Returns
        -------

        """
        # Calculate masking (masked file or detectors)
        if isinstance(beam_center, str):
            # beam center mask XML file: apply mask
            apply_mask(flood_ws, mask=beam_center)  # data_ws reference shall not be invalidated here!
        elif mask_angles is not None and self._instrument == CG3:
            # Mask angle
            # Mask wing detector right top/bottom corners
            if self._is_wing_detector is False:
                # main detector: mask wing
                component = 'wing_detector'
            else:
                # wing detector: mask main
                component = 'detector1'
            apply_mask(flood_ws, Components=component)
            MaskAngle(Workspace=flood_ws, MaxAngle=mask_angles[1], Angle="TwoTheta")
        else:
            # calculate beam center mask from beam center workspace
            # Mask the new beam center by 65 mm (Lisa's magic number)
            masking = list(circular_mask_from_beam_center(flood_ws, self._beam_center_radius))
            # Mask
            apply_mask(flood_ws, mask=masking)  # data_ws reference shall not be invalidated here!

        # Set uncertainties
        # output: masked are zero intensity and zero error
        masked_flood_ws = set_init_uncertainties(flood_ws)

        return masked_flood_ws


def debug_tools():
    # Export sensitivities calculated in file for quick review
    # Export 2D view
    png_name = output_file.split('.')[0] + '.png'
    export_detector_view(sensitivities, png_name)
    # Export to hdf5 for a quick review
    # sens_h5 = h5py.File('{}.h5'.format(output_file.split('.')[0]), 'w')
    sens_group = sens_h5.create_group('Sensitivities')
    sens_group.create_dataset('sensitivities', data=sensitivities)
    sens_group.create_dataset('sensitivities error', data=sensitivities_error)
    sens_h5.close()


def export_detector_view(ws, png_name):
    """Export detector view to a PNG file

    This method is for debugging purpose

    Parameters
    ----------
    ws : ~mantid.api.MatrixWorkspace
        Workspace to plot in detector view
    png_name : str
        Path of the output PNG file

    Returns
    -------
    None

    """
    if isinstance(ws, np.ndarray):
        vec_y = ws
    else:
        vec_y = ws.extractY().transpose()
    vec_y = vec_y.reshape((192, 256)).transpose()

    plt.imshow(vec_y)

    plt.savefig(png_name)

    return
