import numpy as np
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

PREPARE_DATA = {CG2: drtsans.mono.gpsans.api.prepare_data,
                CG3: drtsans.mono.biosans.api.prepare_data,
                EQSANS: drtsans.tof.eqsans.api.prepare_data}

FIND_BEAM_CENTER = {CG2: drtsans.mono.gpsans.find_beam_center,
                    CG3: drtsans.mono.biosans.find_beam_center,
                    EQSANS: drtsans.tof.eqsans.find_beam_center}

CENTER_DETECTOR = {CG2: drtsans.mono.gpsans.center_detector,
                   CG3: drtsans.mono.biosans.center_detector,
                   EQSANS: drtsans.tof.eqsans.center_detector}

CALCULATE_TRANSMISSION = {CG2: drtsans.mono.gpsans.calculate_transmission,
                          CG3: drtsans.mono.biosans.calculate_transmission,
                          EQSANS: drtsans.tof.eqsans.calculate_transmission}

APPLY_TRANSMISSION = {CG2: drtsans.mono.gpsans.apply_transmission_correction,
                      CG3: drtsans.mono.biosans.apply_transmission_correction,
                      EQSANS: drtsans.tof.eqsans.apply_transmission_correction}

SOLID_ANGLE_CORRECTION = {
    CG2: drtsans.mono.gpsans.solid_angle_correction,
    CG3: drtsans.mono.biosans.solid_angle_correction,
}


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
        if instrument not in [CG2, CG3, EQSANS]:
            raise RuntimeError('Instrument {} is not supported'.format(instrument))
        self._instrument = instrument

        # flood runs
        self._flood_runs = None
        # direct beam (center) runs
        self._direct_beam_center_runs = None
        # mask
        self._default_mask = None
        self._extra_mask_dict = dict()
        self._beam_center_radius = 65  # mm

        # Transmission correction (BIOSANS)
        self._transmission_runs = None
        self._transmission_flood_runs = None
        self._theta_dep_correction = False

        # Dark current
        self._dark_current_runs = None

        # Apply solid angle correction or not?
        self._solid_angle_correction = False

        # BIOSANS special
        if self._instrument == CG3:
            self._is_wing_detector = is_wing_detector
        else:
            self._is_wing_detector = False

        # Mask angle for wing detector
        self._wing_det_mask_angle = None
        # Mask angle on main detector
        self._main_det_mask_angle = None

    def set_solid_angle_correction_flag(self, apply_correction):
        """Set the flag to apply solid angle correction

        Parameters
        ----------
        apply_correction : bool
            Flag for applying

        Returns
        -------

        """
        self._solid_angle_correction = apply_correction

        if self._is_wing_detector:
            # solid angle correction on BIOSANS Wing detector cannot be trusted
            self._solid_angle_correction = False

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

    def set_dark_current_runs(self, dark_current_runs):
        """Set dark current runs

        Parameters
        ----------
        dark_current_runs : ~list or ~tuple or int
            Dark current run(s)'s run number(s)

        Returns
        -------
        None

        """
        if dark_current_runs is None:
            self._dark_current_runs = None
        elif isinstance(dark_current_runs, int):
            self._dark_current_runs = [dark_current_runs]
        else:
            self._dark_current_runs = list(dark_current_runs)

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

    def set_masks(self, default_mask, pixels, wing_det_mask_angle=None, main_det_mask_angle=None):
        """Set masks

        Parameters
        ----------
        default_mask : str or None
            mask file name
        pixels : str or None
            pixels to mask.  Example: '1-8,249-256'
        wing_det_mask_angle : float or None
            angle to mask (MaskAngle) to (BIOSANS) wing detector
        main_det_mask_angle : float or None
            angle to mask (MaskAngle) to main detector

        Returns
        -------
        None

        """
        # default mask (XML) file name
        if default_mask is not None:
            self._default_mask = default_mask

        # pixels to mask
        if pixels is not None:
            self._extra_mask_dict[PIXEL] = pixels

        # angles to mask (BIOSANS)
        if wing_det_mask_angle is not None:
            self._wing_det_mask_angle = wing_det_mask_angle
        if main_det_mask_angle is not None:
            self._main_det_mask_angle = main_det_mask_angle

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

    def set_transmission_correction(self, transmission_flood_runs, transmission_beam_run):
        """Set transmission beam run and transmission flood runs

        Parameters
        ----------
        transmission_flood_runs : ~list

        transmission_beam_run : ~list


        Returns
        -------

        """
        if isinstance(transmission_beam_run, int):
            self._transmission_runs = [transmission_beam_run]
        else:
            self._transmission_runs = list(transmission_beam_run)

        if isinstance(transmission_flood_runs, int):
            self._transmission_flood_runs = [transmission_flood_runs]
        else:
            self._transmission_flood_runs = list(transmission_flood_runs)

    def set_theta_dependent_correction_flag(self, flag):
        """Set the flag to do theta dep with transmission correction

        Parameters
        ----------
        flag : bool
            True to do the correction

        Returns
        -------

        """
        self._theta_dep_correction = flag

    def _prepare_flood_data(self, index, beam_center, dark_current_run):
        """Prepare flood data including
        (1) load
        (2) mask: default, pixels
        (3) center detector
        (4) optionally solid angle correction
        (5) optionally dark current correction
        (6) normalization

        Parameters
        ----------
        index
        beam_center
        dark_current_run

        Returns
        -------

        """

        # Prepare data
        # get right prepare_data method specified to instrument type
        prepare_data = PREPARE_DATA[self._instrument]

        bio_specials = dict()
        if self._instrument == CG3:
            bio_specials['center_y_wing'] = beam_center[2]

        print('DARK CURRENT RUN : {}.    type: {}'.format(dark_current_run, type(dark_current_run)))

        # Load data with masking: returning to a list of workspace references
        # processing includes: load, mask, normalize by monitor
        flood_ws = prepare_data(data='{}_{}'.format(self._instrument, self._flood_runs[index]),
                                mask=self._default_mask,
                                btp=self._extra_mask_dict,
                                center_x=beam_center[0],
                                center_y=beam_center[1],
                                dark_current=dark_current_run,
                                overwrite_instrument=False,
                                flux_method='monitor',
                                solid_angle=self._solid_angle_correction,
                                **bio_specials)

        return flood_ws

    @staticmethod
    def _get_masked_detectors(workspace):
        """Get the detector masking information

        Parameters
        ----------
        workspace

        Returns
        -------
        numpy.ndarray
            (N, 1) bool array, True for being masked

        """
        masked_array = workspace.extractE() < 1E-5

        masked_bis = np.where(masked_array)[0]

        print('IMPORTANT ....................  Masked workspace indexes: {},'
              'shape = {},  Original masked array: dtype = {}, shape = {}'
              ''.format(masked_bis, masked_bis.shape, masked_array.dtype, masked_array.shape))

        return masked_array

    @staticmethod
    def _set_mask_value(flood_workspace, use_moving_detector_method, det_mask_array):
        """Set masked pixels' values to NaN or -infinity

        Parameters
        ----------
        flood_workspace
        use_moving_detector_method
        det_mask_array

        Returns
        -------

        """
        # Complete mask array
        total_mask_array = flood_workspace.extractE() < 1E-6

        # Sanity check for detector-patch case
        if not use_moving_detector_method:
            assert total_mask_array.shape == det_mask_array.shape, '{} <> {}'.format(total_mask_array.shape,
                                                                                     det_mask_array.shape)
            assert total_mask_array.dtype == det_mask_array.dtype, 'dtype wrong'

        import time

        time_start = time.time()

        num_spec = flood_workspace.getNumberHistograms()
        for i in range(num_spec):
            if total_mask_array[i] and use_moving_detector_method:
                # only total.. from center mask.  moving detector: set to NaN
                flood_workspace.dataY(i)[0] = np.nan
                flood_workspace.dataE(i)[0] = np.nan
            elif total_mask_array[i] and not use_moving_detector_method and det_mask_array[i]:
                # both masked: pixel and default mask for detector-patch algorithm: set to NaN
                flood_workspace.dataY(i)[0] = np.nan
                flood_workspace.dataE(i)[0] = np.nan
            elif total_mask_array[i]:
                # only total mask.  so it is from beam center. detector-patch algorithm: set to -inf
                flood_workspace.dataY(i)[0] = -np.NINF
                flood_workspace.dataE(i)[0] = -np.NINF
            elif not total_mask_array[i] and not use_moving_detector_method and det_mask_array[i]:
                # logic error: mask missing
                raise RuntimeError('Mask missing')
        # END-FOR

        time_stop = time.time()
        print('Set NaN/Infinity: time = {}'.format(time_stop - time_start))
        print('Number of infinities = {}'.format(len(np.where(np.isinf(flood_workspace.extractY()))[0])))
        print('Number of NaNs       = {}'.format(len(np.where(np.isnan(flood_workspace.extractY()))[0])))

        return flood_workspace

    def execute(self, use_moving_detector_method, min_threshold, max_threshold, output_nexus_name):
        """Main workflow method to calculate sensitivities correction

        Parameters
        ----------
        use_moving_detector_method : bool
            Flag to use 'moving detectors' method; Otherwise, use 'detector patch' method
        min_threshold : float
            minimum threshold of normalized count for GOOD pixels
        max_threshold : float
            maximum threshold of normalized count for GOOD pixels
        output_nexus_name : str
            path to the output processed NeXus file

        Returns
        -------
        None

        """
        # Number of pair of workspaces to process
        num_workspaces_set = len(self._flood_runs)

        # Load beam center runs and calculate beam centers
        beam_centers = [self._calculate_beam_center(i)
                        for i in range(num_workspaces_set)]

        # Set default value to dark current runs
        if self._dark_current_runs is None:
            self._dark_current_runs = [None] * num_workspaces_set

        # Load and process flood data with (1) mask (2) center detector and (3) solid angle correction
        flood_workspaces = [self._prepare_flood_data(i, beam_centers[i], self._dark_current_runs[i])
                            for i in range(num_workspaces_set)]

        # Retrieve masked detectors
        if not use_moving_detector_method:
            det_mask_list = [self._get_masked_detectors(flood_workspaces[i])
                             for i in range(num_workspaces_set)]
        else:
            det_mask_list = [None] * num_workspaces_set

        # Mask beam centers
        flood_workspaces = [self._mask_beam_center(flood_workspaces[i], beam_centers[i])
                            for i in range(num_workspaces_set)]

        # Transmission correction as an option
        if self._transmission_runs is not None and not self._is_wing_detector:
            # Must have transmission run specified and cannot be wing detector (of CG3)
            flood_workspaces = [self._apply_transmission_correction(flood_workspaces[i],
                                                                    self._transmission_runs[i],
                                                                    self._transmission_flood_runs[i])
                                for i in range(num_workspaces_set)]

        # Set the masked pixels' counts to nan and -infinity
        flood_workspaces = [self._set_mask_value(flood_workspaces[i], use_moving_detector_method,
                                                 det_mask_list[i])
                            for i in range(num_workspaces_set)]

        print('Preparation of data is over....')
        print('Number of infinities = {}'.format(len(np.where(np.isinf(flood_workspaces[0].extractY()))[0])))
        print('Number of NaNs       = {}'.format(len(np.where(np.isnan(flood_workspaces[0].extractY()))[0])))

        # DEBUG OUTPUT
        for i in range(len(flood_workspaces)):
            workspace = flood_workspaces[i]
            output_file = 'Step_JustBefore_{}.nxs'.format(str(workspace))
            debug_output(workspace, output_file, note='Step 3 Center Masked {}'.format(str(workspace)))

        # Decide algorithm to prepare sensitivities
        if self._instrument in [CG2, CG3] and use_moving_detector_method is True:
            # Prepare by using moving detector algorithm
            from drtsans.sensitivity_correction_moving_detectors import calculate_sensitivity_correction

            # Calculate sensitivities for each file
            sens_ws = calculate_sensitivity_correction(flood_workspaces,
                                                       threshold_min=min_threshold,
                                                       threshold_max=1.5)

        else:
            # Prepare by Use the sensitivity patch method
            from drtsans.sensitivity_correction_patch import calculate_sensitivity_correction

            if self._instrument == CG3:
                polynomial_order = 3
            else:
                polynomial_order = 2

            # working on 1 and only 1
            sens_ws = calculate_sensitivity_correction(flood_workspaces[0],
                                                       min_threshold=min_threshold,
                                                       max_threshold=max_threshold,
                                                       poly_order=polynomial_order,
                                                       min_detectors_per_tube=50)

        # Export
        SaveNexusProcessed(InputWorkspace=sens_ws, Filename=output_nexus_name)

    def _calculate_beam_center(self, index):
        """Find beam centers for all flood runs

        Beam center run shall be
        (1) masked properly (default mask + top/bottom)
        (2) NOT corrected by solid angle

        Parameters
        ----------
        index : int
            beam center run index mapped to flood run

        Returns
        -------
        ~tuple
            beam center as xc, yc and possible wc for BIOSANS

        """
        if self._direct_beam_center_runs is None and self._instrument == CG3 \
                and self._wing_det_mask_angle is not None:
            # CG3, flood run as direct beam center and mask angle is defined
            # In this case, run shall be loaded to another workspace for masking
            beam_center_run = self._flood_runs[index]
        elif self._direct_beam_center_runs is None:
            raise RuntimeError('Beam center runs must be given for {}'.format(self._instrument))
        else:
            # Direct beam run is specified
            beam_center_run = self._direct_beam_center_runs[index]

        # Prepare data
        # Only applied for BIOSANS with mask_angle case!!! and GPSANS moving detector
        # It is not necessary for EQSANS because data won't be modified at all!
        prepare_data = PREPARE_DATA[self._instrument]
        beam_center_workspace = prepare_data(data='{}_{}'.format(self._instrument, beam_center_run),
                                             mask=self._default_mask,
                                             btp=self._extra_mask_dict,
                                             overwrite_instrument=False,
                                             flux_method='monitor',
                                             solid_angle=False,
                                             output_workspace='BC_{}_{}'.format(self._instrument,
                                                                                beam_center_run))
        # Mask angle for CG3: apply mask on angle
        if self._instrument == CG3 and self._wing_det_mask_angle is not None:
            # mask wing detector
            apply_mask(beam_center_workspace, Components='wing_detector')
            # mask 2-theta angle on main detector
            MaskAngle(Workspace=beam_center_workspace, MinAngle=self._wing_det_mask_angle, Angle="TwoTheta")
        # END-IF

        # Find detector center
        find_beam_center = FIND_BEAM_CENTER[self._instrument]
        beam_center = find_beam_center(beam_center_workspace)

        print('DEBUG: {} beam centers: {}'.format(index, beam_center))
        debug_output(beam_center_workspace, output_file='BeamCenter_{}.nxs'.format(beam_center_run),
                     note='{} beam centers: {}'.format(index, beam_center))

        return beam_center

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

    def _mask_beam_center(self, flood_ws, beam_center):
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
        elif self._main_det_mask_angle is not None and self._instrument == CG3:
            # CG3 special: Mask 2-theta angle
            # Mask wing detector right top/bottom corners
            if self._is_wing_detector is False:
                # main detector: mask wing
                component = 'wing_detector'
            else:
                # wing detector: mask main
                component = 'detector1'
            apply_mask(flood_ws, Components=component)
            # mask 2theta
            MaskAngle(Workspace=flood_ws, MaxAngle=self._main_det_mask_angle, Angle="TwoTheta")
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

    def _apply_transmission_correction(self, flood_ws, transmission_beam_run, transmission_flood_run):
        """Calculate and pply transmission correction

        Parameters
        ----------
        flood_ws : MarixWorkspace
            Flood run workspace to
            transmission correct workspace
        transmission_beam_run : int
            run number for transmission beam run
        transmission_flood_run : int
            run number for transmission flood run

        Returns
        -------
        MatrixWorkspace
            Flood workspace with transmission corrected

        """
        prepare_data = PREPARE_DATA[self._instrument]

        # Load, mask default and pixels, and normalize
        transmission_workspace = prepare_data(data='{}_{}'.format(self._instrument, transmission_beam_run),
                                              mask=self._default_mask, btp=self._extra_mask_dict,
                                              overwrite_instrument=False,
                                              flux_method='time',
                                              solid_angle=False,
                                              output_workspace='TRANS_{}_{}'.format(self._instrument,
                                                                                    transmission_beam_run))
        # Apply mask
        if self._instrument == CG3:
            apply_mask(transmission_workspace, Components='wing_detector')
            MaskAngle(Workspace=transmission_workspace, MinAngle=2 * self._main_det_mask_angle, Angle="TwoTheta")

        # Load, mask default and pixels, normalize
        transmission_flood_ws = prepare_data(data='{}_{}'.format(self._instrument, transmission_flood_run),
                                             mask=self._default_mask, btp=self._extra_mask_dict,
                                             overwrite_instrument=False,
                                             flux_method='time',
                                             solid_angle=self._solid_angle_correction,
                                             output_workspace='TRANS_{}_{}'.format(self._instrument,
                                                                                   transmission_flood_run))
        # Apply mask
        if self._instrument == CG3:
            apply_mask(transmission_flood_ws, Components='wing_detector')
            MaskAngle(Workspace=transmission_flood_ws, MinAngle=2 * self._main_det_mask_angle, Angle="TwoTheta")

        # Zero-Angle Transmission Co-efficients
        calculate_transmission = CALCULATE_TRANSMISSION[self._instrument]
        transmission_corr_ws = calculate_transmission(transmission_flood_ws, transmission_workspace)
        average_zero_angle = np.mean(transmission_corr_ws.readY(0))
        average_zero_angle_error = np.linalg.norm(transmission_corr_ws.readE(0))
        print("\tTransmission Coefficient is....{:.3f} +/- {:.3f}"
              "".format(average_zero_angle, average_zero_angle_error))

        # Apply calculated transmission
        apply_transmission_correction = APPLY_TRANSMISSION[self._instrument]
        flood_ws = apply_transmission_correction(flood_ws, trans_workspace=transmission_corr_ws,
                                                 theta_dependent=self._theta_dep_correction)

        return flood_ws


def debug_output(workspace, output_file, note=''):
    """

    Parameters
    ----------
    workspace : numpy.ndarray
        data to plot
    output_file : str
        output file name as reference

    Returns
    -------

    """
    import h5py
    # from mantid.simpleapi import GeneratePythonScript

    # Save Nexus
    # script_text = GeneratePythonScript(workspace)
    print(note)
    # print(script_text.strip())
    SaveNexusProcessed(InputWorkspace=workspace,
                       Filename=output_file)

    data = workspace.extractY()
    data_error = workspace.extractE()

    # Export sensitivities calculated in file for quick review
    # Export 2D view
    png_name = output_file.split('.')[0] + '.png'
    export_detector_view(data, png_name)
    # Export to hdf5 for a quick review
    sens_h5 = h5py.File('{}.h5'.format(output_file.split('.')[0]), 'w')
    sens_group = sens_h5.create_group('Data')
    sens_group.create_dataset('Data', data=data)
    if data_error is not None:
        sens_group.create_dataset('Data error', data=data_error)
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
    from matplotlib import pyplot as plt
    if isinstance(ws, np.ndarray):
        vec_y = ws
    else:
        vec_y = ws.extractY().transpose()
    try:
        vec_y = vec_y.reshape((192, 256)).transpose()
    except ValueError:
        # exception in case not GPSANS. return
        return

    plt.imshow(vec_y)

    plt.savefig(png_name)

    return
