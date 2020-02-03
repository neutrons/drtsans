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
        # Prepare data
        # get right prepare_data method specified to instrument type
        prepare_data = PREPARE_DATA[self._instrument]

        # Load data with masking: returning to a list of workspace references
        # processing includes: load, mask, normalize by monitor
        flood_workspaces = [prepare_data(data='{}_{}'.format(self._instrument, self._flood_runs[i]),
                                         mask=self._default_mask,
                                         btp=self._extra_mask_dict,
                                         overwrite_instrument=False,
                                         flux_method='monitor')
                            for i in range(len(self._flood_runs))]

        # DEBUG OUTPUT
        for i in range(len(flood_workspaces)):
            workspace = flood_workspaces[i]
            output_file = 'Step_1_{}.nxs'.format(str(workspace))
            debug_output(workspace, output_file, note='Step 1 {}'.format(str(workspace)))

        # Calculate beam centers
        beam_centers = [self._calculate_beam_center(i, flood_workspaces[i]) for i in range(len(flood_workspaces))]

        # DEBUG OUTPUT
        for i in range(len(beam_centers)):
            print('Index {}  Beam Center = {}'.format(i, beam_centers[i]))

        # Center detectors
        self._center_detectors(flood_workspaces, beam_centers)

        # DEBUG OUTPUT
        for i in range(len(flood_workspaces)):
            workspace = flood_workspaces[i]
            output_file = 'Step_2_Centred_{}.nxs'.format(str(workspace))
            debug_output(workspace, output_file, note='Step 2 Detector Centered {}'.format(str(workspace)))

        # Mask beam center
        # BIOSANS - MASK_ANGLES = 1.5, 100.  # 1.5, 57.0   # None
        # self._flood_mask_angle = MASK_ANGLES[1]
        flood_workspaces = [self._mask_beam_center(flood_workspace, beam_centers)
                            for flood_workspace in flood_workspaces]

        # DEBUG OUTPUT
        for i in range(len(flood_workspaces)):
            workspace = flood_workspaces[i]
            output_file = 'Step_3_Masked_{}.nxs'.format(str(workspace))
            debug_output(workspace, output_file, note='Step 3 Center Masked {}'.format(str(workspace)))

        # Solid angle correction
        flood_workspaces = [self._apply_solid_angle_correction(workspace)
                            for workspace in flood_workspaces]

        # Transmission correction
        if self._transmission_runs is not None and self._is_wing_detector is False:
            # calculate transmission corrections
            trans_corr_ws_list = [self._calculate_transmission_correction(
                transmission_run=self._transmission_runs[i],
                transmission_flood_run=self._transmission_flood_runs[i])
                for i in range(len(self._transmission_runs))]
            # apply
            flood_workspaces = [self._apply_transmission_correction(flood_ws=flood_workspaces[i],
                                                                    transmission_corr_ws=trans_corr_ws_list[i],
                                                                    is_theta_dep_corr=self._theta_dep_correction)
                                for i in range(len(flood_workspaces))]

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

            # working on 1 and only 1
            sens_ws = calculate_sensitivity_correction(flood_workspaces[0],
                                                       min_threshold=min_threshold,
                                                       max_threshold=max_threshold)

        # Export
        SaveNexusProcessed(InputWorkspace=sens_ws, Filename=output_nexus_name)

    def _calculate_beam_center(self, index, flood_ws):
        """Find beam centers for all flood runs

        Parameters
        ----------
        flood_ws : ~MatrixWorkspace
            flood workspaces that might be used for beam center runs

        Returns
        -------
        ~tuple
            beam center as xc, yc and possible wc for BIOSANS

        """
        # Determine the beam center runs
        beam_center_workspace = None

        if self._direct_beam_center_runs is None and self._instrument == CG3 \
                and self._wing_det_mask_angle is not None:
            # CG3, flood run as direct beam center and mask angle is defined
            # In this case, run shall be loaded to another workspace for masking
            beam_center_run = self._flood_runs[index]
        elif self._direct_beam_center_runs is None:
            # Use flood run's workspace to find beam center without changing it
            beam_center_workspace = flood_ws
            beam_center_run = None
        else:
            # Direct beam run is specified
            beam_center_run = self._direct_beam_center_runs[index]

        # Prepare data
        # Only applied for BIOSANS with mask_angle case!!! and GPSANS moving detector
        # It is not necessary for EQSANS because data won't be modified at all!
        if beam_center_workspace is None:
            prepare_data = PREPARE_DATA[self._instrument]
            beam_center_workspace = prepare_data(data='{}_{}'.format(self._instrument, beam_center_run),
                                                 mask=self._default_mask, btp=self._extra_mask_dict,
                                                 overwrite_instrument=False,
                                                 flux_method='monitor',
                                                 output_workspace='BC_{}_{}'.format(self._instrument,
                                                                                    beam_center_run))

            # Mask angle for CG3: apply mask on angle
            if self._instrument == CG3 and self._wing_det_mask_angle is not None:
                # mask wing detector
                apply_mask(beam_center_workspace, Components='wing_detector')
                # mask 2-theta angle on main detector
                MaskAngle(Workspace=beam_center_workspace, MinAngle=self._wing_det_mask_angle, Angle="TwoTheta")
            # END-IF
        # END-IF

        # Find detector center
        find_beam_center = FIND_BEAM_CENTER[self._instrument]
        beam_center = find_beam_center(beam_center_workspace)

        print('DEBUG: {} beam centers: {}'.format(index, beam_center))

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

    def _calculate_transmission_correction(self, transmission_run, transmission_flood_run):
        """

        Returns
        -------
        Workspace
            processed transmission workspace

        """
        prepare_data = PREPARE_DATA[self._instrument]

        # Load, mask default and pixels, and normalize
        transmission_workspace = prepare_data(data='{}_{}'.format(self._instrument, transmission_run),
                                              mask=self._default_mask, btp=self._extra_mask_dict,
                                              overwrite_instrument=False,
                                              flux_method='time',
                                              output_workspace='TM_{}_{}'.format(self._instrument,
                                                                                 transmission_run))
        # Apply mask
        if self._instrument == CG3:
            apply_mask(transmission_workspace, Components='wing_detector')
            MaskAngle(Workspace=transmission_workspace, MinAngle=2 * self._main_det_mask_angle, Angle="TwoTheta")

        # Load, mask default and pixels, normalize
        transmission_flood_ws = prepare_data(data='{}_{}'.format(self._instrument, transmission_flood_run),
                                             mask=self._default_mask, btp=self._extra_mask_dict,
                                             overwrite_instrument=False,
                                             flux_method='time',
                                             output_workspace='TM_{}_{}'.format(self._instrument,
                                                                                transmission_flood_run))
        # Apply mask
        if self._instrument == CG3:
            apply_mask(transmission_flood_ws, Components='wing_detector')
            MaskAngle(Workspace=transmission_flood_ws, MinAngle=2 * self._main_det_mask_angle, Angle="TwoTheta")

        # Zero-Angle Transmission Co-efficients
        calculate_transmission = CALCULATE_TRANSMISSION[self._instrument]
        ws_tr = calculate_transmission(transmission_flood_ws, transmission_workspace)
        average_zero_angle = np.mean(ws_tr.readY(0))
        average_zero_angle_error = np.linalg.norm(ws_tr.readE(0))
        print("\tTransmission Coefficient is....{:.3f} +/- {:.3f}"
              "".format(average_zero_angle, average_zero_angle_error))

        return ws_tr

    def _apply_solid_angle_correction(self, flood_ws):
        """Apply solid angle correction

        Parameters
        ----------
        flood_ws

        Returns
        -------
        MatrixWorkspace

        """
        # Transmission correction is only acted on wing detector
        if self._is_wing_detector:
            return flood_ws
        # Flag is off
        if self._solid_angle_correction is False:
            return flood_ws

        # apply solid angle correction
        solid_angle_correction = SOLID_ANGLE_CORRECTION[self._instrument]
        flood_ws = solid_angle_correction(flood_ws)

        return flood_ws

    def _apply_transmission_correction(self, flood_ws, transmission_corr_ws, is_theta_dep_corr):
        """Apply transmission correction

        Parameters
        ----------
        flood_ws
        transmission_corr_ws : workspace
            transmission correct workspace

        Returns
        -------

        """
        apply_transmission_correction = APPLY_TRANSMISSION[self._instrument]

        flood_ws = apply_transmission_correction(flood_ws, trans_workspace=transmission_corr_ws,
                                                 theta_dependent=is_theta_dep_corr)

        return flood_ws


def debug_output(workspace, output_file, note=''):
    """

    Parameters
    ----------
    data : numpy.ndarray
        data to plot
    output_file : str
        output file name as reference

    Returns
    -------

    """
    import h5py
    from mantid.simpleapi import GeneratePythonScript

    # Save Nexus
    script_text = GeneratePythonScript(workspace)
    print(note)
    print(script_text.strip())
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
    import numpy as np
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
