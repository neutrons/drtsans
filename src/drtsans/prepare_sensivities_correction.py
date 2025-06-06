# standard imports
import os
from typing import Optional

# third party imports
import h5py
from mantid.api import mtd
from mantid.kernel import logger
from mantid.simpleapi import (
    SaveNexusProcessed,
    Integration,
    MaskDetectors,
    CreateWorkspace,
)
import numpy as np

# drtsans imports
from drtsans.load import load_events
from drtsans.geometry import panel_names
from drtsans.instruments import extract_run_number
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
import drtsans.mono.gpsans
import drtsans.mono.biosans
from drtsans.process_uncertainties import set_init_uncertainties
from drtsans.sensitivity_correction_moving_detectors import (
    calculate_sensitivity_correction as calculate_sensitivity_correction_moving,
)
from drtsans.sensitivity_correction_patch import (
    calculate_sensitivity_correction as calculate_sensitivity_correction_patch,
)
import drtsans.tof.eqsans


# Constants
CG2 = "CG2"
CG3 = "CG3"
EQSANS = "EQSANS"
PIXEL = "Pixel"
MOVING_DETECTORS = "Moving Detectors"
PATCHING_DETECTORS = "Patching Detectors"

# As this script is a wrapper to handle script prepare_sensitivity
# (e.g. https://github.com/neutrons/drtsans/blob/next/scripts%2Fprepare_sensitivities_biosans.py),
# by which user just needs to set up instrument name but not is required to import the right modules
# (for example drtsans.mono.gpsans or drtsans.tof.eqsans),
# therefore it has to import the correct ones according instrument name in string.
# Using dictionary with instrument name as key is solution for it.

# prepare data  in .../api.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fmono%2Fgpsans%2Fapi.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fmono%2Fbiosans%2Fapi.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Ftof%2Feqsans%2Fapi.py
PREPARE_DATA = {
    CG2: drtsans.mono.gpsans.api.prepare_data,
    CG3: drtsans.mono.biosans.api.prepare_data,
    EQSANS: drtsans.tof.eqsans.api.prepare_data,
}

# Find beam center in .../find_beam_center.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fmono%2Fbiosans%2Fbeam_finder.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fbeam_finder.py
FIND_BEAM_CENTER = {
    CG2: drtsans.mono.gpsans.find_beam_center,
    CG3: drtsans.mono.biosans.find_beam_center,
    EQSANS: drtsans.tof.eqsans.find_beam_center,
}

# Center detector in
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fmono%2Fbiosans%2Fbeam_finder.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fbeam_finder.py
CENTER_DETECTOR = {
    CG2: drtsans.mono.gpsans.center_detector,
    CG3: drtsans.mono.biosans.center_detector,
    EQSANS: drtsans.tof.eqsans.center_detector,
}

# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fsolid_angle.py
SOLID_ANGLE_CORRECTION = {
    CG2: drtsans.mono.gpsans.solid_angle_correction,
    CG3: drtsans.mono.biosans.solid_angle_correction,
    EQSANS: drtsans.solid_angle_correction,
}

# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fsensitivity_correction_moving_detectors.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans%2Fsensitivity_correction_patch.py
CALCULATE_SENSITIVITY_CORRECTION = {
    MOVING_DETECTORS: calculate_sensitivity_correction_moving,
    PATCHING_DETECTORS: calculate_sensitivity_correction_patch,
}


class PrepareSensitivityCorrection(object):
    """Workflow (container) class to prepare sensitivities correction file

    It tries to manage the various configuration and approaches for instrument scientists to prepare
    sensitivities file for EQSANS, GPSANS and BIOSANS.

    """

    @staticmethod
    def _get_masked_detectors(workspace):
        """Get the detector masking information

        Parameters
        ----------
        workspace : ~mantid.api.MatrixWorkspace
            Workspace to get masked detectors' masking status

        Returns
        -------
        numpy.ndarray
            (N, 1) bool array, True for being masked

        """
        # The masked pixels after `set_uncertainties()` will have zero uncertainty.
        # Thus, it is an efficient way to identify them by check uncertainties (E) close to zero
        masked_array = workspace.extractE() < 1e-5
        return masked_array

    @staticmethod
    def _set_mask_value(flood_workspace, det_mask_array, use_moving_detector_method=True):
        """Set masked pixels' values to NaN or -infinity according to mask type and sensitivity correction
        algorithm

        Parameters
        ----------
        flood_workspace :
            Flood data workspace space to have masked pixels' value set
        det_mask_array : numpy.ndarray
            Array to indicate pixel to be masked or not
        use_moving_detector_method : bool
            True for calculating sensitivities by moving detector algorithm;
            otherwise for detector patching algorithm

        Returns
        -------

        """
        # Complete mask array.  Flood workspace has been processed by set_uncertainties.  Therefore all the masked
        # pixels' uncertainties are zero, which is different from other pixels
        total_mask_array = flood_workspace.extractE() < 1e-6

        # Loop through each detector pixel to check its masking state to determine whether its value shall be
        # set to NaN, -infinity or not changed (i.e., for pixels without mask)
        num_spec = flood_workspace.getNumberHistograms()
        problematic_pixels = list()
        for i in range(num_spec):
            if total_mask_array[i][0] and use_moving_detector_method:
                # Moving detector algorithm.  Any masked detector pixel is set to NaN
                flood_workspace.dataY(i)[0] = np.nan
                flood_workspace.dataE(i)[0] = np.nan
            elif total_mask_array[i][0] and not use_moving_detector_method and det_mask_array[i][0]:
                # Patch detector method: Masked as the bad pixels and thus set to NaN
                flood_workspace.dataY(i)[0] = np.nan
                flood_workspace.dataE(i)[0] = np.nan
            elif total_mask_array[i][0]:
                # Patch detector method: Pixels that have not been masked as bad pixels, but have been
                # identified as needing to have values set by the patch applied. To identify them, the
                # value is set to -INF.
                flood_workspace.dataY(i)[0] = -np.inf
                flood_workspace.dataE(i)[0] = -np.inf
            elif not total_mask_array[i][0] and not use_moving_detector_method and det_mask_array[i][0]:
                # Logic error: impossible case
                problematic_pixels.append(i)
        # END-FOR

        # Array
        if len(problematic_pixels) > 0:
            raise RuntimeError(
                f"Impossible case: pixels {problematic_pixels} has local detector mask is on, but total mask is off"
            )

        logger.debug(
            "Patch detector method: Pixels that have not been masked as bad pixels, but have been"
            "identified as needing to have values set by the patch applied. To identify them, the"
            "value is set to -INF."
        )
        logger.notice("Number of infinities = {}".format(len(np.where(np.isinf(flood_workspace.extractY()))[0])))

        logger.debug("Moving/Patch detector algorithm.  Any masked detector pixel is set to NaN")
        logger.notice("Number of NaNs       = {}".format(len(np.where(np.isnan(flood_workspace.extractY()))[0])))

        return flood_workspace

    @staticmethod
    def sum_input_runs(flood_workspaces):
        """Do NaN sum to all input flood workspaces

        Parameters
        ----------
        flood_workspaces

        Returns
        -------

        """
        from mantid.simpleapi import CloneWorkspace

        # Do NaN sum
        y_list = [ws.extractY() for ws in flood_workspaces]
        y_matrix = np.array(y_list)

        nan_sum_matrix = np.nansum(y_matrix, axis=0)

        # clone a workspace
        cloned = CloneWorkspace(InputWorkspace=flood_workspaces[0], OutputWorkspace="FloodSum")

        for iws in range(cloned.getNumberHistograms()):
            cloned.dataY(iws)[0] = nan_sum_matrix[iws][0]

        # output
        SaveNexusProcessed(InputWorkspace=cloned, Filename="SummedFlood.nxs")

    def __init__(self, instrument, component="detector1"):
        """Initialization

        Parameters
        ----------
        instrument : str
            instrument name. One of CG2, CG3, EQSANS
        component : str
            One of the detector panels of the instrument (e.g. detector1, wing_detector, midrange_detector)
        """

        if instrument not in [CG2, CG3, EQSANS]:
            raise RuntimeError("Instrument {} is not supported".format(instrument))
        self._instrument = instrument

        # validate the component as part of the instrument
        if component not in panel_names(instrument):
            raise RuntimeError(f"Component {component} is not part of {instrument}")
        self._component = component

        # flood runs
        self._flood_runs = None
        # direct beam (center) runs
        self._direct_beam_center_runs = None
        # mask
        self._default_mask = None
        self._extra_mask_dict = dict()
        self._beam_center_radius = None  # mm

        # flag to enforce to use IDF XML in NeXus file; otherwise, it may use IDF from Mantid library
        self._enforce_use_nexus_idf = False

        # Transmission correction
        self._transmission_reference_runs = None
        self._transmission_flood_runs = None

        # Dark current
        self._dark_current_runs = None

        # Pixel calibration default
        self._apply_calibration = False

        # Scale Component default
        self._scale_components: Optional[dict] = None

        # Apply solid angle correction or not?
        self._solid_angle_correction = False

    @property
    def beam_center_radius(self) -> float:
        return self._beam_center_radius

    @property
    def scale_components(self) -> Optional[dict]:
        return self._scale_components

    @scale_components.setter
    def scale_components(self, scalings=None):
        """
        Setter method for the `scale_component` property.

        Adjust pixel positions, heights and widths with Mantid algorithm ScaleInstrumentComponent

        Parameters
        ----------
        scalings : dict or None, optional
            dictionary of scaling triads. For instance,
            `scalings={"detector1": [1.0, 2.0, 1.0], "wing_detector":[0.5, 1.0, 0.5]}`
            doubles the height of pixels in "detector1" (Y-axis is the vertical axis in the detector coordinate
            system), and halves the width of pixels in "wing_detector". No changes for "midrange detector".

        Raises
        ------
        TypeError
            If `scalings` is not a list.
        AssertionError
            If the length of `scalings` is not equal to 3.
        """
        if scalings is not None:
            assert isinstance(scalings, dict)
            self._scale_components = scalings

    def set_pixel_calibration_flag(self, apply_calibration):
        """Set the flag to apply pixel calibrations.

        Parameters
        ----------
        apply_calibration : bool, str
            Flag for applying the pixel calibration. No calibration (False), Default (True), Calibration file (str)
        """
        self._apply_calibration = apply_calibration

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
        if isinstance(flood_runs, (int, str)):
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
        elif isinstance(dark_current_runs, (int, str)):
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
        if isinstance(direct_beam_runs, (int, str)):
            # defined as a single value
            self._direct_beam_center_runs = [direct_beam_runs]
        else:
            # shall be a list or tuple
            self._direct_beam_center_runs = list(direct_beam_runs)

    def set_masks(self, default_mask, pixels):
        """Set masks

        Parameters
        ----------
        default_mask : str or ~mantid.api.MaskWorkspace or :py:obj:`list` or None
            Mask to be applied. If :py:obj:`list`, it is a list of
            detector ID's.
            mask file name
        pixels : str or None
            pixels to mask.  Example: '1-8,249-256'

        Returns
        -------
        None

        """
        # default mask file or list of detector IDS or MaskWorkspace
        if default_mask is not None:
            self._default_mask = default_mask

        # pixels to mask
        if pixels is not None:
            self._extra_mask_dict[PIXEL] = pixels

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

    def _get_beam_center_run(self, index):
        r"""
        Input beam center associated to a particular flood run.

        Parameters
        ----------
        index : int
            beam center run index mapped to flood run index

        Returns
        -------
        int, str
            beam center run number (as an int) or absolute file path (as a string)
        """
        if self._direct_beam_center_runs is None:
            raise RuntimeError("Beam center runs must be given for {}".format(self._instrument))
        return self._direct_beam_center_runs[index]

    def _prepare_data_opts(self, beam_center):
        r"""
        Set additional options for function prepare_data

        Parameters
        ----------
        beam_center : list
            beam center for each of the double panels

        Returns
        -------
        dict
        """
        opts = dict(center_x=beam_center[0], center_y=beam_center[1], flux_method=None)
        if self._instrument in [CG2, CG3]:
            opts["flux_method"] = "monitor"
            opts["overwrite_instrument"] = False
            opts["enforce_use_nexus_idf"] = self._enforce_use_nexus_idf
        return opts

    def _get_beam_center_workspace(self, beam_center_run):
        r"""
        Load and prepare the beam center run with customary corrections

        Parameters
        ----------
        beam_center_run : str
            instrument name plus run number or absolute file path

        Returns
        -------
        ~mantid.api.Workspace2D
        """
        # elucidate the name for the output workspace
        if os.path.exists(beam_center_run):
            output_workspace = f"BC_{self._instrument}_{extract_run_number(beam_center_run)}"
        else:
            output_workspace = f"BC_{beam_center_run}"

        beam_center_workspace = PREPARE_DATA[self._instrument](
            data=beam_center_run,
            pixel_calibration=self._apply_calibration,
            mask=self._default_mask,
            btp=self._extra_mask_dict,
            solid_angle=False,
            output_workspace=output_workspace,
            **self._prepare_data_opts(beam_center=[0.0, 0.0]),
        )

        return beam_center_workspace

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
            Beam center as xc, yc and possible wc for BIOSANS

        """
        # Process the input run chosen to calculate the beam center
        beam_center_run = self._get_beam_center_run(index)
        if isinstance(beam_center_run, str) and beam_center_run.isdigit():
            beam_center_run = int(beam_center_run)  # convert to int if possible

        if isinstance(beam_center_run, int):
            beam_center_run = "{}_{}".format(self._instrument, beam_center_run)
        elif isinstance(beam_center_run, str):
            assert os.path.exists(beam_center_run), f"Bean center run {beam_center_run} cannot be found"
        else:
            raise RuntimeError(f"Beam center run {beam_center_run} is not recognized")

        beam_center_workspace = self._get_beam_center_workspace(beam_center_run)
        beam_center = FIND_BEAM_CENTER[self._instrument](beam_center_workspace)
        return beam_center[:-1]  # last item is the fit results, useless here

    def _prepare_flood_data(self, flood_run, beam_center, dark_current_run):
        """Prepare flood data including
        (1) load
        (2) mask: default, pixels
        (3) center detector
        (4) optionally solid angle correction
        (5) optionally dark current correction
        (6) normalization

        Parameters
        ----------
        flood_run: int, str
            flood run number of flood file path

        beam_center

        dark_current_run

        enforce_use_nexus_idf: bool
            flag to enforce to use IDF XML in NeXus file; otherwise, it may use IDF from Mantid library

        Returns
        -------

        """
        if dark_current_run is not None:
            if isinstance(dark_current_run, str) and os.path.exists(dark_current_run):
                pass  # dark current run (given) is a data file: do nothing
            else:
                # dark current is a run number either as an integer or a string cast from integer
                dark_current_run = "{}_{}".format(self._instrument, dark_current_run)

        # Load data with masking: returning to a list of workspace references
        # processing includes: load, mask, normalize by monitor
        if isinstance(flood_run, str) and flood_run.isdigit():  # run number passed on as a string
            flood_run = int(flood_run)  # unequivocally convert run numbers to integers
        if isinstance(flood_run, int):
            flood_run = "{}_{}".format(self._instrument, flood_run)
        elif isinstance(flood_run, str):  # if still a string, it must be an absolute file path
            assert os.path.exists(flood_run)
        else:
            raise ValueError(f"Flood run {flood_run} is not recognized")

        flood_ws = PREPARE_DATA[self._instrument](
            data=flood_run,
            scale_components=self.scale_components,
            pixel_calibration=self._apply_calibration,
            mask=self._default_mask,
            btp=self._extra_mask_dict,
            dark_current=dark_current_run,
            solid_angle=self._solid_angle_correction,
            **self._prepare_data_opts(beam_center),
        )

        # Integrate all the wavelength bins if necessary
        if flood_ws.blocksize() != 1:
            flood_ws = Integration(InputWorkspace=flood_ws, OutputWorkspace=str(flood_ws))

        return flood_ws

    def _mask_beam_center(self, flood_ws, beam_center):
        """Mask beam center

        Mask beam center with two algorithms, tried in sequence:
        1. if beam center mask file is present, use it. It's assumed it will mask the beam center pixels and
           all detector panels but the one of interest.
        1. Mask all detector pixels within a distance of property "beam_center_radius" from the beam axis.

        Parameters
        ----------
        flood_ws : ~mantid.api.MatrixWorkspace
            Mantid workspace for flood data
        beam_center : tuple or str
            if tuple, beam centers (xc, yc, wc) / (xc, yc); str: beam center masks file
        Returns
        -------

        """
        if isinstance(beam_center, str):
            apply_mask(flood_ws, mask=beam_center)
        else:
            masking = list(circular_mask_from_beam_center(flood_ws, self._beam_center_radius))
            apply_mask(flood_ws, mask=masking)  # data_ws reference shall not be invalidated here!

        # Set uncertainties
        # output: masked are zero intensity and zero error
        masked_flood_ws = set_init_uncertainties(flood_ws)

        return masked_flood_ws

    def _export_sensitivity(self, sensitivity_ws, output_nexus_name, parent_flood_run):
        """Process and export sensitivities to a processed NeXus file

        Parameters
        ----------
        sensitivity_ws :  ~mantid.api.MatrixWorkspace
            MatrixWorkspace containing sensitivity and error
        output_nexus_name : str
            Output NeXus file  name
        parent_flood_run : int
            Flood run number to create parent workspace for sensitivity workspace

        Returns
        -------
        None

        """
        # Create a new workspace for output
        instrument_name = {CG2: "GPSANS", CG3: "BIOSANS", EQSANS: "EQSANS_"}[self._instrument]
        if isinstance(parent_flood_run, int):
            event_nexus = "{}{}".format(instrument_name, parent_flood_run)
        else:
            # must be a nexus file already
            event_nexus = parent_flood_run
            assert os.path.exists(event_nexus)

        parent_ws = load_events(run=event_nexus, MetaDataOnly=True, LoadNexusInstrumentXML=self._enforce_use_nexus_idf)

        # Create new sensitivity workspace
        new_sens_name = "{}_new".format(str(sensitivity_ws))
        new_sensitivity_ws = CreateWorkspace(
            DataX=sensitivity_ws.extractX().flatten(),
            DataY=sensitivity_ws.extractY().flatten(),
            DataE=sensitivity_ws.extractE().flatten(),
            NSpec=parent_ws.getNumberHistograms(),
            ParentWorkspace=parent_ws,
            OutputWorkspace=new_sens_name,
        )

        # Mask detectors
        mask_ws_indexes = list()
        for iws in range(new_sensitivity_ws.getNumberHistograms()):
            # get the workspace with -infinity or NaN for masking
            if np.isnan(new_sensitivity_ws.readY(iws)[0]) or np.isinf(new_sensitivity_ws.readY(iws)[0]):
                mask_ws_indexes.append(iws)
        MaskDetectors(Workspace=new_sensitivity_ws, WorkspaceIndexList=mask_ws_indexes)

        # Set all the mask values to NaN
        new_sensitivity_ws = mtd[new_sens_name]
        for iws in mask_ws_indexes:
            new_sensitivity_ws.dataY(iws)[0] = np.nan

        # Save
        SaveNexusProcessed(InputWorkspace=new_sensitivity_ws, Filename=output_nexus_name)

    def _apply_transmission_correction(self, *args):
        r"""Calculate and apply transmission correction

        Warnings
        --------
        This method is implemented only for CG3

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError("Transmission correction is not implemented")

    def execute(
        self,
        use_moving_detector_method,
        min_threshold,
        max_threshold,
        output_nexus_name,
        enforce_use_nexus_idf=False,
        debug_mode=False,
    ):
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
        enforce_use_nexus_idf: bool
            flag to enforce to use IDF XML in NeXus file; otherwise, it may use IDF from Mantid library
        debug_mode: bool
            flag for debugging mode

        Returns
        -------
        None

        """
        # Number of pair of workspaces to process
        num_workspaces_set = len(self._flood_runs)

        self._enforce_use_nexus_idf = enforce_use_nexus_idf

        # Load beam center runs and calculate beam centers
        beam_centers = list()
        for i in range(num_workspaces_set):
            beam_center_i = self._calculate_beam_center(i)
            beam_centers.append(beam_center_i)
            logger.notice("Calculated beam center ({}-th) = {}".format(i, beam_center_i))

        # Set default value to dark current runs
        if self._dark_current_runs is None:
            self._dark_current_runs = [None] * num_workspaces_set

        # Load and process flood data with (1) mask (2) center detector and (3) solid angle correction
        flood_workspaces = list()
        for i in range(num_workspaces_set):
            flood_ws_i = self._prepare_flood_data(self._flood_runs[i], beam_centers[i], self._dark_current_runs[i])
            flood_workspaces.append(flood_ws_i)
            logger.notice(f"Load {i}-th flood run {self._flood_runs[i]} to {flood_ws_i}")

        # Retrieve masked detectors before masking the beam center. These are termed "bad pixels"
        if not use_moving_detector_method:
            bad_pixels_list = list()
            for i in range(num_workspaces_set):
                bad_pixels_list.append(self._get_masked_detectors(flood_workspaces[i]))
        else:
            bad_pixels_list = [None] * num_workspaces_set

        # Mask beam centers
        for i in range(num_workspaces_set):
            flood_workspaces[i] = self._mask_beam_center(flood_workspaces[i], beam_centers[i])

        # Transmission correction as an option
        if self._transmission_reference_runs is not None:
            for i in range(num_workspaces_set):
                flood_workspaces[i] = self._apply_transmission_correction(
                    flood_workspaces[i],
                    self._transmission_reference_runs[i],
                    self._transmission_flood_runs[i],
                    beam_centers[i],
                )

        # Set the masked pixels' counts to nan and -infinity
        for i in range(num_workspaces_set):
            flood_workspaces[i] = self._set_mask_value(
                flood_workspaces[i], bad_pixels_list[i], use_moving_detector_method
            )

        info = "Preparation of data is over.\n"
        for fws in flood_workspaces:
            info += (
                f"{str(fws)}: Number of infinities = {len(np.where(np.isinf(fws.extractY()))[0])},"
                f"Number of NaNs = {len(np.where(np.isnan(fws.extractY()))[0])}\n"
            )
        logger.notice(info)

        # Debug output
        if debug_mode:
            for flood_ws in flood_workspaces:
                SaveNexusProcessed(InputWorkspace=flood_ws, Filename=f"{str(flood_ws)}_flood.nxs")

        # Decide algorithm to prepare sensitivities
        if self._instrument in [CG2, CG3] and use_moving_detector_method is True:
            if debug_mode:
                # nan.sum all the input flood runs to check the coverage of summed counts
                self.sum_input_runs(flood_workspaces)

            # Prepare by using moving detector algorithm
            calculate_sensitivity_correction = CALCULATE_SENSITIVITY_CORRECTION[MOVING_DETECTORS]

            # Calculate sensitivities for each file
            sens_ws = calculate_sensitivity_correction(
                flood_workspaces,
                threshold_min=min_threshold,
                threshold_max=max_threshold,
            )

        else:
            # Prepare by using the sensitivity patch method for a single detector (image)
            # Such as GPSANS, BIOSANS Main detector, BIOSANS wing detector, EQSANS
            calculate_sensitivity_correction = CALCULATE_SENSITIVITY_CORRECTION[PATCHING_DETECTORS]

            # Default polynomial order: CG3 uses order 3.  Others use order 2.
            if self._instrument == CG3:
                polynomial_order = 3
            else:
                polynomial_order = 2

            # This only processes a single image, even for the Bio-SANS.
            # Each detector on the Bio-SANS must be treated independently
            sens_ws = calculate_sensitivity_correction(
                flood_workspaces[0],
                min_threshold=min_threshold,
                max_threshold=max_threshold,
                poly_order=polynomial_order,
                min_detectors_per_tube=50,
                component_name=self._component,
            )

        # Export
        self._export_sensitivity(sens_ws, output_nexus_name, self._flood_runs[0])


def debug_output(workspace, output_file):
    """Exporting a workspace to NeXus file and HDF5 for debugging purpose

    Parameters
    ----------
    workspace : numpy.ndarray
        data to plot
    output_file : str
        output file name as reference

    Returns
    -------
    None

    """
    # Save Nexus
    SaveNexusProcessed(InputWorkspace=workspace, Filename=output_file)

    data = workspace.extractY()
    data_error = workspace.extractE()

    # Export sensitivities calculated in file for quick review
    # Export to hdf5 for a quick review
    sens_h5 = h5py.File("{}.h5".format(output_file.split(".")[0]), "w")
    sens_group = sens_h5.create_group("Data")
    sens_group.create_dataset("Data", data=data)
    if data_error is not None:
        sens_group.create_dataset("Data error", data=data_error)
    sens_h5.close()
