# local imports
from drtsans.mask_utils import circular_mask_from_beam_center, apply_mask
from drtsans.mono.biosans.api import prepare_data
from drtsans.mono.biosans import apply_transmission_correction, calculate_transmission
from drtsans.mono.spice_data import SpiceRun
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection as PrepareBase
from drtsans.process_uncertainties import set_init_uncertainties

# third party imports
from mantid.kernel import logger
from mantid.simpleapi import MaskAngle
import numpy as np

# standard imports
import os
from typing import Union


# Constants
PIXEL = "Pixel"


class PrepareSensitivityCorrection(PrepareBase):
    def __int__(self):
        super().__init__()

    def set_masks(self, default_mask, pixels, wing_det_mask_angle=None, main_det_mask_angle=None):
        """Set masks

        Parameters
        ----------
        default_mask : str or ~mantid.api.MaskWorkspace or :py:obj:`list` or None
            Mask to be applied. If :py:obj:`list`, it is a list of
            detector ID's. If `None`, it is expected that `maskbtp` is not empty.
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

        super().set_masks(default_mask, pixels)

        if wing_det_mask_angle is not None:
            self._wing_det_mask_angle = wing_det_mask_angle
        if main_det_mask_angle is not None:
            self._main_det_mask_angle = main_det_mask_angle

    def _get_beam_center_run(self, index):
        r"""
        Input beam center associated to a particular flood run.

        if direct beam center runs are not defined, then the flood run itself is used as direct beam center run

        Parameters
        ----------
        index : int
            beam center run index mapped to flood run index

        Returns
        -------
        int, str
            beam center run number (as an int) or absolute file path (as a string)
        """
        if self._direct_beam_center_runs is None and self._wing_det_mask_angle is not None:
            return self._flood_runs[index]
        else:
            return super()._get_beam_center_run(index)

    def _prepare_data_opts(self, beam_center):
        r"""Set additional options for function prepare_data

        Parameters
        ----------
        beam_center : List[float]
            beam center for the main detector (X,Y), the wing detector (Y) and optionally the midrange detector (Y)

        Returns
        -------
        dict
        """
        opts = super()._prepare_data_opts(beam_center)
        if len(beam_center) > 2:
            opts["center_y_wing"] = beam_center[2]
        if len(beam_center) > 3:
            opts["center_y_midrange"] = beam_center[3]
        return opts

    def _get_beam_center_workspace(self, beam_center_run):
        r"""
        Load and prepare the beam center run with customary corrections. Also masks the curved detectors.

        Parameters
        ----------
        beam_center_run: str
            instrument name plus run number or absolute file path

        Returns
        -------
        ~mantid.api.Workspace2D
        """
        beam_center_workspace = super()._get_beam_center_workspace(beam_center_run)

        if self._wing_det_mask_angle is not None:
            # mask wing detector
            btp = dict(Components="wing_detector")
            apply_mask(beam_center_workspace, **btp)
            # mask 2-theta angle on main detector
            MaskAngle(
                Workspace=beam_center_workspace,
                MinAngle=self._wing_det_mask_angle,
                Angle="TwoTheta",
            )

        return beam_center_workspace

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

        # TODO (jose borreguero): suspicious "if" block.
        if self._main_det_mask_angle is not None:
            # Mask 2-theta angle
            # Mask wing detector right top/bottom corners
            if self._component == "wing_detector":
                component_to_mask = "detector1"
            else:
                component_to_mask = "wing_detector"
            apply_mask(flood_ws, Components=component_to_mask)
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

    def set_transmission_correction(self, transmission_flood_runs, transmission_reference_runs, beam_trap_factor=2):
        """Set transmission beam run and transmission flood runs

        Parameters
        ----------
        transmission_flood_runs : int or tuple or list
            transmission flood runs
        transmission_reference_runs : int or tuple or list
            transmission reference runs
        beam_trap_factor : float, int
            factor to beam trap size for masking angle

        Returns
        -------

        """

        def format_run_or_runs(run_s):
            """Format input run or runs to list of run or file names

            Parameters
            ----------
            run_s: int or str or ~list
                a run, a NeXus file name, or a list of runs

            Returns
            -------
            ~list

            """
            if isinstance(run_s, (int, str)):
                # an integer or string as run number
                run_list = [run_s]
            else:
                # a sequence, tuple or list
                run_list = list(run_s)

            return run_list

        # transmission reference
        self._transmission_reference_runs = format_run_or_runs(transmission_reference_runs)

        # transmission flood
        self._transmission_flood_runs = format_run_or_runs(transmission_flood_runs)

        # if isinstance(transmission_reference_run, int):
        #     # a run number
        #     self._transmission_reference_runs = [transmission_reference_run]
        # else:
        #     self._transmission_reference_runs = list(transmission_reference_run)
        #
        # if isinstance(transmission_flood_runs, int):
        #     self._transmission_flood_runs = [transmission_flood_runs]
        # else:
        #     self._transmission_flood_runs = list(transmission_flood_runs)

        # Set the beam trap factor for transmission reference and flood run to mask angle
        self._biosans_beam_trap_factor = beam_trap_factor

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

    def _apply_transmission_correction(self, flood_ws, transmission_beam_run, transmission_flood_run, beam_center):
        """Calculate and apply transmission correction

        Parameters
        ----------
        flood_ws : MarixWorkspace
            Flood run workspace to transmission correct workspace
        transmission_beam_run : int or str
            run number for transmission beam run
        transmission_flood_run : int or str
            run number for transmission flood run
        beam_center : ~tuple
            detector center

        Returns
        -------
        MatrixWorkspace
            Flood workspace with transmission corrected

        """
        if isinstance(transmission_beam_run, str) and os.path.exists(transmission_beam_run):
            sans_data = transmission_beam_run
        elif (
            isinstance(transmission_beam_run, str)
            and transmission_beam_run.isdigit()
            or isinstance(transmission_beam_run, int)
        ):
            sans_data = "{}_{}".format(self._instrument, transmission_beam_run)
        else:
            raise TypeError(
                f"Transmission run {transmission_beam_run} of type {type(transmission_beam_run)} "
                f"is not supported to load a NeXus run from it"
            )

        transmission_workspace = prepare_data(
            data=sans_data,
            pixel_calibration=self._apply_calibration,
            mask=self._default_mask,
            btp=self._extra_mask_dict,
            solid_angle=False,
            output_workspace="TRANS_{}_{}".format(self._instrument, transmission_beam_run),
            **self._prepare_data_opts(beam_center),
        )

        # Apply mask
        apply_mask(transmission_workspace, Components="wing_detector")
        MaskAngle(
            Workspace=transmission_workspace,
            MinAngle=self._biosans_beam_trap_factor * self._main_det_mask_angle,
            Angle="TwoTheta",
        )

        if not os.path.exists(transmission_flood_run):
            # given run number: form to CG3_XXX
            mtd_trans_run = "{}_{}".format(self._instrument, transmission_flood_run)
        else:
            # already a file path
            mtd_trans_run = transmission_flood_run

        transmission_flood_ws = prepare_data(
            data=mtd_trans_run,
            pixel_calibration=self._apply_calibration,
            mask=self._default_mask,
            btp=self._extra_mask_dict,
            solid_angle=False,
            output_workspace="TRANS_{}_{}".format(self._instrument, transmission_flood_run),
            **self._prepare_data_opts(beam_center),
        )

        # Apply mask
        apply_mask(transmission_flood_ws, Components="wing_detector")
        MaskAngle(
            Workspace=transmission_flood_ws,
            MinAngle=self._biosans_beam_trap_factor * self._main_det_mask_angle,
            Angle="TwoTheta",
        )

        # Zero-Angle Transmission Co-efficients
        transmission_corr_ws = calculate_transmission(transmission_flood_ws, transmission_workspace)
        average_zero_angle = np.mean(transmission_corr_ws.readY(0))
        average_zero_angle_error = np.linalg.norm(transmission_corr_ws.readE(0))
        logger.notice(
            f"Transmission Coefficient is {average_zero_angle:.3f} +/- "
            f"{average_zero_angle_error:.3f}."
            f"Transmission flood {str(transmission_flood_ws)} and "
            f"transmission {str(transmission_workspace)}"
        )

        # Apply calculated transmission
        flood_ws = apply_transmission_correction(
            flood_ws,
            trans_workspace=transmission_corr_ws,
            theta_dependent=self._theta_dep_correction,
        )

        return flood_ws


def prepare_spice_sensitivities_correction(
    component: str,
    flood_run: SpiceRun,
    direct_beam_run: SpiceRun,
    dark_current_run: SpiceRun,
    apply_solid_angle_correction: bool,
    transmission_flood_run: SpiceRun,
    transmission_reference_run: SpiceRun,
    beam_trap_size_factor: float,
    apply_theta_dependent_correction: bool,
    universal_mask_file: Union[str, None],
    pixels_to_mask: str,
    beam_center_mask_radius: float,
    main_detector_mask_angle: float,
    wing_detector_mask_angle: float,
    min_count_threshold: float,
    max_count_threshold: float,
    sensitivity_file_name: str,
    nexus_dir: str = None,
):
    """Prepare sensitivities from SPICE files

    Parameters
    ----------
    component: str
        Detector panel to use. One of "detector1" or "wing_detector"
    flood_run: SpiceRun
        flood run
    direct_beam_run: SpiceRun or None
        direct beam run
    dark_current_run: SpiceRun
        dark current run
    apply_solid_angle_correction: bool
        Flag to apply solid angle correction to flood run
    transmission_flood_run: SpiceRun
        transmission flood run
    transmission_reference_run: SpiceRun
        transmission reference run
    beam_trap_size_factor: float
        size factor of beam trap given by user
    apply_theta_dependent_correction: bool
        Flag to apply theta dependent correction to transmission run
    universal_mask_file: str
        path to mask file applied to all the runs
    pixels_to_mask: str
        lists of pixels (IDs) to mask
    beam_center_mask_radius: float
        radius of mask for beam center in mm
    main_detector_mask_angle: float
        angle for main detector mask
    wing_detector_mask_angle: float
        angle for wing detector mask
    min_count_threshold: float
        minimum normalized count threshold as a good pixel
    max_count_threshold: float
        maximum normalized count threshold as a good pixel
    sensitivity_file_name: str
        output file name with full path
    nexus_dir: str or None
        directory for nexus file.  None for default.

    """

    CG3 = "CG3"

    # Set up sensitivities preparation configurations
    if component not in ["detector1", "wing_detector"]:
        raise ValueError(f"Unknown component {component}. Must be one of 'detector1' or 'wing_detector'")

    preparer = PrepareSensitivityCorrection(CG3, component=component)
    # Load flood runs
    preparer.set_flood_runs([flood_run.unique_nexus_name(nexus_dir, True)])

    # Process beam center runs
    if direct_beam_run is not None:
        preparer.set_direct_beam_runs([direct_beam_run.unique_nexus_name(nexus_dir, True)])

    # Set extra masks
    preparer.set_masks(
        universal_mask_file,
        pixels_to_mask,
        wing_det_mask_angle=wing_detector_mask_angle,
        main_det_mask_angle=main_detector_mask_angle,
    )

    # Set beam center radius
    if beam_center_mask_radius is not None:
        preparer.set_beam_center_radius(beam_center_mask_radius)

    # Transmission
    if transmission_reference_run is not None:
        trans_flood_file = transmission_flood_run.unique_nexus_name(nexus_dir, True)
        trans_ref_file = transmission_reference_run.unique_nexus_name(nexus_dir, True)
        preparer.set_transmission_correction(
            transmission_flood_runs=[trans_flood_file],
            transmission_reference_runs=[trans_ref_file],
            beam_trap_factor=beam_trap_size_factor,
        )
        preparer.set_theta_dependent_correction_flag(apply_theta_dependent_correction)

    # Dark runs
    if dark_current_run is not None:
        dark_curr_file = dark_current_run.unique_nexus_name(nexus_dir, True)
        preparer.set_dark_current_runs([dark_curr_file])

    # solid angle correction
    preparer.set_solid_angle_correction_flag(apply_solid_angle_correction)

    # Run
    moving_detector = False

    preparer.execute(
        moving_detector,
        min_count_threshold,
        max_count_threshold,
        sensitivity_file_name,
        enforce_use_nexus_idf=True,
        debug_mode=True,
    )
