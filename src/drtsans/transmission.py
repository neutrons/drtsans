import numpy as np
import sys

from mantid import mtd
from mantid.kernel import logger

# https://docs.mantidproject.org/nightly/algorithms/ApplyTransmissionCorrection-v1.html
# https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
# https://docs.mantidproject.org/nightly/algorithms/GroupDetectors-v1.html
# https://docs.mantidproject.org/nightly/algorithms/RebinToWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html
from mantid.simpleapi import (
    ApplyTransmissionCorrection,
    Divide,
    GroupDetectors,
    RebinToWorkspace,
    ReplaceSpecialValues,
)

r""" links to drtsans imports
circular_mask_from_beam_center, masked_detectors available at:
    <https://github.com/neutrons/drtsans/blob/next/src/drtsans/mask_utils.py>
beam_radius <https://github.com/neutrons/drtsans/blob/next/src/drtsans/geometry.py>
"""  # noqa: E501
from drtsans.mask_utils import circular_mask_from_beam_center, masked_detectors

# Symbols to be exported
__all__ = [
    "apply_transmission_correction",
    "calculate_transmission",
    "TransmissionErrorToleranceError",
    "TransmissionNanError",
]


class TransmissionErrorToleranceError(Exception):
    """Exception raised when the transmission error is larger than the transmission error tolerance"""

    pass


class TransmissionNanError(Exception):
    """Exception raised when all transmission values are NaN"""

    pass


def calculate_transmission(
    input_sample, input_reference, radius, radius_unit="mm", transmission_error_tolerance=None, output_workspace=None
):
    """
    Calculate the raw transmission coefficients at zero scattering angle
    from already prepared sample and reference data.

    For EQ-SANS, one additional step fitting the returned raw values is
    necessary. Use `eqsans.calculate_transmission` instead.

    **Mantid algorithms used:**
    :ref:`Divide <algm-Divide-v1>`
    :ref:`GroupDetectors <algm-GroupDetectors-v2>`
    :ref:`RebinToWorkspace <algm-RebinToWorkspace-v1>`


    Parameters
    ----------
    input_sample: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        Prepared sample workspace (possibly obtained with an attenuated beam)
    input_reference: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        Prepared direct beam workspace (possibly obtained with an attenuated beam)
    radius: float
        Radius around the beam center for pixel integration, in millimeters.
        If None, radius will be obtained or calculated using `input_reference` workspace.
    radius_unit: str
        Either 'mm' or 'm', and only used in conjunction with option `radius`.
    transmission_error_tolerance: float | None
        Maximum relative error for transmission. If None, the error will not be checked.
    output_workspace: str
        Name of the output workspace containing the raw transmission values.
        If None, a hidden random name will be provided.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        Workspace containing the raw transmission values

    Raises
    ------
    TransmissionNanError
        If all transmission values are NaN.
    TransmissionToleranceError
        If there are insufficient statistics to calculate the transmission correction.
    """
    if output_workspace is None:
        output_workspace = mtd.unique_hidden_name()

    if radius is None:
        logger.information("Calculating beam radius from sample logs")
        raise NotImplementedError("beam radius must be specified")
    else:
        radius = float(radius) if radius_unit == "mm" else 1.0e3 * radius  # to mm
        logger.information("beam radius is (mm) {}".format(radius))
    if radius <= 0.0:
        raise ValueError("Encountered negative beam radius={}mm".format(radius))

    # Find the identity of the detector pixels falling within the beam area
    detector_ids = circular_mask_from_beam_center(input_reference, radius, unit="mm")
    if not detector_ids:
        raise RuntimeError("No pixels in beam with radius of {:.2f} mm".format(radius))

    # Warn when masking many pixels around the beam center
    warning_message = (
        "Warning: More than half of the detectors within a radius of {:.2f} mm ".format(radius)
        + "from the beam center are masked in the input {0}"
    )
    for run, workspace in dict(sample=input_sample, reference=input_reference).items():
        if len(masked_detectors(workspace, detector_ids)) > len(detector_ids) / 2:
            sys.stderr.write(warning_message.format(run))

    # Add the intensities of the detector pixels within the beam area
    sample_intensity_workspace = GroupDetectors(
        InputWorkspace=input_sample,
        DetectorList=detector_ids,
        OutputWorkspace=mtd.unique_hidden_name(),
    )
    reference_intensity_workspace = GroupDetectors(
        InputWorkspace=input_reference,
        DetectorList=detector_ids,
        OutputWorkspace=mtd.unique_hidden_name(),
    )

    # If the reference workspace used a different wavelength binning than that of the sample workspace, a rebinning
    # step is necessary prior to dividing sample intensities by the reference intensities.
    reference_intensity_workspace = RebinToWorkspace(
        WorkspaceToRebin=reference_intensity_workspace,
        WorkspaceToMatch=sample_intensity_workspace,
        OutputWorkspace=reference_intensity_workspace.name(),
    )

    # RebinToWorkspace may spill some intensity in the reference workspace in the region of wavelengths
    # corresponding to the gap between the lead and skip pulses. We have to harmonize the gap of the
    # reference workspace to that of the sample workspace
    gap_indexes = np.where(sample_intensity_workspace.dataY(0) == 0.0)
    reference_intensity_workspace.dataY(0)[gap_indexes] = 0.0

    # calculate zero angle transmission coefficient(s)
    zero_angle_transmission_workspace = Divide(
        LHSWorkspace=sample_intensity_workspace,
        RHSWorkspace=reference_intensity_workspace,
        OutputWorkspace=output_workspace,
    )

    # Notify of incorrect calculation of zero angle transmission
    # Will happen if the beam centers have been totally masked
    if bool(np.all(np.isnan(zero_angle_transmission_workspace.readY(0)))) is True:
        raise TransmissionNanError("Transmission at zero-angle is NaN")

    non_gap_indexes = np.isfinite(zero_angle_transmission_workspace.readY(0))

    if transmission_error_tolerance:
        # Verify that errors are below the transmission error tolerance
        transmission_relative_error = (
            zero_angle_transmission_workspace.readE(0)[non_gap_indexes]
            / zero_angle_transmission_workspace.readY(0)[non_gap_indexes]
        )
        if not np.all(transmission_relative_error < transmission_error_tolerance):
            i_max = np.argmax(transmission_relative_error)
            rel_error = transmission_relative_error[i_max]
            transmission = zero_angle_transmission_workspace.readY(0)[non_gap_indexes][i_max]
            abs_error = zero_angle_transmission_workspace.readE(0)[non_gap_indexes][i_max]
            raise TransmissionErrorToleranceError(
                f"transmission_error / transmission_value ({abs_error:.4f} / {transmission:.4f} = "
                f"{rel_error:.4f}) > transmission_relative_error_tolerance "
                f"({transmission_error_tolerance:.4f})"
            )

    # Notify of average transmission value
    average_zero_angle_transmission = np.mean(zero_angle_transmission_workspace.readY(0)[non_gap_indexes])
    average_zero_angle_transmission_error = np.linalg.norm(zero_angle_transmission_workspace.readE(0)[non_gap_indexes])
    message = "Average zero angle transmission = {0} +/- {1}".format(
        average_zero_angle_transmission, average_zero_angle_transmission_error
    )
    logger.warning(message)

    # A bit of clean up
    sample_intensity_workspace.delete()
    reference_intensity_workspace.delete()

    return zero_angle_transmission_workspace


def apply_transmission_correction(
    input_workspace,
    trans_workspace=None,
    trans_value=None,
    trans_error=0.0,
    theta_dependent=True,
    output_workspace=None,
):
    r"""
    Correct the intensities with transmission coefficient(s).

    **Mantid algorithms used:**
    :ref:`ApplyTransmissionCorrection <algm-ApplyTransmissionCorrection-v1>`
    :ref:`ReplaceSpecialValues <algm-ReplaceSpecialValues-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Input workspace to correct its intensities
    trans_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace containing the transmission coefficient(s). The result of applying `calculate_transmission`
        to the input workspace. If :py:obj:`None`, `trans_value` will be used.
    trans_value: float
        A single transmission coefficient to correct the intensities. If :py:obj:`None`,
        `trans_workspace` will be used.
    trans_error: float
        Error associated to `trans_value`.
    theta_dependent : bool
        Flag to do theta dependent correction
    output_workspace: str
        Name of the workspace containing the corrected intensities. If :py:obj:`None`, the `input_workspace`
        will be overwritten.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    # kwargs is a list of options to be passed on to Mantid algorithm ApplyTransmissionCorrection
    kwargs = dict(
        InputWorkspace=input_workspace,
        ThetaDependent=theta_dependent,
        OutputWorkspace=output_workspace,
    )

    if trans_workspace is not None:
        # EQ-SANS transmissions in skip-frame mode have transmission values of zero in the wavelength gap.
        # Need to be replaced with one to avoid division of intensities by zero.
        clean_trans_workspace = ReplaceSpecialValues(
            InputWorkspace=trans_workspace,
            SmallNumberThreshold=1.0e-6,
            SmallNumberValue=1.0,
            OutputWorkspace=mtd.unique_hidden_name(),
        )
        kwargs["TransmissionWorkspace"] = clean_trans_workspace
    elif trans_value is not None:  # we are passing a single value for the transmission
        kwargs.update(dict(TransmissionValue=trans_value, TransmissionError=trans_error))
    else:  # we neither passed a transmission workspace nor a single transmission value
        raise RuntimeError("Provide either trans_workspace or trans_value")

    ApplyTransmissionCorrection(**kwargs)

    if trans_workspace is not None:
        clean_trans_workspace.delete()

    return mtd[output_workspace]
