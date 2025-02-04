from mantid.kernel import logger
from drtsans.mono.geometry import beam_radius
from drtsans.transmission import calculate_transmission as raw_calculate_transmission
from drtsans.transmission import apply_transmission_correction

# Symbols to be exported to the eqsans namespace
__all__ = ["calculate_transmission", "apply_transmission_correction"]


def calculate_transmission(
    input_sample,
    input_reference,
    radius=None,
    radius_unit="mm",
    transmission_error_tolerance=None,
    output_workspace=None,
):
    """
    Calculate the raw transmission coefficients at zero scattering angle
    from already prepared sample and reference data.

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
        Maximum relative error for transmission
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
        If all transmission values are NaN
    TransmissionToleranceError
        If there is insufficient statistics to calculate the transmission correction
    """
    if radius is None:
        logger.information("Calculating beam radius from sample logs")
        radius = beam_radius(input_reference, unit="mm")

    zero_angle_transmission_workspace = raw_calculate_transmission(
        input_sample,
        input_reference,
        radius,
        radius_unit=radius_unit,
        transmission_error_tolerance=transmission_error_tolerance,
        output_workspace=output_workspace,
    )

    return zero_angle_transmission_workspace
