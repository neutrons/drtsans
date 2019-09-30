from mantid.simpleapi import (mtd, CreateSingleValuedWorkspace, GroupDetectors, Divide, RebinToWorkspace,
                              DeleteWorkspace, Multiply, DeleteWorkspace)
from drtsans.settings import unique_workspace_dundername
from drtsans.geometry import masked_detectors
from drtsans import transmission
from drtsans.mask_utils import circular_mask_from_beam_center
from drtsans.settings import unique_workspace_dundername as uwd


def empty_beam_intensity(empty_beam_workspace, beam_radius=None, unit='mm', roi=None,
                         attenuator_coefficient=1.0, attenuator_error=0.0, output_workspace=None):
    r"""
    Calculate the intensity impinging on the detector, taking into account attenuation.

    It is assumed the center of the detector has been moved to coincide with the center of the beam.

    Parameters
    ----------
    empty_beam_workspace: str, MatrixWorkspace, EventsWorkspace
        Attenuated intensities collected at the detector with and empty beam.
    beam_radius: float
        Radius of the beam at the detector. If None, it will be estimated with the sample and source apertures.
    unit: str
        Units for the beam radius, either meters ('m') or mili-miters ('mm').
    roi: file path, MaskWorkspace, list
        Region of interest where to collect intensities. If `list`, it is a list of detector ID's. This option
        overrides beam radius.
    attenuator_coefficient: float
        Fraction of the neutrons allowed to pass through the attenuator. Assumed wavelength independent.
    attenuator_error: float
        Estimated error for the attenuator coefficient.
    output_workspace: str
        Name of the workspace containing the calculated intensity. If None, a random hidden name will be
        automatically provided.

    Returns
    -------
    MatrixWorkspace, EventsWorkspace
        A one-spectrum workspace containing either:
        - all events collected within the beam radius if empty_beam_workspace is an EventsWorkspace.
        - intensity versus wavelength spectrum, if empty_beam_workspace is a MatrixWorkspace for a TOF instrument.
        - intensity spectrum with one bin only, if empty_beam_workspace is a MatrixWorkspace for a MONO instrument.
    """
    if output_workspace is None:
        output_workspace = unique_workspace_dundername()

    # Obtain the beam radius from the logs or calculate from the source and sample apertures
    if beam_radius is None:
        beam_radius = transmission.beam_radius(empty_beam_workspace, unit=unit)

    det_ids = circular_mask_from_beam_center(empty_beam_workspace, beam_radius, unit=unit)

    # Avoid pitfall of having too many pixels masked within the beam radius
    if len(masked_detectors(empty_beam_workspace, det_ids)) > len(det_ids) / 2:
        msg = 'More than half of the detectors within a radius of {:.2f} {} '.format(beam_radius, unit) + \
              'from the beam center are masked in the empty beam workspace'
        raise RuntimeError(msg.format)

    # Integrate the intensity within the beam radius, then divide by the attenuation factor
    GroupDetectors(InputWorkspace=empty_beam_workspace, DetectorList=det_ids, OutputWorkspace=output_workspace)
    chi = CreateSingleValuedWorkspace(DataValue=attenuator_coefficient, ErrorValue=attenuator_error)
    Divide(LHSWorkspace=output_workspace, RHSWorkspace=chi, OutputWorkspace=output_workspace)

    return mtd[output_workspace]


def empty_beam_scaling(input_workspace, empty_beam_workspace, beam_radius=None, unit='mm',
                       attenuator_coefficient=1.0, attenuator_error=0.0, output_workspace=None):
    r"""
    Normalize input workspace by the intensity impinging on the detector for an empty beam run,
    taking into account attenuation.

    Parameters
    ----------
    input_workspace: str, MatrixWorkspace, EventsWorkspace
        Workspace to be normalized
    empty_beam_workspace: str, MatrixWorkspace, EventsWorkspace
        Attenuated intensities collected at the detector with and empty beam.
    beam_radius: float
        Radius of the beam at the detector. If None, it will be estimated with the sample and source apertures.
    unit: str
        Units for the beam radius, either meters ('m') or mili-miters ('mm').
    attenuator_coefficient: float
        Fraction of the neutrons allowed to pass through the attenuator. Assumed wavelength independent.
    attenuator_error: float
        Estimated error for the attenuator coefficient.
    output_workspace: str
        Name of the normalized workspace. If ``None``, then the name of ``input_workspace`` will be used,
        thus overwriting ``input_workspace``.

    Returns
    -------
        MatrixWorkspace, EventsWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    # Calculate the intensity impinging on the detector, taking into account attenuation.
    beam_intensity = unique_workspace_dundername()  # temporary workspace
    empty_beam_intensity(empty_beam_workspace, beam_radius=beam_radius, unit=unit,
                         attenuator_coefficient=attenuator_coefficient, attenuator_error=attenuator_error,
                         output_workspace=beam_intensity)

    # Divide the sample intensity by the empty beam intensity
    RebinToWorkspace(WorkspaceToRebin=beam_intensity, WorkspaceToMatch=input_workspace, OutputWorkspace=beam_intensity)
    Divide(LHSWorkspace=input_workspace, RHSWorkspace=beam_intensity, OutputWorkspace=output_workspace)

    DeleteWorkspace(Workspace=beam_intensity)  # the temporary workspace is not needed anymore
    return str(output_workspace)

def standard_sample_scaling(input_workspace, f, f_std, output_workspace=None):
    r"""
    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace to be normalized
    f: ~mantid.api.WorkspaceSingleValue
    f_std: mantid.api.SingleValueWorkspace
    output_workspace: ~mantid.api.MatrixWorkspace
    Returns
    -------
        ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    scaling_factor = Divide(LHSWorkspace=f_std, RHSWorkspace=f, OutputWorkspace=uwd())
    output_workspace = Multiply(LHSWorkspace=input_workspace, RHSWorkspace=scaling_factor)
    DeleteWorkspace(Workspace=scaling_factor)
    return output_workspace