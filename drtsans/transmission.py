import numpy as np

from mantid import mtd
from mantid.kernel import logger
# https://docs.mantidproject.org/nightly/algorithms/ApplyTransmissionCorrection-v1.html
# https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html
# https://docs.mantidproject.org/nightly/algorithms/GroupDetectors-v1.html
# https://docs.mantidproject.org/nightly/algorithms/RebinToWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/ReplaceSpecialValues-v1.html
from mantid.simpleapi import (ApplyTransmissionCorrection, Divide, GroupDetectors, ReplaceSpecialValues,
                              RebinToWorkspace)

from drtsans.mono.geometry import beam_radius
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.geometry import masked_detectors
from drtsans.mask_utils import circular_mask_from_beam_center

# Symbols to be exported to the drtsans namespace
__all__ = ['apply_transmission_correction']


def calculate_transmission(input_sample, input_reference,
                           radius=None, radius_unit='mm',
                           output_workspace=None):
    """
    Calculate the raw transmission coefficients at zero scattering angle
    from already prepared sample and reference data.

    For EQ-SANS, one additional step fitting the returned raw values is
    necessary. Use `eqsans.calculate_transmission` instead.

    Parameters
    ----------
    input_sample: MatrixWorkspace
        Prepared sample workspace (possibly obtained with an attenuated beam)
    input_reference: MatrixWorkspace
        Prepared direct beam workspace (possibly obtained with an attenuated
         beam)
    radius: float
        Radius around the bean center for pixel integration, in millimeters.
        If None, radius will be obtained or calculated using `input_reference`.
    radius_unit: str
        Either 'mm' or 'm', and only used in conjunction with option `radius`.
    output_workspace: str
        Name of the output workspace containing the raw transmission values.
        If None, a hidden random name will be provided.

    Returns
    -------
    MatrixWorkspace
        Workspace containing the raw transmission values
    """
    if output_workspace is None:
        output_workspace = uwd()

    if radius is None:
        logger.information('Calculating beam radius from sample logs')
        radius = beam_radius(input_reference, unit='mm')
    else:

        radius = float(radius) if radius_unit == 'mm' else 1.e3 * radius  # to mm
    if radius <= 0.:
        raise ValueError('Encountered negative beam radius={}mm'.format(radius))

    det_ids = circular_mask_from_beam_center(input_reference, radius, unit='mm')
    if not det_ids:
        raise RuntimeError('No pixels in beam with radius of {:.2f} mm'.format(radius))

    # Avoid pitfall of masking many pixels around the beam center
    msg = 'More than half of the detectors ' +\
          'within a radius of {:.2f} mm '.format(radius) +\
          'from the beam center are masked in the input {0}'
    for k, v in dict(sample=input_sample, reference=input_reference).items():
        if len(masked_detectors(v, det_ids)) > len(det_ids) / 2:
            raise RuntimeError(msg.format(k))

    # by default it sums all the grouped detectors
    gis = GroupDetectors(InputWorkspace=input_sample, DetectorList=det_ids,
                         OutputWorkspace=uwd())
    gir = GroupDetectors(InputWorkspace=input_reference, DetectorList=det_ids,
                         OutputWorkspace=uwd())
    gir = RebinToWorkspace(WorkspaceToRebin=gir, WorkspaceToMatch=gis,
                           OutputWorkspace=gir.name())

    # calculate zero angle transmission coefficient(s)
    zat = Divide(LHSWorkspace=gis, RHSWorkspace=gir,
                 OutputWorkspace=output_workspace)
    av_t, av_e = np.mean(zat.dataY(0)), np.linalg.norm(zat.dataE(0))
    message = 'Average zero angle transmission = {0} +/- {1}'
    logger.notice(message.format(av_t, av_e))
    gis.delete()
    gir.delete()
    return zat


def apply_transmission_correction(input_workspace, trans_workspace=None, trans_value=None, trans_error=0.0,
                                  theta_dependent=True, output_workspace=None):
    r"""
    Correct the intensities with transmission coefficient(s).

    **Mantid algorithms used:**
    :ref:`ApplyTransmissionCorrection <algm-ApplyTransmissionCorrection-v1>`,

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
    output_workspace: str
        Name of the workspace containing the corrected intensities. If :py:obj:`None`, the `input_workspace`
        will be overwritten.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    kwargs = dict(InputWorkspace=input_workspace, ThetaDependent=theta_dependent, OutputWorkspace=output_workspace)
    if trans_workspace is not None:
        # EQ-SANS transmissions in skip-frame mode have transmission values of zero in the wavelength gap.
        # Need to be replaced with one to avoid division of intensities by zero.
        clean_trans_workspace = ReplaceSpecialValues(InputWorkspace=trans_workspace, SmallNumberThreshold=1.0e-6,
                                                     SmallNumberValue=1.0, OutputWorkspace=uwd())
        kwargs['TransmissionWorkspace'] = clean_trans_workspace
    elif trans_value is not None:
        kwargs.update(dict(TransmissionValue=trans_value, TransmissionError=trans_error))
    else:
        raise RuntimeError('Provide either trans_workspace or trans_value')
    ApplyTransmissionCorrection(**kwargs)
    if trans_workspace is not None:
        clean_trans_workspace.delete()
    return mtd[output_workspace]
