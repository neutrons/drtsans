from mantid.simpleapi import (MaskBTP, FindCenterOfMassPosition, MoveInstrumentComponent, Integration)
from mantid.kernel import logger
from drtsans import mask_utils
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.mask_utils import mask_spectra_with_special_values


__all__ = ['center_detector', 'find_beam_center']


def find_beam_center(input_ws, method='center_of_mass', mask=None, **kwargs):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    **Mantid algorithms used:**
        :ref:`FindCenterOfMassPosition <algm-FindCenterOfMassPosition-v2>`,
        :ref:`Integration <algm-Integration-v1>`,
        :ref:`MaskBTP <algm-MaskBTP-v1>`

    Parameters
    ----------
    input_workspace: str, Workspace
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)
    mask: str, ``MaskWorkspace``
        Use a mask in conjuction with `method` to find the beam center
    kwargs: dict
        Parameters to be passed to the method to calculate the center or to MaskBTP or 'panel'
        'panel' is either 'front' or 'back' to mask a whole panel
    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters)
    """
    if method != 'center_of_mass':
        raise NotImplementedError()  # (f'{method} is not implemented')

    # flatten the workspace
    flat_ws = Integration(InputWorkspace=input_ws, OutputWorkspace=uwd())
    mask_spectra_with_special_values(flat_ws)

    # parameters to be passed to masking
    mask_dict = {}
    if mask:
        mask_dict['mask'] = mask

    # mask panel for EQSANS
    panel = kwargs.pop('panel', None)
    if panel:
        MaskBTP(Workspace=flat_ws, Components=panel + '-panel')
    for key in ['Components', 'Bank', 'Tube', 'Pixel']:
        value = kwargs.pop(key, None)
        if value is not None:
            mask_dict[key] = value

    if mask_dict:
        mask_workspace = mask_utils.apply_mask(flat_ws, **mask_dict)
        mask_workspace.delete()  # we don't need the mask workspace so keep it clean

    # find center of mass position
    center = FindCenterOfMassPosition(InputWorkspace=flat_ws.name(), **kwargs)
    center_x, center_y = center
    logger.information("Found beam position: X={:.3} m, Y={:.3} m.".format(
        center_x, center_y))
    return center_x, center_y


def center_detector(input_workspace, center_x, center_y, component='detector1'):
    """Center the detector. Move the `component` by (-center_x, -centery)
    from the current position (relative motion).

    **Mantid algorithms used:**
    :ref:`MoveInstrumentComponent <algm-MoveInstrumentComponent-v1>`,

    Parameters
    ----------
    input_workspace : Workspace2D, str
        The workspace to be centered
    center_x : float
        in meters
    center_y : float
        in meters
    component : string
        name of the detector to be centered
    """

    MoveInstrumentComponent(Workspace=input_workspace,
                            ComponentName=component,
                            X=-center_x,
                            Y=-center_y,
                            RelativePosition=True)
