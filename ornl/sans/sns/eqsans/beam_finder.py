import numpy as np
from mantid.api import mtd
from mantid.simpleapi import (Integration, FindCenterOfMassPosition,
                              MoveInstrumentComponent)
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import (detector_name, sample_detector_distance)
from ornl.sans.sns.eqsans.mask import apply_mask


__all__ = ['center_detector', 'find_beam_center']


def center_detector(input_workspace, mask=None, method='center_of_mass',
                    x=None, y=None, unit='m', relative=False,
                    move_detector=True, **kwargs):
    r"""
    Move the detector on the XY plane to center the beam location

    An estimation of the center of the detector is carried out
    according to ``method``, unless ``x`` and ``y`` absolute coordinates
    are provided (``relative=False``) in which case ``x`` and ``y`` are used.
    If ``x`` and ``y`` are a translation (``relative=False``) then both
    ``method`` and the translation will be applied.

    case 1: ``center_detector(ws)`` position detector after a center of mass estimation.
    The intersection between the neutron beam and the detector array will have
    coordinates (0, 0, Z).

    case 2: ``center_detector(ws, x=x0, y=y0)`` position detector at absolute
    coordinates (x0,y0)

    case 3: ``enter_detector(ws, method=None, x=x0, y=y0, relative=True)`` translate the
    detector by (x0, y0)

    case 4: ``center_detector(ws, x=x0, y=y0, relative=True)`` position detector after a
    center of mass estimation followed by a translation (x0, y0)

    Parameters
    ----------
    input_workspace: str, ~mantid.api.Workspace
        Input workspace containing the instrument
    mask: str, ``MaskWorkspace``
        Use a mask in conjuction with `method` to find the beam center
    method: str
        Method to estimate the center of the beam. `None` for no method
    x: float
        Final position or translation along the X-axis
    y: float
        Final position or translation along the Y-axis
    unit: str
        units of `x` and `y`. Either meters 'm' or mili-meters 'mm'
    relative: Bool
        Values of `x` and `y` are either absolute coordinates or a
        translation.
    move_detector: bool
        Only calculate the final position if this is False
    kwargs: dict
        Parameters to be passed to FindCenterOfMassPosition

    Returns
    -------
    numpy.ndarray
        Coordinates of the detector center
    """
    ws = mtd[str(input_workspace)]
    i = ws.getInstrument()
    rs = i.getComponentByName(detector_name(i)).getPos()
    rf = np.zeros(3)

    # Use `method` when we are not passing absolute coordinates
    abs_xy = x is not None and y is not None and relative is False
    if method is not None and abs_xy is False:
        method_to_alg = dict(center_of_mass=FindCenterOfMassPosition)
        ws_flattened = Integration(InputWorkspace=ws, OutputWorkspace=uwd())
        if mask is not None:
            fmsk = apply_mask(ws_flattened, mask=mask)
            fmsk.delete()
        # (t_x, t_y) is the intersection point of the neutron beam with the
        # detector
        t_x, t_y = list(method_to_alg[method](InputWorkspace=ws_flattened, **kwargs))
        ws_flattened.delete()
        rs = np.array([-t_x, -t_y, rs[-1]])
        rf = rs

    if x is not None and y is not None:
        t_x = x if unit == 'm' else x / 1.e3
        t_y = y if unit == 'm' else y / 1.e3
        if relative is True:
            rf = rs + np.array([t_x, t_y, 0.0])
        else:
            rf = np.array([t_x, t_y, rs[-1]])

    # Recalculate distance from sample to detector
    if move_detector is True:
        MoveInstrumentComponent(ws, X=rf[0], Y=rf[1], Z=rf[2],
                                ComponentName=detector_name(ws),
                                RelativePosition=False)
        sdd = sample_detector_distance(ws, unit='mm', search_logs=False)
        SampleLogs(ws).insert('sample-detector-distance', sdd, unit='mm')

    return rf


def find_beam_center(input_workspace, method='center_of_mass', mask=None, **kwargs):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    Parameters
    ----------
    input_workspace: str, ~mantid.api.Workspace
    method: str
        Method to calculate the beam center( only 'center_of_mass' is
        implemented)
    mask: str, ``MaskWorkspace``
        Path to mask file, or ``MaskWorkspace`` object
    kwargs: dict
        Parameters to be passed to the method to calculate the center

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center (units in meters)
    """
    # `r` below is the position that the center of the detector would adopt
    # such that the coordinates of the neutron beam impinging on the
    # detector has coordinates (0, 0, z)
    if method != 'center_of_mass':
        raise NotImplementedError(f'{method} is not implemented')
    r = center_detector(input_workspace, mask=mask, method=method,
                        move_detector=False, **kwargs)
    return -r[0], -r[1]
