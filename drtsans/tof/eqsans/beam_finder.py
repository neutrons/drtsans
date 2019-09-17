import numpy as np
from mantid.api import mtd
from mantid.simpleapi import Integration, FindCenterOfMassPosition, MoveInstrumentComponent
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import (detector_name, sample_detector_distance)
from drtsans.mask_utils import mask_spectra_with_special_values
from drtsans.tof.eqsans.mask import apply_mask

__all__ = ['center_detector', 'find_beam_center']


def center_detector(input_workspace, mask=None, x=None, y=None, unit='m', relative=False, method='center_of_mass',
                    move_detector=True, **method_kwargs):
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
    x: float
        Final position or translation along the X-axis. Units must be those of option `unit`, which
        defaults to meters.
    y: float
        Final position or translation along the Y-axis. Units must be those of option `unit`, which
        defaults to meters.
    unit: str
        units of `x` and `y`. Either meters 'm' or mili-meters 'mm'. Default is meters
    relative: Bool
        Values of `x` and `y` are either absolute coordinates or a
        translation.
    method: str
        Method to estimate the center of the beam. Default is `center_of_mass` and use `None` for no method.
            Method     -->     Mantid algorithm
        center_of_mass --> FindCenterOfMassPosition
    move_detector: bool
        Only calculate the final position if this is False
    method_kwargs: dict
        Parameters to be passed to the selected method

    Returns
    -------
    numpy.ndarray
        Final position of the detector's center, always in meters.
    """
    if method not in (None, 'center_of_mass'):
        raise RuntimeError('Not implemented method')
    unit_to_meters = dict(m=1., mm=1.e-3)
    method_to_algorithm = dict(center_of_mass=FindCenterOfMassPosition)  # in case we add more methods later
    method_to_algorithm_options = {'center_of_mass': dict(DirectBeam=True)}  # default options of method
    workspace = mtd[str(input_workspace)]
    instrument = workspace.getInstrument()
    starting_position = instrument.getComponentByName(detector_name(instrument)).getPos()
    final_position = np.zeros(3)  # this has to be determined

    # Use `method` is necessary when we are not passing absolute coordinates
    xy_are_absolute_coordinates = x is not None and y is not None and relative is False
    if method is not None and xy_are_absolute_coordinates is False:
        workspace_flattened = Integration(InputWorkspace=workspace, OutputWorkspace=uwd())
        if mask is not None:
            mask_workspace = apply_mask(workspace_flattened, mask=mask)
            mask_workspace.delete()  # we don't need the mask workspace so keep it clean
        mask_spectra_with_special_values(workspace_flattened)
        algorithm = method_to_algorithm[method]
        algorithm_options = method_to_algorithm_options[method]
        algorithm_options.update(method_kwargs)
        # (t_x, t_y) is the intersection point of the neutron beam with the detector
        t_x, t_y = list(algorithm(InputWorkspace=workspace_flattened, **algorithm_options))
        workspace_flattened.delete()
        starting_position = np.array([-t_x, -t_y, starting_position[-1]])
        final_position = starting_position

    if x is not None and y is not None:
        t_x = x * unit_to_meters[unit]
        t_y = y * unit_to_meters[unit]
        if relative is True:
            final_position = starting_position + np.array([t_x, t_y, 0.0])
        else:
            final_position = np.array([t_x, t_y, starting_position[-1]])

    # Recalculate distance from sample to detector
    if move_detector is True:
        final_xyz = dict(zip(('X', 'Y', 'Z'), final_position))
        MoveInstrumentComponent(workspace, ComponentName=detector_name(workspace), RelativePosition=False, **final_xyz)
        sdd = sample_detector_distance(workspace, unit='mm', search_logs=False)
        SampleLogs(workspace).insert('sample-detector-distance', sdd, unit='mm')

    return final_position


def find_beam_center(input_workspace, mask=None, method='center_of_mass', **method_kwargs):
    r"""
    Calculate absolute coordinates of beam impinging on the detector.
    Usually employed for a direct beam run (no sample and not sample holder).

    Parameters
    ----------
    input_workspace: str, ~mantid.api.Workspace
    mask: str, ``MaskWorkspace``
        Path to mask file, or ``MaskWorkspace`` object
    mask: str, MaskWorkspace
        Path to mask file, or MaskWorkspace object
    method_kwargs: dict
        Additional keyword arguments to be passed to the method to calculate the center
            Method     -->     Mantid algorithm
        center_of_mass --> FindCenterOfMassPosition

    Returns
    -------
    tuple
        (X, Y) coordinates of the beam center, always in meters.
    """
    # `r` below is the position that the center of the detector would adopt
    # such that the coordinates of the neutron beam impinging on the
    # detector has coordinates (0, 0, z)
    if method != 'center_of_mass':
        raise NotImplementedError('{} is not implemented'.format(method))
    detector_coordinates = center_detector(input_workspace, mask=mask, method=method, move_detector=False,
                                           **method_kwargs)
    return -detector_coordinates[0], -detector_coordinates[1]
