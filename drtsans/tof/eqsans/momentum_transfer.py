import numpy as np
from collections import namedtuple
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/momentum_transfer.py
import drtsans.momentum_transfer
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/resolution.py
import drtsans.resolution
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/geometry.py
from drtsans.tof.eqsans.geometry import (sample_aperture_diameter, source_aperture_diameter,
                                         source_aperture_sample_distance)
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/geometry.py
from drtsans import geometry as sans_geometry
from drtsans.geometry import pixel_size
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py
from drtsans.samplelogs import SampleLogs
from mantid.kernel import logger

__all__ = ['convert_to_q', 'split_by_frame']


def eqsans_resolution(*args, **kwargs):
    r"""
    Function to compute the resolution for EQSANS. Some parts are calculated in drtsans.resolution

    Parameters
    ----------
    args:
        Either |Q| or Qx and Qy
    kwargs: dict
        The following parameters are required

        - mode: one of 'scalar', 'azimuthal'
        - instrument_parameters: an InstrumentSetupParameters object containg information about L1, apertures, etc.
        - pixel_info: dict containing arrays for two_theta, azimuthal, L2 (per pixel)
        - wavelength: array of wavenegths (same shape as Q)
        - delta_wavelength: array of wavenegth widths (same shape as Q)

    Returns
    -------
    ~np.array or list of arrays
        Resolution along the components of momentum transfer from the input
    """
    mode = kwargs.get('mode')
    instrument_parameters = kwargs.get('instrument_parameters')
    pixel_info = kwargs.get('pixel_info')
    wavelength = kwargs.get('wavelength')
    delta_wavelength = kwargs.get('delta_wavelength')

    L1 = instrument_parameters.l1
    samp_det_distance = pixel_info.l2.reshape(-1, 1)

    # get the general sigma_geom^2 (not mode/ instrument dependent)
    sigma_geom = drtsans.resolution.calculate_sigma_geometry(mode,
                                                             wavelength,
                                                             delta_wavelength,
                                                             pixel_info,
                                                             instrument_parameters)
    # get the moderator part of the resolution (only EQSANS)
    moderator_part = (3.9650 * 0.001 * moderator_time_uncertainty(wavelength) /
                      (wavelength*(L1+samp_det_distance)))**2

    # return resolution according to formulas 10.5 and 10.6 in the master document
    if mode == 'scalar':
        q = args[0]
        return np.sqrt(sigma_geom + np.square(q) * (np.square(delta_wavelength / wavelength) + moderator_part) / 12.)
    if mode == 'azimuthal':
        qx = args[0]
        qy = args[1]
        return [np.sqrt(sigma_geom[0] + np.square(qx) * (np.square(delta_wavelength / wavelength)
                        + moderator_part) / 12.),
                np.sqrt(sigma_geom[1] + np.square(qy) * (np.square(delta_wavelength / wavelength)
                        + moderator_part) / 12.)]

    # should not get here
    raise NotImplementedError('The mode you selected is not yet implemented')


def convert_to_q(ws, mode, resolution_function=eqsans_resolution, **kwargs):
    r"""
    Convert a workspace with units of wavelength into a
    series of arrays: intensity, error, q (or q components),
    delta q (or delta q components), and wavelength

    Using the scattering angle as :math:`2\theta` and azimuthan angle as
    :math:`\phi`,the calculaion of momentum transfer is:

    - 'scalar' mode:

    .. math:: |Q| = \frac{4\pi}{\lambda}\sin\theta

    - 'azimuthal' mode:

    .. math::

       Q_x=\frac{4\pi}{\lambda}\sin\theta\cos\phi

       Q_y=\frac{4\pi}{\lambda}\sin\theta\sin\phi

    - 'crystallographic' mode:

    .. math::

       Q_x=\frac{2\pi}{\lambda}\sin(2\theta)\cos\phi

       Q_y=\frac{2\pi}{\lambda}\sin(2\theta)\sin\phi

       Qz_=\frac{2\pi}{\lambda}(\cos(2\theta)-1)


    It calls drtsans.momentum_transfer.convert_to_q

    Parameters
    ----------

    ws:  str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Workspace in units of wavelength
    mode: str
        Available options are 'scalar', 'azimuthal', and 'crystallographic'
    resolution_function:
        Function to calculate resolution
    kwargs:
        Parameters to be passed to the resolution function

    Returns
    -------
    ~collections.namedtuple
       A namedtuple with fields for

      - intensity
      - error
      - mod_q (:math:`|Q|`) or qx, qy (:math:`Q_x, Q_y`) or qx, qy, qz (:math:`Q_x, Q_y, Q_z`) (depending on the mode)
      - delta_q or delta_qx, delta_qy or delta_qx, delta_qy, delta_qz - the resolution along the q components
      - wavelength

    """
    # check if one wants to override pixel sizes
    pixel_sizes = kwargs.get('pixel_sizes', None)
    # get the InstrumentSetupParameters
    instrument_setup = retrieve_instrument_setup(ws, pixel_sizes)
    return drtsans.momentum_transfer.convert_to_q(ws, mode, resolution_function,
                                                  instrument_parameters=instrument_setup, **kwargs)


def retrieve_instrument_setup(ws):
    """ Get instrument parameter including L1, L2, source aperture diameter and sample aperture radius
    :param ws:
    :return: MomentumTransferResolutionParameters instance
    """
    # Retrieve L1 and L2 from instrument geometry
    l1 = source_aperture_sample_distance(ws, unit='m')
    l2 = sans_geometry.sample_detector_distance(ws, unit='m',
                                                search_logs=False)
    r1 = 0.5 * source_aperture_diameter(ws, unit='m')
    r2 = 0.5 * sample_aperture_diameter(ws, unit='m')

    # Set up the parameter class
    setup_params = drtsans.resolution.InstrumentSetupParameters(l1=l1,
                                                                sample_det_center_dist=l2,
                                                                source_aperture_radius=r1,
                                                                sample_aperture_radius=r2)
    return setup_params


def moderator_time_uncertainty(wl):
    """ Relative Q uncertainty due to emission time jitter
    in the neutron moderator.
    :param wl: float (or ndarray) wavelength [Angstrom]
    :return: float or ~np.array (same shape to wave_length_array) of emission error time
    """
    original_float = False
    if isinstance(wl, float):
        wl = np.array([wl])
        original_float = True
    # init output array to zeros
    time_error = np.zeros_like(wl)

    # formula for lambda>2
    mask = wl > 2.0
    time_error[mask] = 0.0148 * wl[mask]**3 \
        - 0.5233 * wl[mask]**2 \
        + 6.4797 * wl[mask] + 231.99

    # formula for lambda<=2
    mask = wl <= 2.0
    time_error[mask] = 392.31 * wl[mask]**6 \
        - 3169.3 * wl[mask]**5 \
        + 10445 * wl[mask]**4 \
        - 17872 * wl[mask]**3 \
        + 16509 * wl[mask]**2 \
        - 7448.4 * wl[mask] + 1280.5

    # clean up memory
    del mask

    # convert from zero-dimensional nparray to float
    if original_float:
        time_error = float(time_error)

    return time_error


def split_by_frame(input_workspace, *args, **kwargs):
    r""" Split the output of convert_to_q into a list of
    outputs, one for each frame.

    Parameters
    ----------

    input_workspace:  str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Workspace in units of wavelength
    kwargs: ~collections.namedtuple or dict
        Outputs from convert_to_q

    Returns
    -------

    list
        A list with namedtuples
    """
    # transform args to dictionary
    try:
        kwargs = args[0]._asdict()
    except AttributeError:
        pass
    keys = kwargs.keys()
    # get the wavelengths for each frame
    sl = SampleLogs(input_workspace)
    frames = []
    if bool(sl.is_frame_skipping.value) is True:
        logger.information("This is a frame skipping data set.")
        frame1_wavelength_min = sl.wavelength_skip_min.value
        frame1_wavelength_max = sl.wavelength_skip_max.value
        frame2_wavelength_min = sl.wavelength_lead_min.value
        frame2_wavelength_max = sl.wavelength_lead_max.value
        frames.append((frame1_wavelength_min, frame1_wavelength_max))
        frames.append((frame2_wavelength_min, frame2_wavelength_max))
    else:
        frame1_wavelength_min = sl.wavelength_min.value
        frame1_wavelength_max = sl.wavelength_max.value
        frames.append((frame1_wavelength_min, frame1_wavelength_max))

    # create the output list
    output_list = []
    for wl_min, wl_max in frames:
        output = dict()
        wavelength = kwargs['wavelength']
        # use only data between the given wavelengths
        keep = np.logical_and(np.greater_equal(wavelength, wl_min), np.less_equal(wavelength, wl_max))
        for k in keys:
            output[k] = kwargs[k][keep]
        # create the namedtuple from a dictionary
        output_list.append(namedtuple('frame', output)(**output))
    return output_list
