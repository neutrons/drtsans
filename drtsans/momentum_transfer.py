import numpy as np

from mantid.simpleapi import mtd
from drtsans.dataobjects import IQazimuthal, IQcrystal, IQmod, getDataType, DataType
from drtsans.geometry import pixel_centers
from drtsans.settings import namedtuplefy, unpack_v3d
from drtsans.detector import Detector

__all__ = ['convert_to_q', 'convert_to_subpixels']


def convert_to_subpixels(i_q_azimuthal, workspace, mode, Nx=2, Ny=2, detector='detector1'):
    r"""
    Creates the subpixels arrays for the intensity, error, momentum transfer and resolution

    Parameters
    ----------
    i_q_azimuthal: ~drtsans.dataobjects.IQazimuthal
        The named tuple containing the data and momentum transfer in azimuthal mode
    workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        workspace from which contains information about the number of detector pixels
    mode: str
        Available options are 'scalar', 'azimuthal'
    Nx: int
        number of subpixels in the horizontal direction
    Ny: int
        number of subpixels in the vertical direction
    detector: str
        either 'detector1' or 'wing_detector' - used to get the number of pixels

    Returns
    -------
    ~collections.namedtuple
        either ~IQmod or ~IQazimuthal
    """
    if getDataType(i_q_azimuthal) != DataType.IQ_AZIMUTHAL:
        raise RuntimeError("This method requires an input of IQazimuthal type")
    if mode not in ['scalar', 'azimuthal']:
        raise ValueError("mode must be one of 'scalar', 'azimuthal'")
    ws = mtd[str(workspace)]
    detector_information = Detector(ws, detector)
    # assume vertical tubes
    x_pixels = detector_information.n_tubes
    y_pixels = detector_information.n_pixels_per_tube
    # get the step in qx and qy
    qx_max = np.max(i_q_azimuthal.qx)
    qx_min = np.min(i_q_azimuthal.qx)
    d_qx = (qx_max - qx_min)/((x_pixels - 1) * Nx)
    qy_max = np.max(i_q_azimuthal.qy)
    qy_min = np.min(i_q_azimuthal.qy)
    d_qy = (qy_max - qy_min)/((y_pixels - 1) * Ny)

    # create the arrays
    new_qx = []
    new_qy = []
    # only Qx and Qy change
    for i in range(Nx):
        for j in range(Ny):
            new_qx.append(i_q_azimuthal.qx+(i - (Nx - 1) / 2) * d_qx)
            new_qy.append(i_q_azimuthal.qy+(j - (Ny - 1) / 2) * d_qy)
    new_qx = np.concatenate(new_qx)
    new_qy = np.concatenate(new_qy)
    # the rest are just Nx*Ny copies of the original
    new_intensity = np.concatenate([i_q_azimuthal.intensity] * (Nx * Ny))
    new_error = np.concatenate([i_q_azimuthal.error] * (Nx * Ny))
    if i_q_azimuthal.delta_qx is not None:
        new_delta_qx = np.concatenate([i_q_azimuthal.delta_qx] * (Nx * Ny))
    else:
        new_delta_qx = None
    if i_q_azimuthal.delta_qy is not None:
        new_delta_qy = np.concatenate([i_q_azimuthal.delta_qy] * (Nx * Ny))
    else:
        new_delta_qy = None
    if i_q_azimuthal.wavelength is not None:
        new_wavelength = np.concatenate([i_q_azimuthal.wavelength] * (Nx * Ny))
    else:
        new_wavelength = None

    # create the output
    if mode == 'azimuthal':
        iq = IQazimuthal(new_intensity, new_error, new_qx, new_qy,
                         delta_qx=new_delta_qx, delta_qy=new_delta_qy,
                         wavelength=new_wavelength)
    else:
        # transform to |Q|
        new_q = np.sqrt(new_qx**2 + new_qy**2)
        if new_delta_qx is not None and new_delta_qy is not None:
            new_dq = np.sqrt(new_delta_qx**2 + new_delta_qy**2)
        else:
            new_dq = None
        iq = IQmod(new_intensity, new_error, new_q, delta_mod_q=new_dq, wavelength=new_wavelength)
    return iq


def convert_to_q(ws, mode, resolution_function=None, **kwargs):
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

       Q_x=-\frac{4\pi}{\lambda}\sin\theta\cos\phi

       Q_y=\frac{4\pi}{\lambda}\sin\theta\sin\phi

    - 'crystallographic' mode:

    .. math::

       Q_x=\frac{2\pi}{\lambda}\sin(2\theta)\cos\phi

       Q_y=\frac{2\pi}{\lambda}\sin(2\theta)\sin\phi

       Qz_=\frac{2\pi}{\lambda}(\cos(2\theta)-1)

    Note the minus sign in :math:`Q_x` in the azimuthal mode, so it increases
    to the right when looking at the detector.

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

    # check that the workspace is in units of wavelength
    if not ws:
        raise RuntimeError('Workspace cannot be None')
    wsh = mtd[str(ws)]
    if wsh.getAxis(0).getUnit().unitID() != 'Wavelength':
        raise RuntimeError('Input workspace {} for calculate Q and resolution must be in unit Wavelength but not {}'
                           ''.format(wsh, wsh.getAxis(0).getUnit().unitID()))

    # switch according to mode
    if mode == 'scalar':
        return _convert_to_q_scalar(wsh, resolution_function, **kwargs)
    if mode == 'azimuthal':
        return _convert_to_q_azimuthal(wsh, resolution_function, **kwargs)
    if mode == 'crystallographic':
        return _convert_to_q_crystal(wsh, resolution_function, **kwargs)
    raise NotImplementedError('The mode you selected is not yet implemented')


def _convert_to_q_scalar(ws, resolution_function, **kwargs):
    r"""
    Convert to scalar momentum transfer

    **Mantid algorithms used:**
        :ref:`ConvertUnits <algm-ConvertUnits-v1>`,
        :ref:`ConvertToPointData <algm-ConvertToPointData-v1>`,
        :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`

    Parameters
    ----------

    ws:  ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Workspace in units of wavelength
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
      - mod_q
      - delta_q
      - wavelength

    """
    # get the data
    lam = ws.extractX()
    delta_lam = lam[:, 1:] - lam[:, :-1]
    lam = (lam[:, 1:] + lam[:, :-1]) * 0.5  # wavelengths at the bin centers
    intensity = ws.extractY()
    error = ws.extractE()

    # get geometry info from the original workspace for resolution
    info = pixel_info(ws)  # polar coordinates for each pixel
    number_of_bins = lam.shape[1]
    two_theta = np.repeat(info.two_theta, number_of_bins).reshape(-1, number_of_bins)
    mod_q = 4. * np.pi * np.sin(two_theta * 0.5) / lam

    # calculate the  resolution for each pixel
    if resolution_function is not None:
        delta_q = resolution_function(mod_q, mode='scalar',
                                      pixel_info=info,
                                      wavelength=lam,
                                      delta_wavelength=delta_lam,
                                      **kwargs)
    else:
        delta_q = mod_q * 0.0

    # Has the user requested subpixels?
    keep = info.keep.astype(bool)  # valid spectra indexes (unmasked, not monitor)
    n_horizontal, n_vertical = kwargs.get('n_horizontal', 1), kwargs.get('n_vertical', 1)
    [lam, ] = _filter_and_replicate([lam, ], keep, n_horizontal, n_vertical)
    if n_horizontal * n_vertical > 1:
        # Calculate modulus Q for each subpixel
        subpixel_polar_coords = subpixel_info(ws, n_horizontal, n_vertical)
        two_theta = np.repeat(subpixel_polar_coords.two_theta, number_of_bins).reshape(-1, number_of_bins)
        mod_q = (4. * np.pi * np.sin(two_theta * 0.5) / lam.reshape(len(two_theta), number_of_bins)).ravel()
    else:
        # retain only those pixels that are unmasked or not monitor
        [mod_q, ] = _filter_and_replicate([mod_q, ], keep, n_horizontal=1, n_vertical=1)

    # retain only those pixels that are unmasked or not monitor
    # if user requested subpixels, then replicate pixel quantities for each subpixel
    intensity, error, delta_q = _filter_and_replicate([intensity, error, delta_q], keep, n_horizontal, n_vertical)

    return IQmod(intensity=intensity, error=error, mod_q=mod_q, delta_mod_q=delta_q, wavelength=lam)


def _convert_to_q_azimuthal(ws, resolution_function, **kwargs):
    r"""
    Convert to 2D momentum transfer in azimuthal convention

    **Mantid algorithms used:**
        :ref:`ConvertUnits <algm-ConvertUnits-v1>`,
        :ref:`ConvertToPointData <algm-ConvertToPointData-v1>`,
        :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`

    Parameters
    ----------

    ws:  ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Workspace in units of wavelength
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
      - qx
      - qy
      - delta_qx
      - delta_qy
      - wavelength

    """
    # get the data
    lam = ws.extractX()
    delta_lam = lam[:, 1:] - lam[:, :-1]
    lam = (lam[:, 1:] + lam[:, :-1]) * 0.5
    intensity = ws.extractY()
    error = ws.extractE()

    # get geometry info from the original workspace for resolution
    info = pixel_info(ws)
    number_of_bins = lam.shape[1]
    two_theta = np.repeat(info.two_theta, number_of_bins).reshape(-1, number_of_bins)
    mod_q = 4. * np.pi * np.sin(two_theta * 0.5) / lam
    azimuthal = np.repeat(info.azimuthal, number_of_bins).reshape(-1, number_of_bins)
    qx = -mod_q * np.cos(azimuthal)  # note the convention for the left handed reference frame
    qy = mod_q * np.sin(azimuthal)

    # calculate the resolution for pixels
    if resolution_function is not None:
        delta_qx, delta_qy = resolution_function(qx, qy, mode='azimuthal',
                                                 pixel_info=info,
                                                 wavelength=lam,
                                                 delta_wavelength=delta_lam,
                                                 **kwargs)
    else:
        delta_qx = mod_q * 0.0
        delta_qy = delta_qx

    # Has the user requested subpixels?
    keep = info.keep.astype(bool)  # valid spectra indexes (unmasked, not monitor)
    n_horizontal, n_vertical = kwargs.get('n_horizontal', 1), kwargs.get('n_vertical', 1)
    [lam, ] = _filter_and_replicate([lam, ], keep, n_horizontal, n_vertical)
    if n_horizontal * n_vertical > 1:
        # Calculate modulus Q for each subpixel
        subpixel_polar_coords = subpixel_info(ws, n_horizontal, n_vertical)
        two_theta = np.repeat(subpixel_polar_coords.two_theta, number_of_bins).reshape(-1, number_of_bins)
        mod_q = (4. * np.pi * np.sin(two_theta * 0.5) / lam.reshape(len(two_theta), number_of_bins)).ravel()
        azimuthal = np.repeat(subpixel_polar_coords.azimuthal, number_of_bins)
        qx = -mod_q * np.cos(azimuthal)  # note the convention for the left handed reference frame
        qy = mod_q * np.sin(azimuthal)
    else:
        # retain only those pixels that are unmasked or not monitor
        qx, qy = _filter_and_replicate([qx, qy], keep, n_horizontal=1, n_vertical=1)

    # retain only those pixels that are unmasked or not monitor
    # if user requested subpixels, then replicate pixel quantities for each subpixel
    intensity, error, delta_qx, delta_qy = _filter_and_replicate([intensity, error, delta_qx, delta_qy],
                                                                  keep, n_horizontal, n_vertical)

    return IQazimuthal(intensity=intensity, error=error, qx=qx, qy=qy, delta_qx=delta_qx, delta_qy=delta_qy,
                       wavelength=lam)


def _convert_to_q_crystal(ws, resolution_function, **kwargs):
    r"""
    Convert to 3D momentum transfer in crystallographic convention

    **Mantid algorithms used:**
        :ref:`ConvertUnits <algm-ConvertUnits-v1>`,
        :ref:`ConvertToPointData <algm-ConvertToPointData-v1>`,
        :ref:`DeleteWorkspace <algm-DeleteWorkspace-v1>`

    Parameters
    ----------

    ws:  ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Workspace in units of wavelength
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
      - qx
      - qy
      - qz
      - delta_qx
      - delta_qy
      - delta_qz
      - wavelength

    """
    # get the data
    lam = ws.extractX()
    delta_lam = lam[:, 1:] - lam[:, :-1]
    lam = (lam[:, 1:] + lam[:, :-1]) * 0.5
    intensity = ws.extractY()
    error = ws.extractE()

    # get geometry info from the original workspace for resolution
    info = pixel_info(ws)
    number_of_bins = lam.shape[1]
    two_theta = np.repeat(info.two_theta, number_of_bins).reshape(-1, number_of_bins)
    azimuthal = np.repeat(info.azimuthal, number_of_bins).reshape(-1, number_of_bins)
    temp = 2. * np.pi / lam
    qx = temp * np.sin(two_theta) * np.cos(azimuthal)
    qy = temp * np.sin(two_theta) * np.sin(azimuthal)
    qz = temp * (np.cos(two_theta) - 1.)

    # calculate the  resolution for each pixel
    if resolution_function is not None:
        delta_qx, delta_qy, delta_qz = resolution_function(qx, qy, qz, mode='crystallographic',
                                                           pixel_info=info,
                                                           wavelength=lam,
                                                           delta_wavelength=delta_lam,
                                                           **kwargs)
    else:
        delta_qx = lam * 0.0
        delta_qy = delta_qx
        delta_qz = delta_qx

    # Has the user requested subpixels?
    keep = info.keep.astype(bool)  # valid spectra indexes (unmasked, not monitor)
    n_horizontal, n_vertical = kwargs.get('n_horizontal', 1), kwargs.get('n_vertical', 1)
    [lam, ] = _filter_and_replicate([lam, ], keep, n_horizontal, n_vertical)
    if n_horizontal * n_vertical > 1:
        # Calculate modulus Q for each subpixel
        subpixel_polar_coords = subpixel_info(ws, n_horizontal, n_vertical)
        two_theta = np.repeat(subpixel_polar_coords.two_theta, number_of_bins).reshape(-1, number_of_bins)
        azimuthal = np.repeat(subpixel_polar_coords.azimuthal, number_of_bins).reshape(-1, number_of_bins)
        temp = 2. * np.pi / lam.reshape(len(two_theta), number_of_bins)
        qx = (temp * np.sin(two_theta) * np.cos(azimuthal)).ravel()
        qy = (temp * np.sin(two_theta) * np.sin(azimuthal)).ravel()
        qz = (temp * (np.cos(two_theta) - 1.)).ravel()
    else:
        # retain only those pixels that are unmasked or not monitor
        [qx, qy, qz] = _filter_and_replicate([qx, qy, qz], keep, n_horizontal=1, n_vertical=1)

    # retain only those pixels that are unmasked or not monitor
    # if user requested subpixels, then replicate pixel quantities for each subpixel
    intensity, error = _filter_and_replicate([intensity, error], keep, n_horizontal, n_vertical)
    delta_qx, delta_qy, delta_qz = _filter_and_replicate([delta_qx, delta_qy, delta_qz],
                                                         keep, n_horizontal, n_vertical)

    return IQcrystal(intensity=intensity, error=error, qx=qx, qy=qy, qz=qz,
                     delta_qx=delta_qx, delta_qy=delta_qy, delta_qz=delta_qz, wavelength=lam)


def custom_pixel_info(input_workspace, n_horizontal=1, n_vertical=1):
    if n_horizontal == 1 and n_vertical == 1:
        return pixel_info(input_workspace)  # takes less than 1 second for GPSANS detector
    else:
        return subpixel_info(input_workspace, n_horizontal, n_vertical)  # takes about 13 seconds for GPSANS detector


def _masked_or_monitor(spec_info, idx):
    r"""
    Helper function to check if a spectra is valid

    Parameters
    ----------
    spec_info: ~mantid.api.SpectrumInfo
        SpectrumInfo from a workspace
    idx: int
        index

    Returns
    -------
    bool
        True if spectrum has no detectors, the detector is a monitor, or the spectrum is masked
        False otherwise
    """
    return spec_info.isMonitor(idx) or spec_info.isMasked(idx) or not spec_info.hasDetectors(idx)


@namedtuplefy
def pixel_info(input_workspace):
    r"""
    Helper function to extract two theta angle, azimuthal angle, l2, and a "keep" flag for unmasked pixel detectors.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~maFntid.api.MatrixWorkspace
        Name or reference to a Mantid workspace

    Returns
    -------
    ~collections.namedtuple
        A namedtuple with fields for two_theta, azimuthal, l2, keep
    """
    ws = mtd[str(input_workspace)]
    spec_info = ws.spectrumInfo()

    info = [[np.nan, np.nan, np.nan, False] if _masked_or_monitor(spec_info, i) else
            [spec_info.twoTheta(i), spec_info.azimuthal(i), spec_info.l2(i), True]
            for i in range(ws.getNumberHistograms())]
    info = np.array(info)
    return dict(two_theta=info[:, 0], azimuthal=info[:, 1], l2=info[:, 2], keep=info[:, 3])


@namedtuplefy
def subpixel_info(input_workspace, n_horizontal, n_vertical):
    r"""
    Calculate the two theta angle, azimuthal angle, l2 for the subpixels of each unmasked pixel detector.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace
        Name or reference to a Mantid workspace
    n_horizontal: int
        Number of subpixels along the horizontal direction (on the XZ plane)
    n_vertical: int
        Number of subpixels along the vertical direction (along the Y axis)

    Returns
    -------
    ~collections.namedtuple
        A namedtuple with the following fields:
        two_theta, numpy.ndarray of scattering angles for each subpixel of each valid pixel detector.
        azimuthal, numpy.ndarray of angles on the XY plane  for each subpixel of each valid pixel detector.
        l2, numpy.ndarray for the distance between the sample and each subpixel of each valid pixel detector.
        keep, numpy.ndarray for the workspace indexes of the valid pixel detectors.
    """
    workspace = mtd[str(input_workspace)]
    number_spectra = workspace.getNumberHistograms()
    spectrum_info = workspace.spectrumInfo()
    component_info = workspace.componentInfo()

    # Find valid workspace indexes
    def valid_index(idx):
        return not(spectrum_info.isMonitor(idx) or spectrum_info.isMasked(idx) or not spectrum_info.hasDetectors(idx))
    valid_indexes = [idx for idx in range(number_spectra) if valid_index(idx) is True]

    # Find the componentInfo indexes corresponding to the valid workspace indexes
    get_spectrum_definition = spectrum_info.getSpectrumDefinition
    info_indexes = [get_spectrum_definition(idx)[0][0] for idx in valid_indexes]

    # Find the position of the pixel centers
    pixel_positions = pixel_centers(input_workspace, info_indexes)
    sample_position = component_info.samplePosition()
    pixel_positions -= sample_position

    # Fractional coordinates of subpixels are the same for all pixels.
    x_fractional = np.linspace(0.5, -0.5, n_horizontal, endpoint=False) - 1 / (2 * n_horizontal)
    y_fractional = np.linspace(-0.5, 0.5, n_vertical, endpoint=False) + 1 / (2 * n_vertical)
    z_fractional = np.array([0])
    # xyz_fractional.shape = (number_subpixels, 3)
    xyz_fractional = np.array(np.meshgrid(x_fractional, y_fractional, z_fractional)).T.reshape(-1, 3)

    # Find the dimensions of each pixel. Due to calibrations, each pixel will have its own size
    last_info_index = info_indexes[-1]
    nominal_pixel_dimensions = np.array(component_info.shape(last_info_index).getBoundingBox().width())
    scale_factors = np.array([unpack_v3d(component_info.scaleFactor, i) for i in info_indexes])
    pixel_dimensions = scale_factors * nominal_pixel_dimensions  # shape = (number_pixels, 3)

    # Internal coordinates of subpixels are different for each pixel, because each pixel has its own size
    # xyz_internal.shape = (number_pixels, number_subpixels, 3)
    xyz_internal = pixel_dimensions[:, None, :] * xyz_fractional[None, :, :]

    # We obtain coordinates of subpixels in the laboratory frame of reference
    #  by translating each subpixel by the pixel center coordinates
    xyz = xyz_internal + pixel_positions[:, None, :]  # shape = (number_pixels, number_subpixels, 3)
    xyz = xyz.reshape((len(info_indexes) * n_horizontal * n_vertical, 3))

    l2 = np.linalg.norm(xyz, axis=1)  # shape = (number_subpixels * number_pixels, )
    two_theta = np.arccos(xyz[:, 2] / l2)  # cos(two_theta) = z / l2
    azimuthal = np.arctan2(xyz[:, 1], xyz[:, 0])  # numpy.arctan2(y, x) is quadrant-aware

    return dict(two_theta=two_theta, azimuthal=azimuthal, l2=l2, keep=valid_indexes)


def _filter_and_replicate(arrays, unmasked_indexes, n_horizontal=1, n_vertical=1):
    r"""
    Retain only items in the array corresponding to unmasked pixels, and replicate
    these values for each subpixel if subpixels are requested.

    For every array in ```arrays```, the same filter and replication steps are applied.
    Parameters
    ----------
    arrays: list
        List of arrays to be filtered.
    unmasked_indexes: numpy.ndarray
        List of array indexes corresponding to spectra with unmasked pixels.
    n_horizontal: int
        Number of subpixels along the horizontal direction (on the XZ plane)
    n_vertical: int
        Number of subpixels along the vertical direction (along the Y axis)

    Returns
    -------
    list
        List of filtered and replicated arrays
    """
    # It's assumed that arrays are of the shape (...,1) so that reshape(-1) will eliminate the last dimension
    return [np.repeat(array[unmasked_indexes, :], n_horizontal * n_vertical, axis=0).reshape(-1) for array in arrays]
