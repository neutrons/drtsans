from mantid.simpleapi import mtd
from drtsans.dataobjects import IQazimuthal, IQmod
from drtsans.settings import namedtuplefy
import numpy as np


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

       Q_x=\frac{4\pi}{\lambda}\sin\theta\cos\phi

       Q_y=\frac{4\pi}{\lambda}\sin\theta\sin\phi

    - 'crystallographic' mode:

    .. math::

       Q_x=\frac{2\pi}{\lambda}\sin(2\theta)\cos\phi

       Q_y=\frac{2\pi}{\lambda}\sin(2\theta)\sin\phi

       Qz_=\frac{2\pi}{\lambda}(\cos(2\theta)-1)

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
    lam = (lam[:, 1:] + lam[:, :-1]) * 0.5
    intensity = ws.extractY()
    error = ws.extractE()

    # get geometry info from the original workspace for resolution
    info = pixel_info(ws)
    number_of_bins = lam.shape[1]
    two_theta = np.repeat(info.two_theta, number_of_bins).reshape(-1, number_of_bins)
    mod_q = 4. * np.pi * np.sin(two_theta * 0.5) / lam

    # calculate the  resolution
    if resolution_function is not None:
        delta_q = resolution_function(mod_q, mode='scalar',
                                      pixel_info=info,
                                      wavelength=lam,
                                      delta_wavelength=delta_lam,
                                      **kwargs)
    else:
        delta_q = mod_q * 0.0
    # account for masking and monitors
    keep = info.keep.astype(bool)
    lam = lam[keep, :].reshape(-1)
    intensity = intensity[keep, :].reshape(-1)
    error = error[keep, :].reshape(-1)
    mod_q = mod_q[keep, :].reshape(-1)
    delta_q = delta_q[keep, :].reshape(-1)

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
    qx = mod_q * np.cos(azimuthal)
    qy = mod_q * np.sin(azimuthal)

    # calculate the  resolution
    if resolution_function is not None:
        delta_qx, delta_qy = resolution_function(qx, qy, mode='azimuthal',
                                                 pixel_info=info,
                                                 wavelength=lam,
                                                 delta_wavelength=delta_lam,
                                                 **kwargs)
    else:
        delta_qx = mod_q * 0.0
        delta_qy = delta_qx
    # account for masking and monitors
    keep = info.keep.astype(bool)
    lam = lam[keep, :].reshape(-1)
    intensity = intensity[keep, :].reshape(-1)
    error = error[keep, :].reshape(-1)
    qx = qx[keep, :].reshape(-1)
    qy = qy[keep, :].reshape(-1)
    delta_qx = delta_qx[keep, :].reshape(-1)
    delta_qy = delta_qy[keep, :].reshape(-1)

    return IQazimuthal(intensity=intensity, error=error, qx=qx, qy=qy,
                       delta_qx=delta_qx, delta_qy=delta_qy, wavelength=lam)


@namedtuplefy
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

    # calculate the  resolution
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
    # account for masking and monitors
    keep = info.keep.astype(bool)
    lam = lam[keep, :].reshape(-1)
    intensity = intensity[keep, :].reshape(-1)
    error = error[keep, :].reshape(-1)
    qx = qx[keep, :].reshape(-1)
    qy = qy[keep, :].reshape(-1)
    qz = qz[keep, :].reshape(-1)
    delta_qx = delta_qx[keep, :].reshape(-1)
    delta_qy = delta_qy[keep, :].reshape(-1)
    delta_qz = delta_qz[keep, :].reshape(-1)

    return dict(intensity=intensity, error=error, qx=qx, qy=qy, qz=qz,
                delta_qx=delta_qx, delta_qy=delta_qy, delta_qz=delta_qz, wavelength=lam)


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
def pixel_info(ws):
    r"""
    Helper function to extract two theta angle, azimuthal angle, l2, and a "keep" flag.
    The "keep" is false if the spectrum has no detectors, the detector is a monitor,
    or if the spectrum has been masked

    Parameters
    ----------
    ws: ~mantid.api.IEventWorkspace, ~mantid.api.MatrixWorkspace

    Returns
    -------
    ~collections.namedtuple
        A namedtuple with fields for two_theta, azimuthal, l2, keep
    """
    spec_info = ws.spectrumInfo()
    info = [[np.nan, np.nan, np.nan, False] if _masked_or_monitor(spec_info, i) else
            [spec_info.twoTheta(i), spec_info.azimuthal(i), spec_info.l2(i), True]
            for i in range(ws.getNumberHistograms())]
    info = np.array(info)
    return dict(two_theta=info[:, 0], azimuthal=info[:, 1], l2=info[:, 2], keep=info[:, 3])
