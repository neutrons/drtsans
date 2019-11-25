from drtsans.dataobjects import DataType, getDataType
import numpy as np


def _toQmodAndPhi(data):
    '''This function returns the values of qmod and phi that are parallel
    to the original data array. It assumes that the data is IQazimuthal

    Parameters
    ==========
    data: ~drtsans.dataobjects.Azimuthal

    Results
    =======
    tuple
        ```(qmod, phi)``` with the same dimensionality as the data.intensity
        with Q in angstrom and phi in degrees'''
    if not getDataType(data) == DataType.IQ_AZIMUTHAL:
        raise RuntimeError('Calculating qmod and phi only works for IQazimuthal')

    # reshape the qx and qy if intensity array is 2d
    if len(data.intensity.shape) == 2 and len(data.qx.shape) == 1 and len(data.qy.shape) == 1:
        qx = np.tile(data.qx, (data.qy.shape[0], 1))
        qy = np.tile(data.qy, (data.qx.shape[0], 1)).transpose()
    else:
        qx = data.qx
        qy = data.qy

    # calculate q-scalar
    q = np.sqrt(np.square(qx) + np.square(qy))

    # phi is expected to be positive so use cyclical nature of trig functions
    phi = np.arctan2(qy, qx)
    phi[phi < 0.] += 2. * np.pi

    return q, np.rad2deg(phi)

def _binInQAndPhi(data, q_min, q_delta, q_max, phi_delta):
    q, phi = _toQmodAndPhi(data)

    # the bonus two steps is to get the end-point in the array
    q_bins = np.arange(q_min, q_max + 2. * q_delta, q_delta, dtype=float)
    # additional half-circle is to pick up things that are symmetric near phi=0
    phi_bins = np.arange(0., 540. + 2. * phi_delta, phi_delta, dtype=float)

    # create output data
    intensity = np.zeros((phi_bins.size-1, q_bins.size-1), dtype=float)
    error = np.zeros((phi_bins.size-1, q_bins.size-1), dtype=float)

    # do the binning
    for q_val, phi_val, i_val, e_val in zip(q.ravel(), phi.ravel(), data.intensity.ravel(), data.error.ravel()):
        q_index = q_bins.searchsorted(q_val, side='right')
        if q_index >= q_bins.size or q_index == 0:
            print('FAILED TO FIND Q_BIN FOR {} from {}'.format(q_val, q_bins))
            continue

        phi_index = phi_bins.searchsorted(phi_val, side='right')
        if phi_index >= phi_bins.size or q_index == 0:
            print('FAILED TO FIND PHI_BIN FOR {}'.format(phi_val))
            continue

        if i_val == 700.:
            print(i_val, '\t',
                  'Q[{}] {} <= {} < {}'.format(q_index-1, q_bins[q_index-1], q_val, q_bins[q_index]),
                  '\t'
                  'phi[{}] {} <= {} < {}'.format(phi_index-1, phi_bins[phi_index-1], phi_val, phi_bins[phi_index]))

        intensity[phi_index - 1, q_index - 1] += i_val
        error[phi_index - 1, q_index - 1] += e_val

    # bins that didn't accumulate counts are nan
    mask = error == 0.
    intensity[mask] = np.nan
    error[mask] = np.nan

    return intensity, error, phi_bins, q_bins
