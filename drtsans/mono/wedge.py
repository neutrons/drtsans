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
        ```(qmod, phi)``` with the same dimensionality as the data.intensity'''
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

    return q, phi
