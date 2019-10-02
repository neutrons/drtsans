from drtsans.tof.eqsans.momentum_transfer import calculate_q_dq as tof_cal_q_dq
from drtsans.mono.momentum_transfer import calculate_q_dq as mono_cal_q_dq
from mantid.api import AnalysisDataService
import collections
import numpy as np


# List of instrument names (from mantid instrument name) for HFIR/mono/ SANS beamlines
MonoSANS = ['GP-SANS', 'BIO-SANS']
# List of instrument names (from mantid instrument name) for SNS/TOF SANS beamlines
TofSANS = ['EQ-SANS']

# Define named tuple for Q
# np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
# Q, Qx, dQx, Qy, dQy
SansMomentumTransfer = collections.namedtuple('SansMomentumTransfer', 'q qx qy dq dqx dqy')


def calculate_q_dq(ws, pixel_sizes=None, instrument_type=None):
    """Calculate momentum transfer and momentum transfer resolution

    Parameters
    ----------
    ws  : list of workspaces (references or names)
    pixel_sizes : dictionary for detector pixel sizes
    instrument_type: extra flag ("mono" or "tof") for testing purpose (Generic workspace)

    Returns
    -------
    SansMomentumTransfer
        Momentum transfer set including 1D, 2D and resolutions
    """
    if isinstance(ws, list):
        raise RuntimeError('Calculate Q dQ cannot accept workspace {} ({}) as a list.'
                           ''.format(ws, type(ws)))

    # Identify the instrument
    ws = AnalysisDataService.retrieve(str(ws))

    if instrument_type == 'tof' or ws.getInstrument().getName() in TofSANS:
        # TOF workspace (EQ-SANS or Generic)
        qx_matrix, qy_matrix, dqx_matrix, dqy_matrix = tof_cal_q_dq(ws, pixel_sizes)
    elif instrument_type == 'mono' or ws.getInstrument.getName() in MonoSANS:
        # HFIR mono wavelength reactor
        qx_matrix, qy_matrix, dqx_matrix, dqy_matrix = mono_cal_q_dq(ws, pixel_sizes)
    else:
        # Unsupported case
        raise RuntimeError('Instrument {} is not supported without instrument type ({})is not specified'
                           'or not supported'
                           ''.format(ws.getInstrument().getName(), instrument_type))

    # Calculate Q
    q_array = np.sqrt(qx_matrix**2 + dqy_matrix**2)
    dq_array = np.sqrt(dqx_matrix**2 + dqy_matrix**2)

    q_set = SansMomentumTransfer(q_array, qx_matrix, qy_matrix, dq_array, dqx_matrix, dqy_matrix)

    return q_set
