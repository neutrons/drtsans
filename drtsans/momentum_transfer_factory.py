from drtsans.tof.eqsans.momentum_transfer import calculate_q_dq as tof_cal_q_dq
from drtsans.mono.momentum_transfer import calculate_q_dq as mono_cal_q_dq
from mantid.api import AnalysisDataService


# List of instrument names (from mantid instrument name) for HFIR/mono/ SANS beamlines
MonoSANS = ['GP-SANS', 'BIO-SANS']
# List of instrument names (from mantid instrument name) for SNS/TOF SANS beamlines
TofSANS = ['EQ-SANS']


def calculate_q_dq(ws, pixel_sizes=None, instrument_type=None):
    """Calculate momentum transfer and momentum transfer resolution

    Parameters
    ----------
    ws  : list of workspaces (references or names)
    pixel_sizes : dictionary for detector pixel sizes
    instrument_type: extra flag ("mono" or "tof") for testing purpose (Generic workspace)

    Returns
    -------
        2D arrays for Q, Qx, dQx, Qy, dQy
    """
    if not isinstance(ws, list):
        raise RuntimeError('Calculate Q dQ only accept workspace(s) as list but not {}'.format(type(ws)))

    q_list = list()

    for ws_i in ws:
        # Identify the instrument
        ws_i = AnalysisDataService.retrieve(str(ws_i))

        if instrument_type == 'tof' or ws_i.getInstrument().getName() == TofSANS:
            # TOF workspace (EQ-SANS or Generic)
            q_set = tof_cal_q_dq(ws, pixel_sizes)
        elif instrument_type == 'mono' or ws_i.getInstrument.getName() in MonoSANS:
            # HFIR mono wavelength reactor
            q_set = mono_cal_q_dq(ws, pixel_sizes)
        else:
            # Unsupported case
            raise RuntimeError('Instrument {} is not supported without instrument type ({})is not specified'
                               'or not supported'
                               ''.format(ws_i.getInstrument().getName(), instrument_type))

        q_list.append(q_set)
    # END-FOR

    return q_list

