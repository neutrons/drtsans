from multiprocessing import Pool

import pytest

import mantid
from mantid import mtd
from mantid.kernel import logger
from mantid.simpleapi import CloneWorkspace, LoadHFIRSANS


# @pytest.mark.skip(reason="It doesn't pass on the build servers. "
#                          "XML lib incompatibility.")
def skip_test_momentum_tranfer_serial(biosans_f):

    ws = LoadHFIRSANS(
        # Filename='/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/'\
        #      'biosans/BioSANS_exp440_scan0022_0006.xml',
        Filename=biosans_f['anisotropic'],
        OutputWorkspace="ws")



    mt = MomentumTransfer(ws)
    assert mt.qx.shape == mt.qy.shape == mt.dqx.shape == mt.dqy.shape == \
        mt.i.shape == mt.i_sigma.shape == (256*192, )

    table_iq = mt.q2d()
    assert isinstance(table_iq, mantid.dataobjects.TableWorkspace)

    _, ws = mt.bin_into_q2d()
    assert ws.extractY().shape == (256, 192)
    assert ws.extractX().shape == (256, 193)

    _, ws = mt.bin_into_q1d()
    assert ws.extractY().shape == (1, 100)
    assert ws.extractX().shape == (1, 101)

    _, ws_iq_feature = mt.bin_wedge_into_q1d(phi_0=0, phi_aperture=30)
    ws_iq_feature_i = ws_iq_feature.extractY()
    assert ws_iq_feature_i.shape == (1, 100)
    assert ws_iq_feature.extractX().shape == (1, 101)
    ws_iq_feature = CloneWorkspace(ws_iq_feature)

    _, ws_iq_non_feature = mt.bin_wedge_into_q1d(phi_0=0 + 90, phi_aperture=30)
    ws_iq_non_feature_i = ws_iq_non_feature.extractY()
    assert ws_iq_non_feature_i.shape == (1, 100)
    assert ws_iq_non_feature.extractX().shape == (1, 101)
    ws_iq_non_feature = CloneWorkspace(ws_iq_non_feature)

    # Wedge with feature has more counts than that without it
    assert (ws_iq_feature_i.sum() - ws_iq_non_feature_i.sum()) > 200


def bin_in_parallel(params):
    """Runs mantid algorithms in parallel

    Parameters
    ----------
    params : list
        The parameters of this function have to be passed in a list.
        `map` does not support multiple algorithms

    Returns
    -------
    list
        Mantid loses `ws.name()`. We passing a list of ws_name amd WSs.
    """

    component_name, out_ws_prefix = params
    ws = mtd["ws_data_raw"]
    mt = MomentumTransfer(ws,
                          component_name=component_name,
                          out_ws_prefix=out_ws_prefix)

    table_iq = mt.q2d()
    ws_q2d_name, ws_q2d = mt.bin_into_q2d()
    ws_q1d_name, ws_q1d = mt.bin_into_q1d()

    return ((table_iq.name(), table_iq), (ws_q2d_name, ws_q2d), (ws_q1d_name,
                                                                 ws_q1d))


@pytest.mark.skip(reason="Only works as standalone test.")
def skip_test_momentum_tranfer_parallel(biosans_f):
    '''
    Note that we are using `pathos`. That's the only way to serialize
    Mantid WSs back to python.
    '''

    ws = LoadHFIRSANS(
        # Filename='/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/'
        # 'biosans/BioSANS_exp440_scan0022_0006.xml',
        Filename=biosans_f['anisotropic'],
        OutputWorkspace="ws_data_raw")

    parameters = [('detector1', 'main'), ('wing_detector', 'wing')]

    with Pool(processes=2) as pool:
        return_list = pool.map(bin_in_parallel, parameters)

    for ret in return_list:
        for ws_name, ws in ret:
            logger.information("Testing WS: {}".format(ws_name))
            if ws_name.startswith('wing'):
                y_shape = (256, 160)
                x_shape = (256, 161)
            else:
                y_shape = (256, 192)
                x_shape = (256, 193)

            if ws_name.endswith('qxqy'):
                assert ws.extractY().shape == y_shape
                assert ws.extractX().shape == x_shape

            if ws_name.endswith('iq'):
                assert ws.extractY().shape == (1, 100)
                assert ws.extractX().shape == (1, 101)
