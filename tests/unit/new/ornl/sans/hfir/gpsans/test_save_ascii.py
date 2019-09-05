import re
import tempfile

import pytest
from pytest import approx

from mantid.simpleapi import LoadHFIRSANS
from ornl.sans.hfir.momentum_transfer import MomentumTransfer
from ornl.sans.save_ascii import save_ascii_1D, save_xml_1D
from ornl.settings import unique_workspace_name


def numbers_in_line(line, numbers):
    xyz = [float(s) for s in re.findall(r'\d+\.\d+', line)]
    return all(
        [x == approx(n, rel=1.E-03) or x < 0.02 for x, n in zip(xyz, numbers)])


@pytest.mark.skip(reason='skip test until issue #140 resolved')
def test_save_ascii(gpsans_f):

    ws = LoadHFIRSANS(Filename=gpsans_f['sample_transmission'],
                      OutputWorkspace=unique_workspace_name())

    mt = MomentumTransfer(ws)
    _, ws_iqxqy = mt.bin_into_q2d()
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = mt.bin_into_q1d()
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_1D(ws_iq, 'Test GPSANS', tmp.name)
        output_lines = tmp.readlines()
        numbers = (0.158771, 0.000000, 140.542940, 0.025318)
        assert numbers_in_line(output_lines[101], numbers) is True

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_xml_1D(ws_iq, 'Test GPSANS', tmp.name)
        output_lines = tmp.readlines()
        numbers = (0.158594, 0.235702, 0.0253123)

    # with tempfile.NamedTemporaryFile('r+') as tmp:
    #     save_ascii_2D(ws_iqxqy, ws_dqx, ws_dqy, 'Test 2D GPSANS', tmp.name)
    #     output_lines = tmp.readlines()
    #     numbers = (0.137103, 0.081288, 0.0, 1.0, 0.020047, 0.002688)
    #     assert numbers_in_line(output_lines[48900], numbers) is True


if __name__ == '__main__':
    pytest.main()
