import pytest
from pytest import approx
import re
import tempfile
from mantid.simpleapi import LoadHFIRSANS
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.settings import unique_workspace_name
from ornl.sans.save_ascii import save_ascii_1D, save_ascii_2D, save_xml_1D


def numbers_in_line(line, numbers):
    xyz = [float(s) for s in re.findall(r'\d+\.\d+', line)]
    return all([x == approx(n, rel=1.E-03) or x < 0.02
                for x, n in zip(xyz, numbers)])


def test_save_ascii(biosans_sensitivity_dataset):

    ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood'],
        OutputWorkspace=unique_workspace_name())

    wss_name_ws = bin_into_q2d(ws)
    assert len(wss_name_ws) == 3

    ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
    assert ws_iqxqy.extractY().shape == (256, 192)
    assert ws_iqxqy.extractX().shape == (256, 193)

    _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy)
    assert ws_iq.extractY().shape == (1, 100)
    assert ws_iq.extractX().shape == (1, 101)

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_1D(ws_iq, 'Test BioSANS', tmp.name)
        output_lines = tmp.readlines()
        numbers = (0.112749, 107368.964140, 2815.765996, 0.021596)
        assert numbers_in_line(output_lines[101], numbers) is True

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_xml_1D(ws_iq, 'Test BioSANS', tmp.name)
        output_lines = tmp.readlines()
        numbers = (0.112749, 121.167, 3.17761, 0.0215961)
        assert numbers_in_line(output_lines[110], numbers) is True

    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_ascii_2D(ws_iqxqy, ws_dqx, ws_dqy, 'Test BioSANS 2D', tmp.name)
        output_lines = tmp.readlines()
        numbers = (0.077098, 0.081494, 73.000000, 8.544004, 0.015055, 0.001741)
        assert numbers_in_line(output_lines[48388], numbers) is True


if __name__ == '__main__':
    pytest.main()
