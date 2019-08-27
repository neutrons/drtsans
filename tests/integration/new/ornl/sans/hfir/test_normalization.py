import pytest
import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.hfir.normalisation import time, monitor

x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
z = (x + 5 * y) * 10 + 100


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x.tolist(),
                          'intensities': z.tolist()}],
                         indirect=True)
def test_normalization_by_time(generic_workspace):
    ws = generic_workspace
    I_sam = z  # choose sample data from g or z
    t_sam = 5  # any value
    SampleLogs(ws).insert('timer', t_sam, 'Second')
    I_samnorm = I_sam / t_sam
    ws_samnorm = time(ws)
    assert np.allclose(ws_samnorm.extractY().ravel(), I_samnorm.ravel())


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x.tolist(),
                          'intensities': z.tolist()}],
                         indirect=True)
def test_normalization_by_monitor(generic_workspace):
    ws = generic_workspace
    I_sam = z  # choose sample data from g or z
    factor_is = 10**8
    flux_sam = 5 * 10**8
    SampleLogs(ws).insert('monitor', flux_sam)
    I_samnorm = factor_is / flux_sam * I_sam
    ws_samnorm = monitor(ws)
    assert np.allclose(ws_samnorm.extractY().ravel(), I_samnorm.ravel())


if __name__ == '__main__':
    pytest.main()