import pytest
import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.hfir.normalisation import time, monitor

# sample data for integration tests
x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
I_sam = (x + 5 * y) * 10 + 100


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x,
                          'intensities': I_sam}],
                         indirect=True)
def test_normalization_by_time(generic_workspace):
    ws = generic_workspace
    t_sam = 5  # value from test
    SampleLogs(ws).insert('timer', t_sam, 'Second')
    # Calculate following equation 6.1 in master document
    I_samnorm = I_sam / t_sam
    # Calculate using API function
    ws_samnorm = time(ws)
    # Check results from both caluclations  match
    assert np.allclose(ws_samnorm.extractY().ravel(), I_samnorm.ravel())


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x,
                          'intensities': I_sam}],
                         indirect=True)
def test_normalization_by_monitor(generic_workspace):
    ws = generic_workspace
    factor_is = 10**8  # value from test
    flux_sam = 5 * 10**8  # value from test
    SampleLogs(ws).insert('monitor', flux_sam)
    # Calculate following equation 6.4 in master document
    I_samnorm = factor_is / flux_sam * I_sam
    # Calcuate using API function
    ws_samnorm = monitor(ws)
    # Check results match
    assert np.allclose(ws_samnorm.extractY().ravel(), I_samnorm.ravel())


if __name__ == '__main__':
    pytest.main()
