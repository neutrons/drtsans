from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
import ornl.sans.sns.eqsans.dark_current_subtract as dcs


def test_compute_log_ratio(eqsans_w):
    r1 = eqsans_w['data'].run()  # run logs for data workspace
    r2 = eqsans_w['darkcurrent'].run()
    assert_almost_equal(dcs.compute_log_ratio(r1, r2, 'duration'), 0.2, 2)
    assert_almost_equal(dcs.compute_log_ratio(r1, r2, 'proton_charge'), 0.2, 2)
    with pytest.raises(RuntimeError):
        dcs.compute_log_ratio(r1, r2, 'non_existent_log')


def test_duration_ratio(eqsans_w):
    r1 = eqsans_w['data'].run()  # run logs for data workspace
    r2 = eqsans_w['darkcurrent'].run()
    assert_almost_equal(dcs.duration_ratio(r1, r2), 0.2, 2)
    assert_almost_equal(dcs.duration_ratio(r1, r2, 'proton_charge'), 0.2, 2)
    assert_almost_equal(dcs.duration_ratio(r1, r2, 'non_existent_log'), 1, 2)


def test_init():
    alg = dcs.EQSANSDarkCurrentSubtract()
    alg.initialize()


def test_PyExec(eqsans_w):
    alg = dcs.EQSANSDarkCurrentSubtract()
    alg.initialize()
    alg.setProperties(dict(Data=eqsans_w['data'],
                           DarkCurrent=eqsans_w['darkcurrent']))
    alg.PyExec()
    w = alg.getProperty('OutputWorkspace').value
    assert_almost_equal(w.readY(23426)[0], 32.44, 2)


if __name__ == '__main__':
    pytest.main()