from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
import ornl.sans.sns.eqsans.dark_current_subtract as dcs
import mantid.simpleapi as mtds


def test_compute_log_ratio(porasil_slice1m):
    s = porasil_slice1m.w['s'].run()
    dc = porasil_slice1m.w['dc'].run()
    assert_almost_equal(dcs.compute_log_ratio(s, dc, 'duration'), 0.08, 2)
    assert_almost_equal(dcs.compute_log_ratio(s, dc, 'proton_charge'), 0.08, 2)
    with pytest.raises(RuntimeError):
        dcs.compute_log_ratio(s, dc, 'non_existent_log')


def test_duration_ratio(porasil_slice1m):
    s = porasil_slice1m.w['s'].run()
    dc = porasil_slice1m.w['dc'].run()
    assert_almost_equal(dcs.duration_ratio(s, dc), 0.08, 2)
    assert_almost_equal(dcs.duration_ratio(s, dc, 'proton_charge'), 0.08, 2)
    assert_almost_equal(dcs.duration_ratio(s, dc, 'non_existent_log'), 1, 2)


def test_subtract_pixelcount_dark(porasil_slice1m):
    s = porasil_slice1m.w['s']
    dc = porasil_slice1m.w['dc']
    w = dcs.subtract_pixelcount_dark(s, dc)
    id = int(s.extractY().argmax())  # detector with max number of counts
    assert w.getSpectrum(id).getNumberEvents() == 124578 + 67
    assert_almost_equal(w.readY(id)[0], 124578.0 - 0.0859 * 67, 0)


def test_subtract_isotropic_dark(porasil_slice1m):
    s = porasil_slice1m.w['s']
    dc = porasil_slice1m.w['dc']
    # Assume we are working in Wavelength units
    s = mtds.ConvertUnits(s, Target='Wavelength', Emode='Elastic')
    s = mtds.Rebin(s, Params=[0.1, 0.01, 5.0], PreserveEvents=True)
    dc = mtds.ConvertUnits(dc, Target='Wavelength', Emode='Elastic')
    dc = mtds.Rebin(dc, Params=[0.1, 0.01, 5.0], PreserveEvents=True)
    w = dcs.subtract_isotropic_dark(s, dc)
    assert w.getNumberEvents() == 77024214
    assert_almost_equal(w.getSpectrum(0).getWeights().sum(), 7.97455, 5)


def test_init():
    alg = dcs.EQSANSDarkCurrentSubtract()
    alg.initialize()


def test_PyExec(porasil_slice1m):
    alg = dcs.EQSANSDarkCurrentSubtract()
    alg.initialize()
    s = porasil_slice1m.w['s']

    # Pixel Count
    dc = porasil_slice1m.w['dc']
    alg.setProperties(dict(Data=s, DarkCurrent=dc, Method='PixelCount'))
    alg.PyExec()
    w = alg.getProperty('OutputWorkspace').value
    id = int(s.extractY().argmax())  # detector with max number of counts
    assert w.getSpectrum(id).getNumberEvents() == 124578 + 67
    assert_almost_equal(w.readY(id)[0], 124578.0 - 0.0859 * 67, 0)

    # Isotropic
    s = mtds.ConvertUnits(s, Target='Wavelength', Emode='Elastic')
    s = mtds.Rebin(s, Params=[0.1, 0.01, 5.0], PreserveEvents=True)
    dc = porasil_slice1m.w['dc']
    dc = mtds.ConvertUnits(dc, Target='Wavelength', Emode='Elastic')
    dc = mtds.Rebin(dc, Params=[0.1, 0.01, 5.0], PreserveEvents=True)
    alg.setProperties(dict(Data=s, DarkCurrent=dc, Method='Isotropic'))
    alg.PyExec()
    w = alg.getProperty('OutputWorkspace').value
    assert w.getNumberEvents() == 77024214
    assert_almost_equal(w.getSpectrum(0).getWeights().sum(), 7.97455, 5)


if __name__ == '__main__':
    pytest.main()