from __future__ import (absolute_import, division, print_function)

import pytest
from os.path import join as pjn
from numpy.testing import assert_almost_equal
from mantid.simpleapi import Load
from ornl.sans.sns.eqsans.correct_frame import transmitted_bands
from ornl.sans.sns.eqsans.transmission import fit_band, fit_raw


def test_fit_band(refd):
    raw = Load(pjn(refd.new.eqsans, 'test_transmission',
                   'raw_transmission.nxs'))
    bands = transmitted_bands(raw)
    fb = fit_band(raw, bands.lead, 'name=UserFunction,Formula=a*x+b')
    assert_almost_equal(fb.mfit.OutputChi2overDoF, 32, decimal=0)


def test_fit_raw(refd):
    raw = Load(pjn(refd.new.eqsans, 'test_transmission',
                   'raw_transmission.nxs'))
    fitted = fit_raw(raw, 'fitted_transmission')
    ws = fitted.fit
    assert ws.name() == 'fitted_transmission'
    assert_almost_equal(fitted.lead_mfit.OutputChi2overDoF, 32, decimal=0)
    assert_almost_equal(fitted.skip_mfit.OutputChi2overDoF, 16, decimal=0)


if __name__ == '__main__':
    pytest.main()
