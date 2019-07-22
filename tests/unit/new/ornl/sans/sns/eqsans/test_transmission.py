from __future__ import (absolute_import, division, print_function)

import pytest
from os.path import join as pjn
from numpy.testing import assert_almost_equal
from mantid.simpleapi import Load
from ornl.sans.sns.eqsans.correct_frame import transmitted_bands
from ornl.sans.sns.eqsans.transmission import fit_band, fit_raw, beam_radius
from ornl.sans.sns.eqsans.geometry import insert_aperture_logs


@pytest.mark.offline
def test_fit_band(refd):
    raw = Load(pjn(refd.new.eqsans, 'test_transmission',
                   'raw_transmission.nxs'))
    bands = transmitted_bands(raw)
    fb = fit_band(raw, bands.lead, 'name=UserFunction,Formula=a*x+b')
    assert_almost_equal(fb.mfit.OutputChi2overDoF, 1.1, decimal=1)


@pytest.mark.offline
def test_fit_raw(refd):
    raw = Load(pjn(refd.new.eqsans, 'test_transmission',
                   'raw_transmission.nxs'))
    fitted = fit_raw(raw, 'fitted_transmission')
    ws = fitted.transmission
    assert ws.name() == 'fitted_transmission'
    assert_almost_equal(fitted.lead_mfit.OutputChi2overDoF, 1.1, decimal=1)
    assert_almost_equal(fitted.skip_mfit.OutputChi2overDoF, 1.1, decimal=1)


@pytest.mark.offline
def test_beam_radius(refd):
    data_dir = pjn(refd.new.eqsans, 'test_transmission')
    sample = Load(pjn(data_dir, 'sample.nxs'))
    insert_aperture_logs(sample)  # source and sample aperture diameters
    assert_almost_equal(beam_radius(sample), 10.4, decimal=1)


if __name__ == '__main__':
    pytest.main()
