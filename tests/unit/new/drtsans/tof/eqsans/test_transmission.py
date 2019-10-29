import pytest
from os.path import join as pjn
from numpy.testing import assert_almost_equal
from mantid.simpleapi import LoadNexus
from drtsans.settings import (namedtuplefy, unique_workspace_dundername as uwd)
from drtsans.tof.eqsans.correct_frame import transmitted_bands
from drtsans.tof.eqsans.transmission import fit_band, fit_raw_transmission, beam_radius


@pytest.fixture(scope='module')
@namedtuplefy
def trasmission_data(reference_dir):
    data_dir = pjn(reference_dir.new.eqsans, 'test_transmission')
    a = LoadNexus(pjn(data_dir, 'raw_transmission.nxs'))
    b = LoadNexus(pjn(data_dir, 'sample.nxs'))
    c = LoadNexus(pjn(data_dir, 'raw_transmission_skip.nxs'))
    d = LoadNexus(pjn(data_dir, 'sample_skip.nxs'))
    return dict(data_dir=data_dir, raw=a, sample=b, raw_skip=c, sample_skip=d)


def test_beam_radius(trasmission_data):
    assert_almost_equal(beam_radius(trasmission_data.sample), 11.1, decimal=1)
    assert_almost_equal(beam_radius(trasmission_data.sample_skip), 13.5,
                        decimal=1)


def test_fit_band(trasmission_data):
    # Non-skip mode
    bands = transmitted_bands(trasmission_data.raw)
    fb = fit_band(trasmission_data.raw, bands.lead)
    assert_almost_equal(fb.mantid_fit_output.OutputChi2overDoF, 1.1, decimal=1)
    # Frame-skipping mode
    bands = transmitted_bands(trasmission_data.raw_skip)
    fb = fit_band(trasmission_data.raw_skip, bands.lead)
    assert_almost_equal(fb.mantid_fit_output.OutputChi2overDoF, 1.1, decimal=1)
    fb = fit_band(trasmission_data.raw_skip, bands.skip)
    assert_almost_equal(fb.mantid_fit_output.OutputChi2overDoF, 3.6, decimal=0)


def test_fit_raw(trasmission_data):
    fitted = fit_raw_transmission(trasmission_data.raw, output_workspace=uwd())
    assert_almost_equal(fitted.lead_mantid_fit.OutputChi2overDoF, 1.1, decimal=1)
    fitted = fit_raw_transmission(trasmission_data.raw_skip, output_workspace=uwd())
    assert_almost_equal(fitted.lead_mantid_fit.OutputChi2overDoF, 1.1, decimal=1)
    assert_almost_equal(fitted.skip_mantid_fit.OutputChi2overDoF, 3.6, decimal=1)


if __name__ == '__main__':
    pytest.main([__file__])
