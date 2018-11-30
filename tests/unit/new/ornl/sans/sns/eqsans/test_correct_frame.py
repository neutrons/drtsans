from __future__ import (absolute_import, division, print_function)

import pytest
from pytest import approx
from numpy.testing import assert_almost_equal
from mantid.simpleapi import Load
from ornl.sans.sns.eqsans import correct_frame as cf
from ornl.settings import amend_config
from ornl.sans.geometry import source_detector_distance


def test_transmitted_bands():
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        ws = Load(Filename='EQSANS_86217')
        bands = cf.transmitted_bands(ws)
        assert_almost_equal((bands.lead.min, bands.lead.max),
                            (2.48, 6.78), decimal=2)
        assert_almost_equal((bands.skip.min, bands.skip.max),
                            (10.90, 15.23), decimal=2)


def test_transmitted_bands_clipped():
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        ws = Load(Filename='EQSANS_86217')
        sdd = source_detector_distance(ws, units='m')
        bands_0 = cf.transmitted_bands_clipped(ws, sdd, 0.0, 0.0)
        lwc, hwc = (0.139, 0.560)  # expected clippings
        # With no interior clipping
        bands = cf.transmitted_bands_clipped(ws, sdd, 500, 2000,
                                             interior_clip=False)
        # Check clippings for the lead pulse
        b1_0, b2_0 = bands_0.lead.min, bands_0.lead.max
        b1, b2 = bands.lead.min, bands.lead.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0), 0.01)
        # Check clippings for the skip pulse
        b1_0, b2_0 = bands_0.skip.min, bands_0.skip.max
        b1, b2 = bands.skip.min, bands.skip.max
        assert (b1, b2) == approx((b1_0, b2_0 - hwc), 0.01)
        # With interior clipping
        bands = cf.transmitted_bands_clipped(ws, sdd, 500, 2000,
                                             interior_clip=True)
        b1_0, b2_0 = bands_0.lead.min, bands_0.lead.max
        b1, b2 = bands.lead.min, bands.lead.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0 - hwc), 0.01)
        b1_0, b2_0 = bands_0.skip.min, bands_0.skip.max
        b1, b2 = bands.skip.min, bands.skip.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0 - hwc), 0.01)


if __name__ == '__main__':
    pytest.main()
