from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
from mantid.simpleapi import Load
from ornl.sans.sns.eqsans import correct_frame as cf
from ornl.settings import amend_config


def test_transmitted_bands():
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        ws = Load(Filename='EQSANS_86217')
        bands = cf.transmitted_bands(ws)
        assert_almost_equal((bands.lead.min, bands.lead.max),
                            (2.48, 6.78), decimal=2)
        assert_almost_equal((bands.skip.min, bands.skip.max),
                            (10.90, 15.23), decimal=2)


if __name__ == '__main__':
    pytest.main()
