from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
from os.path import join as pjn
from mantid.simpleapi import Load, CompareWorkspaces
from ornl.settings import amend_config
from ornl.sans.transmission import (calculate_transmission,
                                    _apply_transmission_mantid)
from ornl.sans.sns.eqsans.transmission import beam_radius, fit_raw
from ornl.sans.sns.eqsans.geometry import insert_aperture_logs


def test_transmission(refd):
    data_dir = pjn(refd.new.eqsans, 'test_transmission')
    cmp_dir = pjn(data_dir, 'compare')

    def quick_compare(tentative, asset):
        _ws = Load(pjn(cmp_dir, asset))
        return CompareWorkspaces(tentative, _ws)

    with amend_config({'instrumentName': 'EQSANS'}):
        # Prepare data
        sample = Load(pjn(data_dir, 'sample.nxs'))
        insert_aperture_logs(sample)  # source and sample aperture diameters
        reference = Load(pjn(data_dir, 'direct_beam.nxs'))

        # Calculate raw transmission
        raw_transmission = calculate_transmission(sample, reference,
                                                  beam_radius(sample))
        quick_compare(raw_transmission, 'raw_transmission.nxs')

        # Fit the raw transmission
        fitted = fit_raw(raw_transmission, 'fitted_transm')
        quick_compare(fitted.transmission, 'fitted_transmission.nxs')
        assert_almost_equal(fitted.lead_mfit.OutputChi2overDoF, 1.1, decimal=1)

        # Apply the fitted transmission
        corr = _apply_transmission_mantid(sample, trans_ws=fitted.transmission,
                                          theta_dependent=True)
        quick_compare(corr, 'sample_corrected.nxs')


if __name__ == '__main__':
    pytest.main()
