from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
from os.path import join as pjn
from mantid.simpleapi import Load, CompareWorkspaces
from ornl.settings import amend_config
from ornl.sans.transmission import zero_angle_transmission, apply_transmission
from ornl.sans.sns.eqsans.transmission import beam_radius, fit_raw
from ornl.sans.sns.eqsans.geometry import insert_aperture_logs


def test_transmission(refd):
    data_dir = pjn(refd.new.eqsans, 'test_transmission')
    cmp_dir = pjn(data_dir, 'compare')

    def quick_compare(tentative, asset):
        _ws = Load(pjn(cmp_dir, asset))
        return CompareWorkspaces(tentative, _ws)

    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        # Prepare data
        sample = Load(pjn(data_dir, 'sample.nxs'))
        insert_aperture_logs(sample)  # source and sample aperture diameters
        reference = Load(pjn(data_dir, 'direct_beam.nxs'))

        # Calculate raw transmission
        raw = zero_angle_transmission(sample, reference,
                                      beam_radius(sample),
                                      'raw_transm',
                                      delete_temp_wss=False)
        assert raw.detids == [22912, 22913, 22914, 22915, 22916, 23169, 23170,
                              23680, 23681, 23682, 23683, 23684, 23937, 23938]
        quick_compare(raw.transmission, 'raw_transmission.nxs')
        quick_compare(raw.reference, 'reference_grouped.nxs')
        quick_compare(raw.sample, 'sample_grouped.nxs')

        # Fit the raw transmission
        fitted = fit_raw(raw.transmission, 'fitted_transm')
        quick_compare(fitted.transmission, 'fitted_transmission.nxs')
        assert_almost_equal(fitted.lead_mfit.OutputChi2overDoF, 1.1, decimal=1)

        # Apply the fitted transmission
        corr = apply_transmission(sample, 'sample_corr',
                                  trans_ws=fitted.transmission,
                                  theta_dependent=True)
        quick_compare(corr.OutputWorkspace, 'sample_corrected.nxs')
        quick_compare(corr.TransmissionWorkspace, 'theta_transmission.nxs')


if __name__ == '__main__':
    pytest.main()
