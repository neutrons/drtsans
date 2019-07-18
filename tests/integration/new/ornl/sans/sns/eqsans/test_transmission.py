from __future__ import (absolute_import, division, print_function)

import pytest
from numpy.testing import assert_almost_equal
from os.path import join as pjn

from mantid.simpleapi import LoadNexus, CompareWorkspaces, SaveNexus

from ornl.settings import namedtuplefy, unique_workspace_dundername as uwd
from ornl.sans.sns.eqsans.geometry import insert_aperture_logs
from ornl.sans.sns.eqsans import (calculate_transmission,
                                  apply_transmission_correction)


@pytest.fixture(scope='module')
@namedtuplefy
def trans_fix(refd):
    data_dir = pjn(refd.new.eqsans, 'test_transmission')
    cmp_dir = pjn(data_dir, 'compare')

    def quick_compare(tentative, asset):
        r"""asset: str, name of golden standard nexus file"""
        ws = LoadNexus(pjn(cmp_dir, asset), OutputWorkspace=uwd())
        return CompareWorkspaces(tentative, ws).Result

    a = LoadNexus(pjn(data_dir, 'sample.nxs'))
    insert_aperture_logs(a)  # source and sample aperture diameters
    b = LoadNexus(pjn(data_dir, 'direct_beam.nxs'))
    c = LoadNexus(pjn(data_dir, 'sample_skip.nxs'))
    d = LoadNexus(pjn(data_dir, 'direct_beam_skip.nxs'))
    return dict(data_dir=data_dir, sample=a, reference=b, sample_skip=c,
                reference_skip=d, compare=quick_compare)


def test_calculate_transmission(trans_fix):
    # Raw transmission values
    raw = calculate_transmission(trans_fix.sample_skip,
                                 trans_fix.reference_skip, radius=50,
                                 fit_func=None)
    SaveNexus(raw, '/tmp/raw_transmission_skip.nxs')
    raw = calculate_transmission(trans_fix.sample, trans_fix.reference,
                                 fit_func=None)
    assert trans_fix.compare(raw, 'raw_transmission.nxs')
    # Fitted transmission values
    fitted = calculate_transmission(trans_fix.sample, trans_fix.reference)
    trans_fix.compare(fitted, 'fitted_transmission.nxs')
    assert_almost_equal(fitted.lead_mfit.OutputChi2overDoF, 1.1, decimal=1)


def test_apply_transmission(trans_fix):
    corr = apply_transmission_correction(trans_fix.sample, trans_fix.reference,
                                         output_workspace=uwd())
    trans_fix.compare(corr, 'sample_corrected.nxs')


if __name__ == '__main__':
    pytest.main()
