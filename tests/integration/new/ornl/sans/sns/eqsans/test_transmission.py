from __future__ import (absolute_import, division, print_function)

import pytest
from os.path import join as pjn
from mantid.simpleapi import LoadNexus, CompareWorkspaces, SumSpectra

from ornl.settings import (amend_config, namedtuplefy,
                           unique_workspace_dundername as uwd)
from ornl.sans.sns.eqsans.geometry import insert_aperture_logs
from ornl.sans.sns.eqsans.api import prepare_data
from ornl.sans.sns.eqsans import (calculate_transmission,
                                  apply_transmission_correction)


@pytest.fixture(scope='module')
@namedtuplefy
def trans_fix(reference_dir):
    data_dir = pjn(reference_dir.new.eqsans, 'test_transmission')
    cmp_dir = pjn(data_dir, 'compare')

    def quick_compare(tentative, asset):
        r"""asset: str, name of golden standard nexus file"""
        ws = LoadNexus(pjn(cmp_dir, asset), OutputWorkspace=uwd())
        return CompareWorkspaces(tentative, ws, Tolerance=1.E-6).Result

    a = LoadNexus(pjn(data_dir, 'sample.nxs'))
    insert_aperture_logs(a)  # source and sample aperture diameters
    b = LoadNexus(pjn(data_dir, 'direct_beam.nxs'))
    c = LoadNexus(pjn(data_dir, 'sample_skip.nxs'))
    d = LoadNexus(pjn(data_dir, 'direct_beam_skip.nxs'))
    return dict(data_dir=data_dir, sample=a, reference=b, sample_skip=c,
                reference_skip=d, compare=quick_compare)


def test_masked_beam_center(reference_dir, trans_fix):
    r"""Test for an exception raised when the beam centers are masked"""
    mask = pjn(trans_fix.data_dir, 'beam_center_masked.xml')
    with amend_config(data_dir=reference_dir.new.eqsans):
        s = prepare_data("EQSANS_88975", mask=mask, output_workspace=uwd())
        d = prepare_data("EQSANS_88973", mask=mask, output_workspace=uwd())
    with pytest.raises(RuntimeError, match=r'More than half of the detectors'):
        calculate_transmission(s, d, output_workspace=())
    s.delete()
    d.delete()


def test_calculate_raw_transmission(trans_fix):
    raw = calculate_transmission(trans_fix.sample, trans_fix.reference,
                                 fit_func=None)
    assert trans_fix.compare(raw, 'raw_transmission.nxs')
    # big radius because detector is not centered
    raw = calculate_transmission(trans_fix.sample_skip,
                                 trans_fix.reference_skip, radius=50,
                                 fit_func=None)
    assert trans_fix.compare(raw, 'raw_transmission_skip.nxs')


def test_calculate_fitted_transmission(trans_fix):
    fitted = calculate_transmission(trans_fix.sample, trans_fix.reference)
    assert trans_fix.compare(fitted, 'fitted_transmission.nxs')
    # big radius because detector is not centered
    fitted = calculate_transmission(trans_fix.sample_skip,
                                    trans_fix.reference_skip, radius=50)
    assert trans_fix.compare(fitted, 'fitted_transmission_skip.nxs')


def test_apply_transmission(trans_fix):
    trans = calculate_transmission(trans_fix.sample, trans_fix.reference)
    corr = apply_transmission_correction(trans_fix.sample, trans,
                                         output_workspace=uwd())
    corr = SumSpectra(corr, OutputWorkspace=corr.name())
    trans_fix.compare(corr, 'sample_corrected.nxs')
    # big radius because detector is not centered
    trans = calculate_transmission(trans_fix.sample_skip,
                                   trans_fix.reference_skip, radius=50)
    corr = apply_transmission_correction(trans_fix.sample_skip, trans,
                                         output_workspace=uwd())
    corr = SumSpectra(corr, OutputWorkspace=corr.name())
    trans_fix.compare(corr, 'sample_corrected_skip.nxs')


if __name__ == '__main__':
    pytest.main()
