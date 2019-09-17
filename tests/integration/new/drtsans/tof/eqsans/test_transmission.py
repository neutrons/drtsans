import pytest
import json
import numpy as np
from pathlib import Path
from os.path import join as pjn
from mantid.simpleapi import LoadNexus, CompareWorkspaces, SumSpectra

from drtsans.settings import (amend_config, namedtuplefy,
                              unique_workspace_dundername as uwd)
from drtsans.tof.eqsans.geometry import insert_aperture_logs
from drtsans.tof.eqsans.api import prepare_data
from drtsans.tof.eqsans import (calculate_transmission, apply_transmission_correction, find_beam_center,
                                center_detector)


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


@pytest.mark.parametrize('workspace_with_instrument',
                         [dict(name='EQSANS', Nx=10, Ny=10, dx=1e-3, dy=1e-3, zc=0.01)], indirect=True)
def test_calculate_transmission_single_bin(reference_dir, workspace_with_instrument):
    r"""
    This implements Issue #175, addressing master document, section 7.3
    dev - Jose Borreguero <borreguerojm@ornl.gov>
    SME - Changwoo Do <doc1@ornl.gov>
    - List of Mantid algorithms employed:
        FindCenterOfMassPosition <https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html>
    """
    # Gather reference and sample data, then check
    wavelength_bin = [2.0, 5.0]  # some arbitrary wavelength bin
    with open(Path(reference_dir.new.eqsans) / 'test_transmission' / 'transmission_single_bin.json') as json_file:
        data = json.load(json_file)
    reference_counts = np.array(data['reference_data_block'], dtype=float).reshape((10, 10, 1))
    reference_workspace = workspace_with_instrument(axis_values=wavelength_bin, intensities=reference_counts,
                                                    view='array')
    assert reference_workspace.readY(42)[0] == 9
    sample_counts = np.array(data['sample_data_block'], dtype=float).reshape((10, 10, 1))
    sample_workspace = workspace_with_instrument(axis_values=wavelength_bin, intensities=sample_counts, view='array')
    assert sample_workspace.readY(42)[0] == 8

    # Find the beam center and then center reference and sample
    beam_center = find_beam_center(reference_workspace)
    assert beam_center == pytest.approx([-0.00007, 0.00015], abs=1e-5)
    center_detector(reference_workspace, x=-beam_center[0], y=-beam_center[1], unit='m')
    center_detector(sample_workspace, x=-beam_center[0], y=-beam_center[1], unit='m')

    # Calculate transmission, then check
    transmission = calculate_transmission(sample_workspace, reference_workspace,
                                          radius=2.0, radius_unit='mm', fit_func=None)
    assert transmission.readY(0)[0] == pytest.approx(0.7888, abs=1e-04)


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
