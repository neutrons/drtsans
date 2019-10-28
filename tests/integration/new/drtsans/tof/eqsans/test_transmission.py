from math import nan
import numpy as np
import pytest
from os.path import join as pjn
r"""
Hyperlinks to Mantid algorithms
CompareWorkspaces <https://docs.mantidproject.org/nightly/algorithms/CompareWorkspaces-v1.html>
LoadNexus <https://docs.mantidproject.org/nightly/algorithms/LoadNexus-v1.html>
SumSpectra <https://docs.mantidproject.org/nightly/algorithms/SumSpectra-v1.html>
"""
from mantid.simpleapi import CompareWorkspaces, LoadNexus, SumSpectra
from mantid.api import mtd
r"""
Hyperlinks to drtsans functions
amend_config, namedtuplefy, unique_workspace_dundername available at:
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
insert_aperture_logs <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/geometry.py>
prepare_data <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/api.py>
calculate_transmission <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/transmission.py>
apply_transmission_correction <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/transmission.py>
find_beam_center, center_detector <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/beam_finder.py>
"""  # noqa: E501
from drtsans.settings import amend_config, namedtuplefy, unique_workspace_dundername
from drtsans.tof.eqsans.geometry import insert_aperture_logs
from drtsans.tof.eqsans.api import prepare_data
from drtsans.tof.eqsans import (calculate_transmission, apply_transmission_correction, find_beam_center,
                                center_detector)


@pytest.fixture(scope='module')
@namedtuplefy  # converts a dictionary into a namedtuple for more pythonic dereferencing of dictionary keys
def test_data_9a_part_1():
    r"""
    Data for the test that calculates the raw transmission at zero scattering angle starting with sample and
    empty-beam datasets, integrated on wavelength, over a predefined region of interest.

    Dropbox link to the test data:
        <https://www.dropbox.com/s/8ttwjfa1u0q4cyq/calculate_transmission_test.xlsx>

    Returns
    -------
    dict
    """
    return dict(wavelength_range=[2.0, 5.0],  # data integrated over this wavelength range
                emtpy_reference=[[1, 2, 3, 1, 1, 2, 2, 3, 0, 0],  # uncertainties by taking the square root
                                 [3, 2, 2, 4, 5, 4, 6, 3, 2, 1],
                                 [1, 4, 6, 9, 13, 15, 5, 8, 3, 0],
                                 [7, 3, 8, 19, 25, 18, 65, 12, 4, 1],
                                 [2, 5, 9, 28, 79, 201, 41, 16, 2, 5],
                                 [0, 7, 11, 23, 128, 97, 50, 17, 3, 2],
                                 [3, 3, 9, 20, 27, 23, 18, 7, 4, 3],
                                 [1, 2, 5, 9, 9, 15, 4, nan, nan, 1],
                                 [2, 4, 2, 4, 3, 4, 1, nan, nan, 1],
                                 [4, 0, 1, 3, 1, 2, 0, 0, 2, 0, ]],
                sample=[[4, 5, 4, 3, 1, 5, 5, 7, 3, 3],
                        [4, 4, 3, 5, 7, 5, 7, 3, 2, 4],
                        [4, 5, 5, 9, 12, 12, 7, 8, 5, 3],
                        [9, 3, 10, 15, 20, 18, 49, 13, 4, 2],
                        [3, 5, 9, 23, 60, 155, 35, 13, 3, 8],
                        [2, 9, 13, 20, 96, 76, 41, 15, 6, 3],
                        [3, 5, 10, 15, 22, 20, 17, 10, 5, 6],
                        [1, 6, 7, 10, 8, 13, 5, nan, nan, 4],
                        [2, 6, 3, 7, 6, 4, 3, nan, nan, 1],
                        [3, 4, 1, 4, 5, 6, 3, 3, 3, 1]],
                radius=2.0,
                radius_unit='mm',
                transmission=0.7888,
                transmission_uncertainty=0.04071,  # 0.0419,
                precision=1.e-04  # precision when comparing test data with drtsans calculations
                )


@pytest.mark.parametrize('workspace_with_instrument',
                         [dict(name='EQSANS', Nx=10, Ny=10, dx=1.e-3, dy=1.e-3, zc=1.0)], indirect=True)
def test_calculate_transmission_single_bin(test_data_9a_part_1, reference_dir, workspace_with_instrument):
    r"""
    This test calculates the raw transmission at zero scattering angle starting with sample and empty-beam datasets,
    integrated on wavelength, over a predefined region of interest.

    This test implements Issue #175, addressing section 7.3 of the master document.
    Dropbox links to the test:
        <https://www.dropbox.com/s/1lejukntcx3g8p1/calculate_transmission_test.pdf>
        <https://www.dropbox.com/s/8ttwjfa1u0q4cyq/calculate_transmission_test.xlsx>

    Test introduces a detector array with the following properties:
        - 10 tubes.
        - 10 pixels per tube.
        - square pixels of size 1 mili-meter.
        - distance from sample of 1 meter (this is irrelevant for the test, though).

    dev - Jose Borreguero <borreguerojm@ornl.gov>
    SME - Changwoo Do <doc1@ornl.gov>

    - List of Mantid algorithms employed:
        FindCenterOfMassPosition <https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html>
    """
    data = test_data_9a_part_1  # let's cut-down on the verbosity

    # Load the empty reference into a workspace
    reference_counts = np.array(data.emtpy_reference, dtype=float).reshape((10, 10, 1))
    reference_workspace = unique_workspace_dundername()  # some temporary random name for the workspace
    workspace_with_instrument(axis_values=data.wavelength_range,
                              intensities=reference_counts, uncertainties=np.sqrt(reference_counts),
                              view='array', output_workspace=reference_workspace)
    assert mtd[reference_workspace].readY(44)[0] == 128

    # Load the sample into a workspace
    sample_counts = np.array(data.sample, dtype=float).reshape((10, 10, 1))
    sample_workspace = unique_workspace_dundername()  # some temporary random name for the workspace
    workspace_with_instrument(axis_values=data.wavelength_range,
                              intensities=sample_counts,  uncertainties=np.sqrt(sample_counts),
                              view='array', output_workspace=sample_workspace)
    assert mtd[sample_workspace].readY(44)[0] == 96

    # Find the beam center using the empty reference, and then center both reference and sample runs
    # We pass centering options to the underlying Mantid algorithm finding the center of the beam
    # <FindCenterOfMassPosition://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v1.html>
    beam_center = find_beam_center(reference_workspace,
                                   centering_options={'BeamRadius': data.radius, 'Tolerance': 0.1 * data.radius})
    center_detector(reference_workspace, *beam_center)
    center_detector(sample_workspace, *beam_center)

    # Calculate raw (no fitting) transmission at zero angle using drtsans
    transmission = calculate_transmission(sample_workspace, reference_workspace,
                                          radius=data.radius, radius_unit=data.radius_unit, fit_function=None)

    # Verify transmission and associated uncertainty
    assert transmission.readY(0)[0] == pytest.approx(data.transmission, abs=data.precision)
    assert transmission.readE(0)[0] == pytest.approx(data.transmission_uncertainty, abs=data.precision)


@pytest.fixture(scope='module')
@namedtuplefy
def trans_fix(reference_dir):
    data_dir = pjn(reference_dir.new.eqsans, 'test_transmission')
    cmp_dir = pjn(data_dir, 'compare')

    def quick_compare(tentative, asset):
        r"""asset: str, name of golden standard nexus file"""
        ws = LoadNexus(pjn(cmp_dir, asset), OutputWorkspace=unique_workspace_dundername())
        return CompareWorkspaces(tentative, ws, Tolerance=1.E-6).Result

    a = LoadNexus(pjn(data_dir, 'sample.nxs'))
    insert_aperture_logs(a)  # source and sample aperture diameters
    b = LoadNexus(pjn(data_dir, 'direct_beam.nxs'))
    c = LoadNexus(pjn(data_dir, 'sample_skip.nxs'))
    d = LoadNexus(pjn(data_dir, 'direct_beam_skip.nxs'))
    return dict(data_dir=data_dir, sample=a, reference=b, sample_skip=c,
                reference_skip=d, compare=quick_compare)


def test_masked_beam_center(reference_dir, trans_fix):
    r"""
    (this test was written previously to the testset with the instrument team)
    Test for an exception raised when the beam centers are masked
    """
    mask = pjn(trans_fix.data_dir, 'beam_center_masked.xml')
    with amend_config(data_dir=reference_dir.new.eqsans):
        s = prepare_data("EQSANS_88975", mask=mask, output_workspace=unique_workspace_dundername())
        d = prepare_data("EQSANS_88973", mask=mask, output_workspace=unique_workspace_dundername())
    with pytest.raises(RuntimeError, match=r'More than half of the detectors'):
        calculate_transmission(s, d, output_workspace=())
    s.delete()
    d.delete()


def test_calculate_raw_transmission(trans_fix):
    r"""
    (this test was written previously to the testset with the instrument team)
    """
    raw = calculate_transmission(trans_fix.sample, trans_fix.reference,
                                 fit_function=None)
    assert trans_fix.compare(raw, 'raw_transmission.nxs')
    # big radius because detector is not centered
    raw = calculate_transmission(trans_fix.sample_skip,
                                 trans_fix.reference_skip, radius=50,
                                 fit_function=None)
    assert trans_fix.compare(raw, 'raw_transmission_skip.nxs')


def test_calculate_fitted_transmission(trans_fix):
    r"""
    (this test was written previously to the testset with the instrument team)
    """
    #175 fixture for data, instead of loading from a json file
    fitted = calculate_transmission(trans_fix.sample, trans_fix.reference)
    assert trans_fix.compare(fitted, 'fitted_transmission.nxs')
    # big radius because detector is not centered
    fitted = calculate_transmission(trans_fix.sample_skip, trans_fix.reference_skip, radius=50)
    assert trans_fix.compare(fitted, 'fitted_transmission_skip.nxs')


def test_apply_transmission(trans_fix):
    r"""
    (this test was written previously to the testset with the instrument team)
    """
    trans = calculate_transmission(trans_fix.sample, trans_fix.reference)
    corr = apply_transmission_correction(trans_fix.sample, trans, output_workspace=unique_workspace_dundername())
    corr = SumSpectra(corr, OutputWorkspace=corr.name())
    trans_fix.compare(corr, 'sample_corrected.nxs')
    # big radius because detector is not centered
    trans = calculate_transmission(trans_fix.sample_skip, trans_fix.reference_skip, radius=50)
    corr = apply_transmission_correction(trans_fix.sample_skip, trans, output_workspace=unique_workspace_dundername())
    corr = SumSpectra(corr, OutputWorkspace=corr.name())
    trans_fix.compare(corr, 'sample_corrected_skip.nxs')


if __name__ == '__main__':
    pytest.main([__file__])
