import pytest
import numpy as np

# CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
from mantid.simpleapi import CreateWorkspace

# unique_workspace_dundername within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py> # noqa: 501
# normalize_dark_current <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/dark_current.py>  # noqa: E501
# SampleLogs within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans.dark_current import normalize_dark_current


@pytest.fixture(scope='module')
def data_test_16a():
    r"""
    Input and expected output taken from the intro to issue #174
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/174#dark-current-normalized-at-eq-sans>
    """
    l_min, l_step, l_max = 2.5, 0.1, 6.0
    return dict(t_frame=1. / 60. * 1000000.,  # time duration of a frame, in microseconds
                t_low=500,  # sample measurement tof range, in microseconds
                t_high=2000,
                t_sam=100,  # sample measurement duration, in seconds
                t_dc=3600,  # dark current duration, in seconds
                l_min=l_min,  # minimum wavelength
                l_max=l_max,
                l_step=l_step,  # wavelength bin width
                n_pixels=25,  # 25 pixels in the detector
                I_dc=[[40., 45., 50., 55., 60.],
                      [65., 70., 75., 80., 85.],
                      [90., 95., 100., 105., 110.],
                      [115., 120., 125., 130., 135.],
                      [140., 145., 150., 155., 160.]],
                I_dc_norm=[[0.027, 0.0304,  0.0337, 0.0371, 0.0405],
                           [0.0438, 0.0472, 0.0506, 0.054,  0.0573],
                           [0.0607, 0.0641, 0.0675, 0.0708, 0.0742],
                           [0.0776, 0.081,  0.0843, 0.0877, 0.0911],
                           [0.0944, 0.0978, 0.1012, 0.1046, 0.1079]],
                precision=1.e-4  # precision to compare reduction framework to test results
                )


def test_normalize_dark_current(data_test_16a):
    """Test of dark current normalization
    For details see https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/156
    and also https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/174

    dev - Jose Borreguero <borreguerojm@ornl.gov>, Steven Hahn <hahnse@ornl.gov>, Jiao Lin <linjiao@ornl.gov>
    SME - Changwoo Do <doc1@ornl.gov>

    **Mantid algorithms used:**
    :ref:`CreateWorkspace <algm-CreateWorkspace-v1>`,
    <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>

    **drtsans functions used:**
    ~drtsans.settings.unique_workspace_dundername
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
    ~drtsans.samplelogs.SampleLogs
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
    ~drtsans.tof.eqsans.dark_current.normalize_dark_current
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/dark_current.py>
    """
    wavelength_bin_boundaries = np.arange(data_test_16a['l_min'],
                                          data_test_16a['l_max'] + data_test_16a['l_step'] / 2.,
                                          data_test_16a['l_step'])

    # The dark current normalization method `normalize_dark_current` requires a sample run from
    # which to retrieve the wavelength bins and other info, such as the clipping values for the TOF frame.
    # The intensities of this run are not relevant for the purposes of normalizing the dark current,
    # so we set them to zero
    data_workspace = CreateWorkspace(DataX=wavelength_bin_boundaries, UnitX='Wavelength',
                                     DataY=np.zeros((data_test_16a['n_pixels'], len(wavelength_bin_boundaries) - 1)),
                                     NSpec=data_test_16a['n_pixels'], OutputWorkspace=unique_workspace_dundername())

    # The intensity in a detector pixel is the same for all wavelength bins by construction in the test. Thus,
    # we repeat the one given value per detector pixel to be the same for all wavelength bins.
    intensities_list = np.array(data_test_16a['I_dc']).flatten()
    intensities_list = np.repeat(intensities_list[:, np.newaxis], len(wavelength_bin_boundaries) - 1, axis=1)

    # The test does not sum the dark current intensities over all wavelength channels, but the reduction framework
    # does as it follows the master document. In order to match the test results we have to first divide the
    # dark current intensities by the number of wavelength channels
    intensities_list /= len(wavelength_bin_boundaries) - 1

    # The dark current workspace now becomes:
    dark_workspace = CreateWorkspace(DataX=wavelength_bin_boundaries, UnitX='Wavelength', DataY=intensities_list,
                                     NSpec=data_test_16a['n_pixels'], OutputWorkspace=unique_workspace_dundername())

    # Initialize the sample logs. In the reduction framework this would have happened after loading the events file
    # and converting to wavelength
    data_sample_logs = SampleLogs(data_workspace)
    data_sample_logs.insert('tof_frame_width', data_test_16a['t_frame'])
    data_sample_logs.insert('tof_frame_width_clipped', data_test_16a['t_frame'] - data_test_16a['t_low'] - data_test_16a['t_high'])  # noqa: E501
    data_sample_logs.insert('wavelength_min', data_test_16a['l_min'], unit='Angstrom')
    data_sample_logs.insert('wavelength_max', data_test_16a['l_max'], unit='Angstrom')
    data_sample_logs.insert('wavelength_lead_min', data_test_16a['l_min'], unit='Angstrom')
    data_sample_logs.insert('wavelength_lead_max', data_test_16a['l_max'], unit='Angstrom')
    data_sample_logs.insert('is_frame_skipping', False)

    # Initialize the dark current logs. Only the duration of the run is necessary, which is recorded by the data
    # acquisition software.
    dark_sample_log = SampleLogs(dark_workspace)
    dark_sample_log.insert('duration', data_test_16a['t_dc'])

    # Normalization obtained by calling the method id drtsans. Know that this normalization does not take
    # into account the duration of the sample run. The reason is that we can reuse this normalized dark current
    # with sample runs of different duration. The duration of the sample run is incorporated just prior to dark
    # current subtraction.
    dark_normalized = normalize_dark_current(dark_workspace, data_workspace,
                                             output_workspace=unique_workspace_dundername())

    # Compare to test result. Notice that test output, data_test_16a['I_dc_norm'], shows the
    # intensities for one wavelength channel, for all pixels. By construction, the intensities were set the same
    # for all wavelength channels. Thus, we select the first wavelength channel in our computed intensities ([0])
    #
    # Notice also that at this moment we incorporate the duration of the sample run ('t_sam')
    computed_intensities = data_test_16a['t_sam'] * np.transpose(dark_normalized.extractY())[0]
    test_intensities = np.array(data_test_16a['I_dc_norm']).flatten()
    assert computed_intensities == pytest.approx(test_intensities, abs=data_test_16a['precision'])


if __name__ == '__main__':
    pytest.main()
