import pytest
import itertools

# CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
from mantid.simpleapi import CreateWorkspace

# unique_workspace_dundername within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py> # noqa: 501
# SampleLogs within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
# time, monitor within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/normalisation.py>
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans.mono.normalisation import time, monitor


@pytest.fixture(scope='module')
def data_test_16a():
    return dict(t_sam=5,  # duration of the sample run
                wavelength_bin=[2.0, 3.0],  # some arbitrary wavelength bin to aid in the creation of a workspace
                precision=1e-04,  # desired precision for comparisons,
                n_pixels=25,
                I_sam=[[40.,  45.,  50.,  55.,  60.],
                       [65.,  70.,  75.,  80.,  85.],
                       [90.,  95., 100., 105., 110.],
                       [115., 120., 125., 130., 135.],
                       [140., 145., 150., 155., 160.]],
                I_sam_err=[[6.3246, 6.7082, 7.0711, 7.4162, 7.746],
                           [8.0623, 8.3666, 8.6603, 8.9443, 9.2195],
                           [9.4868, 9.7468, 10., 10.247, 10.4881],
                           [10.7238, 10.9545, 11.1803, 11.4018, 11.619],
                           [11.8322, 12.0416, 12.2474, 12.4499, 12.6491]],
                I_samnorm=[[8.,  9., 10., 11., 12.],
                           [13., 14., 15., 16., 17.],
                           [18., 19., 20., 21., 22.],
                           [23., 24., 25., 26., 27.],
                           [28., 29., 30., 31., 32.]],
                I_samnorm_err=[[1.26492, 1.34164, 1.41422, 1.48324, 1.5492],
                               [1.61246, 1.67332, 1.73206, 1.78886, 1.8439],
                               [1.89736, 1.94936, 2., 2.0494, 2.09762],
                               [2.14476, 2.1909, 2.23606, 2.28036, 2.3238],
                               [2.36644, 2.40832, 2.44948, 2.48998, 2.52982]],
                flux_sam=5e8,  # flux at the monitor
                )


def test_normalization_by_time(data_test_16a):
    r"""
    Normalize sample intensities by the duration of the run.

    devs - Steven Hahn <hahnse@ornl.gov>,
           Jiao Lin <linjiao@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>
    SME  - Changwoo Do <doc1@ornl.gov>

    **Mantid algorithms used:**
    :ref:`CreateWorkspace <algm-CreateWorkspace-v1>`,
    <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>

    **drtsans functions used:**
    ~drtsans.settings.unique_workspace_dundername
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
    ~drtsans.samplelogs.SampleLogs
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
    ~drtsans.mono.normalisation.time
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/normalisation.py>
    """
    # Create a Mantid workspace with the sample intensities
    intensities_list = list(itertools.chain(*data_test_16a['I_sam']))
    errors_list = list(itertools.chain(*data_test_16a['I_sam_err']))
    ws = CreateWorkspace(DataX=data_test_16a['wavelength_bin'],
                         DataY=intensities_list, DataE=errors_list,
                         NSpec=data_test_16a['n_pixels'],
                         OutputWorkspace=unique_workspace_dundername())
    # Insert the duration of the run as a metadata item
    SampleLogs(ws).insert('timer', data_test_16a['t_sam'], 'Second')
    # Carry out the normalization
    ws_samnorm = time(ws)
    # Compare normalized intensities to those of the test
    intensities_list = list(itertools.chain(*data_test_16a['I_samnorm']))
    assert ws_samnorm.extractY() == pytest.approx(intensities_list, abs=data_test_16a['precision'])
    # Compare normalized errors to those of the test
    errors_list = list(itertools.chain(*data_test_16a['I_samnorm_err']))
    assert ws_samnorm.extractE() == pytest.approx(errors_list, abs=data_test_16a['precision'])
    ws_samnorm.delete()  # some clean up


def test_normalization_by_monitor(data_test_16a):
    r"""
    Normalize sample intensities by flux at monitor

    devs - Steven Hahn <hahnse@ornl.gov>,
           Jiao Lin <linjiao@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>
    SME  - Changwoo Do <doc1@ornl.gov>

    **Mantid algorithms used:**
    :ref:`CreateWorkspace <algm-CreateWorkspace-v1>`,
    <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>

    **drtsans functions used:**
    ~drtsans.settings.unique_workspace_dundername
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
    ~drtsans.samplelogs.SampleLogs
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
    ~drtsans.mono.normalisation.monitor
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/normalisation.py>
    """
    # Create a Mantid workspace with the sample intensities
    intensities_list = list(itertools.chain(*data_test_16a['I_sam']))
    ws = CreateWorkspace(DataX=data_test_16a['wavelength_bin'],
                         DataY=intensities_list,
                         NSpec=data_test_16a['n_pixels'],
                         OutputWorkspace=unique_workspace_dundername())
    # Insert the flux at the monitor as a metadata item
    SampleLogs(ws).insert('monitor', data_test_16a['flux_sam'])
    ws_samnorm = monitor(ws)
    # Compare normalized intensities to those of the test
    intensities_list = list(itertools.chain(*data_test_16a['I_samnorm']))
    assert ws_samnorm.extractY() == pytest.approx(intensities_list, abs=data_test_16a['precision'])
    ws_samnorm.delete()  # some clean up


if __name__ == '__main__':
    pytest.main()
