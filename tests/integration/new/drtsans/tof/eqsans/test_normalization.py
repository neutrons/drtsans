import pytest
import numpy as np

# CreateWorkspace <https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html>
from mantid.simpleapi import CreateWorkspace

# unique_workspace_dundername within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py> # noqa: 501
# SampleLogs within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
# time, monitor within <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/mono/normalisation.py>
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import (load_events, transform_to_wavelength,
                                normalise_by_time, normalise_by_monitor)


def test_normalise_by_time(reference_dir):
    w = load_events('EQSANS_68168', data_dir=reference_dir.new.eqsans)
    d = SampleLogs(w).duration.value
    w = transform_to_wavelength(w)
    y, e = sum(w.readY(42)), sum(w.readE(42))
    w = normalise_by_time(w)
    assert (sum(w.readY(42)), sum(w.readE(42))) == pytest.approx((y / d, e / d))
    assert SampleLogs(w).normalizing_duration.value == 'duration'
    w.delete()


x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
sigma, mu = 0.3, 0  # gaussian image
d = np.sqrt(x * x + y * y)
g = np.exp(-((d - mu) ** 2 / (2 * sigma ** 2)))


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x,
                           'intensities': g}],
                         indirect=True)
def test_normalization_by_time(generic_workspace):
    ws = generic_workspace
    I_sam = g  # choose sample data from g or z

    t_sam = 5  # any value
    SampleLogs(ws).insert('timer', t_sam, 'Second')
    I_samnorm = I_sam / t_sam
    ws_samnorm = normalise_by_time(ws, log_key='timer')
    assert np.allclose(ws_samnorm.extractY().ravel(), I_samnorm.ravel())


@pytest.fixture(scope='module')
def data_test_16a_by_monitor():
    return dict(precision=1e-04,  # desired precision for comparisons,
                n_pixels=25,
                wavelength_bins=[2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0],
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
                fm=[5, 5, 4, 4, 3, 3, 3, 3, 2, 2],   # flux to monitor ratio
                phi=[20, 40, 30, 25, 20, 10, 5, 5, 5, 5],  # monitor spectrum
                I_samnorm=[[[0.4, 0.45, 0.5, 0.55, 0.6],
                            [0.65, 0.7, 0.75, 0.8, 0.85],
                            [0.9, 0.95, 1., 1.05, 1.1],
                            [1.15, 1.2, 1.25, 1.3, 1.35],
                            [1.4, 1.45, 1.5, 1.55, 1.6]],
                           [[0.2, 0.225, 0.25, 0.275, 0.3],
                            [0.325, 0.35, 0.375, 0.4, 0.425],
                            [0.45, 0.475, 0.5, 0.525, 0.55],
                            [0.575, 0.6, 0.625, 0.65, 0.675],
                            [0.7, 0.725, 0.75, 0.775, 0.8]],
                           [[0.3333, 0.375, 0.4167, 0.4583, 0.5],
                            [0.5417, 0.5833, 0.625, 0.6667, 0.7083],
                            [0.75, 0.7917, 0.8333, 0.875, 0.9167],
                            [0.9583, 1., 1.0417, 1.0833, 1.125],
                            [1.1667, 1.2083, 1.25, 1.2917, 1.3333]],
                           [[0.4, 0.45, 0.5, 0.55, 0.6],
                            [0.65, 0.7, 0.75, 0.8, 0.85],
                            [0.9, 0.95, 1.,  1.05, 1.1],
                            [1.15, 1.2, 1.25, 1.3, 1.35],
                            [1.4, 1.45, 1.5, 1.55, 1.6]],
                           [[0.6667, 0.75, 0.8333, 0.9167, 1.],
                            [1.0833, 1.1667, 1.25, 1.3333, 1.4167],
                            [1.5, 1.5833, 1.6667, 1.75, 1.8333],
                            [1.9167, 2., 2.0833, 2.1667, 2.25],
                            [2.3333, 2.4167, 2.5, 2.5833, 2.6667]],
                           [[1.3333, 1.5, 1.6667, 1.8333, 2.],
                            [2.1667, 2.3333, 2.5, 2.6667, 2.8333],
                            [3., 3.1667, 3.3333, 3.5,  3.6667],
                            [3.8333, 4., 4.1667, 4.3333, 4.5],
                            [4.6667, 4.8333, 5., 5.1667, 5.3333]],
                           [[2.6667, 3., 3.3333, 3.6667, 4.],
                            [4.3333, 4.6667, 5., 5.3333, 5.6667],
                            [6., 6.3333, 6.6667, 7., 7.3333],
                            [7.6667, 8., 8.3333, 8.6667, 9.],
                            [9.3333, 9.6667, 10., 10.3333, 10.6667]],
                           [[2.6667, 3., 3.3333, 3.6667, 4.],
                            [4.3333, 4.6667, 5., 5.3333, 5.6667],
                            [6., 6.3333, 6.6667, 7., 7.3333],
                            [7.6667, 8., 8.3333, 8.6667, 9.],
                            [9.3333, 9.6667, 10., 10.3333, 10.6667]],
                           [[4., 4.5, 5., 5.5, 6.],
                            [6.5, 7., 7.5, 8., 8.5],
                            [9., 9.5, 10., 10.5, 11.],
                            [11.5, 12., 12.5, 13., 13.5],
                            [14., 14.5, 15., 15.5, 16.]],
                           [[4., 4.5, 5., 5.5, 6.],
                            [6.5, 7., 7.5, 8., 8.5],
                            [9., 9.5, 10., 10.5, 11.],
                            [11.5, 12., 12.5, 13., 13.5],
                            [14., 14.5, 15., 15.5, 16.]]]
                )


def test_normalization_by_monitor_spectrum(data_test_16a_by_monitor):
    intensities_list = np.array(data_test_16a_by_monitor['I_sam']).flatten()
    errors_list = np.array(data_test_16a_by_monitor['I_sam_err']).flatten()
    # The intensity in a detector pixel is the same for all wavelength bins
    intensities_list = np.repeat(intensities_list[:, np.newaxis], len(data_test_16a_by_monitor['wavelength_bins']) - 1,
                                 axis=1)
    errors_list = np.repeat(errors_list[:, np.newaxis], len(data_test_16a_by_monitor['wavelength_bins']) - 1, axis=1)
    data_workspace = CreateWorkspace(DataX=data_test_16a_by_monitor['wavelength_bins'],
                                     DataY=intensities_list,
                                     DataE=errors_list,
                                     NSpec=data_test_16a_by_monitor['n_pixels'],
                                     OutputWorkspace=unique_workspace_dundername())
    SampleLogs(data_workspace).insert('is_frame_skipping', False)
    # In the reduction framework, counts at the monitor will be stored in a Mantid workspace
    monitor_workspace = CreateWorkspace(DataX=data_test_16a_by_monitor['wavelength_bins'],
                                        DataY=data_test_16a_by_monitor['fm'],
                                        DataE=np.zeros(len(data_test_16a_by_monitor['wavelength_bins']) - 1),
                                        NSpec=1,
                                        OutputWorkspace=unique_workspace_dundername())
    # In the reduction framework, the flux-to-monitor file will be loaded to a Mantid workspace
    flux_to_monitor_workspace = CreateWorkspace(DataX=data_test_16a_by_monitor['wavelength_bins'],
                                                DataY=data_test_16a_by_monitor['phi'],
                                                DataE=np.zeros(len(data_test_16a_by_monitor['wavelength_bins']) - 1),
                                                NSpec=1,
                                                OutputWorkspace=unique_workspace_dundername())
    # Carry out the normalization
    data_workspace = normalise_by_monitor(data_workspace, flux_to_monitor_workspace, monitor_workspace)
    # Compare to test data. Notice that data_test_16a_by_monitor['I_samnorm'] has shape (10, 5, 5) but
    # data_workspace.extractY() has shape (25, 10). A transpose operation is necessary
    test_intensities = np.transpose(np.array(data_test_16a_by_monitor['I_samnorm']).reshape((10, 25)))
    assert data_workspace.extractY() == pytest.approx(test_intensities, abs=data_test_16a_by_monitor['precision'])


if __name__ == '__main__':
    pytest.main([__file__])
