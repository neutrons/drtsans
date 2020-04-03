from pathlib import Path
import pytest
import numpy as np
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.samplelogs import SampleLogs
from drtsans.mono.gpsans.load import load_histogram
from drtsans.mono.load import load_and_split
from drtsans.geometry import sample_detector_distance


def test_load_histogram(reference_dir):
    filename = str(Path(reference_dir.new.gpsans) / 'CG2_exp206_scan0016_0001.xml')

    name = uwd()
    w = load_histogram(filename, wavelength=4.2, wavelength_spread=0.42, output_workspace=name)
    assert w.name() == name
    z = w.getInstrument().getComponentByName('detector1').getPos()[-1]
    assert z == 6.055
    max_id = np.argmax(w.extractY()[2:].flatten())
    tube_id, pixel_id = max_id // 256, max_id % 256
    assert (tube_id, pixel_id) == (90, 120)
    assert SampleLogs(w)['sample-detector-distance'].value == pytest.approx(6055.0, abs=0.1)

    # override the the sample-to-detector-distance. Should not affect the detector with maximum intensity
    w = load_histogram(filename, wavelength=4.2, wavelength_spread=0.42,
                       sample_to_detector_distance=1.042, unit='m',  output_workspace=name)
    z = w.getInstrument().getComponentByName('detector1').getPos()[-1]
    assert z == 1.042
    max_id = np.argmax(w.extractY()[2:].flatten())
    assert (max_id // 256, max_id % 256) == (tube_id, pixel_id)
    assert SampleLogs(w)['sample-detector-distance'].value == pytest.approx(1042.0, abs=0.1)


def test_load_and_split(reference_dir):
    # Check that is fails with missing required parameters
    filename = str(Path(reference_dir.new.gpsans) / 'CG2_9177.nxs.h5')
    print('Try to load {}'.format(filename))
    with pytest.raises(ValueError) as excinfo:
        load_and_split(filename, data_dir=reference_dir.new.gpsans,
                       sample_to_si_name='CG2:CS:SampleToSi',
                       si_nominal_distance=0.)
    assert "Must provide with time_interval or log_name and log_value_interval" == str(excinfo.value)

    filtered_ws = load_and_split(filename, data_dir=reference_dir.new.biosans, time_interval=50,
                                 sample_to_si_name='CG2:CS:SampleToSi',
                                 si_nominal_distance=0.)

    # suppose to get 3 output workspaces
    assert filtered_ws.size() == 3

    # Get workspaces
    sliced_ws_list = [filtered_ws.getItem(i) for i in range(3)]

    # Verify
    verify_geometry_meta(sliced_ws_list, 19., 83., -83. * 1E-3)


def test_load_and_split_overwrite_ssd(reference_dir):
    """Overwrite sample-silicon-window distance

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Check that is fails with missing required parameters
    filename = str(Path(reference_dir.new.gpsans) / 'CG2_9177.nxs.h5')
    print('Try to load {}'.format(filename))

    with pytest.raises(ValueError) as excinfo:
        load_and_split(filename, data_dir=reference_dir.new.gpsans,
                       sample_to_si_name='CG2:CS:SampleToSi',
                       si_nominal_distance=0.)
    assert "Must provide with time_interval or log_name and log_value_interval" == str(excinfo.value)

    filtered_ws = load_and_split(filename, data_dir=reference_dir.new.biosans, time_interval=50,
                                 sample_to_si_name='CG2:CS:SampleToSi',
                                 si_nominal_distance=0.,
                                 sample_to_si_value=103.)

    # suppose to get 3 output workspaces
    assert filtered_ws.size() == 3

    # Get workspaces
    sliced_ws_list = [filtered_ws.getItem(i) for i in range(3)]

    # Verify
    # SDD = 19. m
    # Sample to Silicon window is changed from 83 to 103 by 20 mm
    # SDD is changed to 19.020 meter
    verify_geometry_meta(sliced_ws_list,
                         expected_sample_detector_distance=19.020,
                         expected_sample_si_distance=103.,
                         expected_sample_position_z=-103 * 1E-3)


def verify_geometry_meta(workspace_list, expected_sample_detector_distance,
                         expected_sample_si_distance,
                         expected_sample_position_z):
    """Assuming there are 2 workspaces in the group

    Parameters
    ----------
    workspace_list: ~list
        list of workspaces
    expected_sample_detector_distance: float
        .. ...  unit = meter
    expected_sample_si_distance: float
        ... ... unit = millimeter
    expected_sample_position_z: float
        ... ...

    Returns
    -------

    """
    for index, workspace in enumerate(workspace_list):
        print('[TEST] workspace {}: {}'.format(index, str(workspace)))

        # check SDD: unit meter
        sdd = sample_detector_distance(workspace, unit='m', search_logs=False)
        assert sdd == pytest.approx(expected_sample_detector_distance, 1E-4)

        # check sample silicon window distance: unit millimeter
        swd = SampleLogs(workspace)['CG2:CS:SampleToSi'].value
        assert swd == pytest.approx(expected_sample_si_distance, 1E-4)

        # sample position: unit meter
        sample_pos_z = workspace.getInstrument().getSample().getPos()[2]
        assert sample_pos_z == pytest.approx(expected_sample_position_z, 1E-4)


if __name__ == '__main__':
    pytest.main([__file__])
