# Test load GPSANS and BIOSANS data
import pytest
import os
from drtsans.mono.load import load_events
from drtsans.mono.meta_data import get_sample_detector_offset
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import sample_detector_distance
from drtsans.load import move_instrument


def test_load_gpsans():
    """Test load GPSANS data

    Returns
    -------

    """
    nexus_file_name = '/HFIR/CG2/IPTS-23801/nexus/CG2_7116.nxs.h5'
    if not os.path.exists(nexus_file_name):
        pytest.skip('Skip due to NeXus file {} is not accessible.'.format(nexus_file_name))

    # Load data
    ws = load_events(nexus_file_name, output_workspace='gptest01', overwrite_instrument=True,
                     detector_offset=0, sample_offset=0)

    sample_offset, detector_offset = get_sample_detector_offset(ws, 0.0, None, None)

    assert sample_offset == 0.
    assert detector_offset == 0.


def test_load_gpsans_overwrite_meta():
    """Test load GPSANS data with overwriting sample/detector position related meta data

    Returns
    -------

    """

    assert True


def test_load_biosans():
    """Test load BIOSANS data

    Returns
    -------

    """
    # Decide to skip data or not
    nexus_file_name = '/HFIR/CG3/IPTS-23782/nexus/CG3_4829.nxs.h5'
    if not os.path.exists(nexus_file_name):
        pytest.skip('Skip due to NeXus file {} is not accessible.'.format(nexus_file_name))

    # Load data
    ws = load_events(nexus_file_name, output_workspace='biotest01', overwrite_instrument=True,
                     detector_offset=0, sample_offset=0)

    # Check current instrument setup and meta data (sample logs)
    logs = SampleLogs(ws)
    print('[TEST INFO] SampleToSi = {} mm'.format(logs.find_log_with_units('CG3:CS:SampleToSi', unit='mm')))
    print('[TEST INFO] Sample to detector distance = {} /{} meter'
          ''.format(sample_detector_distance(ws),
                    sample_detector_distance(ws, log_key='sample_detector_distance', search_logs=True)))

    # Calculate offset without any overwriting
    sample_offset, detector_offset = get_sample_detector_offset(ws, 'CS3:CS:SampleToSi', 71. * 1E-3)
    print('[TEST INFO] Sample offset = {}, Detector offset = {}'
          ''.format(sample_offset, detector_offset))

    # Move sample and detector
    ws = move_instrument(ws, sample_offset, detector_offset)

    # Verify
    new_sample_det_distance = sample_detector_distance(ws)
    print('[TEST INFO] Sample detector distance after moving = {}'.format(new_sample_det_distance))

    assert 1 == 4


def test_load_biosans_overwrite_meta():
    """Test load BIOSANS data with overwriting sample/detector position related meta data

    Returns
    -------

    """
    assert True


if __name__ == '__main__':
    pytest.main([__file__])
