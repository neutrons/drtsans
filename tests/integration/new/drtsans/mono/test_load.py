# Test load GPSANS and BIOSANS data
import pytest
import os
# import numpy as np
from drtsans.mono.load import load_events
from drtsans.mono.meta_data import get_sample_detector_offset
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import sample_detector_distance
from drtsans.load import move_instrument
from mantid.simpleapi import AddSampleLogMultiple


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

    # Check current instrument setup and meta data (sample logs)
    logs = SampleLogs(ws)
    print('[TEST INFO] SampleToSi = {} mm'.format(logs.find_log_with_units('CG2:CS:SampleToSi', unit='mm')))
    raw_sample_det_distance = sample_detector_distance(ws, unit='m', search_logs=False)
    print('[TEST INFO] Sample to detector distance = {} /{} meter'
          ''.format(raw_sample_det_distance,
                    sample_detector_distance(ws, unit='m', log_key='sample_detector_distance', search_logs=True)))

    # sample and detector offsets can only be retrieved from a loaded workspace
    # This is a technical debt
    sample_offset, detector_offset = get_sample_detector_offset(ws, 'CG2:CS:SampleToSi', 0.)

    assert sample_offset == -0.088
    assert detector_offset == -0.088

    # Move instrument
    # Move sample and detector
    ws = move_instrument(ws, sample_offset, detector_offset)

    # Verify
    new_sample_det_distance = sample_detector_distance(ws, unit='m', search_logs=False)
    print('[TEST INFO] Sample detector distance after moving = {} meter'.format(new_sample_det_distance))
    print('[TEST INFO] Sample position = {}'.format(ws.getInstrument().getSample().getPos()))

    assert new_sample_det_distance == raw_sample_det_distance


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
    raw_sample_det_distance = sample_detector_distance(ws)
    print('[TEST INFO] Sample to detector distance = {} /{} meter'
          ''.format(raw_sample_det_distance,
                    sample_detector_distance(ws, log_key='sample_detector_distance', search_logs=True)))

    # Calculate offset without any overwriting
    # sample and detector offsets can only be retrieved from a loaded workspace
    # This is a technical debt
    sample_offset, detector_offset = get_sample_detector_offset(ws, 'CG3:CS:SampleToSi', 71. * 1E-3)
    print('[TEST INFO] Sample offset = {}, Detector offset = {}'
          ''.format(sample_offset, detector_offset))

    assert sample_offset == pytest.approx(0., 1E-12)
    assert detector_offset == pytest.approx(0., 1E-12)

    # In this file, SampleToSi is 71 mm. Thus there won't be any offset for sample and detector
    # assert np.testing.assert_allclose(sample_pos)

    # Second test on SampleToSi distance other than 71.00 mm
    # reset sample log CG3:CS:SampleToSi
    test_sample_si_distance = 74.21
    AddSampleLogMultiple(Workspace=ws, LogNames='{}'.format('CG3:CS:SampleToSi'),
                         LogValues='{}'.format(test_sample_si_distance),
                         LogUnits='mm')

    # Calculate offsets
    # Calculate offset without any overwriting
    sample_offset, detector_offset = get_sample_detector_offset(ws, 'CG3:CS:SampleToSi', 71. * 1E-3)
    print('[TEST INFO 2] Sample offset = {}, Detector offset = {}'
          ''.format(sample_offset, detector_offset))

    # Move sample and detector
    ws = move_instrument(ws, sample_offset, detector_offset)

    # Verify
    # re-calculate the sample detector distance
    new_sample_det_distance = sample_detector_distance(ws, unit='m', search_logs=False)
    # check whether the sample detector distance from meta data in workspace is consistent
    meta_sample_det_distance = sample_detector_distance(ws, unit='mm', search_logs=True)
    print('[TEST INFO 2] Sample detector distance after moving = {} meter'.format(new_sample_det_distance))
    print('[TEST INFO 2] Sample position = {}'.format(ws.getInstrument().getSample().getPos()))

    assert 1 == 3

    assert new_sample_det_distance == raw_sample_det_distance * 1E-3
    assert new_sample_det_distance == meta_sample_det_distance * 1E-3


def test_load_biosans_overwrite_meta():
    """Test load BIOSANS data with overwriting sample/detector position related meta data

    Returns
    -------

    """
    # Decide to skip data or not
    nexus_file_name = '/HFIR/CG3/IPTS-23782/nexus/CG3_4829.nxs.h5'
    if not os.path.exists(nexus_file_name):
        pytest.skip('Skip due to NeXus file {} is not accessible.'.format(nexus_file_name))

    # Load data
    ws = load_events(nexus_file_name, output_workspace='biotest02', overwrite_instrument=True,
                     detector_offset=0, sample_offset=0)

    # Check current instrument setup and meta data (sample logs)
    logs = SampleLogs(ws)
    print('[TEST INFO] SampleToSi = {} mm'.format(logs.find_log_with_units('CG3:CS:SampleToSi', unit='mm')))
    raw_sample_det_distance = sample_detector_distance(ws)
    print('[TEST INFO] Sample to detector distance = {} /{} meter'
          ''.format(raw_sample_det_distance,
                    sample_detector_distance(ws, log_key='sample_detector_distance', search_logs=True)))

    # Calculate offset with overwriting to sample-detector-distance
    sample_offset, detector_offset = get_sample_detector_offset(ws, 'CG3:CS:SampleToSi', 71. * 1E-3,
                                                                overwrite_sample_detector_distance=7.1234)
    print('[TEST INFO] Sample offset = {}, Detector offset = {}'
          ''.format(sample_offset, detector_offset))

    # Move sample and detector
    ws = move_instrument(ws, sample_offset, detector_offset)

    # Verify
    new_sample_det_distance = sample_detector_distance(ws, unit='m', search_logs=False)
    print('[TEST INFO 2] Sample detector distance after moving = {} meter'.format(new_sample_det_distance))
    print('[TEST INFO 2] Sample position = {}'.format(ws.getInstrument().getSample().getPos()))

    assert new_sample_det_distance == pytest.approx(7.1234, 1E-7)


if __name__ == '__main__':
    pytest.main([__file__])
