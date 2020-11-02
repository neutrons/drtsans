import pytest
import os
import numpy as np
from mantid.simpleapi import LoadHFIRSANS, mtd


def test_benchmark_spice(reference_dir):
    """Test the benchmark SPICE file that is created to expose and verify the pixel mapping issue
    between old BIOSANS IDF and new BIOSANS IDF
    """
    # Access the test spice file
    spice_name = os.path.join(reference_dir.new.biosans, 'BioSANS_exp549_scan0020_0001_benchmark.xml')
    assert os.path.exists(spice_name)

    # Load data
    spice_ws = LoadHFIRSANS(Filename=spice_name, OutputWorkspace='CG3_5490020001_Benchmark')
    assert spice_ws

    # Test geometry
    shift_pixels = 127  # this is from the benchmark setup

    # main detector: 192 tubes
    main_det_tuple = list()
    for itube in range(192):
        det_id = 2 + itube * 256 + shift_pixels
        # get position X
        det_pos_x = spice_ws.getDetector(det_id).getPos().X()
        det_pos_y = spice_ws.getDetector(det_id).getPos().Y()
        # get count
        count = int(spice_ws.readY(det_id)[0])
        main_det_tuple.append((det_pos_x, det_pos_y, count))

    # Verify the positions
    # sort
    main_det_tuple.sort(reverse=True)
    # split
    pos_x_list, pos_y_list, count_list = zip(*main_det_tuple)

    # x position: shall be linear decreasing with constant step: from positive X to negative X
    pos_x_array = np.array(pos_x_list)
    pixel_distance_array = pos_x_array[1:] - pos_x_array[:-1]
    assert pixel_distance_array.mean() == pytest.approx(-0.0055)
    assert pixel_distance_array.std() < 5E-17

    # y position: shall be same
    pos_y_array = np.array(pos_y_list)
    assert pos_y_array.std() < 1E-17
    assert pos_y_array.mean() == pytest.approx(-0.00215)

    # counts
    assert np.allclose(np.array(count_list), np.arange(1, 192 + 1))

    # Wing detector: 160 tubes
    main_det_tuple = list()
    for itube in range(160):
        det_id = 2 + (192 + itube) * 256 + shift_pixels
        # get position X
        det_pos_x = spice_ws.getDetector(det_id).getPos().X()
        det_pos_y = spice_ws.getDetector(det_id).getPos().Y()
        # get count
        count = int(spice_ws.readY(det_id)[0])
        main_det_tuple.append((det_pos_x, det_pos_y, count))

    # Verify the positions
    # sort
    main_det_tuple.sort(reverse=True)
    # split
    pos_x_list, pos_y_list, count_list = zip(*main_det_tuple)

    # x position: shall be linear decreasing with constant step: from positive X to negative X
    pos_x_array = np.array(pos_x_list)
    pixel_distance_array = pos_x_array[1:] - pos_x_array[:-1]
    # always to the right but since it is curved not constant
    assert (len(pixel_distance_array[pixel_distance_array >= 0])) == 0

    # y position: shall be same
    pos_y_array = np.array(pos_y_list)
    assert pos_y_array.std() < 1E-17
    assert pos_y_array.mean() == pytest.approx(-0.00215)

    # counts
    assert np.allclose(np.array(count_list), np.arange(1, 160 + 1) * 1000)



if __name__ == '__main__':
    pytest.main([__file__])
