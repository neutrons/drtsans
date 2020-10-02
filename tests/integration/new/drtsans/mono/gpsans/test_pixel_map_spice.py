import pytest
import os
import numpy as np
from tempfile import mkdtemp
from drtsans.mono.bar_scan_pixel_calibration import generate_pixel_map_legacy
from mantid.simpleapi import LoadNexusProcessed


def test_pixel_calibration(reference_dir, cleanfile):
    """

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Set and clean output
    test_output_dir = '/tmp/test_barscan/'
    # TODO cleanfile(test_output_dir)

    expected_calib_nexus = os.path.join(reference_dir.new.gpsans,
                                        'calibrations/CG2_Pixel_Calibration_Expected_111.nxs')
    assert os.path.exists(expected_calib_nexus)
    # test_calib_nexus = os.path.join(test_output_dir, 'runs_1_111/tables/Test_CG2_Pixel_Calibration.nxs')

    test_calib_nexus = generate_pixel_map_legacy()
    print(f'Calibraton file {test_calib_nexus} of type {type(test_calib_nexus)}')
    assert os.path.exists(test_calib_nexus)

    # Compare 2 NeXus file
    compare_pixel_calibration_files(test_calib_nexus, expected_calib_nexus)


def compare_pixel_calibration_files(test_file_name, gold_file_name):
    """Compare 2 calibration file by the pixel calibration value

    Algorithm:
    1. Load both processed NeXus files
    2. Compare workspaces' Y value

    Parameters
    ----------
    test_file_name
    gold_file_name

    Returns
    -------

    """
    # Load calibration file to Mantid workspace
    test_calib_ws = LoadNexusProcessed(Filename=test_file_name)
    gold_calib_ws = LoadNexusProcessed(Filename=gold_file_name)

    # Get calibration values
    test_cal_values = np.array(test_calib_ws.column(1)).flatten()
    gold_cal_values = np.array(gold_calib_ws.column(1)).flatten()

    np.testing.assert_allclose(test_cal_values, gold_cal_values)


if __name__ == '__main__':
    pytest.main(__file__)
