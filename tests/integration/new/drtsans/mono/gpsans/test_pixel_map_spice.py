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

    # template CG2 nexus file for IDF
    TEMPLATE_EVENT_NEXUS = '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/hfir/gpsans/CG2_9177.nxs.h5'

    # Set template event nexus
    template_event_nexus = TEMPLATE_EVENT_NEXUS
    assert os.path.exists(template_event_nexus), f'Template nexus {template_event_nexus} cannot be found.'

    # First and last pt for the barscan: Set by user
    # -------------------------------------------------------------------------------------------------------
    # IPTS 828 Exp 280.  (/HFIR/CG2/IPTS-828/exp280/Datafiles)
    root_dir = '/HFIR/CG2/'
    ipts = 828
    exp_number = 280
    scan_number = 5
    first_pt = 1
    last_pt = 111

    flood_ipts = 828
    flood_exp = 280
    flood_scan = 4
    flood_pt = 1

    # Calculate pixel calibration file
    test_calib_nexus = generate_pixel_map_legacy(ipts, exp_number, scan_number, range(first_pt, last_pt + 1),
                                                 flood_ipts, flood_exp, flood_scan, flood_pt,
                                                 root_dir, template_event_nexus,
                                                 test_output_dir)
    print(f'Calibraton file {test_calib_nexus} of type {type(test_calib_nexus)}')
    assert os.path.exists(test_calib_nexus)

    # Get expected data file
    expected_calib_nexus = os.path.join(reference_dir.new.gpsans,
                                        'calibrations/CG2_Pixel_Calibration_Expected_111.nxs')
    assert os.path.exists(expected_calib_nexus), f'Gold result (file) {expected_calib_nexus} cannot be found.'

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
