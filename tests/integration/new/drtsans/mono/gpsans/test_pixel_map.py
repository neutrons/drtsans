import pytest
import os
from tempfile import mkdtemp
from drtsans.mono.bar_scan_pixel_calibration import generate_pixel_map_legacy


def test_pixel_calibration(reference_dir, cleanfile):
    """

    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    expected_calib_nexus = os.path.join(reference_dir.new.gpsans, 'calibrations/CG2_Pixel_Calibration_Expected_111.nxs')
    assert os.path.exists(expected_calib_nexus)
    test_calib_nexus = '/tmp/test_pixel_calibration/runs_1_111/tables/Test_CG2_Pixel_Calibration.nxs'



if __name__ == '__main__':
    pytest.main(__file__)
