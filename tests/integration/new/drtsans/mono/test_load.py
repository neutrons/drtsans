# Test load GPSANS and BIOSANS data
import pytest
import os


def test_load_gpsans():
    """Test load GPSANS data

    Returns
    -------

    """
    nexus_file_name = '/HFIR/CG2/IPTS-23801/nexus/CG2_7116.nxs.h5'
    if not os.path.exists(nexus_file_name):
        pytest.skip('Skip due to NeXus file {} is not accessible.'.format(nexus_file_name))

    assert True


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
    nexus_file_name = '/HFIR/CG3/IPTS-23782/nexus/CG3_4829.nxs.h5'
    if not os.path.exists(nexus_file_name):
        pytest.skip('Skip due to NeXus file {} is not accessible.'.format(nexus_file_name))

    assert True


def test_load_biosans_overwrite_meta():
    """Test load BIOSANS data with overwriting sample/detector position related meta data

    Returns
    -------

    """
    assert True


if __name__ == '__main__':
    pytest.main([__file__])

