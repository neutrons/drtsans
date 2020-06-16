import pytest
import os
from drtsans.h5_buffer import HDFNode
import h5py


def test_hello_world():
    assert 'Hello World!'


def test_copy_nexus(reference_dir, cleanfile):
    """Test creating event NeXus file, loading it and compare to the original event NeXus.

    Test data: GPSANS run 9166

    Returns
    -------

    """
    # Get the source file
    source_nexus = os.path.join(reference_dir, '.nxs.h5')

    # Duplicate the source file to the temporary directory
    # TODO - this will be replaced by tempfile for future
    output_dir = '/tmp/'
    cleanfile(output_dir)
    target_nexus = os.path.join(output_dir, 'duplicated_cg2_nxs.h5')

    # Load the source
    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = HDFNode(nexus_h5, None)

    # Duplicate
    target_file = h5py.File(target_nexus, 'w')
    source_root.write(target_file)
    target_file.close()

    # Load source file to workspace

    # Load the duplicated

    # Compare counts on each pixel

    # Compare pixels' positions

    # Compare meta data


if __name__ == '__main__':
    pytest.main(__file__)
