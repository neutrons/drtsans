import pytest
import numpy as np
import os
from drtsans.h5_buffer import HDFNode
from drtsans.load import load_events
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
    test_nexus_name = 'CG2_9177.nxs.h5'
    source_nexus = os.path.join(reference_dir.new.gpsans, test_nexus_name)
    assert os.path.exists(source_nexus), f'Test data {source_nexus} does not exist'

    # Duplicate the source file to the temporary directory
    # TODO - this will be replaced by tempfile for future
    output_dir = '/tmp/nexus'
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir('/tmp/nexus')
    target_nexus = os.path.join(output_dir, 'CG2_9177.nxs.h5')

    # Load the source
    nexus_h5 = h5py.File(source_nexus, 'r')
    source_root = HDFNode(nexus_h5, None)

    # Duplicate
    target_file = h5py.File(target_nexus, 'w')
    source_root.write(target_file)
    target_file.close()

    # Load source file to workspace
    source_ws = load_events(test_nexus_name, output_workspace='cg2_source')

    # Load the duplicated
    target_ws = load_events(target_nexus, output_workspace='cg2_duplicated')

    # Compare counts on each pixel
    source_y = source_ws.extractY()
    target_y = target_ws.extractY()
    np.testing.assert_allclose(source_y, target_y)

    # Compare pixels' positions
    num_hist = source_ws.getNumberHistograms()
    for iws in range(0, num_hist, 100):
        source_det_i_pos = source_ws.getInstrument().getDetector(iws).getPos()
        target_det_i_pos = target_ws.getInstrument().getDetector(iws).getPos()
        np.testing.assert_allclose(source_det_i_pos, target_det_i_pos,
                                   err_msg=f'Mismatch is detected at Detector {iws}')
    # Check source position
    source_moderator_pos = source_ws.getInstrument().getSource().getPos()
    target_moderator_pos = target_ws.getInstrument().getSource().getPos()
    np.testing.assert_allclose(source_moderator_pos, target_moderator_pos,
                               err_msg=f'Mismatch is detected at neutron source position')

    # Compare meta data
    assert len(source_ws.getRun().getProperties()) == len(target_ws.getRun().getProperties()), 'Meta data mismatch'


if __name__ == '__main__':
    pytest.main(__file__)
