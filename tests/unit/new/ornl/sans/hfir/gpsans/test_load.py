from os.path import join as pj
import pytest
import re
import numpy as np
from ornl.path import exists
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.hfir.gpsans.load import (permute_contents, serve_data_file,
                                        load_histogram)


@pytest.fixture(scope='module')
def tfs(refd):
    return dict(original=pj(refd.new.gpsans, 'CG2_exp206_scan0016_0001.xml'),
                idf=pj(refd.new.gpsans, 'instrument',
                       'GPSANS_Definition_2019_2100.xml'))


def test_permute_contents(tfs):
    xml = permute_contents(tfs['original'])
    # location in the XML file of the fourth pixel containing two counts
    assert [m.regs[0][0] for m in list(re.finditer('2\t', xml))][4] == 29800


def test_serve_data_file(tfs):
    with serve_data_file(tfs['original']) as data_file:
        assert data_file == tfs['original']
    with serve_data_file(tfs['original'], tfs['idf']) as data_file:
        assert exists(data_file) is True
        assert data_file != tfs['original']
    assert exists(data_file) is False


def test_load_histogram(tfs):
    name = uwd()
    w = load_histogram(tfs['original'], wavelength=4.2, wavelength_spread=0.42,
                       sample_to_detector_distance=1042, idf=None,
                       output_workspace=name)
    assert w.name() == name
    z = w.getInstrument().getComponentByName('detector1').getPos()[-1]
    assert z == 1.042
    max_id = np.argmax(w.extractY()[2:].flatten())
    tube_id, pixel_id = max_id // 256, max_id % 256
    assert (tube_id, pixel_id) == (92, 120)
    # Load with the new IDF
    w = load_histogram(tfs['original'], wavelength=4.2, wavelength_spread=0.42,
                       sample_to_detector_distance=1042, idf=tfs['idf'],
                       output_workspace=name)
    z = w.getInstrument().getComponentByName('detector1').getPos()[-1]
    assert z == 1.042
    max_id = np.argmax(w.extractY()[2:].flatten())
    assert max_id // 256 == 90
    assert max_id % 256 == pixel_id


if __name__ == '__main__':
    pytest.main()
