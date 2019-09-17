from pathlib import Path
import pytest
import numpy as np
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.samplelogs import SampleLogs
from drtsans.hfir.gpsans.load import load_histogram


def test_load_histogram(reference_dir):
    filename = str(Path(reference_dir.new.gpsans) / 'CG2_exp206_scan0016_0001.xml')

    name = uwd()
    w = load_histogram(filename, wavelength=4.2, wavelength_spread=0.42, output_workspace=name)
    assert w.name() == name
    z = w.getInstrument().getComponentByName('detector1').getPos()[-1]
    assert z == 6.055
    max_id = np.argmax(w.extractY()[2:].flatten())
    tube_id, pixel_id = max_id // 256, max_id % 256
    assert (tube_id, pixel_id) == (90, 120)
    assert SampleLogs(w)['sample-detector-distance'].value == pytest.approx(6055.0, abs=0.1)

    # override the the sample-to-detector-distance. Should not affect the detector with maximum intensity
    w = load_histogram(filename, wavelength=4.2, wavelength_spread=0.42,
                       sample_to_detector_distance=1.042, unit='m',  output_workspace=name)
    z = w.getInstrument().getComponentByName('detector1').getPos()[-1]
    assert z == 1.042
    max_id = np.argmax(w.extractY()[2:].flatten())
    assert (max_id // 256, max_id % 256) == (tube_id, pixel_id)
    assert SampleLogs(w)['sample-detector-distance'].value == pytest.approx(1042.0, abs=0.1)


if __name__ == '__main__':
    pytest.main()
