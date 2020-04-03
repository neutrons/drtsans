from pathlib import Path
import pytest
import numpy as np
from drtsans.settings import unique_workspace_dundername as uwd
from drtsans.samplelogs import SampleLogs
from drtsans.mono.gpsans.load import load_histogram
from drtsans.mono.load import load_and_split


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


def test_load_and_split(reference_dir):
    # Check that is fails with missing required parameters
    with pytest.raises(ValueError) as excinfo:
        load_and_split('CG2_9177 ', data_dir=reference_dir.new.biosans,
                       sample_to_si_name='CG2:CS:SampleToSi',
                       si_nominal_distance=0.)
    assert "Must provide with time_interval or log_name and log_value_interval" == str(excinfo.value)

    filtered_ws = load_and_split('CG2_9177.', data_dir=reference_dir.new.biosans, time_interval=1000,
                                 sample_to_si_name='CG2:CS:SampleToSi',
                                 si_nominal_distance=0.)

    print(filtered_ws.size())
    assert filtered_ws.size() == 2


if __name__ == '__main__':
    pytest.main([__file__])
