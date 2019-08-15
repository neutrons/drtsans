import pytest

from mantid import mtd
from ornl.sans.hfir.biosans import load
from ornl.sans.samplelogs import SampleLogs


def test_api_load(biosans_f):
    ws = load(filename=biosans_f['beamcenter'])
    assert ws.name() == "BioSANS_exp402_scan0006_0001"

    # check logs
    sl = SampleLogs(ws)
    wavelength_log = sl.find_log_with_units('wavelength', 'Angstrom')
    assert round(wavelength_log) == 6.0
    wavelength_spread_log = sl.find_log_with_units(
        'wavelength-spread', 'Angstrom')
    assert round(wavelength_spread_log, 1) == 0.8
    wavelength_spread_ratio_log = sl.single_value('wavelength-spread-ratio')
    assert wavelength_spread_ratio_log == pytest.approx(0.1323, abs=1e-3)

    ws_name = "xpto"
    ws = load(filename=biosans_f['beamcenter'], output_workspace=ws_name)
    assert ws.name() == ws_name
    assert ws_name in mtd.getObjectNames()

    ws_name = "xptoxpto"
    ws = load(filename=biosans_f['beamcenter'], output_workspace=ws_name,
              wavelength=12, wavelength_spread=1,
              sample_to_detector_distance=9)
    assert ws_name in mtd.getObjectNames()

    # check logs when some parameters don't come directly from the metadata
    sl = SampleLogs(ws)
    wavelength_log = sl.find_log_with_units('wavelength', 'Angstrom')
    assert round(wavelength_log) == 12.0
    wavelength_spread_log = sl.find_log_with_units(
        'wavelength-spread', 'Angstrom')
    assert wavelength_spread_log == 1.0
    wavelength_spread_ratio_log = sl.single_value('wavelength-spread-ratio')
    assert wavelength_spread_ratio_log == pytest.approx(
        wavelength_spread_log / wavelength_log, abs=1e-3)
