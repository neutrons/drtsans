import pytest

from mantid import mtd

from drtsans.mono.biosans import load_histogram, load_events, transform_to_wavelength
from drtsans.samplelogs import SampleLogs


def test_load_events(reference_dir):
    # default workspace name is file hint
    events_workspace = load_events('CG3_961.nxs.h5', data_dir=reference_dir.new.biosans,
                                   overwrite_instrument=True)

    assert events_workspace.name() == 'BIOSANS_961'

    sample_logs = SampleLogs(events_workspace)
    assert sample_logs.monitor.value == 19173627

    x, y, z = events_workspace.getInstrument().getComponentByName('detector1').getPos()
    assert -1000 * x == pytest.approx(sample_logs.single_value('detector_trans_Readback'), abs=0.001)
    assert z == pytest.approx(sample_logs.single_value('sample_detector_distance'), abs=0.001)


def test_transform_to_wavelength(reference_dir):
    workspace = load_events('CG3_961.nxs.h5', data_dir=reference_dir.new.biosans)
    workspace = transform_to_wavelength(workspace)
    assert workspace.getAxis(0).getUnit().caption() == 'wavelength'


def test_api_load(biosans_f):
    ws = load_histogram(filename=biosans_f['beamcenter'])
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
    ws = load_histogram(
        filename=biosans_f['beamcenter'], output_workspace=ws_name)
    assert ws.name() == ws_name
    assert ws_name in mtd.getObjectNames()

    ws_name = "xptoxpto"
    ws = load_histogram(filename=biosans_f['beamcenter'], output_workspace=ws_name,
                        wavelength=12, wavelength_spread=1, sample_det_cent=9)
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


if __name__ == '__main__':
    pytest.main([__file__])
