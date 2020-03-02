import pytest

from mantid import mtd
from mantid.simpleapi import GroupWorkspaces

from drtsans.mono.biosans import load_histogram, load_events, transform_to_wavelength, merge_data
from drtsans.samplelogs import SampleLogs


def test_load_events(reference_dir):
    events_workspace = load_events('CG3_961', data_dir=reference_dir.new.biosans, overwrite_instrument=True)
    assert events_workspace.name() == 'BIOSANS_961'

    sample_logs = SampleLogs(events_workspace)
    assert sample_logs.monitor.value == 19173627

    x, y, z = events_workspace.getInstrument().getComponentByName('detector1').getPos()
    assert -1000 * x == pytest.approx(sample_logs.single_value('detector_trans_Readback'), abs=0.001)
    assert z == pytest.approx(sample_logs.single_value('sample_detector_distance'), abs=0.001)


def test_transform_to_wavelength(reference_dir):
    workspace = load_events('CG3_961.nxs.h5', data_dir=reference_dir.new.biosans)
    workspace = transform_to_wavelength(workspace)
    assert workspace.getAxis(0).getUnit().caption() == 'Wavelength'


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


def test_merge_data(reference_dir):
    # Merge the same file twice
    workspace1 = load_events('CG3_961.nxs.h5', data_dir=reference_dir.new.biosans, output_workspace='workspace1')

    with pytest.raises(ValueError) as excinfo:
        merge_data('workspace1', "merged")
    assert "is not a Workspace2D" in str(excinfo.value)  # Should complain about wrong workspace type

    workspace1 = transform_to_wavelength(workspace1)
    workspace2 = load_events('CG3_960.nxs.h5', data_dir=reference_dir.new.biosans, output_workspace='workspace2')
    workspace2 = transform_to_wavelength(workspace2)

    sample_logs1 = SampleLogs(workspace1)
    sample_logs2 = SampleLogs(workspace2)

    merged_workspaces = merge_data([workspace1, workspace2],  output_workspace="merged")

    merged_sample_logs = SampleLogs(merged_workspaces)

    # Check monitor and duration increase as the sum
    assert sample_logs1.monitor.value == 19173627
    assert sample_logs2.monitor.value == 1039
    assert merged_sample_logs.monitor.value == 19173627 + 1039
    assert sample_logs1.duration.value == pytest.approx(1809.4842529296875, abs=1e-11)
    assert sample_logs2.duration.value == pytest.approx(0.0833325386047363, abs=1e-11)
    assert merged_sample_logs.duration.value == pytest.approx(1809.4842529296875 + 0.08333253860473633, abs=1e-11)

    # Check Time Series properties increase length
    assert sample_logs1.wavelength.size() == 692
    assert sample_logs2.wavelength.size() == 2
    assert merged_sample_logs.wavelength.size() == 692 + 2

    # Check integrated intensity increases as the total sum
    assert mtd[str(workspace1)].extractY().sum() == 11067715
    assert mtd[str(workspace2)].extractY().sum() == 1
    assert mtd[str(merged_workspaces)].extractY().sum() == 11067715 + 1

    # Test different input formats
    # List of workspace names
    merged_workspaces_2 = merge_data(["workspace1", "workspace2"],  output_workspace="merged2")
    assert SampleLogs(merged_workspaces_2).duration.value == pytest.approx(1809.4842529296875 + 0.08333253860473633,
                                                                           abs=1e-11)

    # Comma separated list of workspace space
    merged_workspaces_3 = merge_data("workspace1, workspace2",  output_workspace="merged3")
    assert SampleLogs(merged_workspaces_3).duration.value == pytest.approx(1809.4842529296875 + 0.08333253860473633,
                                                                           abs=1e-11)

    # Workspace group
    ws_group = GroupWorkspaces('workspace1, workspace2')
    merged_workspaces_4 = merge_data(ws_group,  output_workspace="merged4")
    assert SampleLogs(merged_workspaces_4).duration.value == pytest.approx(1809.4842529296875 + 0.08333253860473633,
                                                                           abs=1e-11)


if __name__ == '__main__':
    pytest.main([__file__])
