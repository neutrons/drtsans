import pytest
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
from mantid.simpleapi import CreateWorkspace
from drtsans.mono.gpsans.api import (
    load_all_files,
    prepare_data_workspaces,
    process_single_configuration,
)
from drtsans.mono.gpsans import reduction_parameters
from drtsans.samplelogs import SampleLogs


def test_load_all_files_simple(datarepo_dir):
    reduction_input = {
        "instrumentName": "CG2",
        "iptsNumber": "23801",
        "sample": {"runNumber": "1338", "transmission": {"runNumber": ""}},
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "emptyTransmission": {"runNumber": ""},
        "beamCenter": {"runNumber": "1338"},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "GPSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.gpsans)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].getHistory()
    assert history.size() == 8
    assert history.getAlgorithm(0).name() == "LoadEventAsWorkspace2D"
    assert history.getAlgorithm(0).getProperty("Filename").value.endswith("sans/hfir/gpsans/CG2_1338.nxs.h5")
    assert history.getAlgorithm(1).name() == "MoveInstrumentComponent"  # moderator
    assert history.getAlgorithm(2).name() == "MoveInstrumentComponent"  # sample-position
    assert history.getAlgorithm(3).name() == "MoveInstrumentComponent"  # detector1
    assert history.getAlgorithm(4).name() == "AddSampleLogMultiple"  # CG2:CS:SampleToSi
    assert history.getAlgorithm(5).name() == "AddSampleLogMultiple"  # sample_detector_distance
    assert history.getAlgorithm(6).name() == "SetUncertainties"
    assert history.getAlgorithm(7).name() == "AddSampleLogMultiple"  # sample_offset

    assert loaded.background is None
    assert loaded.background_transmission is None
    assert str(loaded.center) == str(loaded.sample[0])
    assert loaded.empty is None
    assert loaded.sample_transmission is None
    assert loaded.blocked_beam is None
    assert loaded.dark_current is None
    assert loaded.sensitivity is None
    assert loaded.mask is None


@pytest.mark.parametrize("generic_workspace", [{"name": "ws_raw_histo"}], indirect=True)
def test_prepare_data_workspaces_simple(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    output = prepare_data_workspaces(ws, solid_angle=False)
    clean_workspace(output)
    # this should make a clone of the workspace
    assert ws is not output
    # and change the workspace name automatically
    assert ws.name() == "ws_raw_histo"
    assert output.name() == "ws_processed_histo"

    # This should really do nothing except clone the workspace
    history = output.getHistory()
    assert history.size() == 3
    # first 2 algorithms is creating generic_workspace
    assert history.getAlgorithm(2).name() == "CloneWorkspace"

    # the ws name should change to what is set
    output2 = prepare_data_workspaces(ws, output_workspace_name="foobar", solid_angle=False)
    clean_workspace(output2)
    assert output2.name() == "foobar"


def test_prepare_data_workspaces_center(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    output = prepare_data_workspaces(ws, center_x=0.111, center_y=0.123, solid_angle=False)
    clean_workspace(output)

    # this should make a clone of the workspace
    assert ws is not output
    # and change the workspace name automatically
    assert output.name() == "GenericSANS_processed_histo"

    history = output.getHistory()
    assert history.size() == 4
    # first 2 algorithms is creating generic_workspace
    assert history.getAlgorithm(2).name() == "CloneWorkspace"
    # This is the call to MoveInstrumentComponent
    alg3 = history.getAlgorithm(3)
    assert alg3.name() == "MoveInstrumentComponent"
    assert alg3.getPropertyValue("ComponentName") == "detector1"
    assert alg3.getPropertyValue("RelativePosition") == "1"
    assert alg3.getPropertyValue("X") == "-0.111"
    assert alg3.getPropertyValue("Y") == "-0.123"


def test_prepare_data_workspaces_dark_current(clean_workspace, temp_workspace_name):
    # Create dark current workspace, insert the duration of the dark
    # current run as one of the log entries in the dark current
    # workspace.
    dark_current_workspace = temp_workspace_name()  # arbitrary name for the dark current workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.full(2, 100.0),
        DataE=np.full(2, 10.0),
        NSpec=2,
        OutputWorkspace=dark_current_workspace,
    )
    SampleLogs(dark_current_workspace).insert("duration", 3600.0, "second")

    # Create a sample run workspace.
    data_workspace = temp_workspace_name()  # arbitrary name for the sample workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.array([1.0, 2.0]),
        DataE=np.array([1.0, np.sqrt(2)]),
        NSpec=2,
        OutputWorkspace=data_workspace,
    )
    # Insert the duration of the sample run. The log key must be the
    # same as that used for the dark current, which turns out to be
    # 'duration'
    SampleLogs(data_workspace).insert("duration", 36.0, "second")

    output = prepare_data_workspaces(data_workspace, dark_current=dark_current_workspace, solid_angle=False)
    clean_workspace(output)

    assert output.getHistory().size() == 8

    assert_almost_equal(output.extractY(), [[0], [1]])
    assert_almost_equal(output.extractE(), [[1.00498756], [1.41774469]])


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 2], [3, 4]]}], indirect=True)
def test_prepare_data_workspaces_flux_method(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)
    SampleLogs(ws).insert("duration", 2.0)
    SampleLogs(ws).insert("monitor", 2e9)

    # No normalization
    output = prepare_data_workspaces(ws, flux_method=None, solid_angle=False)
    clean_workspace(output)
    assert output.getHistory().size() == 3
    assert_almost_equal(output.extractY(), [[1], [2], [3], [4]])

    # Normalize by time
    output = prepare_data_workspaces(ws, flux_method="time", solid_angle=False)
    assert output.getHistory().size() == 5
    assert_almost_equal(output.extractY(), [[0.5], [1], [1.5], [2]])

    # Normalize by monitor, should scale by 1e8/(monitor counts)
    output = prepare_data_workspaces(ws, flux_method="monitor", solid_angle=False)
    assert output.getHistory().size() == 5
    assert_almost_equal(output.extractY(), [[0.05], [0.1], [0.15], [0.2]])


def test_prepare_data_workspaces_apply_mask(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)

    # mask_ws
    output = prepare_data_workspaces(ws, mask_ws=[0, 2], solid_angle=False)
    clean_workspace(output)
    history = output.getHistory()
    assert history.size() == 4
    alg3 = history.getAlgorithm(3)
    assert alg3.name() == "MaskDetectors"
    assert alg3.getPropertyValue("DetectorList") == "0,2"


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 1], [1, 1]]}], indirect=True)
def test_prepare_data_workspaces_solid_angle(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # No normalization
    output = prepare_data_workspaces(ws, solid_angle=True)
    clean_workspace(output)
    # CreateWorkspace, LoadInstrument, CloneWorkspace, CloneWorkspace,
    # ClearMaskFlag, SolidAngle, Divide, ReplaceSpecialValues
    assert output.getHistory().size() == 8
    assert_almost_equal(output.extractY(), [[25.6259267], [25.6259267], [25.6259267], [25.6259267]])


def test_prepare_data_workspaces_sensitivity(clean_workspace, temp_workspace_name):
    # Create dark current workspace, insert the duration of the dark
    # current run as one of the log entries in the dark current
    # workspace.
    sensitivity_workspace = temp_workspace_name()  # arbitrary name for the dark current workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.full(2, 2.0),
        DataE=np.full(2, np.sqrt(2)),
        NSpec=2,
        OutputWorkspace=sensitivity_workspace,
    )

    # Create a sample run workspace.
    data_workspace = temp_workspace_name()  # arbitrary name for the sample workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.array([1.0, 2.0]),
        DataE=np.array([1.0, np.sqrt(2)]),
        NSpec=2,
        OutputWorkspace=data_workspace,
    )

    output = prepare_data_workspaces(data_workspace, sensitivity_workspace=sensitivity_workspace, solid_angle=False)
    clean_workspace(output)

    assert output.getHistory().size() == 6

    assert_almost_equal(output.extractY(), [[0.5], [1.0]])
    assert_almost_equal(output.extractE(), [[0.6123724], [1.0]])


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 2], [3, 4]]}], indirect=True)
def test_process_single_configuration_thickness_absolute_scale(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)

    # This should only run prepare_data_workspaces,
    # normalize_by_thickness and scale by absolute_scale
    # The output result should be scaled by y_out = y_in * absolute_scale / thickness

    output = process_single_configuration(ws, solid_angle=False, debug=False)
    clean_workspace(output)

    # CreateWorkspace, LoadInstrument, CloneWorkspace,
    # CreateSingleValuedWorkspace, Divide,
    # CreateSingleValuedWorkspace, Multiply
    assert output.getHistory().size() == 7

    assert_equal(output.extractY(), [[1], [2], [3], [4]])

    output = process_single_configuration(ws, solid_angle=False, absolute_scale=1.5)
    assert_equal(output.extractY(), [[1.5], [3], [4.5], [6]])

    output = process_single_configuration(ws, solid_angle=False, thickness=0.1)
    assert_equal(output.extractY(), [[10], [20], [30], [40]])

    output = process_single_configuration(ws, solid_angle=False, absolute_scale=1.5, thickness=0.1)
    assert_equal(output.extractY(), [[15], [30], [45], [60]])


if __name__ == "__main__":
    pytest.main([__file__])
