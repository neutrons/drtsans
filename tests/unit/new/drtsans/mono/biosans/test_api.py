import pytest
import os
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
from mantid.simpleapi import CreateWorkspace, LoadHFIRSANS, LoadNexusProcessed
from drtsans.mono.biosans.api import (
    load_all_files,
    prepare_data_workspaces,
    prepare_data,
    process_single_configuration,
)
from drtsans.mono.biosans import reduction_parameters
from drtsans.mono.biosans.simulated_events import update_idf
from drtsans.load import load_events
from drtsans.samplelogs import SampleLogs
from drtsans.settings import unique_workspace_dundername as uwd

# standard imports
from unittest.mock import patch as mock_patch


@pytest.mark.skipif(
    not os.path.exists("/HFIR/CG3/IPTS-23782/nexus/CG3_960.nxs.h5"),
    reason="Required data is not available",
)
def test_load_all_files_simple():
    reduction_input = {
        "instrumentName": "CG3",
        "iptsNumber": "23782",
        "sample": {"runNumber": "960", "transmission": {"runNumber": ""}},
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "emptyTransmission": {"runNumber": ""},
        "beamCenter": {"runNumber": "960"},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "BIOSANS", validate=False)
    loaded = load_all_files(reduction_input)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].getHistory()
    assert history.size() == 9
    assert history.getAlgorithm(0).name() == "LoadEventNexus"
    assert history.getAlgorithm(0).getProperty("Filename").value == "/HFIR/CG3/IPTS-23782/nexus/CG3_960.nxs.h5"
    assert history.getAlgorithm(1).name() == "MoveInstrumentComponent"  # moderator
    assert history.getAlgorithm(2).name() == "MoveInstrumentComponent"  # sample-position
    assert history.getAlgorithm(3).name() == "MoveInstrumentComponent"  # detector1
    assert history.getAlgorithm(4).name() == "AddSampleLogMultiple"  # CG2:CS:SampleToSi
    assert history.getAlgorithm(5).name() == "AddSampleLogMultiple"  # sample_detector_distance
    assert history.getAlgorithm(6).name() == "HFIRSANS2Wavelength"
    assert history.getAlgorithm(7).name() == "SetUncertainties"
    assert history.getAlgorithm(8).name() == "AddSampleLogMultiple"  # sample_offset

    assert loaded.background is None
    assert loaded.background_transmission is None
    assert str(loaded.center) == str(loaded.sample[0])
    assert loaded.empty is None
    assert loaded.sample_transmission is None
    assert loaded.blocked_beam is None
    assert loaded.dark_current_main is None
    assert loaded.dark_current_wing is None
    # assert loaded.dark_current_midrange is None, not implemented
    assert loaded.sensitivity_main is None
    assert loaded.sensitivity_wing is None
    # assert loaded.sensitivity_midrange is None , not implemented
    assert loaded.mask is None


@pytest.mark.parametrize("generic_workspace", [{"name": "ws_raw_histo"}], indirect=True)
def test_prepare_data_workspaces_simple(generic_workspace):
    ws = generic_workspace  # friendly name

    output = prepare_data_workspaces(ws)
    # this should make a clone of the workspace
    assert ws is not output
    # and change the workspace name automatically
    assert ws.name() == "ws_raw_histo"
    assert output.name() == "ws_processed_histo"

    output2 = prepare_data_workspaces(ws, output_workspace="foobar")
    # the ws name should change to what is set
    assert ws.name() == "ws_raw_histo"
    assert output2.name() == "foobar"


def test_prepare_data_workspaces_center(biosans_f):
    ws = LoadHFIRSANS(Filename=biosans_f["beamcenter"])

    output = prepare_data_workspaces(ws, center_x=0.111, center_y=0.123, center_y_wing=0.222, solid_angle=False)

    # this should make a clone of the workspace
    assert ws is not output
    # and change the workspace name automatically
    assert output.name() == "ws_processed_histo"

    history = output.getHistory()
    assert history.size() == 4
    # There are 2 calls to MoveInstrumentComponent
    alg2 = history.getAlgorithm(2)
    assert alg2.name() == "MoveInstrumentComponent"
    assert alg2.getPropertyValue("ComponentName") == "detector1"
    assert alg2.getPropertyValue("RelativePosition") == "1"
    assert alg2.getPropertyValue("X") == "-0.111"
    assert alg2.getPropertyValue("Y") == "-0.123"

    alg3 = history.getAlgorithm(3)
    assert alg3.name() == "MoveInstrumentComponent"
    assert alg3.getPropertyValue("ComponentName") == "wing_detector"
    assert alg3.getPropertyValue("RelativePosition") == "1"
    assert alg3.getPropertyValue("X") == "0"
    assert alg3.getPropertyValue("Y") == "-0.222"


def test_prepare_data_workspaces_center_midrange_success(reference_dir):
    # similar test to test_prepare_data_workspaces_center
    # load the file and add the midrange detector
    # generate the output workspace from prepare_data_workspaces
    # with center parameters for main, wing and midrange detectors
    # check the algorithm history to ensure instrument components were moved with the requested coordinates

    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.new.biosans, overwrite_instrument=True)
    assert ws.getInstrument().getComponentByName("midrange_detector") is None
    # add the midrange detector
    ws = update_idf(ws)
    assert ws.getInstrument().getComponentByName("midrange_detector")

    output = prepare_data_workspaces(
        ws, center_x=0.111, center_y=0.123, center_y_wing=0.222, center_y_midrange=0.112, solid_angle=False
    )

    # this should make a clone of the workspace
    assert ws is not output
    # and the output workspace name is changed automatically from above
    assert output.name() == ws.name() + "_processed_histo"

    history = output.getHistory()
    assert history.size() == ws.getHistory().size() + 3 + 1
    # There are: 1 call to CloneWorkspace and 3 calls to MoveInstrumentComponent

    alg2 = history.getAlgorithm(4)
    assert alg2.name() == "MoveInstrumentComponent"
    assert alg2.getPropertyValue("ComponentName") == "detector1"
    assert alg2.getPropertyValue("RelativePosition") == "1"
    assert alg2.getPropertyValue("X") == "-0.111"
    assert alg2.getPropertyValue("Y") == "-0.123"

    alg3 = history.getAlgorithm(5)
    assert alg3.name() == "MoveInstrumentComponent"
    assert alg3.getPropertyValue("ComponentName") == "wing_detector"
    assert alg3.getPropertyValue("RelativePosition") == "1"
    assert alg3.getPropertyValue("X") == "0"
    assert alg3.getPropertyValue("Y") == "-0.222"

    alg4 = history.getAlgorithm(6)
    assert alg4.name() == "MoveInstrumentComponent"
    assert alg4.getPropertyValue("ComponentName") == "midrange_detector"
    assert alg4.getPropertyValue("RelativePosition") == "1"
    assert alg4.getPropertyValue("X") == "0"
    assert alg4.getPropertyValue("Y") == "-0.112"


def test_prepare_data_workspaces_center_midrange_failure(reference_dir):
    # similar test to test_prepare_data_workspaces_center_midrange_success
    # midrange center is required, but not passed
    # results to failure to move the instrument components

    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.new.biosans, overwrite_instrument=True)
    assert ws.getInstrument().getComponentByName("midrange_detector") is None
    # add the midrange detector
    ws = update_idf(ws)
    assert ws.getInstrument().getComponentByName("midrange_detector")

    output = prepare_data_workspaces(
        ws, center_x=0.111, center_y=0.123, center_y_wing=0.222, center_y_midrange=None, solid_angle=False
    )

    # this should make a clone of the workspace
    assert ws is not output
    # and the output workspace name is changed automatically from above
    assert output.name() == ws.name() + "_processed_histo"

    history = output.getHistory()
    assert history.size() == ws.getHistory().size() + 1
    # There is only 1 call to CloneWorkspace
    assert history.getAlgorithm(history.size() - 1).name() == "CloneWorkspace"


def test_prepare_data_center(reference_dir):
    # similar test to test_prepare_data_workspaces_center
    # load the file
    # generate the output workspace from prepare_data with center parameters for main and wing detectors
    # check the algorithm history to ensure instrument components were moved with the requested coordinates
    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.new.biosans, overwrite_instrument=True)

    output = prepare_data(str(ws), center_x=0.111, center_y=0.123, center_y_wing=0.222, solid_angle=False)

    # this should make a clone of the workspace
    assert ws is not output

    history = output.getHistory()
    # There are 2 calls to MoveInstrumentComponent
    alg2 = history.getAlgorithm(6)
    assert alg2.name() == "MoveInstrumentComponent"
    assert alg2.getPropertyValue("ComponentName") == "detector1"
    assert alg2.getPropertyValue("RelativePosition") == "1"
    assert alg2.getPropertyValue("X") == "-0.111"
    assert alg2.getPropertyValue("Y") == "-0.123"

    alg3 = history.getAlgorithm(7)
    assert alg3.name() == "MoveInstrumentComponent"
    assert alg3.getPropertyValue("ComponentName") == "wing_detector"
    assert alg3.getPropertyValue("RelativePosition") == "1"
    assert alg3.getPropertyValue("X") == "0"
    assert alg3.getPropertyValue("Y") == "-0.222"


@mock_patch("drtsans.load.__monitor_counts")
@mock_patch("drtsans.load.LoadEventNexus")
def test_prepare_data_center_midrange_success(mock_LoadEventNexus, mock_monitor_counts, reference_dir):
    # similar test to test_prepare_data_workspaces_center
    # load the file with mock patch
    # generate the output workspace from prepare_data with center parameters for main, wing and midrange detectors
    # check the algorithm history to ensure instrument components were moved with the requested coordinates

    output_workspace = "CG3_92300"
    synthetics_datasets = os.path.join(reference_dir.new.biosans, "synthetic_dataset")
    synthetics_data_path = os.path.join(synthetics_datasets, f"{output_workspace}.nxs.h5")

    mock_LoadEventNexus.return_value = LoadNexusProcessed(synthetics_data_path, OutputWorkspace=output_workspace)
    mock_monitor_counts.return_value = 42
    previous_history = mock_LoadEventNexus.return_value.getHistory().size()

    output = prepare_data(
        synthetics_data_path,
        data_dir=synthetics_datasets,
        output_workspace=output_workspace,
        center_x=0.111,
        center_y=0.123,
        center_y_wing=0.222,
        center_y_midrange=0.112,
        solid_angle=False,
    )

    history = output.getHistory()
    assert history.size() == previous_history + 5 + 3

    # There are 3 calls to MoveInstrumentComponent
    alg2 = history.getAlgorithm(14)
    assert alg2.name() == "MoveInstrumentComponent"
    assert alg2.getPropertyValue("ComponentName") == "detector1"
    assert alg2.getPropertyValue("RelativePosition") == "1"
    assert alg2.getPropertyValue("X") == "-0.111"
    assert alg2.getPropertyValue("Y") == "-0.123"

    alg3 = history.getAlgorithm(15)
    assert alg3.name() == "MoveInstrumentComponent"
    assert alg3.getPropertyValue("ComponentName") == "wing_detector"
    assert alg3.getPropertyValue("RelativePosition") == "1"
    assert alg3.getPropertyValue("X") == "0"
    assert alg3.getPropertyValue("Y") == "-0.222"

    alg4 = history.getAlgorithm(16)
    assert alg4.name() == "MoveInstrumentComponent"
    assert alg4.getPropertyValue("ComponentName") == "midrange_detector"
    assert alg4.getPropertyValue("RelativePosition") == "1"
    assert alg4.getPropertyValue("X") == "0"
    assert alg4.getPropertyValue("Y") == "-0.112"


@mock_patch("drtsans.load.__monitor_counts")
@mock_patch("drtsans.load.LoadEventNexus")
def test_prepare_data_center_midrange_failure(mock_LoadEventNexus, mock_monitor_counts, reference_dir):
    # similar test to test_prepare_data_center_midrange_success
    # midrange center is required, but not passed
    # results to failure to move the instrument components

    output_workspace = "CG3_92300"
    synthetics_datasets = os.path.join(reference_dir.new.biosans, "synthetic_dataset")
    synthetics_data_path = os.path.join(synthetics_datasets, f"{output_workspace}.nxs.h5")

    mock_LoadEventNexus.return_value = LoadNexusProcessed(synthetics_data_path, OutputWorkspace=output_workspace)
    mock_monitor_counts.return_value = 42
    previous_history = mock_LoadEventNexus.return_value.getHistory().size()

    output = prepare_data(
        synthetics_data_path,
        data_dir=synthetics_datasets,
        output_workspace=output_workspace,
        center_x=0.111,
        center_y=0.123,
        center_y_wing=0.222,
        solid_angle=False,
    )

    history = output.getHistory()
    assert history.size() == previous_history + 5


def test_prepare_data_workspaces_dark_current():
    # Create dark current workspace, insert the duration of the dark
    # current run as one of the log entries in the dark current
    # workspace.
    dark_current_workspace = uwd()  # arbitrary name for the dark current workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.full(2, 100.0),
        DataE=np.full(2, 10.0),
        NSpec=2,
        OutputWorkspace=dark_current_workspace,
    )
    SampleLogs(dark_current_workspace).insert("duration", 3600.0, "second")

    # Create a sample run workspace.
    data_workspace = uwd()  # arbitrary name for the sample workspace
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

    assert output.getHistory().size() == 8

    assert_almost_equal(output.extractY(), [[0], [1]])
    assert_almost_equal(output.extractE(), [[1.00498756], [1.41774469]])


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 2], [3, 4]]}], indirect=True)
def test_prepare_data_workspaces_flux_method(generic_workspace):
    ws = generic_workspace  # friendly name
    SampleLogs(ws).insert("duration", 2.0)
    SampleLogs(ws).insert("monitor", 2e9)

    # No normalization
    output = prepare_data_workspaces(ws, flux_method=None, solid_angle=False)
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


def test_prepare_data_workspaces_apply_mask(generic_workspace):
    ws = generic_workspace

    # mask_ws
    output = prepare_data_workspaces(ws, mask_ws=[0, 2], solid_angle=False)
    history = output.getHistory()
    assert history.size() == 4
    alg3 = history.getAlgorithm(3)
    assert alg3.name() == "MaskDetectors"
    assert alg3.getPropertyValue("DetectorList") == "0,2"


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 1], [1, 1]]}], indirect=True)
def test_prepare_data_workspaces_solid_angle(generic_workspace):
    ws = generic_workspace  # friendly name

    # No normalization
    output = prepare_data_workspaces(ws, solid_angle=True)
    # CreateWorkspace, LoadInstrument, CloneWorkspace, CloneWorkspace,
    # ClearMaskFlag, SolidAngle, Divide, ReplaceSpecialValues
    assert output.getHistory().size() == 8
    assert_almost_equal(output.extractY(), [[25.6259267], [25.6259267], [25.6259267], [25.6259267]])


def test_prepare_data_workspaces_sensitivity():
    # Create dark current workspace, insert the duration of the dark
    # current run as one of the log entries in the dark current
    # workspace.
    sensitivity_workspace = uwd()  # arbitrary name for the dark current workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.full(2, 2.0),
        DataE=np.full(2, np.sqrt(2)),
        NSpec=2,
        OutputWorkspace=sensitivity_workspace,
    )

    # Create a sample run workspace.
    data_workspace = uwd()  # arbitrary name for the sample workspace
    CreateWorkspace(
        DataX=[2.5, 3.5],
        DataY=np.array([1.0, 2.0]),
        DataE=np.array([1.0, np.sqrt(2)]),
        NSpec=2,
        OutputWorkspace=data_workspace,
    )

    output = prepare_data_workspaces(data_workspace, sensitivity_workspace=sensitivity_workspace, solid_angle=False)

    assert output.getHistory().size() == 6

    assert_almost_equal(output.extractY(), [[0.5], [1.0]])
    assert_almost_equal(output.extractE(), [[0.6123724], [1.0]])


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 2], [3, 4]]}], indirect=True)
def test_process_single_configuration_thickness_absolute_scale(generic_workspace):
    ws = generic_workspace

    # This should only run prepare_data_workspaces,
    # normalize_by_thickness and scale by absolute_scale
    # The output result should be scaled by y_out = y_in * absolute_scale / thickness

    output, trans = process_single_configuration(ws, solid_angle=False)

    # CreateWorkspace, LoadInstrument, CloneWorkspace,
    # CreateSingleValuedWorkspace, Divide,
    # CreateSingleValuedWorkspace, Multiply
    assert output.getHistory().size() == 7

    assert_equal(output.extractY(), [[1], [2], [3], [4]])
    assert not trans["sample"]
    assert not trans["background"]

    output, _ = process_single_configuration(ws, solid_angle=False, absolute_scale=1.5)
    assert_equal(output.extractY(), [[1.5], [3], [4.5], [6]])

    output, _ = process_single_configuration(ws, solid_angle=False, thickness=0.1)
    assert_equal(output.extractY(), [[10], [20], [30], [40]])

    output, _ = process_single_configuration(ws, solid_angle=False, absolute_scale=1.5, thickness=0.1)
    assert_equal(output.extractY(), [[15], [30], [45], [60]])


if __name__ == "__main__":
    pytest.main([__file__])
