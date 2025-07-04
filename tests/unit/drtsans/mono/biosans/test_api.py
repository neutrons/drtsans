import pytest
from collections import namedtuple
import os
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
from mantid.simpleapi import CreateWorkspace, LoadHFIRSANS, LoadNexusProcessed
from drtsans.mono.biosans import reduction_parameters
from drtsans.mono.biosans.simulated_events import update_idf
from drtsans.load import load_events
import drtsans.plots.api

# from drtsans.plots.api import plot_IQmod, plot_IQazimuthal
from drtsans.samplelogs import SampleLogs
from os.path import join as path_join

from drtsans.dataobjects import IQazimuthal, IQmod, testing, I1DAnnular
from drtsans.mono.biosans.api import (
    load_all_files,
    prepare_data_workspaces,
    prepare_data,
    process_single_configuration,
    file_has_midrange_detector,
    save_i1d_all,
    plot_reduction_output,
    check_overlap_stitch_configuration,
    get_sample_detectordata,
)

# standard imports
from unittest.mock import patch as mock_patch


def test_load_all_files_simple(datarepo_dir):
    reduction_input = {
        "instrumentName": "CG3",
        "iptsNumber": "24665",
        "sample": {"runNumber": "1322", "transmission": {"runNumber": ""}},
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "emptyTransmission": {"runNumber": ""},
        "beamCenter": {"runNumber": "1322"},
        "configuration": {
            "useDefaultMask": False,
            "blockedBeamRunNumber": "",
        },
    }

    reduction_input = reduction_parameters(reduction_input, "BIOSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.biosans)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].getHistory()
    assert history.size() == 8
    assert history.getAlgorithm(0).name() == "LoadEventAsWorkspace2D"
    assert history.getAlgorithm(0).getProperty("Filename").value.endswith("sans/hfir/biosans/CG3_1322.nxs.h5")
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
    assert loaded.dark_current_main is None
    assert loaded.dark_current_wing is None
    # assert loaded.dark_current_midrange is None, not implemented
    assert loaded.sensitivity_main is None
    assert loaded.sensitivity_wing is None
    # assert loaded.sensitivity_midrange is None , not implemented
    assert loaded.mask is None


@pytest.mark.parametrize("generic_workspace", [{"name": "ws_raw_histo"}], indirect=True)
def test_prepare_data_workspaces_simple(clean_workspace, generic_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    output = prepare_data_workspaces(ws)
    clean_workspace(output)
    # this should make a clone of the workspace
    assert ws is not output
    # and change the workspace name automatically
    assert ws.name() == "ws_raw_histo"
    assert output.name() == "ws_processed_histo"

    output2 = prepare_data_workspaces(ws, output_workspace=clean_workspace("foobar"))
    # the ws name should change to what is set
    assert ws.name() == "ws_raw_histo"
    assert output2.name() == "foobar"


@pytest.mark.datarepo
def test_prepare_data_workspaces_center(biosans_f, clean_workspace):
    ws = LoadHFIRSANS(Filename=biosans_f["beamcenter"])
    clean_workspace(ws)

    output = prepare_data_workspaces(ws, center_x=0.111, center_y=0.123, center_y_wing=0.222, solid_angle=False)
    clean_workspace(output)

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


@pytest.mark.mount_eqsans
def test_prepare_data_workspaces_center_midrange_success(reference_dir, has_sns_mount):
    if not has_sns_mount:
        pytest.skip("Skip test if SNS mount is not found.")

    # similar test to test_prepare_data_workspaces_center
    # load the file and add the midrange detector
    # generate the output workspace from prepare_data_workspaces
    # with center parameters for main, wing and midrange detectors
    # check the algorithm history to ensure instrument components were moved with the requested coordinates

    # NOTE: this data is too large for datarepo, and should be run manually on analysis cluster
    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.biosans, overwrite_instrument=True)
    assert ws.getInstrument().getComponentByName("midrange_detector") is None
    # add the midrange detector
    ws = update_idf(ws)
    assert ws.getInstrument().getComponentByName("midrange_detector")

    # this should make a clone of the workspace
    output = prepare_data_workspaces(
        ws, center_x=0.111, center_y=0.123, center_y_wing=0.222, center_y_midrange=0.112, solid_angle=False
    )

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


@pytest.mark.mount_eqsans
def test_prepare_data_workspaces_center_midrange_failure(reference_dir, has_sns_mount):
    if not has_sns_mount:
        pytest.skip("Skip test if SNS mount is not found.")

    # similar test to test_prepare_data_workspaces_center_midrange_success
    # midrange center is required, but not passed
    # results to failure to move the instrument components

    # NOTE: this data is too large for datarepo, and should be run manually on analysis cluster
    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.biosans, overwrite_instrument=True)
    assert ws.getInstrument().getComponentByName("midrange_detector") is None
    # add the midrange detector
    ws = update_idf(ws)
    assert ws.getInstrument().getComponentByName("midrange_detector")

    # this should make a clone of the workspace
    output = prepare_data_workspaces(
        ws, center_x=0.111, center_y=0.123, center_y_wing=0.222, center_y_midrange=None, solid_angle=False
    )

    assert ws is not output
    # and the output workspace name is changed automatically from above
    assert output.name() == ws.name() + "_processed_histo"

    history = output.getHistory()
    assert history.size() == ws.getHistory().size() + 1
    # There is only 1 call to CloneWorkspace
    assert history.getAlgorithm(history.size() - 1).name() == "CloneWorkspace"


@pytest.mark.mount_eqsans
def test_prepare_data_workspaces_apply_mask_detectors_str(reference_dir, has_sns_mount):
    if not has_sns_mount:
        pytest.skip("Skip test if SNS mount is not found.")

    # load the file and add the midrange detector
    # generate the output workspace from prepare_data_workspaces
    # mask detector with a comma separated detector name-string

    # NOTE: this data is too large for datarepo, and should be run manually on analysis cluster
    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.biosans, overwrite_instrument=True)
    assert ws.getInstrument().getComponentByName("midrange_detector") is None
    # add the midrange detector
    ws = update_idf(ws)
    assert ws.getInstrument().getComponentByName("midrange_detector")

    # mask_ws
    mask_detectors = "detector1,wing_detector,midrange_detector"
    output = prepare_data_workspaces(ws, mask_detector=mask_detectors, solid_angle=False)
    history = output.getHistory()
    assert history.size() == 5
    alg = history.lastAlgorithm()
    assert alg.name() == "MaskDetectors"
    assert alg.getPropertyValue("ComponentList") == mask_detectors


@pytest.mark.mount_eqsans
def test_prepare_data_workspaces_apply_mask_detectors_lst(reference_dir, has_sns_mount):
    if not has_sns_mount:
        pytest.skip("Skip test if SNS mount is not found.")

    # load the file and add the midrange detector
    # generate the output workspace from prepare_data_workspaces
    # mask detector with a list of detector names

    # NOTE: this data is too large for datarepo, and should be run manually on analysis cluster
    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.biosans, overwrite_instrument=True)
    assert ws.getInstrument().getComponentByName("midrange_detector") is None
    # add the midrange detector
    ws = update_idf(ws)
    assert ws.getInstrument().getComponentByName("midrange_detector")

    # mask_ws
    mask_detectors = ["detector1", "midrange_detector"]
    output = prepare_data_workspaces(ws, mask_detector=mask_detectors, solid_angle=False)
    history = output.getHistory()
    assert history.size() == 5
    alg = history.lastAlgorithm()
    assert alg.name() == "MaskDetectors"
    assert alg.getPropertyValue("ComponentList") == ",".join(mask_detectors)


@pytest.mark.mount_eqsans
def test_prepare_data_center(reference_dir, has_sns_mount):
    if not has_sns_mount:
        pytest.skip("Skip test if SNS mount is not found.")

    # similar test to test_prepare_data_workspaces_center
    # load the file
    # generate the output workspace from prepare_data with center parameters for main and wing detectors
    # check the algorithm history to ensure instrument components were moved with the requested coordinates

    # NOTE: this data is too large for datarepo, and should be run manually on analysis cluster
    ws = load_events("CG3_957.nxs.h5", data_dir=reference_dir.biosans, overwrite_instrument=True)

    # this should make a clone of the workspace
    output = prepare_data(str(ws), center_x=0.111, center_y=0.123, center_y_wing=0.222, solid_angle=False)

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


@pytest.mark.datarepo
@mock_patch("drtsans.load.__monitor_counts")
@mock_patch("drtsans.load.LoadEventAsWorkspace2D")
def test_prepare_data_center_midrange_success(
    mock_LoadEventAsWorkspace2D, mock_monitor_counts, datarepo_dir, clean_workspace
):
    # similar test to test_prepare_data_workspaces_center
    # load the file with mock patch
    # generate the output workspace from prepare_data with center parameters for main, wing and midrange detectors
    # check the algorithm history to ensure instrument components were moved with the requested coordinates

    output_workspace = clean_workspace("CG3_92300")
    synthetics_datasets = os.path.join(datarepo_dir.biosans, "synthetic_dataset")
    synthetics_data_path = os.path.join(synthetics_datasets, f"{output_workspace}.nxs.h5")

    mock_LoadEventAsWorkspace2D.return_value = LoadNexusProcessed(
        synthetics_data_path, OutputWorkspace=output_workspace
    )
    mock_monitor_counts.return_value = 42
    previous_history = mock_LoadEventAsWorkspace2D.return_value.getHistory().size()

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
    alg2 = history.getAlgorithm(history.size() - 4)
    assert alg2.name() == "MoveInstrumentComponent"
    assert alg2.getPropertyValue("ComponentName") == "detector1"
    assert alg2.getPropertyValue("RelativePosition") == "1"
    assert alg2.getPropertyValue("X") == "-0.111"
    assert alg2.getPropertyValue("Y") == "-0.123"

    alg3 = history.getAlgorithm(history.size() - 3)
    assert alg3.name() == "MoveInstrumentComponent"
    assert alg3.getPropertyValue("ComponentName") == "wing_detector"
    assert alg3.getPropertyValue("RelativePosition") == "1"
    assert alg3.getPropertyValue("X") == "0"
    assert alg3.getPropertyValue("Y") == "-0.222"

    alg4 = history.getAlgorithm(history.size() - 2)
    assert alg4.name() == "MoveInstrumentComponent"
    assert alg4.getPropertyValue("ComponentName") == "midrange_detector"
    assert alg4.getPropertyValue("RelativePosition") == "1"
    assert alg4.getPropertyValue("X") == "0"
    assert alg4.getPropertyValue("Y") == "-0.112"


@pytest.mark.datarepo
@mock_patch("drtsans.load.__monitor_counts")
@mock_patch("drtsans.load.LoadEventAsWorkspace2D")
def test_prepare_data_center_midrange_failure(
    mock_LoadEventAsWorkspace2D, mock_monitor_counts, datarepo_dir, clean_workspace
):
    # similar test to test_prepare_data_center_midrange_success
    # midrange center is required, but not passed
    # results to failure to move the instrument components

    output_workspace = clean_workspace("CG3_92300")
    synthetics_datasets = os.path.join(datarepo_dir.biosans, "synthetic_dataset")
    synthetics_data_path = os.path.join(synthetics_datasets, f"{output_workspace}.nxs.h5")

    mock_LoadEventAsWorkspace2D.return_value = LoadNexusProcessed(
        synthetics_data_path, OutputWorkspace=output_workspace
    )
    mock_monitor_counts.return_value = 42
    previous_history = mock_LoadEventAsWorkspace2D.return_value.getHistory().size()

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


@pytest.mark.datarepo
@mock_patch("drtsans.load.__monitor_counts")
@mock_patch("drtsans.load.LoadEventAsWorkspace2D")
def test_prepare_data_apply_mask_detectors_lst(
    mock_LoadEventAsWorkspace2D, mock_monitor_counts, datarepo_dir, clean_workspace
):
    # load the file with mock patch
    # generate the output workspace from prepare_data
    # mask detector with a list of detector names

    output_workspace = clean_workspace("CG3_92300")
    synthetics_datasets = os.path.join(datarepo_dir.biosans, "synthetic_dataset")
    synthetics_data_path = os.path.join(synthetics_datasets, f"{output_workspace}.nxs.h5")

    mock_LoadEventAsWorkspace2D.return_value = LoadNexusProcessed(
        synthetics_data_path, OutputWorkspace=output_workspace
    )
    mock_monitor_counts.return_value = 42

    # mask_ws
    mask_detectors = ["wing_detector", "midrange_detector"]
    output = prepare_data(
        synthetics_data_path,
        data_dir=synthetics_datasets,
        output_workspace=output_workspace,
        mask_detector=mask_detectors,
        solid_angle=False,
    )

    history = output.getHistory()

    alg = history.getAlgorithm(history.size() - 2)
    assert alg.name() == "MaskDetectors"
    assert alg.getPropertyValue("ComponentList") == ",".join(mask_detectors)


@pytest.mark.datarepo
@mock_patch("drtsans.load.__monitor_counts")
@mock_patch("drtsans.load.LoadEventAsWorkspace2D")
def test_prepare_data_apply_mask_detectors_str(
    mock_LoadEventAsWorkspace2D, mock_monitor_counts, datarepo_dir, clean_workspace
):
    # load the file with mock patch
    # generate the output workspace from prepare_data
    # mask detector with a detector name

    output_workspace = clean_workspace("CG3_92300")
    synthetics_datasets = os.path.join(datarepo_dir.biosans, "synthetic_dataset")
    synthetics_data_path = os.path.join(synthetics_datasets, f"{output_workspace}.nxs.h5")

    mock_LoadEventAsWorkspace2D.return_value = LoadNexusProcessed(
        synthetics_data_path, OutputWorkspace=output_workspace
    )
    mock_monitor_counts.return_value = 42

    # mask_ws
    mask_detectors = "midrange_detector"
    output = prepare_data(
        synthetics_data_path,
        data_dir=synthetics_datasets,
        output_workspace=output_workspace,
        mask_detector=mask_detectors,
        solid_angle=False,
    )

    history = output.getHistory()

    alg = history.getAlgorithm(history.size() - 2)
    assert alg.name() == "MaskDetectors"
    assert alg.getPropertyValue("ComponentList") == mask_detectors


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
    clean_workspace(str(ws))
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
    clean_workspace(str(ws))

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
    clean_workspace(str(ws))

    # No normalization
    output = prepare_data_workspaces(ws, solid_angle=True)
    clean_workspace(output)
    # CreateWorkspace, LoadInstrument, CloneWorkspace, CloneWorkspace,
    # ClearMaskFlag, SolidAngle, Divide, ReplaceSpecialValues
    assert output.getHistory().size() == 8
    assert_almost_equal(output.extractY(), [[25.6259267], [25.6259267], [25.6259267], [25.6259267]])


def test_prepare_data_workspaces_sensitivity(temp_workspace_name, clean_workspace):
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
    clean_workspace(str(ws))

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
    clean_workspace(output)
    assert_equal(output.extractY(), [[1.5], [3], [4.5], [6]])

    output, _ = process_single_configuration(ws, solid_angle=False, thickness=0.1)
    assert_equal(output.extractY(), [[10], [20], [30], [40]])

    output, _ = process_single_configuration(ws, solid_angle=False, absolute_scale=1.5, thickness=0.1)
    assert_equal(output.extractY(), [[15], [30], [45], [60]])


@pytest.mark.skipif(
    not os.path.exists("/HFIR/HB2B/shared/autoreduce/"),
    reason="Skip test if HFIR mount is down.",
)
@pytest.mark.mount_eqsans
def test_has_midrange_detector(has_sns_mount):
    """Unit test for helper function has_midrange_detector."""
    if not has_sns_mount:
        pytest.skip("Skip test if SNS mount is not found.")

    reduction_input = {
        "instrumentName": "CG3",
        "iptsNumber": "23782,23783",
        "sample": {"runNumber": "960", "transmission": {"runNumber": ""}},
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "emptyTransmission": {"runNumber": ""},
        "beamCenter": {"runNumber": "960"},
        "configuration": {"useDefaultMask": False},
    }
    rst = file_has_midrange_detector(
        sample=reduction_input["sample"]["runNumber"],
        instrument_name=reduction_input["instrumentName"],
        ipts=reduction_input["iptsNumber"],
        directory=None,
    )
    assert not rst


@pytest.mark.parametrize(
    "iqmod_dummy, output_files",
    [
        # case when "1DQbinType" is "scalar" or "annular": 1 intensity profile per detector
        (
            [IQmod(intensity=[], error=[], mod_q=[], delta_mod_q=[])],
            ["output_1D_main.txt", "output_1D_midrange.txt", "output_1D_wing.txt", "output_1D_combined.txt"],
        ),
        # case when "1DQbinType" is "wedge": 2 intensity profiles per detector
        (
            [
                IQmod(intensity=[], error=[], mod_q=[], delta_mod_q=[]),
                IQmod(intensity=[], error=[], mod_q=[], delta_mod_q=[]),
            ],
            [
                "output_1D_main_wedge_0.txt",
                "output_1D_midrange_wedge_0.txt",
                "output_1D_wing_wedge_0.txt",
                "output_1D_combined_wedge_0.txt",
                "output_1D_main_wedge_1.txt",
                "output_1D_midrange_wedge_1.txt",
                "output_1D_wing_wedge_1.txt",
                "output_1D_combined_wedge_1.txt",
            ],
        ),
        # case with I1DAnnular
        (
            [I1DAnnular(intensity=[], error=[], phi=[])],
            ["output_1D_main.txt", "output_1D_midrange.txt", "output_1D_wing.txt", "output_1D_combined.txt"],
        ),
    ],
)
def test_save_i1d_all(tmp_path, iqmod_dummy, output_files):
    output_dir = path_join(tmp_path, "1D")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    save_i1d_all(iqmod_dummy, iqmod_dummy, iqmod_dummy, iqmod_dummy, "output", tmp_path, "", True)

    for filename in output_files:
        filepath = path_join(output_dir, filename)
        assert os.path.exists(filepath)


def test_plot_reduction_output(monkeypatch):
    """Unit test for helper function plot_reduction_output."""
    # Mock the plot_i1d function
    plot_i1d_counter = 0

    def mock_plot_i1d(*args, **kwargs):
        nonlocal plot_i1d_counter
        plot_i1d_counter += 1
        return "mock_plot_i1d"

    monkeypatch.setattr("drtsans.plots.api.plot_i1d", mock_plot_i1d)

    # Mock the plot_IQ function
    plot_IQazimuthal_counter = 0

    def mock_plot_IQazimuthal(*args, **kwargs):
        nonlocal plot_IQazimuthal_counter
        plot_IQazimuthal_counter += 1
        return "mock_plot_IQazimuthal"

    monkeypatch.setattr(drtsans.plots.api, "plot_IQazimuthal", mock_plot_IQazimuthal)
    monkeypatch.setattr(drtsans.plots.api, "plot_IQazimuthal", mock_plot_IQazimuthal)

    # Mock allow_overwrite
    allow_overwrite_counter = 0

    def mock_allow_overwrite(*args, **kwargs):
        nonlocal allow_overwrite_counter
        allow_overwrite_counter += 1
        return True

    monkeypatch.setattr(drtsans.path, "allow_overwrite", mock_allow_overwrite)

    # Case 1: with midrange detector
    reduction_input = {
        "instrumentName": "CG3",
        "iptsNumber": "23782",
        "sample": {"runNumber": "960", "transmission": {"runNumber": ""}},
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "emptyTransmission": {"runNumber": ""},
        "beamCenter": {"runNumber": "960"},
        "outputFileName": "tmp",
        "configuration": {
            "useDefaultMask": False,
            "blockedBeamRunNumber": "",
            "darkMainFileName": None,
            "darkWingFileName": None,
            "outputDir": "/tmp",
            "1DQbinType": "wedge",
            "wedges": "mock_wedges",
            "QminMain": 0.1,
            "QmaxMain": 0.2,
            "QminWing": 0.3,
            "QmaxWing": 0.4,
            "QminMidrange": 0.5,
            "QmaxMidrange": 0.6,
        },
        "has_midrange_detector": True,
    }
    # create a namedtuple to mock the reduction output
    output_1 = namedtuple(
        "ReductionOutput",
        [
            # plot_IQazimuthal
            "I2D_main",
            "I2D_wing",
            "I2D_midrange",
            # plot_IQmod
            "I1D_main",
            "I1D_wing",
            "I1D_midrange",
            "I1D_combined",
        ],
    )
    output_1.I1D_main = np.array([1, 2, 3])
    output_1.I1D_wing = np.array([4, 5, 6])
    output_1.I1D_midrange = np.array([7, 8, 9])
    output_1.I1D_combined = np.array([10, 11, 12])
    reduction_output = [output_1]
    plot_reduction_output(reduction_output, reduction_input)

    assert plot_IQazimuthal_counter == 3
    assert plot_i1d_counter == 3
    assert allow_overwrite_counter == 2


@pytest.mark.parametrize(
    "qbintype1d, qmin_name, qmax_name",
    [
        ("scalar", "overlapStitchQmin", "overlapStitchQmax"),
        ("wedge", "wedge1overlapStitchQmin", "wedge1overlapStitchQmax"),
        ("wedge", "wedge2overlapStitchQmin", "wedge2overlapStitchQmax"),
    ],
)
@pytest.mark.parametrize(
    "throws_error, has_midrange, ignore_midrange, reference_detector, qmin_value, qmax_value",
    [
        (True, True, False, "main", [0.01], [0.015]),  # too few Qmin/Qmax for run with midrange detector
        (
            True,
            True,
            True,
            "main",
            [0.01, 0.02],
            [0.015, 0.025],
        ),  # too many Qmin/Qmax for overlapStitchIgnoreMidrange = True
        (
            True,
            False,
            False,
            "main",
            [0.01, 0.02],
            [0.015, 0.025],
        ),  # too many Qmin/Qmax for run without midrange detector
        (True, False, False, "main", [0.01, 0.02, 0.03], [0.015, 0.025, 0.035]),  # lists too long
        (True, False, False, "midrange", None, None),  # midrange selected but not run with
        # valid inputs:
        (False, False, False, "main", None, None),
        (False, False, False, "wing", None, None),
        (False, True, False, "main", None, None),
        (False, False, False, "main", [], []),
        (False, True, False, "midrange", None, None),
        (False, True, False, "main", [], []),
        (False, False, False, "main", [0.01], [0.015]),  # run without midrange detector
        (False, True, False, "main", [0.01, 0.02], [0.015, 0.025]),  # run with midrange detector
        (False, True, True, "main", [0.01], [0.015]),  # run with midrange detector, but it is excluded from stitching
    ],
)
def test_check_overlap_stitch_configuration(
    qbintype1d,
    qmin_name,
    qmax_name,
    throws_error,
    has_midrange,
    ignore_midrange,
    reference_detector,
    qmin_value,
    qmax_value,
):
    """Unit test for helper function check_overlap_stitch_configuration."""
    reduction_config = {}
    reduction_config["1DQbinType"] = qbintype1d
    reduction_config["overlapStitchIgnoreMidrange"] = ignore_midrange
    reduction_config["overlapStitchReferenceDetector"] = reference_detector

    # first, set all overlap stitch parameters to None to only test the ones given by qmin_name and qmax_name
    reduction_config["overlapStitchQmin"] = None
    reduction_config["overlapStitchQmax"] = None
    reduction_config["wedge1overlapStitchQmin"] = None
    reduction_config["wedge1overlapStitchQmax"] = None
    reduction_config["wedge2overlapStitchQmin"] = None
    reduction_config["wedge2overlapStitchQmax"] = None

    # set the overlap stitch parameters given by qmin_name and qmax_name
    reduction_config[qmin_name] = qmin_value
    reduction_config[qmax_name] = qmax_value
    reduction_input = {}
    reduction_input["has_midrange_detector"] = has_midrange
    reduction_input["configuration"] = reduction_config
    if throws_error:
        with pytest.raises(ValueError) as error_info:
            check_overlap_stitch_configuration(reduction_input)
        assert "overlapStitch" in str(error_info.value)
    else:
        check_overlap_stitch_configuration(reduction_input)


@pytest.mark.parametrize("has_midrange", [True, False])
@pytest.mark.parametrize(
    "iqmod1d, iqmod2d",
    [
        # case when "1DQbinType" is "scalar" or "annular": one 1D intensity profile per detector
        (
            [
                IQmod(
                    intensity=[1.0, 1.0, 1.0],
                    error=[0.1, 0.1, 0.1],
                    mod_q=[0.1, 0.2, 0.3],
                    delta_mod_q=[0.01, 0.01, 0.01],
                )
            ],
            IQazimuthal(
                intensity=[1.0, 1.0, 1.0],
                error=[0.1, 0.1, 0.1],
                qx=[0.1, 0.2, 0.3],
                qy=[0.1, 0.2, 0.3],
                delta_qx=[0.01, 0.01, 0.01],
                delta_qy=[0.01, 0.01, 0.01],
            ),
        ),
        # case when "1DQbinType" is "wedge": two 1D intensity profiles per detector
        (
            [
                IQmod(
                    intensity=[1.0, 1.0, 1.0],
                    error=[0.1, 0.1, 0.1],
                    mod_q=[0.1, 0.2, 0.3],
                    delta_mod_q=[0.01, 0.01, 0.01],
                ),
                IQmod(
                    intensity=[2.0, 2.0, 2.0],
                    error=[0.1, 0.1, 0.1],
                    mod_q=[0.1, 0.2, 0.3],
                    delta_mod_q=[0.01, 0.01, 0.01],
                ),
            ],
            IQazimuthal(
                intensity=[1.0, 1.0, 1.0],
                error=[0.1, 0.1, 0.1],
                qx=[0.1, 0.2, 0.3],
                qy=[0.1, 0.2, 0.3],
                delta_qx=[0.01, 0.01, 0.01],
                delta_qy=[0.01, 0.01, 0.01],
            ),
        ),
    ],
)
def test_get_sample_detectordata(has_midrange, iqmod1d, iqmod2d):
    """Test of helper function get_sample_detectordata"""
    # create detector data
    detector_data = get_sample_detectordata(
        iqmod1d, iqmod1d, iqmod1d, iqmod1d, iqmod2d, iqmod2d, iqmod2d, has_midrange
    )

    # test the contents of detector data
    testing.assert_allclose(detector_data["combined"]["i1d"][0], iqmod1d[-1])
    testing.assert_allclose(detector_data["main_0"]["i1d"][0], iqmod1d[0])
    testing.assert_allclose(detector_data["wing_0"]["i1d"][0], iqmod1d[0])
    testing.assert_allclose(detector_data["main_0"]["iqxqy"], iqmod2d)
    testing.assert_allclose(detector_data["wing_0"]["iqxqy"], iqmod2d)
    if has_midrange:
        testing.assert_allclose(detector_data["midrange_0"]["i1d"][0], iqmod1d[0])
        testing.assert_allclose(detector_data["midrange_0"]["iqxqy"], iqmod2d)
    else:
        assert "midrange_0" not in detector_data

    if len(iqmod1d) == 2:  # case when "1DQbinType" is "wedge"
        testing.assert_allclose(detector_data["main_1"]["i1d"][0], iqmod1d[1])
        testing.assert_allclose(detector_data["wing_1"]["i1d"][0], iqmod1d[1])
        if has_midrange:
            testing.assert_allclose(detector_data["midrange_1"]["i1d"][0], iqmod1d[1])
        else:
            assert "midrange_1" not in detector_data


if __name__ == "__main__":
    pytest.main([__file__])
