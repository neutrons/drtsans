import os
import pytest
from collections import namedtuple
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
from mantid.simpleapi import CreateWorkspace, mtd

from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (
    load_all_files,
    prepare_data_workspaces,
    pre_process_single_configuration,
)
from drtsans.samplelogs import SampleLogs


ws_mon_pair = namedtuple("ws_mon_pair", ["data", "monitor"])


def test_load_all_files_simple_interval(datarepo_dir):
    """This test should test options FilterByTimeStart and FilterByTimeStop that can be passed to Mantid algorithms
    LoadEventNexus or LoadEvenAsWorkspace2D.
    """
    # set some filenames to /bin/true so that it will pass the validation without being replaced with default values
    specs = {
        "iptsNumber": 22747,
        "sample": {"runNumber": 105428, "thickness": 1.0},
        "beamCenter": {"runNumber": 105428},
        "outputFileName": "test",
        "dataDirectories": datarepo_dir.eqsans,
        "configuration": {
            "outputDir": "/tmp",
            "useDefaultMask": False,
            "maskFileName": "/bin/true",
            "darkFileName": "/bin/true",
            "sensitivityFileName": "/bin/true",
            "sampleOffset": "0",
            "LogQBinsPerDecade": 10,
            "normalization": "Total charge",
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": "5.18325",
            "instrumentConfigurationDir": os.path.join(datarepo_dir.eqsans, "instrument_configuration"),
        },
    }
    reduction_input = reduction_parameters(specs, instrument_name="EQSANS")

    # replace the /bin/true filenames with None
    reduction_input["configuration"]["maskFileName"] = None
    reduction_input["configuration"]["darkFileName"] = None
    reduction_input["configuration"]["sensitivityFileName"] = None

    start_time = 0.0
    end_time = 1.0
    reduction_input["sample"]["loadOptions"] = {"FilterByTimeStart": start_time, "FilterByTimeStop": end_time}

    loaded = load_all_files(reduction_input)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].data.getHistory()

    assert history.size() == 10
    assert history.getAlgorithm(0).name() == "LoadEventNexus"
    assert history.getAlgorithm(0).getProperty("Filename").value.endswith("sns/eqsans/EQSANS_105428.nxs.h5")
    assert history.getAlgorithm(2).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(3).name() == "SetInstrumentParameter"
    assert history.getAlgorithm(5).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(6).name() == "ConvertUnits"
    assert history.getAlgorithm(7).name() == "Rebin"
    assert history.getAlgorithm(8).name() == "SetUncertainties"
    assert history.getAlgorithm(9).name() == "AddSampleLogMultiple"

    assert loaded.background.data is None
    assert loaded.background_transmission.data is None
    assert loaded.empty.data is None
    assert loaded.sample_transmission.data is None
    assert loaded.dark_current.data is None
    assert loaded.sensitivity is None
    assert loaded.mask is None

    # Verify that if something is changed that it gets applied correctly on reload, use default mask as test
    # First check the current value
    assert not loaded.sample[0].data.detectorInfo().isMasked(1)

    # check interval
    w = loaded.sample[0].data
    assert int(w.extractY().sum()) == 773

    # Change reduction input and rerun load_all_files
    reduction_input["configuration"]["useDefaultMask"] = True
    reduction_input["configuration"]["defaultMask"] = "{'Pixel':'1'}"
    loaded = load_all_files(reduction_input)

    # Check that the value has changed
    assert loaded.sample[0].data.detectorInfo().isMasked(1)


def test_load_all_files_simple(datarepo_dir):
    # set some filenames to /bin/true so that it will pass the validation without being replaced with default values
    specs = {
        "iptsNumber": 22747,
        "sample": {"runNumber": 105428, "thickness": 1.0},
        "beamCenter": {"runNumber": 105428},
        "outputFileName": "test",
        "dataDirectories": datarepo_dir.eqsans,
        "configuration": {
            "outputDir": "/tmp",
            "useDefaultMask": False,
            "maskFileName": "/bin/true",
            "darkFileName": "/bin/true",
            "sensitivityFileName": "/bin/true",
            "sampleOffset": "0",
            "LogQBinsPerDecade": 10,
            "normalization": "Total charge",
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": "5.18325",
            "instrumentConfigurationDir": os.path.join(datarepo_dir.eqsans, "instrument_configuration"),
        },
    }
    reduction_input = reduction_parameters(specs, instrument_name="EQSANS")

    # replace the /bin/true filenames with None
    reduction_input["configuration"]["maskFileName"] = None
    reduction_input["configuration"]["darkFileName"] = None
    reduction_input["configuration"]["sensitivityFileName"] = None

    loaded = load_all_files(reduction_input)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].data.getHistory()

    assert history.size() == 10
    assert history.getAlgorithm(0).name() == "LoadEventNexus"
    assert history.getAlgorithm(0).getProperty("Filename").value.endswith("sns/eqsans/EQSANS_105428.nxs.h5")
    assert history.getAlgorithm(2).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(3).name() == "SetInstrumentParameter"
    assert history.getAlgorithm(5).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(6).name() == "ConvertUnits"
    assert history.getAlgorithm(7).name() == "Rebin"
    assert history.getAlgorithm(8).name() == "SetUncertainties"
    assert history.getAlgorithm(9).name() == "AddSampleLogMultiple"

    assert loaded.background.data is None
    assert loaded.background_transmission.data is None
    assert loaded.empty.data is None
    assert loaded.sample_transmission.data is None
    assert loaded.dark_current.data is None
    assert loaded.sensitivity is None
    assert loaded.mask is None

    # Verify that if something is changed that it gets applied correctly on reload, use default mask as test
    # First check the current value
    assert not loaded.sample[0].data.detectorInfo().isMasked(1)

    # Change reduction input and rerun load_all_files
    reduction_input["configuration"]["useDefaultMask"] = True
    reduction_input["configuration"]["defaultMask"] = "{'Pixel':'1'}"
    loaded = load_all_files(reduction_input)

    # Check that the value has changed
    assert loaded.sample[0].data.detectorInfo().isMasked(1)


@pytest.mark.parametrize("generic_workspace", [{"name": "ws_raw_histo"}], indirect=True)
def test_prepare_data_workspaces_simple(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    output = prepare_data_workspaces(ws_mon_pair(data=ws, monitor=None))
    clean_workspace(output)
    # this should make a clone of the workspace
    assert ws is not output
    # and change the workspace name automatically
    assert ws.name() == "ws_raw_histo"
    assert output.name() == "ws_processed_histo"

    output2 = prepare_data_workspaces(ws_mon_pair(data=ws, monitor=None), output_workspace=clean_workspace("foobar"))
    # the ws name should change to what is set
    assert ws.name() == "ws_raw_histo"
    assert output2.name() == "foobar"


def test_prepare_data_workspaces_dark_current(temp_workspace_name, clean_workspace):
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

    # Chech that this fail for the correct reason, I didn't add the all the required logs
    with pytest.raises(AttributeError) as excinfo:
        prepare_data_workspaces(
            ws_mon_pair(data=mtd[data_workspace], monitor=None),
            dark_current=ws_mon_pair(data=mtd[dark_current_workspace], monitor=None),
            solid_angle=False,
        )
    assert str(excinfo.value) == '"tof_frame_width_clipped" not found in sample logs'
    clean_workspace(data_workspace + "_processed_histo")


@pytest.mark.parametrize("generic_workspace", [{"intensities": [[1, 2], [3, 4]]}], indirect=True)
def test_prepare_data_workspaces_flux_method(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)
    SampleLogs(ws).insert("duration", 2.0)
    SampleLogs(ws).insert("monitor", 2e9)

    # No normalization
    output = prepare_data_workspaces(ws_mon_pair(data=ws, monitor=None), flux_method=None, solid_angle=False)
    clean_workspace(output)
    assert output.getHistory().size() == 3
    assert_almost_equal(output.extractY(), [[1], [2], [3], [4]])

    # Normalize by time
    output = prepare_data_workspaces(ws_mon_pair(data=ws, monitor=None), flux_method="time", solid_angle=False)
    assert output.getHistory().size() == 4
    assert_almost_equal(output.extractY(), [[0.5], [1], [1.5], [2]])

    # need to add test for proton charge and monitor


def test_prepare_data_workspaces_apply_mask(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)

    # mask_ws
    output = prepare_data_workspaces(ws_mon_pair(data=ws, monitor=None), mask_ws=[0, 2], solid_angle=False)
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
    output = prepare_data_workspaces(ws_mon_pair(data=ws, monitor=None), solid_angle=True)
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

    output = prepare_data_workspaces(
        ws_mon_pair(data=mtd[data_workspace], monitor=None),
        sensitivity_workspace=sensitivity_workspace,
        solid_angle=False,
    )
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

    output = pre_process_single_configuration(
        ws_mon_pair(data=ws, monitor=None),
        bkg_ws_raw=ws_mon_pair(data=None, monitor=None),
        solid_angle=False,
    )
    clean_workspace(output)

    # CreateWorkspace, LoadInstrument, CloneWorkspace,
    # CreateSingleValuedWorkspace, Divide,
    # CreateSingleValuedWorkspace, Multiply
    assert output.getHistory().size() == 7

    assert_equal(output.extractY(), [[1], [2], [3], [4]])

    output = pre_process_single_configuration(
        ws_mon_pair(data=ws, monitor=None),
        bkg_ws_raw=ws_mon_pair(data=None, monitor=None),
        solid_angle=False,
        absolute_scale=1.5,
    )
    assert_equal(output.extractY(), [[1.5], [3], [4.5], [6]])

    output = pre_process_single_configuration(
        ws_mon_pair(data=ws, monitor=None),
        bkg_ws_raw=ws_mon_pair(data=None, monitor=None),
        solid_angle=False,
        thickness=0.1,
    )
    assert_equal(output.extractY(), [[10], [20], [30], [40]])

    output = pre_process_single_configuration(
        ws_mon_pair(data=ws, monitor=None),
        bkg_ws_raw=ws_mon_pair(data=None, monitor=None),
        solid_angle=False,
        absolute_scale=1.5,
        thickness=0.1,
    )
    assert_equal(output.extractY(), [[15], [30], [45], [60]])


if __name__ == "__main__":
    pytest.main([__file__])
