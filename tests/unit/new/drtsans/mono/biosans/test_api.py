import pytest
import os
import numpy as np
from numpy.testing import assert_almost_equal
from mantid.simpleapi import CreateWorkspace, LoadHFIRSANS
from drtsans.mono.biosans.api import load_all_files, prepare_data_workspaces
from drtsans.samplelogs import SampleLogs
from drtsans.settings import unique_workspace_dundername as uwd


@pytest.mark.skipif(not os.path.exists('/HFIR/CG3/IPTS-23782/nexus/CG3_960.nxs.h5'),
                    reason="Required data is not available")
def test_load_all_files_simple():
    reduction_input = {
        "instrumentName": "CG3",
        "iptsNumber": "23782",
        "runNumber": "960",
        "transmission": {"runNumber": ""},
        "background": {"runNumber": "",
                       "transmission": {"runNumber": ""}},
        "emptyTrans": {"runNumber": ""},
        "beamCenter": {"runNumber": "960"},
        "configuration": {
            "useMaskFileName": False,
            "useDefaultMask": False,
            "useBlockedBeam": False,
            "useDarkFileName": False,
            "useSensitivityFileName": False,
            "wavelength": "",
            "wavelengthSpread": ""}
    }
    loaded = load_all_files(reduction_input)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].getHistory()
    assert history.size() == 4
    assert history.getAlgorithm(0).name() == "LoadEventNexus"
    assert history.getAlgorithm(0).getProperty("Filename").value == '/HFIR/CG3/IPTS-23782/nexus/CG3_960.nxs.h5'
    assert history.getAlgorithm(1).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(2).name() == "HFIRSANS2Wavelength"
    assert history.getAlgorithm(3).name() == "SetUncertainties"

    assert loaded.background is None
    assert loaded.background_transmission is None
    assert loaded.center is not None
    assert loaded.empty is None
    assert loaded.sample_transmission is None
    assert loaded.blocked_beam is None
    assert loaded.dark_current is None
    assert loaded.sensitivity is None
    assert loaded.mask is None


@pytest.mark.parametrize('generic_workspace', [{'name': 'ws_raw_histo'}], indirect=True)
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
    ws = LoadHFIRSANS(Filename=biosans_f['beamcenter'])

    output = prepare_data_workspaces(ws,
                                     center_x=0.111,
                                     center_y=0.123,
                                     center_y_wing=0.222,
                                     solid_angle=False)

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


def test_prepare_data_workspaces_dark_current():
    # Create dark current workspace, insert the duration of the dark
    # current run as one of the log entries in the dark current
    # workspace.
    dark_current_workspace = uwd()  # arbitrary name for the dark current workspace
    CreateWorkspace(DataX=[2.5, 3.5],
                    DataY=np.full(2, 100.),
                    DataE=np.full(2, 10.),
                    NSpec=2,
                    OutputWorkspace=dark_current_workspace)
    SampleLogs(dark_current_workspace).insert('duration', 3600., 'second')

    # Create a sample run workspace.
    data_workspace = uwd()  # arbitrary name for the sample workspace
    CreateWorkspace(DataX=[2.5, 3.5],
                    DataY=np.array([1., 2.]),
                    DataE=np.array([1., np.sqrt(2)]),
                    NSpec=2,
                    OutputWorkspace=data_workspace)
    # Insert the duration of the sample run. The log key must be the
    # same as that used for the dark current, which turns out to be
    # 'duration'
    SampleLogs(data_workspace).insert('duration', 36., 'second')

    output = prepare_data_workspaces(data_workspace,
                                     dark_current=dark_current_workspace,
                                     solid_angle=False)

    assert output.getHistory().size() == 8

    assert_almost_equal(output.extractY(), [[0], [1]])
    assert_almost_equal(output.extractE(), [[1.00498756], [1.41774469]])


@pytest.mark.parametrize('generic_workspace', [{'intensities': [[1, 2], [3, 4]]}], indirect=True)
def test_prepare_data_workspaces_flux_method(generic_workspace):
    ws = generic_workspace  # friendly name
    SampleLogs(ws).insert("duration", 2.)
    SampleLogs(ws).insert("monitor", 2e9)

    # No normalization
    output = prepare_data_workspaces(ws,
                                     flux_method=None,
                                     solid_angle=False)
    assert output.getHistory().size() == 3
    assert_almost_equal(output.extractY(), [[1], [2], [3], [4]])

    # Normalize by time
    output = prepare_data_workspaces(ws,
                                     flux_method="time",
                                     solid_angle=False)
    assert output.getHistory().size() == 5
    assert_almost_equal(output.extractY(), [[0.5], [1], [1.5], [2]])

    # Normalize by monitor, should scale by 1e8/(monitor counts)
    output = prepare_data_workspaces(ws,
                                     flux_method="monitor",
                                     solid_angle=False)
    assert output.getHistory().size() == 5
    assert_almost_equal(output.extractY(), [[0.05], [0.1], [0.15], [0.2]])


def test_prepare_data_workspaces_apply_mask(generic_workspace):
    ws = generic_workspace

    # mask_ws
    output = prepare_data_workspaces(ws,
                                     mask_ws=[0, 2],
                                     solid_angle=False)
    history = output.getHistory()
    assert history.size() == 4
    alg3 = history.getAlgorithm(3)
    assert alg3.name() == "MaskDetectors"
    assert alg3.getPropertyValue("DetectorList") == "0,2"


@pytest.mark.parametrize('generic_workspace', [{'intensities': [[1, 1], [1, 1]]}], indirect=True)
def test_prepare_data_workspaces_solid_angle(generic_workspace):
    ws = generic_workspace  # friendly name

    # No normalization
    output = prepare_data_workspaces(ws,
                                     solid_angle=True)
    # CreateWorkspace, LoadInstrument, CloneWorkspace, CloneWorkspace,
    # ClearMaskFlag, SolidAngle, Divide, ReplaceSpecialValues
    assert output.getHistory().size() == 8
    assert_almost_equal(output.extractY(), [[25.6259267], [25.6259267], [25.6259267], [25.6259267]])


def test_prepare_data_workspaces_sensitivity():
    # Create dark current workspace, insert the duration of the dark
    # current run as one of the log entries in the dark current
    # workspace.
    sensitivity_workspace = uwd()  # arbitrary name for the dark current workspace
    CreateWorkspace(DataX=[2.5, 3.5],
                    DataY=np.full(2, 2.),
                    DataE=np.full(2, np.sqrt(2)),
                    NSpec=2,
                    OutputWorkspace=sensitivity_workspace)

    # Create a sample run workspace.
    data_workspace = uwd()  # arbitrary name for the sample workspace
    CreateWorkspace(DataX=[2.5, 3.5],
                    DataY=np.array([1., 2.]),
                    DataE=np.array([1., np.sqrt(2)]),
                    NSpec=2,
                    OutputWorkspace=data_workspace)

    output = prepare_data_workspaces(data_workspace,
                                     sensitivity_workspace=sensitivity_workspace,
                                     solid_angle=False)

    assert output.getHistory().size() == 6

    assert_almost_equal(output.extractY(), [[0.5], [1.0]])
    assert_almost_equal(output.extractE(), [[0.6123724], [1.0]])


if __name__ == '__main__':
    pytest.main([__file__])
