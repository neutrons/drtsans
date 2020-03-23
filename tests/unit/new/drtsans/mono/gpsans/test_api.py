import pytest
import os
from numpy.testing import assert_almost_equal
from drtsans.mono.gpsans.api import load_all_files, prepare_data_workspaces
from drtsans.samplelogs import SampleLogs


@pytest.mark.skipif(not os.path.exists('/HFIR/CG2/IPTS-23801/nexus/CG2_1338.nxs.h5'),
                    reason="Required data is not available")
def test_load_all_files_simple():
    reduction_input = {
        "instrumentName": "CG2",
        "iptsNumber": "23801",
        "runNumber": "1338",
        "transmission": {"runNumber": ""},
        "background": {"runNumber": "",
                       "transmission": {"runNumber": ""}},
        "emptyTrans": {"runNumber": ""},
        "beamCenter": {"runNumber": "1338"},
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
    assert history.getAlgorithm(0).getProperty("Filename").value == '/HFIR/CG2/IPTS-23801/nexus/CG2_1338.nxs.h5'
    assert history.getAlgorithm(1).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(2).name() == "HFIRSANS2Wavelength"
    assert history.getAlgorithm(3).name() == "SetUncertainties"

    assert loaded.background is None
    assert loaded.background_transmisson is None
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

    output = prepare_data_workspaces(ws, solid_angle=False)
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
    output2 = prepare_data_workspaces(ws, output_workspace="foobar", solid_angle=False)
    assert output2.name() == "foobar"


def test_prepare_data_workspaces_center(generic_workspace):
    ws = generic_workspace  # friendly name

    output = prepare_data_workspaces(ws,
                                     center_x=0.111,
                                     center_y=0.123,
                                     solid_angle=False)

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


if __name__ == '__main__':
    pytest.main([__file__])
