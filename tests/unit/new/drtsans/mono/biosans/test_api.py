import pytest
import os
from drtsans.mono.biosans.api import load_all_files, prepare_data_workspaces


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


if __name__ == '__main__':
    pytest.main([__file__])
