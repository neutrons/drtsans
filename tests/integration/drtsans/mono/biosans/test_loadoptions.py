"""Integration tests for loadOptions time filtering in BIOSANS."""

import pytest
from mantid.simpleapi import SumSpectra
from drtsans.mono.biosans import load_all_files, reduction_parameters


@pytest.mark.datarepo
def test_load_with_time_filtering(datarepo_dir):
    """Test that FilterByTimeStart and FilterByTimeStop work correctly in loadOptions."""
    reduction_input = {
        "instrumentName": "CG3",
        "sample": {
            "runNumber": "1322",
            "transmission": {"runNumber": ""},
            "loadOptions": {
                "FilterByTimeStart": 0.0,
                "FilterByTimeStop": 10.0,
            },
        },
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "beamCenter": {"runNumber": "1322"},
        "emptyTransmission": {"runNumber": ""},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "BIOSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.biosans)

    # Verify sample workspace was loaded
    assert loaded.sample is not None
    assert len(loaded.sample) == 1

    # Check that workspace was loaded with time filtering
    ws = loaded.sample[0]
    assert ws is not None
    assert ws.getNumberHistograms() > 0

    # Verify time filtering was applied by checking event count
    ws_summed = SumSpectra(ws)
    assert ws_summed.dataY(0)[0] == 2283  # number of events in the first 10 seconds


@pytest.mark.datarepo
def test_load_without_time_filtering_baseline(datarepo_dir):
    """Test loading without time filtering as baseline."""
    reduction_input = {
        "instrumentName": "CG3",
        "sample": {
            "runNumber": "1322",
            "transmission": {"runNumber": ""},
            "loadOptions": {},
        },
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "beamCenter": {"runNumber": "1322"},
        "emptyTransmission": {"runNumber": ""},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "BIOSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.biosans)

    # Verify sample workspace was loaded
    assert loaded.sample is not None
    assert len(loaded.sample) == 1

    ws = loaded.sample[0]
    assert ws is not None
    assert ws.getNumberHistograms() > 0


@pytest.mark.datarepo
def test_load_with_partial_time_range(datarepo_dir):
    """Test loading with only FilterByTimeStart specified."""
    reduction_input = {
        "instrumentName": "CG3",
        "sample": {
            "runNumber": "1322",
            "transmission": {"runNumber": ""},
            "loadOptions": {"FilterByTimeStart": 5.0},
        },
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "beamCenter": {"runNumber": "1322"},
        "emptyTransmission": {"runNumber": ""},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "BIOSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.biosans)

    # Verify sample workspace was loaded successfully with partial time filter
    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    ws = loaded.sample[0]
    assert ws is not None
    assert ws.getNumberHistograms() > 0
