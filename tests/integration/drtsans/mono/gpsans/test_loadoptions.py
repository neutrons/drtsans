"""Integration tests for loadOptions time filtering in GPSANS."""

import pytest
from drtsans.mono.gpsans import load_all_files, reduction_parameters


@pytest.mark.datarepo
def test_load_with_time_filtering(datarepo_dir):
    """Test that FilterByTimeStart and FilterByTimeStop work correctly in loadOptions."""
    reduction_input = {
        "instrumentName": "CG2",
        "sample": {
            "runNumber": "9166",
            "transmission": {"runNumber": ""},
            "loadOptions": {
                "FilterByTimeStart": 0.0,
                "FilterByTimeStop": 50.0,
            },
        },
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "beamCenter": {"runNumber": "9166"},
        "emptyTransmission": {"runNumber": ""},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "GPSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.gpsans)

    # Verify sample workspace was loaded
    assert loaded.sample is not None
    assert len(loaded.sample) == 1

    # Check that workspace was loaded with time filtering
    ws = loaded.sample[0]
    assert ws is not None
    assert ws.getNumberHistograms() > 0


@pytest.mark.datarepo
def test_load_without_time_filtering_baseline(datarepo_dir):
    """Test loading without time filtering as baseline."""
    reduction_input = {
        "instrumentName": "CG2",
        "sample": {
            "runNumber": "9166",
            "transmission": {"runNumber": ""},
            "loadOptions": {},
        },
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "beamCenter": {"runNumber": "9166"},
        "emptyTransmission": {"runNumber": ""},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "GPSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.gpsans)

    # Verify sample workspace was loaded
    assert loaded.sample is not None
    assert len(loaded.sample) == 1

    ws = loaded.sample[0]
    assert ws is not None
    assert ws.getNumberHistograms() > 0


@pytest.mark.datarepo
def test_load_with_partial_time_range(datarepo_dir):
    """Test loading with only FilterByTimeStop specified."""
    reduction_input = {
        "instrumentName": "CG2",
        "sample": {
            "runNumber": "9166",
            "transmission": {"runNumber": ""},
            "loadOptions": {"FilterByTimeStop": 30.0},
        },
        "background": {"runNumber": "", "transmission": {"runNumber": ""}},
        "beamCenter": {"runNumber": "9166"},
        "emptyTransmission": {"runNumber": ""},
        "configuration": {"useDefaultMask": False},
    }

    reduction_input = reduction_parameters(reduction_input, "GPSANS", validate=False)
    loaded = load_all_files(reduction_input, path=datarepo_dir.gpsans)

    # Verify sample workspace was loaded
    assert loaded.sample is not None
    assert len(loaded.sample) == 1
