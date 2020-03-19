import pytest
import os
from drtsans.tof.eqsans.api import load_all_files


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-22747/nexus/EQSANS_105428.nxs.h5'),
                    reason="Required data is not available")
def test_load_all_files_simple():
    reduction_input = {
        "instrumentName": "EQSANS",
        "iptsNumber": "22747",
        "runNumber": "105428",
        "thickness": "1",
        "transmission": {"runNumber": ""},
        "background": {"runNumber": "",
                       "transmission": {"runNumber": ""}},
        "empty": {"runNumber": ""},
        "beamCenter": {"runNumber": "105428"},
        "configuration": {
            "useMaskFileName": False,
            "useDefaultMask": False,
            "useDarkFileName": False,
            "useSensitivityFileName": False,
            "wavelengthStep": "",
            "useDetectorOffset": False,
            "useSampleOffset": True,
            "sampleOffset": "0",
            "useTOFcuts": False,
            "normalization": "Total charge",
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": "5.18325",
        }
    }
    loaded = load_all_files(reduction_input)

    assert loaded.sample is not None
    assert len(loaded.sample) == 1
    history = loaded.sample[0].data.getHistory()
    assert history.size() == 7
    assert history.getAlgorithm(0).name() == "LoadEventNexus"
    assert history.getAlgorithm(0).getProperty("Filename").value == '/SNS/EQSANS/IPTS-22747/nexus/EQSANS_105428.nxs.h5'
    assert history.getAlgorithm(1).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(2).name() == "MoveInstrumentComponent"
    assert history.getAlgorithm(3).name() == "ConvertUnits"
    assert history.getAlgorithm(4).name() == "Rebin"
    assert history.getAlgorithm(5).name() == "SetUncertainties"
    assert history.getAlgorithm(6).name() == "AddSampleLogMultiple"

    assert loaded.background.data is None
    assert loaded.background_transmission.data is None
    assert loaded.empty.data is None
    assert loaded.sample_transmission.data is None
    assert loaded.dark_current.data is None
    assert loaded.sensitivity is None
    assert loaded.mask is None


if __name__ == '__main__':
    pytest.main([__file__])
