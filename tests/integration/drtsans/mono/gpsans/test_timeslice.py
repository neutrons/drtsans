import pytest
import os
from mantid.simpleapi import mtd
from drtsans.mono.gpsans import (
    load_all_files,
    reduction_parameters,
    reduce_single_configuration,
)
from drtsans.samplelogs import SampleLogs
from drtsans.api import NoDataProcessedError


@pytest.mark.datarepo
def test_timeslice(datarepo_dir, temp_directory):
    reduction_input = {
        "instrumentName": "GPSANS",
        "iptsNumber": "24727",
        "sample": {"runNumber": "9166", "thickness": "0.1", "transmission": {"runNumber": "9178"}},
        "outputFileName": "Al4",
        "background": {"runNumber": "9165", "transmission": {"runNumber": "9177"}},
        "beamCenter": {"runNumber": "9177"},
        "emptyTransmission": {"runNumber": "9177"},
        "configuration": {
            "outputDir": "",
            "sampleApertureSize": "",
            "sourceApertureDiameter": "",
            "maskFileName": "",
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-10,247-256'}"],
            "blockedBeamRunNumber": "",
            "darkFileName": "",
            "sensitivityFileName": "sens_c486_noBar.nxs",
            "absoluteScaleMethod": "direct_beam",
            "DBScalingBeamRadius": "40",
            "StandardAbsoluteScale": "1",
            "normalization": "Monitor",
            "sampleOffset": "",
            "useSolidAngleCorrection": True,
            "useThetaDepTransCorrection": True,
            "mmRadiusForTransmission": "40",
            "numQxQyBins": "180",
            "1DQbinType": "scalar",
            "QbinType": "log",
            "LogQBinsPerDecade": "33",
            "useLogQBinsEvenDecade": True,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "numQBins": "",
            "AnnularAngleBin": "",
            "Qmin": "",
            "Qmax": "",
            "useErrorWeighting": False,
            "useMaskBackTubes": False,
            "wavelength": "1.23",
            "wavelengthSpread": "",
            "useTimeSlice": True,
            "timeSliceInterval": 50,
            "logSliceName": "",
            "useLogSlice": False,
            "logSliceInterval": "",
            "smearingPixelSizeX": "1.2345",
            "smearingPixelSizeY": "2.3456",
            "sampleToSi": "234.56",
            "sampleDetectorDistance": "32.11",
        },
    }

    # set output directory
    reduction_input["configuration"]["outputDir"] = temp_directory(prefix="trans_slice")
    reduction_input["dataDirectories"] = datarepo_dir
    reduction_input["configuration"]["sensitivityFileName"] = os.path.join(
        datarepo_dir.gpsans, "overwrite_gold_04282020", "sens_c486_noBar.nxs"
    )

    # check inputs
    reduction_input = reduction_parameters(reduction_input, validate=True)
    # load files
    loaded = load_all_files(
        reduction_input,
        prefix="GP_TEST_LOAD",
        load_params=None,
        path=datarepo_dir.gpsans,
    )
    assert len(loaded.sample) == 7  # 300.048/50 rounded up

    output = reduce_single_configuration(loaded, reduction_input)
    assert len(output) == 7  # nothing skipped, same length as loaded

    # set every second workspace monitor counts to zero
    for ws in loaded.sample[::2]:
        SampleLogs(ws).insert("monitor", 0.0)

    output = reduce_single_configuration(loaded, reduction_input)
    assert len(output) == 3  # 4 out of 7 skipped

    # set every workspace monitor counts to zero
    for ws in loaded.sample:
        SampleLogs(ws).insert("monitor", 0.0)

    with pytest.raises(NoDataProcessedError) as e:
        reduce_single_configuration(loaded, reduction_input)
    assert "No data was processed. Check the input data." in str(e.value)
    mtd.clear()
