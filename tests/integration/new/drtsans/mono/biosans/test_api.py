from math import isclose
import pytest
import os
import threading
from mantid.simpleapi import mtd, DeleteWorkspace
from mantid.api import AnalysisDataService
from drtsans.mono.transmission import calculate_transmission
from drtsans.mono.biosans.api import load_all_files, reduce_single_configuration

# need to hardcode these worspace removals, is there a way to derive them?
# is there a safe way to just ping for a list of workspaces generated by this user?
workspaces = [
    "_bkgd_trans",
    "_empty",
    "_filter",
    "_info",
    "_load_tmp",
    "_load_tmp_monitors",
    "_main_sensitivity",
    "_sample_trans",
    "_wing_sensitivity",
    "barscan_BIOSANS_detector1_20210302",
    "barscan_BIOSANS_wing_detector_20210302",
    "processed_data_main_0",
    "processed_data_main_1",
    "processed_data_main_10",
    "processed_data_main_2",
    "processed_data_main_3",
    "processed_data_main_4",
    "processed_data_main_5",
    "processed_data_main_6",
    "processed_data_main_7",
    "processed_data_main_8",
    "processed_data_main_9",
    "processed_data_wing_0",
    "processed_data_wing_1",
    "processed_data_wing_10",
    "processed_data_wing_2",
    "processed_data_wing_3",
    "processed_data_wing_4",
    "processed_data_wing_5",
    "processed_data_wing_6",
    "processed_data_wing_7",
    "processed_data_wing_8",
    "processed_data_wing_9",
    "TOFCorrectWS",
    "tubewidth_BIOSANS_detector1_20200613",
    "tubewidth_BIOSANS_detector1_20210304",
    "tubewidth_BIOSANS_wing_detector_20200613",
    "tubewidth_BIOSANS_wing_detector_20210304",
]


@pytest.mark.skipif(
    not os.path.exists("/HFIR/HB2B/shared/autoreduce/"),
    reason="Skip test on build server",
)
def test_reduce_single_configuration_slice_transmission_false(generatecleanfile):
    reduction_input = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "BIOSANS",
        "iptsNumber": "24666",
        "dataDirectories": None,
        "sample": {
            "runNumber": "8375",
            "thickness": 0.2,
            "transmission": {"runNumber": "8379", "value": None},
        },
        "background": {
            "runNumber": "8374",
            "transmission": {"runNumber": "8378", "value": None},
        },
        "emptyTransmission": {"runNumber": "8381", "value": None},
        "beamCenter": {"runNumber": "8381"},
        "outputFileName": "r8375_AgBeh_15m18Aqa",
        "configuration": {
            "wavelength": None,
            "wavelengthSpread": None,
            "useTimeSlice": False,
            "useTimeSliceTransmission": False,
            "timeSliceInterval": 60.0,
            "useLogSlice": False,
            "logSliceName": None,
            "logSliceInterval": None,
            "sampleOffset": None,
            "sampleApertureSize": 12.0,
            "sampleDetectorDistance": None,
            "sampleToSi": None,
            "sourceApertureDiameter": None,
            "usePixelCalibration": True,
            "maskFileName": None,
            "useDefaultMask": True,
            "defaultMask": [
                {"Pixel": "1-18,239-256"},
                {"Bank": "18-24,42-48"},
                {"Bank": "49", "Tube": "1"},
                {"Bank": "88", "Tube": "4"},
            ],
            "useMaskBackTubes": False,
            "darkMainFileName": "CG3_8331.nxs.h5",
            "darkWingFileName": "CG3_8331.nxs.h5",
            "normalization": "Monitor",
            "normalizationResortToTime": False,
            "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8367m7p0_bsSVP.nxs",
            "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8369w1p4_bsSVP.nxs",
            "useSolidAngleCorrection": True,
            "blockedBeamRunNumber": None,
            "useThetaDepTransCorrection": True,
            "DBScalingBeamRadius": 40.0,
            "mmRadiusForTransmission": None,
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": 2.286e-09,
            "numMainQxQyBins": 100,
            "numWingQxQyBins": 100,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numMainQBins": None,
            "numWingQBins": None,
            "LogQBinsPerDecadeMain": 25,
            "LogQBinsPerDecadeWing": 25,
            "useLogQBinsDecadeCenter": False,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": None,
            "WedgeMaxAngles": None,
            "autoWedgeQmin": 0.003,
            "autoWedgeQmax": 0.04,
            "autoWedgeQdelta": 0.01,
            "autoWedgeAzimuthalDelta": 1.0,
            "autoWedgePeakWidth": 0.5,
            "autoWedgeBackgroundWidth": 1.0,
            "autoWedgeSignalToNoiseMin": 2.0,
            "AnnularAngleBin": 1.0,
            "useErrorWeighting": False,
            "smearingPixelSizeX": None,
            "smearingPixelSizeY": None,
            "useSubpixels": False,
            "subpixelsX": None,
            "subpixelsY": None,
            "QminMain": 0.0009,
            "QmaxMain": 0.016,
            "QminWing": 0.009,
            "QmaxWing": 0.3,
            "overlapStitchQmin": [0.0105],
            "overlapStitchQmax": [0.0145],
            "wedge1QminMain": 0.003,
            "wedge1QmaxMain": 0.0425,
            "wedge1QminWing": 0.02,
            "wedge1QmaxWing": 0.45,
            "wedge1overlapStitchQmin": 0.025,
            "wedge1overlapStitchQmax": 0.04,
            "wedge2QminMain": 0.003,
            "wedge2QmaxMain": 0.0425,
            "wedge2QminWing": 0.03,
            "wedge2QmaxWing": 0.45,
            "wedge2overlapStitchQmin": 0.03,
            "wedge2overlapStitchQmax": 0.04,
            "wedges": None,
            "symmetric_wedges": True,
        },
        "logslice_data": {},
    }
    reduction_input["configuration"]["outputDir"] = generatecleanfile(
        prefix="trans_slice_false"
    )

    prefix = "sans-backend-test" + str(threading.get_ident()) + "_"
    loaded = load_all_files(reduction_input, prefix)
    _ = reduce_single_configuration(loaded, reduction_input)
    # just need a couple components from reduce
    # but the whole thing needs to be run then a few components pulled
    transmission = calculate_transmission(
        mtd["_sample_trans"],  # pull relevant transmission
        mtd["_empty"],  # pull relevant
    )
    transmission_val = transmission.extractY()[0][0]
    clean_all_ws(prefix)
    assert isclose(transmission_val, 0.5734218305525239)  # provided by s6v
    del _


def clean_all_ws(prefix):
    for workspace in workspaces:
        remove_ws(workspace)
    for object_name in mtd.getObjectNames():
        if object_name.startswith(prefix):
            remove_ws(object_name)


def remove_ws(workspace):
    if AnalysisDataService.doesExist(workspace):
        DeleteWorkspace(workspace)


def test_reduce_single_configuration_slice_transmission_true(generatecleanfile):
    reduction_input = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "BIOSANS",
        "iptsNumber": "24666",
        "dataDirectories": None,
        "sample": {
            "runNumber": "8361",
            "thickness": 0.1,
            "transmission": {"runNumber": "8361", "value": None},
        },
        "background": {
            "runNumber": "8359",
            "transmission": {"runNumber": "8359", "value": None},
        },
        "emptyTransmission": {"runNumber": "8364", "value": None},
        "beamCenter": {"runNumber": "8373"},
        "outputFileName": "r8361_PorB3_15m",
        "configuration": {
            "wavelength": None,
            "wavelengthSpread": None,
            "useTimeSlice": True,
            "useTimeSliceTransmission": True,
            "timeSliceInterval": 60.0,
            "useLogSlice": False,
            "logSliceName": None,
            "logSliceInterval": None,
            "sampleOffset": None,
            "sampleApertureSize": 14.0,
            "sampleDetectorDistance": None,
            "sampleToSi": None,
            "sourceApertureDiameter": None,
            "usePixelCalibration": True,
            "maskFileName": None,
            "useDefaultMask": True,
            "defaultMask": [
                {"Pixel": "1-18,239-256"},
                {"Bank": "18-24,42-48"},
                {"Bank": "49", "Tube": "1"},
                {"Bank": "88", "Tube": "4"},
            ],
            "useMaskBackTubes": False,
            "darkMainFileName": "CG3_8331.nxs.h5",
            "darkWingFileName": "CG3_8331.nxs.h5",
            "normalization": "Monitor",
            "normalizationResortToTime": False,
            "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8367m7p0_bsSVP.nxs",
            "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle490/Sens_f8369w1p4_bsSVP.nxs",
            "useSolidAngleCorrection": True,
            "blockedBeamRunNumber": None,
            "useThetaDepTransCorrection": True,
            "DBScalingBeamRadius": 40.0,
            "mmRadiusForTransmission": None,
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": 2.094e-10,
            "numMainQxQyBins": 100,
            "numWingQxQyBins": 100,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numMainQBins": None,
            "numWingQBins": None,
            "LogQBinsPerDecadeMain": 25,
            "LogQBinsPerDecadeWing": 25,
            "useLogQBinsDecadeCenter": False,
            "useLogQBinsEvenDecade": False,
            "WedgeMinAngles": None,
            "WedgeMaxAngles": None,
            "autoWedgeQmin": 0.003,
            "autoWedgeQmax": 0.04,
            "autoWedgeQdelta": 0.01,
            "autoWedgeAzimuthalDelta": 1.0,
            "autoWedgePeakWidth": 0.5,
            "autoWedgeBackgroundWidth": 1.0,
            "autoWedgeSignalToNoiseMin": 2.0,
            "AnnularAngleBin": 1.0,
            "useErrorWeighting": False,
            "smearingPixelSizeX": None,
            "smearingPixelSizeY": None,
            "useSubpixels": False,
            "subpixelsX": None,
            "subpixelsY": None,
            "QminMain": 0.003,
            "QmaxMain": 0.045,
            "QminWing": 0.03,
            "QmaxWing": 0.9,
            "overlapStitchQmin": [0.0325],
            "overlapStitchQmax": [0.0425],
            "wedge1QminMain": 0.003,
            "wedge1QmaxMain": 0.0425,
            "wedge1QminWing": 0.02,
            "wedge1QmaxWing": 0.45,
            "wedge1overlapStitchQmin": 0.025,
            "wedge1overlapStitchQmax": 0.04,
            "wedge2QminMain": 0.003,
            "wedge2QmaxMain": 0.0425,
            "wedge2QminWing": 0.03,
            "wedge2QmaxWing": 0.45,
            "wedge2overlapStitchQmin": 0.03,
            "wedge2overlapStitchQmax": 0.04,
            "wedges": None,
            "symmetric_wedges": True,
        },
        "logslice_data": {},
    }
    reduction_input["configuration"]["outputDir"] = generatecleanfile(
        prefix="trans_slice_true"
    )
    prefix = "sans-backend-test" + str(threading.get_ident()) + "_"
    loaded = load_all_files(reduction_input, prefix)
    _ = reduce_single_configuration(loaded, reduction_input)
    # just need a couple components from reduce
    # but the whole thing needs to be run then a few components pulled
    transmission = calculate_transmission(
        mtd["_sample_trans"],  # pull relevant transmission
        mtd["_empty"],  # pull relevant
    )

    transmission_val = transmission.extractY()[0][0]
    clean_all_ws(prefix)
    assert isclose(
        transmission_val, 0.7526460467895154  # from above config using older workflow
    )
    del _


if __name__ == "__main__":
    pytest.main(__file__)