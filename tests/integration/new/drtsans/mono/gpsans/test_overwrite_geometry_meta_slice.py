import json
import os
from drtsans.mono.gpsans import load_all_files, reduce_single_configuration


def setup_json_config():

    reduction_input = json.loads(""" {
        "instrumentName": "CG2",
        "iptsNumber": "22474",
        "runNumber": "7872",
        "thickness": "0.05",
        "outputFilename": "CG2_7872time",
        "transmission": {
            "runNumber": "7871",
            "value": ""
        },
        "background": {
            "runNumber": "7853",
            "transmission": {
                "runNumber": "7854",
                "value": ""
            }
        },
        "beamCenter": {
            "runNumber": "7854"
        },
        "emptyTrans": {
            "runNumber": "7854"
        },
        "configuration": {
            "outputDir": "/HFIR/CG2/IPTS-22474/shared/KCL_time",
            "sampleApertureSize": "",
            "maskFileName": "/HFIR/CG2/IPTS-22474/shared/trap12.nxs",
            "useMaskFileName": true,
            "useDefaultMask": true,
            "DefaultMask":["{'Pixel':'1-10,247-256'}"],
            "useBlockedBeam": true,
            "BlockBeamRunNumber":"7845",
            "useDarkFileName": false,
            "darkFilename": "",
            "useSensitivityFileName": true,
            "sensitivityFileName": "/HFIR/CG2/shared/drt_sensitivity/sens_c486_noBar.nxs",
            "UseBarScan": true,
            "BarScanFileName": "",
            "absoluteScaleMethod": "direct_beam",
            "DBScalingBeamRadius": "40",
            "StandardAbsoluteScale": "1",
            "normalization": "Time",
            "sampleOffset": "",
            "useSampleOffset": false,
            "useSolidAngleCorrection": true,
            "useThetaDepTransCorrection": true,
            "mmRadiusForTransmission": "",
            "numQxQyBins": "150",
            "1DQbinType": "scalar",
            "QbinType": "log",
            "EvenDecades": true,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "numQBins": "33",
            "AnnularAngleBin": "5",
            "Qmin": "",
            "Qmax": "",
            "useErrorWeighting": true,
            "useMaskBackTubes": false,
            "wavelength": "",
            "wavelengthSpread": "",
             "timeslice": true,
            "timesliceinterval": "60",
            "logslicename": "",
            "logslice": false,
            "logsliceinterval": ""
        }
    }""")

    return reduction_input


def reduce_data(reduction_input):
    # chekcing if output directory exists, if it doesn't, creates the folder
    output_dir = reduction_input["configuration"]["outputDir"]
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    loaded = load_all_files(reduction_input)

    out = reduce_single_configuration(loaded, reduction_input, prefix='')

    assert out
