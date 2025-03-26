# standard library imports
import os

# third party imports
import pytest

# local imports
from drtsans.mono.biosans import load_all_files, reduce_single_configuration, plot_reduction_output
from drtsans.redparams import reduction_parameters


REDUCTION_INPUT = {
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
        "sampleOffset": None,
        "sampleApertureSize": 14.0,
        "useDefaultMask": True,
        "defaultMask": [
            {"Pixel": "1-18,239-256"},
            {"Bank": "18-24,42-48"},
            {"Bank": "49", "Tube": "1"},
            {"Bank": "88", "Tube": "4"},
        ],
        "useMaskBackTubes": False,
        "normalization": "Monitor",
        "normalizationResortToTime": False,
        "sensitivityMainFileName": "Sens_f6368m4p0_bsSVP.nxs",
        "sensitivityWingFileName": "Sens_f6380w1p4_bsSVP.nxs",
        "useSolidAngleCorrection": True,
        "blockedBeamRunNumber": None,
        "useThetaDepTransCorrection": True,
        "DBScalingBeamRadius": 40.0,
        "mmRadiusForTransmission": None,
        "absoluteScaleMethod": "standard",
        "StandardAbsoluteScale": 2.094e-10,
        "numMainQxQyBins": 100,
        "numWingQxQyBins": 100,
        "numMidrangeQxQyBins": 100,
        "1DQbinType": "annular",
        "QbinType": "linear",
        "AnnularAngleBin": 1.0,
        "numMainQBins": 100,
        "numWingQBins": 100,
        "numMidrangeQBins": 100,
        "useErrorWeighting": False,
        "QminMain": 0.003,
        "QmaxMain": 0.045,
        "QminWing": 0.03,
        "QmaxWing": 0.9,
        "QminMidrange": 0.03,
        "QmaxMidrange": 0.9,
        "overlapStitchQmin": [0.0325],
        "overlapStitchQmax": [0.0425],
    },
}


@pytest.mark.datarepo
def test_reduce_single_configuration(datarepo_dir, temp_directory):
    reduction_input = REDUCTION_INPUT.copy()
    output_dir = temp_directory(prefix="test_annular_binning")
    reduction_input["configuration"]["outputDir"] = output_dir
    reduction_input["dataDirectories"] = datarepo_dir
    reduction_input["configuration"]["sensitivityMainFileName"] = os.path.join(
        datarepo_dir.biosans, "Sens_f6368m4p0_bsSVP.nxs"
    )
    reduction_input["configuration"]["sensitivityWingFileName"] = os.path.join(
        datarepo_dir.biosans, "Sens_f6380w1p4_bsSVP.nxs"
    )
    reduction_input = reduction_parameters(parameters_particular=reduction_input, validate=True)

    load_params = {}
    loaded = load_all_files(
        reduction_input,
        prefix="BIOSANS_TEST_LOAD",
        load_params=load_params,
        path=datarepo_dir.biosans,
    )
    assert len(loaded.sample) == 1
    output = reduce_single_configuration(loaded, reduction_input)
    assert len(output) == 1
    plot_reduction_output(output, reduction_input, loglog=False)
    assert os.path.exists(os.path.join(output_dir, "1D", "r8361_PorB3_15m_1D.png"))
    del output
