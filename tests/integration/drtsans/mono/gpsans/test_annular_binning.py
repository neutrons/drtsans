import os

import pytest

from drtsans.mono.gpsans import (
    load_all_files,
    plot_reduction_output,
    reduce_single_configuration,
)
from drtsans.redparams import reduction_parameters


REDUCTION_INPUT = {
    "instrumentName": "GPSANS",
    "iptsNumber": "24727",
    "sample": {"runNumber": "9166", "thickness": "0.1", "transmission": {"runNumber": "9178"}},
    "outputFileName": "Al4",
    "background": {"runNumber": "9165", "transmission": {"runNumber": "9177"}},
    "beamCenter": {"runNumber": "9177"},
    "emptyTransmission": {"runNumber": "9177"},
    "configuration": {
        "outputDir": "",
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
        "useSolidAngleCorrection": True,
        "useThetaDepTransCorrection": True,
        "mmRadiusForTransmission": "40",
        "numQxQyBins": "180",
        "1DQbinType": "annular",
        "AnnularAngleBin": "1.0",
        "QbinType": "linear",
        "numQBins": "100",
        "Qmin": "",
        "Qmax": "",
        "useErrorWeighting": False,
        "useMaskBackTubes": False,
        "wavelength": "1.23",
        "wavelengthSpread": "",
        "sampleToSi": "234.56",
        "sampleDetectorDistance": "32.11",
    },
}


@pytest.mark.datarepo
def test_reduce_single_configuration(datarepo_dir, temp_directory):
    reduction_input = REDUCTION_INPUT.copy()
    output_dir = temp_directory(prefix="test_annular_binning")
    reduction_input["configuration"]["outputDir"] = output_dir
    reduction_input["dataDirectories"] = datarepo_dir
    reduction_input["configuration"]["sensitivityFileName"] = os.path.join(
        datarepo_dir.gpsans, "overwrite_gold_04282020", "sens_c486_noBar.nxs"
    )
    reduction_input = reduction_parameters(parameters_particular=reduction_input, validate=True)

    loaded = load_all_files(
        reduction_input,
        prefix="GP_TEST_LOAD",
        load_params=None,
        path=datarepo_dir.gpsans,
    )
    assert len(loaded.sample) == 1
    output = reduce_single_configuration(loaded, reduction_input)
    assert len(output) == 1
    plot_reduction_output(output, reduction_input, loglog=False, close_figures=True)
    assert os.path.exists(os.path.join(output_dir, "1D", "Al4_1D.png"))
    del output
