import os
from pathlib import Path

import pytest

from drtsans.redparams import reduction_parameters
from drtsans.tof.eqsans import load_all_files, reduce_single_configuration, plot_reduction_output


REDUCTION_INPUT = {
    "instrumentName": "EQSANS",
    "iptsNumber": "27799",
    "outputFileName": "agbe_125707_test1",
    "sample": {"runNumber": "125707", "thickness": 1, "transmission": {"runNumber": "", "value": "1"}},
    "schemaStamp": "2021-08-06T17:32:34.528330",
    "dataDirectories": "",
    "emptyTransmission": {"runNumber": "125701", "value": ""},
    "background": {"runNumber": "", "transmission": {"runNumber": "", "value": ""}},
    "beamCenter": {"runNumber": "125701"},
    "configuration": {
        "1DQbinType": "annular",
        "AnnularAngleBin": 1.0,
        "QbinType": "linear",
        "Qmax": "",
        "Qmin": "",
        "StandardAbsoluteScale": 1,
        "absoluteScaleMethod": "standard",
        "beamFluxFileName": "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample",
        "cutTOFmax": 2000.0,
        "cutTOFmin": 500.0,
        "darkFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2021B_mp/EQSANS_124667.nxs.h5",
        "defaultMask": "",
        "detectorOffset": "80",
        "fitInelasticIncoh": True,
        "elasticReference": {
            "runNumber": "124680",
            "thickness": "1.0",
            "transmission": {"runNumber": "", "value": "1.0"},
        },
        "elasticReferenceBkgd": {"runNumber": "", "transmission": {"runNumber": "", "value": "0.9"}},
        "fluxMonitorRatioFile": "",
        "maskFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2021B_mp/beamstop_mask_4m_ext.nxs",
        "mmRadiusForTransmission": 25,
        "normalization": "Total charge",
        "numQBins": 100,
        "numQxQyBins": 80,
        "sampleApertureSize": 10,
        "sampleOffset": "314.5",
        "useDefaultMask": True,
        "wavelengthStep": 0.1,
        "wavelengthStepType": "constant Delta lambda",
    },
}


@pytest.mark.datarepo
def test_reduce_single_configuration(datarepo_dir, temp_directory):
    """Test reduction workflow with annular binning

    Note: applies elastic reference normalization and incoherent inelastic correction
    """
    reduction_input = REDUCTION_INPUT.copy()
    output_dir = temp_directory(prefix="test_annular_binning")
    reduction_input["configuration"]["outputDir"] = output_dir

    datadir = Path(datarepo_dir.eqsans) / "test_corrections"
    reduction_input["dataDirectories"] = str(datadir)
    reduction_input["configuration"]["instrumentConfigurationDir"] = str(
        Path(datarepo_dir.eqsans) / "instrument_configuration"
    )
    reduction_input["configuration"]["maskFileName"] = str(Path(datadir) / "beamstop_mask_4m_ext.nxs")
    reduction_input["configuration"]["beamFluxFileName"] = str(
        Path(datarepo_dir.eqsans) / "instrument_configuration" / "bl6_flux_at_sample.dat"
    )
    # Set darkFileName so that it will pass the validator, later set to None
    reduction_input["configuration"]["darkFileName"] = "/bin/true"
    reduction_input["configuration"]["sensitivityFileName"] = str(
        Path(datadir) / "Sensitivity_patched_thinPMMA_4m_124972.nxs"
    )

    reduction_input = reduction_parameters(parameters_particular=reduction_input, validate=True)

    reduction_input["configuration"]["darkFileName"] = None
    loaded = load_all_files(
        reduction_input,
        prefix="EQ_TEST_LOAD",
        load_params=None,
    )
    assert len(loaded.sample) == 1
    output = reduce_single_configuration(loaded, reduction_input)
    assert len(output) == 1
    assert os.path.exists(os.path.join(output_dir, "agbe_125707_test1_Iphi.dat"))
    plot_reduction_output(output, reduction_input)
    assert os.path.exists(os.path.join(output_dir, "agbe_125707_test1_Iphi.png"))
    del output
