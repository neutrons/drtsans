import pytest
import os
import tempfile
import json
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (load_all_files, reduce_single_configuration,plot_reduction_output)  # noqa E402
from drtsans.tof.eqsans import validate_reduction_parameters


@pytest.mark.skipif(not os.path.exists('/SNS/EQSANS/IPTS-26015/nexus/EQSANS_115363.nxs.h5'),
                    reason='Required test data not available')
def test_parse_json():
    """Test the JSON to dictionary
    """
    # Specify JSON string
    json_str = """
        "instrumentName": "EQSANS",
        "iptsNumber": "26015",
        "sample": {
            "runNumber": "115363",
            "thickness": "1.0"
        },
        "background": {
            "runNumber": "115361",
            "transmission": {
                "runNumber": "115357",
                "value": ""
            }
        },
        "emptyTransmission": {
            "runNumber": "115356",
            "value": ""
        },
        "beamCenter": {
            "runNumber": "115356"
        },
        "outputFileName": "test_wavelength_step",
        "configuration": {
            "cutTOFmax": "1500",
            "wavelengthStepType": "constant Delta lambda/lambda",
            "sampleApertureSize": "10",
            "fluxMonitorRatioFile": ("/SNS/EQSANS/"
                                     "IPTS-24769/shared/EQSANS_110943.out"),
            "sensitivityFileName": ("/SNS/EQSANS/shared/NeXusFiles/"
                                    "EQSANS/2020A_mp/"
                                    "Sensitivity_patched_thinPMMA_2o5m_113514_mantid.nxs"),
            "numQBins": "100",
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "5",
            "useSliceIDxAsSuffix": true
        }
        "fitInelasticIncoh": true,
        "elasticReference": 260159121,
        "selectMinIncoh": true
    """

    # Parse
    reduction_input = json.loads(json_str)
    # Validate
    with pytest.raises(RuntimeError):
        # TODO - expect to fail as elastic reference run 260159121 does not exist
        input_config = validate_reduction_parameters(reduction_input)

    # Respecify to use a valid run
    json_str.replace('260159121', '26015')
    # Parse
    reduction_input = json.loads(json_str)
    # Validate
    with pytest.raises(RuntimeError):
        # TODO - expect to fail as elastic reference run 260159121 does not exist
        input_config = validate_reduction_parameters(reduction_input)

    # Specify expected output
    expected_configuration = {
        "instrumentName": "EQSANS",
        "iptsNumber": 26015,
        "sample": {
            "runNumber": 115363,
            "thickness": 1.0
        },
        "background": {
            "runNumber": 115361,
            "transmission": {
                "runNumber": 115357,
                "value": ""
            }
        },
        "emptyTransmission": {
            "runNumber": 115356,
            "value": ""
        },
        "beamCenter": {
            "runNumber": 115356
        },
        "outputFileName": "test_wavelength_step",
        "configuration": {
            "cutTOFmax": "1500",
            "wavelengthStepType": "constant Delta lambda/lambda",
            "sampleApertureSize": "10",
            "fluxMonitorRatioFile": ("/SNS/EQSANS/"
                                     "IPTS-24769/shared/EQSANS_110943.out"),
            "sensitivityFileName": ("/SNS/EQSANS/shared/NeXusFiles/"
                                    "EQSANS/2020A_mp/"
                                    "Sensitivity_patched_thinPMMA_2o5m_113514_mantid.nxs"),
            "numQBins": 100,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "5",
            "useSliceIDxAsSuffix": True
        },

        "fitInelasticIncoh": True,
        "elasticReference": 26015,
        "selectMinIncoh": True
    }

    # input config shall be same as expected config
    assert input_config == expected_configuration


if __name__ == "__main__":
    pytest.main([__file__])