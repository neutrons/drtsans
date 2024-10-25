import glob
import json
import os
from pathlib import Path

import numpy as np
import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import DeleteWorkspace, mtd

from drtsans.redparms import ReductionParameterError
from drtsans.tof.eqsans import reduction_parameters
from drtsans.tof.eqsans.api import (
    load_all_files,
    reduce_single_configuration,
)
from drtsans.tof.eqsans.correction_api import parse_correction_config


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "elastic_reference_run, fitInelasticIncoh, "
    + "skip_elastic, skip_inelastic, expected_do_elastic, expected_do_inelastic",
    [
        ("92160", True, False, False, True, True),  # both
        (None, False, False, False, False, False),  # none
        (None, True, False, False, False, True),  # inelastic only
        ("92160", False, False, False, True, False),  # elastic only
        ("92160", True, True, False, False, True),  # skip elastic
        ("92160", True, False, True, True, False),  # skip inelastic
        ("92160", True, True, True, False, False),  # skip both
    ],
    ids=["both", "none", "elastic_only", "inelastic_only", "skip_elastic", "skip_inelastic", "skip_both"],
)
def test_parse_json(
    elastic_reference_run,
    fitInelasticIncoh,
    skip_elastic,
    skip_inelastic,
    expected_do_elastic,
    expected_do_inelastic,
    datarepo_dir,
):
    """Test the JSON to dictionary"""

    # Specify JSON input
    reduction_input = {
        "instrumentName": "EQSANS",
        "iptsNumber": "26015",
        "sample": {"runNumber": "105428", "thickness": "1.0"},
        "background": {
            "runNumber": "104088",
            "transmission": {"runNumber": "101595", "value": ""},
        },
        "emptyTransmission": {"runNumber": "101595", "value": ""},
        "beamCenter": {"runNumber": "92353"},
        "outputFileName": "test_wavelength_step",
        "configuration": {
            "outputDir": "/path/to/nowhere",
            "cutTOFmax": "1500",
            "wavelengthStepType": "constant Delta lambda/lambda",
            "sampleApertureSize": "10",
            "fluxMonitorRatioFile": None,
            "sensitivityFileName": "/bin/true",
            "numQBins": "100",
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "5",
            "useSliceIDxAsSuffix": True,
            "fitInelasticIncoh": fitInelasticIncoh,
            "elasticReference": {
                "runNumber": elastic_reference_run,
                "thickness": "1.0",
                "transmission": {"runNumber": None, "value": "0.89"},
            },
            "elasticReferenceBkgd": {
                "runNumber": "",
                "transmission": {"runNumber": "", "value": "0.9"},
            },
            "selectMinIncoh": True,
            "maskFileName": "/bin/true",
            "darkFileName": "/bin/true",
        },
    }

    # Validate
    with amend_config(data_dir=datarepo_dir.eqsans):
        input_config = reduction_parameters(reduction_input)

    # Check that inelastic incoherence config items were parsed
    assert input_config["configuration"].get("fitInelasticIncoh") == fitInelasticIncoh
    assert input_config["configuration"]["elasticReference"].get("runNumber") == elastic_reference_run
    assert input_config["configuration"].get("selectMinIncoh")

    # Parse
    correction = parse_correction_config(input_config, skip_elastic, skip_inelastic)
    assert correction.do_elastic_correction == expected_do_elastic
    assert correction.do_inelastic_correction == expected_do_inelastic
    if expected_do_elastic:
        assert correction.elastic_reference
        assert correction.elastic_reference.run_number == "92160"
        assert correction.elastic_reference.thickness == 1.0
        assert correction.elastic_reference.transmission_value == 0.89
        assert correction.elastic_reference.background_run_number is None


@pytest.mark.datarepo
def test_parse_invalid_json(datarepo_dir):
    """Test the JSON to dictionary"""

    invalid_run_num = "260159121"
    valid_run_num = "85550"
    # Specify JSON input
    reduction_input = {
        "instrumentName": "EQSANS",
        "iptsNumber": "26015",
        "sample": {"runNumber": "85550", "thickness": "1.0"},
        "background": {
            "runNumber": "86217",
            "transmission": {"runNumber": "88565", "value": ""},
        },
        "emptyTransmission": {"runNumber": "88901", "value": ""},
        "beamCenter": {"runNumber": "88973"},
        "outputFileName": "test_wavelength_step",
        "configuration": {
            "outputDir": "/path/to/nowhere",
            "cutTOFmax": "1500",
            "wavelengthStepType": "constant Delta lambda/lambda",
            "sampleApertureSize": "10",
            "fluxMonitorRatioFile": None,
            "sensitivityFileName": "/bin/true",
            "numQBins": "100",
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "5",
            "useSliceIDxAsSuffix": True,
            "fitInelasticIncoh": True,
            "elasticReference": {
                "runNumber": invalid_run_num,
                "thickness": "1.0",
                "transmission": {"runNumber": valid_run_num, "value": "0.9"},
            },
            "elasticReferenceBkgd": {
                "runNumber": valid_run_num,
                "transmission": {"runNumber": valid_run_num, "value": "0.9"},
            },
            "selectMinIncoh": True,
            "maskFileName": "/bin/true",
            "darkFileName": "/bin/true",
        },
    }

    # Validate
    with pytest.raises(ReductionParameterError) as excinfo:
        # expect to fail as elastic reference run 260159121 does not exist
        with amend_config(data_dir=datarepo_dir.eqsans):
            reduction_parameters(reduction_input)

    assert "Cannot find events file associated to 260159121" in str(excinfo.value)

    # Respecify to use a valid run
    # json_str.replace('260159121', '26015')
    reduction_input["configuration"]["elasticReference"]["runNumber"] = valid_run_num
    # Defaults and Validate
    with amend_config(data_dir=datarepo_dir.eqsans):
        input_config = reduction_parameters(reduction_input)

    # Check that inelastic incoherence config items were parsed
    assert input_config["configuration"].get("fitInelasticIncoh")
    assert input_config["configuration"]["elasticReference"].get("runNumber") == valid_run_num
    assert input_config["configuration"].get("selectMinIncoh")


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "fitInelasticIncoh, elastic_reference_run",
    [
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ],
)
def test_incoherence_correction_elastic_normalization(
    fitInelasticIncoh, elastic_reference_run, datarepo_dir, temp_directory
):
    """Test incoherence correction with elastic correction"""

    # Set up the configuration dict
    config_json_file = os.path.join(datarepo_dir.eqsans, "test_corrections/agbe_125707_test1.json")
    assert os.path.exists(config_json_file), f"Test JSON file {config_json_file} does not exist."
    with open(config_json_file, "r") as config_json:
        configuration = json.load(config_json)
    assert isinstance(configuration, dict)

    # Create temp output directory
    test_dir = temp_directory()
    # test_dir = "test_output"
    base_name = "EQSANS_125707"
    corrections = (fitInelasticIncoh, elastic_reference_run)
    if corrections == (False, False):
        correction_case = "no_correction"
    elif corrections == (True, False):
        correction_case = "inelastic_correction"
    elif corrections == (False, True):
        correction_case = "elastic_correction"
    elif corrections == (True, True):
        correction_case = "elastic_inelastic_correction"

    output_dir = Path(test_dir) / f"{base_name}_{correction_case}"
    output_dir.mkdir(parents=True, exist_ok=True)

    configuration["configuration"]["outputDir"] = test_dir
    configuration["outputFileName"] = f"{base_name}_{correction_case}"
    configuration["dataDirectories"] = os.path.join(datarepo_dir.eqsans, "test_corrections")
    configuration["configuration"]["outputWavelengthDependentProfile"] = True
    configuration["configuration"]["maskFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_corrections", "beamstop_mask_4m_ext.nxs"
    )
    configuration["configuration"][
        "darkFileName"
    ] = "/bin/true"  # so that it will pass the validator, later set to None
    configuration["configuration"]["sensitivityFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_corrections", "Sensitivity_patched_thinPMMA_4m_124972.nxs"
    )
    configuration["configuration"]["instrumentConfigurationDir"] = os.path.join(
        datarepo_dir.eqsans, "instrument_configuration"
    )
    configuration["configuration"]["beamFluxFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_normalization", "beam_profile_flux.txt"
    )
    # lower number of bins to make I(Q) comparison less flaky
    configuration["configuration"]["numQBins"] = 80
    configuration["configuration"]["numQxQyBins"] = 40

    # override individual correction settings
    if not elastic_reference_run:
        configuration["configuration"]["elasticReference"]["runNumber"] = None
    if not fitInelasticIncoh:
        configuration["configuration"]["fitInelasticIncoh"] = False

    # validate and clean configuration
    input_config = reduction_parameters(configuration)
    input_config["configuration"]["darkFileName"] = None
    loaded = load_all_files(input_config)

    # check loaded JSON file
    if elastic_reference_run:
        assert loaded.elastic_reference.data
    assert loaded.elastic_reference_background.data is None

    # Reduce
    reduction_output = reduce_single_configuration(
        loaded_ws=loaded,
        reduction_input=input_config,
        not_apply_incoherence_correction=not fitInelasticIncoh,
        not_apply_elastic_correction=not elastic_reference_run,
    )
    assert reduction_output

    # Check output result
    iq1d_base_name = f"{base_name}_{correction_case}_Iq.dat"
    test_iq1d_file = os.path.join(test_dir, iq1d_base_name)
    assert os.path.exists(test_iq1d_file), f"Expected test result {test_iq1d_file} does not exist"

    reference_data_dir = os.path.join(datarepo_dir.eqsans, "test_corrections", correction_case)
    np.testing.assert_allclose(
        np.loadtxt(test_iq1d_file),
        np.loadtxt(os.path.join(reference_data_dir, iq1d_base_name)),
    )

    # Check 2D output result
    iq2d_base_name = f"{base_name}_{correction_case}_Iqxqy.dat"
    test_iq2d_file = os.path.join(test_dir, iq2d_base_name)
    assert os.path.exists(test_iq2d_file), f"Expected test result {test_iq2d_file} does not exist"

    np.testing.assert_allclose(
        np.loadtxt(test_iq2d_file, skiprows=4),
        np.loadtxt(os.path.join(reference_data_dir, iq2d_base_name), skiprows=4),
    )

    # check that the wavelength dependent profiles are created
    number_of_wavelengths = 31
    output_dir = os.path.join(test_dir, f"{base_name}_{correction_case}", "slice_0", "frame_0")
    if correction_case in ["elastic_correction", "elastic_inelastic_correction"]:
        assert len(glob.glob(os.path.join(output_dir, "IQ_*_before_k_correction.dat"))) == number_of_wavelengths
        assert len(glob.glob(os.path.join(output_dir, "IQ_*_after_k_correction.dat"))) == number_of_wavelengths
    elif correction_case in ["inelastic_correction", "elastic_inelastic_correction"]:
        assert len(glob.glob(os.path.join(output_dir, "IQ_*_before_b_correction.dat"))) == number_of_wavelengths
        assert len(glob.glob(os.path.join(output_dir, "IQ_*_after_b_correction.dat"))) == number_of_wavelengths
    else:
        pass

    # check the k factor file, if elastic correction is enabled
    if elastic_reference_run:
        k_base_name = f"{base_name}_{correction_case}_elastic_k1d_{base_name}.dat"
        test_k_file = os.path.join(output_dir, k_base_name)
        assert os.path.exists(test_k_file), f"Expected test result {test_k_file} does not exist"
        np.testing.assert_allclose(
            np.loadtxt(test_k_file, delimiter=",", skiprows=1),
            np.loadtxt(os.path.join(reference_data_dir, k_base_name), delimiter=",", skiprows=1),
        )

    # check the b factor file, if inelastic correction is enabled
    if fitInelasticIncoh:
        b_base_name = f"{base_name}_{correction_case}_inelastic_b1d_{base_name}.dat"
        test_b_file = os.path.join(output_dir, b_base_name)
        assert os.path.exists(test_b_file), f"Expected test result {test_b_file} does not exist"
        np.testing.assert_allclose(
            np.loadtxt(test_b_file, delimiter=",", skiprows=1),
            np.loadtxt(os.path.join(reference_data_dir, b_base_name), delimiter=",", skiprows=1),
        )

    # cleanup
    # NOTE: loaded is not a dict that is iterable, so we have to delete the
    #       leftover workspace explicitly
    # _empty:	37.277037 MB
    # _EQSANS_124667_raw_histo:	106.769917 MB
    # _EQSANS_124680_raw_histo:	40.050957 MB
    # _EQSANS_125701_raw_events:	41.005101 MB
    # _EQSANS_125701_raw_histo:	37.276829 MB
    # _EQSANS_125707_raw_histo:	38.326685 MB
    # _mask:	27.092797 MB
    # _sensitivity:	30.03614 MB
    # processed_data_main:	38.327517 MB
    # processed_elastic_ref:	40.051789 MB
    DeleteWorkspace("_empty")
    DeleteWorkspace("_mask")
    DeleteWorkspace("_sensitivity")
    DeleteWorkspace("processed_data_main")
    if elastic_reference_run:
        DeleteWorkspace("processed_elastic_ref")
    for ws in mtd.getObjectNames():
        if str(ws).startswith("_EQSANS_"):
            DeleteWorkspace(ws)


@pytest.mark.datarepo
def test_incoherence_correction_elastic_normalization_weighted(datarepo_dir, temp_directory):
    """Test incoherence correction with elastic correction"""

    # Set up the configuration dict
    config_json_file = os.path.join(datarepo_dir.eqsans, "test_corrections/porsil_29024_abs_inel.json")
    assert os.path.exists(config_json_file), f"Test JSON file {config_json_file} does not exist."
    with open(config_json_file, "r") as config_json:
        configuration = json.load(config_json)
    assert isinstance(configuration, dict)

    # Create temp output directory
    test_dir = temp_directory()

    def run_reduction_and_compare(config, expected_result_basename):
        with amend_config(data_dir=datarepo_dir.eqsans):
            # validate and clean configuration
            input_config = reduction_parameters(config)
            input_config["configuration"]["darkFileName"] = None
            loaded = load_all_files(input_config)

            # Reduce
            reduction_output = reduce_single_configuration(
                loaded, input_config, not_apply_incoherence_correction=False
            )
        assert reduction_output

        test_iq1d_file = os.path.join(test_dir, config["outputFileName"] + "_Iq.dat")
        gold_iq1d_file = os.path.join(datarepo_dir.eqsans, "test_corrections", expected_result_basename + "_Iq.dat")
        # compare
        np.testing.assert_allclose(np.loadtxt(test_iq1d_file), np.loadtxt(gold_iq1d_file))

        test_iq2d_file = os.path.join(test_dir, config["outputFileName"] + "_Iqxqy.dat")
        gold_iq2d_file = os.path.join(datarepo_dir.eqsans, "test_corrections", expected_result_basename + "_Iqxqy.dat")
        # compare
        np.testing.assert_allclose(np.loadtxt(test_iq2d_file, skiprows=4), np.loadtxt(gold_iq2d_file, skiprows=4))

        DeleteWorkspace("_empty")
        DeleteWorkspace("_mask")
        DeleteWorkspace("_sensitivity")
        DeleteWorkspace("processed_data_main")
        for ws in mtd.getObjectNames():
            if str(ws).startswith("_EQSANS_"):
                DeleteWorkspace(ws)

    # Run without intensity weighted correction
    base_name = "EQSANS_132078"
    assert os.path.exists(test_dir), f"Output dir {test_dir} does not exit"
    configuration["configuration"]["outputDir"] = test_dir
    configuration["outputFileName"] = base_name
    configuration["dataDirectories"] = os.path.join(datarepo_dir.eqsans, "test_corrections")
    configuration["configuration"][
        "darkFileName"
    ] = "/bin/true"  # so that it will pass the validator, later set to None
    configuration["configuration"]["sensitivityFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_corrections", "Sensitivity_patched_thinPMMA_4m_129610.nxs"
    )
    configuration["configuration"]["instrumentConfigurationDir"] = os.path.join(
        datarepo_dir.eqsans, "instrument_configuration"
    )
    configuration["configuration"]["maskFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_corrections", "EQSANS_132078_mask.nxs"
    )
    configuration["configuration"]["beamFluxFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_normalization", "beam_profile_flux.txt"
    )
    run_reduction_and_compare(configuration, "EQSANS_132078")

    # Run with weighted
    base_name = "EQSANS_132078_weighted"
    configuration["outputFileName"] = base_name
    configuration["configuration"]["incohfit_intensityweighted"] = True
    configuration["configuration"]["incohfit_factor"] = None
    configuration["configuration"]["incohfit_qmin"] = None
    configuration["configuration"]["incohfit_qmax"] = None
    run_reduction_and_compare(configuration, "EQSANS_132078_weighted")

    # Run with weighted and factor
    base_name = "EQSANS_132078_weighted_factor"
    configuration["outputFileName"] = base_name
    configuration["configuration"]["incohfit_intensityweighted"] = True
    configuration["configuration"]["incohfit_factor"] = 10
    configuration["configuration"]["incohfit_qmin"] = None
    configuration["configuration"]["incohfit_qmax"] = None
    run_reduction_and_compare(configuration, "EQSANS_132078_weighted_factor")

    # Run with weighted and manual q range
    # q-range is set to be the same as what the factor calculation finds
    base_name = "EQSANS_132078_weighted_qrange"
    configuration["outputFileName"] = base_name
    configuration["configuration"]["incohfit_intensityweighted"] = True
    configuration["configuration"]["incohfit_factor"] = None
    configuration["configuration"]["incohfit_qmin"] = 0.065
    configuration["configuration"]["incohfit_qmax"] = 0.224
    run_reduction_and_compare(configuration, "EQSANS_132078_weighted_factor")

    print(f"Output directory: {test_dir}")


@pytest.mark.mount_eqsans
def test_incoherence_correction_elastic_normalization_slices_frames(has_sns_mount, datarepo_dir, temp_directory):
    """Test incoherence correction with elastic correction with time slicing and frame mode"""
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    # Set up the configuration dict
    config_json_file = os.path.join(datarepo_dir.eqsans, "test_corrections/test_136908.json")
    assert os.path.exists(config_json_file), f"Test JSON file {config_json_file} does not exist."
    with open(config_json_file, "r") as config_json:
        configuration = json.load(config_json)
    assert isinstance(configuration, dict)

    # Create temp output directory
    test_dir = temp_directory()
    base_name = "EQSANS_136908"

    assert os.path.exists(test_dir), f"Output dir {test_dir} does not exit"
    configuration["configuration"]["outputDir"] = test_dir
    configuration["outputFileName"] = base_name

    # validate and clean configuration
    input_config = reduction_parameters(configuration)
    input_config["configuration"]["darkFileName"] = None
    loaded = load_all_files(input_config)

    # check loaded JSON file
    assert loaded.elastic_reference.data
    assert loaded.elastic_reference_background.data is None

    # Reduce
    reduction_output = reduce_single_configuration(loaded, input_config, not_apply_incoherence_correction=False)
    assert reduction_output
    print(f"Output directory: {test_dir}")

    # check that the wavelength dependent profiles are created in subdirectories for slices and frames
    for islice in range(3):
        for iframe in range(2):
            if iframe == 0:
                number_of_wavelengths = 29
            else:
                number_of_wavelengths = 28
            output_dir = os.path.join(test_dir, base_name, f"slice_{islice}", f"frame_{iframe}")
            # before k correction
            assert len(glob.glob(os.path.join(output_dir, "IQ_*_before_k_correction.dat"))) == number_of_wavelengths
            # after k correction
            assert len(glob.glob(os.path.join(output_dir, "IQ_*_after_k_correction.dat"))) == number_of_wavelengths
            # before b correction
            assert len(glob.glob(os.path.join(output_dir, "IQ_*_before_b_correction.dat"))) == number_of_wavelengths
            # after b correction
            assert len(glob.glob(os.path.join(output_dir, "IQ_*_after_b_correction.dat"))) == number_of_wavelengths

    # cleanup
    DeleteWorkspace("_empty")
    DeleteWorkspace("_mask")
    DeleteWorkspace("_sensitivity")
    DeleteWorkspace("_bkgd_trans")
    DeleteWorkspace("_sample_trans")
    DeleteWorkspace("TOFCorrectWS")
    DeleteWorkspace("processed_data_main")
    DeleteWorkspace("processed_elastic_ref")
    for ws in mtd.getObjectNames():
        if str(ws).startswith("_EQSANS_") and mtd.doesExist(str(ws)):
            DeleteWorkspace(ws)


if __name__ == "__main__":
    pytest.main([__file__])
