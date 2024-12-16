import glob
import json
import os

import numpy as np
import pytest
from mantid.kernel import amend_config
from mantid.simpleapi import DeleteWorkspace, mtd

from drtsans.redparams import ReductionParameterError
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
        ("92160", [True], False, False, True, [True, True]),
        (None, [False], False, False, False, [False, False]),
        (None, [True], False, False, False, [True, True]),
        ("92160", [False], False, False, True, [False, False]),
        ("92160", [True], True, False, False, [True, True]),
        ("92160", [True], False, True, True, [False, False]),
        ("92160", [True], True, True, False, [False, False]),
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
            "fitInelasticIncoh": [True],
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
    outputFileName = f"{base_name}_{correction_case}"

    configuration["configuration"]["outputDir"] = test_dir
    configuration["outputFileName"] = outputFileName
    configuration["dataDirectories"] = os.path.join(datarepo_dir.eqsans, "test_corrections")
    configuration["configuration"]["outputWavelengthDependentProfile"] = True
    configuration["configuration"]["maskFileName"] = os.path.join(
        datarepo_dir.eqsans, "test_corrections", "beamstop_mask_4m_ext.nxs"
    )
    # Set darkFileName so that it will pass the validator, later set to None
    configuration["configuration"]["darkFileName"] = "/bin/true"
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
        configuration["configuration"]["fitInelasticIncoh"] = [False]

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

    ### Check output results

    elastic_output_dir = os.path.join(test_dir, "info", "elastic_norm", f"{outputFileName}", "slice_0", "frame_0")
    inelastic_output_dir = os.path.join(test_dir, "info", "inelastic_incoh", f"{outputFileName}", "slice_0", "frame_0")

    # Check empty subdirectories not created for no correction case
    if correction_case == "no_correction":
        assert not os.path.exists(elastic_output_dir)
        assert not os.path.exists(inelastic_output_dir)

    # Check 1D output result
    iq1d_base_name = f"{outputFileName}_Iq.dat"
    test_iq1d_file = os.path.join(test_dir, iq1d_base_name)
    assert os.path.exists(test_iq1d_file), f"Expected test result {test_iq1d_file} does not exist"

    reference_data_dir = os.path.join(datarepo_dir.eqsans, "test_corrections", correction_case)
    np.testing.assert_allclose(
        np.loadtxt(test_iq1d_file),
        np.loadtxt(os.path.join(reference_data_dir, iq1d_base_name)),
    )

    # Check 2D output result
    iq2d_base_name = f"{outputFileName}_Iqxqy.dat"
    test_iq2d_file = os.path.join(test_dir, iq2d_base_name)
    assert os.path.exists(test_iq2d_file), f"Expected test result {test_iq2d_file} does not exist"

    np.testing.assert_allclose(
        np.loadtxt(test_iq2d_file, skiprows=4),
        np.loadtxt(os.path.join(reference_data_dir, iq2d_base_name), skiprows=4),
    )

    # Check that the wavelength dependent profiles are created
    number_of_wavelengths = 31
    if correction_case in ["elastic_correction", "elastic_inelastic_correction"]:
        assert (
            len(glob.glob(os.path.join(elastic_output_dir, "IQ_*_before_k_correction.dat"))) == number_of_wavelengths
        )
        assert len(glob.glob(os.path.join(elastic_output_dir, "IQ_*_after_k_correction.dat"))) == number_of_wavelengths
    elif correction_case in ["inelastic_correction", "elastic_inelastic_correction"]:
        assert (
            len(glob.glob(os.path.join(inelastic_output_dir, "IQ_*_before_b_correction.dat"))) == number_of_wavelengths
        )
        assert (
            len(glob.glob(os.path.join(inelastic_output_dir, "IQ_*_after_b_correction.dat"))) == number_of_wavelengths
        )
    else:
        pass

    # check the k factor file, if elastic correction is enabled
    if elastic_reference_run:
        k_base_name = f"{outputFileName}_elastic_k1d_{base_name}.dat"
        test_k_file = os.path.join(elastic_output_dir, k_base_name)
        assert os.path.exists(test_k_file), f"Expected test result {test_k_file} does not exist"
        np.testing.assert_allclose(
            np.loadtxt(test_k_file, delimiter=",", skiprows=1),
            np.loadtxt(os.path.join(reference_data_dir, k_base_name), delimiter=",", skiprows=1),
        )

    # check the b factor file, if inelastic correction is enabled
    if fitInelasticIncoh:
        b_base_name = f"{outputFileName}_inelastic_b1d_{base_name}.dat"
        test_b_file = os.path.join(inelastic_output_dir, b_base_name)
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
@pytest.mark.parametrize(
    "base_name, expected_result_basename, weighted, qmin, qmax, factor",
    [
        ("EQSANS_132078", "EQSANS_132078", False, None, None, None),
        ("EQSANS_132078_weighted", "EQSANS_132078_weighted", True, None, None, None),
        ("EQSANS_132078_weighted_factor", "EQSANS_132078_weighted_factor", True, None, None, 10),
        ("EQSANS_132078_weighted_qrange", "EQSANS_132078_weighted_factor", True, 0.065, 0.224, None),
    ],
    ids=["unweighted", "weighted", "weighted_factor", "weighted_qrange"],
)
def test_incoherence_correction_elastic_normalization_weighted(
    base_name, expected_result_basename, weighted, qmin, qmax, factor, datarepo_dir, temp_directory
):
    """Test incoherence correction with elastic correction"""

    def _run_reduction_and_compare(config, expected_result_basename):
        """Run reduction and compare the output with the expected result"""
        with amend_config(data_dir=datarepo_dir.eqsans):
            # validate and clean configuration
            input_config = reduction_parameters(config)
            input_config["configuration"]["darkFileName"] = None
            loaded = load_all_files(input_config)

            # Reduce
            reduction_output = reduce_single_configuration(loaded, input_config)
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

    # Create temp output directory
    test_dir = temp_directory()
    assert os.path.exists(test_dir), f"Output dir {test_dir} does not exit"

    # Set up the configuration dict
    config_json_file = os.path.join(datarepo_dir.eqsans, "test_corrections/porsil_29024_abs_inel.json")
    assert os.path.exists(config_json_file), f"Test JSON file {config_json_file} does not exist."
    with open(config_json_file, "r") as config_json:
        configuration = json.load(config_json)
    assert isinstance(configuration, dict)

    # Override common configuration values
    configuration["configuration"]["outputDir"] = test_dir
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

    # Override individual correction settings
    configuration["outputFileName"] = base_name
    configuration["configuration"]["incohfit_intensityweighted"] = weighted
    configuration["configuration"]["incohfit_qmin"] = qmin
    configuration["configuration"]["incohfit_qmax"] = qmax
    configuration["configuration"]["incohfit_factor"] = factor
    _run_reduction_and_compare(configuration, expected_result_basename)

    print(f"Output directory: {test_dir}")


@pytest.mark.mount_eqsans
@pytest.mark.parametrize(
    "fitInelasticIncoh, incohfit_intensityweighted, incohfit_qmin, incohfit_qmax, incohfit_factor",
    [
        # Original configuration
        (True, False, 0.04, 0.08, None),
        # Multiple fitInelasticIncoh values
        ([True, False], False, 0.04, 0.08, None),
        # Multiple incohfit_intensityweighted values
        (True, [False, True], 0.04, 0.08, None),
        # Multiple incohfit_qmin/qmax values
        (True, False, [0.04, 0.06], [0.08, 0.1], None),
        # Multiple incohfit_factor values
        (True, False, 0.04, 0.08, [None, 10]),
    ],
)
def test_incoherence_correction_elastic_normalization_slices_frames(
    fitInelasticIncoh,
    incohfit_intensityweighted,
    incohfit_qmin,
    incohfit_qmax,
    incohfit_factor,
    has_sns_mount,
    datarepo_dir,
    temp_directory,
):
    """Test incoherence correction with elastic correction with time slicing and frame mode"""
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    # Set up the configuration dict
    config_json_file = os.path.join(datarepo_dir.eqsans, "test_corrections/test_136908.json")
    assert os.path.exists(config_json_file), f"Test JSON file {config_json_file} does not exist."
    with open(config_json_file, "r") as config_json:
        configuration = json.load(config_json)
    assert isinstance(configuration, dict)

    # Override configuration values
    configuration["configuration"]["fitInelasticIncoh"] = fitInelasticIncoh
    configuration["configuration"]["incohfit_intensityweighted"] = incohfit_intensityweighted
    configuration["configuration"]["incohfit_qmin"] = incohfit_qmin
    configuration["configuration"]["incohfit_qmax"] = incohfit_qmax
    configuration["configuration"]["incohfit_factor"] = incohfit_factor

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
    reduction_output = reduce_single_configuration(loaded, input_config)
    assert reduction_output
    print(f"Output directory: {test_dir}")

    # check that the wavelength dependent profiles are created in subdirectories for slices and frames
    for islice in range(3):
        for iframe in range(2):
            # 29 wavelengths for the first frame, 28 for the second frame
            # unless inelastic correction is disabled for the frame
            if isinstance(fitInelasticIncoh, list) and fitInelasticIncoh[iframe] is False:
                num_wavelengths_k = 29 - iframe
                num_wavelengths_b = 0
            else:
                num_wavelengths_k = num_wavelengths_b = 29 - iframe
            elastic_output_dir = os.path.join(
                test_dir, "info", "elastic_norm", f"{base_name}", f"slice_{islice}", f"frame_{iframe}"
            )
            inelastic_output_dir = os.path.join(
                test_dir, "info", "inelastic_incoh", f"{base_name}", f"slice_{islice}", f"frame_{iframe}"
            )
            # before k correction
            assert (
                len(glob.glob(os.path.join(elastic_output_dir, "IQ_*_before_k_correction.dat"))) == num_wavelengths_k
            )
            # after k correction
            assert len(glob.glob(os.path.join(elastic_output_dir, "IQ_*_after_k_correction.dat"))) == num_wavelengths_k
            # before b correction
            assert (
                len(glob.glob(os.path.join(inelastic_output_dir, "IQ_*_before_b_correction.dat"))) == num_wavelengths_b
            )
            # after b correction
            assert (
                len(glob.glob(os.path.join(inelastic_output_dir, "IQ_*_after_b_correction.dat"))) == num_wavelengths_b
            )

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
