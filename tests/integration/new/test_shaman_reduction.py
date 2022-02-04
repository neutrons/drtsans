import json
import os
import pytest
import subprocess
from tempfile import NamedTemporaryFile
import time

from drtsans import scriptsdir
from drtsans.tof.eqsans import reduction_parameters, update_reduction_parameters

# this should point to the root directory of the code repository
ROOT_DIR = os.path.abspath(os.path.join(__file__, "../../../../"))
# specific output files for given basenames
FILES = {
    "EQSANS_88980": [
        "EQSANS_88980_bkgd_88974_trans.txt",
        "EQSANS_88980_bkgd_88974_raw_trans.txt",
        "EQSANS_88980_trans.txt",
        "EQSANS_88980_raw_trans.txt",
        "EQSANS_88980.nxs",
        "EQSANS_88980_frame_0_Iqxqy.dat",
        "EQSANS_88980_frame_0_Iq.dat",
        "EQSANS_88980_frame_1_Iqxqy.dat",
        "EQSANS_88980_frame_1_Iq.dat",
        "EQSANS_88980_reduction_log.hdf",
        "EQSANS_88980_0_Iqxqy.png",
        "EQSANS_88980_0_Iq.png",
        "EQSANS_88980_1_Iqxqy.png",
        "EQSANS_88980_1_Iq.png",
    ],
    "EQSANS_112300": [
        "EQSANS_112300_Iq.dat",
        "EQSANS_112300_Iq.png",
        "EQSANS_112300_Iqxqy.dat",
        "EQSANS_112300_Iqxqy.png",
        "EQSANS_112300_bkgd_112296_raw_trans.txt",
        "EQSANS_112300_bkgd_112296_trans.txt",
        "EQSANS_112300.nxs",
        "EQSANS_112300_raw_trans.txt",
        "EQSANS_112300_reduction_log.hdf",
        "EQSANS_112300_trans.txt",
    ],
    "CG3_4822": [],
    "CG3_4822_wedge": [
        "1D/CG3_4822_wedge_1D_main_wedge_0.txt",
        "1D/CG3_4822_wedge_1D_main_wedge_1.txt",
        "1D/CG3_4822_wedge_1D_wing_wedge_0.txt",
        "1D/CG3_4822_wedge_1D_wing_wedge_1.txt",
        "1D/CG3_4822_wedge_1D_both_wedge_0.txt",
        "1D/CG3_4822_wedge_1D_both_wedge_1.txt",
        "1D/CG3_4822_wedge_1D_wedge_0.png",
        "1D/CG3_4822_wedge_1D_wedge_1.png",
        "2D/CG3_4822_wedge_2D_main.dat",
        "2D/CG3_4822_wedge_2D_wing.dat",
        "2D/CG3_4822_wedge_2D_main.png",
        "2D/CG3_4822_wedge_2D_wing.png",
    ],
    "CG3_5532": [],  # auto-wedge
    "CG2_8944": [],
}


def write_configfile(input_json_file, basename, tmpdir, config_overrides=dict()):
    """
        Create a new json configuration file with a better place for the output files
        and a standardized basename
    def test_eqsans
        Parameters
        ----------
        input_json_file: str
            json file to modify. Expected to be in ``ROOT_DIR/scripts``
        basename: str
            The value of ``outputFileName``. This makes it significantly easier to find the output files.
        config_overrides: dict
            if it is specified, then additional configuration will be set specific to instruments

        Returns
        -------
        str, str
            The name of the output directory and the name of the (re)configured json file
    """
    # get the full path to a configuration file
    input_json_file = os.path.join(ROOT_DIR, "scripts", input_json_file)
    assert os.path.exists(input_json_file), 'Could not find "{}"'.format(
        input_json_file
    )

    # temporary directory is always readable
    outputdir = tmpdir
    output_json_file = NamedTemporaryFile(
        prefix=os.path.basename(input_json_file).replace(".json", "_"),
        suffix=".json",
        delete=False,
    ).name
    print("Reconfigured json file set to {}".format(output_json_file))

    # read the existing file in
    with open(input_json_file, "r") as handle:
        json_params = json.load(handle)

    # get the cannonical runnumber name to standardize it
    json_params["outputFileName"] = basename
    print('Output basename set to "{}"'.format(json_params["outputFileName"]))

    # update the output directory
    json_params["configuration"]["outputDir"] = outputdir
    print('Output directory set to "{}"'.format(outputdir))

    # set additional values or overwrite the values from script
    for conf_name, value in config_overrides.items():
        json_params["configuration"][conf_name] = value

    # write out the new file
    with open(output_json_file, "w") as handle:
        json.dump(json_params, handle)
    return outputdir, output_json_file


def run_reduction(python_script, json_file):
    python_script = os.path.join(scriptsdir, python_script)
    assert os.path.exists(python_script), 'Could not find "{}"'.format(python_script)

    # run the script
    cmd = "python3 {} {}".format(python_script, json_file)
    print('Running "{}"'.format(cmd))
    start = time.clock()
    proc = subprocess.Popen(cmd, shell=True, universal_newlines=True)
    proc.communicate()

    # 0 = ok
    # 42 = non-empty error log file. This is the case for developers running the test on
    # the console without a special configuration of their logging
    returncode = proc.returncode
    assert returncode in [0, 42]
    print(python_script, "took", time.clock() - start, "seconds")


def check_and_cleanup(outputdir, basename):
    # verify that the output files were created and cleanup
    err_msg = ""
    num_found = 0
    for f in FILES[basename]:
        filename = os.path.join(outputdir, f)
        if not os.path.isfile(filename):
            err_msg += '"{}" does not exist\n'.format(filename)
        else:
            num_found += 1
            os.remove(filename)
    if err_msg != "":
        raise AssertionError(
            f"Only {num_found} expected files are found out of {len(FILES[basename])}:\n{err_msg}"
        )

    # remove the output and error logs if they were created, neither is required to exist
    for ext in [".out", ".err", "_reduction_log.hdf"]:
        logname = os.path.join(outputdir, basename + ext)
        if os.path.isfile(logname):
            os.remove(logname)

    # verify that there isn't anything else with the same prefix
    extra_files = [
        os.path.join(outputdir, name)
        for name in os.listdir(outputdir)
        if basename in name
    ]
    nonempty_files = []
    for filename in extra_files:
        if os.path.isfile(filename) and os.stat(filename).st_size == 0:
            os.remove(filename)
        else:
            nonempty_files.append(filename)
    assert not nonempty_files, "Found extra files {}".format(",".join(nonempty_files))


def check_for_required_files(filenames):
    """Skip the test if a required file is missing"""
    for filename in filenames:
        if not os.path.exists(filename):
            pytest.skip('"{}" is not available'.format(filename))


specs_eqsans = {
    "EQSANS_88980": {
        "iptsNumber": 19800,
        "sample": {
            "runNumber": 88980,
            "thickness": 0.1,
            "transmission": {"runNumber": 88980},
        },
        "background": {"runNumber": 88978, "transmission": {"runNumber": 88974}},
        "beamCenter": {"runNumber": 88973},
        "emptyTransmission": {"runNumber": 88973},
        "configuration": {
            "sampleApertureSize": 30,
            "darkFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5",
            "StandardAbsoluteScale": 0.0208641883,
            "sampleOffset": 0,
        },
    },
    "EQSANS_112300": {
        "iptsNumber": 24769,
        "sample": {
            "runNumber": 112300,
            "thickness": 1.0,
            "transmission": {"runNumber": 112297},
        },
        "background": {"runNumber": 112299, "transmission": {"runNumber": 112296}},
        "beamCenter": {"runNumber": 112295},
        "emptyTransmission": {"runNumber": 112295},
        "configuration": {
            "sampleApertureSize": 10,
            "darkFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2019B_mp/EQSANS_108764.nxs.h5",
            "StandardAbsoluteScale": 1,
            "sampleOffset": 314.5,
        },
    },
}


@pytest.mark.skip(reason="Shaman is deprecated.")
@pytest.mark.parametrize(
    "run_config, basename",
    [
        (specs_eqsans["EQSANS_88980"], "EQSANS_88980"),
        (specs_eqsans["EQSANS_112300"], "EQSANS_112300"),
    ],
    ids=["88980", "112300"],
)
def test_eqsans(run_config, basename, generatecleanfile, cleanfile):
    common_config = {
        "configuration": {
            "maskFileName": "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs",
            "useDefaultMask": True,
            "normalization": "Total charge",
            "fluxMonitorRatioFile": "/SNS/EQSANS/IPTS-24769/shared/EQSANS_110943.out",
            "beamFluxFileName": "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample",
            "absoluteScaleMethod": "standard",
            "detectorOffset": 0,
            "mmRadiusForTransmission": 25,
            "numQxQyBins": 80,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "numQBins": 120,
            "AnnularAngleBin": 5,
            "wavelengthStepType": "constant Delta lambda",
            "wavelengthStep": 0.1,
        }
    }
    input_config = reduction_parameters(
        common_config, "EQSANS", validate=False
    )  # defaults and common options
    input_config = update_reduction_parameters(input_config, run_config, validate=False)
    output_dir = generatecleanfile()
    amendments = {
        "outputFileName": basename,
        "configuration": {"outputDir": output_dir},
    }
    input_config = update_reduction_parameters(
        input_config, amendments, validate=True
    )  # final changes and validation

    json_file = NamedTemporaryFile(suffix=".json", delete=False).name
    with open(json_file, "w") as handle:
        json.dump(input_config, handle)
    cleanfile(json_file)
    run_reduction("eqsans_reduction.py", json_file)
    check_and_cleanup(output_dir, basename)


@pytest.mark.skip(reason="Shaman is deprecated.")
@pytest.mark.parametrize(
    "configfile, basename, required",
    [
        (
            "biosans_reduction.json",
            "CG3_4822",
            ("/HFIR/CG3/IPTS-23782/nexus/CG3_4822.nxs.h5",),
        ),
        (
            "biosans_wedge_reduction.json",
            "CG3_4822_wedge",
            ("/HFIR/CG3/IPTS-23782/nexus/CG3_4822.nxs.h5",),
        ),
        (
            "biosans_autowedge_reduction.json",
            "CG3_5532",
            ("/HFIR/CG3/IPTS-21089/nexus/CG3_5532.nxs.h5",),
        ),
    ],
    ids=["4822", "4822_wedge", "5532_autowedge"],
)
def test_biosans(configfile, basename, required, generatecleanfile, reference_dir):
    check_for_required_files(required)

    # Set for sensitivities files
    config_overwrite_dict = dict()
    sensitivity_dir = os.path.join(reference_dir.new.biosans, "overwrite_gold_04282020")

    #  '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4829m7p0_TDC_SAC.h5'
    config_overwrite_dict["sensitivityMainFileName"] = os.path.join(
        sensitivity_dir, "sens_f4829m7p0_TDC_SAC.h5"
    )

    # '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4835w3p2_TDC_SAC.h5'
    config_overwrite_dict["sensitivityWingFileName"] = os.path.join(
        sensitivity_dir, "sens_f4835w3p2_TDC_SAC.h5"
    )

    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(
        configfile, basename, generatecleanfile(), config_overwrite_dict
    )

    run_reduction("biosans_reduction.py", json_file)

    check_and_cleanup(outputdir, basename)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


@pytest.mark.skip(reason="Shaman is deprecated.")
@pytest.mark.parametrize(
    "configfile, basename, required",
    [
        (
            "gpsans_reduction.json",
            "CG2_8944",
            ("/HFIR/CG2/IPTS-20775/nexus/CG2_8944.nxs.h5",),
        )
    ],
    ids=["8944"],
)
def test_gpsans(configfile, basename, required, generatecleanfile):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename, generatecleanfile())

    run_reduction("gpsans_reduction.py", json_file)

    check_and_cleanup(outputdir, basename)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


if __name__ == "__main__":
    pytest.main([__file__])
