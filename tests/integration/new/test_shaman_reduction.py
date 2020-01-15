import json
import os
import pytest
import subprocess
from tempfile import gettempdir, NamedTemporaryFile
import time

# this should point to the root directory of the code repository
ROOT_DIR = os.path.abspath(os.path.join(__file__, '../../../../'))
# specific extensions for given basenames
EXTENSIONS = {'EQSANS_88980': ['.png', '_bkg_88974_trans.txt', '_frame_1_Iq.txt', '_frame_2_Iq.txt',
                               '_frame_2_Iqxqy.txt', '_qxqy1.png', '_qxqy2.png', '_trans.txt'],
              'CG3_1433': ['_merged_Iq.json', '_merged_Iq.png']}
# double for loop to simplify creating the mix of extensions for CG3_1433
for iq_type in '_Iq', '_Iqxqy', '_wing_Iq', '_wing_Iqxqy':
    for ext in 'json', 'png', 'txt':
        EXTENSIONS['CG3_1433'].append('{}.{}'.format(iq_type, ext))


def write_configfile(input_json_file, basename):
    '''
    Create a new json configuration file with a better place for the output files
    and a standardized basename

    Parameters
    ----------
    input_json_file: str
        json file to modify. Expected to be in ``ROOT_DIR/scripts``
    basename: str
        The value of ``outputFilename``. This makes it significantly easier to find the output files.

    Returns
    -------
    str, str
        The name of the output directory and the name of the (re)configured json file
    '''
    # get the full path to a configuration file
    input_json_file = os.path.join(ROOT_DIR, 'scripts', input_json_file)
    assert os.path.exists(input_json_file), 'Could not find "{}"'.format(input_json_file)

    # temporary directory is always readable
    outputdir = gettempdir()
    output_json_file = NamedTemporaryFile(prefix=os.path.basename(input_json_file).replace('.json', '_'),
                                          suffix='.json', delete=False).name
    print('Reconfigured json file set to {}'.format(output_json_file))

    # read the existing file in
    with open(input_json_file, 'r') as handle:
        json_params = json.load(handle)

    # get the cannonical runnumber name to standardize it
    json_params['outputFilename'] = basename
    print('Output basename set to "{}"'.format(json_params['outputFilename']))

    # update the output directory
    json_params["configuration"]["outputDir"] = outputdir
    print('Output directory set to "{}"'.format(outputdir))

    # write out the new file
    with open(output_json_file, 'w') as handle:
        json.dump(json_params, handle)

    return outputdir, output_json_file


def run_reduction(pythonscript, json_file):
    # determine python script with full path
    scriptdir = os.path.join(os.path.abspath(os.path.curdir), 'scripts')
    pythonscript = os.path.join(scriptdir, pythonscript)
    assert os.path.exists(pythonscript), 'Could not find "{}"'.format(pythonscript)

    # run the script
    cmd = 'python {} {}'.format(pythonscript, json_file)
    print('Running "{}"'.format(cmd))
    start = time.clock()
    proc = subprocess.Popen(cmd,
                            shell=True,
                            universal_newlines=True)
    proc.communicate()

    # 0 = ok
    # 42 = non-empty error log file. This is the case for developers running the test on
    # the console without a special configuration of their logging
    returncode = proc.returncode
    assert returncode in [0, 42]
    print(pythonscript, 'took', time.clock() - start, 'seconds')


def check_and_cleanup(outputdir, basename):
    # verify that the output files were created and cleanup
    for extension in EXTENSIONS[basename]:
        filename = os.path.join(outputdir, basename + extension)
        assert os.path.isfile(filename), '"{}" does not exist'.format(filename)
        os.remove(filename)

    # remove the output and error logs if they were created, neither is required to exist
    for ext in ['.out', '.err']:
        logname = os.path.join(outputdir, basename + ext)
        if os.path.isfile(logname):
            os.remove(logname)


@pytest.mark.parametrize('configfile, basename',
                         [('reduction.json', 'EQSANS_88980')],
                         ids=['88980'])
def test_eqsans(configfile, basename):
    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename)

    run_reduction('eqsans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename)


@pytest.mark.parametrize('configfile, basename',
                         [('biosans_reduction.json', 'CG3_1433')],
                         ids=['1433'])
def test_biosans(configfile, basename):
    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename)

    run_reduction('biosans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename)


@pytest.mark.parametrize('configfile, basename',
                         [('gpsans_reduction.json', 'CG2_1481')],
                         ids=['1481'])
def test_gpsans(configfile, basename):
    # test taken from EXAMPLES_GPSANS5 Al4 main detector

    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename)

    run_reduction('gpsans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename)


if __name__ == '__main__':
    pytest.main([__file__])
