import json
import os
import pytest
import subprocess
from tempfile import gettempdir, NamedTemporaryFile
import time

# this should point to the root directory of the code repository
ROOT_DIR = os.path.abspath(os.path.join(__file__, '../../../../'))
# specific extensions for given basenames
EXTENSIONS = {'EQSANS_88980': ['_bkgd_88974_trans.txt',
                               '_bkgd_88974_raw_trans.txt',
                               '_trans.txt',
                               '_raw_trans.txt',
                               '.nxs',
                               '_wl_frame_0_2D.dat',
                               '_1D_wl_frame_0.dat',
                               '_wl_frame_1_2D.dat',
                               '_1D_wl_frame_1.dat',
                               '_reduction_log.hdf',
                               '_0_2D.png',
                               '_0_1D.png',
                               '_1_2D.png',
                               '_1_1D.png'],
              'EQSANS_112300': ['_1D.dat',
                                '_1D.png',
                                '_2D.dat',
                                '_2D.png',
                                '_bkgd_112296_raw_trans.txt',
                                '_bkgd_112296_trans.txt',
                                '.nxs',
                                '_raw_trans.txt',
                                '_reduction_log.hdf',
                                '_trans.txt'],
              'CG3_1433': ['_merged_Iq.json', '_merged_Iq.png'],
              'CG2_1481': ['_Iq.json', '_Iq.png', '_Iq.txt', '_Iqxqy.json', '_Iqxqy.png', '_Iqxqy.txt']}
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
    cmd = 'python3 {} {}'.format(pythonscript, json_file)
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


def check_for_required_files(filenames):
    '''Skip the test if a required file is missing'''
    for filename in filenames:
        if not os.path.exists(filename):
            pytest.skip('"{}" is not available'.format(filename))


@pytest.mark.parametrize('configfile, basename, required',
                         [('reduction.json', 'EQSANS_88980',
                           ('/SNS/EQSANS/IPTS-19800/nexus/EQSANS_88980.nxs.h5', )),
                          ('eqsans_reduction.json', 'EQSANS_112300',
                           ('/SNS/EQSANS/IPTS-24769/nexus/EQSANS_112300.nxs.h5', ))],
                         ids=['88980', '112300'])
def test_eqsans(configfile, basename, required):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename)

    run_reduction('eqsans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


@pytest.mark.parametrize('configfile, basename, required',
                         [('biosans_reduction.json', 'CG3_1433',
                           ('/HFIR/CG3/IPTS-24665/nexus/CG3_1433.nxs.h5', ))],
                         ids=['1433'])
def test_biosans(configfile, basename, required):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename)

    run_reduction('biosans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


@pytest.mark.parametrize('configfile, basename, required',
                         [('gpsans_reduction.json', 'CG2_1481',
                           ('/HFIR/CG2/IPTS-23801/nexus/CG2_1481.nxs.h5', ))],
                         ids=['1481'])
def test_gpsans(configfile, basename, required):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    outputdir, json_file = write_configfile(configfile, basename)

    run_reduction('gpsans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


if __name__ == '__main__':
    pytest.main([__file__])
