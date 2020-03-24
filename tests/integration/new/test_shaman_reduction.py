import json
import os
import pytest
import subprocess
from tempfile import NamedTemporaryFile
import time

# this should point to the root directory of the code repository
ROOT_DIR = os.path.abspath(os.path.join(__file__, '../../../../'))
# specific extensions for given basenames
EXTENSIONS = {'EQSANS_88980': ['_bkgd_88974_trans.txt',
                               '_bkgd_88974_raw_trans.txt',
                               '_trans.txt',
                               '_raw_trans.txt',
                               '.nxs',
                               '_frame_0_Iqxqy.dat',
                               '_frame_0_Iq.dat',
                               '_frame_1_Iqxqy.dat',
                               '_frame_1_Iq.dat',
                               '_reduction_log.hdf',
                               '_0_Iqxqy.png',
                               '_0_Iq.png',
                               '_1_Iqxqy.png',
                               '_1_Iq.png'],
              'EQSANS_112300': ['_frame_0_Iq.dat',
                                '_Iq.png',
                                '_frame_0_Iqxqy.dat',
                                '_Iqxqy.png',
                                '_bkgd_112296_raw_trans.txt',
                                '_bkgd_112296_trans.txt',
                                '.nxs',
                                '_raw_trans.txt',
                                '_reduction_log.hdf',
                                '_trans.txt'],
              'CG3_4822': [],
              'CG3_5532': [],  # auto-wedge
              'CG2_8944': []}


def write_configfile(input_json_file, basename, tmpdir):
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
    outputdir = str(tmpdir)
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


def unique_basename(basename):
    return os.path.basename(NamedTemporaryFile(prefix=basename+'_', delete=False).name)


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


def check_and_cleanup(outputdir, basename_short, basename):
    # verify that the output files were created and cleanup
    for extension in EXTENSIONS[basename_short]:
        filename = os.path.join(outputdir, basename + extension)
        assert os.path.isfile(filename), '"{}" does not exist'.format(filename)
        os.remove(filename)

    # remove the output and error logs if they were created, neither is required to exist
    for ext in ['.out', '.err', '_reduction_log.hdf']:
        logname = os.path.join(outputdir, basename + ext)
        if os.path.isfile(logname):
            os.remove(logname)

    # verify that there isn't anything else with the same prefix
    extra_files = [os.path.join(outputdir, name) for name in os.listdir(outputdir) if basename in name]
    nonempty_files = []
    for filename in extra_files:
        if os.path.isfile(filename) and os.stat(filename).st_size == 0:
            os.remove(filename)
        else:
            nonempty_files.append(filename)
    assert not nonempty_files, 'Found extra files {}'.format(','.join(nonempty_files))


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
def test_eqsans(configfile, basename, required, tmpdir):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    basename_long = unique_basename(basename)
    outputdir, json_file = write_configfile(configfile, basename_long, tmpdir)
    run_reduction('eqsans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename, basename_long)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


@pytest.mark.parametrize('configfile, basename, required',
                         [('biosans_reduction.json', 'CG3_4822', ('/HFIR/CG3/IPTS-23782/nexus/CG3_4822.nxs.h5', )),
                          ('biosans_autowedge_reduction.json', 'CG3_5532',
                           ('/HFIR/CG3/IPTS-21089/nexus/CG3_5532.nxs.h5', ))],
                         ids=['4822', '5532_wedge'])
def test_biosans(configfile, basename, required, tmpdir):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    basename_long = unique_basename(basename)
    outputdir, json_file = write_configfile(configfile, basename_long, tmpdir)

    run_reduction('biosans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename, basename_long)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


@pytest.mark.parametrize('configfile, basename, required',
                         [('gpsans_reduction.json', 'CG2_8944',
                           ('/HFIR/CG2/IPTS-20775/nexus/CG2_8944.nxs.h5', ))],
                         ids=['8944'])
def test_gpsans(configfile, basename, required, tmpdir):
    check_for_required_files(required)

    # modify the config file and get the output directory and the full path to the new configuration file
    basename_long = unique_basename(basename)
    outputdir, json_file = write_configfile(configfile, basename_long, tmpdir)

    run_reduction('gpsans_reduction.py', json_file)

    check_and_cleanup(outputdir, basename, basename_long)

    # delete the modified configuration file
    if os.path.isfile(json_file):
        os.remove(json_file)


if __name__ == '__main__':
    pytest.main([__file__])
