import json
import os
import pytest
import subprocess
from tempfile import gettempdir

# this should point to the root directory of the code repository
ROOT_DIR = os.path.abspath(os.path.join(__file__, '../../../../'))
# extensions that every reduction should output
STD_EXTENSIONS = ['.err', '.out']
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
        Full path to the json file to modify
    basename: str
        The value of ``outputFilename``. This makes it significantly easier to find the output files.

    Returns
    -------
    str, str
        The name of the output directory and the name of the (re)configured json file
    '''
    outputdir = gettempdir()
    output_json_file = os.path.join(outputdir, os.path.basename(input_json_file))

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


@pytest.mark.parametrize('configfile, basename',
                         [('reduction.json', 'EQSANS_88980'),
                          ('biosans_reduction.json', 'CG3_1433')],
                         ids=['EQSANS', 'BIOSANS'])
def test_eqsans(configfile, basename):
    # get the full path to a configuration file
    json_file = os.path.join(ROOT_DIR, 'scripts', configfile)
    assert os.path.exists(json_file), 'Could not find "{}"'.format(json_file)

    outputdir, json_file = write_configfile(json_file, basename)

    # run the script
    cmd = 'process_reduction.py {}'.format(json_file)
    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            universal_newlines=True)
    proc.communicate()

    # 0 = ok
    # 255 = non-empty error log file. This is the case for developers running the test on
    # the console without a special configuration of their logging
    returncode = proc.returncode
    assert returncode in [0, 255]

    print('standard extensions:', STD_EXTENSIONS)
    print('other extensions:', EXTENSIONS[basename])
    # verify that the output files were created and cleanup
    for extension in STD_EXTENSIONS + EXTENSIONS[basename]:
        filename = os.path.join(outputdir, basename + extension)
        assert os.path.isfile(filename), '"{}" does not exist'.format(filename)
        os.remove(filename)


if __name__ == '__main__':
    pytest.main([__file__])
