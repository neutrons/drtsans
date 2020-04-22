import pytest
import json
from mantid.simpleapi import mtd
from drtsans.mono.gpsans import load_all_files
# from drtsans.mono.gpsans import load_and_split


def gpsans_load_all_files(data_dir, json_file, output_dir):
    """Standard reduction workflow

    Parameters
    ----------'
    data_dir
    json_file
    output_dir

    Returns
    -------

    """
    # USER Input here with scan numbers etc.
    samples = ['9166']
    samples_trans = ['9178']
    sample_thick = ['0.1']
    bkgd = ['9165']
    bkgd_trans = ['9177']

    # Sample names for output
    sample_names = ["Al4"]

    # Import JSON
    with open(json_file) as f:
        reduction_input = json.load(f)

    # set output directory
    reduction_input["configuration"]["outputDir"] = output_dir

    reduction_input["runNumber"] = samples[0]
    reduction_input["transmission"]["runNumber"] = samples_trans[0]
    reduction_input["background"]["runNumber"] = bkgd[0]
    reduction_input["background"]["transmission"]["runNumber"] = bkgd_trans[0]
    reduction_input["outputFilename"] = sample_names[0]
    reduction_input["thickness"] = sample_thick[0]
    loaded = load_all_files(reduction_input,
                            prefix='GP_TEST_LOAD',
                            load_params=None,
                            path=data_dir)
    print('Type of loaded = {}'.format(loaded))

    for ws_name in mtd.getObjectNames():
        print(ws_name)


if __name__ == '__main__':
    pytest.main([__file__])
