# Integration test for overwriting instrument geometry related meta data for GP-SANS
# Tests are
import pytest
import os, json
from drtsans.mono.gpsans import load_all_files, reduce_single_configuration, plot_reduction_output
import time


def reduce_gpsans_data(json_file, output_dir):
    # USER Input here with scan numbers etc.
    samples = ['9166', '9167', '9176']
    samples_trans = ['9178', '9179', '9188']
    sample_thick = ['0.1'] * 3
    bkgd = ['9165', '9165', '9165']
    bkgd_trans = ['9177', '9177', '9177']

    sample_names = ["Al4", "PorasilC3", "PTMA-15"]
    # Test JSON
    # json_file = '/HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test2.json'

    # Import JSON
    with open(json_file) as f:
        reduction_input = json.load(f)

    output_dir = reduction_input["configuration"]["outputDir"]
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    start_time = time.time()
    for i in range(len(samples)):
        reduction_input["runNumber"] = samples[i]
        reduction_input["transmission"]["runNumber"] = samples_trans[i]
        reduction_input["background"]["runNumber"] = bkgd[i]
        reduction_input["background"]["transmission"]["runNumber"] = bkgd_trans[i]
        reduction_input["outputFilename"] = sample_names[i]
        reduction_input["thickness"] = sample_thick[i]
        loaded = load_all_files(reduction_input)
        out = reduce_single_configuration(loaded, reduction_input)
        plot_reduction_output(out, reduction_input, loglog=False)

    end_time = time.time()
    print(end_time - start_time)


def test_no_overwrite():
    """Test reduce 3 sets of data without overwriting

    Returns
    -------

    """
    json_file = '/HFIR/CG2/shared/UserAcceptance/overwrite_meta/gpsans_reduction_test2.json'
    reduce_gpsans_data(json_file, '/tmp/meta_overwrite/test1')


if __name__ == '__main__':
    pytest.main([__file__])

