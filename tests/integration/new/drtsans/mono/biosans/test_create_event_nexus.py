"""
Integration test to create event nexus file
"""
import pytest
import numpy as np
import os
from drtsans.load import load_events
import json
import h5py
from drtsans.mono.biosans import (load_all_files, reduce_single_configuration,
                                  reduction_parameters, validate_reduction_parameters)
from mantid.simpleapi import LoadEventNexus, Rebin, ConvertToMatrixWorkspace, mtd, LoadHFIRSANS
from drtsans.mono.biosans.cg3_spice_to_nexus import generate_event_nexus
from drtsans.mono.convert_xml_to_nexus import EventNexusConverter
from tempfile import mkdtemp
from matplotlib import pyplot as plt


def test_convert_spice_to_nexus(reference_dir, cleanfile):
    """Test to convert SPICE to NeXus

    Parameters
    ----------
    reference_dir
    cleanfile

    Returns
    -------

    """
    # Specify the test data
    # FIXME in this stage, using data in /HFIR/CG3/; data in reference_dir will be used
    spice_data_file = os.path.join(reference_dir.new.biosans, 'BioSANS_exp327_scan0014_0001.xml')
    template_nexus_file = os.path.join(reference_dir.new.biosans, 'CG3_5705.nxs.h5')
    assert os.path.exists(spice_data_file), f'SPICE file {spice_data_file} cannot be located'
    assert os.path.exists(template_nexus_file), f'Template NeXus file {template_nexus_file} cannot be located'

    # Specify the output directory
    if False:
        # FIXME - enable this section after this test is passed
        output_dir = mkdtemp(prefix='cg3spice2nexus')
        cleanfile(output_dir)
    else:
        # FIXME - remove this section after this test is passed
        import shutil
        output_dir = '/tmp/cg3spice2nexus'
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)

    # output file name
    out_nexus_file = os.path.join(output_dir, 'CG3_32700140001.nxs.h5')

    # set DAS meta data log map
    meta_map = {
        'CG3:CS:SampleToSi': ('sample_to_flange', 'mm'),  # same
        'sample_detector_distance': ('sdd', 'm'),  # same
        'wavelength': ('lambda', 'angstroms'),  # angstroms -> A
        'wavelength_spread': ('dlambda', 'fraction'),  # fraction -> None
        'source_aperture_diameter': ('source_aperture_size', 'mm'),  # same
        'sample_aperture_diameter': ('sample_aperture_size', 'mm'),  # same
        'detector_trans_Readback': ('detector_trans', 'mm'),  # same
        'source_distance': ('ource_aperture_sample_aperture_distance', 'm'),  # same. source-aperture-sample-aperture
        'beamtrap_diameter': ('beamtrap_diameter', 'mm'),  # not there
        'ww_rot_Readback': ('det_west_wing_rot', 'degree')
    }

    # init convert
    converter = EventNexusConverter('CG3', 'CG3')
    converter.load_idf(template_nexus_file)
    converter.load_sans_xml(spice_data_file, meta_map)
    converter.generate_event_nexus(out_nexus_file, num_banks=88)

    # Check: file existence
    os.path.exists(out_nexus_file), f'Output file {out_nexus_file} cannot be located'

    # Check instrument node against the original one
    test_nexus_h5 = h5py.File(out_nexus_file, 'r')
    test_idf = test_nexus_h5['entry']['instrument']['instrument_xml']['data'][0]
    expected_nexus_h5 = h5py.File(template_nexus_file, 'r')
    expected_idf = expected_nexus_h5['entry']['instrument']['instrument_xml']['data'][0]
    assert test_idf == expected_idf
    test_nexus_h5.close()
    expected_nexus_h5.close()

    # Load test data
    test_ws_name = 'TestFaked_CG3_32700140001'
    LoadEventNexus(Filename=out_nexus_file, OutputWorkspace=test_ws_name,
                   NumberOfBins=1, LoadNexusInstrumentXML=True)
    ConvertToMatrixWorkspace(InputWorkspace=test_ws_name, OutputWorkspace=test_ws_name)
    test_nexus_ws = mtd[test_ws_name]

    # Load template event nexus
    LoadEventNexus(Filename=template_nexus_file, OutputWorkspace='cg3template',
                   NumberOfBins=1, LoadNexusInstrumentXML=True)
    template_ws = mtd['cg3template']

    # Check number of histograms
    assert test_nexus_ws.getNumberHistograms() == template_ws.getNumberHistograms()

    # Compare units of required DAS logs
    for das_log_name in ['CG3:CS:SampleToSi', 'wavelength', 'wavelength_spread', 'source_aperture_diameter',
                         'sample_aperture_diameter', 'detector_trans_Readback', 'sample_detector_distance',
                         'detector_trans_Readback', 'www_rot_Readback']:
        template_unit = template_ws.run().getProperty(das_log_name).units
        test_unit = test_nexus_ws.run().getProperty(das_log_name).units
        assert template_unit == test_unit, f'DAS log {das_log_name} unit does not match'

    # Check instrument by comparing pixel position
    # Run 9711: detector_trans_Readback = 0.002 mm (to negative X direction)
    # Exp315 Scan 5  Run 60: detector trans = 0.001 mm
    # Thus all pixels of from-SPICE data shall have a postive 1 mm shift
    # Both data have different SDD.  Thus all the pixels will have a constant shift along Z direction
    diff_x = 0.001
    diff_z_list = list()
    for iws in range(0, template_ws.getNumberHistograms(), 10):
        test_pixel_pos = test_nexus_ws.getDetector(iws).getPos()
        expected_pixel_pos = template_ws.getDetector(iws).getPos()
        # constant difference at x
        assert test_pixel_pos[0] - diff_x == pytest.approx(expected_pixel_pos[0], abs=1e-7)
        # y shall be exactly same
        assert test_pixel_pos[1] == pytest.approx(expected_pixel_pos[1], abs=1e-7)
        # z shall have constant difference
        diff_z_list.append(test_pixel_pos[2] - expected_pixel_pos[2])
    # shift along Z-axis shall be a constant
    assert np.array(diff_z_list).std() < 1E-12

    # Load original SPICE file
    spice_ws_name = os.path.basename(spice_data_file).split('.')[0]
    spice_ws_name = f'CG3IntTestSpice_{spice_ws_name}'
    LoadHFIRSANS(Filename=spice_data_file, OutputWorkspace=spice_ws_name)
    spice_ws = mtd[spice_ws_name]

    # compare histograms
    for iws in range(0, test_nexus_ws.getNumberHistograms()):
        assert test_nexus_ws.readY(iws)[0] == pytest.approx(spice_ws.readY(iws + 2)[0], abs=1E-3)

    # compare DAS logs (partial)
    for log_name in ['wavelength', 'source_aperture_diameter', 'sample_aperture_diameter']:
        nexus_log_value = test_nexus_ws.run().getProperty(log_name).value.mean()
        spice_log_value = spice_ws.run().getProperty(log_name).value
        assert nexus_log_value == pytest.approx(spice_log_value, 1e-7)


def test_duplicate_event_nexus(reference_dir, cleanfile):
    """Test duplicating an HDF5/NeXus in 2 different approaches in order to verify EventNexusWriter

    Verification is to load both of the generated Event NeXus to do a comparison

    Test data: BIOSANS run 5709

    Returns
    -------

    """
    # Get the source file
    source_nexus_file = 'CG3_5709.nxs.h5'
    source_nexus_file = os.path.join(reference_dir.new.biosans, source_nexus_file)
    assert os.path.exists(source_nexus_file), f'Test data {source_nexus_file} does not exist'

    # Duplicate the source file to the temporary directory
    output_dir = mkdtemp(prefix='dupcg3nexus')
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    product_dup_nexus = os.path.join(output_dir, 'CG3_5709_product.nxs.h5')

    # Duplicate with both approach
    logs_white_list = ['CG3:CS:SampleToSi', 'sample_detector_distance',
                       'wavelength', 'wavelength_spread',
                       'source_aperture_diameter', 'sample_aperture_diameter',
                       'detector_trans_Readback', 'ww_rot_Readback',
                       'source_aperture_sample_aperture_distance']
    generate_event_nexus(source_nexus_file, product_dup_nexus, logs_white_list)

    # Load source file to workspace
    target_ws = load_events(product_dup_nexus, output_workspace='cg3_product', NumberOfBins=1,
                            LoadNexusInstrumentXML=True)
    source_ws = load_events(source_nexus_file, output_workspace='cg3_source', NumberOfBins=1,
                            LoadNexusInstrumentXML=True)

    # Compare pixels' positions
    num_hist = source_ws.getNumberHistograms()
    for iws in range(0, num_hist):
        source_det_i_pos = source_ws.getInstrument().getDetector(iws).getPos()
        target_det_i_pos = target_ws.getInstrument().getDetector(iws).getPos()
        np.testing.assert_allclose(source_det_i_pos, target_det_i_pos,
                                   err_msg=f'Mismatch is detected at Detector {iws}')
    # Check source position
    source_moderator_pos = source_ws.getInstrument().getSource().getPos()
    target_moderator_pos = target_ws.getInstrument().getSource().getPos()
    print(f'source moderator @ {source_moderator_pos}; re-generated moderator @ {target_moderator_pos}')
    np.testing.assert_allclose(source_moderator_pos, target_moderator_pos,
                               err_msg=f'Mismatch is detected at neutron source position')

    # Compare counts on each pixel
    source_y = source_ws.extractY()
    target_y = target_ws.extractY()
    np.testing.assert_allclose(source_y, target_y)


def verify_histogram(source_nexus, test_nexus):
    """Check whether two NeXus files can render out same result

    Parameters
    ----------
    source_nexus: str
        source/gold nexus file name
    test_nexus: str
        nexus file to test

    Returns
    -------

    """
    # Compare in HDF5 level
    source_h5 = h5py.File(source_nexus, 'r')
    target_h5 = h5py.File(test_nexus, 'r')

    # Check all the banks
    error_hdf = ''
    for bank_id in range(1, 89):
        source_event_ids = source_h5['entry'][f'bank{bank_id}_events']['event_id'][()]
        target_event_ids = target_h5['entry'][f'bank{bank_id}_events']['event_id'][()]
        if source_event_ids.shape != target_event_ids.shape:
            error_hdf += f'Bank {bank_id}  {source_event_ids.shape} vs {target_event_ids.shape}\n'

    # close nexus files as hdf5
    source_h5.close()
    target_h5.close()

    # Report error
    if len(error_hdf) > 0:
        print(error_hdf)
        print(f'source: {source_nexus}')
        print(f'target: {test_nexus}')
        raise AssertionError(error_hdf)

    # Compare with NeXus
    # Load NeXus file
    src_ws = LoadEventNexus(Filename=source_nexus, OutputWorkspace='gold', NumberOfBins=500)
    test_ws = LoadEventNexus(Filename=test_nexus, OutputWorkspace='test', NumberOfBins=1)

    # Compare with raw counts
    error_message = ''
    for i in range(src_ws.getNumberHistograms()):
        if src_ws.readY(i).sum() != test_ws.readY(i)[0]:
            error_message += f'Workspace-index {i} / detector ID {src_ws.getDetector(i).getID()}/' \
                             f'{test_ws.getDetector(i).getID()}: Expected counts = {src_ws.readY(i)},' \
                             f'Actual counts = {test_ws.readY(i)}\n'

    # report error
    if len(error_message) > 0:
        print(error_message)
        print(f'source: {source_nexus}')
        print(f'target: {test_nexus}')
        raise AssertionError(error_message)

    # A more tricky situation: Rebin throws away events
    src_ws = Rebin(InputWorkspace=src_ws, Params='-20000,40000,20000', PreserveEvents=False)

    # Compare counts
    error_message = ''
    for i in range(src_ws.getNumberHistograms()):
        if src_ws.readY(i)[0] != test_ws.readY(i)[0]:
            error_message += f'Workspace-index {i} / detector ID {src_ws.getDetector(i).getID()}/' \
                             f'{test_ws.getDetector(i).getID()}: Expected counts = {src_ws.readY(i)},' \
                             f'Actual counts = {test_ws.readY(i)}\n'

    # write the error message to disk
    report_file_name = os.path.basename(source_nexus).split('.')[0] + '_error_log.txt'
    with open(report_file_name, 'w') as report_file:
        report_file.write(error_message)
        report_file.write(f'source: {source_nexus}\n')
        report_file.write(f'target: {test_nexus}\n')


def test_reduction(reference_dir, cleanfile):
    """Test generate (partially copy) an event Nexus file by
    verifying reduction result between raw and generated event nexus file

    Testing is modified from mono.gpsans.test_overwrite_geometry_meta_data.test_no_overwrite()

    Returns
    -------

    """
    # Set up test
    json_str = generate_testing_json(os.path.join(reference_dir.new.biosans, 'overwrite_gold_04282020'), None, None)

    # Create output directory
    output_dir = mkdtemp(prefix='nexuscg3reduction')
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Run
    reduce_biosans_data(reference_dir.new.biosans, json_str, output_dir, prefix='BioMetaRaw')

    # Get result files
    sample_names = ['csmb_ecoli1h_n2']
    gold_path = os.path.join(reference_dir.new.biosans, 'overwrite_gold_20200714/test1/')

    # Verify
    try:
        # there are some minor difference because Rebin throws away some events
        verify_reduction_results(sample_names, output_dir, gold_path, title='Raw (no overwriting)', prefix='test1',
                                 rel_tol=3e-6)
    except AssertionError as a_error:
        raise AssertionError(f'Reduction result does not match {a_error}')


def generate_testing_json(sens_nxs_dir, sample_to_si_window_distance, sample_to_detector_distance):
    """Generating testing JSON

    Parameters
    ----------
    sens_nxs_dir: str
        directory path to sensitivity files
    sample_to_si_window_distance: float or None
        sample to silicon window distance to overwrite the EPICS value with unit millimeter
    sample_to_detector_distance: float or None
        sample to silicon window distance to overwrite the EPICS value with unit meter

    Returns
    -------
    str
        JSON string
    """
    specs = {
        "iptsNumber": "23782",
        "sample": {"runNumber": "4822", "thickness": "0.1", "transmission": {"runNumber": "4822"}},
        "outputFileName": "CG3_4822",
        "background": {"runNumber": "4821", "transmission": {"runNumber": "4821"}},
        "beamCenter": {"runNumber": "1322"},
        "emptyTransmission": {"runNumber": "5705"},
        "configuration": {
            "outputDir": "/HFIR/CG3/shared/UserAcceptance/override/test1",
            "sampleApertureSize": "14",
            "useDefaultMask": True,
            "defaultMask": ["{'Pixel':'1-12,244-256'}", "{'Bank':'21-24,45-48'}"],
            "darkMainFileName": "CG3_1383.nxs.h5",
            "darkWingFileName": "CG3_1383.nxs.h5",
            "sensitivityMainFileName": "/HFIR/CG3/shared/Cycle486/sens_f4829m7p0_TDC_SAC.h5",
            "sensitivityWingFileName": "/HFIR/CG3/shared/Cycle486/sens_f4835w3p2_TDC_SAC.h5",
            "absoluteScaleMethod": "standard",
            "StandardAbsoluteScale": "0.0055e-8",
            "numMainQxQyBins": "100",
            "numWingQxQyBins": "100",
            "1DQbinType": "scalar",
            "QbinType": "log",
            "useLogQBinsEvenDecade": False,
            "LogQBinsPerDecadeMain": 20,
            "LogQBinsPerDecadeWing": 25,
            "WedgeMinAngles": "-30, 60",
            "WedgeMaxAngles": "30, 120",
            "AnnularAngleBin": "1",
            "QminMain": "0.003",
            "QminWing": "0.003",
            "overlapStitchQmin": "0.075",
            "overlapStitchQmax": "0.095",
            "useTimeSlice": False,
            "timeSliceInterval": "200",
            "usePixelCalibration": False,
            "useSubpixels": False
        }
    }
    reduction_input = reduction_parameters(specs, 'BIOSANS', validate=False)  # add defaults and defer validation
    reduction_config = reduction_input['configuration']  # a handy shortcut

    #  '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4829m7p0_TDC_SAC.h5'
    main_sens = os.path.join(sens_nxs_dir, 'sens_f4829m7p0_TDC_SAC.h5')
    reduction_config['sensitivityMainFileName'] = main_sens

    # '/SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/meta_overwrite/biosans/sens_f4835w3p2_TDC_SAC.h5'
    wing_sens = os.path.join(sens_nxs_dir, 'sens_f4835w3p2_TDC_SAC.h5')
    reduction_config['sensitivityWingFileName'] = wing_sens

    if sample_to_si_window_distance is not None:
        reduction_config['sampleToSi'] = sample_to_si_window_distance

    if sample_to_detector_distance is not None:
        reduction_config['sampleDetectorDistance'] = sample_to_detector_distance

    return json.dumps(reduction_input)  # return a string representation


def reduce_biosans_data(nexus_dir, json_str, output_dir, prefix):
    """Reduce BIOSANS runs

    Parameters
    ----------
    nexus_dir: str
        path to NeXus files
    json_str: str
        configuration json
    output_dir: str
        output directory
    prefix: str
        prefix for output workspaces' names

    Returns
    -------

    """
    # Set up (testing) runs
    sample_names = ['csmb_ecoli1h_n2']
    sample = '5709'
    samples_tran = sample
    background = '5715'
    backgrounds_trans = background

    # Replace by duplicates
    logs_white_list = ['CG3:CS:SampleToSi', 'sample_detector_distance',
                       'wavelength', 'wavelength_spread',
                       'source_aperture_diameter', 'sample_aperture_diameter',
                       'detector_trans_Readback', 'ww_rot_Readback',
                       'source_aperture_sample_aperture_distance']

    # generate sample
    source_sample_nexus = os.path.join(nexus_dir, f'CG3_{sample}.nxs.h5')
    os.path.exists(source_sample_nexus), f'Source sample NeXus {source_sample_nexus} does not exist'
    test_sample_nexus = os.path.join(output_dir, f'CG3_{sample}.nxs.h5')
    generate_event_nexus(source_sample_nexus, test_sample_nexus, logs_white_list)
    verify_histogram(source_sample_nexus, test_sample_nexus)

    # generate background
    source_bkgd_nexus = os.path.join(nexus_dir, f'CG3_{background}.nxs.h5')
    os.path.exists(source_bkgd_nexus), f'Source background NeXus {source_bkgd_nexus} does not exist'
    test_bkgd_nexus = os.path.join(output_dir, f'CG3_{background}.nxs.h5')
    generate_event_nexus(source_bkgd_nexus, test_bkgd_nexus, logs_white_list)

    # sample trans
    if samples_tran == sample:
        samples_tran = test_sample_nexus

    # background
    if backgrounds_trans == background:
        backgrounds_trans = test_bkgd_nexus

    # checking if output directory exists, if it doesn't, creates the folder
    for subfolder in ['1D', '2D']:
        output_folder = os.path.join(output_dir, subfolder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    # Load JSON for configuration
    # start with a fresh set of reduction parameters for every sample, because the reduction "pollutes"
    # the reduction parameters dictionary with additions not allowed by the schema
    reduction_input = json.loads(json_str)
    reduction_input["dataDirectories"] = nexus_dir
    reduction_input["configuration"]["outputDir"] = output_dir
    reduction_input["sample"]["runNumber"] = test_sample_nexus
    reduction_input["sample"]["transmission"]["runNumber"] = samples_tran   # set to test
    reduction_input["background"]["runNumber"] = test_bkgd_nexus  # backgrounds[0]
    reduction_input["background"]["transmission"]["runNumber"] = backgrounds_trans
    reduction_input["outputFileName"] = sample_names[0]

    # beam center: convert and reset
    beam_center_run = reduction_input['beamCenter']['runNumber']
    source_bc_nexus = os.path.join(nexus_dir, f'CG3_{beam_center_run}.nxs.h5')
    os.path.exists(source_bc_nexus), f'Source background NeXus {source_bc_nexus} does not exist'
    test_bc_nexus = os.path.join(output_dir, f'CG3_{beam_center_run}.nxs.h5')
    # specific log: there is no source_aperture_diameter in run 1322
    das_1322_list = ['CG3:CS:SampleToSi', 'sample_detector_distance', 'wavelength', 'wavelength_spread',
                     'detector_trans_Readback', 'ww_rot_Readback',
                     'source_aperture_sample_aperture_distance']
    generate_event_nexus(source_bc_nexus, test_bc_nexus, das_1322_list)
    reduction_input['beamCenter']['runNumber'] = test_bc_nexus

    # always check after updating the parameters
    reduction_input = validate_reduction_parameters(reduction_input)
    loaded = load_all_files(reduction_input,
                            path=nexus_dir,
                            prefix=prefix)

    # reduce
    reduce_single_configuration(loaded, reduction_input)


def verify_reduction_results(sample_names, output_dir, gold_path, title, prefix, rel_tol=1e-7):

    unmatched_errors = ''

    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_test1_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
        # compare
        title_i = '{}: {}'.format(sample_name, title)
        try:
            compare_reduced_iq(output_log_file, gold_log_file, title_i, prefix, rel_tol=rel_tol)
        except AssertionError as unmatched_error:
            unmatched_errors = 'Testing output {} is different from gold result {}:\n{}' \
                               ''.format(output_log_file, gold_log_file, unmatched_error)
    # END-FOR

    # raise error for all
    if unmatched_errors != '':
        raise AssertionError(unmatched_errors)


def compare_reduced_iq(test_log_file, gold_log_file, title, prefix, rel_tol=1e-7):
    """

    Parameters
    ----------
    test_log_file: str
        Absolute
    gold_log_file: str
    title: str
        plot title
    prefix: str
        file name prefix
    rel_tol: float
        relative tolerance of the difference between expected value and actual value

    Returns
    -------

    """
    log_errors = list()

    for is_main_detector in [True, False]:
        assert os.path.exists(gold_log_file)
        vec_q_a, vec_i_a = get_iq1d(test_log_file, is_main=is_main_detector)
        vec_q_b, vec_i_b = get_iq1d(gold_log_file, is_main=is_main_detector)

        try:
            np.testing.assert_allclose(vec_q_a, vec_q_b, rtol=rel_tol)
            np.testing.assert_allclose(vec_i_a, vec_i_b, rtol=rel_tol)
            log_errors.append(None)
        except AssertionError as assert_err:
            if is_main_detector:
                flag = 'Main_detector'
            else:
                flag = 'Wing_detector'
            log_errors.append(f'{flag}: {assert_err}')
            plt.cla()
            plt.plot(vec_q_a, vec_i_a, color='red', label='{} Corrected'.format(flag))
            plt.plot(vec_q_b, vec_i_b, color='black', label='{} Before being corrected'.format(flag))
            plt.yscale('log')
            plt.title(title)
            plt.legend()
            out_name = prefix + '_' + os.path.basename(test_log_file).split('.')[0] + '_{}.png'.format(flag)
            plt.savefig(out_name)
    # END-FOR

    # Report
    if not (log_errors[0] is None and log_errors[1] is None):
        error_message = 'Main: {}; Wing: {}'.format(log_errors[0], log_errors[1])
        raise AssertionError(error_message)


def get_iq1d(log_file_name, is_main=True):
    """

    Parameters
    ----------
    log_file_name: str
        output log file's name
    is_main: bool
        for main or wing

    Returns
    -------

    """
    # Open file and entry
    assert os.path.exists(log_file_name)
    log_h5 = h5py.File(log_file_name, 'r')
    try:
        if is_main:
            iq1d_entry = log_h5['main_0']['I(Q)']
        else:
            iq1d_entry = log_h5['wing_0']['I(Q)']
    except KeyError:
        if is_main:
            iq1d_entry = log_h5['_slice_1']['main_0']['I(Q)']
        else:
            iq1d_entry = log_h5['_slice_1']['wing_0']['I(Q)']

    # Get data with a copy
    vec_q = np.copy(iq1d_entry['Q'][()])
    vec_i = np.copy(iq1d_entry['I'][()])

    # close file
    log_h5.close()

    return vec_q, vec_i


if __name__ == '__main__':
    pytest.main(__file__)
