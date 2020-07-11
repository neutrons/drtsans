"""
Integration test to create event nexus file
"""
import pytest
import numpy as np
import os
from drtsans.load import load_events
import h5py
from drtsans.mono.gpsans import (load_all_files, plot_reduction_output, reduce_single_configuration,
                                 reduction_parameters, update_reduction_parameters)
from drtsans.files.hdf5_rw import GroupNode, DataSetNode
from drtsans.files.event_nexus_nodes import InstrumentNode, DasLogNode, BankNode, MonitorNode
from drtsans.files.event_nexus_rw import generate_events_from_histogram
from drtsans.files.event_nexus_rw import generate_monitor_events_from_count
from drtsans.files.event_nexus_rw import init_event_nexus, parse_event_nexus, EventNeXusWriter
from mantid.simpleapi import LoadEventNexus


def test_duplicate_event_nexus(reference_dir, cleanfile):
    """Test duplicating an HDF5/NeXus in 2 different approaches in order to verify EventNexusWriter

    Verification is to load both of the generated Event NeXus to do a comparison

    Test data: BIOSANS run 5709

    Returns
    -------

    """
    # Get the source file
    source_nexus_file = 'CG3_5709.nxs.h5'
    source_nexus_file = os.path.join(reference_dir.new.bioans, source_nexus_file)
    assert os.path.exists(source_nexus_file), f'Test data {source_nexus_file} does not exist'

    # Duplicate the source file to the temporary directory
    # TODO - this will be replaced by tempfile for future
    output_dir = '/tmp/dupcg3nexus'
    cleanfile(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    prototype_dup_nexus = os.path.join(output_dir, 'CG2_5709_prototype.nxs.h5')
    product_dup_nexus = os.path.join(output_dir, 'CG2_5709_product.nxs.h5')

    # Duplicate with both approach
    generate_event_nexus_prototype(source_nexus_file, prototype_dup_nexus)
    generate_event_nexus(source_nexus_file, product_dup_nexus)

    # Load source file to workspace
    target_ws = load_events(product_dup_nexus, output_workspace='cg3_product', NumberOfBins=2)

    # Load the duplicated
    prototype_ws = load_events(prototype_dup_nexus, output_workspace='cg3_prototype', NumberOfBins=2)

    # Compare pixels' positions
    num_hist = prototype_ws.getNumberHistograms()
    for iws in range(0, num_hist, 100):
        source_det_i_pos = prototype_ws.getInstrument().getDetector(iws).getPos()
        target_det_i_pos = target_ws.getInstrument().getDetector(iws).getPos()
        np.testing.assert_allclose(source_det_i_pos, target_det_i_pos,
                                   err_msg=f'Mismatch is detected at Detector {iws}')
    # Check source position
    source_moderator_pos = prototype_ws.getInstrument().getSource().getPos()
    target_moderator_pos = target_ws.getInstrument().getSource().getPos()
    np.testing.assert_allclose(source_moderator_pos, target_moderator_pos,
                               err_msg=f'Mismatch is detected at neutron source position')

    # Compare counts on each pixel
    source_y = prototype_ws.extractY()
    target_y = target_ws.extractY()
    np.testing.assert_allclose(source_y, target_y)

    # Compare meta data
    assert len(prototype_ws.getRun().getProperties()) == len(target_ws.getRun().getProperties()), 'Meta data mismatch'


def generate_event_nexus(source_nexus, target_nexus):
    """Generate event NeXus properly

    Parameters
    ----------
    source_nexus
    target_nexus

    Returns
    -------

    """
    cg3_num_banks = 88

    # Import essential experimental data from source event nexus file
    nexus_contents = parse_event_nexus(source_nexus, num_banks=48)
    # Generate event nexus writer
    event_nexus_writer = EventNeXusWriter(beam_line='CG2', instrument_name='CG2')

    # set instrument: 88 banks (2 detectors)
    event_nexus_writer.set_instrument_info(cg3_num_banks,  nexus_contents[0])

    # set counts: 88 banks (2 detectors)
    for bank_id in range(1, cg3_num_banks + 1):
        event_nexus_writer.set_bank_histogram(bank_id, nexus_contents[1][bank_id])

    # set meta
    for das_log in nexus_contents[5].values():
        event_nexus_writer.set_meta_data(das_log)

    # time
    start_time = nexus_contents[3]
    end_time = nexus_contents[4]

    # Write file
    event_nexus_writer.generate_event_nexus(target_nexus, start_time, end_time, nexus_contents[2])


def generate_event_nexus_prototype(source_nexus, target_nexus):
    """Generate event NeXus using white list.

    This serves as the prototype to create event nexus from SANS histogram raw data

    White list
    Entry attributes: {'NX_class': b'NXentry'}

    White List Node: /entry/monitor1
    White List Node: /entry/proton_charge

    White List Node: /entry/duration
    White List Node: /entry/start_time
    White List Node: /entry/end_time

    White List Node: /entry/experiment_identifier
    White List Node: /entry/experiment_title
    White List Node: /entry/title
    White List Node: /entry/notes
    White List Node: /entry/raw_frames
    White List Node: /entry/run_number

    White List Node: /entry/total_counts
    White List Node: /entry/total_other_counts
    White List Node: /entry/total_pulses

    Parameters
    ----------
    source_nexus: str
        source event NeXus file name
    target_nexus

    Returns
    -------

    """
    # parse nexus information
    nexus_contents = parse_event_nexus(source_nexus, num_banks=48)

    # Create new nexus file structure
    target_nexus_root = init_event_nexus()

    target_entry_node = target_nexus_root.get_child('entry', is_short_name=True)

    # set instrument node
    set_instrument_node(nexus_contents[0], target_entry_node)

    # set DAS logs
    set_das_log_node(nexus_contents[5], nexus_contents[3], target_entry_node)

    # Add node on the white list
    entry_level_white_list = [
        ('/entry/start_time', nexus_contents[3]),
        ('/entry/end_time', nexus_contents[4])
    ]
    for child_node_name, child_value in entry_level_white_list:
        child_node = DataSetNode(child_node_name)
        child_node.set_string_value(child_value)
        target_entry_node.set_child(child_node)

    # set Bank 1 - 88 (2 detectors)
    max_pulse_time_array = None
    for bank_id in range(1, 88 + 1):
        bank_node_i = set_single_bank_node(nexus_contents[1][bank_id], target_entry_node, bank_id=bank_id,
                                           run_start_time=nexus_contents[3])
        event_time_zeros = bank_node_i.get_child('event_time_zero', is_short_name=True).value
        if max_pulse_time_array is None or event_time_zeros.shape[0] > max_pulse_time_array.shape[0]:
            max_pulse_time_array = event_time_zeros

    # Set monitor node
    set_monitor_node(nexus_contents[2],  nexus_contents[3], target_entry_node, max_pulse_time_array)

    # write
    target_nexus_root.write(target_nexus)


def set_monitor_node(monitor_counts, run_start_time, target_entry_node, event_time_zeros):
    """

    Parameters
    ----------
    monitor_counts: float, int
    target_entry_node
    event_time_zeros: ~numpy.ndarray
        event time zeros
    run_start_time: str, Bytes
        run start time

    Returns
    -------

    """
    # Generate a monitor node
    target_monitor_node = MonitorNode('/entry/monitor1', 'monitor1')

    tof_min = 0.
    tof_max = 10000.
    monitor_events = generate_monitor_events_from_count(monitor_counts, event_time_zeros, tof_min, tof_max)

    target_monitor_node.set_monitor_events(event_index_array=monitor_events.event_index,
                                           event_time_offset_array=monitor_events.event_time_offset,
                                           run_start_time=run_start_time,
                                           event_time_zero_array=event_time_zeros)

    target_entry_node.set_child(target_monitor_node)


def set_single_bank_node(bank_histogram, target_entry_node, bank_id, run_start_time):
    """Test writing bank 9 from histogram

    Parameters
    ----------
    bank_histogram: TofHistogram
        HDF5 file entry
    target_entry_node: GroupNode
        Target (output) group node for /entry/
    bank_id: int
        bank ID (from 1 to 48)
    run_start_time: str, Bytes
        run start time

    Returns
    -------
    BankNode
        newly generated bank node

    """
    # generate events
    nexus_events = generate_events_from_histogram(bank_histogram, 10.)

    try:
        run_start_time = np.string_(run_start_time).decode()
    except AttributeError:
        pass

    # Create bank node for bank
    bank_node = BankNode(name=f'/entry/bank{bank_id}_events', bank_name=f'bank{bank_id}')
    bank_node.set_events(nexus_events.event_id, nexus_events.event_index,
                         nexus_events.event_time_offset, run_start_time,
                         nexus_events.event_time_zero)

    # Link with parent
    target_entry_node.set_child(bank_node)

    return bank_node


def set_instrument_node(xml_idf, target_entry_node):
    """Set instrument node

    Parameters
    ----------
    xml_idf:  str
        IDF content
    target_entry_node

    Returns
    -------

    """
    # Create new instrument node
    instrument_node = InstrumentNode()
    target_entry_node.set_child(instrument_node)

    # Set values
    instrument_node.set_idf(xml_idf, idf_type=b'text/xml', description=b'XML contents of the instrument IDF')
    instrument_node.set_instrument_info(target_station_number=1, beam_line=b'CG2', name=b'CG2', short_name=b'CG2')


def set_das_log_node(das_log_dict, run_start_time, target_entry_node):
    """Set DAS log node in a mixed way

    Parameters
    ----------
    das_log_dict: dict
        das log dictionary containing DasLog objects
    run_start_time: str
        run start time
    target_entry_node: GroupNode
        target node

    Returns
    -------

    """
    target_logs_node = GroupNode('/entry/DASlogs')
    target_entry_node.set_child(target_logs_node)
    # add attribute
    target_logs_node.add_attributes({'NX_class': 'NXcollection'})

    for log_name in das_log_dict:
        set_single_log_node(target_logs_node, das_log_dict[log_name], run_start_time)


def set_single_log_node(log_collection_node, das_log, start_time):
    """

    Parameters
    ----------
    log_collection_node
    das_log: DasLog
    start_time: str

    Returns
    -------

    """
    # Set up a DAS log node
    das_log_node = DasLogNode(log_name=f'/entry/DASlogs/{das_log.name}',
                              log_times=das_log.times,
                              log_values=das_log.values,
                              start_time=start_time,
                              log_unit=das_log.unit)

    if das_log.device is not None:
        if das_log.device.target is None:
            device_target = None
        else:
            device_target = das_log.device.target
        das_log_node.set_device_info(device_id=das_log.device.id,
                                     device_name=das_log.device.name,
                                     target=device_target)

    # append to parent node
    log_collection_node.set_child(das_log_node)


def set_sdd_node(log_collection_node, source_h5):
    # Get times and value for /entry/DASlogs/sample_detector_distance
    ssd_entry = source_h5['entry']['DASlogs']['sample_detector_distance']
    ssd_times = ssd_entry['time'].value
    ssd_start_time = ssd_entry['time'].attrs['start']
    ssd_value = ssd_entry['value'].value
    ssd_value_unit = ssd_entry['value'].attrs['units']

    # Set up a DAS log node
    ssd_test_node = DasLogNode(log_name='/entry/DASlogs/sample_detector_distance',
                               log_times=ssd_times, log_values=ssd_value,
                               start_time=ssd_start_time, log_unit=ssd_value_unit)

    ssd_test_node.set_device_info(device_id=13, device_name=b'Mot-Galil3',
                                  target=b'/entry/DASlogs/CG2:CS:SampleToDetRBV')

    # append to parent node
    log_collection_node.set_child(ssd_test_node)


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
    # Load NeXus file
    src_ws = LoadEventNexus(Filename=source_nexus, OutputWorkspace='gold', NumberOfBins=1)
    test_ws = LoadEventNexus(Filename=test_nexus, OutputWorkspace='test', NumberOfBins=1)

    # Compare counts
    error_message = ''
    for i in range(src_ws.getNumberHistograms()):
        if src_ws.readY(i)[0] != test_ws.readY(i)[0]:
            error_message += f'Workspace-index {i} / detector ID {src_ws.getDetector(i).getID()}/' \
                             f'{test_ws.getDetector(i).getID()}: Expected counts = {src_ws.readY(i)},' \
                             f'Actual counts = {test_ws.readY(i)}\n'
    if error_message != '':
        print(error_message)
        raise AssertionError(error_message)


def not_now_test_reduction(reference_dir, cleanfile):
    """Test generate (partially copy) an event Nexus file by
    verifying reduction result between raw and generated event nexus file

    Testing is modified from mono.gpsans.test_overwrite_geometry_meta_data.test_no_overwrite()

    Returns
    -------

    """
    assert reference_dir
    assert cleanfile
    raise NotImplementedError('Next Step')


def verify_reduction_results(sample_names, output_dir, gold_path, title, prefix):

    unmatched_errors = ''

    for sample_name in sample_names:
        # output log file name
        output_log_file = os.path.join(output_dir, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(output_log_file), 'Output {} cannot be found'.format(output_log_file)
        # gold file
        gold_log_file = os.path.join(gold_path, '{}_reduction_log.hdf'.format(sample_name))
        assert os.path.exists(gold_path), 'Gold file {} cannot be found'.format(gold_log_file)
        # compare
        title_i = '{}: {}'.format(sample_name, title)
        try:
            compare_reduced_iq(output_log_file, gold_log_file, title_i, prefix)
        except AssertionError as unmatched_error:
            unmatched_errors = 'Testing output {} is different from gold result {}:\n{}' \
                               ''.format(output_log_file, gold_log_file, unmatched_error)
    # END-FOR

    # raise error for all
    if unmatched_errors != '':
        raise AssertionError(unmatched_errors)


def compare_reduced_iq(test_log_file, gold_log_file, title, prefix):
    """Compare I(Q) from reduced file and gold file

    Parameters
    ----------
    test_log_file: str
        name of drtsans result log file
    gold_log_file
    title: str
        title of output figure
    prefix: str
        prefix of output file

    Returns
    -------

    """
    # Plot main
    test_q_vec, test_intensity_vec = get_iq1d(test_log_file)
    gold_q_vec, gold_intensity_vec = get_iq1d(gold_log_file)

    # Verify result
    try:
        np.testing.assert_allclose(test_q_vec, test_q_vec, atol=1E-4)
        np.testing.assert_allclose(test_intensity_vec, gold_intensity_vec, atol=1E-7)
    except AssertionError as assert_err:
        from matplotlib import pyplot as plt
        plt.cla()
        plt.plot(test_q_vec, test_intensity_vec, color='red', label='Corrected')
        plt.plot(gold_q_vec, gold_intensity_vec, color='black', label='Before being corrected')
        plt.legend()
        plt.title(title)
        plt.yscale('log')
        out_name = prefix + '_' + os.path.basename(test_log_file).split('.')[0] + '.png'
        plt.savefig(out_name)

        raise assert_err


def get_iq1d(log_file_name):
    """Get I(Q) from output SANS log file

    Parameters
    ----------
    log_file_name: str
        log file name

    Returns
    -------
    tuple
        numpy 1D array for Q, numpy 1D array for intensity

    """
    # Open file and entry
    log_h5 = h5py.File(log_file_name, 'r')

    if '_slice_1' in log_h5:
        data_entry = log_h5['_slice_1']['main']
    else:
        data_entry = log_h5['main']

    # Get data
    iq1d_entry = data_entry['I(Q)']

    # Get data with a copy
    vec_q = np.copy(iq1d_entry['Q'][()])
    vec_i = np.copy(iq1d_entry['I'][()])

    # close file
    log_h5.close()

    return vec_q, vec_i


if __name__ == '__main__':
    pytest.main(__file__)
