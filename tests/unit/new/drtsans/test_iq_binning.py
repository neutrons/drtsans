import numpy as np
import pytest
# https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html
from mantid.simpleapi import LoadEmptyInstrument, AddSampleLog  # AddTimeSeriesLog, Rebin, ConvertUnits,
from drtsans.iq import bin_iq_into_linear_q1d, bin_iq_into_logarithm_q1d, IofQCalculator
from drtsans.momentum_transfer_factory import calculate_q_dq

# This test implements issue #169 to verify
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/tree/169_bin_q1d
# DEV - Wenduo Zhou <petersonpf@ornl.gov> and Joe Osborn <osbornjd@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>, Shuo, Lisa

# this array is a copy & paste from Lisa's PDF for counts on detectors in the detector view
det_view_counts = np.array([[1,      0,   1,   1,   1,   0,   1,   1, 1],
                            [1,     87, 130, 208, 121, 107, 279, 186, 0],
                            [1,    172, 250, 500, 479, 512, 488, 208, 1],
                            [2,    147, 600,   5,  10,   5, 700, 189, 3],
                            [1,    239, 562,  10,   0,  10, 800, 217, 1],
                            [3, np.nan, 550,   5,  10,   5, 689, 228, 1],
                            [1,    178, 567, 503, 469, 489, 499, 156, 0],
                            [0,    108, 350, 389, 409, 253, 192, 209, 1],
                            [1,     80, 193, 120, 148, 108, 201, 210, 1],
                            [1,      3,   2,   0,   1,   1,   0,   1, 1]], dtype=np.float)

# Define some constants
sdd = 5.  # meter
x_pixel_size = 5.5 * 1.E-3  # meter
y_pixel_size = 4.1 * 1.E-3  # meter
wavelength = 6  # angstroms
delta_lambda = 0.15  # angstroms
x_beam_center = 5*x_pixel_size - 0.02749  # meter
y_beam_center = 5.5*y_pixel_size - 0.02059  # meter
R1 = 0.02  # source aperture radius
R2 = 0.007  # sample aperture radius


# Make a mantid workspace for the intensity
@pytest.mark.parametrize('generic_IDF',
                         [{'name': 'GenericSANS',
                           'l1': -15.,
                           'Nx': det_view_counts.shape[1],  # 9
                           'Ny': det_view_counts.shape[0],  # 10
                           'dx': x_pixel_size,
                           'dy': y_pixel_size,
                           'xc': x_beam_center,
                           'yc': y_beam_center,
                           'zc': sdd}],
                         indirect=True)
def test_binning_1d(generic_IDF):
    """ Test binning I(Q) in 1D.
    This test shall verify
    1. Q for each pixel
    2. I(Q)

    Parameters
    ----------
    generic_IDF : fixture of workspace

    Returns
    -------

    """
    with open(r'/tmp/GenericSANS_Definition.xml', 'w') as tmp:
        tmp.write(generic_IDF)
        tmp.close()
    workspace = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS',
                                    OutputWorkspace='GenericMonoSANS')
    workspace.getAxis(0).setUnit('Wavelength')

    # Set value
    wave_length_range = [wavelength - 0.5 * delta_lambda,  wavelength + 0.5*delta_lambda]
    det_counts = det_view_counts[:, ::-1]
    det_counts = det_counts.flatten('F')
    det_counts_error = np.sqrt(det_counts)

    # User golden data
    # expected_q_array = golden_q_array[:, ::-1].flatten('F')

    # assume that the TOF is already frame corrected
    for i in range(workspace.getNumberHistograms()):
        workspace.dataX(i)[:] = wave_length_range  # microseconds
        workspace.dataY(i)[0] = det_counts[i]
        workspace.dataE(i)[0] = det_counts_error[i]
    # END-FOR

    # # Create a single workspace workspace
    # workspace = workspace_with_instrument(axis_values=[wavelength - 0.5*delta_lambda,
    #                                                    wavelength + 0.5*delta_lambda],
    #                                       intensities=det_view_counts)

    print(det_view_counts.shape)
    print(det_view_counts[0])
    print(x_beam_center, y_beam_center)
    # for i in range(90):
    #    det_id = i   # 80 + i  # i * 10
    #    print('Det {}: {} = {} (count)'.format(det_id, workspace.getDetector(det_id).getPos(),
    #                                        workspace.readY(det_id)[0]))

    # Verify detector position
    det_80_pos = workspace.getDetector(80).getPos()
    assert abs(det_80_pos[0] - (5.5 - 27.49) * 0.001) < 1E-6, 'Pixel 80 X-distance to center. Expected to be {}.' \
                                                              'Test is {}'.format((5.5 - 27.49)*0.001, det_80_pos[0])

    assert abs(det_80_pos[1] - (4.1 - 20.59) * 0.001) < 1E-6, 'Pixel 80 Y-distance to center. Expected to be {}.' \
                                                              'Test is {}'.format((4.1 - 20.58)*0.001, det_80_pos[1])

    # Add sample logs for Q resolution
    AddSampleLog(Workspace=workspace, LogName='wavelength', LogText='{}'.format(wavelength),
                 LogType='Number', LogUnit='A')
    AddSampleLog(Workspace=workspace, LogName='wavelength-spread', LogText='{}'.format(delta_lambda),
                 LogType='Number', LogUnit='A')
    AddSampleLog(Workspace=workspace, LogName='source-aperture-diameter', LogText='{}'.format(R1*2.*1000),
                 LogType='Number', LogUnit='mm')
    AddSampleLog(Workspace=workspace, LogName='sample-aperture-diameter', LogText='{}'.format(R2*2.*1000),
                 LogType='Number', LogUnit='mm')

    # Calculate Q and dQ
    q_dq = calculate_q_dq(workspace, instrument_type='mono')
    q_array = q_dq.q
    print('Max Q = {} @ {}'.format(np.max(q_array), np.argmax(q_array)))
    assert q_array.shape == (90, 1)
    # assert abs(q_array[80] - expected_q_array[80]) < 10000000
    # assert abs(np.max(q_array) - 0.013434*0.5) < 0.0001

    # Test the linear binning
    # Can pass the instrument type or the function will grab it from the workspace
    result = bin_iq_into_linear_q1d(workspace, bins=10, q_min=0, q_max=None, instrument='mono')

    # Test instrument geometry
    # assert that returned workspace has some binning scheme as from Lisa's PDF
    assert result is not None

    # Test the logarithm binning
    log_result = bin_iq_into_logarithm_q1d(workspace, bins_per_decade=33,
                                           q_min=0.001, q_max=1.0, instrument='mono')

    assert log_result is not None

    return


def test_generate_log_bins():
    """
    Unit test for the method to generate logarithm bins
    Returns
    -------
    None
    """
    q_min = 0.001
    q_max = 1.
    step_per_decade = 33
    bin_centers, bin_edges = IofQCalculator.determine_log_bin_edges(q_min, q_max, step_per_decade)

    # Verify: bin size, min and max
    assert bin_edges.shape[0] == bin_centers.shape[0] + 1
    assert bin_centers.shape[0] == 100
    assert abs(bin_centers[0] - q_min) < 1.E-12
    assert abs(bin_centers[99] - q_max) < 1.E-12

    return


if __name__ == "__main__":
    pytest.main()
