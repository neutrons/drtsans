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

# From @lzi: qbinning_1D.xlsx
det_view_q = np.array([[0.0057557289, 0.0048832439, 0.0041493828, 0.003639012, 0.00345271, 0.0036405748,
                        0.0041521236, 0.0048867371, 0.0057596805],
                       [0.0052854782, 0.0043190126, 0.0034677269, 0.0028372756, 0.0025940164, 0.0028392797,
                        0.0034710061, 0.0043229619, 0.0052897813],
                       [0.00492126, 0.0038647564, 0.0028822853, 0.0020814781, 0.0017353191, 0.0020842091,
                        0.0028862297, 0.0038691694, 0.0049258812],
                       [0.0046878562, 0.0035627952, 0.0024626644, 0.0014455746, 0.0008766204, 0.0014495043,
                        0.0024672797, 0.0035675817, 0.0046927073],
                       [0.0046052016, 0.0034533152, 0.0023014311, 0.0011495873, 1.80845385701797E-05,
                        0.0011545249, 0.0023063691, 0.0034582533, 0.0046101396],
                       [0.0046812885, 0.0035541491, 0.002450139, 0.0014241316, 0.0008407902,
                        0.0014281203, 0.0024547779, 0.0035589472, 0.0046861464],
                       [0.0049087405, 0.0038488017, 0.0028608564, 0.002051702, 0.0016994889,
                        0.0020544727, 0.0028648303, 0.0038532329, 0.0049133735],
                       [0.0052679865, 0.0042975887, 0.0034410067, 0.002804555, 0.0025581863,
                        0.0028065826, 0.0034443113, 0.0043015577, 0.0052723038],
                       [0.0057343077, 0.0048579767, 0.0041196168, 0.0036050342, 0.00341688,
                        0.0036066118, 0.0041223774, 0.0048614881, 0.0057382741],
                       [0.006283909, 0.0054959302, 0.0048555755, 0.0044273783, 0.0042755688,
                        0.004428663, 0.0048579179, 0.0054990343, 0.0062875287]], dtype=np.float)

# Linear binning
linear_bin_centers = np.ndarray([0.000314376, 0.000943129, 0.001571882, 0.002200635, 0.002829388, 0.003458141,
                                 0.004086894, 0.004715647, 0.005344399, 0.005973152])
linear_bin_right_bound = np.array([0.000628753, 0.001257506, 0.001886259, 0.002515011, 0.003143764, 0.003772517,
                                   0.00440127, 0.005030023, 0.005658776, 0.006287529])
linear_bin_left_bound = np.array([0, 0.000628753, 0.001257506, 0.001886259, 0.002515011, 0.003143764, 0.003772517,
                                  0.00440127 , 0.005030023, 0.005658776])

# Weighted
golden_linear_bin_iq = np.array([0 , 10, 7.4717068928, 438.5442497648, 266.6526914842, 4.7169313371,
                                 1.7008577113, 1.1211574604, 0.9333333333, 1.2], dtype=np.float)
golden_linear_bin_sigmaq = np.array([1, 1.5811388301, 1.1159231524, 6.0452753574, 6.1719722651,
                                     0.4856403678, 0.3764812646, 0.3529490773, 0.3651483717, 0.5477225575],
                                    dtype=np.float)

# Logarithm binning
log_bin_centers = np.array([6.58E-04, 7.05E-04, 7.56E-04, 8.11E-04, 8.70E-04, 9.33E-04, 1.00E-03, 1.07E-03,
                            1.15E-03,
                            1.23E-03, 1.32E-03, 1.42E-03, 1.52E-03, 1.63E-03, 1.75E-03, 1.87E-03, 2.01E-03,
                            2.15E-03,
                            2.31E-03, 2.48E-03, 2.66E-03, 2.85E-03, 3.05E-03, 3.27E-03, 3.51E-03, 3.76E-03,
                            4.04E-03, 4.33E-03,
                            4.64E-03, 4.98E-03, 5.34E-03, 5.72E-03, 6.14E-03], dtype=np.float)


golden_log_weighted_intensity = np.array([0, 0, 0, 0, 10, 0, 0, 10, 10, 0, 0, 5, 0, 0, 473.9472574, 0, 495.9012097,
                                          505.9288538, 660.20558, 567.5506439, 121, 262.7392886, 0, 0, 12.16131385,
                                          2.472004695, 1.989894845, 126.6649782, 1.44, 0.931190118, 0.75, 1.487164355,
                                          1])

golden_log_weighted_sigma = np.array([1, 1, 1, 1, 2.236067977, 1, 1, 3.16227766, 2.236067977, 1, 1, 1.118033989, 1, 1,
                                      15.39394779, 1, 15.74644737, 15.90485545, 18.16873111, 10.65411323, 11,
                                      5.730829877, 1, 1, 0.967205087, 0.703136501, 0.705318163, 4.253821445,
                                      0.489897949, 0.364728885, 0.433012702, 0.704074891, 0.707106781], dtype=np.float)


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


def prepare_test_input_arrays():
    """
    Take the arrays provided by IS in 2D Excel cells and transform to drtsans supported data structure.

    Returns
    -------

    """
    # Transform from 2D array to 1D
    iq_array = det_view_counts.flatten().astype(float)
    sigma_q_array = np.sqrt(iq_array)
    q_array = det_view_q.flatten()
    dq_array = q_array * 0.001  # No test, fake now

    # Set zero counts uncertainties correct to sigma Q
    zero_counts_index = np.where(iq_array < 0.00001)
    sigma_q_array[zero_counts_index] = 1.

    # Remove the nans from the I(Q)
    bad_indices = np.isnan(iq_array)
    good_indices = ~bad_indices

    q_array = q_array[good_indices]
    dq_array = dq_array[good_indices]
    iq_array = iq_array[good_indices]
    sigma_q_array = sigma_q_array[good_indices]

    return q_array, dq_array, iq_array, sigma_q_array


def test_linear_binning():
    """

    Returns
    -------

    """
    # Define target Q range
    q_min = 0
    q_max = 0.00629
    assert abs(q_max - np.max(det_view_q)) < 1E-5
    bins = 10

    # Prepare inputs
    q_array, dq_array, iq_array, sigma_q_array = prepare_test_input_arrays()

    # Bin
    bin_centers, bin_edges = IofQCalculator.determine_linear_bin_edges(q_min, q_max, bins)
    assert bin_centers.shape == (10, )
    assert bin_edges[0] == 0.
    assert abs(bin_centers[9] - 0.0059731523) < 0.00005

    binned_q = IofQCalculator.weighted_binning(q_array, dq_array, iq_array, sigma_q_array, bin_centers, bin_edges)

    # Test for Q bins
    assert binned_q.q.shape == (10, )
    assert pytest.approx(binned_q.q[0], q_max/bins * 0.5, 1E-5)

    # Test for I(Q)
    for i in range(10):
        print('Q[{}]: I = {}, gold I = {}, diff = {}'.format(i, binned_q.i[i], golden_linear_bin_iq[i],
                                                             binned_q.i[i] - golden_linear_bin_iq[i]))
        # assert abs(binned_q.i[i] - golden_linear_bin_iq[i]) < 1E-5
        assert pytest.approx(binned_q.sigma[i], golden_linear_bin_sigmaq[i], 1E-5)
    assert False

    return


def test_log_binning():
    """
    Unit test for the method to generate logarithm bins
    Returns
    -------
    None
    """
    # Get test Q, dQ, I, sigmaI
    q_array, dq_array, iq_array, sigma_q_array = prepare_test_input_arrays()

    # Set logarithm binning
    q_min = 0.001
    q_max = 1.
    step_per_decade = 33
    bin_centers, bin_edges = IofQCalculator.determine_log_bin_edges(q_min, q_max, step_per_decade)

    # Verify: bin size, min and max
    assert bin_edges.shape[0] == bin_centers.shape[0] + 1
    assert bin_centers.shape[0] == 100
    assert abs(bin_centers[0] - q_min) < 1.E-12
    assert abs(bin_centers[99] - q_max) < 1.E-12

    # Bin
    binned_q = IofQCalculator.weighted_binning(q_array, dq_array, iq_array, sigma_q_array, bin_centers, bin_edges)

    # # Test for I(Q)
    num_test_points = golden_log_weighted_intensity.shape[0]
    for i in range(num_test_points):
        print('Q[{}]: I = {}, sigmaI = {}'.format(i, binned_q.i[i], binned_q.sigma[i]))
        # assert pytest.approx(binned_q.i[i], golden_log_weighted_iq[i], 1E-5)
        # assert pytest.approx(binned_q.sigma[i], golden_log_weighted_sigma[i], 1E-5)
    # END-FOR

    return


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
def skip_test_binning_1d_workflow(generic_IDF):
    """ Test the workflow to bin I(Q) in 1D from momentum transfer calculation to calculated binned
    Q, dQ, I(Q), sigma_I(Q)
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


def assign_bins(bin_edges, data_points, det_counts):
    """
    Check the value of each data points and assign them with proper bin indexes
    Parameters
    ----------
    bin_edges
    data_points
    det_counts: ndarray
        detector counts
    Returns
    -------

    """
    import bisect

    bin_index_list = [-1] * data_points.shape[0]

    for i in range(data_points.shape[0]):
        bin_index = bisect.bisect_left(bin_edges, data_points[i])
        if data_points[i] < bin_edges[bin_index]:
            bin_index -= 1
            if bin_index < 0:
                raise NotImplementedError('Implementation error')
        bin_index_list[i] = bin_index
        print('{}, {}, {}, {}, {}, {}'.format(i, det_counts[i], data_points[i], bin_index, bin_edges[bin_index], bin_edges[bin_index+1]))
    # END-FOR

    # count
    counts = 0
    for i in range(bin_edges.shape[0] - 1):
        print('{}-th bin: count = {}'.format(i, bin_index_list.count(i)))
        counts +=  bin_index_list.count(i)
    print('sum = {}'.format(counts))
    for i in range(bin_edges.shape[0] - 1):
        print('{}, {}'.format(i, bin_index_list.count(i)))

    for i in range(data_points.shape[0]):
        print('{}'.format(det_counts[i]))
    # END-FOR
    return


if __name__ == "__main__":
    pytest.main()
