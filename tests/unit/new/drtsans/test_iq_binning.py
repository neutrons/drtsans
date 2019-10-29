import numpy as np
import pytest
from drtsans.iq import bin_iq_into_linear_q1d, bin_iq_into_logarithm_q1d, BinningMethod,\
    determine_1d_linear_bins, determine_1d_log_bins, do_1d_weighted_binning, do_1d_no_weight_binning
import bisect

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
# linear_bin_centers = np.ndarray([0.000314376, 0.000943129, 0.001571882, 0.002200635, 0.002829388, 0.003458141,
#                                  0.004086894, 0.004715647, 0.005344399, 0.005973152], dtype=np.float)
# linear_bin_right_bound = np.array([0.000628753, 0.001257506, 0.001886259, 0.002515011, 0.003143764, 0.003772517,
#                                    0.00440127, 0.005030023, 0.005658776, 0.006287529], dtype=np.float)
# linear_bin_left_bound = np.array([0, 0.000628753, 0.001257506, 0.001886259, 0.002515011, 0.003143764, 0.003772517,
#                                   0.00440127 , 0.005030023, 0.005658776], dtype=np.float)

# Linear-weighted
gold_linear_bin_iq = np.array([0, 10, 7.460646, 575.3626, 242.9668, 4.561801, 4.243206, 1.023098, 0.75, 1.],
                              dtype=np.float)
gold_linear_bin_sigmaq = np.array([1, 1.58113883, 1.115096873, 7.585265905, 4.929166353, 0.57082652, 0.571314703,
                                   0.270330222, 0.433012702, 0.5], dtype=np.float)

# Linear-no weight TODO FIXME - gold data shall be filled
gold_linear_bin_no_weight_iq = np.array([],
                                        dtype=np.float)  # binned I(Q)
gold_linear_bin_no_weight_sq = np.array([],
                                        dtype=np.float)  # binned sigma I(Q)

# Logarithm binning
log_bin_centers = np.array([6.58E-04, 7.05E-04, 7.56E-04, 8.11E-04, 8.70E-04, 9.33E-04, 1.00E-03, 1.07E-03,
                            1.15E-03,
                            1.23E-03, 1.32E-03, 1.42E-03, 1.52E-03, 1.63E-03, 1.75E-03, 1.87E-03, 2.01E-03,
                            2.15E-03,
                            2.31E-03, 2.48E-03, 2.66E-03, 2.85E-03, 3.05E-03, 3.27E-03, 3.51E-03, 3.76E-03,
                            4.04E-03, 4.33E-03,
                            4.64E-03, 4.98E-03, 5.34E-03, 5.72E-03, 6.14E-03], dtype=np.float)

# Log-weighted: NOTE: not the full set but first a few with counts
gold_log_bin_weighted_iq = np.array([0, 0, 0, 0, 10, 0, 0, 10, 10, 0, 0, 5, 0, 0, 473.9472574, 0, 495.9012097,
                                     505.9288538, 660.20558, 567.5506439, 121, 262.7392886, 0, 0, 12.16131385,
                                     2.472004695, 1.989894845, 126.6649782, 1.44, 0.931190118, 0.75, 1.487164355,
                                     1],
                                    dtype=np.float)  # binned I(Q)
gold_log_bin_weighted_sq = np.array([1, 1, 1, 1, 2.236067977, 1, 1, 3.16227766, 2.236067977, 1, 1, 1.118033989, 1, 1,
                                     15.39394779, 1, 15.74644737, 15.90485545, 18.16873111, 10.65411323, 11,
                                     5.730829877, 1, 1, 0.967205087, 0.703136501, 0.705318163, 4.253821445,
                                     0.489897949, 0.364728885, 0.433012702, 0.704074891, 0.707106781],
                                    dtype=np.float)  # binned sigma I(Q)

# Log-no-weight: NOTE: no the full set but first a few with counts
gold_log_bin_no_weight_iq = np.array([],
                                     dtype=np.float)  # binned I(Q)
gold_log_bin_no_weight_sq = np.array([],
                                     dtype=np.float)  # binned Sigma(Q)


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

    # Correct zero counts uncertainties (Sigma(Q)) to 1
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
    """ Test binning I(Q) with linear binning
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

    # Calculate bin centers and bin edges by Qmin, Qmax and number of bins
    bin_centers, bin_edges = determine_1d_linear_bins(q_min, q_max, bins)

    # Test
    assert bin_centers.shape == (10, ), 'Bin shape incorrect'
    assert bin_edges[0] == 0., 'Bin edge shall start from 0'
    assert abs(bin_centers[9] - 0.0059731523) < 0.00005, 'Last bin center is too far away from gold data'

    # Get assignment of each Qi to Qk
    bins_dict = assign_bins(bin_edges, q_array, iq_array, bin_centers)
    i_weighted_array, sigma_i_weighted_array = bin_weighted_prototype(bin_centers, iq_array, bins_dict)

    # Do weighted binning
    binned_q = do_1d_weighted_binning(q_array, dq_array, iq_array, sigma_q_array, bin_centers, bin_edges)

    # Test for Q bins
    assert binned_q.q.shape == (10, )
    assert pytest.approx(binned_q.q[0], q_max/bins * 0.5, 1E-5)

    # Test for I(Q)
    assert np.allclose(i_weighted_array, binned_q.i, 1e-6)
    assert np.allclose(sigma_i_weighted_array, binned_q.sigma, 1e-6)
    assert np.sqrt(np.sum((i_weighted_array - binned_q.i)**2)) < 1e-6

    # Test no-weight binning
    i_noweight_array, sigma_i_noweight_array = do_no_weight_prototype(bin_centers, iq_array, bins_dict)
    no_weight_iq = do_1d_no_weight_binning(q_array, dq_array, iq_array, sigma_q_array, bin_centers, bin_edges)

    # verify
    assert np.allclose(i_noweight_array, no_weight_iq.i, 1e-6)
    assert np.allclose(sigma_i_noweight_array, no_weight_iq.sigma, 1e-6)
    assert np.sqrt(np.sum((i_noweight_array - no_weight_iq.i)**2)) < 1e-6

    # Test to go through wrapper method
    wiq = bin_iq_into_linear_q1d(iq_array, sigma_q_array, q_array, dq_array, bins, q_min, q_max,
                                 BinningMethod.WEIGHTED)
    assert np.allclose(wiq.i, binned_q.i, 1e-6)

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
    q_min = q_array.min()
    q_max = 1.
    step_per_decade = 33
    bin_centers, bin_edges = determine_1d_log_bins(q_min, q_max, step_per_decade)

    for ibin in range(bin_centers.shape[0]):
        print('{}\t{}\t{}\t{}\t'.format(ibin, bin_centers[ibin], bin_edges[ibin],
                                        bin_edges[ibin + 1]))

    # Verify: bin size, min and max are all on the power of 10
    # assert bin_edges.shape[0] == bin_centers.shape[0] + 1
    # assert bin_centers.shape[0] == 100
    # assert abs(bin_centers[0] - q_min) < 1.E-12
    # assert abs(bin_centers[99] - q_max) < 1.E-12

    # Test bins
    bin_assignment_dict = assign_bins(bin_edges, q_array, iq_array, bin_centers)

    # Bin with weighted binning algorithm
    binned_q = do_1d_weighted_binning(q_array, dq_array, iq_array, sigma_q_array, bin_centers, bin_edges)
    # do the weighted binning
    weighted_i_q, weighted_i_sigma_q = bin_weighted_prototype(bin_centers, iq_array, bin_assignment_dict)

    # verify: NaN first
    assert np.allclose(np.isnan(weighted_i_q), np.isnan(binned_q.i)), 'NaN shall be same'

    # verify
    if not np.allclose(weighted_i_q, binned_q.i, 1e-12, True):
        print('Different of log binning I(Q)')
        for k in range(binned_q.i.shape[0]):
            if np.isnan(weighted_i_q[k]) and np.isnan(binned_q.i[k]):
                # print('{}\t{}\t{}\t{}'.format(k, 'NaN', 'NaN', 0))
                pass
            elif not np.isnan(weighted_i_q[k]) and not np.isnan(binned_q.i[k]):
                diff_k = abs(weighted_i_q[k] - binned_q.i[k])
                if abs(diff_k) < 1e-10:
                    print('{}\t{}\t{}\t{}'.format(k, weighted_i_q[k], binned_q.i[k],
                                                  weighted_i_q[k] - binned_q.i[k]))
                else:
                    print('{}\t{}\t{}\t{}\tError'.format(k, weighted_i_q[k], binned_q.i[k],
                                                         weighted_i_q[k] - binned_q.i[k]))
            else:
                pass
                # print('{}\t{}\t{}\t{}\tError'.format(k, weighted_i_q[k], binned_q.i[k],
                #                                      weighted_i_q[k] - binned_q.i[k]))
    # END-FOR

    assert np.allclose(weighted_i_q[~np.isnan(weighted_i_q)],
                       binned_q.i[~np.isnan(binned_q.i)],
                       1e-10, True), 'Binned I(Q) does not match'
    assert np.allclose(weighted_i_sigma_q, binned_q.sigma, 1e-12, True), 'Binned sigma_I(Q) does not match'

    # Test no-weight binning
    no_weight_binned_iq = do_1d_no_weight_binning(q_array, dq_array, iq_array,
                                                  sigma_q_array, bin_centers, bin_edges)
    # do no weight binning
    noweight_i_q, noweight_sigma_q = do_no_weight_prototype(bin_centers, iq_array, bin_assignment_dict)
    # verify: NaN first
    assert np.allclose(np.isnan(noweight_i_q), np.isnan(no_weight_binned_iq.i)), 'No-weight binning NaN shall be same'
    assert np.allclose(np.isnan(noweight_sigma_q), np.isnan(no_weight_binned_iq.sigma)), 'No-weight binning ' \
                                                                                         'Sigma(I) NaN shall be same'

    # verify: non-NaN value
    assert np.allclose(noweight_i_q[~np.isnan(noweight_i_q)],
                       no_weight_binned_iq.i[~np.isnan(no_weight_binned_iq.i)],
                       1e-10, True), 'No-weight binned I(Q) does not match'
    assert np.allclose(noweight_sigma_q[~np.isnan(noweight_sigma_q)],
                       no_weight_binned_iq.sigma[~np.isnan(no_weight_binned_iq.sigma)],
                       1e-10, True), 'No-weight binned sigma_I(Q) does not match'

    # Test to go through wrapper method
    wiq = bin_iq_into_logarithm_q1d(iq_array, sigma_q_array, q_array, dq_array, step_per_decade,
                                    q_min, q_max, BinningMethod.WEIGHTED)
    assert wiq

    # Note: disable the check due to a different algorithm (Lisa vs William) to generate log bins
    # assert np.allclose(wiq.i[~np.isnan(wiq.i)],
    #                    binned_q.i[~np.isnan(binned_q.i)],
    #                    1e-10, True)

    return


def assign_bins(bin_edges, data_points, det_counts, bin_centers):
    """Check the value of each data points and assign them with proper bin indexes

    Parameters
    ----------
    bin_edges
    data_points
    det_counts: ndarray
        detector counts
    bin_centers: ndarray
        bin centers
    Returns
    -------
    dict
        dictionary of Q[k] to detector Q[i]'s list
    """
    # Print bin edges
    print('Bin edges:')
    for i in range(bin_centers.shape[0]):
        print('{}\t{}\t{}'.format(i, bin_edges[i], bin_edges[i + 1]))
    # END-FOR

    # Register binning
    bin_index_list = [-1] * data_points.shape[0]
    bins_dict = dict()  # record [Qi] associated with each Qk
    i_raw_dict = dict()
    for i in range(bin_edges.shape[0] - 1):
        bins_dict[i] = set()
        i_raw_dict[i] = 0
    # END-FOR

    print('Assign bins:')
    for i in range(data_points.shape[0]):
        bin_index = bisect.bisect_left(bin_edges, data_points[i])
        if data_points[i] < bin_edges[bin_index]:
            bin_index -= 1
            if bin_index < 0:
                raise NotImplementedError('Implementation error')
        bin_index_list[i] = bin_index
        print('{}\t{}\t{}\t{}\t({}, {})'
              ''.format(i, det_counts[i], data_points[i], bin_index, bin_edges[bin_index], bin_edges[bin_index+1]))
        # register
        bins_dict[bin_index].add(i)
        if det_counts[i] > 0:
            i_raw_dict[bin_index] += 1
    # END-FOR

    # Binning assignment demo 2
    for k in range(bin_centers.shape[0]):
        print('{}\t{}'.format(k, sorted(list(bins_dict[k]))))
    # END-FOR(k)

    return bins_dict


def bin_weighted_prototype(bin_centers, det_counts, bins_dict):
    """Do weighted binning

    Parameters
    ----------
    bin_centers
    det_counts
    bins_dict

    Returns
    -------
    ndarray, ndarray
    """

    # Do weighted binning
    k_for_print = [7]
    i_raw_array = np.zeros_like(bin_centers, dtype=float)
    w_array = np.zeros_like(bin_centers, dtype=float)

    for k in range(bin_centers.shape[0]):
        i_k_raw = 0.
        w_k = 0.
        for i_det in bins_dict[k]:
            if np.isnan(det_counts[i_det]):
                # ignore NaN
                diff_i_k_raw = 0
                diff_w_k = 0
            elif det_counts[i_det] < 1E-20:
                # zero
                diff_i_k_raw = 0
                diff_w_k = 1
            else:
                # non-zero
                diff_i_k_raw = 1
                diff_w_k = 1./det_counts[i_det]
            # END-IF-ELSE
            i_k_raw += diff_i_k_raw
            w_k += diff_w_k

            if k in k_for_print:
                print('k = {}: i = {}, I(Q) = {}, I_raw_k += {}, w_k += {}'
                      ''.format(k, i_det, det_counts[i_det], diff_i_k_raw, diff_w_k))
            # END-IF
        # END-FOR

        # register
        i_raw_array[k] = i_k_raw
        w_array[k] = w_k
    # END-FOR

    # Final touch for weighted binning
    i_weighted_array = i_raw_array / w_array
    sigma_i_weighed_array = 1 / np.sqrt(w_array)

    # Do the summation
    print('Weighted')
    for i in range(bin_centers.shape[0]):
        print('{}\t{}\t{}\t{}\t{}'
              ''.format(i, bin_centers[i], i_raw_array[i], i_weighted_array[i], sigma_i_weighed_array[i]))

    return i_weighted_array, sigma_i_weighed_array


def do_no_weight_prototype(bin_centers, det_counts, bins_dict):
    """Do weighted binning

    Parameters
    ----------
    bin_centers
    det_counts
    bins_dict

    Returns
    -------
    ndarray, ndarray
    """
    # Do weighted binning
    k_for_print = [7]
    i_sum_array = np.zeros_like(bin_centers, dtype=float)
    i_q_array = np.zeros_like(bin_centers, dtype=float)
    sigma_iq_array = np.zeros_like(bin_centers, dtype=float)

    for k in range(bin_centers.shape[0]):
        i_k_raw = 0.
        sigma_sq_k = 0.
        num_counts = 0
        for i_det in bins_dict[k]:
            if np.isnan(det_counts[i_det]):
                # ignore NaN
                diff_i_k_raw = 0
                diff_sigma2 = 0
            elif det_counts[i_det] < 1E-20:
                # zero
                diff_i_k_raw = 0
                diff_sigma2 = 1
                num_counts += 1
            else:
                # non-zero
                diff_i_k_raw = det_counts[i_det]
                diff_sigma2 = det_counts[i_det]
                num_counts += 1
            # END-IF-ELSE
            i_k_raw += diff_i_k_raw
            sigma_sq_k += diff_sigma2

            if k in k_for_print:
                print('k = {}: i = {}, I(Q) = {}, I_sum_k += {}, sigma^2_k += {}'
                      ''.format(k, i_det, det_counts[i_det], diff_i_k_raw, diff_sigma2))
            # END-IF
        # END-FOR

        # register
        i_sum_array[k] = i_k_raw
        if num_counts == 0:
            # zero counts: NaN
            i_q_array[k] = np.nan
            sigma_iq_array[k] = np.nan
        else:
            # non-zero counts
            i_q_array[k] = i_k_raw / num_counts
            sigma_iq_array[k] = np.sqrt(sigma_sq_k) / num_counts
    # END-FOR

    # Output
    print('No-Weight')
    for i in range(bin_centers.shape[0]):
        print('{}\t{}\t{}\t{}\t{}'
              ''.format(i, bin_centers[i], i_sum_array[i], i_q_array[i], sigma_iq_array[i]))

    return i_q_array, sigma_iq_array


if __name__ == "__main__":
    pytest.main([__file__])
