import numpy as np
import pytest


# 2D "detector view" I(i, j) array
# After flattened, it goes from top left to top right, and then move down to next row
i_det_array = np.array([
    [1, 2, 3, 4, 4, 5, 6, 7, 8],
    [1, 0, 1, 1, 1, 0, 1, 1, 1],
    [1, 87, 130, 208, 121, 107, 279, 186, 0],
    [1, 172, 250, 500, 479, 512, 499, 208, 1],
    [2, 147, 600, 5, 10, 5, 700, 189, 3],
    [1, 239, 562, 10, 0, 10, 800, 217, 1],
    [3, np.nan, 550, 5, 10, 5, 689, 228, 1],
    [1, 178, 567, 503, 469, 489, 499, 156, 0],
    [0, 108, 350, 389, 409, 253, 192, 209, 1],
    [1, 80, 193, 120, 148, 108, 201, 210, 1],
    [1, 3, 2, 0, 1, 1, 0, 1, 1]
    ])


# 2D "detector view" Qx(i, j) array
# After flattened, it goes from top left to top right, and then move down to next row
qx_det_array = np.array([
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492],
        [-0.0092152492, -0.0069114662, -0.0046076581, -0.0023038332, 0, 0.0023038332, 0.0046076581, 0.0069114662,
         0.0092152492]])

# 2D "detector view" Qy(i, j) array
# After flattened, it goes from top left to top right, and then move down to next row
qy_det_array = np.array([
    [-0.006869579, -0.006869579, -0.006869579, -0.006869579, -0.006869579, -0.006869579, -0.006869579,
     -0.006869579, -0.006869579],
    [-0.0051521964, -0.0051521964, -0.0051521964, -0.0051521964, -0.0051521964, -0.0051521964, -0.0051521964,
     -0.0051521964, -0.0051521964],
    [-0.0034348033, -0.0034348033, -0.0034348033, -0.0034348033, -0.0034348033, -0.0034348033, -0.0034348033,
     -0.0034348033, -0.0034348033],
    [-0.0017174034, -0.0017174034, -0.0017174034, -0.0017174034, -0.0017174034, -0.0017174034, -0.0017174034,
     -0.0017174034, -0.0017174034],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0.0017174034, 0.0017174034, 0.0017174034, 0.0017174034, 0.0017174034, 0.0017174034, 0.0017174034,
     0.0017174034, 0.0017174034],
    [0.0034348033, 0.0034348033, 0.0034348033, 0.0034348033, 0.0034348033, 0.0034348033, 0.0034348033,
     0.0034348033, 0.0034348033],
    [0.0051521964, 0.0051521964, 0.0051521964, 0.0051521964, 0.0051521964, 0.0051521964, 0.0051521964,
     0.0051521964, 0.0051521964],
    [0.006869579, 0.006869579, 0.006869579, 0.006869579, 0.006869579, 0.006869579, 0.006869579,
     0.006869579, 0.006869579],
    [0.0085869477, 0.0085869477, 0.0085869477, 0.0085869477, 0.0085869477, 0.0085869477, 0.0085869477,
     0.0085869477, 0.0085869477]
    ])

i_final_det_array = np.array([
    [1, 2, 3, 4, 5, 6, 7, 8],
    [1, 0, 1, 1, 0, 1, 1, 1],
    [1, 87, 130, 152.9969604863, 107, 279, 186, 0],
    [1, 172, 250, 489.2747701736, 512, 499, 208, 1],
    [1.3333333333, 182.0362694301, 580.3786574871, 2.1428571429, 6.6666666667, 746.6666666667, 202.0344827586, 1.5],
    [3, np.nan, 550, 6.6666666667, 5, 689, 228, 1],
    [1, 178, 567, 485.4053497942, 489, 499, 156, 0],
    [0, 108, 350, 398.7493734336, 253, 192, 209, 1],
    [1, 5.7831325301, 3.958974359, 1.4887671845, 1.9816513761, 0.995049505, 1.990521327, 1]
    ])



def test_create_2d_bins():
    """Tester for method to generate 2D binning

    Returns
    -------

    """

    return


def assign_2d_bin_is(qx_array, qy_array, num_x_bins, num_y_bins):
    """Assign each pixel with a proper bin number

    Note: this method is exactly the same as Lisa's and is used to verify np.histogram2d()

    Note: 'is' in the method name stands for Instrument scientist (Lisa)

    Parameters
    ----------
    qx_array: ndarray (N, M)
        Qx
    qy_array: ndarray (N, M)
        Qy
    num_x_bins: integer
        number of bins on Qx array
    num_y_bins: integer
        number of bins on Qy array
    Returns
    -------

    """
    # Calculate delta(q_x) and delta(q_y)
    qx_bin_size = (np.max(qx_array) - np.min(qx_array)) / (num_x_bins - 1)
    qy_bin_size = (np.max(qy_array) - np.min(qy_array)) / (num_y_bins - 1)

    # Assign
    qx_index_array = np.ceil(qx_array / qx_bin_size + 0.5 * (1 + (np.max(qx_array) - np.min(qx_array)) / qx_bin_size))
    qy_index_array = np.ceil(qy_array / qy_bin_size + 0.5 * (1 + (np.max(qy_array) - np.min(qy_array)) / qy_bin_size))

    return qx_index_array, qy_index_array


def bin_iq_2d_is(i_array, qx_array, qy_array, num_x_bins, num_y_bins):
    """Bin I(Qx, Qy) for 2-dimensional Q

    Note: this method implementes Lisa's algorithm and is used to compare np.histogram2d()

    Parameters
    ----------
    i_array
    qx_array
    qy_array
    num_x_bins
    num_y_bins

    Returns
    -------

    """
    # Assign Qx and Qy
    qx_index_array, qy_index_array = assign_2d_bin_is(qx_array, qy_array, num_x_bins, num_y_bins)
    qx_index_array = qx_index_array.astype('int32')
    qy_index_array = qy_index_array.astype(int)

    # Bin (Qx', Qy')
    # TODO - continue from here
    raise NotImplementedError('ASAP as next: {} {} {}'.format(i_array, qx_index_array, qy_index_array))


if __name__ == "__main__":
    pytest.main()
