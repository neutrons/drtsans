import numpy as np
# from drtsans.iq import IofQ


# This test implements issue #169 to verify
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/tree/169_bin_q1d
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>


def generate_test_data(dimension):
    """Generate test data

    Test data including I(Q), Q and dQ depending on 1D or 2D

    Parameters
    ----------
    dimension : int
        1D data or 2D

    Returns
    -------

    """
    if dimension == 1:
        pass
    elif dimension == 2:
        pass
    else:
        raise RuntimeError('Dimension equal to {} is not supported'.format(dimension))

    return


def test_1d_bin_linear_no_wt():
    """Test '1D_bin_linear_no_sub_wt'

    Returns
    -------

    """
    q_min = 0.
    q_max = 0.010
    num_bins = 10

    # Calculate and verify Q bins
    binned_data = IofQ

    # Calculate and verify
    # I(0.0035) =		68.92857
    assert binned_data.i[3] == 68.92857
    # di(0.0035)		2.218889
    assert binned_data.sigma[3] == 2.218889
    # sigma_Q(0.0035) = 		3.722E-05
    assert binned_data.dq[3] == 3.722E-05

    return


def test_1d_bin_log_no_wt():
    """Test '1D_bin_log_no_sub_no_wt'

    Returns
    -------

    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    num_bins = 10  # 10 steps per decade

    # Verify bin edges and bin center

    # Binned I(Q) no-weight
    binned_data = IofQ

    # Verify: 2 I(Q) in bin: Q(3, 2, 3.1), Q(3, 2, 3.2)
    # I(0.0022) = 70.00000
    assert 1 ==1
    # dI(0.0022) = 5.9160797831
    assert 1 == 1
    # sigma_Q(0.0022) = 2.529E-05
    assert 1 == 1

    return


def test_2d_bin_no_wt():
    """Test '2D_bin_no_sub_no_wt'

    Returns
    -------

    """

    return


def determine_2d_bins():
    """

    Returns
    -------

    """
    print('Hello World')

    min_q = -0.006134
    min_q_near = - 0.003254

    qx_min_edge = min_q - 0.5 * (min_q_near - min_q)
    qx_max_edge = 0.0068250907

    qx_step = (qx_max_edge - qx_min_edge) / 5

    print(qx_min_edge)
    print(qx_step)
    print(qx_step - 0.002880)

    return

determine_2d_bins()


