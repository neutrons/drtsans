import numpy as np
from drtsans.dataobjects import IQmod, IQazimuthal
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/iq.py
from drtsans.iq import determine_1d_log_bins, BinningMethod, bin_intensity_into_q1d,\
    select_i_of_q_by_wedge
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/tests/unit/new/drtsans/i_of_q_binning_tests_data.py
from tests.unit.new.drtsans.i_of_q_binning_tests_data import generate_test_data
import pytest

# This module supports testing data for issue #239.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/247

# DEV - Wenduo Zhou <zhouw@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data


def test_1d_bin_log_wedge_no_wt():
    """Test the methods to select I(Qx, Qy) by wedge angles and do the binning

    The test data comes from example in '1D_bin_log_wedge_no_sub_no_wt' from eqsans_tof_q_binning_tests_R5.xlsx.
    File location: https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/uploads/5423db9b77dfd4911bf799c247530865/
                   eqsans_tof_q_binning_tests_R5.xlsx

    Returns
    -------

    """
    """Test '1D_bin_log_wedge_no_sub_no_wt'
    """
    # Define Q range
    q_min = 0.001  # Edge
    q_max = 0.010  # Edge
    step_per_decade = 10  # 10 steps per decade

    min_wedge_angle = -45.
    max_wedge_angle = 45

    # Get data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)
    log_bins = determine_1d_log_bins(q_min, q_max, step_per_decade)

    # Test high level method
    # Define input data
    test_i_q = IQazimuthal(intensity=intensities, error=sigmas, qx=qx_array, qy=qy_array,
                           delta_qx=dqx_array, delta_qy=dqy_array)
    # Select I(Q) inside wedge
    wedge_i_of_q = select_i_of_q_by_wedge(test_i_q, min_wedge_angle, max_wedge_angle)

    # Convert from Q2D to Q1d because the test case expects to I(Q)
    # Re-construct a IQmod  by
    # (1) Q = sqrt(Qx**2 + Qy**2)
    # (2) sigmaQ = sqrt(sigmaQx**2 + sigmaQy**2)
    test_i_q1d = IQmod(intensity=wedge_i_of_q.intensity,
                       error=wedge_i_of_q.error,
                       mod_q=np.sqrt(wedge_i_of_q.qx**2 + wedge_i_of_q.qy**2),
                       delta_mod_q=np.sqrt(wedge_i_of_q.delta_qx**2 + wedge_i_of_q.delta_qy**2))

    binned_iq = bin_intensity_into_q1d(test_i_q1d, log_bins, bin_method=BinningMethod.NOWEIGHT)

    # Verification
    assert binned_iq.intensity[7] == pytest.approx(67.7, abs=1E-10)
    assert binned_iq.error[7] == pytest.approx(2.601922, abs=1E-5)
    assert binned_iq.delta_mod_q[7] == pytest.approx(0.0115, abs=1E-4), 'Q resolution (Q[7] = {}) is ' \
                                                                        'incorrect comparing to {}.' \
                                                                        ''.format(binned_iq.delta_mod_q[7],
                                                                                  0.0115)


if __name__ == '__main__':
    pytest.main([__file__])
