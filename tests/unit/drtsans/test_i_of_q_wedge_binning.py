import pytest

from drtsans.dataobjects import IQazimuthal, q_azimuthal_to_q_modulo
from drtsans.iq import (
    BinningMethod,
    bin_intensity_into_q1d,
    determine_1d_log_bins,
    select_i_of_q_by_wedge,
)
from tests.unit.drtsans.i_of_q_binning_tests_data import generate_test_data

# This module supports testing data for issue #239.
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/247

# DEV - Wenduo Zhou <zhouw@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>

# All tests data are generated in tests.unit.new.drtsans.i_of_q_binning_tests_data


def test_1d_bin_log_wedge_no_wt():
    """Test the methods to select I(Qx, Qy) by wedge angles and do the binning"""
    # Define Q range
    q_min = 0.001  # center min
    q_max = 0.007890  # center max
    step_per_decade = 10  # 10 steps per decade

    min_wedge_angle = -45.0
    max_wedge_angle = 45

    # Get data
    intensities, sigmas, qx_array, dqx_array, qy_array, dqy_array = generate_test_data(2, True)
    log_bins = determine_1d_log_bins(q_min, q_max, True, step_per_decade)

    # Test high level method
    # Define input data
    test_i_q = IQazimuthal(
        intensity=intensities,
        error=sigmas,
        qx=qx_array,
        qy=qy_array,
        delta_qx=dqx_array,
        delta_qy=dqy_array,
    )
    # Select I(Q) inside wedge
    wedge_i_of_q = select_i_of_q_by_wedge(test_i_q, min_wedge_angle, max_wedge_angle)

    # Convert from Q2D to Q1d because the test case expects to I(Q)
    test_i_q1d = q_azimuthal_to_q_modulo(wedge_i_of_q)

    binned_iq = bin_intensity_into_q1d(test_i_q1d, log_bins, bin_method=BinningMethod.NOWEIGHT)

    # Verification
    # bin index = 7, bin center = 0.005586, bin edges = (0.00631, 0.007033)
    assert binned_iq.intensity[7] == pytest.approx(70.0, abs=1e-10)
    assert binned_iq.error[7] == pytest.approx(3.741657387, abs=1e-5)
    assert binned_iq.delta_mod_q[7] == pytest.approx(0.011460, abs=1e-4), (
        "Q resolution (Q[7] = {}) is incorrect comparing to {}.".format(binned_iq.delta_mod_q[7], 0.0115)
    )


if __name__ == "__main__":
    pytest.main([__file__])
