# Main method in this module implement step 2 of
# wavelength dependent inelastic incoherent scattering correction
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/-/issues/689

from drtsans.dataobjects import IQmod, IQazimuthal
from drtsans.dataobjects import verify_same_q1d
import numpy as np


def normalize_by_elastic_reference(i_of_q, ref_i_of_q):

    # check i_of_q and ref_i_of_q shall have same binning
    if not verify_same_q1d(i_of_q, ref_i_of_q):
        raise RuntimeError('blabla')

    # Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process
    q_min, q_max = determine_q_range(i_of_q)

    # Find scale factor K(lambda) , that minimizes sum_q |refI(q, lambda_ref) - K(lambda) refI(q, lambda)|^2
    scale_factor_k_vec = calculate_scale_factor(ref_i_of_q, q_min, q_max)

    # Normalize input I(q, lambda)
    i_of_q = normalize_intensity(i_of_q, scale_factor_k_vec)

    return i_of_q, scale_factor_k_vec, sigma_scale_factor_k_vec


def determine_q_range(iqmod):
    # Determine q_min and q_max  that exist in all I(q, lambda) for the fitting (minimization) process
    return 0., 0.


def calculate_scale_factor(i_of_q, q_min, q_max):
    # i_of_q: reference I(Q)

    # init structure
    scale_factor_k_vec = np.ndarray(shape=(n_lambda,), dtype='float')

    for i_lambda in range(n_lambda):
        # get wavelength
        lambda_i = wavelength(i_of_q, i_lambda)

        # calculate P and S
        # P(lambda_i) = sum_{k=1,N} I_k^{lambda_ref} * I_k^{lambda_i}
        p_lambda_i = 0.
        s_lambda_i = 0.
        for k in range(index_qmin_lambda_i, index_qmax_lambda_i):
            intensity_k_ref = intensity(i_of_q, k, i_ref_lambda)
            intensity_k_i = intensity(i_of_q, k, i_lambda_i)
            p_lambda_i += intensity_k_ref * intensity_k_i
            s_lambda_i += intensity_k_i * intensity_k_i

        # finally
        k_lambda_i = p_lambda_i / s_lambda_i
    # END-FOR

    return scale_factor_k_vec


def normalize_intensity(i_of_q, scale_factor_vec):
    return i_of_q








