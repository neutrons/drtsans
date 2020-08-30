# Test drtsans.auto_wedge
import pytest
import numpy as np
from drtsans.auto_wedge import _create_fit_function, _set_function_param_value, _calculate_function
from mantid.simpleapi import FlatBackground, Gaussian


def test_create_functions():
    """Test create a set of functions from a string
    """
    # set input
    input_function_str = 'name=FlatBackground,A0=26452.031557632155;name=Gaussian,Height=383890.0871952363,' \
                         'PeakCentre=106.28166682679631,Sigma=11.468820017817437;name=Gaussian,' \
                         'Height=685867.5517201225,PeakCentre=276.6831781469825,Sigma=27.242175478846256'

    # execute
    test_functions = _create_fit_function(input_function_str)

    # verify
    assert len(test_functions) == 3, f'There shall be 3 functions in the list but not {len(test_functions)}'
    # background
    background = test_functions[0]
    assert background.name == 'FlatBackground'
    assert background.A0 == pytest.approx(26452.031557632155, 1E-12)
    # gaussian 1
    g1 = test_functions[1]
    assert g1.name == 'Gaussian'
    assert g1.Sigma == pytest.approx(11.468820017817437, 1E-12)
    # gaussian 2
    g2 = test_functions[2]
    assert g2.name == 'Gaussian'
    assert g2.Height == pytest.approx(685867.5517201225, 1E-9)


def test_set_param_value():
    """Test setting parameter values
    """
    # Create input functions
    test_funcs = [FlatBackground(), Gaussian(), Gaussian]

    # Input dictionary
    param_dict = {'f0.A0': 3.2,
                  'f1.PeakCentre': 1.0,
                  'f1.Sigma': 0.42,
                  'f1.Height': 10.11,
                  'f2.Sigma': 1.24,
                  'f2.Height': 5.11,
                  'f2.PeakCentre': 7.2}

    # Execute
    for param_name, param_value in param_dict.items():
        _set_function_param_value(test_funcs, param_name, param_value)

    # Verify
    assert test_funcs[0].A0 == pytest.approx(3.2, 0.05)
    assert test_funcs[1].PeakCentre == pytest.approx(1.0, 0.05)
    assert test_funcs[1].Sigma == pytest.approx(0.42, 0.005)
    assert test_funcs[2].Sigma == pytest.approx(1.24, 0.005)
    assert test_funcs[2].Height == pytest.approx(5.11, 0.005)


def test_calculate_functions():
    """Test calculate a combo function
    """
    # Create test function
    functions = [None] * 3
    functions[0] = FlatBackground(A0=10.)
    functions[1] = Gaussian(PeakCentre=5., Height=40., Sigma=1.)
    functions[2] = Gaussian(PeakCentre=95., Height=80., Sigma=2.)

    # Set up X
    vec_x = np.arange(100)

    # Calculate
    vec_y = _calculate_function(functions, vec_x)

    # Verify
    assert vec_x[5] == pytest.approx(5., 1E-6)
    assert vec_y[5] == pytest.approx(50., 1E-6)

    assert vec_x[95] == pytest.approx(95., 1E-6)
    assert vec_y[95] == pytest.approx(90., 1E-6)


if __name__ == '__main__':
    pytest.main(__file__)
