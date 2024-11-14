import pytest

from drtsans.tof.eqsans.correction_api import listify_incohfit_parameter


@pytest.mark.parametrize(
    "parameter, expected",
    [
        (10, [10.0, 10.0]),
        (None, [None, None]),
        (10, [10.0, 10.0]),
        (10, [10, 10]),
        ([1.0], [1.0, 1.0]),
        ([2.0, 3.0], [2.0, 3.0]),
        ([4.0, 5], [4.0, 5.0]),
        (9.0, [9.0, 9.0]),
        ([3, 4], [3, 4]),
        ([6], [6, 6]),
        (7, [7, 7]),
        ([True], [True, True]),
        ([False, True], [False, True]),
        ([True, False], [True, False]),
        ([False], [False, False]),
        (True, [True, True]),
        (False, [False, False]),
    ],
)
def test_listify_incohfit_parameter(parameter, expected):
    if expected == "error":
        with pytest.raises(ValueError):
            listify_incohfit_parameter(parameter)
    else:
        assert listify_incohfit_parameter(parameter) == expected
