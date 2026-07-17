import numpy as np
import pytest

from drtsans.dataobjects import IQmod
from drtsans.tof.eqsans.correction_api import (
    bypass_correction_for_single_wavelength_bin,
    CorrectionConfiguration,
    listify_incohfit_parameter,
)


def _make_iq1d(wavelengths):
    n = len(wavelengths)
    return IQmod(
        intensity=np.ones(n),
        error=np.ones(n),
        mod_q=np.linspace(0.01, 0.1, n),
        delta_mod_q=None,
        wavelength=np.array(wavelengths),
    )


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


def test_bypass_correction_for_single_wavelength_bin_both_requested(mocker):
    mock_logger = mocker.patch("drtsans.tof.eqsans.correction_api.logger")
    correction_setup = CorrectionConfiguration(
        do_elastic_correction=True,
        do_inelastic_correction=[True, False],
    )
    iq1d = _make_iq1d([3.0])

    bypass_correction_for_single_wavelength_bin(iq1d, correction_setup, frameskip_frame=0)

    assert correction_setup.do_elastic_correction is False
    assert correction_setup.do_inelastic_correction == [False, False]
    assert mock_logger.warning.call_count == 2


def test_bypass_correction_for_single_wavelength_bin_only_inelastic_requested(mocker):
    mock_logger = mocker.patch("drtsans.tof.eqsans.correction_api.logger")
    correction_setup = CorrectionConfiguration(
        do_elastic_correction=False,
        do_inelastic_correction=[True, True],
    )
    iq1d = _make_iq1d([3.0, 3.0, 3.0])  # single unique wavelength, repeated across pixels

    bypass_correction_for_single_wavelength_bin(iq1d, correction_setup, frameskip_frame=1)

    assert correction_setup.do_elastic_correction is False
    assert correction_setup.do_inelastic_correction == [True, False]
    mock_logger.warning.assert_called_once()


def test_bypass_correction_for_single_wavelength_bin_multiple_bins_no_change(mocker):
    mock_logger = mocker.patch("drtsans.tof.eqsans.correction_api.logger")
    correction_setup = CorrectionConfiguration(
        do_elastic_correction=True,
        do_inelastic_correction=[True, True],
    )
    iq1d = _make_iq1d([1.0, 2.0, 3.0])

    bypass_correction_for_single_wavelength_bin(iq1d, correction_setup, frameskip_frame=0)

    assert correction_setup.do_elastic_correction is True
    assert correction_setup.do_inelastic_correction == [True, True]
    mock_logger.warning.assert_not_called()
