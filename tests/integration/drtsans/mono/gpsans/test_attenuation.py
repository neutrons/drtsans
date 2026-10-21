from datetime import date
import os

import pytest
from mantid.simpleapi import LoadEventNexus

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/mono/gpsans/attenuation.py
from drtsans.mono.gpsans import attenuation_factor
from drtsans.mono.gpsans.attenuation import (
    _attenuator_name,
    _load_default_attenuation_coefficients,
    _run_start_time,
)
from drtsans.samplelogs import SampleLogs


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "run_number, effective_date, expected_value, expected_error",
    [
        # empty cell, 2026-06-07, before the 2026-09-24 calibration
        (170056, date(2020, 2, 26), 0.00021825564301391177, 2.682365493341721e-05),
        # empty cell, 2026-10-18, after the 2026-09-24 calibration
        (176202, date(2026, 9, 24), 0.00018121045456012452, 2.7952234707187197e-05),
    ],
)
def test_attenuation_factor_effective_block(
    datarepo_dir, clean_workspace, run_number, effective_date, expected_value, expected_error
):
    """Attenuation factor of the x10k attenuator at 4.76 Å, from the coefficients block effective for the run.

    Expected values are A * exp(-B * wavelength) + C, with uncertainty propagated from the independent
    uncertainties of A, B and C, evaluated at the logged wavelength with the coefficients of the effective block.
    """
    filename = os.path.join(datarepo_dir.gpsans, f"CG2_{run_number}.nxs.h5")
    workspace = LoadEventNexus(Filename=filename, MetaDataOnly=True, OutputWorkspace=f"CG2_{run_number}")
    clean_workspace(workspace)

    assert _attenuator_name(workspace) == "x10k"
    assert SampleLogs(workspace).single_value("wavelength") == pytest.approx(4.76, abs=0.01)
    assert _load_default_attenuation_coefficients(_run_start_time(workspace)).effective_date == effective_date

    value, error = attenuation_factor(workspace)
    assert value == pytest.approx(expected_value, rel=1e-6)
    assert error == pytest.approx(expected_error, rel=1e-6)
