#!/usr/bin/env python
from datetime import date
import importlib.resources
import textwrap

import numpy as np
import pytest

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/mono/gpsans/attenuation.py
from drtsans.mono.gpsans import attenuation_factor
from drtsans.mono.gpsans.attenuation import (
    _ATTENUATOR_NAMES,
    _NO_ATTENUATION,
    _attenuation_factor,
    _attenuator_name,
    _load_custom_attenuation_coefficients,
    _load_default_attenuation_coefficients,
    _load_attenuation_coefficients,
    _parse_timestamped_formula_blocks,
    _run_start_time,
    _select_effective_block,
)

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py
from drtsans.samplelogs import SampleLogs


def test_attenuation_factor(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # Test input and expected values provided by Lisa Debeer-Schmitt, 2020-02-26
    wavelength = 4.75
    attenuator = 6  # x2k
    expected_value = 0.001037270673420313
    expected_error = 7.005200329552345e-05

    # Add sample logs
    SampleLogs(ws).insert("wavelength", wavelength, "Angstrom")
    SampleLogs(ws).insert("attenuator", attenuator)
    SampleLogs(ws).insert("run_start", "2000-01-01T00:00:00")

    value, error = attenuation_factor(ws)
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


def test_attenuation_factor_open_close(generic_workspace, clean_workspace):
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # add wavelength
    SampleLogs(ws).insert("wavelength", 1.54, "Angstrom")

    # add Undefined attenuator
    attenuator = 0  # Undefined
    SampleLogs(ws).insert("attenuator", attenuator)
    assert attenuation_factor(ws) == (1, 0)

    # add Close attenuator
    attenuator = 1  # Close
    SampleLogs(ws).insert("attenuator", attenuator)
    assert attenuation_factor(ws) == (1, 0)

    # add Open attenuator
    attenuator = 2  # Open
    SampleLogs(ws).insert("attenuator", attenuator)
    assert attenuation_factor(ws) == (1, 0)


def test_attenuation_factor_missing_logs(generic_workspace, clean_workspace):
    """This test that correct error messages are return if the required
    attenuator or wavelenght logs is missing
    """
    ws = generic_workspace  # friendly name
    clean_workspace(ws)

    # Missing attenuator and wavelength log
    with pytest.raises(RuntimeError) as excinfo:
        attenuation_factor(ws)
    assert "attenuator" in str(excinfo.value)  # Should complain about missing attenuator

    # Add in attenuator so only missing wavelength log
    SampleLogs(ws).insert("attenuator", 4)
    SampleLogs(ws).insert("run_start", "2000-01-01T00:00:00")
    with pytest.raises(RuntimeError) as excinfo:
        attenuation_factor(ws)
    assert "wavelength" in str(excinfo.value)  # Should complain about missing wavelength


# Coefficients (A, A error, B, B error, C, C error) provided by Lisa Debeer-Schmitt, 2020-02-26
REFERENCE_COEFFICIENTS = {
    "x3": (
        0.3733459538730628,
        0.008609717163831113,
        0.08056544906925872,
        0.008241433507695071,
        0.0724341919138054,
        0.01125160779959418,
    ),
    "x30": (
        0.11696573514650677,
        0.006304060228295941,
        0.25014801934427583,
        0.01012469612884642,
        0.003696051816711061,
        0.0003197928933191539,
    ),
    "x300": (
        0.028719247985112162,
        0.0019190523738874328,
        0.3884993528348815,
        0.010600714703684273,
        0.00017081815634872129,
        9.055642884664314e-06,
    ),
    "x2k": (
        0.015510737042113254,
        0.0008301527045697745,
        0.5840982399579384,
        0.010252064767405953,
        6.966839283167031e-05,
        2.260503164358648e-06,
    ),
    "x10k": (
        0.00563013075327734,
        0.0005203265715819975,
        0.6961581698084675,
        0.01938010154115584,
        1.3123049075167468e-05,
        1.5266828654446554e-06,
    ),
    "x100k": (
        0.1439135754790426,
        0.005573924841205431,
        0.30824770207752383,
        0.011967728090637404,
        0.006739099792400909,
        0.0007076026868930973,
    ),
}


def _write_formula_file(path, formula="A * exp(-B * wavelength) + C", effective=None, attenuators=None):
    attenuators = attenuators or {
        "x2k": {
            "A": (0.02, 0.001),
            "B": (0.5, 0.01),
            "C": (1.0e-4, 2.0e-6),
        }
    }
    lines = []
    if effective is not None:
        lines.append(f"[effective {effective}]")
    lines.append(f"formula = {formula}")
    lines.append("")
    for attenuator, parameters in attenuators.items():
        lines.append(f"[attenuator {attenuator}]")
        for parameter, (value, error) in parameters.items():
            lines.append(f"{parameter} = {value}, {error}")
        lines.append("")
    path.write_text("\n".join(lines))
    return path


def _coefficients_as_tuples(formula_block):
    return {
        attenuator: tuple(
            value
            for parameter in formula_block.parameters
            for value in (
                parameter_values[parameter].value,
                parameter_values[parameter].error,
            )
        )
        for attenuator, parameter_values in formula_block.coefficients.items()
    }


def test_default_coefficients_file():
    default_file = importlib.resources.files("drtsans.configuration") / "GPSANS_attenuation_coefficients.txt"
    with importlib.resources.as_file(default_file) as filename:
        coefficients = _parse_timestamped_formula_blocks(filename)[0]
    assert coefficients.effective_date == date(1990, 1, 1)
    assert coefficients.formula == "A * exp(-B * wavelength) + C"
    assert coefficients.parameters == ("A", "B", "C")
    assert _coefficients_as_tuples(coefficients) == REFERENCE_COEFFICIENTS
    assert set(coefficients.coefficients) == set(_ATTENUATOR_NAMES.values()) - _NO_ATTENUATION
    # the default file is used when no file is given
    assert _coefficients_as_tuples(_load_attenuation_coefficients()) == REFERENCE_COEFFICIENTS
    assert (
        _coefficients_as_tuples(_load_default_attenuation_coefficients("2000-01-01T00:00:00"))
        == REFERENCE_COEFFICIENTS
    )


def test_attenuation_factor_custom_file(generic_workspace, clean_workspace, tmp_path):
    ws = generic_workspace
    clean_workspace(ws)
    wavelength = 4.75
    SampleLogs(ws).insert("wavelength", wavelength, "Angstrom")
    SampleLogs(ws).insert("attenuator", 6)  # x2k
    SampleLogs(ws).insert("run_start", "2000-01-01T00:00:00")

    # Custom file with x2k coefficients differing from the default ones
    custom_coefficients = (0.02, 0.001, 0.5, 0.01, 1.0e-4, 2.0e-6)
    coefficients_file = tmp_path / "custom_coefficients.txt"
    _write_formula_file(
        coefficients_file,
        attenuators={
            "x2k": {
                "A": custom_coefficients[0:2],
                "B": custom_coefficients[2:4],
                "C": custom_coefficients[4:6],
            }
        },
    )

    value, error = attenuation_factor(ws, coefficients_file)
    expected_value, expected_error = _attenuation_factor(*custom_coefficients, wavelength)
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)
    # the result differs from the one obtained with the default coefficients
    assert value != pytest.approx(attenuation_factor(ws)[0])

    # attenuator missing from the custom file
    SampleLogs(ws).insert("attenuator", 7)  # x10k
    with pytest.raises(ValueError, match="x10k not found"):
        attenuation_factor(ws, coefficients_file)


@pytest.mark.parametrize(
    "contents, message",
    [
        ("formula = A\n\n[attenuator x2k]\nA = 0.1", "parameter lines must have the form"),
        ("formula = A\n\n[attenuator x2k]\nA = half, 0.01", "could not convert string to float"),
        ("formula = A\n\n[attenuator x2k]\nA = nan, 0.01", "finite numbers"),
        ("formula = A\n\n[attenuator x2k]\nA = 0.1, inf", "finite numbers"),
    ],
)
def test_malformed_coefficients_file(tmp_path, contents, message):
    coefficients_file = tmp_path / "malformed_coefficients.txt"
    coefficients_file.write_text(contents)
    with pytest.raises(ValueError) as excinfo:
        _load_custom_attenuation_coefficients(coefficients_file)
    assert message in str(excinfo.value)


def test_repeated_attenuator_in_coefficients_file(tmp_path):
    coefficients_file = tmp_path / "repeated_coefficients.txt"
    coefficients_file.write_text("formula = A\n\n[attenuator x2k]\nA = 0.1, 0.01\n\n[attenuator x2k]\nA = 0.2, 0.01\n")
    with pytest.raises(ValueError) as excinfo:
        _load_custom_attenuation_coefficients(coefficients_file)
    assert "attenuator x2k is repeated" in str(excinfo.value)


def test_attenuation_factor_nonexistent_file(generic_workspace, clean_workspace, tmp_path):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", 6)  # x2k
    SampleLogs(ws).insert("run_start", "2000-01-01T00:00:00")
    with pytest.raises(FileNotFoundError):
        attenuation_factor(ws, tmp_path / "nonexistent.txt")
    # the file is not read when the beam is not attenuated
    SampleLogs(ws).insert("attenuator", 2)  # Open
    assert attenuation_factor(ws, tmp_path / "nonexistent.txt") == (1, 0)


def test_attenuation_factor_unknown_attenuator(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", 9)
    with pytest.raises(ValueError, match="Unknown attenuator log value 9"):
        attenuation_factor(ws)


def test_timestamp_selection_uses_most_recent_effective_block(tmp_path):
    coefficients_file = tmp_path / "timestamped.txt"
    coefficients_file.write_text(
        textwrap.dedent(
            """
            [effective 1990-01-01]
            formula = A
            [attenuator x2k]
            A = 1.0, 0.1

            [effective 2026-05-01]
            formula = A
            [attenuator x2k]
            A = 2.0, 0.2
            """
        )
    )
    blocks = _parse_timestamped_formula_blocks(coefficients_file)
    assert _select_effective_block(blocks, "2026-04-30T23:59:59").coefficients["x2k"]["A"].value == pytest.approx(1.0)
    assert _select_effective_block(blocks, "2026-05-01T00:00:00").coefficients["x2k"]["A"].value == pytest.approx(2.0)


def test_run_start_time_uses_start_time_before_run_start_and_run_begin(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("run_begin", "2025-01-01T00:00:00")
    SampleLogs(ws).insert("run_start", "2026-01-01T00:00:00")
    SampleLogs(ws).insert("start_time", "2027-01-01T00:00:00")
    assert _run_start_time(ws).year == 2027


def test_run_start_time_uses_run_start_before_run_begin(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("run_begin", "2025-01-01T00:00:00")
    SampleLogs(ws).insert("run_start", "2026-01-01T00:00:00")
    assert _run_start_time(ws).year == 2026


def test_run_start_time_uses_run_begin_fallback(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("run_begin", "2025-01-01T00:00:00")
    assert _run_start_time(ws).year == 2025


def test_default_coefficients_missing_run_timestamp_raises(generic_workspace, clean_workspace):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", 6)
    with pytest.raises(RuntimeError, match="start_time, run_start or run_begin"):
        attenuation_factor(ws)


def test_default_coefficients_before_earliest_effective_date_raises():
    with pytest.raises(ValueError, match="earlier than the earliest attenuation calibration"):
        _load_default_attenuation_coefficients("1989-12-31T23:59:59")


def test_alternate_formula_value_and_uncertainty(tmp_path):
    coefficients_file = tmp_path / "alternate.txt"
    _write_formula_file(
        coefficients_file,
        formula="A + B * wavelength + C * exp(-D * wavelength)",
        attenuators={
            "x2k": {
                "A": (1.0, 0.1),
                "B": (2.0, 0.2),
                "C": (3.0, 0.3),
                "D": (0.5, 0.05),
            }
        },
    )
    formula_block = _load_custom_attenuation_coefficients(coefficients_file)
    value, error = formula_block.compiled.value_function(1.0, 2.0, 3.0, 0.5, 4.0), None
    expected_value = 1.0 + 2.0 * 4.0 + 3.0 * np.exp(-0.5 * 4.0)
    expected_error = np.sqrt(
        0.1**2 + (4.0 * 0.2) ** 2 + (np.exp(-0.5 * 4.0) * 0.3) ** 2 + ((-4.0 * 3.0 * np.exp(-0.5 * 4.0)) * 0.05) ** 2
    )
    evaluated_value, error = attenuation_factor_from_block(formula_block, "x2k", 4.0)
    assert value == pytest.approx(expected_value)
    assert evaluated_value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


def attenuation_factor_from_block(formula_block, attenuator_name, wavelength):
    from drtsans.mono.gpsans.attenuation import _evaluate_formula_with_error

    return _evaluate_formula_with_error(formula_block, attenuator_name, wavelength)


def test_formula_without_wavelength_is_wavelength_independent(generic_workspace, clean_workspace, tmp_path):
    coefficients_file = tmp_path / "constant.txt"
    _write_formula_file(
        coefficients_file,
        formula="A + B",
        attenuators={"x2k": {"A": (0.2, 0.01), "B": (0.3, 0.02)}},
    )
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("attenuator", 6)
    value, error = attenuation_factor(ws, coefficients_file)
    assert value == pytest.approx(0.5)
    assert error == pytest.approx(np.sqrt(0.01**2 + 0.02**2))


def test_additional_supported_functions(tmp_path):
    coefficients_file = tmp_path / "additional-functions.txt"
    _write_formula_file(
        coefficients_file,
        formula="log10(A) + sinh(B) + tanh(C) + abs(D)",
        attenuators={"x2k": {"A": (100.0, 1.0), "B": (0.5, 0.01), "C": (0.25, 0.02), "D": (-0.3, 0.03)}},
    )
    formula_block = _load_custom_attenuation_coefficients(coefficients_file)
    value, error = attenuation_factor_from_block(formula_block, "x2k", 4.0)
    expected_value = np.log10(100.0) + np.sinh(0.5) + np.tanh(0.25) + abs(-0.3)
    expected_error = np.sqrt(
        (1.0 / (100.0 * np.log(10.0)) * 1.0) ** 2
        + (np.cosh(0.5) * 0.01) ** 2
        + ((1.0 / np.cosh(0.25) ** 2) * 0.02) ** 2
        + ((-1.0) * 0.03) ** 2
    )
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


@pytest.mark.parametrize(
    "contents, message",
    [
        ("formula = A + D\n[attenuator x2k]\nA = 0.1, 0.01\n", "D is not a parameter"),
        ("formula = A\n[attenuator x2k]\nA = 0.1, 0.01\nB = 0.2, 0.02\n", "parameters not used"),
        ("formula = A\n[attenuator x2k]\nA = 0.1, 0.01\nA = 0.2, 0.02\n", "parameter A is repeated"),
        (
            "formula = A + B\n[attenuator x2k]\nA = 0.1, 0.01\nB = 0.2, 0.02\n[attenuator x30]\nA = 0.3, 0.03\n",
            "does not define the same parameters",
        ),
        ("formula = exp()\n[attenuator x2k]\nA = 0.1, 0.01\n", "invalid call to exp"),
        ("x2k,0.1,0.01,0.5,0.01,0.001,0.0001\n", "old comma-separated attenuation coefficients"),
    ],
)
def test_invalid_formula_blocks(tmp_path, contents, message):
    coefficients_file = tmp_path / "invalid.txt"
    coefficients_file.write_text(contents)
    with pytest.raises(ValueError, match=message):
        _load_custom_attenuation_coefficients(coefficients_file)


def test_duplicate_effective_date_raises(tmp_path):
    coefficients_file = tmp_path / "duplicate-effective.txt"
    coefficients_file.write_text(
        textwrap.dedent(
            """
            [effective 1990-01-01]
            formula = A
            [attenuator x2k]
            A = 1.0, 0.1

            [effective 1990-01-01]
            formula = A
            [attenuator x2k]
            A = 2.0, 0.2
            """
        )
    )
    with pytest.raises(ValueError, match="effective date 1990-01-01 is repeated"):
        _parse_timestamped_formula_blocks(coefficients_file)


@pytest.mark.parametrize(
    "log_value, name",
    [
        (0, "Undefined"),
        (1, "Close"),
        (2, "Open"),
        (3, "x3"),
        (4, "x30"),
        (5, "x300"),
        (6, "x2k"),
        (7, "x10k"),
        (8, "x100k"),
    ],
)
def test_attenuator_name(generic_workspace, clean_workspace, log_value, name):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("attenuator", log_value)
    assert _attenuator_name(ws) == name


@pytest.mark.parametrize("log_value", [9, 5.5, 52.9971])
def test_attenuator_name_invalid(generic_workspace, clean_workspace, log_value):
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("attenuator", log_value)
    with pytest.raises(ValueError, match="Unknown attenuator log value"):
        _attenuator_name(ws)


def test_attenuator_name_negative(generic_workspace, clean_workspace):
    """A negative log value, the x2k stage position (mm) of a run converted from a SPICE file, is Undefined"""
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", -195.999101)
    assert _attenuator_name(ws) == "Undefined"
    assert attenuation_factor(ws) == (1, 0)


@pytest.mark.parametrize("log_value", [0.5, 1.5, 2.5])
def test_attenuator_name_non_integer_below_three(generic_workspace, clean_workspace, log_value):
    """A non-integer log value below 3, the average of Undefined, Close and Open entries, is Undefined"""
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("wavelength", 4.75, "Angstrom")
    SampleLogs(ws).insert("attenuator", log_value)
    assert _attenuator_name(ws) == "Undefined"
    assert attenuation_factor(ws) == (1, 0)


if __name__ == "__main__":
    pytest.main([__file__])
