import math

import numpy as np
import pytest

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/expressions.py
from drtsans.expressions import (
    SUPPORTED_CONSTANTS,
    SUPPORTED_FUNCTIONS,
    ExpressionError,
    parse_expression,
    reserved_names,
)


def _evaluate(text, parameters, wavelength=4.0):
    """Parse ``text`` with fitted parameters and ``wavelength``, and evaluate value and uncertainty."""
    expression = parse_expression(text, variables=("wavelength",), parameters=tuple(parameters), with_uncertainty=True)
    values = [value for value, _ in parameters.values()]
    errors = [error for _, error in parameters.values()]
    return expression.evaluate_with_uncertainty(values, errors, [wavelength])


def test_value_and_uncertainty():
    parameters = {"A": (1.0, 0.1), "B": (2.0, 0.2), "C": (3.0, 0.3), "D": (0.5, 0.05)}
    value, error = _evaluate("A + B * wavelength + C * exp(-D * wavelength)", parameters)
    expected_value = 1.0 + 2.0 * 4.0 + 3.0 * np.exp(-0.5 * 4.0)
    expected_error = np.sqrt(
        0.1**2 + (4.0 * 0.2) ** 2 + (np.exp(-0.5 * 4.0) * 0.3) ** 2 + ((-4.0 * 3.0 * np.exp(-0.5 * 4.0)) * 0.05) ** 2
    )
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


def test_value_function_argument_order():
    expression = parse_expression("A + B * wavelength", variables=("wavelength",), parameters=("A", "B"))
    # parameters first, in declaration order, then variables
    assert expression.value_function(1.0, 2.0, 4.0) == pytest.approx(9.0)
    assert expression.derivative_functions == {}


def test_additional_supported_functions():
    parameters = {"A": (100.0, 1.0), "B": (0.5, 0.01), "C": (0.25, 0.02), "D": (-0.3, 0.03)}
    value, error = _evaluate("log10(A) + sinh(B) + tanh(C) + abs(D)", parameters)
    expected_value = np.log10(100.0) + np.sinh(0.5) + np.tanh(0.25) + abs(-0.3)
    expected_error = np.sqrt(
        (1.0 / (100.0 * np.log(10.0)) * 1.0) ** 2
        + (np.cosh(0.5) * 0.01) ** 2
        + ((1.0 / np.cosh(0.25) ** 2) * 0.02) ** 2
        + ((-1.0) * 0.03) ** 2
    )
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


def test_erf_value_and_uncertainty():
    value, error = _evaluate("erf(A)", {"A": (0.5, 0.1)})
    assert value == pytest.approx(math.erf(0.5))
    assert error == pytest.approx((2.0 / math.sqrt(math.pi)) * math.exp(-(0.5**2)) * 0.1)


@pytest.mark.parametrize("name", sorted(SUPPORTED_FUNCTIONS))
def test_every_supported_function_evaluates(name):
    expression = parse_expression(f"{name}(x)", variables=("x",))
    assert math.isfinite(expression(0.5))


@pytest.mark.parametrize(
    "text, parameters, non_finite",
    [
        ("sqrt(A)", {"A": (-1.0, 0.1)}, "value"),
        ("sqrt(A)", {"A": (0.0, 0.1)}, "uncertainty"),
        ("A / B", {"A": (1.0, 0.1), "B": (0.0, 0.1)}, "value"),
        ("A ** 0.5", {"A": (-1.0, 0.1)}, "value"),
        ("A ** 0.5", {"A": (0.0, 0.1)}, "uncertainty"),
    ],
)
def test_non_finite_results_are_returned(text, parameters, non_finite):
    # callers decide how to report non-finite results, so evaluation returns them instead of raising
    value, error = _evaluate(text, parameters)
    if non_finite == "value":
        assert not math.isfinite(value)
    else:
        assert math.isfinite(value)
        assert not math.isfinite(error)


def test_caller_chosen_variable_name():
    for variable in ("x", "wavelength"):
        expression = parse_expression(f"0.95 - 0.01*({variable} - 16)", variables=(variable,))
        assert expression(6.0) == pytest.approx(1.05)
        assert expression.uses(variable)
    with pytest.raises(ExpressionError, match="wavelength is not an allowed name") as exc_info:
        parse_expression("0.95 - 0.01*(wavelength - 16)", variables=("x",))
    assert exc_info.value.unknown_name == "wavelength"


def test_constant_expression_ignores_variable():
    expression = parse_expression("2 * pi", variables=("x",))
    assert not expression.uses("x")
    assert expression(1.0) == pytest.approx(2.0 * math.pi)
    assert expression(10.0) == pytest.approx(2.0 * math.pi)


def test_surrounding_whitespace_is_ignored():
    assert parse_expression("  0.9 + x \n", variables=("x",))(0.05) == pytest.approx(0.95)


@pytest.mark.parametrize("text, expected", [("+x", 2.0), ("-x", -2.0), ("-(-x)", 2.0), ("2 * -x", -4.0)])
def test_unary_operators(text, expected):
    assert parse_expression(text, variables=("x",))(2.0) == pytest.approx(expected)


@pytest.mark.parametrize("value, expected", [(0.9, 0.9), (-1, -1.0), (2, 2.0)])
def test_number_input(value, expected):
    assert parse_expression(value, variables=("x",))(6.0) == pytest.approx(expected)


@pytest.mark.parametrize("value", [True, False, None, [0.9], {"x": 1}])
def test_non_string_non_number_input_is_rejected(value):
    with pytest.raises(ExpressionError, match="expressions must be strings or numbers"):
        parse_expression(value, variables=("x",))


@pytest.mark.parametrize(
    "text, message",
    [
        ("x =", "invalid formula syntax"),
        ("'0.9'", "formula constants must be numbers"),
        ("x // 2", "unsupported formula operator"),
        ("not x", "unsupported unary formula operator"),
        ("gamma(x)", "formula calls are limited to supported functions"),
        ("__import__('os')", "formula calls are limited to supported functions"),
        ("np.exp(x)", "formula calls are limited to supported functions"),
        ("exp(x=1)", "formula functions do not accept keyword arguments"),
        ("exp()", "invalid call to exp"),
        ("x.real", "formula contains unsupported syntax"),
        ("x[0]", "formula contains unsupported syntax"),
        ("x < 1", "formula contains unsupported syntax"),
        ("lambda: x", "formula contains unsupported syntax"),
    ],
)
def test_unsupported_syntax(text, message):
    with pytest.raises(ExpressionError, match=message) as exc_info:
        parse_expression(text, variables=("x",))
    assert exc_info.value.unknown_name is None


def test_unknown_name_is_reported():
    with pytest.raises(ExpressionError) as exc_info:
        parse_expression("A + D * wavelength", variables=("wavelength",), parameters=("A",))
    assert str(exc_info.value) == "D is not an allowed name"
    assert exc_info.value.unknown_name == "D"


def test_reserved_names():
    names = reserved_names(("wavelength",))
    assert "wavelength" in names
    assert set(SUPPORTED_FUNCTIONS) <= names
    assert set(SUPPORTED_CONSTANTS) <= names
    assert "x" not in names


@pytest.mark.parametrize(
    "variables, parameters, message",
    [
        (("x",), ("pi",), "pi is a reserved name"),
        (("exp",), (), "exp is a reserved name"),
        (("x",), ("1A",), "is not a valid variable or parameter name"),
        (("x",), ("x",), "names are repeated"),
    ],
)
def test_invalid_declared_names(variables, parameters, message):
    with pytest.raises(ExpressionError, match=message):
        parse_expression("1", variables=variables, parameters=parameters)


def test_evaluation_argument_errors():
    with_parameters = parse_expression("A * x", variables=("x",), parameters=("A",))
    with pytest.raises(TypeError, match="evaluate_with_uncertainty"):
        with_parameters(1.0)
    with pytest.raises(TypeError, match="without uncertainty support"):
        with_parameters.evaluate_with_uncertainty([1.0], [0.1], [1.0])

    without_parameters = parse_expression("x", variables=("x",))
    with pytest.raises(TypeError, match="expected 1 variable values"):
        without_parameters(1.0, 2.0)

    with_uncertainty = parse_expression("A * x", variables=("x",), parameters=("A",), with_uncertainty=True)
    with pytest.raises(TypeError, match="expected 1 parameter values"):
        with_uncertainty.evaluate_with_uncertainty([1.0, 2.0], [0.1], [1.0])
    with pytest.raises(TypeError, match="expected 1 variable values, got 2"):
        with_uncertainty.evaluate_with_uncertainty([1.0], [0.1], [1.0, 2.0])
