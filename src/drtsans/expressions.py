"""Restricted parsing of mathematical expressions read from configuration text.

Expressions are parsed with :py:mod:`ast` and translated node by node into SymPy, so only numbers, declared
names, arithmetic operators and a fixed list of functions are accepted. Nothing in the expression text is passed to
:py:func:`eval`, and invalid expressions are rejected when they are parsed rather than when they are evaluated.

Callers choose the names of the independent variables (for example ``wavelength`` or ``x``) and, optionally, the
names of fitted parameters whose uncertainties are propagated to first order.
"""

import ast
from dataclasses import dataclass, field
import math
from typing import Callable, Dict, FrozenSet, Optional, Sequence, Tuple, Union
import unicodedata

import numpy as np
import sympy

__all__ = [
    "Expression",
    "ExpressionError",
    "SUPPORTED_CONSTANTS",
    "SUPPORTED_FUNCTIONS",
    "parse_expression",
    "reserved_names",
]


def _log10(argument):
    return sympy.log(argument, 10)


SUPPORTED_FUNCTIONS = {
    "abs": sympy.Abs,
    "acos": sympy.acos,
    "asin": sympy.asin,
    "atan": sympy.atan,
    "cos": sympy.cos,
    "cosh": sympy.cosh,
    "erf": sympy.erf,
    "exp": sympy.exp,
    "log": sympy.log,
    "log10": _log10,
    "sin": sympy.sin,
    "sinh": sympy.sinh,
    "sqrt": sympy.sqrt,
    "tan": sympy.tan,
    "tanh": sympy.tanh,
}
SUPPORTED_CONSTANTS = {"pi": sympy.pi}
_LAMBDIFY_MODULES = [{"erf": math.erf}, "numpy"]
# Largest decimal exponent of a 64-bit float. Constant powers beyond it cannot be evaluated as floats, and SymPy
# would compute them exactly while parsing, which takes unbounded time and memory.
_MAX_DECIMAL_EXPONENT = 308


class ExpressionError(ValueError):
    """Report an invalid expression.

    The message carries no file or line context; callers add it before logging and re-raising.

    Parameters
    ----------
    message
        Human-readable validation failure.
    unknown_name
        Name that is neither a declared variable, a declared parameter, nor a supported constant, when that is the
        cause of the failure. Callers use it to build a message specific to their configuration format.
    """

    def __init__(self, message: str, unknown_name: Optional[str] = None):
        super().__init__(message)
        self.unknown_name = unknown_name


@dataclass(frozen=True)
class Expression:
    """Represent a validated expression ready for numerical evaluation.

    Parameters
    ----------
    text
        Expression as given by the caller, converted to a string.
    expression
        Symbolic SymPy expression.
    variables
        Ordered names of the independent variables.
    parameters
        Ordered names of the fitted parameters. Empty when the expression has no parameters.
    value_function
        Lambdified callable that evaluates ``expression``. Arguments are the parameter values followed by the
        variable values, in the order of ``parameters`` and ``variables``.
    derivative_functions
        Lambdified partial derivatives keyed by parameter name, with the same arguments as ``value_function``.
        Empty unless the expression was parsed with ``with_uncertainty=True``.
    """

    text: str
    expression: sympy.Expr
    variables: Tuple[str, ...]
    parameters: Tuple[str, ...]
    value_function: Callable
    derivative_functions: Dict[str, Callable] = field(default_factory=dict)

    def uses(self, name: str) -> bool:
        """Return whether the expression references a variable or parameter.

        Parameters
        ----------
        name
            Variable or parameter name.

        Returns
        -------
        bool
            :py:obj:`True` if ``name`` appears in the expression.
        """
        return name in {symbol.name for symbol in self.expression.free_symbols}

    def __call__(self, *variable_values: float) -> float:
        """Evaluate an expression that has no parameters.

        Parameters
        ----------
        variable_values
            One value per variable, in the order of ``variables``.

        Returns
        -------
        float
            Expression value. Division by zero or the fractional power of a negative number yields ``inf`` or
            ``nan`` instead of raising, so callers must check the result.

        Raises
        ------
        TypeError
            If the expression has parameters, or the number of values differs from the number of variables.
        """
        if self.parameters:
            raise TypeError("expressions with parameters are evaluated with evaluate_with_uncertainty")
        if len(variable_values) != len(self.variables):
            raise TypeError(f"expected {len(self.variables)} variable values, got {len(variable_values)}")
        return _evaluate(self.value_function, variable_values)

    def evaluate_with_uncertainty(
        self,
        parameter_values: Sequence[float],
        parameter_errors: Sequence[float],
        variable_values: Sequence[float],
    ) -> Tuple[float, float]:
        """Evaluate the expression and its uncertainty propagated from the parameter uncertainties.

        The parameters are treated as independent, so the variance is
        :math:`\\sum_p (\\partial f / \\partial p)^2 \\sigma_p^2`. Variables are taken as exact.

        Parameters
        ----------
        parameter_values
            One value per parameter, in the order of ``parameters``.
        parameter_errors
            One-sigma uncertainty per parameter, in the order of ``parameters``.
        variable_values
            One value per variable, in the order of ``variables``.

        Returns
        -------
        tuple of float
            Expression value and its uncertainty. Either may be ``inf`` or ``nan``; callers must check them.

        Raises
        ------
        TypeError
            If the expression was parsed without uncertainty support, or the numbers of values do not match.
        """
        if self.parameters and not self.derivative_functions:
            raise TypeError("expression was parsed without uncertainty support")
        if len(parameter_values) != len(self.parameters) or len(parameter_errors) != len(self.parameters):
            raise TypeError(f"expected {len(self.parameters)} parameter values and uncertainties")
        if len(variable_values) != len(self.variables):
            raise TypeError(f"expected {len(self.variables)} variable values, got {len(variable_values)}")
        arguments = list(parameter_values) + list(variable_values)
        value = _evaluate(self.value_function, arguments)
        variance = 0.0
        for parameter, error in zip(self.parameters, parameter_errors):
            derivative = _evaluate(self.derivative_functions[parameter], arguments)
            variance += (derivative * error) ** 2
        return value, float(np.sqrt(variance))


def reserved_names(variables: Sequence[str]) -> FrozenSet[str]:
    """Return the names that cannot be used for parameters.

    Parameters
    ----------
    variables
        Names of the independent variables.

    Returns
    -------
    frozenset of str
        Variable names, supported function names and supported constant names.
    """
    return frozenset({*variables, *SUPPORTED_FUNCTIONS, *SUPPORTED_CONSTANTS})


def parse_expression(
    value: Union[str, int, float],
    variables: Sequence[str],
    parameters: Sequence[str] = (),
    with_uncertainty: bool = False,
) -> Expression:
    """Parse an expression and compile it for numerical evaluation.

    Supported syntax is numbers, the names in ``variables`` and ``parameters``, the operators ``+``, ``-``, ``*``,
    ``/`` and ``**``, unary ``+`` and ``-``, parentheses, the functions in :py:data:`SUPPORTED_FUNCTIONS` and the
    constant ``pi``.

    Parameters
    ----------
    value
        Expression text, or a plain number for a constant expression.
    variables
        Names of the independent variables, for example ``("wavelength",)``. Variables are optional in the
        expression.
    parameters
        Names of fitted parameters. Parameters are optional in the expression.
    with_uncertainty
        If :py:obj:`True`, also compile the partial derivatives with respect to each parameter, which
        :py:meth:`Expression.evaluate_with_uncertainty` requires.

    Returns
    -------
    Expression
        Validated expression with its numerical callables.

    Raises
    ------
    ExpressionError
        If ``value`` is not a string or number, uses unsupported syntax, references undeclared names, has a
        complex value, or has a constant power beyond the range of 64-bit floats, or if a variable or parameter
        name is invalid or conflicts with a reserved name.
    """
    variables = tuple(variables)
    parameters = tuple(parameters)
    _validate_names(variables, parameters)

    variable_symbols = {name: sympy.Symbol(name, real=True) for name in variables}
    parameter_symbols = {name: sympy.Symbol(name, real=True) for name in parameters}
    symbols = {**parameter_symbols, **variable_symbols}

    # bool is a subclass of int, but True or False in a configuration is a mistake, not a number
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ExpressionError(f"expressions must be strings or numbers, not {type(value).__name__}")
    if isinstance(value, str):
        try:
            parsed = ast.parse(value.strip(), mode="eval")
        except SyntaxError as error:
            raise ExpressionError(f"invalid formula syntax: {error.msg}") from None
        expression = _ast_to_sympy(parsed.body, symbols)
    else:
        expression = sympy.Float(value) if isinstance(value, float) else sympy.Integer(value)

    # SymPy simplifies constant parts while parsing, so sqrt(-1) becomes I and 1/0 or log(0) becomes complex infinity
    if expression.has(sympy.I):
        raise ExpressionError("formula has a complex value")
    if expression.has(sympy.zoo):
        raise ExpressionError("formula has an infinite value of undefined sign, for example from a division by zero")

    arguments = [parameter_symbols[name] for name in parameters] + [variable_symbols[name] for name in variables]
    value_function = sympy.lambdify(arguments, expression, modules=_LAMBDIFY_MODULES)
    derivative_functions = {}
    if with_uncertainty:
        derivative_functions = {
            name: sympy.lambdify(arguments, sympy.diff(expression, parameter_symbols[name]), modules=_LAMBDIFY_MODULES)
            for name in parameters
        }
    return Expression(
        text=str(value),
        expression=expression,
        variables=variables,
        parameters=parameters,
        value_function=value_function,
        derivative_functions=derivative_functions,
    )


def _validate_names(variables: Tuple[str, ...], parameters: Tuple[str, ...]) -> None:
    """Check that variable and parameter names are distinct identifiers that do not shadow reserved names.

    Parameters
    ----------
    variables
        Names of the independent variables.
    parameters
        Names of the fitted parameters.

    Raises
    ------
    ExpressionError
        If a name is not an identifier, is not in Unicode NFKC normal form, is repeated, or is a supported function
        or constant name.
    """
    names = variables + parameters
    for name in names:
        # Python converts identifiers to NFKC form when parsing, so a name such as fullwidth "\uff2b" would be
        # read as "K" in the expression and no longer match its declaration
        if not name.isidentifier() or unicodedata.normalize("NFKC", name) != name:
            raise ExpressionError(f"{name!r} is not a valid variable or parameter name")
        if name in SUPPORTED_FUNCTIONS or name in SUPPORTED_CONSTANTS:
            raise ExpressionError(f"{name} is a reserved name")
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise ExpressionError(f"variable and parameter names are repeated: {repeated}")


def _evaluate(function: Callable, arguments: Sequence[float]) -> float:
    """Evaluate a lambdified function with NumPy scalar arguments.

    Parameters
    ----------
    function
        Lambdified callable.
    arguments
        Positional argument values.

    Returns
    -------
    float
        Function value.
    """
    # Lambdify emits plain operators for division and powers. With NumPy scalars, division by zero and fractional
    # powers of negative numbers yield inf or nan, which callers report, instead of raising ZeroDivisionError or
    # returning a complex number as Python floats do.
    with np.errstate(all="ignore"):
        return float(function(*[np.float64(argument) for argument in arguments]))


def _check_power_range(base: sympy.Expr, exponent: sympy.Expr) -> None:
    """Reject a constant power whose magnitude is beyond the range of 64-bit floats.

    Parameters
    ----------
    base
        Constant base of the power.
    exponent
        Constant exponent of the power.

    Raises
    ------
    ExpressionError
        If :math:`|\\text{exponent} \\cdot \\log_{10}|\\text{base}||` exceeds the largest decimal exponent of a
        64-bit float.
    """
    magnitude = float(sympy.Abs(base))
    if magnitude == 0.0 or not math.isfinite(magnitude):
        return
    decimal_exponent = float(sympy.Abs(exponent)) * abs(math.log10(magnitude))
    if decimal_exponent > _MAX_DECIMAL_EXPONENT:
        raise ExpressionError("formula constant is out of range")


def _ast_to_sympy(node: ast.AST, symbols: Dict[str, sympy.Symbol]) -> sympy.Expr:
    """Translate an allowed Python expression AST into a SymPy expression.

    Parameters
    ----------
    node
        AST node from ``ast.parse(..., mode="eval")``.
    symbols
        Declared variable and parameter names mapped to SymPy symbols.

    Returns
    -------
    sympy.Expr
        Symbolic expression equivalent to ``node``.

    Raises
    ------
    ExpressionError
        If ``node`` contains unsupported syntax, unsupported function calls, undeclared names, or a constant power
        beyond the range of 64-bit floats.
    """
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return sympy.Float(node.value) if isinstance(node.value, float) else sympy.Integer(node.value)
        raise ExpressionError("formula constants must be numbers")
    if isinstance(node, ast.Name):
        if node.id in symbols:
            return symbols[node.id]
        if node.id in SUPPORTED_CONSTANTS:
            return SUPPORTED_CONSTANTS[node.id]
        raise ExpressionError(f"{node.id} is not an allowed name", unknown_name=node.id)
    if isinstance(node, ast.BinOp):
        left = _ast_to_sympy(node.left, symbols)
        right = _ast_to_sympy(node.right, symbols)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.Pow):
            if left.is_number and right.is_number:
                _check_power_range(left, right)
            return left**right
        raise ExpressionError("unsupported formula operator")
    if isinstance(node, ast.UnaryOp):
        operand = _ast_to_sympy(node.operand, symbols)
        if isinstance(node.op, ast.UAdd):
            return operand
        if isinstance(node.op, ast.USub):
            return -operand
        raise ExpressionError("unsupported unary formula operator")
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in SUPPORTED_FUNCTIONS:
            raise ExpressionError("formula calls are limited to supported functions")
        if node.keywords:
            raise ExpressionError("formula functions do not accept keyword arguments")
        arguments = [_ast_to_sympy(argument, symbols) for argument in node.args]
        try:
            return SUPPORTED_FUNCTIONS[node.func.id](*arguments)
        except (TypeError, ValueError) as error:
            raise ExpressionError(f"invalid call to {node.func.id}: {error}") from None
    raise ExpressionError("formula contains unsupported syntax")
