import ast
from dataclasses import dataclass, field
from datetime import date, datetime
import importlib.resources
import keyword
import math
from pathlib import Path
from typing import Dict, Optional, Tuple, Union
from zoneinfo import ZoneInfo

from dateutil.parser import parse as parse_date
from mantid.api import MatrixWorkspace
from mantid.kernel import logger
import numpy as np
import sympy

from drtsans.samplelogs import SampleLogs

__all__ = ["attenuation_factor"]

_ATTENUATOR_NAMES = {
    0: "Undefined",
    1: "Close",
    2: "Open",
    3: "x3",
    4: "x30",
    5: "x300",
    6: "x2k",
    7: "x10k",
    8: "x100k",
}
_NO_ATTENUATION = {"Undefined", "Close", "Open"}
_DEFAULT_COEFFICIENTS_FILE = "GPSANS_attenuation_coefficients.txt"
_INSTRUMENT_TIME_ZONE = ZoneInfo("America/New_York")


def _log10(argument):
    return sympy.log(argument, 10)


_SUPPORTED_FUNCTIONS = {
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
_SUPPORTED_CONSTANTS = {"pi": sympy.pi}
_RESERVED_NAMES = frozenset({"wavelength", *_SUPPORTED_FUNCTIONS, *_SUPPORTED_CONSTANTS})


@dataclass(frozen=True)
class _ParameterValue:
    """Store one fitted parameter from an attenuation coefficients file.

    Parameters
    ----------
    value
        Best-fit parameter value.
    error
        One-sigma uncertainty of ``value``.
    line_number
        Source line where the parameter was declared, used in validation diagnostics.
    """

    value: float
    error: float
    line_number: int


@dataclass
class _FormulaBlockBuilder:
    """Accumulate raw parser state for one attenuation formula block.

    Parameters
    ----------
    source
        Path or label of the file being parsed, used in error messages.
    effective_date
        Effective date for a default-file block. Custom blocks leave this unset.
    effective_line_number
        Source line for ``effective_date``.
    formula
        Raw formula string from the ``formula =`` line.
    formula_line_number
        Source line for ``formula``.
    coefficients
        Attenuator names mapped to raw parameter declarations.
    attenuator_line_numbers
        Source lines for attenuator block headers.
    """

    source: Union[str, Path]
    effective_date: Optional[date] = None
    effective_line_number: Optional[int] = None
    formula: Optional[str] = None
    formula_line_number: Optional[int] = None
    coefficients: Dict[str, Dict[str, _ParameterValue]] = field(default_factory=dict)
    attenuator_line_numbers: Dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class _CompiledFormula:
    """Represent a validated attenuation expression ready for numerical evaluation.

    Parameters
    ----------
    expression
        Symbolic SymPy expression for the attenuation factor.
    parameters
        Ordered fitted parameter names declared by the first attenuator block.
    uses_wavelength
        Whether ``expression`` references the reserved variable ``wavelength``.
    value_function
        Lambdified callable that evaluates ``expression``.
    derivative_functions
        Lambdified partial derivatives keyed by fitted parameter name.
    """

    expression: sympy.Expr
    parameters: Tuple[str, ...]
    uses_wavelength: bool
    value_function: object
    derivative_functions: Dict[str, object]


@dataclass(frozen=True)
class _FormulaBlock:
    """Store a validated attenuation formula block.

    Parameters
    ----------
    formula
        Raw formula string as written in the coefficients file.
    parameters
        Ordered fitted parameter names. Every attenuator in the block defines these names.
    coefficients
        Attenuator names mapped to parameter values and uncertainties.
    compiled
        Compiled formula and derivative callables.
    effective_date
        Effective date for default-file blocks. Custom-file blocks leave this unset.
    """

    formula: str
    parameters: Tuple[str, ...]
    coefficients: Dict[str, Dict[str, _ParameterValue]]
    compiled: _CompiledFormula
    effective_date: Optional[date] = None

    def coefficients_for_log(self) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Return attenuation coefficients in the nested shape written to reduction logs.

        Returns
        -------
        dict
            Mapping of attenuator name to parameter name to ``{"value": value, "error": error}``.
        """
        return {
            attenuator: {
                parameter: {"value": parameter_value.value, "error": parameter_value.error}
                for parameter, parameter_value in parameters.items()
            }
            for attenuator, parameters in self.coefficients.items()
        }


def attenuation_factor(
    input_workspace: Union[str, MatrixWorkspace], coefficients_file: Optional[Union[str, Path]] = None
) -> Tuple[float, float]:
    r"""Return the attenuation factor and uncertainty for a GPSANS workspace.

    The attenuation formula and fitted parameter values are read from ``coefficients_file``. If
    ``coefficients_file`` is :py:obj:`None`, the timestamped packaged default file is selected by the
    workspace ``start_time`` sample log, with ``run_start`` and ``run_begin`` as fallbacks. Custom files are
    untimestamped formula blocks and always override the default file.

    The formula may depend on wavelength, in Å, through the optional variable ``wavelength``. Formula
    uncertainty is propagated from the independent fitted-parameter uncertainties. See
    :ref:`user.corrections.attenuation` for the file format.

    If the attenuator is one of Undefined, Close or Open then a
    attenuation factor of 1 with uncertainty 0 is returned. A negative log value, as found in runs converted from
    SPICE files, and a non-integer log value between 0 and 3 are taken as Undefined, with a warning.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Input workspace for which to calculate attenuation factor
    coefficients_file: str, ~pathlib.Path
        Path to the custom file of attenuator fit coefficients. If :py:obj:`None`, the timestamped default file
        in the ``drtsans.configuration`` package is used.

    Returns
    -------
    float, float
        attenuation_factor, attenuation_factor_uncertainty

    Raises
    ------
    RuntimeError
        If the workspace has no "attenuator" sample log, no run timestamp for default coefficients, or no
        "wavelength" sample log when the selected formula depends on wavelength.
    FileNotFoundError
        If the coefficients file does not exist
    ValueError
        If the attenuator log value is unknown, if the attenuator is missing from the coefficients file,
        or if the coefficients file is malformed.
    """
    attenuator_name = _attenuator_name(input_workspace)

    if attenuator_name in _NO_ATTENUATION:
        return 1, 0

    if coefficients_file is None:
        formula_block = _load_default_attenuation_coefficients(_run_start_time(input_workspace))
    else:
        formula_block = _load_custom_attenuation_coefficients(coefficients_file)
    return _attenuator_transmission(input_workspace, attenuator_name, formula_block, coefficients_file)


def _attenuator_transmission(
    input_workspace: Union[str, MatrixWorkspace],
    attenuator_name: str,
    coefficients: _FormulaBlock,
    coefficients_file: Optional[Union[str, Path]] = None,
) -> Tuple[float, float]:
    """Evaluate the transmitted fraction for one attenuator.

    Parameters
    ----------
    input_workspace
        Workspace containing the ``wavelength`` sample log when the selected formula uses wavelength.
    attenuator_name
        Attenuator name produced by :func:`_attenuator_name`.
    coefficients
        Validated attenuation formula block selected for the reduction.
    coefficients_file
        Source file path used only for error messages. ``None`` means the packaged default file.

    Returns
    -------
    tuple of float
        Attenuation factor and propagated uncertainty.

    Raises
    ------
    ValueError
        If ``attenuator_name`` is not present in ``coefficients``.
    RuntimeError
        If the formula uses wavelength and the workspace has no ``wavelength`` sample log.
    """
    if attenuator_name in _NO_ATTENUATION:
        return 1, 0

    if attenuator_name not in coefficients.coefficients:
        if coefficients_file is None:
            coefficients_file = _DEFAULT_COEFFICIENTS_FILE
        message = f"Attenuator {attenuator_name} not found in the attenuation coefficients file {coefficients_file}"
        logger.error(message)
        raise ValueError(message)

    wavelength = (
        SampleLogs(input_workspace).single_value("wavelength") if coefficients.compiled.uses_wavelength else 0.0
    )
    return _evaluate_formula_with_error(coefficients, attenuator_name, wavelength)


def _attenuator_name(input_workspace: Union[str, MatrixWorkspace]) -> str:
    """Return the attenuation-file name for the workspace attenuator.

    Parameters
    ----------
    input_workspace
        Workspace containing the ``attenuator`` sample log.

    Returns
    -------
    str
        Attenuator name used in the coefficients file.

    Raises
    ------
    RuntimeError
        If the workspace has no ``attenuator`` sample log.
    ValueError
        If the log value is not one of the supported GPSANS attenuator indices.
    """
    attenuator = SampleLogs(input_workspace).single_value("attenuator")
    if attenuator < 0:
        logger.warning(
            f"Negative attenuator log value {attenuator}, probably an attenuator stage position (mm) "
            "from a SPICE file. The attenuator is taken as Undefined and no attenuation correction is applied"
        )
        return _ATTENUATOR_NAMES[0]
    if attenuator < 3 and attenuator not in _ATTENUATOR_NAMES:
        logger.warning(
            f"Non-integer attenuator log value {attenuator}, probably the attenuator moved between the Undefined, "
            "Close and Open positions during the run. The attenuator is taken as Undefined and no attenuation "
            "correction is applied"
        )
        return _ATTENUATOR_NAMES[0]
    if attenuator not in _ATTENUATOR_NAMES:
        message = (
            f"Unknown attenuator log value {attenuator}. Valid values are 0 to 8. "
            "Runs converted from SPICE files may instead store the attenuator stage position (mm)"
        )
        logger.error(message)
        raise ValueError(message)
    return _ATTENUATOR_NAMES[attenuator]


def _load_attenuation_coefficients(
    coefficients_file: Optional[Union[str, Path]] = None,
) -> _FormulaBlock:
    """Load attenuation coefficients through the backward-compatible private entry point.

    Parameters
    ----------
    coefficients_file
        Custom coefficients file to read. If ``None``, the current packaged default block is selected using
        today's instrument-local date.

    Returns
    -------
    _FormulaBlock
        Validated attenuation formula block.

    Notes
    -----
    Runtime attenuation calculations should call :func:`_load_default_attenuation_coefficients` with a run
    timestamp for the packaged default file. This wrapper exists for tests and metadata callers that only need the
    current packaged block.
    """
    if coefficients_file is None:
        return _load_default_attenuation_coefficients(datetime.now(tz=_INSTRUMENT_TIME_ZONE))
    return _load_custom_attenuation_coefficients(coefficients_file)


def _load_default_attenuation_coefficients(run_start: Union[str, datetime]) -> _FormulaBlock:
    """Load the packaged attenuation history and select the block effective for a run.

    Parameters
    ----------
    run_start
        Run timestamp used to select the most recent effective block not later than the run date.

    Returns
    -------
    _FormulaBlock
        Selected default attenuation formula block.

    Raises
    ------
    ValueError
        If the packaged file is malformed or ``run_start`` predates the earliest effective block.
    """
    with importlib.resources.as_file(
        importlib.resources.files("drtsans.configuration") / _DEFAULT_COEFFICIENTS_FILE
    ) as default_file:
        blocks = _parse_timestamped_formula_blocks(default_file)
    return _select_effective_block(blocks, run_start)


def _load_custom_attenuation_coefficients(coefficients_file: Union[str, Path]) -> _FormulaBlock:
    """Load a custom untimestamped attenuation coefficients file.

    Parameters
    ----------
    coefficients_file
        Path to a custom attenuation coefficients file.

    Returns
    -------
    _FormulaBlock
        Validated custom attenuation formula block.

    Raises
    ------
    FileNotFoundError
        If ``coefficients_file`` does not exist.
    ValueError
        If the file is malformed or contains a default-file-only effective header.
    """
    return _parse_formula_block(coefficients_file, allow_effective=False)


def _parse_timestamped_formula_blocks(coefficients_file: Union[str, Path]) -> Tuple[_FormulaBlock, ...]:
    """Parse all effective-dated formula blocks in the packaged coefficients file.

    Parameters
    ----------
    coefficients_file
        Path to the packaged attenuation coefficients history.

    Returns
    -------
    tuple of _FormulaBlock
        Validated formula blocks, one per effective date.

    Raises
    ------
    ValueError
        If the file has malformed blocks or duplicate effective dates.
    """
    builders = []
    current_builder = None
    current_attenuator = None
    effective_dates = {}

    with open(coefficients_file, "r") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            effective_date = _parse_effective_header(line, coefficients_file, line_number)
            if effective_date is not None:
                if effective_date in effective_dates:
                    _raise_parse_error(
                        coefficients_file,
                        line_number,
                        f"effective date {effective_date.isoformat()} is repeated",
                    )
                effective_dates[effective_date] = line_number
                current_builder = _FormulaBlockBuilder(
                    source=coefficients_file,
                    effective_date=effective_date,
                    effective_line_number=line_number,
                )
                builders.append(current_builder)
                current_attenuator = None
                continue

            if current_builder is None:
                _raise_parse_error(
                    coefficients_file,
                    line_number,
                    "default attenuation coefficients must start with [effective YYYY-MM-DD]",
                )
            current_attenuator = _parse_block_line(current_builder, current_attenuator, line, line_number)

    if not builders:
        _raise_parse_error(coefficients_file, None, "no [effective YYYY-MM-DD] blocks found")
    return tuple(_validate_formula_block(builder) for builder in builders)


def _parse_formula_block(coefficients_file: Union[str, Path], allow_effective: bool = False) -> _FormulaBlock:
    """Parse a single attenuation formula block.

    Parameters
    ----------
    coefficients_file
        Path to the file containing one formula block.
    allow_effective
        Whether to accept an ``[effective YYYY-MM-DD]`` header in the block.

    Returns
    -------
    _FormulaBlock
        Validated formula block.

    Raises
    ------
    ValueError
        If the file is malformed or contains an effective header when ``allow_effective`` is ``False``.
    """
    builder = _FormulaBlockBuilder(source=coefficients_file)
    current_attenuator = None

    with open(coefficients_file, "r") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            effective_date = _parse_effective_header(line, coefficients_file, line_number)
            if effective_date is not None:
                if not allow_effective:
                    _raise_parse_error(
                        coefficients_file,
                        line_number,
                        "custom attenuation coefficients files must not contain an [effective YYYY-MM-DD] section",
                    )
                if builder.effective_date is not None:
                    _raise_parse_error(coefficients_file, line_number, "effective date is repeated")
                builder.effective_date = effective_date
                builder.effective_line_number = line_number
                current_attenuator = None
                continue

            current_attenuator = _parse_block_line(builder, current_attenuator, line, line_number)

    return _validate_formula_block(builder)


def _parse_block_line(
    builder: _FormulaBlockBuilder, current_attenuator: Optional[str], line: str, line_number: int
) -> Optional[str]:
    """Parse one meaningful line within a formula block.

    Parameters
    ----------
    builder
        Mutable parser state for the current formula block.
    current_attenuator
        Attenuator section currently receiving parameter lines, or ``None`` before the first section.
    line
        Stripped non-comment input line.
    line_number
        One-based source line number.

    Returns
    -------
    str or None
        Updated current attenuator name.

    Raises
    ------
    ValueError
        If the line is not valid in its current parser context.
    """
    if line.startswith("[") and line.endswith("]"):
        return _parse_attenuator_block(builder, line, line_number)
    if line.lower().startswith("formula"):
        if current_attenuator is not None:
            _raise_parse_error(builder.source, line_number, "formula must be declared before attenuator blocks")
        name, value = _split_assignment(builder.source, line, line_number)
        if name != "formula":
            _raise_parse_error(builder.source, line_number, "expected formula = <expression>")
        if builder.formula is not None:
            _raise_parse_error(builder.source, line_number, "formula is repeated")
        if not value:
            _raise_parse_error(builder.source, line_number, "formula is empty")
        builder.formula = value
        builder.formula_line_number = line_number
        return current_attenuator
    if current_attenuator is None:
        if "," in line:
            _raise_parse_error(
                builder.source,
                line_number,
                "old comma-separated attenuation coefficients are no longer supported; "
                "use formula = ... and [attenuator name] parameter blocks",
            )
        _raise_parse_error(builder.source, line_number, "expected formula or [attenuator name] block")

    parameter, parameter_value = _parse_parameter_line(builder.source, line, line_number)
    parameters = builder.coefficients[current_attenuator]
    if parameter in parameters:
        _raise_parse_error(builder.source, line_number, f"parameter {parameter} is repeated for {current_attenuator}")
    parameters[parameter] = parameter_value
    return current_attenuator


def _parse_attenuator_block(builder: _FormulaBlockBuilder, line: str, line_number: int) -> str:
    """Parse an attenuator section header.

    Parameters
    ----------
    builder
        Mutable parser state for the current formula block.
    line
        Section header line, expected to have the form ``[attenuator name]``.
    line_number
        One-based source line number.

    Returns
    -------
    str
        Attenuator name declared by the section.

    Raises
    ------
    ValueError
        If the section is malformed or repeats an attenuator name.
    """
    header = line[1:-1].strip()
    prefix = "attenuator "
    if not header.lower().startswith(prefix):
        _raise_parse_error(builder.source, line_number, f"unknown section [{header}]")
    name = header[len(prefix) :].strip()
    if not name:
        _raise_parse_error(builder.source, line_number, "empty attenuator name")
    if name in builder.coefficients:
        _raise_parse_error(builder.source, line_number, f"attenuator {name} is repeated")
    builder.coefficients[name] = {}
    builder.attenuator_line_numbers[name] = line_number
    return name


def _parse_parameter_line(source: Union[str, Path], line: str, line_number: int) -> Tuple[str, _ParameterValue]:
    """Parse one fitted-parameter declaration.

    Parameters
    ----------
    source
        Path or label of the file being parsed.
    line
        Parameter declaration with the form ``parameter = value, uncertainty``.
    line_number
        One-based source line number.

    Returns
    -------
    tuple
        Parameter name and parsed value/error record.

    Raises
    ------
    ValueError
        If the parameter name is invalid, reserved, duplicated by the caller, or has non-finite values.
    """
    parameter, values = _split_assignment(source, line, line_number)
    if not parameter:
        _raise_parse_error(source, line_number, "empty parameter name")
    if not parameter.isidentifier() or keyword.iskeyword(parameter):
        _raise_parse_error(source, line_number, f"{parameter} is not a valid parameter name")
    if parameter in _RESERVED_NAMES:
        _raise_parse_error(source, line_number, f"{parameter} is a reserved name")
    fields = [field.strip() for field in values.split(",")]
    if len(fields) != 2:
        _raise_parse_error(source, line_number, "parameter lines must have the form parameter = value, uncertainty")
    try:
        value = float(fields[0])
        error = float(fields[1])
    except ValueError as error:
        _raise_parse_error(source, line_number, str(error))
    if not math.isfinite(value) or not math.isfinite(error):
        _raise_parse_error(source, line_number, "parameter values and uncertainties must be finite numbers")
    return parameter, _ParameterValue(value=value, error=error, line_number=line_number)


def _validate_formula_block(builder: _FormulaBlockBuilder) -> _FormulaBlock:
    """Validate parser state and compile it into a formula block.

    Parameters
    ----------
    builder
        Raw parser state for one formula block.

    Returns
    -------
    _FormulaBlock
        Validated formula block with compiled value and derivative callables.

    Raises
    ------
    ValueError
        If the block is missing required content, attenuator parameter sets differ, or the formula and declared
        parameters are inconsistent.
    """
    if builder.formula is None:
        _raise_parse_error(builder.source, None, "attenuation coefficients block is missing a formula")
    if not builder.coefficients:
        _raise_parse_error(
            builder.source, builder.formula_line_number, "attenuation coefficients block has no attenuators"
        )

    first_attenuator = next(iter(builder.coefficients))
    parameters = tuple(builder.coefficients[first_attenuator])
    if not parameters:
        _raise_parse_error(
            builder.source,
            builder.attenuator_line_numbers[first_attenuator],
            f"attenuator {first_attenuator} has no parameters",
        )

    parameter_set = set(parameters)
    for attenuator, attenuator_parameters in builder.coefficients.items():
        current_set = set(attenuator_parameters)
        if current_set != parameter_set:
            missing = sorted(parameter_set - current_set)
            extra = sorted(current_set - parameter_set)
            details = []
            if missing:
                details.append(f"missing {', '.join(missing)}")
            if extra:
                details.append(f"extra {', '.join(extra)}")
            _raise_parse_error(
                builder.source,
                builder.attenuator_line_numbers[attenuator],
                f"attenuator {attenuator} does not define the same parameters as {first_attenuator}: "
                + "; ".join(details),
            )

    compiled = _compile_formula(builder.formula, parameters, builder.source, builder.formula_line_number)
    formula_parameters = _formula_parameter_names(compiled.expression, parameters)
    missing = sorted(formula_parameters - parameter_set)
    extra = sorted(parameter_set - formula_parameters)
    if missing:
        _raise_parse_error(
            builder.source, builder.formula_line_number, f"formula parameters missing from block: {missing}"
        )
    if extra:
        _raise_parse_error(builder.source, builder.formula_line_number, f"parameters not used by formula: {extra}")

    return _FormulaBlock(
        formula=builder.formula,
        parameters=parameters,
        coefficients=builder.coefficients,
        compiled=compiled,
        effective_date=builder.effective_date,
    )


def _compile_formula(
    formula: str,
    parameters: Tuple[str, ...],
    source: Optional[Union[str, Path]] = None,
    line_number: Optional[int] = None,
) -> _CompiledFormula:
    """Compile a validated attenuation formula string.

    Parameters
    ----------
    formula
        Formula expression from the coefficients file.
    parameters
        Ordered fitted parameter names declared by the first attenuator block.
    source
        Optional file path or label for diagnostics.
    line_number
        Optional source line number for diagnostics.

    Returns
    -------
    _CompiledFormula
        SymPy expression and numerical callables for the formula and its partial derivatives.

    Raises
    ------
    ValueError
        If the expression uses unsupported syntax or names not declared by the attenuator blocks.
    """
    try:
        parsed = ast.parse(formula, mode="eval")
    except SyntaxError as error:
        _raise_parse_error(source, line_number, f"invalid formula syntax: {error.msg}")

    symbols = {parameter: sympy.Symbol(parameter, real=True) for parameter in parameters}
    wavelength_symbol = sympy.Symbol("wavelength", real=True)
    expression = _ast_to_sympy(parsed.body, symbols, wavelength_symbol, source, line_number)
    uses_wavelength = wavelength_symbol in expression.free_symbols
    unknown_symbols = {str(symbol) for symbol in expression.free_symbols} - set(parameters) - {"wavelength"}
    if unknown_symbols:
        _raise_parse_error(source, line_number, f"unknown formula symbols: {sorted(unknown_symbols)}")

    arguments = [symbols[parameter] for parameter in parameters] + [wavelength_symbol]
    value_function = sympy.lambdify(arguments, expression, modules="numpy")
    derivative_functions = {
        parameter: sympy.lambdify(arguments, sympy.diff(expression, symbols[parameter]), modules="numpy")
        for parameter in parameters
    }
    return _CompiledFormula(
        expression=expression,
        parameters=parameters,
        uses_wavelength=uses_wavelength,
        value_function=value_function,
        derivative_functions=derivative_functions,
    )


def _evaluate_formula_with_error(
    formula_block: _FormulaBlock, attenuator_name: str, wavelength: float
) -> Tuple[float, float]:
    """Evaluate an attenuator formula and its propagated uncertainty.

    Parameters
    ----------
    formula_block
        Validated formula block containing coefficients and derivative callables.
    attenuator_name
        Name of the attenuator to evaluate.
    wavelength
        Wavelength in Å. Ignored by formulas that do not use ``wavelength``.

    Returns
    -------
    tuple of float
        Formula value and uncertainty from independent-parameter propagation.
    """
    parameter_values = formula_block.coefficients[attenuator_name]
    values = [parameter_values[parameter].value for parameter in formula_block.parameters]
    arguments = values + [wavelength]
    value = float(formula_block.compiled.value_function(*arguments))
    variance = 0.0
    for parameter in formula_block.parameters:
        derivative = float(formula_block.compiled.derivative_functions[parameter](*arguments))
        variance += (derivative * parameter_values[parameter].error) ** 2
    return value, float(np.sqrt(variance))


def _select_effective_block(blocks: Tuple[_FormulaBlock, ...], run_start: Union[str, datetime]) -> _FormulaBlock:
    """Select the default block effective for a run timestamp.

    Parameters
    ----------
    blocks
        Effective-dated formula blocks parsed from the packaged coefficients file.
    run_start
        Run timestamp to compare against block effective dates.

    Returns
    -------
    _FormulaBlock
        Most recent block whose effective date is not later than the run date.

    Raises
    ------
    ValueError
        If no block is effective for ``run_start``.
    """
    if not blocks:
        raise ValueError("No attenuation coefficient blocks are available")
    run_date = _as_instrument_date(run_start)
    sorted_blocks = sorted(blocks, key=lambda block: block.effective_date)
    selected = None
    for block in sorted_blocks:
        if block.effective_date <= run_date:
            selected = block
        else:
            break
    if selected is None:
        earliest = sorted_blocks[0].effective_date.isoformat()
        message = (
            f"Run timestamp {run_date.isoformat()} is earlier than the earliest attenuation calibration {earliest}"
        )
        logger.error(message)
        raise ValueError(message)
    return selected


def _run_start_time(input_workspace: Union[str, MatrixWorkspace]) -> datetime:
    """Read the run timestamp used for default attenuation calibration selection.

    Parameters
    ----------
    input_workspace
        Workspace whose sample logs contain ``start_time``, ``run_start`` or ``run_begin``.

    Returns
    -------
    datetime
        Run timestamp. Naive timestamps are interpreted later as instrument-local time.

    Raises
    ------
    RuntimeError
        If none of ``start_time``, ``run_start`` or ``run_begin`` is present.
    """
    logs = SampleLogs(input_workspace)
    for log_name in ("start_time", "run_start", "run_begin"):
        try:
            return _as_datetime(logs.single_value(log_name))
        except RuntimeError:
            continue
    message = (
        "Workspace must contain a start_time, run_start or run_begin sample log to select attenuation coefficients"
    )
    logger.error(message)
    raise RuntimeError(message)


def _parse_effective_header(line: str, source: Union[str, Path], line_number: int) -> Optional[date]:
    """Parse an optional default-file effective-date header.

    Parameters
    ----------
    line
        Stripped input line.
    source
        Path or label of the file being parsed.
    line_number
        One-based source line number.

    Returns
    -------
    datetime.date or None
        Effective date when ``line`` is an ``[effective YYYY-MM-DD]`` header; otherwise ``None``.

    Raises
    ------
    ValueError
        If the line is an effective header with an invalid date.
    """
    if not (line.startswith("[") and line.endswith("]")):
        return None
    header = line[1:-1].strip()
    prefix = "effective "
    if not header.lower().startswith(prefix):
        return None
    effective_text = header[len(prefix) :].strip()
    try:
        effective = date.fromisoformat(effective_text)
    except ValueError:
        _raise_parse_error(source, line_number, f"invalid effective date {effective_text!r}")
    return effective


def _split_assignment(source: Optional[Union[str, Path]], line: str, line_number: Optional[int]) -> Tuple[str, str]:
    """Split a configuration assignment into stripped name and value text.

    Parameters
    ----------
    source
        Optional file path or label for diagnostics.
    line
        Input line expected to contain ``=``.
    line_number
        Optional one-based source line number.

    Returns
    -------
    tuple of str
        Assignment name and value.

    Raises
    ------
    ValueError
        If ``line`` does not contain an assignment separator.
    """
    if "=" not in line:
        _raise_parse_error(source, line_number, "expected assignment containing '='")
    name, value = line.split("=", 1)
    return name.strip(), value.strip()


def _ast_to_sympy(
    node: ast.AST,
    symbols: Dict[str, sympy.Symbol],
    wavelength_symbol: sympy.Symbol,
    source: Optional[Union[str, Path]],
    line_number: Optional[int],
) -> sympy.Expr:
    """Translate an allowed Python expression AST into a SymPy expression.

    Parameters
    ----------
    node
        AST node from ``ast.parse(..., mode="eval")``.
    symbols
        Allowed fitted-parameter names mapped to SymPy symbols.
    wavelength_symbol
        SymPy symbol for the reserved ``wavelength`` variable.
    source
        Optional file path or label for diagnostics.
    line_number
        Optional source line number for diagnostics.

    Returns
    -------
    sympy.Expr
        Symbolic expression equivalent to ``node``.

    Raises
    ------
    ValueError
        If ``node`` contains unsupported syntax, unsupported function calls, or undeclared names.
    """
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return sympy.Float(node.value) if isinstance(node.value, float) else sympy.Integer(node.value)
        _raise_parse_error(source, line_number, "formula constants must be numbers")
    if isinstance(node, ast.Name):
        if node.id in symbols:
            return symbols[node.id]
        if node.id == "wavelength":
            return wavelength_symbol
        if node.id in _SUPPORTED_CONSTANTS:
            return _SUPPORTED_CONSTANTS[node.id]
        _raise_parse_error(source, line_number, f"{node.id} is not a parameter of any attenuator block")
    if isinstance(node, ast.BinOp):
        left = _ast_to_sympy(node.left, symbols, wavelength_symbol, source, line_number)
        right = _ast_to_sympy(node.right, symbols, wavelength_symbol, source, line_number)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.Pow):
            return left**right
        _raise_parse_error(source, line_number, "unsupported formula operator")
    if isinstance(node, ast.UnaryOp):
        operand = _ast_to_sympy(node.operand, symbols, wavelength_symbol, source, line_number)
        if isinstance(node.op, ast.UAdd):
            return operand
        if isinstance(node.op, ast.USub):
            return -operand
        _raise_parse_error(source, line_number, "unsupported unary formula operator")
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in _SUPPORTED_FUNCTIONS:
            _raise_parse_error(source, line_number, "formula calls are limited to supported functions")
        if node.keywords:
            _raise_parse_error(source, line_number, "formula functions do not accept keyword arguments")
        arguments = [
            _ast_to_sympy(argument, symbols, wavelength_symbol, source, line_number) for argument in node.args
        ]
        try:
            return _SUPPORTED_FUNCTIONS[node.func.id](*arguments)
        except (TypeError, ValueError) as error:
            _raise_parse_error(source, line_number, f"invalid call to {node.func.id}: {error}")
    _raise_parse_error(source, line_number, "formula contains unsupported syntax")


def _formula_parameter_names(expression: sympy.Expr, parameters: Tuple[str, ...]) -> set:
    """Return declared parameter names that appear in a compiled expression.

    Parameters
    ----------
    expression
        SymPy expression compiled from the formula string.
    parameters
        Declared fitted parameter names.

    Returns
    -------
    set
        Subset of ``parameters`` referenced by ``expression``.
    """
    free_symbol_names = {symbol.name for symbol in expression.free_symbols}
    return set(parameters) & free_symbol_names


def _as_instrument_date(timestamp: Union[str, datetime]) -> date:
    """Convert a timestamp to an instrument-local calendar date.

    Parameters
    ----------
    timestamp
        String or datetime timestamp from a workspace sample log.

    Returns
    -------
    datetime.date
        Date in the GPSANS instrument time zone.
    """
    return _as_datetime(timestamp).astimezone(_INSTRUMENT_TIME_ZONE).date()


def _as_datetime(timestamp: Union[str, datetime]) -> datetime:
    """Normalize a timestamp string or datetime object.

    Parameters
    ----------
    timestamp
        Timestamp value to normalize.

    Returns
    -------
    datetime.datetime
        Timezone-aware timestamp. Naive timestamps are interpreted in the GPSANS instrument time zone.
    """
    if isinstance(timestamp, datetime):
        result = timestamp
    else:
        result = parse_date(str(timestamp))
    if result.tzinfo is None:
        result = result.replace(tzinfo=_INSTRUMENT_TIME_ZONE)
    return result


def _raise_parse_error(source: Optional[Union[str, Path]], line_number: Optional[int], message: str) -> None:
    """Log and raise a formatted attenuation coefficients parse error.

    Parameters
    ----------
    source
        Optional file path or label for diagnostics.
    line_number
        Optional one-based source line number.
    message
        Human-readable validation failure.

    Raises
    ------
    ValueError
        Always raised with the formatted message.
    """
    location = f"line {line_number} in " if line_number is not None else ""
    source_text = str(source) if source is not None else "attenuation formula"
    full_message = f"Invalid {location}attenuation coefficients file {source_text}: {message}"
    logger.error(full_message)
    raise ValueError(full_message)


def _attenuation_factor(A, A_e, B, B_e, C, C_e, wavelength):
    """Evaluate the legacy GPSANS attenuation formula.

    Parameters
    ----------
    A, B, C
        Legacy fitted coefficients for ``A * exp(-B * wavelength) + C``.
    A_e, B_e, C_e
        One-sigma uncertainties for ``A``, ``B`` and ``C``.
    wavelength
        Wavelength in Å.

    Returns
    -------
    tuple of float
        Formula value and propagated uncertainty.

    Notes
    -----
    This helper is kept for tests and regression calculations while runtime evaluation uses formula blocks.
    """
    scale = A * np.exp(-B * wavelength) + C
    scale_error_amp = np.exp(-B * wavelength)
    scale_error_exp_const = A * np.exp(-B * wavelength) * (-wavelength)
    scale_error_bkgd = 1
    scale_error = np.sqrt(
        (scale_error_amp * A_e) ** 2 + (scale_error_exp_const * B_e) ** 2 + (scale_error_bkgd * C_e) ** 2
    )
    return scale, scale_error
