import importlib.resources
import math
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from mantid.api import MatrixWorkspace
from mantid.kernel import logger
import numpy as np

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py
from drtsans.samplelogs import SampleLogs

# Functions exposed to the general user (public) API
__all__ = ["attenuation_factor"]

# Value of the "attenuator" sample log mapped to the attenuator name used in the coefficients file
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

# Attenuator names for which the beam is not attenuated
_NO_ATTENUATION = {"Undefined", "Close", "Open"}

_DEFAULT_COEFFICIENTS_FILE = "GPSANS_attenuation_coefficients.txt"


def attenuation_factor(
    input_workspace: Union[str, MatrixWorkspace], coefficients_file: Optional[Union[str, Path]] = None
) -> Tuple[float, float]:
    r"""This gets the wavelength and attenuator value from the workspace
    logs then calculates the attenuation factor based on the fitted
    parameters for each different attenuator based on the equation

    .. math::

       attenuation factor = A e^{-B \lambda} + C

    with the wavelength :math:`\lambda` in Å, :math:`B` in 1/Å, and :math:`A`, :math:`C` dimensionless.
    The fitted parameters are read from ``coefficients_file``. See :ref:`user.corrections.attenuation`
    for the file format.

    The attenuation scale factor and the uncertainty for this is
    returned.

    The attenuator value pulled from the logs is mapped to the
    attenuator name by:

    0: "Undefined"
    1: "Close"
    2: "Open",
    3: "x3",
    4: "x30",
    5: "x300",
    6: "x2k",
    7: "x10k",
    8: "x100k"

    If the attenuator is one of Undefined, Close or Open then a
    attenuation factor of 1 with uncertainty 0 is returned. A negative log value, as found in runs converted from
    SPICE files, is taken as Undefined, with a warning.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Input workspace for which to calculate attenuation factor
    coefficients_file: str, ~pathlib.Path
        Path to the file of attenuator fit coefficients. If :py:obj:`None`, the file
        ``GPSANS_attenuation_coefficients.txt`` in the ``drtsans.configuration`` package is used.

    Returns
    -------
    float, float
        attenuation_factor, attenuation_factor_uncertainty

    Raises
    ------
    RuntimeError
        If the workspace has no "attenuator" sample log, or no "wavelength" sample log when the beam is attenuated
    FileNotFoundError
        If the coefficients file does not exist
    ValueError
        If the attenuator log value is unknown, if the attenuator is missing from the coefficients file,
        or if the coefficients file is malformed.
    """
    attenuator_name = _attenuator_name(input_workspace)

    if attenuator_name in _NO_ATTENUATION:
        # return scale factor of 1 and error 0, without reading the coefficients file
        return 1, 0

    coefficients = _load_attenuation_coefficients(coefficients_file)
    return _attenuator_transmission(input_workspace, attenuator_name, coefficients, coefficients_file)


def _attenuator_transmission(
    input_workspace: Union[str, MatrixWorkspace],
    attenuator_name: str,
    coefficients: Dict[str, Tuple[float, ...]],
    coefficients_file: Optional[Union[str, Path]] = None,
) -> Tuple[float, float]:
    """Transmitted fraction of an attenuator, and its uncertainty, at the wavelength of the workspace.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace containing the "wavelength" sample log (Å)
    attenuator_name
        Attenuator name, as returned by :py:func:`_attenuator_name`
    coefficients
        Attenuation coefficients, as returned by :py:func:`_load_attenuation_coefficients`
    coefficients_file
        Path to the file the coefficients were read from, only used in error messages. If :py:obj:`None`,
        the default file is assumed.

    Returns
    -------
    float, float
        Transmitted fraction and its uncertainty. (1, 0) for the attenuators "Undefined", "Close" and "Open".

    Raises
    ------
    ValueError
        If the attenuator is missing from the coefficients
    """
    if attenuator_name in _NO_ATTENUATION:
        return 1, 0

    if attenuator_name not in coefficients:
        if coefficients_file is None:
            coefficients_file = _DEFAULT_COEFFICIENTS_FILE
        message = f"Attenuator {attenuator_name} not found in the attenuation coefficients file {coefficients_file}"
        logger.error(message)
        raise ValueError(message)

    wavelength = SampleLogs(input_workspace).single_value("wavelength")
    return _attenuation_factor(*coefficients[attenuator_name], wavelength)


def _attenuator_name(input_workspace: Union[str, MatrixWorkspace]) -> str:
    """Name of the attenuator given by the "attenuator" sample log of the workspace.

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Workspace containing the "attenuator" sample log

    Returns
    -------
    str
        One of "Undefined", "Close", "Open", "x3", "x30", "x300", "x2k", "x10k", "x100k". A negative log value
        is taken as "Undefined", with a warning.

    Raises
    ------
    RuntimeError
        If the workspace has no "attenuator" sample log
    ValueError
        If the log value is positive and not an integer from 0 to 8, for instance when the attenuator changed
        during the run
    """
    attenuator = SampleLogs(input_workspace).single_value("attenuator")
    if attenuator < 0:
        # Runs converted from SPICE files store the attenuator stage position (mm) instead of the attenuator index
        logger.warning(
            f"Negative attenuator log value {attenuator}, probably an attenuator stage position (mm) "
            "from a SPICE file. The attenuator is taken as Undefined and no attenuation correction is applied"
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
) -> Dict[str, Tuple[float, ...]]:
    """Read the attenuator fit coefficients from a comma-separated file.

    Lines that are blank or start with ``#`` are ignored. Every other line must contain the attenuator name
    followed by six finite numbers: A, A error, B, B error, C, C error. Each attenuator name must be unique.

    Parameters
    ----------
    coefficients_file
        Path to the attenuation coefficients file. If :py:obj:`None`, the file
        ``GPSANS_attenuation_coefficients.txt`` in the ``drtsans.configuration`` package is used.

    Returns
    -------
    dict
        Attenuator name mapped to the tuple (A, A error, B, B error, C, C error), in the order of the file

    Raises
    ------
    FileNotFoundError
        If the coefficients file does not exist
    ValueError
        If a line does not have seven fields, if the attenuator name is empty or repeated, or if any of the last six
        fields is not a finite number
    """
    if coefficients_file is None:
        with importlib.resources.as_file(
            importlib.resources.files("drtsans.configuration") / _DEFAULT_COEFFICIENTS_FILE
        ) as default_file:
            return _load_attenuation_coefficients(default_file)

    coefficients = {}
    with open(coefficients_file, "r") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in line.split(",")]
            try:
                if len(fields) != 7:
                    raise ValueError(f"expected 7 comma-separated fields, found {len(fields)}")
                name = fields[0]
                if not name:
                    raise ValueError("empty attenuator name")
                if name in coefficients:
                    raise ValueError(f"attenuator {name} is repeated")
                values = tuple(float(field) for field in fields[1:])
                if not all(math.isfinite(value) for value in values):
                    raise ValueError("coefficients must be finite numbers")
                coefficients[name] = values
            except ValueError as error:
                message = f"Invalid line {line_number} in attenuation coefficients file {coefficients_file}: {error}"
                logger.error(message)
                raise ValueError(message) from error
    return coefficients


def _attenuation_factor(A, A_e, B, B_e, C, C_e, wavelength):
    """
    This calculates the function
        A * exp(-B * λ) + C
    along with the uncertainty
    """
    scale = A * np.exp(-B * wavelength) + C
    scale_error_Amp = np.exp(-B * wavelength)
    scale_error_exp_const = A * np.exp(-B * wavelength) * (-wavelength)
    scale_error_bkgd = 1
    scale_error = np.sqrt(
        (scale_error_Amp * A_e) ** 2 + (scale_error_exp_const * B_e) ** 2 + (scale_error_bkgd * C_e) ** 2
    )
    return scale, scale_error
