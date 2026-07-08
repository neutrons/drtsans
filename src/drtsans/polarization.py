# standard library imports
from collections import namedtuple
from enum import StrEnum
from dataclasses import dataclass
from math import fsum
from typing import ClassVar, Generator, List, Optional, Union

# third party imports
import h5py
from mantid.api import AnalysisDataService
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import CreateSingleValuedWorkspace, DeleteLog, DeleteWorkspace, logger, mtd, RenameWorkspace
import numpy as np
import sympy as sp

# drtsans imports
from drtsans.api import _set_uncertainty_from_numpy
from drtsans.path import abspath
from drtsans.samplelogs import SampleLogs
from drtsans.type_hints import MantidWorkspace


class PolarizationLevel(StrEnum):
    NONE = "none"  # also 0
    HALF = "half"  # also 1
    FULL = "full"  # also 2

    @classmethod
    def from_int(cls, level: int) -> "PolarizationLevel":
        r"""
        Convert an integer polarization mode to a PolarizationLevel enum.

        Parameters
        ----------
        level : int
            Integer representation of the polarization mode:
            - 0: OFF (unpolarized)
            - 1: HALF (polarizer is active)
            - 2: FULL (both polarizer and analyzer are active)

        Raises
        ------
        ValueError
            If the input integer does not correspond to a valid polarization mode.
        """
        if level == 0:
            return cls.NONE
        elif level == 1:
            return cls.HALF
        elif level == 2:
            return cls.FULL
        else:
            raise ValueError(f"Invalid polarization mode integer: {level}. Must be 0, 1, or 2.")

    @classmethod
    def get(cls, source: Union[str, EventWorkspace]) -> "PolarizationLevel":
        r"""
        Find if a run is polarized, and to what degree.

        Parameters
        ----------
        source : str, EventWorkspace
            Either an instrument plus run number identifier (e.g., CG2_1235), a file name in the current directory,
            an absolute file path, the name of an events workspace, or an EventWorkspace object.

        Raises
        ------
        TypeError
            If the input is neither a string nor an EventWorkspace.
        """
        mode = 0
        # Determines polarization mode from workspace or file metadata
        if isinstance(source, EventWorkspace):
            sample_logs = SampleLogs(source)
            if PV_POLARIZER in sample_logs:
                mode += int(sample_logs.single_value(PV_POLARIZER))
                if PV_ANALYZER in sample_logs:
                    mode += int(sample_logs.single_value(PV_ANALYZER))
        elif isinstance(source, str):
            # case 1: `source` is the name of a workspace in the AnalysisDataService
            if AnalysisDataService.doesExist(source):
                if isinstance(mtd[source], EventWorkspace):
                    sample_logs = SampleLogs(source)
                    if PV_POLARIZER in sample_logs:
                        mode += int(sample_logs.single_value(PV_POLARIZER))
                        if PV_ANALYZER in sample_logs:
                            mode += int(sample_logs.single_value(PV_ANALYZER))
                else:
                    raise TypeError(f"The workspace '{source}' is not an EventWorkspace")
            else:
                # case 2: `source` is a file path
                filepath = abspath(source)
                with h5py.File(filepath, "r") as file_handle:
                    dataset = file_handle.get(f"/entry/DASlogs/{PV_POLARIZER}")
                    if dataset is not None:  # log not found, assume unpolarized
                        mode += int(dataset["value"][()])
                        dataset = file_handle.get(f"/entry/DASlogs/{PV_ANALYZER}")
                        if dataset is not None:
                            mode += int(dataset["value"][()])
        else:
            raise TypeError(f"{source} must be either a string or an EventWorkspace object")
        return cls.from_int(mode)


class PolarizationCrossSection(StrEnum):
    """Enumerate the possible spin cross-section states based on flipper and analyzer status."""

    NONE = "none"  # polarizer and analyzer disengaged
    OFF = "off"  # polarizer flipper off, analyzer disengaged
    ON = "on"  # polarizer flipper on, analyzer disengaged
    OFF_OFF = "off_off"  # polarizer flipper off, analyzer at Zero state
    OFF_ON = "off_on"  # polarizer flipper off, analyzer at Pi state
    ON_OFF = "on_off"  # polarizer flipper on, analyzer at Zero state
    ON_ON = "on_on"  # polarizer flipper on, analyzer at Pi state

    @classmethod
    def get(cls, workspace: MantidWorkspace) -> "PolarizationCrossSection":
        """Retrieve the polarization cross-section from the sample logs of the given workspace."""
        return cls(str(SampleLogs(workspace).single_value(cls.logname)))

    def log(self, workspace: MantidWorkspace):
        """Insert the polarization cross-section into the sample logs of the given workspace."""
        SampleLogs(workspace).insert(name=self.__class__.logname, value=self.value)

    @property
    def level(self) -> PolarizationLevel:
        if self == PolarizationCrossSection.NONE:
            return PolarizationLevel.NONE
        if self in {PolarizationCrossSection.OFF, PolarizationCrossSection.ON}:
            return PolarizationLevel.HALF
        return PolarizationLevel.FULL


# Add class variable after the enum definition
PolarizationCrossSection.logname = "polarization: cross-section"


class PolarizationState(StrEnum):
    """Enumerate the possible spin states based on the upstream and downstream neutrons"""

    NONE = "none"  # no polarization
    UP = "up"  # half polarization, upstream spin up
    DOWN = "down"  # half polarization, upstream spin down
    UP_UP = "up_up"  # full polarization, upstream spin up, downstream spin up
    UP_DOWN = "up_down"  # full polarization, upstream spin up, downstream spin down
    DOWN_UP = "down_up"  # full polarization, upstream spin down, downstream spin up
    DOWN_DOWN = "down_down"  # full polarization, upstream spin down, downstream spin down

    @classmethod
    def get(cls, workspace: MantidWorkspace) -> "PolarizationState":
        """Retrieve the polarization state from the sample logs of the given workspace."""
        return cls(str(SampleLogs(workspace).single_value(cls.logname)))

    def log(self, workspace: MantidWorkspace):
        """Insert the polarization state into the sample logs of the given workspace."""
        SampleLogs(workspace).insert(name=self.__class__.logname, value=self.value)

    @property
    def level(self) -> PolarizationLevel:
        if self == PolarizationState.NONE:
            return PolarizationLevel.NONE
        if self in {PolarizationState.UP, PolarizationState.DOWN}:
            return PolarizationLevel.HALF
        return PolarizationLevel.FULL


# Add class variable after the enum definition
PolarizationState.logname = "polarization: state"


__all__ = [
    "PV_POLARIZER",
    "PV_POLARIZER_FLIPPER",
    "PV_POLARIZER_VETO",
    "PV_ANALYZER",
    "PV_ANALYZER_FLIPPER",
    "PV_ANALYZER_VETO",
    "half_polarization",
    "SimulatedPolarizationLogs",
]

# Names of processing variables related to polarization, stored in the sample logs of the Nexus Event file
# For the moment, we assume these PVs are consistent across instruments.
PV_POLARIZER = "Polarizer"
PV_POLARIZER_FLIPPER = "PolarizerFlipper"
PV_POLARIZER_VETO = "PolarizerVeto"
PV_ANALYZER = "Analyzer"
PV_ANALYZER_FLIPPER = "AnalyzerFlipper"
PV_ANALYZER_VETO = "AnalyzerVeto"


def polarized_sample(reduction_parameters: dict) -> bool:
    r"""
    Determine if the sample run involves polarized neutrons.

    This function checks for polarization under `configuration['polarization']['level']`.
    If the settings are missing, it examines the sample metadata to determine the polarization level.

    Parameters
    ----------
    reduction_parameters : dict
        Dictionary of reduction configuration parameters. It can be either the full reduction input containing
        the "configuration" key or just the reduction configuration itself.

    Returns
    -------
    bool
        True if the sample is polarized (polarization level is not NONE), False otherwise.

    Raises
    ------
    ValueError
        If multiple sample runs are provided and the first sample represents a polarized run.
        Currently, we can't reduce polarization for summed datasets.

    Notes
    -----
    The function modifies the reduction configuration by setting the polarization level in
    `configuration['polarization']['level']` if not previously specified.
    """

    if "configuration" in reduction_parameters:
        reduction_config = reduction_parameters["configuration"]
        directories = reduction_parameters.get("dataDirectories", None)
    else:
        reduction_config = reduction_parameters
        directories = None

    if reduction_config.get("polarization", {}).get("level", None) is None:
        if reduction_config == reduction_parameters:
            logger.warning("Unable to resolve polarization level. Setting to NONE by default.")
            reduction_config["polarization"] = {"level": str(PolarizationLevel.NONE)}
        else:
            sample = reduction_parameters["sample"]["runNumber"].strip()
            multiple_samples = len(sample.split(",")) > 1
            if multiple_samples:
                sample = sample.split(",")[0]  # inquire from the first sample
            sample_filepath = abspath(
                sample,
                instrument=reduction_parameters["instrumentName"],
                ipts=reduction_parameters["iptsNumber"],
                directory=directories,
                search_archive=True,
            )
            level = PolarizationLevel.get(sample_filepath)
            if multiple_samples and level != PolarizationLevel.NONE:
                raise ValueError("Can't do polarization reduction on summed data sets")
            reduction_config["polarization"] = {"level": str(level)}

    return reduction_config["polarization"]["level"] != PolarizationLevel.NONE


def _calc_flipping_ratio(polarization):
    """Calculates the flipping ratio from the polarization state

    Parameters
    ----------
    polarization: str, ~mantid.api.MatrixWorkspace
        Polarization state

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The ratio of flipping
    """
    value = polarization.extractY()[0]
    uncertainty = polarization.extractE()[0]
    if len(value) == 1 and len(uncertainty) == 1:
        value, uncertainty = value[0], uncertainty[0]

    uncertainty = 2.0 * uncertainty / np.square(1 - value)
    value = (1 + value) / (1 - value)

    if isinstance(value, float):  # create a single value workspace
        return CreateSingleValuedWorkspace(
            DataValue=value,
            ErrorValue=uncertainty,
            OutputWorkspace=mtd.unique_hidden_name(),
            EnableLogging=False,
        )
    else:
        # if it is an array should call CreateWorkspace(EnableLogging=False)
        raise NotImplementedError("Somebody needs to create an output from {} (type={})".format(value, type(value)))


def _calc_half_polarization_up(flipper_off, flipper_on, efficiency, flipping_ratio):
    """This calculates the spin up workspace

    Parameters
    ----------
    flipper_off_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper off measurement
    flipper_on_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper on measurement
    efficiency: ~mantid.api.MatrixWorkspace
        Flipper efficiency
    flipping_ratio: ~mantid.api.MatrixWorkspace
        The ratio of flipping

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The spin up workspace
    """
    __spin_up = flipper_off + (flipper_off - flipper_on) / (efficiency * (flipping_ratio - 1.0))

    e = efficiency.extractY()[0]
    F = flipping_ratio.extractY()[0]

    # the uncertainties (numerically) aren't correct because of build-up of numerical errors.
    # Recalculate them based on the proper equation based on a hand calculation of the partial derivatives
    m0_part = np.square(flipper_off.extractE()[0][0] * (1 + 1 / (e * (F - 1))))
    m1_part = np.square(flipper_on.extractE()[0][0] * (1 / (e * (F - 1))))

    mixed = np.square((flipper_off.extractY()[0][0] - flipper_on.extractY()[0][0]) / (e * (F - 1)))
    mixed *= np.square(efficiency.extractE()[0][0] / e) + np.square(flipping_ratio.extractE()[0][0] / (F - 1))
    sup_err = np.sqrt(m0_part + m1_part + mixed)

    # set the uncertainty in the workspace
    __spin_up = _set_uncertainty_from_numpy(__spin_up, sup_err)

    return __spin_up


def _calc_half_polarization_down(flipper_off, flipper_on, efficiency, flipping_ratio):
    """This calculates the spin down workspace

    Parameters
    ----------
    flipper_off_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper off measurement
    flipper_on_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper on measurement
    efficiency: ~mantid.api.MatrixWorkspace
        Flipper efficiency
    flipping_ratio: ~mantid.api.MatrixWorkspace
        The ratio of flipping

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The spin down workspace
    """
    __spin_down = flipper_off - (flipper_off - flipper_on) / (efficiency * (1.0 - 1.0 / flipping_ratio))

    e = efficiency.extractY()[0]
    F = flipping_ratio.extractY()[0]

    # the uncertainties (numerically) aren't correct because of build-up of numerical errors.
    # Recalculate them based on the proper equation based on a hand calculation of the partial derivatives
    m0_part = np.square(flipper_off.extractE()[0][0] * (1 - 1 / (e * (1 - 1 / F))))
    m1_part = np.square(flipper_on.extractE()[0][0] * (1 / (e * (1 - 1 / F))))
    mixed = np.square((flipper_off.extractY()[0][0] - flipper_on.extractY()[0][0]) / (e * (1 - 1 / F)))
    mixed *= np.square(efficiency.extractE()[0][0] / e) + np.square(
        flipping_ratio.extractE()[0][0] / (F * F * (1 - 1 / F))
    )

    sdn_err = np.sqrt(m0_part + m1_part + mixed)

    # set the uncertainty in the workspace
    __spin_down = _set_uncertainty_from_numpy(__spin_down, sdn_err)

    return __spin_down


def half_polarization(
    flipper_off_workspace,
    flipper_on_workspace,
    polarization,
    efficiency,
    spin_up_workspace=None,
    spin_down_workspace=None,
):
    """Calculate the spin up/down workspaces from flipper on/off.

    **Mantid algorithms used:**
    :ref:`RenameWorkspace <algm-RenameWorkspace-v1>`

    Parameters
    ----------
    flipper_off_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper off measurement
    flipper_on_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper on measurement
    polarization: str, ~mantid.api.MatrixWorkspace
        Polarization state
    efficiency: str, ~mantid.api.MatrixWorkspace
        Flipper efficiency
    spin_up_workspace: str
        Name of the resulting spin up workspace. If :py:obj:`None`, then
        ``flipper_off_workspace`` will be overwritten.
    spin_down_workspace: str
        Name of the resulting spin down workspace. If :py:obj:`None`, then
        ``flipper_on_workspace`` will be overwritten.

    Returns
    -------
    py:obj:`tuple` of 2 ~mantid.api.MatrixWorkspace
    """
    if spin_up_workspace is None:
        spin_up_workspace = str(flipper_off_workspace)
    if spin_down_workspace is None:
        spin_down_workspace = str(flipper_on_workspace)

    # this is denoted as "F" in the master document
    flipping_ratio = _calc_flipping_ratio(polarization)

    __spin_up = _calc_half_polarization_up(flipper_off_workspace, flipper_on_workspace, efficiency, flipping_ratio)
    __spin_down = _calc_half_polarization_down(flipper_off_workspace, flipper_on_workspace, efficiency, flipping_ratio)

    spin_up_workspace = RenameWorkspace(InputWorkspace=__spin_up, OutputWorkspace=spin_up_workspace)
    spin_down_workspace = RenameWorkspace(InputWorkspace=__spin_down, OutputWorkspace=spin_down_workspace)

    DeleteWorkspace(flipping_ratio)

    return (spin_up_workspace, spin_down_workspace)


class PolarizationDecoder:
    """
    Base class for converting measured device cross-sections into spin-state cross-sections.

    Parameters
    ----------
    reduction_config : dict
        Reduction configuration containing optional ``polarization.polarizer`` entries.
        The polarizer ``polarization`` and flipper ``efficiency`` values may be numeric
        constants or expressions in wavelength symbol ``x``.

    Attributes
    ----------
    p : callable
        Wavelength-dependent incident-beam polarization, bounded in ``[-1, 1]``.
    e : callable
        Wavelength-dependent polarizer flipper efficiency, bounded in ``[0, 1]``.
    """

    _x = sp.Symbol("x")  # wavelength symbol used in all sympy expressions

    @classmethod
    def _sympify(cls, value) -> sp.Expr:
        """Parse a config value (string or number) into a sympy expression in x (wavelength)."""
        return sp.sympify(str(value))

    @classmethod
    def _lambdify(cls, expr: sp.Expr):
        """Convert a sympy expression in x into a numpy-callable function."""
        return sp.lambdify(cls._x, expr, "numpy")

    @staticmethod
    def _validate_polarization(name: str, value: float):
        """
        Validate that a polarization value is physically meaningful.

        Parameters
        ----------
        name : str
            Human-readable name used in the exception message.
        value : float
            Value to validate.

        Raises
        ------
        ValueError
            If ``value`` is zero or not in the interval ``[-1, 1]``.
        """
        if value == 0 or not -1 <= value <= 1:
            raise ValueError(f"{name} must be non-zero and in the interval [-1, 1].")

    @staticmethod
    def _validate_efficiency(name: str, value: float):
        """
        Validate that an efficiency value is physically meaningful.

        Parameters
        ----------
        name : str
            Human-readable name used in the exception message.
        value : float
            Value to validate.

        Raises
        ------
        ValueError
            If ``value`` is not in the interval ``(0, 1]``.
        """
        if not 0 < value <= 1:
            raise ValueError(f"{name} must be in the interval (0, 1].")

    def __init__(self, reduction_config: dict):
        """
        Load polarizer properties from the reduction configuration.

        Parameters
        ----------
        reduction_config : dict
            Reduction configuration containing optional ``polarization.polarizer``
            entries for ``polarization`` and ``efficiency``. Missing values default
            to perfect polarization and perfect flipper efficiency.
        """
        polarizer = reduction_config.get("polarization", {}).get("polarizer", {})
        self.p = self._lambdify(self._sympify(polarizer.get("polarization", "1")))
        self.e = self._lambdify(self._sympify(polarizer.get("efficiency", "1")))

    def decode(self, device_cross_sections):
        """
        Decode measured device cross-sections into spin-state cross-sections.

        Parameters
        ----------
        device_cross_sections : list
            Workspaces corresponding to measured polarizer/analyzer device states.

        Raises
        ------
        NotImplementedError
            Always raised by the base class. Subclasses implement the instrument
            geometry for half or full polarization.
        """
        raise NotImplementedError("Subclasses must implement the decode method.")


class HalfPolarizationDecoder(PolarizationDecoder):
    """
    Decoder for half-polarization data.

    Half polarization uses two measured device cross-sections: polarizer flipper
    off (``S0``) and on (``S1``). The decoder applies the inverse of the 2 x 2
    encoding matrix from Section 9.1 of the master requirements document to
    recover spin-up and spin-down cross-sections.
    """

    def decoding_matrix(self, wavelength: float) -> np.ndarray:
        """
        Return the half-polarization decoding matrix at a wavelength.

        Parameters
        ----------
        wavelength : float
            Neutron wavelength in Angstrom.

        Returns
        -------
        numpy.ndarray
            A 2 x 2 matrix transforming measured ``[S0, S1]`` device
            cross-sections into spin-state ``[S_up, S_down]`` cross-sections.

        Raises
        ------
        ValueError
            If polarizer polarization is outside ``[-1, 1]`` or flipper
            efficiency is outside ``[0, 1]``.
        """
        p = self.p(wavelength)
        self._validate_polarization("Polarization", p)
        e = self.e(wavelength)
        self._validate_efficiency("Flipper efficiency", e)
        dl = (1 - p) / (2 * e * p)  # spin-down leakage
        ul = (1 + p) / (2 * e * p)  # spin-up leakage
        return np.array(
            [
                [1 + dl, -dl],
                [1 - ul, ul],
            ]
        )

    def decode(self, device_cross_sections: list[MantidWorkspace]):
        """
        Decode two half-polarization device cross-sections into spin states.

        Parameters
        ----------
        device_cross_sections : list of MantidWorkspace
            Workspaces tagged with ``PolarizationCrossSection.OFF`` and
            ``PolarizationCrossSection.ON`` sample logs.

        Returns
        -------
        list of MantidWorkspace
            Decoded spin-up and spin-down workspaces, tagged with
            ``PolarizationState.UP`` and ``PolarizationState.DOWN``.

        Raises
        ------
        ValueError
            If exactly two device cross-section workspaces are not supplied.
        """
        if len(device_cross_sections) != 2:
            raise ValueError(
                f"Half polarization requires exactly 2 device cross-sections, got {len(device_cross_sections)}."
            )

        # Find out which cross-section is S⁰ and which is S¹
        by_state = {PolarizationCrossSection.get(ws): ws for ws in device_cross_sections}
        s0 = by_state[PolarizationCrossSection.OFF]  # flipper off → S⁰
        s1 = by_state[PolarizationCrossSection.ON]  # flipper on  → S¹
        spinor = np.array([s0, s1], dtype=object)

        # Find the decoding matrix
        wavelength = SampleLogs(s0).single_value("wavelength")
        M = self.decoding_matrix(wavelength).astype(object)

        # Find the spin states from the device cross-sections, then inject the spin state as a sample-log
        spin_cross_sections = list(M @ spinor)
        for ws, state in zip(spin_cross_sections, [PolarizationState.UP, PolarizationState.DOWN]):
            DeleteLog(Workspace=ws, Name=PolarizationCrossSection.logname)  # delete the device cross-section samplelog
            state.log(ws)  # inject the spin state samplelog
        return spin_cross_sections


class FullPolarizationDecoder(PolarizationDecoder):
    """
    Decoder for full-polarization data.

    Full polarization uses four measured device cross-sections formed from the
    polarizer flipper state and the analyzer state. The decoder numerically
    inverts the wavelength-dependent 4 x 4 encoding matrix from Section 9.2 of
    the master requirements document to recover the four spin-state
    cross-sections.
    """

    def __init__(self, reduction_config: dict):
        """
        Load polarizer and analyzer properties from the reduction configuration.

        Parameters
        ----------
        reduction_config : dict
            Reduction configuration containing optional ``polarization.polarizer``
            and ``polarization.analyzer`` entries. Analyzer ``polarizationZero``
            and ``polarizationPi`` values may be numeric constants or expressions
            in wavelength symbol ``x``. Missing values default to 1 and -1,
            respectively.
        """
        super().__init__(reduction_config)
        analyzer = reduction_config.get("polarization", {}).get("analyzer", {})
        self.p_0 = self._lambdify(self._sympify(analyzer.get("polarizationZero", "1")))
        self.p_pi = self._lambdify(self._sympify(analyzer.get("polarizationPi", "-1")))

    def encoding_matrix(self, wavelength: float) -> np.ndarray:
        """
        Return the full-polarization encoding matrix at a wavelength.

        The matrix maps spin-state cross-sections ordered as
        ``[S_up_up, S_up_down, S_down_up, S_down_down]`` to measured device
        cross-sections ordered as ``[S00, S10, S0pi, S1pi]``.

        Parameters
        ----------
        wavelength : float
            Neutron wavelength in Angstrom.

        Returns
        -------
        numpy.ndarray
            A 4 x 4 encoding matrix built from polarizer polarization, flipper
            efficiency, and analyzer zero/pi-state polarizations.

        Raises
        ------
        ValueError
            If any polarization value is outside ``[-1, 1]`` or efficiency
            value is outside ``[0, 1]``.
        """
        p = self.p(wavelength)
        self._validate_polarization("Polarizer polarization", p)
        e = self.e(wavelength)
        self._validate_efficiency("Flipper efficiency", e)
        p_0 = self.p_0(wavelength)
        self._validate_polarization("Analyzer zero-state polarization", p_0)
        p_pi = self.p_pi(wavelength)
        self._validate_polarization("Analyzer pi-state polarization", p_pi)

        # Eq. 9.3 of the Master document converts signed polarization to the ratio terms used by Eqs. 9.16-9.19.
        # For the pi analyzer state, the Master ratio multiplies spin-down transmission terms,
        # so its ratio is reciprocal to the signed up/down polarization convention used here.
        polarizer_up_fraction = (1 + p) / 2
        polarizer_down_fraction = (1 - p) / 2

        flipper_on_up_fraction = e * polarizer_down_fraction + (1 - e) * polarizer_up_fraction
        flipper_on_down_fraction = e * polarizer_up_fraction + (1 - e) * polarizer_down_fraction

        analyzer_0_up_fraction = (1 + p_0) / 2
        analyzer_0_down_fraction = (1 - p_0) / 2
        analyzer_pi_up_fraction = (1 + p_pi) / 2
        analyzer_pi_down_fraction = (1 - p_pi) / 2

        return np.array(
            [
                [
                    polarizer_up_fraction * analyzer_0_up_fraction,
                    polarizer_up_fraction * analyzer_0_down_fraction,
                    polarizer_down_fraction * analyzer_0_up_fraction,
                    polarizer_down_fraction * analyzer_0_down_fraction,
                ],
                [
                    flipper_on_up_fraction * analyzer_0_up_fraction,
                    flipper_on_up_fraction * analyzer_0_down_fraction,
                    flipper_on_down_fraction * analyzer_0_up_fraction,
                    flipper_on_down_fraction * analyzer_0_down_fraction,
                ],
                [
                    polarizer_up_fraction * analyzer_pi_up_fraction,
                    polarizer_up_fraction * analyzer_pi_down_fraction,
                    polarizer_down_fraction * analyzer_pi_up_fraction,
                    polarizer_down_fraction * analyzer_pi_down_fraction,
                ],
                [
                    flipper_on_up_fraction * analyzer_pi_up_fraction,
                    flipper_on_up_fraction * analyzer_pi_down_fraction,
                    flipper_on_down_fraction * analyzer_pi_up_fraction,
                    flipper_on_down_fraction * analyzer_pi_down_fraction,
                ],
            ]
        )

    def decode(self, device_cross_sections: list[MantidWorkspace]):
        """
        Decode four full-polarization device cross-sections into spin states.

        Parameters
        ----------
        device_cross_sections : list of MantidWorkspace
            Workspaces tagged with ``OFF_OFF``, ``ON_OFF``, ``OFF_ON``, and
            ``ON_ON`` device cross-section sample logs.

        Returns
        -------
        list of MantidWorkspace
            Decoded spin-state workspaces ordered as ``UP_UP``, ``UP_DOWN``,
            ``DOWN_UP``, and ``DOWN_DOWN``. Device cross-section logs are removed
            and replaced with polarization state logs.

        Raises
        ------
        ValueError
            If exactly four device cross-section workspaces are not supplied.
        numpy.linalg.LinAlgError
            If the encoding matrix is singular at the workspace wavelength.
        """
        if len(device_cross_sections) != 4:
            raise ValueError(
                f"Full polarization requires exactly 4 device cross-sections, got {len(device_cross_sections)}."
            )

        # Find out which cross-section is S00, S10, S0pi, and S1pi.
        by_state = {PolarizationCrossSection.get(ws): ws for ws in device_cross_sections}
        s00 = by_state[PolarizationCrossSection.OFF_OFF]
        s10 = by_state[PolarizationCrossSection.ON_OFF]
        s0pi = by_state[PolarizationCrossSection.OFF_ON]
        s1pi = by_state[PolarizationCrossSection.ON_ON]
        spinor = np.array([s00, s10, s0pi, s1pi], dtype=object)

        wavelength = SampleLogs(s00).single_value("wavelength")
        M = np.linalg.inv(self.encoding_matrix(wavelength)).astype(object)

        spin_cross_sections = list(M @ spinor)
        for ws, state in zip(
            spin_cross_sections,
            [
                PolarizationState.UP_UP,
                PolarizationState.UP_DOWN,
                PolarizationState.DOWN_UP,
                PolarizationState.DOWN_DOWN,
            ],
        ):
            DeleteLog(Workspace=ws, Name=PolarizationCrossSection.logname)
            state.log(ws)
        return spin_cross_sections


def polarization_decoder(device_cross_sections, reduction_config):
    decoders = {"half": HalfPolarizationDecoder, "full": FullPolarizationDecoder}
    polarization_level = reduction_config["polarization"]["level"]
    decoder = decoders[polarization_level](reduction_config)
    return decoder.decode(device_cross_sections)


# A simple way to encode the name and specifications for one of the time generators methods of class SimulatedLogs
# Example: polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 1.0, "alive_duration": 0.2})
TimesGeneratorSpecs = namedtuple("TimesGeneratorSpecs", ["name", "kwargs"])


@dataclass
class SimulatedPolarizationLogs:
    """A simulated log for testing purposes."""

    polarizer: int = 0
    polarizer_flipper: Optional[TimesGeneratorSpecs] = None
    polarizer_veto: Optional[TimesGeneratorSpecs] = None
    analyzer: int = 0
    analyzer_flipper: Optional[TimesGeneratorSpecs] = None
    analyzer_veto: Optional[TimesGeneratorSpecs] = None

    # Class variables
    flipper_generators: ClassVar[List[str]] = ["heartbeat", "binary_pulse", "cycled_intervals"]
    veto_generators: ClassVar[List[str]] = ["binary_pulse"]  # available time-generators for veto intervals

    def __post_init__(self):
        # Validate input polarizer and analyzer flipper generators
        for flipper, device in zip([self.polarizer_flipper, self.analyzer_flipper], ["polarizer", "analyzer"]):
            if flipper and flipper.name not in self.flipper_generators:
                raise ValueError(
                    f"The {device} flipper generator must be one of {self.flipper_generators}, got '{flipper.name}'"
                )
        # Validate input polarizer and analyzer veto generators
        for veto, device in zip([self.polarizer_veto, self.analyzer_veto], ["polarizer", "analyzer"]):
            if veto and veto.name not in self.veto_generators:
                raise ValueError(
                    f"The {device} veto generator must be one of {self.veto_generators}, got '{veto.name}'"
                )

    def heartbeat(
        self, interval: float, dead_time: Optional[float] = 0.0, upper_bound: Optional[float] = None
    ) -> Generator[float, None, None]:
        """
        Generate a sequence of timestamps at regular intervals, starting at or later than dead_time.

        Parameters
        ----------
        interval : float
            The time interval between consecutive timestamps, in seconds.
        dead_time: float
            The initial time period, in seconds, during which no times are generated. Defaults to 0.0.
        upper_bound : float, optional
            The maximum time value to generate, in seconds. If None, the generator will continue indefinitely.

        Yields
        ------
        float
            The next timestamp in the sequence.

        Examples
        --------
        >>> log = SimulatedPolarizationLogs()
        >>> list(log.heartbeat(interval=1.0, upper_bound=5.0))
        [0, 1.0, 2.0, 3.0, 4.0, 5.0]
        """
        elapsed = 0
        while elapsed <= upper_bound if upper_bound is not None else True:
            if elapsed >= dead_time:
                yield elapsed
            elapsed += interval

    def cycled_intervals(
        self,
        intervals: List[float],
        dead_time: Optional[float] = 0.0,
        upper_bound: Optional[float] = None,
    ) -> Generator[float, None, None]:
        """
        Generate a sequence of timestamps by repeatedly cycling through a list of time intervals.

        Starting from time zero, each successive timestamp is computed by adding the next interval
        in ``intervals`` to the current elapsed time. Once all intervals have been consumed, the
        cycle restarts from the first interval. Timestamps are only yielded when at or after
        ``dead_time``.

        Parameters
        ----------
        intervals : list of float
            A list of time intervals, in seconds, to cycle through repeatedly.
        dead_time : float, optional
            The initial time period, in seconds, during which no times are generated. Defaults to 0.0.
        upper_bound : float, optional
            The maximum time value to generate, in seconds. If :py:obj:`None`, the generator
            will continue indefinitely.

        Yields
        ------
        float
            The next timestamp in the cycled intervals sequence.

        Examples
        --------
        >>> log = SimulatedPolarizationLogs()
        >>> list(log.cycled_intervals(intervals=[1.0, 2.0, 1.5], dead_time=0.0, upper_bound=10.0))
        [0.0, 1.0, 3.0, 4.5, 5.5, 7.5, 9.0]
        """
        accumulated, elapsed, i = [], 0.0, 0
        while elapsed <= upper_bound if upper_bound is not None else True:
            if elapsed >= dead_time:
                yield elapsed
            if i < len(intervals):
                accumulated.append(intervals[i])
                i += 1
            else:
                accumulated.append(intervals[0])
                i = 1
            # timestamps are computed with math.fsum over all accumulated intervals rather than by
            # incremental addition, because summing many small floats incrementally causes rounding errors to
            # drift (e.g. yielding 0.19999999 instead of 0.2).
            elapsed = fsum(accumulated)

    def binary_pulse(
        self,
        interval: float,
        alive_duration: float,
        dead_time: Optional[float] = 0.0,
        upper_bound: Optional[float] = None,
    ) -> Generator[float, None, None]:
        """
        Generate a sequence of timestamps with a binary pulse pattern, starting at or later than dead_time.

        The timestamps alternate between the start and end of a veto period, which is centered within each interval.

        Parameters
        ----------
        interval : float
            The time interval between consecutive pulses, in seconds.
        alive_duration : float
            The duration of the period above the zero baseline, in seconds. Must be less than `interval`.
        dead_time: float
            The initial time period, in seconds, during which no times are generated. Defaults to 0.0.
        upper_bound : float, optional
            The maximum time value to generate, in seconds. If None, the generator will continue indefinitely.

        Yields
        ------
        float
            The next timestamp in the binary pulse sequence.

        Examples
        --------
        >>> log = SimulatedPolarizationLogs()
        >>> list(log.binary_pulse(interval=3.0, alive_duration=1.0, upper_bound=10))
        [0, 2.5, 3.5, 5.5, 6.5, 8.5, 9.5]
        """
        if not alive_duration < interval:
            raise ValueError("Veto duration must be less than the interval")
        elapsed, latest_pulse_time, veto_half, continue_while = 0.0, interval, alive_duration / 2, True
        if dead_time == 0.0:
            yield elapsed
        while continue_while:
            for elapsed in [latest_pulse_time - veto_half, latest_pulse_time + veto_half]:
                if (upper_bound is None) or (elapsed <= upper_bound):
                    if elapsed >= dead_time:
                        yield elapsed
                else:
                    continue_while = False  # exit the outer while loop
                    break  # exit the immediate `for` loop
            latest_pulse_time += interval

    def times_generator(self, pv_name: str, **options: dict) -> Optional[Generator[float, None, None]]:
        """
        Create a generator that yields time points

        This method selects the appropriate time generator function (e.g., `heartbeat` or `binary_pulse`)
        based on the PV name and its associated generator specifications. Additional options can be passed
        to override or extend the generator specifications.

        Parameters
        ----------
        pv_name : str
            The name of the process variable (e.g., 'PolarizerFlipper', 'PolarizerVeto', etc.).
        **options : dict
            Additional keyword arguments to override or extend the generator's default arguments.

        Raises
        ------
        KeyError
            If the provided PV name does not match any known process variable.
        AttributeError
            If the generator function associated with the PV name is not found.

        Examples
        --------
        >>> logs = SimulatedPolarizationLogs(
        ...     polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 1.0}),
        ...     polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 1.0, "alive_duration": 0.4})
        ... )
        >>> list(logs.times_generator(PV_POLARIZER_FLIPPER, upper_bound=6.3))
        [0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        >>> list(logs.times_generator(PV_POLARIZER_VETO, upper_bound=6.3))
        [0.0, 0.8, 1.2, 1.8, 2.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2]
        >>> logs.times_generator(PV_ANALYZER_FLIPPER)
        None
        >>> logs.times_generator(PV_ANALYZER_VETO)
        None
        """
        # conversion between PV name and class field
        converter = {
            PV_POLARIZER_FLIPPER: self.polarizer_flipper,
            PV_POLARIZER_VETO: self.polarizer_veto,
            PV_ANALYZER_FLIPPER: self.analyzer_flipper,
            PV_ANALYZER_VETO: self.analyzer_veto,
        }
        specs_field = converter[pv_name]
        if specs_field is None:
            return None
        generator_function = getattr(self, specs_field.name)
        kwargs = {**specs_field.kwargs, **options}
        return generator_function(**kwargs)

    def inject(self, input_workspace: MantidWorkspace):
        """
        Injects simulated log data into a Mantid workspace.

        This method adds polarizer and analyzer values as single-valued logs, and generates time-series logs
        for flippers and veto process variables based on their associated time-generator specifications.
        Values for flipper and veto time-series are either 0 or 1, and always start with 0 for simplicity.

        Parameters
        ----------
        input_workspace : MantidWorkspace
            The Mantid workspace into which the simulated logs will be injected.

        Raises
        ------
        AttributeError
            If the workspace does not contain required sample log entries like `run_start` or `duration`.

        Examples
        --------
        >>> workspace = CreateSingleValuedWorkspace(OutputWorkspace="example")
        >>> sample_logs = SampleLogs(workspace)
        >>> sample_logs.insert("start_time", "2023-10-01T00:00:00")
        >>> sample_logs.insert("duration", 300)
        >>> logs = SimulatedPolarizationLogs(
        ...     polarizer=1,
        ...     polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 1.0}),
        ...     polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 2.0, "alive_duration": 0.5})
        ... )
        >>> logs.inject(workspace)
        """
        sample_logs = SampleLogs(input_workspace)

        # Retrieve the run start time and duration from the sample logs, handling potential attribute differences.
        try:
            run_start = sample_logs.run_start.value
        except AttributeError:
            run_start = sample_logs.start_time.value
        duration: float = sample_logs.duration.value  # in seconds

        # insert polarizer and analyzer types
        sample_logs.insert(name=PV_POLARIZER, value=self.polarizer)
        sample_logs.insert(name=PV_ANALYZER, value=self.analyzer)

        # insert the time series
        for pv in [PV_POLARIZER_FLIPPER, PV_POLARIZER_VETO, PV_ANALYZER_FLIPPER, PV_ANALYZER_VETO]:
            times = self.times_generator(pv, upper_bound=duration)
            if times is None:  # no time generator specifications for this PV
                continue
            times = list(times)  # run the generator to get all times
            values = [i % 2 for i in range(len(times))]  # Alternating zeros and ones
            sample_logs.insert_time_series(name=pv, start_time=run_start, elapsed_times=times, values=values)
