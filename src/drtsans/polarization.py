# standard library imports
from dataclasses import dataclass
from typing import ClassVar, Generator, List, Optional, Union

# third party imports
from mantid.dataobjects import EventWorkspace
from mantid.simpleapi import CreateSingleValuedWorkspace, DeleteWorkspace, mtd, RenameWorkspace
import numpy as np

# drtsans imports
from drtsans.api import _set_uncertainty_from_numpy


__all__ = ["half_polarization"]


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


@dataclass
class SimulatedLogs:
    """A simulated log for testing purposes."""
    polarizer: int = 0
    polarizer_flipper: Optional[str] = None
    polarizer_veto: Optional[str] = None
    analyzer: int = 0
    analyzer_flipper: Optional[str] = None
    analyzer_veto: Optional[str] = None

    # Class variables
    flipper_generators: ClassVar[List[str]] = ["heartbeat"]
    veto_generators: ClassVar[List[str]] = ["binary_pulse"]

    def __post_init__(self):
        # Validate input polarizer and analyzer flipper generators
        for flipper, device in zip([self.polarizer_flipper, self.analyzer_flipper], ["polarizer", "analyzer"]):
            if flipper and flipper not in self.flipper_generators:
                raise ValueError(
                    f"The {device} flipper generator must be one of {self.flipper_generators}, got '{flipper}'"
                )
        # Validate input polarizer and analyzer veto generators
        for veto, device in zip([self.polarizer_veto, self.analyzer_veto], ["polarizer", "analyzer"]):
            if veto and veto not in self.veto_generators:
                raise ValueError(
                    f"The {device} veto generator must be one of {self.veto_generators}, got '{veto}'"
                )

    def heartbeat(self, interval: float, dead_time: Optional[float] = 0.0, upper_bound: Optional[float] = None) -> Generator[float, None, None]:
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
        >>> log = SimulatedLogs()
        >>> list(log.heartbeat(interval=1.0, upper_bound=5.0))
        [0, 1.0, 2.0, 3.0, 4.0, 5.0]
        """
        elapsed = 0
        while elapsed <= upper_bound if upper_bound is not None else True:
            if elapsed >= dead_time:
                yield elapsed
            elapsed += interval

    def binary_pulse(self, interval: float, veto_duration: float, dead_time: Optional[float] = 0.0, upper_bound: Optional[float] = None) -> Generator[float, None, None]:
        """
        Generate a sequence of timestamps with a binary pulse pattern, starting at or later than dead_time.

        The timestamps alternate between the start and end of a veto period, which is centered within each interval.

        Parameters
        ----------
        interval : float
            The time interval between consecutive pulses, in seconds.
        veto_duration : float
            The duration of the veto period, in seconds. Must be less than the interval.
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
        >>> log = SimulatedLogs()
        >>> list(log.binary_pulse(interval=3.0, veto_duration=1.0, upper_bound=10))
        [0, 2.5, 3.5, 5.5, 6.5, 8.5, 9.5]
        """
        assert veto_duration < interval, "Veto duration must be less than the interval"
        elapsed, intervals, veto_half, continue_while = 0.0, interval, veto_duration / 2, True
        yield elapsed
        while continue_while:
            for elapsed in [intervals - veto_half, intervals + veto_half]:
                if (upper_bound is None) or (elapsed <= upper_bound):
                    if elapsed >= dead_time:
                        yield elapsed
                else:
                    continue_while = False  # exit the outer while loop
                    break  # exit the inmediate `for` loop
            intervals += interval
