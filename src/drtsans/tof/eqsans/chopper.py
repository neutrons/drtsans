r"""
This module provides class `EQSANSDiskChopperSet` representing the set of disk choppers.
Prior to 2026, EQSANS had four disk choppers (two of them paired as a double chopper).
Starting in 2026, EQSANS has six disk choppers (three double choppers).
The main goal of the module is to find the set of neutron wavelength
bands transmitted by the chopper set, given definite choppers settings such as aperture and starting phase.
"""

from dataclasses import dataclass, field

import numpy as np
from drtsans.chopper import DiskChopper, DiskChopperConfiguration
from drtsans.samplelogs import SampleLogs
from drtsans.frame_mode import FrameMode
from drtsans.path import exists
from mantid.api import Run
from mantid.simpleapi import LoadNexusProcessed, mtd

from drtsans.wavelength import Wbands


EQSANS_SIX_CHOPPERS_LOGS = ["Speed5", "Phase5", "Speed6", "Phase6"]


@dataclass(frozen=True)
class EQSANSFourChoppersConfiguration(DiskChopperConfiguration):
    """Configuration for the set of four disk choppers in EQSANS.

    This is the configuration used at EQSANS prior to the installation of the
    additional two choppers in 2026.
    """

    n_choppers: int = 4
    aperture: list[float] = field(default_factory=lambda: [129.605, 179.989, 230.010, 230.007])
    to_source: list[float] = field(default_factory=lambda: [5.700, 7.800, 9.497, 9.507])
    offsets: dict[FrameMode, list[float]] = field(
        default_factory=lambda: {
            FrameMode.not_skip: [9507.0, 9471.0, 9829.7, 9584.3],
            FrameMode.skip: [19024.0, 18820.0, 19714.0, 19360.0],
        }
    )


@dataclass(frozen=True)
class EQSANSSixChoppersConfiguration(DiskChopperConfiguration):
    """Configuration for the set of six disk choppers in EQSANS.

    This is the configuration used at EQSANS after the installation of the
    additional two choppers in 2026, resulting in three double choppers.
    """

    n_choppers: int = 6
    aperture: list[float] = field(default_factory=lambda: [129.600, 180.000, 230.010, 230.007, 129.600, 180.000])
    to_source: list[float] = field(default_factory=lambda: [5.7178, 7.7998, 9.4998, 9.5058, 5.7238, 7.8058])
    offsets: dict[FrameMode, list[float]] = field(
        default_factory=lambda: {
            # TODO: Update these offsets when provided by the instrument scientists (EWM 15020)
            FrameMode.not_skip: [9507.0, 9471.0, 9829.7, 9584.3, 0.0, 0.0],
            FrameMode.skip: [19024.0, 18820.0, 19714.0, 19360.0, 0.0, 0.0],
        }
    )


class EQSANSDiskChopperSet:
    r"""
    Set of disks choppers installed in EQSANS.

    Parameters
    ----------
    other: file name, workspace, Run object, run number
        Load the chopper settings from this object.
    """

    #: Neutrons of a given wavelength :math:`\lambda` emitted from the moderator follow a distribution of delayed
    #: emission times that depends on the wavelength, and is characterized by function
    #: :math:`FWHM(\lambda) \simeq pulsewidth \cdot \lambda`.
    #: This is the default :math:`pulsewidth` in micro-sec/Angstrom.
    _pulse_width = 20

    #: The number of wavelength bands transmitted by a disk chopper is determined by the slowest emitted neutron,
    #: expressed as the maximum wavelength. This is the default cut-off maximum wavelength, in Angstroms.
    _cutoff_wl = 35

    def __init__(self, other):
        # Load choppers settings from the logs
        if isinstance(other, Run) or str(other) in mtd:
            sample_logs = SampleLogs(other)
        elif exists(other):
            ws = LoadNexusProcessed(other)
            sample_logs = SampleLogs(ws)
        else:
            raise RuntimeError("{} is not a valid file name, workspace, Run object or run number".format(other))

        # Get the chopper configuration (4 or 6 choppers)
        self.chopper_config: DiskChopperConfiguration = self.get_chopper_configuration(sample_logs)

        self._n_choppers = self.chopper_config.n_choppers
        self._aperture = self.chopper_config.aperture
        self._to_source = self.chopper_config.to_source
        self._offsets = self.chopper_config.offsets

        self._choppers = list()
        for chopper_index in range(self._n_choppers):
            aperture = self._aperture[chopper_index]
            to_source = self._to_source[chopper_index]
            speed = sample_logs["Speed{}".format(1 + chopper_index)].value.mean()
            sensor_phase = sample_logs["Phase{}".format(1 + chopper_index)].value.mean()
            ch = DiskChopper(to_source, aperture, speed, sensor_phase)
            ch.pulse_width = self._pulse_width
            ch.cutoff_wl = self._cutoff_wl
            self._choppers.append(ch)

        # Determine period and if frame skipping mode from the first chopper
        ch = self._choppers[0]
        # example of frame skipping: chopper speed 30 Hz, pulse frequency 60 Hz: abs(30 - 60) / 2 = 15
        condition = abs(ch.speed - sample_logs.frequency.value.mean()) / 2 > 1
        self.frame_mode = FrameMode.skip if condition else FrameMode.not_skip

        # Select appropriate offsets, based on the frame-skip mode.
        for chopper_index in range(self._n_choppers):
            ch = self._choppers[chopper_index]
            ch.offset = self._offsets[self.frame_mode][chopper_index]

    def transmission_bands(self, cutoff_wl: float = None, delay: float = 0, pulsed: bool = False) -> Wbands:
        r"""
        Wavelength bands transmitted by the chopper apertures. The number of bands is determined by the
        slowest neutrons emitted from the moderator.

        Parameters
        ----------
        cutoff_wl: float
            maximum wavelength of incoming neutrons. Discard slower neutrons when finding the transmission bands.
        delay: float
            Additional time-of-flight to include in the calculations. For instance, this could be a multiple
            of the the pulse period.
        pulsed: bool
            Include a correction due to delayed emission of neutrons from the moderator. See
            :const:`~drtsans.tof.eqsans.chopper.EQSANSDiskChopperSet._pulse_width` for a
            more detailed explanation.

        Returns
        -------
        ~drtsans.wavelength.Wbands
            Set of wavelength bands transmitted by the chopper.
        """
        if cutoff_wl is None:
            cutoff_wl = self._cutoff_wl
        # Filter out the choppers with zero speed, which do not contribute to the transmission bands
        moving_choppers = [ch for ch in self._choppers if not np.isclose(ch.speed, 0.0)]
        if not moving_choppers:
            return Wbands()
        # Transmission bands of the first chopper
        wb = moving_choppers[0].transmission_bands(cutoff_wl, delay, pulsed)
        # Find the common transmitted bands between the first chopper
        # and the ensuing choppers
        for ch in moving_choppers[1:]:
            wb_other = ch.transmission_bands(cutoff_wl, delay, pulsed)
            wb *= wb_other
        # We end up with the transmission bands of the chopper set
        return wb

    def get_chopper_configuration(self, sample_logs: SampleLogs) -> DiskChopperConfiguration:
        r"""
        Get the chopper configuration (number of choppers, apertures, distances to source, offsets)
        based on the sample logs.

        Parameters
        ----------
        sample_logs: ~drtsans.samplelogs.SampleLogs
            Sample logs containing the run information.

        Returns
        -------
        DiskChopperConfiguration
            Configuration of the disk choppers.
        """
        has_six_choppers = all(log in sample_logs for log in EQSANS_SIX_CHOPPERS_LOGS)
        if has_six_choppers:
            return EQSANSSixChoppersConfiguration()
        else:
            return EQSANSFourChoppersConfiguration()

    @property
    def period(self):
        return self._choppers[0].period

    def __getitem__(self, item):
        return self._choppers[item]

    @property
    def pulse_width(self):
        return self._pulse_width
