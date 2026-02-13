r"""
This module provides class `EQSANSDiskChopperSet` representing the set of disk choppers.
Prior to 2026, EQSANS had four disk choppers (two of them paired as a double chopper).
Starting in 2026, EQSANS has six disk choppers (three double choppers).
The main goal of the module is to find the set of neutron wavelength
bands transmitted by the chopper set, given definite choppers settings such as aperture and starting phase.
"""

import importlib
import json

import numpy as np
from drtsans.chopper import DiskChopper, DiskChopperSetConfiguration
from drtsans.samplelogs import SampleLogs
from drtsans.frame_mode import FrameMode
from drtsans.path import exists
from mantid.api import Run
from mantid.simpleapi import LoadNexusProcessed, mtd

from drtsans.wavelength import Wbands


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
        self.chopper_config: DiskChopperSetConfiguration = self.get_chopper_configuration(sample_logs.start_time.value)

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

    def get_chopper_configuration(self, start_time: str) -> DiskChopperSetConfiguration:
        r"""
        Get the chopper configuration (number of choppers, apertures, distances to source, offsets)
        from the JSON configuration file based on the log "start_time".

        Parameters
        ----------
        start_time
            String representing the run start time in the format: "YYYY-MM-DDThh:mm:ssZ"

        Returns
        -------
        DiskChopperSetConfiguration
            Configuration of the disk choppers.
        """
        # Get daystamp from sample logs (format: YYYYMMDD)
        start_time_str = start_time[0:10]  # "YYYY-MM-DD"
        daystamp = int(start_time_str.replace("-", ""))  # Convert to YYYYMMDD integer

        # Load configuration from JSON file
        with importlib.resources.open_text("drtsans.configuration", "EQSANS_chopper_configurations.json") as file:
            configs = json.load(file)

        return DiskChopperSetConfiguration.from_json(configs, daystamp)

    @property
    def period(self):
        return self._choppers[0].period

    def __getitem__(self, item):
        return self._choppers[item]

    @property
    def pulse_width(self):
        return self._pulse_width

    @property
    def _n_choppers(self):
        return self.chopper_config.n_choppers

    @property
    def _aperture(self):
        return self.chopper_config.aperture

    @property
    def _to_source(self):
        return self.chopper_config.to_source

    @property
    def _offsets(self):
        return self.chopper_config.offsets
