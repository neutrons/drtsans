from __future__ import (absolute_import, division, print_function)

from ornl.sans.chopper import DiskChopper
from ornl.sans.samplelogs import SampleLogs
from enum import Enum


class FrameMode(Enum):
    r"""
    Selects if instrument operating in frame-skipping mode
    """
    not_skip = 0
    skip = 1


class EQSANSDiskChopperSet(object):

    _pulse_width = 20  # micro-sec/Angstrom
    # Discard neutrons with wavelengths above this cutoff, in Angstroms
    _cutoff_wl = 35
    _n_choppers = 4
    # Opening angles
    _aperture = [129.605, 179.989, 230.010, 230.007]
    # Distance to neutron source (moderator), in meters
    _to_source = [5.700, 7.800, 9.497, 9.507]
    # Phase offsets, in micro-seconds.
    _offsets = {FrameMode.not_skip: [9507., 9471., 9829.7, 9584.3],
                FrameMode.skip: [19024., 18820., 19714., 19360.]}

    def __init__(self, other):
        r"""
        Customized disk chopper object for EQSANS

        Parameters
        ----------
        chopper_index: int
            From zero to three, since we have four choppers
        other: file name, workspace, Run object, run number
            Load the sample logs from this object to determine chopper settings
        """

        # Load choppers settings from the logs
        sl = SampleLogs(other)
        self._choppers = list()
        for chopper_index in range(self._n_choppers):
            aperture = self._aperture[chopper_index]
            to_source = self._to_source[chopper_index]
            speed = sl['Speed{}'.format(1 + chopper_index)].value.mean()
            sensor_phase = sl['Phase{}'.format(1 + chopper_index)].value.mean()
            ch = DiskChopper(to_source, aperture, speed, sensor_phase)
            ch.pulse_width = self._pulse_width
            ch.cutoff_wl = self._cutoff_wl
            self._choppers.append(ch)

        # Determine period and if frame skipping mode from the first chopper
        ch = self._choppers[0]
        self.period = ch.period
        condition = abs(ch.speed - sl.frequency.value.mean()) / 2 > 1
        self.frame_mode = FrameMode.skip if condition else FrameMode.not_skip

        # Select appropriate offsets
        for chopper_index in range(self._n_choppers):
            ch = self._choppers[chopper_index]
            ch.offset = self._offsets[self.frame_mode][chopper_index]

    def transmission_bands(self, cutoff_wl=None, delay=0, pulsed=False):
        r"""
        Wavelength bands transmitted by the chopper set

        Parameters
        ----------
        cutoff_wl: float
            maximum wavelength of incoming neutrons
        delay: float
            additional time-of-flight to include in the calculations
        pulsed: bool
            correction due to delayed emission of neutrons from the moderator

        Returns
        -------
        Wbands
        """
        ch = self._choppers[0]
        wb = ch.transmission_bands(cutoff_wl, delay, pulsed)
        for ch in self._choppers[1:]:
            wb *= ch.transmission_bands(cutoff_wl, delay, pulsed)
        return wb
