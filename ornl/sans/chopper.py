from __future__ import (absolute_import, division, print_function)

from ornl.sans.wavelength import Wband, Wbands


class DiskChopper(object):

    _pulse_width = 20  # micro-sec/Angstrom
    _cutoff_wl = 35  # maximum wavelength of incoming neutrons, in Angstroms

    def __init__(self, to_source, aperture, speed, sensor_phase, offset=0):
        r"""

        Parameters
        ----------
        to_source: float
            Distance to the neutron source (moderator) in meters
        aperture: float
            Opening window, in degrees
        speed: float
            rotational frequency, in Hz
        nominal_phase: float
            phase of the installed sensor
        offset: float
            phase between the installed sensor and the middle of the
            opening transmission window
        """
        self.to_source = to_source
        self.aperture = aperture
        self.speed = speed
        self.sensor_phase = sensor_phase
        self.offset = offset

    @property
    def pulse_width(self):
        r"""
        Neutrons of a given wavelength :math:`\lambda` emitted from the
        moderator have a distribution of delayed times characterized by
        a :match:`FWHM(\lambda) \simeq pulsewidth \cdot \lambda`.
        """
        return self._pulse_width

    @pulse_width.setter
    def pulse_width(self, value):
        r"""
        Override default pulse width, defined at the DiskChopper class level
        """
        self._pulse_width = value

    @property
    def cutoff_wl(self):
        r"""
        Discard neutrons emitted from the moderator with a wavelength above
        this quantity
        """
        return self._cutoff_wl

    @cutoff_wl.setter
    def cutoff_wl(self, value):
        r"""
        Override default cutoff wavelength, defined at the DiskChopper
        class level
        """
        self._cutoff_wl = value

    @property
    def phase(self):
        r"""
        Time needed by the chopper for the middle of the transmisison
        window to cut accross the axis defined by the neutron beam
        """
        return self.sensor_phase - self.offset

    @property
    def period(self):
        r"""
        Time required by the chopper for a full spin, in microseconds
        """
        return 1.0e6 / self.speed

    @property
    def transmission_duration(self):
        r"""
        Time for the transmission window to spin the whole aperture,
        in microseconds
        """
        return self.period * (self.aperture / 360.0)

    @property
    def opening_phase(self):
        r"""
        Time required by the opening edge of the transmission to hit
        the beam, in microseconds
        """
        return self.phase - 0.5 * self.transmission_duration

    @property
    def closing_phase(self):
        r"""
        Time required by the closing edge of the transmission to hit
        the beam, in microseconds
        """
        return self.phase + 0.5 * self.transmission_duration

    @property
    def rewind(self):
        r"""
        Rewind the chopper to obtain the minimum positive time
        for the closing edge of the transmission window to hit the
        axis defined by the neutron beam.

        Returns
        -------
        float
            Phase of the opening edge when the chopper is rewound.
            Can return a negative value
        """
        t_closing = self.closing_phase
        while t_closing < 0:
            t_closing += self.period
        return t_closing - self.transmission_duration

    def wavelength(self, tof, delay=0, pulsed=False):
        r"""
        Convert time of flight of the arriving neutron wavelength

        Parameters
        ----------
        tof: float
            time of flight, in microseconds
        delay: float
            Additional time-of-flight to include in the calculations
        pulsed: bool
            Include the correction due to delayed emission of neutrons
            from the moderator

        Returns
        -------
        float
            Neutron wavelength (in Angstroms). Returns zero for negative
            input `tof`

        """
        sigma = 3.9560346e-03  # plank constant divided by neutron mass
        loc = self.to_source
        if pulsed is True:
            loc += sigma * self._pulse_width
        wl = sigma * (tof + delay) / loc
        if wl < 0:
            wl = 0
        return wl

    def tof(self, wavelength, delay=0, pulsed=False):
        r"""
        Convert wavelength of arriving neutron to time of flight

        Parameters
        ----------
        wavelength: float
            wavelength of the arriving neutron, in microseconds
        delay: float
            Additional time-of-flight to include in the calculations
        pulsed: bool
            Include the correction due to delayed emission of neutrons
            from the moderator

        Returns
        -------
        float
            time of flight (in Angstroms)
        """
        sigma = 3.9560346e-03  # plank constant divided by neutron mass
        loc = self.to_source
        if pulsed is True:
            loc += sigma * self._pulse_width
        return wavelength * loc / sigma - delay

    def transmission_bands(self, cutoff_wl=None, delay=0, pulsed=False):
        r"""
        Wavelength bands transmitted by the chopper

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
            Set of transmitted wavelength bands
        """
        if cutoff_wl is None:
            cutoff_wl = self.cutoff_wl
        wb = Wbands()
        t_opening = self.rewind
        # shortest wavelength, obtained with pulsed correction if needed
        opening_wl = self.wavelength(t_opening, delay, pulsed)
        while opening_wl < cutoff_wl:
            # slowest wavelength, obtained with no pulse correction
            t = t_opening + self.transmission_duration
            closing_wl = self.wavelength(t, delay, False)
            if closing_wl > cutoff_wl:
                closing_wl = cutoff_wl
            wb += Wband(opening_wl, closing_wl)
            t_opening += self.period
            opening_wl = self.wavelength(t_opening, delay, pulsed)
        return wb
