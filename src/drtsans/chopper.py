r"""
This module provides class `DiskChopper` representing a rotating chopper with an aperture of certain width. The
main goal is to find the set of neutron wavelength bands transmitted by the chopper, given definite chopper
settings such as aperture and starting phase.
"""

from dataclasses import dataclass, field
from typing import Any
from drtsans.frame_mode import FrameMode
from drtsans.type_hints import EmissionDelay
from drtsans.wavelength import BROGLIE_NEUTRON, Wband, Wbands, from_tof as wavelength_from_tof


class DiskChopperSetConfigurationParsingError(Exception):
    """Raised when there is an error parsing the disk chopper set configuration from JSON."""


@dataclass(frozen=True)
class DiskChopperSetConfiguration:
    """Configuration for a set of disk choppers.

    Attributes
    ----------
    n_choppers: int
        Number of single disk choppers in the set.
    aperture: list[float]
        List of transmission aperture widths (in degrees) for each chopper.
    to_source: list[float]
        List of distances to the neutron source (in meters) for each chopper.
    offsets: dict[FrameMode, list[float]]
        Dictionary mapping frame modes to lists of offset phases (in micro seconds) for each chopper.
        These values are required to calibrate the value reported in the metadata. The combination on
        the reported phase and this offset is the time (starting from the current pulse) at which the
        middle of the choppers apertures will intersect with the neutron beam axis.
    """

    n_choppers: int
    aperture: list[float] = field(default_factory=list)
    to_source: list[float] = field(default_factory=list)
    offsets: dict[FrameMode, list[float]] = field(default_factory=dict)

    @classmethod
    def from_json(cls, json_config: Any, target_daystamp: int) -> "DiskChopperSetConfiguration":
        """Get chopper configuration from JSON object based on daystamp.

        Selects the configuration with the largest daystamp that is less than or equal to
        the target daystamp.

        Parameters
        ----------
        json_config
            JSON configuration
        target_daystamp
            8-digit integer whose digits are to be understood as YYYYMMDD (e.g., 20260209)

        Returns
        -------
        DiskChopperSetConfiguration
            The configuration object matching the target daystamp

        Raises
        ------
        DiskChopperSetConfigurationParsingError
             If there is an error parsing the JSON object or if no valid configuration is found for the target daystamp
        """
        # Find all entries that have the required keys
        required = {
            "n_choppers",
            "aperture",
            "to_source",
            "offsets",
            "daystamp",
        }
        valid_configs = [entry for entry in json_config if required.issubset(set([str(v) for v in entry.keys()]))]
        if not valid_configs:
            raise DiskChopperSetConfigurationParsingError("No valid configuration entries found")

        # Find all configurations with daystamp <= target_daystamp
        valid_configs = [cfg for cfg in valid_configs if cfg.get("daystamp", 0) <= target_daystamp]
        if not valid_configs:
            raise DiskChopperSetConfigurationParsingError(
                f"No valid configuration found on or before daystamp {target_daystamp}"
            )

        # Select the configuration with the largest daystamp
        selected = max(valid_configs, key=lambda x: x.get("daystamp", 0))

        # Parse and validate n_choppers
        try:
            n_choppers = int(selected["n_choppers"])
            if n_choppers <= 0:
                raise ValueError("n_choppers must be a positive integer")
        except (ValueError, TypeError) as e:
            raise DiskChopperSetConfigurationParsingError(f"Invalid n_choppers value '{selected['n_choppers']}': {e}")

        # Parse and validate aperture
        try:
            aperture = [float(x) for x in selected["aperture"]]
            if len(aperture) != n_choppers:
                raise ValueError(f"aperture list length ({len(aperture)}) does not match n_choppers ({n_choppers})")
        except (ValueError, TypeError) as e:
            raise DiskChopperSetConfigurationParsingError(f"Invalid aperture value '{selected['aperture']}': {e}")

        # Parse and validate to_source
        try:
            to_source = [float(x) for x in selected["to_source"]]
            if len(to_source) != n_choppers:
                raise ValueError(f"to_source list length ({len(to_source)}) does not match n_choppers ({n_choppers})")
        except (ValueError, TypeError) as e:
            raise DiskChopperSetConfigurationParsingError(f"Invalid to_source value '{selected['to_source']}': {e}")

        # Parse offsets from strings to FrameMode enums and validate
        offsets_dict = {}
        for mode_str, values in selected["offsets"].items():
            try:
                mode = FrameMode[mode_str]  # Convert string to FrameMode enum
                offset_values = [float(x) for x in values]
                if len(offset_values) != n_choppers:
                    raise ValueError(
                        f"offsets list length ({len(offset_values)}) for mode '{mode_str}' "
                        f"does not match n_choppers ({n_choppers})"
                    )
                offsets_dict[mode] = offset_values
            except KeyError:
                raise DiskChopperSetConfigurationParsingError(f"Invalid frame mode '{mode_str}' in offsets")
            except (ValueError, TypeError) as e:
                raise DiskChopperSetConfigurationParsingError(f"Invalid offsets value for mode '{mode_str}': {e}")

        return cls(
            n_choppers=n_choppers,
            aperture=aperture,
            to_source=to_source,
            offsets=offsets_dict,
        )


class DiskChopper:
    r"""
    Rotating disk chopper with an aperture of a certain width letting neutrons through.

    The angular position of the middle of the chopper aperture at the moment a neutron pulse happens is given
    by a metadata entry (the sensor phase) and an additional angle (the offset). The offset server to calibrate
    the value reported by the metadata.

    Parameters
    ----------
    to_source: float
        Distance to the neutron source (moderator) in meters
    aperture: float
        Width of the opening window letting neutrons through, in degrees
    speed: float
        rotational frequency, in Hz
    sensor_phase: float
        phase reported by the installed sensor in the metadata. It's the time for the chopper to rotate by
        and angle created by the following three points: (1) the center of the chopper, (2) the middle of
        the aperture, and (3) the point of intersection of the chopper and the pulse prompt neutrons. Units
        are micro seconds
    offset: float
        Additional phase difference to be added to the `sensor_phase` due to miscalibrations. The offset calibrates
        the value `sensor_phase` reported by the metadata. Units are in micro seconds.
    """

    #: The number of wavelength bands transmitted by a disk chopper is determined by the slowest emitted neutron,
    #: expressed as the maximum wavelength. This is the default cut-off maximum wavelength, in Angstroms.
    _cutoff_wl = 35

    def __init__(self, to_source, aperture, speed, sensor_phase, offset=0):
        self.to_source = to_source
        self.aperture = float(aperture)
        self.speed = float(speed)
        self.sensor_phase = float(sensor_phase)
        self.offset = float(offset)

    @property
    def cutoff_wl(self):
        r"""
        Discard neutrons transmitted by the disk chopper having a wavelength above this quantity. This
        property can override the default cutoff wavelength :const:`~drtsans.chopper.DiskChopper._cutoff_wl`.
        """
        return self._cutoff_wl

    @cutoff_wl.setter
    def cutoff_wl(self, value):
        r"""
        Override default cutoff wavelength :const:`~drtsans.chopper.DiskChopper._cutoff_wl`.
        """
        self._cutoff_wl = value

    @property
    def phase(self):
        r"""
        Time (starting from the current pulse) when the middle of the chopper aperture will
        intersect with the neutron beam axis, in micro seconds.
        """
        return self.sensor_phase - self.offset

    @property
    def period(self):
        r"""
        Time span required by the chopper for a full spin, in micro seconds.
        """
        return 1.0e6 / self.speed

    @property
    def transmission_duration(self):
        r"""
        Time span taking the chopper to spin an angle equal to its aperture, in micro seconds.
        """
        return self.period * (self.aperture / 360.0)

    @property
    def opening_phase(self):
        r"""
        Time (starting from the current pulse) when the opening edge of the chopper aperture will
        intersect with the neutron beam axis, in micro seconds.
        """
        return self.phase - 0.5 * self.transmission_duration

    @property
    def closing_phase(self):
        r"""
        Time (starting from the current pulse) when the closing edge of the chopper aperture will
        intersect with the neutron beam axis, in micro seconds.
        """
        return self.phase + 0.5 * self.transmission_duration

    @property
    def rewind(self):
        r"""
        Spin the chopper backwards until the chopper aperture intersects with the neutron beam axis.
        At this point, the :const:`~drtsans.chopper.DiskChopper.opening_phase` will be negative, and the
        :const:`~drtsans.chopper.DiskChopper.closing_phase` will be positive.

        Returns
        -------
        float
            Opening phase, in micro seconds. The opening phase will be negative (most likely) or zero.
        """
        t_closing = self.closing_phase
        while t_closing < 0:
            t_closing += self.period
        return t_closing - self.transmission_duration

    def wavelength(self, tof, delay=0, emission_delay: EmissionDelay = None):
        r"""
        Convert time-of-flight to neutron wavelength, for a neutron that has traveled the distance from the
        moderator to the chopper.

        Solves the implicit equation

        .. math::

            \frac{D}{v(\lambda)} + t_0(\lambda) = t_m + d

        where :math:`D` is the distance from moderator to chopper, :math:`v(\lambda) = h/(m\lambda)` is
        the neutron velocity, and :math:`t_0(\lambda)` is the delayed emission time from the moderator.
        When no emission delay is provided the equation is explicit:
        :math:`\lambda = \frac{h}{m} \frac{t_m + d}{D}`.

        Parameters
        ----------
        tof: float
            time of flight, in micro seconds
        delay: float
            Additional time-of-flight to include in the calculations. For instance, this could be a multiple
            of the the pulse period.
        emission_delay: callable, optional
            Function returning the delayed emission time (in microseconds) for a neutron of a given
            wavelength (in Angstroms). If :py:obj:`None`, no emission-time correction is applied.

        Returns
        -------
        float
            Neutron wavelength (in Angstroms). Returns zero for negative `tof`.
        """
        return wavelength_from_tof(tof, delay, distance=self.to_source, emission_delay=emission_delay)

    def tof(self, wavelength, delay=0, emission_delay: EmissionDelay = None):
        r"""
        Convert wavelength to *measured* time-of-flight, for a neutron that has traveled the distance from the
        moderator to the chopper.

        The measured time of flight :math:`t_m` is the sum of the actual flight time and the delayed emission
        time :math:`t_0(\lambda)` from the moderator, minus any additional delay :math:`d`:

        .. math::

            t_m = \frac{D}{v} + t_0(\lambda) - d = \frac{D \lambda}{h/m} + t_0(\lambda) - d

        where :math:`D` is the distance from moderator to chopper and :math:`v = h/(m\lambda)` is the neutron
        velocity.

        Parameters
        ----------
        wavelength: float
            wavelength of the neutron, in Angstroms.
        delay: float
            Additional time-of-flight to include in the calculations. For instance, this could be a multiple
            of the pulse period.
        emission_delay: callable, optional
            Function returning the delayed emission time (in microseconds) for a neutron of a given
            wavelength (in Angstroms). If :py:obj:`None`, no emission-time correction is applied.

        Returns
        -------
        float
            time-of-flight, in micro seconds.
        """
        velocity = BROGLIE_NEUTRON / wavelength  # neutron velocity, in meters/microsecond
        t0 = emission_delay(wavelength) if emission_delay is not None else 0.0
        return self.to_source / velocity + t0 - delay

    def transmission_bands(self, cutoff_wl=None, delay=0, emission_delay: EmissionDelay = None):
        r"""
        Wavelength bands transmitted by the chopper aperture. The number of bands is determined by the
        slowest neutrons emitted from the moderator.

        A neutron of wavelength :math:`\lambda` reaches the chopper at time
        :math:`t_0(\lambda) + D / v(\lambda)`, and is transmitted when that time falls inside the
        aperture's opening window. Both edges of each band therefore carry the emission-delay
        correction, which shifts the band towards shorter wavelengths by an amount inversely
        proportional to the chopper's distance to the source. Omitting ``emission_delay`` yields
        instead the purely geometric band, which is what the data acquisition system phases the
        choppers against.

        Parameters
        ----------
        cutoff_wl: float
            maximum wavelength of incoming neutrons. Discard slower neutrons when finding the transmission bands.
        delay: float
            Additional time-of-flight to include in the calculations. For instance, this could be a multiple
            of the the pulse period.
        emission_delay: callable, optional
            Function returning the delayed emission time (in microseconds) for a neutron of a given
            wavelength (in Angstroms). If :py:obj:`None`, no emission-time correction is applied.

        Returns
        -------
        ~drtsans.wavelength.Wbands
            Set of wavelength bands transmitted by the chopper.
        """
        if cutoff_wl is None:
            cutoff_wl = self.cutoff_wl
        wb = Wbands()
        t_opening = self.rewind
        # shortest wavelength, corrected for emission delay if provided
        opening_wl = self.wavelength(t_opening, delay, emission_delay)
        while opening_wl < cutoff_wl:
            # slowest wavelength, corrected for emission delay if provided
            t = t_opening + self.transmission_duration
            closing_wl = self.wavelength(t, delay, emission_delay)
            if closing_wl > cutoff_wl:
                closing_wl = cutoff_wl
            wb += Wband(opening_wl, closing_wl)
            t_opening += self.period
            opening_wl = self.wavelength(t_opening, delay, emission_delay)
        return wb
