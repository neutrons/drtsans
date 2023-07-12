# package imports
from drtsans.beam_finder import _calculate_neutron_drop
from drtsans.geometry import (
    panel_names,
    spectrum_info_ranges,
    get_pixel_masks,
    get_xyz,
    get_solid_angles,
)
from drtsans.samplelogs import SampleLogs

# third party imports
from mantid.api import mtd
from mantid.api import Workspace as MantidWorkspace
from mantid.kernel import DateAndTime
import numpy as np

# standard imports
import math
from typing import Callable, Optional, List, Union


def twotheta_to_xyz(twotheta_cross_section: Callable) -> Callable:
    r"""
    Decorator to pass pixel coordinates (x, y, z) to a cross section that instead accepts
    scattering angles as arguments

    Parameters
    ----------
    twotheta_cross_section
        Cross section function that accepts scattering angles as arguments
        twotheta_cross_section(two_thetas: np.ndarray) -> np.ndarray

    Returns
    -------
    Cross section function that uses the spatial pixel coordinates (x, y, z) as arguments

    """

    def xy_cross_section(x, y, z):
        two_thetas = np.arctan(np.sqrt(x**2 + y**2) / z)
        return twotheta_cross_section(two_thetas)

    return xy_cross_section


def insert_events(
    input_workspace: Union[str, MantidWorkspace],
    xy_cross_section: Callable,
    center_x: float = 0.0,
    center_y: float = 0.0,
    components: Optional[Union[str, List[str]]] = None,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
    gravity_correction: bool = True,
) -> None:
    r"""
    Insert events into one (or more) of the detector panels.

    Event TOF are randomly selected between 1000 and 16665 microseconds. All events are assigned to the first pulse.


    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    cross_section
        A function that takes the pixels' coordinates and returns the raw neutron counts.
        cross_section(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray
    center_x
        the X-coordinate at the point of intersection of the beam axis with the main detector when the
        main detector is centered at (0, 0, Z). Default: 0.0
    center_y
        the Y-coordinate at the point of intersection of the beam axis with the main detector when the
        main detector is centered at (0, 0, Z). Default: 0.0
    components
        One (or more) of the double-panels in the instrument (e.g. "detector1", "wing_detector"). If None, all
        components are used.
    component_efficiencies
        The efficiency of each component. If a single value is provided, it is used for all components. It is
        a multiplicative factor of the pixel's cross section.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the pixel cross section.
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the pixel.
    gravity_correction
        If True, the Y-component coordinate of each pixel is shifted "upwards" to correct the gravity drop in the
        neutron's path after scattering from the sample. This shift is applied before ``xy_cross_section``
        is called. There's only one value of the gravity drop for each component. Log entry "wavelength" is
        required for this correction.
    """
    workspace = mtd[str(input_workspace)]

    if isinstance(components, str):
        components = [components]

    # validate components exist:
    if components is not None:
        for component in components:
            if component not in panel_names(input_workspace):
                raise ValueError(f"Component {component} not found in workspace {str(input_workspace)}")
    else:
        components = panel_names(input_workspace)  # all components

    # validate efficiencies:
    if isinstance(component_efficiencies, (int, float)):
        efficiencies = [component_efficiencies] * len(components)
    else:
        assert len(components) == len(component_efficiencies)
        efficiencies = component_efficiencies

    for efficiency, component in zip(efficiencies, components):
        x, y, z = get_xyz(workspace, component)  # pixes coordinates, in meters
        # pixels coordinates (x, y) must be centered as (x-center_x, y-center_y) before calculating the cross-section
        x -= center_x
        y -= center_y
        if gravity_correction:
            wavelength = np.mean(SampleLogs(input_workspace).wavelength.value)
            position = workspace.getInstrument().getComponentByName(component).getPos().norm()
            y += _calculate_neutron_drop(position, wavelength)

        cross_sections = efficiency * xy_cross_section(x, y, z)  # simulated intensities

        if solid_angle_correction:
            solid_angles = get_solid_angles(workspace, component, back_panel_attenuation=back_panel_attenuation)
            cross_sections *= solid_angles
        neutron_counts = cross_sections.astype(int)
        mask = get_pixel_masks(workspace, component)
        neutron_counts[mask] = 0  # don't insert counts for masked pixels
        #
        # Insert as many events as neutron counts
        # generate random TOF values for each event
        TOF_MIN = 1000.0  # a thousand microseconds, naively selected
        TOF_MAX = 16665.0  # inverse of 60 Hz, in microseconds
        tofs = np.sort(np.random.randint(TOF_MIN, TOF_MAX + 1, max(neutron_counts))).astype(float)
        # we'll use one pulse time for all events
        pulse_time = DateAndTime(SampleLogs(workspace).start_time.value)

        first, next_to_last = spectrum_info_ranges(input_workspace, component)
        for idx, neutron_count in zip(range(first, next_to_last), neutron_counts):
            if neutron_count <= 0:
                continue
            spectrum = workspace.getSpectrum(idx)
            for i_tof in range(neutron_count):
                spectrum.addEventQuickly(tofs[i_tof], pulse_time)


def insert_beam_spot(
    input_workspace: Union[str, MantidWorkspace],
    center_x: float = 0.0,
    center_y: float = 0.0,
    diameter: float = 0.03,
    max_counts_in_pixel: int = 1000,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
    gravity_correction: bool = True,
) -> None:
    r"""
    Insert events into the main detector panel, simulating a beam spot.

    Neutron counts are simulated as a 2D Gaussian distribution, with the maximum counts at the center of the beam spot.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    center_x
        X-coordinate of the center of the beam spot, in meters. Default: 0.0
    center_y
        Y-coordinate of the center of the beam spot, in meters. Default: 0.0
    diameter
        Diameter of the beam spot, in meters
    max_counts_in_pixel
        Maximum number of counts in a single pixel, asumming the pixel is located at the center of the beam spot,
        on the front panel, and a nominal distance of 1 meter from the sample.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the intensity
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the intensity
    gravity_correction
        If True, the Y-component coordinate of each pixel is shifted "upwards" to correct the gravity drop in the
        neutron's path after scattering from the sample. This shift is applied before ``xy_cross_section``
        is called. There's only one value of the gravity drop for each component. Log entry "wavelength" is
        required for this correction.
    """

    def _xy_gaussian(x, y, _):
        return max_counts_in_pixel * np.exp(-(x**2 + y**2) / (2 * diameter**2))

    insert_events(
        input_workspace,
        _xy_gaussian,
        center_x=center_x,
        center_y=center_y,
        components="detector1",
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
        gravity_correction=gravity_correction,
    )


def insert_background(
    input_workspace: Union[str, MantidWorkspace],
    components: Optional[Union[str, List[str]]] = None,
    flavor: str = "gaussian noise",
    **flavor_kwargs: dict,
) -> None:
    r"""
    Insert background events into one (or more) of the detector panels.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    components
        One (or more) of the double-panels in the instrument (e.g. "detector1", "wing_detector"). If None, all
        components are used.
    flavor
        The type of background to insert. Supported flavors are:
         "gaussian noise", event counts randomly selected from a Gaussian distribution
         "flat noise", event counts randomly selected from a flat distribution
    flavor_kwargs
        Keyword arguments passed to the background generator function.
        For "gaussian noise", the following arguments are supported:
            mean: float
                Mean of the Gaussian distribution
            stddev: float
                Standard deviation of the Gaussian distribution
        For "flat noise", the following arguments are supported:
            min_counts: float
                Minimum number of neutrons counts (bigger or equal to 0)
            max_counts: float
                Maximum number of neutrons counts
    """

    def _gaussian_noise(x, y, _):
        assert len(x) == len(y)
        noise = np.random.normal(loc=flavor_kwargs["mean"], scale=flavor_kwargs["stddev"], size=len(x))
        noise[noise < 0] = 0.0
        return noise

    def _flat_noise(x, y, _):
        assert len(x) == len(y)
        assert flavor_kwargs["min_counts"] >= 0 and flavor_kwargs["min_counts"] < flavor_kwargs["max_counts"]
        return np.random.randint(flavor_kwargs["min_counts"], flavor_kwargs["max_counts"] + 1, len(x)).astype(float)

    event_generator = {"gaussian noise": _gaussian_noise, "flat noise": _flat_noise}[flavor]

    insert_events(
        input_workspace,
        event_generator,
        components=components,
        component_efficiencies=1.0,
        back_panel_attenuation=1.0,
        solid_angle_correction=False,
        gravity_correction=False,
    )


def insert_events_isotropic(
    input_workspace: Union[str, MantidWorkspace],
    center_x: float = 0.0,
    center_y: float = 0.0,
    components: Optional[Union[str, List[str]]] = None,
    max_counts_in_pixel: int = 1000,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
) -> None:
    r"""
    Isotropic scattering, all directions scatter equally.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    center_x
        X-coordinate of the center of the beam spot, in meters. Default: 0.0
    center_y
        Y-coordinate of the center of the beam spot, in meters. Default: 0.0
    components
        One (or more) of the double-panels in the instrument (e.g. "detector1", "wing_detector"). If None, all
        components are used.
    max_counts_in_pixel
        Maximum number of counts in a single pixel, asumming the pixel is located at the center of the beam spot,
        on the front panel, and a nominal distance of 1 meter from the sample.
    component_efficiencies
        The efficiency of each component. If a single value is provided, it is used for all components. It is
        a multiplicative factor of the pixel cross section.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the intensity
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the intensity
    """

    @twotheta_to_xyz
    def _twotheta_isotropic(twothetas):
        return np.repeat(1.0 * max_counts_in_pixel, len(twothetas))

    insert_events(
        input_workspace,
        _twotheta_isotropic,
        center_x=center_x,
        center_y=center_y,
        components=components,
        component_efficiencies=component_efficiencies,
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
        gravity_correction=False,
    )


def insert_events_ring(
    input_workspace: Union[str, MantidWorkspace],
    center_x: float = 0.0,
    center_y: float = 0.0,
    components: Optional[Union[str, List[str]]] = None,
    max_counts_in_pixel: int = 1000,
    twotheta_center: float = 6.0,
    twotheta_dev: float = 1.0,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
    gravity_correction: bool = True,
) -> None:
    r"""
    Insert events into the detector pixels of the designated components according to a gaussian function dependent
    on the two-theta scattering angle. The intensity pattern on the detectors is a single ring.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    center_x
        X-coordinate of the center of the beam spot, in meters. Default: 0.0
    center_y
        Y-coordinate of the center of the beam spot, in meters. Default: 0.0
    components
        One (or more) of the double-panels in the instrument (e.g. "detector1", "wing_detector"). If None, all
        components are used.
    max_counts_in_pixel
        Maximum number of counts in a single pixel, asumming the pixel is located at the center of the beam spot,
        on the front panel, and a nominal distance of 1 meter from the sample.
    twotheta_center
        the center of the gaussian function, in degrees. Default: 6.0. Must be a positive number.
    twotheta_dev
        the standard deviation of the gaussian function, in degrees. Default: 1.0. Must be a positive number.
    component_efficiencies
        The efficiency of each component. If a single value is provided, it is used for all components. It is
        a multiplicative factor of the pixel cross section.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the intensity
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the intensity
    gravity_correction
        If True, the Y-component coordinate of each pixel is shifted "upwards" to correct the gravity drop in the
        neutron's path after scattering from the sample. This shift is applied before ``xy_cross_section``
        is called. There's only one value of the gravity drop for each component. Log entry "wavelength" is
        required for this correction.
    """

    @twotheta_to_xyz
    def _twotheta_gaussian(twothetas):
        exponent = -((twothetas - np.radians(twotheta_center)) ** 2) / np.radians(twotheta_dev) ** 2
        return max_counts_in_pixel * np.exp(exponent)

    return insert_events(
        input_workspace,
        _twotheta_gaussian,
        center_x=center_x,
        center_y=center_y,
        components=components,
        component_efficiencies=component_efficiencies,
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
        gravity_correction=gravity_correction,
    )


def insert_events_sin_squared(
    input_workspace: Union[str, MantidWorkspace],
    center_x: float = 0.0,
    center_y: float = 0.0,
    components: Optional[Union[str, List[str]]] = None,
    max_counts_in_pixel: int = 1000,
    period: float = 16.0,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
    gravity_correction: bool = True,
) -> None:
    r"""
    Insert events into the detector pixels of the designated components according to a squared sinusoidal function
    dependent on the two-theta scattering angle. The intensity pattern on the detectors is a series of rings.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    center_x
        X-coordinate of the center of the beam spot, in meters. Default: 0.0
    center_y
        Y-coordinate of the center of the beam spot, in meters. Default: 0.0
    components
        One (or more) of the double-panels in the instrument (e.g. "detector1", "wing_detector"). If None, all
        components are used.
    max_counts_in_pixel
        Maximum number of counts in a single pixel, asumming the pixel is located at the center of the beam spot,
        on the front panel, and a nominal distance of 1 meter from the sample.
    period : float
        The period of the sin function, in degrees. Default: 16.0.
    component_efficiencies
        The efficiency of each component. If a single value is provided, it is used for all components. It is
        a multiplicative factor of the pixel cross section.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the intensity
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the intensity
    gravity_correction
        If True, the Y-component coordinate of each pixel is shifted "upwards" to correct the gravity drop in the
        neutron's path after scattering from the sample. This shift is applied before ``xy_cross_section``
        is called. There's only one value of the gravity drop for each component. Log entry "wavelength" is
        required for this correction.
    """

    @twotheta_to_xyz
    def _twotheta_sin_squared(twothetas):
        return max_counts_in_pixel * np.sin(2 * math.pi * twothetas / np.radians(period)) ** 2

    return insert_events(
        input_workspace,
        _twotheta_sin_squared,
        center_x=center_x,
        center_y=center_y,
        components=components,
        component_efficiencies=component_efficiencies,
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
        gravity_correction=gravity_correction,
    )
