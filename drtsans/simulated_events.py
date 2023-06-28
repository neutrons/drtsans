# package imports
from drtsans.geometry import (
    panel_names,
    spectrum_info_ranges,
    get_pixel_masks,
    get_xy,
    get_twothetas,
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


def insert_events(input_workspace: Union[str, MantidWorkspace], component: str, neutron_counts: np.ndarray) -> None:
    r"""
    Insert events into one of the detector panels.

    Event TOF are randomly selected between 1000 and 16665 microseconds. All events are assigned to the first pulse.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    component
        One of the double-panels in the instrument (e.g. "detector1", "wing_detector")
    neutron_counts
        Number of events to insert into each of the detector pixels of the components. The length of this array
        should be equal to the number of pixels in the component.
    """
    workspace = mtd[str(input_workspace)]
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


def insert_xy_events(
    input_workspace: Union[str, MantidWorkspace],
    xy_cross_section: Callable,
    components: Optional[Union[str, List[str]]] = None,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
) -> None:
    r"""
    Insert events into one (or more) of the detector panels.

    Event TOF are randomly selected between 1000 and 16665 microseconds. All events are assigned to the first pulse.


    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    xy_cross_section
        A function that takes the position of each pixel and returns a neutron count.
    components
        One (or more) of the double-panels in the instrument (e.g. "detector1", "wing_detector"). If None, all
        components are used.
    component_efficiencies
        The efficiency of each component. If a single value is provided, it is used for all components. It is
        a multiplicative factor of the pixel cross section.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the pixel cross section.
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the pixel.
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
        x, y = get_xy(workspace, component)  # pixes coordinates, in meters
        cross_sections = efficiency * xy_cross_section(x, y)  # simulated intensities
        if solid_angle_correction:
            solid_angles = get_solid_angles(workspace, component, back_panel_attenuation=back_panel_attenuation)
            cross_sections *= solid_angles
        neutron_counts = cross_sections.astype(int)
        mask = get_pixel_masks(workspace, component)
        neutron_counts[mask] = 0  # don't insert counts for masked pixels
        insert_events(input_workspace, component, neutron_counts)


def insert_twotheta_events(
    input_workspace: Union[str, MantidWorkspace],
    twotheta_cross_section: Callable,
    twotheta_units: str = "degrees",
    components: Optional[Union[str, List[str]]] = None,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
) -> None:
    r"""
    Insert events into the detector pixels of the designated components according to the given two-theta
    dependent cross section.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    twotheta_cross_section
        A function returning the number of neutrons scattered into a detector pixel facing the sample and
        placed at a unit distance from the sample. The function must take a single argument, a
        numpy array of scattering angles.
    twotheta_units
        the units of the scattering angles passed to the `theta_cross_section` function. Default: "degrees". Can be
        "degrees" or "radians".
    components : str, list
        the component(s) to insert events into. Default: "midrange_detector". Can be one component (passed on as a
        string) or a list of components (passed on as a list of strings).
    component_efficiencies
        The efficiency of each component. If a single value is provided, it is used for all components. It is
        a multiplicative factor of the pixel cross section.
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the pixel.
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
        two_thetas = get_twothetas(workspace, component, twotheta_units)
        cross_sections = efficiency * twotheta_cross_section(two_thetas)
        if solid_angle_correction:
            solid_angles = get_solid_angles(workspace, component, back_panel_attenuation=back_panel_attenuation)
            cross_sections *= solid_angles
        neutron_counts = cross_sections.astype(int)
        mask = get_pixel_masks(workspace, component)
        neutron_counts[mask] = 0  # don't insert counts for masked pixels
        insert_events(input_workspace, component, neutron_counts)


def insert_beam_spot(
    input_workspace: Union[str, MantidWorkspace],
    center_x: float,
    center_y: float,
    diameter: float = 0.03,
    max_counts_in_pixel: int = 1000,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
) -> None:
    r"""
    Insert events into the main detector panel, simulating a beam spot.

    Neutron counts are simulated as a 2D Gaussian distribution, with the maximum counts at the center of the beam spot.

    Parameters
    ----------
    input_workspace
        Mantid workspace or name of the workspace
    center_x
        X-coordinate of the center of the beam spot, in meters
    center_y
        Y-coordinate of the center of the beam spot, in meters
    diameter
        Diameter of the beam spot, in meters
    max_counts_in_pixel
        Maximum number of counts in a single pixel, asumming the pixel is located at the center of the beam spot,
        on the front panel, and a nominal distance of 1 meter from the sample.
    back_panel_attenuation
        The attenuation of the back panel. This is a multiplicative factor of the intensity
    solid_angle_correction
        If True, the solid angle of each pixel is taken into account. This is a multiplicative factor of the intensity
    """

    def _xy_gaussian(x, y):
        return max_counts_in_pixel * np.exp(-((x - center_x) ** 2 + (y - center_y) ** 2) / (2 * diameter**2))

    insert_xy_events(
        input_workspace,
        _xy_gaussian,
        components="detector1",
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
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

    def _gaussian_noise(x, y):
        assert len(x) == len(y)
        noise = np.random.normal(loc=flavor_kwargs["mean"], scale=flavor_kwargs["stddev"], size=len(x))
        noise[noise < 0] = 0.0
        return noise

    def _flat_noise(x, y):
        assert len(x) == len(y)
        assert flavor_kwargs["min_counts"] >= 0
        return np.random.randint(flavor_kwargs["min_counts"], flavor_kwargs["max_counts"] + 1, len(x)).astype(float)

    event_generator = {"gaussian noise": _gaussian_noise, "flat noise": _flat_noise}[flavor]

    insert_xy_events(
        input_workspace,
        event_generator,
        components=components,
        component_efficiencies=1.0,
        back_panel_attenuation=1.0,
        solid_angle_correction=False,
    )


def insert_events_isotropic(
    input_workspace: Union[str, MantidWorkspace],
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

    def _twotheta_isotropic(twothetas):
        return np.repeat(1.0 * max_counts_in_pixel, len(twothetas))

    insert_twotheta_events(
        input_workspace,
        _twotheta_isotropic,
        components=components,
        component_efficiencies=component_efficiencies,
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
    )


def insert_events_ring(
    input_workspace: Union[str, MantidWorkspace],
    components: Optional[Union[str, List[str]]] = None,
    max_counts_in_pixel: int = 1000,
    twotheta_center: float = 6.0,
    twotheta_dev: float = 1.0,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
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
    """

    def _twotheta_gaussian(twothetas):
        return max_counts_in_pixel * np.exp(-((twothetas - twotheta_center) ** 2) / twotheta_dev**2)

    return insert_twotheta_events(
        input_workspace,
        _twotheta_gaussian,
        twotheta_units="degrees",
        components=components,
        component_efficiencies=component_efficiencies,
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
    )


def insert_events_sin_squared(
    input_workspace: Union[str, MantidWorkspace],
    components: Optional[Union[str, List[str]]] = None,
    max_counts_in_pixel: int = 1000,
    period: float = 16.0,
    component_efficiencies: Union[int, float, List] = 1.0,
    back_panel_attenuation: float = 0.5,
    solid_angle_correction: bool = True,
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
    """

    def _twotheta_sin_squared(twothetas):
        return max_counts_in_pixel * np.sin(2 * math.pi * twothetas / period) ** 2

    return insert_twotheta_events(
        input_workspace,
        _twotheta_sin_squared,
        twotheta_units="degrees",
        components=components,
        component_efficiencies=component_efficiencies,
        back_panel_attenuation=back_panel_attenuation,
        solid_angle_correction=solid_angle_correction,
    )
