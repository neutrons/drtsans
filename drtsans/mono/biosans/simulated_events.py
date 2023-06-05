# package imports
from drtsans.instruments import fetch_idf
from drtsans.samplelogs import SampleLogs
from drtsans.mono.biosans.geometry import (
    get_twothetas,
    get_solid_angles,
    get_position_south_detector,
    has_midrange_detector,
    info_ids,
)

# third party imports
from mantid.api import mtd
from mantid.kernel import DateAndTime
from mantid.simpleapi import LoadInstrument, MoveInstrumentComponent
import numpy as np

# standard library imports
import math
import tempfile


def update_idf(input_workspace):
    r"""
    Download the latest BIOSANS IDF from the Mantid GitHub repository and load it into a workspace.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to load the instrument into.

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with the instrument loaded.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        idf = fetch_idf("BIOSANS_Definition.xml", output_directory=temp_dir)
        LoadInstrument(Workspace=input_workspace, Filename=idf, RewriteSpectraMap=True)
    return mtd[str(input_workspace)]


def insert_twotheta_events(
    input_workspace,
    twotheta_cross_section,
    twotheta_units="degrees",
    components="midrange_detector",
    efficiencies=1.0,
    solid_angle_correction=True,
):
    r"""
    Insert events into the detector pixels of the designated components according to the given two-theta
    dependent cross section.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to insert events into.

    twotheta_cross_section : function
        A function returning the number of neutrons scattered into a BIOSANS detector pixel facing the sample and
        placed at a unit distance from the sample. The function must take a single argument, a
        numpy array of scattering angles.

    twotheta_units : str
        the units of the scattering angles passed to the `theta_cross_section` function. Default: "degrees". Can be
        "degrees" or "radians".

    components : str, list
        the component(s) to insert events into. Default: "midrange_detector". Can be one component (passed on as a
        string) or a list of components (passed on as a list of strings).

    efficiencies : float, list
        the probability that a neutron hitting a pixel in the component will be detected. Default: 1.0. Can be a
        single value (passed on as a float) or a list of values (passed on as a list of floats).

    solid_angle_correction : bool
        correct the number of neutrons returned by the `theta_cross_section` function by the solid angle of the
        pixel, cos(vertical_angle) / L2^2 = R / L2^3, where L2 is the distance from the sample to the pixel center
        and R is the projection of L2 onto the horizontal plane. Thus, this correction equals 1 for a pixel detector
        positioned on the horizontal plane at a unit distance from the sample.

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with events inserted.
    """
    workspace = mtd[str(input_workspace)]
    if "midrange_detector" in components and not has_midrange_detector(workspace):
        workspace = update_idf(workspace)

    if isinstance(components, str):
        components = [components]  # promote to a list
    if isinstance(efficiencies, (int, float)):
        efficiencies = [efficiencies]
    assert len(components) == len(efficiencies)

    tof_min = 1000.0  # a thousand microseconds, naively selected
    tof_max = 16665.0  # inverse of 60 Hz, in microseconds
    pulse_time = DateAndTime(SampleLogs(workspace).start_time.value)  # add one second to the start time
    for efficiency, component in zip(efficiencies, components):
        two_thetas = get_twothetas(workspace, component, twotheta_units)
        cross_sections = efficiency * twotheta_cross_section(two_thetas)
        if solid_angle_correction:
            solid_angles = get_solid_angles(workspace, component)
            cross_sections *= solid_angles
        neutron_counts = cross_sections.astype(int)
        tofs = np.sort(np.random.randint(tof_min, tof_max + 1, max(neutron_counts))).astype(float)
        first, next_to_last = info_ids[component]["spectrum_info_range"]
        for idx, neutron_count in zip(range(first, next_to_last), neutron_counts):
            spectrum = workspace.getSpectrum(idx)
            for i_tof in range(neutron_count):
                spectrum.addEventQuickly(tofs[i_tof], pulse_time)
    return workspace


def insert_events_isotropic(input_workspace, counts_per_pixel=1000, components="midrange_detector", efficiencies=1.0):
    r"""
    Isotropic scattering, all directions scatter equally.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to insert events into.

    twotheta_center : float
        the center of the gaussian function, in degrees. Default: 6.0. Must be a positive number.

    twotheta_dev : float
        the standard deviation of the gaussian function, in degrees. Default: 1.0. Must be a positive number.

    max_counts_in_pixel : int
        The maximum number of counts to insert into a single pixel of optimal efficiency and situated in the front
        panel of the component, at 1 meter from the sample. Default: 1000.

    components : str, list
        the component(s) to insert events into. Default: "midrange_detector". Can be one component (passed on as a
        string) or a list of components (passed on as a list of strings).

    efficiencies : float, list
        the probability that a neutron hitting a pixel in the component will be detected. Default: 1.0. Can be a
        single value (passed on as a float) or a list of values (passed on as a list of floats).

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with events inserted.
    """

    def _twotheta_isotropic(twothetas):
        return np.repeat(1.0 * counts_per_pixel, len(twothetas))

    return insert_twotheta_events(
        input_workspace, _twotheta_isotropic, components=components, efficiencies=efficiencies
    )


def insert_events_ring(
    input_workspace,
    twotheta_center=6.0,
    twotheta_dev=1.0,
    max_counts_in_pixel=1000,
    components="midrange_detector",
    efficiencies=1.0,
):
    r"""
    Insert events into the detector pixels of the designated components according to a gaussian function dependent
    on the two-theta scattering angle. The intensity pattern on the detectors is a single ring.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to insert events into.

    twotheta_center : float
        the center of the gaussian function, in degrees. Default: 6.0. Must be a positive number.

    twotheta_dev : float
        the standard deviation of the gaussian function, in degrees. Default: 1.0. Must be a positive number.

    max_counts_in_pixel : int
        The maximum number of counts to insert into a single pixel of optimal efficiency and situated in the front
        panel of the component, at 1 meter from the sample. Default: 1000.

    components : str, list
        the component(s) to insert events into. Default: "midrange_detector". Can be one component (passed on as a
        string) or a list of components (passed on as a list of strings).

    efficiencies : float, list
        the probability that a neutron hitting a pixel in the component will be detected. Default: 1.0. Can be a
        single value (passed on as a float) or a list of values (passed on as a list of floats).

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with events inserted.
    """

    def _twotheta_gaussian(twothetas):
        return max_counts_in_pixel * np.exp(-((twothetas - twotheta_center) ** 2) / twotheta_dev**2)

    return insert_twotheta_events(
        input_workspace, _twotheta_gaussian, twotheta_units="degrees", components=components, efficiencies=efficiencies
    )


def insert_events_sin_squared(
    input_workspace, period=16.0, max_counts_in_pixel=1000, components="midrange_detector", efficiencies=1.0
):
    r"""
    Insert events into the detector pixels of the designated components according to a squared sinusoidal function
    dependent on the two-theta scattering angle. The intensity pattern on the detectors is a series of rings.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to insert events into.

    period : float
        The period of the sin function, in degrees. Default: 16.0.

    max_counts_in_pixel : int
        The maximum number of counts to insert into a single pixel of optimal efficiency and situated in the front
        panel of the component, at 1 meter from the sample. Default: 1000.

    components : str, list
        the component(s) to insert events into. Default: "midrange_detector". Can be one component (passed on as a
        string) or a list of components (passed on as a list of strings).

    efficiencies : float, list
        the probability that a neutron hitting a pixel in the component will be detected. Default: 1.0. Can be a
        single value (passed on as a float) or a list of values (passed on as a list of floats).

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with events inserted.
    """

    def _twotheta_sin_squared(twothetas):
        return max_counts_in_pixel * np.sin(2 * math.pi * twothetas / period) ** 2

    return insert_twotheta_events(
        input_workspace,
        _twotheta_sin_squared,
        twotheta_units="degrees",
        components=components,
        efficiencies=efficiencies,
    )


def insert_events_offset_center(
    input_workspace, detector1_shift_x, detector1_shift_y, max_counts=10000, bin_center_smear=0.05
):
    r"""
    Insert events into the South detector and shift its center on the XY plane. The intensity pattern is a smeared
    bright spot in the center of the detector.

    The intensity at each pixel is rescaled by its solid angle, by the efficiency of the component, and whether
    the pixel is located in the front or back panel of the component.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to insert events into.
    detector1_shift_x : float
        translate the South detector center by this amount in the X direction, in meters.
    detector1_shift_y : float
        translate the South detector center by this amount in the Y direction, in meters.
    max_counts : int
        The maximum number of counts to insert into a single pixel of optimal efficiency and situated in the front
        panel of the component, at 1 meter from the sample. Default: 10000.
    bin_center_smear : float
        the width of the bright spot. Default: 0.05.

    Returns
    -------

    """
    twotheta_dev = float(np.degrees(bin_center_smear / get_position_south_detector(input_workspace)))
    insert_events_ring(
        input_workspace,
        twotheta_center=0.0,
        twotheta_dev=twotheta_dev,
        max_counts_in_pixel=max_counts,
        components="detector1",
        efficiencies=1.0,
    )
    MoveInstrumentComponent(
        Workspace=input_workspace,
        ComponentName="detector1",
        X=detector1_shift_x,
        Y=detector1_shift_y,
        RelativePosition=True,
    )
