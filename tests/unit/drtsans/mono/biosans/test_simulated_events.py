# standard imports
from pathlib import Path
from typing import Any, Callable, List, Union
from unittest.mock import patch as mock_patch

# third party imports
from mantid.dataobjects import EventWorkspace
from mantid.kernel import DateAndTime
from mantid.simpleapi import (
    AddSampleLog,
    DeleteWorkspace,
    LoadEmptyInstrument,
    LoadNexusProcessed,
    mtd,
    Plus,
    SaveNexus,
)
import numpy as np
import pytest

# drtsans imports
from drtsans.instruments import empty_instrument_workspace
from drtsans.mono.load import transform_to_wavelength
from drtsans.mono.biosans import (
    load_all_files,
    reduce_single_configuration,
    reduction_parameters,
    update_reduction_parameters,
)
from drtsans.mono.biosans.geometry import (
    set_angle_midrange_detector,
    set_angle_wing_detector,
    set_position_south_detector,
)
from drtsans.mono.biosans.simulated_events import update_idf
from drtsans.simulated_events import insert_background, insert_beam_spot, insert_events_isotropic, insert_events_ring
from drtsans.samplelogs import SampleLogs
from drtsans.settings import namedtuplefy


def test_update_idf(fetch_idf, temp_workspace_name):
    idf_old = fetch_idf("BIOSANS_Definition_2019_2023.xml")
    workspace = LoadEmptyInstrument(InstrumentName="BIOSANS", Filename=idf_old, OutputWorkspace=temp_workspace_name())
    assert workspace.getInstrument().getComponentByName("midrange_detector") is None
    workspace = update_idf(workspace)
    assert workspace.getInstrument().getComponentByName("midrange_detector")


def construct_file_name(run_number: Union[str, int]) -> str:
    r"""Convention in this module to construct a file name from the run number"""
    return f"CG3_{run_number}.nxs"


def missing_files(root_path: Path, run_numbers: List[str], file_names: List[str] = list()) -> bool:
    r"""
    Find if the necessary files for the reduction are missing

    Parameters
    ----------
    root_path
        Directory containing all requested files
    run_numbers
        list of run numbers for which corresponding Nexus files are expected within ``root_path``
    file_names
        list of file names expected within ``root_path``

    Returns
    -------
    True if one or more file is missing or if ``root_path`` is not a directory
    """
    root_path = Path(root_path)  # just in case :)
    if not root_path.is_dir():
        return True
    for run_number in run_numbers:
        if not Path(root_path / construct_file_name(run_number)).is_file():
            return True
    for name in file_names:
        if not Path(root_path / name).is_file():
            return True
    return False


def create_six_rings_pattern(config: dict, metadata: dict):
    r"""
    Create a set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of six rings distributed in the detectors of BIOSANS,
    plus additional files to simulate the:
    - transmission of the sample
    - background run (run without the sample)
    - background transmission run
    - empty transmission run

    Parameters
    ----------
    config
        configuration for the reduction workflow
    metadata
        configuration for the instrument, run, and neutrons
    """
    # verify directory exists or can be created if it does not exist
    root_path = Path(config["dataDirectories"])
    if not root_path.is_dir():
        root_path.mkdir(parents=True)

    # we're going to be reusing some intensity patterns more than once. For instance, we'll reuse a
    # flat noise pattern for the sample run as well as the transmission run. Thus,
    # cache these intensity patterns
    events_cache = {}

    def save_and_delete(
        input_workspace: Union[str, EventWorkspace], run_number: Union[int, str] = None, filename: str = None
    ):
        r"""Save `input_workspace` to a Nexus file associated to `run_number`, then delete `input_workspace`"""
        if run_number:
            filename = construct_file_name(run_number)
        SaveNexus(InputWorkspace=str(input_workspace), Filename=str(root_path / filename))
        DeleteWorkspace(str(input_workspace))

    def common_empty_workspace(run_number: Union[int, str] = None, events=True) -> str:
        r"""Create an empty workspace with the appropriate configuration for the instrument, run, and neutrons"""
        # Create an empty events workspace for BIOSANS
        workspace_name = mtd.unique_hidden_name()
        workspace_events = empty_instrument_workspace(
            workspace_name, instrument_name="CG3", event_workspace=events, monitors_have_spectra=False
        )
        workspace_events.getAxis(0).setUnit("TOF")
        workspace_events.getAxis(1).setUnit("Label")
        # Adjust instrument and log instrument configuration
        set_position_south_detector(workspace_name, distance=metadata["sample_detector_distance"])
        set_angle_wing_detector(workspace_name, angle=metadata["angle_wing_detector"])
        set_angle_midrange_detector(workspace_name, angle=metadata["angle_midrange_detector"])
        for log_name in ["sample_aperture_diameter", "source_aperture_diameter", "CG3:CS:SampleToSi"]:
            AddSampleLog(
                Workspace=workspace_name,
                LogType="Number Series",
                NumberType="Double",
                LogName=log_name,
                LogText=str(metadata[log_name]),
                LogUnit="mm",
            )
        # Logs for the run configuration
        sample_logs = SampleLogs(workspace_name)
        if run_number:
            sample_logs.insert("run_number", str(run_number))
        sample_logs.insert("start_time", metadata["start_time"])
        sample_logs.insert("end_time", metadata["end_time"])
        sample_logs.insert("duration", metadata["duration"], unit="second")
        sample_logs.insert("attenuator", metadata["attenuator"])
        # Logs for the neutrons configuration
        sample_logs.insert("monitor", metadata["monitor"])
        for log_name, unit in [("wavelength", "A"), ("wavelength_spread", "A")]:
            AddSampleLog(
                Workspace=workspace_name,
                LogType="Number Series",
                NumberType="Double",
                LogName=log_name,
                LogText=str(metadata[log_name]),
                LogUnit=unit,
            )
        # log the proton charge time-series
        pulse_times = list(np.arange(0.0, metadata["duration"], metadata["pulse_period"] / 1000.0))  # in seconds
        values = [metadata["pulse_charge"]] * len(pulse_times)
        sample_logs.insert_time_series(
            "proton_charge", pulse_times, values, start_time=metadata["start_time"], unit="picoCoulombs"
        )
        return workspace_name

    def _insert_isotropic_events(input_workspace: Union[str, EventWorkspace]):
        r"""
        Insert events impinging an isotropic intensity pattern on the BIOSANS detectors into `input_workspace`
        `max_counts_in_pixel` and `component_efficiencies` are chosen so that each detector pixel receives
        about one neutron per pulse (thus 60 neutrons per pixel in a 1 second runs). This amounts to 130MB file size.
        """
        if events_cache.get("isotropic", None) is None:
            events_cache["isotropic"] = common_empty_workspace()  # initialize the cache
            # Insert events at every prompt-pulse time
            for pulse_time in metadata["pulse_times"]:
                insert_events_isotropic(
                    events_cache["isotropic"],
                    center_x=metadata["CenterX"],
                    center_y=metadata["CenterY"],
                    max_counts_in_pixel=111,
                    components=["detector1", "midrange_detector", "wing_detector"],
                    component_efficiencies=[1.0, 0.25, 0.025],
                    back_panel_attenuation=1.0,  # back-panel tubes receive same intensity as the front-panel ones
                    solid_angle_correction=True,
                    pulse_time=pulse_time,
                    lambda_distribution=metadata.get("wavelength_distribution", None),
                )
        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["isotropic"], OutputWorkspace=input_workspace)

    def _insert_flat_noise(input_workspace: Union[str, EventWorkspace]):
        r"""Insert noise events, at most one event per pixel in any given pulse. Thus at most 60 events per pixel
        for a 1 second run"""
        if events_cache.get("flat_noise", None) is None:
            events_cache["flat_noise"] = common_empty_workspace()
            for pulse_time in metadata["pulse_times"]:
                insert_background(
                    events_cache["flat_noise"],
                    flavor="flat noise",
                    flavor_kwargs=dict(min_counts=0, max_counts=1),
                    pulse_time=pulse_time,
                    lambda_distribution=metadata.get("wavelength_distribution", None),
                )
        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["flat_noise"], OutputWorkspace=input_workspace)

    def _insert_beam_spot(input_workspace: Union[str, EventWorkspace], max_counts_in_pixel: int):
        if events_cache.get("beam_spot", None) is None:
            events_cache["beam_spot"] = common_empty_workspace()
            for pulse_time in metadata["pulse_times"]:
                insert_beam_spot(
                    events_cache["beam_spot"],
                    center_x=metadata["CenterX"],
                    center_y=metadata["CenterY"],
                    diameter=0.015,
                    max_counts_in_pixel=1000,
                    pulse_time=pulse_time,
                    lambda_distribution=metadata.get("wavelength_distribution", None),
                )

        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["beam_spot"], OutputWorkspace=input_workspace)

    # SAMPLE RUN (six time-resolved rings, one isotropic scattering, one flat noise)
    ws_sample = common_empty_workspace(run_number=config["sample"]["runNumber"])
    pulse_count = len(metadata["pulse_times"])
    two_theta_at_max_list = list(metadata["two_theta_at_max"]) * int(pulse_count / len(metadata["two_theta_at_max"]))
    assert len(two_theta_at_max_list) == len(metadata["pulse_times"])
    target_detector_lists = metadata["target_detectors"] * int(pulse_count / len(metadata["target_detectors"]))
    assert len(target_detector_lists) == len(metadata["pulse_times"])

    def _efficiencies(detectors: List[str]):
        r"""tune down the efficiencies the midrange and wing pixels so that we don't have too many events there"""
        effs = {"detector1": 1.0, "midrange_detector": 0.3, "wing_detector": 0.05}
        return [effs[detector] for detector in detectors]

    # Every prompt-pulse, we create a ring pattern with a particular scatteting angle (`two_theta)
    # on selected dectector panels (`target_detectors`)
    for pulse_time, two_theta, target_detectors in zip(
        metadata["pulse_times"], two_theta_at_max_list, target_detector_lists
    ):
        insert_events_ring(
            ws_sample,
            twotheta_center=two_theta,  # in degrees
            twotheta_dev=0.3,
            max_counts_in_pixel=400,
            center_x=metadata["CenterX"],
            center_y=metadata["CenterY"],
            components=target_detectors,
            component_efficiencies=_efficiencies(target_detectors),
            back_panel_attenuation=1.0,  # ignore shading of the back-panel tubes due to the front-panel tubes
            solid_angle_correction=True,
            gravity_correction=True,
            pulse_time=pulse_time,
            lambda_distribution=metadata.get("wavelength_distribution", None),
        )
    _insert_isotropic_events(ws_sample)
    _insert_flat_noise(ws_sample)
    save_and_delete(ws_sample, run_number=config["sample"]["runNumber"])

    # BACKGROUND RUN (one isotropic scattering, one flat noise)
    ws_background = common_empty_workspace(run_number=config["background"]["runNumber"])
    _insert_isotropic_events(ws_background)
    _insert_flat_noise(ws_background)
    save_and_delete(ws_background, run_number=config["background"]["runNumber"])

    # EMPTY TRANSMISSION and BEAM CENTER RUN (one beam spot and one flat noise)
    ws_beam_center = common_empty_workspace(run_number=config["emptyTransmission"]["runNumber"])
    max_allowed_count = 500  # max counts in any pixel of the beam spot
    _insert_beam_spot(ws_beam_center, max_allowed_count)
    _insert_flat_noise(ws_beam_center)
    save_and_delete(ws_beam_center, run_number=config["emptyTransmission"]["runNumber"])
    # we're making the beam center run the same as the empty transmission, so not export again

    # SAMPLE TRANSMISSION AND BACKGROUND TRANSMISSION RUN (one beam spot and one flat noise)
    ws_transmission = common_empty_workspace(run_number=config["sample"]["transmission"]["runNumber"])
    _insert_beam_spot(ws_transmission, int(max_allowed_count * metadata["transmission"]))
    _insert_flat_noise(ws_transmission)
    save_and_delete(ws_transmission, run_number=config["sample"]["transmission"]["runNumber"])
    # we're making the background transmission the same run number than the sample transmission, so no export again


@pytest.mark.datarepo
@pytest.fixture(scope="module")
@namedtuplefy
def six_rings_pattern(datarepo_dir) -> dict:
    r"""
    A set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of six time-resolved rings in the detectors of BIOSANS.
    Read :func:`create_ring_pattern` for more details.

    Files are stored in the data repository under subdirectory :code:`/biosans/simulated_events/six_rings_pattern`

    Returns
    -------
    After decorator `namedtuplefy` is applied to this function, the function returns a namedtuple with the
    following fields:
    - config: a dictionary with the reduction parameters
    - metadata: a dictionary with additional information not encoded in the reduction parameters
    """
    assert Path(datarepo_dir.biosans).is_dir(), f"Data repository {datarepo_dir.biosans} not found"
    datadir = Path(datarepo_dir.biosans) / "simulated_events" / "six_rings_pattern"
    # configuration for the reduction workflow
    config = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "BIOSANS",
        "iptsNumber": "12345",
        "dataDirectories": f"{str(datadir)}",
        "sample": {
            "runNumber": "92310",
            "thickness": 0.2,
            "transmission": {"runNumber": "92330", "value": None},
        },
        "background": {
            "runNumber": "92320",
            "transmission": {"runNumber": "92330", "value": None},
        },
        "emptyTransmission": {"runNumber": "92300", "value": None},
        "beamCenter": {
            "runNumber": "92300",
            "method": "center_of_mass",
            "com_centering_options": {"IntegrationRadius": 0.07},
        },
        "outputFileName": None,
        "configuration": {
            "outputDir": None,
            "wavelength": 2.5,
            "wavelengthSpread": 0.2,
            "useTimeSlice": True,
            "timeSliceOffset": 0.0,
            "timeSlicePeriod": 6.0 / 60,
            "timeSliceInterval": 3.0 / 60,
            "sampleApertureSize": 12.0,
            "usePixelCalibration": False,
            "useDefaultMask": True,
            "defaultMask": [
                {"Pixel": "1-18,239-256"},
            ],
            "darkMainFileName": None,
            "darkWingFileName": None,
            "darkMidrangeFileName": None,
            "normalization": "Monitor",
            "sensitivityMainFileName": None,
            "sensitivityWingFileName": None,
            "sensitivityMidrangeFileName": None,
            "useSolidAngleCorrection": True,
            "useThetaDepTransCorrection": False,
            "DBScalingBeamRadius": 40.0,
            "StandardAbsoluteScale": 1.0,
            "numMainQxQyBins": 100,
            "numMidrangeQxQyBins": 100,
            "numWingQxQyBins": 100,
            "numMainQBins": 20,
            "numMidrangeQBins": 50,
            "numWingQBins": 50,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "QminMain": 0.07,
            "QmaxMain": 0.2,
            "QminMidrange": 0.07,
            "QmaxMidrange": 0.4,
            "QminWing": 0.27,
            "QmaxWing": 0.7,
            "overlapStitchQmin": [0.1, 0.27],
            "overlapStitchQmax": [0.17, 0.4],
            "useErrorWeighting": False,
        },
    }
    wavelength = config["configuration"]["wavelength"]
    # configuration for the instrument, run, and neutrons
    metadata = {
        # instrument configuration
        "sample_detector_distance": 8.0,  # from sample to the main detector, in meters
        "angle_midrange_detector": 2.5,  # degrees away from the beam axis
        "angle_wing_detector": 5.0,  # degrees away from the beam axis
        "sample_aperture_diameter": -1.0,  # in mili-meters
        "source_aperture_diameter": 10.0,  # in mili-meters
        "CG3:CS:SampleToSi": 65.0,  # in mili-meters
        "attenuator": 2,
        "CenterX": -0.008,  # beam spot center, in meters
        "CenterY": -0.016,
        # run configuration
        "start_time": "2023-08-01 00:00:00",
        "end_time": "2023-08-01 00:00:01",  # run finishes just one second after it starts
        "duration": 1.0,  # in seconds, 60 pulses in total
        "transmission": 0.75,  # the only property of the sample
        # neutrons configuration
        "monitor": 4200,  # counts in the monitor
        "wavelength": wavelength,
        "wavelength_spread": config["configuration"]["wavelengthSpread"],
        "wavelength_distribution": lambda events_count: np.repeat(wavelength, events_count),  # monochromatic
        "pulse_period": 1000.0 / 60.0,  # in mili-seconds
        "pulse_charge": 17529040.0,  # in pico-Coulombs
        # scattering angles (in degrees) for the ring patterns
        "two_theta_at_max": np.array(
            [
                2.8,
                7.0,
                11.0,
                3.4,
                8.5,
                13.0,
            ]
        ),
        "target_detectors": [
            ["detector1", "midrange_detector"],  # receives a ring intensity pattern at two_theta 2.8
            ["midrange_detector", "wing_detector"],  # receive a ring intensity pattern at two_theta 7.0
            ["wing_detector"],  # receive a ring intensity pattern at two_theta 11.0
            ["detector1", "midrange_detector"],  # receives a ring intensity pattern at two_theta 3.4
            ["midrange_detector", "wing_detector"],  # receive a ring intensity pattern at two_theta 8.0
            ["wing_detector"],  # receive a ring intensity pattern at two_theta 13.0
        ],
    }

    # Insert a list of prompt-pulse DateAndTime instances in the `metadata` dictionary
    first_pulse_time, end_time = DateAndTime(metadata["start_time"]), DateAndTime(metadata["end_time"])
    pulse_period = metadata["pulse_period"] * 1000000  # from mili-seconds to nano-seconds
    current_pulse_time, pulse_count, metadata["pulse_times"] = first_pulse_time, 0, list()
    while current_pulse_time < end_time:
        metadata["pulse_times"].append(current_pulse_time)
        pulse_count += 1
        current_pulse_time = first_pulse_time + round(pulse_count * pulse_period)

    # Four our monochromatic wavelength distribution with wavelength `W`, the scattering ring patterns have
    # intensity maxima at scattering angles stored in `metadata[two_theta_at_max]`.
    # The corresponding momentum transfer modulus Q are
    # `Q = 2 * k_i * sin(two_theta_at_max / 2)`, with `k_i = 2 * pi / W` the incident wave vector.
    theta_at_max = np.radians(metadata["two_theta_at_max"]) / 2.0
    metadata["Q_at_max_I"] = 2 * (2 * np.pi / wavelength) * np.sin(theta_at_max)  # store these values

    # Generate the simulated Nexus event files if they're not in the database
    if missing_files(
        datadir,
        run_numbers=["92300", "92310", "92320", "92330"],
    ):
        create_six_rings_pattern(config, metadata)

    return dict(config=config, metadata=metadata)


def _mock_LoadEventAsWorkspace2D(*_, **kwargs):
    # Substitute LoadEventAsWorkspace2D with LoadNexusProcessed because our synthetic files were created with SaveNexus
    ws = LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])
    ws = transform_to_wavelength(ws)
    return ws


def _mock_LoadEventNexus(*_, **kwargs):
    # Substitute LoadEventNexus with LoadNexusProcessed because our synthetic files were created with SaveNexus
    ws = LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])
    return ws


@pytest.mark.datarepo
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
@mock_patch("drtsans.load.LoadEventAsWorkspace2D", new=_mock_LoadEventAsWorkspace2D)
def test_split_six_rings(six_rings_pattern: dict, temp_directory: Callable[[Any], str]):
    r"""
    Scattering from a simulated sample leaving an intensity pattern in the shape of six time-resolved rings
    in the detectors of BIOSANS.
    - The duration of the simulated run is 1 second. Thus, there are 60 prompt-pulses in the run
    - We split the 1 second run in chunks of duration 1/20 seconds, thus we have 20 chunks. In the first chunk,
      we simulate an intensity pattern of three rings at scattering angles  2.8, 7.0, 11.0. In the second chunk,
      we simulate an intensity pattern of three rings at scattering angles 3.4, 8.5, 13.0. Subsequent chunks
      alternate between the [2.8, 7.0, 11.0] and [3.4, 8.5, 13.0] patterns.
    - We use reduction parameters  `"timeSlicePeriod": 6.0 / 60, "timeSliceInterval": 3.0 / 60,`. This means we'll
      collect and reduce together all the [2.8, 7.0, 11.0] chunks. Similarly for the [3.4, 8.5, 13.0] chunks.
    - The ring at two_theta=2.8 is detected by `detector1` and `midrange_detector`. The ring at two_theta=7.0 is
      detected by `midrange_detector` and `wing_detector`. This means we'll be able to stitch together the intensity
      profiles obtained at each detector. Similarly, the ring at two_theta=3.4 is is detected by `detector1` and
      `midrange_detector`, and the ring at two_theta=8.5 is detected by `midrange_detector` and `wing_detector`.
    """
    # pad missing parameters with default values from schema EQSANS.json
    config = reduction_parameters(six_rings_pattern.config, "BIOSANS", validate=False)

    # insert parameters customized for the current test
    sample_run_number = config["sample"]["runNumber"]
    amendments = {
        "outputFileName": f"CG2_{sample_run_number}",
        "configuration": {
            "outputDir": temp_directory(prefix="testBIOSANSSplitSixRings_"),
        },
    }
    config = update_reduction_parameters(config, amendments, validate=True)
    metadata = six_rings_pattern.metadata

    def _mock_monitor_split_and_log(
        monitor: str, monitor_group: str, sample_group: str, is_mono: bool, filter_events: dict
    ):
        r"""
        Divide the total number of monitor counts by the number of time slices, then insert the resulting value
         as log entry 'monitor' in the corresponding splited sample workspaces

        Parameters
        ----------
        sample_group
            name of the WorkspaceGroup containing the splited workspaces for the sample
        """
        splitted_workspaces_count = mtd[sample_group].getNumberOfEntries()
        monitor_count_per_slice = metadata["monitor"] / splitted_workspaces_count
        for n in range(splitted_workspaces_count):
            SampleLogs(mtd[sample_group].getItem(n)).insert("monitor", monitor_count_per_slice)

    # load all necessary files
    with mock_patch("drtsans.load._monitor_split_and_log", side_effect=_mock_monitor_split_and_log):
        loaded = load_all_files(config, path=config["dataDirectories"])

    # do the actual reduction
    reduction_output = reduce_single_configuration(loaded, config)

    ##############
    # Find if the peaks of the different I(Q) signal lie near the expected Q values
    ##############

    def fetch_profile(time_slice="first", profile="I1D_main"):
        slice_index = 0 if time_slice == "first" else 1
        iqmod = getattr(reduction_output[slice_index], profile)[0]
        intensities = iqmod.intensity
        intensities[np.isnan(intensities)] = 0.0
        return iqmod.mod_q, intensities

    # I(Q) for the main detector, first time slice
    q_values, intensities = fetch_profile(time_slice="first", profile="I1D_main")
    expected = metadata["Q_at_max_I"][0]  # expected Q-value for the peak of I(Q) profile I1D_main
    q_at_max = q_values[np.argmax(intensities)]
    assert abs(expected - q_at_max) < 0.1 * expected

    # I(Q) for the midrange detector, first time slice
    q_values, intensities = fetch_profile(time_slice="first", profile="I1D_midrange")
    # first peak
    expected = metadata["Q_at_max_I"][0]  # expected Q-value for the first peak of I(Q) profile I1D_midrange
    q_at_max = q_values[np.argmax(intensities[0:20])]  # the first peak is within the first 20 Q-bins
    assert abs(expected - q_at_max) < 0.1 * expected
    # second peak
    expected = metadata["Q_at_max_I"][1]  # expected Q-value for the second peak of I(Q) profile I1D_midrange
    q_at_max = q_values[20 + np.argmax(intensities[20:])]  # second peak located after the first 20 Q-bins
    assert abs(expected - q_at_max) < 0.1 * expected

    # I(Q) for the wing detector, first time slice
    q_values, intensities = fetch_profile(time_slice="first", profile="I1D_wing")
    # first peak
    expected = metadata["Q_at_max_I"][1]  # expected Q-value for the first peak of I(Q) profile I1D_wing
    q_at_max = q_values[np.argmax(intensities[0:20])]
    assert abs(expected - q_at_max) < 0.1 * expected
    # second peak
    expected = metadata["Q_at_max_I"][2]  # expected Q-value for the second peak of I(Q) profile I1D_wing
    q_at_max = q_values[20 + np.argmax(intensities[20:])]
    assert abs(expected - q_at_max) < 0.1 * expected

    # I(Q) for the main detector, second time slice
    q_values, intensities = fetch_profile(time_slice="second", profile="I1D_main")
    expected = metadata["Q_at_max_I"][3]
    q_at_max = q_values[np.argmax(intensities)]
    assert abs(expected - q_at_max) < 0.1 * expected

    # I(Q) for the midrange detector, second time slice
    q_values, intensities = fetch_profile(time_slice="second", profile="I1D_midrange")
    # first peak
    expected = metadata["Q_at_max_I"][3]
    q_at_max = q_values[np.argmax(intensities[0:20])]
    assert abs(expected - q_at_max) < 0.1 * expected
    # second peak
    expected = metadata["Q_at_max_I"][4]
    q_at_max = q_values[20 + np.argmax(intensities[20:])]
    assert abs(expected - q_at_max) < 0.1 * expected

    # I(Q) for the wing detector, second time slice
    q_values, intensities = fetch_profile(time_slice="second", profile="I1D_wing")
    # first peak
    expected = metadata["Q_at_max_I"][4]
    q_at_max = q_values[np.argmax(intensities[0:20])]
    assert abs(expected - q_at_max) < 0.1 * expected
    # second peak
    expected = metadata["Q_at_max_I"][5]
    q_at_max = q_values[20 + np.argmax(intensities[20:])]
    assert abs(expected - q_at_max) < 0.1 * expected


if __name__ == "__main__":
    pytest.main([__file__])
