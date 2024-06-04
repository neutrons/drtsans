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
    LoadNexusProcessed,
    MoveInstrumentComponent,
    mtd,
    Plus,
    SaveNexus,
)
import numpy as np
import pytest

# drtsans imports
from drtsans.dataobjects import IQmod
from drtsans.instruments import empty_instrument_workspace
from drtsans.mono.load import transform_to_wavelength
from drtsans.mono.gpsans import (
    load_all_files,
    reduce_single_configuration,
    reduction_parameters,
    update_reduction_parameters,
)
from drtsans.samplelogs import SampleLogs
from drtsans.settings import namedtuplefy
from drtsans.simulated_events import insert_background, insert_beam_spot, insert_events_isotropic, insert_events_ring


def construct_file_name(run_number: Union[str, int]) -> str:
    r"""Convention in this module to construct a file name from the run number"""
    return f"CG2_{run_number}.nxs"


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


def create_three_rings_pattern(config: dict, metadata: dict):
    r"""
    Create a set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of three rings in the main detector of GPSANS, plus additional files to simulate the:
    - transmission of the sample
    - background run (run without the sample)
    - background transmission run
    - empty transmission run

    Parameters
    ----------
    config
        Reduction configuration, usually determined by the instrument scientist and the user
    metadata
        additional configuration, usually stored in the logs of the Nexus files, or encoded in the instrument geometry
        and the neutron pulse.

    """
    # verify directory exists or can be created if it does not exist
    root_path = Path(config["dataDirectories"][0])
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
        r"""Create an empty workspace with the correct instrument definition file and metadata"""
        # Create an empty events workspace for GPSANS with the detector at a certain position
        workspace_name = mtd.unique_hidden_name()
        workspace_events = empty_instrument_workspace(
            workspace_name, instrument_name="CG2", event_workspace=events, monitors_have_spectra=False
        )
        workspace_events.getAxis(0).setUnit("TOF")
        workspace_events.getAxis(1).setUnit("Label")
        MoveInstrumentComponent(
            Workspace=workspace_name,
            ComponentName="detector1",
            Z=metadata["sample_detector_distance"],
            RelativePosition=False,
        )
        # Insert sample logs
        sample_logs = SampleLogs(workspace_name)
        if run_number:
            sample_logs.insert("run_number", str(run_number))
        sample_logs.insert("attenuator", metadata["attenuator"])
        sample_logs.insert("duration", metadata["duration"], unit="second")
        sample_logs.insert("start_time", metadata["start_time"])
        sample_logs.insert("end_time", metadata["end_time"])
        sample_logs.insert("monitor", metadata["monitor"])

        for config_name, log_name, unit in [
            ("sampleApertureSize", "sample_aperture_diameter", "mm"),
            ("sourceApertureDiameter", "source_aperture_diameter", "mm"),
            ("wavelength", "wavelength", "A"),
            ("wavelengthSpread", "wavelength_spread", "A"),
        ]:
            AddSampleLog(
                Workspace=workspace_name,
                LogType="Number Series",
                NumberType="Double",
                LogName=log_name,
                LogText=str(config["configuration"][config_name]),
                LogUnit=unit,
            )

        for meta_name, log_name, unit in [
            ("sample_detector_distance", "sample_detector_distance", "m"),
            ("sample_to_Si", "CG2:CS:SampleToSi", "mm"),
        ]:
            AddSampleLog(
                Workspace=workspace_name,
                LogType="Number Series",
                NumberType="Double",
                LogName=log_name,
                LogText=str(metadata[meta_name]),
                LogUnit=unit,
            )

        # insert proton charge
        times = list(np.arange(0.0, metadata["duration"], metadata["pulse_period"] / 1000.0))  # in seconds
        values = [metadata["pulse_charge"]] * len(times)
        sample_logs.insert_time_series(
            "proton_charge", times, values, start_time=metadata["start_time"], unit="picoCoulombs"
        )
        return workspace_name

    def _insert_isotropic_events(input_workspace: Union[str, EventWorkspace]):
        if events_cache.get("isotropic", None) is None:
            events_cache["isotropic"] = common_empty_workspace()
            for pulse_time in metadata["pulse_times"]:
                insert_events_isotropic(
                    events_cache["isotropic"],
                    center_x=metadata["CenterX"],
                    center_y=metadata["CenterY"],
                    max_counts_in_pixel=25,
                    back_panel_attenuation=1.0,
                    solid_angle_correction=True,
                    pulse_time=pulse_time,
                    lambda_distribution=metadata.get("wavelength_distribution", None),
                )
        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["isotropic"], OutputWorkspace=input_workspace)

    def _insert_flat_noise(input_workspace: Union[str, EventWorkspace]):
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
                    max_counts_in_pixel=max_counts_in_pixel,
                    pulse_time=pulse_time,
                    lambda_distribution=metadata.get("wavelength_distribution", None),
                )

        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["beam_spot"], OutputWorkspace=input_workspace)

    # SAMPLE RUN (three time-resolved rings, one isotropic scattering, one flat noise)
    ws_sample = common_empty_workspace(run_number=config["sample"]["runNumber"])
    for pulse_time, two_theta in zip(
        metadata["pulse_times"], np.tile(metadata["two_theta_at_max"], len(metadata["pulse_times"]))
    ):
        insert_events_ring(
            ws_sample,
            twotheta_center=two_theta,  # in degrees
            twotheta_dev=0.5,
            max_counts_in_pixel=50,
            center_x=metadata["CenterX"],
            center_y=metadata["CenterY"],
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
    transmission_value = 0.75  # expected ratio between the transmission (with the sample) and the empty-transmission
    _insert_beam_spot(ws_transmission, int(max_allowed_count * transmission_value))
    _insert_flat_noise(ws_transmission)
    save_and_delete(ws_transmission, run_number=config["sample"]["transmission"]["runNumber"])
    # we're making the background transmission the same run number than the sample transmission, so no export again

    return None


@pytest.mark.datarepo
@pytest.fixture(scope="module")
@namedtuplefy
def three_rings_pattern(datarepo_dir) -> dict:
    r"""
    A set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of three time-resolved rings in the main detector of GPSANS.
    Read :func:`create_ring_pattern` for more details.

    Files are stored in the data repository under subdirectory :code:`/gpsans/simulated_events/three_rings_pattern`

    Returns
    -------
    After decorator `namedtuplefy` is applied to this function, the function returns a namedtuple with the
    following fields:
    - config: a dictionary with the reduction parameters
    - metadata: a dictionary with additional information not encoded in the reduction parameters
    """
    assert Path(datarepo_dir.gpsans).is_dir(), f"Data repository {datarepo_dir.gpsans} not found"
    datadir = Path(datarepo_dir.gpsans) / "simulated_events" / "three_rings_pattern"
    # configuration for the reduction workflow
    config = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "GPSANS",
        "iptsNumber": "12345",
        "dataDirectories": [str(datadir)],
        "sample": {"runNumber": "92310", "thickness": "0.1", "transmission": {"runNumber": "92330", "value": ""}},
        "background": {"runNumber": "92320", "transmission": {"runNumber": "92330", "value": ""}},
        "emptyTransmission": {"runNumber": "92300", "value": ""},
        "beamCenter": {
            "runNumber": "92300",
            "method": "gaussian",
            "gaussian_centering_options": {"theta": {"value": 0.0, "vary": False}},
        },
        "outputFileName": None,
        "configuration": {
            "outputDir": None,
            "wavelength": 2.5,
            "wavelengthSpread": 0.2,
            "sampleOffset": "",
            "useDetectorOffset": False,
            "sampleApertureSize": 4.0,  # in mm
            "sourceApertureDiameter": 40.0,  # in mm
            "usePixelCalibration": False,
            "maskFileName": "",
            "useDefaultMask": True,
            "defaultMask": [{"Pixel": "1-10,247-256"}],
            "useMaskBackTubes": False,
            "darkFileName": None,
            "normalization": "Monitor",
            "sensitivityFileName": None,
            "useSolidAngleCorrection": True,
            "useThetaDepTransCorrection": False,
            "DBScalingBeamRadius": "40.0",
            "mmRadiusForTransmission": "",
            "absoluteScaleMethod": "direct_beam",
            "StandardAbsoluteScale": "1.0",
            "numQxQyBins": 192,
            "1DQbinType": "scalar",
            "QbinType": "log",
            "numQBins": 100,
            "LogQBinsPerDecade": "",
            "useLogQBinsDecadeCenter": "",
            "useLogQBinsEvenDecade": False,
            "Qmin": "",
            "Qmax": "",
            "useErrorWeighting": False,
            "smearingPixelSizeX": "",
            "smearingPixelSizeY": "",
            "useSubpixels": False,
        },
    }
    wavelength = config["configuration"]["wavelength"]
    # configuration for the neutrons, instruments and experiment
    metadata = {
        "sample_detector_distance": 5.0,  # in meters
        "sample_to_Si": 1.0,  # in mm
        "monitor": 3.0e9,  # assume all files have the same number of counts
        "attenuator": 2,
        "wavelength_distribution": lambda events_count: np.repeat(wavelength, events_count),  # monochromatic
        "pulse_period": 1000.0 / 60.0,  # in mili-seconds
        "start_time": "2023-08-01 00:00:00",
        "end_time": "2023-08-01 00:00:01",  # 60 pulses in total
        "pulse_charge": 17529040.0,  # in pico-Coulombs
        "duration": 1.0,  # in seconds,
        "CenterX": 0.025239,
        "CenterY": 0.0170801,
        "two_theta_at_max": np.array(
            [1.0, 3.0, 5.0]
        ),  # in degrees, scattering angles with maximum scattered intensity
    }

    # create the list of pulse times
    pulse_start, end_time = [DateAndTime(metadata[time]) for time in ["start_time", "end_time"]]
    pulse_period = metadata["pulse_period"] * 1000000  # from mili-seconds to nano-seconds
    pulse_current, pulse_count, metadata["pulse_times"] = pulse_start, 0, list()
    while pulse_current < end_time:
        metadata["pulse_times"].append(pulse_current)
        pulse_count += 1
        pulse_current = pulse_start + round(pulse_count * pulse_period)

    # Four our monochromatic wavelength distribution with wavelength `W`, the scattering ring pattern has a maximum at
    # scattering angle `two_theta_at_max`, thus the momentum transfer modulus Q corresponding to this maximum
    # is Q = 2 * k_i * sin(two_theta_at_max / 2), with k_i = 2 * pi / W
    theta_at_max = np.radians(metadata["two_theta_at_max"]) / 2.0
    metadata["Q_at_max_I"] = 2 * (2 * np.pi / wavelength) * np.sin(theta_at_max)

    if missing_files(
        datadir,
        run_numbers=["92300", "92310", "92320", "92330"],
    ):
        create_three_rings_pattern(config, metadata)

    return dict(config=config, metadata=metadata)


def _mock_LoadEventAsWorkspace2D(*_, **kwargs):
    ws = LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])
    ws = transform_to_wavelength(ws)
    return ws


def _mock_LoadEventNexus(*_, **kwargs):
    ws = LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])
    return ws


@pytest.mark.datarepo
@mock_patch("drtsans.load.LoadEventAsWorkspace2D", new=_mock_LoadEventAsWorkspace2D)
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
def test_split_three_rings(three_rings_pattern: dict, temp_directory: Callable[[Any], str]):
    r"""
    Reduce a single configuration that simulates scattering from a sample which imprints an intensity pattern in the
    shape of three time-resolved rings in the main detector of GPSANS.

    Parameters
    ----------
    three_rings_pattern
        Fixture evaluating to a reduction configuration, containing all the necessary run numbers
        and reduction parameters
    temp_directory
        Fixture evaluating to a function. When invoked, returns the path to a directory which will be erased upon
        completion of the test, irrespective of the test being successful or not.
    """
    # pad missing parameters with default values from schema EQSANS.json
    config = reduction_parameters(three_rings_pattern.config, "GPSANS", validate=False)

    # insert parameters customized for the current test
    sample_run_number = config["sample"]["runNumber"]
    amendments = {
        "outputFileName": f"CG2_{sample_run_number}",
        "configuration": {
            "outputDir": temp_directory(prefix="testGPSANSSplitThreeRings_"),
            "useTimeSlice": True,
            "timeSliceInterval": 1.0 / 60,
            "timeSliceOffset": 0.0,
            "timeSlicePeriod": 3.0 / 60,
        },
    }
    config = update_reduction_parameters(config, amendments, validate=True)
    metadata = three_rings_pattern.metadata

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

    # find the Q-modulus where the scattering intensity is maximum, and compare to what's expected
    minimum_peak_intensity = 800.0  # all three peaks have a maximum intensity bigger than this number
    for peak_index, q_at_max_i in enumerate(metadata["Q_at_max_I"]):
        i_vs_qmod: IQmod = reduction_output[peak_index].I1D_main[0]  # 1D intensity profile
        closest_index = np.argmin(np.abs(i_vs_qmod.mod_q - q_at_max_i))
        assert i_vs_qmod.intensity[closest_index] > minimum_peak_intensity


if __name__ == "__main__":
    pytest.main([__file__])
