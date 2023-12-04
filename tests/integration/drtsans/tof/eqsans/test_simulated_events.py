# standard library imports
from pathlib import Path
from typing import Any, Callable, List, Tuple, Union
from unittest.mock import patch as mock_patch

# third party imports
from mantid.simpleapi import (
    AddSampleLog,
    ConvertUnits,
    DeleteWorkspace,
    LoadNexusProcessed,
    MoveInstrumentComponent,
    mtd,
    Plus,
    Rebin,
    SaveNexus,
)
from mantid.dataobjects import EventWorkspace
from mantid.kernel import DateAndTime
import numpy as np
from numpy.testing import assert_allclose
import pytest

# local imports
from drtsans.dataobjects import IQmod
from drtsans.instruments import empty_instrument_workspace
from drtsans.load import load_events as generic_load_events
from drtsans.load import load_and_split as generic_load_and_split
from drtsans.samplelogs import SampleLogs
from drtsans.settings import amend_config, namedtuplefy, unique_workspace_dundername
from drtsans.simulated_events import insert_background, insert_beam_spot, insert_events_isotropic, insert_events_ring
from drtsans.tof.eqsans import (
    load_all_files,
    reduce_single_configuration,
    reduction_parameters,
    update_reduction_parameters,
)


def _histogram_all_events(
    input_workspace: Union[str, EventWorkspace],
    units: str = "Wavelength",
    binning: str = "0.1,0.1,30",  # first bin boundary, bin width, last bin boundary
) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Convert events from time-of-flight to selected units, then histogram each spectrum and add all histograms, thus
    ending with just one histogram for all events.

    Parameters
    ----------
    input_workspace
    output_workspace
    units
    binning

    Returns
    -------
    Tuple with two items, representing the histogram's bin boundaries and the histogram intensities, respectively
    """
    import mantid.simpleapi as mapi

    temp_workspace = unique_workspace_dundername()  # a name not taken by any other already existing workspace
    mapi.ConvertUnits(InputWorkspace=input_workspace, OutputWorkspace=temp_workspace, Target=units)
    mapi.Rebin(InputWorkspace=temp_workspace, OutputWorkspace=temp_workspace, Params=binning, PreserveEvents=False)
    mapi.SumSpectra(InputWorkspace=temp_workspace, OutputWorkspace=temp_workspace)
    bins, intensities = mapi.mtd[temp_workspace].extractX()[0], mapi.mtd[temp_workspace].extractY()[0]
    DeleteWorkspace(temp_workspace)
    return bins, intensities


def construct_file_name(run_number: Union[str, int]) -> str:
    r"""Convention in this module to construct a file name from the run number"""
    return f"EQSANS_{run_number}.nxs"


def missing_files(root_path: Path, run_numbers: List[str], file_names: List[str]) -> bool:
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


def create_ring_pattern(config: dict, metadata: dict):
    r"""
    Create a set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of a ring in the main detector of EQSANS, plus additional files to simulate the:
    - transmission of the sample
    - background run (run without the sample)
    - background transmission run
    - empty transmission run
    - dark run

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

    # create only one flat noise workspace for all events
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
        # Create an empty events workspace for EQSANS with the detector at a certain position
        workspace_name = unique_workspace_dundername()
        workspace_events = empty_instrument_workspace(
            workspace_name, instrument_name="EQSANS", event_workspace=events, monitors_have_spectra=False
        )
        workspace_events.getAxis(0).setUnit("TOF")
        workspace_events.getAxis(1).setUnit("Label")
        MoveInstrumentComponent(
            Workspace=workspace_name, ComponentName="detector1", Z=metadata["Z"], RelativePosition=False
        )
        # Insert sample logs
        sample_logs = SampleLogs(workspace_name)
        if run_number:
            sample_logs.insert("run_number", str(run_number))
        sample_logs.insert("duration", metadata["duration"], unit="second")
        sample_logs.insert("start_time", metadata["start_time"])
        sample_logs.insert("end_time", metadata["end_time"])
        AddSampleLog(
            Workspace=workspace_name,
            LogType="Number Series",
            NumberType="Double",
            LogName="wavelength",
            LogText=str(metadata["wavelength"]),
            LogUnit="A",
        )
        for slit in ["vBeamSlit", "vBeamSlit2", "vBeamSlit3"]:
            AddSampleLog(
                Workspace=workspace_name,
                LogType="Number Series",
                NumberType="Int",
                LogName=slit,
                LogText=str(metadata[slit]),
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
            insert_events_isotropic(
                events_cache["isotropic"],
                center_x=metadata["CenterX"],
                center_y=metadata["CenterY"],
                max_counts_in_pixel=100,
                back_panel_attenuation=1.0,
                solid_angle_correction=True,
                pulse_time=metadata.get("pulse_time", None),
                lambda_distribution=metadata.get("wavelength_distribution", None),
            )
        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["isotropic"], OutputWorkspace=input_workspace)

    def _insert_flat_noise(input_workspace: Union[str, EventWorkspace]):
        if events_cache.get("flat_noise", None) is None:
            events_cache["flat_noise"] = common_empty_workspace()
            insert_background(
                events_cache["flat_noise"],
                flavor="flat noise",
                flavor_kwargs=dict(min_counts=0, max_counts=2),
                pulse_time=metadata.get("pulse_time", None),
                lambda_distribution=metadata.get("wavelength_distribution", None),
            )
        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["flat_noise"], OutputWorkspace=input_workspace)

    def _insert_beam_spot(input_workspace: Union[str, EventWorkspace], max_counts_in_pixel: int):
        if events_cache.get("beam_spot", None) is None:
            events_cache["beam_spot"] = common_empty_workspace()
            insert_beam_spot(
                events_cache["beam_spot"],
                center_x=metadata["CenterX"],
                center_y=metadata["CenterY"],
                diameter=0.015,
                max_counts_in_pixel=max_counts_in_pixel,
                pulse_time=metadata.get("pulse_time", None),
                lambda_distribution=metadata.get("wavelength_distribution", None),
            )
        Plus(LHSWorkspace=input_workspace, RHSWorkspace=events_cache["beam_spot"], OutputWorkspace=input_workspace)

    # SAMPLE RUN (one ring, one isotropic scattering, one flat noise)
    ws_sample = common_empty_workspace(run_number=config["sample"]["runNumber"])
    insert_events_ring(
        ws_sample,
        twotheta_center=metadata["two_theta_at_max"],  # in degrees
        twotheta_dev=0.5,
        max_counts_in_pixel=200,
        center_x=metadata["CenterX"],
        center_y=metadata["CenterY"],
        back_panel_attenuation=1.0,  # ignore shading of the back-panel tubes due to the front-panel tubes
        solid_angle_correction=True,
        gravity_correction=True,
        pulse_time=metadata.get("pulse_time", None),
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
    max_allowed_count = 1000  # max counts in any pixel of the beam spot
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


@pytest.fixture(scope="module")
@namedtuplefy
def ring_pattern(datarepo_dir) -> dict:
    r"""
    A set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of a ring in the main detector of EQSANS. Read :func:`create_ring_pattern` for more details.

    Files are stored in the data repository under subdirectory :code:`/eqsans/simulated_events/ring_pattern`

    Returns
    -------
    After decorator `namedtuplefy` is applied to this function, the function returns a namedtuple with the
    following fields:
    - config: a dictionary with the reduction parameters
    - metadata: a dictionary with additional information not encoded in the reduction parameters
    """
    assert Path(datarepo_dir.eqsans).is_dir(), f"Data repository {datarepo_dir.eqsans} not found"
    datadir = Path(datarepo_dir.eqsans) / "simulated_events" / "ring_pattern"
    config = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "EQSANS",
        "dataDirectories": [str(datadir)],
        "iptsNumber": "12345",
        "sample": {
            "runNumber": "92310",
            "thickness": 0.1,
            "transmission": {"runNumber": "92330", "value": None},
        },
        "background": {"runNumber": "92320", "transmission": {"runNumber": "92330", "value": None}},
        "emptyTransmission": {"runNumber": "92300", "value": None},
        "beamCenter": {
            "runNumber": "92300",
            "method": "center_of_mass",
            "com_centering_options": {"IntegrationRadius": 0.07},
        },
        "configuration": {
            "instrumentConfigurationDir": str(Path(datarepo_dir.eqsans) / "instrument_configuration"),
            "maskFileName": str(Path(datadir) / "beamstop60_mask_4m.nxs"),
            "beamFluxFileName": str(Path(datarepo_dir.eqsans) / "instrument_configuration" / "bl6_flux_at_sample.dat"),
            "darkFileName": str(Path(datadir) / "empty_dark_file.nxs"),
            "sensitivityFileName": str(Path(datarepo_dir.eqsans) / "sensitivities" / "default_sensitivity.nxs"),
            "fluxMonitorRatioFile": None,
            "cutTOFmin": 0.0,
            "cutTOFmax": 0.0,
            "wavelengthStepType": "constant Delta lambda",
            "wavelengthStep": 0.1,
            "useDefaultMask": True,
            "sampleApertureSize": 30,
            "StandardAbsoluteScale": 1.03,
            "sampleOffset": 0,
            "normalization": "Total charge",
            "absoluteScaleMethod": "standard",
            "detectorOffset": 0,
            "mmRadiusForTransmission": 25,
            "numQxQyBins": 80,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "useErrorWeighting": False,
            "numQBins": 120,
            "AnnularAngleBin": 5,
        },
    }
    wavelength = 2.5  # wavelength of the incident beam neutrons, in Angstroms
    metadata = {
        "wavelength": wavelength,
        "wavelength_min": wavelength - 1.0,
        "wavelength_max": wavelength + 1.0,
        "wavelength_distribution": lambda events_count: np.repeat(wavelength, events_count),  # monochromatic
        "Z": 5.0,  # distance along the nominal Z-axis from the origin of coords. to the main detector, in meters
        "start_time": "2023-08-01 00:00:00",
        "end_time": "2023-08-01 00:01:00",
        "pulse_period": 1000.0 / 60.0,  # in mili-seconds
        "pulse_charge": 17529040.0,  # in pico-Coulombs
        "duration": 60.0,  # in seconds,
        "CenterX": 0.025239,
        "CenterY": 0.0170801,
        "vBeamSlit": 5,
        "vBeamSlit2": 5,
        "vBeamSlit3": 5,
        "two_theta_at_max": 3.0,  # in degrees, scattering angle with maximum scattered intensity
    }
    # Four our monochromatic wavelength distribution with wavelength `W`, the scattering ring pattern has a maximum at
    # scattering angle `two_theta_at_max`, thus the momentum transfer modulus Q corresponding to this maximum
    # is Q = 2 * k_i * sin(two_theta_at_max / 2), with k_i = 2 * pi / W
    theta_at_max = metadata["two_theta_at_max"] * np.pi / 360.0  # in radians
    metadata["Q_at_max_I"] = 2 * (2 * np.pi / metadata["wavelength"]) * np.sin(theta_at_max)

    if missing_files(
        datadir,
        run_numbers=["92300", "92310", "92320", "92330"],
        file_names=["beamstop60_mask_4m.nxs", "empty_dark_file.nxs"],
    ):
        create_ring_pattern(config, metadata)
    return dict(config=config, metadata=metadata)


def create_three_rings_pattern(config: dict, metadata: dict):
    r"""
    Create a set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of a ring in the main detector of EQSANS, plus additional files to simulate the:
    - transmission of the sample
    - background run (run without the sample)
    - background transmission run
    - empty transmission run
    - dark run

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

    # create only one flat noise workspace for all events
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
        # Create an empty events workspace for EQSANS with the detector at a certain position
        workspace_name = unique_workspace_dundername()
        workspace_events = empty_instrument_workspace(
            workspace_name, instrument_name="EQSANS", event_workspace=events, monitors_have_spectra=False
        )
        workspace_events.getAxis(0).setUnit("TOF")
        workspace_events.getAxis(1).setUnit("Label")
        MoveInstrumentComponent(
            Workspace=workspace_name, ComponentName="detector1", Z=metadata["Z"], RelativePosition=False
        )
        # Insert sample logs
        sample_logs = SampleLogs(workspace_name)
        if run_number:
            sample_logs.insert("run_number", str(run_number))
        sample_logs.insert("duration", metadata["duration"], unit="second")
        sample_logs.insert("start_time", metadata["start_time"])
        sample_logs.insert("end_time", metadata["end_time"])
        AddSampleLog(
            Workspace=workspace_name,
            LogType="Number Series",
            NumberType="Double",
            LogName="wavelength",
            LogText=str(metadata["wavelength"]),
            LogUnit="A",
        )
        for slit in ["vBeamSlit", "vBeamSlit2", "vBeamSlit3"]:
            AddSampleLog(
                Workspace=workspace_name,
                LogType="Number Series",
                NumberType="Int",
                LogName=slit,
                LogText=str(metadata[slit]),
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
                    max_counts_in_pixel=50,
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

    # SAMPLE RUN (one ring, one isotropic scattering, one flat noise)
    ws_sample = common_empty_workspace(run_number=config["sample"]["runNumber"])

    for pulse_time, two_theta in zip(
        metadata["pulse_times"], np.tile(metadata["two_theta_at_max"], len(metadata["pulse_times"]))
    ):
        insert_events_ring(
            ws_sample,
            twotheta_center=two_theta,  # in degrees
            twotheta_dev=0.5,
            max_counts_in_pixel=100,
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
    _, intensities = _histogram_all_events(ws_transmission, units="Wavelength", binning="0.0,1.0,5.0")
    save_and_delete(ws_transmission, run_number=config["sample"]["transmission"]["runNumber"])
    # we're making the background transmission the same run number than the sample transmission, so no export again

    return None


@pytest.fixture(scope="module")
@namedtuplefy
def three_rings_pattern(datarepo_dir) -> dict:
    r"""
    A set of processed event Nexus files simulating the scattering from a sample that leaves an intensity
    pattern in the shape of a ring in the main detector of EQSANS. Read :func:`create_ring_pattern` for more details.

    Files are stored in the data repository under subdirectory :code:`/eqsans/simulated_events/ring_pattern`

    Returns
    -------
    After decorator `namedtuplefy` is applied to this function, the function returns a namedtuple with the
    following fields:
    - config: a dictionary with the reduction parameters
    - metadata: a dictionary with additional information not encoded in the reduction parameters
    """
    assert Path(datarepo_dir.eqsans).is_dir(), f"Data repository {datarepo_dir.eqsans} not found"
    datadir = Path(datarepo_dir.eqsans) / "simulated_events" / "three_rings_pattern"
    config = {
        "schemaStamp": "2020-04-15T21:09:52.745905",
        "instrumentName": "EQSANS",
        "dataDirectories": [str(datadir)],
        "iptsNumber": "12345",
        "sample": {
            "runNumber": "92310",
            "thickness": 0.1,
            "transmission": {"runNumber": "92330", "value": None},
        },
        "background": {"runNumber": "92320", "transmission": {"runNumber": "92330", "value": None}},
        "emptyTransmission": {"runNumber": "92300", "value": None},
        "beamCenter": {
            "runNumber": "92300",
            "method": "center_of_mass",
            "com_centering_options": {"IntegrationRadius": 0.07},
        },
        "configuration": {
            "instrumentConfigurationDir": str(Path(datarepo_dir.eqsans) / "instrument_configuration"),
            "maskFileName": str(Path(datadir) / "beamstop60_mask_4m.nxs"),
            "beamFluxFileName": str(Path(datarepo_dir.eqsans) / "instrument_configuration" / "bl6_flux_at_sample.dat"),
            "darkFileName": str(Path(datadir) / "empty_dark_file.nxs"),
            "sensitivityFileName": str(Path(datarepo_dir.eqsans) / "sensitivities" / "default_sensitivity.nxs"),
            "fluxMonitorRatioFile": None,
            "cutTOFmin": 0.0,
            "cutTOFmax": 0.0,
            "wavelengthStepType": "constant Delta lambda",
            "wavelengthStep": 0.1,
            "useDefaultMask": True,
            "sampleApertureSize": 30,
            "StandardAbsoluteScale": 1.03,
            "sampleOffset": 0,
            "normalization": "Total charge",
            "absoluteScaleMethod": "standard",
            "detectorOffset": 0,
            "mmRadiusForTransmission": 25,
            "numQxQyBins": 80,
            "1DQbinType": "scalar",
            "QbinType": "linear",
            "useErrorWeighting": False,
            "numQBins": 120,
            "AnnularAngleBin": 5,
        },
    }
    wavelength = 2.5  # wavelength of the incident beam neutrons, in Angstroms
    metadata = {
        "wavelength": wavelength,
        "wavelength_min": wavelength - 1.0,
        "wavelength_max": wavelength + 1.0,
        "wavelength_distribution": lambda events_count: np.repeat(wavelength, events_count),  # monochromatic
        "Z": 5.0,  # distance along the nominal Z-axis from the origin of coords. to the main detector, in meters
        "pulse_period": 1000.0 / 60.0,  # in mili-seconds
        "start_time": "2023-08-01 00:00:00",
        "end_time": "2023-08-01 00:00:01",  # 60 pulses in total
        "pulse_charge": 17529040.0,  # in pico-Coulombs
        "duration": 1.0,  # in seconds,
        "CenterX": 0.025239,
        "CenterY": 0.0170801,
        "vBeamSlit": 5,
        "vBeamSlit2": 5,
        "vBeamSlit3": 5,
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
    metadata["Q_at_max_I"] = 2 * (2 * np.pi / metadata["wavelength"]) * np.sin(theta_at_max)

    if missing_files(
        datadir,
        run_numbers=["92300", "92310", "92320", "92330"],
        file_names=["beamstop60_mask_4m.nxs", "empty_dark_file.nxs"],
    ):
        create_three_rings_pattern(config, metadata)
    return dict(config=config, metadata=metadata)


def _mock_LoadEventNexus(*args, **kwargs):
    # Substitute LoadEventNexus with LoadNexusProcessed because our synthetic files were created with SaveNexus
    return LoadNexusProcessed(Filename=kwargs["Filename"], OutputWorkspace=kwargs["OutputWorkspace"])


def _mock_eqsans_load_events(*args, **kwargs):
    output_workspace = generic_load_events(*args, **kwargs)
    SampleLogs(output_workspace).insert("is_frame_skipping", 0)
    return output_workspace


@pytest.mark.datarepo
@mock_patch("drtsans.tof.eqsans.api.load_events", new=_mock_eqsans_load_events)
@mock_patch("drtsans.tof.eqsans.load.load_events", new=_mock_eqsans_load_events)
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
def test_reduce_ring(ring_pattern: dict, temp_directory: Callable[[Any], str]):
    r"""
    Reduce a single configuration that simulates scattering from a sample which imprints an intensity pattern in the
    shape of a ring in the main detector of EQSANS.

    Parameters
    ----------
    ring_pattern
        Fixture evaluating to a reduction configuration, containing all the necessary run numbers
        and reduction parameters
    temp_directory
        Fixture evaluating to a function. When invoked, returns the path to a directory which will be erased upon
        completion of the test, irrespective of the test being successful or not.
    """
    # pad missing parameters with default values from schema EQSANS.json
    config = reduction_parameters(ring_pattern.config, "EQSANS", validate=False)

    # insert parameters customized for the current test
    sample_run_number = config["sample"]["runNumber"]
    amendments = {
        "outputFileName": f"EQSANS_{sample_run_number}",
        "configuration": {
            "outputDir": temp_directory(prefix="testReduceRing_"),
        },
    }
    config = update_reduction_parameters(config, amendments, validate=True)
    config["configuration"]["darkFileName"] = None  # force no dark field correction
    metadata = ring_pattern.metadata

    def _mock_transform_to_wavelength(*args, **kwargs):
        r"""This mock requires information from fixture ring_pattern, hence we define it inside the body of the
        test"""
        input_workspace = args[0]
        output_workspace = kwargs.get("output_workspace", str(input_workspace))
        ConvertUnits(
            InputWorkspace=input_workspace,
            Target="Wavelength",
            Emode="Elastic",
            OutputWorkspace=output_workspace,
        )
        w_min, w_max = metadata["wavelength_min"], metadata["wavelength_max"]
        Rebin(
            InputWorkspace=output_workspace,
            Params=(w_min, kwargs["bin_width"], w_max),
            PreserveEvents=False,
            OutputWorkspace=output_workspace,
        )
        # log wavelength relevant info
        sample_logs = SampleLogs(output_workspace)
        sample_logs.insert("wavelength_min", w_min, unit="Angstrom")
        sample_logs.insert("wavelength_max", w_max, unit="Angstrom")
        sample_logs.insert("wavelength_lead_min", w_min, unit="Angstrom")
        sample_logs.insert("wavelength_lead_max", w_max, unit="Angstrom")
        # log time-of-flight relevant info (needed for dark-current correction)
        low_tof_clip, high_tof_clip = kwargs["low_tof_clip"], kwargs["high_tof_clip"]
        pulse_period = metadata["pulse_period"]
        sample_logs.insert("low_tof_clip", low_tof_clip, unit="ms"),
        sample_logs.insert("high_tof_clip", high_tof_clip, unit="ms"),
        sample_logs.insert("tof_frame_width", pulse_period, unit="ms")
        tof_width_clipped = pulse_period - low_tof_clip - high_tof_clip
        sample_logs.insert("tof_frame_width_clipped", tof_width_clipped, unit="ms")
        return mtd[output_workspace], None

    # load all necessary files
    with mock_patch("drtsans.tof.eqsans.load.transform_to_wavelength", side_effect=_mock_transform_to_wavelength):
        with amend_config(data_dir=config["dataDirectories"]):
            loaded = load_all_files(config)

    # do the actual reduction
    reduction_output = reduce_single_configuration(loaded, config)

    # find the Q-modulus where the scattering intensity is maximum, and compare to what's expected
    i_vs_qmod: IQmod = reduction_output[0].I1D_main[0]
    q_at_max_i = i_vs_qmod.mod_q[i_vs_qmod.intensity.argmax()]  # Q with maximum scattering intensity
    assert_allclose(q_at_max_i, metadata["Q_at_max_I"], atol=0.002)


@pytest.mark.datarepo
@mock_patch("drtsans.tof.eqsans.api.load_events", new=_mock_eqsans_load_events)
@mock_patch("drtsans.tof.eqsans.load.load_events", new=_mock_eqsans_load_events)
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
def test_reduce_three_rings(three_rings_pattern: dict, temp_directory: Callable[[Any], str]):
    r"""
    Reduce a single configuration that simulates scattering from a sample which imprints an intensity pattern in the
    shape of a ring in the main detector of EQSANS.

    Parameters
    ----------
    ring_pattern
        Fixture evaluating to a reduction configuration, containing all the necessary run numbers
        and reduction parameters
    temp_directory
        Fixture evaluating to a function. When invoked, returns the path to a directory which will be erased upon
        completion of the test, irrespective of the test being successful or not.
    """
    # pad missing parameters with default values from schema EQSANS.json
    config = reduction_parameters(three_rings_pattern.config, "EQSANS", validate=False)

    # insert parameters customized for the current test
    sample_run_number = config["sample"]["runNumber"]
    amendments = {
        "outputFileName": f"EQSANS_{sample_run_number}",
        "configuration": {
            "outputDir": temp_directory(prefix="testReduceThreeRings_"),
        },
    }
    config = update_reduction_parameters(config, amendments, validate=True)
    config["configuration"]["darkFileName"] = None  # force no dark field correction
    metadata = three_rings_pattern.metadata

    def _mock_transform_to_wavelength(*args, **kwargs):
        r"""This mock requires information from fixture ring_pattern, hence we define it inside the body of the
        test"""
        input_workspace = args[0]
        output_workspace = kwargs.get("output_workspace", str(input_workspace))
        ConvertUnits(
            InputWorkspace=input_workspace,
            Target="Wavelength",
            Emode="Elastic",
            OutputWorkspace=output_workspace,
        )
        w_min, w_max = metadata["wavelength_min"], metadata["wavelength_max"]
        Rebin(
            InputWorkspace=output_workspace,
            Params=(w_min, kwargs["bin_width"], w_max),
            PreserveEvents=False,
            OutputWorkspace=output_workspace,
        )

        # add small random intensities to each histogram bin to avoid having bins with no intensities,
        # necessary when the neutron flux is monochromatic
        workspace = mtd[output_workspace]
        for i in range(workspace.getNumberHistograms()):
            intensities = workspace.readY(i)
            int_rand_max = 0.01 * max(intensities)
            int_rand = np.random.uniform(low=0.0, high=int_rand_max, size=len(intensities))
            workspace.dataY(i)[:] += int_rand
            workspace.dataE(i)[:] += np.sqrt(int_rand)

        # log wavelength relevant info
        sample_logs = SampleLogs(output_workspace)
        sample_logs.insert("wavelength_min", w_min, unit="Angstrom")
        sample_logs.insert("wavelength_max", w_max, unit="Angstrom")
        sample_logs.insert("wavelength_lead_min", w_min, unit="Angstrom")
        sample_logs.insert("wavelength_lead_max", w_max, unit="Angstrom")
        # log time-of-flight relevant info (needed for dark-current correction)
        low_tof_clip, high_tof_clip = kwargs["low_tof_clip"], kwargs["high_tof_clip"]
        pulse_period = metadata["pulse_period"]
        sample_logs.insert("low_tof_clip", low_tof_clip, unit="ms"),
        sample_logs.insert("high_tof_clip", high_tof_clip, unit="ms"),
        sample_logs.insert("tof_frame_width", pulse_period, unit="ms")
        tof_width_clipped = pulse_period - low_tof_clip - high_tof_clip
        sample_logs.insert("tof_frame_width_clipped", tof_width_clipped, unit="ms")
        return mtd[output_workspace], None

    # load all necessary files
    with mock_patch("drtsans.tof.eqsans.load.transform_to_wavelength", side_effect=_mock_transform_to_wavelength):
        with amend_config(data_dir=config["dataDirectories"]):
            loaded = load_all_files(config)

    # do the actual reduction
    reduction_output = reduce_single_configuration(loaded, config)

    # find the Q-modulus where the scattering intensity is maximum, and compare to what's expected
    i_vs_qmod: IQmod = reduction_output[0].I1D_main[0]
    q_values, minimum_peak_intensity = i_vs_qmod.mod_q, 100.0
    for q_at_max_i in metadata["Q_at_max_I"]:
        closest_index = np.argmin(np.abs(q_values - q_at_max_i))
        assert i_vs_qmod.intensity[closest_index] > minimum_peak_intensity


def _mock_eqsans_load_and_split(*args, **kwargs):
    workspace_group = generic_load_and_split(*args, **kwargs)
    for workspace in workspace_group:
        SampleLogs(workspace).insert("is_frame_skipping", 0)
    return workspace_group


@pytest.mark.datarepo
@mock_patch("drtsans.tof.eqsans.load.load_and_split", new=_mock_eqsans_load_and_split)
@mock_patch("drtsans.tof.eqsans.api.load_events", new=_mock_eqsans_load_events)
@mock_patch("drtsans.tof.eqsans.load.load_events", new=_mock_eqsans_load_events)
@mock_patch("drtsans.load.LoadEventNexus", new=_mock_LoadEventNexus)
def test_split_three_rings(three_rings_pattern: dict, temp_directory: Callable[[Any], str]):
    r"""
    Reduce a single configuration that simulates scattering from a sample which imprints an intensity pattern in the
    shape of a ring in the main detector of EQSANS.

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
    config = reduction_parameters(three_rings_pattern.config, "EQSANS", validate=False)

    # insert parameters customized for the current test
    sample_run_number = config["sample"]["runNumber"]
    amendments = {
        "outputFileName": f"EQSANS_{sample_run_number}",
        "configuration": {
            "outputDir": temp_directory(prefix="testSplitThreeRings_"),
            "useTimeSlice": True,
            "timeSliceInterval": 1.0 / 60,
            "timeSliceOffset": 0.0,
            "timeSlicePeriod": 3.0 / 60,
        },
    }
    config = update_reduction_parameters(config, amendments, validate=True)
    config["configuration"]["darkFileName"] = None  # force no dark field correction
    metadata = three_rings_pattern.metadata

    def _mock_transform_to_wavelength(*args, **kwargs):
        r"""This mock requires information from fixture ring_pattern, hence we define it inside the body of the
        test"""
        input_workspace = args[0]
        output_workspace = kwargs.get("output_workspace", str(input_workspace))
        ConvertUnits(
            InputWorkspace=input_workspace,
            Target="Wavelength",
            Emode="Elastic",
            OutputWorkspace=output_workspace,
        )
        w_min, w_max = metadata["wavelength_min"], metadata["wavelength_max"]
        Rebin(
            InputWorkspace=output_workspace,
            Params=(w_min, kwargs["bin_width"], w_max),
            PreserveEvents=False,
            OutputWorkspace=output_workspace,
        )

        # add small random intensities to each histogram bin to avoid having bins with no intensities,
        # necessary when the neutron flux is monochromatic
        workspace = mtd[output_workspace]
        for i in range(workspace.getNumberHistograms()):
            intensities = workspace.readY(i)
            int_rand_max = 0.01 * max(intensities)
            int_rand = np.random.uniform(low=0.0, high=int_rand_max, size=len(intensities))
            workspace.dataY(i)[:] += int_rand
            workspace.dataE(i)[:] += np.sqrt(int_rand)

        # log wavelength relevant info
        sample_logs = SampleLogs(output_workspace)
        sample_logs.insert("wavelength_min", w_min, unit="Angstrom")
        sample_logs.insert("wavelength_max", w_max, unit="Angstrom")
        sample_logs.insert("wavelength_lead_min", w_min, unit="Angstrom")
        sample_logs.insert("wavelength_lead_max", w_max, unit="Angstrom")
        # log time-of-flight relevant info (needed for dark-current correction)
        low_tof_clip, high_tof_clip = kwargs["low_tof_clip"], kwargs["high_tof_clip"]
        pulse_period = metadata["pulse_period"]
        sample_logs.insert("low_tof_clip", low_tof_clip, unit="ms"),
        sample_logs.insert("high_tof_clip", high_tof_clip, unit="ms"),
        sample_logs.insert("tof_frame_width", pulse_period, unit="ms")
        tof_width_clipped = pulse_period - low_tof_clip - high_tof_clip
        sample_logs.insert("tof_frame_width_clipped", tof_width_clipped, unit="ms")
        return mtd[output_workspace], None

    # load all necessary files
    with mock_patch("drtsans.tof.eqsans.load.transform_to_wavelength", side_effect=_mock_transform_to_wavelength):
        with amend_config(data_dir=config["dataDirectories"]):
            loaded = load_all_files(config)

    # do the actual reduction
    reduction_output = reduce_single_configuration(loaded, config)

    # find the Q-modulus where the scattering intensity is maximum, and compare to what's expected
    minimum_peak_intensity = 300
    for peak_index, q_at_max_i in enumerate(metadata["Q_at_max_I"]):
        i_vs_qmod: IQmod = reduction_output[peak_index].I1D_main[0]  # 1D intensity profile
        closest_index = np.argmin(np.abs(i_vs_qmod.mod_q - q_at_max_i))
        assert i_vs_qmod.intensity[closest_index] > minimum_peak_intensity


if __name__ == "__main__":
    pytest.main([__file__])
