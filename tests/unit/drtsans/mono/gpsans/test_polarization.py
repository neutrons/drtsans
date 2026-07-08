import h5py
import numpy as np
from mantid.kernel import amend_config
from mantid.simpleapi import CreateSingleValuedWorkspace, CreateWorkspace, DeleteWorkspace, mtd
from numpy.testing import assert_equal, assert_array_almost_equal
import pytest

from drtsans.instruments import empty_instrument_workspace
from drtsans.polarization import (
    FullPolarizationDecoder,
    HalfPolarizationDecoder,
    PolarizationDecoder,
    PolarizationLevel,
    PolarizationCrossSection,
    PolarizationState,
    polarized_sample,
    _calc_flipping_ratio,
    half_polarization,
    SimulatedPolarizationLogs,
    TimesGeneratorSpecs,
    PV_ANALYZER,
    PV_ANALYZER_FLIPPER,
    PV_ANALYZER_VETO,
    PV_POLARIZER,
    PV_POLARIZER_FLIPPER,
    PV_POLARIZER_VETO,
)
from drtsans.samplelogs import SampleLogs


class TestPolarizationLevel:
    def test_from_int(self):
        assert PolarizationLevel.from_int(0) == PolarizationLevel.NONE
        assert PolarizationLevel.from_int(0) == "none"  # a StrEnum instance also compares equal to its string value
        assert PolarizationLevel.from_int(1) == PolarizationLevel.HALF
        assert PolarizationLevel.from_int(1) == "half"
        assert PolarizationLevel.from_int(2) == PolarizationLevel.FULL
        assert PolarizationLevel.from_int(2) == "full"
        with pytest.raises(ValueError) as excinfo:
            PolarizationLevel.from_int(3)
        assert "Invalid polarization mode integer: 3. Must be 0, 1, or 2." in str(excinfo.value)

    def test_get(self, tmp_path):
        # case "event Nexus file" with no polarization metadata
        nexus_file = tmp_path / "CG2_12345.nxs.h5"
        with h5py.File(nexus_file, "w") as write_handle:
            write_handle.require_group("/entry/DASlogs")
        assert PolarizationLevel.get(str(nexus_file)) == PolarizationLevel.NONE

        # case "event Nexus file" with half polarization metadata
        with h5py.File(nexus_file, "w") as write_handle:
            group = write_handle.require_group(f"/entry/DASlogs/{PV_POLARIZER}")
            group.create_dataset("value", data=0)
        assert PolarizationLevel.get(str(nexus_file)) == PolarizationLevel.NONE
        with h5py.File(nexus_file, "r+") as write_handle:
            group = write_handle[f"/entry/DASlogs/{PV_POLARIZER}"]
            group["value"][...] = 1
        assert PolarizationLevel.get(str(nexus_file)) == PolarizationLevel.HALF

        # case "event Nexus file" with full polarization metadata
        with h5py.File(nexus_file, "r+") as write_handle:
            group = write_handle.require_group(f"/entry/DASlogs/{PV_ANALYZER}")
            group.create_dataset("value", data=0)
        assert PolarizationLevel.get(str(nexus_file)) == PolarizationLevel.HALF
        with h5py.File(nexus_file, "r+") as write_handle:
            group = write_handle[f"/entry/DASlogs/{PV_ANALYZER}"]
            group["value"][...] = 1
        assert PolarizationLevel.get(str(nexus_file)) == PolarizationLevel.FULL

        # case workspace that is not an EventWorkspace
        ws = CreateSingleValuedWorkspace(DataValue=0.95, ErrorValue=0.01, OutputWorkspace=mtd.unique_hidden_name())
        with pytest.raises(TypeError) as excinfo:
            PolarizationLevel.get(ws)  # passing workspace object
        assert f"{str(ws)} must be either a string or an EventWorkspace object" in str(excinfo.value)
        with pytest.raises(TypeError) as excinfo:
            PolarizationLevel.get(str(ws))  # passing workspace name
        assert f"The workspace '{ws}' is not an EventWorkspace" in str(excinfo.value)

        # case EventWorkspace with no polarization metadata
        ws = empty_instrument_workspace(str(ws), instrument_name="GPSANS", event_workspace=True)
        assert PolarizationLevel.get(ws) == PolarizationLevel.NONE
        assert PolarizationLevel.get(str(ws)) == PolarizationLevel.NONE

        # case EventWorkspace with half polarization metadata
        SampleLogs(ws).insert(PV_POLARIZER, 0)
        assert PolarizationLevel.get(ws) == PolarizationLevel.NONE
        assert PolarizationLevel.get(str(ws)) == PolarizationLevel.NONE
        SampleLogs(ws).insert(PV_POLARIZER, 1)
        assert PolarizationLevel.get(ws) == PolarizationLevel.HALF
        assert PolarizationLevel.get(str(ws)) == PolarizationLevel.HALF

        # case EventWorkspace with full polarization metadata
        SampleLogs(ws).insert(PV_ANALYZER, 0)
        assert PolarizationLevel.get(ws) == PolarizationLevel.HALF
        assert PolarizationLevel.get(str(ws)) == PolarizationLevel.HALF
        SampleLogs(ws).insert(PV_ANALYZER, 1)
        assert PolarizationLevel.get(ws) == PolarizationLevel.FULL
        assert PolarizationLevel.get(str(ws)) == PolarizationLevel.FULL


class TestPolarizationCrossSection:
    def test_values(self):
        values = sorted([member.value for member in PolarizationCrossSection])
        assert values == ["none", "off", "off_off", "off_on", "on", "on_off", "on_on"]

    def test_log_get(self):
        ws = CreateSingleValuedWorkspace(DataValue=0.95, ErrorValue=0.01, OutputWorkspace=mtd.unique_hidden_name())
        for cross_section in PolarizationCrossSection:
            cross_section.log(ws)
            assert PolarizationCrossSection.get(ws) == str(cross_section)  # StrEnum can compare to its string value

    def test_level(self):
        assert PolarizationCrossSection.NONE.level == "none"
        for cross_section in ["on", "off"]:
            assert PolarizationCrossSection(cross_section).level == "half"
        for cross_section in ["off_off", "off_on", "on_off", "on_on"]:
            assert PolarizationCrossSection(cross_section).level == "full"


class TestPolarizationState:
    def test_values(self):
        values = sorted([member.value for member in PolarizationState])
        assert values == ["down", "down_down", "down_up", "none", "up", "up_down", "up_up"]

    def test_log_get(self):
        ws = CreateSingleValuedWorkspace(DataValue=0.95, ErrorValue=0.01, OutputWorkspace=mtd.unique_hidden_name())
        for cross_section in PolarizationState:
            cross_section.log(ws)
            assert PolarizationState.get(ws) == str(cross_section)  # StrEnum can compare to its string value

    def test_level(self):
        assert PolarizationState.NONE.level == "none"
        for cross_section in ["up", "down"]:
            assert PolarizationState(cross_section).level == "half"
        for cross_section in ["down_down", "down_up", "up_down", "up_up"]:
            assert PolarizationState(cross_section).level == "full"


def test_polarized_sample(tmp_path):
    # Test when polarization level is already set to NONE in reduction config.
    reduction_config = {"polarization": {"level": "none"}}
    result = polarized_sample(reduction_config)
    assert result is False

    # case: polarization level already set to HALF
    reduction_config = {"polarization": {"level": "half"}}
    assert polarized_sample(reduction_config) is True

    # case: polarization level already set to FULL
    reduction_config = {"polarization": {"level": "full", "extra_key": "extra_value"}}
    assert polarized_sample(reduction_config) is True
    assert reduction_config["polarization"]["extra_key"] == "extra_value"

    with amend_config(data_dir=str(tmp_path)):
        # case: single sample run with no polarization metatadata
        nexus_file = tmp_path / "CG2_12345.nxs.h5"
        with h5py.File(nexus_file, "w") as write_handle:
            write_handle.require_group("/entry/DASlogs")
        reduction_input = {
            "configuration": {},
            "sample": {"runNumber": "12345"},
            "instrumentName": "CG2",
            "iptsNumber": "1234",
        }
        assert polarized_sample(reduction_input) is False
        assert reduction_input["configuration"]["polarization"]["level"] == "none"

        # case: single sample run with half polarization
        with h5py.File(nexus_file, "a") as write_handle:
            group = write_handle.require_group(f"/entry/DASlogs/{PV_POLARIZER}")
            group.create_dataset("value", data=1)
        reduction_input["configuration"]["polarization"] = {"extra_key": "extra_value"}  # clear "level""
        assert polarized_sample(reduction_input) is True
        assert reduction_input["configuration"]["polarization"]["level"] == "half"
        assert "extra_key" not in reduction_input["configuration"]["polarization"]  # deleted obsolete key

        # case: single sample run with full polarization
        with h5py.File(nexus_file, "a") as write_handle:
            group = write_handle.require_group(f"/entry/DASlogs/{PV_ANALYZER}")
            group.create_dataset("value", data=1)
        reduction_input["configuration"]["polarization"] = {}  # clear the polarization entry
        assert polarized_sample(reduction_input) is True
        assert reduction_input["configuration"]["polarization"]["level"] == "full"

        # case: multiple sample runs when all are unpolarized
        for run_num in ["12349", "12350", "12351"]:
            nexus_file = tmp_path / f"CG2_{run_num}.nxs.h5"
            with h5py.File(nexus_file, "w") as write_handle:
                write_handle.require_group("/entry/DASlogs")
        reduction_input = {
            "sample": {"runNumber": "12349, 12350, 12351"},
            "instrumentName": "CG2",
            "iptsNumber": "1234",
            "configuration": {"polarization": {}},
        }
        assert polarized_sample(reduction_input) is False
        assert reduction_input["configuration"]["polarization"]["level"] == "none"

        # case: multiple sample runs with polarization raises ValueError
        nexus_file = tmp_path / "CG2_12349.nxs.h5"
        with h5py.File(nexus_file, "a") as write_handle:
            group = write_handle.require_group(f"/entry/DASlogs/{PV_POLARIZER}")
            group.create_dataset("value", data=1)
        reduction_input["configuration"]["polarization"] = {}  # clear the polarization entry
        with pytest.raises(ValueError) as excinfo:
            polarized_sample(reduction_input)
        assert "Can't do polarization reduction on summed data sets" in str(excinfo.value)


def test_flipping_ratio(temp_workspace_name, clean_workspace):
    """Test for the calculation of the flipping ratio in section 9.1. This was used to determine
    that the uncertainty needs to be calculated separately because of accumulation of numeric
    error.

    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Lisa DeBeer-Schmitt <debeerschmlm@ornl.gov>
          Mike Fitzsimmons <fitzsimmonsm@ornl.gov>
    """
    # this is called "P" in the document
    polarization = CreateSingleValuedWorkspace(DataValue=0.95, ErrorValue=0.01, OutputWorkspace=temp_workspace_name())
    flipping_ratio_expected = (1.0 + 0.95) / (1.0 - 0.95)
    flipping_ratio_err_expected = (2.0 * polarization.readE(0)[0]) / 0.0025  # denominator is (1-p)^2

    flipping_ratio = _calc_flipping_ratio(polarization)
    clean_workspace(flipping_ratio)

    assert flipping_ratio.extractY() == flipping_ratio_expected
    assert flipping_ratio.extractE() == pytest.approx(flipping_ratio_err_expected)


def test_half_polarization(temp_workspace_name):
    """Test the calculation and application of the half polarization from section 9.1
    and requires reading section 9.0 for definition variables.

    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Lisa DeBeer-Schmitt <debeerschmlm@ornl.gov>
          Mike Fitzsimmons <fitzsimmonsm@ornl.gov>
    """
    # this is called "P" in the document
    polarization = CreateSingleValuedWorkspace(DataValue=0.95, ErrorValue=0.01, OutputWorkspace=temp_workspace_name())
    # this is called "e" in the document
    efficiency = CreateSingleValuedWorkspace(DataValue=0.998, ErrorValue=0.001, OutputWorkspace=temp_workspace_name())

    # values for the measured flipper off (M0) and flipper on (M1)
    M0 = CreateSingleValuedWorkspace(DataValue=10000, ErrorValue=100, OutputWorkspace=temp_workspace_name())
    M1 = CreateSingleValuedWorkspace(DataValue=8100, ErrorValue=90, OutputWorkspace=temp_workspace_name())

    # expected results
    SpinUpExp = CreateSingleValuedWorkspace(
        DataValue=10050.100,
        ErrorValue=103.2046,
        OutputWorkspace=temp_workspace_name(),  # was 103.205
    )
    SpinDownExp = CreateSingleValuedWorkspace(
        DataValue=8046.0925, ErrorValue=93.2163, OutputWorkspace=temp_workspace_name()
    )

    # do the calculation
    SpinUp, SpinDown = half_polarization(M0, M1, polarization, efficiency)

    # compare with what was expected
    assert SpinUp.extractY()[0][0] == pytest.approx(SpinUpExp.extractY()[0][0])
    assert SpinUp.extractE()[0][0] == pytest.approx(SpinUpExp.extractE()[0][0])
    assert SpinDown.extractY()[0][0] == pytest.approx(SpinDownExp.extractY()[0][0])
    assert SpinDown.extractE()[0][0] == pytest.approx(SpinDownExp.extractE()[0][0])


class TestSimulatedLogs:
    def test_valid_flipper_generators(self):
        log = SimulatedPolarizationLogs(
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {}),
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {}),
        )
        assert log.polarizer_flipper.name == "heartbeat"
        assert log.analyzer_flipper.name == "heartbeat"

    def test_invalid_flipper_generators(self):
        with pytest.raises(ValueError) as excinfo:
            SimulatedPolarizationLogs(polarizer_flipper=TimesGeneratorSpecs("invalid_generator", {}))
        assert (
            "The polarizer flipper generator must be one of ['heartbeat', 'binary_pulse', 'cycled_intervals']"
            in str(excinfo.value)
        )
        with pytest.raises(ValueError) as excinfo:
            SimulatedPolarizationLogs(analyzer_flipper=TimesGeneratorSpecs("invalid_generator", {}))
        assert (
            "The analyzer flipper generator must be one of ['heartbeat', 'binary_pulse', 'cycled_intervals']"
            in str(excinfo.value)
        )

    def test_valid_veto_generators(self):
        log = SimulatedPolarizationLogs(
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {}),
        )
        assert log.polarizer_veto.name == "binary_pulse"
        assert log.analyzer_veto.name == "binary_pulse"

    def test_invalid_veto_generators(self):
        with pytest.raises(ValueError) as excinfo:
            SimulatedPolarizationLogs(polarizer_veto=TimesGeneratorSpecs("invalid_veto", {}))
        assert "polarizer veto generator must be one of ['binary_pulse']" in str(excinfo.value)
        with pytest.raises(ValueError) as excinfo:
            SimulatedPolarizationLogs(analyzer_veto=TimesGeneratorSpecs("invalid_veto", {}))
        assert "analyzer veto generator must be one of ['binary_pulse']" in str(excinfo.value)

    def test_heartbeat_generator(self):
        times = SimulatedPolarizationLogs().heartbeat(interval=1.0, upper_bound=10)
        assert_array_almost_equal(list(times), [0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        times = SimulatedPolarizationLogs().heartbeat(interval=1.0, dead_time=3.5, upper_bound=10)
        assert_array_almost_equal(list(times), [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])

    def test_binary_pulse_generator(self):
        times = SimulatedPolarizationLogs().binary_pulse(interval=3.0, alive_duration=1.0, upper_bound=10)
        assert_array_almost_equal(list(times), [0, 2.5, 3.5, 5.5, 6.5, 8.5, 9.5], decimal=2)
        times = SimulatedPolarizationLogs().binary_pulse(
            interval=3.0, alive_duration=1.0, dead_time=2.7, upper_bound=10
        )
        assert_array_almost_equal(list(times), [3.5, 5.5, 6.5, 8.5, 9.5], decimal=2)

    def test_cycled_intervals_generator(self):
        times = SimulatedPolarizationLogs().cycled_intervals(
            intervals=[2.0 / 60, 1.0 / 60], dead_time=0, upper_bound=0.2
        )
        assert list(times)[-1] == 0.2  # makes sure rounding errors don't accumulate
        times = SimulatedPolarizationLogs().cycled_intervals(intervals=[1.0, 2.0], upper_bound=9.0)
        assert_array_almost_equal(list(times), [0.0, 1.0, 3.0, 4.0, 6.0, 7.0, 9.0])
        times = SimulatedPolarizationLogs().cycled_intervals(intervals=[1.0, 2.0], dead_time=3.5, upper_bound=9.0)
        assert_array_almost_equal(list(times), [4.0, 6.0, 7.0, 9.0])

    def test_times_generator(self):
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 1.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 1.0, "alive_duration": 0.2}),
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 2.0}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 2.0, "alive_duration": 0.4}),
        )
        times = logs.times_generator(PV_POLARIZER_FLIPPER, upper_bound=6.0)
        assert_array_almost_equal(list(times), [0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], decimal=2)
        times = logs.times_generator(PV_POLARIZER_VETO, upper_bound=6.0)
        assert_array_almost_equal(list(times), [0.0, 0.9, 1.1, 1.9, 2.1, 2.9, 3.1, 3.9, 4.1, 4.9, 5.1, 5.9], decimal=2)
        times = logs.times_generator(PV_ANALYZER_FLIPPER, upper_bound=6.0)
        assert_array_almost_equal(list(times), [0, 2.0, 4.0, 6.0], decimal=2)
        times = logs.times_generator(PV_ANALYZER_VETO, upper_bound=6.0)
        assert_array_almost_equal(list(times), [0.0, 1.8, 2.2, 3.8, 4.2, 5.8], decimal=2)

    def test_inject(self):
        # create a workspace with required sample logs
        workspace = CreateSingleValuedWorkspace(OutputWorkspace=mtd.unique_hidden_name())
        sample_logs = SampleLogs(workspace)
        sample_logs.insert("start_time", "2023-10-01T00:00:00")
        sample_logs.insert("duration", 300)  # 5 minutes
        # inject the simulated logs. Notice there are no longs for the analyzer veto
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 60.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "alive_duration": 1.0}),
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 120}),
            analyzer_veto=None,
        )
        logs.inject(workspace)
        # check polarizer and analyzer values
        assert sample_logs[PV_POLARIZER].value == 1
        assert sample_logs[PV_ANALYZER].value == 2
        # check time-series values
        assert_equal(sample_logs[PV_POLARIZER_FLIPPER].value, [0, 1, 0, 1, 0, 1])
        assert_equal(sample_logs[PV_POLARIZER_VETO].value, [0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
        assert_equal(sample_logs[PV_ANALYZER_FLIPPER].value, [0, 1, 0])
        # check last time in the time-series times
        assert "T00:05:00.0" in str(sample_logs[PV_POLARIZER_FLIPPER].times[-1])
        assert "T00:04:59.5" in str(sample_logs[PV_POLARIZER_VETO].times[-1])
        assert "T00:04:00.0" in str(sample_logs[PV_ANALYZER_FLIPPER].times[-1])
        assert (PV_ANALYZER_VETO in sample_logs) is False


@pytest.fixture
def half_pol_workspaces():
    """Two Workspace2D objects tagged with OFF/ON cross-section logs and wavelength=6."""
    names, ws_list = [], []
    for cross_section, y_val in [
        (PolarizationCrossSection.OFF, 10.0),
        (PolarizationCrossSection.ON, 8.0),
    ]:
        name = mtd.unique_hidden_name()
        ws = CreateWorkspace(DataX=[5.0, 7.0], DataY=[y_val], DataE=[0.1], OutputWorkspace=name)
        SampleLogs(ws).insert("wavelength", 6.0)
        cross_section.log(ws)
        names.append(name)
        ws_list.append(ws)
    yield ws_list
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def full_pol_workspaces():
    """Four Workspace2D objects tagged with full-polarization cross-section logs and wavelength=6."""
    config = {
        "polarization": {
            "polarizer": {"polarization": str(1.0 / 3.0), "efficiency": "1.0"},
            "analyzer": {"polarizationZero": "0.5", "polarizationPi": "0.25"},
        }
    }
    spin_values = np.array([13.0, 17.0, 19.0, 23.0])
    device_values = FullPolarizationDecoder(config).encoding_matrix(wavelength=6.0) @ spin_values

    names, ws_list = [], []
    for cross_section, y_val in zip(
        [
            PolarizationCrossSection.OFF_OFF,
            PolarizationCrossSection.ON_OFF,
            PolarizationCrossSection.OFF_ON,
            PolarizationCrossSection.ON_ON,
        ],
        device_values,
    ):
        name = mtd.unique_hidden_name()
        ws = CreateWorkspace(DataX=[5.0, 7.0], DataY=[y_val], DataE=[0.1], OutputWorkspace=name)
        SampleLogs(ws).insert("wavelength", 6.0)
        cross_section.log(ws)
        names.append(name)
        ws_list.append(ws)
    yield ws_list, spin_values, config
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


class TestPolarizationDecoder:
    def test_constant_polarization_and_efficiency(self):
        config = {"polarization": {"polarizer": {"polarization": "0.9", "efficiency": "0.8"}}}
        decoder = PolarizationDecoder(config)
        assert decoder.p(6.0) == pytest.approx(0.9)
        assert decoder.e(6.0) == pytest.approx(0.8)

    def test_linear_polarization(self):
        config = {"polarization": {"polarizer": {"polarization": "0.95 - 0.01*(x - 16)", "efficiency": "1"}}}
        decoder = PolarizationDecoder(config)
        assert decoder.p(16.0) == pytest.approx(0.95)
        assert decoder.p(6.0) == pytest.approx(0.95 - 0.01 * (6.0 - 16.0))

    def test_numeric_value_coercion(self):
        config = {"polarization": {"polarizer": {"polarization": 0.9, "efficiency": 0.8}}}
        decoder = PolarizationDecoder(config)
        assert decoder.p(10.0) == pytest.approx(0.9)
        assert decoder.e(10.0) == pytest.approx(0.8)

    def test_defaults_are_unity(self):
        decoder = PolarizationDecoder({})
        for wavelength in [5.0, 10.0, 16.0]:
            assert decoder.p(wavelength) == pytest.approx(1.0)
            assert decoder.e(wavelength) == pytest.approx(1.0)

    def test_decode_not_implemented(self):
        with pytest.raises(NotImplementedError):
            PolarizationDecoder({}).decode([])

    @pytest.mark.parametrize("value", [-1.0, 1.0])
    def test_polarization_includes_signed_bounds(self, value):
        PolarizationDecoder._validate_polarization("Polarization", value)

    @pytest.mark.parametrize("value", [-1.1, 0.0, 1.1])
    def test_polarization_rejects_zero_and_outside_signed_bounds(self, value):
        with pytest.raises(ValueError, match=r"Polarization must be non-zero and in the interval \[-1, 1\]"):
            PolarizationDecoder._validate_polarization("Polarization", value)

    def test_efficiency_includes_upper_bound(self):
        PolarizationDecoder._validate_efficiency("Flipper efficiency", 1.0)

    @pytest.mark.parametrize("value", [-0.1, 0.0, 1.1])
    def test_efficiency_rejects_zero_and_outside_unit_bounds(self, value):
        with pytest.raises(ValueError, match=r"Flipper efficiency must be in the interval \(0, 1\]"):
            PolarizationDecoder._validate_efficiency("Flipper efficiency", value)


class TestHalfPolarizationDecoder:
    def _make_decoder(self, polarization, efficiency):
        config = {"polarization": {"polarizer": {"polarization": str(polarization), "efficiency": str(efficiency)}}}
        return HalfPolarizationDecoder(config)

    # --- Group A: decoding_matrix (pure numpy) ---

    def test_identity_at_perfect_polarization(self):
        """P=1, e=1 → perfect instrument → decoding matrix is the identity."""
        decoder = self._make_decoder(polarization=1.0, efficiency=1.0)
        M = decoder.decoding_matrix(wavelength=6.0)
        np.testing.assert_array_almost_equal(M, np.eye(2))

    def test_matrix_shape(self):
        decoder = self._make_decoder(polarization=0.9, efficiency=0.95)
        assert decoder.decoding_matrix(wavelength=6.0).shape == (2, 2)

    def test_matrix_values_known_case(self):
        """P=1/3, e=1 → dl=1, ul=2 → M = [[2,-1],[-1,2]]."""
        decoder = self._make_decoder(polarization=1.0 / 3.0, efficiency=1.0)
        M = decoder.decoding_matrix(wavelength=6.0)
        np.testing.assert_array_almost_equal(M, [[2, -1], [-1, 2]])

    def test_polarization_below_negative_one_raises(self):
        decoder = self._make_decoder(polarization=-1.1, efficiency=1.0)
        with pytest.raises(ValueError, match="Polarization must be non-zero and in the interval"):
            decoder.decoding_matrix(wavelength=6.0)

    def test_efficiency_below_zero_raises(self):
        decoder = self._make_decoder(polarization=0.9, efficiency=-0.1)
        with pytest.raises(ValueError, match=r"Flipper efficiency must be in the interval \(0, 1\]"):
            decoder.decoding_matrix(wavelength=6.0)

    def test_efficiency_zero_raises(self):
        decoder = self._make_decoder(polarization=0.9, efficiency=0.0)
        with pytest.raises(ValueError, match=r"Flipper efficiency must be in the interval \(0, 1\]"):
            decoder.decoding_matrix(wavelength=6.0)

    def test_efficiency_above_one_raises(self):
        decoder = self._make_decoder(polarization=0.9, efficiency=1.1)
        with pytest.raises(ValueError, match=r"Flipper efficiency must be in the interval \(0, 1\]"):
            decoder.decoding_matrix(wavelength=6.0)

    def test_negative_polarization_is_in_valid_range(self):
        decoder = self._make_decoder(polarization=-0.5, efficiency=0.8)
        assert decoder.decoding_matrix(wavelength=6.0).shape == (2, 2)

    def test_matrix_inverts_encoding(self):
        """Decoding matrix M_dec is the inverse of the physical encoding matrix M_enc."""

        def encoding_matrix(p, e):
            # Row 0: flipper off — fraction of S↑ and S↓ that pass the polarizer
            # Row 1: flipper on  — flipper flips spin with efficiency e before analysis
            return np.array(
                [
                    [(1 + p) / 2, (1 - p) / 2],
                    [
                        (1 - e) * (1 + p) / 2 + e * (1 - p) / 2,
                        (1 - e) * (1 - p) / 2 + e * (1 + p) / 2,
                    ],
                ]
            )

        for p, e in [(0.9, 0.95), (0.5, 0.8), (1.0 / 3.0, 1.0), (0.95, 0.998)]:
            decoder = self._make_decoder(polarization=p, efficiency=e)
            M_dec = decoder.decoding_matrix(wavelength=6.0)
            M_enc = encoding_matrix(p, e)
            np.testing.assert_array_almost_equal(M_dec @ M_enc, np.eye(2), decimal=10)

    # --- Group B: decode (Mantid workspaces) ---

    def test_returns_two_workspaces(self, half_pol_workspaces, clean_workspace):
        decoder = self._make_decoder(polarization=0.9, efficiency=0.95)
        result = decoder.decode(half_pol_workspaces)
        for ws in result:
            clean_workspace(ws)
        assert len(result) == 2

    def test_output_logs_have_polarization_state(self, half_pol_workspaces, clean_workspace):
        decoder = self._make_decoder(polarization=0.9, efficiency=0.95)
        result = decoder.decode(half_pol_workspaces)
        for ws in result:
            clean_workspace(ws)
        assert PolarizationState.get(result[0]) == PolarizationState.UP
        assert PolarizationState.get(result[1]) == PolarizationState.DOWN

    def test_output_logs_lack_cross_section(self, half_pol_workspaces, clean_workspace):
        decoder = self._make_decoder(polarization=0.9, efficiency=0.95)
        result = decoder.decode(half_pol_workspaces)
        for ws in result:
            clean_workspace(ws)
        for ws in result:
            assert PolarizationCrossSection.logname not in SampleLogs(ws)

    def test_wrong_count_raises(self, half_pol_workspaces):
        decoder = self._make_decoder(polarization=0.9, efficiency=0.95)
        with pytest.raises(ValueError, match="exactly 2 device cross-sections"):
            decoder.decode([half_pol_workspaces[0]])
        with pytest.raises(ValueError, match="exactly 2 device cross-sections"):
            decoder.decode(half_pol_workspaces + half_pol_workspaces)

    def test_order_invariant(self, half_pol_workspaces, clean_workspace):
        """Passing [s0, s1] or [s1, s0] yields the same decoded spin states."""
        decoder = self._make_decoder(polarization=0.9, efficiency=0.95)
        result_normal = decoder.decode(half_pol_workspaces)
        result_reversed = decoder.decode(list(reversed(half_pol_workspaces)))
        for ws in result_normal + result_reversed:
            clean_workspace(ws)
        np.testing.assert_array_almost_equal(result_normal[0].readY(0), result_reversed[0].readY(0))
        np.testing.assert_array_almost_equal(result_normal[1].readY(0), result_reversed[1].readY(0))

    def test_intensity_at_perfect_polarization(self, half_pol_workspaces, clean_workspace):
        """P=1, e=1: identity decoding → S↑ = S⁰ and S↓ = S¹."""
        decoder = self._make_decoder(polarization=1.0, efficiency=1.0)
        ws_off, ws_on = half_pol_workspaces
        result = decoder.decode(half_pol_workspaces)
        for ws in result:
            clean_workspace(ws)
        np.testing.assert_array_almost_equal(result[0].readY(0), ws_off.readY(0))
        np.testing.assert_array_almost_equal(result[1].readY(0), ws_on.readY(0))

    def test_intensity_known_case(self, half_pol_workspaces, clean_workspace):
        """P=1/3, e=1: M=[[2,-1],[-1,2]] → S↑=2*10-8=12, S↓=-10+2*8=6."""
        decoder = self._make_decoder(polarization=1.0 / 3.0, efficiency=1.0)
        result = decoder.decode(half_pol_workspaces)
        for ws in result:
            clean_workspace(ws)
        np.testing.assert_array_almost_equal(result[0].readY(0), [12.0])
        np.testing.assert_array_almost_equal(result[1].readY(0), [6.0])


class TestFullPolarizationDecoder:
    def test_loads_polarizer_and_analyzer_values(self):
        config = {
            "polarization": {
                "polarizer": {"polarization": "0.9", "efficiency": "0.8"},
                "analyzer": {"polarizationZero": "0.7", "polarizationPi": "0.6"},
            }
        }
        decoder = FullPolarizationDecoder(config)

        assert decoder.p(6.0) == pytest.approx(0.9)
        assert decoder.e(6.0) == pytest.approx(0.8)
        assert decoder.p_0(6.0) == pytest.approx(0.7)
        assert decoder.p_pi(6.0) == pytest.approx(0.6)

    def test_loads_wavelength_dependent_analyzer_values(self):
        config = {
            "polarization": {
                "analyzer": {
                    "polarizationZero": "0.95 - 0.01*(x - 16)",
                    "polarizationPi": "0.90 - 0.02*(x - 16)",
                }
            }
        }
        decoder = FullPolarizationDecoder(config)

        assert decoder.p_0(16.0) == pytest.approx(0.95)
        assert decoder.p_0(6.0) == pytest.approx(0.95 - 0.01 * (6.0 - 16.0))
        assert decoder.p_pi(16.0) == pytest.approx(0.90)
        assert decoder.p_pi(6.0) == pytest.approx(0.90 - 0.02 * (6.0 - 16.0))

    def test_analyzer_defaults_select_ideal_opposite_states(self):
        decoder = FullPolarizationDecoder({})

        for wavelength in [5.0, 10.0, 16.0]:
            assert decoder.p_0(wavelength) == pytest.approx(1.0)
            assert decoder.p_pi(wavelength) == pytest.approx(-1.0)

    def test_encoding_matrix_shape(self):
        config = {
            "polarization": {
                "polarizer": {"polarization": "0.9", "efficiency": "0.95"},
                "analyzer": {"polarizationZero": "0.5", "polarizationPi": "0.25"},
            }
        }
        decoder = FullPolarizationDecoder(config)

        assert decoder.encoding_matrix(wavelength=6.0).shape == (4, 4)

    def test_encoding_matrix_values_known_case(self):
        config = {
            "polarization": {
                "polarizer": {"polarization": str(1.0 / 3.0), "efficiency": "1.0"},
                "analyzer": {"polarizationZero": "0.5", "polarizationPi": "0.25"},
            }
        }
        decoder = FullPolarizationDecoder(config)

        expected = np.array(
            [
                [1.0 / 2.0, 1.0 / 6.0, 1.0 / 4.0, 1.0 / 12.0],
                [1.0 / 4.0, 1.0 / 12.0, 1.0 / 2.0, 1.0 / 6.0],
                [5.0 / 12.0, 1.0 / 4.0, 5.0 / 24.0, 1.0 / 8.0],
                [5.0 / 24.0, 1.0 / 8.0, 5.0 / 12.0, 1.0 / 4.0],
            ]
        )

        np.testing.assert_array_almost_equal(decoder.encoding_matrix(wavelength=6.0), expected)

    def test_encoding_matrix_allows_perfect_polarizer(self):
        config = {
            "polarization": {
                "polarizer": {"polarization": "1.0", "efficiency": "1.0"},
                "analyzer": {"polarizationZero": "0.5", "polarizationPi": "0.25"},
            }
        }
        decoder = FullPolarizationDecoder(config)

        expected = np.array(
            [
                [3.0 / 4.0, 1.0 / 4.0, 0.0, 0.0],
                [0.0, 0.0, 3.0 / 4.0, 1.0 / 4.0],
                [5.0 / 8.0, 3.0 / 8.0, 0.0, 0.0],
                [0.0, 0.0, 5.0 / 8.0, 3.0 / 8.0],
            ]
        )

        np.testing.assert_array_almost_equal(decoder.encoding_matrix(wavelength=6.0), expected)

    def test_encoding_matrix_defaults_are_finite_and_select_distinct_spin_states(self):
        decoder = FullPolarizationDecoder({})

        expected = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )

        np.testing.assert_array_almost_equal(decoder.encoding_matrix(wavelength=6.0), expected)

    def test_encoding_matrix_allows_signed_polarization_values(self):
        config = {
            "polarization": {
                "polarizer": {"polarization": "-0.5", "efficiency": "0.8"},
                "analyzer": {"polarizationZero": "0.25", "polarizationPi": "-0.25"},
            }
        }
        decoder = FullPolarizationDecoder(config)

        matrix = decoder.encoding_matrix(wavelength=6.0)
        assert matrix.shape == (4, 4)
        assert np.all(np.isfinite(matrix))

    @pytest.mark.parametrize(
        "polarization_config, match",
        [
            ({"polarizer": {"polarization": "-1.1"}}, "Polarizer polarization"),
            ({"polarizer": {"polarization": "0"}}, "Polarizer polarization"),
            ({"polarizer": {"polarization": "1.1"}}, "Polarizer polarization"),
            ({"polarizer": {"efficiency": "-0.1"}}, "Flipper efficiency"),
            ({"polarizer": {"efficiency": "0"}}, "Flipper efficiency"),
            ({"polarizer": {"efficiency": "1.1"}}, "Flipper efficiency"),
            ({"analyzer": {"polarizationZero": "-1.1"}}, "Analyzer zero-state polarization"),
            ({"analyzer": {"polarizationZero": "0"}}, "Analyzer zero-state polarization"),
            ({"analyzer": {"polarizationPi": "0"}}, "Analyzer pi-state polarization"),
            ({"analyzer": {"polarizationPi": "1.1"}}, "Analyzer pi-state polarization"),
        ],
    )
    def test_encoding_matrix_rejects_out_of_range_values(self, polarization_config, match):
        config = {"polarization": polarization_config}
        decoder = FullPolarizationDecoder(config)

        with pytest.raises(ValueError, match=match):
            decoder.encoding_matrix(wavelength=6.0)

    def test_decode_returns_four_workspaces(self, full_pol_workspaces, clean_workspace):
        device_cross_sections, _, config = full_pol_workspaces
        decoder = FullPolarizationDecoder(config)

        result = decoder.decode(device_cross_sections)
        for ws in result:
            clean_workspace(ws)

        assert len(result) == 4

    def test_decode_output_logs_have_polarization_state(self, full_pol_workspaces, clean_workspace):
        device_cross_sections, _, config = full_pol_workspaces
        decoder = FullPolarizationDecoder(config)

        result = decoder.decode(device_cross_sections)
        for ws in result:
            clean_workspace(ws)

        assert PolarizationState.get(result[0]) == PolarizationState.UP_UP
        assert PolarizationState.get(result[1]) == PolarizationState.UP_DOWN
        assert PolarizationState.get(result[2]) == PolarizationState.DOWN_UP
        assert PolarizationState.get(result[3]) == PolarizationState.DOWN_DOWN
        for ws in result:
            assert PolarizationCrossSection.logname not in SampleLogs(ws)

    def test_decode_wrong_count_raises(self, full_pol_workspaces):
        device_cross_sections, _, config = full_pol_workspaces
        decoder = FullPolarizationDecoder(config)

        with pytest.raises(ValueError, match="exactly 4 device cross-sections"):
            decoder.decode(device_cross_sections[:3])
        with pytest.raises(ValueError, match="exactly 4 device cross-sections"):
            decoder.decode(device_cross_sections + device_cross_sections)

    def test_decode_order_invariant(self, full_pol_workspaces, clean_workspace):
        device_cross_sections, _, config = full_pol_workspaces
        decoder = FullPolarizationDecoder(config)

        result_normal = decoder.decode(device_cross_sections)
        result_reversed = decoder.decode(list(reversed(device_cross_sections)))
        for ws in result_normal + result_reversed:
            clean_workspace(ws)

        for normal, reversed_ in zip(result_normal, result_reversed):
            np.testing.assert_array_almost_equal(normal.readY(0), reversed_.readY(0))

    def test_decode_intensity_known_case(self, full_pol_workspaces, clean_workspace):
        device_cross_sections, spin_values, config = full_pol_workspaces
        decoder = FullPolarizationDecoder(config)

        result = decoder.decode(device_cross_sections)
        for ws in result:
            clean_workspace(ws)

        for ws, expected in zip(result, spin_values):
            np.testing.assert_array_almost_equal(ws.readY(0), [expected])


if __name__ == "__main__":
    pytest.main([__file__])
