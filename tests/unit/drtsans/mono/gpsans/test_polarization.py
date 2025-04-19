# third-party imports
from mantid.simpleapi import CreateSingleValuedWorkspace, mtd
from numpy.testing import assert_equal, assert_array_almost_equal
import pytest

# drtsans imports
from drtsans import half_polarization
from drtsans.polarization import (
    _calc_flipping_ratio,
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


def test_flipping_ratio(temp_workspace_name, clean_workspace):
    """Test for the calclation of the flipping ratio in section 9.1. This was used to determine
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
        assert "The polarizer flipper generator must be one of ['heartbeat']" in str(excinfo.value)
        with pytest.raises(ValueError) as excinfo:
            SimulatedPolarizationLogs(analyzer_flipper=TimesGeneratorSpecs("invalid_generator", {}))
        assert "The analyzer flipper generator must be one of ['heartbeat']" in str(excinfo.value)

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
        times = SimulatedPolarizationLogs().binary_pulse(interval=3.0, veto_duration=1.0, upper_bound=10)
        assert_array_almost_equal(list(times), [0, 2.5, 3.5, 5.5, 6.5, 8.5, 9.5], decimal=2)
        times = SimulatedPolarizationLogs().binary_pulse(
            interval=3.0, veto_duration=1.0, dead_time=2.7, upper_bound=10
        )
        assert_array_almost_equal(list(times), [3.5, 5.5, 6.5, 8.5, 9.5], decimal=2)

    def test_times_generator(self):
        logs = SimulatedPolarizationLogs(
            polarizer=1,
            polarizer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 1.0}),
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 1.0, "veto_duration": 0.2}),
            analyzer=2,
            analyzer_flipper=TimesGeneratorSpecs("heartbeat", {"interval": 2.0}),
            analyzer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 2.0, "veto_duration": 0.4}),
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
            polarizer_veto=TimesGeneratorSpecs("binary_pulse", {"interval": 60.0, "veto_duration": 1.0}),
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


if __name__ == "__main__":
    pytest.main([__file__])
