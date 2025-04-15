# third-party imports
from mantid.simpleapi import CreateSingleValuedWorkspace
import pytest

# drtsans imports
from drtsans import half_polarization
from drtsans.polarization import _calc_flipping_ratio, SimulatedLogs


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
        log = SimulatedLogs(polarizer_flipper="heartbeat", analyzer_flipper="heartbeat")
        assert log.polarizer_flipper == "heartbeat"
        assert log.analyzer_flipper == "heartbeat"

    def test_invalid_flipper_generators(self):
        with pytest.raises(ValueError) as excinfo:
            SimulatedLogs(polarizer_flipper="invalid_generator")
        assert "The polarizer flipper generator must be one of ['heartbeat']" in str(excinfo.value)
        with pytest.raises(ValueError) as excinfo:
            SimulatedLogs(analyzer_flipper="invalid_generator")
        assert "The analyzer flipper generator must be one of ['heartbeat']" in str(excinfo.value)




    def test_valid_veto_generators(self):
        log = SimulatedLogs(polarizer_veto="binary_pulse", analyzer_veto="binary_pulse")
        assert log.polarizer_veto == "binary_pulse"
        assert log.analyzer_veto == "binary_pulse"

    def test_invalid_veto_generators(self):
        with pytest.raises(ValueError, match="polarizer veto generator must be one of ['binary_pulse']"):
            SimulatedLogs(polarizer_veto="invalid_veto")
        with pytest.raises(ValueError, match="analyzer veto generator must be one of ['binary_pulse']"):
            SimulatedLogs(analyzer_veto="invalid_veto")

    def test_heartbeat_generator(self):
        log = SimulatedLogs()
        generator = log.heartbeat(interval=1.0)
        assert next(generator) == 0
        assert next(generator) == 1.0
        assert next(generator) == 2.0

    def test_binary_pulse_generator(self):
        log = SimulatedLogs()
        generator = log.binary_pulse(interval=2.0, duration=1.0)
        assert next(generator) == 0
        assert next(generator) == 1.5
        assert next(generator) == 2.5
        assert next(generator) == 3.5


if __name__ == "__main__":
    pytest.main([__file__])
