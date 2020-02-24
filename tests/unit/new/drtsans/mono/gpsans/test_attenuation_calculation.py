#!/usr/bin/env python
import pytest
from drtsans.mono.gpsans import attenuation_factor
from mantid.simpleapi import AddSampleLog, CreateSampleWorkspace


def test_attenuation_factor():
    # Test input and expected values provided by Lisa Debeer-Schmitt
    wavelength = 4.75
    attenuator = 6  # x2k
    expected_value = 0.0010404484927383326
    expected_error = 0.0004406323142152028

    # Create workspace and add sample logs
    ws_test_attenuation_factor = CreateSampleWorkspace()
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='wavelength',
                 LogText='{}'.format(wavelength),
                 LogType='Number Series',
                 LogUnit='A')
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='attenuator',
                 LogText='{}'.format(attenuator),
                 LogType='Number Series')

    value, error = attenuation_factor(ws_test_attenuation_factor)
    assert value == pytest.approx(expected_value)
    assert error == pytest.approx(expected_error)


def test_attenuation_factor_open_close():
    # Create workspace and add sample logs
    ws_test_attenuation_factor = CreateSampleWorkspace()
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='wavelength',
                 LogText='1.54',
                 LogType='Number Series',
                 LogUnit='A')

    attenuator = 0  # Undefined
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='attenuator',
                 LogText='{}'.format(attenuator),
                 LogType='Number Series')
    assert attenuation_factor(ws_test_attenuation_factor) == (1, 0)

    attenuator = 1  # Close
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='attenuator',
                 LogText='{}'.format(attenuator),
                 LogType='Number Series')
    assert attenuation_factor(ws_test_attenuation_factor) == (1, 0)

    attenuator = 2  # Open
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='attenuator',
                 LogText='{}'.format(attenuator),
                 LogType='Number Series')
    assert attenuation_factor(ws_test_attenuation_factor) == (1, 0)


def test_attenuation_factor_missing_logs():
    # Create workspace
    ws_test_attenuation_factor = CreateSampleWorkspace()

    # Missing attenuator log
    with pytest.raises(RuntimeError):
        attenuation_factor(ws_test_attenuation_factor)

    # Add attenuator log
    AddSampleLog(Workspace=ws_test_attenuation_factor,
                 LogName='attenuator',
                 LogText='4',
                 LogType='Number Series')

    # Missing wavelength log
    with pytest.raises(RuntimeError):
        attenuation_factor(ws_test_attenuation_factor)


if __name__ == '__main__':
    pytest.main([__file__])
