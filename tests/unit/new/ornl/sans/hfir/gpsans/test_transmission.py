#!/usr/bin/env python
from __future__ import print_function

import pytest


def _instrument_geometry(gpsans_f, workspace_to_center):
    '''
    Finds the beamcenter and places the instrument in the right position.
    '''

    from ornl.sans.hfir.gpsans import beam_finder
    from mantid.simpleapi import (
        MoveInstrumentComponent, LoadHFIRSANS)

    __beamcenter = LoadHFIRSANS(Filename=gpsans_f['beamcenter_off_setted'])
    x, y = beam_finder.direct_beam_center(__beamcenter)

    MoveInstrumentComponent(
        Workspace=workspace_to_center, ComponentName='detector1', X=-x, Y=-y)


@pytest.mark.offline
def test_calculate_transmission(gpsans_f):
    '''

    '''
    from ornl.sans.transmission import (zero_angle_transmission,
                                        calculate_radius_from_input_ws)
    from mantid.simpleapi import LoadHFIRSANS

    input_sample_ws = LoadHFIRSANS(Filename=gpsans_f['sample_scattering'])

    _instrument_geometry(gpsans_f, input_sample_ws)

    input_reference_ws = LoadHFIRSANS(Filename=gpsans_f['sample_transmission'])

    _instrument_geometry(gpsans_f, input_reference_ws)
    radius = calculate_radius_from_input_ws(input_sample_ws)

    calculated_transmission = zero_angle_transmission(
        input_sample_ws, input_reference_ws, radius)
    assert calculated_transmission.readY(
        0)[0] == pytest.approx(0.0014, abs=1e-4)
    assert calculated_transmission.readE(
        0)[0] == pytest.approx(0.0035, abs=1e-4)


@pytest.mark.offline
def test_apply_transmission_with_ws(gpsans_f):
    '''
    '''

    from ornl.sans.transmission import apply_transmission_mantid
    from mantid.simpleapi import LoadHFIRSANS, CreateWorkspace

    trans_value = 0.5191
    trans_ws = CreateWorkspace(
        DataX=[3.8, 4.2],
        DataY=[trans_value],
        DataE=[0.0141],
        UnitX="Wavelength"
    )

    ws_sample = LoadHFIRSANS(Filename=gpsans_f['sample_scattering'])

    ws_sample_corrected = apply_transmission_mantid(
        ws_sample, trans_ws=trans_ws, theta_dependent=False)

    assert ws_sample.readY(9100)[0] == pytest.approx(3.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(3.0 / trans_value, abs=1e-3)


@pytest.mark.offline
def test_apply_transmission_with_values(gpsans_f):
    '''
    '''

    from ornl.sans.transmission import apply_transmission_mantid
    from mantid.simpleapi import LoadHFIRSANS

    trans_value = 0.5191
    trans_error = 0.0141

    ws_sample = LoadHFIRSANS(Filename=gpsans_f['sample_scattering'])

    ws_sample_corrected = apply_transmission_mantid(
        ws_sample, trans_value=trans_value, trans_error=trans_error,
        theta_dependent=False)

    assert ws_sample.readY(9100)[0] == pytest.approx(3.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(3.0 / trans_value, abs=1e-3)


@pytest.mark.offline
def test_apply_transmission_correction(gpsans_f):
    '''
    This is the function that users / scientists use 90% of the time
    '''
    from ornl.sans.transmission import apply_transmission_correction
    from mantid.simpleapi import LoadHFIRSANS

    input_sample_ws = LoadHFIRSANS(Filename=gpsans_f['sample_scattering'])
    input_reference_ws = LoadHFIRSANS(Filename=gpsans_f['sample_transmission'])

    input_sample_corrected_ws = apply_transmission_correction(
        input_sample_ws, input_reference_ws, theta_dependent=False)

    # I'm not sure if this is correct!
    assert input_sample_corrected_ws.readY(9100)[0] == \
        pytest.approx(22.45, abs=1e-2)


if __name__ == '__main__':
    pytest.main()
