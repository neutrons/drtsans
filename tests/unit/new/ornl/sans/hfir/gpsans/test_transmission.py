#!/usr/bin/env python
from __future__ import print_function

import pytest


@pytest.fixture(scope='module')
def sample_scattering_sum_ws(gpsans_full_dataset):
    '''
    Sums the contents of gpsans_full_dataset['sample_scattering_list'] into
    a single WS
    Runs once per module
    '''

    from mantid.simpleapi import LoadHFIRSANS

    filename_list = gpsans_full_dataset['sample_scattering_list']

    __acc = LoadHFIRSANS(Filename=filename_list[0])
    for filename in filename_list[1:]:
        __ws = LoadHFIRSANS(Filename=filename)
        __acc += __ws

    return __acc


def _instrument_geometry(gpsans_full_dataset, workspace_to_center):
    '''
    Finds the beamcenter and places the instrument in the right position.
    '''

    from ornl.sans.hfir.gpsans import beam_finder
    from mantid.simpleapi import MoveInstrumentComponent, LoadHFIRSANS

    __beamcenter = LoadHFIRSANS(
        Filename=gpsans_full_dataset['beamcenter'])
    x, y = beam_finder.direct_beam_center(__beamcenter)

    MoveInstrumentComponent(
        Workspace=workspace_to_center, ComponentName='detector1', X=-x, Y=-y)


@pytest.mark.offline
def test_calculate_transmission(gpsans_full_dataset, sample_scattering_sum_ws):
    '''

    '''
    from ornl.sans.transmission import (zero_angle_transmission,
                                        calculate_radius_from_input_ws)
    from mantid.simpleapi import LoadHFIRSANS

    input_sample_ws = sample_scattering_sum_ws

    _instrument_geometry(gpsans_full_dataset, input_sample_ws)

    input_reference_ws = LoadHFIRSANS(
        Filename=gpsans_full_dataset['sample_transmission'])

    _instrument_geometry(gpsans_full_dataset, input_reference_ws)
    radius = calculate_radius_from_input_ws(input_sample_ws)

    calculated_transmission = zero_angle_transmission(
        input_sample_ws, input_reference_ws, radius)
    assert calculated_transmission.readY(
        0)[0] == pytest.approx(0.0087, abs=1e-4)
    assert calculated_transmission.readE(
        0)[0] == pytest.approx(0.0113, abs=1e-4)


@pytest.mark.offline
def test_apply_transmission_with_ws(gpsans_full_dataset, sample_scattering_sum_ws):
    '''
    '''

    from ornl.sans.transmission import apply_transmission_mantid
    from mantid.simpleapi import CreateWorkspace

    trans_value = 0.5191
    trans_ws = CreateWorkspace(
        DataX=[3.8, 4.2],
        DataY=[trans_value],
        DataE=[0.0141],
        UnitX="Wavelength"
    )

    ws_sample = sample_scattering_sum_ws

    ws_sample_corrected = apply_transmission_mantid(
        ws_sample, trans_ws=trans_ws, theta_dependent=False)

    assert ws_sample.readY(9100)[0] == pytest.approx(20.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(20.0 / trans_value, abs=1e-3)


@pytest.mark.offline
def test_apply_transmission_with_values(gpsans_full_dataset, sample_scattering_sum_ws):
    '''
    '''

    from ornl.sans.transmission import apply_transmission_mantid

    trans_value = 0.5191
    trans_error = 0.0141

    ws_sample = sample_scattering_sum_ws

    ws_sample_corrected = apply_transmission_mantid(
        ws_sample, trans_value=trans_value, trans_error=trans_error,
        theta_dependent=False)

    assert ws_sample.readY(9100)[0] == pytest.approx(20.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(20.0 / trans_value, abs=1e-3)


@pytest.mark.offline
def test_apply_transmission_correction(gpsans_full_dataset, sample_scattering_sum_ws):
    '''
    This is the function that users / scientists use 90% of the time
    '''
    from ornl.sans.transmission import apply_transmission_correction
    from mantid.simpleapi import LoadHFIRSANS

    input_sample_ws = sample_scattering_sum_ws

    input_reference_ws = LoadHFIRSANS(
        Filename=gpsans_full_dataset['sample_transmission'])

    input_sample_corrected_ws = apply_transmission_correction(
        input_sample_ws, input_reference_ws, theta_dependent=False)

    # I'm not sure if this is correct!
    assert input_sample_corrected_ws.readY(9100)[0] == \
        pytest.approx(14.86, abs=1e-2)


if __name__ == '__main__':
    pytest.main()
