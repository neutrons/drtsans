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


@pytest.fixture(scope='module')
def dataset_center(gpsans_full_dataset):
    '''
    Finds the beamcenter and places the instrument in the right position.
    '''
    from ornl.sans.hfir.gpsans import beam_finder
    from mantid.simpleapi import LoadHFIRSANS

    __beamcenter = LoadHFIRSANS(
        Filename=gpsans_full_dataset['beamcenter'])
    x, y = beam_finder.direct_beam_center(__beamcenter)
    return x, y


@pytest.mark.offline
def test_calculate_transmission(gpsans_full_dataset, sample_scattering_sum_ws,
                                dataset_center):
    '''

    '''
    from ornl.sans.transmission import (
        calculate_transmission_value, calculate_transmission_ws)
    from mantid.simpleapi import LoadHFIRSANS, MoveInstrumentComponent

    x, y = dataset_center[0], dataset_center[1]
    input_sample_ws = sample_scattering_sum_ws
    MoveInstrumentComponent(
        Workspace=input_sample_ws, ComponentName='detector1', X=-x, Y=-y)

    input_reference_ws = LoadHFIRSANS(
        Filename=gpsans_full_dataset['sample_transmission'])

    MoveInstrumentComponent(
        Workspace=input_reference_ws, ComponentName='detector1', X=-x, Y=-y)

    calculated_transmission_ws = calculate_transmission_ws(
        input_sample_ws, input_reference_ws)

    calculated_transmission_value = calculate_transmission_value(
        input_sample_ws, input_reference_ws)

    assert calculated_transmission_value[0] == pytest.approx(0.08224, abs=1e-4)
    assert calculated_transmission_value[1] == pytest.approx(0.01267, abs=1e-4)

    assert calculated_transmission_ws.readY(0)[0] == \
        pytest.approx(0.08224, abs=1e-4)
    assert calculated_transmission_ws.readE(0)[0] == \
        pytest.approx(0.01267, abs=1e-4)


@pytest.mark.offline
def test_apply_transmission_with_ws(gpsans_full_dataset,
                                    sample_scattering_sum_ws, dataset_center):
    '''
    '''

    from ornl.sans.transmission import _apply_transmission_mantid
    from mantid.simpleapi import CreateWorkspace, MoveInstrumentComponent

    trans_value = 0.5191
    trans_ws = CreateWorkspace(
        DataX=[3.8, 4.2],
        DataY=[trans_value],
        DataE=[0.0141],
        UnitX="Wavelength"
    )

    ws_sample = sample_scattering_sum_ws
    x, y = dataset_center[0], dataset_center[1]
    MoveInstrumentComponent(
        Workspace=ws_sample, ComponentName='detector1', X=-x, Y=-y)

    ws_sample_corrected = _apply_transmission_mantid(
        ws_sample, trans_ws=trans_ws, theta_dependent=False)

    assert ws_sample.readY(9100)[0] == pytest.approx(20.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(20.0 / trans_value, abs=1e-3)


@pytest.mark.offline
def test_apply_transmission_with_values(gpsans_full_dataset, dataset_center,
                                        sample_scattering_sum_ws):
    '''
    '''

    from ornl.sans.transmission import _apply_transmission_mantid
    from mantid.simpleapi import MoveInstrumentComponent

    trans_value = 0.5191
    trans_error = 0.0141

    ws_sample = sample_scattering_sum_ws
    x, y = dataset_center[0], dataset_center[1]
    MoveInstrumentComponent(
        Workspace=ws_sample, ComponentName='detector1', X=-x, Y=-y)

    ws_sample_corrected = _apply_transmission_mantid(
        ws_sample, trans_value=trans_value, trans_error=trans_error,
        theta_dependent=False)

    assert ws_sample.readY(9100)[0] == pytest.approx(20.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(20.0 / trans_value, abs=1e-3)


@pytest.mark.offline
def test_apply_transmission_correction_ws(gpsans_full_dataset, dataset_center,
                                          sample_scattering_sum_ws):
    '''
    This is the function that users / scientists use 90% of the time
    '''
    from ornl.sans.transmission import apply_transmission_correction_ws
    from mantid.simpleapi import MoveInstrumentComponent, CreateWorkspace

    input_sample_ws = sample_scattering_sum_ws
    x, y = dataset_center[0], dataset_center[1]
    MoveInstrumentComponent(
        Workspace=input_sample_ws, ComponentName='detector1', X=-x, Y=-y)

    transmission_ws = CreateWorkspace(
        DataX=[0, 1], DataY=[0.08224400871459694], DataE=[0.012671053121947698]
    )

    input_sample_corrected_ws = apply_transmission_correction_ws(
        input_sample_ws, transmission_ws, theta_dependent=False)

    assert input_sample_corrected_ws.readY(9100)[0] == \
        pytest.approx(243.178, abs=1e-2)
    assert input_sample_corrected_ws.readE(9100)[0] == \
        pytest.approx(58.312, abs=1e-2)


@pytest.mark.offline
def test_apply_transmission_correction_value(gpsans_full_dataset,
                                             sample_scattering_sum_ws):
    '''
    This is the function that users / scientists use 90% of the time
    '''
    from ornl.sans.transmission import apply_transmission_correction_value

    input_sample_ws = sample_scattering_sum_ws

    # Zero angle transmission values as the input_reference_ws above
    trans_value = 0.08224400871459694
    trans_error = 0.012671053121947698

    input_sample_corrected_ws = apply_transmission_correction_value(
        input_sample_ws, trans_value=trans_value, trans_error=trans_error,
        theta_dependent=False)

    # Note the corrected values are the same as above
    assert input_sample_corrected_ws.readY(9100)[0] == \
        pytest.approx(243.178, abs=1e-2)
    assert input_sample_corrected_ws.readE(9100)[0] == \
        pytest.approx(58.312, abs=1e-2)


if __name__ == '__main__':
    pytest.main()
