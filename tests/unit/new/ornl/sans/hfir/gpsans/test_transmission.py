import pytest
from pytest import approx
from mantid.simpleapi import (LoadHFIRSANS, MoveInstrumentComponent,
                              CreateWorkspace)
from ornl.settings import unique_workspace_dundername as uwd
from ornl.sans.transmission import (calculate_transmission,
                                    apply_transmission_correction)


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


def test_calculate_transmission(gpsans_full_dataset, sample_scattering_sum_ws,
                                dataset_center):
    x, y = dataset_center[0], dataset_center[1]
    input_sample_ws = sample_scattering_sum_ws
    MoveInstrumentComponent(
        Workspace=input_sample_ws, ComponentName='detector1', X=-x, Y=-y)

    input_reference_ws = LoadHFIRSANS(
        Filename=gpsans_full_dataset['sample_transmission'])

    MoveInstrumentComponent(
        Workspace=input_reference_ws, ComponentName='detector1', X=-x, Y=-y)

    trans = calculate_transmission(input_sample_ws, input_reference_ws)

    assert trans.readY(0)[0] == approx(0.08224, abs=1e-4)
    assert trans.readE(0)[0] == approx(0.01267, abs=1e-4)


def test_apply_transmission_correction(gpsans_full_dataset,
                                       sample_scattering_sum_ws,
                                       dataset_center):
    trans_value = 0.5191
    trans_ws = CreateWorkspace(
        DataX=[3.8, 4.2],
        DataY=[trans_value],
        DataE=[0.0141],
        UnitX="Wavelength"
    )

    ws = sample_scattering_sum_ws
    x, y = dataset_center[0], dataset_center[1]
    MoveInstrumentComponent(Workspace=ws, ComponentName='detector1',
                            X=-x, Y=-y)

    ws_c = apply_transmission_correction(ws, trans_workspace=trans_ws,
                                         theta_dependent=False,
                                         output_workspace=uwd())

    assert ws.readY(9100)[0] == approx(20.0, abs=1e-3)
    assert ws_c.readY(9100)[0] == approx(20.0 / trans_value, abs=1e-3)


def test_apply_transmission_with_values(gpsans_full_dataset, dataset_center,
                                        sample_scattering_sum_ws):
    trans_value = 0.5191
    trans_error = 0.0141

    ws = sample_scattering_sum_ws
    x, y = dataset_center[0], dataset_center[1]
    MoveInstrumentComponent(Workspace=ws, ComponentName='detector1',
                            X=-x, Y=-y)

    ws_c = apply_transmission_correction(ws, trans_value=trans_value,
                                         trans_error=trans_error,
                                         theta_dependent=False,
                                         output_workspace=uwd())

    assert ws.readY(9100)[0] == approx(20.0, abs=1e-3)
    assert ws_c.readY(9100)[0] == approx(20.0 / trans_value, abs=1e-3)


def test_apply_transmission_correction_ws(gpsans_full_dataset, dataset_center,
                                          sample_scattering_sum_ws):
    ws = sample_scattering_sum_ws
    x, y = dataset_center[0], dataset_center[1]
    MoveInstrumentComponent(Workspace=ws, ComponentName='detector1',
                            X=-x, Y=-y)

    trans_ws = CreateWorkspace(DataX=[0, 1], DataY=[0.08224400871459694],
                               DataE=[0.012671053121947698])

    ws_c = apply_transmission_correction(ws, trans_ws, theta_dependent=False,
                                         output_workspace=uwd())

    assert ws_c.readY(9100)[0] == approx(243.178, abs=1e-2)
    assert ws_c.readE(9100)[0] == approx(58.312, abs=1e-2)


def test_apply_transmission_correction_value(gpsans_full_dataset,
                                             sample_scattering_sum_ws):
    ws = sample_scattering_sum_ws
    # Zero angle transmission values
    trans_value = 0.08224400871459694
    trans_error = 0.012671053121947698

    ws_c = apply_transmission_correction(ws, trans_value=trans_value,
                                         trans_error=trans_error,
                                         theta_dependent=False,
                                         output_workspace=uwd())
    # Note the corrected values are the same as above
    assert ws_c.readY(9100)[0] == approx(243.178, abs=1e-2)
    assert ws_c.readE(9100)[0] == approx(58.312, abs=1e-2)


if __name__ == '__main__':
    pytest.main()
