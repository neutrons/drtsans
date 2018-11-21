#!/usr/bin/env python
from __future__ import print_function

import pytest


def _instrument_geometry(gpsans_f, workspace_to_center):
    '''
    Finds the beamcenter and places the instrument in the right position.
    '''

    from ornl.sans.hfir.gpsans import beam_finder
    from mantid import mtd
    from mantid.simpleapi import (
        MoveInstrumentComponent, LoadSpice2D)

    ws_name = "__beamcenter"
    LoadSpice2D(Filename=gpsans_f['beamcenter_off_setted'],
                OutputWorkspace=ws_name)
    ws = mtd[ws_name]
    x, y = beam_finder.direct_beam_center(ws)

    MoveInstrumentComponent(
        Workspace=workspace_to_center, ComponentName='detector1', X=-x, Y=-y)


def test_calculate_transmission(gpsans_f):
    '''

    '''
    from ornl.sans.transmission import (zero_angle_transmission,
                                        calculate_radius_from_input_ws)
    from mantid.simpleapi import LoadSpice2D
    from mantid import mtd

    input_sample_ws_mame = 'input_sample_ws_name'
    LoadSpice2D(Filename=gpsans_f['sample_scattering'],
                OutputWorkspace=input_sample_ws_mame)
    input_sample_ws = mtd[input_sample_ws_mame]
    _instrument_geometry(gpsans_f, input_sample_ws)
    input_reference_ws_name = 'input_reference_ws_name'
    LoadSpice2D(Filename=gpsans_f['sample_transmission'],
                OutputWorkspace=input_reference_ws_name)
    input_reference_ws = mtd[input_reference_ws_name]
    _instrument_geometry(gpsans_f, input_reference_ws)
    output_ws_name = "__test_out"

    radius = calculate_radius_from_input_ws(input_sample_ws)
    calculated_transmission = zero_angle_transmission(
        input_sample_ws, input_reference_ws, radius,
        output_ws_name, delete_temp_wss=True).transmission
    assert calculated_transmission.readY(
        0)[0] == pytest.approx(0.0014, abs=1e-4)
    assert calculated_transmission.readE(
        0)[0] == pytest.approx(0.0035, abs=1e-4)


def test_apply_transmission_with_ws(gpsans_f):
    '''
    '''

    from ornl.sans.transmission import apply_transmission
    from mantid.simpleapi import LoadSpice2D, CreateWorkspace
    from mantid import mtd

    trans_ws_name = "_transmission"
    CreateWorkspace(
        OutputWorkspace=trans_ws_name,
        DataX=[3.8, 4.2],
        DataY=[0.5191],
        DataE=[0.0141],
        UnitX="Wavelength"
    )
    trans_ws = mtd[trans_ws_name]

    ws_sample_name = 'ws_sample'
    LoadSpice2D(Filename=gpsans_f['sample_scattering'],
                OutputWorkspace=ws_sample_name)
    ws_sample = mtd[ws_sample_name]

    ws_sample_corrected_name = 'ws_sample_corrected_name'
    apply_transmission(ws_sample, ws_sample_corrected_name, trans_ws=trans_ws)
    ws_sample_corrected = mtd[ws_sample_corrected_name]

    assert ws_sample.readY(9100)[0] == pytest.approx(3.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(5.8557131, abs=1e-3)


def test_apply_transmission_with_values(gpsans_f):
    '''
    '''

    from ornl.sans.transmission import apply_transmission
    from mantid.simpleapi import LoadSpice2D
    from mantid import mtd

    trans_value = 0.5191
    trans_error = 0.0141

    ws_sample_name = 'ws_sample'
    LoadSpice2D(Filename=gpsans_f['sample_scattering'],
                OutputWorkspace=ws_sample_name)
    ws_sample = mtd[ws_sample_name]

    ws_sample_corrected_name = 'ws_sample_corrected_name'
    apply_transmission(ws_sample, ws_sample_corrected_name,
                       trans_value=trans_value, trans_error=trans_error)
    ws_sample_corrected = mtd[ws_sample_corrected_name]

    assert ws_sample.readY(9100)[0] == pytest.approx(3.0, abs=1e-3)
    assert ws_sample_corrected.readY(
        9100)[0] == pytest.approx(5.8557131, abs=1e-3)


if __name__ == '__main__':
    pytest.main()
