#!/usr/bin/env python
from __future__ import print_function

from mantid.api import mtd
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
    from ornl.sans.transmission import calculate_transmission
    from mantid.simpleapi import LoadSpice2D

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

    calculated_transmission = calculate_transmission(
        input_sample_ws, input_reference_ws, output_ws_name,
        radius=None, delete_temp_wss=True)

    assert calculated_transmission.readY(
        0)[0] == pytest.approx(0.5191, abs=1e-4)
    assert calculated_transmission.readE(
        0)[0] == pytest.approx(0.0141, abs=1e-4)


def test_apply_transmission(gpsans_f):
    '''
    '''

    # TODO
    pass
