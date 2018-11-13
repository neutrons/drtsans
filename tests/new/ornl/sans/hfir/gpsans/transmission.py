#!/usr/bin/env python
from __future__ import print_function

import pytest


def _instrument_geometry(gpsans_f):
    '''
    Finds the beamcenter and places the instrument in the right position.
    '''

    from ornl.sans.hfir.gpsans import beam_finder
    from mantid import mtd
    from mantid.simpleapi import (
        MoveInstrumentComponent)

    x, y = beam_finder.direct_beam_center(gpsans_f['beamcenter'])
    ws = mtd['ws']  # This comes from direct_beam_center function
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)


def test_calculate_transmission(gpsans_f):
    '''

    '''
    from ornl.sans.transmission import calculate_transmission
    from mantid.simpleapi import (
        LoadSpice2D)

    _instrument_geometry(gpsans_f)

    # calculate_transmission(input_sample_ws, input_reference_ws, output_ws,
    #                           radius=None, delete_temp_wss=True)

    input_sample_ws = LoadSpice2D(gpsans_f['sample_scattering'])
    input_reference_ws = LoadSpice2D(gpsans_f['sample_transmission'])
    output_ws = "__test_out"

    calculate_transmission(input_sample_ws, input_reference_ws, output_ws,
                           radius=None, delete_temp_wss=True)


def test_apply_transmission(gpsans_f):
    '''
    '''
    _instrument_geometry(gpsans_f)
    # TODO
