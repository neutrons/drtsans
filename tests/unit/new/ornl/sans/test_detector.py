#!/usr/bin/env python
from __future__ import print_function

from mantid.simpleapi import LoadEmptyInstrument
from ornl.sans.detector import Component
from ornl.settings import unique_workspace_name


def test_detector_biosans():

    ws = LoadEmptyInstrument(InstrumentName='biosans',
                             OutputWorkspace=unique_workspace_name())

    d = Component(ws, "detector1")
    assert 192 == d.dim_x
    assert 256 == d.dim_y
    assert 192*256 == d.dims
    assert 2 == d.first_index

    d = Component(ws, "wing_detector")
    assert 160 == d.dim_x
    assert 256 == d.dim_y
    assert 160*256 == d.dims
    assert 192*256 + 2 == d.first_index


def test_detector_gpsans():

    ws = LoadEmptyInstrument(InstrumentName='cg2',
                             OutputWorkspace=unique_workspace_name())

    d = Component(ws, "detector1")
    assert 192 == d.dim_x
    assert 256 == d.dim_y
    assert 192*256 == d.dims
    assert 2 == d.first_index


def test_detector_eqsans():

    ws = LoadEmptyInstrument(InstrumentName='eqsans',
                             OutputWorkspace=unique_workspace_name())

    d = Component(ws, "detector1")
    assert 192 == d.dim_x
    assert 256 == d.dim_y
    assert 192*256 == d.dims
    assert 1 == d.first_index
