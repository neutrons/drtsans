#!/usr/bin/env python
from __future__ import print_function

from mantid.simpleapi import LoadHFIRSANS
from ornl.sans.detector import Detector
from ornl.settings import unique_workspace_name


def test_detector_biosans(biosans_f):

    ws = LoadHFIRSANS(
        Filename=biosans_f['beamcenter'],
        OutputWorkspace=unique_workspace_name())
    det = Detector(ws)
    assert 2 == det.get_number_of_monitors()
    assert (192, 256) == det.get_detector_dimensions("detector1")
    assert (160, 256) == det.get_detector_dimensions("wing_detector")


def test_detector_cache_biosans(biosans_f):

    import time
    import math

    ws = LoadHFIRSANS(
        Filename=biosans_f['beamcenter'],
        OutputWorkspace=unique_workspace_name())

    det = Detector(ws)

    t_mons_1_start = time.time()
    det.get_number_of_monitors()
    t_mons_1_stop = time.time()
    t_mons_1_elapsed = t_mons_1_stop - t_mons_1_start

    t_mons_2_start = time.time()
    det.get_number_of_monitors()
    t_mons_2_stop = time.time()
    t_mons_2_elapsed = t_mons_2_stop - t_mons_2_start

    # The 2nd time should be at least 2 orders of magnitude faster
    assert abs(int(math.log10(t_mons_1_elapsed - t_mons_2_elapsed))) >= 2

    t_det_1_start = time.time()
    det.get_detector_dimensions("detector1")
    t_det_1_stop = time.time()

    t_det_2_start = time.time()
    det.get_detector_dimensions("detector1")
    t_det_2_stop = time.time()
    assert abs(int(math.log10(t_det_1_stop - t_det_1_start) -
                             (t_det_2_stop - t_det_2_start))) >= 4


def test_detector_gpsans(gpsans_f):

    ws = LoadHFIRSANS(
        Filename=gpsans_f['beamcenter'],
        OutputWorkspace=unique_workspace_name())
    det = Detector(ws)
    assert 2 == det.get_number_of_monitors()
    assert (192, 256) == det.get_detector_dimensions("detector1")
