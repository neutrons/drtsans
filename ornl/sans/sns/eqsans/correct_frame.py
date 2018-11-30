from __future__ import (absolute_import, division, print_function)

import numpy as np

from mantid.simpleapi import ConvertUnits, Rebin

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans.chopper import EQSANSDiskChopperSet
from ornl.sans.frame_mode import FrameMode
from ornl.sans import wavelength as wlg
from ornl.settings import namedtuplefy
from ornl.sans.geometry import source_detector_distance


@namedtuplefy
def transmitted_bands(ws):
    r"""
    Wavelength bands of the lead and skipped pulses transmitted by
    the choppers

    Parameters
    ----------
    ws: MatrixWorkspace
        Input Workspace containing all necessary info

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - lead, WBand object for the wavelength band of the lead pulse
        - skipped, Wband for the skipped pulse. None if not operating in
            the skipped frame mode
    """
    sl = SampleLogs(ws)
    pulse_period = 1.e6 / sl.frequency.value.mean()  # 10^6/60 micro-seconds
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers
    # Wavelength band of neutrons from the leading pulse transmitted
    # by the chopper system
    lead_band = ch.transmission_bands(pulsed=True)[0]
    # Wavelength from the previous, skipped pulse.
    skip_band = ch.transmission_bands(delay=pulse_period, pulsed=True)[0]\
        if ch.frame_mode == FrameMode.skip else None
    return dict(lead=lead_band, skip=skip_band)


@namedtuplefy
def limiting_tofs(ws, sdd):
    r"""
    Minimimum and maximum TOF's for neutrons coming from the lead and skip
    pulses, assuming a detector were placed at `distance`.

    Parameters
    ----------
    ws: MatrixWorkspace
            Input Workspace containing all necessary info
    sdd: float
        Distance from source to detector, in meters

    Returns
    -------
    namedtuple
        Fields of the namedtuple
        lead: tuple, tof_min and tof_max for the neutrons of the lead pulse
        skip: tuple, tof_min and tof_max for the neutrons of the skip
            pulse. `None` if the frame mode is 'non-skipping'
    """
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers
    bands = transmitted_bands(ws)
    lead = (wlg.tof(bands.lead.min, sdd, ch.pulse_width),
            wlg.tof(bands.lead.max, sdd))
    skip = None if ch.frame_mode == FrameMode.not_skip else \
        (wlg.tof(bands.skip.min, sdd, ch.pulse_width),
         wlg.tof(bands.skip.max, sdd))
    return dict(lead=lead, skip=skip)


@namedtuplefy
def transmitted_bands_clipped(ws, sdd, ltc, htc, interior_clip=False):
    r"""
    Wavelength bands of the lead and skipped pulses transmitted by
    the choppers taking into account the TOF clippings for neutrons
    arriving at the detector, assuming a detector were placed at `distance`.

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    ltc: float
        trim neutrons of the leading pulse with a TOF smaller than the
        minimal TOF plus this value. Units in micro-seconds
    sdd: float
        Distance from source to detector, in meters
    htc: float
        trim neutrons of the leading pulse with TOF bigger than the maximal
        TOF minus this value.  Units in micro-seconds
    interior_clip: False
        If True, also trim slow neutrons from the lead pulse (using `htc`) and
        fast neutrons from the skip pulse (using `ltc`)
    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - lead, WBand object for the wavelength band of the lead pulse
        - skipped, Wband for the skipped pulse. None if not operating in
            the skipped frame mode
    """
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers
    lwc = wlg.from_tof(ltc, sdd, ch.pulse_width)  # low wavelength clip
    hwc = wlg.from_tof(htc, sdd)  # high wavelength clip
    bands = transmitted_bands(ws)
    if ch.frame_mode == FrameMode.not_skip:
        lead = wlg.Wband(bands.lead.min + lwc, bands.lead.max - hwc)
        skip = None
    else:
        l_i = bands.lead.max - hwc if interior_clip is True else bands.lead.max
        l_a = bands.skip.min + lwc if interior_clip is True else bands.skip.min
        lead = wlg.Wband(bands.lead.min + lwc, l_i)
        skip = wlg.Wband(l_a, bands.skip.max - hwc)
    return dict(lead=lead, skip=skip)


def correct_frame(ws, s2c):
    r"""
    Assign the correct TOF to each event.

    TOFS may need to be shifted one (or more) frame widths and if working
    in frame skipped mode, one pulse period.

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    s2c: float
        Distance from source to detecting component (detector1, monitor) in
        meters
    ltc: float
        trim neutrons of the leading pulse with a TOF smaller than the
        minimal TOF plus this value
    htc: float
        trim neutrons f the leading pulse with TOF bigger than the maxima
         TOF minus this value.

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - ws: EventsWorkspace, the input workspace
        - lead_band: WBand wavelength band for the lead pulse
        - skip_band: WBand wavelength band for the skip pulse. `None` if not
            working in frame-skipping mode
    """
    sl = SampleLogs(ws)
    pulse_period = 1.e6 / sl.frequency.value.mean()  # 10^6/60 micro-seconds
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers
    # The TOF values recorded are never be bigger than the frame width,
    # which is also the choppers' rotational period.
    frame_width = ch.period  # either 10^6/60 or 10^6/30 micro-seconds

    # Find how many frame widths elapsed from the time the neutrons of the
    # lead pulse were emitted and the time the neutrons arrived to the
    # detector bank. This time must be added to the stored TOF values
    tof_min, tof_max = limiting_tofs(ws, s2c).lead
    frames_offset_time = frame_width * int(tof_min / frame_width)
    for i in range(ws.getNumberHistograms()):
        sp = ws.getSpectrum(i)
        tofs = sp.getTofs()
        if len(tofs) == 0:
            continue  # no events found in this detector
        tofs += frames_offset_time  # shift times to the correct frame
        # TOF values smaller than that of the fastest neutrons have been
        # 'folded' by the data acquisition system. They must be shifted
        tofs[np.where(tofs < tof_min)] += frame_width
        # Events from the skipped pulse are delayed by one pulse period
        if ch.frame_mode == FrameMode.skip:
            tofs[np.where(tofs > tof_min + pulse_period)] += pulse_period
        pulse_times = list(sp.getPulseTimes())
        sp.clear(False)
        for tof, pt in zip(tofs, pulse_times):
            sp.addEventQuickly(tof, pt)
    # Amend the logs
    sl.is_frame_skipping = True if ch.frame_mode == FrameMode.skip else False


def correct_detector_frame(ws):
    correct_frame(ws, source_detector_distance(ws, units='m'))


def convert_to_wavelength(ws, bands, bin_width, out_ws, events=False):
    r"""
    Convert a time-of-fligth events workspace to a wavelength workpsace

    Parameters
    ----------
    ws: EventsWorkspace
        Input workspace, after TOF frame correction
    bands: namedtuple
        Output of running `transmitted_bands_clipped` on the workspace
    bin_width: float
        Bin width in Angstroms
    out_ws: str
        Name of the output workspace
    events: bool
        Do we preserve events?
    Returns
    -------
    MatrixWorspace
        EventsWorkspace or Matrix2DWorkspace
    """
    fm = (EQSANSDiskChopperSet(ws).frame_mode == FrameMode.skip)  # frame skip?

    # Convert to Wavelength and rebin
    _ws = ConvertUnits(ws, Target='Wavelength', Emode='Elastic')
    w_min = bands.lead.min
    w_max = bands.lead.max if fm is False else bands.skip.max
    _ws = Rebin(_ws, Params=[w_min, bin_width, w_max],
                PreserveEvents=events, OutputWorkspace=out_ws)

    # Discard neutrons in between bands.lead.max and bands.skip.min
    if fm is True:
        to_zero = np.where((_ws.dataX(0) > bands.lead.max) &
                           (_ws.dataX(0) < bands.skip.min))[0]
        for i in range(_ws.getNumberHistograms()):
            _ws.dataY(i)[to_zero] = 0.0
            _ws.dataE(i)[to_zero] = 0.0

    # Insert bands information in the logs
    sl = SampleLogs(_ws)
    sl.wavelength_min, sl.wavelength_max = w_min, w_max
    if fm is True:
        sl.wavelength_lead_min = bands.lead.min
        sl.wavelength_lead_max = bands.lead.max
        sl.wavelength_skip_min = bands.skip.min
        sl.wavelength_skip_max = bands.skip.max

    return _ws
