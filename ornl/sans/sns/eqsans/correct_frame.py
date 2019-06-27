from __future__ import (absolute_import, division, print_function)

import numpy as np

from mantid.simpleapi import (ConvertUnits, Rebin, EQSANSCorrectFrame)

from ornl.settings import (optional_output_workspace,
                           unique_workspace_dundername as uwd)
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
    the choppers of the workspace.

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
    pulse_period = 1.e6 / sl.single_value('frequency')  # 10^6/60 micro-seconds
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers
    # Wavelength band of neutrons from the leading pulse transmitted
    # by the chopper system
    lead_band = ch.transmission_bands(pulsed=True)[0]
    # Wavelength from the previous, skipped pulse.
    skip_band = ch.transmission_bands(delay=pulse_period, pulsed=True)[0]\
        if ch.frame_mode == FrameMode.skip else None
    return dict(lead=lead_band, skip=skip_band)


@namedtuplefy
def clipped_bands_from_logs(ws):
    r"""
    Retrieve the wavelength bands over which we expect non-zero intensity. We
    inspect the log.

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
    lead = wlg.Wband(sl.wavelength_lead_min.value,
                     sl.wavelength_lead_max.value)
    if bool(sl.is_frame_skipping.value) is True:
        skip = wlg.Wband(sl.wavelength_skip_min.value,
                         sl.wavelength_skip_max.value)
    else:
        skip = None
    return dict(lead=lead, skip=skip)


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
def transmitted_bands_clipped(ws, sdd, low_tof_clip, high_tof_clip,
                              interior_clip=False):
    r"""
    Wavelength bands of the lead and skipped pulses transmitted by
    the choppers taking into account the TOF clippings for neutrons
    arriving at the detector, assuming a detector were placed at `distance`.

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    sdd: float
        Distance from source to detector, in meters
    low_tof_clip: float
        trim neutrons of the leading pulse with a TOF smaller than the
        minimal TOF plus this value. Units in micro-seconds
    high_tof_clip: float
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
    lwc = wlg.from_tof(low_tof_clip, sdd, ch.pulse_width)  # low wavel. clip
    hwc = wlg.from_tof(high_tof_clip, sdd)  # high wavelength clip
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


def log_tof_structure(ws, low_tof_clip, high_tof_clip, interior_clip=False):
    r"""
    Append to the logs relevant information about the time of flight
    frame and structure

    Insert clipping times

    Parameters
    ----------
    ws: MatrixWorkspace
        InputWorkspace
    sdd: float
        Distance from source to detector, in meters
    low_tof_clip: float
        trim neutrons of the leading pulse with a TOF smaller than the
        minimal TOF plus this value. Units in micro-seconds
    high_tof_clip: float
        trim neutrons of the leading pulse with TOF bigger than the maximal
        TOF minus this value.  Units in micro-seconds
    interior_clip: False
        If True, also trim slow neutrons from the lead pulse (using `htc`) and
        fast neutrons from the skip pulse (using `ltc`)
    """
    ch = EQSANSDiskChopperSet(ws)
    sl = SampleLogs(ws)
    sl.insert('tof_frame_width', ch.period, unit='ms')
    clip_times = 1 if interior_clip is False else 2
    tof_width_clipped = ch.period - clip_times * (low_tof_clip + high_tof_clip)
    sl.insert('tof_frame_width_clipped', tof_width_clipped, unit='ms')


def correct_frame(ws, source_to_component_distance):
    r"""
    Assign the correct TOF to each event.

    TOFS may need to be shifted one (or more) frame widths and if working
    in frame skipped mode, one pulse period.

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    source_to_component_distance: float
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

    tof_min, tof_max = limiting_tofs(ws, source_to_component_distance).lead
    EQSANSCorrectFrame(ws,
                       PulsePeriod=pulse_period,
                       MinTOF=tof_min,
                       FrameWidth=frame_width,
                       FrameSkipping=(ch.frame_mode is FrameMode.skip))
    # Amend the logs
    fr_skip = 1 if ch.frame_mode == FrameMode.skip else 0
    sl.insert('is_frame_skipping', fr_skip)


def correct_detector_frame(ws):
    correct_frame(ws, source_detector_distance(ws, units='m'))


def band_gap_indexes(ws, bands):
    r"""
    Convention to define the indexes of the gap band.

    For runs in skipped frame mode, there is a wavelength band in between
    the bands of the lead and skipped pulse. This range has zero neutron
    counts. This function returns the indexes of `ws.dataY(0)` array
    corresponding to the band gap.

    Parameters
    ----------
    ws: MatrixWorkspace
        Input workspace with units of Wavelength on the X-axis
    bands: namedtuple
        Output of running `transmitted_bands_clipped` on the workspace
    Returns
    -------
    list
        Indexes of array `ws.dataY(i)` where intensity is zero. Empty list if
        working frame mode is not skipped-mode
    """
    if bands.skip is None:
        return list()
    else:
        return (np.where((ws.dataX(0) > bands.lead.max) &
                         (ws.dataX(0) < bands.skip.min))[0]).tolist()


@optional_output_workspace
def convert_to_wavelength(ws, bands, bin_width, events=False):
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
    events: bool
        Do we preserve events?
    Returns
    -------
    MatrixWorspace
        EventsWorkspace or Matrix2DWorkspace
    """
    fm = (EQSANSDiskChopperSet(ws).frame_mode == FrameMode.skip)  # frame skip?

    # Convert to Wavelength and rebin
    _ws = ConvertUnits(ws, Target='Wavelength', Emode='Elastic',
                       OutputWorkspace=uwd())
    w_min = bands.lead.min
    w_max = bands.lead.max if fm is False else bands.skip.max
    _ws = Rebin(_ws, Params=[w_min, bin_width, w_max],
                PreserveEvents=events, OutputWorkspace=_ws.name())

    # Discard neutrons in between bands.lead.max and bands.skip.min
    if fm is True:
        to_zero = band_gap_indexes(_ws, bands)
        for i in range(_ws.getNumberHistograms()):
            _ws.dataY(i)[to_zero] = 0.0
            _ws.dataE(i)[to_zero] = 1.0

    # Insert bands information in the logs
    sl = SampleLogs(_ws)
    sl.wavelength_min, sl.wavelength_max = w_min, w_max
    sl.insert('wavelength_lead_min', bands.lead.min, unit='Angstrom')
    sl.insert('wavelength_lead_max', bands.lead.max, unit='Angstrom')
    if fm is True:
        sl.insert('wavelength_skip_min', bands.skip.min, unit='Angstrom')
        sl.insert('wavelength_skip_max', bands.skip.max, unit='Angstrom')
    return _ws
