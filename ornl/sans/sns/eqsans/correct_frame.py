from __future__ import (absolute_import, division, print_function)

import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans.chopper import EQSANSDiskChopperSet
from ornl.sans.frame_mode import FrameMode
from ornl.sans import wavelength as wlg
from ornl.settings import namedtuplefy
from ornl.sans.geometry import sample_source_distance, source_detector_distance


def replace_tofs(ws):
    for i in range(ws.getNumberHistograms()):
        sp = ws.getSpectrum(i)
        tofs = sp.getTofs()
        pts = np.asarray(sp.getPulseTimes())
        sp.clear(False)
        for tof, pt in zip(tofs, pts):
            sp.addEventQuickly(tof, pt)
    return ws


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
def transmitted_bands_clipped(ws, ltc, htc, ltc2=None, htc2=None):
    r"""
    Wavelength bands of the lead and skipped pulses transmitted by
    the choppers taking into account the TOF clipping cutoffs

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    ltc: float
        trim neutrons of the leading pulse with a TOF smaller than the
        minimal TOF plus this value
    htc: float
        trim neutrons f the leading pulse with TOF bigger than the maxima
         TOF minus this value.
    ltc2: float
        as `ltc` but for neutrons of the skipped pulse. If None, then we use
        `ltc`
    htc2: float
        as `htc` but for neutrons of the skipped pulse. If None, then we use
        `htc`

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - lead, WBand object for the wavelength band of the lead pulse
        - skipped, Wband for the skipped pulse. None if not operating in
            the skipped frame mode
    """
    sdd = source_detector_distance(ws)
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers

    def clip_band(band, lc, hc):
        tof_min = wlg.tof(band.min, sdd, ch.pulse_width) + lc
        tof_max = wlg.tof(band.max, sdd) - hc
        return wlg.Wband(wlg.from_tof(tof_min, sdd, ch.pulse_width),
                         wlg.from_tof(tof_max, sdd))

    bands = transmitted_bands(ws)
    lead_band = clip_band(bands.lead, ltc, htc)
    if bands.skip is None:
        skip_band = None
    else:
        if ltc2 is None:
            ltc2 = ltc
        if htc2 is None:
            htc2 = htc
        skip_band = clip_band(bands.skip, ltc2, htc2)

    return dict(lead=lead_band, skip=skip_band)


@namedtuplefy
def correct_frame(ws, ltc=0, htc=0, ltc2=None, htc2=None,
                  component='detector1'):
    r"""
    Assign the correct TOF to each event.

    TOFS may need to be shifted one (or more) frame widths and if working
    in frame skipped mode, one pulse period.

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    component: str
        Either 'detector1' or 'monitor1' to correct time of flights in the
        detector bank or in the monitor
    ltc: float
        trim neutrons of the leading pulse with a TOF smaller than the
        minimal TOF plus this value
    htc: float
        trim neutrons f the leading pulse with TOF bigger than the maxima
         TOF minus this value.
    ltc2: float
        as `ltc` but for neutrons of the skipped pulse. If None, then we use
        `ltc`
    htc2: float
        as `htc` but for neutrons of the skipped pulse. If None, then we use
        `htc`

    Returns
    -------
    namedtuple
        Fields of the namedtuple:
        - ws: EventsWorkspace, the input workspace
        - lead_band: WBand wavelength band for the lead pulse
        - skip_band: WBand wavelength band for the skip pulse. `None` if not
            working in frame-skipping mode
    """
    if ltc2 is None:
        ltc2 = ltc
    if htc2 is None:
        htc2 = htc

    sl = SampleLogs(ws)
    pulse_period = 1.e6 / sl.frequency.value.mean()  # 10^6/60 micro-seconds
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers

    # The TOF values recorded are never be bigger than the frame width,
    # which is also the choppers' rotational period.
    first_chopper = ch[0]
    frame_width = first_chopper.period

    instrument = ws.getInstrument()
    if component == 'detector1':  # Find distances from source to detectors
        sample = instrument.getSample()
        s2s = sample_source_distance(ws, units='m')
        s2c_i = [s2s + instrument.getDetector(i).getDistance(sample)
                 for i in range(ws.getNumberHistograms())]
        s2c = min(s2c_i)  # nominal distance from source to component
    elif component == 'monitor1':  # Find distance from source to monitor
        source = instrument.getSource()
        s2c = instrument.getDetector(-1).getDistance(source)
        s2c_i = [s2c, ]

    bands = transmitted_bands(ws)

    def tof_min_max(band, distance):
        return (wlg.tof(band.min, distance, ch.pulse_width),
                wlg.tof(band.max, s2c))

    # Find how many frame widths elapsed from the time the neutrons of the
    # lead pulse were emitted and the time the neutrons arrived to the
    # detector bank. This time must be added to the stored TOF values
    tof_min, tof_max = tof_min_max(bands.lead, s2c)
    frames_offset_time = frame_width * int((tof_min + tof_max)
                                           / (2 * frame_width))

    for i in range(ws.getNumberHistograms()):
        sp = ws.getSpectrum(i)
        tofs = sp.getTofs()
        if len(tofs) == 0:
            continue  # no events found in this detector
        tofs += frames_offset_time  # shift times to the correct frame
        d = s2c_i[i]  # distance from source to detector
        # TOF boundaries for the leading pulse
        tof_min, tof_max = tof_min_max(bands.lead, d)
        # TOF values smaller than that of the fastest neutrons have been
        # 'folded' by the data acquisition system. They must be shifted
        tofs[np.where(tofs < tof_min)] += frame_width
        if ch.frame_mode == FrameMode.skip:
            # Neutrons slower than the slowest neutron of the leading pulse
            # originated in the skipped pulse. Their TOFs must be
            # increased by a pulse period
            tofs[np.where(tofs > tof_max)] += pulse_period
        else:
            # Neutrons stored as slower than the slowest neutron are actually
            # fast neutrons whose TOFs have been shifted by a frame width
            tofs[np.where(tofs > tof_max)] -= frame_width
        # Discard events in the leading pulse outside the valid TOF range
        valid_range = np.where((tofs > tof_min + ltc) &
                               (tofs < tof_max - htc))[0]
        if ch.frame_mode == FrameMode.skip:
            # Discard events in the skipped pulse utside the valid TOF range
            tof_min, tof_max = tof_min_max(bands.skip, d)
            r12 = np.where((tofs > tof_min + ltc2) &
                           (tofs < tof_max - htc2))[0]
            valid_range = np.append(valid_range, r12)
        # reset the event list with only the valid TOF's
        valid_tofs = tofs[valid_range]
        valid_pts = np.asarray(sp.getPulseTimes())[valid_range]
        sp.clear(False)
        for tof, pt in zip(valid_tofs, valid_pts):
            sp.addEventQuickly(tof, pt)

    return dict(ws=ws, lead_band=bands.lead, skip_band=bands.skip)
