from __future__ import (absolute_import, division, print_function)

import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans.chopper import EQSANSDiskChopperSet
from ornl.sans.frame_mode import FrameMode
from ornl.sans import wavelength as wlg
from ornl.sans.core_utils import namedtuplefy


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
def correct_frame(ws, ltc=0, htc=0, ltc2=None, htc2=None):
    r"""
    Assign the correct TOF to each event.

    TOFS may need to be shifted one (or more) frame widths and if working
    in frame skipped mode, one pulse period.

    Parameters
    ----------
    ws: EventsWorkspace
        Data workspace
    ltc: float
        trim neutrons with TOF smaller than this value
    htc: float
        trim neutrons with TOF bigger than the pulse period minus this value.
    ltc2: float
        trim neutrons in the skipped pulse with TOF smaller than this.
        If None, then we use ltc
    htc2: float
        trim neutrons in the skipped pulse with TOF bigger than this.
        If None, then we use htc.


    Returns
    -------
    dict:
    """
    if ltc2 is None:
        ltc2 = ltc
    if htc2 is None:
        htc2 = htc
    sl = SampleLogs(ws)
    pulse_period = 1.e6 / sl.frequency.value.mean()  # 10^6/60 micro-seconds
    ch = EQSANSDiskChopperSet(ws)  # object representing the four choppers
    # Wavelength band of neutrons from the leading pulse transmitted
    # by the chopper system
    lead_band = ch.transmission_bands(pulsed=True)[0]
    # Wavelength from the previous, skipped pulse.
    skip_band = ch.transmission_bands(delay=pulse_period, pulsed=True)[0]\
        if ch.frame_mode == FrameMode.skip else None

    # The TOF values recorded are never be bigger than the frame width,
    # which is also the choppers' rotational period.
    first_chopper = ch[0]
    frame_width = first_chopper.period

    # Find distances from source to detectors
    instrument = ws.getInstrument()
    sample = instrument.getSample()
    source = instrument.getSource()
    s2s = sample.getDistance(source)  # distance from source to sample
    s2d_i = [s2s + instrument.getDetector(i).getDistance(sample)
             for i in range(ws.getNumberHistograms())]  # source to detectors
    s2d = min(s2d_i)  # nominal distance from source to detector bank

    # Find how many frame widths elapsed from the time the neutrons of the
    # lead pulse were emitted and the time the neutrons arrived to the
    # detector bank. This time must be added to the stored TOF values
    tof_min = wlg.tof(lead_band.min, s2d, ch.pulse_width)
    tof_max = wlg.tof(lead_band.max, s2d)
    frames_offset_time = frame_width * int((tof_min + tof_max) / (2 * frame_width))

    # Correct the TOF for each event
    for i in range(ws.getNumberHistograms()):
        sp = ws.getSpectrum(i)
        tofs = sp.getTofs()
        if len(tofs) == 0:
            continue  # no events found in this detector
        tofs += frames_offset_time  # shift times to the correct frame
        d = s2d_i[i]  # distance from source to detector
        # TOF boundaries for the leading pulse
        tof_min = wlg.tof(lead_band.min, d, ch.pulse_width)
        tof_max = wlg.tof(lead_band.max, d)
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
            # Discard events in the skipped pulse within the tof cutoffs
            tof_min = wlg.tof(skip_band.min, d, ch.pulse_width)
            tof_max = wlg.tof(skip_band.max, d)
            r12 = np.where((tofs > tof_min + ltc2) &
                           (tofs < tof_max - htc2))[0]
            valid_range = np.append(valid_range, r12)
        valid_tofs = tofs[valid_range]
        valid_pts = np.asarray(sp.getPulseTimes())[valid_range]
        # reset the event list with only the valid TOF's
        sp.clear(False)
        for tof, pt in zip(valid_tofs, valid_pts):
            sp.addEventQuickly(tof, pt)
    return dict(ws=ws, lead_band=lead_band, skip_band=skip_band)