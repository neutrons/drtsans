# flake8: noqa
# https://docs.mantidproject.org/nightly/algorithms/Stitch1D-v3.html
import os, numpy as np
from mantid.simpleapi import WorkspaceFactory, Stitch1D

def stitch(q1,y1,erry1,errq1, q2,y2,erry2,errq2, startoverlap, stopoverlap):
    r"""
    The algorithm calculates the stitched spectrum out of two spectrum, one
    for low q, another for high q.

    **Mantid algorithms used:**
    :ref:`Stitch1D <algm-Stitch1D-v3>`,

    Parameters
    ----------
    q1,y1,erry1,errq1: array
        numpy arrays of left spectrum
    q2,y2,erry2,errq2: array
        numpy arrays of right spectrum
    startoverlap: float
        Starting point of overlap for q axis
    stopoverlap: float
        Stopping point of overlap for q axis

    Returns
    -------
    qout,yout,erryout,errqout, scale
    """
    # stitch using Mantid Stitch1D algorithm.
    # first interpolate the data in the overlap region
    newq = np.linspace(startoverlap, stopoverlap, 100)
    ws1 = _toWS(q1, y1, newq)
    ws2 = _toWS(q2, y2, newq)
    # then stitch
    _, scale = Stitch1D(LHSWorkspace=ws1, RHSWorkspace=ws2, StartOverlap=startoverlap, EndOverlap=stopoverlap)
    # put together the output data. it is the concatenation of low_q_spectrum[q<startoverlap] and high_q_spectrum[q>startoverlap]
    lowqrange = q1<startoverlap
    highqrange = q2>startoverlap
    qout = np.concatenate((q1[lowqrange], q2[highqrange]))
    yout = np.concatenate((y1[lowqrange], scale*y2[highqrange]))
    erryout = np.concatenate((erry1[lowqrange], scale*erry2[highqrange]))
    errqout = np.concatenate((errq1[lowqrange], errq2[highqrange]))
    return qout, yout, erryout, errqout, scale


def _toWS(x,y, newx):
    newy = np.interp(newx, x, y)
    nbins = newx.size
    ws = WorkspaceFactory.create(
        "Workspace2D", NVectors=1, XLength=nbins, YLength=nbins)
    ws.setX(0, newx)
    ws.setY(0, newy)
    return ws
