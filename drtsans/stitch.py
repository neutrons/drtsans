import os, numpy as np
from mantid.simpleapi import WorkspaceFactory, Stitch1D

def stitch(q1,y1,erry1,errq1, q2,y2,erry2,errq2, startoverlap, stopoverlap):
    newq = np.linspace(startoverlap, stopoverlap, 100)
    ws1 = _toWS(q1, y1, newq)
    ws2 = _toWS(q2, y2, newq)
    # stitch
    _, scale = Stitch1D(LHSWorkspace=ws1, RHSWorkspace=ws2, StartOverlap=startoverlap, EndOverlap=stopoverlap)
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
