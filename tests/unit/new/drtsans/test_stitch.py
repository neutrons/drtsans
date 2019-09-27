import os, numpy as np
import pytest
from mantid.simpleapi import Load, CreateWorkspace, WorkspaceFactory, Stitch1D


def test_stitch(reference_dir):
    def loadData(filename):
        datadir = os.path.join(reference_dir.new.eqsans, "test_stitch")
        data = os.path.join(datadir, filename)
        x,y,erry,errx = np.loadtxt(data, skiprows=2).T
        return x,y,erry,errx
    def toWS(x,y, newx):
        newy = np.interp(newx, x, y)
        nbins = newx.size
        ws = WorkspaceFactory.create(
            "Workspace2D", NVectors=1, XLength=nbins, YLength=nbins)
        ws.setX(0, newx)
        ws.setY(0, newy)
        return ws
    # load data
    (x1, y1, erry1, errx1) = loadData("sMCM_cc_4m_Iq.txt")
    (x2, y2, erry2, errx2) = loadData("sMCM_cc_hq_Iq.txt")
    startoverlap, stopoverlap = 0.04, 0.08
    newx = np.linspace(startoverlap, stopoverlap, 100)
    ws1 = toWS(x1, y1, newx)
    ws2 = toWS(x2, y2, newx)
    # stitch
    _, scale = Stitch1D(LHSWorkspace=ws1, RHSWorkspace=ws2, StartOverlap=startoverlap, EndOverlap=stopoverlap)
    lowqrange = x1<startoverlap
    highqrange = x2>startoverlap
    xout = np.concatenate((x1[lowqrange], x2[highqrange]))
    yout = np.concatenate((y1[lowqrange], scale*y2[highqrange]))
    # plot
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    expected_x,expected_y,expected_err,expected_errx = loadData("sMCM_cc_stitched.txt")
    ax.plot(expected_x, expected_y, '+', label='expected')
    # ax.errorbar(stitched)
    ax.plot(x1, y1, 'v', label='low q')
    ax.plot(x2, y2*scale, 'o', label='high q scaled')
    ax.plot(xout, yout, label='stitched')
    ax.legend()
    plt.savefig("stitch.pdf")
    # check
    assert np.allclose(xout, expected_x)
    assert np.allclose(yout, expected_y, rtol=2e-4)
    return


if __name__ == '__main__':
    pytest.main()
 
