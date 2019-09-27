import os, numpy as np
import pytest
from mantid.simpleapi import Load, CreateWorkspace, WorkspaceFactory, Stitch1D
from drtsans.stitch import stitch

headless = 'DISPLAY' not in os.environ or not os.environ['DISPLAY']


def test_stitch(reference_dir):
    # function to load test data
    def loadData(filename):
        datadir = os.path.join(reference_dir.new.eqsans, "test_stitch")
        data = os.path.join(datadir, filename)
        x,y,erry,errx = np.loadtxt(data, skiprows=2).T
        return x,y,erry,errx
    # load test data
    (x1, y1, erry1, errx1) = loadData("sMCM_cc_4m_Iq.txt")
    (x2, y2, erry2, errx2) = loadData("sMCM_cc_hq_Iq.txt")
    # parameters for stitching
    startoverlap, stopoverlap = 0.04, 0.08
    # stitch
    xout, yout, erryout, errxout, scale = stitch(x1,y1,erry1,errx1, x2,y2,erry2,errx2, startoverlap, stopoverlap)
    if not headless:
        # plot
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots()
        expected_x,expected_y,expected_err,expected_errx = loadData("sMCM_cc_stitched.txt")
        ax.plot(expected_x, expected_y, '+', label='expected')
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
 
