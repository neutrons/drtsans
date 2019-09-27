# flake8: noqa
import os
import numpy as np
import pytest
from drtsans.stitch import stitch

headless = 'DISPLAY' not in os.environ or not os.environ['DISPLAY']


def test_stitch(reference_dir):
    # function to load test data
    def loadData(filename):
        datadir = os.path.join(reference_dir.new.eqsans, "test_stitch")
        data = os.path.join(datadir, filename)
        q,y,erry,errq = np.loadtxt(data, skiprows=2).T
        return q,y,erry,errq
    # load test data
    (q1, y1, erry1, errq1) = loadData("sMCM_cc_4m_Iq.txt")
    (q2, y2, erry2, errq2) = loadData("sMCM_cc_hq_Iq.txt")
    # parameters for stitching
    startoverlap, stopoverlap = 0.04, 0.08
    # stitch
    qout, yout, erryout, errqout, scale = stitch(q1,y1,erry1,errq1, q2,y2,erry2,errq2, startoverlap, stopoverlap)
    # expected results
    expected_q,expected_y,expected_err,expected_errq = loadData("sMCM_cc_stitched.txt")
    if not headless:
        # plot
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots()
        ax.loglog(expected_q, expected_y, '+', label='expected')
        ax.loglog(q1, y1, 'v', mfc='none', label='low q')
        ax.loglog(q2, y2*scale, 'o', mfc='none', label='high q scaled')
        ax.loglog(qout, yout, label='stitched')
        ax.plot([startoverlap, startoverlap], [0, np.max(yout)*10], 'k-')
        ax.plot([stopoverlap, stopoverlap], [0, np.max(yout)*10], 'k-')
        ax.set_ylim(0, np.max(yout)*1.1)
        ax.legend()
        plt.savefig("stitch.png")
    # check
    np.testing.assert_almost_equal(qout, expected_q)
    np.testing.assert_almost_equal(yout, expected_y, decimal=4)
    np.testing.assert_almost_equal(erryout, expected_err, decimal=4)
    np.testing.assert_almost_equal(errqout, expected_errq, decimal=7)
    return


if __name__ == '__main__':
    pytest.main()
 
