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
        ax.plot(expected_q, expected_y, '+', label='expected')
        ax.plot(q1, y1, 'v', label='low q')
        ax.plot(q2, y2*scale, 'o', label='high q scaled')
        ax.plot(qout, yout, label='stitched')
        ax.legend()
        plt.savefig("stitch.pdf")
    # check
    assert np.allclose(qout, expected_q)
    assert np.allclose(yout, expected_y, rtol=2e-4)
    assert np.allclose(erryout, expected_err, rtol=2e-4)
    assert np.allclose(errqout, expected_errq, rtol=2e-4)
    return


if __name__ == '__main__':
    pytest.main()
 
