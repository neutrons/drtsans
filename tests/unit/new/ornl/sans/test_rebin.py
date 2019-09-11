from __future__ import (absolute_import, division, print_function)

import pytest


@pytest.fixture
def fake_events():
    def _(wavelengths):
        from mantid.simpleapi import CreateSampleWorkspace
        import mantid.kernel
        ws = CreateSampleWorkspace(
            WorkspaceType='Event',
            Function='Flat background',
            NumBanks=1,
            BankPixelWidth=1,
            NumEvents=0,
            XUnit='Wavelength',
            XMax=10,
            BinWidth=1)
        sp0 = ws.getSpectrum(0)
        date = mantid.kernel.DateAndTime("2019-09-09T00:00")
        for lam in wavelengths:
            sp0.addEventQuickly(lam, date)
        return ws
    return _


def test_linear(fake_events):
    """Test linear binning by creating a workspace with fake events at specific wavelengths,
    binning into a histogram, and test against expected output.

    For details see https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/208
    """
    # create a workspace with fake events
    ws = fake_events([2.57, 3.05, 2.76, 3.13, 2.84])
    # binning the events to linear bins 2.5, 2.6, ..., 3.2
    from mantid.simpleapi import Rebin
    ws = Rebin(InputWorkspace=ws, Params='2.5, .1, 3.2')  # start, step, end of bin edges
    # verify
    import numpy as np
    assert np.allclose(ws.readX(0), np.arange(2.5, 3.21, .1))
    assert np.allclose(ws.readY(0), [1., 0., 1., 1., 0., 1., 1.])
    return


if __name__ == '__main__':
    pytest.main()
