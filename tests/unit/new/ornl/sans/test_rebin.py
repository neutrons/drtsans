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
    ws = fake_events([1., 2.])
    from mantid.simpleapi import Rebin
    ws = Rebin(InputWorkspace=ws, Params='.1')
    return


if __name__ == '__main__':
    pytest.main()
