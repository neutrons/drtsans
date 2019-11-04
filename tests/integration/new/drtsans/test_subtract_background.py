from drtsans import subtract_background
# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py
from drtsans.settings import unique_workspace_dundername as uwd
# https://docs.mantidproject.org/nightly/algorithms/CompareWorkspaces-v1.html
# https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
from mantid.simpleapi import CompareWorkspaces, CreateWorkspace, DeleteWorkspace
import numpy as np
import pytest


class _Data1D(object):
    '''This is a factory class for generating the 1d example data

    The equations are taken from the spreadsheet supplied for 1d data
    '''
    scale_factor = 0.92
    Sig_scale = 0.005 * scale_factor

    Q_Scale = np.linspace(.01, .99, 99)  # 99 numbers from 0.01 to 0.99
    I_Background_1d = np.power(Q_Scale * 10, -4) + .57  # power-law maximum (Porod-scattering)
    Sig_background_1d = I_Background_1d * 0.02

    I_Data_1d = scale_factor * I_Background_1d + 0.2 * (0.5 - np.absolute(0.5 - Q_Scale))
    Sig_data_1d = 0.01 * I_Data_1d

    I_output_1d = I_Data_1d - scale_factor * I_Background_1d
    Sig_output_1d = np.sqrt(np.power(Sig_data_1d, 2) + np.power(scale_factor * Sig_background_1d, 2)
                            + np.power(Sig_scale * I_Background_1d, 2))

    def create(self, datatype, y=None, e=None):
        '''This function creates data based on supplied information. There are pre-defined ``datatype``
        of ``data`` and ``background``. All others are custom.
        '''
        if datatype == 'data':
            y = self.I_Data_1d
            e = self.Sig_data_1d
        elif datatype == 'background':
            y = self.I_Background_1d
            e = self.Sig_background_1d
        elif datatype == 'output':
            y = self.I_output_1d
            e = self.Sig_output_1d
        else:
            raise RuntimeError('Unknown data type="{}"'.format(datatype))

        # create a workspace with the correct signal and uncertainties and random name
        return CreateWorkspace(DataX=self.Q_Scale, DataY=y, DataE=e,
                               UnitX='momentumtransfer', OutputWorkspace=uwd())


# -------------------- 1d tests
def test_data_not_background_1d():
    '''This tests that the ``data`` is not equal to the ``background``

    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Ken Littrell <littrellkc@ornl.gov>
    '''
    factory = _Data1D()

    data = factory.create('data')
    background = factory.create('background')

    assert not CompareWorkspaces(data, background).Result

    DeleteWorkspace(data)
    DeleteWorkspace(background)


def test_subtract_background_1d():
    '''This tests that ``data - scale * background`` and its uncertainties gives the expected result.

    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Ken Littrell <littrellkc@ornl.gov>
    '''
    factory = _Data1D()

    # create workspaces with the input data
    data = factory.create('data')
    background = factory.create('background')
    expected = factory.create('output')

    # do the calculation using the framework in-place
    observed = subtract_background(data, background, scale=factory.scale_factor, scale_error=factory.Sig_scale)

    # check the results
    np.testing.assert_equal(observed.extractX(), expected.extractX())
    np.testing.assert_equal(observed.extractY(), expected.extractY())
    np.testing.assert_almost_equal(observed.extractE(), expected.extractE())  # sqrts aren't quite the same

    # cleanup workspaces that were created
    for wksp in [data, background, expected]:
        DeleteWorkspace(wksp)


# -------------------- 2d tests
def test_subtract_background_2d_linearized():
    pass


def test_subtract_background_2d():
    pass


if __name__ == '__main__':
    pytest.main([__file__])
