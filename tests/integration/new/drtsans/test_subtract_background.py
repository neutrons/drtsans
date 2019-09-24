# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py
from drtsans.settings import unique_workspace_dundername as uwd
# https://docs.mantidproject.org/nightly/algorithms/CompareWorkspaces-v1.html
# https://docs.mantidproject.org/nightly/algorithms/CreateWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html
from mantid.simpleapi import CompareWorkspaces, CreateWorkspace, DeleteWorkspace
import numpy as np
import pytest

# these equations are taken from the spreadsheet supplied
scale_factor = 0.92
Sig_scale = 0.005 * scale_factor

Q_Scale = np.linspace(.01, .99, 99)  # 99 numbers from 0.01 to 0.99
I_Background = np.power(Q_Scale * 10, -4) + .57  # power-law maximum (Porod-scattering)
Sig_background = I_Background * 0.02

I_Data = scale_factor * I_Background + 0.2 * (0.5 - np.absolute(0.5 - Q_Scale))
Sig_data = 0.01 * I_Data

I_output = I_Data - scale_factor * I_Background
Sig_output = np.sqrt(np.power(Sig_data, 2) + np.power(scale_factor * Sig_background, 2)
                     + np.power(Sig_scale * I_Background, 2))

Q_dummy = Q_Scale.copy()
I_dummy = I_Background.copy()
Sig_Dummy1 = Sig_background.copy()


def create_workspace(datatype, y=None, e=None, x=Q_Scale):
    '''This function creates data based on supplied information. There are pre-defined ``datatype``
    of ``data`` and ``background``. All others are custom.
    '''
    if datatype == 'data':
        y = I_Data
        e = Sig_data
    elif datatype == 'background':
        y = I_Background
        e = Sig_background
    elif datatype == 'custom':
        if y is None:
            raise RuntimeError('Must supply signal')
        if e is None:
            e = np.sqrt(y)
    else:
        raise RuntimeError('Unknown data type={}'.format(datatype))

    # create a workspace with the correct signal and uncertainties and random name
    return CreateWorkspace(DataX=x, DataY=y, DataE=e,
                           UnitX='momentumtransfer', OutputWorkspace=uwd())


def test_data_not_background():
    '''This tests that the ``data`` is not equal to the ``background``

    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Ken Littrell <littrellkc@ornl.gov>
    '''
    data = create_workspace('data')
    background = create_workspace('background')

    assert not CompareWorkspaces(data, background).Result

    DeleteWorkspace(data)
    DeleteWorkspace(background)


def test_subtract_background():
    '''This tests that ``data - scale * background`` and its uncertainties gives the expected result.

    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Ken Littrell <littrellkc@ornl.gov>
    '''
    # create workspaces with the input data
    data = create_workspace('data')
    background = create_workspace('background')
    scale = create_workspace('custom', scale_factor, Sig_scale, x=42.)  # only workspaces can carry errors
    expected = create_workspace('custom', I_output, Sig_output)

    # do the calculation using the framework
    observed = data - scale * background

    # check the results
    np.testing.assert_equal(observed.extractY(), expected.extractY())
    np.testing.assert_almost_equal(observed.extractE(), expected.extractE())  # sqrts aren't quite the same

    # cleanup workspaces that were created
    DeleteWorkspace(data)
    DeleteWorkspace(background)
    DeleteWorkspace(scale)
    DeleteWorkspace(expected)
    DeleteWorkspace(observed)


if __name__ == '__main__':
    pytest.main()
