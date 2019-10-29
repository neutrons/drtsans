# https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py
from drtsans.settings import unique_workspace_dundername as uwd
# https://docs.mantidproject.org/nightly/algorithms/CloneWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/Minus-v1.html
# https://docs.mantidproject.org/nightly/algorithms/RebinToWorkspace-v1.html
# https://docs.mantidproject.org/nightly/algorithms/RenameWorkspace-v1.html
# mtd is https://docs.mantidproject.org/nightly/api/python/mantid/api/AnalysisDataServiceImpl.html
from mantid.simpleapi import CloneWorkspace, CreateSingleValuedWorkspace, Minus, mtd, RebinToWorkspace, RenameWorkspace
import numpy as np

__all__ = ['half_polarization', 'subtract_background']


def subtract_background(input_workspace, background, scale=1.0, scale_error=0.0,
                        output_workspace=None):
    r"""
    Subtract a prepared background from a prepared sample.

    Perform a rebin if sample and background have different binning.

    **Mantid algorithms used:**
    :ref:`CloneWorkspace <algm-CloneWorkspace-v1>`
    :ref:`CreateSingleValuedWorkspace <algm-CreateSingleValuedWorkspace-v1>`
    :ref:`Minus <algm-Minus-v1>`
    :ref:`Multiply <algm-Multiply-v1>`
    :ref:`RebinToWorkspace <algm-RebinToWorkspace-v1>`

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        Sample workspace.
    background: str, ~mantid.api.MatrixWorkspace
        Background workspace.
    scale: float
        Rescale background intensities by this multiplicative factor before
        subtraction from the sample.
    scale_error: float
        Uncertainty in scale factor
    output_workspace: str
        Name of the sample corrected by the background. If :py:obj:`None`, then
        ``input_workspace`` will be overwritten.

    Returns
    -------
    ~mantid.api.MatrixWorkspace
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)

    # get handle on input_workspace
    input_workspace = mtd[str(input_workspace)]

    workspaces_to_delete = []

    # make the background match the input_workspace binning if possible
    if input_workspace.isHistogramData() and background.isHistogramData():
        # rebin the background to match the input_workspace
        # this will do nothing if the x-axis is already the same
        # put it in the output_workspace to save space in memory
        background_rebinned = RebinToWorkspace(WorkspaceToRebin=background,
                                               WorkspaceToMatch=input_workspace,
                                               OutputWorkspace=uwd())
    else:
        # subtract will check that the x-axis is identical
        # otherwise the subtraction will fail
        background_rebinned = CloneWorkspace(InputWorkspace=background, OutputWorkspace=uwd())
    workspaces_to_delete.append(str(background_rebinned))

    # need to create a special object if the uncertainty in the scale was specified
    if scale_error != 0.:
        scale = CreateSingleValuedWorkspace(DataValue=scale, ErrorValue=scale_error,
                                            OutputWorkspace=uwd())
        workspaces_to_delete.append(str(scale))

    background_rebinned *= scale

    Minus(LHSWorkspace=input_workspace,
          RHSWorkspace=background_rebinned,
          OutputWorkspace=output_workspace)

    for name in workspaces_to_delete:
        if mtd.doesExist(name):
            mtd.remove(name)

    return mtd[output_workspace]


def _calc_flipping_ratio(polarization):
    '''Calculates the flipping ratio from the polarization state

    Parameters
    ----------
    polarization: str, ~mantid.api.MatrixWorkspace
        Polarization state

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The ratio of flipping
    '''
    value = polarization.extractY()[0]
    uncertainty = polarization.extractE()[0]
    if len(value) == 1 and len(uncertainty) == 1:
        value, uncertainty = value[0], uncertainty[0]

    uncertainty = 2. * uncertainty / np.square(1 - value)
    value = (1 + value) / (1 - value)

    if isinstance(value, float):  # create a single value workspace
        return CreateSingleValuedWorkspace(DataValue=value, ErrorValue=uncertainty,
                                           OutputWorkspace=uwd(), EnableLogging=False)
    else:
        # if it is an array should call CreateWorkspace(EnableLogging=False)
        raise NotImplementedError('Somebody needs to create an output from {} (type={})'.format(value, type(value)))


def _set_uncertainty_from_numpy(wksp, uncertainty):
    '''Inner function to override the uncertainties on a workspace. This is a detail that
    is easy to mess up'''
    # TODO add support for more workspace types / dimensions
    # this works well for single value workspaces
    for i in range(uncertainty.shape[0]):
        wksp.setE(i, np.array([uncertainty[i]]))

    return wksp


def _calc_half_polarization_up(flipper_off, flipper_on, efficiency, flipping_ratio):
    '''This calculates the spin up workspace

    Parameters
    ----------
    flipper_off_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper off measurement
    flipper_on_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper on measurement
    efficiency: ~mantid.api.MatrixWorkspace
        Flipper efficiency
    flipping_ratio: ~mantid.api.MatrixWorkspace
        The ratio of flipping

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The spin up workspace
    '''
    __spin_up = flipper_off + (flipper_off - flipper_on) / (efficiency * (flipping_ratio - 1.))

    e = efficiency.extractY()[0]
    F = flipping_ratio.extractY()[0]

    # the uncertainties (numerically) aren't correct because of build-up of numerical errors.
    # Recalculate them based on the proper equation based on a hand calculation of the partial derivatives
    m0_part = np.square(flipper_off.extractE()[0][0] * (1 + 1/(e*(F-1))))
    m1_part = np.square(flipper_on.extractE()[0][0] * (1 / (e*(F-1))))

    mixed = np.square((flipper_off.extractY()[0][0] - flipper_on.extractY()[0][0])/(e * (F-1)))
    mixed *= (np.square(efficiency.extractE()[0][0] / e)
              + np.square(flipping_ratio.extractE()[0][0] / (F-1)))
    sup_err = np.sqrt(m0_part + m1_part + mixed)

    # set the uncertainty in the workspace
    __spin_up = _set_uncertainty_from_numpy(__spin_up, sup_err)

    return __spin_up


def _calc_half_polarization_down(flipper_off, flipper_on, efficiency, flipping_ratio):
    '''This calculates the spin down workspace

    Parameters
    ----------
    flipper_off_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper off measurement
    flipper_on_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper on measurement
    efficiency: ~mantid.api.MatrixWorkspace
        Flipper efficiency
    flipping_ratio: ~mantid.api.MatrixWorkspace
        The ratio of flipping

    Returns
    -------
    ~mantid.api.MatrixWorkspace
        The spin down workspace
    '''
    __spin_down = flipper_off - (flipper_off - flipper_on) / (efficiency * (1. - 1./flipping_ratio))

    e = efficiency.extractY()[0]
    F = flipping_ratio.extractY()[0]

    # the uncertainties (numerically) aren't correct because of build-up of numerical errors.
    # Recalculate them based on the proper equation based on a hand calculation of the partial derivatives
    m0_part = np.square(flipper_off.extractE()[0][0] * (1 - 1/(e*(1 - 1/F))))
    m1_part = np.square(flipper_on.extractE()[0][0] * (1 / (e*(1-1/F))))
    mixed = np.square((flipper_off.extractY()[0][0] - flipper_on.extractY()[0][0])/(e * (1-1/F)))
    mixed *= (np.square(efficiency.extractE()[0][0] / e)
              + np.square(flipping_ratio.extractE()[0][0] / (F*F*(1-1/F))))

    sdn_err = np.sqrt(m0_part + m1_part + mixed)

    # set the uncertainty in the workspace
    __spin_down = _set_uncertainty_from_numpy(__spin_down, sdn_err)

    return __spin_down


def half_polarization(flipper_off_workspace, flipper_on_workspace, polarization, efficiency,
                      spin_up_workspace=None, spin_down_workspace=None):
    '''Calculate the spin up/down workspaces from flipper on/off.

    **Mantid algorithms used:**
    :ref:`RenameWorkspace <algm-RenameWorkspace-v1>`

    Parameters
    ----------
    flipper_off_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper off measurement
    flipper_on_workspace: str, ~mantid.api.MatrixWorkspace
        Flipper on measurement
    polarization: str, ~mantid.api.MatrixWorkspace
        Polarization state
    efficiency: str, ~mantid.api.MatrixWorkspace
        Flipper efficiency
    spin_up_workspace: str
        Name of the resulting spin up workspace. If :py:obj:`None`, then
        ``flipper_off_workspace`` will be overwritten.
    spin_down_workspace: str
        Name of the resulting spin down workspace. If :py:obj:`None`, then
        ``flipper_on_workspace`` will be overwritten.

    Returns
    -------
    py:obj:`tuple` of 2 ~mantid.api.MatrixWorkspace
    '''
    if spin_up_workspace is None:
        spin_up_workspace = str(flipper_off_workspace)
    if spin_down_workspace is None:
        spin_down_workspace = str(flipper_on_workspace)

    # this is denoted as "F" in the master document
    flipping_ratio = _calc_flipping_ratio(polarization)

    __spin_up = _calc_half_polarization_up(flipper_off_workspace, flipper_on_workspace, efficiency, flipping_ratio)
    __spin_down = _calc_half_polarization_down(flipper_off_workspace, flipper_on_workspace, efficiency, flipping_ratio)

    spin_up_workspace = RenameWorkspace(InputWorkspace=__spin_up, OutputWorkspace=spin_up_workspace)
    spin_down_workspace = RenameWorkspace(InputWorkspace=__spin_down, OutputWorkspace=spin_down_workspace)

    return (spin_up_workspace, spin_down_workspace)
