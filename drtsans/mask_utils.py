import numpy as np
from mantid.api import mtd
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadMask, MaskDetectors, MaskBTP, ExtractMask, MaskSpectra)
from drtsans.settings import unique_workspace_dundername as uwd


def mask_as_numpy_array(w, invert=False):
    """Return mask in pixels as numpy array of bool. Items are True if masked
    :param w: input workspace
    :param invert: invert the array, items are True if unmasked
    :return: numpy.ndarray(bool)
    """
    mask = [w.getDetector(i).isMasked()
            for i in range(w.getNumberHistograms())]
    mask = np.asarray(mask)
    return mask if invert is False else np.invert(mask)


def masked_indexes(w, invert=False):
    """List of masked workspaces indexes
    :param w: input workspace
    :param invert: Return list of unmasked workspace indexes if True
    :return: numpy.ndarray(bool)
    """
    mask = mask_as_numpy_array(w, invert=invert)
    return np.where(mask)[0]


def apply_mask(w, mask=None, output_workspace=None, **btp):
    r"""
    Apply a mask to a workspace.

    The function accepts a path to a mask file or a MaskWorkspace,
    plus options for algorithm MaskBTP.

    Parameters
    ----------
    w: Workspace
        Workspace to be masked
    mask: mask file path, MaskWorkspace, list
        Mask to be applied. If `list`, it is a list of detector ID's. If
        `None`, it is expected that `maskbtp` is not empty.
    btp: dict
        Options to Mantid algorithm MaskBTP. Will be used if `mask=None`

    Returns
    -------
    MaskWorkspace
        Combination of mask and MaskBTP
    """
    w = str(w)
    if output_workspace is None:
        output_workspace = uwd()
    instrument = mtd[w].getInstrument().getName()
    if mask is not None:
        if isinstance(mask, str):
            wm = LoadMask(Instrument=instrument, InputFile=mask,
                          RefWorkspace=w, OutputWorkspace=uwd())
            MaskDetectors(Workspace=w, MaskedWorkspace=wm)
            wm.delete()  # delete temporary workspace
        elif isinstance(mask, MaskWorkspace):
            MaskDetectors(Workspace=w, MaskedWorkspace=mask)
        elif isinstance(mask, list):
            MaskDetectors(Workspace=w, DetectorList=mask)
    if bool(btp):
        MaskBTP(Workspace=w, **btp)
    return ExtractMask(InputWorkspace=w,
                       OutputWorkspace=output_workspace).OutputWorkspace


def mask_spectra_with_special_values(input_workspace, output_workspace=None):
    r"""
    Mask spectra in a workspace containing non-finite values.

    Non-finite values are evaluated with `numpy.isfinite`

    Parameters
    ----------
    input_workspace: str, MatrixWorkspace
    special_values: list
        List of string representations for special `float` values. The special value can be obtained by applying
        `float` to the string, e.g. float('nan').
    output_workspace : str
        Name of the normalised workspace. If None, the name of the input
        workspace is chosen (the input workspace is overwritten).

    Returns
    -------
    list
        Workspace indexes masked. Returns zero if no spectra are masked.
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    workspace = mtd[str(input_workspace)]
    intensities = workspace.extractY()
    non_finite_indexes = np.argwhere(np.isfinite(np.sum(intensities, axis=-1)) == False)  # noqa: E712
    non_finite_indexes = non_finite_indexes.flatten().tolist()
    if len(non_finite_indexes) > 0:
        MaskSpectra(InputWorkspace=input_workspace, InputWorkspaceIndexType='WorkspaceIndex',
                    InputWorkspaceIndexSet=non_finite_indexes, OutputWorkspace=output_workspace)
    return len(non_finite_indexes)
