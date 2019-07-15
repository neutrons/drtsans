from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadMask, MaskDetectors, MaskBTP, ExtractMask)
from ornl.settings import unique_workspace_dundername as uwd


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
    mask: mask file path, MaskWorkspace
        Mask to be applied. If `None`, it is expected that `maskbtp`
        is not empty
    btp: dict
        Options to Mantid algorithm MaskBTP. Will be used if `mask=None`

    Returns
    -------
    MaskWorkspace
        Combination of mask and MaskBTP
    """
    if output_workspace is None:
        output_workspace = uwd()
    instrument = w.getInstrument().getName()
    w = str(w)
    if mask is not None:
        if isinstance(mask, str):
            wm = LoadMask(Instrument=instrument, InputFile=mask,
                          RefWorkspace=w, OutputWorkspace=uwd())
            MaskDetectors(Workspace=w, MaskedWorkspace=wm)
            wm.delete()  # delete temporary workspace
        elif isinstance(mask, MaskWorkspace):
            MaskDetectors(Workspace=w, MaskedWorkspace=mask)
    if bool(btp):
        MaskBTP(Workspace=w, **btp)
    return ExtractMask(InputWorkspace=w,
                       OutputWorkspace=output_workspace).OutputWorkspace
