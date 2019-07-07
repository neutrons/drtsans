from __future__ import (absolute_import, division, print_function)

import numpy as np
from mantid.dataobjects import MaskWorkspace
from mantid.simpleapi import (LoadMask, MaskDetectors, MaskBTP, CloneWorkspace)
from ornl.settings import (optional_output_workspace,
                           unique_workspace_dundername as uwd)


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


@optional_output_workspace
def apply_mask(w, mask=None, **btp):
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
    instrument = w.getInstrument().name()
    if isinstance(mask, str):
        wm = LoadMask(Instrument=instrument, InputFile=mask, RefWorkspace=w,
                      OutputWorkspace=uwd())
    elif isinstance(mask, MaskWorkspace):
        wm = CloneWorkspace(mask, OutputWorkspace=uwd())
    else:
        raise RuntimeError('mask not understood')
    if bool(btp):
        wm2 = MaskBTP(instrument=instrument, **btp, OutputWorkspace=uwd())
        wm += wm2
    MaskDetectors(Workspace=w, MaskedWorkspace=wm)
    return wm
