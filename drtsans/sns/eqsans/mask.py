from mantid.simpleapi import MaskBTP
from ornl.sans import mask_utils

__all__ = ['apply_mask']


def apply_mask(w, panel=None, mask=None, output_workspace=None, **btp):
    r"""
    Apply a mask or region-of-interest to a workspace.

    The function accepts a path to a mask file, a MaskWorkspace, or options
    to algorithm :ref:`MaskBTP <algm-MaskBTP-v1>`.

    Parameters
    ----------
    w: ~mantid.api.IEventWorkspace
        Workspace to be masked
    panel: str
        Either 'front' or 'back' to mask a whole panel
    mask: mask file path, ``MaskWorkspace``, :py:obj:`list`
        Additional mask to be applied. If :py:obj:`list`, it is a list of
        detector ID's.
    btp: dict
        Options to Mantid algorithm :ref:`MaskBTP <algm-MaskBTP-v1>`. Will be used if ``mask=None``

    Returns
    -------
    MaskWorkspace
        Combination of panel, mask, and :ref:`MaskBTP <algm-MaskBTP-v1>` masks
    """
    if panel:
        MaskBTP(Workspace=w, instrument='EQ-SANS', Components=panel + '-panel')
    # let apply_mask handle selecting output workspace's name
    return mask_utils.apply_mask(w, mask=mask,
                                 output_workspace=output_workspace, **btp)
