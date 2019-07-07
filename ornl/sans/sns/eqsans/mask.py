from ornl.sans import mask_utils

__all__ = ['apply_mask']


def apply_mask(w, panel=None, mask=None, **btp):
    r"""
    Apply a mask or region-of-interest to a workspace.

    The function accepts a path to a mask file, a MaskWorkspace, or options
    to algorithm MaskBTP.

    Parameters
    ----------
    w: Workspace
        Workspace to be masked
    panel: str
        Either 'front' or 'back' to mask a whole panel
    mask: mask file path, MaskWorkspace
        Mask to be applied
    btp: dict
        Options to Mantid algorithm MaskBTP. Will be used if `mask=None`
    """
    comp = dict(front='')
    if panel:
        btp.update(Components=panel + '-panel')
        mask_utils.apply_mask(w, **btp)
    else:
        apply_mask(w, mask=mask, **btp)
