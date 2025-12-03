from mantid.simpleapi import logger, mtd

r"""
Hyperlinks to drtsans functions
exists, registered_workspace <https://github.com/neutrons/drtsans/blob/next/src/drtsans/path.py>
"""  # noqa: E501
from drtsans import subtract_background
from drtsans.path import registered_workspace
from drtsans.tof.eqsans.dark_current import subtract_dark_current  # noqa E402
from drtsans.tof.eqsans.normalization import normalize_by_flux  # noqa E402

__all__ = ["subtract_blocked_beam"]


def subtract_blocked_beam(
    input_workspace, blocked_beam=None, flux_method=None, flux=None, dark_current=None, output_workspace=None
):
    r"""
    Subtracts a blocked beam background from the input workspace, with optional
    dark current subtraction and flux normalization.

    Parameters
    ----------
    input_workspace : Workspace or str
        The workspace from which the blocked beam background will be subtracted.
    blocked_beam : namedtuple
        (~mantid.dataobjects.Workspace2D, ~mantid.dataobjects.Workspace2D)
        An object containing the blocked beam data. If None or its data is None, subtraction is skipped.
    flux_method : str, optional
        The method used for flux normalization. If "monitor", blocked beam subtraction is skipped.
    flux : Workspace or str, optional
        The workspace or value used for flux normalization.
    dark_current : namedtuple
        (~mantid.dataobjects.Workspace2D, ~mantid.dataobjects.Workspace2D)
        An object containing dark current data. If provided, dark current is subtracted from the blocked beam.
    output_workspace : str, optional
        The name of the output workspace. If None, defaults to the name of the input workspace.
    Returns
    -------
    None
    """
    if blocked_beam is None or blocked_beam.data is None:
        return

    if flux_method == "monitor":
        logger.warning(
            "Blocked beam run was supplied but subtraction is not compatible with monitor flux normalization. "
            "Skipping blocked beam subtraction."
        )
        return

    if output_workspace is None:
        output_workspace = str(input_workspace)

    bb_ws_name = str(blocked_beam.data).replace("_raw_histo", "_processed_histo")
    if not registered_workspace(bb_ws_name):
        mtd[str(blocked_beam.data)].clone(OutputWorkspace=bb_ws_name)

        if dark_current is not None and dark_current.data is not None:
            subtract_dark_current(bb_ws_name, dark_current.data)

        normalize_by_flux(bb_ws_name, flux, method=flux_method)

    subtract_background(input_workspace, mtd[bb_ws_name], output_workspace=output_workspace)
