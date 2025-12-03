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
    input_workspace=None, blocked_beam=None, flux_method=None, flux=None, dark_current=None, output_workspace=None
):
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
