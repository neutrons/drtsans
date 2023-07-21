from drtsans.instruments import fetch_idf, empty_instrument_workspace
from drtsans.mono.biosans.geometry import has_midrange_detector
from drtsans.geometry import spectrum_info_ranges

# third party imports
from mantid.api import mtd
from mantid.dataobjects import Workspace2D
from mantid.simpleapi import CopyLogs

# standard imports
from typing import Union


def clone_component_intensities(
    input_workspace: Union[str, Workspace2D],
    output_workspace: str = None,
    input_component: str = "wing_detector",
    output_component: str = "midrange_detector",
    copy_logs: bool = True,
) -> Workspace2D:
    r"""
    Clone intensities into the detector pixels of a component from another component.

    Parameters
    ----------
    input_workspace
        The workspace to copy intensities from.
    output_workspace
        If :py:obj:`None`, the name of the input_workspace will be used, thus overwriting the input workspace.
    input_component
        The component to copy intensities from.
    output_component
        The component to copy intensities to.
    copy_logs
        If True, the sample logs will be copied to the output workspace

    Returns
    -------
    The workspace with intensities inserted.

    Raises
    ------
    ValueError
        If ``input_component=="midrange_detector" but ``input_workspace`` lacks the midrange detector.
    ValueError
        If ``output_component=="midrange_detector" but ``output_workspace`` lacks the midrange detector.
    """
    if input_component == "midrange_detector" and not has_midrange_detector(input_workspace):
        raise ValueError("No midrange detector in input workspace")

    if output_workspace is None:
        output_workspace = str(input_workspace)
    else:
        empty_instrument_workspace(
            output_workspace=output_workspace,
            filename=fetch_idf("BIOSANS_Definition.xml"),
            event_workspace=False,
            monitors_have_spectra=(mtd[str(input_workspace)].getSpectrum(0).getDetectorIDs()[0] < 0),
        )
        if copy_logs:
            CopyLogs(input_workspace, output_workspace)

    if output_component == "midrange_detector" and not has_midrange_detector(output_workspace):
        raise ValueError("No midrange detector in output workspace")

    in_first, in_last = spectrum_info_ranges(input_workspace, input_component)
    out_first, out_last = spectrum_info_ranges(output_workspace, output_component)
    assert in_last - in_first >= out_last - out_first, "Not enough pixels to clone"

    readers = [getattr(mtd[str(input_workspace)], reader) for reader in ("readX", "readY", "readE")]
    writers = [getattr(mtd[str(output_workspace)], writer) for writer in ("dataX", "dataY", "dataE")]
    for in_wi, out_wi in zip(range(in_first, in_last), range(out_first, out_last)):
        for reader, writer in zip(readers, writers):
            writer(out_wi)[:] = reader(in_wi)

    return mtd[str(output_workspace)]
