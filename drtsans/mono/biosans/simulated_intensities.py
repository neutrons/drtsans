from drtsans.instruments import fetch_idf, empty_instrument_workspace
from drtsans.mono.biosans.geometry import has_midrange_detector
from drtsans.geometry import spectrum_info_ranges
from drtsans.path import registered_workspace

# third party imports
from mantid.api import mtd
from mantid.dataobjects import Workspace2D
from mantid.simpleapi import CopyLogs, LoadInstrument

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
    elif registered_workspace(output_workspace):
        pass
    else:
        empty_instrument_workspace(
            output_workspace=output_workspace,
            filename=fetch_idf("BIOSANS_Definition.xml"),
            event_workspace=False,
            monitors_have_spectra=(mtd[str(input_workspace)].getSpectrum(0).getDetectorIDs()[0] < 0),
        )
        # copy the units from the input workspace
        origin, target = mtd[str(input_workspace)], mtd[str(output_workspace)]
        origin_unit = origin.getAxis(0).getUnit().unitID()
        target.getAxis(0).setUnit(origin_unit)
        target.setYUnit(origin.YUnit())
        if copy_logs:
            CopyLogs(input_workspace, output_workspace)
        # This load forces the panels (detector1, wing_detector, midrange_detector) to moved according to the logs
        LoadInstrument(
            Workspace=output_workspace, Filename=fetch_idf("BIOSANS_Definition.xml"), RewriteSpectraMap=False
        )

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


def insert_midrange_detector(
    input_workspace: Union[str, Workspace2D],
    output_workspace: str = None,
) -> Workspace2D:
    r"""
    Insert the midrange detector into a Workspace2D lacking the midrange detector.

    Parameters
    ----------
    input_workspace
        The workspace to copy intensities from.
    output_workspace
        If :py:obj:`None`, the name of the input_workspace will be used, thus overwriting the input workspace.

    Returns
    -------
    The workspace with the midrange detector.
    """
    if output_workspace is None:
        output_workspace = str(input_workspace)
    clone_component_intensities(
        input_workspace,
        output_workspace=output_workspace,
        copy_logs=True,
        input_component="detector1",
        output_component="detector1",
    )
    clone_component_intensities(
        input_workspace,
        output_workspace=output_workspace,
        input_component="wing_detector",
        output_component="wing_detector",
    )
    return mtd[str(output_workspace)]
