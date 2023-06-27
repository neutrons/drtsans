from drtsans.instruments import fetch_idf, empty_instrument_workspace
from drtsans.mono.biosans.geometry import has_midrange_detector
from drtsans.tubecollection import TubeCollection

# third party imports
from mantid.api import mtd
from mantid.simpleapi import CopyLogs


def clone_component_intensities(
    input_workspace,
    output_workspace,
    input_component="wing_detector",
    output_component="midrange_detector",
    copy_logs=True,
):
    r"""
    Clone intensities into the detector pixels of a component from another component.

    Parameters
    ----------
    input_workspace : str, ~Mantid.api.Workspace
        The workspace to copy intensities from.

    output_workspace : str
        If :py:obj:`None`, the name of the input_workspace will be used, thus overwriting the input workspace.

    input_component : str
        The component to copy intensities from. Default: "wing_detector".

    output_component : str
        The component to copy intensities to. Default: "midrange_detector".

    copy_logs : bool
        If True, the sample logs will be copied to the output workspace

    Returns
    -------
    ~Mantid.api.Workspace
        The workspace with intensities inserted.
    """
    input_workspace = str(input_workspace)
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

    if input_component == "midrange_detector" and not has_midrange_detector(input_workspace):
        raise ValueError("No midrange detector in input workspace")

    input_collection = TubeCollection(input_workspace, input_component)
    output_collection = TubeCollection(output_workspace, output_component)
    # verify that the number of pixels and tubes is greater or equal in input_component
    if len(input_collection) < len(output_collection):
        raise ValueError(f"Not enough tubes to copy {len(input_collection)} < {len(output_collection)}.")
    if len(input_collection[0]) < len(output_collection[0]):
        raise ValueError(f"Not enough tube pixels to copy {len(input_collection[0])} < {len(output_collection[0])}.")

    # copy data tube for tube, pixel for pixel
    for input_tube, output_tube in zip(input_collection, output_collection):
        input_pixel_indexes = input_tube.spectrum_info_index
        output_pixel_indexes = output_tube.spectrum_info_index
        output_dataX = mtd[str(output_workspace)].dataX
        output_dataY = mtd[str(output_workspace)].dataY
        output_dataE = mtd[str(output_workspace)].dataE
        input_readX = mtd[str(input_workspace)].readX
        input_readY = mtd[str(input_workspace)].readY
        input_readE = mtd[str(input_workspace)].readE
        for input_pixel, output_pixel in zip(input_pixel_indexes, output_pixel_indexes):
            output_dataX(output_pixel)[:] = input_readX(input_pixel)
            output_dataY(output_pixel)[:] = input_readY(input_pixel)
            output_dataE(output_pixel)[:] = input_readE(input_pixel)

    return mtd[str(output_workspace)]
