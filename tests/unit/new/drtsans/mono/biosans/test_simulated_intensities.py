# local imports
from drtsans.instruments import fetch_idf
from drtsans.mono.biosans.simulated_intensities import clone_component_intensities
from drtsans.tubecollection import TubeCollection

# third party imports
from mantid.api import mtd
from mantid.simpleapi import DeleteWorkspace, LoadEmptyInstrument
import pytest

# standard imports


@pytest.mark.parametrize("output_workspace_param", ["test_clone_output", None])
def test_clone_component_intensities(reference_dir, output_workspace_param):
    LoadEmptyInstrument(
        InstrumentName="BIOSANS",
        Filename=fetch_idf("BIOSANS_Definition.xml"),
        OutputWorkspace="test_clone_input",
    )
    input_workspace = mtd["test_clone_input"]

    # set intensities in wing detector pixels corresponding to first and last pixel in midrange detector
    wing_collection = TubeCollection(input_workspace, "wing_detector")
    midrange_collection = TubeCollection(input_workspace, "midrange_detector")
    wing_pixel_first = wing_collection[0][0].spectrum_info_index
    wing_pixel_last = wing_collection[len(midrange_collection) - 1][-1].spectrum_info_index
    input_workspace.dataY(wing_pixel_first)[:] = 44.0
    input_workspace.dataY(wing_pixel_last)[:] = 44.0

    output_workspace = clone_component_intensities(
        input_workspace, output_workspace_param, input_component="wing_detector", output_component="midrange_detector"
    )

    # assert
    midrange_pixel_first = midrange_collection[0][0].spectrum_info_index
    midrange_pixel_last = midrange_collection[-1][-1].spectrum_info_index
    assert output_workspace.readY(midrange_pixel_first) == 44.0
    assert output_workspace.readY(midrange_pixel_last) == 44.0

    if str(input_workspace) != str(output_workspace):
        DeleteWorkspace(output_workspace)
    DeleteWorkspace(input_workspace)


if __name__ == "__main__":
    pytest.main([__file__])
