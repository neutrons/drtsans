# local imports
from drtsans.instruments import fetch_idf
from drtsans.mono.biosans.simulated_intensities import clone_component_intensities
from drtsans.geometry import spectrum_info_ranges

# third party imports
from mantid.simpleapi import LoadEmptyInstrument
import pytest


@pytest.mark.parametrize("reuse", [False, True])
def test_clone_component_intensities(reference_dir, reuse, temp_workspace_name):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS",
        Filename=fetch_idf("BIOSANS_Definition.xml"),
        OutputWorkspace=temp_workspace_name(),
    )

    # set intensities in wing detector pixels corresponding to first and last pixel in midrange detector
    golden_value = 42
    in_first, in_last = spectrum_info_ranges(workspace, "wing_detector")
    out_first, out_last = spectrum_info_ranges(workspace, "midrange_detector")
    workspace.dataY(in_first)[:] = golden_value
    workspace.dataY(in_first + (out_last - out_first) - 1)[:] = golden_value

    output_workspace = temp_workspace_name() if reuse is False else None
    output_workspace = clone_component_intensities(
        workspace,
        output_workspace=output_workspace,
        input_component="wing_detector",
        output_component="midrange_detector",
    )

    # assert
    assert output_workspace.readY(out_first) == golden_value
    assert output_workspace.readY(out_last - 1) == golden_value


if __name__ == "__main__":
    pytest.main([__file__])
