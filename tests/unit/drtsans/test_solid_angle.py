import numpy as np
from os.path import join
import pytest
from mantid.simpleapi import (
    CompareWorkspaces,
    Load,
    LoadEmptyInstrument,
    MoveInstrumentComponent,
)
from drtsans import calculate_solid_angle, solid_angle_correction


@pytest.mark.datarepo
def test_solid_angle(datarepo_dir, clean_workspace):
    pytest.skip("This test fails, defect written up in EWM Defect 2841")
    # Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName="eqsans")
    clean_workspace(wsInput)

    # Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName="detector1", RelativePosition="0", Z="5")

    # Apply solid angle correction
    wsOutput = solid_angle_correction(wsInput, detector_type="VerticalTube")
    clean_workspace(wsOutput)

    # Let's do some validation
    assert wsOutput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(datarepo_dir.eqsans, "test_solid_angle.nxs"))
    result, messages = CompareWorkspaces(wsOutput, reference_workspace)
    clean_workspace(messages)
    assert result


@pytest.mark.datarepo
def test_solid_angle_optional_output(datarepo_dir, clean_workspace):
    pytest.skip("This test fails, defect written up in EWM Defect 2841")
    # Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName="eqsans")
    clean_workspace(wsInput)

    # Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName="detector1", RelativePosition="0", Z="5")

    # Apply solid angle correction
    wsOutput = solid_angle_correction(
        wsInput, detector_type="VerticalTube", output_workspace=clean_workspace("wsOutput")
    )

    # Let's do some validation
    assert wsOutput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(datarepo_dir.eqsans, "test_solid_angle.nxs"))
    result, messages = CompareWorkspaces(wsOutput, reference_workspace)
    clean_workspace(messages)
    assert result


@pytest.mark.datarepo
def test_solid_angle_input_output(datarepo_dir, clean_workspace):
    pytest.skip("This test fails, defect written up in EWM Defect 2841")
    # Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName="eqsans")
    clean_workspace(wsInput)

    # Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName="detector1", RelativePosition="0", Z="5")

    # Apply solid angle correction
    wsInput = solid_angle_correction(wsInput, detector_type="VerticalTube", output_workspace="wsInput")

    # Let's do some validation
    assert wsInput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(datarepo_dir.eqsans, "test_solid_angle.nxs"))
    result, messages = CompareWorkspaces(wsInput, reference_workspace)
    clean_workspace(messages)
    assert result


@pytest.mark.parametrize(
    "instrument, num_monitor",
    [("BIOSANS", 2), ("GPSANS", 2), ("EQSANS", 1)],
    ids=["BIOSANS", "GPSANS", "EQSANS"],
)
def test_solid_angle_calculation(instrument, num_monitor, clean_workspace):
    ws_empty = LoadEmptyInstrument(InstrumentName=instrument)
    clean_workspace(ws_empty)
    ws = calculate_solid_angle(ws_empty)
    clean_workspace(ws)

    values = ws.extractY()
    # monitors have zero solid angle
    np.testing.assert_equal(values[:1], 0.0)
    # everything other than monitors should have a positive solid angle
    assert np.all(values[num_monitor:] > 0.0)


if __name__ == "__main__":
    pytest.main([__file__])
