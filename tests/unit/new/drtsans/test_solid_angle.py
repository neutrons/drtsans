from os.path import join
import pytest
from mantid.simpleapi import (CompareWorkspaces, Load,
                              LoadEmptyInstrument,
                              MoveInstrumentComponent)
from drtsans import solid_angle_correction


def test_solid_angle(reference_dir):
    # Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName='eqsans')

    # Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName='detector1',
                            RelativePosition='0', Z='5')

    # Apply solid angle correction
    wsOutput = solid_angle_correction(wsInput, detector_type='VerticalTube')

    # Let's do some validation
    assert wsOutput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(reference_dir.new.eqsans,
                                             'test_solid_angle.nxs'))
    assert CompareWorkspaces(wsOutput, reference_workspace)


def test_solid_angle_optional_output(reference_dir):
    # Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName='eqsans')

    # Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName='detector1',
                            RelativePosition='0', Z='5')

    # Apply solid angle correction
    wsOutput = solid_angle_correction(wsInput, detector_type='VerticalTube',
                                      output_workspace='wsOutput')

    # Let's do some validation
    assert wsOutput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(reference_dir.new.eqsans,
                                             'test_solid_angle.nxs'))
    assert CompareWorkspaces(wsOutput, reference_workspace)


def test_solid_angle_input_output(reference_dir):
    # Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName='eqsans')

    # Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName='detector1',
                            RelativePosition='0', Z='5')

    # Apply solid angle correction
    wsInput = solid_angle_correction(wsInput, detector_type='VerticalTube',
                                     output_workspace='wsInput')

    # Let's do some validation
    assert wsInput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(reference_dir.new.eqsans,
                                             'test_solid_angle.nxs'))
    assert CompareWorkspaces(wsInput, reference_workspace)


if __name__ == '__main__':
    pytest.main()
