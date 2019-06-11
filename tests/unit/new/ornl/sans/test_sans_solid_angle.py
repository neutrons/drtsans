from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from os.path import join
import pytest
from mantid.simpleapi import (CompareWorkspaces, Load,
                              LoadEmptyInstrument,
                              MoveInstrumentComponent)
from ornl.sans import solid_angle_correction as sac


def test_sans_solid_angle(refd):
    #Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName='eqsans')

    #Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName='detector1', RelativePosition='0', Z='5')

    #Apply solid angle correction
    wsOutput = sac.solid_angle_correction(wsInput, detector_type='Normal')

    #Let's do some validation
    assert wsOutput.getNumberHistograms(), 49153
    reference_workspace = Load(Filename=join(refd.new.eqsans, 'test_sans_solid_angle.nxs'))
    assert CompareWorkspaces(wsOutput, reference_workspace)


if __name__ == '__main__':
    pytest.main()
