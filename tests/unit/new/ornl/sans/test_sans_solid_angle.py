from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pytest
from mantid.simpleapi import LoadEmptyInstrument, MoveInstrumentComponent

from ornl.sans import solid_angle_correction


@pytest.fixture(scope='module')
def test_sans_solid_angle():
    instrument = 'eqsans'

    #Load empty instrument
    wsInput = LoadEmptyInstrument(InstrumentName=instrument)

    #Move the detector 5m away from the sample
    MoveInstrumentComponent(Workspace=wsInput, ComponentName='detector1', RelativePosition='0', Z='5')

    #Solid Angle correction
    wsOutput = solid_angle_correction(wsInput)

    #Let's do some validation
    assert wsOutput.getNumberHistograms(), 49153
    #SA is bigger in the center: 8.99095e-07 vs 9.4172e-07
    correctionEdge = wsOutput.y(48896)[0];
    correctionCenter = wsOutput.y(25984)[0];
    assert correctionCenter > correctionEdge;

