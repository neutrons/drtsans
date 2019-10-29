import pytest

# https://docs.mantidproject.org/nightly/algorithms/LoadHFIRSANS-v1.html
# https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
# https://docs.mantidproject.org/nightly/algorithms/MoveInstrumentComponent-v1.htm
from mantid.simpleapi import (FindCenterOfMassPosition, LoadHFIRSANS,
                              MoveInstrumentComponent)
from drtsans.mono.biosans import beam_finder
from drtsans.samplelogs import SampleLogs

# Note for testing beam center: The FindCenterOfMassPosition algorithm
# needs some pixels outside the given examples, to be able to perform
# the iterative approach. So, for a 2x2 example, we need a 4x4 instrument
# to pad with zeros on each side


def _verify_pixel(wksp, index, position, counts):
    assert wksp.readY(index)[0] == counts, 'wksp_index={}'.format(index)
    assert wksp.detectorInfo().position(index) == pytest.approx(position), 'wksp_index={}'.format(index)


@pytest.mark.parametrize('generic_workspace',
                         [{'intensities': [[0,   0,    0,    0,    0,    0,    0,    0,   0, 0],
                                           [0, 173,  261,  327,  334,  307,  226,  146,  52, 0],
                                           [0, 351,  548,  740,  879,  958,  862,  503, 287, 0],
                                           [0, 317,  526,  663,  881,  896,  839,  583, 354, 0],
                                           [0, 685, 1173, 1380, 1726, 1635, 1383, 1022, 589, 0],
                                           [0, 300,  474,  592,  696,  783,  746,  541, 416, 0],
                                           [0, 327,  505,  716,  817,  842,  798,  674, 471, 0],
                                           [0, 134,  201,  269,  310,  292,  226,  148,  74, 0],
                                           [0,  51,   58,  101,  154,  138,  168,  104,  64, 0],
                                           [0,   0,    0,    0,    0,    0,    0,   0,    0, 0]],
                           'xc': 0.53075, 'dx': 0.0055,  # meters
                           'yc': 0.54675, 'dy': 0.0045,  # meters
                           'zc': 15.5}],  # meters
                         indirect=True)
def test_beam_finder_excel(generic_workspace):
    r"""
    Testing section 3.1 in the master document
    Functions to test: drtsans.mono.biosans.find_beam_center

    Underlying Mantid algorithms:
      FindCenterOfMassPosition https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
    This is different detector position from :func:`test_beam_finder_excel2`, but otherwise the same
    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Venky Pingali <pingalis@ornl.gov>
    """
    ws = generic_workspace  # friendly name
    SampleLogs(ws).insert('wavelength', 6.0, 'Angstrom')

    # verify the data is in correctly
    _verify_pixel(ws, 11, [0.5500, 0.5310, 15.5], 173)
    _verify_pixel(ws, 21, [0.5445, 0.5310, 15.5], 351)
    _verify_pixel(ws, 31, [0.5390, 0.5310, 15.5], 317)
    _verify_pixel(ws, 81, [0.5115, 0.5310, 15.5], 51)
    _verify_pixel(ws, 88, [0.5115, 0.5625, 15.5], 64)

    # run the function to calculate the beam center
    x, y, y_gravity = beam_finder.find_beam_center(ws, sample_det_cent_main_detector=15.5,
                                                   sample_det_cent_wing_detector=1.13)

    # within .1mm
    assert x == pytest.approx(0.5331, abs=0.0001)
    assert y == pytest.approx(0.5468, abs=0.0001)
    assert y_gravity == pytest.approx(0.54675 + 0.002694, abs=0.0001)


@pytest.mark.parametrize('generic_workspace',
                         [{'intensities': [[0,   0,    0,    0,    0,    0,    0,    0,   0, 0],
                                           [0, 173,  261,  327,  334,  307,  226,  146,  52, 0],
                                           [0, 351,  548,  740,  879,  958,  862,  503, 287, 0],
                                           [0, 317,  526,  663,  881,  896,  839,  583, 354, 0],
                                           [0, 685, 1173, 1380, 1726, 1635, 1383, 1022, 589, 0],
                                           [0, 300,  474,  592,  696,  783,  746,  541, 416, 0],
                                           [0, 327,  505,  716,  817,  842,  798,  674, 471, 0],
                                           [0, 134,  201,  269,  310,  292,  226,  148,  74, 0],
                                           [0,  51,   58,  101,  154,  138,  168,  104,  64, 0],
                                           [0,   0,    0,    0,    0,    0,    0,   0,    0, 0]],
                           'xc': 0.01075, 'dx': 0.0055,  # meters
                           'yc': -0.01025, 'dy': 0.0045,  # meters
                           'zc': 15.5}],  # meters
                         indirect=True)
def test_beam_finder_excel2(generic_workspace):
    r"""
    Testing section 3.1 in the master document
    Functions to test: drtsans.mono.biosans.find_beam_center

    Underlying Mantid algorithms:
      FindCenterOfMassPosition https://docs.mantidproject.org/nightly/algorithms/FindCenterOfMassPosition-v2.html
    This is different detector position from :func:`test_beam_finder_excel`, but otherwise the same
    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Venky Pingali <pingalis@ornl.gov>
    """
    ws = generic_workspace  # friendly name
    SampleLogs(ws).insert('wavelength', 6.0, 'Angstrom')

    # verify the data is in correctly
    _verify_pixel(ws, 11, [0.0300, -0.0260, 15.5], 173)
    _verify_pixel(ws, 21, [0.0245, -0.0260, 15.5], 351)
    _verify_pixel(ws, 31, [0.0190, -0.0260, 15.5], 317)
    _verify_pixel(ws, 81, [-0.0085, -0.0260, 15.5], 51)
    _verify_pixel(ws, 88, [-0.0085,  0.0055, 15.5], 64)

    # run the function to calculate the beam center
    x, y, y_gravity = beam_finder.find_beam_center(ws, sample_det_cent_main_detector=15.5,
                                                   sample_det_cent_wing_detector=1.13)

    # within .1mm
    assert x == pytest.approx(0.0131, abs=0.0001)
    assert y == pytest.approx(-0.0102, abs=0.0001)
    assert y_gravity == pytest.approx(-0.0102 + 0.002694, abs=0.0001)


def test_beam_finder(biosans_f):
    '''
    Test with the new beam finder

    1. Find the beamcenter x,y
    2. Move detector1 x,y according to beamcenter x,y
    3. Find gravity
    4. Move wing_detector y according to Gravity
    '''

    ws = LoadHFIRSANS(Filename=biosans_f['beamcenter'])

    # 0.00144037741238 -0.0243732351545 -0.0267
    x, y, y_gravity = beam_finder.find_beam_center(ws)

    assert x == pytest.approx(0.0014, abs=1e-3)
    assert y == pytest.approx(-0.0243, abs=1e-3)
    assert y_gravity == pytest.approx(-0.0216, abs=1e-3)

    instrument = ws.getInstrument()
    pos_main = instrument.getComponentByName("detector1").getPos()
    pos_wing = instrument.getComponentByName("wing_detector").getPos()

    # Let's center the instrument and get the new center:
    # Move down and left
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='detector1', X=-x, Y=-y)

    pos_main_1 = instrument.getComponentByName("detector1").getPos()
    pos_wing_1 = instrument.getComponentByName("wing_detector").getPos()

    assert pos_main[0] - pos_main_1[0] == x
    assert pos_main[1] - pos_main_1[1] == y
    assert pos_wing == pos_wing_1  # we did not touch it

    # Now let's correct the wing detector for the gravity drop
    # Relative movement up words
    MoveInstrumentComponent(
        Workspace=ws, ComponentName='wing_detector', X=-x, Y=-y_gravity)

    pos_main_2 = instrument.getComponentByName("detector1").getPos()
    pos_wing_2 = instrument.getComponentByName("wing_detector").getPos()

    assert pos_main_1 == pos_main_2
    assert pos_wing_2[1] == pos_main_2[1] + (abs(y_gravity) - abs(y))

    # After the re-centring we should be at (0,0)
    # Note that to give the same results we need to enter the center
    # estimates as the previous results!
    center = FindCenterOfMassPosition(InputWorkspace=ws,
                                      CenterX=-x, CenterY=-y)
    x1, y1 = center
    # Tolerance 1e-3 == millimeters
    assert x1 == pytest.approx(0.0, abs=1e-3)
    assert y1 == pytest.approx(0.0, abs=1e-3)

    # let's the test our wrap function. The results should be the same.
    x2, y2, y_gravity2 = beam_finder.find_beam_center(ws, centering_options=dict(CenterX=-x, CenterY=-y))

    assert x2 == pytest.approx(0.0, abs=1e-3) == x1
    assert y2 == pytest.approx(0.0, abs=1e-3) == y1
    assert y_gravity2 == pytest.approx(0.0 + y_gravity - y, abs=1e-3)
    assert abs(y_gravity2) > abs(y2)


if __name__ == '__main__':
    pytest.main([__file__])
