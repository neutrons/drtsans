# flake8: noqa
import unittest
import math
from functools import reduce

from mantid.simpleapi import mtd
from reduction_workflow.instruments.sans.sns_command_interface import *
from reduction_workflow.instruments.sans.hfir_command_interface import *

import mock_api as eqsans


class EQSANS_api(unittest.TestCase):
    """
        Simple example scripts for EQSANS reduction,
        comparing the old and new APIs
    """

    def check_iq(self, iq_ws):
        """ Compare an I(q) result with reference data """
        q_ref = [0.00627899344457, 0.00446576953342, -0.000642819111528,
                 -0.101534608008, -0.310757393963, -0.131130075581,
                 -0.0323358349459, 0.00226873035507, 0.00467903730804,
                 0.00476547751556, 0.00454621777737, 0.00461420993593,
                 0.00444299343831, 0.00444274263385, 0.00442913410634,
                 0.0044129942791, 0.00423460966621, 0.00427546254601,
                 0.00441012526601, 0.00433252773473, 0.00432164034684,
                 0.00439791288994, 0.00469528641197, 0.00486253958713,
                 0.00480024141719, 0.00468809133506, 0.00481356988653,
                 0.00483214945498, 0.00483136210702, 0.00487216167692,
                 0.00492088466747, 0.00488815974658, 0.00514516516681,
                 0.00614152997835, 0.0085220782664, 0.0116563229355,
                 0.0140532844266, 0.0147721414837, 0.0142834882465,
                 0.0140842231454, 0.0140740183002, 0.0140196950919,
                 0.0141709266154, 0.0147495317946, 0.0145935952665,
                 0.0121393978471, 0.008514735452, 0.00597619704824,
                 0.00494511669476, 0.00451223808871, 0.00436704422346,
                 0.00423573250156, 0.00416240398741, 0.0040986674837,
                 0.00403915181736, 0.00401955881536, 0.00391070469289,
                 0.00388059685075, 0.00387805144206, 0.00378120272082,
                 0.00377189967075, 0.00362619922524, 0.00360440134901,
                 0.00343862502064, 0.00331973424475, 0.00317669581008,
                 0.00304771024322, 0.00287494022365, 0.00277546500525,
                 0.00284655272675, 0.00294610358984, 0.00305379983798,
                 0.00291370978474, 0.00312756430255, 0.00336163490614,
                 0.00347958253506, 0.00368613734536, 0.00409759160975,
                 0.00415979771501, 0.00424967221252, 0.00412307025487,
                 0.00404350132103, 0.00443294240065, 0.0039460980164,
                 0.00401991710877, 0.004356062666, 0.00430232102261,
                 0.00377806951331, 0.00347046410843, 0.00290652092552,
                 0.00254292463803, 0.00204285999082, 0.00145608832153,
                 0.00250844350263, 0.00123810738814, 0.000881609990444,
                 0.000446780891197, 9.50473596521e-05, 0.00034358746158, 0.0]

        iq = iq_ws.readY(0)
        diff = [math.fabs(q_ref[i] - iq[i]) < 0.0001 for i in range(7, 100)]
        output = reduce(lambda x, y: x and y, diff)
        if not output:
            print(iq)
        return output

    def test_simple_reduction_new_api(self):
        """ Simple reduction example using the new API """
        # Find beam center
        # This can be done like this:
        #     x, y = eqsans.find_beam_center()
        x, y = 96.29, 126.15

        ws = eqsans.load_events("EQSANS_104088",
                                beam_center_x=x, beam_center_y=y)
        ws = eqsans.prepare_data(ws, tubes=False)

        # Find transmission beam center, or use the one we have
        # Apply transmission
        ws = eqsans.apply_transmission(ws, 1.0, 0)

        # Now the background
        ws_bck = eqsans.load_events("EQSANS_105428",
                                    beam_center_x=x, beam_center_y=y)
        ws_bck = eqsans.prepare_data(ws_bck, tubes=False)

        # Find transmission beam center, or use the one we have
        # Apply transmission
        ws_bck = eqsans.apply_transmission(ws_bck, 1.0, 0)
        ws = eqsans.subtract_background(ws, ws_bck)

        iq = eqsans.iq(ws)
        self.assertTrue(self.check_iq(iq))

    def skip_test_simple_reduction_old_api(self):
        """ Simple reduction example using the old API """
        EQSANS()
        SetBeamCenter(96.29, 126.15)
        AppendDataFile("EQSANS_104088")
        UseConfig(False)
        UseConfigTOFTailsCutoff(False)
        UseConfigMask(False)
        TotalChargeNormalization(normalize_to_beam=False)
        SetTransmission(1.0, 0.0, False)
        Background("EQSANS_105428")
        Resolution(10)
        Reduce1D()

        self.assertTrue(self.check_iq(mtd['EQSANS_104088_Iq']))


if __name__ == '__main__':
    unittest.main()
