# flake8: noqa
from __future__ import (absolute_import, division, print_function)
import unittest

from mantid.simpleapi import mtd, ExtractMask
from reduction_workflow.instruments.sans.sns_command_interface import *
from reduction_workflow.instruments.sans.hfir_command_interface import *

import mock_api as eqsans


class EQSANS_api(unittest.TestCase):
    """
        More complete example scripts for EQSANS reduction,
        comparing the old and new APIs
    """

    def check_iq(self, iq_frame1_ws, iq_frame2_ws):
        """ Verify the reduced output """
        iq1 = iq_frame1_ws.readY(0)
        iq2 = iq_frame2_ws.readY(0)
        return (iq1[-1]-0.00875232) < 0.00001 and (iq2[-1]+2.40273e-05) < 0.00001

    def test_old_api(self):
        """ Real-life EQSANS reduction script """
        mask60_ws4m = Load(Filename="/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs")
        ws604m, masked60_detectors4m = ExtractMask(InputWorkspace=mask60_ws4m, OutputWorkspace="__edited_mask604m")
        detector_ids604m = [int(i) for i in masked60_detectors4m]
        Mask_BS604m = detector_ids604m

        EQSANS()
        SolidAngle(detector_tubes=True)
        DarkCurrent("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5")
        TotalChargeNormalization()
        SetAbsoluteScale(0.0208641883)
        AzimuthalAverage(n_bins=100, n_subpix=1, log_binning=True)
        IQxQy(nbins=75)
        MaskDetectors(Mask_BS604m)
        OutputPath("/tmp")
        UseConfigTOFTailsCutoff(True)
        UseConfigMask(True)
        ReductionSingleton().reduction_properties["SampleOffset"] = 340
        ReductionSingleton().reduction_properties["DetectorOffset"] = 0
        Resolution(sample_aperture_diameter=10)
        PerformFlightPathCorrection(True)
        DirectBeamCenter("88973")
        SensitivityCorrection("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017A_mp/Sensitivity_patched_thinPMMA_4m_79165_event.nxs", min_sensitivity=0.5, max_sensitivity=1.5, use_sample_dc=True)
        DivideByThickness(0.1)
        DirectBeamTransmission("88975", "88973", beam_radius=5)
        ThetaDependentTransmission(True)
        AppendDataFile(["88980"])
        CombineTransmissionFits(False)
        Background(["88979"])
        BckDirectBeamTransmission("88974", "88973", beam_radius=5)
        BckThetaDependentTransmission(True)
        BckCombineTransmissionFits(False)
        #SaveIq(process='None')
        Reduce()

        return self.check_iq(mtd['88980_frame1_Iq'], mtd['88980_frame2_Iq'])

    def test_new_api(self):
        """ Same example using the new API """

        mask60_ws4m = Load(Filename="/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs")
        ws604m, masked60_detectors4m = ExtractMask(InputWorkspace=mask60_ws4m, OutputWorkspace="__edited_mask604m")
        detector_ids604m = [int(i) for i in masked60_detectors4m]
        Mask_BS604m = detector_ids604m

        beam_flux_file = "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample"
        absolute_scale = 0.0208641883
        sample_thickness = 0.1  # mm

        # The following remains to be added
        """
        DarkCurrent("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5")

        MaskDetectors(Mask_BS604m)
        UseConfigMask(True)

        SensitivityCorrection("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017A_mp/Sensitivity_patched_thinPMMA_4m_79165_event.nxs", min_sensitivity=0.5, max_sensitivity=1.5, use_sample_dc=True)

        DirectBeamTransmission("EQSANS_88975", "EQSANS_88973", beam_radius=5)
        ThetaDependentTransmission(True)
        CombineTransmissionFits(False)

        BckDirectBeamTransmission("EQSANS_88974", "EQSANS_88973", beam_radius=5)
        BckThetaDependentTransmission(True)
        BckCombineTransmissionFits(False)
        """

        # Find beam center
        x, y = eqsans.find_beam_center("EQSANS_88973")

        ws = eqsans.load_events("EQSANS_88980",
                             sample_offset=340)
        ws = eqsans.prepare_data(ws, normalize_to_monitor=False,
                              beam_profile=beam_flux_file)

        # Apply transmission
        ws = eqsans.apply_transmission(ws, x, y)

        # Now the background
        ws_bck = eqsans.load_events("EQSANS_88979",
                                 sample_offset=340)
        ws_bck = eqsans.prepare_data(ws_bck)

        # Find transmission beam center, or use the one we have
        # Apply transmission
        ws_bck = eqsans.apply_transmission(ws_bck)
        ws = eqsans.subtract_background(ws, ws_bck)

        ws *= absolute_scale
        ws /= sample_thickness

        iq = eqsans.iq(ws, number_of_bins=100, log_binning=True,
                    sample_aperture=10.0)
        iqxqy = eqsans.iqxqy(ws, number_of_bins=75, log_binning=False)


if __name__ == '__main__':
    unittest.main()
