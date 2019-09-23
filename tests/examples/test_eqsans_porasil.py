# flake8: noqa
import unittest
import numpy as np

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
        return np.fabs(iq1[-1]-0.00875232) < 0.00001 and np.fabs(iq2[-1]+2.40273e-05) < 0.00001

    def test_old_api(self):
        """ Real-life EQSANS reduction script """
        mtd.clear()
        mask60_ws4m = Load(Filename="/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs")
        ws604m, masked60_detectors4m = ExtractMask(InputWorkspace=mask60_ws4m, OutputWorkspace="__edited_mask604m")
        detector_ids604m = [int(i) for i in masked60_detectors4m]

        EQSANS()
        SolidAngle(detector_tubes=True)
        DarkCurrent("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5")
        TotalChargeNormalization()
        SetAbsoluteScale(0.0208641883)
        AzimuthalAverage(n_bins=100, n_subpix=1, log_binning=True)
        IQxQy(nbins=75)
        MaskDetectors(detector_ids604m)
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

        ref = np.loadtxt('data/ref_porasil_complete_frame1.txt')
        iq = mtd['88980_frame1_Iq'].readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

        ref = np.loadtxt('data/ref_porasil_complete_frame2.txt')
        iq = mtd['88980_frame2_Iq'].readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

    def test_old_api_load_only(self):
        """ Real-life EQSANS reduction script: first part only """
        mtd.clear()
        EQSANS()
        TotalChargeNormalization()
        SetBeamCenter(90.927, 131.474)
        AzimuthalAverage(n_bins=100, n_subpix=1, log_binning=True)
        OutputPath("/tmp")
        UseConfigTOFTailsCutoff(True)
        UseConfigMask(True)
        ReductionSingleton().reduction_properties["SampleOffset"] = 340
        ReductionSingleton().reduction_properties["DetectorOffset"] = 0
        Resolution(sample_aperture_diameter=10)
        PerformFlightPathCorrection(True)
        AppendDataFile(["88980"])
        SetTransmission(1.0, 0.0, False)
        SaveIq(process='None')
        Reduce()

        ref = np.loadtxt('data/ref_load_only_frame1.txt')
        iq = mtd['88980_frame1_Iq'].readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

        ref = np.loadtxt('data/ref_load_only_frame2.txt')
        iq = mtd['88980_frame2_Iq'].readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

    def test_new_api_load_only(self):
        """ Real-life EQSANS reduction script: first part only """
        mtd.clear()
        beam_flux_file = "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample"
        x = 90.927
        y = 131.474

        ws = eqsans.load_events("EQSANS_88980",
                                beam_center_x=x, beam_center_y=y,
                                use_config_tof_cuts=True,
                                use_config=True,
                                use_config_mask=True,
                                correct_for_flight_path=True,
                                sample_offset=340)
        ws = eqsans.prepare_data(ws, normalize_to_monitor=False,
                                 beam_profile=beam_flux_file,
                                 tubes=False)

        _, iq_f2 = eqsans.iq(ws, number_of_bins=100, log_binning=True,
                             sample_aperture=10.0)

        ref = np.loadtxt('data/ref_load_only_frame2.txt')
        q = iq_f2.readX(0)
        iq = iq_f2.readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

    def test_old_api_no_transmission(self):
        """ Real-life EQSANS reduction script """
        mtd.clear()
        EQSANS()
        SolidAngle(detector_tubes=True)
        DarkCurrent("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5")
        TotalChargeNormalization()
        SetAbsoluteScale(0.0208641883)
        AzimuthalAverage(n_bins=100, n_subpix=1, log_binning=True)
        OutputPath("/tmp")
        UseConfigTOFTailsCutoff(True)
        UseConfigMask(True)
        ReductionSingleton().reduction_properties["SampleOffset"] = 340
        ReductionSingleton().reduction_properties["DetectorOffset"] = 0
        Resolution(sample_aperture_diameter=10)
        PerformFlightPathCorrection(True)
        DirectBeamCenter("88973")
        DivideByThickness(0.1)
        AppendDataFile(["88980"])
        SetTransmission(1.0, 0.0, False)
        Reduce()

        ref = np.loadtxt('data/ref_no_transmission_frame1.txt')
        iq = mtd['88980_frame1_Iq'].readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

        ref = np.loadtxt('data/ref_no_transmission_frame2.txt')
        iq = mtd['88980_frame2_Iq'].readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        self.assertTrue(np.sum(residuals) < 0.01)

    def test_new_api_no_transmission(self):
        mtd.clear()
        beam_flux_file = "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample"
        absolute_scale = 0.0208641883
        sample_thickness = 0.1  # mm
        dark_data_file = "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5"

        # Find beam center
        # The beam center found is slightly different because
        # the new API simply called the beam finder with default options,
        # which include the extra loading options.
        x, y = eqsans.find_beam_center("EQSANS_88973")
        self.assertAlmostEqual(x, 90.927, delta=0.1)
        self.assertAlmostEqual(y, 131.474, delta=0.1)

        ws = eqsans.load_events("EQSANS_88980",
                                beam_center_x=x, beam_center_y=y,
                                use_config_tof_cuts=True,
                                use_config=True,
                                use_config_mask=True,
                                correct_for_flight_path=True,
                                sample_offset=340)
        ws = eqsans.prepare_data(ws, normalize_to_monitor=False,
                                 beam_profile=beam_flux_file,
                                 dark_data_file=dark_data_file)

        ws *= absolute_scale
        ws /= sample_thickness

        _, iq_f2 = eqsans.iq(ws, number_of_bins=100, log_binning=True,
                             sample_aperture=10.0)

        ref = np.loadtxt('data/ref_no_transmission_frame2.txt')
        iq = iq_f2.readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2
        # Skip the garbage points at the lowest and highest bins
        self.assertTrue(np.sum(residuals[15:-5])/len(residuals[15:-5]) < 1.0)

    def test_new_api(self):
        """ Same example using the new API """
        mtd.clear()
        mask60_ws4m = Load(Filename="/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs")
        ws604m, masked60_detectors4m = ExtractMask(InputWorkspace=mask60_ws4m, OutputWorkspace="__edited_mask604m")
        detector_ids604m = [int(i) for i in masked60_detectors4m]

        beam_flux_file = "/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample"
        absolute_scale = 0.0208641883
        sample_thickness = 0.1  # mm
        dark_data_file = "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_86275.nxs.h5"
        sensitivity_file = "/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017A_mp/Sensitivity_patched_thinPMMA_4m_79165_event.nxs"

        # Find beam center
        x, y = eqsans.find_beam_center("EQSANS_88973")

        ws = eqsans.load_events("EQSANS_88980",
                                beam_center_x=x, beam_center_y=y,
                                use_config_tof_cuts=True,
                                use_config=True,
                                use_config_mask=True,
                                correct_for_flight_path=True,
                                sample_offset=340)
        ws = eqsans.prepare_data(ws, normalize_to_monitor=False,
                                 beam_center_x=x, beam_center_y=y,
                                 sensitivity_file=sensitivity_file,
                                 beam_profile=beam_flux_file,
                                 masked_detector_list=detector_ids604m,
                                 dark_data_file=dark_data_file)

        # Apply transmission
        ws = eqsans.compute_and_apply_transmission(ws,
                                                   sample_data_file="/SNS/EQSANS/IPTS-19800/nexus/EQSANS_88975.nxs.h5",
                                                   empty_data_file="/SNS/EQSANS/IPTS-19800/nexus/EQSANS_88973.nxs.h5",
                                                   dark_data_file=dark_data_file,
                                                   beam_center_x=x, beam_center_y=y)

        # Now the background
        ws_bck = eqsans.load_events("EQSANS_88979",
                                    beam_center_x=x, beam_center_y=y,
                                    use_config_tof_cuts=True,
                                    use_config=True,
                                    use_config_mask=True,
                                    correct_for_flight_path=True,
                                    sample_offset=340)
        ws_bck = eqsans.prepare_data(ws_bck, normalize_to_monitor=False,
                                     beam_center_x=x, beam_center_y=y,
                                     sensitivity_file=sensitivity_file,
                                     beam_profile=beam_flux_file,
                                     masked_detector_list=detector_ids604m,
                                     dark_data_file=dark_data_file)

        # Apply transmission
        ws_bck = eqsans.compute_and_apply_transmission(ws_bck,
                                                       sample_data_file="/SNS/EQSANS/IPTS-19800/nexus/EQSANS_88974.nxs.h5",
                                                       empty_data_file="/SNS/EQSANS/IPTS-19800/nexus/EQSANS_88973.nxs.h5",
                                                       dark_data_file=dark_data_file,
                                                       beam_center_x=x, beam_center_y=y)
        ws = eqsans.subtract_background(ws, ws_bck)

        ws *= absolute_scale
        ws /= sample_thickness

        _, iq_f2 = eqsans.iq(ws, number_of_bins=100, log_binning=True,
                       sample_aperture=10.0)

        ref = np.loadtxt('data/ref_porasil_complete_frame2.txt')
        iq = iq_f2.readY(0)
        pos = ref.T[2]>0
        residuals = (ref.T[1][pos] - iq[pos])**2 / ref.T[2][pos]**2

        # Skip the garbage points at the lowest and highest bins
        self.assertTrue(np.sum(residuals[15:-5])/len(residuals[15:-5]) < 0.1)

if __name__ == '__main__':
    unittest.main()
