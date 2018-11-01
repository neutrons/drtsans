"""EQSANS reduction script
"""
from __future__ import (absolute_import, division, print_function)

import pytest
import os
from os.path import join as pjn

from mantid.kernel import ConfigService
from mantid.simpleapi import (Load, ExtractMask, LoadEventNexus)

import reduction_workflow.instruments.sans.sns_command_interface as sns_cli
import reduction_workflow.instruments.sans.hfir_command_interface as hfir_cli
import reduction_workflow.command_interface as main_cli

shared_dir = '/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp'
ipts_dir = '/SNS/EQSANS/IPTS-20196/nexus'
output_dir = '/tmp/sans_rewrite/tests/integration/' \
             'legacy/ornl/sans/eqsans/test_reduction'
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

def test_reduction():
    # Set instrument
    config = ConfigService.Instance()
    previous_instrument = config['instrumentName']
    config['instrumentName'] = 'EQSANS'
    try:
        # This is the 1.3m mask
        mask60_ws1m = Load(Filename=pjn(shared_dir, 'beamstop60_mask_4m.nxs'))
        ws601m, masked60_detectors1m =\
            ExtractMask(mask60_ws1m, OutputWorkspace="__edited_mask601m")
        mask_bs601m = [int(i) for i in masked60_detectors1m]
        #
        start_time, time_window = 0, 600
        end_time = start_time + time_window
        FSca1m = '92164'
        work_space_filename = '{}_{:04d}'.format(FSca1m, end_time)
        full_file_name = pjn(ipts_dir, 'EQSANS_{}.nxs.h5'.format(FSca1m))
        LoadEventNexus(Filename=full_file_name,
                       OutputWorkspace=work_space_filename,
                       FilterByTimeStart=str(start_time),
                       FilterByTimeStop=str(end_time))
        #
        sns_cli.EQSANS()
        sns_cli.TotalChargeNormalization()
        sns_cli.AzimuthalAverage(n_bins=200, n_subpix=1, log_binning=True)
        sns_cli.UseConfigTOFTailsCutoff(True)
        sns_cli.UseConfigMask(True)
        sns_cli.Resolution(sample_aperture_diameter=10)
        sns_cli.PerformFlightPathCorrection(True)
        sns_cli.CombineTransmissionFits(False)
        sns_cli.BckCombineTransmissionFits(False)
        sns_cli.SaveIq(process='None')

        hfir_cli.SolidAngle(detector_tubes=True)
        hfir_cli.DarkCurrent(pjn(shared_dir, 'EQSANS_89157.nxs.h5'))
        hfir_cli.SetAbsoluteScale(1.00)
        hfir_cli.IQxQy(nbins=75)
        hfir_cli.MaskDetectors(mask_bs601m)
        hfir_cli.DirectBeamCenter('92160')
        sensitivity_file = 'Sensitivity_patched_thinPMMA_1o3m_92313_event.nxs'
        hfir_cli.SensitivityCorrection(pjn(shared_dir, sensitivity_file),
                                       min_sensitivity=0.5,
                                       max_sensitivity=1.5,
                                       use_sample_dc=True)
        # thickness in cm...by using a value of 1.0, the data is uncorrected
        # for the sample thickness. Divide the reduced data by the thickness
        # in cm to get into absolute units
        hfir_cli.DivideByThickness(0.1)
        hfir_cli.DirectBeamTransmission('92162', '92160', beam_radius=3)
        hfir_cli.ThetaDependentTransmission(True)
        hfir_cli.Background('92163')
        hfir_cli.BckDirectBeamTransmission('92161', '92160', beam_radius=3)
        hfir_cli.BckThetaDependentTransmission(True)

        reduction_pr = main_cli.ReductionSingleton().reduction_properties
        reduction_pr['SampleOffset'] = 34
        reduction_pr['DetectorOffset'] = 0

        main_cli.OutputPath(output_dir)
        main_cli.AppendDataFile([work_space_filename])
        """
        #
        # Original ordering
        #
        sns_cli.EQSANS()
        hfir_cli.SolidAngle(detector_tubes=True)
        hfir_cli.DarkCurrent(pjn(shared_dir, 'EQSANS_89157.nxs.h5'))
        sns_cli.TotalChargeNormalization()
        hfir_cli.SetAbsoluteScale(1.00)
        sns_cli.AzimuthalAverage(n_bins=200, n_subpix=1, log_binning=True)
        hfir_cli.IQxQy(nbins=75)
        hfir_cli.MaskDetectors(Mask_BS601m)
        main_cli.OutputPath(output_dir)
        sns_cli.UseConfigTOFTailsCutoff(True)
        sns_cli.UseConfigMask(True)
        reduction_pr = main_cli.ReductionSingleton().reduction_properties
        reduction_pr['SampleOffset'] = 34
        reduction_pr['DetectorOffset'] = 0
        sns_cli.Resolution(sample_aperture_diameter=10)
        sns_cli.PerformFlightPathCorrection(True)
        hfir_cli.DirectBeamCenter('92160')
        sensitivity_file = 'Sensitivity_patched_thinPMMA_1o3m_92313_event.nxs'
        hfir_cli.SensitivityCorrection(pjn(shared_dir, sensitivity_file),
                                       min_sensitivity=0.5,
                                       max_sensitivity=1.5,
                                       use_sample_dc=True)
        # thickness in cm...by using a value of 1.0, the data is uncorrected
        # for the sample thickness. Divide the reduced data by the thickness
        # in cm to get into absolute units
        hfir_cli.DivideByThickness(0.1)
        hfir_cli.DirectBeamTransmission('92162', '92160', beam_radius=3)
        hfir_cli.ThetaDependentTransmission(True)
        main_cli.AppendDataFile([WorkSpaceFilename])
        sns_cli.CombineTransmissionFits(False)
        hfir_cli.Background('92163')
        hfir_cli.BckDirectBeamTransmission('92161', '92160', beam_radius=3)
        hfir_cli.BckThetaDependentTransmission(True)
        sns_cli.BckCombineTransmissionFits(False)
        sns_cli.SaveIq(process='None')
        """
        main_cli.Reduce()
    finally:
        config['instrumentName'] = previous_instrument


if __name__ == '__main__':
    pytest.main()