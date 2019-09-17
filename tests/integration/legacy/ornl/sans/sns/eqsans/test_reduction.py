"""EQSANS reduction script
"""

import shutil
import tempfile
from os.path import join as pjn

import pytest

import reduction_workflow.command_interface as main_cli
import reduction_workflow.instruments.sans.hfir_command_interface as hfir_cli
import reduction_workflow.instruments.sans.sns_command_interface as sns_cli
from mantid.kernel import ConfigService
from mantid.simpleapi import (CompareWorkspaces, ExtractMask, Load,
                              LoadEventNexus)

shared_dir = '/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp'
ipts_dir = '/SNS/EQSANS/IPTS-20196/nexus'

output_dir = tempfile.mkdtemp()


@pytest.mark.skip(reason="Mantid master is failing for EQSANS reduction.")
def test_reduction(reference_dir):

    # Set specific configuration
    config = ConfigService.Instance()
    previous_instrument = config['instrumentName']
    config['instrumentName'] = 'EQSANS'
    previous_archive = config['datasearch.searcharchive']
    config['datasearch.searcharchive'] = 'hfir,sns'

    # Let's rock and roll!
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

        main_cli.Reduce()

        # Check reduction agrees with reference data
        ref_data_dir = pjn(reference_dir.legacy.eqsans, 'test_reduction')
        for file in ('92164_0600_Iq.txt',
                     '92164_0600_Iqxy.dat',
                     '92164_0600_transmission_fit.txt',
                     '92164_0600_transmission.txt'):
            outw = Load(Filename=pjn(output_dir, file))
            refw = Load(Filename=pjn(ref_data_dir, file))
            cmp, mesg = CompareWorkspaces(outw, refw, Tolerance=1e-04)
            assert cmp is True
    finally:
        # clean up your mess, baby
        shutil.rmtree(output_dir)
        config['instrumentName'] = previous_instrument
        config['datasearch.searcharchive'] = previous_archive


if __name__ == '__main__':
    pytest.main()

"""Original script in /SNS/EQSANS/IPTS-20196/shared/heller/test/porasil_slice1m.py
# EQSANS reduction script
# Script automatically generated on Sun Mar 19 06:21:58 2017

import mantid
from mantid.simpleapi import *
from reduction_workflow.instruments.sans.sns_command_interface import *

config = ConfigService.Instance()
config['instrumentName']='EQSANS'

# This is the 1.3m mask
mask60_ws1m = Load(Filename="/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/beamstop60_mask_4m.nxs")
ws601m, masked60_detectors1m = ExtractMask(InputWorkspace=mask60_ws1m, OutputWorkspace="__edited_mask601m")
detector_ids601m = [int(i) for i in masked60_detectors1m]
Mask_BS601m = detector_ids601m

# The 1.3m reduction

start_time = 0
time_window = 600
end_time = start_time + time_window

FSca1m = "92164"
WorkSpaceFilename = FSca1m + "_%04d" %(end_time)
FullFileName = "/SNS/EQSANS/IPTS-20196/nexus/EQSANS_"+FSca1m+".nxs.h5"
LoadEventNexus(Filename=FullFileName,OutputWorkspace=WorkSpaceFilename,FilterByTimeStart=str(start_time),FilterByTimeStop=str(end_time))

EQSANS()
SolidAngle(detector_tubes=True)
DarkCurrent("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/EQSANS_89157.nxs.h5")
TotalChargeNormalization()
SetAbsoluteScale(1.00)
AzimuthalAverage(n_bins=200, n_subpix=1, log_binning=True)
IQxQy(nbins=75)
MaskDetectors(Mask_BS601m)
OutputPath("/SNS/EQSANS/IPTS-20196/shared/heller/test")
UseConfigTOFTailsCutoff(True)
UseConfigMask(True)
ReductionSingleton().reduction_properties["SampleOffset"] = 340
ReductionSingleton().reduction_properties["DetectorOffset"] = 0
Resolution(sample_aperture_diameter=10)
PerformFlightPathCorrection(True)
DirectBeamCenter("92160")
SensitivityCorrection("/SNS/EQSANS/shared/NeXusFiles/EQSANS/2017B_mp/Sensitivity_patched_thinPMMA_1o3m_92313_event.nxs", min_sensitivity=0.5, max_sensitivity=1.5, use_sample_dc=True)
#the thickness is in cm...by using a value of 1.0, the data is essentially uncorrected for the sample thickness
#divide the reduced data by the thickness in cm to get into absolute units
DivideByThickness(0.1)
DirectBeamTransmission("92162", "92160", beam_radius=3)
ThetaDependentTransmission(True)
AppendDataFile([WorkSpaceFilename])
CombineTransmissionFits(False)
Background("92163")
BckDirectBeamTransmission("92161", "92160", beam_radius=3)
BckThetaDependentTransmission(True)
BckCombineTransmissionFits(False)
SaveIq(process='None')
Reduce()
"""  # noqa: E501
