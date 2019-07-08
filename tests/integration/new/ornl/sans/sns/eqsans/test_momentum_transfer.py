#!/usr/bin/env python
from __future__ import print_function

import numpy as np

from mantid import mtd
from mantid.simpleapi import (AddSampleLog, ConfigService, CreateWorkspace,
                              ExtractSpectra, MoveInstrumentComponent, Rebin)
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.sans.sns.eqsans import normalisation
from ornl.settings import unique_workspace_name
from ornl.sans.sns.eqsans import load_events, transform_to_wavelength
# from ornl.sans.sns.eqsans import center_detector
# from ornl.sans.sns.eqsans import geometry


def legacy_reduction():

    from reduction_workflow.command_interface import AppendDataFile, Reduce
    from reduction_workflow.instruments.sans import (
        sns_command_interface as eqsans)
    import tempfile

    configI = ConfigService.Instance()
    configI["facilityName"] = 'SNS'

    eqsans.EQSANS()
    AppendDataFile('EQSANS_68200')
    eqsans.UseConfig(False)
    eqsans.UseConfigTOFTailsCutoff(False)
    eqsans.UseConfigMask(False)
    eqsans.SetBeamCenter(90, 132.5)
    # eqsans.SetTransmission(1.0, 1.0)
    eqsans.SetTOFTailsCutoff(low_cut=500, high_cut=2000)
    eqsans.SetWavelengthStep(step=0.1)
    eqsans.Resolution(sample_aperture_diameter=10.0)
    eqsans.NoTransmission()
    eqsans.NoNormalization()
    eqsans.NoSolidAngle()
    eqsans.TotalChargeNormalization()
    eqsans.OutputPath(tempfile.gettempdir())
    eqsans.AzimuthalAverage(n_bins=110)
    eqsans.IQxQy(nbins=110, log_binning=False)
    Reduce()

    return mtd['EQSANS_68200']


def test_momentum_tranfer_serial():
    '''

    '''

    ws = load_events('EQSANS_68200', detector_offset=0,
                     sample_offset=0)

    ws = transform_to_wavelength(ws, bin_width=0.1,
                                 low_tof_clip=500,
                                 high_tof_clip=2000)
    
    # center_detector will do this when jose fixes the Z issue
    # Center the beam: xyz: 0.0165214,0.0150392,4.00951
    instrument = ws.getInstrument()
    component = instrument.getComponentByName("detector1")
    z = component.getPos()[2]
    MoveInstrumentComponent(Workspace='ws', ComponentName='detector1',
                            X=-0.025, Y=-0.016, Z=z, RelativePosition=False)

    flux_ws = normalisation.load_beam_flux_file(
        '/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample',
        out_ws='flux_ws', ws_reference=ws)

    ws = normalisation.proton_charge_and_flux(ws, flux_ws, "ws")

    # This is not working! slit4 missing
    # geometry.sample_aperture_diameter(ws)
    # To date this does not add keyword to logs
    # geometry.source_aperture_diameter(ws)

    # this temporary. Waiting for jose to finish the issues above
    AddSampleLog(Workspace=ws, LogName='sample-aperture-diameter',
                 LogText='10.', LogType='Number', LogUnit='mm')
    AddSampleLog(Workspace=ws, LogName='source-aperture-diameter',
                 LogText='20.', LogType='Number', LogUnit='mm')

    rebin_start, rebin_end, rebin_step = 2.6, 5.6, 0.2

    ws = Rebin(
        InputWorkspace=ws,
        OutputWorkspace="ws_rebin",
        Params="{:.2f},{:.2f},{:.2f}".format(
            rebin_start, rebin_step, rebin_end)
    )

    bins = np.arange(rebin_start, rebin_end, rebin_step)

    for index, bin_start in enumerate(bins):

        ws_extracted = ExtractSpectra(
            InputWorkspace=ws,
            OutputWorkspace="slice_{:}".format(index),
            XMin=bin_start, XMax=bin_start+rebin_step)

        wavelength_mean = bin_start + rebin_step/2

        AddSampleLog(Workspace=ws_extracted, LogName='wavelength',
                     LogText="{:.2f}".format(wavelength_mean),
                     LogType='Number', LogUnit='Angstrom')

        AddSampleLog(Workspace=ws_extracted, LogName='wavelength-spread',
                     LogText='0.2', LogType='Number', LogUnit='Angstrom')
        ws_prefix = "ws_{}".format(index)

        # 2D
        wss_name_ws = bin_into_q2d(ws_extracted, out_ws_prefix=ws_prefix)
        assert len(wss_name_ws) == 3

        ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
        assert ws_iqxqy.extractY().shape == (256, 192)
        assert ws_iqxqy.extractX().shape == (256, 193)
        # 1D
        _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy,
                                out_ws_prefix=ws_prefix)
        assert ws_iq.extractY().shape == (1, 100)
        assert ws_iq.extractX().shape == (1, 101)

    # Calculate an average I(q)

    mean_x = np.mean([mtd["ws_{}_iq".format(index)].extractX() for
                      index, _ in enumerate(bins)], axis=0)
    mean_y = np.mean([mtd["ws_{}_iq".format(index)].extractY() for
                      index, _ in enumerate(bins)], axis=0)
    mean_e = np.sqrt(np.sum([(mtd["ws_{}_iq".format(index)].extractE())**2 for
                             index, _ in enumerate(bins)], axis=0))

    mean_iq = CreateWorkspace(
        DataX=mean_x,
        DataY=mean_y,
        DataE=mean_e,
        # Dx=dq_bin_centers,  # bin centers!!
        NSpec=1,
        UnitX='MomentumTransfer',
        YUnitLabel='Counts',
        OutputWorkspace="mean_iq"
    )

    assert mean_iq is not None
