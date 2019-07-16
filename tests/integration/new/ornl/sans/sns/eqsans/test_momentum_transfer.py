#!/usr/bin/env python
from __future__ import print_function

import numpy as np

from mantid import mtd
from mantid.simpleapi import (AddSampleLog, ConfigService, CreateWorkspace,
                              ExtractSpectra, MoveInstrumentComponent, Rebin,
                              ReplaceSpecialValues)
from ornl.sans.momentum_transfer import MomentumTransfer
from ornl.sans.sns.eqsans import normalisation
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
    eqsans.SetTOFTailsCutoff(low_cut=500, high_cut=2000)
    eqsans.SetWavelengthStep(step=0.1)
    eqsans.Resolution(sample_aperture_diameter=10.0)
    eqsans.NoTransmission()
    eqsans.NoNormalization()
    eqsans.NoSolidAngle()
    eqsans.TotalChargeNormalization()
    eqsans.OutputPath(tempfile.gettempdir())
    eqsans.AzimuthalAverage(n_bins=100)
    eqsans.IQxQy(nbins=100, log_binning=False)
    Reduce()

    return mtd['EQSANS_68200_iq']


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
    MoveInstrumentComponent(Workspace=ws, ComponentName='detector1',
                            X=-0.025, Y=-0.016, Z=z, RelativePosition=False)

    flux_ws = normalisation.load_beam_flux_file(
        '/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample',
        output_workspace='flux_ws', ws_reference=ws)

    ws = normalisation.normalise_by_proton_charge_and_flux(ws, flux_ws, "ws")

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
        wss_name_ws = bin_into_q2d(
            ws_extracted, out_ws_prefix=ws_prefix, bins=[100, 100])
        assert len(wss_name_ws) == 3

        ws_iqxqy, ws_dqx, ws_dqy = [ws[1] for ws in wss_name_ws]
        # assert ws_iqxqy.extractY().shape == (256, 192)
        # assert ws_iqxqy.extractX().shape == (256, 193)
        # 1D
        _, ws_iq = bin_into_q1d(ws_iqxqy, ws_dqx, ws_dqy,
                                out_ws_prefix=ws_prefix,
                                # The same binning as legacy
                                bins=np.linspace(0.00119358, 0.460477, 100))
        # assert ws_iq.extractY().shape == (1, 100)
        # assert ws_iq.extractX().shape == (1, 101)

    all_bins = np.array([mtd["ws_{}_iq".format(index)].extractX() for
                         index, _ in enumerate(bins)])
    q_max = all_bins.max()
    q_min = all_bins.min()

    [Rebin(
        InputWorkspace=mtd["ws_{}_iq".format(index)],
        OutputWorkspace="ws_rebin_{}_iq".format(index),
        Params='{},0.004,{}'.format(q_min, q_max)) for
        index, _ in enumerate(bins)]

    # Calculate an average I(q)
    q = [mtd["ws_rebin_{}_iq".format(index)].extractX().ravel() for
         index, _ in enumerate(bins)]
    i = [mtd["ws_rebin_{}_iq".format(index)].extractY().ravel() for
         index, _ in enumerate(bins)]
    dq = [mtd["ws_rebin_{}_iq".format(index)].extractE().ravel() for
          index, _ in enumerate(bins)]

    q = np.ma.masked_invalid(q)
    i = np.ma.masked_invalid(i)
    dq = np.ma.masked_invalid(dq)

    # Calculate an average I(q)
    mean_q = np.mean(q, axis=0)
    mean_i = np.mean(i, axis=0)
    mean_dq = np.sqrt(np.sum(dq**2, axis=0))

    mean_iq = CreateWorkspace(
        DataX=mean_q.filled(np.nan),
        DataY=mean_i.filled(np.nan),
        DataE=mean_dq.filled(np.nan),
        # Dx=dq_bin_centers,  # bin centers!!
        NSpec=1,
        UnitX='MomentumTransfer',
        YUnitLabel='Counts',
        OutputWorkspace="mean_rebin_iq"
    )

    iq_old = legacy_reduction()

    iq_new = ReplaceSpecialValues(
        InputWorkspace=mean_iq, NaNValue=0, InfinityValue=0)

    max_old = iq_old.readY(0).max()
    max_new = iq_new.readY(0).max()

    iq_old_final = iq_old * max_new/max_old

    assert iq_old_final is not None

    # plt.plot(iq_new.extractY().ravel())
    # plt.plot(iq_old_final.extractY().ravel())
    # plt.show()
