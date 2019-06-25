#!/usr/bin/env python
from __future__ import print_function

import numpy as np

from mantid import mtd
from mantid.simpleapi import (AddSampleLog, ConfigService, CreateWorkspace,
                              ExtractSpectra, MoveInstrumentComponent, Rebin)
from ornl.sans.momentum_transfer import bin_into_q1d, bin_into_q2d
from ornl.sans.sns.eqsans import reduce
from ornl.settings import unique_workspace_name


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
    eqsans.SetTransmission(1.0, 1.0)
    eqsans.SetTOFTailsCutoff(low_cut=500, high_cut=2000)
    eqsans.SetWavelengthStep(step=0.1)
    # eqsans.NoTransmission()
    eqsans.NoSolidAngle()
    eqsans.OutputPath(tempfile.gettempdir())
    Reduce()

    return mtd['EQSANS_68200']


def test_momentum_tranfer_serial():
    '''

    '''

    ws = reduce.load_w('EQSANS_68200', low_tof_clip=500,
                       high_tof_clip=2000, dw=0.1)
    # Center the beam: xyz: 0.0165214,0.0150392,4.00951
    MoveInstrumentComponent(Workspace='ws', ComponentName='detector1',
                            X=-0.025, Y=-0.016, RelativePosition=False)

    # this temporary. Waiting for jose to finish `eqsans.load_events`
    # https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/issues/66
    AddSampleLog(Workspace=ws, LogName='sample-aperture-diameter',
                 LogText='20.', LogType='Number', LogUnit='mm')
    AddSampleLog(Workspace=ws, LogName='source-aperture-diameter',
                 LogText='20.', LogType='Number', LogUnit='mm')

    rebin_start, rebin_end, rebin_step = 2.6, 5.6, 0.2

    ws = Rebin(
        InputWorkspace=ws,
        OutputWorkspace=unique_workspace_name(suffix="_rebin"),
        Params="{:.2f},{:.2f},{:.2f}".format(
            rebin_start, rebin_step, rebin_end)
    )

    bins = np.arange(rebin_start, rebin_end, rebin_step)

    for index, bin_start in enumerate(bins):

        ws_extracted = ExtractSpectra(
            InputWorkspace=ws,
            OutputWorkspace=unique_workspace_name(suffix="_{:}".format(index)),
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
