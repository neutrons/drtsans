#!/usr/bin/env python
from __future__ import print_function

import numpy as np

import mantid
from mantid import mtd
from mantid.simpleapi import (AddSampleLog, ConfigService, ExtractSpectra,
                              Rebin)
from ornl.sans.sns.eqsans import (load_events, normalisation,
                                  transform_to_wavelength)
from ornl.sans.sns.eqsans.momentum_transfer import MomentumTransfer

from ornl.sans.sns.eqsans import center_detector
from ornl.sans.sns.eqsans import geometry


def legacy_reduction():

    from reduction_workflow.command_interface import AppendDataFile, Reduce
    from reduction_workflow.instruments.sans import (sns_command_interface as
                                                     eqsans)
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

    ws = load_events('EQSANS_68200', detector_offset=0, sample_offset=0)

    ws = transform_to_wavelength(ws,
                                 bin_width=0.1,
                                 low_tof_clip=500,
                                 high_tof_clip=2000)

    center_detector(ws, x=-0.025, y=-0.016, units='m')

    flux_ws = normalisation.load_beam_flux_file(
        '/SNS/EQSANS/shared/instrument_configuration/bl6_flux_at_sample',
        output_workspace='flux_ws',
        ws_reference=ws)

    ws = normalisation.normalise_by_proton_charge_and_flux(ws, flux_ws, "ws")

    # geometry.sample_aperture_diameter is not working: slit4 missing
    # We hard code the sample_aperture_diameter instead
    AddSampleLog(Workspace=ws,
                 LogName='sample-aperture-diameter',
                 LogText='10.',
                 LogType='Number',
                 LogUnit='mm')

    geometry.source_aperture_diameter(ws)

    rebin_start, rebin_end, rebin_step = 2.6, 5.6, 0.2

    ws = Rebin(InputWorkspace=ws,
               OutputWorkspace="ws_rebin",
               Params="{:.2f},{:.2f},{:.2f}".format(rebin_start, rebin_step,
                                                    rebin_end))

    bins = np.arange(rebin_start, rebin_end, rebin_step)

    mt_sum = MomentumTransfer()

    for index, bin_start in enumerate(bins):
        total_pixels_in_detector = 256 * 192
        ws_extracted = ExtractSpectra(InputWorkspace=ws,
                                      XMin=bin_start,
                                      XMax=bin_start + rebin_step)

        wavelength_mean = bin_start + rebin_step / 2
        AddSampleLog(Workspace=ws_extracted,
                     LogName='wavelength',
                     LogText="{:.2f}".format(wavelength_mean),
                     LogType='Number',
                     LogUnit='Angstrom')
        AddSampleLog(Workspace=ws_extracted,
                     LogName='wavelength-spread',
                     LogText='0.2',
                     LogType='Number',
                     LogUnit='Angstrom')

        mt_extracted = MomentumTransfer(ws_extracted)
        mt_sum += mt_extracted

        assert mt_extracted.qx.shape == mt_extracted.qy.shape == \
            mt_extracted.dqx.shape == mt_extracted.dqy.shape == \
            mt_extracted.i.shape == mt_extracted.i_sigma.shape == \
            (total_pixels_in_detector, )

        assert mt_sum.qx.shape == mt_sum.qy.shape == mt_sum.dqx.shape == \
            mt_sum.dqy.shape == mt_sum.i.shape == mt_sum.i_sigma.shape == \
            (total_pixels_in_detector + (index * total_pixels_in_detector),)

    ws_sum_table = mt_sum.q2d()
    assert isinstance(ws_sum_table, mantid.dataobjects.TableWorkspace)

    _, ws_sum_q2d = mt_sum.bin_into_q2d()
    assert ws_sum_q2d.extractY().shape == (256, 192)
    assert ws_sum_q2d.extractX().shape == (256, 193)

    _, ws_sum_q1d = mt_sum.bin_into_q1d()
    assert ws_sum_q1d.extractY().shape == (1, 100)
    assert ws_sum_q1d.extractX().shape == (1, 101)

    ws_iq_old = legacy_reduction()
    ws_iq_new = ws_sum_q1d

    max_old = ws_iq_old.readY(0).max()
    max_new = ws_iq_new.readY(0).max()

    ws_iq_old = ws_iq_old * max_new / max_old

    ws_iq_new = Rebin(InputWorkspace=ws_iq_new,
                      Params="0.00119358, .00463923, 0.457006261051")
    ws_iq_old = Rebin(InputWorkspace=ws_iq_old,
                      Params="0.00119358, .00463923, 0.457006261051")

    assert np.allclose(ws_iq_new.extractY(), ws_iq_old.extractY(), rtol=1)

    # plt.loglog(ws_iq_new.extractY().ravel(), label = 'New')
    # plt.loglog(ws_iq_old.extractY().ravel(), label = 'Old')
    # plt.legend()
    # plt.show()
