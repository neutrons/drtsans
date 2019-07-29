#!/usr/bin/env python
from __future__ import print_function

import os

import numpy as np

import mantid
from mantid import mtd
from mantid.simpleapi import AddSampleLog, ConfigService, ExtractSpectra, Rebin
from ornl.sans.sns.eqsans import (center_detector, geometry, load_events,
                                  normalisation, transform_to_wavelength,
                                  prepare_data)
from ornl.sans.sns.eqsans.momentum_transfer import (MomentumTransfer, iq,
                                                    iqxqy,
                                                    prepare_momentum_transfer)


def legacy_reduction(refd):

    from reduction_workflow.command_interface import AppendDataFile, Reduce
    from reduction_workflow.instruments.sans import (sns_command_interface as
                                                     eqsans)
    import tempfile

    configI = ConfigService.Instance()
    configI["facilityName"] = 'SNS'

    eqsans.EQSANS()
    AppendDataFile(os.path.join(refd.new.eqsans, 'EQSANS_68200_event.nxs'))
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

    return mtd['EQSANS_68200_event_iq']


def test_momentum_tranfer_serial(refd):

    ws = load_events(os.path.join(refd.new.eqsans, 'EQSANS_68200_event.nxs'),
                     detector_offset=0,
                     sample_offset=0)

    ws = transform_to_wavelength(ws,
                                 bin_width=0.1,
                                 low_tof_clip=500,
                                 high_tof_clip=2000)

    center_detector(ws, x=-0.025, y=-0.016, unit='m')

    flux_ws = normalisation.load_beam_flux_file(os.path.join(
        refd.new.eqsans, 'test_normalisation', 'beam_profile_flux.txt'),
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

    ws_iq_old = legacy_reduction(refd)
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


def test_api(refd):

    ws = load_events(os.path.join(refd.new.eqsans, 'EQSANS_68200_event.nxs'),
                     detector_offset=0,
                     sample_offset=0)

    ws = transform_to_wavelength(ws,
                                 bin_width=0.1,
                                 low_tof_clip=500,
                                 high_tof_clip=2000)

    center_detector(ws, x=-0.025, y=-0.016, unit='m')

    flux_ws = normalisation.load_beam_flux_file(os.path.join(
        refd.new.eqsans, 'test_normalisation', 'beam_profile_flux.txt'),
        output_workspace='flux_ws',
        ws_reference=ws)

    ws = normalisation.normalise_by_proton_charge_and_flux(ws, flux_ws, "ws")

    #
    table_ws = prepare_momentum_transfer(
        ws, wavelength_binning=[2.6, 0.2, 5.6])
    assert mtd.doesExist(ws.name() + "_table")

    table_ws = table_ws[0]

    iq_ws = iq(table_ws)
    assert mtd.doesExist(ws.name() + "_iq")

    _ = iq(table_ws, bins=100, log_binning=True)
    iq_ws = mtd[ws.name() + "_iq"]
    # check if it is log binning:
    q = iq_ws.readX(0)
    q_divide = q[1:] / q[:-1]
    assert np.allclose(q_divide, q_divide[0])

    _ = iqxqy(table_ws)
    assert mtd.doesExist(ws.name() + "_iqxqy")


def test_api_frame_skipping(refd):

    db_ws = load_events(os.path.join(refd.new.eqsans, "EQSANS_88973.nxs.h5"))
    center = center_detector(db_ws)

    ws = prepare_data(os.path.join(refd.new.eqsans, "EQSANS_88980.nxs.h5"),
                      x_center=-center[0], y_center=-center[1],
                      # use_config_tof_cuts=True,
                      # use_config=True,
                      # use_config_mask=True,
                      # correct_for_flight_path=True,
                      # flux=beam_flux_file,
                      # mask=detector_ids604m,
                      # dark_current=dark_data_file,
                      # sensitivity_file_path=sensitivity_file,
                      sample_offset=340)

    ws_frame1, ws_frame2 = prepare_momentum_transfer(ws)

    ws_frame1_iq = iq(ws_frame1, bins=150, log_binning=True)
    ws_frame2_iq = iq(ws_frame2, bins=150, log_binning=True)

    _ = iqxqy(ws_frame1)
    _ = iqxqy(ws_frame2)

    assert not np.allclose(ws_frame1_iq.extractX(), ws_frame2_iq.extractX())

    assert ws_frame1_iq is not None
    assert ws_frame2_iq is not None
