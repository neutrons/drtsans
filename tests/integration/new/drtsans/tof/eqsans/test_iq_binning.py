import os

# import numpy as np

# import mantid
from mantid import mtd
from mantid.simpleapi import AddSampleLog, ConfigService, Rebin  # ExtractSpectra MaskAngle,
from drtsans.tof.eqsans import (center_detector, geometry, load_events, normalization, transform_to_wavelength)
from drtsans.iq import bin_iq_into_linear_q1d, BinningMethod
from drtsans.convert_to_q import convert_to_q


# Integration test on I(Q) binning algorithms for EQ-SANS


def legacy_reduction(reference_dir):

    from reduction_workflow.command_interface import AppendDataFile, Reduce
    from reduction_workflow.instruments.sans import (sns_command_interface as
                                                     eqsans)
    import tempfile

    configI = ConfigService.Instance()
    configI["facilityName"] = 'SNS'

    eqsans.EQSANS()
    AppendDataFile(os.path.join(reference_dir.new.eqsans, 'EQSANS_68200_event.nxs'))
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


def test_iq_binning_serial(reference_dir):
    """
    Integration test workflow for binning I(Q)
    Parameters
    ----------
    reference_dir

    Returns
    -------

    """
    # Load event data
    ws = load_events(os.path.join(reference_dir.new.eqsans, 'EQSANS_68200_event.nxs'),
                     detector_offset=0,
                     sample_offset=0)

    # Convert to wave length
    ws = transform_to_wavelength(ws, bin_width=0.1, low_tof_clip=500, high_tof_clip=2000)

    # Calibration in the next few steps
    center_detector(ws, center_x=0.025, center_y=0.016)

    flux_ws = normalization.load_beam_flux_file(os.path.join(
        reference_dir.new.eqsans, 'test_normalization', 'beam_profile_flux.txt'),
        output_workspace='flux_ws', data_workspace=ws)

    ws = normalization.normalize_by_proton_charge_and_flux(ws, flux_ws, "ws")

    # Prepare to calculate Q, dQ and bin I(Q)
    # NOTE: geometry.sample_aperture_diameter is not working: slit4 missing in EQSANS_68200_event.nxs
    # We hard code the sample_aperture_diameter instead
    AddSampleLog(Workspace=ws,
                 LogName='sample-aperture-diameter',
                 LogText='10.',
                 LogType='Number',
                 LogUnit='mm')
    geometry.source_aperture_diameter(ws)

    # Rebin the workspace: rebin to a coarse binning of wave length matching to final resolution I(Q)
    rebin_start, rebin_end, rebin_step = 2.6, 5.6, 0.2
    ws = Rebin(InputWorkspace=ws,
               OutputWorkspace="ws_rebin",
               Params="{:.2f},{:.2f},{:.2f}".format(rebin_start, rebin_step,
                                                    rebin_end))

    # Calculate Q and dQ
    r = convert_to_q(ws, mode='scalar', resolution_function=None)
    iq_array = r[0]
    sigma_iq_array = r[1]
    q_array = r[2]
    dq_array = r[3]

    # Bin I(Q)
    final_q_min = 0
    i_of_q = bin_intensity_into_q1d(iq_array, sigma_iq_array, q_array, dq_array, bins=10, q_min=final_q_min,
                                    linear_binning=True, bin_method=BinningMethod.WEIGHTED)
    assert i_of_q

    # TODO - continue from here
    # assert i_of_q.iq.shape == (256, 192)
    # assert i_of_q.q.shape == (256, 193)

    # bins = np.arange(rebin_start, rebin_end, rebin_step)

    # TODO - The following test shall be removed with new API
    # mt_sum = sfer()

    # for index, bin_start in enumerate(bins):
    #     total_pixels_in_detector = 256 * 192
    #     ws_extracted = ExtractSpectra(InputWorkspace=ws,
    #                                   XMin=bin_start,
    #                                   XMax=bin_start + rebin_step)

    #     wavelength_mean = bin_start + rebin_step / 2
    #     AddSampleLog(Workspace=ws_extracted,
    #                  LogName='wavelength',
    #                  LogText="{:.2f}".format(wavelength_mean),
    #                  LogType='Number',
    #                  LogUnit='Angstrom')
    #     AddSampleLog(Workspace=ws_extracted,
    #                  LogName='wavelength-spread',
    #                  LogText='0.2',
    #                  LogType='Number',
    #                  LogUnit='Angstrom')

    #     mt_extracted = MomentumTransfer(ws_extracted)
    #     mt_sum += mt_extracted

    #     assert mt_extracted.qx.shape == mt_extracted.qy.shape == \
    #         mt_extracted.dqx.shape == mt_extracted.dqy.shape == \
    #         mt_extracted.i.shape == mt_extracted.i_sigma.shape == \
    #         (total_pixels_in_detector, )

    #     assert mt_sum.qx.shape == mt_sum.qy.shape == mt_sum.dqx.shape == \
    #         mt_sum.dqy.shape == mt_sum.i.shape == mt_sum.i_sigma.shape == \
    #         (total_pixels_in_detector + (index * total_pixels_in_detector),)

    # ws_sum_table = mt_sum.q2d()
    # assert isinstance(ws_sum_table, mantid.dataobjects.TableWorkspace)

    # _, ws_sum_q2d = mt_sum.bin_into_q2d()
    # assert ws_sum_q2d.extractY().shape == (256, 192)
    # assert ws_sum_q2d.extractX().shape == (256, 193)

    # _, ws_sum_q1d = mt_sum.bin_into_q1d()
    # assert ws_sum_q1d.extractY().shape == (1, 100)
    # assert ws_sum_q1d.extractX().shape == (1, 101)

    # ws_iq_old = legacy_reduction(reference_dir)
    # ws_iq_new = ws_sum_q1d

    # max_old = ws_iq_old.readY(0).max()
    # max_new = ws_iq_new.readY(0).max()

    # ws_iq_old = ws_iq_old * max_new / max_old

    # ws_iq_new = Rebin(InputWorkspace=ws_iq_new,
    #                   Params="0.00119358, .00463923, 0.457006261051")
    # ws_iq_old = Rebin(InputWorkspace=ws_iq_old,
    #                   Params="0.00119358, .00463923, 0.457006261051")

    # assert np.allclose(ws_iq_new.extractY(), ws_iq_old.extractY(), rtol=1)

    # plt.loglog(ws_iq_new.extractY().ravel(), label = 'New')
    # plt.loglog(ws_iq_old.extractY().ravel(), label = 'Old')
    # plt.legend()
    # plt.show()


def skip_test_api(reference_dir):

    ws = load_events(os.path.join(reference_dir.new.eqsans, 'EQSANS_68200_event.nxs'),
                     detector_offset=0,
                     sample_offset=0)
    assert ws

    # ws = transform_to_wavelength(ws,
    #                              bin_width=0.1,
    #                              low_tof_clip=500,
    #                              high_tof_clip=2000)

    # center_detector(ws, x=-0.025, y=-0.016, unit='m')

    # flux_ws = normalization.load_beam_flux_file(os.path.join(
    #     reference_dir.new.eqsans, 'test_normalization', 'beam_profile_flux.txt'),
    #     output_workspace='flux_ws',
    #     ws_reference=ws)

    # ws = normalization.normalize_by_proton_charge_and_flux(ws, flux_ws, "ws")

    # #
    # table_ws = prepare_momentum_transfer(
    #     ws, wavelength_binning=[2.6, 0.2, 5.6])
    # assert mtd.doesExist(ws.name() + "_table")

    # table_ws = table_ws[0]

    # assert mtd.doesExist(ws.name() + "_iq")

    # iq_ws = cal_iq(table_ws, bins=100, log_binning=True)
    # assert mtd[ws.name() + "_iq"] is not None

    # # check if it is log binning:
    # q = iq_ws.readX(0)
    # q_divide = q[1:] / q[:-1]
    # assert np.allclose(q_divide, q_divide[0])

    # # TODO: this need some Anisotopic data to tes. See gpsans test.
    # # TODO: Ask BL Scientists.
    # # Test Wedge
    # iq_wedge_ws = iq_wedge(table_ws)
    # assert iq_wedge_ws

    # # Tests annulus
    # iq_annular_ws_1 = iq_annular(table_ws, q_min=0.0001, q_max=0.1,
    #                              bins=20, suffix="_annular_iq_1")
    # assert iq_annular_ws_1
    # y_ws_1 = iq_annular_ws_1.extractY().ravel()
    # iq_annular_ws_2 = iq_annular(table_ws, q_min=0.1, q_max=0.3,
    #                              bins=40, suffix="_annular_iq_2")
    # assert iq_annular_ws_2
    # y_ws_2 = iq_annular_ws_2.extractY().ravel()
    # assert abs(np.log10(y_ws_1.sum())-np.log10(y_ws_2.sum())) > 2  # orders magnitude

    # iqxqy_ws = iqxqy(table_ws)
    # assert iqxqy_ws is not None
    # assert mtd.doesExist(ws.name() + "_iqxqy")

    # # test 2d with log_binning=True
    # iqxqy_ws_log = iqxqy(table_ws, log_binning=True)
    # qx = iqxqy_ws_log.readX(0)
    # # Q log bining in 2D may have different negative and positive bin widths.
    # # The Qx must be thus divided in two: the beggining (negative) and at the
    # # end (positive). Let's pick a few bins (e.g. 30) where we now the division
    # # with the anterior bin is constant.
    # qx_divide_negative = qx[1:31] / qx[:30]  # beginning off Qx
    # qx_divide_positive = qx[-30:] / qx[-31:-1]  # end of qx
    # assert np.allclose(qx_divide_negative, qx_divide_negative[0])
    # assert np.allclose(qx_divide_positive, qx_divide_positive[0])

    # # Test Mask now
    # MaskAngle(Workspace=ws, MaxAngle=0.4, Angle="TwoTheta")

    # table_ws_masked = prepare_momentum_transfer(
    #     ws, wavelength_binning=[2.6, 0.2, 5.6], prefix='mask')

    # table_ws_masked = table_ws_masked[0]

    # iq_ws_masked = cal_iq(table_ws_masked, bins=100, log_binning=True)
    # assert np.all(
    #     iq_ws_masked.extractY()[0][:5] < iq_ws.extractY()[0][:5])

    # iqxqy_ws_masked = iqxqy(table_ws_masked)
    # # Masked values around beam center
    # assert iqxqy_ws_masked.readY(51)[46] == 0
    # assert iqxqy_ws_masked.readY(51)[47] == 0


def skip_test_api_frame_skipping(reference_dir):

    db_ws = load_events(os.path.join(reference_dir.new.eqsans, "EQSANS_88973.nxs.h5"))
    assert db_ws

    # center = center_detector(db_ws)

    # ws = prepare_data(os.path.join(reference_dir.new.eqsans, "EQSANS_88980.nxs.h5"),
    #                   x_center=-center[0], y_center=-center[1],
    #                   # use_config_tof_cuts=True,
    #                   # use_config=True,
    #                   # use_config_mask=True,
    #                   # correct_for_flight_path=True,
    #                   # flux=beam_flux_file,
    #                   # mask=detector_ids604m,
    #                   # dark_current=dark_data_file,
    #                   # sensitivity_file_path=sensitivity_file,
    #                   sample_offset=340)

    # ws_frame1, ws_frame2 = prepare_momentum_transfer(ws)

    # ws_frame1_iq = cal_iq(ws_frame1, bins=150, log_binning=True)
    # ws_frame2_iq = cal_iq(ws_frame2, bins=150, log_binning=True)
    # assert ws_frame1_iq is not None
    # assert ws_frame2_iq is not None
    # assert not np.allclose(ws_frame1_iq.extractX(), ws_frame2_iq.extractX())

    # ws_frame1_iqxqy = iqxqy(ws_frame1)
    # ws_frame2_iqxqy = iqxqy(ws_frame2)
    # assert ws_frame1_iqxqy is not None
    # assert ws_frame2_iqxqy is not None
