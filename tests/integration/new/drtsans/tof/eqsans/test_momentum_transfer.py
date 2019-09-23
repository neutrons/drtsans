import os

import numpy as np

# import mantid
from mantid import mtd
from mantid.simpleapi import AddSampleLog, ConfigService, MaskAngle, Rebin  # ExtractSpectra
from drtsans.tof.eqsans import (center_detector, geometry, load_events, normalisation, prepare_data,
                                transform_to_wavelength)


# This integration test is to verify the workflow to calculate Q and dQ
# q_resolution_per_pixel() is tested in test_resolution
# DEV - Wenduo Zhou <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>


def test_momentum_transfer_serial(reference_dir):

    print(mtd)
    print(mantid.simpleapi.ConfigService)
    print(mantid.simpleapi.MaskAngle)
    print(drtsans.tof.eqsans.prepare_data)

    ws = load_events(os.path.join(reference_dir.new.eqsans, 'EQSANS_68200_event.nxs'),
                     detector_offset=0,
                     sample_offset=0)

    ws = transform_to_wavelength(ws,
                                 bin_width=0.1,
                                 low_tof_clip=500,
                                 high_tof_clip=2000)

    center_detector(ws, x=-0.025, y=-0.016, unit='m')

    flux_ws = normalisation.load_beam_flux_file(os.path.join(
        reference_dir.new.eqsans, 'test_normalisation', 'beam_profile_flux.txt'),
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
    print(bins)

    # TODO - calculate q and dq

    # Check: total_pixels_in_detector = 256 * 192 with number of wave length

    return

