from __future__ import absolute_import, division, print_function

import tempfile

import pytest


'''
For every flood:
    LoadData
    Prepare Flood:
        - Geometry
        - SA
        - DC subtraction
        - Monitor Normalization
        - Transmission

    sensitivity = Calculate sensitivity

    Mask beamstop

Join (average) all the sensitivities in one single file
Save file as nexus

'''


@pytest.mark.offline
def test_sensitivity_procedural(biosans_sensitivity_dataset):

    from mantid.simpleapi import (CalculateSensitivity, LoadHFIRSANS, LoadMask,
                                  MaskDetectors, MoveInstrumentComponent,
                                  RenameWorkspace, ReplaceSpecialValues,
                                  SANSMaskDTP, SANSSolidAngle, SaveNexus)
    from ornl.sans.hfir.biosans.beam_finder import direct_beam_center
    from ornl.sans.hfir.dark_current import subtract_normalised_dark
    from ornl.sans.hfir.normalisation import time
    from ornl.sans.sensitivity import inf_value_to_mask, interpolate_mask
    from ornl.sans.transmission import (apply_transmission_correction_value,
                                        calculate_transmission)
    # Load the files into WS
    dark_current_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['dark_current'])
    flood_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood'])
    flood_beamcenter_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood_beamcenter'])
    empty_transmission_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['empty_transmission'])
    flood_mask_ws = LoadMask(
        Instrument='BioSANS',
        InputFile=biosans_sensitivity_dataset['flood_mask'])

    # This is for the main detector, let's mask the wing detector
    SANSMaskDTP(InputWorkspace=dark_current_ws, Detector="2")
    SANSMaskDTP(InputWorkspace=flood_ws, Detector="2")
    SANSMaskDTP(InputWorkspace=flood_beamcenter_ws, Detector="2")
    SANSMaskDTP(InputWorkspace=empty_transmission_ws, Detector="2")

    # Let's correct the data first
    ###########################################################################
    # DC normalisation
    dark_current_norm_ws = time(dark_current_ws)
    # DC Subtraction
    flood_dc_corrected_ws = subtract_normalised_dark(
        flood_ws, dark_current_norm_ws)
    flood_dc_corrected_ws = RenameWorkspace(
        InputWorkspace=flood_dc_corrected_ws,
        OutputWorkspace="flood_dc_corrected_ws")

    empty_transmission_dc_corrected_ws = subtract_normalised_dark(
        empty_transmission_ws, dark_current_norm_ws)
    empty_transmission_dc_corrected_ws = RenameWorkspace(
        InputWorkspace=empty_transmission_dc_corrected_ws,
        OutputWorkspace="empty_transmission_dc_corrected_ws")

    flood_beamcenter_dc_corrected_ws = subtract_normalised_dark(
        flood_beamcenter_ws, dark_current_norm_ws)
    flood_beamcenter_dc_corrected_ws = RenameWorkspace(
        InputWorkspace=flood_beamcenter_dc_corrected_ws,
        OutputWorkspace="flood_beamcenter_dc_corrected_ws")

    ###########################################################################
    # Find the beam center
    x, y, y_gravity = direct_beam_center(
        flood_beamcenter_dc_corrected_ws, tolerance=0.00125, direct_beam=True,
        beam_radius=0.0155)
    assert x == pytest.approx(0.011, abs=1e-3)
    assert y == pytest.approx(-0.012, abs=1e-3)

    ###########################################################################
    # Get the right geometry
    MoveInstrumentComponent(
        Workspace=flood_dc_corrected_ws, ComponentName='detector1', X=-x, Y=-y)
    MoveInstrumentComponent(
        Workspace=empty_transmission_dc_corrected_ws,
        ComponentName='detector1', X=-x, Y=-y)

    ###########################################################################
    # Normalization (In the original script they use time normalization)
    flood_dc_time_corrected_ws = time(flood_dc_corrected_ws)
    flood_dc_time_corrected_ws = RenameWorkspace(
        InputWorkspace=flood_dc_time_corrected_ws,
        OutputWorkspace="flood_dc_time_corrected_ws")

    empty_transmission_dc_time_corrected_ws = time(
        empty_transmission_dc_corrected_ws)
    empty_transmission_dc_time_corrected_ws = RenameWorkspace(
        InputWorkspace=empty_transmission_dc_time_corrected_ws,
        OutputWorkspace="empty_transmission_dc_time_corrected_ws")

    ###########################################################################
    # Solid Angle correction
    solid_angle_ws = SANSSolidAngle(
        InputWorkspace=flood_dc_time_corrected_ws, Type='Tube')

    flood_dc_time_sa_corrected_ws = flood_dc_time_corrected_ws / solid_angle_ws

    flood_dc_time_sa_corrected_ws = ReplaceSpecialValues(
        InputWorkspace=flood_dc_time_sa_corrected_ws,
        NaNValue=0, InfinityValue=0)

    # No need to recalculate solid angle.
    # Probably we don't need it for the transmission.
    empty_transmission_time_sa_corrected_ws =\
        empty_transmission_dc_time_corrected_ws / solid_angle_ws

    empty_transmission_time_sa_corrected_ws = ReplaceSpecialValues(
        InputWorkspace=empty_transmission_time_sa_corrected_ws,
        NaNValue=0, InfinityValue=0)

    ###########################################################################
    # This is only to get transmission from the flat measurement
    # The value will be used in the wing detector
    calculated_transmission_value,  calculated_transmission_error = \
        calculate_transmission(flood_dc_time_sa_corrected_ws,
                               empty_transmission_time_sa_corrected_ws)

    ###########################################################################
    # Sensitivity calculation
    sensitivity_ws = CalculateSensitivity(
        InputWorkspace=flood_dc_time_sa_corrected_ws, MinThreshold=0.5,
        MaxThreshold=1.5)

    ###########################################################################
    # Load and mask sensitivity according to the beamstop
    MaskDetectors(Workspace=sensitivity_ws, MaskedWorkspace=flood_mask_ws)

    ###########################################################################
    # Let's interpolate the masks
    sensitivity_interpolated_ws = interpolate_mask(
        sensitivity_ws, polynomial_degree=1)
    sensitivity_interpolated_ws = RenameWorkspace(
        InputWorkspace=sensitivity_interpolated_ws,
        OutputWorkspace="sensitivity_interpolated_ws")

    # assert sensitivity_interpolated_ws.readY(21640)[0] == Property.EMPTY_DBL
    # inf_value_to_mask(sensitivity_interpolated_ws)

    assert sensitivity_interpolated_ws.readY(
        21640)[0] == pytest.approx(1.0, abs=1e-03)
    assert sensitivity_interpolated_ws.readE(
        21640)[0] == pytest.approx(0.0, abs=1e-03)
    assert sensitivity_interpolated_ws.detectorInfo().isMasked(21640)

    ###########################################################################
    # Save

    with tempfile.NamedTemporaryFile(
            delete=True, prefix="sensitivivity_cg3_main_",
            suffix=".nxs") as temp_file:
        SaveNexus(InputWorkspace=sensitivity_interpolated_ws,
                  Filename=temp_file.name,
                  Title='CG3 sensitivivity')
        print("Saved NeXus sensitivity file to: {}".format(temp_file.name))

###############################################################################
# Let's do the wing now!
###############################################################################
    # Because we used Masks the wing detector is Zeroed. Load again!

    # Load the files into WS
    dark_current_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['dark_current'])
    flood_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['flood'])

    # Let's mask the main detector
    SANSMaskDTP(InputWorkspace=dark_current_ws, Detector="1")
    SANSMaskDTP(InputWorkspace=flood_ws, Detector="1")

    # DC normalisation
    dark_current_norm_ws = time(dark_current_ws)
    # DC Subtraction
    flood_dc_corrected_ws = subtract_normalised_dark(
        flood_ws, dark_current_norm_ws)
    flood_dc_corrected_ws = RenameWorkspace(
        InputWorkspace=flood_dc_corrected_ws,
        OutputWorkspace="flood_dc_corrected_ws")

    # Geometry with gravity correction
    MoveInstrumentComponent(
        Workspace=flood_dc_corrected_ws,
        ComponentName='detector1', X=-x, Y=-y_gravity)

    # Normalization (In the original script they use time normalization)
    flood_dc_time_corrected_ws = time(flood_dc_corrected_ws)
    flood_dc_time_corrected_ws = RenameWorkspace(
        InputWorkspace=flood_dc_time_corrected_ws,
        OutputWorkspace="flood_dc_time_corrected_ws")

    # No solid angle correction

    # Apply Transmission Correction
    flood_dc_time_trans_corrected_wing_ws = \
        apply_transmission_correction_value(
            flood_dc_time_corrected_ws, calculated_transmission_value,
            calculated_transmission_error, theta_dependent=False)
    flood_dc_time_trans_corrected_wing_ws = RenameWorkspace(
        InputWorkspace=flood_dc_time_trans_corrected_wing_ws,
        OutputWorkspace="flood_dc_time_trans_corrected_wing_ws")

    # Sensitivity calculation
    sensitivity_ws = CalculateSensitivity(
        InputWorkspace=flood_dc_time_trans_corrected_wing_ws, MinThreshold=0.5,
        MaxThreshold=1.5)

    inf_value_to_mask(sensitivity_ws)
    assert sensitivity_ws.readY(52867)[0] == 1
    assert sensitivity_ws.readE(52867)[0] == 0
    assert sensitivity_ws.detectorInfo().isMasked(52867)


@pytest.mark.offline
def test_sensitivity_detector(biosans_sensitivity_dataset):

    from ornl.sans.sensitivity import Detector
    from mantid.simpleapi import LoadHFIRSANS
    import numpy as np

    # Load the files into WS
    dark_current_ws = LoadHFIRSANS(
        Filename=biosans_sensitivity_dataset['dark_current'])

    d = Detector(dark_current_ws, "detector1")
    assert d.first_det_id == 3
    assert d.last_det_id == (
        d.n_pixels_per_tube*d.n_tubes + d.first_det_id - 1)
    assert d.detector_id_to_ws_index[d.first_det_id] == 2
    
    # Get WS indices for the first tube
    start_ws_index, stop_ws_index = next(d.tube_ws_indices)
    assert start_ws_index == d.detector_id_to_ws_index[d.first_det_id]
    assert stop_ws_index == d.detector_id_to_ws_index[
        d.first_det_id+d.n_pixels_per_tube]

    #
    y, e = d.ws_values(start_ws_index, stop_ws_index)
    assert len(y) == d.n_pixels_per_tube
    assert len(e) == d.n_pixels_per_tube

    start_ws_index, stop_ws_index = next(d.tube_ws_indices)
    assert start_ws_index == d.detector_id_to_ws_index[
        d.first_det_id+d.n_pixels_per_tube]
    y2, _ = d.ws_values(start_ws_index, stop_ws_index)
    assert np.array_equal(y, y)
    assert not np.array_equal(y, y2)
    