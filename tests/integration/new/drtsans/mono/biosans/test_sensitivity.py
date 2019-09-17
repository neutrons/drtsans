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

    from mantid.simpleapi import (CalculateEfficiency, LoadHFIRSANS, LoadMask,
                                  MaskDetectors, MoveInstrumentComponent,
                                  RenameWorkspace, ReplaceSpecialValues,
                                  MaskBTP, SolidAngle, SaveNexus)
    from drtsans.mono.biosans.beam_finder import find_beam_center
    from drtsans.mono.dark_current import subtract_normalised_dark
    from drtsans.mono.normalisation import time
    from drtsans.sensitivity import inf_value_to_mask, interpolate_mask
    from drtsans.transmission import (apply_transmission_correction, calculate_transmission)
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
    MaskBTP(Workspace=dark_current_ws, Components='wing_detector')
    MaskBTP(Workspace=flood_ws, Components='wing_detector')
    MaskBTP(Workspace=flood_beamcenter_ws, Components='wing_detector')
    MaskBTP(Workspace=empty_transmission_ws, Components='wing_detector')

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
    x, y, y_gravity = find_beam_center(
        flood_beamcenter_dc_corrected_ws, Tolerance=0.00125, DirectBeam=True,
        BeamRadius=0.0155)
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
    solid_angle_ws = SolidAngle(
        InputWorkspace=flood_dc_time_corrected_ws, Method='VerticalTube')

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
    trans = calculate_transmission(flood_dc_time_sa_corrected_ws,
                                   empty_transmission_time_sa_corrected_ws)
    calculated_transmission_value = trans.dataY(0)[0]
    calculated_transmission_error = trans.dataE(0)[0]
    ###########################################################################
    # Sensitivity calculation
    sensitivity_ws = CalculateEfficiency(
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
    MaskBTP(Workspace=dark_current_ws, Components='detector1')
    MaskBTP(Workspace=flood_ws, Components="detector1")

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
        apply_transmission_correction(
            flood_dc_time_corrected_ws,
            trans_value=calculated_transmission_value,
            trans_error=calculated_transmission_error,
            theta_dependent=False)
    flood_dc_time_trans_corrected_wing_ws = RenameWorkspace(
        InputWorkspace=flood_dc_time_trans_corrected_wing_ws,
        OutputWorkspace="flood_dc_time_trans_corrected_wing_ws")

    # Sensitivity calculation
    sensitivity_ws = CalculateEfficiency(
        InputWorkspace=flood_dc_time_trans_corrected_wing_ws, MinThreshold=0.5,
        MaxThreshold=1.5)

    inf_value_to_mask(sensitivity_ws)
    assert sensitivity_ws.readY(52867)[0] == 1
    assert sensitivity_ws.readE(52867)[0] == 0
    assert sensitivity_ws.detectorInfo().isMasked(52867)


@pytest.mark.offline
def test_sensitivity_detector(biosans_sensitivity_dataset):
    '''This tests the drtsans.sensitivity.Detector with data from BioSANS
    '''

    from drtsans.sensitivity import Detector
    from mantid.simpleapi import LoadHFIRSANS, MaskBTP
    from mantid.kernel import Property
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
    d.next_tube()
    start_ws_index, stop_ws_index = d.get_current_ws_indices()
    assert start_ws_index == d.detector_id_to_ws_index[d.first_det_id]
    assert stop_ws_index == d.detector_id_to_ws_index[
        d.first_det_id+d.n_pixels_per_tube]
    ws_indices_range = d.get_current_ws_indices_range()
    assert len(ws_indices_range) == d.n_pixels_per_tube
    assert start_ws_index == ws_indices_range[0]
    assert stop_ws_index-1 == ws_indices_range[-1]

    # Now the data for the same tube
    y, e = d.get_ws_data()
    assert len(y) == d.n_pixels_per_tube
    assert len(e) == d.n_pixels_per_tube

    # move to the 2d tube
    d.next_tube()
    start_ws_index, stop_ws_index = d.get_current_ws_indices()

    assert start_ws_index == d.detector_id_to_ws_index[
        d.first_det_id+d.n_pixels_per_tube]
    y2, _ = d.get_ws_data()
    assert np.array_equal(y, y)
    assert not np.array_equal(y, y2)

    # move to the 3rd tube
    d.next_tube()
    # Mask 3rd tube
    MaskBTP(Workspace=dark_current_ws, Tube="3")
    pixels_masked = d.get_pixels_masked()
    # All pixels should be masked
    assert np.count_nonzero(pixels_masked) == 256

    # move to the 4th tube
    d.next_tube()
    pixels_masked = d.get_pixels_masked()
    # None of pixels should be masked
    assert np.count_nonzero(pixels_masked) == 0

    # Get the infinite pixels now
    pixels_infinite = d.get_pixels_infinite()
    # None of pixels should be infinite
    assert np.count_nonzero(pixels_infinite) == 0
    # Let's mock infinite pixels. Set a tube to infinite
    start_ws_index, stop_ws_index = d.get_current_ws_indices()
    for ws_idx in range(start_ws_index, stop_ws_index):
        dark_current_ws.setY(ws_idx, np.array([Property.EMPTY_DBL]))
        dark_current_ws.setE(ws_idx, np.array([1]))
    # Get the infinite pixels now
    pixels_infinite = d.get_pixels_infinite()
    # All of pixels should be infinite
    assert np.count_nonzero(pixels_infinite) == 256
    # Get the coordinates Y of this tube:
    y_coordinates = d.get_y_coordinates()
    assert len(y_coordinates) == d.n_pixels_per_tube

    # # move to the 5th tube
    d.next_tube()
    y_coordinates_2 = d.get_y_coordinates()
    # Arrays must be equal
    assert np.array_equal(y_coordinates, y_coordinates_2)
