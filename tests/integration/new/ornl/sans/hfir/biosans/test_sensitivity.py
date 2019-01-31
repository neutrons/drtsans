from __future__ import absolute_import, division, print_function

import tempfile

import pytest

from mantid.simpleapi import (CalculateSensitivity, SANSMaskDTP,
                              LoadHFIRSANS, LoadMask, ReplaceSpecialValues,
                              MaskDetectors, MoveInstrumentComponent,
                              SANSSolidAngle, SaveNexus)
from ornl.sans.hfir.biosans.beam_finder import direct_beam_center
from ornl.sans.hfir.dark_current import subtract_normalised_dark
from ornl.sans.hfir.normalisation import time
from ornl.sans.transmission import calculate_transmission
from ornl.sans.sensitivity import interpolate_mask

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
    ############################################################################
    # DC normalisation
    dark_current_norm_ws = time(dark_current_ws)
    # DC Subtraction
    flood_dc_corrected_ws = subtract_normalised_dark(
        flood_ws, dark_current_norm_ws)
    empty_transmission_dc_corrected_ws = subtract_normalised_dark(
        empty_transmission_ws, dark_current_norm_ws)
    flood_beamcenter_dc_corrected_ws = subtract_normalised_dark(
        flood_beamcenter_ws, dark_current_norm_ws)

    ############################################################################
    # Find the beam center
    x, y, y_gravity = direct_beam_center(
        flood_beamcenter_dc_corrected_ws, tolerance=0.00125, direct_beam=True,
        beam_radius=0.0155)
    assert x == pytest.approx(0.011, abs=1e-3)
    assert y == pytest.approx(-0.012, abs=1e-3)

    ############################################################################
    # Get the right geometry
    MoveInstrumentComponent(
        Workspace=flood_ws, ComponentName='detector1', X=-x, Y=-y)
    MoveInstrumentComponent(
        Workspace=empty_transmission_ws, ComponentName='detector1', X=-x, Y=-y)

    ############################################################################
    # Normalization (In the original script they use time normalization)
    flood_dc_time_corrected_ws = time(flood_dc_corrected_ws)
    empty_transmission_dc_time_corrected_ws = time(
        empty_transmission_dc_corrected_ws)

    ############################################################################
    # Solid Angle correction
    solid_angle_ws = SANSSolidAngle(
        InputWorkspace=flood_dc_time_corrected_ws, Type='Tube')
    flood_dc_time_sa_corrected_ws = flood_dc_time_corrected_ws / solid_angle_ws
    # No need to recalculate solid angle.
    # Probably we don't need it for transmission.
    empty_transmission_time_sa_corrected_ws =\
        empty_transmission_dc_time_corrected_ws / solid_angle_ws
    
    empty_transmission_time_sa_corrected_ws = ReplaceSpecialValues(
        InputWorkspace=empty_transmission_time_sa_corrected_ws,
        NaNValue=0, InfinityValue=0)

    ############################################################################
    # This is only to get transmission from the flat measurement
    # The value will be used in the wing detector
    calculated_transmission_value,  calculated_transmission_error = \
        calculate_transmission(flood_dc_time_sa_corrected_ws,
                               empty_transmission_time_sa_corrected_ws)

    ############################################################################
    # Sensitivity calculation
    sensitivity_ws = CalculateSensitivity(
        InputWorkspace=flood_dc_time_sa_corrected_ws, MinSensitivity=0.1,
        MaxSensitivity=2.0)

    ############################################################################
    # Load and mask sensitivity according to the beamstop
    MaskDetectors(Workspace=sensitivity_ws, MaskedWorkspace=flood_mask_ws)

    ############################################################################
    # Let's interpolate the masks
    sensitivity_interpolated_ws = interpolate_mask(
        sensitivity_ws, polynomial_degree=2)

    ############################################################################
    # Save

    with tempfile.NamedTemporaryFile(
            delete=False, prefix="sensitivivity_biosans_",
            suffix=".nxs") as temp_file:
        SaveNexus(InputWorkspace=sensitivity_interpolated_ws,
                  Filename=temp_file.name,
                  Title='CG2 exp206 sensitivivity')
        print("Saved NeXus sensitivity file to: {}".format(temp_file.name))
