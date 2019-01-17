from __future__ import absolute_import, division, print_function

import tempfile

import pytest

from mantid import mtd
from mantid.simpleapi import (CalculateSensitivity, ClearMaskFlag, SANSMaskDTP,
                              LoadEmptyInstrument, LoadHFIRSANS, LoadMask,
                              MaskDetectors, MoveInstrumentComponent,
                              ReplaceSpecialValues, SANSSolidAngle, SaveNexus)
from ornl.sans.hfir.biosans.beam_finder import direct_beam_center
from ornl.sans.hfir.dark_current import subtract_normalized_dark
from ornl.sans.hfir.normalisation import time
from ornl.sans.transmission import (apply_transmission,
                                    calculate_radius_from_input_ws,
                                    zero_angle_transmission)
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
    flood_dc_corrected_ws = subtract_normalized_dark(
        flood_ws, dark_current_norm_ws)
    empty_transmission_dc_corrected_ws = subtract_normalized_dark(
        empty_transmission_ws, dark_current_norm_ws)
    flood_beamcenter_dc_corrected_ws = subtract_normalized_dark(
        flood_beamcenter_ws, dark_current_norm_ws)
    ############################################################################
    # Normalization (In the original script they use time normalization)
    flood_dc_time_corrected_ws = time(flood_dc_corrected_ws)
    empty_transmission_time_corrected_ws = time(
        empty_transmission_dc_corrected_ws)
    flood_beamcenter_time_corrected_ws = time(flood_beamcenter_dc_corrected_ws)

    ############################################################################
    # Find the beam center
    # First mask de wing detector (number 2)
    x, y, y_gravity = direct_beam_center(
        flood_beamcenter_time_corrected_ws, tolerance=0.00125, direct_beam=True, 
        beam_radius=0.0155)
    print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))
    # Since the Loader translates the X the values are almost the same
    # Test if it's right!
    assert x == pytest.approx(0.011, abs=1e-3)
    assert y == pytest.approx(-0.012, abs=1e-3)
    # Get the right geometry
    MoveInstrumentComponent(
        Workspace=flood_ws, ComponentName='detector1', X=-x, Y=-y)
    ############################################################################    
    # Solid Angle correction
    solid_angle_ws = SANSSolidAngle(
        InputWorkspace=flood_dc_corrected_ws, Type='Tube')
    flood_dc_time_sa_corrected_ws = flood_dc_time_corrected_ws / solid_angle_ws
    # Do we need solid angle? Solid angle is the same (?)
    empty_transmission_time_sa_corrected_ws =\
        empty_transmission_time_corrected_ws / solid_angle_ws    
    ############################################################################    
    # Transmission correction
    radius = calculate_radius_from_input_ws(
        empty_transmission_time_sa_corrected_ws)
    # This returns a WS with the sensitivity value + error
    calculated_transmission_ws = zero_angle_transmission(
        flood_dc_time_sa_corrected_ws, empty_transmission_time_sa_corrected_ws,
        radius, "x?????", delete_temp_wss=True).transmission
    flood_dc_time_sa_trans_corrected_name = \
        'flood_dc_time_sa_trans_corrected_ws'
    apply_transmission(flood_dc_time_sa_corrected_ws,
                       flood_dc_time_sa_trans_corrected_name,
                       trans_ws=calculated_transmission_ws,
                       theta_dependent=False)
    flood_all_corrected = mtd[flood_dc_time_sa_trans_corrected_name]
    flood_all_corrected = ReplaceSpecialValues(
        InputWorkspace=flood_all_corrected, NaNValue=0, InfinityValue=0)

    ############################################################################    
    # Sensitivity calculation
    sensitivity_ws = CalculateSensitivity(InputWorkspace=flood_all_corrected,
                                          MinSensitivity=0.1,
                                          MaxSensitivity=1.7)
    ############################################################################    
    # Load and mask sensitivity according to the beamstop
    MaskDetectors(Workspace=sensitivity_ws, MaskedWorkspace=flood_mask_ws)
    ############################################################################    
    # Let's interpolate the masks
    #interpolate_mask(sensitivity_ws, flood_mask_ws, polynomial_degree=2)

    with tempfile.NamedTemporaryFile(
            delete=False, prefix="sensitivivity_biosans_",
            suffix=".nxs") as temp_file:
        SaveNexus(InputWorkspace=sensitivity_ws,
                  Filename=temp_file.name,
                  Title='CG2 exp206 sensitivivity')
        print("Saved NeXus sensitivity file to: {}".format(temp_file.name))
