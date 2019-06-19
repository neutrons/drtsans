from __future__ import absolute_import, division, print_function

import tempfile

import pytest

from mantid import mtd
from mantid.kernel import Property
from mantid.simpleapi import (CalculateEfficiency, ClearMaskFlag,
                              LoadEmptyInstrument, LoadHFIRSANS, LoadMask,
                              MaskDetectors, MoveInstrumentComponent,
                              ReplaceSpecialValues, SANSSolidAngle, SaveNexus)
from ornl.sans.hfir.gpsans.beam_finder import direct_beam_center
from ornl.sans.hfir.normalisation import monitor
from ornl.sans.sensitivity import inf_value_to_mask


'''
For every flood:
    LoadData
    Prepare Flood:
        - Geometry
        - SA
        - DC subtraction
        - Monitor Normalization

    sensitivity = Calculate sensitivity

    Mask beamstop

Join (average) all the sensitivities in one single file
Save file as nexus

'''


@pytest.mark.skip(reason="current CalculateEfficiency algorithm not aplicable")
@pytest.mark.offline
def test_sensitivity_procedural(gpsans_sensitivity_dataset):
    dark_current_ws = LoadHFIRSANS(
        Filename=gpsans_sensitivity_dataset['dark_current'])

    for trans in [0, 200, 400]:
        # Get the file names
        flood_file = gpsans_sensitivity_dataset["flood_trans_{}".format(trans)]
        flood_beamcenter_file = gpsans_sensitivity_dataset[
            "flood_trans_{}_beamcenter".format(trans)]
        flood_mask_file = gpsans_sensitivity_dataset[
            "flood_trans_{}_mask".format(trans)]
        #
        # Load the files into WS
        flood_ws = LoadHFIRSANS(Filename=flood_file)
        flood_beamcenter_ws = LoadHFIRSANS(Filename=flood_beamcenter_file)
        flood_mask_ws = LoadMask(Instrument='CG2', InputFile=flood_mask_file,
                                 RefWorkspace=flood_ws.name())
        #
        # Find the beam center
        x, y = direct_beam_center(flood_beamcenter_ws)
        print("Beam center found = ({:.3}, {:.3}) meters.".format(x, y))

        # Since the Loader translates the X the values are almost the same
        # Test if it's right!
        assert x == pytest.approx(-0.021, abs=1e-3)
        assert y == pytest.approx(-0.03, abs=1e-2)
        #
        # Get the right geometry
        MoveInstrumentComponent(Workspace=flood_ws.name(),
                                ComponentName='detector1', X=-x, Y=-y)
        #
        # DC Subtraction
        flood_dc_corrected_ws = flood_ws - dark_current_ws
        #
        # Solid Angle correction
        solid_angle_ws = SANSSolidAngle(
            InputWorkspace=flood_dc_corrected_ws.name(), Type='Tube')
        flood_dc_sa_corrected_ws = flood_dc_corrected_ws / solid_angle_ws
        #
        # Monitor Normalization
        flood_dc_sa_mon_corrected_ws = monitor(flood_dc_sa_corrected_ws)
        #
        # Sensitivity
        sensitivity_ws_name = "sensitivity_{}".format(trans)
        CalculateEfficiency(
            InputWorkspace=flood_dc_sa_mon_corrected_ws.name(),
            OutputWorkspace=sensitivity_ws_name,
            MinThreshold=0.5, MaxThreshold=1.5)
        #
        # Load and mask sensitivity according to the beamstop
        MaskDetectors(
            Workspace=sensitivity_ws_name, MaskedWorkspace=flood_mask_ws)
        ####################################
        # This is for the next step: the averaging of all sensitivities
        # WSs with 0 in the Mask, 1 elsewhere
        zero_ws_name = "zero_{}".format(trans)
        LoadEmptyInstrument(InstrumentName='cg2', OutputWorkspace=zero_ws_name)
        MaskDetectors(
            Workspace=zero_ws_name, MaskedWorkspace=sensitivity_ws_name)
        ClearMaskFlag(Workspace=zero_ws_name, ComponentName='detector1')
        ClearMaskFlag(Workspace=sensitivity_ws_name, ComponentName='detector1')
        assert mtd[sensitivity_ws_name].readY(31874)[0] == Property.EMPTY_DBL

        inf_value_to_mask(mtd[sensitivity_ws_name])
        assert mtd[sensitivity_ws_name].readY(31874)[0] == 1
        assert mtd[sensitivity_ws_name].readE(31874)[0] == 0
        assert mtd[sensitivity_ws_name].detectorInfo().isMasked(31874)

    # Sum the zeros ws: Max = 3, Min = 2
    zeros_summed_ws = mtd['zero_0'] + mtd['zero_200'] + mtd['zero_400']
    sensitivities_summed_ws = mtd['sensitivity_0'] + \
        mtd['sensitivity_200'] + mtd['sensitivity_400']

    # Divide
    sensitivity_final_ws = sensitivities_summed_ws / zeros_summed_ws
    sensitivity_final_ws = ReplaceSpecialValues(
        InputWorkspace=sensitivity_final_ws, NaNValue=0, InfinityValue=0)

    # this is a bad pixel!
    assert sensitivity_final_ws.detectorInfo().isMasked(31874)

    with tempfile.NamedTemporaryFile(delete=True, prefix="sensitivivity_cg2_",
                                     suffix=".nxs") as temp_file:
        SaveNexus(InputWorkspace=sensitivity_final_ws,
                  Filename=temp_file.name,
                  Title='CG2 exp206 sensitivivity')
        print("Saved NeXus sensitivity file to: {}".format(temp_file.name))
