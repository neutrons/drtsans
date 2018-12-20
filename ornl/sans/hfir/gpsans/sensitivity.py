# -*- coding: utf-8 -*-
# flake8: noqa

from mantid.simpleapi import (
    LoadHFIRSANS, SANSMaskDTP, FindCenterOfMassPosition)


def sensitibity_file(files):

    Load(
    Filename='/home/rhf/Documents/SANS/GPSANS/20181018-SensitivityFiles/CG2_exp130_scan0023_0001.xml\
    ,/home/rhf/Documents/SANS/GPSANS/20181018-SensitivityFiles/CG2_exp130_scan0021_0001.xml,\
    /home/rhf/Documents/SANS/GPSANS/20181018-SensitivityFiles/CG2_exp130_scan0019_0001.xml',
    OutputWorkspace='MultiFiles')

    #
    # Draw the masks by hand
    #


    # WSs with 0 in the Mask, 1 elsewhere
    LoadEmptyInstrument(InstrumentName='cg2', OutputWorkspace='zero_19')
    MaskDetectors(Workspace='zero_19', MaskedWorkspace='CG2_exp130_scan0019_0001')
    ClearMaskFlag(Workspace='zero_19', ComponentName='detector1')

    LoadEmptyInstrument(InstrumentName='cg2', OutputWorkspace='zero_21')
    MaskDetectors(Workspace='zero_21', MaskedWorkspace='CG2_exp130_scan0021_0001')
    ClearMaskFlag(Workspace='zero_21', ComponentName='detector1')

    LoadEmptyInstrument(InstrumentName='cg2', OutputWorkspace='zero_23')
    MaskDetectors(Workspace='zero_23', MaskedWorkspace='CG2_exp130_scan0023_0001')
    ClearMaskFlag(Workspace='zero_23', ComponentName='detector1')

    # Sum the zeros ws: Max = 3, Min = 2
    Plus(LHSWorkspace='zero_19', RHSWorkspace='zero_21', OutputWorkspace='zero')
    Plus(LHSWorkspace='zero', RHSWorkspace='zero_23', OutputWorkspace='zero')

    # Sum the data
    ClearMaskFlag(Workspace='CG2_exp130_scan0019_0001', ComponentName='detector1')
    ClearMaskFlag(Workspace='CG2_exp130_scan0021_0001', ComponentName='detector1')
    ClearMaskFlag(Workspace='CG2_exp130_scan0023_0001', ComponentName='detector1')

    Plus(LHSWorkspace='CG2_exp130_scan0019_0001', RHSWorkspace='CG2_exp130_scan0021_0001', OutputWorkspace='data')
    Plus(LHSWorkspace='data', RHSWorkspace='CG2_exp130_scan0023_0001', OutputWorkspace='data')

    # Divide
    Divide(LHSWorkspace='data', RHSWorkspace='zero', OutputWorkspace='final', WarnOnZeroDivide=False)

    #####
