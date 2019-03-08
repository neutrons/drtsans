from __future__ import (absolute_import, division, print_function)

import pytest
import numpy as np

from mantid.simpleapi import (Load, EQSANSLoad, SumSpectra, RebinToWorkspace,
                              MoveInstrumentComponent, ConvertUnits,
                              CloneWorkspace)
from mantid.api import AnalysisDataService
from ornl.settings import amend_config
from ornl.sans.sns.eqsans.correct_frame import correct_detector_frame
from ornl.sans.sns.eqsans import correct_frame as cf
from ornl.sans.geometry import source_detector_distance


trials = dict(  # configurations with frame skipped
    skip_1m=('EQSANS_86217', 0.02, 1.3),
    skip_4m=('EQSANS_92353', 0.02, 4.0),
    skip_5m=('EQSANS_85550', 0.02, 5.0),
    # configurations with no frame skipped
    nonskip_1m=('EQSANS_101595', 0.02, 1.3),
    nonskip_4m=('EQSANS_88565', 0.02, 4.0),
    nonskip_8m=('EQSANS_88901', 0.02, 8.0)
    )


def compare_to_eqsans_load2(ws, wo, dl, s2d, ltc, htc):
    eq_out = EQSANSLoad(InputWorkspace=wo.name(),
                        OutputWorkspace='from_EQSANSLoad',
                        NoBeamCenter=False,
                        UseConfigBeam=False,
                        UseConfigTOFCuts=False,
                        LowTOFCut=ltc,
                        highTOFCut=htc,
                        SkipTOFCorrection=False,
                        WavelengthStep=dl,
                        UseConfigMask=False,
                        UseConfig=False,
                        CorrectForFlightPath=False,
                        PreserveEvents=False,
                        SampleDetectorDistance=s2d * 1e3,
                        SampleDetectorDistanceOffset=0.0,
                        SampleOffset=0.0,
                        DetectorOffset=0.0,
                        LoadMonitors=False)
    sm = SumSpectra(eq_out.OutputWorkspace)
    ws_l = SumSpectra(ws)
    ws_l = RebinToWorkspace(ws_l, sm)
    non_zero = np.where((sm.dataY(0) > 0) & (ws_l.dataY(0) > 0))[0]
    a = sm.dataY(0)[non_zero]
    b = ws_l.dataY(0)[non_zero]
    assert abs(1.0 - np.mean(b / a)) < 0.1
    AnalysisDataService.remove(sm.name())
    AnalysisDataService.remove(ws_l.name())
    AnalysisDataService.remove(eq_out.OutputWorkspace.name())


def compare_to_eqsans_load(ws, wo, dl, s2d):
    eq_out = EQSANSLoad(InputWorkspace=wo.name(),
                        OutputWorkspace='from_EQSANSLoad',
                        NoBeamCenter=False,
                        UseConfigBeam=False,
                        UseConfigTOFCuts=False,
                        LowTOFCut=500,
                        highTOFCut=2000,
                        SkipTOFCorrection=False,
                        WavelengthStep=dl,
                        UseConfigMask=False,
                        UseConfig=False,
                        CorrectForFlightPath=False,
                        PreserveEvents=False,
                        SampleDetectorDistance=s2d * 1e3,
                        SampleDetectorDistanceOffset=0.0,
                        SampleOffset=0.0,
                        DetectorOffset=0.0,
                        LoadMonitors=False)
    try:
        sm = SumSpectra(eq_out.OutputWorkspace)
        # wave_info = eq_out.OutputMessage.split('Wavelength range:')[-1]
        # lambdas = [float(x) for x in re.findall(r'\d+\.\d+', wave_info)]
        # compare intensities summed up over all detectors
        ws_l = ConvertUnits(ws, Target='Wavelength', Emode='Elastic')
        ws_l = RebinToWorkspace(ws_l, sm, PreserveEvents=False)
        ws_l = SumSpectra(ws_l)
        non_zero = np.where((sm.dataY(0) > 0) & (ws_l.dataY(0) > 0))[0]
        a = sm.dataY(0)[non_zero]
        b = ws_l.dataY(0)[non_zero]
        assert abs(1.0 - np.mean(b / a)) < 0.1
    finally:
        AnalysisDataService.remove(sm.name())
        AnalysisDataService.remove(ws_l.name())
        AnalysisDataService.remove(eq_out.OutputWorkspace.name())


def test_correct_detector_frame():
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        for run_number, wavelength_bin, sdd in trials.values():
            wo = Load(Filename=run_number, OutputWorkspace='original')
            ws = CloneWorkspace(wo, OutputWorkspace='copy')
            MoveInstrumentComponent(ws, ComponentName='detector1', Z=sdd)
            correct_detector_frame(ws)
            compare_to_eqsans_load(ws, wo, wavelength_bin, sdd)
            AnalysisDataService.remove(ws.name())
            AnalysisDataService.remove(wo.name())


def test_convert_to_wavelength():
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        for run_number, wavelength_bin, sadd in trials.values():
            wo = Load(Filename=run_number, OutputWorkspace='original')
            ws = CloneWorkspace(wo, OutputWorkspace='copy')
            MoveInstrumentComponent(ws, ComponentName='detector1', Z=sadd)
            cf.correct_detector_frame(ws)
            sodd = source_detector_distance(ws, units='m')
            bands = cf.transmitted_bands_clipped(ws, sodd, 500, 2000,
                                                 interior_clip=True)
            ws = cf.convert_to_wavelength(ws, bands, wavelength_bin, 'copy')
            compare_to_eqsans_load2(ws, wo, wavelength_bin, sadd, 500, 2000)


if __name__ == '__main__':
    pytest.main()
