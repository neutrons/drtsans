from __future__ import (absolute_import, division, print_function)

import pytest
import numpy as np
from numpy.testing import assert_allclose
import re

from mantid.kernel import ConfigService
from mantid.simpleapi import (Load, EQSANSLoad, SumSpectra, RebinToWorkspace,
                              MoveInstrumentComponent, ConvertUnits,
                              CloneWorkspace)
from mantid.api import AnalysisDataService
from ornl.settings import amend_config
from ornl.sans.sns.eqsans.correct_frame import correct_frame


def compare_to_eqsans_load(alg_out, wo, dl, s2d, ltc, htc):
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
    try:
        sm = SumSpectra(eq_out.OutputWorkspace)
        wave_info = eq_out.OutputMessage.split('Wavelength range:')[-1]
        lambdas = [float(x) for x in re.findall(r'\d+\.\d+', wave_info)]
        # compare wavelengths of the lead pulse
        assert_allclose((alg_out.lead_band.min, alg_out.lead_band.max),
                        lambdas[0:2], rtol=0.15)
        # compare wavelengths of the skipped pulse
        if alg_out.skip_band is not None:
            assert_allclose((alg_out.skip_band.min, alg_out.skip_band.max),
                            lambdas[2:4], rtol=0.15)
        # compare intensities summed up over all detectors
        ws_l = ConvertUnits(alg_out.ws, Target='Wavelength', Emode='Elastic')
        ws_l = RebinToWorkspace(ws_l, sm, PreserveEvents=False)
        ws_l = SumSpectra(ws_l)

        non_zero = np.where(sm.dataY(0) > 0.0)
        a = sm.dataY(0)[non_zero]
        b = ws_l.dataY(0)[non_zero]
        assert abs(1.0 - np.mean(b / a)) < 0.05
    finally:
        AnalysisDataService.remove(sm.name())
        AnalysisDataService.remove(ws_l.name())
        AnalysisDataService.remove(eq_out.OutputWorkspace.name())


def test_correct_frame():
    trials = dict(  # configurations with frame skipped
                  skip_1m=('EQSANS_86217', 200, 1000, 0.02, 1.3),
                  skip_4m=('EQSANS_92353', 200, 200, 0.02, 4.0),
                  skip_5m=('EQSANS_85550', 200, 1500, 0.02, 5.0),
                  # configurations with no frame skipped
                  nonskip_1m=('EQSANS_101595', 200, 1000, 0.02, 1.3),
                  nonskip_4m=('EQSANS_88565', 200, 1000, 0.02, 4.0),
                  nonskip_8m=('EQSANS_88901', 200, 1500, 0.02, 8.0))

    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        for run_number, low_tof_cut, high_tof_cut, wavelength_bin,\
                source_to_detector_distance in trials.values():
            wo = Load(Filename=run_number, OutputWorkspace='original')
            ws = CloneWorkspace(wo, OutputWorkspace='copy')
            MoveInstrumentComponent(ws, ComponentName='detector1',
                                    Z=source_to_detector_distance)
            alg_out = correct_frame(ws, low_tof_cut, high_tof_cut)
            compare_to_eqsans_load(alg_out, wo, wavelength_bin,
                                   source_to_detector_distance,
                                   low_tof_cut, high_tof_cut)
            AnalysisDataService.remove(alg_out.ws.name())
            AnalysisDataService.remove(wo.name())


if __name__ == '__main__':
    pytest.main()
