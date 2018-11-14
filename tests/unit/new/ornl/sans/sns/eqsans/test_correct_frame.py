from __future__ import (absolute_import, division, print_function)

import pytest
import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose
import re

from mantid.kernel import ConfigService
from mantid.simpleapi import (Load, EQSANSLoad, SumSpectra, RebinToWorkspace,
                              MoveInstrumentComponent, ConvertUnits,
                              CloneWorkspace, SaveNexus)
from mantid.api import AnalysisDataService
from ornl.sans.sns.eqsans.correct_frame import correct_frame#, replace_tofs


def compare_to_eqsans_load(alg_out, wo, dl, s2d, ltc, htc, run_number):
    if run_number == '92144':
        print('hello')
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
                        lambdas[0:2], rtol=0.1)
        #print((alg_out.lead_band.min, alg_out.lead_band.max), lambdas[0:2])
        # compare wavelengths of the skipped pulse
        if alg_out.skip_band is not None:
            assert_allclose((alg_out.skip_band.min, alg_out.skip_band.max),
                            lambdas[2:4], rtol=0.1)
            #print((alg_out.skip_band.min, alg_out.skip_band.max), lambdas[2:4])
        # compare intensities summed up over all detectors
        ws_l = ConvertUnits(alg_out.ws, Target='Wavelength', Emode='Elastic')
        ws_l = RebinToWorkspace(ws_l, sm,
                                PreserveEvents=False)  # Params=[sm.dataX(0)[0], dl, sm.dataX(0)[-1]])
        ws_l = SumSpectra(ws_l)
        SaveNexus(ws_l, '/tmp/junk_ws.nxs')
        SaveNexus(sm, '/tmp/junk_sm.nxs')
        non_zero = np.where(sm.dataY(0) > 0.0)
        a = sm.dataY(0)[non_zero]
        b = ws_l.dataY(0)[non_zero]
        #print(abs(1.0 - np.mean(b / a)))
        assert abs(1.0 - np.mean(b / a)) < 0.05
    finally:
        AnalysisDataService.remove(sm.name())
        AnalysisDataService.remove(ws_l.name())
        AnalysisDataService.remove(eq_out.OutputWorkspace.name())


def test_correct_frame():
    config = ConfigService.Instance()
    previous_instrument = config['instrumentName']
    config['instrumentName'] = 'EQSANS'
    previous_archive = config['datasearch.searcharchive']
    config['datasearch.searcharchive'] = 'on'

    trials = dict(skip=('92353', 200, 200, 0.02, 4.0),
                  # configurations with no frame skipped
                  porasil_1m=('92164', 200, 1000, 0.02, 1.3),
                  porasil_4m=('92149', 200, 1000, 0.02, 4.0),
                  porasil_8m=('92144', 200, 1500, 0.02, 8.0))
    #trials = dict(porasil_1m=('92164', 200, 1000, 0.02, 1.3))
    try:
        for run_number, low_tof_cut, high_tof_cut, wavelength_bin,\
            source_to_detector_distance in trials.values():
            print('#########\n\nrun_number =', run_number, '\n\n#########')
            wo = Load(Filename=run_number, OutputWorkspace='original')
            ws = CloneWorkspace(wo, OutputWorkspace='copy')
            MoveInstrumentComponent(ws, ComponentName='detector1',
                                    Z=source_to_detector_distance)
            alg_out = correct_frame(ws, low_tof_cut, high_tof_cut)
            compare_to_eqsans_load(alg_out, wo, wavelength_bin,
                                   source_to_detector_distance,
                                   low_tof_cut, high_tof_cut,
                                   run_number)
            AnalysisDataService.remove(alg_out.ws.name())
            AnalysisDataService.remove(wo.name())
    finally:
        config['instrumentName'] = previous_instrument
        config['datasearch.searcharchive'] = previous_archive


if __name__ == '__main__':
    #pytest.main()
    test_correct_frame()