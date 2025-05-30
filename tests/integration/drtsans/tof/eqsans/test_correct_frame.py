import pytest
from pytest import approx
import numpy as np
from os.path import join as path_join

from mantid.simpleapi import (
    Load,
    DeleteWorkspace,
    EQSANSLoad,
    SumSpectra,
    RebinToWorkspace,
    MoveInstrumentComponent,
    ConvertUnits,
    CloneWorkspace,
    LoadNexusMonitors,
    LoadEventNexus,
)
from mantid.api import AnalysisDataService
from mantid.kernel import amend_config
from drtsans.tof.eqsans.correct_frame import (
    correct_detector_frame,
    correct_monitor_frame,
    smash_monitor_spikes,
)
from drtsans.tof.eqsans.load import load_events_monitor
from drtsans.tof.eqsans import correct_frame
from drtsans.geometry import source_detector_distance


# (run, wavelength_bin, source-detector-distance, monitor-min-TOF)
trials = dict(
    # configurations with frame skipped
    skip_1m=("EQSANS_86217.nxs.h5", 0.02, 1.3, -1),
    skip_4m=("EQSANS_92353.nxs.h5", 0.02, 4.0, -1),
    skip_5m=("EQSANS_85550.nxs.h5", 0.02, 5.0, -1),
    # configurations with no frame skipped
    nonskip_1m=("EQSANS_101595.nxs.h5", 0.02, 1.3, 15388),
    nonskip_4m=("EQSANS_88565.nxs.h5", 0.02, 4.0, 25565),
    nonskip_8m=("EQSANS_88901.nxs.h5", 0.02, 8.0, 30680),
)


def compare_to_eqsans_load(ws, wo, dl, s2d, ltc, htc):
    r"""Compare wavelength-dependent intensities  after
    TOF correction with EQSANSCorrectFrame and EQSANSLoad.
    For comparison, sum intensites over all detectors.
    """
    eq_out = EQSANSLoad(
        InputWorkspace=wo.name(),
        OutputWorkspace=AnalysisDataService.unique_hidden_name(),
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
        LoadMonitors=False,
    )
    try:
        sm = SumSpectra(eq_out.OutputWorkspace)
        # wave_info = eq_out.OutputMessage.split('Wavelength range:')[-1]
        # lambdas = [float(x) for x in re.findall(r'\d+\.\d+', wave_info)]
        # compare intensities summed up over all detectors
        ws_l = ConvertUnits(ws, Target="Wavelength", Emode="Elastic")
        ws_l = RebinToWorkspace(ws_l, sm, PreserveEvents=False)
        ws_l = SumSpectra(ws_l)
        non_zero = np.where((sm.dataY(0) > 0) & (ws_l.dataY(0) > 0))[0]
        a = sm.dataY(0)[non_zero]
        b = ws_l.dataY(0)[non_zero]
        assert abs(1.0 - np.mean(b / a)) < 0.1
    finally:
        for w in (eq_out.OutputWorkspace, ws_l, sm):
            if AnalysisDataService.doesExist(w.name()):
                AnalysisDataService.remove(w.name())


@pytest.mark.datarepo
def test_correct_detector_frame(datarepo_dir, temp_workspace_name):
    for v in trials.values():
        run_number, wavelength_bin, sdd = v[0:3]
        wo = LoadEventNexus(path_join(datarepo_dir.eqsans, run_number), OutputWorkspace=temp_workspace_name())
        ws = LoadEventNexus(path_join(datarepo_dir.eqsans, run_number), OutputWorkspace=temp_workspace_name())
        MoveInstrumentComponent(ws, ComponentName="detector1", Z=sdd)
        correct_detector_frame(ws, path_to_pixel=False)
        compare_to_eqsans_load(ws, wo, wavelength_bin, sdd, 500, 2000)


@pytest.mark.datarepo
def test_smash_monitor_spikes(datarepo_dir):
    # Smash two spikes
    w = load_events_monitor("EQSANS_88565.nxs.h5", data_dir=datarepo_dir.eqsans)
    w = smash_monitor_spikes(w)
    assert max(w.dataY(0)) < 1e3
    DeleteWorkspace(w)

    # Monitor data is crap
    w = load_events_monitor("EQSANS_101595.nxs.h5", data_dir=datarepo_dir.eqsans)
    with pytest.raises(RuntimeError, match="Monitor spectrum is flat"):
        smash_monitor_spikes(w)
    DeleteWorkspace(w)


@pytest.mark.datarepo
def test_correct_monitor_frame(datarepo_dir):
    for k, v in trials.items():
        with amend_config(data_dir=datarepo_dir.eqsans):
            w = LoadNexusMonitors(v[0], LoadOnly="Events", OutputWorkspace=AnalysisDataService.unique_hidden_name())
        if not bool(k.find("skip")):  # run in skip frame mode
            with pytest.raises(RuntimeError, match="cannot correct monitor"):
                correct_monitor_frame(w)
        else:
            correct_monitor_frame(w)
            assert w.getSpectrum(0).getTofMin() == approx(v[3], abs=1)


@pytest.mark.datarepo
def test_convert_to_wavelength(datarepo_dir):
    with amend_config(data_dir=datarepo_dir.eqsans):
        for v in trials.values():
            run_number, wavelength_bin, sadd = v[0:3]
            wo = Load(Filename=run_number, OutputWorkspace=AnalysisDataService.unique_hidden_name())
            ws = CloneWorkspace(wo, OutputWorkspace=AnalysisDataService.unique_hidden_name())
            MoveInstrumentComponent(ws, ComponentName="detector1", Z=sadd)
            correct_frame.correct_detector_frame(ws, path_to_pixel=False)
            sodd = source_detector_distance(ws, unit="m")
            bands = correct_frame.transmitted_bands_clipped(ws, sodd, 500, 2000, interior_clip=True)
            ws = correct_frame.convert_to_wavelength(ws, bands, wavelength_bin)
            compare_to_eqsans_load(ws, wo, wavelength_bin, sadd, 500, 2000)


if __name__ == "__main__":
    pytest.main([__file__])
