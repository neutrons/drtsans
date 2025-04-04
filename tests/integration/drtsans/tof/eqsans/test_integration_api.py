"""
Test top-level API
"""

from os.path import join as pj
import pytest
from pytest import approx
import numpy as np
from mantid.dataobjects import EventWorkspace

# https://docs.mantidproject.org/nightly/algorithms/CompareWorkspaces-v1.html
# https://docs.mantidproject.org/nightly/algorithms/LoadNexus-v1.html
# https://docs.mantidproject.org/nightly/algorithms/SumSpectra-v1.html
from mantid.simpleapi import SumSpectra, mtd, LoadNexus, CompareWorkspaces
from mantid.simpleapi import DeleteWorkspace
from mantid.kernel import amend_config

# public API
from drtsans.tof import eqsans
from drtsans import solid_angle_correction

# protected API
from drtsans.settings import namedtuplefy
from drtsans.samplelogs import SampleLogs
from drtsans.geometry import main_detector_name

keys = (
    "run",
    "num_events",
    "sdd",
    "ssd",
    "min_tof",
    "max_tof",
    "skip_frame",
    "w_min",
    "w_max",
    "flux_normalized",
)

values = (
    ("EQSANS_86217", 508339, 1300, 14122, 9595, 59633, True, 2.61, 14.72, 540),
    ("EQSANS_92353", 262291, 4000, 14122, 11288, 61309, True, 2.59, 12.98, 431),
    ("EQSANS_85550", 270022, 5000, 14122, 11914, 61930, True, 2.59, 12.43, 599),
    ("EQSANS_101595", 289989, 1300, 14122, 7657, 24384, False, 2.11, 5.65, 150),
    ("EQSANS_88565", 19362, 4000, 14122, 45486, 62172, False, 10.02, 13.2, 743),
    ("EQSANS_88901", 340431, 8000, 14122, 67202, 83868, False, 11.99, 14.62, 56013),
)

run_sets = [{k: v for k, v in zip(keys, value)} for value in values]


@pytest.fixture(scope="module")
def flux_file(datarepo_dir):
    return pj(datarepo_dir.eqsans, "test_normalization", "beam_profile_flux.txt")


@pytest.fixture(scope="module", params=run_sets, ids=[item["run"] for item in run_sets])
@namedtuplefy
def run_infoset(datarepo_dir, request):
    run_set = request.param
    run = run_set["run"]
    tof_workspace = mtd.unique_hidden_name()
    with amend_config(data_dir=datarepo_dir.eqsans):
        eqsans.load_events(run, output_workspace=tof_workspace)
    wavelength_workspace, bands = eqsans.transform_to_wavelength(
        tof_workspace,
        low_tof_clip=500,
        high_tof_clip=2000,
        output_workspace=mtd.unique_hidden_name(),
    )
    wavelength_workspace = eqsans.set_init_uncertainties(wavelength_workspace)
    return {**run_set, **dict(ws=tof_workspace, wl=wavelength_workspace)}


class TestLoadEvents(object):
    @pytest.mark.datarepo
    def test_loading_file(self, run_infoset):
        tof_workspace = mtd[run_infoset.ws]
        assert isinstance(tof_workspace, EventWorkspace)
        assert tof_workspace.getNumberEvents() == run_infoset.num_events

    @pytest.mark.datarepo
    def test_geometry(self, run_infoset):
        ws = mtd[run_infoset.ws]
        # assert distance of detector1 same as that in detectorZ of the logs
        instrument = ws.getInstrument()
        det = instrument.getComponentByName(main_detector_name(instrument))
        d1 = det.getDistance(instrument.getSample())
        assert run_infoset.sdd == pytest.approx(d1 * 1000, abs=1)
        # Check logs
        sl = SampleLogs(ws)
        assert run_infoset.ssd == approx(sl.single_value("source-sample-distance"), abs=1)
        assert run_infoset.sdd == approx(sl.single_value("sample-detector-distance"), abs=1)

    @pytest.mark.datarepo
    def test_tofs(self, run_infoset):
        ws = mtd[run_infoset.ws]
        assert ws.getTofMin() == pytest.approx(run_infoset.min_tof, abs=1)
        assert ws.getTofMax() == pytest.approx(run_infoset.max_tof, abs=1)
        assert bool(SampleLogs(ws).is_frame_skipping.value) == run_infoset.skip_frame

    @pytest.mark.datarepo
    def test_offsets(self, datarepo_dir):
        with amend_config(data_dir=datarepo_dir.eqsans):
            workspace = eqsans.load_events(
                "EQSANS_86217",
                output_workspace=mtd.unique_hidden_name(),
                detector_offset=42,
                sample_offset=-24,
            )
        sample_logs = SampleLogs(workspace)
        source_sample_distance = sample_logs.single_value("source-sample-distance")
        source_detector_distance = sample_logs.single_value("sample-detector-distance")
        assert 14098.0 == pytest.approx(source_sample_distance, abs=0.1)
        assert 1366.0 == pytest.approx(source_detector_distance, abs=0.1)


@pytest.mark.datarepo
def test_transform_to_wavelength(run_infoset):
    ws, bands = eqsans.transform_to_wavelength(
        run_infoset.ws,
        low_tof_clip=500,
        high_tof_clip=2000,
        output_workspace=mtd.unique_hidden_name(),
    )
    ws = eqsans.set_init_uncertainties(ws)
    sl = SampleLogs(ws)
    assert sl.wavelength_min.value == approx(run_infoset.w_min, abs=0.05)
    assert sl.wavelength_max.value == approx(run_infoset.w_max, abs=0.05)
    # assert zero uncertainty assignment
    for i in range(ws.getNumberHistograms()):
        zci = np.where(ws.readY(i) == 0)[0]  # zero count indices
        np.testing.assert_equal(ws.readE(i)[zci], np.ones(len(zci)))


@pytest.mark.datarepo
def test_normalize_by_flux(run_infoset, flux_file):
    data_workspace = run_infoset.wl
    normalized_data_workspace = eqsans.normalize_by_flux(
        data_workspace,
        flux_file,
        method="proton charge",
        output_workspace=mtd.unique_hidden_name(),
    )
    normalized_data_workspace = SumSpectra(normalized_data_workspace)
    assert 1.0e6 * np.average(normalized_data_workspace.readY(0)) == approx(run_infoset.flux_normalized, abs=1.0)
    # clean up
    DeleteWorkspace(normalized_data_workspace)


@pytest.mark.datarepo
def test_subtract_background(datarepo_dir):
    data_dir = pj(datarepo_dir.eqsans, "test_subtract_background")
    ws = LoadNexus(pj(data_dir, "sample.nxs"), OutputWorkspace=mtd.unique_hidden_name())
    ws_name = ws.name()
    wb = LoadNexus(pj(data_dir, "background.nxs"), OutputWorkspace=mtd.unique_hidden_name())
    ws_wb = eqsans.subtract_background(ws, wb, scale=0.42)
    assert ws_wb.name() == ws_name
    assert max(ws_wb.readY(0)) < 1.0e-09


@pytest.mark.datarepo
def test_prepare_monitors(datarepo_dir):
    data_dirs = [
        datarepo_dir.eqsans,
        pj(datarepo_dir.eqsans, "test_integration_api"),
    ]
    with amend_config(data_dir=data_dirs):
        # Raises for a run in skip frame mode
        with pytest.raises(RuntimeError, match="cannot correct monitor"):
            eqsans.prepare_monitors("EQSANS_92353")
        # Compare for a run in non-skip frame mode
        w = eqsans.prepare_monitors("EQSANS_88565")
        assert w.name() == "EQSANS_88565_monitors"
        sl = SampleLogs(w)
        assert sl.wavelength_min.value == approx(9.9, abs=0.1)
        assert sl.wavelength_max.value == approx(13.6, abs=0.1)
        v = LoadNexus(
            "EQSANS_88565_monitors_wav.nxs",
            OutputWorkspace=mtd.unique_hidden_name(),
        )
        assert CompareWorkspaces(w, v, Tolerance=1e-3, ToleranceRelErr=True).Result is True
    # cleanup
    DeleteWorkspace(w)
    DeleteWorkspace(v)
    # NOTE:
    # somehow EQSANS_92353_monitors:	134.563688 MB
    # is left in ADS
    DeleteWorkspace("EQSANS_92353_monitors")


@pytest.mark.datarepo
def test_solid_angle(run_infoset):
    ws2 = solid_angle_correction(
        run_infoset.ws,
        output_workspace=mtd.unique_hidden_name(),
        detector_type="VerticalTube",
    )
    assert isinstance(ws2, EventWorkspace)
    assert ws2.getNumberEvents() == run_infoset.num_events


@pytest.mark.parametrize(
    "name", ["save_ascii_1D", "save_cansas_xml_1D", "save_cansas_nx", "save_nist_dat", "save_nexus"]
)
def test_api_contents(name):
    assert name in dir(eqsans)


if __name__ == "__main__":
    pytest.main([__file__])
